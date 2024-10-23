# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import shutil
import tempfile
from unittest.mock import Mock, patch, call, ANY

import pandas as pd
from ncbi.datasets import GenomeApi, ApiClient, ApiException
from q2_types.feature_data import DNAFASTAFormat
from q2_types.genome_data import (LociDirectoryFormat,
                                  ProteinsDirectoryFormat)
from qiime2.plugin.testing import TestPluginBase
from rescript.ncbi import _default_ranks
from rescript.ncbi_datasets import (_get_assembly_descriptors,
                                    _fetch_and_extract_dataset,
                                    _fetch_taxonomy, get_ncbi_genomes)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestNCBIDatasets(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.fake_acc_ids = ['AC_12.1', 'AC_23.2']
        self.fake_assembly_ids = ['GCF_123', 'GCF_234']
        self.fake_tax_ids = ['1234', '2345']

        zip_path = os.path.join(self.temp_dir.name, 'ncbi-dataset')
        self.zipped_data = f"{zip_path}.zip"
        shutil.make_archive(
            zip_path, 'zip', self.get_data_path('ncbi-dataset')
        )

    def generate_fake_response(self, token, genome_count=2, is_error=False):
        if is_error:
            msg = [Mock(error=Mock(message='Something went wrong.'))]
            assemblies = None
        else:
            msg = None
            assemblies = [
                Mock(assembly=Mock(
                    assembly_accession=self.fake_assembly_ids[i],
                    org=Mock(tax_id=self.fake_tax_ids[i])
                )) for i in range(genome_count)
            ]
        return Mock(
            next_page_token=token, assemblies=assemblies, messages=msg
        )

    def test_get_assembly_descriptors_one_page(self):
        with patch.object(GenomeApi, 'assembly_descriptors_by_taxon') as p:
            p.return_value = self.generate_fake_response(None, 2)
            api_instance = GenomeApi(ApiClient())

            obs = _get_assembly_descriptors(
                api_instance, ['complete_chromosome'], 'refseq',
                True, 10, 'some taxon', True
            )

            exp = {
                k: v for k, v
                in zip(self.fake_assembly_ids[:2], self.fake_tax_ids[:2])
            }
            self.assertDictEqual(obs, exp)
            p.assert_called_once_with(
                taxon='some taxon', page_size=10,
                filters_assembly_source='refseq',
                filters_assembly_level=['complete_chromosome'],
                filters_reference_only=True,
                filters_has_annotation=True,
                tax_exact_match=True, page_token=''
            )

    def test_get_assembly_descriptors_many_pages(self):
        with patch.object(GenomeApi, 'assembly_descriptors_by_taxon') as p:
            p.side_effect = [
                self.generate_fake_response('token1', 1),
                self.generate_fake_response('token2', 1),
                self.generate_fake_response(None, 1)
            ]
            api_instance = GenomeApi(ApiClient())

            obs = _get_assembly_descriptors(
                api_instance, ['complete_chromosome'], 'refseq',
                True, 10, 'some taxon', True
            )

            exp = {'GCF_123': '1234', }
            self.assertDictEqual(obs, exp)
            p.assert_has_calls(
                [call(taxon='some taxon', page_size=10,
                      filters_assembly_source='refseq',
                      filters_assembly_level=['complete_chromosome'],
                      filters_reference_only=True,
                      filters_has_annotation=True,
                      tax_exact_match=True, page_token=i)
                 for i in ('', 'token1', 'token2')]
            )

    def test_get_assembly_descriptors_with_api_error(self):
        with patch.object(GenomeApi, 'assembly_descriptors_by_taxon') as p:
            p.side_effect = ApiException('Santa is not real.')
            api_instance = GenomeApi(ApiClient())

            with self.assertRaisesRegex(
                    Exception, r'error while calling NCBI.*Santa is not real.'
            ):
                _get_assembly_descriptors(
                    api_instance, ['complete_chromosome'], 'refseq',
                    True, 10, 'some taxon', True
                )

    def test_get_assembly_descriptors_no_result(self):
        with patch.object(GenomeApi, 'assembly_descriptors_by_taxon') as p:
            p.return_value = self.generate_fake_response(None, is_error=True)
            api_instance = GenomeApi(ApiClient())

            with self.assertRaisesRegex(
                    Exception, ".*Please update your query.*"
            ):
                _get_assembly_descriptors(
                    api_instance, ['complete_chromosome'], 'refseq',
                    True, 10, 'invalid taxon', True
                )

    @patch('tempfile.TemporaryDirectory')
    def test_fetch_and_extract_dataset(self, p):
        test_temp_dir = MockTempDir()
        p.return_value = test_temp_dir

        genomes = DNAFASTAFormat()
        loci = LociDirectoryFormat()
        proteins = ProteinsDirectoryFormat()
        with open(self.zipped_data, 'rb') as fin:
            fake_response = Mock(data=fin.read())

            obs_accessions = \
                _fetch_and_extract_dataset(
                    fake_response, genomes, loci, proteins,
                    only_genomic=True
                )

        exp_accessions = {'GCA_000008865.2': ['BA000007.3'], }
        self.assertDictEqual(obs_accessions, exp_accessions)

        for f in [f'{x}_loci.gff' for x in exp_accessions.keys()]:
            self.assertTrue(os.path.isfile(os.path.join(str(loci), f)))
        for f in [f'{x}_proteins.fasta' for x in exp_accessions.keys()]:
            self.assertTrue(os.path.isfile(os.path.join(str(proteins), f)))

    @patch('tempfile.TemporaryDirectory')
    def test_fetch_and_extract_dataset_all_chromosomes(self, p):
        test_temp_dir = MockTempDir()
        p.return_value = test_temp_dir

        genomes = DNAFASTAFormat()
        loci = LociDirectoryFormat()
        proteins = ProteinsDirectoryFormat()

        with open(self.zipped_data, 'rb') as fin:
            fake_response = Mock(data=fin.read())

            obs_accessions = \
                _fetch_and_extract_dataset(
                    fake_response, genomes, loci, proteins,
                    only_genomic=False
                )

        exp_accessions = {
            'GCA_000008865.2': ['BA000007.3', 'AB011548.2', 'AB011549.2'],
        }
        self.assertDictEqual(obs_accessions, exp_accessions)

        for f in [f'{x}_loci.gff' for x in exp_accessions.keys()]:
            self.assertTrue(os.path.isfile(os.path.join(str(loci), f)))
        for f in [f'{x}_proteins.fasta' for x in exp_accessions.keys()]:
            self.assertTrue(os.path.isfile(os.path.join(str(proteins), f)))

    @patch('rescript.ncbi_datasets.get_taxonomies')
    def test_fetch_taxonomy(self, p):
        p.return_value = ({'GCF_123': 'this;is;some;taxonomy',
                           'GCF_234': 'that;is;another;taxonomy'}, [])

        obs_taxa = _fetch_taxonomy(
            self.fake_assembly_ids,
            self.fake_tax_ids,
            pd.Series(
                {'GCF_123': ['AC_12.1'], 'GCF_234': ['AC_23.2']},
                name="assembly_id"
            ).explode(),
            _default_ranks,
            True
        )

        exp_taxa = pd.DataFrame(
            {'Taxon': ['this;is;some;taxonomy', 'that;is;another;taxonomy']},
            index=self.fake_acc_ids
        )
        exp_taxa.index.name = 'Feature ID'
        pd.testing.assert_frame_equal(obs_taxa, exp_taxa)
        p.assert_called_once_with(
            taxids={'GCF_123': '1234', 'GCF_234': '2345'},
            ranks=_default_ranks, rank_propagation=True,
            logging_level='INFO', n_jobs=2, request_lock=ANY
        )

    @patch('rescript.ncbi_datasets.get_taxonomies')
    def test_fetch_taxonomy_bad_accs(self, p):
        p.return_value = ({}, ['ACC1', 'ACC2'])

        with self.assertRaisesRegex(
            Exception, r'Invalid taxonomy.*\: ACC1, ACC2. Please check.*'
        ):
            _fetch_taxonomy(
                self.fake_assembly_ids,
                self.fake_tax_ids,
                pd.Series(),
                _default_ranks,
                True
            )

    # just test that everything works together
    def test_get_ncbi_genomes(self):
        with patch.object(GenomeApi, 'assembly_descriptors_by_taxon') as p1, \
                patch.object(GenomeApi, 'download_assembly_package') as p2, \
                patch('rescript.ncbi_datasets.get_taxonomies') as p3, \
                open(self.zipped_data, 'rb') as fin:
            p1.return_value = self.generate_fake_response(None, 2)
            fake_response = Mock(data=fin.read())
            p2.return_value = fake_response
            p3.return_value = ({'AC_12.1': 'this;is;some;taxonomy',
                                'AC_23.2': 'that;is;another;taxonomy'}, [])

            result = get_ncbi_genomes(
                taxon='some taxon', assembly_source='refseq',
                assembly_levels=['scaffold', 'contig'],
                only_reference=True, page_size=12, tax_exact_match=True,
                only_genomic=True
            )
            obs_genomes, obs_loci, obs_proteins, obs_taxa = result

            p1.assert_called_once_with(
                taxon='some taxon', page_size=12,
                filters_assembly_source='refseq',
                filters_assembly_level=['scaffold', 'contig'],
                filters_reference_only=True,
                filters_has_annotation=True,
                tax_exact_match=True, page_token=''
            )
            p2.assert_called_once_with(
                self.fake_assembly_ids, exclude_sequence=False,
                include_annotation_type=['PROT_FASTA', 'GENOME_GFF'],
                _preload_content=False
            )
            p3.assert_called_once_with(
                taxids={
                    k: v for k, v
                    in zip(self.fake_assembly_ids, self.fake_tax_ids)
                }, ranks=_default_ranks, rank_propagation=True,
                logging_level='INFO', n_jobs=2, request_lock=ANY
            )
            exp_ids = ('GCA_000008865.2',)
            self.assertIsInstance(obs_genomes, DNAFASTAFormat)
            self.assertIsInstance(obs_loci, LociDirectoryFormat)
            self.assertIsInstance(obs_proteins, ProteinsDirectoryFormat)
            for f in [f'{x}_loci.gff' for x in exp_ids]:
                self.assertTrue(os.path.isfile(os.path.join(str(obs_loci), f)))
            for f in [f'{x}_proteins.fasta' for x in exp_ids]:
                self.assertTrue(
                    os.path.isfile(os.path.join(str(obs_proteins), f)))
