# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_eukaryome import (_assemble_rrna_url,
                                    _retrieve_data_from_eukaryome,
                                    _make_fasta_str, _make_taxonomy_df,
                                    _process_eukaryome_data,
                                    _extract_zip_euk_seq_file,
                                    _extract_7zip_euk_seq_file
                                    )
from q2_types.feature_data import (TSVTaxonomyFormat,
                                   DNAFASTAFormat,
                                   DNAIterator)
from unittest.mock import patch
import pandas as pd
from pandas import DataFrame
import os
import tempfile
import zipfile
import py7zr


class TestGetEukaryome(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.eukaryome_fasta = self.get_data_path('eukaryome-seqs.fasta')
        self.eukaryome_seqs = DNAFASTAFormat(
                            self.eukaryome_fasta,
                            # self.get_data_path('eukaryome-seqs.fasta'),
                            mode='r')
        self.eukaryome_parsed_fasta = \
            self.get_data_path('eukaryome-parsed-seqs.fasta')
        self.eukaryome_parsed_seqs = DNAFASTAFormat(
                            self.eukaryome_parsed_fasta,
                            mode='r')
        self.eukaryome_parsed_tax = TSVTaxonomyFormat(
                            self.get_data_path('eukaryome-parsed-taxa.tsv'),
                            mode='r')

    def test_assemble_rrna_url_01(self):
        obs_url = _assemble_rrna_url(
                            rrna_gene='SSU',
                            version='1.9.3')
        exp_url = (
                 'https://sisu.ut.ee/wp-content/uploads/sites/643/'
                 'General_EUK_SSU_v1.9.3.zip')
        self.assertEqual(obs_url, exp_url)

    def test_assemble_rrna_url_02(self):
        obs_url = _assemble_rrna_url(
                            rrna_gene='longread',
                            version='1.9.2')
        exp_url = (
                 'https://sisu.ut.ee/wp-content/uploads/sites/643/'
                 'General_EUK_longread_v1.9.2.zip')
        self.assertEqual(obs_url, exp_url)

    def test_extract_zip_euk_seq_file(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            zfn = os.path.join(tmpdirname, 'eukaryome-pseqs.zip')
            with zipfile.ZipFile(zfn, "w", zipfile.ZIP_DEFLATED) as myzip:
                myzip.write(self.eukaryome_parsed_fasta,
                            arcname='eukaryome-pseqs.fasta')

            obs_uzs = _extract_zip_euk_seq_file(
                        tmpdirname, 'eukaryome-pseqs', zfn)
            exp_seqs_dict = {seq.metadata['id'].split(';', 1)[0]: str(seq)
                             for seq in self.eukaryome_parsed_seqs.view(
                                DNAIterator)}
            obs_seqs_dict = {seq.metadata['id']: str(seq)
                             for seq in obs_uzs}
            self.assertEqual(exp_seqs_dict, obs_seqs_dict)

    def test_extract_7zip_euk_seq_file(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            zfn = os.path.join(tmpdirname, 'eukaryome-pseqs.zip')
            with zipfile.ZipFile(zfn, "w", zipfile.ZIP_DEFLATED) as myzip:
                szfn = os.path.join(tmpdirname, 'eukaryome-pseqs.7z')
                with py7zr.SevenZipFile(szfn, 'w') as my7z:
                    my7z.write(self.eukaryome_parsed_fasta,
                               arcname='eukaryome-pseqs.fasta')
                myzip.write(szfn, arcname='eukaryome-pseqs.7z')

            obs_uzs = _extract_7zip_euk_seq_file(
                        tmpdirname, 'eukaryome-pseqs', zfn)
            exp_seqs_dict = {seq.metadata['id'].split(';', 1)[0]: str(seq)
                             for seq in self.eukaryome_parsed_seqs.view(
                                 DNAIterator)}
            obs_seqs_dict = {seq.metadata['id']: str(seq)
                             for seq in obs_uzs}
            self.assertEqual(exp_seqs_dict, obs_seqs_dict)

    def test_make_fasta_str(self):
        seqid = 'Seq01'
        seqstr = 'ATGCCGTA'
        exp_fasta_str = '>Seq01\nATGCCGTA\n'
        obs_fasta_str = _make_fasta_str(seqid, seqstr)
        self.assertEqual(obs_fasta_str, exp_fasta_str)

    def test_make_taxonomy_df(self):
        tax_data = [['ID1', 'p__Rotifera;c__Bdelloidea'],
                    ['ID2', 'p__Rotifera;c__Monogononta']]
        exp_df = pd.DataFrame(tax_data,
                              columns=['Feature ID', 'Taxon'])
        exp_df.set_index('Feature ID', inplace=True)
        tax_dict = {'ID1': 'p__Rotifera;c__Bdelloidea',
                    'ID2': 'p__Rotifera;c__Monogononta'}
        obs_df = _make_taxonomy_df(tax_dict)
        self.assertEqual(True, obs_df.equals(exp_df))

    def test_process_eukaryome_data(self):
        euk_data = self.eukaryome_seqs.view(DNAIterator)
        obs_seqs, obs_tax = _process_eukaryome_data(euk_data)
        # tax
        otdf = obs_tax.view(DataFrame)
        etdf = self.eukaryome_parsed_tax.view(DataFrame)
        self.assertEqual(True, otdf.equals(etdf))
        # seq
        exp_seqs_dict = {seq.metadata['id'].split(';', 1)[0]: str(seq)
                         for seq in self.eukaryome_seqs.view(DNAIterator)}
        obs_seqs_dict = {seq.metadata['id']: str(seq)
                         for seq in obs_seqs.view(DNAIterator)}
        self.assertEqual(exp_seqs_dict, obs_seqs_dict)

    def test_get_eukaryome(self):
        def _makey_fakey(fake_seq, fake_tax):
            return (self.eukaryome_parsed_seqs,
                    self.eukaryome_parsed_tax)
        with patch('rescript.get_eukaryome._retrieve_data_from_eukaryome',
                   new=_makey_fakey):
            res = rescript.actions.get_eukaryome_data(rrna_gene=['SSU'])
            self.assertEqual(
                str(res[0].type),
                'Collection[FeatureData[Sequence]]')
            self.assertEqual(
                str(res[1].type),
                'Collection[FeatureData[Taxonomy]]')

    @patch("rescript.get_eukaryome.urlretrieve")
    def test_run_eukaryome_error(self, mock_urlretrieve):
        expected_message = ('Unable to retrieve the following file '
                            'from Eukaryome:\n'
                            'url\n: Simulated network error')
        mock_urlretrieve.side_effect = Exception("Simulated network error")
        with self.assertRaisesRegex(Exception, expected_message):
            _retrieve_data_from_eukaryome("url", "gene")
