# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# from qiime2.plugins import rescript
from qiime2.plugin.testing import TestPluginBase
from rescript.get_pr2 import (_assemble_pr2_urls, _compile_taxonomy_output,
                              _fetch_url, _get_paths, _unzip_fasta,
                              _unzip_taxonomy)
from q2_types.feature_data import (TaxonomyFormat, DNAFASTAFormat,
                                   DNAIterator)
import pandas as pd
from pandas.testing import assert_series_equal
from unittest.mock import patch
import tempfile
import os
import gzip
import shutil


class TestGetPR2(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.pr2_tax_path = self.get_data_path('pr2-tax.tsv')
        self.pr2_tax = TaxonomyFormat(self.pr2_tax_path, mode='r')
        self.pr2_tax_df = self.pr2_tax.view(pd.DataFrame)

        self.pr2_taxf_path = self.get_data_path('pr2-tax-formatted.tsv')
        self.pr2_taxf = TaxonomyFormat(self.pr2_taxf_path, mode='r')
        self.pr2_taxf_df = self.pr2_taxf.view(pd.DataFrame)

        self.pr2_fasta_path = self.get_data_path('pr2-seqs.fasta')
        self.pr2_seqs = DNAFASTAFormat(self.pr2_fasta_path, mode='r')

    # test that appropriate URLs are assembled
    def test_assemble_pr2_urls(self):
        obs_query_urls = _assemble_pr2_urls(version='5.0.0')
        print('obs queries: ', obs_query_urls)

        exp_query_urls = {'fasta': 'https://github.com/pr2database/'
                          'pr2database/releases/download/v5.0.0/'
                          'pr2_version_5.0.0_SSU_mothur.fasta.gz',
                          'tax': 'https://github.com/pr2database/'
                          'pr2database/releases/download/v5.0.0/'
                          'pr2_version_5.0.0_SSU_mothur.tax.gz'}
        print('exp queries: ', exp_query_urls)
        self.assertEqual(obs_query_urls, exp_query_urls)

    def test_compile_taxonomy_output(self):
        obs_tax = _compile_taxonomy_output(self.pr2_tax_df,
                                           ranks=['supergroup',
                                                  'subdivision',
                                                  'genus',
                                                  'species'])
        exp_tax_dict = {
            'AM235536.1.1627_U': 'sgr__Archaeplastida;dvs__Streptophyta_X;'
                                 'g__Uromyrtus;s__Uromyrtus_metrosideros',
            'EU163792.1.1363_U': 'sgr__TSAR;dvs__Ciliophora;'
                                 'g__Trichostomatia_XX;'
                                 's__Trichostomatia_XX_sp.',
            'JQ388892.1.1909_U': 'sgr__Obazoa;dvs__Metazoa;g__Myxobolus;'
                                 's__Myxobolus_musculi'
            }
        exp_tax = pd.Series(exp_tax_dict)
        exp_tax.name = 'Taxon'
        exp_tax.index.name = 'Feature ID'
        assert_series_equal(obs_tax, exp_tax)

    def test_get_paths(self):
        url = ('https://github.com/pr2database/'
               'pr2database/releases/download/v5.0.0/'
               'pr2_version_5.0.0_SSU_mothur.fasta.gz')
        tmpdirname = 'my_tmp_dir'
        obs_inp, obs_outp = _get_paths(url, tmpdirname)

        exp_inp = 'my_tmp_dir/pr2_version_5.0.0_SSU_mothur.fasta.gz'
        exp_outp = 'my_tmp_dir/pr2_version_5.0.0_SSU_mothur.fasta'

        self.assertEqual(obs_inp, exp_inp)
        self.assertEqual(obs_outp, exp_outp)

    def test_unzip_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            pr2_guz = os.path.join(tmpdirname, 'pr2-seqs.fasta')
            pr2_gz = os.path.join(tmpdirname, 'pr2-seqs.fasta.gz')
            with open(self.pr2_fasta_path) as pr2_seq_data:
                with gzip.open(pr2_gz, 'wt') as gzfo:
                    shutil.copyfileobj(pr2_seq_data, gzfo)

                obs_pr2_iter = _unzip_fasta(pr2_gz,
                                            pr2_guz)
                obs_seqs_dict = {seq.metadata['id']: str(seq)
                                 for seq in obs_pr2_iter}

                exp_seqs_dict = {seq.metadata['id']: str(seq)
                                 for seq in self.pr2_seqs.view(
                                     DNAIterator)}
                self.assertEqual(obs_seqs_dict, exp_seqs_dict)

    def test_unzip_taxonomy(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            pr2t_guz = os.path.join(tmpdirname, 'pr2-tax.tsv')
            pr2t_gz = os.path.join(tmpdirname, 'pr2-tax.tsv.gz')
            with open(self.pr2_tax_path) as pr2t_tax_data:
                with gzip.open(pr2t_gz, 'wt') as gzto:
                    shutil.copyfileobj(pr2t_tax_data, gzto)

                obs_tax_df = _unzip_taxonomy(pr2t_gz,
                                             pr2t_guz)

                self.assertEqual(True, obs_tax_df.equals(self.pr2_tax_df))

    @patch("rescript.get_pr2.urlretrieve")
    def test_run_pr2_error(self, mock_urlretrieve):
        expected_message = ("Unable to retrieve the following "
                            "file from PR2:\n"
                            "url\n: Simulated network error")
        mock_urlretrieve.side_effect = Exception("Simulated network error")
        with self.assertRaisesRegex(Exception, expected_message):
            _fetch_url("url", "path")
