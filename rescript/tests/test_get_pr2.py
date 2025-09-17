# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from qiime2.plugin.testing import TestPluginBase
# from qiime2.plugins import rescript
from rescript.get_pr2 import (_assemble_pr2_urls, _compile_taxonomy_output,
                              _fetch_url)
from q2_types.feature_data import (HeaderlessTSVTaxonomyFormat,
                                   TaxonomyFormat,
                                   DNAFASTAFormat)  # , DNAIterator)
import pandas as pd
from pandas.testing import assert_series_equal
from unittest.mock import patch


class TestGetPR2(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.pr2_tax_df = TaxonomyFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/pr2-tax.tsv'),
                            mode="r").view(pd.DataFrame)
        self.pr2_tax = HeaderlessTSVTaxonomyFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/pr2-tax.tsv'),
                            mode='r')
        self.pr2_tax_format = TaxonomyFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/pr2-tax-formatted.tsv'),
                            mode='r')
        self.pr2_seqs = DNAFASTAFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/pr2-seqs.fasta'),
                            mode='r')

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

    @patch("rescript.get_pr2.urlretrieve")
    def test_run_pr2_error(self, mock_urlretrieve):
        expected_message = ("Unable to retrieve the following "
                            "file from PR2:\n"
                            "url\n: Simulated network error")
        mock_urlretrieve.side_effect = Exception("Simulated network error")
        with self.assertRaisesRegex(Exception, expected_message):
            _fetch_url("url", "path")
