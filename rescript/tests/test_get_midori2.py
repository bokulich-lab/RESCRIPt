# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_midori2 import (_assemble_midori2_urls,
                                  _retrieve_data_from_midori2)
from q2_types.feature_data import (TaxonomyFormat,
                                   DNAFASTAFormat,
                                   DNAIterator)
import pandas as pd
from unittest.mock import patch
# from urllib.request import urlopen
# from urllib.error import HTTPError


class TestGetMidori2(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.midori2_tax = TaxonomyFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/midori2-taxa.tsv'),
                            mode='r')
        self.midori2_seqs = DNAFASTAFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/midori2-seqs.fasta'),
                            mode='r')

    # test that appropriate URLs are assembled
    def test_assemble_midori2_urls_01(self):
        fasta_url, tax_url = _assemble_midori2_urls(
                                mito_gene='CO1',
                                version='GenBank263_2024-10-13',
                                ref_seq_type='longest',
                                unspecified_species=True,)
        exp_fasta_url = (
                 'https://www.reference-midori.info/download/'
                 'Databases/GenBank263_2024-10-13/QIIME_sp/longest/MIDORI2_'
                 'LONGEST_NUC_SP_GB263_CO1_QIIME.fasta.gz')
        exp_tax_url = (
                 'https://www.reference-midori.info/download/'
                 'Databases/GenBank263_2024-10-13/QIIME_sp/longest/MIDORI2_'
                 'LONGEST_NUC_SP_GB263_CO1_QIIME.taxon.gz')
        self.assertEqual(fasta_url, exp_fasta_url)
        self.assertEqual(tax_url, exp_tax_url)

    def test_assemble_midori2_urls_02(self):
        fasta_url, tax_url = _assemble_midori2_urls(
                                mito_gene='ND4L',
                                version='GenBank260_2024-04-15',
                                ref_seq_type='uniq',
                                unspecified_species=False,)
        exp_fasta_url = (
                 'https://www.reference-midori.info/download/'
                 'Databases/GenBank260_2024-04-15/QIIME/uniq/MIDORI2_'
                 'UNIQ_NUC_GB260_ND4L_QIIME.fasta.gz')
        exp_tax_url = (
                 'https://www.reference-midori.info/download/'
                 'Databases/GenBank260_2024-04-15/QIIME/uniq/MIDORI2_'
                 'UNIQ_NUC_GB260_ND4L_QIIME.taxon.gz')
        self.assertEqual(fasta_url, exp_fasta_url)
        self.assertEqual(tax_url, exp_tax_url)

    @patch("rescript.get_midori2.urlretrieve")
    def test_run_amrfinder_u_error(self, mock_urlretrieve):
        expected_message = ("Unable to retrieve the following "
                            "file from MIDORI2:\n"
                            "url_seq\n: Simulated network error")
        mock_urlretrieve.side_effect = Exception("Simulated network error")
        with self.assertRaisesRegex(Exception, expected_message):
            _retrieve_data_from_midori2("url_seq", "url_tax")

    def test_get_midori2(self):
        def _makey_fakey(fake_seq, fake_tax):
            return (self.midori2_seqs.view(DNAIterator),
                    self.midori2_tax.view(pd.DataFrame))

        with patch('rescript.get_midori2._retrieve_data_from_midori2',
                   new=_makey_fakey):
            res = rescript.actions.get_midori2_data(mito_gene=['ND4L'])
            self.assertEqual(
                str(res[0].type),
                'Collection[FeatureData[Sequence]]')
            self.assertEqual(
                str(res[1].type),
                'Collection[FeatureData[Taxonomy]]')
