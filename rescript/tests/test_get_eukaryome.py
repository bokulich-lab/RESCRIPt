# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_eukaryome import (_assemble_rrna_url,
                                    _retrieve_data_from_eukaryome)
from q2_types.feature_data import (TSVTaxonomyFormat,
                                   DNAFASTAFormat)
from unittest.mock import patch


class TestGetEukaryome(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.eukaryome_seqs = DNAFASTAFormat(
                            self.get_data_path('eukaryome-seqs.fasta'),
                            mode='r')
        self.eukaryome_parsed_seqs = DNAFASTAFormat(
                            self.get_data_path('eukaryome-parsed-seqs.fasta'),
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
