# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib
# import gzip
# import os
# import shutil
# import tempfile
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_midori2 import (_assemble_midori2_url,
                                  _retrieve_data_from_midori2)
from q2_types.feature_data import (TaxonomyFormat,
                                   DNAFASTAFormat,
                                   DNAIterator)
import pandas as pd
from unittest.mock import patch


class TestGetMidori2(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.midori2_tax = \
            TaxonomyFormat(importlib.resources.files(
                'rescript') / 'tests' / 'data' / 'midori2-taxa.tsv', mode='r'
            )
        self.midori2_seqs = \
            DNAFASTAFormat(importlib.resources.files(
                'rescript') / 'tests' / 'data' / 'midori2-seqs.fasta', mode='r'
            )

    # test that appropriate URLs are assembled
    def test_assemble_midori2_url_01(self):
        url = _assemble_midori2_url(
                    version='273',
                    unspecified_species=True)
        exp_url = (
                 'https://zenodo.org/records/22685701/files/'
                 'MIDORI2_GB273_QIIME_sp.zip')
        self.assertEqual(url, exp_url)

    def test_assemble_midori2_url_02(self):
        url = _assemble_midori2_url(
                    version='265',
                    unspecified_species=False)
        exp_url = (
                 'https://zenodo.org/records/23047760/files/'
                 'MIDORI2_GB265_QIIME.zip')
        self.assertEqual(url, exp_url)

    @patch("rescript.get_midori2.urlretrieve")
    def test_run_midori_error(self, mock_urlretrieve):
        expected_message = ("Unable to retrieve the following "
                            "file from MIDORI2:\n"
                            "url_data\n: Simulated network error")
        mock_urlretrieve.side_effect = Exception("Simulated network error")
        with self.assertRaisesRegex(Exception, expected_message):
            _retrieve_data_from_midori2("url_data")

    def test_get_midori2(self):
        def _makey_fakey(fake_data):
            return ({'seq_res': self.midori2_seqs.view(DNAIterator)},
                    {'tax_res': self.midori2_tax.view(pd.DataFrame)})

        with patch('rescript.get_midori2._retrieve_data_from_midori2',
                   new=_makey_fakey):
            res = rescript.actions.get_midori2_data(version='273')
            self.assertEqual(
                str(res[0].type),
                'Collection[FeatureData[Sequence]]')
            self.assertEqual(
                str(res[1].type),
                'Collection[FeatureData[Taxonomy]]')
