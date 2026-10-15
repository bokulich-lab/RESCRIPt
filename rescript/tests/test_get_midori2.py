# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib
import gzip
import os
import shutil
import tempfile
import zipfile
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

    def test_retrieve_data_from_zip_and_decompression_is_closed(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            mzip_path = os.path.join(tmpdirname, 'midori2-refs.zip')
            with zipfile.ZipFile(mzip_path, mode='w') as zipf:
                fasta_gz = os.path.join(tmpdirname, 'midori2-seqs.fasta.gz')
                taxon_gz = os.path.join(tmpdirname, 'midori2-taxa.taxon.gz')
                for source, destination in [
                        (self.midori2_seqs.path, fasta_gz),
                        (self.midori2_tax.path, taxon_gz)]:
                    with zipf.open(destination, 'w') as zip_stream:
                        with gzip.GzipFile(
                                fileobj=zip_stream, mode="wb") as gz_file:
                            with open(source, 'rb') as original_file:
                                shutil.copyfileobj(original_file, gz_file)

            mzip_path_copy = os.path.join(tmpdirname,
                                          'midori2-refs-copy.zip')

            def copy_download(url, mpc=mzip_path_copy):
                shutil.copyfile(mzip_path, mpc)

            with patch('rescript.get_midori2.urlretrieve',
                       side_effect=copy_download):
                obs_seqs, obs_tax = _retrieve_data_from_midori2(
                                        'midori2-refs.zip')

        obs_seq_dict = {seq.metadata['id']: str(seq)
                        for fid, seq_iter in obs_seqs.items()
                        for seq in seq_iter}

        exp_seq_dict = {seq.metadata['id']: str(seq)
                        for seq in self.midori2_seqs.view(
                            DNAIterator)}

        obs_tax_list = [tax for tax in obs_tax.values()]

        self.assertEqual(obs_seq_dict, exp_seq_dict)
        self.assertTrue(obs_tax_list[0].equals(
                        self.midori2_tax.view(pd.DataFrame)))

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
