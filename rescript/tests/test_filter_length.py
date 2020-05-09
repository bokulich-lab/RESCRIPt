# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAIterator
from qiime2.plugins import rescript

from rescript.filter_length import _seq_length_within_range


import_data = qiime2.Artifact.import_data


# These tests just checks the plugin action and validation steps; individual
# filtering options and edge cases for length filtering tested below.
class TestFilterByTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))
        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))

    def test_filter_seqs_by_taxon(self):
        labels = ['Bacteria', 'Paenibacillus', 'vaginalis', 's__casei']
        min_lens = [270, 295, 295, 250]
        max_lens = [500, 500, 500, 290]
        global_min = 270
        global_max = 500
        filtered, failed = rescript.actions.filter_seqs_by_taxon(
            self.seqs, self.taxa, labels, min_lens, max_lens, global_min,
            global_max)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {'B1', 'B1b', 'C1', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B2', 'B3', 'B1a',
                          'C1a', 'C1b', 'C1c', 'C1d'}
        self.assertEqual(failed_ids, exp_failed_ids)

    # this test makes sure that empty outputs pass
    def test_filter_seqs_by_taxon_no_seqs_pass_filter(self):
        # all seqs are < 300 min_len
        filtered, failed = rescript.actions.filter_seqs_by_taxon(
            self.seqs, self.taxa, ['Bacteria'], [300])
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = set()
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3',
                          'B1a', 'B1b', 'C1', 'C1a', 'C1b', 'C1c', 'C1d', 'C2'}
        self.assertEqual(failed_ids, exp_failed_ids)

    # this test makes sure that empty outputs pass
    def test_filter_seqs_by_taxon_no_failures(self):
        # all seqs are > 100 min_len
        filtered, failed = rescript.actions.filter_seqs_by_taxon(
            self.seqs, self.taxa, ['Bacteria'], [100])
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {
            'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B1a', 'B1b', 'C1',
            'C1a', 'C1b', 'C1c', 'C1d', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = set()
        self.assertEqual(failed_ids, exp_failed_ids)

    def test_filter_seqs_by_taxon_no_filters_error(self):
        with self.assertRaisesRegex(
                ValueError, "No filters were applied.*min_lens, max_lens."):
            rescript.actions.filter_seqs_by_taxon(
                self.seqs, self.taxa, ['Bacteria'])

    def test_filter_seqs_by_taxon_index_mismatch_error(self):
        missing_taxa = import_data(
            'FeatureData[Taxonomy]',
            self.taxa.view(pd.Series).drop(['C1', 'C2']))
        with self.assertRaisesRegex(
                ValueError, "sequences are missing.*C1, C2"):
            rescript.actions.filter_seqs_by_taxon(
                self.seqs, missing_taxa, ['Bacteria'], [1200])

    def test_filter_seqs_by_taxon_min_lens_mismatch(self):
        with self.assertRaisesRegex(
                ValueError, "labels and min_lens must contain"):
            rescript.actions.filter_seqs_by_taxon(
                self.seqs, self.taxa, ['Bacteria'], [1200, 100])

    def test_filter_seqs_by_taxon_max_lens_mismatch(self):
        with self.assertRaisesRegex(
                ValueError, "labels and max_lens must contain"):
            rescript.actions.filter_seqs_by_taxon(
                self.seqs, self.taxa, ['Bacteria', 'Archaea'], None, [300])


class TestSeqLengthWithinRange(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        # fake seq; only length matters here
        self.seq = 'A' * 265
        # List of taxa filter search terms that match sequence's taxonomy.
        # This will always max keys of min/max_lens, as tested above, so we
        # just make that assumption in these tests.
        self.taxahits = ['Bacteria', 'Paenibacillus']

    def test_seq_length_within_range_min_len_false(self):
        # taxid mins filtering criteria fail
        self.assertFalse(
            _seq_length_within_range(
                self.seq, self.taxahits,
                {'Bacteria': 270, 'Paenibacillus': 260}, None, None, None))

    def test_seq_length_within_range_min_len_true(self):
        # taxid mins filtering criteria pass
        self.assertTrue(
            _seq_length_within_range(
                self.seq, self.taxahits,
                {'Bacteria': 260, 'Paenibacillus': 260}, None, None, None))

    def test_seq_length_within_range_max_len_false(self):
        # taxid maxs filtering criteria fail
        self.assertFalse(
            _seq_length_within_range(
                self.seq, self.taxahits, None,
                {'Bacteria': 270, 'Paenibacillus': 260}, None, None))

    def test_seq_length_within_range_max_len_true(self):
        # taxid maxs filtering criteria pass
        self.assertTrue(
            _seq_length_within_range(
                self.seq, self.taxahits, None,
                {'Bacteria': 270, 'Paenibacillus': 270}, None, None))

    def test_seq_length_within_range_hypothetical_true_no_filters(self):
        # never get here, no limits
        self.assertTrue(
            _seq_length_within_range(
                self.seq, self.taxahits, None, None, None, None))

    def test_seq_length_within_range_global_max_true(self):
        # global_max pass
        self.assertTrue(
            _seq_length_within_range(
                self.seq, self.taxahits, None, None, None, 270))

    def test_seq_length_within_range_global_max_false(self):
        # global_max fail
        self.assertFalse(
            _seq_length_within_range(
                self.seq, self.taxahits, None, None, None, 260))

    def test_seq_length_within_range_global_min_true(self):
        # global_max pass
        self.assertTrue(
            _seq_length_within_range(
                self.seq, self.taxahits, None, None, 260, None))

    def test_seq_length_within_range_global_min_false(self):
        # global_max fail
        self.assertFalse(
            _seq_length_within_range(
                self.seq, self.taxahits, None, None, 270, None))


# This method is just a vsearch wrapper with basic validation, so save on tests
class TestFilterGlobally(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        # Lengths: 12 X 291 nt; 4 X 264 nt
        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))

    def test_filter_seqs_by_taxon_no_filters_error(self):
        with self.assertRaisesRegex(
                ValueError, "No filters were applied.*global_min, global_max"):
            rescript.actions.filter_seqs_globally(self.seqs)

    def test_filter_seqs_globally_by_min_length(self):
        # filter out seqs < 270 nt (N = 4)
        filtered, failed = rescript.actions.filter_seqs_globally(
            self.seqs, global_min=270, global_max=None)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3',
                            'B1a', 'B1b', 'C1', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'C1a', 'C1b', 'C1c', 'C1d'}
        self.assertEqual(failed_ids, exp_failed_ids)

    def test_filter_seqs_globally_by_max_length(self):
        # filter out seqs > 280 nt (N = 12)
        filtered, failed = rescript.actions.filter_seqs_globally(
            self.seqs, global_min=None, global_max=280)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {'C1a', 'C1b', 'C1c', 'C1d'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3',
                          'B1a', 'B1b', 'C1', 'C2'}
        self.assertEqual(failed_ids, exp_failed_ids)

    # this test makes sure that empty outputs pass
    def test_filter_seqs_globally_all_filtered_out(self):
        # all seqs are < 270 or > 280
        filtered, failed = rescript.actions.filter_seqs_globally(
            self.seqs, global_min=270, global_max=280)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = set()
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3',
                          'B1a', 'B1b', 'C1', 'C1a', 'C1b', 'C1c', 'C1d', 'C2'}
        self.assertEqual(failed_ids, exp_failed_ids)

    # this test makes sure that empty outputs pass
    def test_filter_seqs_globally_no_failures(self):
        # all seqs are > 100 min_len
        filtered, failed = rescript.actions.filter_seqs_globally(
            self.seqs, global_min=100, global_max=300)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {
            'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B1a', 'B1b', 'C1',
            'C1a', 'C1b', 'C1c', 'C1d', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = set()
        self.assertEqual(failed_ids, exp_failed_ids)
