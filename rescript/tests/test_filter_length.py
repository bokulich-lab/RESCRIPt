# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pandas as pd
import qiime2
import pandas.util.testing as pdt
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAIterator
from qiime2.plugins import rescript

from rescript.filter_length import _seq_length_within_range


import_data = qiime2.Artifact.import_data


# These tests mostly check the plugin action and validation steps; individual
# filtering options and edge cases for length filtering tested below.
class TestFilterByTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))
        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))

    # test that nested taxonomic filtering works.
    # The test seqs contain a few different taxa and some seqs of different
    # lengths; this filters based on multiple criteria, and ensures that nested
    # filtering works when more stringent filters are applied for genus/species
    # level than at kingdom level.
    def test_filter_seqs_length_by_taxon_nested(self):
        # if a sequence matches multiple taxonomic terms in the search, we
        # grab the most stringent: longest minimum length and/or shortest
        # maximum length for filtering.
        labels = ['Bacteria', 'Paenibacillus', 'vaginalis', 's__casei']
        min_lens = [270, 295, 295, 250]
        max_lens = [500, 500, 500, 290]
        global_min = 270
        global_max = 500
        filtered, failed = rescript.actions.filter_seqs_length_by_taxon(
            sequences=self.seqs, taxonomy=self.taxa, labels=labels,
            min_lens=min_lens, max_lens=max_lens, global_min=global_min,
            global_max=global_max)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {'B1', 'B1b', 'C1', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B2', 'B3', 'B1a',
                          'C1a', 'C1b', 'C1c', 'C1d'}
        self.assertEqual(failed_ids, exp_failed_ids)

    # this test also tests that nested filtering works, but covers two other
    # cases:
    # 1. keep any seqs that match NONE of the search terms. This case is only
    #    looking for Lactobacillaceae, so all of the Paenibacillaceae should be
    #    ignored and retained without checking length (no global filters)
    # 2. ensure that the most stringent filter is applied when more specific
    #    labels are searched: stringent family overrides lenient species. In
    #    this case, s__acidilacti is an extraneous label (less stringent filter
    #    than the genus-level filter), but s__damnosus can't get away!
    def test_filter_seqs_by_taxon_nested_keep_taxa_without_label_hit(self):
        labels = ['f__Lactobacillaceae', 's__acidilacti', 's__pseudocasei']
        # all seqs are len=291 nt, except for s__acidilacti and an unknown
        # f__Lactobacillaceae that will be filtered here; the more lenient
        # s__acidilacti filter is ignored. The more stringent s__pseudocasei
        # max_len filter, however, gets applied (as already tested above).
        min_lens = [270, 260, 1]
        max_lens = [500, 500, 280]
        filtered, failed = rescript.actions.filter_seqs_length_by_taxon(
            sequences=self.seqs, taxonomy=self.taxa, labels=labels,
            min_lens=min_lens, max_lens=max_lens, global_min=None,
            global_max=None)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B1a', 'B2',
                            'B3', 'C1', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'B1b', 'C1a', 'C1b', 'C1c', 'C1d'}
        self.assertEqual(failed_ids, exp_failed_ids)

    # this test makes sure that empty outputs pass
    def test_filter_seqs_length_by_taxon_no_seqs_pass_filter(self):
        # all seqs are < 300 min_len
        filtered, failed = rescript.actions.filter_seqs_length_by_taxon(
            sequences=self.seqs, taxonomy=self.taxa, labels=['Bacteria'],
            min_lens=[300])
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = set()
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3',
                          'B1a', 'B1b', 'C1', 'C1a', 'C1b', 'C1c', 'C1d', 'C2'}
        self.assertEqual(failed_ids, exp_failed_ids)

    # this test makes sure that empty outputs pass
    def test_filter_seqs_length_by_taxon_no_failures(self):
        # all seqs are > 100 min_len
        filtered, failed = rescript.actions.filter_seqs_length_by_taxon(
            sequences=self.seqs, taxonomy=self.taxa, labels=['Bacteria'],
            min_lens=[100])
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {
            'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3', 'B1a', 'B1b', 'C1',
            'C1a', 'C1b', 'C1c', 'C1d', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = set()
        self.assertEqual(failed_ids, exp_failed_ids)

    def test_filter_seqs_length_by_taxon_no_filters_error(self):
        with self.assertRaisesRegex(
                ValueError, "No filters were applied.*min_lens, max_lens."):
            rescript.actions.filter_seqs_length_by_taxon(
                sequences=self.seqs, taxonomy=self.taxa, labels=['Bacteria'])

    def test_filter_seqs_length_by_taxon_index_mismatch_error(self):
        missing_taxa = import_data(
            'FeatureData[Taxonomy]',
            self.taxa.view(pd.Series).drop(['C1', 'C2']))
        with self.assertRaisesRegex(ValueError, "IDs are missing.*C1, C2"):
            rescript.actions.filter_seqs_length_by_taxon(
                sequences=self.seqs, taxonomy=missing_taxa,
                labels=['Bacteria'], min_lens=[1200])

    def test_filter_seqs_length_by_taxon_min_lens_mismatch(self):
        with self.assertRaisesRegex(
                ValueError, "labels and min_lens must contain"):
            rescript.actions.filter_seqs_length_by_taxon(
                sequences=self.seqs, taxonomy=self.taxa, labels=['Bacteria'],
                min_lens=[1200, 100])

    def test_filter_seqs_length_by_taxon_max_lens_mismatch(self):
        with self.assertRaisesRegex(
                ValueError, "labels and max_lens must contain"):
            rescript.actions.filter_seqs_length_by_taxon(
                sequences=self.seqs, taxonomy=self.taxa,
                labels=['Bacteria', 'Archaea'], min_lens=None, max_lens=[300])


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
                sequence=self.seq, taxahits=self.taxahits,
                mins={'Bacteria': 270, 'Paenibacillus': 260}, maxs=None,
                global_min=None, global_max=None))

    def test_seq_length_within_range_min_len_true(self):
        # taxid mins filtering criteria pass
        self.assertTrue(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits,
                mins={'Bacteria': 260, 'Paenibacillus': 260}, maxs=None,
                global_min=None, global_max=None))

    def test_seq_length_within_range_max_len_false(self):
        # taxid maxs filtering criteria fail
        self.assertFalse(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits,
                mins=None, maxs={'Bacteria': 270, 'Paenibacillus': 260},
                global_min=None, global_max=None))

    def test_seq_length_within_range_max_len_true(self):
        # taxid maxs filtering criteria pass
        self.assertTrue(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits,
                mins=None, maxs={'Bacteria': 270, 'Paenibacillus': 270},
                global_min=None, global_max=None))

    def test_seq_length_within_range_hypothetical_true_no_filters(self):
        # never get here, no limits
        self.assertTrue(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits, mins=None,
                maxs=None, global_min=None, global_max=None))

    def test_seq_length_within_range_global_max_true(self):
        # global_max pass
        self.assertTrue(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits, mins=None,
                maxs=None, global_min=None, global_max=270))

    def test_seq_length_within_range_global_max_false(self):
        # global_max fail
        self.assertFalse(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits, mins=None,
                maxs=None, global_min=None, global_max=260))

    def test_seq_length_within_range_global_min_true(self):
        # global_max pass
        self.assertTrue(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits, mins=None,
                maxs=None, global_min=260, global_max=None))

    def test_seq_length_within_range_global_min_false(self):
        # global_max fail
        self.assertFalse(
            _seq_length_within_range(
                sequence=self.seq, taxahits=self.taxahits, mins=None,
                maxs=None, global_min=270, global_max=None))


# This method is just a vsearch wrapper with basic validation, so save on tests
class TestFilterGlobally(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        # Lengths: 12 X 291 nt; 4 X 264 nt
        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))

    def test_filter_seqs_length_by_taxon_no_filters_error(self):
        with self.assertRaisesRegex(
                ValueError, "No filters were applied.*global_min, global_max"):
            rescript.actions.filter_seqs_length(self.seqs)

    def test_filter_seqs_length_by_min_length(self):
        # filter out seqs < 270 nt (N = 4)
        filtered, failed = rescript.actions.filter_seqs_length(
            self.seqs, global_min=270, global_max=None)
        filtered_ids = {
            seq.metadata['id'] for seq in filtered.view(DNAIterator)}
        exp_filtered_ids = {'A1', 'A2', 'A3', 'A4', 'A5', 'B1', 'B2', 'B3',
                            'B1a', 'B1b', 'C1', 'C2'}
        self.assertEqual(filtered_ids, exp_filtered_ids)
        failed_ids = {seq.metadata['id'] for seq in failed.view(DNAIterator)}
        exp_failed_ids = {'C1a', 'C1b', 'C1c', 'C1d'}
        self.assertEqual(failed_ids, exp_failed_ids)

    def test_filter_seqs_length_by_max_length(self):
        # filter out seqs > 280 nt (N = 12)
        filtered, failed = rescript.actions.filter_seqs_length(
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
    def test_filter_seqs_length_all_filtered_out(self):
        # all seqs are < 270 or > 280
        filtered, failed = rescript.actions.filter_seqs_length(
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
    def test_filter_seqs_length_no_failures(self):
        # all seqs are > 100 min_len
        filtered, failed = rescript.actions.filter_seqs_length(
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


class TestFilterTaxa(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))

    def test_filter_taxa_invalid(self):
        with self.assertRaisesRegex(ValueError, "No filtering criteria"):
            filtered, = rescript.actions.filter_taxa(self.taxa)

    def test_filter_taxa_by_ids(self):
        ids = pd.Index(['A1', 'B1'], name='Feature ID')
        ids_to_keep = qiime2.Metadata(pd.DataFrame(index=ids))
        exp_taxa = pd.Series(
            ['k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; '
             'f__Paenibacillaceae; g__Paenibacillus; s__chondroitinus',
             'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
             'f__Lactobacillaceae; g__Lactobacillus; s__brevis'],
            index=ids, name='Taxon')
        filtered, = rescript.actions.filter_taxa(
            self.taxa, ids_to_keep=ids_to_keep)
        pdt.assert_series_equal(filtered.view(pd.Series), exp_taxa)

    def test_filter_taxa_by_ids_invalid_ids(self):
        ids = pd.DataFrame(
            index=pd.Index(['A1', 'B1', 'D5'], name='Feature ID'))
        ids_to_keep = qiime2.Metadata(ids)
        with self.assertRaisesRegex(ValueError, "IDs are missing.*D5"):
            filtered, = rescript.actions.filter_taxa(
                self.taxa, ids_to_keep=ids_to_keep)

    def test_filter_taxa_by_include(self):
        ids = pd.Index(['C1', 'C2'], name='Feature ID')
        exp_taxa = pd.Series(
            ['k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
             'f__Lactobacillaceae; g__Pediococcus; s__damnosus'] * 2,
            index=ids, name='Taxon')
        filtered, = rescript.actions.filter_taxa(
            self.taxa, include=['damnosus'])
        pdt.assert_series_equal(filtered.view(pd.Series), exp_taxa)

    # note I slip in a little trick here to test order of operations
    # the include statement is run but effectively useless, as exclusion is
    # subsequently done at the order level
    def test_filter_taxa_by_exclude(self):
        ids = pd.Index(['A1', 'A2'], name='Feature ID')
        exp_taxa = pd.Series(
            ['k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; '
             'f__Paenibacillaceae; g__Paenibacillus; s__chondroitinus'] * 2,
            index=ids, name='Taxon')
        filtered, = rescript.actions.filter_taxa(
            self.taxa, include=['brevis', 'Paenibacillus'],
            exclude=['Lactobacillales', 'alvei'])
        pdt.assert_series_equal(filtered.view(pd.Series), exp_taxa)

    # but now look here: we exclude taxa, like above, but explicitly add them
    # back with ids_to_keep
    def test_filter_taxa_by_complex_query(self):
        ids = pd.Index(['A1'], name='Feature ID')
        ids_to_keep = qiime2.Metadata(pd.DataFrame(index=ids))
        exp_taxa = pd.Series(
            ['k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; '
             'f__Paenibacillaceae; g__Paenibacillus; s__chondroitinus',
             'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
             'f__Lactobacillaceae; g__Lactobacillus; s__brevis'],
            index=ids.union({'B1'}), name='Taxon')
        exp_taxa.index.name = 'Feature ID'
        filtered, = rescript.actions.filter_taxa(
            self.taxa, ids_to_keep=ids_to_keep, include=['brevis'],
            exclude=['o__Bacillales'])
        pdt.assert_series_equal(filtered.view(pd.Series), exp_taxa)

    def test_filter_taxa_fail_all_filtered_out(self):
        with self.assertRaisesRegex(ValueError, "All features were filtered"):
            filtered, = rescript.actions.filter_taxa(
                self.taxa, exclude=['Bacteria'])
