# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
import qiime2
import pandas as pd
import pandas.testing as pdt

from rescript import cross_validate
from ..cross_validate import _evaluate_classifications_stats


import_data = qiime2.Artifact.import_data


class TestPipelines(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        # drop feature C1b because it is missing species level
        self.taxa_series = pd.read_csv(
            self.get_data_path('derep-taxa.tsv'),
            sep='\t', index_col=0).squeeze('columns').drop('C1b')
        self.taxa = import_data('FeatureData[Taxonomy]', self.taxa_series)
        seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))
        self.seqs = import_data(
            'FeatureData[Sequence]', seqs.view(pd.Series).drop('C1b'))

    def test_evaluate_classifications_stats(self):
        # simulate predicted classifications at genus level
        taxa = self.taxa_series.copy().apply(
            lambda x: ';'.join(x.split(';')[:6]))
        taxa = qiime2.Artifact.import_data('FeatureData[Taxonomy]', taxa)
        # first round we just make sure this runs
        _evaluate_classifications_stats([self.taxa], [taxa])
        # now the same but input multiple times to test lists of inputs
        obs = _evaluate_classifications_stats([self.taxa, taxa], [taxa, taxa])
        # now inspect and validate the contents
        # we inspect the second eval results to compare perfect match vs.
        # simulated genus-level classification (when species are expected)
        obs_df = obs.to_dataframe()
        exp = pd.DataFrame({
            'Level': {
                '1': 1.0, '2': 2.0, '3': 3.0, '4': 4.0, '5': 5.0, '6': 6.0,
                '7': 7.0, '8': 1.0, '9': 2.0, '10': 3.0, '11': 4.0,
                '12': 5.0, '13': 6.0},
            'Precision': {
                '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0, '6': 1.0,
                '7': 0.0, '8': 1.0, '9': 1.0, '10': 1.0, '11': 1.0,
                '12': 1.0, '13': 1.0},
            'Recall': {
                '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0, '6': 1.0,
                '7': 0.0, '8': 1.0, '9': 1.0, '10': 1.0, '11': 1.0,
                '12': 1.0, '13': 1.0},
            'F-Measure': {
                '1': 1.0, '2': 1.0, '3': 1.0, '4': 1.0, '5': 1.0, '6': 1.0,
                '7': 0.0, '8': 1.0, '9': 1.0, '10': 1.0, '11': 1.0,
                '12': 1.0, '13': 1.0},
            'Dataset': {'1': '1', '2': '1', '3': '1', '4': '1', '5': '1',
                        '6': '1', '7': '1', '8': '2', '9': '2', '10': '2',
                        '11': '2', '12': '2', '13': '2'}}).sort_index()
        pdt.assert_frame_equal(obs_df.sort_index(), exp, check_names=False)

    def test_evaluate_classifications_mismatch_input_count(self):
        with self.assertRaisesRegex(
                ValueError, "Input must contain an equal number"):
            rescript.actions.evaluate_classifications(
                [self.taxa], [self.taxa, self.taxa])

    def test_evaluate_classifications_expected_id_superset_valid(self):
        taxa = qiime2.Artifact.import_data(
            'FeatureData[Taxonomy]', self.taxa.view(pd.Series).drop('A1'))
        rescript.actions.evaluate_classifications([self.taxa], [taxa])
        self.assertTrue(True)

    def test_evaluate_classifications_observed_id_superset_invalid(self):
        taxa = self.taxa.view(pd.Series)
        taxa['new_garbage'] = 'this;is;most;unexpected'
        taxa = qiime2.Artifact.import_data('FeatureData[Taxonomy]', taxa)
        with self.assertRaisesRegex(
                ValueError, "Indices of pair 1 do not match"):
            rescript.actions.evaluate_classifications([self.taxa], [taxa])


class TestTaxaUtilities(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))
        self.taxa = taxa.view(pd.Series)
        seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))
        self.seqs = seqs.view(pd.Series)

    def test_calculate_per_rank_precision_recall(self):
        # trim the reference taxa at different positions to simulate
        # classification results with underclassification
        warped_taxa = pd.Series([';'.join(t.split(';')[:n]) for n, t in zip(
            [7] * 4 + [5] * 8 + [4] * 4, self.taxa.values)])
        exp = pd.DataFrame(
            {'Level': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7},
             'Precision': {
                0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0},
             'Recall': {
                0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 0.75, 5: 0.25, 6: 0.25},
             'F-Measure': {
                0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 0.8571428571428571, 5: 0.4,
                6: 0.4}})
        obs = cross_validate._calculate_per_rank_precision_recall(
            self.taxa, warped_taxa)
        pdt.assert_frame_equal(exp, obs)
