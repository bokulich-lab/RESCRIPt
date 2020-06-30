# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
import qiime2
import pandas as pd
import pandas.util.testing as pdt

from rescript import evaluate


import_data = qiime2.Artifact.import_data


class TestEvaluateUtilities(TestPluginBase):
    package = 'rescript.tests'

    # test that warning is raised when there are fewer labels than Taxonomies
    # and that missing labels are labeled numerically
    def test_process_labels_warning(self):
        labels = ['a', 'b']
        # this function just looks at lists, not the actual contents type
        dummy_taxonomies = [1, 2, 3, 4, 5]
        with self.assertWarnsRegex(
                UserWarning, "taxonomies and labels are different lengths"):
            new_labels = evaluate._process_labels(labels, dummy_taxonomies)
        self.assertEqual(new_labels, ['a', 'b', 3, 4, 5])

    # test that _process_labels trims labels when too many are given
    def test_process_labels_too_many_labels(self):
        labels = ['a', 'b', 'c', 'd', 'e']
        dummy_taxonomies = [1, 2, 3]
        new_labels = evaluate._process_labels(labels, dummy_taxonomies)
        self.assertEqual(new_labels, ['a', 'b', 'c'])

    def test_process_labels_empty(self):
        dummy_taxonomies = [1, 2, 3]
        new_labels = evaluate._process_labels(None, dummy_taxonomies)
        self.assertEqual(new_labels, [1, 2, 3])


class TestEvaluateTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))

    # this just tests that the pipeline runs, other tests test proper function
    def test_pipeline(self):
        rescript.actions.evaluate_taxonomy([self.taxa], ["name"], "")

    def test_taxonomic_depth(self):
        obs_depths = evaluate._taxonomic_depth(self.taxa.view(pd.Series), "")
        exp_depths = pd.Series({
            'A1': 7,
            'A2': 7,
            'A3': 7,
            'A4': 7,
            'A5': 7,
            'B1': 7,
            'B2': 7,
            'B3': 7,
            'B1a': 7,
            'B1b': 7,
            'C1': 7,
            'C2': 7,
            'C1a': 7,
            'C1b': 5,
            'C1c': 7,
            'C1d': 7})
        pdt.assert_series_equal(obs_depths.sort_index(),
                                exp_depths.sort_index(), check_names=False)

    def test_summarize_taxonomic_depth(self):
        obs = evaluate.summarize_taxonomic_depth(self.taxa.view(pd.Series), "")
        exp = pd.DataFrame({
            'Number of Features Terminating at Depth': {
                1: 0, 2: 0, 3: 0, 4: 0, 5: 1, 6: 0, 7: 15},
            'Proportion of Features Terminating at Depth': {
                1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0625, 6: 0.0, 7: 0.9375},
            'Number of Features Classified at Depth': {
                1: 16, 2: 16, 3: 16, 4: 16, 5: 16, 6: 15, 7: 15},
            'Proportion of Features Classified at Depth': {
                1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 0.9375, 7: 0.9375},
            'Number of Features Unclassified at Depth': {
                1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 1},
            'Proportion of Features Unclassified at Depth': {
                1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0625, 7: 0.0625}})
        print(obs.to_dict())
        pdt.assert_frame_equal(obs.sort_index(),
                               exp.sort_index(), check_names=False)

    def test_taxonomic_entropy(self):
        obs_ent = evaluate._taxonomic_entropy(self.taxa.view(pd.Series), "", 7)
        exp_ent = pd.DataFrame(
            {'Unique Labels': {
                1: 1.0, 2: 1.0, 3: 1.0, 4: 2.0, 5: 2.0, 6: 4.0, 7: 9.0},
             'Taxonomic Entropy': {
                1: 0.0,
                2: 0.0,
                3: 0.0,
                4: 0.6210863745552451,
                5: 0.6210863745552451,
                6: 1.263740679332812,
                7: 2.1006789212792603}})
        print(obs_ent)
        pdt.assert_frame_equal(obs_ent.sort_index(),
                               exp_ent.sort_index(), check_names=False)


# this test class ensures that _evaluate_taxonomy works with real reference
# database taxonomies. Add to this list as we validate additional reference
# databases / taxonomy styles.
# Currently validated taxonomies:
# greengenes
# SILVA
# UNITE
class TestEvaluateRealLiveTaxonomies(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.taxa = import_data(
            'FeatureData[Taxonomy]',
            self.get_data_path('real-taxa-test.tsv')).view(pd.Series)

    def test_evaluate_taxonomy_real_live_taxonomies_no_rank_handle(self):
        exp = pd.DataFrame({
            'Unique Labels': {
                1: 6.0, 2: 8.0, 3: 8.0, 4: 9.0, 5: 9.0, 6: 9.0, 7: 9.0},
            'Taxonomic Entropy': {1: 1.6769877743224173, 2: 2.0431918705451206,
                                  3: 2.0431918705451206, 4: 2.1972245773362196,
                                  5: 2.1972245773362196, 6: 2.1972245773362196,
                                  7: 2.1972245773362196},
            'Number of Features Terminating at Depth': {
                1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 9},
            'Proportion of Features Terminating at Depth': {
                1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 1.0},
            'Number of Features Classified at Depth': {
                1: 9, 2: 9, 3: 9, 4: 9, 5: 9, 6: 9, 7: 9},
            'Proportion of Features Classified at Depth': {
                1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0},
            'Number of Features Unclassified at Depth': {
                1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0},
            'Proportion of Features Unclassified at Depth': {
                1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0}})

        obs = evaluate._evaluate_taxonomy(self.taxa, rank_handle_regex=None)
        pdt.assert_frame_equal(obs, exp, check_names=False)

    def test_evaluate_taxonomy_real_live_taxonomies_ggsilva_rank_handle(self):
        exp = pd.DataFrame({
            'Unique Labels': {
                1: 6.0, 2: 8.0, 3: 8.0, 4: 9.0, 5: 9.0, 6: 9.0, 7: 9.0},
            'Taxonomic Entropy': {1: 1.6769877743224173, 2: 2.0431918705451206,
                                  3: 2.0431918705451206, 4: 2.1972245773362196,
                                  5: 2.1972245773362196, 6: 2.1972245773362196,
                                  7: 2.1972245773362196},
            'Number of Features Terminating at Depth': {
                1: 1, 2: 0, 3: 0, 4: 0, 5: 2, 6: 0, 7: 6},
            'Proportion of Features Terminating at Depth': {
                1: 0.1111111111111111, 2: 0.0, 3: 0.0, 4: 0.0,
                5: 0.2222222222222222, 6: 0.0, 7: 0.6666666666666666},
            'Number of Features Classified at Depth': {
                1: 9, 2: 8, 3: 8, 4: 8, 5: 8, 6: 6, 7: 6},
            'Proportion of Features Classified at Depth': {
                1: 1.0, 2: 0.8888888888888888, 3: 0.8888888888888888,
                4: 0.8888888888888888, 5: 0.8888888888888888,
                6: 0.6666666666666666, 7: 0.6666666666666666},
            'Number of Features Unclassified at Depth': {
                1: 0, 2: 1, 3: 1, 4: 1, 5: 1, 6: 3, 7: 3},
            'Proportion of Features Unclassified at Depth': {
                1: 0.0, 2: 0.1111111111111111, 3: 0.1111111111111111,
                4: 0.1111111111111111, 5: 0.1111111111111111,
                6: 0.3333333333333333, 7: 0.3333333333333333}})

        obs = evaluate._evaluate_taxonomy(
            self.taxa, rank_handle_regex='^[dkpcofgs]__')
        pdt.assert_frame_equal(obs, exp, check_names=False)

    def test_evaluate_taxonomy_real_live_taxonomies_unite_rank_handle(self):
        exp = pd.DataFrame({
            'Unique Labels': {
                1: 6.0, 2: 8.0, 3: 8.0, 4: 9.0, 5: 9.0, 6: 9.0, 7: 9.0},
            'Taxonomic Entropy': {1: 1.6769877743224173, 2: 2.0431918705451206,
                                  3: 2.0431918705451206, 4: 2.1972245773362196,
                                  5: 2.1972245773362196, 6: 2.1972245773362196,
                                  7: 2.1972245773362196},
            'Number of Features Terminating at Depth': {
                1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 7},
            'Proportion of Features Terminating at Depth': {
                1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.1111111111111111,
                7: 0.7777777777777778},
            'Number of Features Classified at Depth': {
                1: 8, 2: 8, 3: 8, 4: 8, 5: 8, 6: 8, 7: 7},
            'Proportion of Features Classified at Depth': {
                1: 0.8888888888888888, 2: 0.8888888888888888,
                3: 0.8888888888888888, 4: 0.8888888888888888,
                5: 0.8888888888888888, 6: 0.8888888888888888,
                7: 0.7777777777777778},
            'Number of Features Unclassified at Depth': {
                1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 2},
            'Proportion of Features Unclassified at Depth': {
                1: 0.1111111111111111, 2: 0.1111111111111111,
                3: 0.1111111111111111, 4: 0.1111111111111111,
                5: 0.1111111111111111, 6: 0.1111111111111111,
                7: 0.2222222222222222}})

        obs = evaluate._evaluate_taxonomy(
            self.taxa, rank_handle_regex='^[dkpcofgs]__unidentified')
        pdt.assert_frame_equal(obs, exp, check_names=False)
