# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
import qiime2
import pandas as pd
import numpy as np
import pandas.util.testing as pdt
from q2_types.feature_data import DNAIterator

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


class TestEvaluateSeqs(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.seqs = import_data('FeatureData[Sequence]',
                                self.get_data_path('derep-test.fasta'))

    # this just tests that the action runs, other tests test proper function
    def test_evaluate_seqs_visualizer(self):
        rescript.actions.evaluate_seqs([self.seqs])

    def test_evaluate_seqs(self):
        obs, lens = evaluate._evaluate_seqs(
            [self.seqs.view(DNAIterator)], ["name"])
        exp = pd.DataFrame({
            'name': {'Length min': 264.0, 'Length 1%': 264.0,
                     'Length 25%': 284.25, 'Length median': 291.0,
                     'Length 75%': 291.0, 'Length 99%': 291.0,
                     'Length max': 291.0, 'N uniques': 6.0,
                     'Sequence Entropy': 1.63}})
        pdt.assert_frame_equal(obs.sort_index(), exp.sort_index(),
                               check_names=False)
        np.testing.assert_array_equal(lens, np.array([[291] * 12 + [264] * 4]))

    def test_evaluate_seqs_low_entropy_multiple_sets(self):
        s1 = ['ACTGATCGTGATGCTGATCGATGCTGATCGATCG',
              'GTGTGTGAGTTATCGTGACGTGTAGCTGACGTAG',
              'ACGTGTACTGTGACTGATGCTGACTGTGGTATAT',
              'ACGAGTCTGAC',
              'ACGTGTACGTGTAGCTGTAGC',
              'CGTTGATGCTGTGATGCTACTGTGACTGATGCGTAGCGTAC']
        s2 = ['ACTGATCGTGATGCTGATCGATGCTGATCGATCG',
              'AAAAAAAAAAAAAAAAAAAAAAA',
              'AAAAAAAAAAAAAAAAAAAAAAAAAAAA']
        s3 = ['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
              'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
              'AAAAAAAAAAAAAAAAAAAAAAAAAAAA']
        obs, lens = evaluate._evaluate_seqs([s1, s2, s3], ['s1', 's2', 's3'])
        exp = pd.DataFrame({
            's1': {'Length min': 11.0, 'Length 1%': 11.5, 'Length 25%': 24.25,
                   'Length median': 34.0, 'Length 75%': 34.0,
                   'Length 99%': 40.65, 'Length max': 41.0, 'N uniques': 6.0,
                   'Sequence Entropy': 1.79},
            's2': {'Length min': 23.0, 'Length 1%': 23.1, 'Length 25%': 25.5,
                   'Length median': 28.0, 'Length 75%': 31.0,
                   'Length 99%': 33.88, 'Length max': 34.0, 'N uniques': 3.0,
                   'Sequence Entropy': 1.1},
            's3': {'Length min': 28.0, 'Length 1%': 28.12, 'Length 25%': 31.0,
                   'Length median': 34.0, 'Length 75%': 34.0,
                   'Length 99%': 34.0, 'Length max': 34.0, 'N uniques': 2.0,
                   'Sequence Entropy': 0.64}})
        exp_lens = [np.array([34, 34, 34, 11, 21, 41]),
                    np.array([34, 23, 28]),
                    np.array([34, 34, 28])]
        pdt.assert_frame_equal(obs.sort_index(), exp.sort_index(),
                               check_names=False)
        for a1, a2 in zip(lens, exp_lens):
            np.testing.assert_array_equal(a1, a2)
