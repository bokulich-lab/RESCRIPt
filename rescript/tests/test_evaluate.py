# ----------------------------------------------------------------------------
# Copyright (c) 2019--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
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


class TestEvaluateTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))

    # this just tests that the pipeline runs, other tests test proper function
    def test_pipeline(self):
        rescript.actions.evaluate_taxonomy([self.taxa], "")

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
        obs_ent = evaluate._taxonomic_entropy(self.taxa.view(pd.Series), "")
        exp_ent = pd.DataFrame(
            {'Unique Labels': {
                1: 1.0, 2: 1.0, 3: 1.0, 4: 2.0, 5: 2.0, 6: 3.0, 7: 8.0},
             'Taxonomic Entropy': {
                1: 0.0,
                2: 0.0,
                3: 0.0,
                4: 0.6210863745552451,
                5: 0.6210863745552451,
                6: 1.0986122886681096,
                7: 1.9913464134109882}})
        pdt.assert_frame_equal(obs_ent.sort_index(),
                               exp_ent.sort_index(), check_names=False)
