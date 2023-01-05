# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
import qiime2
import pandas as pd
import pandas.util.testing as pdt

from rescript import cross_validate


import_data = qiime2.Artifact.import_data


class TestPipelines(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        # drop feature C1b because it is missing species level
        self.taxa_series = pd.read_csv(
            self.get_data_path('derep-taxa.tsv'), sep='\t', index_col=0,
            squeeze=True).drop('C1b')
        self.taxa = import_data('FeatureData[Taxonomy]', self.taxa_series)
        seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))
        self.seqs = import_data(
            'FeatureData[Sequence]', seqs.view(pd.Series).drop('C1b'))
        self.stratified_taxonomy = self.taxa_series.copy().str.replace(
            '; s__brevis', '').str.replace('; s__vaginalis', '').str.replace(
                '; s__pseudocasei', '').sort_index()

    def test_evaluate_cross_validate_k3(self):
        exp, obs, _ = rescript.actions.evaluate_cross_validate(
            self.seqs, self.taxa, k=3)
        # exp_exp (expected ground truth taxonomies)
        # This will equal the original taxonomy except singleton labels will
        # be truncated to reflect stratification.
        exp_exp = self.stratified_taxonomy
        # exp_obs (expected observations)
        exp_obs = pd.Series({
            'A1': palvei,
            'A2': palvei,
            'A3': paeni,
            'A4': palvei,
            'A5': palvei,
            'B1': lcasei,
            'B1a': lcasei,
            'B1b': lacto,
            'B2': lacto,
            'B3': lacto,
            'C1': pdamnosus,
            'C1a': pacidilacti,
            'C1c': pacidilacti,
            'C1d': pacidilacti,
            'C2': pdamnosus}).sort_index()
        pdt.assert_series_equal(
            exp_exp, exp.view(pd.Series).sort_index(), check_names=False)
        pdt.assert_series_equal(
            exp_obs, obs.view(pd.Series).sort_index(), check_names=False)

    def test_evaluate_fit_classifier(self):
        # exp species should equal the input taxonomy when k='disable'
        classifier, evaluation, obs = rescript.actions.evaluate_fit_classifier(
            self.seqs, self.taxa)
        # obs species will equal best possible predictive accuracy.
        exp_obs = pd.Series({
            'A1': palvei,
            'A2': palvei,
            'A3': palvei,
            'A4': palvei,
            'A5': palvei,
            'B1': lcasei,
            'B1a': lcasei,
            'B1b': lacto,
            'B2': lcasei,
            'B3': lcasei,
            'C1': pdamnosus,
            'C1a': pacidilacti,
            'C1c': pacidilacti,
            'C1d': pacidilacti,
            'C2': pdamnosus})
        pdt.assert_series_equal(
            obs.view(pd.Series).sort_index(), exp_obs, check_names=False)

    def test_evaluate_classifications(self):
        # simulate predicted classifications at genus level
        taxa = self.taxa_series.copy().apply(
            lambda x: ';'.join(x.split(';')[:6]))
        taxa = qiime2.Artifact.import_data('FeatureData[Taxonomy]', taxa)
        # first round we just make sure this runs
        rescript.actions.evaluate_classifications([self.taxa], [taxa])
        # now the same but input multiple times to test lists of inputs
        vol, = rescript.actions.evaluate_classifications(
            [self.taxa, taxa], [taxa, taxa])
        # now inspect and validate the contents
        # we inspect the second eval results to compare perfect match vs.
        # simulated genus-level classification (when species are expected)
        vol.export_data(self.temp_dir.name)
        html_path = os.path.join(self.temp_dir.name, 'data.tsv')
        vol = qiime2.Metadata.load(html_path).to_dataframe()
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
        pdt.assert_frame_equal(vol.sort_index(), exp, check_names=False)

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

    def test_validate_even_rank_taxonomy_pass(self):
        taxa = self.taxa.copy().drop('C1b')
        cross_validate._validate_even_rank_taxonomy(taxa)

    def test_validate_even_rank_taxonomy_fail(self):
        with self.assertRaisesRegex(ValueError, "too short: C1b"):
            cross_validate._validate_even_rank_taxonomy(self.taxa)

    def test_validate_indices_match_pass(self):
        cross_validate._validate_indices_match(
            self.taxa.index, self.seqs.index)

    def test_validate_indices_match_fail(self):
        taxa = self.taxa.copy().drop(['A1', 'B1'])
        with self.assertRaisesRegex(ValueError, "one input: A1, B1"):
            cross_validate._validate_indices_match(taxa.index, self.seqs.index)


class TestRelabelStratifiedTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.valid_taxonomies = {
            'k__Bacteria',
            'k__Bacteria; p__Firmicutes',
            'k__Bacteria; p__Firmicutes; c__Bacilli',
            'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales',
            'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
            'f__Lactobacillaceae',
            'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
            'f__Lactobacillaceae; g__Lactobacillus',
            'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
            'f__Lactobacillaceae; g__Lactobacillus; s__casei'}

    def test_relabel_stratified_taxonomy_known_species(self):
        species = ('k__Bacteria; p__Firmicutes; c__Bacilli; '
                   'o__Lactobacillales; f__Lactobacillaceae; '
                   'g__Lactobacillus; s__casei')
        exp = ('k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
               'f__Lactobacillaceae; g__Lactobacillus; s__casei')
        obs = cross_validate._relabel_stratified_taxonomy(
            species, self.valid_taxonomies)
        self.assertEqual(exp, obs)

    def test_relabel_stratified_taxonomy_unknown_species(self):
        species = ('k__Bacteria; p__Firmicutes; c__Bacilli; '
                   'o__Lactobacillales; f__Lactobacillaceae; '
                   'g__Lactobacillus; s__reuteri')
        exp = ('k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; '
               'f__Lactobacillaceae; g__Lactobacillus')
        obs = cross_validate._relabel_stratified_taxonomy(
            species, self.valid_taxonomies)
        self.assertEqual(exp, obs)

    def test_relabel_stratified_taxonomy_unknown_kingdom(self):
        species = 'k__Peanut'
        with self.assertRaisesRegex(RuntimeError, "unknown kingdom"):
            cross_validate._relabel_stratified_taxonomy(
                species, self.valid_taxonomies)


paeni = 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; ' \
        'f__Paenibacillaceae; g__Paenibacillus'
palvei = paeni + '; s__alvei'
lactobacillaceae = 'k__Bacteria; p__Firmicutes; c__Bacilli; ' \
                   'o__Lactobacillales; f__Lactobacillaceae'
lacto = lactobacillaceae + '; g__Lactobacillus'
pedio = lactobacillaceae + '; g__Pediococcus'
lcasei = lacto + '; s__casei'
pdamnosus = pedio + '; s__damnosus'
pacidilacti = pedio + '; s__acidilacti'
