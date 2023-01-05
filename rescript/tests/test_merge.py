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


import_data = qiime2.Artifact.import_data


class TestMergeTaxa(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.merge_taxa = rescript.actions.merge_taxa

        m1 = pd.DataFrame({
            'Taxon': {
                # test merging full species vs. ambiguous level annotations
                '370253': 'k__Bacteria; p__Firmicutes; c__Clostridia; '
                          'o__Clostridiales; f__Ruminococcaceae; '
                          'g__Faecalibacterium; s__prausnitzii',
                # test merging identical annotations
                '2562097': 'k__Bacteria; p__Firmicutes; c__Bacilli; '
                           'o__Bacillales; f__Bacillaceae; g__Bacillus; s__',
                # test merging incomplete vs. full annotations
                '2562091': 'k__Bacteria; p__Actinobacteria; c__Acidimicrobiia;'
                           ' o__Acidimicrobiales; f__Microthrixaceae; g__; '
                           's__',
                # test merging subset vs. superset of same lineage
                '4361279': 'k__Bacteria; p__Proteobacteria; '
                           'c__Betaproteobacteria; o__Burkholderiales; '
                           'f__Oxalobacteraceae; g__; s__',
                '4369464': 'k__Bacteria; p__; c__; o__; f__; g__; s__',
                # test merging unique feature (only found in dataframe a)
                'unique2': 'k__Bacteria; p__; c__; o__; f__; g__; s__blah'},
            'confidence': {
                '370253': 1.0,
                '2562097': 0.99,
                '2562091': 1.0,
                '4361279': 0.9,
                '4369464': 0.8,
                'unique2': 0.77}})
        m1.index.name = 'Feature ID'
        self.m1 = import_data('FeatureData[Taxonomy]', pd.DataFrame(m1))

        # Note: second column is "consensus" in frame b, but "confidence" in a
        m2 = pd.DataFrame({
            'Taxon': {
                '370253': 'k__Bacteria; p__Firmicutes; c__Clostridia; '
                          'o__Clostridiales; f__; g__; s__',
                '2562097': 'k__Bacteria; p__Firmicutes; c__Bacilli; '
                           'o__Bacillales; f__Bacillaceae; g__Bacillus; s__',
                '2562091': 'k__Bacteria; p__Actinobacteria; c__Acidimicrobiia;'
                           ' o__Acidimicrobiales; f__',
                '4361279': 'k__Bacteria; p__Proteobacteria; '
                           'c__Betaproteobacteria; o__Burkholderiales',
                '4369464': 'k__Bacteria; p__Proteobacteria; '
                           'c__Alphaproteobacteria; o__Rhizobiales; '
                           'f__Rhizobiaceae; g__Rhizobium; s__leguminosarum',
                # test merging unique feature (only found in dataframe b)
                'unique1': 'k__Bacteria; p__; c__; o__; f__; g__; s__blah'},
            'Consensus': {
                '370253': 0.95,
                '2562097': 0.9,
                '2562091': 0.85,
                '4361279': 0.87,
                '4369464': 0.99,
                'unique1': 0.75}})
        m2.index.name = 'Feature ID'
        self.m2 = import_data('FeatureData[Taxonomy]', pd.DataFrame(m2))

        # test merging entirely unique dataframe (no overlap with a or b)
        m3 = pd.DataFrame({
            'Taxon': {
                'unique': 'k__Bacteria; p__; c__; o__; f__; g__; s__blah'},
            'confidence': {
                'unique': 0.76}})
        m3.index.name = 'Feature ID'
        self.m3 = import_data('FeatureData[Taxonomy]', pd.DataFrame(m3))

        m4 = pd.DataFrame({
            'Taxon': {
                '370253': 'k__Protozoa; p__; c__; o__; f__; g__; s__blah',
                '2562097': 'k__Bacteria; p__Firmicutes; c__Bacilli; '
                           'o__Bacillales; f__Bacillaceae; g__Bacillus; s__',
                '2562091': 'k__Bacteria; p__Actinobacteria; c__Acidimicrobiia;'
                           ' o__Acidimicrobiales; f__',
                '4361279': 'k__Bacteria; p__Proteobacteria; '
                           'c__Betaproteobacteria; o__Burkholderiales',
                '4369464': 'k__Bacteria; p__Proteobacteria; '
                           'c__Alphaproteobacteria; o__Rhizobiales; '
                           'f__Rhizobiaceae; g__Rhizobium; s__leguminosarum',
                # test merging unique feature (only found in dataframe b)
                'unique1': 'k__Bacteria; p__; c__; o__; f__; g__; s__blah'
            },
            'consensus': {
                '370253': 0.3,
                '2562097': 0.9,
                '2562091': 0.85,
                '4361279': 0.87,
                '4369464': 0.99,
                'unique1': 0.75
            }
        })
        m4.index.name = 'Feature ID'
        self.m4 = import_data('FeatureData[Taxonomy]', pd.DataFrame(m2))

        # same as above, but these are single-column series
        self.s1 = import_data('FeatureData[Taxonomy]', m1['Taxon'])
        self.s2 = import_data('FeatureData[Taxonomy]', m2['Taxon'])
        self.s3 = import_data('FeatureData[Taxonomy]', m3['Taxon'])
        self.s4 = import_data('FeatureData[Taxonomy]', m4['Taxon'])

    def test_merge_taxa_lca(self):
        one_col, = self.merge_taxa([self.s1, self.s2, self.s3], 'lca', '')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                       'o__Acidimicrobiales',
            '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                       'f__Bacillaceae;g__Bacillus;s__',
            '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                      'o__Clostridiales',
            '4361279': 'k__Bacteria;p__Proteobacteria;'
                       'c__Betaproteobacteria;o__Burkholderiales',
            '4369464': 'k__Bacteria',
            'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'}})
        pdt.assert_frame_equal(
            one_col.view(pd.DataFrame), exp, check_names=False)
        multi_col, = self.merge_taxa([self.m1, self.m2, self.m3], 'lca', '')
        pdt.assert_frame_equal(
            multi_col.view(pd.DataFrame), exp, check_names=False)

    def test_merge_taxa_lca_rank_handle(self):
        result, = self.merge_taxa(
            [self.s1, self.s2, self.s3], 'lca')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'Bacteria;Actinobacteria;Acidimicrobiia;'
                       'Acidimicrobiales',
            '2562097': 'Bacteria;Firmicutes;Bacilli;Bacillales;'
                       'Bacillaceae;Bacillus;',
            '370253': 'Bacteria;Firmicutes;Clostridia;Clostridiales',
            '4361279': 'Bacteria;Proteobacteria;'
                       'Betaproteobacteria;Burkholderiales',
            '4369464': 'Bacteria',
            'unique': 'Bacteria;;;;;;blah',
            'unique1': 'Bacteria;;;;;;blah',
            'unique2': 'Bacteria;;;;;;blah'}})
        pdt.assert_frame_equal(
            result.view(pd.DataFrame), exp, check_names=False)

    def test_merge_taxa_lca_rank_handle_plus_new_rank_handle(self):
        result, = self.merge_taxa(
            [self.s1, self.s2, self.s3], 'lca',
            new_rank_handle='greengenes')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                       'o__Acidimicrobiales',
            '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                       'f__Bacillaceae;g__Bacillus;s__',
            '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                      'o__Clostridiales',
            '4361279': 'k__Bacteria;p__Proteobacteria;'
                       'c__Betaproteobacteria;o__Burkholderiales',
            '4369464': 'k__Bacteria',
            'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'}})
        pdt.assert_frame_equal(
            result.view(pd.DataFrame), exp, check_names=False)

    def test_merge_taxa_len(self):
        one_col, = self.merge_taxa([self.s1, self.s2, self.s3], 'len', '')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                       'o__Acidimicrobiales;f__Microthrixaceae;g__;s__',
            '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                       'f__Bacillaceae;g__Bacillus;s__',
            '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                      'o__Clostridiales;f__;g__;s__',
            '4361279': 'k__Bacteria;p__Proteobacteria;c__Betaproteobacteria;'
                       'o__Burkholderiales;f__Oxalobacteraceae;g__;s__',
            '4369464': 'k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;'
                       'o__Rhizobiales;f__Rhizobiaceae;g__Rhizobium;'
                       's__leguminosarum',
            'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'}})
        pdt.assert_frame_equal(
            one_col.view(pd.DataFrame), exp, check_names=False)
        exp = pd.concat([exp, pd.DataFrame({
            'Confidence': {
                '2562091': 1.0, '2562097': np.nan, '370253': np.nan,
                '4361279': 0.9, '4369464': np.nan, 'unique': 0.76,
                'unique1': np.nan, 'unique2': 0.77},
            'Consensus': {
                '2562091': np.nan, '2562097': 0.9, '370253': 0.95,
                '4361279': np.nan, '4369464': 0.99, 'unique': np.nan,
                'unique1': 0.75, 'unique2': np.nan}
        })], axis=1)
        multi_col, = self.merge_taxa([self.m1, self.m2, self.m3], 'len', '')
        multi_col = multi_col.view(pd.DataFrame).apply(
            lambda x: pd.to_numeric(x, errors='ignore'))
        pdt.assert_frame_equal(multi_col, exp, check_names=False)

    def test_merge_taxa_len_rank_handle(self):
        result, = self.merge_taxa(
            [self.m1, self.m2, self.m3], 'len')
        exp = pd.DataFrame({
            'Taxon': {
                '2562091': 'Bacteria;Actinobacteria;Acidimicrobiia;'
                           'Acidimicrobiales;Microthrixaceae;;',
                '2562097': 'Bacteria;Firmicutes;Bacilli;Bacillales;'
                           'Bacillaceae;Bacillus;',
                # note how this expected taxon differs from that in mode=len
                # tests above; this is the benefit of rank_handle, that it
                # can pick the longest lineage after removing ambiguous labels
                '370253': 'Bacteria;Firmicutes;Clostridia;Clostridiales;'
                          'Ruminococcaceae;Faecalibacterium;prausnitzii',
                '4361279': 'Bacteria;Proteobacteria;Betaproteobacteria'
                           ';Burkholderiales;Oxalobacteraceae;;',
                '4369464': 'Bacteria;Proteobacteria;Alphaproteobacteria;'
                           'Rhizobiales;Rhizobiaceae;Rhizobium;leguminosarum',
                'unique': 'Bacteria;;;;;;blah',
                'unique1': 'Bacteria;;;;;;blah',
                'unique2': 'Bacteria;;;;;;blah'},
            'Confidence': {'2562091': 1.0, '2562097': np.nan, '370253': 1.0,
                           '4361279': 0.9, '4369464': np.nan, 'unique': 0.76,
                           'unique1': np.nan, 'unique2': 0.77},
            'Consensus': {'2562091': np.nan, '2562097': 0.9, '370253': np.nan,
                          '4361279': np.nan, '4369464': 0.99, 'unique': np.nan,
                          'unique1': 0.75, 'unique2': np.nan}})
        result = result.view(pd.DataFrame).apply(
            lambda x: pd.to_numeric(x, errors='ignore'))
        pdt.assert_frame_equal(result, exp, check_names=False)

    def test_merge_taxa_len_rank_handle_plus_new_rank_handle(self):
        result, = self.merge_taxa(
            [self.m1, self.m2, self.m3], 'len',
            new_rank_handle='greengenes')
        exp = pd.DataFrame({
            'Taxon': {
                '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                           'o__Acidimicrobiales;f__Microthrixaceae;g__;s__',
                '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;'
                           'o__Bacillales;f__Bacillaceae;g__Bacillus;s__',
                # note how this expected taxon differs from that in mode=len
                # tests above; this is the benefit of rank_handle, that it
                # can pick the longest lineage after removing ambiguous labels
                '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                          'o__Clostridiales;f__Ruminococcaceae;'
                          'g__Faecalibacterium;s__prausnitzii',
                '4361279': 'k__Bacteria;p__Proteobacteria;'
                           'c__Betaproteobacteria;o__Burkholderiales;'
                           'f__Oxalobacteraceae;g__;s__',
                '4369464': 'k__Bacteria;p__Proteobacteria;'
                           'c__Alphaproteobacteria;o__Rhizobiales;'
                           'f__Rhizobiaceae;g__Rhizobium;s__leguminosarum',
                'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
                'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
                'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'},
            'Confidence': {'2562091': 1.0, '2562097': np.nan, '370253': 1.0,
                           '4361279': 0.9, '4369464': np.nan, 'unique': 0.76,
                           'unique1': np.nan, 'unique2': 0.77},
            'Consensus': {'2562091': np.nan, '2562097': 0.9, '370253': np.nan,
                          '4361279': np.nan, '4369464': 0.99, 'unique': np.nan,
                          'unique1': 0.75, 'unique2': np.nan}})
        result = result.view(pd.DataFrame).apply(
            lambda x: pd.to_numeric(x, errors='ignore'))
        pdt.assert_frame_equal(result, exp, check_names=False)

    # score should fail if there is only one column in any dataset
    def test_merge_taxa_score_one_column(self):
        with self.assertRaisesRegex(IndexError, "second column"):
            self.merge_taxa([self.m1, self.s2], 'score', '')

    def test_merge_taxa_score_multi_column(self):
        new, = self.merge_taxa([self.m1, self.m2, self.m3], 'score', '')
        exp = pd.DataFrame({
            'Taxon': {
                '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                           'o__Acidimicrobiales;f__Microthrixaceae;g__;s__',
                '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales'
                           ';f__Bacillaceae;g__Bacillus;s__',
                '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                          'o__Clostridiales;f__Ruminococcaceae;'
                          'g__Faecalibacterium;s__prausnitzii',
                '4361279': 'k__Bacteria;p__Proteobacteria;'
                           'c__Betaproteobacteria;o__Burkholderiales;'
                           'f__Oxalobacteraceae;g__;s__',
                '4369464': 'k__Bacteria;p__Proteobacteria;'
                           'c__Alphaproteobacteria;o__Rhizobiales;'
                           'f__Rhizobiaceae;g__Rhizobium;s__leguminosarum',
                'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
                'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
                'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'},
            'Confidence': {
                '2562091': 1.0, '2562097': 0.99, '370253': 1.0, '4361279': 0.9,
                '4369464': 0, 'unique': 0.76, 'unique1': 0, 'unique2': 0.77},
            'Consensus': {
                '2562091': 0, '2562097': 0, '370253': 0, '4361279': 0,
                '4369464': 0.99, 'unique': 0, 'unique1': 0.75, 'unique2': 0},
            'Score': {
                '2562091': 1.0, '2562097': 0.99, '370253': 1.0, '4361279': 0.9,
                '4369464': 0.99, 'unique': 0.76, 'unique1': 0.75,
                'unique2': 0.77}})
        new = new.view(pd.DataFrame).apply(
            lambda x: pd.to_numeric(x, errors='ignore'))
        pdt.assert_frame_equal(new, exp, check_names=False)

    # compare vs. LCA tests above to see super(set) mode function; the longer
    # taxonomy is selected if it is a superset of the other, otherwise LCA
    def test_merge_taxa_super_lca(self):
        one_col, = self.merge_taxa([self.s1, self.s2, self.s3], 'super', '')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                       'o__Acidimicrobiales;f__Microthrixaceae;g__;s__',
            '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                       'f__Bacillaceae;g__Bacillus;s__',
            '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                      'o__Clostridiales;f__Ruminococcaceae;'
                      'g__Faecalibacterium;s__prausnitzii',
            '4361279': 'k__Bacteria;p__Proteobacteria;'
                       'c__Betaproteobacteria;o__Burkholderiales;'
                       'f__Oxalobacteraceae;g__;s__',
            '4369464': 'k__Bacteria;p__Proteobacteria;'
                       'c__Alphaproteobacteria;o__Rhizobiales;'
                       'f__Rhizobiaceae;g__Rhizobium;s__leguminosarum',
            'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'}})
        pdt.assert_frame_equal(
            one_col.view(pd.DataFrame), exp, check_names=False)
        multi_col, = self.merge_taxa([self.m1, self.m2, self.m3], 'super', '')
        pdt.assert_frame_equal(
            multi_col.view(pd.DataFrame), exp, check_names=False)

    # compare vs. LCA tests above to see super(set) mode function; the longer
    # taxonomy is selected if it is a superset of the other, otherwise LCA
    def test_merge_taxa_super_lca_rank_handle(self):
        result, = self.merge_taxa(
            [self.s1, self.s2, self.s3], 'super')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'Bacteria;Actinobacteria;Acidimicrobiia;'
                       'Acidimicrobiales;Microthrixaceae;;',
            '2562097': 'Bacteria;Firmicutes;Bacilli;Bacillales;'
                       'Bacillaceae;Bacillus;',
            '370253': 'Bacteria;Firmicutes;Clostridia;Clostridiales;'
                      'Ruminococcaceae;Faecalibacterium;prausnitzii',
            '4361279': 'Bacteria;Proteobacteria;'
                       'Betaproteobacteria;Burkholderiales;Oxalobacteraceae;;',
            '4369464': 'Bacteria;Proteobacteria;Alphaproteobacteria;'
                       'Rhizobiales;Rhizobiaceae;Rhizobium;leguminosarum',
            'unique': 'Bacteria;;;;;;blah',
            'unique1': 'Bacteria;;;;;;blah',
            'unique2': 'Bacteria;;;;;;blah'}})
        pdt.assert_frame_equal(
            result.view(pd.DataFrame), exp, check_names=False)
        multi_col, = self.merge_taxa([self.m1, self.m2, self.m3], 'super')
        pdt.assert_frame_equal(
            multi_col.view(pd.DataFrame), exp, check_names=False)

    def test_merge_taxa_super_lca_with_unassigned(self):
        result, = self.merge_taxa([self.s1, self.s4], 'super')
        exp = pd.DataFrame({
            'Taxon': {
                '2562091': 'Bacteria;Actinobacteria;Acidimicrobiia;'
                           'Acidimicrobiales;Microthrixaceae;;',
                '2562097': 'Bacteria;Firmicutes;Bacilli;Bacillales;'
                           'Bacillaceae;Bacillus;',
                '370253':  'Unassigned',
                '4361279': 'Bacteria;Proteobacteria;'
                           'Betaproteobacteria;Burkholderiales;'
                           'Oxalobacteraceae;;',
                '4369464': 'Bacteria;Proteobacteria;Alphaproteobacteria;'
                           'Rhizobiales;Rhizobiaceae;Rhizobium;leguminosarum',
                'unique1': 'Bacteria;;;;;;blah',
                'unique2': 'Bacteria;;;;;;blah'
            }
        })
        pdt.assert_frame_equal(
            result.view(pd.DataFrame), exp, check_names=False)

    # majority is super without superstring collapsing
    def test_merge_taxa_super_lca_majority(self):
        one_col, = self.merge_taxa([self.s1, self.s2, self.s3], 'majority', '')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                       'o__Acidimicrobiales',
            '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                       'f__Bacillaceae;g__Bacillus;s__',
            '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                      'o__Clostridiales',
            '4361279': 'k__Bacteria;p__Proteobacteria;'
                       'c__Betaproteobacteria;o__Burkholderiales;'
                       'f__Oxalobacteraceae;g__;s__',
            '4369464': 'k__Bacteria',
            'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'}})
        pdt.assert_frame_equal(
            one_col.view(pd.DataFrame), exp, check_names=False)
        multi_col, = self.merge_taxa(
            [self.m1, self.m2, self.m3], 'majority', '')
        pdt.assert_frame_equal(
            multi_col.view(pd.DataFrame), exp, check_names=False)

    # note that majority with rank handle behaves the same as super LCA above
    # (same expected outputs) since the super/substrings in the test data
    # consist of empty ranks vs. rank handles.
    def test_merge_taxa_majority_rank_handle(self):
        result, = self.merge_taxa(
            [self.s1, self.s2, self.s3], 'majority')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'Bacteria;Actinobacteria;Acidimicrobiia;'
                       'Acidimicrobiales;Microthrixaceae;;',
            '2562097': 'Bacteria;Firmicutes;Bacilli;Bacillales;'
                       'Bacillaceae;Bacillus;',
            '370253': 'Bacteria;Firmicutes;Clostridia;Clostridiales;'
                      'Ruminococcaceae;Faecalibacterium;prausnitzii',
            '4361279': 'Bacteria;Proteobacteria;'
                       'Betaproteobacteria;Burkholderiales;Oxalobacteraceae;;',
            '4369464': 'Bacteria;Proteobacteria;Alphaproteobacteria;'
                       'Rhizobiales;Rhizobiaceae;Rhizobium;leguminosarum',
            'unique': 'Bacteria;;;;;;blah',
            'unique1': 'Bacteria;;;;;;blah',
            'unique2': 'Bacteria;;;;;;blah'}})
        pdt.assert_frame_equal(
            result.view(pd.DataFrame), exp, check_names=False)
        multi_col, = self.merge_taxa([self.m1, self.m2, self.m3], 'majority')
        pdt.assert_frame_equal(
            multi_col.view(pd.DataFrame), exp, check_names=False)

    # make a more challenging test for majority with label redundancy and
    # no clear majority
    def test_merge_taxa_majority_rank_handle_true_majority(self):
        s4 = pd.DataFrame({
            'Taxon': {
                # ensure majority works (genus level) and at species level
                # LCA works on 3-way disagreement (no majority)
                '370253': 'k__Bacteria; p__Firmicutes; c__Clostridia; '
                          'o__Clostridiales; f__Ruminococcaceae; '
                          'g__Faecalibacterium; s__not_prausnitzii',
                # ensure that empty ranks are ignored by LCA majority
                '2562097': 'k__Bacteria; p__Firmicutes; c__Bacilli; '
                           'o__Bacillales; f__Bacillaceae; g__Bacillus',
                # majority rules
                '2562091': 'k__Bacteria; p__Actinobacteria; c__Acidimicrobiia;'
                           ' o__Acidimicrobiales; f__Microthrixaceae; g__; '
                           's__'}})
        s4.index.name = 'Feature ID'
        s4 = import_data('FeatureData[Taxonomy]', pd.DataFrame(s4))
        result, = self.merge_taxa(
            [self.s1, self.s2, self.s3, s4], 'majority', '')
        exp = pd.DataFrame({'Taxon': {
            '2562091': 'k__Bacteria;p__Actinobacteria;c__Acidimicrobiia;'
                       'o__Acidimicrobiales;f__Microthrixaceae;g__;s__',
            '2562097': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                       'f__Bacillaceae;g__Bacillus;s__',
            '370253': 'k__Bacteria;p__Firmicutes;c__Clostridia;'
                      'o__Clostridiales;f__Ruminococcaceae;'
                      'g__Faecalibacterium',
            '4361279': 'k__Bacteria;p__Proteobacteria;'
                       'c__Betaproteobacteria;o__Burkholderiales;'
                       'f__Oxalobacteraceae;g__;s__',
            '4369464': 'k__Bacteria',
            'unique': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique1': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah',
            'unique2': 'k__Bacteria;p__;c__;o__;f__;g__;s__blah'}})
        pdt.assert_frame_equal(
            result.view(pd.DataFrame), exp, check_names=False)
