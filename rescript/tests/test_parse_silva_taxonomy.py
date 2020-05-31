# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from qiime2.plugin.testing import TestPluginBase
from rescript.parse_silva_taxonomy import (parse_silva_taxonomy,
                                           _keep_allowed_chars, _prep_taxranks,
                                           _prep_taxmap, ALLOWED_RANKS,
                                           SELECTED_RANKS,
                                           _build_base_silva_taxonomy,
                                           _validate_taxrank_taxtree,
                                           _compile_taxonomy_output)
from skbio.tree import TreeNode
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal


class TestParseSilvaTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.bad_char_str = '~`!@#$%^&*;:\"?<>'
        self.good_char_str = '1 am a G00d str1ng.'
        self.taxmap = self.get_data_path('taxmap_slv_ssu_ref_nr_test.txt')
        # note the tax ranks and tree files below need to be consistant,
        # i.e., share 100% of the taxids for later tests to work!
        # Thus the resulting taxmap file is limited to two example
        # accessions with the matching taxonomy. Hence taxmap2.
        # We'll keep taxmap 1 for better test case with more data.
        self.taxranks = self.get_data_path('tax_slv_ssu_test.txt')
        self.taxonomy_tree = self.get_data_path('taxid_tree.tre')
        self.taxmap2 = self.get_data_path('taxmap_test_match_tree.txt')
        self.taxonomy_tree2 = self.get_data_path('taxid_tree_missing_id.tre')

    def test_keep_allowed_chars(self):
        obs_bad = _keep_allowed_chars(self.bad_char_str)
        exp_bad = ''
        self.assertEqual(obs_bad, exp_bad)
        obs_good = _keep_allowed_chars(self.good_char_str)
        exp_good = '1_am_a_G00d_str1ng.'
        self.assertEqual(obs_good, exp_good)

    def test_prep_taxranks(self):
        input_taxranks = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxonomy]', self.taxranks)
        input_taxranks = input_taxranks.view(pd.DataFrame)
        obs_taxranks = _prep_taxranks(input_taxranks)
        obs_taxranks.sort_index(inplace=True)
        dd = {'taxid': ['2', '11084', '42913', '42914', '42915',
                        '11089', '24228', '24229', '42916', '42917',
                        '42918'],
              'taxid_taxonomy': ['Archaea', 'Aenigmarchaeota',
                                 'Aenigmarchaeia', 'Aenigmarchaeales',
                                 'Candidatus_Aenigmarchaeum',
                                 'Deep_Sea_Euryarchaeotic_Group(DSEG)',
                                 'Altiarchaeota', 'Altiarchaeia',
                                 'Altiarchaeales', 'Altiarchaeaceae',
                                 'Candidatus_Altiarchaeum'],
              'taxrank': ['domain', 'phylum', 'class', 'order', 'genus',
                          'class', 'phylum', 'class', 'order', 'family',
                          'genus']}
        exp_taxranks = pd.DataFrame(dd)
        exp_taxranks.set_index('taxid', inplace=True)
        exp_taxranks.sort_index(inplace=True)
        assert_frame_equal(obs_taxranks, exp_taxranks)

    def test_prep_taxmap(self):
        input_taxmap = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxidMap]', self.taxmap)
        input_taxmap = input_taxmap.view(pd.DataFrame)
        obs_taxmap = _prep_taxmap(input_taxmap)
        obs_taxmap.sort_index(inplace=True)
        dm = {'Feature ID': ['A16379.1.1485', 'A45315.1.1521', 'A61579.1.1437',
                             'AAAA02020713.1.1297'],
              'organism_name': ['[Haemophilus]_ducreyi', 'Bacillus_sp.',
                                'Thermopallium_natronophilum',
                                'Oryza_sativa'],
              'taxid': ['3698', '45177', '46692', '46463']}
        exp_taxmap = pd.DataFrame(dm)
        exp_taxmap.set_index('Feature ID', inplace=True)
        exp_taxmap.sort_index(inplace=True)
        assert_frame_equal(obs_taxmap, exp_taxmap)

    def test_build_base_silva_taxonomy(self):
        input_taxranks = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxonomy]', self.taxranks)
        input_taxranks = input_taxranks.view(pd.DataFrame)
        input_taxranks = _prep_taxranks(input_taxranks)
        input_taxtree = qiime2.Artifact.import_data(
                            'Phylogeny[Rooted]', self.taxonomy_tree)
        input_taxtree = input_taxtree.view(TreeNode)
        obs_taxonomy = _build_base_silva_taxonomy(input_taxtree,
                                                  input_taxranks,
                                                  ALLOWED_RANKS)
        obs_taxonomy.sort_index(inplace=True)
        tid = {'taxid': ['2', '11084', '42913', '42914', '42915',
                         '11089', '24228', '24229', '42916', '42917', '42918'],
               'd__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                       'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                       'Archaea'],
               'sk__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea'],
               'k__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                       'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                       'Archaea'],
               'ks__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea'],
               'sp__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea'],
               'p__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                       'Aenigmarchaeota', 'Aenigmarchaeota', 'Aenigmarchaeota',
                       'Altiarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                       'Altiarchaeota', 'Altiarchaeota'],
               'ps__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                        'Altiarchaeota', 'Altiarchaeota', 'Altiarchaeota'],
               'pi__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                        'Altiarchaeota', 'Altiarchaeota', 'Altiarchaeota'],
               'sc__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                        'Altiarchaeota', 'Altiarchaeota', 'Altiarchaeota'],
               'c__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeia', 'Aenigmarchaeia',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia',
                       'Altiarchaeia'],
               'cs__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeia', 'Aenigmarchaeia',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia',
                        'Altiarchaeia'],
               'ci__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeia', 'Aenigmarchaeia',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia',
                        'Altiarchaeia'],
               'so__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeia', 'Aenigmarchaeia',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia',
                        'Altiarchaeia'],
               'o__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeales', 'Aenigmarchaeales',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeales',
                       'Altiarchaeales'],
               'os__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeales', 'Aenigmarchaeales',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeales',
                        'Altiarchaeales'],
               'sf__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeales', 'Aenigmarchaeales',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeales',
                        'Altiarchaeales'],
               'f__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeales', 'Aenigmarchaeales',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeaceae',
                       'Altiarchaeaceae'],
               'fs__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeales', 'Aenigmarchaeales',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeaceae',
                        'Altiarchaeaceae'],
               'g__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeales', 'Candidatus_Aenigmarchaeum',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeaceae',
                       'Candidatus_Altiarchaeum']}
        exp_taxonomy = pd.DataFrame(tid)
        exp_taxonomy.set_index('taxid', inplace=True)
        exp_taxonomy.sort_index(inplace=True)
        assert_frame_equal(obs_taxonomy, exp_taxonomy)

    def test_compile_taxonomy_output_default(self):
        # process taxrank data
        input_taxrank = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxonomy]', self.taxranks)
        input_taxrank = input_taxrank.view(pd.DataFrame)
        input_taxrank = _prep_taxranks(input_taxrank)
        # process taxonomy tree
        input_taxtree = qiime2.Artifact.import_data(
                            'Phylogeny[Rooted]', self.taxonomy_tree)
        input_taxtree = input_taxtree.view(TreeNode)
        # make silva tax
        silva_tax = _build_base_silva_taxonomy(input_taxtree, input_taxrank,
                                               ALLOWED_RANKS)
        # process taxmap data
        input_taxmap = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxidMap]', self.taxmap2)
        input_taxmap = input_taxmap.view(pd.DataFrame)
        input_taxmap = _prep_taxmap(input_taxmap)
        # merge taxmap and taxrank info
        updated_taxmap = pd.merge(input_taxmap, silva_tax, left_on='taxid',
                                  right_index=True)
        # compile_taxonomy
        obs_6r_tax = _compile_taxonomy_output(updated_taxmap,
                                              include_species_labels=False,
                                              selected_ranks=SELECTED_RANKS)
        obs_6r_tax.sort_index(inplace=True)
        # expected 6-rank taxonomy
        t1 = ''.join('d__Archaea; p__Aenigmarchaeota; c__Aenigmarchaeia; '
                     'o__Aenigmarchaeales; f__Aenigmarchaeales; '
                     'g__Candidatus_Aenigmarchaeum')
        t2 = ''.join('d__Archaea; p__Altiarchaeota; c__Altiarchaeia; '
                     'o__Altiarchaeales; f__Altiarchaeaceae; '
                     'g__Candidatus_Altiarchaeum')
        exp_6r_tax = pd.Series([t1, t2],
                               index=['AB600437.1.1389', 'AB301876.1.930'])
        exp_6r_tax.rename('Taxon', inplace=True)
        exp_6r_tax.index.name = 'Feature ID'
        exp_6r_tax.sort_index(inplace=True)
        assert_series_equal(obs_6r_tax, exp_6r_tax)

    def test_validate_taxrank_taxtree_fail(self):
        # taxonomy ranks
        input_taxrank = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxonomy]', self.taxranks)
        input_taxrank = input_taxrank.view(pd.DataFrame)
        input_taxrank = _prep_taxranks(input_taxrank)
        # tree
        input_taxtree = qiime2.Artifact.import_data(
                            'Phylogeny[Rooted]', self.taxonomy_tree2)
        input_taxtree = input_taxtree.view(TreeNode)
        self.assertRaises(ValueError, _validate_taxrank_taxtree,
                          input_taxrank, input_taxtree)

    def test_parse_silva_taxonomy(self):
        # taxrank
        input_taxrank = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxonomy]', self.taxranks)
        input_taxrank = input_taxrank.view(pd.DataFrame)
        # tree
        input_taxtree = qiime2.Artifact.import_data(
                            'Phylogeny[Rooted]', self.taxonomy_tree)
        input_taxtree = input_taxtree.view(TreeNode)
        # taxmap
        input_taxmap = qiime2.Artifact.import_data(
                            'FeatureData[SILVATaxidMap]', self.taxmap2)
        input_taxmap = input_taxmap.view(pd.DataFrame)
        # observed
        obs_res = parse_silva_taxonomy(input_taxtree, input_taxmap,
                                       input_taxrank,
                                       include_species_labels=True)
        obs_res.sort_index(inplace=True)
        # expected:
        t1 = ''.join('d__Archaea; p__Aenigmarchaeota; c__Aenigmarchaeia; '
                     'o__Aenigmarchaeales; f__Aenigmarchaeales; '
                     'g__Candidatus_Aenigmarchaeum; s__uncultured_archaeon')
        t2 = ''.join('d__Archaea; p__Altiarchaeota; c__Altiarchaeia; '
                     'o__Altiarchaeales; f__Altiarchaeaceae; '
                     'g__Candidatus_Altiarchaeum; s__uncultured_archaeon')
        exp_res = pd.Series([t1, t2],
                            index=['AB600437.1.1389', 'AB301876.1.930'])
        exp_res.rename('Taxon', inplace=True)
        exp_res.index.name = 'Feature ID'
        exp_res.sort_index(inplace=True)
        assert_series_equal(obs_res, exp_res)
