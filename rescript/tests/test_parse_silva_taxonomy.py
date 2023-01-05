# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import pkg_resources
from qiime2.plugin.testing import TestPluginBase
from rescript.parse_silva_taxonomy import (parse_silva_taxonomy,
                                           _keep_allowed_chars, _prep_taxranks,
                                           _prep_taxmap, ALLOWED_RANKS,
                                           DEFAULT_RANKS,
                                           _build_base_silva_taxonomy,
                                           _validate_taxrank_taxmap_taxtree,
                                           _compile_taxonomy_output,
                                           _get_clean_organism_name,
                                           _get_terminal_taxon)
from skbio.tree import TreeNode
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal


class TestParseSilvaTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        # taxonomy mapping file 1
        tm_path = pkg_resources.resource_filename('rescript.types.tests',
                                                  'data/silva_taxamap.tsv')
        tm = qiime2.Artifact.import_data('FeatureData[SILVATaxidMap]', tm_path)
        self.taxmap = tm.view(pd.DataFrame)
        # Note: below, the tax ranks and tree files below need to be
        # consistant, i.e., share 100% of the taxids for later tests to work!
        # Thus the resulting taxmap file is limited to two example
        # accessions with the matching taxonomy. Hence taxmap2.
        # We'll keep taxmap 1 for better test case with more data.
        # So self.taxmap2, self.taxtree2, and self.tr must match.
        # silva taxonomy ranks file:
        tr_path = pkg_resources.resource_filename('rescript.types.tests',
                                                  'data/silva_taxa.tsv')
        tr = qiime2.Artifact.import_data('FeatureData[SILVATaxonomy]', tr_path)
        self.taxranks = tr.view(pd.DataFrame)
        # silva taxonomy tree file 1
        tt = qiime2.Artifact.import_data(
            'Phylogeny[Rooted]', self.get_data_path('taxid_tree.tre'))
        self.taxtree = tt.view(TreeNode)
        # taxonomy mapping file 2
        tm2 = qiime2.Artifact.import_data(
            'FeatureData[SILVATaxidMap]',
            self.get_data_path('taxmap_test_match_tree.txt'))
        self.taxmap2 = tm2.view(pd.DataFrame)
        # taxonomy tree file with missing taxid:
        tt2 = qiime2.Artifact.import_data(
            'Phylogeny[Rooted]',
            self.get_data_path('taxid_tree_missing_id.tre'))
        self.taxtree2 = tt2.view(TreeNode)

    def test_keep_allowed_chars(self):
        obs_bad = _keep_allowed_chars('~`!@#$%^&*;:\"?<>')
        exp_bad = ''
        self.assertEqual(obs_bad, exp_bad)
        obs_good = _keep_allowed_chars('1 am a G00d str1ng.')
        exp_good = '1_am_a_G00d_str1ng.'
        self.assertEqual(obs_good, exp_good)
        se_1 = _keep_allowed_chars('Glaciecola sp. #105')
        sx_1 = 'Glaciecola_sp._105'
        self.assertEqual(se_1, sx_1)
        se_2 = _keep_allowed_chars(
            'Acetobacterium sp. enrichment culture clone DhR^2/LM-A07')
        sx_2 = 'Acetobacterium_sp._enrichment_culture_clone_DhR2/LM-A07'
        self.assertEqual(se_2, sx_2)
        se_3 = _keep_allowed_chars('Pseudomonas poae RE*1-1-14')
        sx_3 = 'Pseudomonas_poae_RE1-1-14'
        self.assertEqual(se_3, sx_3)
        se_4 = _keep_allowed_chars('Yersinia pseudotuberculosis PB1/+')
        sx_4 = 'Yersinia_pseudotuberculosis_PB1/'
        self.assertEqual(se_4, sx_4)

    def test_get_terminal_taxon(self):
        term_tax = _get_terminal_taxon(
            'Archaea;Aenigmarchaeota;Aenigmarchaeia;Aenigmarchaeales;')
        exp_term_tax = 'Aenigmarchaeales'
        self.assertEqual(term_tax, exp_term_tax)

    def test_prep_taxranks(self):
        obs_taxranks = _prep_taxranks(self.taxranks)
        obs_taxranks.sort_index(inplace=True)
        dd = {'taxid': ['2', '11084', '42913', '42914', '42915',
                        '11089', '24228', '24229', '42916', '42917'],
              'taxid_taxonomy': ['Archaea', 'Aenigmarchaeota',
                                 'Aenigmarchaeia', 'Aenigmarchaeales',
                                 'Candidatus_Aenigmarchaeum',
                                 'Deep_Sea_Euryarchaeotic_Group(DSEG)',
                                 'Altiarchaeota', 'Altiarchaeia',
                                 'Altiarchaeales', 'Altiarchaeaceae'],
              'taxrank': ['domain', 'phylum', 'class', 'order', 'genus',
                          'class', 'phylum', 'class', 'order', 'family']}
        exp_taxranks = pd.DataFrame(dd)
        exp_taxranks.set_index('taxid', inplace=True)
        exp_taxranks.sort_index(inplace=True)
        assert_frame_equal(obs_taxranks, exp_taxranks)

    def test_get_clean_organism_name(self):
        org_name = 'Acetobacterium sp. enrichment culture clone DhR^2/LM-A07'
        obs_clean_name = _get_clean_organism_name(org_name)
        exp_clean_name = 'Acetobacterium_sp.'
        self.assertEqual(obs_clean_name, exp_clean_name)

    def test_prep_taxmap(self):
        obs_taxmap = _prep_taxmap(self.taxmap)
        obs_taxmap.sort_index(inplace=True)
        dm = {'Feature ID': ['A16379.1.1485', 'A45315.1.1521', 'A61579.1.1437',
                             'AAAA02020713.1.1297', 'AAAA02020714.1.1202',
                             'AAAA02038450.2584.4394', 'AAAA02039541.11.1886',
                             'AAAA02041579.2617.4209', 'AAAA02046270.117.1956',
                             'AAAA02048270.689.2185'],
              'organism_name': ['[Haemophilus]_ducreyi', 'Bacillus_sp.',
                                'Thermopallium_natronophilum',
                                'Oryza_sativa', 'Oryza_sativa', 'Oryza_sativa',
                                'Oryza_sativa', 'Oryza_sativa', 'Oryza_sativa',
                                'Oryza_sativa'],
              'taxid': ['3698', '45177', '46692', '46463', '2852', '10099',
                        '47183', '4432', '4432', '44317']}
        exp_taxmap = pd.DataFrame(dm)
        exp_taxmap.set_index('Feature ID', inplace=True)
        exp_taxmap.sort_index(inplace=True)
        assert_frame_equal(obs_taxmap, exp_taxmap)

    def test_build_base_silva_taxonomy(self):
        input_taxranks = _prep_taxranks(self.taxranks)
        obs_taxonomy = _build_base_silva_taxonomy(self.taxtree,
                                                  input_taxranks,
                                                  ALLOWED_RANKS,
                                                  rank_propagation=True)
        obs_taxonomy.sort_index(inplace=True)
        tid = {'taxid': ['2', '11084', '42913', '42914', '42915',
                         '11089', '24228', '24229', '42916', '42917'],
               'd__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                       'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea'],
               'sk__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea'],
               'k__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                       'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea'],
               'ks__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea'],
               'sp__': ['Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea',
                        'Archaea', 'Archaea', 'Archaea', 'Archaea', 'Archaea'],
               'p__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                       'Aenigmarchaeota', 'Aenigmarchaeota', 'Aenigmarchaeota',
                       'Altiarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                       'Altiarchaeota'],
               'ps__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                        'Altiarchaeota', 'Altiarchaeota'],
               'pi__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                        'Altiarchaeota', 'Altiarchaeota'],
               'sc__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Aenigmarchaeota',
                        'Aenigmarchaeota', 'Altiarchaeota', 'Altiarchaeota',
                        'Altiarchaeota', 'Altiarchaeota'],
               'c__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeia', 'Aenigmarchaeia',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia'],
               'cs__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeia', 'Aenigmarchaeia',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia'],
               'ci__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeia', 'Aenigmarchaeia',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia'],
               'so__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeia', 'Aenigmarchaeia',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeia', 'Altiarchaeia'],
               'o__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeales', 'Aenigmarchaeales',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeales'],
               'os__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeales', 'Aenigmarchaeales',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeales'],
               'sf__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeales', 'Aenigmarchaeales',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeales'],
               'f__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeales', 'Aenigmarchaeales',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeaceae'],
               'fs__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                        'Aenigmarchaeales', 'Aenigmarchaeales',
                        'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                        'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeaceae'],
               'g__': ['Archaea', 'Aenigmarchaeota', 'Aenigmarchaeia',
                       'Aenigmarchaeales', 'Candidatus_Aenigmarchaeum',
                       'Deep_Sea_Euryarchaeotic_Group(DSEG)', 'Altiarchaeota',
                       'Altiarchaeia', 'Altiarchaeales', 'Altiarchaeaceae']}
        exp_taxonomy = pd.DataFrame(tid)
        exp_taxonomy.set_index('taxid', inplace=True)
        exp_taxonomy.sort_index(inplace=True)
        assert_frame_equal(obs_taxonomy, exp_taxonomy)

    def test_compile_taxonomy_output_default(self):
        input_taxrank = _prep_taxranks(self.taxranks)
        silva_tax = _build_base_silva_taxonomy(self.taxtree, input_taxrank,
                                               ALLOWED_RANKS,
                                               rank_propagation=True)
        input_taxmap = _prep_taxmap(self.taxmap2)
        updated_taxmap = pd.merge(input_taxmap, silva_tax, left_on='taxid',
                                  right_index=True)
        obs_6r_tax = _compile_taxonomy_output(updated_taxmap,
                                              ranks=DEFAULT_RANKS,
                                              include_species_labels=False)
        obs_6r_tax.sort_index(inplace=True)
        # expected 6-rank taxonomy
        t1 = ("d__Archaea; p__Aenigmarchaeota; c__Aenigmarchaeia; "
              "o__Aenigmarchaeales; f__Aenigmarchaeales; "
              "g__Candidatus_Aenigmarchaeum")
        exp_6r_tax = pd.Series(t1, index=['AB600437.1.1389'])
        exp_6r_tax.rename('Taxon', inplace=True)
        exp_6r_tax.index.name = 'Feature ID'
        exp_6r_tax.sort_index(inplace=True)
        assert_series_equal(obs_6r_tax, exp_6r_tax)

    def test_validate_taxrank_taxmap_taxtree_pass(self):
        # all files should meet the criteria and return None
        input_taxrank = _prep_taxranks(self.taxranks)
        input_taxmap = _prep_taxmap(self.taxmap2)
        exp = _validate_taxrank_taxmap_taxtree(input_taxrank, input_taxmap,
                                               self.taxtree)
        self.assertEqual(None, exp)

    def test_validate_taxrank_taxmap_taxtree_fail(self):
        # test for missing taxid in tree file
        input_taxrank = _prep_taxranks(self.taxranks)
        input_taxmap = _prep_taxmap(self.taxmap2)
        self.assertRaises(ValueError, _validate_taxrank_taxmap_taxtree,
                          input_taxrank, input_taxmap, self.taxtree2)

    def test_validate_taxrank_taxmap_taxtree_fail2(self):
        # test for case which taxmap file has more taxids
        # than is present in the taxonomy tree.
        input_taxrank = _prep_taxranks(self.taxranks)
        input_taxmap = _prep_taxmap(self.taxmap)
        self.assertRaises(ValueError, _validate_taxrank_taxmap_taxtree,
                          input_taxrank, input_taxmap, self.taxtree)

    def test_parse_silva_taxonomy(self):
        obs_res = parse_silva_taxonomy(self.taxtree, self.taxmap2,
                                       self.taxranks,
                                       include_species_labels=True)
        obs_res.sort_index(inplace=True)
        # expected:
        t1 = ("d__Archaea; p__Aenigmarchaeota; c__Aenigmarchaeia; "
              "o__Aenigmarchaeales; f__Aenigmarchaeales; "
              "g__Candidatus_Aenigmarchaeum; s__uncultured_archaeon")
        exp_res = pd.Series(t1, index=['AB600437.1.1389'])
        exp_res.rename('Taxon', inplace=True)
        exp_res.index.name = 'Feature ID'
        exp_res.sort_index(inplace=True)
        assert_series_equal(obs_res, exp_res)

    def test_parse_silva_taxonomy_no_propagation(self):
        obs_res = parse_silva_taxonomy(self.taxtree, self.taxmap2,
                                       self.taxranks,
                                       include_species_labels=False,
                                       rank_propagation=False)
        obs_res.sort_index(inplace=True)
        # expected:
        t1 = ("d__Archaea; p__Aenigmarchaeota; c__Aenigmarchaeia; "
              "o__Aenigmarchaeales; f__; g__Candidatus_Aenigmarchaeum")
        exp_res = pd.Series(t1, index=['AB600437.1.1389'])
        exp_res.rename('Taxon', inplace=True)
        exp_res.index.name = 'Feature ID'
        exp_res.sort_index(inplace=True)
        assert_series_equal(obs_res, exp_res)

    def test_parse_silva_taxonomy_no_propagation_with_species(self):
        obs_res = parse_silva_taxonomy(self.taxtree, self.taxmap2,
                                       self.taxranks,
                                       include_species_labels=True,
                                       rank_propagation=False)
        obs_res.sort_index(inplace=True)
        # expected:
        t1 = ("d__Archaea; p__Aenigmarchaeota; c__Aenigmarchaeia; "
              "o__Aenigmarchaeales; f__; g__Candidatus_Aenigmarchaeum; "
              "s__uncultured_archaeon")
        exp_res = pd.Series(t1, index=['AB600437.1.1389'])
        exp_res.rename('Taxon', inplace=True)
        exp_res.index.name = 'Feature ID'
        exp_res.sort_index(inplace=True)
        assert_series_equal(obs_res, exp_res)

    def test_parse_silva_taxonomy_no_prop_outoforder_noclass_wsp(self):
        obs_res = parse_silva_taxonomy(self.taxtree, self.taxmap2,
                                       self.taxranks,
                                       include_species_labels=True,
                                       rank_propagation=False,
                                       ranks=['domain', 'phylum', 'genus',
                                              'family', 'order'])
        obs_res.sort_index(inplace=True)
        # expected:
        t1 = ("d__Archaea; p__Aenigmarchaeota; o__Aenigmarchaeales; "
              "f__; g__Candidatus_Aenigmarchaeum; s__uncultured_archaeon")
        exp_res = pd.Series(t1, index=['AB600437.1.1389'])
        exp_res.rename('Taxon', inplace=True)
        exp_res.index.name = 'Feature ID'
        exp_res.sort_index(inplace=True)
        assert_series_equal(obs_res, exp_res)
