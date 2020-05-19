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
                                           _prep_taxmap)
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.tree import Phylogeny, Rooted
from skbio.tree import TreeNode
import rescript
from rescript.types._format import (
    SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat, SILVATaxidMapFormat,
    SILVATaxidMapDirectoryFormat)
from rescript.types._type import SILVATaxonomy, SILVATaxidMap
from io import StringIO
import pandas as pd
from pandas.testing import assert_frame_equal


class TestParseSilvaTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.bad_char_str = '~`!@#$%^&*;:\"?<>'
        self.good_char_str = '1 am a G00d str1ng.'
        self.taxranks = self.get_data_path('tax_slv_ssu_test.txt')
        self.taxmap = self.get_data_path('taxmap_slv_ssu_ref_nr_test.txt')


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
        dd = { 'taxid': ['2', '11084', '42913', '42914', '42915',
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
