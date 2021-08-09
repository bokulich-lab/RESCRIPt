# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2 import Metadata
from q2_types.feature_data import TSVTaxonomyFormat
import pandas as pd

from rescript.edit_taxonomy import edit_taxonomy

import_data = qiime2.Artifact.import_data


class TestReplaceTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        tax_fp = self.get_data_path('escherichia_shigella_taxonomy.txt')
        self.taxonomy = TSVTaxonomyFormat(tax_fp, mode='r').view(pd.Series)

    def test_edit_taxonomy_pass(self):
        replc = self.get_data_path('taxonomy-replacement-pass.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=md_replc_col,
                                   use_regex=False)
        obs_dict = obs_series.to_dict()

        exp_dict = {'Sal01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'Sal02': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'UncultSal': ('d__SUPER_DUPER_BACTERIA; '
                                  'p__Proteobacteria; '
                                  'c__Gammaproteobacteria; '
                                  'o__Enterobacterales; '
                                  'f__Enterobacteriaceae; '
                                  'g__Escherichia-Shigella; '
                                  's__uncultured_Salmonella'),
                    'Esch01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__'),
                    'Shig01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__')}

        self.maxDiff = None
        self.assertEqual(obs_dict, exp_dict)

    def test_edit_taxonomy_pass_full_string(self):
        replc = self.get_data_path('taxonomy-replacement-full-strings.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=md_replc_col,
                                   use_regex=False)
        obs_dict = obs_series.to_dict()

        exp_dict = {'Sal01': ('d__Bacteria; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'Sal02': ('d__Bacteria; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'UncultSal': ('d__Bacteria; '
                                  'p__Proteobacteria; '
                                  'c__Gammaproteobacteria; '
                                  'o__Enterobacterales; '
                                  'f__Enterobacteriaceae; '
                                  'g__Escherichia-Shigella; '
                                  's__uncultured_Salmonella'),
                    'Esch01': ('d__Bacteria; p__Proteobacteria; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; '
                               's__unknwon_Escherichia-Shigella'),
                    'Shig01': ('d__Bacteria; p__Proteobacteria; '
                               'c__Gammaproteobacteria; '
                               'o__Enterobacterales; f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; '
                               's__unknwon_Escherichia-Shigella')}

        self.maxDiff = None
        self.assertEqual(obs_dict, exp_dict)

    def test_edit_taxonomy_pass_regex(self):
        replc = self.get_data_path('taxonomy-replacement-regex.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=md_replc_col,
                                   use_regex=True)
        obs_dict = obs_series.to_dict()

        exp_dict = {'Sal01': ('d__Bacteria; p__LAME-PYHLA; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'Sal02': ('d__Bacteria; p__LAME-PYHLA; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'UncultSal': ('d__Bacteria; '
                                  'p__LAME-PYHLA; '
                                  'c__Gammaproteobacteria; '
                                  'o__Enterobacterales; '
                                  'f__Enterobacteriaceae; '
                                  'g__Escherichia-Shigella; '
                                  's__UNCIVILIZED_Salmonella'),
                    'Esch01': ('d__Bacteria; p__LAME-PYHLA; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__'),
                    'Shig01': ('d__Bacteria; p__LAME-PYHLA; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__')}

        self.maxDiff = None
        self.assertEqual(obs_dict, exp_dict)

    def test_edit_taxonomy_pass_regex_cl(self):
        ss = ['p__.*; c__', 's__uncultured', 'g__Escherichia.*; s__',
              r'g__\[Salmonella\]', 'g__Shigella']
        rs = ['p__LAME-PYHLA; c__', 's__UNCIVILIZED',
              'g__Escherichia-Shigella; s__', 'g__Escherichia-Shigella',
              'g__Escherichia-Shigella']
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   search_strings=ss,
                                   replacement_strings=rs,
                                   use_regex=True)
        obs_dict = obs_series.to_dict()

        exp_dict = {'Sal01': ('d__Bacteria; p__LAME-PYHLA; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'Sal02': ('d__Bacteria; p__LAME-PYHLA; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'UncultSal': ('d__Bacteria; '
                                  'p__LAME-PYHLA; '
                                  'c__Gammaproteobacteria; '
                                  'o__Enterobacterales; '
                                  'f__Enterobacteriaceae; '
                                  'g__Escherichia-Shigella; '
                                  's__UNCIVILIZED_Salmonella'),
                    'Esch01': ('d__Bacteria; p__LAME-PYHLA; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__'),
                    'Shig01': ('d__Bacteria; p__LAME-PYHLA; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__')}

        self.maxDiff = None
        self.assertEqual(obs_dict, exp_dict)

    def test_edit_taxonomy_pass_regex_map_and_cl(self):
        # tests that patterns are handled identically wheater or not they
        # are supplied by the command line or mapping file.
        # will work if using r'<regex>' or `\\` to esccape characters.
        ss = ['p__.*; c__', 's__uncultured', 'g__Escherichia.*; s__',
              r'g__\[Salmonella\]', 'g__Shigella']
        rs = ['p__LAME-PYHLA; c__', 's__UNCIVILIZED',
              'g__Escherichia-Shigella; s__', 'g__Escherichia-Shigella',
              'g__Escherichia-Shigella']
        obs_series_cl = edit_taxonomy(taxonomy=self.taxonomy,
                                      search_strings=ss,
                                      replacement_strings=rs,
                                      use_regex=True)
        obs_dict_cl = obs_series_cl.to_dict()

        replc = self.get_data_path('taxonomy-replacement-regex.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series_map = edit_taxonomy(taxonomy=self.taxonomy,
                                       replacement_map=md_replc_col,
                                       use_regex=True)
        obs_dict_map = obs_series_map.to_dict()

        exp_dict = {'Sal01': ('d__Bacteria; p__LAME-PYHLA; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'Sal02': ('d__Bacteria; p__LAME-PYHLA; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'UncultSal': ('d__Bacteria; '
                                  'p__LAME-PYHLA; '
                                  'c__Gammaproteobacteria; '
                                  'o__Enterobacterales; '
                                  'f__Enterobacteriaceae; '
                                  'g__Escherichia-Shigella; '
                                  's__UNCIVILIZED_Salmonella'),
                    'Esch01': ('d__Bacteria; p__LAME-PYHLA; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__'),
                    'Shig01': ('d__Bacteria; p__LAME-PYHLA; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__')}

        self.maxDiff = None
        self.assertEqual(obs_dict_cl, exp_dict)
        self.assertEqual(obs_dict_map, exp_dict)

    def test_edit_taxonomy_pass_cl(self):
        ss = ['g__[Salmonella]', 'g__Escherichia;', 'g__Shigella',
              'd__Bacteria;']
        rs = ['g__Escherichia-Shigella', 'g__Escherichia-Shigella;',
              'g__Escherichia-Shigella', 'd__SUPER_DUPER_BACTERIA;']
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   search_strings=ss,
                                   replacement_strings=rs,
                                   use_regex=False)
        obs_dict = obs_series.to_dict()

        exp_dict = {'Sal01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'Sal02': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'UncultSal': ('d__SUPER_DUPER_BACTERIA; '
                                  'p__Proteobacteria; '
                                  'c__Gammaproteobacteria; '
                                  'o__Enterobacterales; '
                                  'f__Enterobacteriaceae; '
                                  'g__Escherichia-Shigella; '
                                  's__uncultured_Salmonella'),
                    'Esch01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__'),
                    'Shig01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__')}

        self.maxDiff = None
        self.assertEqual(obs_dict, exp_dict)

    def test_edit_taxonomy_unequal_string_lists_fail_cl(self):
        ss = ['g__[Salmonella]', 'g__Escherichia;', 'g__Shigella']
        rs = ['g__Escherichia-Shigella', 'g__Escherichia-Shigella;',
              'g__Escherichia-Shigella', 'd__SUPER_DUPER_BACTERIA;']
        self.assertRaises(ValueError, edit_taxonomy, taxonomy=self.taxonomy,
                          search_strings=ss, replacement_strings=rs)

    def test_cl_string_precedence_over_map_file(self):
        ss = ['g__[Salmonella]', 'g__Escherichia;', 'g__Shigella',
              'd__Bacteria;']
        rs = ['g__Escherichia-Shigella', 'g__Escherichia-Shigella;',
              'g__Escherichia-Shigella', 'd__SUPER_DUPER_BACTERIA;']
        replc = self.get_data_path('taxonomy-replacement-full-strings.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=md_replc_col,
                                   search_strings=ss,
                                   replacement_strings=rs,
                                   use_regex=False)
        obs_dict = obs_series.to_dict()

        exp_dict = {'Sal01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'Sal02': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                              'c__Gammaproteobacteria; o__Enterobacterales; '
                              'f__Enterobacteriaceae; '
                              'g__Escherichia-Shigella; '
                              's__Salmonella_enterica'),
                    'UncultSal': ('d__SUPER_DUPER_BACTERIA; '
                                  'p__Proteobacteria; '
                                  'c__Gammaproteobacteria; '
                                  'o__Enterobacterales; '
                                  'f__Enterobacteriaceae; '
                                  'g__Escherichia-Shigella; '
                                  's__uncultured_Salmonella'),
                    'Esch01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__'),
                    'Shig01': ('d__SUPER_DUPER_BACTERIA; p__Proteobacteria; '
                               'c__Gammaproteobacteria; o__Enterobacterales; '
                               'f__Enterobacteriaceae; '
                               'g__Escherichia-Shigella; s__')}

        self.maxDiff = None
        self.assertEqual(obs_dict, exp_dict)
