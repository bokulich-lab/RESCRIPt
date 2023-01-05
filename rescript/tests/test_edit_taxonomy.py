# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
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

        # setup taxonomy to be edited
        tax_fp = self.get_data_path('escherichia_shigella_taxonomy.txt')
        self.taxonomy = TSVTaxonomyFormat(tax_fp, mode='r').view(pd.Series)

        # setup full string replacement
        replc = self.get_data_path('taxonomy-replacement-full-strings.txt')
        md_replc = Metadata.load(replc)
        self.md_replc_col = md_replc.get_column('replacements')

        # setup substring replacement
        ss_replc = self.get_data_path('taxonomy-replacement-pass.txt')
        md_ss_replc = Metadata.load(ss_replc)
        self.md_ss_replc_col = md_ss_replc.get_column('replacements')

        # setup substring regex replacement
        ssr_replc = self.get_data_path('taxonomy-replacement-regex.txt')
        md_ssr_replc = Metadata.load(ssr_replc)
        self.md_ssr_replc_col = md_ssr_replc.get_column('replacements')

        # setup reusable dicts
        self.exp_dict_00 = {'Sal01': ('d__SUPER_DUPER_BACTERIA; '
                                      'p__Proteobacteria; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
                                      'f__Enterobacteriaceae; '
                                      'g__Escherichia-Shigella; '
                                      's__Salmonella_enterica'),
                            'Sal02': ('d__SUPER_DUPER_BACTERIA; '
                                      'p__Proteobacteria; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
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
                            'Esch01': ('d__SUPER_DUPER_BACTERIA; '
                                       'p__Proteobacteria; '
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; s__'),
                            'Shig01': ('d__SUPER_DUPER_BACTERIA; '
                                       'p__Proteobacteria; '
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; s__')}

        self.exp_dict_01 = {'Sal01': ('d__Bacteria; p__LAME-PYHLA; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
                                      'f__Enterobacteriaceae; '
                                      'g__Escherichia-Shigella; '
                                      's__Salmonella_enterica'),
                            'Sal02': ('d__Bacteria; p__LAME-PYHLA; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
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
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; s__'),
                            'Shig01': ('d__Bacteria; p__LAME-PYHLA; '
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; s__')}
        self.exp_dict_02 = {'Sal01': ('d__Bacteria; p__Proteobacteria; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
                                      'f__Enterobacteriaceae; '
                                      'g__Escherichia-Shigella; '
                                      's__Salmonella_enterica'),
                            'Sal02': ('d__Bacteria; p__Proteobacteria; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
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
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; '
                                       's__unknown_Escherichia-Shigella'),
                            'Shig01': ('d__Bacteria; p__Proteobacteria; '
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; '
                                       's__unknown_Escherichia-Shigella')}
        self.exp_dict_03 = {'Sal01': ('d__Bacteria; p__Proteobacteria; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
                                      'f__Enterobacteriaceae; '
                                      'g__Escherichia-Shigella; '
                                      's__Salmonella_enterica'),
                            'Sal02': ('d__Bacteria; p__Proteobacteria; '
                                      'c__Gammaproteobacteria; '
                                      'o__Enterobacterales; '
                                      'f__Enterobacteriaceae; '
                                      'g__Escherichia-Shigella; '
                                      's__Salmonella_enterica'),
                            'UncultSal': ('d__Bacteria; '
                                          'p__Proteobacteria; '
                                          'c__Gammaproteobacteria; '
                                          'o__Enterobacterales; '
                                          'f__Enterobacteriaceae; '
                                          'g__Escherichia-Shigella; '
                                          's__UNCIVIL_Salmonella'),
                            'Esch01': ('d__Bacteria; p__Proteobacteria; '
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; '
                                       's__unknown_Escherichia-Shigella'),
                            'Shig01': ('d__Bacteria; p__Proteobacteria; '
                                       'c__Gammaproteobacteria; '
                                       'o__Enterobacterales; '
                                       'f__Enterobacteriaceae; '
                                       'g__Escherichia-Shigella; '
                                       's__unknown_Escherichia-Shigella')}

    def test_edit_taxonomy_pass(self):
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=self.md_ss_replc_col,
                                   use_regex=False)
        obs_dict = obs_series.to_dict()
        self.maxDiff = None
        self.assertEqual(obs_dict, self.exp_dict_00)

    def test_edit_taxonomy_pass_full_string(self):
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=self.md_replc_col,
                                   use_regex=False)
        obs_dict = obs_series.to_dict()
        self.maxDiff = None
        self.assertEqual(obs_dict, self.exp_dict_02)

    def test_edit_taxonomy_pass_regex(self):
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=self.md_ssr_replc_col,
                                   use_regex=True)
        obs_dict = obs_series.to_dict()
        self.maxDiff = None
        self.assertEqual(obs_dict, self.exp_dict_01)

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
        self.maxDiff = None
        self.assertEqual(obs_dict, self.exp_dict_01)

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
        self.maxDiff = None
        self.assertEqual(obs_dict, self.exp_dict_00)

    def test_edit_taxonomy_unequal_string_lists_fail_cl(self):
        ss = ['g__[Salmonella]', 'g__Escherichia;', 'g__Shigella']
        rs = ['g__Escherichia-Shigella', 'g__Escherichia-Shigella;',
              'g__Escherichia-Shigella', 'd__SUPER_DUPER_BACTERIA;']
        self.assertRaises(ValueError, edit_taxonomy, taxonomy=self.taxonomy,
                          search_strings=ss, replacement_strings=rs)

    def test_merge_cl_and_map_file(self):
        ss = ['s__uncultured_Salmonella']
        rs = ['s__UNCIVIL_Salmonella']
        obs_series = edit_taxonomy(taxonomy=self.taxonomy,
                                   replacement_map=self.md_replc_col,
                                   search_strings=ss,
                                   replacement_strings=rs,
                                   use_regex=False)
        obs_dict = obs_series.to_dict()
        self.maxDiff = None
        self.assertEqual(obs_dict, self.exp_dict_03)

    def test_unequal_ranks(self):
        with self.assertWarnsRegex(UserWarning, "Warning: the number of "
                                   "taxonomy ranks in the output "):
            ss = ['d__Bacter.*__']
            rs = ['d__BAD_REPLACEMENT']
            edit_taxonomy(taxonomy=self.taxonomy,
                          search_strings=ss,
                          replacement_strings=rs,
                          use_regex=True)
