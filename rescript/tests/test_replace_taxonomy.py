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

from rescript.replace_taxonomy import replace_taxonomy

import_data = qiime2.Artifact.import_data


class TestReplaceTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        tax_fp = self.get_data_path('escherichia_shigella_taxonomy.txt')
        self.taxonomy = TSVTaxonomyFormat(tax_fp, mode='r').view(pd.Series)

    def test_replace_taxonomy_pass(self):
        replc = self.get_data_path('taxonomy-replacement-pass.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series = replace_taxonomy(self.taxonomy, md_replc_col, True)
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

    def test_replace_taxonomy_pass_full_string(self):
        replc = self.get_data_path('taxonomy-replacement-full-strings.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series = replace_taxonomy(self.taxonomy, md_replc_col, False)
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

    def test_replace_taxonomy_pass_regex(self):
        replc = self.get_data_path('taxonomy-replacement-regex.txt')
        md_replc = Metadata.load(replc)
        md_replc_col = md_replc.get_column('replacements')
        obs_series = replace_taxonomy(self.taxonomy, md_replc_col, True)
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
