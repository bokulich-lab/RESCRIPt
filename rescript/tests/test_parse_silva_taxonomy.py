# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from qiime2.plugin.testing import TestPluginBase
from rescript.parse_silva_taxonomy import parse_silva_taxonomy
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.tree import Phylogeny, Rooted
import rescript
from rescript.types._format import (
    SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat, SILVATaxidMapFormat,
    SILVATaxidMapDirectoryFormat)
from rescript.types._type import SILVATaxonomy, SILVATaxidMap



# class TestParseSilvaTaxonomy(TestPluginBase):
#     package = 'rescript.tests'
#
#     def setUp(self):
#         super().setUp()
#         input_fp = self.get_data_path('cleanseq-test-1.fasta')
#         self.taxonomy_str =
#
#
#     def test_keep_allowed_chars(self):
#         obs
#
