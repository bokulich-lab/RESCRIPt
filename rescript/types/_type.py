# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType
from q2_types.feature_data import FeatureData


SILVATaxonomy = SemanticType(
    'SILVATaxonomy', variant_of=FeatureData.field['type'])
SILVATaxidMap = SemanticType(
    'SILVATaxidMap', variant_of=FeatureData.field['type'])
