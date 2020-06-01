# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._format import (SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat,
                      SILVATaxidMapFormat, SILVATaxidMapDirectoryFormat,
                      RNAFASTAFormat, RNASequencesDirectoryFormat)
from ._type import SILVATaxonomy, SILVATaxidMap, RNASequence

__version__ = '2020.6'

__all__ = ['SILVATaxonomyFormat', 'SILVATaxonomyDirectoryFormat',
           'SILVATaxidMapFormat', 'SILVATaxidMapDirectoryFormat',
           'SILVATaxonomy', 'SILVATaxidMap', 'RNAFASTAFormat',
           'RNASequencesDirectoryFormat', 'RNASequence']
