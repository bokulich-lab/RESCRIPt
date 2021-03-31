# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import DNAFASTAFormat, RNAFASTAFormat

from rescript._utilities import _rna_to_dna


# This exposes the transformer as its own method
def reverse_transcribe(rna_sequences: RNAFASTAFormat) -> DNAFASTAFormat:
    return _rna_to_dna(str(rna_sequences))
