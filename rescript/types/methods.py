# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import RNAIterator, DNAIterator
from rescript._utilities import _rna_to_dna_iterator


# This exposes the transformer as its own method
def reverse_transcribe(rna_sequences: RNAIterator) -> DNAIterator:
    generator = _rna_to_dna_iterator(rna_sequences)
    return DNAIterator(generator)
