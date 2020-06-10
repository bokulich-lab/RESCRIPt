# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import (DNAFASTAFormat, AlignedDNAIterator)


def degap_sequences(aligned_sequences:
                    AlignedDNAIterator) -> DNAFASTAFormat:
    result = DNAFASTAFormat()
    with result.open() as out_fasta:
        for seq in aligned_sequences:
            dg_seq = seq.degap()
            dg_seq.write(out_fasta)
    return result
