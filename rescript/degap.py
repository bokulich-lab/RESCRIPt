# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import (DNAFASTAFormat, AlignedDNAIterator)


def degap_seqs(aligned_sequences:
               AlignedDNAIterator,
               min_length: int = 1) -> DNAFASTAFormat:
    result = DNAFASTAFormat()
    with result.open() as out_fasta:
        for seq in aligned_sequences:
            dg_seq = seq.degap()
            #  If seq is all gaps, then dg_seq will be an empty string
            #  and we'll not write it out.
            if len(dg_seq) >= min_length:
                dg_seq.write(out_fasta)
    return result
