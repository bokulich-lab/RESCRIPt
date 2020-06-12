# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
from q2_types.feature_data import DNAFASTAFormat, DNAIterator


def _filt_seq_with_degenerates(seq, num_degenerates):
    degenerates_in_seq = sum(seq.degenerates())
    return degenerates_in_seq >= num_degenerates


def _filter_homopolymer(seq, homopolymer_length):
    nhl = homopolymer_length - 1  # due to how regex is written
    regex_str = "([ACGTURYSWKMBDHVN])\\1{%s,}" % nhl
    homopolymers = [p for p in re.finditer(regex_str, str(seq))]
    return any(homopolymers)


def cull_seqs(sequences: DNAIterator, num_degenerates: int = 5,
              homopolymer_length: int = 8) -> DNAFASTAFormat:
    result = DNAFASTAFormat()
    with result.open() as out_fasta:
        for seq in sequences:
            degen = _filt_seq_with_degenerates(seq, num_degenerates)
            if not degen:
                poly = _filter_homopolymer(seq, homopolymer_length)
                if not poly:  # if we make it here, write seq to file
                    seq.write(out_fasta)
    return result
