# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
from multiprocessing import cpu_count

from joblib import Parallel, delayed
from q2_types.feature_data import DNAFASTAFormat, DNAIterator


def _filt_seq_with_degenerates(seq, num_degenerates):
    degenerates_in_seq = sum(seq.degenerates())
    return degenerates_in_seq >= num_degenerates


def _filter_homopolymer(seq, homopolymer_length):
    nhl = homopolymer_length - 1  # due to how regex is written
    regex_str = "([ACGTURYSWKMBDHVN])\\1{%s,}" % nhl
    homopolymers = [p for p in re.finditer(regex_str, str(seq))]
    return any(homopolymers)


def _filter_seq(seq_obj, num_degenerates, homopolymer_length):
    '''Process a single sequence'''
    degen = _filt_seq_with_degenerates(seq_obj, num_degenerates)
    if not degen:
        poly = _filter_homopolymer(seq_obj, homopolymer_length)
        if not poly:  # if we make it here, write seq to file
            return seq_obj


def cull_seqs(sequences: DNAIterator, num_degenerates: int = 5,
              homopolymer_length: int = 8) -> DNAFASTAFormat:
    result = DNAFASTAFormat()

    n_jobs = cpu_count()
    parallel = Parallel(n_jobs=n_jobs, backend='loky')
    seqs = parallel(
        delayed(_filter_seq)(
            seq,
            num_degenerates,
            homopolymer_length) for seq in sequences)

    with result.open() as out:
        for seq in seqs:
            seq.write(out) if seq else False

    return result
