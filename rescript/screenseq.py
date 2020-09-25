# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re

from joblib import Parallel, delayed
from q2_types.feature_data import DNAFASTAFormat, DNAIterator


def worker(seq_obj, num_degenerates, homopolymer_length):
    '''Process a single sequence and place it in a queue'''
    degen = _filt_seq_with_degenerates(seq_obj, num_degenerates)
    if not degen:
        poly = _filter_homopolymer(seq_obj, homopolymer_length)
        if not poly:  # if we make it here, write seq to file
            return seq_obj


def writer(out_result, q):
    '''Get a result from the queue and write it out'''
    with out_result.open() as out:
        while True:
            seq_obj = q.get()
            if seq_obj is None:
                break
            seq_obj.write(out)


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

    parallel = Parallel(n_jobs=4, backend='loky')
    seqs = parallel(delayed(worker)(seq, num_degenerates, homopolymer_length) for seq in sequences)

    with result.open() as out:
        for seq in seqs:
            if seq:
                seq.write(out)

    return result
