# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from joblib import Parallel, delayed
from q2_types.feature_data import (DNAFASTAFormat, AlignedDNAIterator)


def worker(seq_obj):
    '''Degap a single sequence and place the result in a queue'''
    dg_seq = seq_obj.degap()
    return dg_seq


def degap_seqs(aligned_sequences:
               AlignedDNAIterator,
               min_length: int = 1) -> DNAFASTAFormat:
    result = DNAFASTAFormat()

    parallel = Parallel(n_jobs=4, backend='loky')
    seqs = parallel(delayed(worker)(seq) for seq in aligned_sequences)

    with result.open() as out:
        for seq in seqs:
            if len(seq) >= min_length:
                seq.write(out)


    return result
