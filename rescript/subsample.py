# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import random
from itertools import tee
from typing import Union
from q2_types.feature_data import DNAIterator


def subsample_fasta(sequences: DNAIterator,
                    subsample_size: Union[float, int] = 0.1,
                    random_seed: int = 42) -> DNAIterator:
    # set seed for reproducibility
    random.seed(random_seed)

    # clone the iterator
    seqs_iter, seqs_iter_clone = tee(sequences)

    # get all ids
    ids = [seq.metadata['id'] for seq in seqs_iter]

    # get a sample of the ids
    if isinstance(subsample_size, float):
        ids_sample = set(random.sample(ids, int(subsample_size * len(ids))))
    else:
        if subsample_size <= len(ids):
            ids_sample = set(random.sample(ids, subsample_size))
        else:
            raise ValueError(f"The requested sample size {subsample_size} "
                             f"is bigger than the total number of "
                             f"sequences {len(ids)}")

    # create new iterator that will yield sequences
    # based on the list generated above
    def sample_from_iterator():
        for seq in seqs_iter_clone:
            id = seq.metadata['id']
            if id in ids_sample:
                ids_sample.remove(id)
                yield seq

    return DNAIterator(sample_from_iterator())
