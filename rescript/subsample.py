# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import random
from q2_types.feature_data import DNAIterator


def subsample_fasta(sequences: DNAIterator,
                    subsample_size: float = 0.1,
                    random_seed: int = 1) -> DNAIterator:
    # set seed for reproducibility
    random.seed(random_seed)

    # create new iterator that will randomly yield sequences
    def sample_from_iterator():
        for seq in sequences:
            if random.random() < subsample_size:
                yield seq

    return DNAIterator(sample_from_iterator())
