import random
from typing import Union

from q2_types.feature_data import (AlignedDNAFASTAFormat, DNAFASTAFormat)

from rescript._utilities import _read_dna_fasta, _read_dna_alignment_fasta


def _read_seqs(sequences: Union[DNAFASTAFormat, AlignedDNAFASTAFormat]):
    return _read_dna_fasta(str(sequences)) if isinstance(
        sequences, DNAFASTAFormat) else _read_dna_alignment_fasta(
        str(sequences))


def subsample_fasta(sequences: AlignedDNAFASTAFormat,
                    subsample_size: Union[float, int] = 0.1,
                    random_seed: int = 42) -> \
        AlignedDNAFASTAFormat:
    # set seed for reproducibility
    random.seed(random_seed)

    # get all ids
    seqs_iter = _read_seqs(sequences)
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

    # write out sampled sequences
    result = DNAFASTAFormat() if isinstance(
        sequences, DNAFASTAFormat) else AlignedDNAFASTAFormat()

    seqs_iter = _read_seqs(sequences)
    with result.open() as out_fasta:
        for seq in seqs_iter:
            id = seq.metadata['id']
            if id in ids_sample:
                seq.write(out_fasta)
                ids_sample.remove(id)

    return result
