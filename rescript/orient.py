# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import tempfile
from multiprocessing import Pool, Manager, Process, cpu_count

from joblib import Parallel, delayed
from q2_types.feature_data import DNAFASTAFormat, DNAIterator

from ._utilities import run_command


def worker(seq_obj, orientation):
    '''Reverse complement a sequence if necessary and place the result in a queue'''
    if orientation == '+':
        return seq_obj
    elif orientation == '-':
        return seq_obj.reverse_complement()


def writer(out_result, q):
    '''Get a result from the queue and write it out'''
    with out_result.open() as out:
        while True:
            seq_obj = q.get()
            if seq_obj is None:
                break
            seq_obj.write(out)


def orient_seqs(sequences: DNAFASTAFormat,
                reference_sequences: DNAFASTAFormat,
                perc_identity: float = 0.9,
                query_cov: float = 0.9,
                threads: int = 1,
                left_justify: bool = False,
                ) -> (DNAFASTAFormat, DNAFASTAFormat):
    matched_temp, notmatched = DNAFASTAFormat(), DNAFASTAFormat()
    # use vsearch to search query seqs against reference database
    # report orientation of query seqs relative to reference seqs.
    with tempfile.NamedTemporaryFile() as out:
        # note: qmask is disabled as DNAFASTAFormat requires all output seqs
        # to be uppercase. Could loop through output seqs to convert to upper
        # but which is faster: disabling masking or looping through with skbio?
        cmd = ['vsearch', '--usearch_global', str(sequences),
               '--matched', str(matched_temp), '--notmatched', str(notmatched),
               '--db', str(reference_sequences), '--id', str(perc_identity),
               '--maxaccepts', '1', '--strand', 'both', '--qmask', 'none',
               '--query_cov', str(query_cov), '--threads', str(threads),
               '--userfields', 'qstrand', '--userout', out.name]
        if left_justify:
            cmd.append('--leftjust')
        run_command(cmd)
        with open(out.name, 'r') as orient:
            orientations = [line.strip() for line in orient]

    # if any query seqs are in reverse orientation, reverse complement
    if '-' in orientations:
        matched = DNAFASTAFormat()

        parallel = Parallel(n_jobs=4, backend='loky')
        seqs = parallel(delayed(worker)(seq, orientation) for seq, orientation in zip(
                     matched_temp.view(DNAIterator), orientations))

        with matched.open() as out:
            for seq in seqs:
                seq.write(out)

    else:
        matched = matched_temp

    return matched, notmatched
