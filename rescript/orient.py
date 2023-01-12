# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import DNAFASTAFormat, DNAIterator

from ._utilities import run_command


def orient_seqs(
    sequences: DNAFASTAFormat,
    reference_sequences: DNAFASTAFormat = None,
    threads: int = 1,
) -> (DNAFASTAFormat, DNAFASTAFormat):
    oriented, notmatched = DNAFASTAFormat(), DNAFASTAFormat()
    if reference_sequences is not None:
        # use vsearch to orient seqs against reference database
        # note: qmask is disabled as DNAFASTAFormat requires all output
        # seqs to be uppercase. Could loop through output seqs to convert
        # to upper but which is faster: disabling masking or looping
        # through with skbio?
        # TODO: add new parameters for orient
        cmd = [
            'vsearch',
            '--orient', str(sequences),
            '--fastaout', str(oriented),
            '--notmatched', str(notmatched),
            '--db', str(reference_sequences),
            '--qmask', 'none',
            '--threads', str(threads),
        ]
        run_command(cmd)

    else:
        oriented = DNAFASTAFormat()
        with oriented.open() as out_fasta:
            for seq in sequences.view(DNAIterator):
                seq.reverse_complement().write(out_fasta)

    return oriented, notmatched
