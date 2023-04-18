# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from warnings import warn

from ._utilities import run_command


def _warn_deprecated(value, default, name):
    if value != default:
        warn(
            f"{name} is deprecated and will be removed from RESCRIPt 2023.5 "
            "(RESCRIPt 2023.2 will be the last release with this parameter).",
            FutureWarning
        )


def orient_seqs(sequences: DNAFASTAFormat,
                reference_sequences: DNAFASTAFormat = None,
                perc_identity: float = 0.9,
                query_cov: float = 0.9,
                threads: int = 1,
                left_justify: bool = False,
                ) -> (DNAFASTAFormat, DNAFASTAFormat):
    matched_temp, notmatched = DNAFASTAFormat(), DNAFASTAFormat()
    if reference_sequences is not None:
        # use vsearch to search query seqs against reference database
        # report orientation of query seqs relative to reference seqs.
        with tempfile.NamedTemporaryFile() as out:
            # note: qmask is disabled as DNAFASTAFormat requires all output
            # seqs to be uppercase. Could loop through output seqs to convert
            # to upper but which is faster: disabling masking or looping
            # through with skbio?
            cmd = ['vsearch',
                   '--usearch_global', str(sequences),
                   '--matched', str(matched_temp),
                   '--notmatched', str(notmatched),
                   '--db', str(reference_sequences),
                   '--id', str(perc_identity),
                   '--maxaccepts', '1',
                   '--strand', 'both',
                   '--qmask', 'none',
                   '--query_cov', str(query_cov),
                   '--threads', str(threads),
                   '--userfields', 'qstrand',
                   '--userout', out.name]
            if left_justify:
                cmd.append('--leftjust')
            run_command(cmd)
            with open(out.name, 'r') as orient:
                orientations = [line.strip() for line in orient]

            # Warn about parameters that will be deprecated
            _warn_deprecated(perc_identity, 0.9, 'perc_identity')
            _warn_deprecated(query_cov, 0.9, 'query_cov')
            _warn_deprecated(left_justify, False, 'left_justify')

        # if any query seqs are in reverse orientation, reverse complement
        if '-' in orientations:
            matched = DNAFASTAFormat()
            with matched.open() as out_fasta:
                for seq, orientation in zip(
                        matched_temp.view(DNAIterator), orientations):
                    if orientation == '+':
                        seq.write(out_fasta)
                    elif orientation == '-':
                        seq.reverse_complement().write(out_fasta)
        else:
            matched = matched_temp
    else:
        matched = DNAFASTAFormat()
        with matched.open() as out_fasta:
            for seq in sequences.view(DNAIterator):
                seq.reverse_complement().write(out_fasta)

    return matched, notmatched
