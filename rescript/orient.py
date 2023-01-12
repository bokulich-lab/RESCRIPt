# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from warnings import warn

from ._utilities import run_command


def _warn_deprecated(value, default, name):
    if value != default:
        warn(
            f"{name} is not available in this current version.",
            DeprecationWarning,
        )


def orient_seqs(sequences: DNAFASTAFormat,
                reference_sequences: DNAFASTAFormat = None,
                perc_identity: float = 0.9,
                query_cov: float = 0.9,
                threads: int = 1,
                left_justify: bool = False,
                ) -> (DNAFASTAFormat, DNAFASTAFormat):
    oriented, notmatched = DNAFASTAFormat(), DNAFASTAFormat()
    if reference_sequences is not None:
        # use vsearch to orient seqs against reference database
        # note: qmask is disabled as DNAFASTAFormat requires all output
        # seqs to be uppercase. Could loop through output seqs to convert
        # to upper but which is faster: disabling masking or looping
        # through with skbio?
        # NOTE: --id, --query_cov and --leftjust are disabled since
        # they are not available with orient
        cmd = [
            'vsearch',
            '--orient', str(sequences),
            '--fastaout', str(oriented),
            '--notmatched', str(notmatched),
            '--db', str(reference_sequences),
            # '--id', str(perc_identity),
            '--qmask', 'none',
            # '--query_cov', str(query_cov),
            '--threads', str(threads),
        ]
        # if left_justify:
        #     cmd.append('--leftjust')
        run_command(cmd)

        # Warn about parameters were deprecated
        _warn_deprecated(perc_identity, 0.9, 'perc_identity')
        _warn_deprecated(query_cov, 0.9, 'query_cov')
        _warn_deprecated(left_justify, False, 'left_justify')

    else:
        oriented = DNAFASTAFormat()
        with oriented.open() as out_fasta:
            for seq in sequences.view(DNAIterator):
                seq.reverse_complement().write(out_fasta)

    return oriented, notmatched
