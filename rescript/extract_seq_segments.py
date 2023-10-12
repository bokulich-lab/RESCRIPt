# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import DNAFASTAFormat

from ._utilities import run_command


def extract_seq_segments(input_sequences: DNAFASTAFormat,
                         reference_segment_sequences: DNAFASTAFormat,
                         perc_identity: float = 0.7,
                         min_seq_len: int = None,
                         threads: int = 1
                         ) -> (DNAFASTAFormat, DNAFASTAFormat):
    # Note: the `--db` is actually comprised of the short amplicon
    # sequences. That is, they are the reference db. The file fed into
    # `--usearch_global` would be the "full-length" sequences, e.g.
    # sequences pulled from GenBank or SILVA, etc... These are actually
    # the 'query' sequences that have there regions extracted.
    # See here for more details:
    # https://www.drive5.com/usearch/manual8.1/utax_user_train.html
    #
    # Warning: I hard set '--strand both', as multiple segments can be
    # extracted from the same sequence, resulting in multiple output
    # seqeunces to have the same ID, causing this action to fail.,
    #
    # WARNING: I had to add `--qmask none` to this command to prevent
    # the output from having mixed-case characters. That is, by
    # default, the `dust` algorithm is run, I think due to the
    # `usearch_global` command.
    extr_seq_segs = DNAFASTAFormat()
    no_matches = DNAFASTAFormat()
    with extr_seq_segs.open() as extr_f, no_matches.open() as no_match_f:
        cmd = ['vsearch',
               '--usearch_global', str(input_sequences),
               '--db', str(reference_segment_sequences),
               '--id', str(perc_identity),
               '--strand', 'plus',
               '--threads', str(threads),
               '--qmask', 'none',
               '--qsegout', extr_f.name,
               '--notmatched', no_match_f.name]
        if min_seq_len:
            cmd += ['--minseqlength', '%i' % min_seq_len]
        run_command(cmd)
    return (extr_seq_segs, no_matches)
