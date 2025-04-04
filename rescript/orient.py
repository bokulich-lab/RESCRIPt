# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import gzip
import tempfile
import pandas as pd
from itertools import zip_longest
from typing import List, Dict, Any

from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt)

from ._utilities import run_command


def _add_optional_parameters(cmd: List[str], **kwargs: Dict[str, Any]) -> None:
    """
    Add optional parameters to a command.

    Parameters
    ----------
    cmd : List[str]
        The command to add optional parameters to.
    kwargs : Dict[str, Any]
        The optional parameters to add to the command.
    """
    for parameter, value in kwargs.items():
        if isinstance(value, bool) and value is True:
            cmd.append(f"--{parameter}")
        elif isinstance(value, str) and value != "":
            cmd.append(f"--{parameter}")
            cmd.append(f"{value}")


def orient_seqs(
    sequences: DNAFASTAFormat,
    reference_sequences: DNAFASTAFormat = None,
    dbmask: str = None,
    relabel: str = None,
    relabel_keep: bool = None,
    relabel_md5: bool = None,
    relabel_self: bool = None,
    relabel_sha1: bool = None,
    sizein: bool = None,
    sizeout: bool = None,
) -> (DNAFASTAFormat, DNAFASTAFormat):
    oriented, notmatched = DNAFASTAFormat(), DNAFASTAFormat()
    if reference_sequences is not None:
        # use vsearch to orient seqs against reference database
        # note: qmask is disabled as DNAFASTAFormat requires all output
        # seqs to be uppercase. Could loop through output seqs to convert
        # to upper but which is faster: disabling masking or looping
        # through with skbio?
        cmd = [
            'vsearch',
            '--orient', str(sequences),
            '--fastaout', str(oriented),
            '--notmatched', str(notmatched),
            '--db', str(reference_sequences),
            '--qmask', 'none',
        ]

        _add_optional_parameters(
            cmd,
            dbmask=dbmask,
            relabel=relabel,
            relabel_keep=relabel_keep,
            relabel_md5=relabel_md5,
            relabel_self=relabel_self,
            relabel_sha1=relabel_sha1,
            sizein=sizein,
            sizeout=sizeout,
        )

        run_command(cmd)

    else:
        oriented = DNAFASTAFormat()
        with oriented.open() as out_fasta:
            for seq in sequences.view(DNAIterator):
                seq.reverse_complement().write(out_fasta)

    return oriented, notmatched


# note: this function is no longer used, but I would like to keep it
# here, as this could be implemented for RC of FASTA seqs, see #224
def _vsearch_revcomp_fastq(seqs_fp, out_fp):
    cmd = [
        'vsearch',
        '--fastx_revcomp', str(seqs_fp),
        '--fastqout', str(out_fp),
    ]
    run_command(cmd)
    run_command(['gzip', str(out_fp)])


# adapted from q2-demux
def read_fastq(filepath):
    fh = gzip.open(filepath, 'rt')
    # vsearch strips off the lane info, so we will do the same for consistency
    for head, seq, qh, q in zip_longest(*[fh] * 4):
        yield (head.strip().split(' ')[0], seq.strip(), qh.strip(), q.strip())


def orient_reads(
    sequences: CasavaOneEightSingleLanePerSampleDirFmt,
    reference_sequences: DNAFASTAFormat = None,
    dbmask: str = None,
) -> (CasavaOneEightSingleLanePerSampleDirFmt,
      CasavaOneEightSingleLanePerSampleDirFmt):
    oriented = CasavaOneEightSingleLanePerSampleDirFmt()
    notmatched = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest
    # inspect manifest to determine whether reads are joined or paired-end
    # the CasavaOneEightSingleLanePerSampleDirFmt creates a manifest with
    # both forward and reverse columns, but the reverse column is empty if
    # the data are joined.
    paired = manifest.reverse.any()

    for _sample in manifest.itertuples():
        fwd_name = os.path.basename(_sample.forward)
        fwd_path_in = str(sequences.path / fwd_name)
        # vsearch outputs are uncompressed, so clip the extension
        fwd_path_out = os.path.splitext(str(oriented.path / fwd_name))[0]
        fwd_notmatched = os.path.splitext(str(notmatched.path / fwd_name))[0]
        if paired:
            rev_name = os.path.basename(_sample.reverse)
            rev_path_in = str(sequences.path / rev_name)
            # vsearch outputs are uncompressed, clip extension
            rev_path_out = os.path.splitext(str(oriented.path / rev_name))[0]
            rev_notmatched = os.path.splitext(
                str(notmatched.path / rev_name))[0]

        # use vsearch to orient reads against reference database
        # this is only done to detect read orientations, outputs are tossed
        # except in the case of pre-joined reads, in which case vsearch
        # does the job (reverse-complement mis-oriented joined reads)
        df = None
        with tempfile.NamedTemporaryFile() as tabbedout:
            cmd = [
                'vsearch',
                '--orient', str(fwd_path_in),
                '--fastqout', str(fwd_path_out),
                '--notmatched', str(fwd_notmatched),
                '--db', str(reference_sequences),
                '--qmask', 'none',
                '--tabbedout', tabbedout.name,
            ]

            _add_optional_parameters(
                cmd,
                dbmask=dbmask,
            )

            run_command(cmd)
            # run_command(['gzip', str(fwd_path_out)])

            # parse the read orientation report from vsearch
            df = pd.read_csv(
                tabbedout.name, sep='\t', index_col=0, header=None)
            df.columns = ['orientation', 'f_hits', 'r_hits']

        # For joined reads, the outputs above should be fine and we are done.
        # But for PE reads the fun is only beginning...

        if paired:
            # for PE reads we want to throw out the outputs that were created
            # by vsearch, as we do not actually want reverse-complemented reads
            # instead we want to use the tabbedout report created in that step
            # to identify misoriented reads, and swap F/R reads.
            # note that opening files in 'w' mode will intentionally overwrite.
            with open(str(fwd_path_out), 'wb') as f_out:
                with open(str(rev_path_out), 'wb') as r_out:
                    with open(str(fwd_notmatched), 'wb') as f_notmatched:
                        with open(str(rev_notmatched), 'wb') as r_notmatched:
                            # read in paired reads (should be in same order)
                            # one fastq sequence pair at a time
                            for _f, _r in zip(read_fastq(str(fwd_path_in)),
                                              read_fastq(str(rev_path_in))):
                                # grab read ID from the header line
                                read_id = _f[0].lstrip('@')
                                orientation = df.loc[read_id, 'orientation']
                                # if correctly oriented, keep
                                if orientation == '+':
                                    f_out.write(
                                        ('\n'.join(_f) + '\n').encode('utf-8'))
                                    r_out.write(
                                        ('\n'.join(_r) + '\n').encode('utf-8'))
                                # if reverse orientation, swap F/R reads
                                elif orientation == '-':
                                    f_out.write(
                                        ('\n'.join(_r) + '\n').encode('utf-8'))
                                    r_out.write(
                                        ('\n'.join(_f) + '\n').encode('utf-8'))
                                # the alternative is no match so we toss...
                                else:
                                    f_notmatched.write(
                                        ('\n'.join(_f) + '\n').encode('utf-8'))
                                    r_notmatched.write(
                                        ('\n'.join(_r) + '\n').encode('utf-8'))
            # gzip all files for that sample (output format expects .gz)
            run_command(['gzip', str(rev_notmatched)])
            run_command(['gzip', str(rev_path_out)])
        run_command(['gzip', str(fwd_notmatched)])
        run_command(['gzip', str(fwd_path_out)])

    return oriented, notmatched
