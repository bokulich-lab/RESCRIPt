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
from skbio import DNA

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


def _vsearch_revcomp_fastq(seqs_fp, out_fp):
    cmd = [
        'vsearch',
        '--fastx_revcomp', str(seqs_fp),
        '--fastqout', str(out_fp),
    ]
    run_command(cmd)
    run_command(['gzip', str(out_fp)])


def read_fastq(filepath):
    # adapted from q2-demux
    fh = gzip.open(filepath, 'rt')
    for header, seq, qh, qual in zip_longest(*[fh] * 4):
        yield (header.strip(), seq.strip(), qh.strip(), qual.strip())


def orient_reads(
    sequences: CasavaOneEightSingleLanePerSampleDirFmt,
    reference_sequences: DNAFASTAFormat = None,
    dbmask: str = None,
) -> (CasavaOneEightSingleLanePerSampleDirFmt,
      CasavaOneEightSingleLanePerSampleDirFmt):
    oriented = CasavaOneEightSingleLanePerSampleDirFmt()
    notmatched = CasavaOneEightSingleLanePerSampleDirFmt()
    manifest = sequences.manifest
    # inspect manifest to determine whether reads are single or paired-end
    # the CasavaOneEightSingleLanePerSampleDirFmt creates a manifest with
    # both forward and reverse columns, but the reverse column is empty if
    # the data are SE.
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

        if reference_sequences is not None:
            # use vsearch to orient reads against reference database
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
                run_command(['gzip', str(fwd_path_out)])

                if paired:
                    # use the tabbedout orientation report to re-orient reverse
                    df = pd.read_csv(
                        tabbedout.name, sep='\t', index_col=0, header=None)
                    df.columns = ['orientation', 'f_hits', 'r_hits']
                    with open(str(rev_path_out), 'wb') as rev_out:
                        with open(str(rev_notmatched), 'wb') as out_notmatched:
                            for _read in read_fastq(str(rev_path_in)):
                                # vsearch strips off the lane info, so reverse
                                # reads should match this pattern
                                _read = list(_read)
                                _read[0] = _read[0].split(' ')[0]
                                head, _seq, qh, q = _read
                                read_id = head.lstrip('@')
                                orientation = df.loc[read_id, 'orientation']
                                if orientation == '+':
                                    rev_out.write(
                                        ('\n'.join(_read) + '\n').encode(
                                            'utf-8'))
                                elif orientation == '-':
                                    rc = str(DNA(_seq).reverse_complement())
                                    # qual score should be simply reversed
                                    rev_out.write(
                                        ('\n'.join([head, rc, qh, q[::-1]]) +
                                         '\n').encode('utf-8'))
                                else:
                                    # the alternative is no match...
                                    out_notmatched.write(
                                        ('\n'.join(_read) +
                                         '\n').encode('utf-8'))
                    run_command(['gzip', str(rev_path_out)])
        # Revcomp mode: if no reference is passed, revcomp all
        else:
            _vsearch_revcomp_fastq(fwd_path_in, fwd_path_out)
            if paired:
                _vsearch_revcomp_fastq(rev_path_in, rev_path_out)

        # Handle notmatched files
        # This notmatched file might only be created when non-empty
        # and does not exist when performing a simple revcomp
        # but a file is required for the output format
        # so we will manufacture any missing files
        notmatched_fps = [str(fwd_notmatched)]
        if paired:
            notmatched_fps += [str(rev_notmatched)]
        for notmatched_fp in notmatched_fps:
            # create empty file if needed
            if not os.path.exists(notmatched_fp):
                open(notmatched_fp, 'w').close()
            run_command(['gzip', str(notmatched_fp)])

    return oriented, notmatched
