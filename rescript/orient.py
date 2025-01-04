# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from typing import List, Dict, Any

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
    threads: int = 1,
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
            '--threads', str(threads),
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
