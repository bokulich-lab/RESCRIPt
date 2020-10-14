from typing import Union

import qiime2
from q2_types.feature_data import (
    AlignedDNAFASTAFormat, DNAFASTAFormat,
    DNAIterator, AlignedDNAIterator)
from skbio import DNA


def _trim_sequence(sequence: DNA,
                   position_start: int,
                   position_end: int) -> DNA:
    return DNA(sequence[position_start:position_end])


def _trim_all_sequences(aligned_sequences: AlignedDNAFASTAFormat,
                        trim_positions: dict) -> AlignedDNAFASTAFormat:
    """
    Trim all sequences within given alignment based on provided positions.

    Arguments:
        aligned_sequences (AlignedDNAFASTAFormat): original, aligned sequences
        trim_positions (dict): dictionary containing positions for trimming

    Returns:
        result (AlignedDNAFASTAFormat): trimmed aligned sequences
    """

    result = AlignedDNAFASTAFormat()
    with result.open() as out_fasta:
        for seq in aligned_sequences.view(DNAIterator):
            seq_trimmed = _trim_sequence(
                seq,
                trim_positions["pos_start"],
                trim_positions["pos_end"])
            seq_trimmed.write(out_fasta)
    return result


def _find_terminal_positions(primer_positions: dict) -> (int, int):
    """
    Identify left- (5') and rightmost (3') trimming position. If both primers
    were used return: (start index of the forward primer, end index of
    the reverse primer). If only forward primer was found return
    (its start index, None). If only reverse primer was found return
    (None, its end index).

    Arguments:
        primer_positions (dict): dictionary containing start and end positions
                                    of forward and reverse primers

    Returns:
        pos_start (int): leftmost (5') position index for trimming
        pos_end (int): rightmost (3') position index for trimming
    """

    pos_start, pos_end = None, None
    # when we found both primers
    if len(primer_positions) == 2:
        pos_start = min([x["pos_start"] for x in primer_positions.values()])
        pos_end = max([x["pos_end"] for x in primer_positions.values()])
    # when only fwd primer was used
    elif "primer_fwd" in primer_positions.keys():
        pos_start = primer_positions["primer_fwd"]["pos_start"]
    # when only rev primer was used
    elif "primer_rev" in primer_positions.keys():
        pos_end = primer_positions["primer_rev"]["pos_end"]
    return pos_start, pos_end


def _locate_primer_positions(
        alignment_with_primers: AlignedDNAFASTAFormat) -> dict:
    """
    Identify position of each primer within the alignment.

    Arguments:
        alignment_with_primers (AlignedDNAFASTAFormat): sequence alignment
                        containing at least one aligned primer

    Returns:
        (dict): dictionary containing trimming positions
                using 0-based indexing
    """

    primers_aligned = dict()
    for aln_seq in alignment_with_primers.view(DNAIterator):
        if aln_seq.metadata["id"] in ["primer_fwd", "primer_rev"]:
            primers_aligned[aln_seq.metadata["id"]] = (str(aln_seq))

    primer_positions = dict()
    for primer_id, primer_seq in primers_aligned.items():
        primer_positions[primer_id] = {
            'pos_start': next(
                (i for i, nt in enumerate(primer_seq) if nt != "-")),
            'pos_end': len(primer_seq) - next(
                (i for i, nt in enumerate(primer_seq[::-1]) if nt != "-"))
        }

    pos_start, pos_end = _find_terminal_positions(primer_positions)

    # TODO: probably some validation like in _prepare_positions
    #  should happen here too

    return {"pos_start": pos_start, "pos_end": pos_end}


def _prepare_positions(position_start: Union[int, None],
                       position_end: Union[int, None],
                       aln_len: int) -> dict:
    """
    Prepare a dictionary with start and end position for sequence trimming.

    Arguments:
        position_start (int, None): leftmost (5') position; 1-based indexing
        position_end (int, None): rightmost (3') position; 1-based indexing
        aln_len (int): length of the alignment

    Returns:
        (dict): dictionary containing trimming positions
                using 0-based indexing
    """

    if not position_start and not position_end:
        raise ValueError("Neither primers nor sequence positions were "
                         "provided - nothing to be trimmed.")
    if position_start:
        if position_end and position_start >= position_end:
            raise ValueError("Start position should be smaller than end "
                             "position while trimming.")
        elif position_start >= aln_len:
            raise ValueError("Start position should be smaller than "
                             "alignment length.")
        else:
            position_start -= 1
    if position_end and position_end > aln_len:
        position_end = aln_len

    return {'pos_start': position_start, 'pos_end': position_end}


def _process_primers(
        primer_fwd: Union[str, None],
        primer_rev: Union[str, None]) -> DNAFASTAFormat:
    """
    Convert provided primers into skbio DNA format. Will reverse complement
    the reverse primer, if provided.

    Arguments:
        primer_fwd (str, None): forward primer
        primer_rev (str, None): reverse primer

    Returns:
        primers_fasta (DNAFASTAFormat): primers in FASTA format
    """

    primers = {
        'primer_fwd': DNA(
            primer_fwd, metadata={
                'id': 'primer_fwd'})
        if primer_fwd else None,
        'primer_rev': DNA(
            primer_rev, metadata={
                'id': 'primer_rev'}).reverse_complement()
        if primer_rev else None
    }

    # save primers in that format to pass them to mafft_add
    primers_fasta = DNAFASTAFormat()
    with primers_fasta.open() as out:
        [primer.write(out) for primer in primers.values() if primer]

    return primers_fasta


def _trim_alignment(expand_alignment_action,
                    aligned_sequences,
                    primer_fwd=None,
                    primer_rev=None,
                    position_start=None,
                    position_end=None) -> AlignedDNAFASTAFormat:
    """
    Trim alignment based on primer alignment or explicitly specified
    positions. When at least one primer sequence is given, primer-based
    trimming will be performed, otherwise position-based trimming is done.

    Arguments:
        expand_alignment_action: qiime action for multiple seq. alignment
        aligned_sequences (AlignedDNAFASTAFormat): alignment to be trimmed
        primer_fwd (str, None): forward primer used to find the leftmost (5')
                            trimming position
        primer_rev (str, None): reverse primer used to find the rightmost (3')
                            trimming position
        position_start (int, None): leftmost (5') trimming position
        position_end (int, None): rightmost (3') trimming position

    Returns:
        result (AlignedFASTAFormat): trimmed alignment

    """

    # if primers are provided use them
    # otherwise fall back to position-based slicing
    if primer_fwd is not None or primer_rev is not None:
        primers_fasta = _process_primers(primer_fwd, primer_rev)
        primers = qiime2.Artifact.import_data(
            'FeatureData[Sequence]', primers_fasta)

        # expand the existing alignment by addition of primers
        alignment_with_primers = expand_alignment_action(
            alignment=aligned_sequences, sequences=primers)[0]

        # find trim positions based on primer positions within alignment
        primer_positions = _locate_primer_positions(alignment_with_primers)
    else:
        # find length of the alignment
        aln_len = len(
            next(
                aligned_sequences.view(AlignedDNAIterator).generator))

        primer_positions = _prepare_positions(
            position_start, position_end, aln_len)

    result = _trim_all_sequences(aligned_sequences, primer_positions)

    return result


def trim_alignment(ctx,
                   aligned_sequences,
                   primer_fwd=None,
                   primer_rev=None,
                   position_start=None,
                   position_end=None):
    expand_alignment_action = ctx.get_action('alignment', 'mafft_add')

    result = _trim_alignment(
        expand_alignment_action,
        aligned_sequences,
        primer_fwd,
        primer_rev,
        position_start,
        position_end)

    return qiime2.Artifact.import_data('FeatureData[AlignedSequence]', result)
