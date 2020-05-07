#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import skbio
import skbio.io
from qiime2.plugin import (Str, Plugin, Choices, List, Citations, Range, Int,
                           Float)
from q2_types.feature_data import (FeatureData, Sequence, DNAFASTAFormat,
                                   DNAIterator)
import qiime2
import rescript


def _filt_seq_with_ambig_bases(seq, num_ambig_bases):
    ambig_bases_in_seq = sum(seq.degenerates())
    if ambig_bases_in_seq >= num_ambig_bases:
        if ambig_bases_in_seq == 0: # handle case for num_ambig_bases==0
            return False
        else:
            return True
    else:
        return False

def _filter_homopolymer(seq, homopolymer_length):
    nhl = homopolymer_length - 1 # due to how regex is written
    if nhl < 2:
        raise ValueError("Homopolymer length must be >= 2!")
    else:
        regex_str = "([ACGTURYSWKMBDHVN])\\1{%s,}" % nhl
        for p in re.finditer(regex_str, str(seq)):
            if len(p.group()) >= homopolymer_length:
                return True
            else:
                continue
        return False

def clean_sequences(sequences: DNAFASTAFormat, num_ambig_bases: int=5,
                    homopolymer_length: int=8) -> DNAFASTAFormat:

    result = DNAFASTAFormat()
    result_fp = str(result)

    with open(result_fp, 'w') as out_fasta:
        for seq in sequences.view(DNAIterator):
            ambig = _filt_seq_with_ambig_bases(seq, num_ambig_bases)
            if ambig == False:
                poly = _filter_homopolymer(seq, homopolymer_length)
                if poly == False: # if we make it here, write seq to file
                    seq.write(out_fasta)
    return result

