# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import skbio
from q2_types.feature_data import (
    DNAFASTAFormat, DNAIterator, AlignedDNAFASTAFormat, RNAFASTAFormat)

from ..plugin_setup import plugin
from ._format import SILVATaxonomyFormat, SILVATaxidMapFormat
from rescript._utilities import (_rna_to_dna, _rna_to_dna_iterator,
                                 _read_fasta, _dna_iterator_to_aligned_fasta)


def _read_dataframe(fh, header=0):
    # Using `dtype=object` and `set_index` to avoid type casting/inference
    # of any columns or the index. E.g., taxids should be left as str.
    df = pd.read_csv(fh, sep='\t', header=header, dtype='str')
    df.set_index(df.columns[0], drop=True, append=False, inplace=True)
    df.index.name = 'id'
    return df


@plugin.register_transformer
def _1(data: pd.DataFrame) -> (SILVATaxonomyFormat):
    ff = SILVATaxonomyFormat()
    with ff.open() as fh:
        data.to_csv(fh, sep='\t', header=False)
    return ff


@plugin.register_transformer
def _2(ff: SILVATaxonomyFormat) -> (pd.DataFrame):
    with ff.open() as fh:
        df = _read_dataframe(fh, header=None)
        return df


@plugin.register_transformer
def _4(data: pd.DataFrame) -> (SILVATaxidMapFormat):
    ff = SILVATaxidMapFormat()
    with ff.open() as fh:
        data.to_csv(fh, sep='\t', header=True)
    return ff


@plugin.register_transformer
def _5(ff: SILVATaxidMapFormat) -> (pd.DataFrame):
    with ff.open() as fh:
        df = _read_dataframe(fh, header=0)
        # normalize column names
        df.columns = ['start', 'stop', 'path', 'organism_name', 'taxid']
        return df


@plugin.register_transformer
def _6(data: RNAFASTAFormat) -> DNAFASTAFormat:
    iterator = _read_fasta(str(data), constructor=skbio.RNA)
    ff = _rna_to_dna(iterator)
    return ff


@plugin.register_transformer
def _7(data: RNAFASTAFormat) -> DNAIterator:
    iterator = _read_fasta(str(data), constructor=skbio.RNA)
    generator = _rna_to_dna_iterator(iterator)
    return DNAIterator(generator)


@plugin.register_transformer
def _8(data: DNAIterator) -> AlignedDNAFASTAFormat:
    return _dna_iterator_to_aligned_fasta(data)
