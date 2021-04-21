# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from q2_types.feature_data import (
    DNAIterator, AlignedDNAFASTAFormat, RNAFASTAFormat, DNAFASTAFormat)

from ..plugin_setup import plugin
from ._format import SILVATaxonomyFormat, SILVATaxidMapFormat
from rescript._utilities import (
    _rna_to_dna, _read_dna_fasta, _dna_iterator_to_aligned_fasta)


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


# TODO: this should probably be moved to q2-types
@plugin.register_transformer
def _6(data: RNAFASTAFormat) -> DNAFASTAFormat:
    return _rna_to_dna(str(data))


# TODO: this should probably be moved to q2-types
@plugin.register_transformer
def _7(data: RNAFASTAFormat) -> DNAIterator:
    converted_dna = _rna_to_dna(str(data))
    generator = _read_dna_fasta(str(converted_dna))
    return DNAIterator(generator)


# TODO: this should probably be moved to q2-types
@plugin.register_transformer
def _8(data: DNAIterator) -> AlignedDNAFASTAFormat:
    return _dna_iterator_to_aligned_fasta(data)
