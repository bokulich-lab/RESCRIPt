# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import shutil
import gzip
from urllib.request import urlretrieve
from pathlib import Path
import pandas as pd
from q2_types.feature_data import (TaxonomyFormat, DNAFASTAFormat,
                                   DNAIterator)

MITO_GENE_LIST = ['A6', 'A8', 'CO1', 'CO2', 'CO3', 'Cytb',
                  'ND1', 'ND2', 'ND3', 'ND4L', 'ND4',
                  'ND5', 'ND6', 'lrRNA', 'srRNA', 'all']


def _assemble_midori2_urls(mito_gene,
                           version='GenBank267_2025-06-19',
                           ref_seq_type='uniq',
                           unspecified_species=False,
                           ):
    print('Building URLs ...')
    if unspecified_species:
        spec = 'QIIME_sp'
        sp = 'SP_'
    else:
        spec = 'QIIME'
        sp = ''

    vernum = version.strip('GenBank').split('_')[0]

    base_url = ('https://www.reference-midori.info/download/Databases/'
                '{full_version}/{spec}/{ref_seq_type}/MIDORI2_'
                '{ref_seq_type_u}_NUC_'
                '{sp}GB{vernum}_{mito_gene}_QIIME.{file_type}.gz')

    base_format_dict = {'full_version': version,
                        'spec': spec,
                        'ref_seq_type': ref_seq_type,
                        'ref_seq_type_u': ref_seq_type.upper(),
                        'sp': sp,
                        'vernum': vernum,
                        'mito_gene': mito_gene}

    fasta_url = base_url.format(**base_format_dict, **{'file_type': 'fasta'})
    tax_url = base_url.format(**base_format_dict, **{'file_type': 'taxon'})

    return fasta_url, tax_url


def _retrieve_data_from_midori2(fasta_url, tax_url):

    print('\nDownloading and processing raw files ... \n')

    with tempfile.TemporaryDirectory() as tmpdirname:
        for url in [fasta_url, tax_url]:
            compressed_fn = Path(url).name
            uncompressed_fn = Path(url).stem
            in_path = os.path.join(tmpdirname, compressed_fn)
            out_path = os.path.join(tmpdirname, uncompressed_fn)

            try:
                print('Retrieving data from {0}'.format(url))
                urlretrieve(url, in_path)
            except Exception as e:
                raise ValueError(
                    'Unable to retrieve the following file '
                    'from MIDORI2:\n{url}\n: {e}'.format(
                                      url=url, e=e))

            print('  Unzipping {0}...\n'.format(in_path))
            with gzip.open(in_path, 'rt') as gz_in:
                with open(out_path, 'w') as gz_out:
                    shutil.copyfileobj(gz_in, gz_out)
                    if out_path.endswith('fasta'):
                        seqs = DNAFASTAFormat(out_path,
                                              mode="r").view(DNAIterator)
                    elif out_path.endswith('taxon'):
                        tax = TaxonomyFormat(out_path,
                                             mode="r").view(pd.DataFrame)
    return seqs, tax


# Note: add ability to trim txids from each rank.
# In the interim, can do this with `edit-taxonomy`
# using the following regex:
#     `--p-search-strings '_\d+(;)|_\d+($)'`.
def get_midori2_data(
    mito_gene: list,
    version: str = 'GenBank267_2025-06-19',
    ref_seq_type: str = 'uniq',
    unspecified_species: bool = False,
        ) -> (DNAIterator, pd.DataFrame):

    seq_res = {}
    tax_res = {}

    if 'all' in mito_gene:
        MITO_GENE_LIST.pop()  # remove 'all' from choices.
        mito_gene = MITO_GENE_LIST

    for gene in mito_gene:
        fasta_url, tax_url = _assemble_midori2_urls(
                                  mito_gene=gene,
                                  version=version,
                                  ref_seq_type=ref_seq_type,
                                  unspecified_species=unspecified_species,
                                  )

        seqs, tax = _retrieve_data_from_midori2(fasta_url, tax_url)
        seq_res[gene + '_seqs'] = seqs
        tax_res[gene + '_taxa'] = tax

    print('\n Saving files...\n')
    return seq_res, tax_res
