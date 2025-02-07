# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import shutil
import gzip
import warnings
from urllib.request import urlretrieve
from urllib.error import HTTPError
from pathlib import Path
import pandas as pd
from collections import OrderedDict
from q2_types.feature_data import (TaxonomyFormat, DNAFASTAFormat,
                                   DNAIterator)

# Note: No need for '5.0.1' '4.14.1'. Documentaiton says to refer to
# the '5.0.0' and '4.14.0' files respectively.
# Also, might need to include extra code for 4.13 and earlier as
# the 16S and 18S files are separated. Note the slight differences
# in naming convention:
#   pr2_version_5.0.0_SSU_mothur.fasta.gz
#  VS
#   pr2_version_4.13.0_16S_mothur.fasta.gz
#   pr2_version_4.13.0_18S_mothur.fasta.gz

_allowed_pr2_ranks = OrderedDict({'domain': 'd__', 'supergroup': 'sgr__',
                                  'division': 'dv__', 'subdivision': 'dvs__',
                                  'class': 'c__', 'order': 'o__',
                                  'family': 'f__', 'genus': 'g__',
                                  'species': 's__'})
_default_pr2_ranks = ['domain', 'supergroup', 'division', 'subdivision',
                      'class', 'order', 'family', 'genus', 'species']


def get_pr2_data(
    version: str = '5.0.0',
    ranks: list = None,
        ) -> (DNAIterator, pd.Series):

    urls = _assemble_pr2_urls(version=version)
    seqs, tax = _retrieve_data_from_pr2(urls, ranks)

    print('\n Saving files...\n')
    return seqs, tax


def _assemble_pr2_urls(version='5.0.0'):
    urls_to_retrieve = []
    base_url = ('https://github.com/pr2database/pr2database/releases/download/'
                'v{ver}/pr2_version_{ver}_SSU_mothur.{ext}.gz')

    for data_type in ['fasta', 'tax']:
        urls_to_retrieve.append(base_url.format(**{'ver': version,
                                                   'ext': data_type}))
    return urls_to_retrieve


def _compile_taxonomy_output(tax, ranks):
    # prepare dataframe with all ranks:
    prefix_list = list(_allowed_pr2_ranks.values())
    tax[prefix_list] = \
        tax['Taxon'].str.strip(';').str.split(';', expand=True)
    tax.drop('Taxon', axis=1, inplace=True)
    # prepend prefixes
    tax.loc[:, prefix_list] = \
        tax.loc[:, prefix_list].apply(lambda x: x.name + x)
    # sort user defined ranks in case provided out of order
    # then only return user specified ranks
    sorted_ranks = [p for r, p in _allowed_pr2_ranks.items() if r in ranks]
    taxonomy = tax.loc[:, sorted_ranks].agg('; '.join, axis=1)
    taxonomy.rename('Taxon', inplace=True)
    return taxonomy


def _retrieve_data_from_pr2(urls_to_retrieve, ranks):
    # Perform check that the `urls_to_retriev` should only
    # contain 2 files, a 'fasta' and 'taxonomy' file.

    if ranks is None:
        ranks = _default_pr2_ranks

    print('\nDownloading and processing raw files ... \n')

    with tempfile.TemporaryDirectory() as tmpdirname:
        for url in urls_to_retrieve:
            compressed_fn = Path(url).name
            uncompressed_fn = Path(url).stem
            in_path = os.path.join(tmpdirname, compressed_fn)
            out_path = os.path.join(tmpdirname, uncompressed_fn)

            try:
                print('Retrieving data from {0}'.format(url))
                urlretrieve(url, in_path)
            except HTTPError:
                msg = ("Unable to retrieve the followng file from PR2:\n"
                       + url)
                warnings.warn(msg, UserWarning)

            print('  Unzipping {0}...\n'.format(in_path))
            with gzip.open(in_path, 'rt') as gz_in:
                with open(out_path, 'w') as gz_out:
                    shutil.copyfileobj(gz_in, gz_out)
                    if out_path.endswith('fasta'):
                        seqs = DNAFASTAFormat(out_path,
                                              mode="r").view(DNAIterator)
                    elif out_path.endswith('tax'):
                        tax = TaxonomyFormat(out_path,
                                             mode="r").view(pd.DataFrame)
                        updated_tax = _compile_taxonomy_output(tax, ranks)
    return seqs, updated_tax
