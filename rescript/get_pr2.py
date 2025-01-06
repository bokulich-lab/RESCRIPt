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
from pandas import DataFrame
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


def get_pr2_data(
    version: str = '5.0.0',
        ) -> (DNAIterator, DataFrame):

    urls = _assemble_pr2_urls(version=version)
    seqs, tax = _retrieve_data_from_pr2(urls)

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


def _retrieve_data_from_pr2(urls_to_retrieve):
    # Perform check that the `urls_to_retriev` should only
    # contain 2 files, a 'fasta' and 'taxonomy' file.

    # seqs = DNAFASTAFormat()
    # tax = HeaderlessTSVTaxonomyFormat()

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
                                             mode="r").view(DataFrame)
    return seqs, tax
