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
import zipfile
import gzip
from urllib.request import urlretrieve
from pathlib import Path
import pandas as pd
from q2_types.feature_data import (TaxonomyFormat, DNAFASTAFormat,
                                   DNAIterator)


ver_id_dict = {
               '273': '22685701',
               '272': '23112349',
               '271': '23110935',
               '270': '23093552',
               '269': '23092095',
               '268': '23090517',
               '267': '23071879',
               '266': '23049405',
               '265': '23047760',
               }


def _assemble_midori2_url(version='273',
                          unspecified_species=False,
                          ):

    print('Building URL ...')
    # unspecified species check
    if unspecified_species:
        sp = '_sp'
    else:
        sp = ''

    base_url = ('https://zenodo.org/records/{id}/'
                'files/MIDORI2_GB{ver}_QIIME{sp}.zip')
    base_format_dict = {'id': ver_id_dict[version],
                        'ver': version,
                        'sp': sp}

    final_url = base_url.format(**base_format_dict)
    return final_url


def _retrieve_data_from_midori2(url):

    print('\nDownloading and processing MIDORI2 data ... \n')

    with tempfile.TemporaryDirectory() as tmpdirname:
        compressed_fn = Path(url).name
        in_path = os.path.join(tmpdirname, compressed_fn)

        try:
            print('Retrieving data from:\n\t{0}'.format(url),
                  '\n\nGrab some coffee, this may take a while.\n\n')
            urlretrieve(url, in_path)
        except Exception as e:
            raise ValueError(
                'Unable to retrieve the following file '
                'from MIDORI2:\n{url}\n: {e}'.format(
                                  url=url, e=e))

        print('Unzipping {0}...\n'.format(in_path))
        seqs_dict = {}
        tax_dict = {}
        with zipfile.ZipFile(in_path, 'r') as zipdir:
            for f_info in zipdir.infolist():
                if f_info.is_dir() or not f_info.filename.endswith('.gz'):
                    continue
                fn = f_info.filename
                with zipdir.open(f_info.filename) as zipf:
                    bn = os.path.basename(fn).rsplit('.', 1)[0]
                    with gzip.open(zipf, 'rb') as gz_in:
                        ofp = os.path.join(tmpdirname, bn)
                        with open(ofp, 'wb') as gz_out:
                            shutil.copyfileobj(gz_in, gz_out)

                    if gz_out.name.endswith('fasta'):
                        seqs_dict[bn] = DNAFASTAFormat(
                                            gz_out.name,
                                            mode="r").view(DNAIterator)
                    elif gz_out.name.endswith('taxon'):
                        tax_dict[bn] = TaxonomyFormat(
                                           gz_out.name,
                                           mode="r").view(pd.DataFrame)
    return seqs_dict, tax_dict


def get_midori2_data(
    version: str = '273',
    unspecified_species: bool = False,
        ) -> (DNAIterator, pd.DataFrame):

    seq_res = {}
    tax_res = {}

    url = _assemble_midori2_url(
                                version=version,
                                unspecified_species=unspecified_species,
                                )

    seq_res, tax_res = _retrieve_data_from_midori2(url)

    print('\n Saving files...\n')
    return seq_res, tax_res
