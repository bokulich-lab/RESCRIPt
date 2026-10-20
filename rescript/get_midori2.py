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
from q2_types.feature_data import (DNASequencesDirectoryFormat,
                                   HeaderlessTSVTaxonomyDirectoryFormat)


VER_ID_DICT = {
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
    base_format_dict = {'id': VER_ID_DICT[version],
                        'ver': version,
                        'sp': sp}

    final_url = base_url.format(**base_format_dict)
    return final_url


def _retrieve_data_from_midori2(url):
    with tempfile.TemporaryDirectory() as tmpdirname:
        compressed_fn = Path(url).name
        in_path = os.path.join(tmpdirname, compressed_fn)
        try:
            print('Retrieving data from:\n{0}'.format(url),
                  '\n\nGrab some coffee, this may take a while.\n')
            urlretrieve(url, in_path)
        except Exception as e:
            raise ValueError(
                'Unable to retrieve the following file '
                'from MIDORI2:\n{url}\n: {e}'.format(
                                  url=url, e=e))

        print('Unzipping {0} ...'.format(compressed_fn))
        seqs_dict = {}
        tax_dict = {}

        with zipfile.ZipFile(in_path, 'r') as zipdir:
            for f_info in zipdir.infolist():
                bn = Path(f_info.filename)
                if f_info.is_dir() or bn.suffix != '.gz':
                    continue
                bn = bn.with_suffix('')
                if bn.suffix == '.fasta':
                    res = DNASequencesDirectoryFormat()
                    res_suffix = 'sequences'
                    res_dict = seqs_dict
                elif bn.suffix == '.taxon':
                    res = HeaderlessTSVTaxonomyDirectoryFormat()
                    res_suffix = 'taxonomy'
                    res_dict = tax_dict
                else:
                    continue
                with zipdir.open(f_info.filename) as zipf:
                    with gzip.open(zipf, 'rb') as gz_in:
                        with res.file.path_maker().open('wb') as out:
                            shutil.copyfileobj(gz_in, out)

                bofn = '_'.join([bn.stem, res_suffix])
                res_dict[bofn] = res
    return seqs_dict, tax_dict


def get_midori2_data(
    version: str = '273',
    unspecified_species: bool = False,
        ) -> (DNASequencesDirectoryFormat,
              HeaderlessTSVTaxonomyDirectoryFormat):

    url = _assemble_midori2_url(
                                version=version,
                                unspecified_species=unspecified_species,
                                )

    seq_res, tax_res = _retrieve_data_from_midori2(url)

    print('\nSaving files...\n')
    return seq_res, tax_res
