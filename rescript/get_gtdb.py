# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import shutil
import tarfile
import gzip
import warnings
import pandas as pd
from urllib.request import urlretrieve
from collections import defaultdict
from urllib.error import HTTPError

from q2_types.feature_data import (TSVTaxonomyFormat, DNAFASTAFormat,
                                   DNAIterator)

# Different versions may have different file names for archaea and
# bacteria. for example 'ar53' and 'bac120' mean that the GTDB phylogeny
# is based on 53 and 120 concatenated proteins (cp), respectively.
# If this changes we can set up a conditional statemnt below.
VERSION_MAP_DICT = {'220.0': {'Archaea': 'ar53', 'Bacteria': 'bac120'},
                    '214.1': {'Archaea': 'ar53', 'Bacteria': 'bac120'},
                    '214.0': {'Archaea': 'ar53', 'Bacteria': 'bac120'},
                    '207.0': {'Archaea': 'ar53', 'Bacteria': 'bac120'},
                    '202.0': {'Archaea': 'ar122', 'Bacteria': 'bac120'}}


def get_gtdb_data(
    version: str = '220.0',
    domain: str = 'Both',
    db_type: str = 'SpeciesReps',
        ) -> (TSVTaxonomyFormat, DNAFASTAFormat):

    queries = _assemble_queries(version=version,
                                db_type=db_type,
                                domain=domain)
    tax_q, seqs_q = _retrieve_data_from_gtdb(queries)

    print('\n Saving files...\n')
    return tax_q, seqs_q


def _assemble_queries(version='220.0',
                      db_type='SpeciesReps',
                      domain='Both'):
    queries = []
    base_url = 'https://data.gtdb.ecogenomic.org/releases/'
    base_version = version.split('.')[0]
    # ^^ Set `base_version` variable becuase number after the decimal is
    # only used for the directory. GTDB trims this off for the actual
    # file names...
    # GTDB v220 started storing the ssu_reps FASTA files
    # as 'fna.gz' instead of their usual 'tar.gz'.
    if version == '220.0':
        stype = 'fna'
    else:
        stype = 'tar'

    if db_type == 'SpeciesReps':
        ver_dom_dict = defaultdict(lambda: defaultdict(dict))

        if domain == 'Both':
            ver_dom_dict[version] = VERSION_MAP_DICT[version]
        else:
            ver_dom_dict[version][domain] = VERSION_MAP_DICT[version][domain]

        full_url = (base_url + 'release{bver}/{ver}/genomic_files_reps/'
                    '{cp}_ssu_reps_r{bver}.{stype}.gz')

        for version, dcp in ver_dom_dict.items():
            for dom, cp in dcp.items():
                queries.append((dom,
                                full_url.format(**{'ver': version,
                                                   'bver': base_version,
                                                   'cp': cp,
                                                   'stype': stype})))
    elif db_type == 'All':
        # Note: GTDB does not maintain separate 'Bacteria' and
        # 'Archaea' files for 'All'. This is only done for
        # the 'SpeciesReps'. Again, account for filename changes
        # to 'fna.gz' in v220, i.e. 'stype'.
        full_url = (base_url + 'release{bver}/{ver}/genomic_files_all/'
                    'ssu_all_r{bver}.{stype}.gz')

        queries.append((db_type,
                        full_url.format(**{'ver': version,
                                           'bver': base_version,
                                           'stype': stype})))
    return queries


def parse_gtdb_taxonomy(tax_str):
    tax = ' '.join(tax_str.split()[0:2])
    return tax


def _get_gtdb_data_path(tmpdirname, url, basename):
    destination = os.path.join(tmpdirname, basename)
    try:
        print('Retrieving data from {0}'.format(url))
        urlretrieve(url, destination)
    except HTTPError:
        msg = ("Unable to retrieve the followng file from GTDB:\n "
               + url)
        warnings.warn(msg, UserWarning)
    return destination


def _extract_seq_file(tmpdirname, uncompressed_fn, destination):
    # GTDB v220 updated the file format for the "ssu_rep" FASTA files
    # to 'fna.gz' instead of `tar.gz`. We now need to check which
    # form the file is in.
    out_path = os.path.join(tmpdirname, uncompressed_fn)

    if tarfile.is_tarfile(destination):
        with tarfile.open(destination, 'r') as tar:
            print('  Untarring {0}...\n'.format(uncompressed_fn))
            tar.extract(member=uncompressed_fn,
                        path=tmpdirname)
            seqs = DNAFASTAFormat(out_path,
                                  mode="r").view(DNAIterator)
    if destination.endswith('fna.gz'):
        with gzip.open(destination, 'rt') as gz_in:
            print('  Unzipping {0}...\n'.format(uncompressed_fn))
            with open(out_path, 'w') as gz_out:
                shutil.copyfileobj(gz_in, gz_out)
                seqs = DNAFASTAFormat(out_path,
                                      mode="r").view(DNAIterator)
    return seqs


def _retrieve_data_from_gtdb(queries):
    proc_seqs = DNAFASTAFormat()
    proc_tax = TSVTaxonomyFormat()
    tax_dict = {}

    print('\nDownloading and processing raw files ... \n')
    with \
        tempfile.TemporaryDirectory() as tmpdirname, \
            proc_seqs.open() as out_fasta, \
            proc_tax.open() as out_tax:

        for domain, url in queries:
            # variable setup
            basename = os.path.basename(url)
            uncompressed_fn = basename.split('.')[0]+'.fna'
            destination = _get_gtdb_data_path(tmpdirname, url, basename)
            seqs = _extract_seq_file(tmpdirname, uncompressed_fn, destination)
            print('  Writing data from \'{0}\'.\n'.format(domain))
            for seq in seqs:
                seq.write(out_fasta)  # write seq to new fasta file
                # add taxonomy to dict:
                tax_dict[seq.metadata['id']] = parse_gtdb_taxonomy(
                         seq.metadata['description'])
        # set up final taxonomy dataframe:
        print('  Sequences processed.')
        print('  Processing taxonomy...')
        parsed_taxonomy_df = pd.DataFrame.from_dict(tax_dict, orient='index')
        parsed_taxonomy_df.index.name = 'Feature ID'
        parsed_taxonomy_df.rename(columns={parsed_taxonomy_df.columns[0]:
                                           'Taxon'}, inplace=True)
        parsed_taxonomy_df.to_csv(out_tax, sep='\t')
        print('  Taxonomy processed.')
        return proc_tax, proc_seqs
