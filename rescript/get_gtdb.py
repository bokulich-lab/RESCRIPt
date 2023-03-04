# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import tarfile
import warnings

import qiime2
from urllib.request import urlretrieve
from collections import defaultdict
from urllib.error import HTTPError
from rescript.get_data import _gzip_decompress

# Different versions may have different file names for archaea and
# bacteria. for example 'ar53' and 'bac120' mean that the GTDB phylogeny
# is based on 53 and 120 concatenated proteins (cp), respectively.
# If this changes we can set up a conditional statemnt below.
VERSION_MAP_DICT = {'207': {'Archaea': 'ar53', 'Bacteria': 'bac120'},
                    '202': {'Archaea': 'ar122', 'Bacteria': 'bac120'}}


def get_gtdb_data(ctx, version='207', domain='Both'):

    ver_dom_dict = defaultdict(lambda: defaultdict(dict))

    # Subset dict if needed, but keep same structure
    # Although we can run the following merge actions on a list of one
    # i.e. 'Archaea', we do not want to confuse anyone when looking
    # at provenance, by running a merge command for no reason.
    if domain == 'Both':
        ver_dom_dict[version] = VERSION_MAP_DICT[version]
    else:
        ver_dom_dict[version][domain] = VERSION_MAP_DICT[version][domain]

    queries = _assemble_queries(ver_dom_dict)
    tax_q, seqs_q = _retrieve_data_from_gtdb(queries)

    if domain == 'Both':
        merge_gtdb_seqs = ctx.get_action('feature_table', 'merge_seqs')
        merge_gtdb_taxonomy = ctx.get_action('feature_table', 'merge_taxa')
        print('\n  Merging taxonomy data...')
        merged_gtdb_tax, = merge_gtdb_taxonomy(data=tax_q)
        print('\n  Merging sequence data...')
        merged_gtdb_seqs, = merge_gtdb_seqs(data=seqs_q)
    else:
        merged_gtdb_tax = tax_q[0]
        merged_gtdb_seqs = seqs_q[0]

    print('\n Saving files...\n')
    return merged_gtdb_tax, merged_gtdb_seqs


def _assemble_queries(ver_dom_dict):
    queries = defaultdict(lambda: list())

    base_seq_url = ('https://data.gtdb.ecogenomic.org/releases/release{ver}/'
                    '{ver}.0/genomic_files_reps/{cp}_ssu_reps_r{ver}.tar.gz')
    base_tax_url = ('https://data.gtdb.ecogenomic.org/releases/release{ver}/'
                    '{ver}.0/{cp}_taxonomy_r{ver}.tsv.gz')

    for version, dcp in ver_dom_dict.items():
        for dom, cp in dcp.items():

            queries['Taxonomy'].append(
                (dom,
                 base_tax_url.format(**{'ver': version, 'cp': cp}),
                 'FeatureData[Taxonomy]',
                 'HeaderlessTSVTaxonomyFormat'))

            queries['Sequence'].append(
               (dom,
                base_seq_url.format(**{'ver': version, 'cp': cp}),
                'FeatureData[Sequence]',
                'DNAFASTAFormat'))
    return queries


def _retrieve_data_from_gtdb(queries):
    '''
    Download data from gtdb, given a list of queries.

    queries: {'Taxonomy':(domain, url, type, format), (...),
              'Sequence':(domain, url, type, format), (...)}
    '''

    tax_results = []
    seq_results = []

    print('\nDownloading and processing raw files ... \n')

    with tempfile.TemporaryDirectory() as tmpdirname:
        for sttype, q_info in queries.items():
            for domain, url, dtype, fmt in q_info:
                print('Retrieving {0} for {1} from {2}'.format(
                      sttype, domain, url))
                # grab url
                bn = os.path.basename(url)
                destination = os.path.join(tmpdirname, bn)
                try:
                    urlretrieve(url, destination)
                except HTTPError:
                    msg = ("Unable to retrieve the followng file from GTDB:\n "
                           + url)
                    warnings.warn(msg, UserWarning)
                # seq files are contained within `tar.gz`
                if tarfile.is_tarfile(destination):
                    seq_results.append(get_tar_data(destination, bn,
                                                    tmpdirname, dtype, fmt))
                else:
                    tax_results.append(get_gzipped_data(destination, bn,
                                                        dtype, fmt))
        return tax_results, seq_results


def get_tar_data(tar_loc, base_fn, tmpdirname, dtype, fmt):
    try:
        untarred_fn = base_fn.split('.')[0]+'.fna'
        print('  Untarring {0}...\n'.format(base_fn))
        with tarfile.open(tar_loc, 'r') as tar:
            tar.extract(member=untarred_fn,
                        path=tmpdirname)
        return qiime2.Artifact.import_data(dtype, os.path.join(
                                       tmpdirname, untarred_fn), fmt)
    except OSError:
        raise OSError(('{0}: either does not exist or can not '
                       'be untarred!'.format(tar_loc)))


def get_gzipped_data(zipped_loc, base_fn, dtype, fmt):
    try:
        unzipped_destination = os.path.splitext(zipped_loc)[0]
        print('  Unzipping {0}...\n'.format(base_fn))
        _gzip_decompress(zipped_loc, unzipped_destination)
        return qiime2.Artifact.import_data(dtype, unzipped_destination, fmt)
    except OSError:
        raise OSError(('{0}: either does not exist or can not '
                       'be unzipped!'.format(unzipped_destination)))
