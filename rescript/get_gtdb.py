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
import gzip
import tarfile

import qiime2
from urllib.request import urlretrieve
from collections import defaultdict
# from urllib.error import HTTPError
# from q2_types.feature_data import DNAFASTAFormat, DNAIterator
# from q2_feature_table import merge_seqs, merge_taxa


# different versions may have different file names for archaea and
# bacteria. for example 'ar53' and 'bac120' mean that the GTDB phylogeny
# is based on 53 and 120 concatenated proteins (cp), respectively.
# If this changes we can set up a conditional statemnt below.
VERSION_DICT = {'207': {'Archaea': 'ar53', 'Bacteria': 'bac120'},
                '202': {'Archaea': 'ar122', 'Bacteria': 'bac120'}}


def get_gtdb_data(ctx, version='207', domain='Both'):

    merge_gtdb_seqs = ctx.get_action('feature_table', 'merge_seqs')
    merge_gtdb_taxonomy = ctx.get_action('feature_table', 'merge_taxa')

    d = defaultdict(lambda: defaultdict(dict))

    if domain == 'Both':
        d[version] = VERSION_DICT[version]
    else:
        d[version][domain] = VERSION_DICT[version][domain]

    # download data from gtdb
    print('\nDownloading and processing raw files ... \n')
    queries = _assemble_gtdb_data_urls(d)
    gtdb_data = _retrieve_data_from_gtdb(queries)
    tax_q, seqs_q = gtdb_data

    print('\n Merging taxonomy data...')
    merged_gtdb_tax, = merge_gtdb_taxonomy(data=tax_q)
    print('\n Merging sequence data...')
    merged_gtdb_seqs, = merge_gtdb_seqs(data=seqs_q)

    print('\n Saving files...')

    return merged_gtdb_tax, merged_gtdb_seqs


def _assemble_gtdb_data_urls(d):
    '''Generate gtdb urls, given database version and reference target.'''

    # Compile URLs
    base_seq_url = ('https://data.gtdb.ecogenomic.org/releases/release%s/'
                    '%s.0/genomic_files_reps/%s_ssu_reps_r%s.tar.gz')
    base_tax_url = ('https://data.gtdb.ecogenomic.org/releases/release%s/'
                    '%s.0/%s_taxonomy_r%s.tsv.gz')

    # taxonomy
    tax_queries = [(domain + '_taxonomy',
                   base_tax_url % (version, version, cp, version),
                   'FeatureData[Taxonomy]', 'HeaderlessTSVTaxonomyFormat')
                   for version, dcp in d.items() for domain, cp in
                   dcp.items()]
    # sequences
    seq_queries = [(domain + '_sequences',
                   base_seq_url % (version, version, cp, version),
                   'FeatureData[Sequence]', 'DNAFASTAFormat')
                   for version, dcp in d.items() for domain, cp in
                   dcp.items()]

    return tax_queries, seq_queries


def _retrieve_data_from_gtdb(queries):
    '''
    Download data from gtdb, given a list of queries.

    queries: list of tuples of (str, str)
        (name, urlpath)
    '''
    tax_results = []
    seq_results = []
    with tempfile.TemporaryDirectory() as tmpdirname:
        for q in queries:
            for name, url, type, format in q:

                print('Retrieving {0} from {1}'.format(name, url))
                # grab url
                bn = os.path.basename(url)
                destination = os.path.join(tmpdirname, bn)
                urlretrieve(url, destination)

                # seq files are contained within `tar.gz``
                if tarfile.is_tarfile(destination):
                    try:
                        untarred_fn = bn.split('.')[0]+'.fna'
                        untarred_destination = os.path.join(
                                               tmpdirname, untarred_fn)
                        print('  Untarring {0}...'.format(bn))
                        with tarfile.open(destination, 'r') as tar:
                            tar.extract(member=untarred_fn,
                                        path=tmpdirname)
                    except OSError:
                        pass

                    seq_results.append(qiime2.Artifact.import_data(
                                                        type,
                                                        untarred_destination,
                                                        format))

                # taxonomy files are contained within ``.gz`
                else:
                    try:
                        unzipped_destination = os.path.splitext(destination)[0]
                        print('  Unzipping {0}...'.format(name))
                        _gzip_decompress(destination, unzipped_destination)
                        destination = unzipped_destination
                    except OSError:
                        pass

                    tax_results.append(qiime2.Artifact.import_data(
                                                                   type,
                                                                   destination,
                                                                   format))

        return tax_results, seq_results


def _gzip_decompress(input_fp, output_fp):
    with gzip.open(input_fp, 'rt') as temp_in:
        with open(output_fp, 'w') as temp_out:
            shutil.copyfileobj(temp_in, temp_out)
