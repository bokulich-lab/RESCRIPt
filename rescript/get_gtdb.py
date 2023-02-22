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
VERSION_DICT = {'207': {'Archaea': 'ar53', 'Bacteria': 'bac120'},
                '202': {'Archaea': 'ar122', 'Bacteria': 'bac120'}}


def get_gtdb_data(ctx, version='207', domain='Both'):

    d = defaultdict(lambda: defaultdict(dict))

    # Subset dict if needed, but keep same structure
    # Although we can run the following merge actions on a list of one
    # i.e. 'Archaea', we do not want to confuse anyone when looking
    # at provenance, by running a merge command for no reason.
    if domain == 'Both':
        d[version] = VERSION_DICT[version]

        queries = _assemble_queries(d)
        tax_q, seqs_q = _retrieve_data_from_gtdb(queries)

        merge_gtdb_seqs = ctx.get_action('feature_table', 'merge_seqs')
        merge_gtdb_taxonomy = ctx.get_action('feature_table', 'merge_taxa')
        print('\n  Merging taxonomy data...')
        merged_gtdb_tax, = merge_gtdb_taxonomy(data=tax_q)
        print('\n  Merging sequence data...')
        merged_gtdb_seqs, = merge_gtdb_seqs(data=seqs_q)

    else:
        d[version][domain] = VERSION_DICT[version][domain]
        queries = _assemble_queries(d)
        tax_q, seqs_q = _retrieve_data_from_gtdb(queries)
        merged_gtdb_tax = tax_q[0]
        merged_gtdb_seqs = seqs_q[0]

    print('\n Saving files...\n')
    return merged_gtdb_tax, merged_gtdb_seqs


def _assemble_queries(d):
    ql = defaultdict(lambda: list())

    base_seq_url = ('https://data.gtdb.ecogenomic.org/releases/release%s/'
                    '%s.0/genomic_files_reps/%s_ssu_reps_r%s.tar.gz')
    base_tax_url = ('https://data.gtdb.ecogenomic.org/releases/release%s/'
                    '%s.0/%s_taxonomy_r%s.tsv.gz')

    for version, dcp in d.items():
        for dom, cp in dcp.items():

            ql['Taxonomy'].append(
                (dom,
                 base_tax_url % (version, version, cp, version),
                 'FeatureData[Taxonomy]',
                 'HeaderlessTSVTaxonomyFormat'))

            ql['Sequence'].append(
               (dom,
                base_seq_url % (version, version, cp, version),
                'FeatureData[Sequence]',
                'DNAFASTAFormat'))
    return ql


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
                # seq files are contained within `tar.gz``
                if tarfile.is_tarfile(destination):
                    try:
                        untarred_fn = bn.split('.')[0]+'.fna'
                        print('  Untarring {0}...\n'.format(bn))
                        with tarfile.open(destination, 'r') as tar:
                            tar.extract(member=untarred_fn,
                                        path=tmpdirname)
                    except OSError:
                        pass
                    seq_results.append(qiime2.Artifact.import_data(
                                        dtype,
                                        os.path.join(tmpdirname, untarred_fn),
                                        fmt))
                # taxonomy files are contained within ``.gz`
                else:
                    try:
                        unzipped_destination = os.path.splitext(destination)[0]
                        print('  Unzipping {0}...\n'.format(bn))
                        _gzip_decompress(destination, unzipped_destination)
                        destination = unzipped_destination
                    except OSError:
                        pass
                    tax_results.append(qiime2.Artifact.import_data(
                                                                dtype,
                                                                destination,
                                                                fmt))
        return tax_results, seq_results
