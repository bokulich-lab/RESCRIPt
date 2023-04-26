# ----------------------------------------------------------------------------
# Copyright (c) 2022-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import shutil
import gzip
# import warnings
# import hashlib
# import pandas as pd
import qiime2
from urllib.request import urlretrieve
# from urllib.error import HTTPError
# from q2_types.feature_data import MixedCaseDNAFASTAFormat, DNAFASTAFormat
# from q2_feature_table import merge_seqs, merge_taxa


def get_rdp_data(ctx,
                 ref_db='Archaea+Bacteria',
                 include_species_labels=False,
                 rank_propagation=True,
                 ranks=None):

    parse_taxonomy = ctx.get_action('rescript', 'parse_rdp_taxonomy')
    merge_rdp_seqs = ctx.get_action('feature_table', 'merge_seqs')
    merge_rdp_taxonomy = ctx.get_action('feature_table', 'merge_taxa')

    # download data from RDP
    print('\nDownloading and processing raw files may take some time... '
          'get some coffee.\n')
    queries = _assemble_rdp_data_urls(ref_db)
    rdp_seq_data = _retrieve_data_from_rdp(queries)

    num_seq_files = len(rdp_seq_data)

    if num_seq_files == 2:
        rdp_tax_df = []
        rdp_seq_sl = []
        for rdp_seqs in rdp_seq_data:
            print('Parsing RDP taxonomy...')
            taxonomy, = parse_taxonomy(
                            rdp_reference_sequences=rdp_seqs,
                            include_species_labels=include_species_labels,
                            rank_propagation=rank_propagation,
                            ranks=ranks)
            # rdp_tax_df.append(taxonomy.view(pd.DataFrame))
            # rdp_seq_sl.append(rdp_seqs.view(pd.Series))
            rdp_tax_df.append(taxonomy)
            rdp_seq_sl.append(rdp_seqs)
        print('Merging taxonomy ...')
        merged_rdp_tax, = merge_rdp_taxonomy(data=rdp_tax_df)
        print('Merging sequences ...')
        merged_rdp_seqs, = merge_rdp_seqs(data=rdp_seq_sl)
        print('Writing artifact files...')
        return merged_rdp_seqs, merged_rdp_tax
    elif num_seq_files == 1:
        rdp_seqs = rdp_seq_data[0]
        rdp_tax, = parse_taxonomy(
                        rdp_reference_sequences=rdp_seqs,
                        include_species_labels=include_species_labels,
                        rank_propagation=rank_propagation,
                        ranks=ranks)
        print('Writing artifact files...')
        return rdp_seqs, rdp_tax
    elif num_seq_files == 0:
        raise ValueError('No RDP data to process!')
    else:
        raise ValueError('There are more than 2 RDP files! There should '
                         'be only 1 or 2 files!')
    # return rdp_seqs, rdp_tax


def _assemble_rdp_data_urls(ref_db):
    '''Generate RDP urls, given database version and reference target.'''
    # Compile URLs
    base_url = 'http://rdp.cme.msu.edu/download/'

    if ref_db == 'Archaea+Bacteria':
        archaea_url = base_url + 'current_Archaea_unaligned.fa.gz'
        bacteria_url = base_url + 'current_Bacteria_unaligned.fa.gz'
        return [('archaeal sequences', archaea_url),
                ('bacterial sequences', bacteria_url)]
    elif ref_db == 'Fungi':
        fungi_url = base_url + 'current_Fungi_unaligned.fa.gz'
        return [('fungal sequences', fungi_url)]
    else:
        raise ValueError('Reference database must be either, '
                         '\'Bacteria+Archaea\' or \'Fungi\'.')


def _retrieve_data_from_rdp(queries):
    '''
    Download data from RDP, given a list of queries.

    queries: list of tuples of (str, str)
        (name, urlpath)
    '''
    results = []
    with tempfile.TemporaryDirectory() as tmpdirname:
        for name, query in queries:
            print('Retrieving {0} from: {1}'.format(name, query))
            # grab url
            destination = os.path.join(tmpdirname, os.path.basename(query))
            urlretrieve(query, destination)

            # gunzip on demand (RDP releases are gzipped)
            try:
                unzipped_destination = os.path.splitext(destination)[0]
                print('  Unzipping {0}...'.format(name))
                _gzip_decompress(destination, unzipped_destination)
                destination = unzipped_destination
            except OSError:
                pass

            # import rdp seqs
            print('Importing RDP sequence data...')
            results.append(qiime2.Artifact.import_data(
                                                'FeatureData[Sequence]',
                                                destination,
                                                'MixedCaseDNAFASTAFormat'))
        print(results)
        return results


def _gzip_decompress(input_fp, output_fp):
    with gzip.open(input_fp, 'rt') as temp_in:
        with open(output_fp, 'w') as temp_out:
            shutil.copyfileobj(temp_in, temp_out)
