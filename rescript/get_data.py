# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import hashlib
import shutil
import gzip
import warnings

import qiime2
from urllib.request import urlretrieve
from urllib.error import HTTPError
from .types._format import RNAFASTAFormat


def get_silva_data(ctx,
                   version='138',
                   target='SSURef_NR99',
                   include_species_labels=False,
                   rank_propagation=True,
                   download_sequences=True):
    # download data from SILVA
    print('Downloading raw files may take some time... get some coffee.')
    queries = _assemble_silva_data_urls(version, target, download_sequences)
    results = _retrieve_data_from_silva(queries)
    # parse taxonomy
    parse_taxonomy = ctx.get_action('rescript', 'parse_silva_taxonomy')
    taxonomy, = parse_taxonomy(
        taxonomy_tree=results['taxonomy tree'],
        taxonomy_map=results['taxonomy map'],
        taxonomy_ranks=results['taxonomy ranks'],
        include_species_labels=include_species_labels,
        rank_propagation=rank_propagation)
    # if skipping sequences, need to output an empty sequence file.
    if not download_sequences:
        results['sequences'] = qiime2.Artifact.import_data(
            'FeatureData[RNASequence]', RNAFASTAFormat())
    return results['sequences'], taxonomy


def _assemble_silva_data_urls(version, target, download_sequences=True):
    '''Generate SILVA urls, given database version and reference target.'''
    # assemble target urls
    ref_map = {'SSURef_NR99': 'ssu_ref_nr',
               'SSURef_Nr99': 'ssu_ref_nr',
               'SSURef': 'ssu_ref',
               'LSURef': 'lsu_ref'}
    # handle silly inconsistencies in filenames between versions
    if target == 'SSURef_NR99' and int(version) < 138:
        target = 'SSURef_Nr99'
    insert = ref_map[target]
    # now compile URLs
    base_url = 'https://www.arb-silva.de/fileadmin/silva_databases/'\
               'release_{0}/Exports/'.format(version)
    base_url_seqs = base_url + 'SILVA_{0}_{1}_tax_silva.fasta.gz'.format(
        version, target)
    base_url_taxmap = '{0}taxonomy/taxmap_slv_{1}_{2}'.format(
        base_url, insert, version)
    # More SILVA release inconsistencies
    if target == 'SSURef' and version == '132':
        base_url_taxmap += '-corrected.txt.gz'
    else:
        base_url_taxmap += '.txt.gz'
    base_url_tax = '{0}taxonomy/tax_slv_{1}_{2}'.format(
        base_url, insert.split('_')[0], version)
    tree_url = base_url_tax + '.tre'
    tax_url = base_url_tax + '.txt'
    if version == '138':
        tree_url += '.gz'
        tax_url += '.gz'

    # download and validate silva files
    queries = [('sequences', base_url_seqs, 'FeatureData[RNASequence]'),
               ('taxonomy map', base_url_taxmap, 'FeatureData[SILVATaxidMap]'),
               ('taxonomy tree', tree_url, 'Phylogeny[Rooted]'),
               ('taxonomy ranks', tax_url, 'FeatureData[SILVATaxonomy]')]

    # optionally skip downloading sequences
    if not download_sequences:
        queries = queries[1:]

    return queries


def _retrieve_data_from_silva(queries):
    '''
    Download data from SILVA, given a list of queries.

    queries: list of tuples of (str, str, str)
        (name, urlpath, QIIME 2 artifact type)
    '''
    results = dict()
    with tempfile.TemporaryDirectory() as tmpdirname:
        for name, query, dtype in queries:
            print('retrieving {0} from: {1}'.format(name, query))
            # grab url
            destination = os.path.join(tmpdirname, os.path.basename(query))
            urlretrieve(query, destination)
            file_md5 = _get_md5(destination)
            # grab expected md5
            # NOTE: SILVA is missing md5s for some files, so we will just skip
            try:
                md5_destination = os.path.join(tmpdirname, 'md5')
                urlretrieve(query + '.md5', md5_destination)
                exp_md5 = _read_silva_md5(md5_destination)
                # validate md5 checksum
                _validate_md5(exp_md5, file_md5, query)
            # if we get an HTTPError just move on, the md5 file does not exist
            except HTTPError:
                msg = ("No md5 file was detected in the SILVA archive for the "
                       "following file. No action is required, but be aware "
                       "that md5-hash validation was not performed for this "
                       "file: " + query)
                warnings.warn(msg, UserWarning)
            # gunzip on demand (SILVA releases are inconsistently gzipped)
            try:
                unzipped_destination = os.path.splitext(destination)[0]
                _gzip_decompress(destination, unzipped_destination)
                destination = unzipped_destination
            except OSError:
                pass
            # import as artifacts
            results[name] = qiime2.Artifact.import_data(dtype, destination)
    return results


def _validate_md5(exp_md5, file_md5, filename):
    if not exp_md5 == file_md5:
        raise ValueError('md5 sums do not match. Manually verify md5 sums '
                         'before proceeding.\nTarget file: {0}\nExpected md5: '
                         '{1}\nObserved md5: {2}\n'.format(
                            filename, exp_md5, file_md5))


# This function is specific for reading the SILVA md5 record files, which are
# a single line txt file in the format "md5  filename"
def _read_silva_md5(file):
    with open(file, 'r') as _md5:
        _md5 = _md5.read().split(' ')[0]
    return _md5


def _get_md5(file, chunksize=8192):
    md5_hash = hashlib.md5()
    with open(file, "rb") as f:
        for chunk in iter(lambda: f.read(chunksize), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def _gzip_decompress(input_fp, output_fp):
    with gzip.open(input_fp, 'rt') as temp_in:
        with open(output_fp, 'w') as temp_out:
            shutil.copyfileobj(temp_in, temp_out)
