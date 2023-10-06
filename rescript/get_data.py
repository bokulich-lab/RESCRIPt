# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
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
import requests
import tarfile

import qiime2
from urllib.request import urlretrieve
from urllib.error import HTTPError
from q2_types.feature_data import RNAFASTAFormat


def _unite_dois_to_urls(DOIs):
    '''Generate UNITE urls, given their DOIs.'''
    # Make DOIs iterable
    DOIs = [DOIs] if isinstance(DOIs, str) else DOIs
    print('Get URLs for these DOIs:', DOIs)
    base_url = 'https://api.plutof.ut.ee/v1/public/dois/'\
               '?format=vnd.api%2Bjson&identifier='
    # Eventual output
    URLs = set()
    # For each DOI, get download URL of file
    for DOI in DOIs:
        query_data = requests.get(base_url + DOI).json()
        # Updates can be made to files in a DOI, so on the advice of the devs,
        # only return the last file uploaded with this -1  vv
        URL = query_data['data'][0]['attributes']['media'][-1]['url']
        URLs.add(URL)
    return URLs


def _unite_get_url(version, taxon_group, singletons=False):
    '''Generate UNITE urls, given database version and reference target.'''
    # Lookup DOIs for various databases, source: https://unite.ut.ee/repository.php
    unite_dois = {
        '9.0': {'fungi':      {False: '10.15156/BIO/2938079', True: '10.15156/BIO/2938080'},
                'eukaryotes': {False: '10.15156/BIO/2938081', True: '10.15156/BIO/2938082'}},
        # Old version 9.0 is not listed here
        '8.3': {'fungi':      {False: '10.15156/BIO/1264708', True: '10.15156/BIO/1264763'},
                'eukaryotes': {False: '10.15156/BIO/1264819', True: '10.15156/BIO/1264861'}},
        '8.2': {'fungi':      {False: '10.15156/BIO/786385', True: '10.15156/BIO/786387'},
                'eukaryotes': {False: '10.15156/BIO/786386', True: '10.15156/BIO/786388'}},
        '8.0': {'fungi':      {False: '', True: '10.15156/BIO/786349'},
                'eukaryotes': {False: '', True: ''}},  # All other 8.0 are in zip files
    }
    # There's got to be a better way! See https://stackoverflow.com/questions/25833613/safe-method-to-get-value-of-nested-dictionary
    try:
        # Check if we have the DOI requested
        target_doi = unite_dois[version][taxon_group][singletons]
    except KeyError as ke:
        print('Unknown DOI for this value: ' + str(ke))
        raise
    return _unite_dois_to_urls(target_doi)

_unite_get_url(version='9.0', taxon_group='fungi')

def _unite_download_targz(url):
    print('Downloading ' + url)
    
    # with tempfile.TemporaryDirectory() as tmp_dir:
    tmp_dir = tempfile.mkdtemp()
    
    response = requests.get(url, stream=True)
    if response.status_code != 200:
        raise ValueError("Failed to download the file from " + url)
    
    tar_file_path = os.path.join(tmp_dir, 'unitefile.tar.gz')
    with open(tar_file_path, 'wb') as f:
        f.write(response.content)
    
    # Extract only the 'developer' subdirectory
    with tarfile.open(tar_file_path, 'r:gz') as tar:
        # Ensure that 'developer' exists in the tar file
        members = [member for member in tar.getmembers() if member.name.startswith('developer')]
        if not members:
            raise ValueError("No 'developer' subdirectory found in the .tar.gz file.")
        
        for member in members:
            member.name = os.path.basename(member.name) # Strip the 'developer' prefix
            tar.extract(member, path=tmp_dir)
    
    return tmp_dir

# Test it by downloading this file
_unite_download_targz('https://files.plutof.ut.ee/public/orig/59/12/591225E8985EFC44B595C79AF5F467421B4D9A95093A0811B13CB4CC13A6DA46.tgz')

    # import as artifacts
    results[name] = qiime2.Artifact.import_data(dtype, destination)



def get_silva_data(ctx,
                   version='138.1',
                   target='SSURef_NR99',
                   include_species_labels=False,
                   rank_propagation=True,
                   ranks=None,
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
        ranks=ranks,
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
               'LSURef_NR99': 'lsu_ref_nr',
               'LSURef': 'lsu_ref'}
    # handle silly inconsistencies in filenames between versions
    if target == 'SSURef_NR99' and float(version) < 138:
        target = 'SSURef_Nr99'
    insert = ref_map[target]

    # Now compile URLs
    base_url = 'https://www.arb-silva.de/fileadmin/silva_databases/'\
               'release_{0}/Exports/'.format(version.replace('.', '_'))
    # ^^ Note the '.format()' above. This handles folder / file naming
    # inconsistency. e.g. folder name for version 138.1 contains the string
    # "138_1", while the filenames within that directory contain the string
    # "138.1". We'll need to modify the variable `version` and replace '.'
    # with '_' for the folder name. May need to update this in the future,
    # if we find more inconsistencies.

    # construct file urls
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

    # add ".gz" for the following versions:
    if version in ['138', '138.1']:
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
