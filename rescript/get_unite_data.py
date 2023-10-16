# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import requests
import tarfile

import qiime2


def _unite_get_doi(version, taxon_group, singletons):
    '''Lookup UNITE DOIs from included list'''
    # Lookup DOIs for databases, see: https://unite.ut.ee/repository.php
    unite_dois = {
        '9.0': {'fungi': {False: '10.15156/BIO/2938079',
                          True: '10.15156/BIO/2938080'},
                'eukaryotes': {False: '10.15156/BIO/2938081',
                               True: '10.15156/BIO/2938082'}},
        # Old version 9.0 is not listed here
        '8.3': {'fungi': {False: '10.15156/BIO/1264708',
                          True: '10.15156/BIO/1264763'},
                'eukaryotes': {False: '10.15156/BIO/1264819',
                               True: '10.15156/BIO/1264861'}},
        '8.2': {'fungi': {False: '10.15156/BIO/786385',
                          True: '10.15156/BIO/786387'},
                'eukaryotes': {False: '10.15156/BIO/786386',
                               True: '10.15156/BIO/786388'}},
        '8.0': {'fungi': {False: '',
                          # All other 8.0 are in zip files
                          True: '10.15156/BIO/786349'},
                'eukaryotes': {False: '',
                               True: ''}},
    }
    try:
        # Check if we have the DOI requested
        doi = unite_dois[version][taxon_group][singletons]
    except KeyError as ke:
        print('Unknown DOI for this value: ' + str(ke))
        raise
    return doi


# Testing
# _unite_get_doi(version='9.0', taxon_group='fungi', singletons=False)


def _unite_get_url(doi):
    '''Query plutof API for UNITE url'''
    print('Get URLs for these DOIs:', doi)
    base_url = 'https://api.plutof.ut.ee/v1/public/dois/'\
               '?format=vnd.api%2Bjson&identifier='
    query_data = requests.get(base_url + doi).json()
    # Updates can be made to files in a DOI, so on the advice of the devs,
    # only return the last (newest) file with this -1  vv
    URL = query_data['data'][0]['attributes']['media'][-1]['url']
    return URL


# Testing
# _unite_get_url(doi='10.15156/BIO/2938079')


# Make tmp_dir for standalone testing
# tmp_dir = tempfile.mkdtemp()
def _unite_get_raw_files(url, download_path):
    '''Download and extract all fasta and txt files'''
    response = requests.get(url, stream=True)
    if response.status_code != 200:
        raise ValueError("Failed to download the file from " + url)
    # Save .tgz file
    unite_file_path = os.path.join(download_path, 'unitefile.tar.gz')
    with open(unite_file_path, 'wb') as f:
        f.write(response.content)
    # Extract only the 'developer' subdirectory
    with tarfile.open(unite_file_path, 'r:gz') as tar:
        # Ensure that 'developer' exists in the tar file
        members = [member for member in tar.getmembers()
                   if member.name.startswith('developer')]
        if not members:
            raise ValueError("No 'developer' subdirectory found")
        for member in members:
            # Strip the 'developer' prefix
            member.name = os.path.basename(member.name)
            tar.extract(member, path=download_path)
    return


# Test it by downloading this file
# _unite_get_raw_files(_unite_get_url(doi='10.15156/BIO/2938079'), tmp_dir)


def _unite_get_qza(cluster_id, download_path):
    '''
    Find and import raw files with matching cluster_id

    Returns: Tuple containing tax_results and seq_results
    '''
    tax_results = []
    seq_results = []
    # Find all files...
    for root, dirs, files in os.walk(download_path):
        # ... with the matching cluster_id
        filtered_files = [file for file in files if cluster_id in file]
        for file in filtered_files:
            fp = os.path.join(root, file)
            if file.endswith('.txt'):
                tax_results.append(
                    qiime2.Artifact.import_data('FeatureData[Taxonomy]', fp))
            elif file.endswith('.fasta'):
                seq_results.append(
                    qiime2.Artifact.import_data('FeatureData[Sequence]', fp))
    return tax_results, seq_results


# Testing
# _unite_get_qza('99', tmp_dir)


def get_unite_data(version, taxon_group, cluster_id='99', singletons=False):
    '''
    Get Qiime2 artifacts for a given version of UNITE

    Returns: Tuple containing tax_results and seq_results
    '''
    doi = _unite_get_doi(version, taxon_group, singletons)
    url = _unite_get_url(doi)
    with tempfile.TemporaryDirectory() as tmpdirname:
        print('Temporary directory:', tmpdirname)
        _unite_get_raw_files(url, tmpdirname)
        return _unite_get_qza(cluster_id, tmpdirname)


# Testing
# example_output=get_unite_data(version='9.0', taxon_group='fungi')
# example_output
