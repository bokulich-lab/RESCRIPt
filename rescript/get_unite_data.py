# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import warnings
import requests
import tarfile

import qiime2
from urllib.request import urlretrieve
from urllib.error import HTTPError
from q2_types.feature_data import DNAFASTAFormat


def _unite_get_doi(version, taxon_group, singletons):
    '''Lookup UNITE DOIs from included list'''
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
        doi = unite_dois[version][taxon_group][singletons]
    except KeyError as ke:
        print('Unknown DOI for this value: ' + str(ke))
        raise
    return doi

# Testing
_unite_get_doi(version='9.0', taxon_group='fungi', singletons=False)

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
_unite_get_url(doi = '10.15156/BIO/2938079')

# with tempfile.TemporaryDirectory() as tmp_dir:
tmp_dir = tempfile.mkdtemp()


def _unite_get_tgz(url, download_path):
    print('Downloading ' + url)
    response = requests.get(url, stream=True)
    if response.status_code != 200:
        raise ValueError("Failed to download the file from " + url)
    # Save .tgz file
    tar_file_path = os.path.join(download_path, 'unitefile.tar.gz')
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
            tar.extract(member, path=download_path)
    
    return download_path

# Test it by downloading this file
_unite_get_tgz(_unite_get_url(doi = '10.15156/BIO/2938079'), tmp_dir)


def get_unite_data(version, taxon_group, singletons=False):
    doi = _unite_get_doi(version, taxon_group, singletons)
    url = _unite_get_url(doi)
    # Eventual output
    results = {'sequences': [], 'taxonomy': []}
    # open temp output files, somehow?
    with tempfile.TemporaryDirectory() as tmpdirname:
        # tmp_dir = tempfile.mkdtemp()
        print('Temporary directory:', tmpdirname)
        _unite_get_tgz(url, download_path=tmpdirname)
        # Find and import artifacts based on suffix
        for root, dirs, files in os.walk(tmpdirname):
            for file in files:
                print(results)
                if file.endswith('.fasta'):
                    fasta_file_name = os.path.join(root, file)
                    print('found fasta: ' + fasta_file_name)
                    results['sequences'].append(qiime2.Artifact.import_data('FeatureData[Sequence]', fasta_file_name))
                elif file.endswith('.txt'):
                    txt_file_name = os.path.join(root, file)
                    print('found txt: ' + txt_file_name)
                    results['taxonomy'].append(qiime2.Artifact.import_data('FeatureData[Taxonomy]', txt_file_name))
    return results

# Testing
example_output = get_unite_data(version='9.0', taxon_group='fungi')

example_output['sequences'][1]
