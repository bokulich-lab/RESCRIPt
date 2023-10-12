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


# import as artifacts
# results[name] = qiime2.Artifact.import_data(dtype, destination)

def get_unite_data(version, taxon_group, singletons=False):
    doi = _unite_get_doi(version, taxon_group, singletons)
    url = _unite_get_url(doi)
    # Eventual output
    results = {'sequences': [], 'taxonomy': []}
    # open temp output files, somehow?
    # with tempfile.TemporaryDirectory() as tmp_dir:
    tmp_dir = tempfile.mkdtemp()
    print('Temporary directory:', tmp_dir)
    _unite_get_tgz(url, download_path=tmp_dir)
    # Find and import artifacts based on suffix
    for root, dirs, files in os.walk(tmp_dir):
        for file in files:
            print(results)
            if file.endswith('.fasta'):
                fasta_file_name = os.path.join(root, file)
                print('found fasta: ' + fasta_file_name)
                with open(fasta_file_name, 'r') as fasta_file:
                    # Read the content of the file and append it as a Python object
                    fasta_content = fasta_file.read()
                    results['sequences'].append(qiime2.Artifact.import_data('FeatureData[RNASequence]', fasta_content))
            elif file.endswith('.txt'):
                txt_file_name = os.path.join(root, file)
                print('found txt: ' + txt_file_name)
                results['taxonomy'].append(qiime2.Artifact.import_data('FeatureData[Taxonomy]', txt_file_name))
    return results

get_unite_data(version='9.0', taxon_group='fungi')


# How do I import data?
# head /tmp/tmpduape95o/sh_refs_qiime_ver9_97_25.07.2023_dev.fasta
fasta_file = open("/tmp/tmpduape95o/sh_refs_qiime_ver9_97_25.07.2023_dev.fasta", 'r')

qiime2.Artifact.import_data('FeatureData[Sequence]', fasta_file)
qiime2.Artifact.import_data('FeatureData[Sequence]', str(fasta_file))
qiime2.Artifact.import_data('FeatureData[RNASequence]', fasta_file)
qiime2.Artifact.import_data('FeatureData[RNASequence]', str(fasta_file))
qiime2.Artifact.import_data('FeatureData[DNASequence]', fasta_file)
qiime2.Artifact.import_data('FeatureData[DNASequence]', str(fasta_file))

fasta_content = fasta_file.read()
qiime2.Artifact.import_data('FeatureData[Sequence]', str(fasta_content))
qiime2.Artifact.import_data('FeatureData[Sequence]', fasta_content)
qiime2.Artifact.import_data('FeatureData[RNASequence]', str(fasta_content))
qiime2.Artifact.import_data('FeatureData[RNASequence]', fasta_content)

qiime2.Artifact.import_data('FeatureData[RNASequence]', "/tmp/tmpduape95o/sh_refs_qiime_ver9_97_25.07.2023_dev.fasta")

with open("/tmp/tmpduape95o/sh_refs_qiime_ver9_97_25.07.2023_dev.fasta", 'r') as fasta_file:
    qiime2.Artifact.import_data('FeatureData[RNASequence]', str(fasta_file))



with open("/tmp/tmpduape95o/sh_taxonomy_qiime_ver9_97_25.07.2023_dev.txt", 'r') as file:
    qiime2.Artifact.import_data('FeatureData[Taxonomy]', file)

with open("/tmp/tmpduape95o/sh_taxonomy_qiime_ver9_97_25.07.2023_dev.txt", 'r') as file:
    # Read the content of the file and append it as a Python object
    txt_content = file.read()
    qiime2.Artifact.import_data('FeatureData[Taxonomy]', txt_content)

qiime2.Artifact.import_data('FeatureData[Taxonomy]', "/tmp/tmpduape95o/sh_taxonomy_qiime_ver9_97_25.07.2023_dev.txt")

