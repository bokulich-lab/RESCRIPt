# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import tarfile
import requests
from requests.exceptions import HTTPError

import qiime2

# Source: https://unite.ut.ee/repository.php
UNITE_DOIS = {
    "9.0": {
        "fungi": {
            False: "10.15156/BIO/2938079",
            True: "10.15156/BIO/2938080",
        },
        "eukaryotes": {
            False: "10.15156/BIO/2938081",
            True: "10.15156/BIO/2938082",
        },
    },
    # Old version 9.0 is not listed here
    "8.3": {
        "fungi": {
            False: "10.15156/BIO/1264708",
            True: "10.15156/BIO/1264763",
        },
        "eukaryotes": {
            False: "10.15156/BIO/1264819",
            True: "10.15156/BIO/1264861",
        },
    },
    "8.2": {
        "fungi": {
            False: "10.15156/BIO/786385",
            True: "10.15156/BIO/786387",
        },
        "eukaryotes": {
            False: "10.15156/BIO/786386",
            True: "10.15156/BIO/786388",
        },
    },
}


def _unite_get_url(
    version: str = None, taxon_group: str = None, singletons: bool = None
) -> str:
    """Get DOI from included list, then query plutof API for UNITE url"""
    # Get matching DOI
    doi = UNITE_DOIS[version][taxon_group][singletons]
    # Build URL
    base_url = (
        "https://api.plutof.ut.ee/v1/public/dois/"
        "?format=vnd.api%2Bjson&identifier="
    )
    query_data = requests.get(base_url + doi).json()
    # Updates can be made to files in a DOI, so on the advice of the devs,
    # only return the last (newest) file with this -1  vv
    URL = query_data["data"][0]["attributes"]["media"][-1]["url"]
    return URL


def _unite_get_tgz(
    url: str = None, download_path: str = None, retries: int = 10
) -> str:
    """Download compressed database"""
    for retry in range(retries):
        # Track downloaded size
        file_size = 0
        try:
            response = requests.get(url, stream=True)
            # Save .tgz file
            unite_file_path = os.path.join(download_path, "unitefile.tar.gz")
            with open(unite_file_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        file_size += len(chunk)
            # Check if the downloaded size matches the expected size
            if file_size == int(response.headers.get("content-length", 0)):
                return unite_file_path  # done!
            else:
                raise ValueError("File download incomplete")
        except HTTPError as e:
            print(
                "Request failed with code "
                + str(e.response.status_code)
                + ", on try "
                + str(retry)
            )
        except ValueError:
            print("File incomplete, on try " + str(retry))


def _unite_get_artifacts(
    tgz_file: str = None, cluster_id: str = None
) -> (None, None):
    """
    Find and import files with matching cluster_id from .tgz

    Returns: Tuple containing tax_results and seq_results
    """
    # Eventual output
    tax_results = []
    seq_results = []
    with tempfile.TemporaryDirectory() as tmpdirname:
        print("Temporary directory:", tmpdirname)
        # Extract from the .tgz file
        with tarfile.open(tgz_file, "r:gz") as tar:
            # Keep only _dev files
            members = [m for m in tar.getmembers() if "_dev" in m.name]
            if not members:
                raise ValueError("No '_dev' files found")
            for member in members:
                # Keep only base name
                member.name = os.path.basename(member.name)
                tar.extract(member, path=tmpdirname)
        # Find and import the raw files...
        for root, dirs, files in os.walk(tmpdirname):
            # ... with the matching cluster_id
            filtered_files = [
                f for f in files if f.split("_")[4] == cluster_id
            ]
            if not filtered_files:
                raise ValueError(
                    "No files found with cluster_id = " + cluster_id
                )
            for file in filtered_files:
                fp = os.path.join(root, file)
                if file.endswith(".txt"):
                    tax_results.append(
                        qiime2.Artifact.import_data(
                            "FeatureData[Taxonomy]",
                            fp,
                            "HeaderlessTSVTaxonomyFormat",
                        )
                    )
                elif file.endswith(".fasta"):
                    seq_results.append(
                        qiime2.Artifact.import_data(
                            "FeatureData[Sequence]",
                            fp,
                            "MixedCaseDNAFASTAFormat",
                        )
                    )
    return tax_results, seq_results


def get_unite_data(
    version: str = None,
    taxon_group: str = None,
    cluster_id: str = "99",
    singletons: bool = False,
) -> (None, None):
    """
    Get Qiime2 artifacts for a given version of UNITE

    Returns: Tuple containing tax_results and seq_results
    """
    url = _unite_get_url(version, taxon_group, singletons)
    with tempfile.TemporaryDirectory() as tmpdirname:
        print("Temporary directory:", tmpdirname)
        tar_file_path = _unite_get_tgz(url, tmpdirname)
        return _unite_get_artifacts(tar_file_path, cluster_id)
