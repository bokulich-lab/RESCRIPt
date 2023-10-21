# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import pkg_resources
import tempfile
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_unite import (
    _unite_get_doi,
    _unite_get_url,
    _unite_get_tgz,
    _unite_get_artifacts,
)

from urllib.request import urlopen
from urllib.error import HTTPError
from unittest.mock import patch

# Global vars for use between functions
obs_dois = []


class TestGetUNITE(TestPluginBase):
    package = "rescript.tests"

    def test_unite_get_doi(self):
        # All supported inputs to _unite_get_doi
        versions = ["9.0", "8.3", "8.2"]
        taxon_groups = ["fungi", "eukaryotes"]
        singletons = [True, False]
        # Eventual output
        obs_dois = []
        # for all combinations
        for v in versions:
            for tg in taxon_groups:
                for s in singletons:
                    # add it to the list
                    obs_dois.append(_unite_get_doi(v, tg, s))
        exp_dois = [
            "10.15156/BIO/2938080",
            "10.15156/BIO/2938079",
            "10.15156/BIO/2938082",
            "10.15156/BIO/2938081",
            "10.15156/BIO/1264763",
            "10.15156/BIO/1264708",
            "10.15156/BIO/1264861",
            "10.15156/BIO/1264819",
            "10.15156/BIO/786387",
            "10.15156/BIO/786385",
            "10.15156/BIO/786388",
            "10.15156/BIO/786386",
        ]
        self.assertEqual(obs_dois, exp_dois)

    # Depends on _unite_get_doi
    def test_unite_get_url(self):
        # All supported inputs to _unite_get_doi
        versions = ["9.0", "8.3", "8.2"]
        taxon_groups = ["fungi", "eukaryotes"]
        singletons = [True, False]
        # for all combinations
        for v in versions:
            for tg in taxon_groups:
                for s in singletons:
                    # First we lookup DOIs from the table
                    doi = _unite_get_doi(v, tg, s)
                    # Then we get and try to open their url
                    try:
                        url = _unite_get_url(doi)
                        urlopen(url)
                    except HTTPError:
                        raise ValueError("Failed to open URL: " + url)

    ##########
    # Alternative versions of these two functions. Instead of using a
    # hardcoded list of DOIs, they populate a global obs_dois list. 
    ##########

    # test by popbulating obs_dois list
    def test_unite_get_doi2(self):
        # All supported inputs to _unite_get_doi
        versions = ["9.0", "8.3", "8.2"]
        taxon_groups = ["fungi", "eukaryotes"]
        singletons = [True, False]
        # for all combinations
        for v in versions:
            for tg in taxon_groups:
                for s in singletons:
                    # add it to the list
                    obs_dois.append(_unite_get_doi(v, tg, s))
        if obs_dois:
            self.assertTrue(True)

    # Test by querying the newly full obs_dois list
    def test_unite_get_url2(self):
        for doi in obs_dois:
            try:
                url = _unite_get_url(doi)
                urlopen(url)
            except HTTPError:
                raise ValueError("Failed to open URL: " + url)

    # Download a single, small file from the API endpoint
    def test_unite_get_raw_files(self):
        # Use a single (unrelated) download URL
        url = "https://files.plutof.ut.ee/doi/C8/E4/C8E4A8E6A7C4C00EACE3499C51E550744A259A98F8FE25993B1C7B9E7D2170B2.zip" # noqa
        with tempfile.TemporaryDirectory() as tmpdirname:
            _unite_get_tgz(url, tmpdirname)
            # If function finished, pass this test
            self.assertTrue(True)

    def test_unite_get_artifacts(self):
        # with patch()
        self.assertTrue(True)
