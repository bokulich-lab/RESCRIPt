# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import tempfile
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_unite import (
    _unite_get_doi,
    _unite_get_url,
    _unite_get_tgz,
    _unite_get_artifacts,
    get_unite_data
)

from urllib.request import urlopen
from urllib.error import HTTPError
from unittest.mock import patch

# Global vars for use between functions
obs_dois = []


class TestGetUNITE(TestPluginBase):
    package = "rescript.tests"

    def setUp(self):
        super().setUp()
        self.unitefile = pkg_resources.resource_filename(
            "rescript.tests", "data/unitefile.tgz"
        )

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
    def test_unite_get_tgz(self):
        # Use a single, small, and unrelated file for testing
        # url = _unite_get_url('10.15156/BIO/587477')
        url = "https://files.plutof.ut.ee/doi/C9/F6/C9F687C997F72F674AA539CB80BF5D5BF6D1F402A2ACF840B20322850D3DFBA4.zip"  # noqa E501
        with tempfile.TemporaryDirectory() as tmpdirname:
            _unite_get_tgz(url, tmpdirname)

    def test_unite_get_artifacts(self):
        # Test on small data/unitefile.tgz with two items inside
        res_one, res_two = _unite_get_artifacts(
            self.unitefile, cluster_id="97"
        )
        self.assertEqual(str(res_one[0].type), "FeatureData[Taxonomy]")
        self.assertEqual(str(res_two[0].type), "FeatureData[Sequence]")
        # test missing files or misspelled cluster_id
        with self.assertRaises(ValueError):
            _unite_get_artifacts(self.unitefile, "nothing")

    # This tests the full get_unite_data pipeline with toy data.
    # All relevant internals are tested elsewhere in this test class, so
    # this just ensures that the full pipeline works.
    # Downloading is mock`ed with patch.
    def test_get_unite_data(self):
        # def _toy_data(give_me_anything_i_shall_ignore_it):
        #     res_one, res_two = _unite_get_artifacts(
        #         self.unitefile, cluster_id="97"
        #     )
        #     return res_one, res_two
        # default
        with patch("rescript.get_unite._unite_get_tgz", new=self.unitefile):
            rescript.actions.get_unite_data(
                version='8.3', taxon_group='fungi', cluster_id = "97"
            )
            self.assertTrue(True)

    def test_get_unite_data2(self):
        # def _toy_data(give_me_anything_i_shall_ignore_it):
        #     res_one, res_two = _unite_get_artifacts(
        #         self.unitefile, cluster_id="97"
        #     )
        #     return res_one, res_two
        # default
        with patch("rescript.get_unite._unite_get_tgz", new=self.unitefile):
            get_unite_data(version="8.3", taxon_group="fungi", cluster_id="97")
            self.assertTrue(True)