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
from rescript.get_unite import (
    UNITE_DOIS,
    _unite_get_url,
    _unite_get_tgz,
    _unite_get_artifacts,
    get_unite_data,
)

from urllib.request import urlopen
from urllib.error import HTTPError
from unittest.mock import patch


class TestGetUNITE(TestPluginBase):
    package = "rescript.tests"

    def setUp(self):
        super().setUp()
        self.unitefile = pkg_resources.resource_filename(
            "rescript.tests", "data/unitefile.tgz"
        )

    # Requires internet access
    def test_unite_get_url(self):
        # for all combinations...
        for v in UNITE_DOIS.keys():
            for tg in UNITE_DOIS[v].keys():
                for s in UNITE_DOIS[v][tg].keys():
                    # ... try to get the URL
                    try:
                        url = _unite_get_url(v, tg, s)
                        urlopen(url)
                    except HTTPError:
                        raise ValueError("No URL for combo: " + v + tg + s)

    # Requires internet access
    def test_unite_get_tgz(self):
        # Download a single, small, unrelated file for testing
        url = "https://files.plutof.ut.ee/doi/C9/F6/C9F687C997F72F674AA539CB80BF5D5BF6D1F402A2ACF840B20322850D3DFBA4.zip"  # noqa E501
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Test working URL
            _unite_get_tgz(url, tmpdirname)

    def test_unite_get_artifacts(self):
        # Test on small data/unitefile.tgz with two items inside
        res_one, res_two = _unite_get_artifacts(
            self.unitefile, cluster_id="97"
        )
        # Column names and one feature from TaxonomyFormat
        self.assertEqual(
            res_one["Taxon"]["SH1140752.08FU_UDB013072_reps"],
            "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Thelephorales;"
            "f__Thelephoraceae;g__Tomentella;s__unidentified",
        )
        self.assertEqual(
            str(type(res_two)),
            "<class 'q2_types.feature_data._transformer.DNAIterator'>",
        )
        # test missing files or misspelled cluster_id
        with self.assertRaises(ValueError):
            _unite_get_artifacts(self.unitefile, "nothing")

    # This tests the function with toy data.
    # All relevant internals are tested elsewhere in this test class.
    # Downloading is mock`ed with patch.
    def test_get_unite_data(self):
        with patch(
            "rescript.get_unite._unite_get_tgz", return_value=self.unitefile
        ):
            get_unite_data(version="8.3", taxon_group="fungi", cluster_id="97")
            self.assertTrue(True)
