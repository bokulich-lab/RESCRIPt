# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
import tempfile
import pandas.core.frame
import q2_types.feature_data
from qiime2.plugin.testing import TestPluginBase
from rescript.get_unite import (
    UNITE_DOIS,
    _unite_get_url,
    _unite_get_tgz,
    _unite_get_artifacts,
    get_unite_data,
)

from urllib.request import urlopen
from unittest.mock import patch, Mock


class TestGetUNITE(TestPluginBase):
    package = "rescript.tests"

    def setUp(self):
        super().setUp()
        self.unitefile = pkg_resources.resource_filename(
            "rescript.tests", "data/unitefile.tgz"
        )
        self.unitefile_no_dev = pkg_resources.resource_filename(
            "rescript.tests", "data/unitefile_no_dev.tgz"
        )

    # Requires internet access
    def test_unite_get_url(self):
        # for all combinations...
        for v in UNITE_DOIS.keys():
            for tg in UNITE_DOIS[v].keys():
                for s in UNITE_DOIS[v][tg].keys():
                    # ... try to get the URL
                    url = _unite_get_url(v, tg, s)
                    urlopen(url)
        self.assertTrue(True)

    def test_unite_get_tgz(self):
        with tempfile.TemporaryDirectory() as tmpdirname:
            # mock the response object
            mock_response = Mock()
            mock_response.iter_content.return_value = [b"mock"]
            mock_response.headers.get.return_value = "4"  # matches content
            # mock successful download
            with patch("requests.get", return_value=mock_response):
                _unite_get_tgz("fakeURL", tmpdirname)
            # real failed download
            with self.assertRaisesRegex(ValueError, "File incomplete on try"):
                _unite_get_tgz("https://files.plutof.ut.ee/nope", tmpdirname)

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
        # test no _dev files found
        with self.assertRaises(ValueError):
            _unite_get_artifacts(self.unitefile_no_dev, cluster_id="97")
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
            res = get_unite_data(
                version="8.3", taxon_group="fungi", cluster_id="97"
            )
            self.assertEqual(len(res), 2)
            self.assertTrue(isinstance(res[0], pandas.core.frame.DataFrame))
            self.assertTrue(
                isinstance(
                    res[1], q2_types.feature_data._transformer.DNAIterator
                )
            )
