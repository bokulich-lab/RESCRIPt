# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import pkg_resources
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_data import (_assemble_silva_data_urls,
                               _retrieve_data_from_silva)
from rescript.plugin_setup import _SILVA_VERSIONS, _SILVA_TARGETS
from urllib.request import urlopen
from urllib.error import HTTPError
from unittest.mock import patch


class TestGetSILVA(TestPluginBase):
    package = 'rescript.tests'

    # test that appropriate URLs are assembled, and those URLs work
    def test_assemble_silva_data_urls(self):
        for version in _SILVA_VERSIONS:
            for target in _SILVA_TARGETS:
                # don't test incompatible version/target settings (these are
                # prevented in the user interface by TypeMap)
                if target == 'LSURef' and version == '138':
                    continue
                obs = _assemble_silva_data_urls(version, target)
                # validate URLs
                for _, u, _ in obs:
                    try:
                        urlopen(u)
                    except HTTPError:
                        raise ValueError('Failed to open URL: ' + u)

    def test_retrieve_data_from_silva(self):
        # we just test download with a single small query here
        # this test just makes sure the method can download and validate
        # we do not check the outputs, since a successful return implies
        # that the contents are valid and imported successfully.
        queries = [
            ('taxa', 'https://www.arb-silva.de/fileadmin/silva_databases/'
                     'release_138/Exports/taxonomy/tax_slv_ssu_138.tre.gz',
             'Phylogeny[Rooted]')]
        _retrieve_data_from_silva(queries)
        self.assertTrue(True)

    # This tests the full get_silva_data pipeline, using mock data and
    # skipping the seqs for the sake of time. All relevant internals are
    # tested elsewhere in this test class, or in TestGetSILVA below, so this
    # just ensures that the full pipeline operates seamlessly, using `mock`
    # to mock data download and slip in fake data in its place.
    def test_get_silva_data(self):

        def _fake_data_on_demand(give_me_anything_i_shall_ignore_it):
            tr = qiime2.Artifact.import_data(
                'FeatureData[SILVATaxonomy]',
                pkg_resources.resource_filename(
                    'rescript.types.tests', 'data/silva_taxa.tsv'))
            tt = qiime2.Artifact.import_data(
                'Phylogeny[Rooted]', self.get_data_path('taxid_tree.tre'))
            tm2 = qiime2.Artifact.import_data(
                'FeatureData[SILVATaxidMap]',
                self.get_data_path('taxmap_test_match_tree.txt'))
            fake_dict = {'taxonomy tree': tt,
                         'taxonomy map': tm2,
                         'taxonomy ranks': tr}
            return fake_dict

        with patch('rescript.get_data._retrieve_data_from_silva',
                   new=_fake_data_on_demand):
            rescript.actions.get_silva_data(
                version='132', target='SSURef_NR99', download_sequences=False)
            self.assertTrue(True)
