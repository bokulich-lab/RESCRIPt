# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase

from rescript.bv_brc import fetch_genomes_bv_brc


class TestPipelines(TestPluginBase):
    package = 'rescript.tests'

    def test_fetch_genomes_bv_brc(self):
        query = "?eq(genome_id,224308.43)"
        query2 = "?eq(taxon_id,224308)"
        fetch_genomes_bv_brc(query2)