# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# import qiime2
# import pkg_resources
from qiime2.plugin.testing import TestPluginBase
# from qiime2.plugins import rescript
from rescript.get_gtdb import (VERSION_DICT, _assemble_gtdb_data_urls)

# from urllib.request import urlopen
# from urllib.error import HTTPError
# from unittest.mock import patch


class TestGetGTDB(TestPluginBase):
    package = 'rescript.tests'

    # test that appropriate URLs are assembled, and those URLs work
    def test_assemble_gtdb_data_urls(self):
        obs_tax_queries, obs_seq_queries = _assemble_gtdb_data_urls(
                                                       VERSION_DICT)
        obs_tax_urls = [tax_q[1] for tax_q in obs_tax_queries]
        obs_seq_urls = [seq_q[1] for seq_q in obs_seq_queries]

        exp_tax_urls = [('https://data.gtdb.ecogenomic.org/releases/'
                         'release207/207.0/ar53_taxonomy_r207.tsv.gz'),
                        ('https://data.gtdb.ecogenomic.org/releases/'
                         'release207/207.0/bac120_taxonomy_r207.tsv.gz'),
                        ('https://data.gtdb.ecogenomic.org/releases/'
                         'release202/202.0/ar122_taxonomy_r202.tsv.gz'),
                        ('https://data.gtdb.ecogenomic.org/releases'
                         '/release202/202.0/bac120_taxonomy_r202.tsv.gz')]
        exp_seq_urls = [('https://data.gtdb.ecogenomic.org/releases'
                         '/release207/207.0/genomic_files_reps/'
                         'ar53_ssu_reps_r207.tar.gz'),
                        ('https://data.gtdb.ecogenomic.org/releases/'
                         'release207/207.0/genomic_files_reps/'
                         'bac120_ssu_reps_r207.tar.gz'),
                        ('https://data.gtdb.ecogenomic.org/releases/'
                         'release202/202.0/genomic_files_reps/'
                         'ar122_ssu_reps_r202.tar.gz'),
                        ('https://data.gtdb.ecogenomic.org/releases/'
                         'release202/202.0/genomic_files_reps/'
                         'bac120_ssu_reps_r202.tar.gz')]
        self.assertEqual(obs_tax_urls, exp_tax_urls)
        self.assertEqual(obs_seq_urls, exp_seq_urls)
