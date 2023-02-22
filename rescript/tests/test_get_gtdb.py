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
from rescript.get_gtdb import (VERSION_DICT, _assemble_queries,
                               _retrieve_data_from_gtdb)

from urllib.request import urlopen
from urllib.error import HTTPError
# from unittest.mock import patch


class TestGetGTDB(TestPluginBase):
    package = 'rescript.tests'

    # test that appropriate URLs are assembled
    def test_assemble_queries(self):
        queries = _assemble_queries(VERSION_DICT)

        obs_tax_urls = [q_info[1] for q_info in queries['Taxonomy']]
        obs_seq_urls = [q_info[1] for q_info in queries['Sequence']]

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

        # test that these URLs work
        for u in obs_tax_urls + obs_seq_urls:
            try:
                urlopen(u)
            except HTTPError:
                raise ValueError('Failed to open URL: ' + u)

    def test_retrieve_data_from_gtdb(self):
        queries = {'Taxonomy': [('Archaea', (
                        'https://data.gtdb.ecogenomic.org/releases/'
                        'release207/207.0/ar53_taxonomy_r207.tsv.gz'),
                        'FeatureData[Taxonomy]',
                        'HeaderlessTSVTaxonomyFormat')],
                   'Sequence': [('Archaea', (
                        'https://data.gtdb.ecogenomic.org/releases'
                        '/release207/207.0/genomic_files_reps/'
                        'ar53_ssu_reps_r207.tar.gz'),
                        'FeatureData[Sequence]', 'DNAFASTAFormat')]}
        _retrieve_data_from_gtdb(queries)
        self.assertTrue(True)  # if we make it here, it worked.
