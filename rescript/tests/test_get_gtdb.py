# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pkg_resources
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from rescript.get_gtdb import _assemble_queries, parse_gtdb_taxonomy
from q2_types.feature_data import (TSVTaxonomyFormat, DNAFASTAFormat)

from urllib.request import urlopen
from urllib.error import HTTPError
from unittest.mock import patch


class TestGetGTDB(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.gtdb_tax = TSVTaxonomyFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/gtdb-tax.tsv'),
                            mode='r')
        self.gtdb_seqs = DNAFASTAFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/gtdb-seqs.fasta'),
                            mode='r')
        self.gtdb_arch_tax = TSVTaxonomyFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/gtdb-taxa-archaea.tsv'),
                            mode='r')
        self.gtdb_arch_seqs = DNAFASTAFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/gtdb-seqs-archaea.fasta'),
                            mode='r')
        self.gtdb_bact_tax = TSVTaxonomyFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/gtdb-taxa-bacteria.tsv'),
                            mode='r')
        self.gtdb_bact_seqs = DNAFASTAFormat(
                            pkg_resources.resource_filename(
                                 'rescript.tests',
                                 'data/gtdb-seqs-bacteria.fasta'),
                            mode='r')

    # test that appropriate URLs are assembled
    def test_assemble_species_rep_queries(self):
        # checking v220 as GTDB updated the file format for the "ssu_rep"
        # FASTA files to 'fna.gz' from their usual 'tar.gz'.
        obs_query_urls = _assemble_queries('220.0', 'SpeciesReps', 'Both')
        print('obs queries: ', obs_query_urls)

        exp_query_urls = [('Archaea',
                           'https://data.gtdb.ecogenomic.org/releases/'
                           'release220/220.0/genomic_files_reps/'
                           'ar53_ssu_reps_r220.fna.gz'),
                          ('Bacteria',
                           'https://data.gtdb.ecogenomic.org/releases/'
                           'release220/220.0/genomic_files_reps/'
                           'bac120_ssu_reps_r220.fna.gz')]
        print('exp queries: ', exp_query_urls)
        self.assertEqual(obs_query_urls, exp_query_urls)

        # test that these URLs work
        for _, u in obs_query_urls:
            try:
                urlopen(u)
            except HTTPError:
                raise ValueError('Failed to open URL: ' + u)

    def test_assemble_species_rep_queries_archaea(self):
        obs_query_urls = _assemble_queries('202.0', 'SpeciesReps', 'Archaea')
        print('obs queries: ', obs_query_urls)

        exp_query_urls = [('Archaea',
                           'https://data.gtdb.ecogenomic.org/releases/'
                           'release202/202.0/genomic_files_reps/'
                           'ar122_ssu_reps_r202.tar.gz')]
        print('exp queries: ', exp_query_urls)
        self.assertEqual(obs_query_urls, exp_query_urls)

        # test that these URLs work
        for _, u in obs_query_urls:
            try:
                urlopen(u)
            except HTTPError:
                raise ValueError('Failed to open URL: ' + u)

    def test_assemble_queries_all(self):
        obs_query_urls = _assemble_queries('207.0', 'All')
        print('obs queries: ', obs_query_urls)

        exp_query_urls = [('All',
                           'https://data.gtdb.ecogenomic.org/releases/'
                           'release207/207.0/genomic_files_all/'
                           'ssu_all_r207.tar.gz')]
        print('exp queries: ', exp_query_urls)
        self.assertEqual(obs_query_urls, exp_query_urls)

        # test that these URLs work
        for _, u in obs_query_urls:
            try:
                urlopen(u)
            except HTTPError:
                raise ValueError('Failed to open URL: ' + u)

    def test_get_gtdb(self):
        def _makey_fakey_both(faking_ignore_this):
            return self.gtdb_tax, self.gtdb_seqs

        def _makey_fakey_arch(faking_ignore_this):
            return self.gtdb_arch_tax, self.gtdb_arch_seqs

        def _makey_fakey_bact(faking_ignore_this):
            return self.gtdb_bact_tax, self.gtdb_bact_seqs

        # default (both domains, version 214)
        with patch('rescript.get_gtdb._retrieve_data_from_gtdb',
                   new=_makey_fakey_both):
            res = rescript.actions.get_gtdb_data(
                version='214.1', db_type='SpeciesReps', domain='Both')
            self.assertEqual(str(res[0].type), 'FeatureData[Taxonomy]')
            self.assertEqual(str(res[1].type), 'FeatureData[Sequence]')

        # just grab archaea domain, and version 202
        with patch('rescript.get_gtdb._retrieve_data_from_gtdb',
                   new=_makey_fakey_arch):
            resa = rescript.actions.get_gtdb_data(
                version='202.0', domain='Archaea')
            self.assertEqual(str(resa[0].type), 'FeatureData[Taxonomy]')
            self.assertEqual(str(resa[1].type), 'FeatureData[Sequence]')

        # just grab bacteria domain, and version 207
        with patch('rescript.get_gtdb._retrieve_data_from_gtdb',
                   new=_makey_fakey_bact):
            resb = rescript.actions.get_gtdb_data(
                version='207.0', domain='Bacteria')
            self.assertEqual(str(resb[0].type), 'FeatureData[Taxonomy]')
            self.assertEqual(str(resb[1].type), 'FeatureData[Sequence]')

        # Non-species rep version 214.1)
        with patch('rescript.get_gtdb._retrieve_data_from_gtdb',
                   new=_makey_fakey_both):
            resc = rescript.actions.get_gtdb_data(
                version='214.1', db_type='All')
            self.assertEqual(str(resc[0].type), 'FeatureData[Taxonomy]')
            self.assertEqual(str(resc[1].type), 'FeatureData[Sequence]')

    def test_parse_gtdb_taxonomy(self):
        tax_in = ('d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Oenococcus;s__Oenococcus oeni '
                  '[locus_tag=NZ_AQVA01000009.1] [location=77871..79431] '
                  '[ssu_len=1561] [contig_len=79790]')
        exp = ('d__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;'
               'f__Lactobacillaceae;g__Oenococcus;s__Oenococcus oeni')
        self.assertEqual(parse_gtdb_taxonomy(tax_in), exp)