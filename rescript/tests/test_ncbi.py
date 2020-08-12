# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pickle
from unittest import mock

from requests.exceptions import (
    HTTPError, ChunkedEncodingError, ConnectionError, ReadTimeout)
from xml.parsers.expat import ExpatError
import qiime2
from qiime2 import Metadata
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from pandas import DataFrame
from q2_types.feature_data import DNAIterator
import rescript.ncbi as ncbi

import_data = qiime2.Artifact.import_data


class TestNCBI(TestPluginBase):
    package = 'rescript.tests'

    def mock_get(self, *args, **kwargs):
        if self.ncbi_exception:
            exception = self.ncbi_exception
            self.ncbi_exception = None
            raise exception
        for response in self.responses:
            if response[0] == ['get', args, kwargs]:
                return response[1]

    def mock_post(self, *args, **kwargs):
        if self.ncbi_exception:
            exception = self.ncbi_exception
            self.ncbi_exception = None
            raise exception
        for response in self.responses:
            if response[0][:2] == ['post', args]:
                req_ids = set(kwargs['data']['id'].split(','))
                res_ids = set(response[0][2]['data']['id'].split(','))
                if req_ids == res_ids:
                    return response[1]

    def setUp(self):
        super().setUp()

        self.get_ncbi_data = rescript.methods.get_ncbi_data

        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('ncbi-seqs.fasta'))
        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('ncbi-taxa.tsv'))
        self.non_standard_taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('ns-ncbi-taxa.tsv'))
        ncbi_responses = self.get_data_path('ncbi-responses.pkl')
        with open(ncbi_responses, 'rb') as fh:
            self.responses = pickle.load(fh)

        self.actual_requests = ncbi.requests
        ncbi.requests = mock.Mock(**{'get.side_effect': self.mock_get,
                                     'post.side_effect': self.mock_post})
        self.ncbi_exception = None

    def tearDown(self):
        super().tearDown()

        ncbi.requests = self.actual_requests

    def test_get_ncbi_data_accession_ids_no_rank_propagation(self):
        df = DataFrame(index=['M59083.2', 'AJ234039.1'])
        df.index.name = 'id'
        md = Metadata(df)

        acc_seq, acc_tax = self.get_ncbi_data(accession_ids=md,
                                              rank_propagation=False)
        acc_seq = {s.metadata['id']: str(s) for s in acc_seq.view(DNAIterator)}
        seqs = {s.metadata['id']: str(s) for s in self.seqs.view(DNAIterator)}
        self.assertEqual(acc_seq, seqs)

        acc_tax = acc_tax.view(DataFrame).to_dict()
        taxa = self.taxa.view(DataFrame).to_dict()
        self.assertEqual(acc_tax, taxa)

    def test_get_ncbi_data_bad_accession(self):
        df = DataFrame(index=['M59083.2', 'not_an_accession', 'AJ234039.1'])
        df.index.name = 'id'
        md = Metadata(df)
        with self.assertLogs(level='WARNING') as log:
            with self.assertRaises(RuntimeError):
                acc_seq, acc_tax = self.get_ncbi_data(accession_ids=md)
        self.assertEqual(log.output,
                         ['WARNING:root:Some IDs have invalid value and were '
                          'omitted. Maximum ID value 18446744073709551615'])

    def test_get_ncbi_data_rank_propagation_nonstandard_ranks(self):
        df = DataFrame(index=['M59083.2', 'AJ234039.1'])
        df.index.name = 'id'
        md = Metadata(df)

        acc_seq, acc_tax = self.get_ncbi_data(
                accession_ids=md,
                ranks=['subkingdom', 'subphylum', 'subclass', 'suborder',
                       'subfamily', 'subgenus', 'subspecies'])
        acc_seq = {s.metadata['id']: str(s) for s in acc_seq.view(DNAIterator)}
        seqs = {s.metadata['id']: str(s) for s in self.seqs.view(DNAIterator)}
        self.assertEqual(acc_seq, seqs)

        acc_tax = acc_tax.view(DataFrame).to_dict()
        taxa = self.non_standard_taxa.view(DataFrame).to_dict()
        self.assertEqual(acc_tax, taxa)

    def test_get_ncbi_data_query(self):
        que_seq, que_tax = self.get_ncbi_data(query='M59083.2 OR AJ234039.1')

        que_seq = {s.metadata['id']: str(s) for s in que_seq.view(DNAIterator)}
        seqs = {s.metadata['id']: str(s) for s in self.seqs.view(DNAIterator)}
        self.assertEqual(que_seq, seqs)

        que_tax = que_tax.view(DataFrame).to_dict()
        taxa = self.taxa.view(DataFrame).to_dict()
        self.assertEqual(que_tax, taxa)

    def test_get_ncbi_data_query_mushroom_one(self):
        seq, tax = self.get_ncbi_data(query='MT345279.1')

        tax = tax.view(DataFrame)
        self.assertEqual(
            tax['Taxon']['MT345279.1'],
            'k__Fungi; p__Basidiomycota; c__Agaricomycetes; o__Boletales; '
            'f__Boletaceae; g__Boletus; s__edulis'
        )

    def test_get_ncbi_data_query_mushroom_two(self):
        seq, tax = self.get_ncbi_data(
            query='MT345279.1',
            ranks=['domain', 'phylum', 'subphylum', 'superfamily', 'family',
                   'subfamily', 'genus', 'species', 'subspecies']
        )

        tax = tax.view(DataFrame)
        self.assertEqual(
            tax['Taxon']['MT345279.1'],
            'd__Plants and Fungi; p__Basidiomycota; ps__Agaricomycotina; '
            'sf__Boletineae; f__Boletaceae; fs__Boletoideae; g__Boletus; '
            's__edulis; ssb__edulis'
        )

    def test_get_ncbi_data_query_mushroom_three(self):
        seq, tax = self.get_ncbi_data(
            query='MT345279.1',
            ranks=['domain', 'phylum', 'subphylum', 'superfamily', 'family',
                   'subfamily', 'genus', 'species', 'subspecies'],
            rank_propagation=False
        )

        tax = tax.view(DataFrame)
        self.assertEqual(
            tax['Taxon']['MT345279.1'],
            'd__Plants and Fungi; p__Basidiomycota; ps__Agaricomycotina; '
            'sf__; f__Boletaceae; fs__Boletoideae; g__Boletus; '
            's__edulis; ssb__'
        )

    def test_ncbi_fails(self):
        exceptions = [ChunkedEncodingError(), ConnectionError(), ReadTimeout(),
                      ExpatError(), RuntimeError('bad record')]
        http_exception = HTTPError()
        http_exception.response = self.responses[0][1]
        http_exception.response.status_code = 429
        exceptions.append(http_exception)

        for exception in exceptions:
            self.ncbi_exception = exception
            with self.assertLogs(level='DEBUG') as log:
                seq, tax = self.get_ncbi_data(query='MT345279.1')
                tax = tax.view(DataFrame)
                self.assertEqual(
                    tax['Taxon']['MT345279.1'],
                    'k__Fungi; p__Basidiomycota; c__Agaricomycetes; '
                    'o__Boletales; f__Boletaceae; g__Boletus; s__edulis'
                )
            self.assertTrue('Retrying' in log.output[0])
