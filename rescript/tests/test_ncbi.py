# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import warnings

import qiime2
from qiime2 import Metadata
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from pandas import DataFrame
from q2_types.feature_data import DNAIterator

import_data = qiime2.Artifact.import_data


class TestNCBI(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.get_ncbi_data = rescript.methods.get_ncbi_data

        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('ncbi-seqs.fasta'))
        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('ncbi-taxa.tsv'))
        self.non_standard_taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('ns-ncbi-taxa.tsv'))

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

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            acc_seq, acc_tax = self.get_ncbi_data(accession_ids=md)
            self.assertGreaterEqual(len(w), 2)
            num_user_warnings = sum(_w.category is UserWarning for _w in w)
            self.assertGreaterEqual(num_user_warnings, 2)

        acc_seq = {s.metadata['id']: str(s) for s in acc_seq.view(DNAIterator)}
        seqs = {s.metadata['id']: str(s) for s in self.seqs.view(DNAIterator)}
        self.assertEqual(acc_seq, seqs)

        acc_tax = acc_tax.view(DataFrame).to_dict()
        taxa = self.taxa.view(DataFrame).to_dict()
        self.assertEqual(acc_tax, taxa)

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
