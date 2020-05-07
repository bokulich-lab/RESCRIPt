#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from rescript.cleanseq import clean_sequences
from q2_types.feature_data import DNAFASTAFormat, DNAIterator


import_data = qiime2.Artifact.import_data


class TestCleanseq(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        input_fp = self.get_data_path('cleanseq-test-1.fasta')
        self.seqs1 = DNAFASTAFormat(input_fp, mode='r')

    def test_cleanseq_default_params(self):
        # Test default params: max_degenerates = 4, homopolymer_length = 8
        obs = clean_sequences(self.seqs1)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Ambig2', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)

    def test_cleanseq_ambig_two_hpoly_seven(self):
        # Test params: max_degenerates = 2, homopolymer_length = 9
        # Keep seq with 2 ambigs, and 7 hpoly
        obs = clean_sequences(self.seqs1, max_degenerates=2,
                              homopolymer_length=7)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Ambig2', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)

    def test_cleanseq_ambig_zero(self):
        # Test params: mmax_degenerates = 0
        # Only seqs w/ no ambig bases returned
        obs = clean_sequences(self.seqs1, max_degenerates=0,
                              homopolymer_length=9)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Hpoly8', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)
