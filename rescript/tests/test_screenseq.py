# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat, DNAIterator, RNAFASTAFormat

from rescript.screenseq import cull_seqs


import_data = qiime2.Artifact.import_data


class TestCleanseq(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        input_fp = self.get_data_path('cleanseq-test-1.fasta')
        self.seqs1 = DNAFASTAFormat(input_fp, mode='r').view(DNAIterator)

    def test_cull_seqs_default_params(self):
        # Test default params: num_degenerates = 5, homopolymer_length = 8
        obs = cull_seqs(self.seqs1)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Ambig2', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)

    def test_cull_seqs_rna_default_params(self):
        # Test default params work with RNA seqs as input
        rna_path = self.get_data_path('cleanseq-test-1-rna.fasta')
        rna_seqs = RNAFASTAFormat(rna_path, mode='r').view(DNAIterator)
        obs = cull_seqs(rna_seqs)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Ambig2', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)

    def test_cull_seqs_degen_seven_hpoly_seven(self):
        # Test params: num_degenerates = 7, homopolymer_length = 7
        # Keep seq with 2 ambigs, and 7 hpoly
        obs = cull_seqs(self.seqs1, num_degenerates=7, homopolymer_length=7)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Ambig2', 'Ambig6', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)

    def test_cull_seqs_degen_seven_hpoly_nine(self):
        # Test params: num_degenerates = 7, homopolymer_length = 9
        # All seqs should pass.
        obs = cull_seqs(self.seqs1, num_degenerates=7, homopolymer_length=9)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Ambig6', 'Ambig2', 'Hpoly8', 'Hpoly8Ambig1', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)

    def test_cull_seqs_degen_one(self):
        # Test params: mnum_degenerates = 2
        # Only seqs without degen bases returned
        obs = cull_seqs(self.seqs1, num_degenerates=1, homopolymer_length=9)
        obs_ids = {seq.metadata['id'] for seq in obs.view(DNAIterator)}
        exp_ids = {'Hpoly8', 'cleanseq'}
        self.assertEqual(obs_ids, exp_ids)
