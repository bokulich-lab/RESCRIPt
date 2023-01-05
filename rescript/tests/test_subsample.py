# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (AlignedDNAFASTAFormat, DNAIterator,
                                   AlignedDNAIterator, DNAFASTAFormat)

from rescript.subsample import subsample_fasta

import_data = qiime2.Artifact.import_data


class TestSubsampleSeq(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        aligned_input_fp = self.get_data_path('trim-test-alignment.fasta')
        unaligned_input_fp = self.get_data_path('trim-test-sequences.fasta')
        self.alignedseqs = AlignedDNAFASTAFormat(
            aligned_input_fp, mode='r').view(AlignedDNAIterator)
        self.seqs = DNAFASTAFormat(
            unaligned_input_fp, mode='r').view(DNAIterator)

    def test_subsample_seqs(self):
        obs = subsample_fasta(
            self.alignedseqs, subsample_size=0.4, random_seed=42)

        obs_seqs = {seq.metadata['id']: str(seq) for seq in obs}
        exp_seqs = {'s2': 'AATTTTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACG'
                          'CCGCGTGGGGGATGA-CGGCCTTCGGGTTGTAAACTCCTTTCGCCAG'
                          'GGACGAAGCGTTTTG-----------',
                    's3': '-----TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATG'
                          'CCGCGTGTGTGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGG'
                          'GGAGGAAAAGCTTATGGTTAAAAAAA',
                    's4': '-----TGGGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATG'
                          'CCGCGTGCGGGAAGANAGGCCTTCGGGTTGTAAACCGCTTTTGTTCG'
                          'GGAAGAAATCCTCTGGGCTAAAAAAA', }
        self.assertEqual(obs_seqs, exp_seqs)

    def test_subsample_seqs_count_unaligned(self):
        obs = subsample_fasta(
            self.seqs, subsample_size=0.2, random_seed=42)

        obs_seqs = {seq.metadata['id']: str(seq) for seq in obs}
        exp_seqs = {'s2': 'AATTTTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACG'
                          'CCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGG'
                          'GACGAAGCGTTTTG', }
        self.assertEqual(obs_seqs, exp_seqs)
