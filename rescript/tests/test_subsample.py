# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (AlignedDNAFASTAFormat, DNAIterator)

from rescript.subsample import subsample_fasta

import_data = qiime2.Artifact.import_data


class TestSubsampleSeq(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        input_fp = self.get_data_path('trim-test-alignment.fasta')
        self.alignedseqs = AlignedDNAFASTAFormat(input_fp, mode='r')

    def test_subsample_seqs_count(self):
        obs = subsample_fasta(
            self.alignedseqs, subsample_size=2, random_seed=42)

        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        exp_seqs = {'s1': '-----TAGGGAATCTTCCACAATGGGTGCAAACCTGATGGAGCAATG'
                          'CCGCGTGAGTGAAGANAGGTCTTCGGATCGTAAAGCTCTGTTGTTAG'
                          'AGAAGAACACGTGCTAGG--------',
                    's5': '---AATGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATG'
                          'CCGCGTGTGTGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGG'
                          'GGAGGAAAAGCTTGTGGTTAA-----', }
        self.assertEqual(obs_seqs, exp_seqs)

    def test_subsample_seqs_fraction(self):
        obs = subsample_fasta(
            self.alignedseqs, subsample_size=0.8, random_seed=42)

        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        exp_seqs = {'s1': '-----TAGGGAATCTTCCACAATGGGTGCAAACCTGATGGAGCAATG'
                          'CCGCGTGAGTGAAGANAGGTCTTCGGATCGTAAAGCTCTGTTGTTAG'
                          'AGAAGAACACGTGCTAGG--------',
                    's2': 'AATTTTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACG'
                          'CCGCGTGGGGGATGA-CGGCCTTCGGGTTGTAAACTCCTTTCGCCAG'
                          'GGACGAAGCGTTTTG-----------',
                    's3': '-----TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATG'
                          'CCGCGTGTGTGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGG'
                          'GGAGGAAAAGCTTATGGTTAAAAAAA',
                    's5': '---AATGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATG'
                          'CCGCGTGTGTGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGG'
                          'GGAGGAAAAGCTTGTGGTTAA-----', }
        self.assertEqual(obs_seqs, exp_seqs)

    def test_subsample_seqs_count_exceeded(self):
        with self.assertRaisesRegex(
                ValueError, '.*sample size 6 is bigger.*sequences 5'):
            subsample_fasta(
                self.alignedseqs, subsample_size=6, random_seed=42)
