# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (AlignedDNAFASTAFormat, AlignedDNAIterator,
                                   DNAIterator)

from rescript.degap import degap_sequences


import_data = qiime2.Artifact.import_data


class TestDegapSeq(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        input_fp = self.get_data_path('degap-test-alignment.fasta')
        self.alignedseqs = AlignedDNAFASTAFormat(
                                input_fp, mode='r').view(
                                AlignedDNAIterator)

    def test_degap_sequences(self):
        #  remove all '-' and '.' chars.
        #  seq 's8' should not be returned as it is all gaps.
        obs = degap_sequences(self.alignedseqs)

        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        exp_seqs = {'s1': ('TAGGGAATCTTCCACAATGGGTGCAAACCTGATGGAGCAATGCCGCGT'
                           'GAGTGAAGANAGGTCTTCGGATCGTAAAGCTCTGTTGTTAGAGAAGAA'
                           'CACGTGCTAGG'),
                    's2': ('AATTTTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGC'
                           'CGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGA'
                           'CGAAGCGTTTTG'),
                    's3': ('TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGT'
                           'GTGTGAAGAAGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAA'
                           'AGCTTATGGTTAAAAAAA'),
                    's4': ('TGGGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATGCCGCGT'
                           'GCGGGAAGANAGGCCTTCGGGTTGTAAACCGCTTTTGTTCGGGAAGAA'
                           'ATCCTCTGGGCTAAAAAAA'),
                    's5': ('AATGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGC'
                           'GTGTGTGAAGAAGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGA'
                           'AAAGCTTGTGGTTAA'),
                    's6': ('G'),
                    's7': ('GG')}
        self.assertEqual(obs_seqs, exp_seqs)
