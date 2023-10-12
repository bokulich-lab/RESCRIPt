# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from rescript.extract_seq_segments import extract_seq_segments

import_data = qiime2.Artifact.import_data


class TestExtractSeqSegments(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        input_seq_fp = self.get_data_path('full-length-query.fasta')
        self.seqs = DNAFASTAFormat(input_seq_fp, mode='r')

        input_segref_fp = self.get_data_path('ref-segs.fasta')
        self.segref = DNAFASTAFormat(input_segref_fp, mode='r')

    def test_extract_seq_segments(self):
        ext_segs, unmat_seqs,  = extract_seq_segments(
                        input_sequences=self.seqs,
                        reference_segment_sequences=self.segref,
                        perc_identity=0.9)

        obs_ext_segs = {seq.metadata['id']: str(seq)
                        for seq in ext_segs.view(DNAIterator)}
        exp_ext_segs = {
                   'fulllen01': ('ATCTTTATGTTTAGAAAAACAAGGGTTTTAATTGAAA'
                                 'AACTAGAAGAAAAAGG'),
                   'fulllen02': ('ATCCTGTTTTTCGAAAACAAACAAAGGTTTGTAAAGA'
                                 'CAGAAGCAGAATACAAAAACAAAAG')}
        self.assertEqual(obs_ext_segs, exp_ext_segs)

        obs_unmat_seqs = {seq.metadata['id']: str(seq)
                          for seq in unmat_seqs.view(DNAIterator)}
        exp_unmat_seqs = {
                   'fulllen03': (
                    'ACCATTTCCCATGCATCATCCTGGTGTAGTATTTGTGTCTATGTCAATTAAAGGGA'
                    'CTAAAAAGTAGTAAAAAATTTAACTAGTAAGCGTAAGGATAATGATATGGACTGTG'
                    'AATGATTCAATAATAGAGATTATTTGCCCATATATGTTTATTTGTAGAGGTATTGT'
                    'ATCTATCGAAATTGGGATAAGATCCAAGGATTTTAGTTCGGATACATTTGTGTTGT'
                    'GAAAGAGTGGAGTGAATGAGAAAGATAGTGAATTTTGTTTGAACTGAATCGCTGAT'
                    'GAAAAAAAAAGGAAGATAAATATTTAGGAAGTAAAATTACCTTTTTATTGGGGATA'
                    'GAGGGACTTGAACCCTCACGATTTCTAAAGTCGACGGATTTTCCTCTTACTATAAA'
                    'TTTCATTGTTGTCGGTATTGACATGTAGAATGGGACTCTCTCTTTATTCTCGTCCG'
                    'ATTAATCAGTTTTTCAAAAGATCTATCAAACTCTGGAATGAATGATTTAATAATTG'
                    'AATTGAATATTCAATTCTTCTTTCAACTTCCATTGGACTGGATTCACAATAACTCT'
                    'AATTCTGAATTTTTCATATTTTAATTATCCATATAATATAATGGATTCGAGTCATG'
                    'ATTAATCGTTTGATTTGATATGTCAGTATGTATACGTATTTATTAGGTATATAGGA'
                    'ACATCTTTTCTGTAATTTTGATAGACGGATTCCAACTACCAACGAAACGTAGTCAA'
                    'CTTCATTCGTTAGAACAGCTTCCATTGAGTCTCTGCACCTATCCGAATTCTAGTTT'
                    'GATAAAACTAAGGTTTATCAAACTAAGGATTTGGCTCAGGATTGCCCAT'),
                   'fulllen04': (
                    'ACCATTTCCCATGCATCATCCTGGTGTAGTGTTTGTGTCTATGTCAATTAAAGGGA'
                    'CTAAAAAGTAGTAAAAAATTTAACTAGTAAGCGTAAGGATAATGATATGGACTGTG'
                    'AATGATTCAATAATAGAGATTATTTTCCCATATATGCTTATTTGTAGAGGTATTGT'
                    'ATCTATCGAAATTGGGATAAGATCCAAGGATTTTAGTTCGGATACATTTGTGTTGT'
                    'GAAAGAGTGGAGTGAATGAGAAAGATAGTGAATTTTGTTTGAACTGAATCGCTGAT'
                    'GAAAAAAAAAGGAGGATAAATATTTAGGAAGTAAAATTACCTTTTTATTGGGGATA'
                    'GAGGGACTTGAACCCTCACGATTTCTAAAGTCGACGGATTTTCCTCTTACTATAAA'
                    'TTTCATTGTTGTCGGTATTGACATGTAGAATGGGACTCTCTCTTTATTCTCGTCCG'
                    'ATTAATCAGTTTTTCAAAAGATCTATCAAACTCTGGAATGAATGATTTAATAATTG'
                    'AATTGAATATTCAATTCTTCTTTCAACTTCCATTGGACTGGATTCACAATAACTCT'
                    'AATTCTGAATTTTTCATATTATAATTATCCATATAATATAATGGATTCGAGTCATG'
                    'ATTAATCGTTTGATTTGATATGTCAGTATGTATACGTATGTATTAGGTATATAGGA'
                    'ACATCTTTTCTGTAATTTTGATAGACGGATTCCAACTACCAACGAAACGTAGTCAA'
                    'CTTCATTCGTTAGAACAGCTTCCATTGAGTCTCTGCACCTATCCTTATTCTAGTTT'
                    'GATAAAACTAAGGTTTATCAAACTAAGGATTTGGCTCAGGATTGCCCAT')}
        self.assertEqual(obs_unmat_seqs, exp_unmat_seqs)
