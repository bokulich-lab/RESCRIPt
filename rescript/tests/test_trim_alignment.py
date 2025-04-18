# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
import skbio
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (
    AlignedDNAFASTAFormat,
    DNAIterator,)

from qiime2.plugins.alignment.methods import mafft_add

from rescript.trim_alignment import (
    _prepare_positions, _process_primers, _locate_primer_positions,
    _trim_all_sequences, _trim_alignment)

import_data = qiime2.Artifact.import_data


class FakeCtx:
    def __init__(self, alignment_dict):
        self.alignments = alignment_dict

    def aln1(self, *args, **kwargs):
        return [self.alignments[1]]

    def aln2(self, *args, **kwargs):
        return [self.alignments[2]]

    def aln3(self, *args, **kwargs):
        return [self.alignments[3]]

    def get_action(self, which_alignment):
        return getattr(self, f"aln{which_alignment}")


class TestExtractAlignmentRegion(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        aligned_silva_seqs_fp = self.get_data_path(
            'small-silva-full-len-alignment.fasta')
        self.aligned_silva_seqs = AlignedDNAFASTAFormat(
             aligned_silva_seqs_fp, mode='r')

        self.v4_primers_dict = {
            "forward": "GTGYCAGCMGCCGCGGTAA",
            "reverse": "GGACTACNVGGGTWTCTAAT"
        }

        silva_alignment_fp = self.get_data_path(
            'small-silva-full-len-alignment.fasta')
        self.silva_alignment = AlignedDNAFASTAFormat(
            silva_alignment_fp, mode='r')
        self.aligned_silva_seqs_art = qiime2.Artifact.import_data(
                                       'FeatureData[AlignedSequence]',
                                       self.silva_alignment)

        silva_alignment_v4_trim_keeplen_nokeep_primer_fp = self.get_data_path(
            'small-silva-v4-trim-keeplength-nokeepprimers.fasta')
        self.silva_v4_trim_keeplen_nokeep_primer = AlignedDNAFASTAFormat(
            silva_alignment_v4_trim_keeplen_nokeep_primer_fp, mode='r')

        silva_alignment_v4_trim_keeplen_wo_primers_fp = self.get_data_path(
            'small-silva-v4-trim-keeplength.fasta')
        self.silva_v4_trim_keeplen = AlignedDNAFASTAFormat(
            silva_alignment_v4_trim_keeplen_wo_primers_fp, mode='r')

        silva_alignment_v4_trim_no_keeplen_wo_primers_fp = self.get_data_path(
            'small-silva-v4-trim-no-keeplength.fasta')
        self.silva_v4_trim_no_keeplen = AlignedDNAFASTAFormat(
            silva_alignment_v4_trim_no_keeplen_wo_primers_fp, mode='r')

        aligned_seqs_fp = self.get_data_path('trim-test-alignment.fasta')
        self.aligned_seqs = qiime2.Artifact.import_data(
            'FeatureData[AlignedSequence]', aligned_seqs_fp)
        self.aligned_seqs_fasta = AlignedDNAFASTAFormat(
            aligned_seqs_fp, mode='r')

        self.primers_dict = {
            "forward": "GGGAATCTTCCACAATGG",
            "reverse": "GTGTTCTTCTCTAACAACAG"
        }

        aligned_with_primers_fp = self.get_data_path(
            'trim-test-alignment-with-primers.fasta')
        self.aligned_with_primers = qiime2.Artifact.import_data(
            'FeatureData[AlignedSequence]', aligned_with_primers_fp)
        self.aligned_with_primers_fasta = AlignedDNAFASTAFormat(
            aligned_with_primers_fp, mode='r')

        self.aligned_mess_fasta = AlignedDNAFASTAFormat(
            self.get_data_path(
                'trim-test-alignment-with-primers-mess.fasta'), mode='r')
        self.aligned_with_fwd_fasta = AlignedDNAFASTAFormat(
            self.get_data_path('trim-test-alignment-fwd.fasta'), mode='r')
        self.aligned_with_rev_fasta = AlignedDNAFASTAFormat(
            self.get_data_path('trim-test-alignment-rev.fasta'), mode='r')
        self.trimmed_fasta = AlignedDNAFASTAFormat(
            self.get_data_path('trim-test-sequences-trimmed.fasta'), mode='r')

        self.fake_ctx = FakeCtx({1: self.aligned_with_primers_fasta,
                                 2: self.aligned_with_fwd_fasta,
                                 3: self.aligned_with_rev_fasta})

        self.exp_seqs_both_primers = {
            's1': ('GGGAATCTTCCACAATGGGTGCAAACCTGATGGAGCAATGCCGCGTGAG'
                   'TGAAGANAGGTCTTCGGATCGTAAAGCTCTGTTGTTAGAGAAGAACAC'),
            's2': ('GGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGG'
                   'GGATGA-CGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCG'),
            's3': ('GGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTG'
                   'TGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAG'),
            's4': ('GGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATGCCGCGTGCG'
                   'GGAAGANAGGCCTTCGGGTTGTAAACCGCTTTTGTTCGGGAAGAAATC'),
            's5': ('GGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTG'
                   'TGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAG'), }
        self.exp_seqs_both_no_primers = {
            's1': ('GGGAATCTTCCACAATGGGTGCAAACCTGATGGAGCAATGCCGCGTGAG'
                   'TGAAGANAGGTCTTCGGATCGTAAAGCTCTGTTGTTAGAGAAGAACAC'),
            's2': ('GGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGG'
                   'GGATGA-CGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCG'),
            's3': ('GGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTG'
                   'TGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAG'),
            's4': ('GGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATGCCGCGTGCG'
                   'GGAAGANAGGCCTTCGGGTTGTAAACCGCTTTTGTTCGGGAAGAAATC'),
            's5': ('GGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTG'
                   'TGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAG'), }
        self.exp_seqs_only_fwd = {
            's1': ('GGGAATCTTCCACAATGGGTGCAAACCTGATGGAGCAATGCCGCGTGAG'
                   'TGAAGANAGGTCTTCGGATCGTAAAGCTCTGTTGTTAGAGAAGAACACG'
                   'TGCTAGG--------'),
            's2': ('GGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGG'
                   'GGATGA-CGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGT'
                   'TTTG-----------'),
            's3': ('GGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTG'
                   'TGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAGC'
                   'TTATGGTTAAAAAAA'),
            's4': ('GGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATGCCGCGTGCG'
                   'GGAAGANAGGCCTTCGGGTTGTAAACCGCTTTTGTTCGGGAAGAAATCC'
                   'TCTGGGCTAAAAAAA'),
            's5': ('GGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTG'
                   'TGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAGC'
                   'TTGTGGTTAA-----'), }
        self.exp_seqs_only_rev = {
            's1': ('-----TAGGGAATCTTCCACAATGGGTGCAAACCTGATGGAGCAATGC'
                   'CGCGTGAGTGAAGANAGGTCTTCGGATCGTAAAGCTCTGTTGTTAGAG'
                   'AAGAACAC'),
            's2': ('AATTTTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGC'
                   'CGCGTGGGGGATGA-CGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGG'
                   'ACGAAGCG'),
            's3': ('-----TGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGC'
                   'CGCGTGTGTGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGG'
                   'AGGAAAAG'),
            's4': ('-----TGGGGAATTTTGGACAATGGGGGCAACCCTGATCCAGCCATGC'
                   'CGCGTGCGGGAAGANAGGCCTTCGGGTTGTAAACCGCTTTTGTTCGGG'
                   'AAGAAATC'),
            's5': ('---AATGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGC'
                   'CGCGTGTGTGAAGA-AGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGG'
                   'AGGAAAAG'), }

    def test_prepare_positions(self):
        pos_start, pos_end, aln_len = 4, 16, 20
        obs_pos = _prepare_positions(pos_start, pos_end, aln_len)
        exp_pos = {"start": 3, "end": 16}
        self.assertDictEqual(obs_pos, exp_pos)

    def test_prepare_positions_no_start(self):
        pos_start, pos_end, aln_len = None, 16, 20
        obs_pos = _prepare_positions(pos_start, pos_end, aln_len)
        exp_pos = {"start": None, "end": 16}
        self.assertDictEqual(obs_pos, exp_pos)

    def test_prepare_positions_no_end(self):
        pos_start, pos_end, aln_len = 4, None, 20
        obs_pos = _prepare_positions(pos_start, pos_end, aln_len)
        exp_pos = {"start": 3, "end": None}
        self.assertDictEqual(obs_pos, exp_pos)

    def test_prepare_positions_end_too_high(self):
        pos_start, pos_end, aln_len = 4, 24, 20
        obs_pos = _prepare_positions(pos_start, pos_end, aln_len)
        exp_pos = {"start": 3, "end": 20}
        self.assertDictEqual(obs_pos, exp_pos)

    def test_prepare_positions_both_none(self):
        with self.assertRaisesRegex(
                ValueError, 'Neither primers nor.*positions were provided'):
            _prepare_positions(None, None, 20)

    def test_prepare_positions_start_too_high1(self):
        with self.assertRaisesRegex(
                ValueError, 'Start.*should be smaller.*end position'):
            _prepare_positions(18, 10, 20)

    def test_prepare_positions_start_too_high2(self):
        with self.assertRaisesRegex(
                ValueError, 'Start.*should be smaller.*alignment'):
            _prepare_positions(22, None, 20)

    def test_process_primers(self):
        obs_primers = _process_primers(
            self.primers_dict["forward"],
            self.primers_dict["reverse"])
        obs_primers = {seq.metadata['id']: str(seq)
                       for seq in obs_primers.view(DNAIterator)}
        exp_primers = {"forward": "GGGAATCTTCCACAATGG",
                       "reverse": "CTGTTGTTAGAGAAGAACAC"}
        self.assertDictEqual(obs_primers, exp_primers)

    def test_process_primers_only_fwd(self):
        obs_primers = _process_primers(self.primers_dict["forward"], None)
        obs_primers = {seq.metadata['id']: str(seq)
                       for seq in obs_primers.view(DNAIterator)}
        exp_primers = {"forward": "GGGAATCTTCCACAATGG"}
        self.assertDictEqual(obs_primers, exp_primers)

    def test_process_primers_only_rev(self):
        obs_primers = _process_primers(None, self.primers_dict["reverse"])
        obs_primers = {seq.metadata['id']: str(seq)
                       for seq in obs_primers.view(DNAIterator)}
        exp_primers = {"reverse": "CTGTTGTTAGAGAAGAACAC"}
        self.assertDictEqual(obs_primers, exp_primers)

    def test_locate_positions(self):
        obs_pos = _locate_primer_positions(self.aligned_with_primers_fasta,
                                           keep_primer_location=True)
        exp_pos = {"start": 7, "end": 104}
        self.assertDictEqual(obs_pos, exp_pos)

    def test_locate_positions_only_fwd(self):
        obs_pos = _locate_primer_positions(self.aligned_with_fwd_fasta,
                                           keep_primer_location=True)
        exp_pos = {"start": 7, "end": None}
        self.assertDictEqual(obs_pos, exp_pos)

    def test_locate_positions_only_rev(self):
        obs_pos = _locate_primer_positions(self.aligned_with_rev_fasta,
                                           keep_primer_location=True)
        exp_pos = {"start": None, "end": 104}
        self.assertDictEqual(obs_pos, exp_pos)

    def test_locate_positions_no_primers(self):
        obs_pos = _locate_primer_positions(self.aligned_seqs_fasta,
                                           keep_primer_location=True)
        exp_pos = {"start": None, "end": None}
        self.assertDictEqual(obs_pos, exp_pos)

    # test locating positions in an alignment where wrong primers were used;
    # here, two random primers were aligned and the reverse not only is
    # aligned with gaps but also its start is upstream to the start of
    # the forward primer - clearly something went wrong with the alignment;
    # user should be warned
    def test_locate_positions_strange_alignment(self):
        with self.assertRaisesRegex(
                ValueError, 'Reverse primer overlaps'):
            _locate_primer_positions(self.aligned_mess_fasta,
                                     keep_primer_location=True)

    # test trimming with both, start and end, positions given
    def test_trim_all_sequences(self):
        obs = _trim_all_sequences(
            self.aligned_seqs,
            {"start": 7, "end": 104})
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_both_primers)

    # test trimming with only end position given
    def test_trim_all_sequences_no_fwd(self):
        obs = _trim_all_sequences(
            self.aligned_seqs,
            {"start": None, "end": 104})
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_only_rev)

    # test trimming with only start position given
    def test_trim_all_sequences_no_rev(self):
        obs = _trim_all_sequences(
            self.aligned_seqs,
            {"start": 7, "end": None})
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_only_fwd)

    # test trimming when both primers are given and keeplength = False
    # tests against expected alignment length
    def test_trim_alignment_keeplen_false(self):
        obs_v4_nokeep_aln = _trim_alignment(
            mafft_add,
            self.aligned_silva_seqs_art,
            self.v4_primers_dict["forward"],
            self.v4_primers_dict["reverse"],
            keeplength=False,
            keep_primer_location=True)

        obs_aln = skbio.io.read(str(obs_v4_nokeep_aln), into=skbio.TabularMSA,
                                constructor=skbio.DNA)
        exp_aln = skbio.io.read(str(self.silva_v4_trim_no_keeplen),
                                into=skbio.TabularMSA,
                                constructor=skbio.DNA)
        self.assertEqual(obs_aln, exp_aln)

    # test trimming when both primers are given and keeplength = True
    # tests against expected alignment length
    def test_trim_alignment_keeplen_true(self):
        obs_v4_keep_aln = _trim_alignment(
            mafft_add,
            self.aligned_silva_seqs_art,
            self.v4_primers_dict["forward"],
            self.v4_primers_dict["reverse"],
            keeplength=True,
            keep_primer_location=True)

        obs_aln = skbio.io.read(str(obs_v4_keep_aln), into=skbio.TabularMSA,
                                constructor=skbio.DNA)
        exp_aln = skbio.io.read(str(self.silva_v4_trim_keeplen),
                                into=skbio.TabularMSA,
                                constructor=skbio.DNA)
        self.assertEqual(obs_aln, exp_aln)

    # test trimming when both primers are given and keeplength = True
    # tests against expected alignment length w/o returning primer
    # sequence in the output.
    def test_trim_alignment_keeplen_true_no_keep_primers_false(self):
        obs_v4_keep_nokeep_primer_aln = _trim_alignment(
            mafft_add,
            self.aligned_silva_seqs_art,
            self.v4_primers_dict["forward"],
            self.v4_primers_dict["reverse"],
            keeplength=True,
            keep_primer_location=False)

        obs_aln = skbio.io.read(str(obs_v4_keep_nokeep_primer_aln),
                                into=skbio.TabularMSA,
                                constructor=skbio.DNA)
        exp_aln = skbio.io.read(str(self.silva_v4_trim_keeplen_nokeep_primer),
                                into=skbio.TabularMSA,
                                constructor=skbio.DNA)
        print(obs_aln, '\n', exp_aln)
        self.assertEqual(obs_aln, exp_aln)

    # test trimming when only fwd primer is given
    def test_trim_alignment_only_fwd(self):
        obs = _trim_alignment(
            self.fake_ctx.get_action(2),
            self.aligned_seqs_fasta,
            self.primers_dict["forward"],
            None,
            keep_primer_location=True)
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_only_fwd)

    # test trimming when only rev primer is given
    def test_trim_alignment_only_rev(self):
        obs = _trim_alignment(
            self.fake_ctx.get_action(3),
            self.aligned_seqs_fasta,
            None,
            self.primers_dict["reverse"],
            keep_primer_location=True)
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_only_rev)

    # test trimming when both positions are given
    def test_trim_alignment_by_positions(self):
        obs = _trim_alignment(
            self.fake_ctx.get_action(1),
            self.aligned_seqs_fasta,
            None, None, 8, 104,
            keep_primer_location=True)
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_both_primers)

    # test trimming when only start position is given
    def test_trim_alignment_by_position_left(self):
        obs = _trim_alignment(
            self.fake_ctx.get_action(2),
            self.aligned_seqs_fasta,
            None, None, 8, None,
            keep_primer_location=True)
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_only_fwd)

    # test trimming when only end position is given
    def test_trim_alignment_by_position_right(self):
        obs = _trim_alignment(
            self.fake_ctx.get_action(3),
            self.aligned_seqs_fasta,
            None, None, None, 104,
            keep_primer_location=True)
        obs_seqs = {seq.metadata['id']: str(seq)
                    for seq in obs.view(DNAIterator)}
        self.assertDictEqual(obs_seqs, self.exp_seqs_only_rev)
