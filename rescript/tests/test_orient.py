# ----------------------------------------------------------------------------
# Copyright (c) 2019-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import os
import gzip
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAIterator
from qiime2.plugins import rescript
from rescript.orient import _add_optional_parameters
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt)
from pandas.testing import assert_index_equal

import_data = qiime2.Artifact.import_data


class TestOrientSeqs(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.seqs = \
            import_data('FeatureData[Sequence]',
                        self.get_data_path('mixed-orientations.fasta'))
        self.ref = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))

        self.rc = \
            import_data('FeatureData[Sequence]',
                        self.get_data_path('mixed-orientations-rc.fasta')
                        )

    def test_reorient_default(self):
        # this test checks that expected IDs AND reoriented seqs are returned
        reoriented, unmatched, = rescript.actions.orient_seqs(
            sequences=self.seqs, reference_sequences=self.ref)
        # all input seqs should pass, minus the two junk seqs. Those that pass
        # are copies of the ref seqs with some mismatches/gaps inserted
        # so the oriented seq ids should match the ref ids.
        reoriented_seqs = {seq.metadata['id']: str(seq)
                           for seq in reoriented.view(DNAIterator)}
        exp_reoriented_seqs = {seq.metadata['id']: str(seq)
                               for seq in self.ref.view(DNAIterator)}
        self.assertEqual(reoriented_seqs.keys(), exp_reoriented_seqs.keys())
        # the oriented seqs should also match the ref seqs, except for those
        # with some mismatches inserted...
        exp_reoriented_seqs['A1'] = (
            'AAAAAAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCG'
            'CGTAGGTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACT'
            'GGCGAGCTAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCG'
            'GGAAGAATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGG'
            'GGAGCAAACAGGATTAGATACCCTGGTAGTC')
        exp_reoriented_seqs['A2'] = (
            'AGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAG'
            'GTGGCCCGTTAAGTGGCTGGTGAAATCCCGGGGCTCAACTCCGGGGCTGCCGGTCAGACTGGCGA'
            'GGAAGAGCACGGTAGGGGCAGATGGAATTCCCGGTGTAGCGGTGGAATGCGTAGATATCGGGAAG'
            'AATACCAGTGGCGAAGGCGTTCTGCTGGGCCGTTGCTGACACTGAGGCGCGACAGCGTGGGGAGC'
            'AAACAGGATTAGATACCCTGGTAGTC')
        exp_reoriented_seqs['A3'] = exp_reoriented_seqs['A3'][:-21] + 'T' * 21
        self.assertEqual(
            set(reoriented_seqs.values()), set(exp_reoriented_seqs.values()))
        # make sure the junk sequences don't match anything and are discarded
        unmatched_ids = {
            seq.metadata['id'] for seq in unmatched.view(DNAIterator)}
        exp_unmatched_ids = {'JUNK', 'MOREJUNK'}
        self.assertEqual(unmatched_ids, exp_unmatched_ids)

    def test_reorient_no_ref(self):
        reoriented, unmatched = rescript.actions.orient_seqs(
            sequences=self.seqs, reference_sequences=None,
            )
        unmatched_ids = {seq.metadata['id']
                         for seq in unmatched.view(DNAIterator)}
        self.assertEqual(unmatched_ids, set([]))
        exp_seqs = [seq for seq in self.rc.view(DNAIterator)]
        test_seqs = [seq for seq in reoriented.view(DNAIterator)]
        for exp, test in zip(*(exp_seqs, test_seqs)):
            self.assertEqual(str(exp), str(test))
            self.assertEqual(exp.metadata['id'], test.metadata['id'])

    def test_add_optional_parameters(self):
        expected = [
            "--dbmask", "dust",
            "--relabel", "new_id",
            "--relabel_keep",
            "--relabel_md5",
            "--relabel_self",
            "--relabel_sha1",
            "--sizein",
            "--sizeout",
        ]
        result = []
        _add_optional_parameters(
            result,
            dbmask="dust",
            relabel="new_id",
            relabel_keep=True,
            relabel_md5=True,
            relabel_self=True,
            relabel_sha1=True,
            sizein=True,
            sizeout=True,
        )
        self.assertEqual(result, expected)


class TestOrientReads(TestPluginBase):
    package = 'rescript.tests'

    # Check that contents of fastq files match per sample
    def _validate_fastqs_equal(self, obs, exp):
        obs = obs.view(CasavaOneEightSingleLanePerSampleDirFmt)
        exp = exp.view(CasavaOneEightSingleLanePerSampleDirFmt)

        # obs and exp manifest indices should be identical
        assert_index_equal(obs.manifest.index, exp.manifest.index)

        # Iterate over each sample, side-by-side
        def _compare_contents(exp, obs, _seq_fp):
            exp_path = str(exp.path / _seq_fp)
            obs_path = str(obs.path / _seq_fp)
            with gzip.open(str(exp_path), 'r') as e:
                with gzip.open(str(obs_path), 'r') as o:
                    self.assertEqual(e.read(), o.read())

        paired = exp.manifest.reverse.any()
        for _id in exp.manifest.itertuples():
            _seq_fp = os.path.basename(_id.forward)
            _compare_contents(exp, obs, _seq_fp)
            if paired:
                _seq_fp = os.path.basename(_id.reverse)
                _compare_contents(exp, obs, _seq_fp)

    def setUp(self):
        super().setUp()
        self.seqs = \
            import_data('SampleData[PairedEndSequencesWithQuality]',
                        self.get_data_path('paired_end_data'))
        self.seqs_joined = \
            import_data('SampleData[JoinedSequencesWithQuality]',
                        self.get_data_path('single_end_data'))
        self.ref = import_data(
            'FeatureData[Sequence]', self.get_data_path('Human-Kneecap.fasta'))

    # this test checks that expected IDs AND reoriented seqs are returned
    # the test data consists of 3 read pairs in the correct orientation
    # and one pair that has R1/R2 switched to simulate mixed-orientation
    # and one pair should have unknown orientation/no similarity to any refseq
    # and thus will be removed.
    def test_reorient_default(self):
        expected = \
            import_data('SampleData[PairedEndSequencesWithQuality]',
                        self.get_data_path('paired_end_data_expected'))
        reoriented, unmatched, = rescript.actions.orient_reads(
            sequences=self.seqs, reference_sequences=self.ref)
        self._validate_fastqs_equal(reoriented, expected)

    # and again but for pre-joined reads
    # in this case one read is in the wrong orientation and will be
    # reoriented. One is garbage with no similarity to refseqs and will be
    # removed. Technically these data are single-end reads but for
    # practical reasons they can masquerade as pre-joined reads just as
    # well for our purposes (since we do not actually care about seq
    # topology or quality info, we just care about orientation!)
    def test_reorient_default_joined(self):
        expected = \
            import_data('SampleData[JoinedSequencesWithQuality]',
                        self.get_data_path('single_end_data_expected'))
        reoriented, unmatched, = rescript.actions.orient_reads(
            sequences=self.seqs_joined, reference_sequences=self.ref)
        self._validate_fastqs_equal(reoriented, expected)
