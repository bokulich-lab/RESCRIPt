# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAIterator
from qiime2.plugins import rescript

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
