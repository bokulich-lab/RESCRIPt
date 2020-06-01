# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pkg_resources

import qiime2
from q2_types.feature_data import DNAFASTAFormat, DNAIterator

from .test_types_formats_transformers import RescriptTypesTestPluginBase
from qiime2.plugins import rescript


class TestReverseTranscribe(RescriptTypesTestPluginBase):

    def setUp(self):
        super().setUp()
        dna_path = pkg_resources.resource_filename(
            'rescript.tests', 'data/derep-test.fasta')
        self.dna_seqs = DNAFASTAFormat(dna_path, mode='r').view(DNAIterator)

    # this is technically already tested as a transformer and is a really
    # shallow utility, but we'll just test again as a stand-alone for fun.
    def test_reverse_transcribe(self):
        rna_seqs = qiime2.Artifact.import_data(
            'FeatureData[RNASequence]',
            self.get_data_path('derep-test-rna.fasta'))
        obs, = rescript.actions.reverse_transcribe(rna_seqs)
        exp = self.dna_seqs
        for observed, expected in zip(obs.view(DNAIterator), exp):
            self.assertEqual(observed, expected)
