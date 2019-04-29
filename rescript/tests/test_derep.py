# ----------------------------------------------------------------------------
# Copyright (c) 2019--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
import qiime2
import pandas as pd
import pandas.util.testing as pdt


import_data = qiime2.Artifact.import_data


class TestDerep(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.dereplicate = rescript.actions.dereplicate

        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))
        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))

    def test_dereplicate(self):
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa)
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; '
                  'f__Paenibacillaceae; g__Paenibacillus; s__chondroitinus',
            'B1': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales;'
                  ' f__Lactobacillaceae; g__Lactobacillus; s__brevis',
            'C1': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales;'
                  ' f__Lactobacillaceae; g__Pediococcus; s__damnosus',
            'B1a': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
                   '; f__Lactobacillaceae; g__Lactobacillus; s__vaginalis',
            'B1b': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
                   '; f__Lactobacillaceae; g__Lactobacillus; s__pseudocasei',
            'C1a': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales'
                   '; f__Lactobacillaceae; g__Pediococcus; s__acidilacti',
            'A3': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; '
                  'f__Paenibacillaceae; g__Paenibacillus; s__alvei',
            'B2': 'k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales;'
                  ' f__Lactobacillaceae; g__Lactobacillus; s__casei'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)
