# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
import qiime2
import pandas as pd
import pandas.testing as pdt

from rescript.dereplicate import _backfill_taxonomy
from rescript._utilities import _return_stripped_taxon_rank_list


import_data = qiime2.Artifact.import_data


class TestDerep(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.dereplicate = rescript.actions.dereplicate

        self.seqs = import_data(
            'FeatureData[Sequence]', self.get_data_path('derep-test.fasta'))
        self.fungi_seqs = import_data(
            'FeatureData[Sequence]',
            self.get_data_path('derep-test-fungi.fasta'))
        self.taxa = import_data(
            'FeatureData[Taxonomy]', self.get_data_path('derep-taxa.tsv'))
        self.fungi_taxa = import_data(
            'FeatureData[Taxonomy]',
            self.get_data_path('derep-test-fungi.tsv'))
        self.seqsnumericids = import_data(
            'FeatureData[Sequence]', self.get_data_path(
                'derep-test-numericIDs.fasta'))
        self.taxanumericids = import_data(
            'FeatureData[Taxonomy]', self.get_data_path(
                'derep-taxa-numericIDs.tsv'))

    def test_dereplicate_uniq(self):
        seqs, taxa, = self.dereplicate(
            self.seqs, self.taxa, mode='uniq', rank_handles=['disable'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__chondroitinus',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__brevis',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__damnosus',
            'B1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__vaginalis',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei',
            'C1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Pediococcus;s__acidilacti',
            'A3': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__alvei',
            'B2': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__casei',
            'C1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)
        # use derep_prefix=True;should still obtain same result if the prefix
        # seqs bear unique taxonomic labels, as seen in this test case
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='uniq',
                                       derep_prefix=True,
                                       rank_handles=['disable'])
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_uniq_99_perc(self):
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='uniq',
                                       perc_identity=0.99,
                                       rank_handles=['disable'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__chondroitinus',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__brevis',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__damnosus',
            'B1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__vaginalis',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei',
            'C1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Pediococcus;s__acidilacti',
            'A3': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__alvei',
            'B2': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__casei',
            'C1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)
        # use derep_prefix=True;should still obtain same result if the prefix
        # seqs bear unique taxonomic labels, as seen in this test case
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='uniq',
                                       perc_identity=0.99, derep_prefix=True,
                                       rank_handles=['disable'])
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_lca(self):
        seqs, taxa, = self.dereplicate(
            self.seqs, self.taxa, mode='lca', rank_handles=['disable'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__damnosus',
            'B1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__vaginalis',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei',
            'C1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_super_lca_majority(self):
        seqs, taxa, = self.dereplicate(
            self.seqs, self.taxa, mode='super')
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__alvei',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__casei',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__damnosus',
            'B1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__vaginalis',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei',
            'C1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Pediococcus;s__acidilacti'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_super_lca_majority_perc99(self):
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='super',
                                       perc_identity=0.99)
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__alvei',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__casei',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__acidilacti',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    # test that LCA taxonomy assignment works when derep_prefix=True
    # here derep_prefix + LCA leads to collapsed C-group seqs + LCA taxonomy
    def test_dereplicate_prefix_lca(self):
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='lca',
                                       derep_prefix=True,
                                       rank_handles=['disable'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae',
            'B1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__vaginalis',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    # test that LCA taxonomy assignment works when derep_prefix=True
    # seqs F05, F06, F07 are identical. Thus having the longest shared prefix.
    # The representative taxonomy should appear in the output as F05, as that
    # has the shorter taxonomy. Thus the LCA is set to that.
    # Here derep_prefix + LCA taxonomy
    def test_dereplicate_prefix_lca_fungi(self):
        seqs, taxa, = self.dereplicate(self.fungi_seqs, self.fungi_taxa,
                                       mode='lca',
                                       derep_prefix=True,
                                       rank_handles=['disable'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'F01': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota',
            'F02': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Debaryomycetaceae;'
                   'g__Candida-Lodderomyces_clade;s__Candida_parapsilosis',
            'F03': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Debaryomycetaceae;'
                   'g__Candida-Lodderomyces_clade;s__Candida_parapsilosis',
            'F04': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Debaryomycetaceae;'
                   'g__Candida-Lodderomyces_clade;s__Candida_parapsilosis',
            'F05': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Saccharomycetaceae;g__Lachancea',
            'F08': 'd__Eukaryota',
            'F09': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Saccharomycetaceae;'
                   'g__Saccharomyces;s__Saccharomyces_cerevisiae'
             }})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_lca_99_perc(self):
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='lca',
                                       perc_identity=0.99,
                                       rank_handles=['disable'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_majority(self):
        seqs, taxa, = self.dereplicate(
            self.seqs, self.taxa, mode='majority')
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__alvei',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__casei',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__damnosus',
            'B1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__vaginalis',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei',
            'C1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Pediococcus;s__acidilacti'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    # test that majority taxonomy assignment works when derep_prefix=True
    # all C-group seqs should be merged, and P. acidilacti is the majority
    def test_dereplicate_prefix_majority(self):
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='majority',
                                       derep_prefix=True)
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__alvei',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__casei',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__acidilacti',
            'B1a': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__vaginalis',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_majority_perc99(self):
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='majority',
                                       perc_identity=0.99)
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__alvei',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__casei',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Pediococcus;s__acidilacti',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_dereplicate_majority_perc98_fungi(self):
        seqs, taxa, = self.dereplicate(self.fungi_seqs, self.fungi_taxa,
                                       mode='majority',
                                       perc_identity=0.98,
                                       rank_handles=['disable'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'F01': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Debaryomycetaceae;'
                   'g__Candida-Lodderomyces_clade;s__Candida_parapsilosis',
            'F05': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Saccharomycetaceae;g__Lachancea;'
                   's__Lachancea_thermotolerans',
            'F09': 'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                   'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                   'o__Saccharomycetales;sf__Saccharomycetaceae;'
                   'g__Saccharomyces;s__Saccharomyces_cerevisiae'
                   }})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    # the above tests check actual derep functionality;this test just makes
    # sure that the same tests/modes above operate on numeric seq IDs, using
    # the same test data above (with numeric IDs).
    # See https://github.com/bokulich-lab/RESCRIPt/issues/49
    def test_dereplicate_numericIDs(self):
        self.dereplicate(self.seqsnumericids, self.taxanumericids, mode='uniq')
        self.assertTrue(True)
        self.dereplicate(self.seqsnumericids, self.taxanumericids, mode='lca')
        self.assertTrue(True)
        self.dereplicate(self.seqsnumericids, self.taxanumericids,
                         mode='majority')
        self.assertTrue(True)

    # Now test with backfilling. These parameters were chosen to set up a
    # variety of backfill levels.
    def test_dereplicate_lca_99_perc_backfill(self):
        # Note it is okay for inconsistent spacing around the ';' for
        # delimiting ranks. See `C1` for example. All that matters
        # is that `;` is delimiting the ranks.
        seqs, taxa, = self.dereplicate(self.seqs, self.taxa, mode='lca',
                                       perc_identity=0.99,
                                       rank_handles=['kingdom', 'phylum',
                                                     'class', 'order',
                                                     'family', 'genus',
                                                     'species'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'A1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                  'f__Paenibacillaceae;g__Paenibacillus;s__',
            'B1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__Lactobacillus;s__',
            'C1': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;'
                  'f__Lactobacillaceae;g__;s__',
            'B1b': 'k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales'
                   ';f__Lactobacillaceae;g__Lactobacillus;s__pseudocasei'}})
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_derep_lca_90_perc_backfill_fungi(self):
        seqs, taxa, = self.dereplicate(self.fungi_seqs, self.fungi_taxa,
                                       mode='lca',
                                       perc_identity=0.90,
                                       rank_handles=[
                                           'domain', 'superkingdom',
                                           'kingdom', 'subkingdom', 'phylum',
                                           'subphylum', 'class', 'order',
                                           'superfamily', 'genus', 'species'])
        exp_taxa = pd.DataFrame({'Taxon': {
            'F01': 'd__Eukaryota;sk__;k__;ks__;p__;ps__;c__;o__;sf__;g__;s__'}
            })
        pdt.assert_frame_equal(taxa.view(pd.DataFrame).sort_index(),
                               exp_taxa.sort_index(), check_names=False)
        pdt.assert_index_equal(seqs.view(pd.Series).sort_index().index,
                               exp_taxa.sort_index().index, check_names=False)

    def test_backfill_taxonomy(self):

        default_rank_handle = "d__; p__; c__; o__; f__; g__; s__"

        def _backfill_series(series, rank_handles=default_rank_handle):
            rank_handles = _return_stripped_taxon_rank_list(rank_handles)
            return series.apply(_backfill_taxonomy,
                                args=([rank_handles, ';']))

        taxa = self.taxa.view(pd.Series).sort_index()
        taxa = taxa.apply(lambda x: ';'.join(
                          _return_stripped_taxon_rank_list(x)))
        exp_taxa = taxa.copy()

        # note: taxonomy is unchanged if rank handles are shorter than taxon
        backfilled_taxa = _backfill_series(taxa, "my;taxonomy;is;too;short")
        pdt.assert_series_equal(backfilled_taxa, exp_taxa, check_names=False)
        # manually backfill to match expected
        exp_taxa.loc['C1b'] += ';g__;s__'
        # backfill with defaults
        backfilled_taxa = _backfill_series(taxa)
        pdt.assert_series_equal(backfilled_taxa, exp_taxa, check_names=False)
        # trim back arbitrarily to backfill again
        trimmed_taxa = backfilled_taxa.apply(
            lambda x: ';'.join(x.split(';')[:3]))
        # manually backfill
        exp_taxa = trimmed_taxa.apply(lambda x: x + ';o__;f__;g__;s__')
        backfilled_taxa = _backfill_series(trimmed_taxa)
        pdt.assert_series_equal(backfilled_taxa, exp_taxa, check_names=False)
        # backfill to root
        # note: taxon labels can never be empty, so this test covers cases
        # where there is no classification beyond root/domain/kingdom
        backfilled_taxa = _backfill_series(trimmed_taxa.apply(lambda x: 'd__'))
        exp_taxa = trimmed_taxa.apply(lambda x: ';'.join(
                                      _return_stripped_taxon_rank_list(
                                          default_rank_handle)))
        pdt.assert_series_equal(backfilled_taxa, exp_taxa, check_names=False)
        # backfill custom labels
        custom_rank_handles = "p;e;a;n;u;t;s"
        exp_taxa = trimmed_taxa.apply(lambda x: x + ';n;u;t;s')
        backfilled_taxa = _backfill_series(trimmed_taxa, custom_rank_handles)
        pdt.assert_series_equal(backfilled_taxa, exp_taxa, check_names=False)

    def test_backfill_long_taxonomy(self):
        default_rank_handle = ('d__; sk__; k__; ks__; p__; ps__; c__; o__; '
                               'sf__; g__; s__; ssb__; for__')
        rank_handles = _return_stripped_taxon_rank_list(default_rank_handle)
        taxa = ['d__Eukaryota; sk__Nucletmycea; k__Fungi; ks__Dikarya; '
                'p__Ascomycota',
                'd__Eukaryota',
                'd__Eukaryota; sk__Nucletmycea; k__Fungi; ks__Dikarya; '
                'p__Ascomycota; ps__Saccharomycotina; c__Saccharomycetes; '
                'o__Saccharomycetales; sf__Debaryomycetaceae; '
                'g__Candida-Lodderomyces_clade; s__Candida_parapsilosis']
        backfilled_taxa = [_backfill_taxonomy(t, rank_handles, ';')
                           for t in taxa]
        exp_taxa = ['d__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                    'p__Ascomycota;ps__;c__;o__;sf__;g__;s__;ssb__;for__',
                    'd__Eukaryota;sk__;k__;ks__;p__;ps__;c__;o__;sf__;g__;s__;'
                    'ssb__;for__',
                    'd__Eukaryota;sk__Nucletmycea;k__Fungi;ks__Dikarya;'
                    'p__Ascomycota;ps__Saccharomycotina;c__Saccharomycetes;'
                    'o__Saccharomycetales;sf__Debaryomycetaceae;'
                    'g__Candida-Lodderomyces_clade;s__Candida_parapsilosis;'
                    'ssb__;for__']
        self.assertEqual(backfilled_taxa, exp_taxa)
