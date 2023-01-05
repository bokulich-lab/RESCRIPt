# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import pandas as pd
import pandas.util.testing as pdt
import numpy as np
import tempfile
import pkg_resources

from qiime2.plugin import ValidationError
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import (
    FeatureData, DNAIterator, AlignedDNAFASTAFormat, AlignedDNAIterator,
    DNAFASTAFormat, RNAFASTAFormat)

from rescript._utilities import _read_fasta, _read_dna_alignment_fasta
from rescript.types import (SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat,
                            SILVATaxidMapFormat, SILVATaxidMapDirectoryFormat,
                            SILVATaxonomy, SILVATaxidMap)


class RescriptTypesTestPluginBase(TestPluginBase):
    package = 'rescript.types.tests'

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.TemporaryDirectory(
            prefix='rescript-test-temp-')

    def tearDown(self):
        self.temp_dir.cleanup()

    def get_data_path(self, filename):
        return pkg_resources.resource_filename(self.package,
                                               'data/%s' % filename)


class TestSILVATypesAndFormats(RescriptTypesTestPluginBase):

    def test_silva_taxonomy_format_validate_positive(self):
        filepath = self.get_data_path('silva_taxa.tsv')
        format = SILVATaxonomyFormat(filepath, mode='r')
        # These should both just succeed
        format.validate('min')
        format.validate('max')

    def test_silva_taxonomy_format_validate_negative_incorrect_row_count(self):
        filepath = self.get_data_path('silva_taxamap.tsv')
        format = SILVATaxonomyFormat(filepath, mode='r')
        with self.assertRaisesRegex(ValidationError, 'Expected.*5 fields'):
            format.validate()

    def test_silva_taxonomy_format_validate_negative_invalid_taxonomy(self):
        filepath = self.get_data_path('silva_taxa.tsv')
        # index_col=1 sets the wrong index, creating an invalid taxonomy
        badtax = pd.read_csv(filepath, sep='\t', index_col=1, header=None)
        junkpath = os.path.join(self.temp_dir.name, 'junk.tsv')
        badtax.to_csv(junkpath, sep='\t', header=False)
        format = SILVATaxonomyFormat(junkpath, mode='r')
        with self.assertRaisesRegex(
                ValidationError, 'not in SILVA taxonomy format'):
            format.validate()

    def test_silva_taxonomy_format_validate_negative_nonnumeric_values(self):
        filepath = self.get_data_path('silva_taxa.tsv')
        badtax = pd.read_csv(filepath, sep='\t', index_col=0, header=None)
        # turn the taxids to a non-numeric column
        badtax[1] = 'blarg'
        junkpath = os.path.join(self.temp_dir.name, 'junk.tsv')
        badtax.to_csv(junkpath, sep='\t', header=False)
        format = SILVATaxonomyFormat(junkpath, mode='r')
        with self.assertRaisesRegex(ValidationError, 'non-numeric value'):
            format.validate()

    def test_silva_taxonomy_format_validate_negative_empty_file(self):
        filepath = self.get_data_path('empty.tsv')
        format = SILVATaxonomyFormat(filepath, mode='r')
        with self.assertRaisesRegex(
                ValidationError, 'must be at least one.*record'):
            format.validate()

    def test_silva_taxidmap_format_validate_positive(self):
        filepath = self.get_data_path('silva_taxamap.tsv')
        format = SILVATaxidMapFormat(filepath, mode='r')
        # These should both just succeed
        format.validate('min')
        format.validate('max')

    def test_silva_taxidmap_format_validate_positive_version132(self):
        filepath = self.get_data_path('silva_taxamap_v132.tsv')
        format = SILVATaxidMapFormat(filepath, mode='r')
        # These should both just succeed
        format.validate('min')
        format.validate('max')

    def test_silva_taxidmap_format_validate_negative_bad_header(self):
        filepath = self.get_data_path('silva_taxa.tsv')
        format = SILVATaxidMapFormat(filepath, mode='r')
        with self.assertRaisesRegex(
                ValidationError, 'Header line does not match SILVA format'):
            format.validate()

    def test_silva_taxidmap_format_validate_negative_row_count(self):
        filepath = self.get_data_path('trash_taxamap.tsv')
        format = SILVATaxidMapFormat(filepath, mode='r')
        with self.assertRaisesRegex(ValidationError, 'Expected.*6 fields'):
            format.validate()

    def test_silva_taxidmap_format_validate_negative_invalid_taxonomy(self):
        filepath = self.get_data_path('silva_taxamap.tsv')
        badtax = pd.read_csv(filepath, sep='\t', index_col=0, header=0)
        # replace the taxonomy with non-semicolon-delimited strings
        badtax['path'] = 'foobar'
        junkpath = os.path.join(self.temp_dir.name, 'junk.tsv')
        badtax.to_csv(junkpath, sep='\t', header=True)
        format = SILVATaxidMapFormat(junkpath, mode='r')
        with self.assertRaisesRegex(
                ValidationError, 'not in SILVA taxonomy format'):
            format.validate()

    def test_silva_taxidmap_format_validate_negative_nonnumeric_values(self):
        filepath = self.get_data_path('silva_taxamap.tsv')
        badtax = pd.read_csv(filepath, sep='\t', index_col=0, header=0)
        # turn the taxids to a non-numeric column
        badtax['taxid'] = 'peanuts'
        junkpath = os.path.join(self.temp_dir.name, 'junk.tsv')
        badtax.to_csv(junkpath, sep='\t', header=True)
        format = SILVATaxidMapFormat(junkpath, mode='r')
        with self.assertRaisesRegex(ValidationError, 'non-numeric value'):
            format.validate()

    def test_silva_taxidmap_format_validate_negative_empty_file(self):
        filepath = self.get_data_path('empty_with_header.tsv')
        format = SILVATaxidMapFormat(filepath, mode='r')
        with self.assertRaisesRegex(
                ValidationError, 'must be at least one.*record'):
            format.validate()


class TestRegistrations(RescriptTypesTestPluginBase):
    def test_silva_taxonomy_semantic_type_registration(self):
        self.assertRegisteredSemanticType(SILVATaxonomy)

    def test_silva_taxid_map_type_registration(self):
        self.assertRegisteredSemanticType(SILVATaxidMap)

    def test_silva_taxonomy_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            FeatureData[SILVATaxonomy], SILVATaxonomyDirectoryFormat)

    def test_silva_taxid_map_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            FeatureData[SILVATaxidMap], SILVATaxidMapDirectoryFormat)


class TestSILVATransformers(RescriptTypesTestPluginBase):

    def setUp(self):
        super().setUp()
        self.tax = pd.DataFrame(
            {1: {'Archaea;': '2'}, 2: {'Archaea;': 'domain'},
             3: {'Archaea;': np.nan}, 4: {'Archaea;': np.nan}})
        self.tax.index.name = 'id'

        self.taxmap = pd.DataFrame({
            'start': {'A16379': '1'}, 'stop': {'A16379': '1485'},
            'path': {'A16379': 'Bacteria;Proteobacteria;Gammaproteobacteria;'
                               'Pasteurellales;Pasteurellaceae;Haemophilus;'},
            'organism_name': {'A16379': '[Haemophilus] ducreyi'},
            'taxid': {'A16379': '3698'}})
        self.taxmap.index.name = 'id'

    def test_pd_dataframe_to_silva_taxonomy_format(self):
        transformer = self.get_transformer(pd.DataFrame, SILVATaxonomyFormat)
        obs = transformer(self.tax)
        obs = pd.read_csv(
            str(obs), sep='\t', header=None, index_col=0, dtype='str')
        obs.index.name = 'id'
        pdt.assert_frame_equal(self.tax, obs, check_dtype=False)

    def test_silva_taxonomy_format_to_pd_dataframe(self):
        _, obs = self.transform_format(
            SILVATaxonomyFormat, pd.DataFrame, 'silva_taxa.tsv')
        pdt.assert_frame_equal(self.tax, obs[:1], check_dtype=False)

    def test_pd_dataframe_to_silva_taxidmap_format(self):
        transformer = self.get_transformer(pd.DataFrame, SILVATaxidMapFormat)
        obs = transformer(self.taxmap)
        obs = pd.read_csv(
            str(obs), sep='\t', header=0, index_col=0, dtype='str')
        obs.index.name = 'id'
        pdt.assert_frame_equal(self.taxmap, obs, check_dtype=False)

    def test_silva_taxidmap_format_to_pd_dataframe(self):
        _, obs = self.transform_format(
            SILVATaxidMapFormat, pd.DataFrame, 'silva_taxamap.tsv')
        pdt.assert_frame_equal(self.taxmap, obs[:1], check_dtype=False)


class TestRNATransformers(RescriptTypesTestPluginBase):

    def setUp(self):
        super().setUp()
        dna_path = pkg_resources.resource_filename(
            'rescript.tests', 'data/derep-test.fasta')
        self.dna_seqs = DNAFASTAFormat(dna_path, mode='r').view(DNAIterator)

    def test_rna_fasta_format_to_dna_fasta_format(self):
        # transform RNA to DNA (reverse transcribe)
        input, obs = self.transform_format(
            RNAFASTAFormat, DNAFASTAFormat, 'derep-test-rna.fasta')
        self.assertIsInstance(obs, DNAFASTAFormat)
        # load expected DNA seqs (already reverse transcribed)
        exp = self.dna_seqs
        # convert to DNAIterator to iterate over seqs, confirm that
        # reverse transcription occurred as expected.
        obs = _read_fasta(str(obs))
        for observed, expected in zip(obs, exp):
            self.assertEqual(observed, expected)

    def test_rna_fasta_format_to_dna_iterator(self):
        input, obs = self.transform_format(
            RNAFASTAFormat, DNAIterator, 'derep-test-rna.fasta')
        exp = self.dna_seqs
        for observed, expected in zip(obs, exp):
            self.assertEqual(observed, expected)


class TestDNAIteratorTransformers(RescriptTypesTestPluginBase):
    def setUp(self):
        super().setUp()
        self.aligned_dna_path = pkg_resources.resource_filename(
            'rescript.tests', 'data/trim-test-alignment.fasta')
        self.aligned_dna_seqs = AlignedDNAFASTAFormat(
            self.aligned_dna_path, mode='r').view(AlignedDNAIterator)

    def test_aligned_dna_iterator_to_dna_fasta(self):
        transformer = self.get_transformer(DNAIterator, AlignedDNAFASTAFormat)
        obs = transformer(self.aligned_dna_seqs)
        obs = _read_dna_alignment_fasta(str(obs))
        exp = _read_dna_alignment_fasta(self.aligned_dna_path)
        for observed, expected in zip(obs, exp):
            self.assertEqual(observed, expected)
