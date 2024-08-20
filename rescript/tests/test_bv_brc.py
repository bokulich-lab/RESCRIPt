# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
from unittest.mock import Mock, patch, mock_open, MagicMock

import pandas as pd
from q2_types.feature_data import (MixedCaseDNAFASTAFormat,
                                   ProteinFASTAFormat,
                                   TSVTaxonomyDirectoryFormat)
from q2_types.genome_data import GenomeSequencesDirectoryFormat
from qiime2.plugin.testing import TestPluginBase

from rescript.bv_brc import fetch_genomes_bv_brc, fetch_metadata_bv_brc, \
    fetch_genome_features_bv_brc, fetch_taxonomy_bv_brc, id_list_handling, \
    error_handling, download_data, json_to_fasta, transform_taxonomy_df, \
    parse_lineage_names_with_ranks, convert_fasta_to_uppercase


class TestIDListHandling(TestPluginBase):
    package = 'rescript.tests'

    def test_error_both_parameters_given(self):
        with self.assertRaisesRegex(ValueError,
                                    "Parameters rql_query and ids can't be "
                                    "used simultaneously."):
            id_list_handling(rql_query="some_query",
                             ids=[1, 2, 3],
                             parameter_name="ids",
                             data_field="id")

    def test_error_neither_parameter_given(self):
        with self.assertRaisesRegex(ValueError,
                                    "At least one of the parameters rql_query "
                                    "and ids has to be given."):
            id_list_handling(rql_query="",
                             ids=[],
                             parameter_name="ids",
                             data_field="id")

    def test_correct_rql_query_generation(self):
        result = id_list_handling(
            rql_query="",
            ids=[1, 2, 3],
            parameter_name="ids",
            data_field="id")
        expected_query = "in(id,(1,2,3))"
        self.assertEqual(result, expected_query)


class TestErrorHandling(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()
        self.response = Mock()

    def test_no_data_found(self):
        self.response.text = "[]"

        with self.assertRaisesRegex(ValueError, "No data"):
            error_handling(self.response, data_type="genome")

    def test_database_error_occurred_undefined_field_object(self):
        self.response.text = ('A Database Error Occured: {"code": 400, '
                              '"msg": "undefined field object in RQL"}')

        with self.assertRaisesRegex(ValueError, "undefined field object"):
            error_handling(self.response, data_type="genome")

    def test_database_error_occurred_undefined_field(self):
        self.response.text = ('A Database Error Occured: {"code": 400, '
                              '"msg": "undefined field"}')

        with self.assertRaisesRegex(ValueError, "undefined field"):
            error_handling(self.response, data_type="genome")

    def test_database_error_occurred_general_error(self):
        self.response.text = ('A Database Error Occured: {"code": 500, "msg": '
                              '"Internal Server Error"}')

        with self.assertRaisesRegex(ValueError, "Internal Server Error"):
            error_handling(self.response, data_type="genome")

    def test_unhandled_response(self):
        self.response.text = "Unexpected error message"

        with self.assertRaisesRegex(ValueError, "Unexpected error"):
            error_handling(self.response, data_type="genome")


class TestDownloadData(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.requests.get')
    @patch('rescript.bv_brc.error_handling')
    def test_download_data_success(self, mock_error_handling,
                                   mock_requests_get):
        # Mock the requests.get response for a successful request
        mock_response = Mock()
        mock_response.status_code = 200
        mock_requests_get.return_value = mock_response

        url = "http://example.com/data"
        data_type = "some_type"

        result = download_data(url, data_type)

        mock_requests_get.assert_called_once_with(url)
        self.assertEqual(result, mock_response)

    @patch('rescript.bv_brc.requests.get')
    @patch('rescript.bv_brc.error_handling')
    def test_download_data_error_400(self, mock_error_handling,
                                     mock_requests_get):
        # Mock the requests.get response for a 400 Bad Request
        mock_response = Mock()
        mock_response.status_code = 400
        mock_requests_get.return_value = mock_response

        url = "http://example.com/data"
        data_type = "some_type"

        download_data(url, data_type)

        mock_requests_get.assert_called_once_with(url)
        mock_error_handling.assert_called_once_with(mock_response, data_type)

    @patch('rescript.bv_brc.requests.get')
    @patch('rescript.bv_brc.error_handling')
    def test_download_data_other_error(self, mock_error_handling,
                                       mock_requests_get):
        # Mock the requests.get response for any other error
        mock_response = Mock()
        mock_response.status_code = 500
        mock_response.text = "Server Error"
        mock_requests_get.return_value = mock_response

        url = "http://example.com/data"
        data_type = "some_type"

        with self.assertRaisesRegex(ValueError, "Server Error"):
            download_data(url, data_type)

        mock_requests_get.assert_called_once_with(url)
        mock_error_handling.assert_not_called()


class TestJsonToFasta(TestPluginBase):
    package = 'rescript.tests'

    def setUp(self):
        super().setUp()

        self.json_input_1 = [
            {
                "genome_id": "genome1",
                "accession": "acc1",
                "description": "desc1",
                "genome_name": "genome_name1",
                "sequence": "ATGC"
            }
        ]

        self.json_input_2 = [
            {
                "genome_id": "genome2",
                "accession": "acc2",
                "description": "desc2",
                "genome_name": "genome_name2",
                "sequence": "CGTA"
            }
        ]

    @patch('rescript.bv_brc.open', new_callable=mock_open)
    def test_json_to_fasta_single_genome(self, mock_file):
        json_to_fasta(self.json_input_1, "/fake/dir")

        # Expected FASTA content
        expected_fasta = ">accn|acc1   desc1   [genome_name1 | genome1]\nATGC"

        # Check if the file was created with the correct path and content
        mock_file.assert_called_once_with("/fake/dir/genome1.fasta", 'w')
        mock_file().write.assert_called_once_with(expected_fasta)

    @patch('rescript.bv_brc.open', new_callable=mock_open)
    def test_json_to_fasta_multiple_genomes(self, mock_file):
        json_to_fasta(self.json_input_1 + self.json_input_2, "/fake/dir")

        # Expected FASTA content
        expected_fasta_genome1 = (">accn|acc1   desc1   [genome_name1 | "
                                  "genome1]\nATGC")
        expected_fasta_genome2 = (">accn|acc2   desc2   [genome_name2 | "
                                  "genome2]\nCGTA")

        # Check if the files were created with the correct path and content
        mock_file().write.assert_any_call(expected_fasta_genome1)
        mock_file().write.assert_any_call(expected_fasta_genome2)

    @patch('rescript.bv_brc.open', new_callable=mock_open)
    def test_json_to_fasta_multiple_sequences_same_genome(self, mock_file):
        json_to_fasta(self.json_input_1 + self.json_input_1, "/fake/dir")

        # Expected FASTA content
        expected_fasta = (">accn|acc1   desc1   [genome_name1 | "
                          "genome1]\nATGC\n"
                          ">accn|acc1   desc1   [genome_name1 | "
                          "genome1]\nATGC")

        # Check if the file was created with the correct path and content
        mock_file.assert_called_once_with("/fake/dir/genome1.fasta", 'w')
        mock_file().write.assert_called_once_with(expected_fasta)


class TestFetchGenomeFeaturesBVBR(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.id_list_handling')
    @patch.object(MixedCaseDNAFASTAFormat, 'open')
    @patch.object(ProteinFASTAFormat, 'open')
    def test_fetch_genome_features_bv_brc(
            self, mock_protein_open, mock_genes_open, mock_id_list_handling,
            mock_download_data
    ):
        # Mock the id_list_handling function
        mock_id_list_handling.return_value = ("in(feature_id, "
                                              "(feature1,feature2))")

        # Mock the responses from download_data
        mock_genes_response = MagicMock()
        mock_genes_response.text = ">gene1\nATGC\n>gene2\nATGC"
        mock_proteins_response = MagicMock()
        mock_proteins_response.text = (
            ">protein1\nMVLSPADKTNVK\n>protein2\nMVLSPADKTNVK")
        mock_download_data.side_effect = [mock_genes_response,
                                          mock_proteins_response]

        # Mock file write actions
        mock_genes_file = MagicMock()
        mock_protein_file = MagicMock()
        mock_genes_open.return_value.__enter__.return_value = mock_genes_file
        mock_protein_open.return_value.__enter__.return_value = (
            mock_protein_file)

        # Call the function
        genes, proteins = fetch_genome_features_bv_brc(
            rql_query="in(feature_id, (feature1,feature2))",
            feature_ids=["feature1", "feature2"]
        )

        # Assertions
        mock_id_list_handling.assert_called_once_with(
            rql_query="in(feature_id, (feature1,feature2))",
            ids=["feature1", "feature2"],
            parameter_name="feature_ids",
            data_field="feature_id"
        )

        mock_download_data.assert_any_call(
            url="https://www.bv-brc.org/api/genome_feature/?in(feature_id, "
                "(feature1,feature2))&http_accept=application/dna+fasta",
            data_type="genome_feature"
        )

        mock_download_data.assert_any_call(
            url="https://www.bv-brc.org/api/genome_feature/?in(feature_id, "
                "(feature1,feature2))&http_accept=application/protein+fasta",
            data_type="genome_feature"
        )

        # Check that the correct data is written to the correct files
        mock_genes_file.write.assert_called_once_with(
            ">gene1\nATGC\n>gene2\nATGC")
        mock_protein_file.write.assert_called_once_with(
            ">protein1\nMVLSPADKTNVK\n>protein2\nMVLSPADKTNVK")

        self.assertIsInstance(genes, MixedCaseDNAFASTAFormat)
        self.assertIsInstance(proteins, ProteinFASTAFormat)

    def test_convert_fasta_to_uppercase(self):
        input_fasta = ">header1\natgca\ngtacg\n>header2\nttgaa\ncctg"
        expected_output = ">header1\nATGCA\nGTACG\n>header2\nTTGAA\nCCTG"

        result = convert_fasta_to_uppercase(input_fasta)
        self.assertEqual(result, expected_output)


class TestFetchGenomesBVBRC(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.json_to_fasta')
    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.id_list_handling')
    def test_fetch_genomes_bv_brc(
            self, mock_id_list_handling, mock_download_data, mock_json_to_fasta
    ):
        # Mock the id_list_handling function
        mock_id_list_handling.return_value = "genome_id=in(genome1,genome2)"

        # Mock the download_data response
        mock_response = MagicMock()
        mock_response.json.return_value = {'genomes': ['genome_data']}
        mock_download_data.return_value = mock_response

        # Call the function
        genomes = fetch_genomes_bv_brc(
            rql_query="genome_id=in(genome1,genome2)",
            genome_ids=["genome1", "genome2"]
        )

        # Assertions
        mock_id_list_handling.assert_called_once_with(
            rql_query="genome_id=in(genome1,genome2)",
            ids=["genome1", "genome2"],
            parameter_name="genome_ids",
            data_field="genome_id"
        )

        mock_download_data.assert_called_once_with(
            url="https://www.bv-brc.org/api/genome_sequence/"
                "?genome_id=in(genome1,genome2)",
            data_type="genome_sequence"
        )

        mock_json_to_fasta.assert_called_once_with(
            {'genomes': ['genome_data']},
            str(genomes)
        )

        self.assertIsInstance(genomes, GenomeSequencesDirectoryFormat)


class TestFetchMetadataBVBR(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.qiime2.Metadata')
    @patch('rescript.bv_brc.pd.read_csv')
    @patch('rescript.bv_brc.download_data')
    def test_fetch_metadata_bv_brc(self, mock_download_data,
                                   mock_read_csv, mock_metadata):
        # Mock the download_data response
        mock_response = MagicMock()
        mock_response.text = (
            "id\tcolumn1\tcolumn2\n1\tdata1\tdata2\n2\tdata3\tdata4")
        mock_download_data.return_value = mock_response

        # Mock the pandas read_csv return value
        mock_df = pd.DataFrame({
            'column1': ['data1', 'data3'],
            'column2': ['data2', 'data4']
        }, index=pd.Index(['1', '2'], name='id'))
        mock_read_csv.return_value = mock_df

        # Mock qiime2.Metadata creation
        mock_metadata_instance = MagicMock()
        mock_metadata.return_value = mock_metadata_instance

        # Call the function
        fetch_metadata_bv_brc(
            data_type="genome",
            rql_query="genome_id=in(1,2)"
        )

        # Assertions
        mock_download_data.assert_called_once_with(
            url="https://www.bv-brc.org/api/genome/"
                "?genome_id=in(1,2)&http_accept=text/tsv",
            data_type="genome"
        )

        mock_read_csv.assert_called_once()
        args, kwargs = mock_read_csv.call_args
        self.assertEqual(kwargs['sep'], '\t')

        self.assertEqual(args[0].getvalue(), "id\tcolumn1\tcolumn2\n1\tdata1"
                                             "\tdata2\n2\tdata3\tdata4")

        mock_metadata.assert_called_once_with(mock_df)


class TestFetchTaxonomyBVBR(TestPluginBase):
    package = 'rescript.tests'

    @patch('pandas.DataFrame.to_csv')
    @patch('rescript.bv_brc.transform_taxonomy_df')
    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.pd.read_csv')
    @patch('rescript.bv_brc.id_list_handling')
    def test_fetch_taxonomy_bv_brc(
            self, mock_id_list_handling, mock_read_csv, mock_download_data,
            mock_transform_taxonomy_df, mock_to_csv
    ):
        # Mock the id_list_handling function
        mock_id_list_handling.return_value = "taxon_id=in(taxon1,taxon2)"

        # Mock the download_data response
        mock_response = MagicMock()
        mock_response.text = (
            "id\trank1\trank2\n1\tdata1\tdata2\n2\tdata3\tdata4")
        mock_download_data.return_value = mock_response

        # Prepare mocks for file output
        with patch('builtins.open', unittest.mock.mock_open()):
            directory = fetch_taxonomy_bv_brc(
                rql_query="taxon_id=in(taxon1,taxon2)",
                ranks=['rank1', 'rank2'],
                taxon_ids=["taxon1", "taxon2"]
            )

            # Assertions
            mock_id_list_handling.assert_called_once_with(
                rql_query="taxon_id=in(taxon1,taxon2)",
                ids=["taxon1", "taxon2"],
                parameter_name="taxon_ids",
                data_field="taxon_id"
            )

            mock_download_data.assert_called_once_with(
                url="https://www.bv-brc.org/api/taxonomy/"
                    "?taxon_id=in(taxon1,taxon2)&http_accept=text/tsv",
                data_type="taxonomy"
            )

            self.assertIsInstance(directory, TSVTaxonomyDirectoryFormat)

    @patch('rescript.bv_brc.parse_lineage_names_with_ranks')
    def test_transform_taxonomy_df(self, mock_parse_lineage_names_with_ranks):
        # Mock the parse_lineage_names_with_ranks function
        mock_parse_lineage_names_with_ranks.side_effect = \
            lambda lineage_names, lineage_ranks, ranks: "Mocked Taxon"

        # Create a sample DataFrame
        df = pd.DataFrame({
            'taxon_id': ['taxon1', 'taxon2'],
            'lineage_names': ['name1;name2', 'name3;name4'],
            'lineage_ranks': ['rank1;rank2', 'rank3;rank4']
        })

        ranks = ['rank1', 'rank2', 'rank3']

        # Call the function
        result_df = transform_taxonomy_df(df, ranks)

        # Expected DataFrame after transformation
        expected_df = pd.DataFrame({
            'Feature ID': ['taxon1', 'taxon2'],
            'Taxon': ['Mocked Taxon', 'Mocked Taxon']
        }).set_index('Feature ID')

        # Assert that the result matches the expected DataFrame
        pd.testing.assert_frame_equal(result_df, expected_df)

    def test_parse_with_missing_ranks(self):
        lineage_names = "Bacteria;Proteobacteria;Enterobacteriaceae"
        lineage_ranks = "kingdom;phylum;family"
        ranks = ['kingdom', 'phylum', 'class', 'order', 'genus', 'species']

        result = parse_lineage_names_with_ranks(lineage_names, lineage_ranks,
                                                ranks)
        expected = "k__Bacteria; p__Proteobacteria; c__; o__; g__; s__"

        self.assertEqual(result, expected)

    def test_parse_with_no_ranks_provided(self):
        lineage_names = (
            "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;"
            "Enterobacteriaceae;Escherichia;coli")
        lineage_ranks = "kingdom;phylum;class;order;family;genus;species"
        ranks = None  # Should fall back to _default_ranks

        result = parse_lineage_names_with_ranks(lineage_names, lineage_ranks,
                                                ranks)
        expected = ("k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; "
                    "o__Enterobacterales; f__Enterobacteriaceae; "
                    "g__Escherichia; s__coli")

        self.assertEqual(result, expected)
