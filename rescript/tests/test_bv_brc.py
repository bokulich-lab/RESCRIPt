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
import qiime2
from qiime2.plugin.testing import TestPluginBase

from rescript.bv_brc import get_bv_brc_genomes, get_bv_brc_metadata, \
    get_bv_brc_genome_features, parameter_validation, \
    error_handling, download_data, create_genome_fasta, \
    create_taxonomy_entry, get_loci, read_tsv_data_with_dtypes, process_loci, \
    get_sequences, get_taxonomy, create_taxonomy, get_genome_sequences


class TestParameterValidation(TestPluginBase):
    package = 'rescript.tests'

    def test_missing_data_type(self):
        # Test when data_type is None
        with self.assertRaisesRegex(ValueError, "data-type"):
            parameter_validation()

    def test_rql_query_and_other_params(self):
        # Test when rql_query is specified with other conflicting parameters
        with self.assertRaisesRegex(ValueError,
                                    "rql_query.*can't.*simultaneously"):
            parameter_validation(rql_query="some_query", ids=[1, 2],
                                 data_field="genome_id", data_type="genome")

    def test_metadata_and_other_params(self):
        # Test when metadata is specified with other conflicting parameters
        with self.assertRaisesRegex(ValueError,
                                    "metadata.*can't.*simultaneously"):
            parameter_validation(metadata="metadata", ids=[1, 2],
                                 data_field="genome_id", data_type="genome")

    def test_ids_without_data_field(self):
        # Test when ids is specified without data_field
        with self.assertRaisesRegex(ValueError, r"ids.*data-field"):
            parameter_validation(ids=[1, 2], data_type="genome")

    def test_no_rql_query_ids_metadata(self):
        # Test when neither rql_query, ids, nor metadata is specified
        with self.assertRaisesRegex(ValueError, "rql-query.*ids.*metadata"):
            parameter_validation(data_type="genome")

    def test_invalid_data_field_for_data_type(self):
        # Test when the data_field is not valid for the given data_type
        with self.assertRaisesRegex(ValueError, "data-field.*permitted"):
            parameter_validation(ids=[1, 2], data_field="invalid_field",
                                 data_type="genome")

    def test_valid_rql_query_generation(self):
        rql_query = parameter_validation(
            ids=["Bacillus subtilis", "Bacteroidales bacterium"],
            data_field="species",
            data_type="genome"
        )
        self.assertEqual(rql_query, "in(species,(%22Bacillus%20subtilis%22,"
                                    "%22Bacteroidales%20bacterium%22))")

    def test_valid_rql_query_with_metadata(self):
        # Create mock metadata objects
        mock_metadata = MagicMock()
        mock_series = MagicMock()

        # Mock the .name attribute of the Series to return "species"
        mock_series.name = "species"

        # Mock the return value of to_series to be the mock_series
        mock_metadata.to_series.return_value = mock_series

        # Mock the ids in the series (mimicking a list of ids)
        mock_series.__iter__.return_value = iter(
            ["Bacillus subtilis", "Bacteroidales bacterium"])

        # Call the function with the mock metadata
        rql_query = parameter_validation(
            metadata=mock_metadata,
            data_type="genome"
        )

        # Assert that the rql_query is correctly generated
        self.assertEqual(rql_query, "in(species,(%22Bacillus%20subtilis%22,"
                                    "%22Bacteroidales%20bacterium%22))")


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

    @patch('rescript.bv_brc.requests.post')
    @patch('rescript.bv_brc.read_tsv_data_with_dtypes')
    @patch('rescript.bv_brc.error_handling')
    def test_download_data_json_batch(self, mock_error_handling,
                                      mock_read_data_with_dtypes, mock_post):
        # Mock the response for the first batch (25,000 entries)
        mock_response_1 = MagicMock()
        mock_response_1.status_code = 200
        mock_response_1.json.return_value = [{"id": i} for i in range(25000)]

        # Mock the response for the second batch (fewer than 25,000 entries)
        mock_response_2 = MagicMock()
        mock_response_2.status_code = 200
        mock_response_2.json.return_value = [{"id": i} for i in
                                             range(25000, 25010)]

        # The first call to requests.post returns 25,000 entries,
        # the second call returns 10 entries.
        mock_post.side_effect = [mock_response_1, mock_response_2]

        # Call the function for JSON
        result = download_data(data_type="genome", query="eq(id,1)",
                               accept="application/json", select=["genome_id"])

        # Check that the result is as expected (25,000 + 10 entries)
        self.assertEqual(len(result), 25010)
        self.assertEqual(result[0], {"id": 0})
        self.assertEqual(result[-1], {"id": 25009})

        # Ensure requests.post is called with the correct parameters
        mock_post.assert_any_call(
            url="https://www.bv-brc.org/api/genome/",
            data="eq(id,1)&limit(25000,0)&select(genome_id)",
            headers={
                'Content-Type': 'application/rqlquery+x-www-form-urlencoded',
                'ACCEPT': 'application/json'}
        )

        mock_post.assert_any_call(
            url="https://www.bv-brc.org/api/genome/",
            data="eq(id,1)&limit(25000,25000)&select(genome_id)",
            headers={
                'Content-Type': 'application/rqlquery+x-www-form-urlencoded',
                'ACCEPT': 'application/json'}
        )

    @patch('rescript.bv_brc.requests.post')
    @patch('rescript.bv_brc.read_tsv_data_with_dtypes')
    @patch('rescript.bv_brc.error_handling')
    def test_download_data_tsv(self, mock_error_handling,
                               mock_read_data_with_dtypes, mock_post):
        # Mock the response for TSV data type
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_post.return_value = mock_response

        # Mock reading TSV data
        mock_df = pd.DataFrame({"id": [1, 2]})
        mock_read_data_with_dtypes.return_value = mock_df

        # Call the function for TSV
        result = download_data(data_type="genome", query="eq(id,1)",
                               accept="text/tsv")

        # Check that the result is a DataFrame and as expected
        pd.testing.assert_frame_equal(result, pd.DataFrame({"id": [1, 2]}))

        # Ensure requests.post is called with the correct parameters
        mock_post.assert_called_with(
            url="https://www.bv-brc.org/api/genome/",
            data="eq(id,1)&limit(25000,0)",
            headers={
                'Content-Type': 'application/rqlquery+x-www-form-urlencoded',
                'ACCEPT': 'text/tsv'}
        )

        # Ensure read_tsv_data_with_dtypes was called
        mock_read_data_with_dtypes.assert_called_once_with(
            response=mock_response, data_type="genome")

    @patch('rescript.bv_brc.requests.post')
    @patch('rescript.bv_brc.error_handling')
    def test_download_data_gff(self, mock_error_handling, mock_post):
        # Mock the response for GFF data type
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "mock_gff_data"
        mock_post.return_value = mock_response

        # Call the function for GFF
        result = download_data(data_type="genome", query="eq(id,1)",
                               accept="application/gff")

        # Check that the result is as expected
        self.assertEqual(result, "mock_gff_data")

        # Ensure requests.post is called with the correct parameters
        mock_post.assert_called_with(
            url="https://www.bv-brc.org/api/genome/",
            data="eq(id,1)&limit(25000,0)",
            headers={
                'Content-Type': 'application/rqlquery+x-www-form-urlencoded',
                'ACCEPT': 'application/gff'}
        )

    @patch('rescript.bv_brc.requests.post')
    @patch('rescript.bv_brc.error_handling')
    def test_download_data_400_error(self, mock_error_handling, mock_post):
        # Mock the response for 400 error
        mock_response = MagicMock()
        mock_response.status_code = 400
        mock_post.return_value = mock_response

        # Call the function and check that error_handling is called
        download_data(data_type="genome", query="eq(id,1)",
                      accept="application/json")

        # Ensure error_handling was called for the 400 response
        mock_error_handling.assert_called_once_with(response=mock_response,
                                                    data_type="genome")

    @patch('rescript.bv_brc.requests.post')
    def test_download_data_unexpected_error(self, mock_post):
        # Mock the response for unexpected error
        mock_response = MagicMock()
        mock_response.status_code = 500
        mock_response.text = "Internal Server Error"
        mock_post.return_value = mock_response

        # Check that the function raises a ValueError for non-200, non-400
        # responses
        with self.assertRaises(ValueError) as context:
            download_data(data_type="genome", query="eq(id,1)",
                          accept="application/json")

        # Ensure the ValueError contains the right message
        self.assertEqual(str(context.exception), "Internal Server Error")


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
        result = create_genome_fasta(self.json_input_1)

        # Expected FASTA content
        expected_fasta = ">acc1   desc1   [genome_name1 | genome1]\nATGC"

        # Check if the file was created with the correct path and content
        mock_file.assert_called_once_with(f"{str(result)}/genome1.fasta", 'w')
        mock_file().write.assert_called_once_with(expected_fasta)

    @patch('rescript.bv_brc.open', new_callable=mock_open)
    def test_json_to_fasta_multiple_genomes(self, mock_file):
        create_genome_fasta(self.json_input_1 + self.json_input_2)

        # Expected FASTA content
        expected_fasta_genome1 = (">acc1   desc1   [genome_name1 | "
                                  "genome1]\nATGC")
        expected_fasta_genome2 = (">acc2   desc2   [genome_name2 | "
                                  "genome2]\nCGTA")

        # Check if the files were created with the correct path and content
        mock_file().write.assert_any_call(expected_fasta_genome1)
        mock_file().write.assert_any_call(expected_fasta_genome2)

    @patch('rescript.bv_brc.open', new_callable=mock_open)
    def test_json_to_fasta_multiple_sequences_same_genome(self, mock_file):
        result = create_genome_fasta(self.json_input_1 + self.json_input_1)

        # Expected FASTA content
        expected_fasta = (">acc1   desc1   [genome_name1 | genome1]\nATGC\n"
                          ">acc1   desc1   [genome_name1 | genome1]\nATGC")

        # Check if the file was created with the correct path and content
        mock_file.assert_called_once_with(f"{str(result)}/genome1.fasta", 'w')
        mock_file().write.assert_called_once_with(expected_fasta)


class TestGetBvBrcGenomes(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.get_taxonomy')
    @patch('rescript.bv_brc.get_genome_sequences')
    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.parameter_validation')
    def test_get_bv_brc_genomes(self, mock_parameter_validation,
                                mock_download_data, mock_get_genome_sequences,
                                mock_get_taxonomy):

        # Call the function
        get_bv_brc_genomes(
            ids_metadata=MagicMock(name='NumericMetadataColumn'),
            rql_query=None,
            data_field="mock_field",
            ids=["id1", "id2"],
            ranks=["rank1", "rank2"],
            rank_propagation=True
        )


class TestGetBvBrcMetadata(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.parameter_validation')
    def test_get_bv_brc_metadata(self, mock_parameter_validation,
                                 mock_download_data):
        # Mock the return value of parameter_validation
        mock_parameter_validation.return_value = 'rql(query)'

        # Mock the return value of download_data
        mock_df = pd.DataFrame({
            'id': ['id1', 'id2'],
            'feature': ['value1', 'value2'],
            'empty_field': [' ', None]
        }).set_index('id')
        mock_download_data.return_value = mock_df

        # Call the function
        result_metadata = get_bv_brc_metadata(
            ids_metadata=None,
            data_type='genome',
            data_field='genome_id',
            ids=['id1', 'id2']
        )

        # Assertions on the returned qiime2.Metadata
        self.assertIsInstance(result_metadata, qiime2.Metadata)

        # Extract the DataFrame from the result and check its values
        result_df = result_metadata.to_dataframe()

        # Check that the DataFrame's index and columns are correct
        self.assertEqual(result_df.index.tolist(), ['id1', 'id2'])
        self.assertIn('feature', result_df.columns)

        # Check that empty/space-only fields are replaced with NaN
        self.assertTrue(pd.isna(result_df.loc['id1', 'empty_field']))
        self.assertTrue(pd.isna(result_df.loc['id2', 'empty_field']))


class TestGetBvBrcGenomeFeatures(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.get_loci')
    @patch('rescript.bv_brc.get_taxonomy')
    @patch('rescript.bv_brc.get_sequences')
    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.parameter_validation')
    def test_get_bv_brc_genome_features(self, mock_parameter_validation,
                                        mock_download_data, mock_get_sequences,
                                        mock_get_taxonomy, mock_get_loci):
        # Mocked return values for the external functions
        mock_parameter_validation.return_value = "mocked_rql_query"
        mock_download_data.return_value = [{"genome_id": "genome1"}]
        mock_get_sequences.return_value = (MagicMock(), MagicMock())
        mock_get_taxonomy.return_value = MagicMock()
        mock_get_loci.return_value = MagicMock()

        # Call the function
        get_bv_brc_genome_features(
            ids_metadata=None,
            rql_query=None,
            data_field="mock_field",
            ids=["id1", "id2"],
            ranks=["rank1", "rank2"],
            rank_propagation=True
        )


class TestCreateTaxonomyEntry(TestPluginBase):
    package = 'rescript.tests'

    def test_basic_functionality(self):
        lineage_names = ("Animalia;Chordata;Mammalia;Primates;Hominidae;Homo;"
                         "sapiens")
        lineage_ranks = "kingdom;phylum;class;order;family;genus;species"
        expected_result = ("k__Animalia; p__Chordata; c__Mammalia; "
                           "o__Primates; f__Hominidae; g__Homo; s__sapiens")

        result = create_taxonomy_entry(lineage_names, lineage_ranks)
        self.assertEqual(result, expected_result)

    def test_with_missing_ranks(self):
        lineage_names = "Animalia;Chordata;Hominidae;sapiens"
        lineage_ranks = "kingdom;phylum;family;species"
        expected_result = ("k__Animalia; p__Chordata; c__; o__; f__Hominidae; "
                           "g__; s__sapiens")

        result = create_taxonomy_entry(lineage_names, lineage_ranks,
                                       rank_propagation=False)
        self.assertEqual(result, expected_result)

    def test_rank_propagation(self):
        lineage_names = "Animalia;Chordata;Mammalia;Homo"
        lineage_ranks = "kingdom;phylum;class;genus"
        expected_result = ("k__Animalia; p__Chordata; c__Mammalia; "
                           "o__Mammalia; f__Mammalia; g__Homo; s__Homo")

        result = create_taxonomy_entry(lineage_names, lineage_ranks,
                                       rank_propagation=True)
        self.assertEqual(result, expected_result)

    def test_genus_species_split(self):
        lineage_names = ("Animalia;Chordata;Mammalia;Primates;"
                         "Hominidae;Homo sapiens")
        lineage_ranks = "kingdom;phylum;class;order;family;species"
        expected_result = ("k__Animalia; p__Chordata; c__Mammalia;"
                           " o__Primates; f__Hominidae; g__Homo; s__sapiens")

        result = create_taxonomy_entry(lineage_names, lineage_ranks)
        self.assertEqual(result, expected_result)

    def test_genus_only_split(self):
        lineage_names = ("Animalia;Chordata;Mammalia;Primates;Hominidae;"
                         "Homo;Homo sapiens")
        lineage_ranks = "kingdom;phylum;class;order;family;genus;species"
        expected_result = ("k__Animalia; p__Chordata; c__Mammalia; "
                           "o__Primates; f__Hominidae; g__Homo; s__sapiens")

        result = create_taxonomy_entry(lineage_names, lineage_ranks)
        self.assertEqual(result, expected_result)

    def test_no_species_in_ranks(self):
        lineage_names = ("Animalia;Chordata;Mammalia;Primates;"
                         "Hominidae;Homo sapiens")
        lineage_ranks = "kingdom;phylum;class;order;family;genus"
        expected_result = ("k__Animalia; p__Chordata; c__Mammalia; "
                           "o__Primates; f__Hominidae; g__Homo sapiens")

        result = create_taxonomy_entry(lineage_names, lineage_ranks,
                                       ranks=['kingdom', 'phylum',
                                              'class', 'order',
                                              'family', 'genus'])
        self.assertEqual(result, expected_result)

    def test_custom_ranks(self):
        lineage_names = "Metazoa"
        lineage_ranks = "superkingdom"
        expected_result = "sk__Metazoa"

        result = create_taxonomy_entry(lineage_names, lineage_ranks,
                                       ranks=['superkingdom'])
        self.assertEqual(result, expected_result)


class TestGetTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.create_taxonomy')
    @patch('rescript.bv_brc.download_data')
    def test_get_taxonomy(self, mock_download_data, mock_create_taxonomy):
        # Mock response_sequences (list of dicts)
        response_sequences = [
            {"taxon_id": "taxon1", "feature_id": "feature1"},
            {"taxon_id": "taxon2", "feature_id": "feature2"},
            {"taxon_id": "taxon3", "feature_id": "feature3"},
            {"taxon_id": "taxon1", "feature_id": "feature4"},
        ]

        # Mock the download_data and create_taxonomy functions
        mock_download_data.return_value = MagicMock()
        mock_create_taxonomy.return_value = MagicMock()

        # Define test parameters
        ranks = ["kingdom", "phylum", "class"]
        rank_propagation = True
        accession_name = "feature_id"

        # Call the function
        get_taxonomy(
            response_sequences=response_sequences,
            ranks=ranks,
            rank_propagation=rank_propagation,
            accession_name=accession_name
        )

        # Check that the taxon_ids are extracted correctly and that duplicates
        # are removed
        expected_taxon_ids = {"taxon1", "taxon2", "taxon3"}
        extracted_taxon_ids = {str(entry["taxon_id"]) for entry in
                               response_sequences}
        self.assertEqual(extracted_taxon_ids, expected_taxon_ids)


class TestCreateTaxonomy(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.TSVTaxonomyFormat')
    @patch('rescript.bv_brc.create_taxonomy_entry')
    def test_create_taxonomy(self, mock_create_taxonomy_entry,
                             mock_TSVTaxonomyFormat):
        # Mock the input data for taxonomy_bvbrc DataFrame
        taxonomy_bvbrc = pd.DataFrame({
            'taxon_id': ['taxon1', 'taxon2'],
            'lineage_names': [['Bacteria', 'Proteobacteria'],
                              ['Bacteria', 'Firmicutes']],
            'lineage_ranks': [['domain', 'phylum'], ['domain', 'phylum']]
        })

        # Mock the create_taxonomy_entry to return a fake taxonomy string
        mock_create_taxonomy_entry.side_effect = \
            lambda lineage_names, lineage_ranks, rank_propagation, ranks: (
                ";".join(lineage_names))

        # Mock response_sequences (list of dicts)
        response_sequences = [
            {"taxon_id": "taxon1", "feature_id": "feature1"},
            {"taxon_id": "taxon2", "feature_id": "feature2"}
        ]

        # Mock TSVTaxonomyFormat to return a file-like object
        mock_taxonomy_file = MagicMock()
        (mock_TSVTaxonomyFormat.return_value.open.
         return_value.__enter__).return_value = mock_taxonomy_file

        # Call the create_taxonomy function
        create_taxonomy(
            taxonomy_bvbrc=taxonomy_bvbrc,
            response_sequences=response_sequences,
            ranks=['domain', 'phylum'],
            rank_propagation=True,
            accession_name="feature_id"
        )

        # Ensure that the correct data was written to the file
        written_data = "\t".join(["Feature ID", "Taxon"]) + "\n" + \
                       "\t".join(
                           ["feature1", "Bacteria;Proteobacteria"]) + "\n" + \
                       "\t".join(["feature2", "Bacteria;Firmicutes"]) + "\n"

        # Check if the 'write' method was actually called with the right data
        mock_taxonomy_file.write.assert_called_once_with(written_data)


class TestReadDataWithDtypes(TestPluginBase):
    package = 'rescript.tests'

    def test_read_data_with_dtypes(self):
        # Mock response with a TSV file containing antibiotics data
        mock_response = MagicMock()
        mock_response.text = """_version_\tantibiotic_name\tcas_id\tdescription
1.0\tPenicillin\tCAS1234\tAntibiotic description
2.0\tAmoxicillin\tCAS5678\tAntibiotic description 2
"""

        # Call the function with mock data and data_type "antibiotics"
        df = read_tsv_data_with_dtypes(mock_response, "antibiotics")

        # Expected DataFrame output
        expected_data = {
            "_version_": ["1.0", "2.0"],
            "antibiotic_name": ["Penicillin", "Amoxicillin"],
            "cas_id": ["CAS1234", "CAS5678"],
            "description": ["Antibiotic description",
                            "Antibiotic description 2"]
        }
        expected_df = pd.DataFrame(expected_data)

        # Check if the DataFrame matches the expected output
        pd.testing.assert_frame_equal(df, expected_df)

    def test_no_data_raises_value_error(self):
        # Mock response with only the column headers
        mock_response = MagicMock()
        mock_response.text = """_version_\tantibiotic_name\tcas_id"""

        # Assert that a ValueError is raised when no data rows are present
        with self.assertRaises(ValueError):
            read_tsv_data_with_dtypes(mock_response, "antibiotics")


class TestGetLoci(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.process_loci')
    @patch('builtins.open', new_callable=mock_open)
    @patch('rescript.bv_brc.LociDirectoryFormat')
    def test_get_loci(self, mock_loci_dir_format, mock_open_file,
                      mock_process_loci, mock_download_data):
        # Set up the mocks
        mock_download_data.return_value = "mock_gff_data"
        mock_process_loci.return_value = "processed_gff_data"

        # Create a mock directory format
        mock_loci_dir = MagicMock()
        mock_loci_dir_format.return_value = mock_loci_dir

        # Mock response_sequences with genome_id information
        response_sequences = [
            {"genome_id": "genome1"},
            {"genome_id": "genome2"},
            {"genome_id": "genome1"}
        ]

        # Call the function with response_sequences
        get_loci(response_sequences)

        # Check that download_data is called for each unique genome_id
        mock_download_data.assert_any_call(data_type="genome_feature",
                                           query="eq(genome_id,genome1)",
                                           accept="application/gff")
        mock_download_data.assert_any_call(data_type="genome_feature",
                                           query="eq(genome_id,genome2)",
                                           accept="application/gff")

        # Ensure it is only called twice (for unique genome_ids)
        self.assertEqual(mock_download_data.call_count, 2)

        # Check that process_loci is called with the right data
        mock_process_loci.assert_any_call(gff_string="mock_gff_data")

        # Check that open was called with the correct paths and the data was
        # written
        mock_open_file.assert_any_call(
            mock_loci_dir.__str__() + '/genome1.gff', 'w')
        mock_open_file.assert_any_call(
            mock_loci_dir.__str__() + '/genome2.gff', 'w')

        # Ensure the processed data was written to the file
        mock_open_file().write.assert_any_call("processed_gff_data")

    def test_process_loci(self):
        # Input GFF string with both headers and data lines
        input_data = """##gff-version 3
##sequence-region accn|NC_000001.11 1 1000
accn|NC_000001.11\tRefSeq\tregion\t1\t1000\t.\t+\t.\tID=region0;
accn|NC_000002.11\tRefSeq\tgene\t1\t1000\t.\t+\t.\tID=gene0;"""

        # Expected output after processing (removing "accn|")
        expected_output = """##gff-version 3
##sequence-region accn|NC_000001.11 1 1000
NC_000001.11\tRefSeq\tregion\t1\t1000\t.\t+\t.\tID=region0;
NC_000002.11\tRefSeq\tgene\t1\t1000\t.\t+\t.\tID=gene0;"""

        # Call the function with the input data
        result = process_loci(input_data)

        # Check if the result matches the expected output
        self.assertEqual(result, expected_output)


class TestGetSequences(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.download_data')
    @patch('builtins.open', new_callable=mock_open)
    def test_get_sequences(self, mock_open, mock_download_data):
        # Mock data for genome_features
        genome_features = [
            {
                "genome_id": "genome1",
                "feature_id": "feature1",
                "na_sequence_md5": "md5na1",
                "aa_sequence_md5": "md5aa1"
            },
            {
                "genome_id": "genome2",
                "feature_id": "feature2",
                "na_sequence_md5": "md5na2"
            }
        ]

        # Mock the return value of download_data
        mock_download_data.return_value = [
            {"md5": "md5na1", "sequence": "atgc"},
            {"md5": "md5aa1", "sequence": "MKV"},
            {"md5": "md5na2", "sequence": "gtca"}
        ]

        # Call the function
        get_sequences(genome_features)

        # Assertions to check that download_data was called correctly
        mock_download_data.assert_called_once_with(
            data_type="feature_sequence",
            query=unittest.mock.ANY,
            accept="application/json",
            select=["md5", "sequence"]
        )

        # Check that genome1.fasta and genome2.fasta are in the paths used
        # for file opening
        open_calls = [call[0][0] for call in mock_open.call_args_list]

        self.assertTrue(any('genome1.fasta' in call for call in open_calls))
        self.assertTrue(any('genome2.fasta' in call for call in open_calls))

        # Check if the correct sequences were written to the genes file
        mock_open().write.assert_any_call('>feature1\nATGC\n')
        mock_open().write.assert_any_call('>feature2\nGTCA\n')

        # Check if the correct sequences were written to the proteins file
        mock_open().write.assert_any_call('>feature1\nMKV\n')


class TestGetGenomeSequences(TestPluginBase):
    package = 'rescript.tests'

    @patch('rescript.bv_brc.download_data')
    @patch('rescript.bv_brc.create_genome_fasta')
    def test_get_genome_sequences(self, mock_create_genome_fasta,
                                  mock_download_data):
        # Sample response_genomes to be used as input
        response_genomes = [
            {'genome_id': '12345'},
            {'genome_id': '67890'}
        ]

        # Call the function
        get_genome_sequences(response_genomes)

        # Assert that download_data was called with the correct arguments
        mock_download_data.assert_called_once_with(
            data_type="genome_sequence",
            query=unittest.mock.ANY,
            accept="application/json",
            select=["accession", "description", "genome_name", "genome_id",
                    "sequence"]
        )
