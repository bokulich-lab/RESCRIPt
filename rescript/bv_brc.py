# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from io import StringIO
import os
import qiime2
import pandas as pd
import requests
from q2_types.feature_data import MixedCaseDNAFASTAFormat, ProteinFASTAFormat, TSVTaxonomyDirectoryFormat
from q2_types.genome_data import GenomeSequencesDirectoryFormat

from rescript.ncbi import _allowed_ranks, _default_ranks
import json


def fetch_genomes_bv_brc(
        rql_query: str = None,
        genome_ids: list = None
) -> GenomeSequencesDirectoryFormat:

    # Parameter validation
    rql_query = id_list_handling(rql_query=rql_query,
                                 ids=genome_ids,
                                 parameter_name="genome_ids",
                                 data_field="genome_id"
                                 )

    # Define output format
    genomes = GenomeSequencesDirectoryFormat()

    # Get requests response
    response = download_data(
        url=f"https://www.bv-brc.org/api/genome_sequence/?{rql_query}",
        data_type="genome_sequence",
    )

    # Transform
    json_to_fasta(response.json(), str(genomes))

    return genomes


def fetch_metadata_bv_brc(data_type: str, rql_query: str) -> qiime2.Metadata:

    # Get requests response
    response = download_data(
        url=f"https://www.bv-brc.org/api/{data_type}/?{rql_query}&http_accept=text/tsv",
        data_type=data_type
    )

    tsv_data = StringIO(response.text)
    metadata = pd.read_csv(tsv_data, sep='\t')
    metadata.index.name = "id"
    metadata.index = metadata.index.astype(str)

    return qiime2.Metadata(metadata)


def fetch_taxonomy_bv_brc(
        rql_query: str,
        ranks: list = None,
        taxon_ids: list = None,
) -> TSVTaxonomyDirectoryFormat:

    # Parameter validation
    rql_query = id_list_handling(rql_query=rql_query,
                                 ids=taxon_ids,
                                 parameter_name="taxon_ids",
                                 data_field="taxon_id"
                                 )

    # Define output format
    directory = TSVTaxonomyDirectoryFormat()

    # Get requests response
    response = download_data(
        url=f"https://www.bv-brc.org/api/taxonomy/?{rql_query}&http_accept=text/tsv",
        data_type="taxonomy"
    )

    # Convert to data frame
    tsv_data = StringIO(response.text)
    metadata = pd.read_csv(tsv_data, sep='\t')

    # Transform metadata to TSVTaxonomyFormat
    taxonomy = transform_taxonomy_df(df=metadata, ranks=ranks)
    taxonomy.to_csv(os.path.join(str(directory), "taxonomy.tsv"), sep="\t")

    return directory


def parse_lineage_names_with_ranks(lineage_names, lineage_ranks, ranks):
    # Set ranks to default if no list is specified
    if not ranks:
        ranks = _default_ranks

    # Split the lineage names and ranks by ';'
    lineage_split = lineage_names.split(';')
    rank_split = lineage_ranks.split(';')

    # Dictionary to map taxonomic ranks to their prefixes for the specified ranks
    rank_to_prefix = {key: _allowed_ranks[key] for key in ranks if key in ranks}

    # Initialize the list for the parsed lineage
    parsed_lineage = []

    # Loop over each rank and assign the corresponding prefix and name
    for rank, name in zip(rank_split, lineage_split):
        prefix = rank_to_prefix.get(rank, None)
        if prefix:
            parsed_lineage.append(f"{prefix}{name}")
        else:
            pass

    # Ensure all taxonomic levels are covered (fill in missing levels with just the
    # prefix)
    final_lineage = []
    for required_prefix in rank_to_prefix.values():
        # Check if any parsed_lineage item starts with the required prefix
        match = next(
            (item for item in parsed_lineage if item.startswith(required_prefix)), None)
        if match:
            final_lineage.append(match)
        else:
            final_lineage.append(required_prefix)

    # Join the parsed lineage names with '; '
    return '; '.join(final_lineage)


def transform_taxonomy_df(df, ranks):
    # Apply the transformation
    df['Taxon'] = df.apply(
        lambda row: parse_lineage_names_with_ranks(lineage_names=row['lineage_names'],
                                                   lineage_ranks=row['lineage_ranks'],
                                                   ranks=ranks), axis=1)

    # Rename columns and set index
    df = df.rename(columns={'taxon_id': 'Feature ID'})
    df = df[['Feature ID', 'Taxon']]
    df = df.set_index('Feature ID')
    return df


def fetch_genome_features_bv_brc(
        rql_query: str = None,
        feature_ids: list = None,
) -> (MixedCaseDNAFASTAFormat, ProteinFASTAFormat):

    # Parameter validation
    rql_query = id_list_handling(rql_query=rql_query,
                                 ids=feature_ids,
                                 parameter_name="feature_ids",
                                 data_field="feature_id")

    # Define output formats
    genes = MixedCaseDNAFASTAFormat()
    proteins = ProteinFASTAFormat()

    # Construct URLs for genes and proteins downloads
    base_url = "https://www.bv-brc.org/api/genome_feature/?"
    genes_url = base_url + f"{rql_query}&http_accept=application/dna+fasta"
    proteins_url = base_url + f"{rql_query}&http_accept=application/protein+fasta"

    # Get requests response for genes and proteins
    response_genes = download_data(url=genes_url, data_type="genome_feature")
    response_proteins = download_data(url=proteins_url, data_type="genome_feature")

    # Save genes and proteins as FASTA files
    fasta_genes = response_genes.text
    with genes.open() as file:
        file.write(fasta_genes)

    fasta_proteins = response_proteins.text
    with proteins.open() as file:
        file.write(fasta_proteins)

    return genes, proteins


def json_to_fasta(json, output_dir):
    # Dictionary to hold sequences grouped by genome_id
    fasta_files = {}

    # Loop over all entries in dict
    for entry in json:
        genome_id = entry['genome_id']
        if genome_id not in fasta_files:
            fasta_files[genome_id] = []

        # Construct FASTA format to be identical to BV-BRC FASTA headers
        header = (f">accn|{entry['accession']}   {entry['description']}   "
                  f"[{entry['genome_name']} | {genome_id}]")

        fasta_files[genome_id].append(f"{header}\n{entry['sequence'].upper()}")

    # Write each genome_id's sequences to a separate FASTA file
    for genome_id, sequences in fasta_files.items():
        fasta_content = "\n".join(sequences)
        fasta_filename = os.path.join(output_dir, f"{genome_id}.fasta")

        with open(fasta_filename, 'w') as fasta_file:
            fasta_file.write(fasta_content)


def download_data(url, data_type):
    # Get requests response
    response = requests.get(url)

    # If response is correct return it
    if response.status_code == 200:
        return response

    # Error handling if response incorrect
    elif response.status_code == 400:
        error_handling(response, data_type)
    else:
        raise ValueError(response.text)


def error_handling(response, data_type):
    # No data found for query or incorrect RQL query
    if response.text == "[]":
        raise ValueError("No data could be retrieved. Either because of an "
                         "incorrect RQL query or because no data exists for the "
                         "query.")

    elif response.text.startswith("A Database Error Occured:"):

        # Parse the response dict
        json_str = response.text[response.text.find('{'):]
        response_dict = json.loads(json_str)

        # Incorrect RQL operator
        if response_dict['msg'].startswith("undefined field object"):
            raise ValueError(
                f"Error code {response_dict['code']}: {response_dict['msg']}. "
                f"Incorrect RQL query operator."
            )

        # Incorrect field for data type
        elif response_dict['msg'].startswith("undefined field"):
            raise ValueError(
                f"Error code {response_dict['code']}: {response_dict['msg']}. \n"
                f"Allowed fields for data type {data_type}: \n{data_fields[data_type]}"
            )

        else:
            raise ValueError(
                f"Error code {response_dict['code']}: {response_dict['msg']}."
            )

    else:
        raise ValueError(response.text)


def id_list_handling(rql_query: str, ids: list, parameter_name: str, data_field: str):
    # Error if rql_query and ids parameters are given
    if rql_query and ids:
        raise ValueError(f"Parameters rql_query and {parameter_name} can't be used "
                         "simultaneously.")

    # Error if rql_query and ids parameters are not given
    elif not rql_query and not ids:
        raise ValueError("At least one of the parameters rql_query and "
                         f"{parameter_name} has to be given.")

    # construct the RQL queries
    elif ids:
        rql_query = f"in({data_field},({','.join(map(str, ids))}))"

    return rql_query


data_fields = {
    "antibiotics": [
        "_version_",
        "antibiotic_name",
        "atc_classification",
        "canonical_smiles",
        "cas_id",
        "date_inserted",
        "date_modified",
        "description",
        "drugbank_interactions",
        "inchi_key",
        "isomeric_smiles",
        "mechanism_of_action",
        "molecular_formula",
        "molecular_weight",
        "pharmacological_classes",
        "pharmacology",
        "pubchem_cid",
        "pubchem_cid_i",
        "synonyms"
    ],
    "enzyme_class_ref": [
        "_version_",
        "date_inserted",
        "date_modified",
        "ec_description",
        "ec_number",
        "go"
    ],
    "epitope": [
        "_version_",
        "assay_results",
        "bcell_assays",
        "comments",
        "date_inserted",
        "date_modified",
        "end",
        "epitope_id",
        "epitope_sequence",
        "epitope_type",
        "host_name",
        "mhc_assays",
        "organism",
        "protein_accession",
        "protein_id",
        "protein_name",
        "start",
        "taxon_id",
        "taxon_lineage_ids",
        "taxon_lineage_names",
        "tcell_assays",
        "total_assays"
    ],
    "epitope_assay": [
        "_version_",
        "assay_group",
        "assay_id",
        "assay_measurement",
        "assay_measurement_unit",
        "assay_method",
        "assay_result",
        "assay_type",
        "authors",
        "date_inserted",
        "date_modified",
        "end",
        "epitope_id",
        "epitope_sequence",
        "epitope_type",
        "host_name",
        "host_taxon_id",
        "mhc_allele",
        "mhc_allele_class",
        "organism",
        "pdb_id",
        "pmid",
        "protein_accession",
        "protein_id",
        "protein_name",
        "start",
        "taxon_id",
        "taxon_lineage_ids",
        "taxon_lineage_names",
        "title"
    ],
    "experiment": [
        "_version_",
        "additional_data",
        "additional_metadata",
        "biosets",
        "date_inserted",
        "date_modified",
        "detection_instrument",
        "doi",
        "exp_description",
        "exp_id",
        "exp_name",
        "exp_poc",
        "exp_protocol",
        "exp_title",
        "exp_type",
        "experimenters",
        "genome_id",
        "measurement_technique",
        "organism",
        "pmid",
        "public_identifier",
        "public_repository",
        "samples",
        "strain",
        "study_description",
        "study_institution",
        "study_name",
        "study_pi",
        "study_title",
        "taxon_id",
        "taxon_lineage_ids",
        "treatment_amount",
        "treatment_duration",
        "treatment_name",
        "treatment_type"
    ],
    "bioset": [
        "_version_",
        "additional_data",
        "additional_metadata",
        "analysis_group_1",
        "analysis_group_2",
        "analysis_method",
        "bioset_criteria",
        "bioset_description",
        "bioset_id",
        "bioset_name",
        "bioset_result",
        "bioset_type",
        "date_inserted",
        "date_modified",
        "entity_count",
        "entity_type",
        "exp_id",
        "exp_name",
        "exp_title",
        "exp_type",
        "genome_id",
        "organism",
        "protocol",
        "result_type",
        "strain",
        "study_description",
        "study_institution",
        "study_name",
        "study_pi",
        "study_title",
        "taxon_id",
        "taxon_lineage_ids",
        "treatment_amount",
        "treatment_duration",
        "treatment_name",
        "treatment_type"
    ],
    "bioset_result": [
        "_version_",
        "bioset_description",
        "bioset_id",
        "bioset_name",
        "bioset_type",
        "counts",
        "date_inserted",
        "date_modified",
        "entity_id",
        "entity_name",
        "entity_type",
        "exp_id",
        "exp_name",
        "exp_title",
        "exp_type",
        "feature_id",
        "fpkm",
        "gene",
        "gene_id",
        "genome_id",
        "id",
        "locus_tag",
        "log2_fc",
        "organism",
        "other_ids",
        "other_value",
        "p_value",
        "patric_id",
        "product",
        "protein_id",
        "result_type",
        "strain",
        "taxon_id",
        "tpm",
        "treatment_amount",
        "treatment_duration",
        "treatment_name",
        "treatment_type",
        "uniprot_id",
        "z_score"
    ],
    "gene_ontology_ref": [
        "_version_",
        "date_inserted",
        "date_modified",
        "definition",
        "go_id",
        "go_name",
        "ontology"
    ],
    "genome": [
        "_version_",
        "additional_metadata",
        "altitude",
        "antimicrobial_resistance",
        "antimicrobial_resistance_evidence",
        "assembly_accession",
        "assembly_method",
        "authors",
        "bioproject_accession",
        "biosample_accession",
        "biovar",
        "body_sample_site",
        "body_sample_subsite",
        "cds",
        "cds_ratio",
        "cell_shape",
        "checkm_completeness",
        "checkm_contamination",
        "chromosomes",
        "clade",
        "class",
        "coarse_consistency",
        "collection_date",
        "collection_year",
        "comments",
        "common_name",
        "completion_date",
        "contig_l50",
        "contig_n50",
        "contigs",
        "core_families",
        "core_family_ratio",
        "culture_collection",
        "date_inserted",
        "date_modified",
        "depth",
        "disease",
        "family",
        "fine_consistency",
        "gc_content",
        "genbank_accessions",
        "genome_id",
        "genome_length",
        "genome_name",
        "genome_quality",
        "genome_quality_flags",
        "genome_status",
        "genus",
        "geographic_group",
        "geographic_location",
        "gram_stain",
        "h1_clade_global",
        "h1_clade_us",
        "h3_clade",
        "h5_clade",
        "h_type",
        "habitat",
        "host_age",
        "host_common_name",
        "host_gender",
        "host_group",
        "host_health",
        "host_name",
        "host_scientific_name",
        "hypothetical_cds",
        "hypothetical_cds_ratio",
        "isolation_comments",
        "isolation_country",
        "isolation_site",
        "isolation_source",
        "kingdom",
        "lab_host",
        "latitude",
        "lineage",
        "longitude",
        "mat_peptide",
        "missing_core_family_ids",
        "mlst",
        "motility",
        "n_type",
        "ncbi_project_id",
        "nearest_genomes",
        "optimal_temperature",
        "order",
        "organism_name",
        "other_clinical",
        "other_environmental",
        "other_names",
        "other_typing",
        "outgroup_genomes",
        "owner",
        "oxygen_requirement",
        "p2_genome_id",
        "partial_cds",
        "partial_cds_ratio",
        "passage",
        "pathovar",
        "patric_cds",
        "ph1n1_like",
        "phenotype",
        "phylum",
        "plasmids",
        "plfam_cds",
        "plfam_cds_ratio",
        "public",
        "publication",
        "reference_genome",
        "refseq_accessions",
        "refseq_cds",
        "refseq_project_id",
        "rrna",
        "salinity",
        "season",
        "segment",
        "segments",
        "sequencing_centers",
        "sequencing_depth",
        "sequencing_platform",
        "sequencing_status",
        "serovar",
        "species",
        "sporulation",
        "sra_accession",
        "state_province",
        "strain",
        "subclade",
        "subtype",
        "superkingdom",
        "taxon_id",
        "taxon_lineage_ids",
        "taxon_lineage_names",
        "temperature_range",
        "trna",
        "type_strain",
        "user_read",
        "user_write"
    ],
    "strain": [
        "1_pb2",
        "2_pb1",
        "3_pa",
        "4_ha",
        "5_np",
        "6_na",
        "7_mp",
        "8_ns",
        "_version_",
        "collection_date",
        "collection_year",
        "date_inserted",
        "date_modified",
        "family",
        "genbank_accessions",
        "genome_ids",
        "genus",
        "geographic_group",
        "h_type",
        "host_common_name",
        "host_group",
        "host_name",
        "id",
        "isolation_country",
        "l",
        "lab_host",
        "m",
        "n_type",
        "other_segments",
        "owner",
        "passage",
        "public",
        "s",
        "season",
        "segment_count",
        "species",
        "status",
        "strain",
        "subtype",
        "taxon_id",
        "taxon_lineage_ids",
        "taxon_lineage_names",
        "user_read",
        "user_write"
    ],
    "genome_amr": [
        "_version_",
        "antibiotic",
        "computational_method",
        "computational_method_performance",
        "computational_method_version",
        "date_inserted",
        "date_modified",
        "evidence",
        "genome_id",
        "genome_name",
        "id",
        "laboratory_typing_method",
        "laboratory_typing_method_version",
        "laboratory_typing_platform",
        "measurement",
        "measurement_sign",
        "measurement_unit",
        "measurement_value",
        "owner",
        "pmid",
        "public",
        "resistant_phenotype",
        "source",
        "taxon_id",
        "testing_standard",
        "testing_standard_year",
        "user_read",
        "user_write",
        "vendor"
    ],
    "feature_sequence": [
        "_version_",
        "date_inserted",
        "date_modified",
        "md5",
        "sequence",
        "sequence_type"
    ],
    "genome_feature": [
        "aa_length",
        "aa_sequence_md5",
        "accession",
        "alt_locus_tag",
        "annotation",
        "brc_id",
        "classifier_round",
        "classifier_score",
        "codon_start",
        "date_inserted",
        "date_modified",
        "end",
        "feature_id",
        "feature_type",
        "figfam_id",
        "gene",
        "gene_id",
        "genome_id",
        "genome_name",
        "go",
        "location",
        "na_length",
        "na_sequence_md5",
        "notes",
        "og_id",
        "owner",
        "p2_feature_id",
        "patric_id",
        "pdb_accession",
        "pgfam_id",
        "plfam_id",
        "product",
        "property",
        "protein_id",
        "public",
        "refseq_locus_tag",
        "segments",
        "sequence_id",
        "sog_id",
        "start",
        "strand",
        "taxon_id",
        "uniprotkb_accession",
        "user_read",
        "user_write"
    ],
    "genome_sequence": [
        "_version_",
        "accession",
        "chromosome",
        "date_inserted",
        "date_modified",
        "description",
        "gc_content",
        "genome_id",
        "genome_name",
        "gi",
        "length",
        "mol_type",
        "owner",
        "p2_sequence_id",
        "plasmid",
        "public",
        "release_date",
        "segment",
        "sequence",
        "sequence_id",
        "sequence_md5",
        "sequence_status",
        "sequence_type",
        "taxon_id",
        "topology",
        "user_read",
        "user_write",
        "version"
    ],
    "id_ref": [
        "_version_",
        "date_inserted",
        "date_modified",
        "id",
        "id_type",
        "id_value",
        "uniprotkb_accession"
    ],
    "misc_niaid_sgc": [
        "_version_",
        "date_inserted",
        "date_modified",
        "gene_symbol_collection",
        "genus",
        "has_clones",
        "has_proteins",
        "selection_criteria",
        "species",
        "strain",
        "target_id",
        "target_status"
    ],
    "pathway": [
        "_version_",
        "accession",
        "alt_locus_tag",
        "annotation",
        "date_inserted",
        "date_modified",
        "ec_description",
        "ec_number",
        "feature_id",
        "gene",
        "genome_ec",
        "genome_id",
        "genome_name",
        "id",
        "owner",
        "pathway_class",
        "pathway_ec",
        "pathway_id",
        "pathway_name",
        "patric_id",
        "product",
        "public",
        "refseq_locus_tag",
        "sequence_id",
        "taxon_id",
        "user_read",
        "user_write"
    ],
    "pathway_ref": [
        "_version_",
        "date_inserted",
        "date_modified",
        "ec_description",
        "ec_number",
        "id",
        "map_location",
        "map_name",
        "map_type",
        "occurrence",
        "pathway_class",
        "pathway_id",
        "pathway_name"
    ],
    "ppi": [
        "_version_",
        "category",
        "date_inserted",
        "date_modified",
        "detection_method",
        "domain_a",
        "domain_b",
        "evidence",
        "feature_id_a",
        "feature_id_b",
        "gene_a",
        "gene_b",
        "genome_id_a",
        "genome_id_b",
        "genome_name_a",
        "genome_name_b",
        "id",
        "interaction_type",
        "interactor_a",
        "interactor_b",
        "interactor_desc_a",
        "interactor_desc_b",
        "interactor_type_a",
        "interactor_type_b",
        "pmid",
        "refseq_locus_tag_a",
        "refseq_locus_tag_b",
        "score",
        "source_db",
        "source_id",
        "taxon_id_a",
        "taxon_id_b"
    ],
    "protein_family_ref": [
        "_version_",
        "date_inserted",
        "date_modified",
        "family_id",
        "family_product",
        "family_type"
    ],
    "sequence_feature": [
        "aa_sequence_md5",
        "aa_variant",
        "additional_metadata",
        "comments",
        "date_inserted",
        "date_modified",
        "end",
        "evidence_code",
        "feature_id",
        "genbank_accession",
        "gene",
        "genome_id",
        "genome_name",
        "id",
        "length",
        "patric_id",
        "product",
        "publication",
        "refseq_locus_tag",
        "segment",
        "segments",
        "sf_category",
        "sf_id",
        "sf_name",
        "sf_sequence",
        "sf_sequence_md5",
        "source",
        "source_aa_sequence",
        "source_id",
        "source_sf_location",
        "source_strain",
        "start",
        "subtype",
        "taxon_id",
        "variant_types"
    ],
    "sequence_feature_vt": [
        "additional_metadata",
        "comments",
        "date_inserted",
        "date_modified",
        "id",
        "sf_category",
        "sf_id",
        "sf_name",
        "sf_sequence",
        "sf_sequence_md5",
        "sfvt_genome_count",
        "sfvt_genome_ids",
        "sfvt_id",
        "sfvt_sequence",
        "sfvt_sequence_md5",
        "sfvt_variations"
    ],
    "sp_gene": [
        "_version_",
        "alt_locus_tag",
        "antibiotics",
        "antibiotics_class",
        "classification",
        "date_inserted",
        "date_modified",
        "e_value",
        "evidence",
        "feature_id",
        "function",
        "gene",
        "genome_id",
        "genome_name",
        "id",
        "identity",
        "organism",
        "owner",
        "patric_id",
        "pmid",
        "product",
        "property",
        "property_source",
        "public",
        "query_coverage",
        "refseq_locus_tag",
        "same_genome",
        "same_genus",
        "same_species",
        "source",
        "source_id",
        "subject_coverage",
        "taxon_id",
        "user_read",
        "user_write"
    ],
    "sp_gene_ref": [
        "_version_",
        "antibiotics",
        "antibiotics_class",
        "assertion",
        "classification",
        "date_inserted",
        "date_modified",
        "function",
        "gene_id",
        "gene_name",
        "genus",
        "gi",
        "id",
        "locus_tag",
        "organism",
        "pmid",
        "product",
        "property",
        "source",
        "source_id",
        "species"
    ],
    "spike_lineage": [
        "_version_",
        "country",
        "date_inserted",
        "date_modified",
        "growth_rate",
        "id",
        "lineage",
        "lineage_count",
        "lineage_of_concern",
        "month",
        "prevalence",
        "region",
        "sequence_features",
        "total_isolates"
    ],
    "spike_variant": [
        "_version_",
        "aa_variant",
        "country",
        "date_inserted",
        "date_modified",
        "growth_rate",
        "id",
        "lineage_count",
        "month",
        "prevalence",
        "region",
        "sequence_features",
        "total_isolates"
    ],
    "structured_assertion": [
        "_version_",
        "comment",
        "date_inserted",
        "date_modified",
        "evidence_code",
        "feature_id",
        "id",
        "owner",
        "patric_id",
        "pmid",
        "property",
        "public",
        "refseq_locus_tag",
        "score",
        "source",
        "user_read",
        "user_write",
        "value"
    ],
    "subsystem": [
        "_version_",
        "active",
        "class",
        "date_inserted",
        "date_modified",
        "feature_id",
        "gene",
        "genome_id",
        "genome_name",
        "id",
        "owner",
        "patric_id",
        "product",
        "public",
        "refseq_locus_tag",
        "role_id",
        "role_name",
        "subclass",
        "subsystem_id",
        "subsystem_name",
        "superclass",
        "taxon_id",
        "user_read",
        "user_write"
    ],
    "subsystem_ref": [
        "_version_",
        "class",
        "date_inserted",
        "date_modified",
        "description",
        "id",
        "notes",
        "pmid",
        "role_id",
        "role_name",
        "subclass",
        "subsystem_id",
        "subsystem_name",
        "superclass"
    ],
    "taxonomy": [
        "_version_",
        "cds_mean",
        "cds_sd",
        "core_families",
        "core_family_ids",
        "description",
        "division",
        "genetic_code",
        "genome_count",
        "genome_length_mean",
        "genome_length_sd",
        "genomes",
        "genomes_f",
        "hypothetical_cds_ratio_mean",
        "hypothetical_cds_ratio_sd",
        "lineage",
        "lineage_ids",
        "lineage_names",
        "lineage_ranks",
        "other_names",
        "parent_id",
        "plfam_cds_ratio_mean",
        "plfam_cds_ratio_sd",
        "taxon_id",
        "taxon_id_i",
        "taxon_name",
        "taxon_rank"
    ],
    "protein_structure": [
        "alignments",
        "authors",
        "date_inserted",
        "date_modified",
        "feature_id",
        "file_path",
        "gene",
        "genome_id",
        "institution",
        "method",
        "organism_name",
        "patric_id",
        "pdb_id",
        "pmid",
        "product",
        "release_date",
        "resolution",
        "sequence",
        "sequence_md5",
        "taxon_id",
        "taxon_lineage_ids",
        "taxon_lineage_names",
        "title",
        "uniprotkb_accession"
    ],
    "protein_feature": [
        "aa_sequence_md5",
        "classification",
        "comments",
        "date_inserted",
        "date_modified",
        "description",
        "e_value",
        "end",
        "evidence",
        "feature_id",
        "feature_type",
        "gene",
        "genome_id",
        "genome_name",
        "id",
        "interpro_description",
        "interpro_id",
        "length",
        "patric_id",
        "product",
        "publication",
        "refseq_locus_tag",
        "score",
        "segments",
        "sequence",
        "source",
        "source_id",
        "start",
        "taxon_id"
    ],
    "surveillance": [
        "additional_metadata",
        "alcohol_or_other_drug_dependence",
        "breastfeeding",
        "chest_imaging_interpretation",
        "chronic_conditions",
        "collection_city",
        "collection_country",
        "collection_date",
        "collection_latitude",
        "collection_longitude",
        "collection_poi",
        "collection_season",
        "collection_state_province",
        "collection_year",
        "collector_institution",
        "collector_name",
        "comments",
        "contact_email_address",
        "contributing_institution",
        "date_inserted",
        "date_modified",
        "daycare_attendance",
        "days_elapsed_to_disease_status",
        "days_elapsed_to_sample_collection",
        "days_elapsed_to_vaccination",
        "diagnosis",
        "dialysis",
        "disease_severity",
        "disease_status",
        "duration_of_exposure",
        "duration_of_treatment",
        "ecmo",
        "education",
        "embargo_end_date",
        "exposure",
        "exposure_type",
        "genome_id",
        "geographic_group",
        "hospitalization_duration",
        "hospitalized",
        "host_age",
        "host_capture_status",
        "host_common_name",
        "host_ethnicity",
        "host_group",
        "host_habitat",
        "host_health",
        "host_height",
        "host_id_type",
        "host_identifier",
        "host_natural_state",
        "host_race",
        "host_sex",
        "host_species",
        "host_weight",
        "human_leukocyte_antigens",
        "id",
        "infections_within_five_years",
        "influenza_like_illness_over_the_past_year",
        "initiation_of_treatment",
        "intensive_care_unit",
        "last_update_date",
        "longitudinal_study",
        "maintenance_medication",
        "nursing_home_residence",
        "onset_hours",
        "other_vaccinations",
        "oxygen_saturation",
        "packs_per_day_for_how_many_years",
        "pathogen_test_interpretation",
        "pathogen_test_result",
        "pathogen_test_type",
        "pathogen_type",
        "post_visit_medications",
        "pre_visit_medications",
        "pregnancy",
        "primary_living_situation",
        "profession",
        "project_identifier",
        "sample_accession",
        "sample_identifier",
        "sample_material",
        "sample_receipt_date",
        "sample_transport_medium",
        "sequence_accession",
        "source_of_vaccine_information",
        "species",
        "strain",
        "submission_date",
        "subtype",
        "sudden_onset",
        "symptoms",
        "taxon_lineage_ids",
        "tobacco_use",
        "travel_history",
        "treatment",
        "treatment_dosage",
        "treatment_type",
        "trimester_of_pregnancy",
        "types_of_allergies",
        "use_of_personal_protective_equipment",
        "vaccination_type",
        "vaccine_dosage",
        "vaccine_lot_number",
        "vaccine_manufacturer",
        "ventilation"
    ],
    "serology": [
        "additional_metadata",
        "collection_city",
        "collection_country",
        "collection_date",
        "collection_state",
        "collection_year",
        "comments",
        "contributing_institution",
        "date_inserted",
        "date_modified",
        "genbank_accession",
        "geographic_group",
        "host_age",
        "host_age_group",
        "host_common_name",
        "host_health",
        "host_identifier",
        "host_sex",
        "host_species",
        "host_type",
        "id",
        "positive_definition",
        "project_identifier",
        "sample_accession",
        "sample_identifier",
        "serotype",
        "strain",
        "taxon_lineage_ids",
        "test_antigen",
        "test_interpretation",
        "test_pathogen",
        "test_result",
        "test_type",
        "virus_identifier"
    ]
}
