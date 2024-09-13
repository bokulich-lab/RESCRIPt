# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from collections import OrderedDict
from io import StringIO
import os
from typing import Union
from urllib.parse import quote

import qiime2
import pandas as pd
import requests
from q2_types.feature_data import TSVTaxonomyFormat
from q2_types.genome_data import (GenomeSequencesDirectoryFormat,
                                  GenesDirectoryFormat,
                                  ProteinsDirectoryFormat, LociDirectoryFormat)

from rescript.ncbi import _allowed_ranks, _default_ranks
import json


def get_bv_brc_metadata(
        ids_metadata: Union[qiime2.NumericMetadataColumn,
        qiime2.CategoricalMetadataColumn] = None,
        data_type: str = None,
        rql_query: str = None,
        data_field: str = None,
        ids: list = None,
) -> qiime2.Metadata:
    # Parameter validation and creation of RQL query
    rql_query = parameter_validation(rql_query=rql_query,
                                     ids=ids,
                                     data_type=data_type,
                                     data_field=data_field,
                                     metadata=ids_metadata)

    # Download metadata as df
    metadata = download_data(data_type=data_type,
                             query=rql_query,
                             accept="text/tsv",
                             )

    # Set index of metadata df
    metadata.index.name = "id"
    metadata.index = metadata.index.astype(str)

    # Replace empty values and values consisting of spaces only with np.nan
    metadata.replace(to_replace=r'^\s*$', value=pd.NA, regex=True,
                     inplace=True)
    metadata.replace([None], pd.NA, inplace=True)
    return qiime2.Metadata(metadata)


def get_bv_brc_genomes(
        ids_metadata: Union[qiime2.NumericMetadataColumn,
                            qiime2.CategoricalMetadataColumn] = None,
        rql_query: str = None,
        data_field: str = None,
        ids: list = None,
        ranks: list = None,
        rank_propagation: bool = True,
) -> (GenomeSequencesDirectoryFormat, TSVTaxonomyFormat):
    # Parameter validation
    rql_query = parameter_validation(rql_query=rql_query,
                                     ids=ids,
                                     data_type="genome",
                                     data_field=data_field,
                                     metadata=ids_metadata
                                     )

    response_genomes = download_data(data_type="genome",
                                     query=rql_query,
                                     accept="application/json",
                                     select=["genome_id", "taxon_id"]
                                     )

    # Get genome sequences and create FASTA files
    genomes = get_genome_sequences(response_genomes=response_genomes)

    # Get taxonomy for sequences
    taxonomy = get_taxonomy(response_sequences=response_genomes,
                            ranks=ranks,
                            rank_propagation=rank_propagation,
                            accession_name="accession")

    return genomes, taxonomy


def get_genome_sequences(response_genomes):
    # Extract genome ids from response (list of dicts)
    genome_ids = set([str(entry['genome_id']) for entry in response_genomes])

    # Fetch the genome sequences for all genome ids
    genome_sequences = download_data(
        data_type="genome_sequence",
        query=f"in(genome_id,({','.join(genome_ids)}))",
        accept="application/json",
        select=["accession", "description", "genome_name",
                "genome_id", "sequence"]
    )

    # Create FASTA files from sequences
    genomes = create_genome_fasta(genome_sequences=genome_sequences)

    return genomes


def get_bv_brc_genome_features(
        ids_metadata: Union[qiime2.NumericMetadataColumn,
                            qiime2.CategoricalMetadataColumn] = None,
        rql_query: str = None,
        data_field: str = None,
        ids: list = None,
        ranks: list = None,
        rank_propagation: bool = True,
) -> (
        GenesDirectoryFormat, ProteinsDirectoryFormat,
        TSVTaxonomyFormat, LociDirectoryFormat):
    # Parameter validation
    rql_query = parameter_validation(rql_query=rql_query,
                                     ids=ids,
                                     data_type="genome_feature",
                                     data_field=data_field,
                                     metadata=ids_metadata
                                     )

    # Download genome_features data object as JSON
    genome_features = download_data(data_type="genome_feature",
                                    query=rql_query,
                                    accept="application/json",
                                    select=["genome_id", "feature_id",
                                            "aa_sequence_md5",
                                            "na_sequence_md5", "taxon_id"]
                                    )

    # Download nucleotide and protein sequences for features
    genes, proteins = get_sequences(genome_features=genome_features)

    # Download taxonomy for feature sequences
    taxonomy = get_taxonomy(response_sequences=genome_features,
                            ranks=ranks,
                            rank_propagation=rank_propagation,
                            accession_name="feature_id"
                            )

    # Download GFF files for genome ids
    loci = get_loci(response_sequences=genome_features)

    return genes, proteins, taxonomy, loci


def get_sequences(genome_features):
    genes = GenesDirectoryFormat()
    proteins = ProteinsDirectoryFormat()

    # Extract md5 values for nucleotide and protein sequences from
    # genome_features JSON
    md5_ids = set([item[key] for item in genome_features for key in
                   ['aa_sequence_md5', 'na_sequence_md5'] if
                   key in item])

    # Download sequences corresponding to the md5 values
    feature_sequences = download_data(data_type="feature_sequence",
                                      query=f"in(md5,({','.join(md5_ids)}))",
                                      accept="application/json",
                                      select=["md5", "sequence"])

    # Create dict with md5 id as keys and sequence as value for faster look up
    lookup_dict = {entry['md5']: entry['sequence'] for entry in
                   feature_sequences}

    fasta_genes = {}
    fasta_proteins = {}

    # Loop through list of dicts
    for entry in genome_features:
        # Extract genome_id
        genome_id = entry["genome_id"]

        if 'na_sequence_md5' in entry:
            # Create FASTA entry for the sequence
            fasta_na = (f">{entry['feature_id']}\n"
                        f"{lookup_dict[entry['na_sequence_md5']].upper()}\n")

            # Add FASTA entry to entries with the same genome id
            fasta_genes[genome_id] = fasta_genes.get(genome_id, "") + fasta_na

        if 'aa_sequence_md5' in entry:
            # Create FASTA entry for the sequence
            fasta_aa = (f">{entry['feature_id']}\n"
                        f"{lookup_dict[entry['aa_sequence_md5']].upper()}\n")

            # Add FASTA entry to entries with the same genome id
            fasta_proteins[genome_id] = (fasta_proteins.get(genome_id, "") +
                                         fasta_aa)

    # Save genes and proteins as FASTA files one file per genome_id
    for genome_id, fasta_sequences in fasta_genes.items():
        with open(os.path.join(str(genes), f"{genome_id}.fasta"),
                  'w') as fasta_file:
            fasta_file.write(fasta_sequences)

    for genome_id, fasta_sequences in fasta_proteins.items():
        with open(os.path.join(str(proteins), f"{genome_id}.fasta"),
                  'w') as fasta_file:
            fasta_file.write(fasta_sequences)

    return genes, proteins


def get_loci(response_sequences):
    # Init loci dir format
    loci = LociDirectoryFormat()

    # Extract genome ids from genome feature JSON
    genome_ids = set([str(entry['genome_id']) for entry in response_sequences])

    # Download GFF files for all genome ids. For every id there is a separate
    # request
    for genome_id in genome_ids:
        gff = download_data(data_type="genome_feature",
                            query=f"eq(genome_id,{genome_id})",
                            accept="application/gff",
                            )

        # Process the loci data and write to file
        with open(os.path.join(str(loci), f"{genome_id}.gff"), 'w') as file:
            file.write(process_loci(gff_string=gff))

    return loci


def process_loci(gff_string):
    # Split the string into lines
    lines = gff_string.splitlines()
    modified_lines = []

    for line in lines:
        if line.startswith("#"):
            # Keep header lines unchanged
            modified_lines.append(line)
        else:
            # Remove the first five characters "accn|"
            modified_lines.append(line[5:])

    # Join the lines back into a single string
    return "\n".join(modified_lines)


def get_taxonomy(response_sequences, ranks, rank_propagation, accession_name):
    # Extract all taxon_ids from response_sequences JSON
    taxon_ids = set([str(entry['taxon_id']) for entry in response_sequences])

    # Download taxonomy for all taxon_ids as df
    taxonomy = download_data(data_type="taxonomy",
                             query=f"in(taxon_id,({','.join(taxon_ids)}))",
                             accept="text/tsv",
                             select=["taxon_id", "lineage_names",
                                     "lineage_ranks"])

    # Transform the df to conform with TSVTaxonomyFormat
    return create_taxonomy(
        taxonomy_bvbrc=taxonomy,
        response_sequences=response_sequences,
        ranks=ranks,
        rank_propagation=rank_propagation,
        accession_name=accession_name,
    )


def create_taxonomy_entry(
        lineage_names,
        lineage_ranks,
        ranks=None,
        rank_propagation=True
):
    # Set ranks to default if no list is specified
    if not ranks:
        ranks = _default_ranks

    # Split the lineage names and ranks by ';'
    lineage_split = lineage_names.split(';')
    rank_split = lineage_ranks.split(';')

    # Initialize the taxonomy dictionary with the allowed ranks
    taxonomy = OrderedDict((r, None) for r in ranks)

    # Loop over each rank and assign the corresponding name
    for rank, name in zip(rank_split, lineage_split):
        if rank in ranks:
            taxonomy[rank] = name

    # Handle genus and species splitting logic
    if 'genus' in ranks and taxonomy.get('species'):
        species = taxonomy.get('species')
        # If genus includes the species name cut it out and put it into species
        if taxonomy.get('genus'):
            if species.startswith(taxonomy['genus'] + ' '):
                species = species[len(taxonomy['genus']) + 1:]
                taxonomy['species'] = species
        # If species includes the genus name cut it out and put it into genus
        elif ' ' in species:
            genus, species = species.split(' ', 1)
            taxonomy['genus'] = genus
            taxonomy['species'] = species

    # Apply rank propagation if enabled. Assigns last higher rank to lower
    # undefined ranks
    if rank_propagation:
        last_label = None
        for rank in taxonomy:
            if taxonomy[rank] is None:
                taxonomy[rank] = last_label
            last_label = taxonomy[rank]

    # Create taxonomy entry
    taxonomy_entry = '; '.join(
        f"{_allowed_ranks.get(rank, '')}{name if name else ''}" for rank, name
        in taxonomy.items()
    )
    return taxonomy_entry


def create_taxonomy(taxonomy_bvbrc, response_sequences,
                    ranks, rank_propagation, accession_name):
    # Create qiime style taxonomy entries for all rows in taxonomy_bvbrc
    taxonomy_bvbrc['Taxon'] = (
        taxonomy_bvbrc.apply(lambda row:
                             create_taxonomy_entry(
                                 lineage_names=row['lineage_names'],
                                 lineage_ranks=row['lineage_ranks'],
                                 rank_propagation=rank_propagation,
                                 ranks=ranks), axis=1))

    taxonomy_df = pd.DataFrame(columns=['Feature ID', 'Taxon'])

    # Loop through each JSON dictionary in the list
    for entry in response_sequences:
        # Get the accession and taxon_id from the JSON dictionary
        accession = entry.get(accession_name)
        taxon_id = str(entry.get('taxon_id'))

        # Look up the corresponding taxon in taxonomy_bvbrc using taxon_id
        taxon_name = taxonomy_bvbrc.loc[
            taxonomy_bvbrc['taxon_id'] == taxon_id, 'Taxon'].values

        # Create a new row as a DataFrame with accession and corresponding
        # taxon
        new_row = pd.DataFrame(
            {'Feature ID': accession, 'Taxon': taxon_name})

        # Append the new row to the taxonomy_df
        taxonomy_df = pd.concat([taxonomy_df, new_row], ignore_index=True)

    # Set index
    taxonomy_df = taxonomy_df.set_index('Feature ID')

    # Write taxonomy_df to taxonomy format
    taxonomy = TSVTaxonomyFormat()

    with taxonomy.open() as file_handle:
        file_handle.write(taxonomy_df.to_csv(sep="\t"))

    return taxonomy


def create_genome_fasta(genome_sequences):
    genomes = GenomeSequencesDirectoryFormat()

    # Dictionary to hold sequences grouped by genome_id
    fasta_entries = {}

    # Loop over all dicts in genome_sequences
    for entry in genome_sequences:
        genome_id = entry['genome_id']
        if genome_id not in fasta_entries:
            fasta_entries[genome_id] = []

        # Construct FASTA format to be identical to BV-BRC FASTA headers
        header = (f">{entry['accession']}   {entry['description']}   "
                  f"[{entry['genome_name']} | {entry['genome_id']}]")

        # Append FASTA entry to fasta_entries
        fasta_entries[genome_id].append(
            f"{header}\n{entry['sequence'].upper()}")

    # Write each genome_id's sequences to a separate FASTA file
    for genome_id, sequences in fasta_entries.items():
        fasta_content = "\n".join(sequences)
        fasta_filename = os.path.join(str(genomes), f"{genome_id}.fasta")

        with open(fasta_filename, 'w') as fasta_file:
            fasta_file.write(fasta_content)

    return genomes


def download_data(data_type, query, accept, select=None):
    # URL and headers for requests
    base_url = "https://www.bv-brc.org/api/"
    url = base_url + data_type + "/"
    headers = {'Content-Type': 'application/rqlquery+x-www-form-urlencoded',
               'ACCEPT': accept}

    results = [] if accept == "application/json" else pd.DataFrame()

    # BV-BRC sets an upper limit of 25000 entries that can be fetched at once
    batch_size = 25000
    start = 0

    while True:
        # Create the data string with current batch start and limit
        data = f"{query}&limit({batch_size},{start})"

        # Add select parameter to only download specified fields
        if select:
            data = data + f"&select({','.join(select)})"

        # POST request
        response = requests.post(url=url, data=data, headers=headers)

        # If the response is successful, process the data
        if response.status_code == 200:
            if accept == "application/json":
                result = response.json()
                results.extend(result)  # Add the current batch to the list
            elif accept == "text/tsv":
                result = read_tsv_data_with_dtypes(response=response,
                                                   data_type=data_type)

                results = pd.concat([results, result], ignore_index=True)
            else:  # application == application/gff
                return response.text

            # If the number of results is less than the batch size, break the
            # loop
            if len(result) < batch_size:
                break

            # Increment the start for the next batch
            start += batch_size

        # Handle errors
        elif response.status_code == 400:
            error_handling(response=response, data_type=data_type)
            break
        else:
            raise ValueError(response.text)

    return results


def read_tsv_data_with_dtypes(response, data_type):
    tsv_data = StringIO(response.text)

    # Read only the header to get a list of the column names
    columns_in_file = pd.read_csv(tsv_data, sep='\t', nrows=0).columns.tolist()

    # Filter the dtype dictionary to include only the columns that
    # exist in the file
    dtype_dict = data_fields_bvbrc.get(data_type, {})
    filtered_dtype_dict = {col: dtype_dict[col] for col in columns_in_file if
                           col in dtype_dict}

    # Move the file pointer to the beginning of the file-like object
    tsv_data.seek(0)

    # Read the entire file with the filtered dtype dictionary
    df = pd.read_csv(tsv_data, sep='\t', dtype=filtered_dtype_dict)

    # Raise value error if no data was retrieved
    if len(df) == 0:
        raise ValueError("No data could be retrieved. Either because of an "
                         "incorrect RQL query or because no data exists for "
                         "the query.")

    return df


def error_handling(response, data_type):
    # No data found for query or incorrect RQL query
    if response.text == "[]":
        raise ValueError("No data could be retrieved. Either because of an "
                         "incorrect RQL query or because no data exists for "
                         "the query.")

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
                f"Error code {response_dict['code']}: {response_dict['msg']}. "
                f"\nAllowed fields for data type {data_type}: "
                f"\n{list(data_fields_bvbrc[data_type].keys())}"
            )

        # Handling any other errors that start with "A Database Error Occured:"
        else:
            raise ValueError(
                f"Error code {response_dict['code']}: {response_dict['msg']}."
            )

    # Handling any other error codes
    else:
        raise ValueError(response.text)


def parameter_validation(rql_query=None,
                         ids=None,
                         data_type=None,
                         data_field=None,
                         metadata=None):
    # Error if any data_type is None
    if not data_type:
        raise ValueError("Parameter 'data-type' has to be specified.")

    local = locals().copy()
    # Error if any other parameter is specified simultaneously with rql_query
    # or metadata
    for parameter_1 in ["rql_query", "metadata"]:
        for parameter_2, value in local.items():
            if (parameter_2 != parameter_1 and
                    parameter_2 != "data_type" and
                    local[parameter_1] is not None and
                    value is not None):
                raise ValueError(
                    f"Parameters '{parameter_1}' and '{parameter_2}' "
                    "can't be used simultaneously.")

    # Error if ids or data_fields is specified without the other
    if (ids and not data_field) or (data_field and not ids):
        raise ValueError(
            "If parameter 'ids' is given, parameter 'data-field' has to be "
            "specified and vice versa.")

    # Error if rql_query, metadata and ids parameters are not given
    if not rql_query and not ids and not metadata:
        raise ValueError("At least one of the parameters 'rql-query', 'ids' "
                         "or 'ids_metadata' has to be specified.")

    # Extract data_field and ids from metadata
    if metadata is not None:
        data_field = metadata.to_series().name
        ids = metadata.to_series()

    # Check if data_field is allowed for specified data_type
    if (data_field is not None and
            data_field not in data_fields_bvbrc[data_type].keys()):
        raise ValueError(
            f"The data-field '{data_field}' is not permitted for the "
            f"data-type '{data_type}'.\nAllowed data fields are: "
            f"{list(data_fields_bvbrc[data_type].keys())}")

    # Construct the RQL query
    if ids is not None or metadata is not None:
        # Join the quoted ids with commas and percent encode the values
        joined_ids = ','.join(quote(f'"{str(id_)}"') for id_ in ids)

        # Final query
        rql_query = f'in({data_field},({joined_ids}))'

    return rql_query


data_fields_bvbrc = {
    "antibiotics": {
        "_version_": str,
        "antibiotic_name": str,
        "atc_classification": str,
        "canonical_smiles": str,
        "cas_id": str,
        "date_inserted": str,
        "date_modified": str,
        "description": str,
        "drugbank_interactions": str,
        "inchi_key": str,
        "isomeric_smiles": str,
        "mechanism_of_action": str,
        "molecular_formula": str,
        "molecular_weight": str,
        "pharmacological_classes": str,
        "pharmacology": str,
        "pubchem_cid": str,
        "pubchem_cid_i": str,
        "synonyms": str
    },
    "enzyme_class_ref": {
        "_version_": str,
        "date_inserted": str,
        "date_modified": str,
        "ec_number": str,
        "go": str
    },
    "epitope": {
        "_version_": str,
        "assay_results": str,
        "bcell_assays": str,
        "comments": str,
        "date_inserted": str,
        "date_modified": str,
        "end": float,
        "epitope_id": str,
        "epitope_sequence": str,
        "epitope_type": str,
        "host_name": str,
        "mhc_assays": str,
        "organism": str,
        "protein_accession": str,
        "protein_id": str,
        "protein_name": str,
        "start": float,
        "taxon_id": str,
        "taxon_lineage_ids": str,
        "taxon_lineage_names": str,
        "tcell_assays": str,
        "total_assays": float
    },
    "epitope_assay": {
        "_version_": str,
        "assay_group": str,
        "assay_id": str,
        "assay_measurement": str,
        "assay_measurement_unit": str,
        "assay_method": str,
        "assay_result": str,
        "assay_type": str,
        "authors": str,
        "date_inserted": str,
        "date_modified": str,
        "end": float,
        "epitope_id": str,
        "epitope_sequence": str,
        "epitope_type": str,
        "host_name": str,
        "host_taxon_id": str,
        "mhc_allele": str,
        "mhc_allele_class": str,
        "organism": str,
        "pdb_id": str,
        "pmid": str,
        "protein_accession": str,
        "protein_id": str,
        "protein_name": str,
        "start": float,
        "taxon_id": str,
        "taxon_lineage_ids": str,
        "taxon_lineage_names": str,
        "title": str
    },
    "experiment": {
        "_version_": str,
        "additional_data": str,
        "additional_metadata": str,
        "biosets": float,
        "date_inserted": str,
        "date_modified": str,
        "detection_instrument": str,
        "doi": str,
        "exp_description": str,
        "exp_id": str,
        "exp_name": str,
        "exp_poc": str,
        "exp_protocol": str,
        "exp_title": str,
        "exp_type": str,
        "experimenters": str,
        "genome_id": str,
        "measurement_technique": str,
        "organism": str,
        "pmid": str,
        "public_identifier": str,
        "public_repository": str,
        "samples": float,
        "strain": str,
        "study_description": str,
        "study_institution": str,
        "study_name": str,
        "study_pi": str,
        "study_title": str,
        "taxon_id": str,
        "taxon_lineage_ids": str,
        "treatment_amount": str,
        "treatment_duration": str,
        "treatment_name": str,
        "treatment_type": str
    },
    "bioset": {
        "_version_": str,
        "additional_data": str,
        "additional_metadata": str,
        "analysis_group_1": str,
        "analysis_group_2": str,
        "analysis_method": str,
        "bioset_criteria": str,
        "bioset_description": str,
        "bioset_id": str,
        "bioset_name": str,
        "bioset_result": str,
        "bioset_type": str,
        "date_inserted": str,
        "date_modified": str,
        "entity_count": str,
        "entity_type": str,
        "exp_id": str,
        "exp_name": str,
        "exp_title": str,
        "exp_type": str,
        "genome_id": str,
        "organism": str,
        "protocol": str,
        "result_type": str,
        "strain": str,
        "study_description": str,
        "study_institution": str,
        "study_name": str,
        "study_pi": str,
        "study_title": str,
        "taxon_id": str,
        "taxon_lineage_ids": str,
        "treatment_amount": str,
        "treatment_duration": str,
        "treatment_name": str,
        "treatment_type": str
    },
    "bioset_result": {
        "_version_": str,
        "bioset_description": str,
        "bioset_id": str,
        "bioset_name": str,
        "bioset_type": str,
        "counts": float,
        "date_inserted": str,
        "date_modified": str,
        "entity_id": str,
        "entity_name": str,
        "entity_type": str,
        "exp_id": str,
        "exp_name": str,
        "exp_title": str,
        "exp_type": str,
        "feature_id": str,
        "fpkm": float,
        "gene": str,
        "gene_id": str,
        "genome_id": str,
        "id": str,
        "locus_tag": str,
        "log2_fc": float,
        "organism": str,
        "other_ids": str,
        "other_value": float,
        "p_value": float,
        "patric_id": str,
        "product": str,
        "protein_id": str,
        "result_type": str,
        "strain": str,
        "taxon_id": str,
        "tpm": float,
        "treatment_amount": str,
        "treatment_duration": str,
        "treatment_name": str,
        "treatment_type": str,
        "uniprot_id": str,
        "z_score": float
    },
    "gene_ontology_ref": {
        "_version_": str,
        "date_inserted": str,
        "date_modified": str,
        "definition": str,
        "go_id": str,
        "go_name": str,
        "ontology": str
    },
    "genome": {
        "_version_": str,
        "additional_metadata": str,
        "altitude": str,
        "antimicrobial_resistance": str,
        "antimicrobial_resistance_evidence": str,
        "assembly_accession": str,
        "assembly_method": str,
        "authors": str,
        "bioproject_accession": str,
        "biosample_accession": str,
        "biovar": str,
        "body_sample_site": str,
        "body_sample_subsite": str,
        "cds": float,
        "cds_ratio": float,
        "cell_shape": str,
        "checkm_completeness": float,
        "checkm_contamination": float,
        "chromosomes": float,
        "clade": str,
        "class": str,
        "coarse_consistency": float,
        "collection_date": str,
        "collection_year": float,
        "comments": str,
        "common_name": str,
        "completion_date": str,
        "contig_l50": float,
        "contig_n50": float,
        "contigs": float,
        "core_families": float,
        "core_family_ratio": float,
        "culture_collection": str,
        "date_inserted": str,
        "date_modified": str,
        "depth": str,
        "disease": str,
        "family": str,
        "fine_consistency": float,
        "gc_content": float,
        "genbank_accessions": str,
        "genome_id": str,
        "genome_length": float,
        "genome_name": str,
        "genome_quality": str,
        "genome_quality_flags": str,
        "genome_status": str,
        "genus": str,
        "geographic_group": str,
        "geographic_location": str,
        "gram_stain": str,
        "h1_clade_global": str,
        "h1_clade_us": str,
        "h3_clade": str,
        "h5_clade": str,
        "h_type": float,
        "habitat": str,
        "host_age": str,
        "host_common_name": str,
        "host_gender": str,
        "host_group": str,
        "host_health": str,
        "host_name": str,
        "host_scientific_name": str,
        "hypothetical_cds": float,
        "hypothetical_cds_ratio": float,
        "isolation_comments": str,
        "isolation_country": str,
        "isolation_site": str,
        "isolation_source": str,
        "kingdom": str,
        "lab_host": str,
        "latitude": str,
        "lineage": str,
        "longitude": str,
        "mat_peptide": float,
        "missing_core_family_ids": str,
        "mlst": str,
        "motility": str,
        "n_type": float,
        "ncbi_project_id": str,
        "nearest_genomes": str,
        "optimal_temperature": str,
        "order": str,
        "organism_name": str,
        "other_clinical": str,
        "other_environmental": str,
        "other_names": str,
        "other_typing": str,
        "outgroup_genomes": str,
        "owner": str,
        "oxygen_requirement": str,
        "p2_genome_id": str,
        "partial_cds": float,
        "partial_cds_ratio": float,
        "passage": str,
        "pathovar": str,
        "patric_cds": float,
        "ph1n1_like": str,
        "phenotype": str,
        "phylum": str,
        "plasmids": float,
        "plfam_cds": float,
        "plfam_cds_ratio": float,
        "public": str,
        "publication": str,
        "reference_genome": str,
        "refseq_accessions": str,
        "refseq_cds": float,
        "refseq_project_id": str,
        "rrna": float,
        "salinity": str,
        "season": str,
        "segment": str,
        "segments": float,
        "sequencing_centers": str,
        "sequencing_depth": str,
        "sequencing_platform": str,
        "sequencing_status": str,
        "serovar": str,
        "species": str,
        "sporulation": str,
        "sra_accession": str,
        "state_province": str,
        "strain": str,
        "subclade": str,
        "subtype": str,
        "superkingdom": str,
        "taxon_id": str,
        "taxon_lineage_ids": str,
        "taxon_lineage_names": str,
        "temperature_range": str,
        "trna": float,
        "type_strain": str,
        "user_read": str,
        "user_write": str
    },
    "strain": {
        "1_pb2": str,
        "2_pb1": str,
        "3_pa": str,
        "4_ha": str,
        "5_np": str,
        "6_na": str,
        "7_mp": str,
        "8_ns": str,
        "_version_": str,
        "collection_date": str,
        "collection_year": float,
        "date_inserted": str,
        "date_modified": str,
        "family": str,
        "genbank_accessions": str,
        "genome_ids": str,
        "genus": str,
        "geographic_group": str,
        "h_type": float,
        "host_common_name": str,
        "host_group": str,
        "host_name": str,
        "id": str,
        "isolation_country": str,
        "l": str,
        "lab_host": str,
        "m": str,
        "n_type": float,
        "other_segments": str,
        "owner": str,
        "passage": str,
        "public": str,
        "s": str,
        "season": str,
        "segment_count": float,
        "species": str,
        "status": str,
        "strain": str,
        "subtype": str,
        "taxon_id": str,
        "taxon_lineage_ids": str,
        "taxon_lineage_names": str,
        "user_read": str,
        "user_write": str
    },
    "genome_amr": {
        "_version_": str,
        "antibiotic": str,
        "computational_method": str,
        "computational_method_performance": str,
        "computational_method_version": str,
        "date_inserted": str,
        "date_modified": str,
        "evidence": str,
        "genome_id": str,
        "genome_name": str,
        "id": str,
        "laboratory_typing_method": str,
        "laboratory_typing_method_version": str,
        "laboratory_typing_platform": str,
        "measurement": str,
        "measurement_sign": str,
        "measurement_unit": str,
        "measurement_value": str,
        "owner": str,
        "pmid": str,
        "public": str,
        "resistant_phenotype": str,
        "source": str,
        "taxon_id": str,
        "testing_standard": str,
        "testing_standard_year": float,
        "user_read": str,
        "user_write": str,
        "vendor": str
    },
    "feature_sequence": {
        "_version_": str,
        "date_inserted": str,
        "date_modified": str,
        "md5": str,
        "sequence": str,
        "sequence_type": str
    },
    "genome_feature": {
        "aa_length": float,
        "aa_sequence_md5": str,
        "accession": str,
        "alt_locus_tag": str,
        "annotation": str,
        "brc_id": str,
        "classifier_round": float,
        "classifier_score": float,
        "codon_start": float,
        "date_inserted": str,
        "date_modified": str,
        "end": float,
        "feature_id": str,
        "feature_type": str,
        "figfam_id": str,
        "gene": str,
        "gene_id": str,
        "genome_id": str,
        "genome_name": str,
        "go": str,
        "location": str,
        "na_length": float,
        "na_sequence_md5": str,
        "notes": str,
        "og_id": str,
        "owner": str,
        "p2_feature_id": str,
        "patric_id": str,
        "pdb_accession": str,
        "pgfam_id": str,
        "plfam_id": str,
        "product": str,
        "property": str,
        "protein_id": str,
        "public": str,
        "refseq_locus_tag": str,
        "segments": str,
        "sequence_id": str,
        "sog_id": str,
        "start": float,
        "strand": str,
        "taxon_id": str,
        "uniprotkb_accession": str,
        "user_read": str,
        "user_write": str
    },
    "genome_sequence": {
        "_version_": str,
        "accession": str,
        "chromosome": str,
        "date_inserted": str,
        "date_modified": str,
        "description": str,
        "gc_content": float,
        "genome_id": str,
        "genome_name": str,
        "gi": float,
        "length": float,
        "mol_type": str,
        "owner": str,
        "p2_sequence_id": str,
        "plasmid": str,
        "public": str,
        "release_date": str,
        "segment": str,
        "sequence": str,
        "sequence_id": str,
        "sequence_md5": str,
        "sequence_status": str,
        "sequence_type": str,
        "taxon_id": str,
        "topology": str,
        "user_read": str,
        "user_write": str,
        "version": str
    },
    "id_ref": {
        "_version_": str,
        "date_inserted": str,
        "date_modified": str,
        "id": str,
        "id_type": str,
        "id_value": str,
        "uniprotkb_accession": str
    },
    "misc_niaid_sgc": {
        "_version_": str,
        "date_inserted": str,
        "date_modified": str,
        "gene_symbol_collection": str,
        "genus": str,
        "has_clones": str,
        "has_proteins": str,
        "selection_criteria": str,
        "species": str,
        "strain": str,
        "target_id": str,
        "target_status": str
    },
    "pathway": {
        "_version_": str,
        "accession": str,
        "alt_locus_tag": str,
        "annotation": str,
        "date_inserted": str,
        "date_modified": str,
        "ec_description": str,
        "ec_number": str,
        "feature_id": str,
        "gene": str,
        "genome_ec": str,
        "genome_id": str,
        "genome_name": str,
        "id": str,
        "owner": str,
        "pathway_class": str,
        "pathway_ec": str,
        "pathway_id": str,
        "pathway_name": str,
        "patric_id": str,
        "product": str,
        "public": str,
        "refseq_locus_tag": str,
        "sequence_id": str,
        "taxon_id": str,
        "user_read": str,
        "user_write": str
    },
    "pathway_ref": {
        "_version_": str,
        "date_inserted": str,
        "date_modified": str,
        "ec_description": str,
        "ec_number": str,
        "id": str,
        "map_location": str,
        "map_name": str,
        "map_type": str,
        "occurrence": float,
        "pathway_class": str,
        "pathway_id": str,
        "pathway_name": str
    },
    "ppi": {
        "_version_": str,
        "category": str,
        "date_inserted": str,
        "date_modified": str,
        "detection_method": str,
        "domain_a": str,
        "domain_b": str,
        "evidence": str,
        "feature_id_a": str,
        "feature_id_b": str,
        "gene_a": str,
        "gene_b": str,
        "genome_id_a": str,
        "genome_id_b": str,
        "genome_name_a": str,
        "genome_name_b": str,
        "id": str,
        "interaction_type": str,
        "interactor_a": str,
        "interactor_b": str,
        "interactor_desc_a": str,
        "interactor_desc_b": str,
        "interactor_type_a": str,
        "interactor_type_b": str,
        "pmid": str,
        "refseq_locus_tag_a": str,
        "refseq_locus_tag_b": str,
        "score": str,
        "source_db": str,
        "source_id": str,
        "taxon_id_a": str,
        "taxon_id_b": str
    },
    "protein_family_ref": {
        "_version_": str,
        "date_inserted": str,
        "date_modified": str,
        "family_id": str,
        "family_product": str,
        "family_type": str
    },
    "sequence_feature": {
        "aa_sequence_md5": str,
        "aa_variant": str,
        "additional_metadata": str,
        "comments": str,
        "date_inserted": str,
        "date_modified": str,
        "end": float,
        "evidence_code": str,
        "feature_id": str,
        "genbank_accession": str,
        "gene": str,
        "genome_id": str,
        "genome_name": str,
        "id": str,
        "length": float,
        "patric_id": str,
        "product": str,
        "publication": str,
        "refseq_locus_tag": str,
        "segment": str,
        "segments": str,
        "sf_category": str,
        "sf_id": str,
        "sf_name": str,
        "sf_sequence": str,
        "sf_sequence_md5": str,
        "source": str,
        "source_aa_sequence": str,
        "source_id": str,
        "source_sf_location": str,
        "source_strain": str,
        "start": float,
        "subtype": str,
        "taxon_id": str,
        "variant_types": str
    },
    "sequence_feature_vt": {
        "additional_metadata": str,
        "comments": str,
        "date_inserted": str,
        "date_modified": str,
        "id": str,
        "sf_category": str,
        "sf_id": str,
        "sf_name": str,
        "sf_sequence": str,
        "sf_sequence_md5": str,
        "sfvt_genome_count": str,
        "sfvt_genome_ids": str,
        "sfvt_id": str,
        "sfvt_sequence": str,
        "sfvt_sequence_md5": str,
        "sfvt_variations": str
    },
    "sp_gene": {
        "_version_": str,
        "alt_locus_tag": str,
        "antibiotics": str,
        "antibiotics_class": str,
        "classification": str,
        "date_inserted": str,
        "date_modified": str,
        "e_value": str,
        "evidence": str,
        "feature_id": str,
        "function": str,
        "gene": str,
        "genome_id": str,
        "genome_name": str,
        "id": str,
        "identity": float,
        "organism": str,
        "owner": str,
        "patric_id": str,
        "pmid": str,
        "product": str,
        "property": str,
        "property_source": str,
        "public": str,
        "query_coverage": float,
        "refseq_locus_tag": str,
        "same_genome": float,
        "same_genus": float,
        "same_species": float,
        "source": str,
        "source_id": str,
        "subject_coverage": float,
        "taxon_id": str,
        "user_read": str,
        "user_write": str
    },
    "sp_gene_ref": {
        "_version_": str,
        "antibiotics": str,
        "antibiotics_class": str,
        "assertion": str,
        "classification": str,
        "date_inserted": str,
        "date_modified": str,
        "function": str,
        "gene_id": str,
        "gene_name": str,
        "genus": str,
        "gi": str,
        "id": str,
        "locus_tag": str,
        "organism": str,
        "pmid": str,
        "product": str,
        "property": str,
        "source": str,
        "source_id": str,
        "species": str
    },
    "spike_lineage": {
        "_version_": str,
        "country": str,
        "date_inserted": str,
        "date_modified": str,
        "growth_rate": float,
        "id": str,
        "lineage": str,
        "lineage_count": float,
        "lineage_of_concern": str,
        "month": str,
        "prevalence": float,
        "region": str,
        "sequence_features": str,
        "total_isolates": float
    },
    "spike_variant": {
        "_version_": str,
        "aa_variant": str,
        "country": str,
        "date_inserted": str,
        "date_modified": str,
        "growth_rate": float,
        "id": str,
        "lineage_count": float,
        "month": str,
        "prevalence": float,
        "region": str,
        "sequence_features": str,
        "total_isolates": float
    },
    "structured_assertion": {
        "_version_": str,
        "comment": str,
        "date_inserted": str,
        "date_modified": str,
        "evidence_code": str,
        "feature_id": str,
        "id": str,
        "owner": str,
        "patric_id": str,
        "pmid": str,
        "property": str,
        "public": str,
        "refseq_locus_tag": str,
        "score": str,
        "source": str,
        "user_read": str,
        "user_write": str,
        "value": str
    },
    "subsystem": {
        "_version_": str,
        "active": str,
        "class": str,
        "date_inserted": str,
        "date_modified": str,
        "feature_id": str,
        "gene": str,
        "genome_id": str,
        "genome_name": str,
        "id": str,
        "owner": str,
        "patric_id": str,
        "product": str,
        "public": str,
        "refseq_locus_tag": str,
        "role_id": str,
        "role_name": str,
        "subclass": str,
        "subsystem_id": str,
        "subsystem_name": str,
        "superclass": str,
        "taxon_id": str,
        "user_read": str,
        "user_write": str
    },
    "subsystem_ref": {
        "_version_": str,
        "class": str,
        "date_inserted": str,
        "date_modified": str,
        "description": str,
        "id": str,
        "notes": str,
        "pmid": str,
        "role_id": str,
        "role_name": str,
        "subclass": str,
        "subsystem_id": str,
        "subsystem_name": str,
        "superclass": str
    },
    "taxonomy": {
        "_version_": str,
        "cds_mean": float,
        "cds_sd": float,
        "core_families": float,
        "core_family_ids": str,
        "description": str,
        "division": str,
        "genetic_code": str,
        "genome_count": float,
        "genome_length_mean": float,
        "genome_length_sd": float,
        "genomes": float,
        "genomes_f": str,
        "hypothetical_cds_ratio_mean": float,
        "hypothetical_cds_ratio_sd": float,
        "lineage": str,
        "lineage_ids": str,
        "lineage_names": str,
        "lineage_ranks": str,
        "other_names": str,
        "parent_id": str,
        "plfam_cds_ratio_mean": float,
        "plfam_cds_ratio_sd": float,
        "taxon_id": str,
        "taxon_id_i": str,
        "taxon_name": str,
        "taxon_rank": str
    },
    "protein_structure": {
        "alignments": str,
        "authors": str,
        "date_inserted": str,
        "date_modified": str,
        "feature_id": str,
        "file_path": str,
        "gene": str,
        "genome_id": str,
        "institution": str,
        "method": str,
        "organism_name": str,
        "patric_id": str,
        "pdb_id": str,
        "pmid": str,
        "product": str,
        "release_date": str,
        "resolution": str,
        "sequence": str,
        "sequence_md5": str,
        "taxon_id": str,
        "taxon_lineage_ids": str,
        "taxon_lineage_names": str,
        "title": str,
        "uniprotkb_accession": str
    },
    "protein_feature": {
        "aa_sequence_md5": str,
        "classification": str,
        "comments": str,
        "date_inserted": str,
        "date_modified": str,
        "description": str,
        "e_value": str,
        "end": float,
        "evidence": str,
        "feature_id": str,
        "feature_type": str,
        "gene": str,
        "genome_id": str,
        "genome_name": str,
        "id": str,
        "interpro_description": str,
        "interpro_id": str,
        "length": float,
        "patric_id": str,
        "product": str,
        "publication": str,
        "refseq_locus_tag": str,
        "score": float,
        "segments": str,
        "sequence": str,
        "source": str,
        "source_id": str,
        "start": float,
        "taxon_id": str
    },
    "surveillance": {
        "additional_metadata": str,
        "alcohol_or_other_drug_dependence": str,
        "breastfeeding": str,
        "chest_imaging_interpretation": str,
        "chronic_conditions": str,
        "collection_city": str,
        "collection_country": str,
        "collection_date": str,
        "collection_latitude": float,
        "collection_longitude": float,
        "collection_poi": str,
        "collection_season": str,
        "collection_state_province": str,
        "collection_year": str,
        "collector_institution": str,
        "collector_name": str,
        "comments": str,
        "contact_email_address": str,
        "contributing_institution": str,
        "date_inserted": str,
        "date_modified": str,
        "daycare_attendance": str,
        "days_elapsed_to_disease_status": str,
        "days_elapsed_to_sample_collection": str,
        "days_elapsed_to_vaccination": str,
        "diagnosis": str,
        "dialysis": str,
        "disease_severity": str,
        "disease_status": str,
        "duration_of_exposure": str,
        "duration_of_treatment": str,
        "ecmo": str,
        "education": str,
        "embargo_end_date": str,
        "exposure": str,
        "exposure_type": str,
        "genome_id": str,
        "geographic_group": str,
        "hospitalization_duration": str,
        "hospitalized": str,
        "host_age": str,
        "host_capture_status": str,
        "host_common_name": str,
        "host_ethnicity": str,
        "host_group": str,
        "host_habitat": str,
        "host_health": str,
        "host_height": str,
        "host_id_type": str,
        "host_identifier": str,
        "host_natural_state": str,
        "host_race": str,
        "host_sex": str,
        "host_species": str,
        "host_weight": str,
        "human_leukocyte_antigens": str,
        "id": str,
        "infections_within_five_years": str,
        "influenza_like_illness_over_the_past_year": str,
        "initiation_of_treatment": str,
        "intensive_care_unit": str,
        "last_update_date": str,
        "longitudinal_study": str,
        "maintenance_medication": str,
        "nursing_home_residence": str,
        "onset_hours": str,
        "other_vaccinations": str,
        "oxygen_saturation": str,
        "packs_per_day_for_how_many_years": str,
        "pathogen_test_interpretation": str,
        "pathogen_test_result": str,
        "pathogen_test_type": str,
        "pathogen_type": str,
        "post_visit_medications": str,
        "pre_visit_medications": str,
        "pregnancy": str,
        "primary_living_situation": str,
        "profession": str,
        "project_identifier": str,
        "sample_accession": str,
        "sample_identifier": str,
        "sample_material": str,
        "sample_receipt_date": str,
        "sample_transport_medium": str,
        "sequence_accession": str,
        "source_of_vaccine_information": str,
        "species": str,
        "strain": str,
        "submission_date": str,
        "subtype": str,
        "sudden_onset": str,
        "symptoms": str,
        "taxon_lineage_ids": str,
        "tobacco_use": str,
        "travel_history": str,
        "treatment": str,
        "treatment_dosage": str,
        "treatment_type": str,
        "trimester_of_pregnancy": str,
        "types_of_allergies": str,
        "use_of_personal_protective_equipment": str,
        "vaccination_type": str,
        "vaccine_dosage": str,
        "vaccine_lot_number": str,
        "vaccine_manufacturer": str,
        "ventilation": str
    },
    "serology": {
        "additional_metadata": str,
        "collection_city": str,
        "collection_country": str,
        "collection_date": str,
        "collection_state": str,
        "collection_year": str,
        "comments": str,
        "contributing_institution": str,
        "date_inserted": str,
        "date_modified": str,
        "genbank_accession": str,
        "geographic_group": str,
        "host_age": str,
        "host_age_group": str,
        "host_common_name": str,
        "host_health": str,
        "host_identifier": str,
        "host_sex": str,
        "host_species": str,
        "host_type": str,
        "id": str,
        "positive_definition": str,
        "project_identifier": str,
        "sample_accession": str,
        "sample_identifier": str,
        "serotype": str,
        "strain": str,
        "taxon_lineage_ids": str,
        "test_antigen": str,
        "test_interpretation": str,
        "test_pathogen": str,
        "test_result": str,
        "test_type": str,
        "virus_identifier": str
    }
}
