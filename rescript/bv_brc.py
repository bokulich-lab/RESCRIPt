# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from io import StringIO

import qiime2
import pandas as pd
import requests
from q2_types.feature_data import MixedCaseDNAFASTAFormat


def json_to_fasta(json: dict):
    fasta_output = []
    for entry in json:
        header = (f">accn|{entry['sequence_id']}   {entry['description']}   "
                  f"[{entry['genome_name']} | {entry['genome_id']}]")
        fasta_output.append(f"{header}\n{entry['sequence']}")
    return "\n".join(fasta_output)


def fetch_genomes_bv_brc(rql_query: str) -> (MixedCaseDNAFASTAFormat, qiime2.Metadata):
    genomes = MixedCaseDNAFASTAFormat()

    # Make the GET request for metadata
    url_metadata = f"https://www.bv-brc.org/api/genome/{rql_query}&http_accept=text/tsv"
    response_metadata = requests.get(url_metadata)

    if response_metadata.status_code == 200:
        # Convert TSV data to dataframe
        tsv_data = StringIO(response_metadata.text)
        metadata = pd.read_csv(tsv_data, sep='\t', index_col="genome_id")

        metadata.index.name = "id"
        metadata.index = metadata.index.astype(str)

        # Extract all genome_ids out of dataframe
        genome_ids = metadata.index.tolist()
    else:
        raise ValueError("Error")

    # Make the GET request for sequences
    url_sequences = (f"https://www.bv-brc.org/api/genome_sequence/"
                     f"?in(genome_id,({','.join(genome_ids)}))")
    response_sequences = requests.get(url_sequences)

    if response_sequences.status_code == 200:
        # Convert JSON to FASTA
        fasta = json_to_fasta(response_sequences.json())

        # Write FASTA format to file
        with genomes.open() as file:
            file.write(fasta)
    else:
        raise ValueError("Error")

    return genomes, qiime2.Metadata(metadata)
