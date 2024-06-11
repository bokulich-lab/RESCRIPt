# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import json
import os
import shutil
import tempfile
import zipfile
from copy import deepcopy
from typing import List, Dict

import ncbi.datasets as nd
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
import pandas as pd
import skbio
from multiprocessing import Manager
from ncbi.datasets import ApiException
from q2_types.feature_data import DNAFASTAFormat
from q2_types.genome_data import (LociDirectoryFormat,
                                  ProteinsDirectoryFormat)

from rescript.ncbi import get_taxonomies, _default_ranks


def _get_assembly_descriptors(
        api_instance, assembly_levels, assembly_source, only_reference,
        page_size, taxon, tax_exact_match
) -> Dict[str, str]:
    assembly_to_taxon = {}
    next_page_token = ''
    while True:
        try:
            genome_summary = api_instance.assembly_descriptors_by_taxon(
                taxon=taxon,
                page_size=page_size,
                filters_assembly_source=assembly_source,
                filters_assembly_level=assembly_levels,
                filters_reference_only=only_reference,
                # TODO: possibly, make this configurable - probably will
                #   require handling some genomes without annotation
                #   (maybe an issue for loci/proteins types?)
                filters_has_annotation=True,
                tax_exact_match=tax_exact_match,
                page_token=next_page_token
            )
        except ApiException as e:
            raise Exception('Unexpected error while calling NCBI Datasets '
                            f'API: {e}.')

        if not genome_summary.assemblies:
            msg = 'The query to the NCBI Dataset API did not return ' \
                  'any genomes. '
            if not genome_summary.messages:
                msg += 'It is possible that your query was too restrictive. ' \
                       'Try to adjust one of the following parameters: ' \
                       '"assembly_source", "assembly_levels", ' \
                       '"only_reference", "taxon" or "tax_exact_match" ' \
                       'and try again.'
                raise Exception(msg)

            msg += 'Please update your query. The following ' \
                   'error was reported:\n ' \
                   f'{genome_summary.messages[0].error.message}'
            raise Exception(msg)

        next_page_token = genome_summary.next_page_token

        genome_assembly = [x.assembly for x in genome_summary.assemblies]
        assembly_to_taxon.update(
            {x.assembly_accession: x.org.tax_id for x in genome_assembly}
        )

        if not next_page_token:
            break
    return assembly_to_taxon


def _extract_accession_ids(seq_reports: List[dict]):
    accessions, assembly_id = {'refseq': [], 'genbank': []}, None
    for rep in seq_reports:
        if not assembly_id:
            assembly_id = rep.get('assemblyAccession')
        accessions['refseq'].append(rep.get('refseqAccession'))
        accessions['genbank'].append(rep.get('genbankAccession'))
    return assembly_id, accessions


def _fetch_and_extract_dataset(
        api_response, genomes, loci, proteins, only_genomic=False
):
    with tempfile.TemporaryDirectory() as tmp:
        result_path = os.path.join(tmp, 'datasets.zip')
        with open(result_path, 'wb') as f:
            f.write(api_response.data)

        with zipfile.ZipFile(result_path, 'r') as zipped:
            zipped.extractall(tmp)

        # find and move all the genome sequences
        genome_seq_fps = glob.glob(
            os.path.join(tmp, 'ncbi_dataset', 'data', '*', '*_genomic.fna')
        )
        acc_to_assembly = {}
        with open(str(genomes), 'a') as fin:
            for f in genome_seq_fps:
                dp = os.path.dirname(f)
                summary_fp = os.path.join(dp, 'sequence_report.jsonl')
                with open(summary_fp, 'r') as jsonl_file:
                    molecules = [json.loads(line) for line in jsonl_file]

                if only_genomic:
                    molecules = [
                        x for x in molecules
                        if x['assignedMoleculeLocationType'] == 'Chromosome'
                    ]

                assembly_id, accession_ids = _extract_accession_ids(molecules)
                molecule_ids = [
                    *accession_ids['refseq'], *accession_ids['genbank']
                ]

                seq = skbio.read(
                    f, format='fasta', constructor=skbio.DNA, lowercase=True
                )
                acc_to_assembly[assembly_id] = []
                for s in seq:
                    _id = s.metadata['id']
                    if s.metadata['id'] in molecule_ids:
                        skbio.io.write(s, format='fasta', into=fin)
                        acc_to_assembly[assembly_id].append(_id)

        # find and move all the gff files
        loci_fps = glob.glob(
            os.path.join(tmp, 'ncbi_dataset', 'data', '*', 'genomic.gff'))
        for f in loci_fps:
            _id = f.split('/')[-2]
            shutil.move(f, os.path.join(str(loci), f'{_id}_loci.gff'))

        # find and move all the protein translations
        protein_fps = glob.glob(
            os.path.join(tmp, 'ncbi_dataset', 'data', '*', 'protein.faa'))
        for f in protein_fps:
            _id = f.split('/')[-2]
            shutil.move(f, os.path.join(
                str(proteins), f'{_id}_proteins.fasta'))

    return acc_to_assembly


def _fetch_taxonomy(
        all_acc_ids: list,
        all_tax_ids: list,
        accession_to_assembly: pd.Series
):
    manager = Manager()
    taxa, bad_accs = get_taxonomies(
        taxids={k: v for k, v in zip(all_acc_ids, all_tax_ids)},
        ranks=_default_ranks, rank_propagation=True,
        logging_level='INFO', n_jobs=2, request_lock=manager.Lock()
    )
    # technically, this should never happen as the taxa accession
    # numbers are extracted from the API response
    if bad_accs:
        raise Exception('Invalid taxonomy accession numbers were found: '
                        f'{", ".join(bad_accs)}. Please check your query '
                        f'and try again.')
    taxa = pd.DataFrame(taxa, index=['Taxon']).T
    taxa = pd.merge(
        taxa, accession_to_assembly,
        left_index=True, right_index=True, how="outer"
    ).set_index('assembly_id')
    taxa.index.name = 'Feature ID'
    return taxa


def _fetch_and_extract_all(api_instance, assemblies, only_genomic):
    genomes = DNAFASTAFormat()
    loci = LociDirectoryFormat()
    proteins = ProteinsDirectoryFormat()
    accession_map = {}

    chunks = [
        assemblies[i:i + 20] for i in range(0, len(assemblies), 20)
    ]
    for assembly_subset in chunks:
        api_response = api_instance.download_assembly_package(
            assembly_subset,
            exclude_sequence=False,
            include_annotation_type=['PROT_FASTA', 'GENOME_GFF'],
            _preload_content=False
        )
        accessions = _fetch_and_extract_dataset(
            api_response, genomes, loci, proteins,
            only_genomic=only_genomic
        )
        accession_map.update(accessions)
    accession_map = pd.Series(accession_map, name='assembly_id')

    return genomes, loci, proteins, accession_map


def get_ncbi_genomes(
        taxon: str,
        assembly_source: str = 'refseq',
        assembly_levels: List[str] = ['complete_genome'],
        only_reference: bool = True,
        only_genomic: bool = False,
        tax_exact_match: bool = False,
        page_size: int = 20,
) -> (DNAFASTAFormat, LociDirectoryFormat,
      ProteinsDirectoryFormat, pd.DataFrame):
    # we use a deepcopy of assembly_levels because the new versions of
    # the NCBI datasets toolkit somehow, magically, update it what
    # causes a conflict in the q2cli resulting in a failure to dump
    # the action's params to provenance
    assembly_descriptors_params = {
        'assembly_levels': deepcopy(assembly_levels),
        'assembly_source': assembly_source,
        'only_reference': only_reference, 'page_size': page_size,
        'taxon': taxon, 'tax_exact_match': tax_exact_match
    }
    with DatasetsApiClient() as api_client:
        api_instance = nd.GenomeApi(api_client)

        assembly_to_taxon = _get_assembly_descriptors(
            api_instance=api_instance, **assembly_descriptors_params
        )

        genomes, loci, proteins, accession_map = _fetch_and_extract_all(
            api_instance, list(assembly_to_taxon.keys()), only_genomic
        )

        taxa = _fetch_taxonomy(
            assembly_to_taxon.keys(),
            assembly_to_taxon.values(),
            accession_map.explode()
        )

    return genomes, loci, proteins, taxa
