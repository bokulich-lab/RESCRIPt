# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import glob
import os
import shutil
import tempfile
import zipfile
from copy import deepcopy
from typing import List

import ncbi.datasets as nd
from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
import pandas as pd
import skbio
from multiprocessing import Manager
from ncbi.datasets import ApiException
from q2_types.feature_data import DNAFASTAFormat, DNAIterator
from q2_types.genome_data import (LociDirectoryFormat,
                                  ProteinsDirectoryFormat)

from rescript.ncbi import get_taxonomies, _default_ranks


def _get_assembly_descriptors(
        api_instance, assembly_levels, assembly_source, only_reference,
        page_size, taxon, tax_exact_match
):
    # all_acc_ids = []
    # all_tax_ids = []
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

        # all_acc_ids.extend([x.assembly_accession for x in genome_assembly])
        # all_tax_ids.extend([x.org.tax_id for x in genome_assembly])
        assembly_to_taxon.update(
            {x.assembly_accession: x.org.tax_id for x in genome_assembly}
        )

        if not next_page_token:
            break
    return assembly_to_taxon


def _extract_accession_ids(path: str):
    assembly_id = '_'.join(os.path.basename(path).split("_")[:2])
    seq = skbio.read(
        path, format='fasta', constructor=skbio.DNA, lowercase=True
    )
    accession_map = {s.metadata['id']: assembly_id for s in seq}
    return accession_map


def _fetch_and_extract_dataset(api_response):
    with tempfile.TemporaryDirectory() as tmp:
        result_path = os.path.join(tmp, 'datasets.zip')
        with open(result_path, 'wb') as f:
            f.write(api_response.data)

        with zipfile.ZipFile(result_path, 'r') as zipped:
            zipped.extractall(tmp)

        # find and move all the genome sequences
        genomes = DNAFASTAFormat()
        genome_seq_fps = glob.glob(
            os.path.join(tmp, 'ncbi_dataset', 'data', '*', '*_genomic.fna')
        )
        accession_to_assembly = {}
        with open(str(genomes), 'a') as fin:
            for f in genome_seq_fps:
                seq = skbio.read(
                    f, format='fasta', constructor=skbio.DNA, lowercase=True
                )
                skbio.io.write(seq, format='fasta', into=fin)
                accession_to_assembly.update(_extract_accession_ids(f))
        accession_to_assembly = pd.Series(
            accession_to_assembly, name='assembly_id'
        )

        # find and move all the gff files
        loci = LociDirectoryFormat()
        loci_fps = glob.glob(
            os.path.join(tmp, 'ncbi_dataset', 'data', '*', 'genomic.gff'))
        for f in loci_fps:
            _id = f.split('/')[-2]
            shutil.move(f, os.path.join(str(loci), f'{_id}_loci.gff'))

        # find and move all the protein translations
        proteins = ProteinsDirectoryFormat()
        protein_fps = glob.glob(
            os.path.join(tmp, 'ncbi_dataset', 'data', '*', 'protein.faa'))
        for f in protein_fps:
            _id = f.split('/')[-2]
            shutil.move(f, os.path.join(
                str(proteins), f'{_id}_proteins.fasta'))
    return genomes, loci, proteins, accession_to_assembly


def _fetch_taxonomy(all_acc_ids, all_tax_ids, accession_to_assembly):
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
    taxa = accession_to_assembly.replace(
        taxa.index, taxa.values, inplace=False
    )
    taxa.name = 'Taxon'
    taxa.index.name = 'Feature ID'
    return taxa


def get_ncbi_genomes(
        taxon: str,
        assembly_source: str = 'refseq',
        assembly_levels: List[str] = ['complete_genome'],
        only_reference: bool = True,
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

        api_response = api_instance.download_assembly_package(
            list(assembly_to_taxon.keys()),
            exclude_sequence=False,
            include_annotation_type=['PROT_FASTA', 'GENOME_GFF'],
            _preload_content=False
        )

        results = _fetch_and_extract_dataset(api_response)
        genomes, loci, proteins, accession_map = results
        taxa = _fetch_taxonomy(
            assembly_to_taxon.keys(),
            assembly_to_taxon.values(),
            accession_map
        )

    return genomes, loci, proteins, taxa
