# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import time
import warnings

import requests
from xmltodict import parse
from pandas import DataFrame
from q2_types.feature_data import DNAIterator
from skbio import DNA
from qiime2 import Metadata
from collections import OrderedDict

_default_ranks = [
    'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
]
_default_allowed_ranks = [
    'domain', 'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum',
    'subphylum', 'infraphylum', 'superclass', 'class', 'subclass',
    'infraclass', 'superorder', 'order', 'suborder', 'superfamily', 'family',
    'subfamily', 'genus', 'species'
]


def get_ncbi_data(
        query: str = None, accession_ids: Metadata = None,
        ranks: list = None, allowed_ranks: list = None,
        disable_rank_propagation: bool = False,
        entrez_delay: float = 0.334) -> (DNAIterator, DataFrame):
    if query is None and accession_ids is None:
        raise ValueError('Query or accession_ids must be supplied')
    if ranks is None:
        ranks = _default_ranks
    if allowed_ranks is None:
        allowed_ranks = _default_allowed_ranks
    if not disable_rank_propagation and not (set(ranks) <= set(allowed_ranks)):
        raise ValueError('Ranks must be a subset of allowed ranks')

    if query:
        seqs, taxids = get_nuc_for_query(query, entrez_delay)

    if accession_ids:
        accs = accession_ids.get_ids()
        if query and seqs:
            accs = accs - seqs.keys()
            if accs:
                acc_seqs, acc_taxids = get_nuc_for_accs(accs, entrez_delay)
                seqs.update(acc_seqs)
                taxids.update(acc_taxids)
        else:
            seqs, taxids = get_nuc_for_accs(accs, entrez_delay)

    taxa = get_taxonomies(
        taxids, ranks, allowed_ranks, disable_rank_propagation, entrez_delay)

    seqs = DNAIterator(DNA(v, metadata={'id': k}) for k, v in seqs.items())
    taxa = DataFrame(taxa, index=['Taxon']).T
    taxa.index.name = 'Feature ID'

    return seqs, taxa


def _get(params, ids=None, entrez_delay=0.):
    if ids:
        assert len(ids) >= 1, "need at least one id"
        epost = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
        data = {'db': params['db'], 'id': ','.join(ids)}
        time.sleep(entrez_delay)
        r = requests.post(epost, data=data)
        r.raise_for_status()
        webenv = parse(r.content)['ePostResult']
        params['WebEnv'] = webenv['WebEnv']
        params['query_key'] = webenv['QueryKey']
        expected_num_records = len(ids)
    else:
        esearch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
        time.sleep(entrez_delay)
        r = requests.get(esearch, params=params)
        r.raise_for_status()
        webenv = parse(r.content)['eSearchResult']
        params = dict(
            db='nuccore', rettype='fasta', retmode='xml',
            WebEnv=webenv['WebEnv'], query_key=webenv['QueryKey']
        )
        expected_num_records = int(webenv['Count'])

    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    # Need to better handle the case of invalid ids
    # Need to better handle intermittent fault where we
    #  don't get the whole data set back
    data = []
    while len(data) < expected_num_records:
        params['retstart'] = len(data)
        time.sleep(entrez_delay)
        r = requests.get(efetch, params=params)
        r.raise_for_status()
        chunk = parse(r.content)
        chunk = list(chunk.values()).pop()
        chunk = list(chunk.values()).pop()
        if expected_num_records == 1:
            return [chunk]
        data.extend(chunk)
    return data


def get_nuc_for_accs(accs, entrez_delay=0.):
    params = dict(
        db='nuccore', rettype='fasta', retmode='xml'
    )
    records = _get(params, accs, entrez_delay)
    seqs = {}
    taxids = {}
    for rec in records:
        seqs[rec['TSeq_accver']] = rec['TSeq_sequence']
        taxids[rec['TSeq_accver']] = rec['TSeq_taxid']
    return seqs, taxids


def get_nuc_for_query(query, entrez_delay=0.):
    params = dict(
        db='nuccore', term=query, usehistory='y', retmax=0
    )
    records = _get(params, entrez_delay=entrez_delay)
    seqs = {}
    taxids = {}
    for rec in records:
        seqs[rec['TSeq_accver']] = rec['TSeq_sequence']
        taxids[rec['TSeq_accver']] = rec['TSeq_taxid']
    return seqs, taxids


def get_taxonomies(taxids, ranks=None, allowed_ranks=None,
                   disable_rank_propagation=False, entrez_delay=0.):
    # download the taxonomies
    params = dict(db='taxonomy')
    ids = list(map(str, taxids.values()))
    records = _get(params, ids, entrez_delay)
    taxa = {}

    # parse the taxonomies
    if disable_rank_propagation:
        allowed_ranks = ranks
    for rec in records:
        taxonomy = OrderedDict([('NCBI_Division', rec['Division'])])
        taxonomy.update((r, None) for r in allowed_ranks)
        for rank in rec['LineageEx']['Taxon']:
            if rank['Rank'] in allowed_ranks:
                taxonomy[rank['Rank']] = rank['ScientificName']
        species = rec['ScientificName']
        if 'genus' in allowed_ranks:
            if taxonomy['genus']:
                if species.startswith(taxonomy['genus'] + ' '):
                    species = species[len(taxonomy['genus']) + 1:]
            elif ' ' in species:
                genus, species = species.split(' ', 1)
                taxonomy['genus'] = genus
        if 'species' in allowed_ranks:
            taxonomy['species'] = species
        last_label = taxonomy['NCBI_Division']
        for rank in taxonomy:
            if taxonomy[rank] is None:
                if disable_rank_propagation:
                    taxonomy[rank] = ''
                else:
                    taxonomy[rank] = last_label
            last_label = taxonomy[rank]
        taxa[rec['TaxId']] = taxonomy

    # return the taxonomies
    missing_accs = []
    missing_taxids = {}
    tax_strings = {}
    for acc, taxid in taxids.items():
        if taxid in taxa:
            taxonomy = taxa[taxid]
            ts = '; '.join([r[0] + '__' + taxonomy[r] for r in ranks])
            tax_strings[acc] = ts
        else:
            missing_accs.append(acc)
            missing_taxids.add(taxid)
    if missing_accs:
        warnings.warn('The following accessions did not have valid taxids: ' +
                      ', '.join(missing_accs) + '. The bad taxids were: ' +
                      ', '.join(missing_taxids), UserWarning)
    return tax_strings
