# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import time
import logging

import requests
from requests.exceptions import (
    HTTPError, ChunkedEncodingError, ConnectionError, ReadTimeout)
from xml.parsers.expat import ExpatError
from xmltodict import parse
from pandas import DataFrame
from q2_types.feature_data import DNAIterator
from skbio import DNA
from qiime2 import Metadata
from collections import OrderedDict

_entrez_params = dict(tool='qiime2-rescript', email='b.kaehler@adfa.edu.au')

_default_ranks = [
    'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'
]

_allowed_ranks = OrderedDict([
    ('domain', 'd__'),
    ('superkingdom', 'sk__'),
    ('kingdom', 'k__'),
    ('subkingdom', 'ks__'),
    ('superphylum', 'sp__'),
    ('phylum', 'p__'),
    ('subphylum', 'ps__'),
    ('infraphylum', 'pi__'),
    ('superclass', 'sc__'),
    ('class', 'c__'),
    ('subclass', 'cs__'),
    ('infraclass', 'ci__'),
    ('cohort', 'co__'),
    ('superorder', 'so__'),
    ('order', 'o__'),
    ('suborder', 'os__'),
    ('infraorder', 'oi__'),
    ('parvorder', 'op__'),
    ('superfamily', 'sf__'),
    ('family', 'f__'),
    ('subfamily', 'fs__'),
    ('tribe', 't__'),
    ('subtribe', 'ts__'),
    ('genus', 'g__'),
    ('subgenus', 'gs__'),
    ('species group', 'ss__'),
    ('species subgroup', 'sgs__'),
    ('species', 's__'),
    ('subspecies', 'ssb__'),
    ('forma', 'f__')
])


def get_ncbi_data(
        query: str = None, accession_ids: Metadata = None,
        ranks: list = None, rank_propagation: bool = True,
        entrez_delay: float = 0.334, logging_level: str = None
        ) -> (DNAIterator, DataFrame):
    if logging_level is not None:
        logging.basicConfig(level=logging_level)

    if ranks is None:
        ranks = _default_ranks

    if query is None and accession_ids is None:
        raise ValueError('Query or accession_ids must be supplied')
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
        taxids, ranks, rank_propagation, entrez_delay)

    seqs = DNAIterator(DNA(v, metadata={'id': k}) for k, v in seqs.items())
    taxa = DataFrame(taxa, index=['Taxon']).T
    taxa.index.name = 'Feature ID'

    return seqs, taxa


def _robustify(http_request, *args):
    max_retries = 10
    backoff_factor = 1
    max_backoff = 120
    status_forcelist = [429, 500, 502, 503, 504]
    exception_forcelist = (
        ChunkedEncodingError, ConnectionError, ExpatError, ReadTimeout)
    for retry in range(max_retries):
        try:
            return http_request(*args)
        except HTTPError as e:
            if e.response.status_code == 400:  # because of missing ids
                return []
            if e.response.status_code not in status_forcelist:
                raise e
            logging.debug('Request failed with code '
                          + str(e.response.status_code) + '. Retrying.')
            last_exception = e
        except RuntimeError as e:
            if str(e) == 'bad record':
                logging.debug('Request failed. Retrying.')
            else:
                raise e
            last_exception = e
        except exception_forcelist as e:
            logging.debug('Request failed with exception:\n' +
                          str(e) + '\nRetrying.')
            last_exception = e
        time.sleep(min(backoff_factor*2**retry, max_backoff))
    raise RuntimeError(
        'Maximum retries (10) exceeded for HTTP request. Persistent trouble '
        'downloading from NCBI. Last exception was:\n' + str(last_exception))


def _epost(params, ids, entrez_delay=0.):
    assert len(ids) >= 1, "need at least one id"
    epost = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
    data = {'db': params['db'], 'id': ','.join(ids)}

    def request(params):
        time.sleep(entrez_delay)
        r = requests.post(
            epost, data=data, params=_entrez_params, timeout=5)
        r.raise_for_status()
        webenv = parse(r.content)['ePostResult']
        if 'ERROR' in webenv:
            if isinstance(webenv['ERROR'], list):
                for error in webenv['ERROR']:
                    logging.warning(error)
            else:
                logging.warning(webenv['ERROR'])
        if 'WebEnv' not in webenv:
            raise ValueError('No data for given ids')
        params = dict(params)
        params['WebEnv'] = webenv['WebEnv']
        params['query_key'] = webenv['QueryKey']
        return params
    return _robustify(request, params)


def _esearch(params, entrez_delay=0.):
    esearch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

    def request(params):
        time.sleep(entrez_delay)
        r = requests.get(esearch, params=params, timeout=5)
        r.raise_for_status()
        webenv = parse(r.content)['eSearchResult']
        if 'WebEnv' not in webenv:
            raise ValueError('No sequences for given query')
        params = dict(
            db='nuccore', rettype='fasta', retmode='xml',
            WebEnv=webenv['WebEnv'],
            query_key=webenv['QueryKey'], **_entrez_params
        )
        expected_num_records = int(webenv['Count'])
        return params, expected_num_records
    return _robustify(request, params)


def _efetch_5000(params, entrez_delay=0.):
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params['retmax'] = 5000

    def request():
        time.sleep(entrez_delay)
        r = requests.get(efetch, params=params, timeout=5)
        r.raise_for_status()
        data = parse(r.content)
        data = list(data.values()).pop()
        data = list(data.values()).pop()
        # check that we got everything
        if isinstance(data, list):
            for rec in data:
                if not isinstance(rec, dict):
                    logging.debug('bad record:\n' + str(rec))
                    raise RuntimeError('bad record')
            else:
                return data
        elif isinstance(data, dict):
            return [data]
        else:
            logging.debug('bad record:\n' + str(data))
            raise RuntimeError('bad record')
    return _robustify(request)


def _large_warning():
    logging.warning(
        'This query could result in more than 100 requests to NCBI. If you '
        'are not running it on the weekend or between 9 pm and 5 am Eastern '
        'Time weekdays, it may result in NCBI blocking your IP address. See '
        'https://www.ncbi.nlm.nih.gov/home/about/policies/ for details.')


def _ungotten_ids(ids, data):
    ungotten = []
    for _id in ids:
        ungotten_id = True
        for datum in data:
            if _id in str(datum):
                ungotten_id = False
                break
        if ungotten_id:
            ungotten.append(_id)
    if len(ungotten) > 10:
        error = ('\nMore than 10 records were missing. The first 10 were ')
    else:
        error = '\nThe following records were missing '
    error += ', '.join(ungotten[:10]) + '.\n'
    return error


def _get_for_ids(params, ids, entrez_delay=0.):
    logging.info('Downloading ' + str(len(ids)) + ' records')
    ids = list(ids)
    data = []
    for chunk in range(0, len(ids), 5000):
        ids_chunk = ids[chunk:chunk+5000]
        chunk_params = _epost(params, ids_chunk, entrez_delay)
        data_chunk = _efetch_5000(chunk_params, entrez_delay)
        if len(ids_chunk) != len(data_chunk):
            error = "Download did not finish."
            error += _ungotten_ids(ids_chunk, data_chunk)
            raise RuntimeError(error)
        data.extend(data_chunk)
        logging.info('got ' + str(len(data)) + ' records')
    return data


def get_nuc_for_accs(accs, entrez_delay=0.):
    if len(accs) > 125000:
        _large_warning()
    params = dict(
        db='nuccore', rettype='fasta', retmode='xml', **_entrez_params
    )
    records = _get_for_ids(params, accs, entrez_delay)
    seqs = {}
    taxids = {}
    for rec in records:
        seqs[rec['TSeq_accver']] = rec['TSeq_sequence']
        taxids[rec['TSeq_accver']] = rec['TSeq_taxid']
    return seqs, taxids


def get_nuc_for_query(query, entrez_delay=0.):
    params = dict(
        db='nuccore', term=query, usehistory='y', retmax=0, **_entrez_params
    )
    params, expected_num_records = _esearch(params, entrez_delay)
    if expected_num_records > 166666:
        _large_warning()
    records = []
    logging.info('Downloading ' + str(expected_num_records) + ' sequences')
    for chunk in range(0, expected_num_records, 5000):
        params['retstart'] = chunk
        data_chunk = _efetch_5000(params, entrez_delay)
        records.extend(data_chunk)
        logging.info('got ' + str(len(records)) + ' sequences')
    if len(records) != expected_num_records:
        raise RuntimeError(
            'Download did not finish. Expected ' + str(expected_num_records) +
            ' sequences, but only got ' + str(len(records)))
    seqs = {}
    taxids = {}
    for rec in records:
        seqs[rec['TSeq_accver']] = rec['TSeq_sequence']
        taxids[rec['TSeq_accver']] = rec['TSeq_taxid']
    return seqs, taxids


def get_taxonomies(
        taxids, ranks=None, rank_propagation=False, entrez_delay=0.):
    # download the taxonomies
    params = dict(db='taxonomy', **_entrez_params)
    ids = set(map(str, taxids.values()))
    records = _get_for_ids(params, ids, entrez_delay)
    taxa = {}

    # parse the taxonomies
    for rec in records:
        if rank_propagation:
            taxonomy = OrderedDict([('NCBI_Division', rec['Division'])])
            taxonomy.update((r, None) for r in _allowed_ranks)
            for rank in rec['LineageEx']['Taxon']:
                if rank['Rank'] in _allowed_ranks:
                    taxonomy[rank['Rank']] = rank['ScientificName']
            species = rec['ScientificName']
            if taxonomy['genus']:
                if species.startswith(taxonomy['genus'] + ' '):
                    species = species[len(taxonomy['genus']) + 1:]
            elif ' ' in species:
                genus, species = species.split(' ', 1)
                taxonomy['genus'] = genus
            taxonomy['species'] = species
            last_label = taxonomy['NCBI_Division']
            for rank in taxonomy:
                if taxonomy[rank] is None:
                    taxonomy[rank] = last_label
                last_label = taxonomy[rank]
        else:
            taxonomy = {r: '' for r in ranks}
            if 'domain' in ranks:
                taxonomy['domain'] = rec['Division']
            elif 'kingdom' in ranks:
                taxonomy['kingdom'] = rec['Division']
            for rank in rec['LineageEx']['Taxon']:
                taxonomy[rank['Rank']] = rank['ScientificName']
            species = rec['ScientificName']
            # if we care about genus and genus is in the species label and
            # we don't already know genus, split it out
            if 'genus' in ranks:
                if taxonomy['genus']:
                    if species.startswith(taxonomy['genus'] + ' '):
                        species = species[len(taxonomy['genus']) + 1:]
                elif ' ' in species:
                    genus, species = species.split(' ', 1)
                    taxonomy['genus'] = genus
            taxonomy['species'] = species

        taxa[rec['TaxId']] = taxonomy

    # return the taxonomies
    missing_accs = []
    missing_taxids = set()
    tax_strings = {}
    for acc, taxid in taxids.items():
        if taxid in taxa:
            taxonomy = taxa[taxid]
            ts = '; '.join([_allowed_ranks[r] + taxonomy[r] for r in ranks])
            tax_strings[acc] = ts
        else:
            missing_accs.append(acc)
            missing_taxids.add(taxid)
    if missing_accs:
        logging.warning(
            'The following accessions did not have valid taxids: ' +
            ', '.join(missing_accs) + '. The bad taxids were: ' +
            ', '.join(missing_taxids), UserWarning)
    return tax_strings
