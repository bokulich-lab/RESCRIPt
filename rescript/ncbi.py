# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import time
import logging
import json

import requests
from requests.exceptions import (
    HTTPError, ChunkedEncodingError, ConnectionError, ReadTimeout)
from xml.parsers.expat import ExpatError
from xmltodict import parse
from pandas import DataFrame
from q2_types.feature_data import DNAIterator, ProteinIterator
from skbio import DNA, Protein
from qiime2 import Metadata
from collections import OrderedDict
from joblib import Parallel, delayed
from multiprocessing import Manager

# This the tool-email pair for this tool that is registered with NCBI.
# If you change it, register a new tool-email pair with NCBI at
# eutilities@ncbi.nlm.nih.gov.
_entrez_params = dict(tool='qiime2-rescript', email='b.kaehler@adfa.edu.au')
# The entrez delay could be reduced to 0.1 with an NCBI API KEY but there is no
# point because this isn't the bottlekneck.
_entrez_delay = 0.334

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
        logging_level: str = None, n_jobs: int = 1
) -> (DNAIterator, DataFrame):
    if ranks is None:
        ranks = _default_ranks
    if query is None and accession_ids is None:
        raise ValueError('Query or accession_ids must be supplied')

    seqs, taxa = _get_ncbi_data(query, accession_ids, ranks, rank_propagation,
                                logging_level, n_jobs, 'nuccore')

    seqs = DNAIterator(DNA(v, metadata={'id': k}) for k, v in seqs.items())
    taxa = DataFrame(taxa, index=['Taxon']).T
    taxa.index.name = 'Feature ID'

    return seqs, taxa


def get_ncbi_data_protein(
        query: str = None, accession_ids: Metadata = None,
        ranks: list = None, rank_propagation: bool = True,
        logging_level: str = None, n_jobs: int = 1
        ) -> (ProteinIterator, DataFrame):
    if ranks is None:
        ranks = _default_ranks
    if query is None and accession_ids is None:
        raise ValueError('Query or accession_ids must be supplied')

    seqs, taxa = _get_ncbi_data(query, accession_ids, ranks, rank_propagation,
                                logging_level, n_jobs, 'protein')

    seqs = ProteinIterator(
        Protein(v, metadata={'id': k}) for k, v in seqs.items())
    taxa = DataFrame(taxa, index=['Taxon']).T
    taxa.index.name = 'Feature ID'

    return seqs, taxa


def _get_ncbi_data(query: str = None, accession_ids: Metadata = None,
                   ranks: list = None, rank_propagation: bool = True,
                   logging_level: str = None, n_jobs: int = 1,
                   db: str = 'nuccore'):
    manager = Manager()
    request_lock = manager.Lock()

    if query:
        seqs, taxids = get_data_for_query(
            query, logging_level, n_jobs, request_lock, _entrez_delay, db)

    if accession_ids:
        accs = accession_ids.get_ids()
        if query and seqs:
            accs = accs - seqs.keys()
            if accs:
                acc_seqs, acc_taxids = get_data_for_accs(
                    accs, logging_level, n_jobs, request_lock,
                    _entrez_delay, db)
                seqs.update(acc_seqs)
                taxids.update(acc_taxids)
        else:
            seqs, taxids = get_data_for_accs(
                accs, logging_level, n_jobs, request_lock, _entrez_delay, db)

    taxa, bad_accs = get_taxonomies(taxids, ranks, rank_propagation,
                                    logging_level, n_jobs, request_lock,
                                    _entrez_delay)
    for acc in bad_accs:
        del seqs[acc]

    return seqs, taxa


def _get_logger(logging_level):
    logger = logging.getLogger('rescript')
    if not hasattr(logger, 'isSetUp'):
        formatter = logging.Formatter(
            '%(levelname)s:%(asctime)-15s:%(processName)s:%(message)s')
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        if logging_level:
            logger.setLevel(logging_level)
        logger.isSetUp = True
    return logger


def _robustify(http_request, logger, *args):
    max_retries = 10
    backoff_factor = 1
    max_backoff = 120
    status_forcelist = [400, 429, 500, 502, 503, 504]
    exception_forcelist = (
        ChunkedEncodingError, ConnectionError, ExpatError, ReadTimeout)
    for retry in range(max_retries):
        try:
            return http_request(*args)
        except HTTPError as e:
            if e.response.status_code == 400:
                logger.debug(
                    'Request failed with code 400. This could be because all '
                    'of the requested accession ids are invalid, or it could '
                    'just be a temporary failure. Retrying.')
            elif e.response.status_code in status_forcelist:
                logger.debug('Request failed with code '
                             + str(e.response.status_code) + '. Retrying.')
            else:
                raise e
            last_exception = e
        except RuntimeError as e:
            if str(e) == 'bad record':
                logger.debug('Request failed. Retrying.')
            else:
                raise e
            last_exception = e
        except exception_forcelist as e:
            logger.debug('Request failed with exception\n' +
                         type(e).__name__ + ': ' + str(e) + '\nRetrying.')
            last_exception = e
        time.sleep(min(backoff_factor*2**retry, max_backoff))
    raise RuntimeError(
        'Maximum retries (10) exceeded for HTTP request. Persistent trouble '
        'downloading from NCBI. Last exception was\n' +
        type(last_exception).__name__ + ': ' + str(last_exception))


def _epost(params, ids, request_lock, logging_level, entrez_delay=0.334):
    assert len(ids) >= 1, "need at least one id"
    epost = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi'
    data = {'db': params['db'], 'id': ','.join(ids)}
    logger = _get_logger(logging_level)

    def request(params):
        request_lock.acquire()
        logger.debug('request lock acquired')
        time.sleep(entrez_delay)
        try:
            r = requests.post(
                epost, data=data, params=_entrez_params, timeout=10,
                stream=True)
        finally:
            request_lock.release()
            logger.debug('request lock released')
        r.raise_for_status()
        webenv = parse(r.content)['ePostResult']
        if 'ERROR' in webenv:
            if isinstance(webenv['ERROR'], list):
                for error in webenv['ERROR']:
                    logger.warning(error)
            else:
                logger.warning(webenv['ERROR'])
        if 'WebEnv' not in webenv:
            raise ValueError('No data for given ids')
        params = dict(params)
        params['WebEnv'] = webenv['WebEnv']
        params['query_key'] = webenv['QueryKey']
        return params
    return _robustify(request, logger, params)


def _esearch(params, logging_level, entrez_delay=0.334):
    esearch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

    def request(params):
        time.sleep(entrez_delay)
        r = requests.get(esearch, params=params, timeout=10)
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
    return _robustify(request, _get_logger(logging_level), params)


def _efetch_5000(params, request_lock, logging_level, entrez_delay=0.334):
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    params['retmax'] = 5000
    logger = _get_logger(logging_level)

    def request():
        request_lock.acquire()
        logger.debug('request lock acquired')
        time.sleep(entrez_delay)
        try:
            r = requests.get(efetch, params=params, timeout=10, stream=True)
        finally:
            request_lock.release()
            logger.debug('request lock released')
        r.raise_for_status()
        data = parse(r.content)
        data = list(data.values()).pop()
        data = list(data.values()).pop()
        # check that we got everything
        if isinstance(data, list):
            for rec in data:
                if not isinstance(rec, dict):
                    logger.debug('bad record:\n' + str(rec))
                    raise RuntimeError('bad record')
            else:
                return data
        elif isinstance(data, dict):
            return [data]
        else:
            logger.debug('bad record:\n' + str(data))
            raise RuntimeError('bad record')
    return _robustify(request, logger)


def _large_warning(logging_level):
    _get_logger(logging_level).warning(
        'This query could result in more than 100 requests to NCBI. If you '
        'are not running it on the weekend or between 9 pm and 5 am Eastern '
        'Time weekdays, it may result in NCBI blocking your IP address. See '
        'https://www.ncbi.nlm.nih.gov/home/about/policies/ for details.')


def _ungotten_ids(ids, data):
    error = ("Partial download. Expected " + str(len(ids)) +
             ' records, but got ' + str(len(data)) + '.')
    id_keys = ['TSeq_accver', 'TaxId']
    gotten = set()
    for record in data:
        for id_key in id_keys:
            if id_key in record:
                gotten.add(record[id_key])
                break
    ungotten = list(set(ids) - gotten)
    if len(ungotten) > 10:
        error += ('\nMore than 10 ids were missing. Ten were: ')
    else:
        error += '\nThe following ids were missing: '
    error += ', '.join(ungotten[:10]) + '.\n'
    return error


def _get_id_chunk(
        ids_chunk, params, request_lock, logging_level, raise_on_partial,
        entrez_delay):
    chunk_params = _epost(
        params, ids_chunk, request_lock, logging_level, entrez_delay)
    data_chunk = _efetch_5000(
        chunk_params, request_lock, logging_level, entrez_delay)
    logger = _get_logger(logging_level)
    if len(ids_chunk) != len(data_chunk):
        error = _ungotten_ids(ids_chunk, data_chunk)
        if raise_on_partial:
            raise RuntimeError(error)
        else:
            logger.warning(error)
    logger.info('got ' + str(len(data_chunk)) + ' records')
    return data_chunk


def _get_for_ids(params, ids, logging_level, n_jobs, request_lock,
                 raise_on_partial, entrez_delay=0.334):
    logger = _get_logger(logging_level)
    logger.info('Downloading ' + str(len(ids)) + ' records')
    ids = list(ids)
    parallel = Parallel(n_jobs=n_jobs, backend='loky')
    chunky = parallel(delayed(_get_id_chunk)(ids[chunk:chunk+5000], params,
                                             request_lock, logging_level,
                                             raise_on_partial, entrez_delay)
                      for chunk in range(0, len(ids), 5000))
    data = [chunk for chunks in chunky for chunk in chunks]

    return data


def get_data_for_accs(accs, logging_level, n_jobs, request_lock,
                      entrez_delay=0.334, db='nuccore'):
    if len(accs) > 125000:
        _large_warning(logging_level)
    params = dict(
        db=db, rettype='fasta', retmode='xml', **_entrez_params
    )
    records = _get_for_ids(params, accs, logging_level, n_jobs, request_lock,
                           True, entrez_delay)
    seqs = {}
    taxids = {}
    for rec in records:
        seqs[rec['TSeq_accver']] = rec['TSeq_sequence']
        taxids[rec['TSeq_accver']] = rec['TSeq_taxid']
    return seqs, taxids


def _get_query_chunk(
        chunk, params, entrez_delay, expected_num_records, request_lock,
        logging_level):
    params['retstart'] = chunk
    data_chunk = _efetch_5000(params, request_lock, logging_level,
                              entrez_delay)
    if len(data_chunk) != min(5000, expected_num_records-chunk):
        logger = _get_logger(logging_level)
        logger.warn(
            'Expected ' + str(min(5000, expected_num_records-chunk)) +
            ' sequences in this chunk, but got ' + str(len(data_chunk)) +
            '. I do not know why, or which sequences are missing.')
    logger = _get_logger(logging_level)
    logger.info('got ' + str(len(data_chunk)) + ' sequences')
    return data_chunk


def get_data_for_query(query, logging_level, n_jobs, request_lock,
                       entrez_delay=0.334, db='nuccore'):
    params = dict(
        db=db, term=query, usehistory='y', retmax=0, **_entrez_params
    )
    params, expected_num_records = _esearch(params, logging_level,
                                            entrez_delay)
    if expected_num_records > 166666:
        _large_warning(logging_level)
    records = []
    logger = _get_logger(logging_level)
    logger.info('Downloading ' + str(expected_num_records) + ' sequences')

    parallel = Parallel(n_jobs=n_jobs, backend='loky')
    chunky = parallel(delayed(_get_query_chunk)(chunk, params, entrez_delay,
                                                expected_num_records,
                                                request_lock, logging_level)
                      for chunk in range(0, expected_num_records, 5000))
    records = [chunk for chunks in chunky for chunk in chunks]

    seqs = {}
    taxids = {}
    for rec in records:
        try:
            seqs[rec['TSeq_accver']] = rec['TSeq_sequence']
            taxids[rec['TSeq_accver']] = rec['TSeq_taxid']
        except KeyError as e:
            if str(e) == "'TSeq_accver'":
                seqs[rec['TSeq_sid']] = rec['TSeq_sequence']
                taxids[rec['TSeq_sid']] = rec['TSeq_taxid']
                logger.warning(
                    'Using ' + rec['TSeq_sid'] + ' as a sequence identifier, '
                    'because it did not come down with an accession version.')
            else:
                raise
    return seqs, taxids


def get_taxonomies(
        taxids, ranks, rank_propagation, logging_level, n_jobs, request_lock,
        entrez_delay=0.334):
    # download the taxonomies
    params = dict(db='taxonomy', **_entrez_params)
    ids = set(map(str, taxids.values()))
    records = _get_for_ids(params, ids, logging_level, n_jobs, request_lock,
                           False, entrez_delay)

    # parse the taxonomies
    taxa = {}
    for rec in records:
        try:
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
            if 'AkaTaxIds' in rec:
                for akaTaxId in rec['AkaTaxIds']:
                    taxa[rec['AkaTaxIds']['TaxId']] = taxonomy
        except (KeyError, TypeError) as e:  # these are those we've seen so far
            logger = _get_logger(logging_level)
            logger.warning(
                'Got exception\n' + type(e).__name__ + ': ' + str(e) +
                '\nfrom taxon\n' + json.dumps(rec, indent=1) +
                '\nSkipping this taxon.')

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
        logger = _get_logger(logging_level)
        logger.warning(
            'The following accessions were deleted from the sequence database '
            'because there was a problem with their taxonomies: ' +
            ', '.join(missing_accs) + '.\nThe problematic taxids were: ' +
            ', '.join(missing_taxids) + '.')
    return tax_strings, missing_accs
