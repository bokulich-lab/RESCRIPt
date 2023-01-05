# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import pandas as pd
from q2_types.feature_data import DNAFASTAFormat, DNAIterator

from ._utilities import run_command

ERROR_FILTER_OPTIONS = (
    'No filters were applied. One or more of the following filter settings '
    'are required: ')


def filter_seqs_length(sequences: DNAFASTAFormat,
                       global_min: int = None,
                       global_max: int = None,
                       threads: int = 1) -> (DNAFASTAFormat, DNAFASTAFormat):
    # Validate filtering options
    if global_min is global_max is None:
        raise ValueError(ERROR_FILTER_OPTIONS + 'global_min, global_max.')
    # Filter with vsearch (for global filtering alone, this should be quick)
    result = DNAFASTAFormat()
    failures = DNAFASTAFormat()
    cmd = ['vsearch', '--fastx_filter', str(sequences), '--fastaout',
           str(result), '--fastaout_discarded', str(failures),
           '--threads', str(threads)]
    if global_min is not None:
        cmd.extend(['--fastq_minlen', str(global_min)])
    if global_max is not None:
        cmd.extend(['--fastq_maxlen', str(global_max)])
    run_command(cmd)
    return result, failures


def filter_seqs_length_by_taxon(sequences: DNAFASTAFormat,
                                taxonomy: pd.Series,
                                labels: str,
                                min_lens: int = None,
                                max_lens: int = None,
                                global_min: int = None,
                                global_max: int = None
                                ) -> (DNAFASTAFormat, DNAFASTAFormat):
    # Validate filtering options
    if min_lens is max_lens is None:
        raise ValueError(ERROR_FILTER_OPTIONS + 'min_lens, max_lens.')

    # validate that all seqIDs are present in taxonomy
    # Note we view as DNAIterator to take a first pass (should take a few
    # seconds) as initial validation before performing length filtering.
    seq_ids = {i.metadata['id'] for i in sequences.view(DNAIterator)}
    _index_is_superset(seq_ids, set(taxonomy.index))

    # set filter options
    mins = maxs = None
    if min_lens is not None:
        if len(labels) != len(min_lens):
            raise ValueError(
                'labels and min_lens must contain the same number of elements')
        else:
            mins = {k: v for k, v in zip(labels, min_lens)}

    if max_lens is not None:
        if len(labels) != len(max_lens):
            raise ValueError(
                'labels and max_lens must contain the same number of elements')
        else:
            maxs = {k: v for k, v in zip(labels, max_lens)}

    # Stream seqs, apply filter(s)
    result = DNAFASTAFormat()
    failures = DNAFASTAFormat()
    with result.open() as out_fasta, failures.open() as out_failed:
        for seq in sequences.view(DNAIterator):
            # taxon is required, we always use taxon-based filtering
            # grab taxon affiliation for seq
            taxon = taxonomy[seq.metadata['id']]
            # search taxon for filter terms
            # NOTE: we find all matching search terms and pass them all to
            # _seq_length_within_range below; that function determines and
            # applies the most stringent matching length thresholds.
            taxahits = [t for t in labels if t in taxon]
            # if there are no taxahits or global filters, just write out
            if not any(taxahits) and global_min is global_max is None:
                seq.write(out_fasta)
            # if there are taxahits or global filters, always check length
            elif _seq_length_within_range(seq, taxahits, mins, maxs,
                                          global_min, global_max):
                seq.write(out_fasta)
            else:
                seq.write(out_failed)
    return result, failures


def filter_taxa(taxonomy: pd.Series, ids_to_keep: qiime2.Metadata = None,
                include: str = None, exclude: str = None) -> pd.Series:
    if include is exclude is ids_to_keep is None:
        raise ValueError('No filtering criteria were applied!')

    ids = taxonomy.index
    print('Input features: ' + str(len(ids)))
    # some initial validation
    if ids_to_keep:
        ids_to_keep = ids_to_keep.ids
        _index_is_superset(set(ids_to_keep), ids)

    if include:
        ids = ids[taxonomy.str.contains('|'.join(include))]

    if exclude:
        bad = taxonomy.index[taxonomy.str.contains('|'.join(exclude))]
        ids = ids.difference(bad)

    # if not using exclude or include, we only want explicit ids_to_keep
    if include is exclude is None:
        ids = ids_to_keep
    # otherwise add back ids_to_keep for explicit inclusion after filtering
    elif ids_to_keep:
        ids = ids.union(ids_to_keep)

    filtered_ids = len(ids)
    print('Output features: ' + str(filtered_ids))

    if filtered_ids == 0:
        raise ValueError("All features were filtered, resulting in an "
                         "empty collection of taxonomies.")

    taxonomy = taxonomy.reindex(ids)
    taxonomy.index.name = 'Feature ID'

    return taxonomy


def _seq_length_within_range(sequence, taxahits, mins, maxs, global_min,
                             global_max):
    '''
    Check is sequence length is within acceptable range.

    sequence: str
        (DNA) sequence
    taxahits: list of str
        taxonomic labels matching filter query.
    mins: dict
        Dict of {str: int} specifying minimum length threshold associated with
        each taxonomic label. Labels must be superset of taxahits.
    maxs: dict
        Dict of {str: int} specifying maximum length threshold associated with
        each taxonomic label. Labels must be superset of taxahits.
    global_min: int
        Minimum length threshold to apply globally.
    global_max: int
        Maximum length threshold to apply globally.

    Return bool
        True if seq within range, False if out of range.
    '''
    seqlen = len(sequence)
    minlen = []
    maxlen = []
    # find all applicable min/max filter thresholds.
    # If a sequence matches multiple taxonomic terms in the search, we find
    # the most stringent threshold: the longest minimum length and/or shortest
    # maximum length for filtering.
    if mins is not None:
        minlen.extend([mins[k] for k in taxahits])
    if maxs is not None:
        maxlen.extend([maxs[k] for k in taxahits])
    if global_min is not None:
        minlen.append(global_min)
    if global_max is not None:
        maxlen.append(global_max)
    # get most stringent thresholds for taxon: min max and max min
    if any(maxlen):
        maxlen = min(maxlen)
    else:
        maxlen = seqlen
    if any(minlen):
        minlen = max(minlen)
    else:
        minlen = 0
    # Finally check seqlen
    return maxlen >= seqlen >= minlen


def _index_is_superset(index1, index2):
    '''
    Validate that index1 is a subset of index2.
    idx1: set
    idx2: set
    '''
    diff = index1.difference(index2)
    if len(diff) > 0:
        raise ValueError('The following IDs are missing from '
                         'the taxonomy: ' + ', '.join(sorted(diff)))
