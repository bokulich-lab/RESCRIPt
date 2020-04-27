#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2019--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import scipy
import qiime2 as q2

from collections import Counter
from itertools import zip_longest

from ._utilities import _taxon_to_list


def evaluate_taxonomy(ctx,
                      taxonomies,
                      rank_handle=None):
    summaries = []
    for n, taxonomy in enumerate(taxonomies, 1):
        summary = _evaluate_taxonomy(taxonomy.view(pd.Series), rank_handle)
        summary['Dataset'] = str(n)
        summaries.append(summary)
    results = pd.concat(summaries).reset_index()
    results.index = pd.Index(
        [str(i) for i in range(1, len(results.index) + 1)], name='id')
    results = q2.Metadata(results)
    volatility = ctx.get_action('longitudinal', 'volatility')
    plots, = volatility(metadata=results,
                        state_column='Level',
                        default_group_column='Dataset',
                        default_metric='Taxonomic Entropy')
    return plots


def _evaluate_taxonomy(taxonomy, rank_handle):
    if rank_handle is None:
        rank_handle = ""
    # Count number of unique taxa and unclassifieds at each level
    summary = summarize_taxonomic_depth(taxonomy, rank_handle=rank_handle)
    # Measure taxonomic entropy at each level
    entropy = _taxonomic_entropy(taxonomy, rank_handle=rank_handle)
    entropy = entropy.merge(summary, left_index=True, right_index=True)
    entropy.index.name = 'Level'
    return entropy


def summarize_taxonomic_depth(taxonomy, rank_handle):
    # measure number of levels in each taxonomy
    depths = _taxonomic_depth(taxonomy, rank_handle=rank_handle)
    # count number + proportion of taxonomies per depth
    depths = depths.value_counts()
    total = remaining = depths.sum()
    proportions = depths / total
    depths = pd.concat([depths, proportions], axis=1)
    depths.columns = ['Number of Features Terminating at Depth',
                      'Proportion of Features Terminating at Depth']
    # iterate over depths up to max depth
    # fill in empty depths (i.e., indicate no termination)
    # calculate cumulative sum of (un)classifieds per depth
    depths_idx = depths.index
    classified = []
    unclassified = []
    for d in range(1, depths_idx.max() + 1):
        if d not in depths_idx:
            depths.loc[d] = [0, 0]
        classified.append(remaining)
        unclassified.append(total - remaining)
        remaining -= depths.loc[d, 'Number of Features Terminating at Depth']
    depths = depths.sort_index()
    classified_prop = [c / total for c in classified]
    unclassified_prop = [u / total for u in unclassified]
    depths['Number of Features Classified at Depth'] = classified
    depths['Proportion of Features Classified at Depth'] = classified_prop
    depths['Number of Features Unclassified at Depth'] = unclassified
    depths['Proportion of Features Unclassified at Depth'] = unclassified_prop
    depths.index.name = 'Level'
    return depths


def _taxonomic_entropy(taxonomy, rank_handle):
    # convert each taxonomy to a list
    taxa_lists = [_taxon_to_list(t, rank_handle=rank_handle)
                  for t in taxonomy.values]
    # taxonomic labels should accumulate at each rank (e.g., to avoid
    # counting identical species names as duplicates when genus is unique)
    taxa_lists = [[';'.join(t[:i]) for i in range(1, len(t) + 1)]
                  for t in taxa_lists]
    # coalesce taxonomies at each rank
    ranks = [i for i in zip_longest(*taxa_lists)]
    # Count unique values at each rank, excluding 'None'
    unique_counts = [[v for k, v in Counter(r).items() if k not in [None, '']]
                     for r in ranks]
    # measure entropy at each rank, excluding 'None'
    entropy = {n: [len(r), scipy.stats.entropy(r)]
               for n, r in enumerate(unique_counts, 1)}
    return pd.DataFrame(entropy,
                        index=['Unique Labels', 'Taxonomic Entropy']).T


def _taxonomic_depth(taxonomy, rank_handle):
    return taxonomy.apply(lambda x: len(set(_taxon_to_list(
        x, rank_handle=rank_handle)) - {None, ''}))
