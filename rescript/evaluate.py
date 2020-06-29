# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import scipy
import warnings
import qiime2 as q2

from collections import Counter
from itertools import zip_longest

from ._utilities import _taxon_to_list


def evaluate_taxonomy(ctx,
                      taxonomies,
                      labels=None,
                      rank_handle_regex=None):
    labels = _process_labels(labels, taxonomies)
    summaries = []
    for name, taxonomy in zip(labels, taxonomies):
        summary = _evaluate_taxonomy(
            taxonomy.view(pd.Series), rank_handle_regex)
        # note: what if n equals an existing name? do we want to error out or
        # just assume a user knows what they are doing? For now I am just
        # adding this info to the warning to make this behavior transparent.
        summary['Dataset'] = str(name)
        summaries.append(summary)
    results = pd.concat(summaries).reset_index()
    # convert index to strings
    results.index = pd.Index(
        [str(i) for i in range(1, len(results.index) + 1)], name='id')
    results = q2.Metadata(results)
    volatility = ctx.get_action('longitudinal', 'volatility')
    plots, = volatility(metadata=results,
                        state_column='Level',
                        default_group_column='Dataset',
                        default_metric='Taxonomic Entropy')
    return plots


def _process_labels(labels, taxonomies):
    if labels is None:
        labels = []
    n_labels, n_taxonomies = len(labels), len(taxonomies)
    # warn if fewer labels than taxonomies
    if n_labels < n_taxonomies:
        _warn_uneven_length()
        labels = [name if name is not None else n for n, (name, t)
                  in enumerate(zip_longest(labels, taxonomies), 1)]
    # trim labels to number of Taxonomies
    elif n_labels > n_taxonomies:
        labels = labels[:n_taxonomies]
    return labels


def _warn_uneven_length():
    msg = ("The lists of input taxonomies and labels are different lengths. "
           "Additional taxonomies will be labeled numerically by their order "
           "in the inputs. Note that if these numbers equal any existing "
           "labels, those data will be grouped in the visualization.")
    warnings.warn(msg, UserWarning)


def _evaluate_taxonomy(taxonomy, rank_handle_regex):
    if rank_handle_regex is None:
        rank_handle_regex = ""
    # Count number of unique taxa and unclassifieds at each level
    summary = summarize_taxonomic_depth(
        taxonomy, rank_handle_regex=rank_handle_regex)
    # Measure taxonomic entropy at each level
    max_depth = summary.index.max()
    entropy = _taxonomic_entropy(
        taxonomy, rank_handle_regex=rank_handle_regex, max_depth=max_depth)
    entropy = entropy.merge(summary, left_index=True, right_index=True)
    entropy.index.name = 'Level'
    return entropy


def summarize_taxonomic_depth(taxonomy, rank_handle_regex):
    # measure number of levels in each taxonomy
    depths = _taxonomic_depth(taxonomy, rank_handle_regex=rank_handle_regex)
    # count number + proportion of taxonomies per depth
    depths = depths.value_counts()
    total = remaining = len(taxonomy)
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
    for d in range(0, depths_idx.max() + 1):
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
    # drop depth 0 on return (indicating presence of empty taxonomies)
    # these unclassifications are counted, but we don't want to plot at
    # level 0 since it is confusing — level 1 is the real root.
    return depths.drop(0)


def _taxonomic_entropy(taxonomy, rank_handle_regex, max_depth):
    # convert each taxonomy to a list
    taxa_lists = [[t for t in _taxon_to_list(v, rank_handle=rank_handle_regex)
                   if t not in [None, '']]
                  for v in taxonomy.values]
    # taxonomic labels should accumulate at each rank (e.g., to avoid
    # counting identical species names as duplicates when genus is unique)
    taxa_lists = [[';'.join(t[:i]) for i in range(1, max_depth + 1)]
                  for t in taxa_lists]
    # coalesce taxonomies at each rank
    ranks = [i for i in zip_longest(*taxa_lists)]
    # Count unique values at each rank, excluding 'None'
    # use Counter instead of set to pass counts to entropy
    unique_counts = [[v for k, v in Counter(r).items()]
                     for r in ranks]
    # measure entropy at each rank, excluding 'None'
    entropy = {n: [len(r), scipy.stats.entropy(r)]
               for n, r in enumerate(unique_counts, 1)}
    return pd.DataFrame(entropy,
                        index=['Unique Labels', 'Taxonomic Entropy']).T


def _taxonomic_depth(taxonomy, rank_handle_regex):
    return taxonomy.apply(lambda x: len([i for i in _taxon_to_list(
        x, rank_handle=rank_handle_regex) if i not in [None, '']]))
