# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2 as q2

from .evaluate import _taxonomic_depth, _process_labels


def _evaluate_classifications_stats(expected_taxonomies,
                                    observed_taxonomies,
                                    labels=None):
    # Validate inputs.
    if len(expected_taxonomies) != len(observed_taxonomies):
        raise ValueError('Expected and Observed Taxonomies do not match. '
                         'Input must contain an equal number of expected and '
                         'observed taxonomies.')
    labels = _process_labels(labels, expected_taxonomies)
    # Do a quick iteration over input pairs to validate indices match.
    # Why loop over twice? So this does not fail midway through computing
    # results on a big batch of inputs if a user makes a mistake.
    expected_taxonomies = [t.view(pd.Series) for t in expected_taxonomies]
    observed_taxonomies = [t.view(pd.Series) for t in observed_taxonomies]
    for n, (t1, t2) in enumerate(zip(
            expected_taxonomies, observed_taxonomies), 1):
        try:
            _validate_index_is_superset(t1.index, t2.index)
        except ValueError:
            raise ValueError(
                'Expected and Observed Taxonomies do not match. Expected '
                'taxonomy must be a superset of observed taxonomies. Indices '
                'of pair {0} do not match.'.format(n))

    results = []
    for t1, t2, name in zip(expected_taxonomies, observed_taxonomies, labels):
        # Align Indices
        expected_taxonomy, observed_taxonomy = t1.align(t2, join='inner')
        # Evaluate classification accuracy
        precision_recall = _calculate_per_rank_precision_recall(
            expected_taxonomy, observed_taxonomy)
        precision_recall['Dataset'] = str(name)
        results.append(precision_recall)
    precision_recall = pd.concat(results)
    # convert index to strings
    precision_recall.index = pd.Index(
        [str(i) for i in range(1, len(precision_recall.index) + 1)], name='id')

    return q2.Metadata(precision_recall)


def evaluate_classifications(ctx,
                             expected_taxonomies,
                             observed_taxonomies,
                             labels=None):
    lineplot = ctx.get_action('vizard', 'lineplot')

    md = _evaluate_classifications_stats(
        expected_taxonomies=expected_taxonomies,
        observed_taxonomies=observed_taxonomies,
        labels=labels)

    plots, = lineplot(metadata=md,
                      x_measure='Level',
                      y_measure='F-Measure',
                      group_by='Dataset',
                      title='RESCRIPt Evaluate Classifications')
    return plots


def _calculate_per_rank_precision_recall(expected_taxonomies,
                                         observed_taxonomies):
    max_depth = max(_taxonomic_depth(expected_taxonomies, "").max(),
                    _taxonomic_depth(observed_taxonomies, "").max())
    precision_recall = []
    for level in range(1, max_depth + 1):
        exp = expected_taxonomies.apply(
            lambda x: ';'.join(x.split(';')[:level]))
        obs = observed_taxonomies.apply(
            lambda x: ';'.join(x.split(';')[:level]))
        p, r, f = _precision_recall_fscore(exp, obs)
        precision_recall.append((level, p, r, f))
    precision_recall = pd.DataFrame(
        precision_recall,
        columns=['Level', 'Precision', 'Recall', 'F-Measure'])
    return precision_recall


# ported from q2_quality_control with permission of nbokulich
# this computes modified precision calculation: underclassifications count as
# false negatives at level L, but not as false positives.
def _precision_recall_fscore(exp, obs, sample_weight=None):
    # precision, recall, fscore, calculated using microaveraging
    if sample_weight is None:
        sample_weight = [1]*len(exp)
    tp, fp, fn = 0, 0, 0
    for e, o, w in zip(exp, obs, sample_weight):
        if o == e:
            # tp for the true class, the rest are tn
            tp += w
        elif e.startswith(o) or o in (
                'Unclassified', 'Unassigned', 'No blast hit', 'other'):
            # fn for the true class
            fn += w
            # no fp for the predicted class, because it was right to some level
            # the rest are tn
        else:
            # fp the the predicted class
            fp += w
            # fn for the true class, the rest are tn
            fn += w

    # avoid divide by zero error. If no true positives, all scores = 0
    if tp == 0:
        return 0, 0, 0
    else:
        p = tp / (tp + fp)
        r = tp / (tp + fn)
        f = 2.*p*r / (p + r)

    return p, r, f


def _validate_index_is_superset(idx1, idx2):
    '''
    match indices of two pd.Index objects.
    idx1: pd.Index or array-like
    idx2: pd.Index or array-like
    '''
    diff = idx2.difference(idx1)
    if len(diff) > 0:
        raise ValueError('The taxonomy IDs must be a superset of the sequence '
                         'IDs. The following feature IDs are missing from the '
                         'sequences: ' + ', '.join(diff))
