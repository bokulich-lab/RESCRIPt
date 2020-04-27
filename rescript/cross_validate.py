#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2019--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2 as q2
from sklearn.model_selection import StratifiedKFold
from collections import defaultdict

from .evaluate import _taxonomic_depth


def cross_validate(ctx,
                   sequences,
                   taxonomy,
                   k=3,
                   random_state=0,
                   reads_per_batch=0,
                   n_jobs=1,
                   confidence=0.7):
    '''
    taxonomy: pd.Series of taxonomy labels
    sequences: pd.Series of sequences
    k: number of kfold cv splits to perform.
    random_state: random state for cv.
    '''
    # Validate inputs
    taxa = taxonomy.view(pd.Series)
    _validate_even_rank_taxonomy(taxa)
    seqs = sequences.view(pd.Series)
    _validate_indices_match(taxa, seqs)

    fit = ctx.get_action('feature_classifier', 'fit_classifier_naive_bayes')
    classify = ctx.get_action('feature_classifier', 'classify_sklearn')

    # Deploy perfect classifier! (no CV, lots of data leakage)
    if k == 'disable':
        classifier, = fit(reference_reads=sequences,
                          reference_taxonomy=taxonomy)
        observed_taxonomy, = classify(reads=sequences,
                                      classifier=classifier,
                                      reads_per_batch=reads_per_batch,
                                      n_jobs=n_jobs,
                                      confidence=confidence,
                                      read_orientation='same')
        return taxonomy, observed_taxonomy

    # Otherwise carry on with some form of CV
    taxonomy = taxa
    sequences = seqs

    # TODO: add option for unstratified (CV without fancy stratification)

    # build training/test strata from taxonomic tree to stratify orphan taxa
    strata = _stratify_taxa(taxonomy, sequences, k)
    strata = zip(*[(s, t) for t, ss in strata for s in ss])
    # split taxonomy into training and test sets
    train_test_data = _generate_train_test_data(
        *strata, taxonomy, sequences, k, random_state)
    # now we perform CV classification
    expected_taxonomies = []
    observed_taxonomies = []
    for train_seqs, test_seqs, train_taxa, test_taxa in train_test_data:
        ref_seqs = q2.Artifact.import_data('FeatureData[Sequence]', train_seqs)
        ref_taxa = q2.Artifact.import_data('FeatureData[Taxonomy]', train_taxa)
        reads = q2.Artifact.import_data('FeatureData[Sequence]', test_seqs)
        # TODO: incorporate different methods? taxonomic weights? params?
        classifier, = fit(reference_reads=ref_seqs,
                          reference_taxonomy=ref_taxa)
        observed_taxonomy, = classify(reads=reads,
                                      classifier=classifier,
                                      reads_per_batch=reads_per_batch,
                                      n_jobs=n_jobs,
                                      confidence=confidence,
                                      read_orientation='same')
        # compile observed + expected taxonomies for evaluation outside of loop
        expected_taxonomies.append(test_taxa)
        observed_taxonomies.append(observed_taxonomy.view(pd.Series))

    # Merge expected/observed taxonomies
    expected_taxonomies = pd.concat(expected_taxonomies)
    expected_taxonomies = q2.Artifact.import_data(
        'FeatureData[Taxonomy]', expected_taxonomies)
    observed_taxonomies = pd.concat(observed_taxonomies)
    observed_taxonomies = q2.Artifact.import_data(
        'FeatureData[Taxonomy]', observed_taxonomies)
    return expected_taxonomies, observed_taxonomies


def evaluate_classifications(ctx, expected_taxonomies, observed_taxonomies):
    volatility = ctx.get_action('longitudinal', 'volatility')
    # Validate inputs.
    if len(expected_taxonomies) != len(observed_taxonomies):
        raise ValueError('Expected and Observed Taxonomies do not match. '
                         'Input must contain an equal number of expected and '
                         'observed taxonomies.')
    # Do a quick iteration over input pairs to validate indices match.
    # Why loop over twice? So this does not fail midway through computing
    # results on a big batch of inputs if a user makes a mistake.
    expected_taxonomies = [t.view(pd.Series) for t in expected_taxonomies]
    observed_taxonomies = [t.view(pd.Series) for t in observed_taxonomies]
    for n, (t1, t2) in enumerate(zip(
            expected_taxonomies, observed_taxonomies), 1):
        # if set(t1.index) != set(t2.index):
        try:
            _validate_indices_match(t1, t2)
        except ValueError:
            raise ValueError(
                'Expected and Observed Taxonomies do not match. Taxonomy row '
                'indices must match in each pair of taxonomies. Indices of '
                'pair {0} do not match.'.format(n))

    results = []
    for n, (t1, t2) in enumerate(zip(
            expected_taxonomies, observed_taxonomies), 1):
        # Align Indices
        expected_taxonomy, observed_taxonomy = t1.align(t2)

        # Evaluate classification accuracy
        precision_recall = _calculate_per_rank_precision_recall(
            expected_taxonomy, observed_taxonomy)
        precision_recall['Dataset'] = str(n)
        results.append(precision_recall)
    precision_recall = pd.concat(results)
    # convert index to strings
    precision_recall.index = pd.Index(
        [str(i) for i in range(1, len(precision_recall.index) + 1)], name='id')
    plots, = volatility(metadata=q2.Metadata(precision_recall),
                        state_column='Level',
                        default_group_column='Dataset',
                        default_metric='F-Measure')
    return plots


def _generate_train_test_data(X, y, taxonomy, sequences, k, random_state):
    '''
    X: tuple of feature IDs
    y: tuple of taxonomy labels
    taxonomy: pd.Series of taxonomy labels
    sequences: pd.Series of sequences
    k: number of kfold cv splits to perform.
    random_state: random state for cv.
    '''
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=random_state)
    for train, test in skf.split(X, y):
        train = {X[i] for i in train}
        test = {X[i] for i in test}
        # subset sequences and taxonomies into training/test sets based on ids
        train_seqs = sequences[train]
        test_seqs = sequences[test]
        train_taxa = taxonomy[train]
        test_taxa = taxonomy[test]
        # Compile set of valid stratified taxonomy labels:
        train_taxonomies = set()
        for t in train_taxa.values:
            t = t.split(';')
            for level in range(1, len(t)+1):
                train_taxonomies.add(';'.join(t[:level]))
        # relabel test taxonomy expected labels using stratified set
        # If a taxonomy in the test set doesn't exist in the training set, trim
        # it until it does
        test_taxa = test_taxa.apply(
            lambda x: _relabel_stratified_taxonomy(x, train_taxonomies))
        yield train_seqs, test_seqs, train_taxa, test_taxa


def _relabel_stratified_taxonomy(taxonomy, valid_taxonomies):
    '''
    Relabel taxonomy label to match nearest valid taxonomic label.
    taxonomy: str
    valid_taxonomies: set
    '''
    t = taxonomy.split(';')
    for level in range(len(t), 0, -1):
        if ';'.join(t[:level]) in valid_taxonomies:
            return ';'.join(t[:level]).strip()
    else:
        raise RuntimeError('unknown kingdom in query set')


def _stratify_taxa(taxonomy, sequences, k):
    tree = _build_tree(taxonomy, sequences)
    strata = _get_strata(tree, k)
    return strata


class Node(object):
    def __init__(self):
        self.tips = set()
        self.children = defaultdict(Node)

    def add(self, sid, taxon):
        self.tips.add(sid)
        if taxon:
            self.children[taxon[0]].add(sid, taxon[1:])

    def __lt__(self, other):
        return len(self.tips) < len(other.tips)


def _build_tree(taxonomy, sequences):
    '''
    taxonomy: pd.Series of taxonomic labels.
    '''
    tree = Node()
    filtered_ids = sequences.index
    for sid, taxon in taxonomy.items():
        if sid not in filtered_ids:
            continue
        tree.add(sid, taxon.split(';'))
    return tree


def _get_strata(tree, k, taxon=[]):
    if not tree.children:
        return [(';'.join(taxon), tree.tips)]
    sorted_children = map(list, map(reversed, tree.children.items()))
    sorted_children = iter(sorted(sorted_children))
    misc = set()
    # aggregate children with small tip sets
    for child, level in sorted_children:
        if len(child.tips) >= k:
            break
        misc.update(child.tips)
    else:  # all the tips are in misc
        return [(';'.join(taxon), misc)]
    # grab the tips from the child that's not in misc
    strata = _get_strata(child, k, taxon+[level])
    # get the rest
    for child, level in sorted_children:
        strata.extend(_get_strata(child, k, taxon+[level]))
    # if there were more than k misc, make a stratum for them
    if len(misc) > k:
        strata = [(';'.join(taxon), misc)] + strata
    else:  # randomly add them to the other strata
        for i, tip in enumerate(misc):
            strata[i % len(strata)][1].add(tip)
    return strata


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


def _validate_even_rank_taxonomy(taxa):
    '''
    Check + raise error if taxonomy does not have 100% even taxonomic levels.
    taxa: pd.Series of taxonomic labels.
    '''
    depths = _taxonomic_depth(taxa, rank_handle="")
    uniq_values = depths.unique()
    if len(uniq_values) > 1:
        max_value = max(uniq_values)
        raise ValueError('Taxonomic label depth is uneven. All taxonomies '
                         'must have the same number of semicolon-delimited '
                         'ranks. The following features are too short: ' +
                         ', '.join(depths[depths < max_value].index))


def _validate_indices_match(series1, series2):
    '''
    match indices of two pd.Series objects.
    '''
    diff = series1.index.symmetric_difference(series2.index)
    if len(diff) > 0:
        raise ValueError('Input features must match. The following features '
                         'are missing from one input: ' + ', '.join(diff))
