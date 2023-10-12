# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2 as q2
import timeit
from warnings import filterwarnings
from sklearn.model_selection import StratifiedKFold

from q2_types.feature_data import DNAFASTAFormat, DNAIterator

from .evaluate import _taxonomic_depth, _process_labels


def evaluate_fit_classifier(ctx,
                            sequences,
                            taxonomy,
                            reads_per_batch='auto',
                            n_jobs=1,
                            confidence=0.7):
    '''
    taxonomy: FeatureData[Taxonomy] artifact of taxonomy labels
    sequences: FeatureData[Sequence] artifact of sequences
    k: number of kfold cv splits to perform.
    '''
    # Validate inputs
    start = timeit.default_timer()
    taxa, seq_ids = _validate_cross_validate_inputs(taxonomy, sequences)
    taxa = taxa.loc[seq_ids]
    taxonomy = q2.Artifact.import_data('FeatureData[Taxonomy]', taxa)
    new_time = _check_time(start, 'Validation')

    fit = ctx.get_action('feature_classifier', 'fit_classifier_naive_bayes')
    classify = ctx.get_action('feature_classifier', 'classify_sklearn')
    _eval = ctx.get_action('rescript', 'evaluate_classifications')

    # Deploy perfect classifier! (no CV, lots of data leakage)
    classifier, = fit(reference_reads=sequences,
                      reference_taxonomy=taxonomy)
    new_time = _check_time(new_time, 'Training')
    observed_taxonomy, = classify(reads=sequences,
                                  classifier=classifier,
                                  reads_per_batch=reads_per_batch,
                                  n_jobs=n_jobs,
                                  confidence=confidence,
                                  read_orientation='same')
    new_time = _check_time(new_time, 'Classification')
    evaluation, = _eval([taxonomy], [observed_taxonomy])
    _check_time(new_time, 'Evaluation')
    _check_time(start, 'Total Runtime')
    return classifier, evaluation, observed_taxonomy


def evaluate_cross_validate(ctx,
                            sequences,
                            taxonomy,
                            k=3,
                            random_state=0,
                            reads_per_batch='auto',
                            n_jobs=1,
                            confidence=0.7):
    '''
    taxonomy: FeatureData[Taxonomy] artifact of taxonomy labels
    sequences: FeatureData[Sequence] artifact of sequences
    k: number of kfold cv splits to perform.
    random_state: random state for cv.
    '''
    # silence impertinent sklearn warnings:
    # 1. classifier version (the classifier is not saved or reused)
    # 2. class sizes (this is handled by taxonomic stratification/relabeling)
    msg = 'The TaxonomicClassifier.*cannot be used with other versions'
    filterwarnings("ignore", message=msg, category=UserWarning)
    filterwarnings(
        "ignore", message='The least populated class', category=UserWarning)
    # Validate inputs
    start = timeit.default_timer()
    taxa, seq_ids = _validate_cross_validate_inputs(taxonomy, sequences)
    new_time = _check_time(start, 'Validation')

    fit = ctx.get_action('feature_classifier', 'fit_classifier_naive_bayes')
    classify = ctx.get_action('feature_classifier', 'classify_sklearn')
    _eval = ctx.get_action('rescript', 'evaluate_classifications')

    # split taxonomy into training and test sets
    train_test_data = _generate_train_test_data(taxa, k, random_state)
    # now we perform CV classification
    expected_taxonomies = []
    observed_taxonomies = []
    for n, (train_taxa, test_taxa) in enumerate(train_test_data):
        train_ids = train_taxa.index
        test_ids = test_taxa.index
        train_seqs, test_seqs = _split_fasta(sequences, train_ids, test_ids)
        ref_taxa = q2.Artifact.import_data('FeatureData[Taxonomy]', train_taxa)
        new_time = _check_time(new_time, 'Fold {0} split'.format(n))
        # TODO: incorporate different methods? taxonomic weights? params?
        classifier, = fit(reference_reads=train_seqs,
                          reference_taxonomy=ref_taxa)
        new_time = _check_time(new_time, 'Fold {0} fit'.format(n))
        observed_taxonomy, = classify(reads=test_seqs,
                                      classifier=classifier,
                                      reads_per_batch=reads_per_batch,
                                      n_jobs=n_jobs,
                                      confidence=confidence,
                                      read_orientation='same')
        # compile observed + expected taxonomies for evaluation outside of loop
        expected_taxonomies.append(test_taxa)
        observed_taxonomies.append(observed_taxonomy.view(pd.Series))
        new_time = _check_time(new_time, 'Fold {0} classify'.format(n))

    # Merge expected/observed taxonomies
    expected_taxonomies = q2.Artifact.import_data(
        'FeatureData[Taxonomy]', pd.concat(expected_taxonomies))
    observed_taxonomies = q2.Artifact.import_data(
        'FeatureData[Taxonomy]', pd.concat(observed_taxonomies))
    evaluation, = _eval([expected_taxonomies], [observed_taxonomies])
    _check_time(new_time, 'Evaluation')
    _check_time(start, 'Total Runtime')
    return expected_taxonomies, observed_taxonomies, evaluation


def _check_time(old_time, name='Time'):
    new_time = timeit.default_timer()
    print('{0}: {1:.2f}s'.format(name, new_time - old_time))
    return new_time


# input validation for cross-validation functions
def _validate_cross_validate_inputs(taxonomy, sequences):
    taxa = taxonomy.view(pd.Series)
    # taxonomies must have even ranks (this is used for confidence estimation
    # with the current NB classifier; we could relax this if we implement other
    # methods later on for CV classification).
    _validate_even_rank_taxonomy(taxa)
    seq_ids = {i.metadata['id'] for i in sequences.view(DNAIterator)}
    _validate_index_is_superset(set(taxa.index), seq_ids)
    return taxa, seq_ids


def _split_fasta(sequences, train_ids, test_ids):
    '''
    Split FeatureData[Sequence] artifact into two, based on two sets of IDs.
    sequences: FeatureData[Sequence] Artifact
    train_ids: set
    test_ids: set
    '''
    train_seqs = DNAFASTAFormat()
    test_seqs = DNAFASTAFormat()
    with train_seqs.open() as _train, test_seqs.open() as _test:
        for s in sequences.view(DNAIterator):
            _id = s.metadata['id']
            if s.metadata['id'] in train_ids:
                _train.write('>%s\n%s\n' % (_id, str(s)))
            elif s.metadata['id'] in test_ids:
                _test.write('>%s\n%s\n' % (_id, str(s)))
    train_seqs = q2.Artifact.import_data('FeatureData[Sequence]', train_seqs)
    test_seqs = q2.Artifact.import_data('FeatureData[Sequence]', test_seqs)
    return train_seqs, test_seqs


def evaluate_classifications(ctx,
                             expected_taxonomies,
                             observed_taxonomies,
                             labels=None):
    volatility = ctx.get_action('longitudinal', 'volatility')
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
    plots, = volatility(metadata=q2.Metadata(precision_recall),
                        state_column='Level',
                        default_group_column='Dataset',
                        default_metric='F-Measure')
    return plots


def _generate_train_test_data(taxonomy, k, random_state):
    '''
    taxonomy: pd.Series of taxonomy labels
    k: number of kfold cv splits to perform.
    random_state: random state for cv.
    '''
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=random_state)
    for train, test in skf.split(taxonomy.index, taxonomy.values):
        # subset sequences and taxonomies into training/test sets based on ids
        train_taxa = taxonomy.iloc[train]
        test_taxa = taxonomy.iloc[test]
        # Compile set of valid stratified taxonomy labels:
        train_taxonomies = _get_valid_taxonomic_labels(train_taxa)
        # relabel test taxonomy expected labels using stratified set
        # If a taxonomy in the test set doesn't exist in the training set, trim
        # it until it does
        test_taxa = test_taxa.apply(
            lambda x: _relabel_stratified_taxonomy(x, train_taxonomies))
        yield train_taxa, test_taxa


def _get_valid_taxonomic_labels(taxonomy):
    valid_labels = {
        ';'.join(t.split(';')[:level]) for t in taxonomy.unique()
        for level in range(1, len(t.split(';'))+1)}
    return valid_labels


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
        raise RuntimeError('unknown kingdom in query set: ' + taxonomy)


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
    depths = taxa.str.count(';')
    uniq_values = depths.unique()
    if len(uniq_values) > 1:
        max_value = max(uniq_values)
        raise ValueError('Taxonomic label depth is uneven. All taxonomies '
                         'must have the same number of semicolon-delimited '
                         'ranks. The following features are too short: ' +
                         ', '.join(depths[depths < max_value].index))


def _validate_indices_match(idx1, idx2):
    '''
    match indices of two pd.Index objects.
    idx1: pd.Index
    idx2: pd.Index or array-like
    '''
    diff = idx1.symmetric_difference(idx2)
    if len(diff) > 0:
        raise ValueError('Input features must match. The following features '
                         'are missing from one input: ' + ', '.join(diff))


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
