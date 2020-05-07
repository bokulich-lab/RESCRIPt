#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Str, Plugin, Choices, List, Citations, Range, Int,
                           Float, Visualization)
from .merge import merge_taxa
from .dereplicate import dereplicate
from .evaluate import evaluate_taxonomy
from .cleanseq import clean_sequences
from .cross_validate import cross_validate, evaluate_classifications
from q2_types.feature_data import FeatureData, Taxonomy, Sequence
from q2_feature_classifier.classifier import (_parameter_descriptions,
                                              _classify_parameters)

import rescript


citations = Citations.load('citations.bib', package='rescript')

plugin = Plugin(
    name='rescript',
    version=rescript.__version__,
    website="https://github.com/nbokulich/RESCRIPt",
    package='rescript',
    description=('Reference sequence annotation and curation pipeline.'),
    short_description=(
        'Pipeline for reference sequence annotation and curation.'),
)


VOLATILITY_PLOT_XAXIS_INTERPRETATION = (
    'The x-axis in these plots represents the taxonomic '
    'levels present in the input taxonomies so are labeled numerically '
    'instead of by rank, but typically for 7-level taxonomies these will '
    'represent: 1 = domain/kingdom, 2 = phylum, 3 = class, 4 = order, '
    '5 = family, 6 = genus, 7 = species.')

rank_handle_description = (
    'Regular expression indicating which taxonomic rank a label '
    'belongs to; this handle is stripped from the label '
    'prior to operating on the taxonomy. The net '
    'effect is that ambiguous or empty levels can be '
    'removed prior to comparison, enabling selection of '
    'taxonomies with more complete taxonomic information. '
    'For example, "^[kpcofgs]__" will recognize greengenes rank '
    'handles. ')

rank_handle_extra_note = (
    'Note that rank_handles are removed but not replaced; use the '
    'new_rank_handle parameter to replace the rank handles.')


plugin.pipelines.register_function(
    function=cross_validate,
    inputs={'sequences': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy]},
    parameters={
        'k': Int % Range(2, None) | Str % Choices(['disable']),
        'random_state': Int % Range(0, None),
        'reads_per_batch': _classify_parameters['reads_per_batch'],
        'n_jobs': _classify_parameters['n_jobs'],
        'confidence': _classify_parameters['confidence']},
    outputs=[('expected_taxonomy', FeatureData[Taxonomy]),
             ('observed_taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'sequences': 'Reference sequences to use for classifier '
                     'training/testing.',
        'taxonomy': 'Reference taxonomy to use for classifier '
                    'training/testing.'},
    parameter_descriptions={
        'k': 'Number of stratified folds. Set to "disable" to disable k-fold '
             'cross-validation. This results in a "perfect" classifier that '
             'knows the correct identity of each input sequence. Such a leaky '
             'classifier indicates the upper limit of classification accuracy '
             'based on sequence information alone, as misclassifications are '
             'an indication of unresolvable kmer profiles.',
        'random_state': 'Seed used by the random number generator.',
        'reads_per_batch': _parameter_descriptions['reads_per_batch'],
        'n_jobs': _parameter_descriptions['n_jobs'],
        'confidence': _parameter_descriptions['confidence']},
    output_descriptions={
        'expected_taxonomy': 'Expected taxonomic label for each input '
                             'sequence. If k-fold CV is enabled, taxonomic '
                             'labels may be truncated due to stratification.',
        'observed_taxonomy': 'Observed taxonomic label for each input '
                             'sequence, predicted by cross-validation.'},
    name=('Evaluate DNA sequence reference database via cross-validated '
          'taxonomic classification.'),
    description=(
        'Evaluate DNA sequence reference database via cross-validated '
        'taxonomic classification. Unique taxonomic labels are truncated to '
        'enable appropriate label stratification. See the cited reference '
        '(Bokulich et al. 2018) for more details.'),
    citations=[citations['bokulich2018optimizing']]
)


plugin.pipelines.register_function(
    function=evaluate_classifications,
    inputs={'expected_taxonomies': List[FeatureData[Taxonomy]],
            'observed_taxonomies': List[FeatureData[Taxonomy]]},
    parameters={},
    outputs=[('evaluation', Visualization)],
    input_descriptions={
        'expected_taxonomies': 'True taxonomic labels for one more more sets '
                               'of features.',
        'observed_taxonomies': 'Predicted classifications of same sets of '
                               'features, input in same order as '
                               'expected_taxonomies.'},
    parameter_descriptions={},
    output_descriptions={
        'evaluation': 'Visualization of classification accuracy results.'},
    name=('Interactively evaluate taxonomic classification accuracy.'),
    description=(
        'Evaluate taxonomic classification accuracy by comparing one or more '
        'sets of true taxonomic labels to the predicted taxonomies for the '
        'same set(s) of features. Output an interactive line plot of '
        'classification accuracy for each pair of expected/observed '
        'taxonomies. ' + VOLATILITY_PLOT_XAXIS_INTERPRETATION),
    citations=[citations['bokulich2018optimizing'],
               citations['bokulich2017q2']]
)


plugin.methods.register_function(
    function=merge_taxa,
    inputs={'data': List[FeatureData[Taxonomy]]},
    parameters={
        'mode': Str % Choices(['len', 'lca', 'score']),
        'rank_handle': Str,
        'new_rank_handle': Str},
    outputs=[('merged_data', FeatureData[Taxonomy])],
    input_descriptions={
        'data': 'Two or more feature taxonomies to be merged.'},
    parameter_descriptions={
        'mode': 'How to merge feature taxonomies: "len" will select the '
                'taxonomy with the most elements (e.g., species level will '
                'beat genus level); "lca" will find the least common ancestor '
                'and report this consensus taxonomy; "score" will select the '
                'taxonomy with the highest score (e.g., confidence or '
                'consensus score). Note that "score" assumes that this score '
                'is always contained as the second column in a feature '
                'taxonomy dataframe.',
        'rank_handle': rank_handle_description + rank_handle_extra_note,
        'new_rank_handle': 'A semicolon-delimited string of rank handles to '
                           'prepend to taxonomic labels at each rank. For '
                           'example, "k__;p__;c__;o__;f__;g__;s__" will '
                           'prepend greengenes-style rank handles. Note that '
                           'merged taxonomies will only contain as many '
                           'levels as there are handles if this parameter is '
                           'used. So "k__;p__" will trim all taxonomies to '
                           'phylum level, even if longer annotations exist.'},
    name='Compare taxonomies and select the longest, highest scoring, or find '
         'the least common ancestor.',
    description='Compare taxonomy annotations and choose the best one. Can '
                'select the longest taxonomy annotation, the highest scoring, '
                'or the least common ancestor. Note: when a tie occurs, the '
                'last taxonomy added takes precedent.',
)


plugin.methods.register_function(
    function=dereplicate,
    inputs={'sequences': FeatureData[Sequence],
            'taxa': FeatureData[Taxonomy]},
    parameters={
        'mode': Str % Choices(['uniq', 'lca', 'majority']),
        'threads': Int % Range(1, 256),
        'perc_identity': Float % Range(0, 1, inclusive_start=False,
                                       inclusive_end=True)},
    outputs=[('dereplicated_sequences', FeatureData[Sequence]),
             ('dereplicated_taxa', FeatureData[Taxonomy])],
    input_descriptions={
        'sequences': 'Sequences to be dereplicated',
        'taxa': 'Taxonomic classifications of sequences to be dereplicated'},
    parameter_descriptions={
        'mode': 'How to handle dereplication when sequences map to distinct '
                'taxonomies. "uniq" will retain all sequences with unique '
                'taxonomic affiliations. "lca" will find the least common '
                'ancestor among all taxa sharing a sequence. "majority" will '
                'find the most common taxonomic label associated with that '
                'sequence; note that in the event of a tie, "majority" will '
                'pick the winner arbitrarily.',
        'threads': 'Number of computation threads to use (1 to 256). The '
                   'number of threads should be lesser or equal to the number '
                   'of available CPU cores.',
        'perc_identity': 'The percent identity at which clustering should be '
                         'performed. This parameter maps to vsearch\'s --id '
                         'parameter.'
    },
    name='Dereplicate features with matching sequences and taxonomies.',
    description=(
        'Dereplicate FASTA format sequences and taxonomies wherever '
        'sequences and taxonomies match; duplicated sequences and taxonomies '
        'are dereplicated using the "mode" parameter to either: retain all '
        'sequences that have unique taxonomic annotations even if the '
        'sequences are duplicates (uniq); or return only dereplicated '
        'sequences labeled by either the least common ancestor (lca) or the '
        'most common taxonomic label associated with sequences in that '
        'cluster (majority).'),
    citations=[citations['rognes2016vsearch']]
)


plugin.pipelines.register_function(
    function=evaluate_taxonomy,
    inputs={'taxonomies': List[FeatureData[Taxonomy]]},
    parameters={'rank_handle': Str},
    outputs=[('taxonomy_stats', Visualization)],
    input_descriptions={
        'taxonomies': 'One or more taxonomies to evaluate.'},
    parameter_descriptions={
        'rank_handle': rank_handle_description,
    },
    name='Compute summary statistics on taxonomy artifact(s).',
    description=(
        'Compute summary statistics on taxonomy artifact(s) and visualize as '
        'interactive lineplots. Summary statistics include the number of '
        'unique labels, taxonomic entropy, and the number of features that '
        'are (un)classified at each taxonomic level. This action is useful '
        'for both reference taxonomies and classification results. ' +
        VOLATILITY_PLOT_XAXIS_INTERPRETATION),
    citations=[citations['bokulich2017q2']]
)

plugin.methods.register_function(
    function=clean_sequences,
    inputs={
        'sequences': FeatureData[Sequence]
        },
    parameters={
        'max_degenerates': Int % Range(0, None),
        'homopolymer_length': Int % Range(2, None)
        },
    outputs=[('cleaned_sequences', FeatureData[Sequence])],
    input_descriptions={
        'sequences': 'Sequences to be removed that contain more than the '
                     'allowed number of degenerate bases and/or long '
                     'homopolymers.'
        },
    parameter_descriptions={
        'max_degenerates': 'Maximum number of degenerate (IUPAC compliant) '
                           'DNA bases allowed in a sequence.',
        'homopolymer_length': 'Sequences with a homopolymer of this length '
                              'or greater will be removed',
    },
    output_descriptions={
        'cleaned_sequences': 'The resulting cleaned sequences.'
        },
    name='Remove sequences with excessive ambiguous base calls and/or '
         'homopolymers.',
    description=(
        'Removes DNA sequences that have the specified number, or more, '
        'of IUPAC compliant ambiguous bases. Remaining sequences are removed '
        'if they contain homopolymers equal to or longer than the specified '
        'length.'
        )
)
