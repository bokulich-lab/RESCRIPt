# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.core.type import TypeMatch
from qiime2.plugin import (Str, Plugin, Choices, List, Citations, Range, Int,
                           Float, Visualization, Bool, TypeMap, Metadata)

from .subsample import subsample_fasta
from .trim_alignment import trim_alignment
from .merge import merge_taxa
from .dereplicate import dereplicate
from .evaluate import evaluate_taxonomy, evaluate_seqs
from .screenseq import cull_seqs
from .degap import degap_seqs
from .parse_silva_taxonomy import (parse_silva_taxonomy, ALLOWED_RANKS,
                                   DEFAULT_RANKS)
from .get_data import get_silva_data
from .cross_validate import (evaluate_cross_validate,
                             evaluate_classifications,
                             evaluate_fit_classifier)
from .filter_length import (filter_seqs_length_by_taxon, filter_seqs_length,
                            filter_taxa)
from .orient import orient_seqs
from q2_types.feature_data import (FeatureData, Taxonomy, Sequence,
                                   AlignedSequence, RNASequence,
                                   AlignedRNASequence, ProteinSequence)
from q2_types.tree import Phylogeny, Rooted
from q2_feature_classifier.classifier import (_parameter_descriptions,
                                              _classify_parameters)
from q2_feature_classifier._taxonomic_classifier import TaxonomicClassifier

import rescript
from rescript._utilities import _rank_handles
from rescript.types._format import (
    SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat, SILVATaxidMapFormat,
    SILVATaxidMapDirectoryFormat)
from rescript.types._type import SILVATaxonomy, SILVATaxidMap
from rescript.types.methods import reverse_transcribe
from rescript.ncbi import (
    get_ncbi_data, _default_ranks, _allowed_ranks, get_ncbi_data_protein)


citations = Citations.load('citations.bib', package='rescript')

plugin = Plugin(
    name='rescript',
    version=rescript.__version__,
    website="https://github.com/nbokulich/RESCRIPt",
    package='rescript',
    description=('Reference sequence annotation and curation pipeline.'),
    short_description=(
        'Pipeline for reference sequence annotation and curation.'),
    citations=[citations['Robeson2020rescript']]
)


SILVA_LICENSE_NOTE = (
    'NOTE: THIS ACTION ACQUIRES DATA FROM THE SILVA DATABASE. SEE '
    'https://www.arb-silva.de/silva-license-information/ FOR MORE INFORMATION '
    'and be aware that earlier versions may be released under a different '
    'license.')

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
    'For example, "^[dkpcofgs]__" will recognize greengenes or silva rank '
    'handles. ')

rank_handle_extra_note = (
    'Note that rank_handles are removed but not replaced; use the '
    'new_rank_handle parameter to replace the rank handles.')

labels_description = (
    'List of labels to use for labeling evaluation results in the resulting '
    'visualization. Inputs are labeled with labels in the order that each '
    'is input. If there are fewer labels than inputs (or no labels), '
    'unnamed inputs are labeled numerically in sequential order. Extra '
    'labels are ignored.')

super_lca_desc = (
    '"super" finds the LCA consensus while giving preference to '
    'majority labels and collapsing substrings into superstrings. '
    'For example, when a more specific taxonomy does not '
    'contradict a less specific taxonomy, the more specific is '
    'chosen. That is, "g__Faecalibacterium; s__prausnitzii", '
    'will be preferred over "g__Faecalibacterium; s__"')

plugin.pipelines.register_function(
    function=evaluate_fit_classifier,
    inputs={'sequences': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy]},
    parameters={
        'reads_per_batch': _classify_parameters['reads_per_batch'],
        'n_jobs': _classify_parameters['n_jobs'],
        'confidence': _classify_parameters['confidence']},
    outputs=[('classifier', TaxonomicClassifier),
             ('evaluation', Visualization),
             ('observed_taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'sequences': 'Reference sequences to use for classifier '
                     'training/testing.',
        'taxonomy': 'Reference taxonomy to use for classifier '
                    'training/testing.'},
    parameter_descriptions={
        'reads_per_batch': _parameter_descriptions['reads_per_batch'],
        'n_jobs': _parameter_descriptions['n_jobs'],
        'confidence': _parameter_descriptions['confidence']},
    output_descriptions={
        'classifier': 'Trained naive Bayes taxonomic classifier.',
        'evaluation': 'Visualization of classification accuracy results.',
        'observed_taxonomy': 'Observed taxonomic label for each input '
                             'sequence, predicted by the trained classifier.'},
    name=('Evaluate and train naive Bayes classifier on reference sequences.'),
    description=(
        'Train a naive Bayes classifier on a set of reference sequences, then '
        'test performance accuracy on this same set of sequences. This '
        'results in a "perfect" classifier that "knows" the correct identity '
        'of each input sequence. Such a leaky classifier indicates the upper '
        'limit of classification accuracy based on sequence information '
        'alone, as misclassifications are an indication of unresolvable kmer '
        'profiles. This test simulates the case where all query sequences '
        'are present in a fully comprehensive reference database. To simulate '
        'more realistic conditions, see `evaluate_cross_validate`. THE '
        'CLASSIFIER OUTPUT BY THIS PIPELINE IS PRODUCTION-READY and can be '
        're-used for classification of other sequences (provided the '
        'reference data are viable), hence THIS PIPELINE IS USEFUL FOR '
        'TRAINING FEATURE CLASSIFIERS AND THEN EVALUATING THEM ON-THE-FLY.'),
    citations=[citations['bokulich2018optimizing']]
)


plugin.pipelines.register_function(
    function=evaluate_cross_validate,
    inputs={'sequences': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy]},
    parameters={
        'k': Int % Range(2, None),
        'random_state': Int % Range(0, None),
        'reads_per_batch': _classify_parameters['reads_per_batch'],
        'n_jobs': _classify_parameters['n_jobs'],
        'confidence': _classify_parameters['confidence']},
    outputs=[('expected_taxonomy', FeatureData[Taxonomy]),
             ('observed_taxonomy', FeatureData[Taxonomy]),
             ('evaluation', Visualization)],
    input_descriptions={
        'sequences': 'Reference sequences to use for classifier '
                     'training/testing.',
        'taxonomy': 'Reference taxonomy to use for classifier '
                    'training/testing.'},
    parameter_descriptions={
        'k': 'Number of stratified folds.',
        'random_state': 'Seed used by the random number generator.',
        'reads_per_batch': _parameter_descriptions['reads_per_batch'],
        'n_jobs': _parameter_descriptions['n_jobs'],
        'confidence': _parameter_descriptions['confidence']},
    output_descriptions={
        'expected_taxonomy': 'Expected taxonomic label for each input '
                             'sequence. Taxonomic labels may be truncated due '
                             'to k-fold CV and stratification.',
        'observed_taxonomy': 'Observed taxonomic label for each input '
                             'sequence, predicted by cross-validation.',
        'evaluation': 'Visualization of cross-validated accuracy results.'},
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
    parameters={'labels': List[Str]},
    outputs=[('evaluation', Visualization)],
    input_descriptions={
        'expected_taxonomies': 'True taxonomic labels for one more more sets '
                               'of features.',
        'observed_taxonomies': 'Predicted classifications of same sets of '
                               'features, input in same order as '
                               'expected_taxonomies.'},
    parameter_descriptions={'labels': labels_description},
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
        'mode': Str % Choices(['len', 'lca', 'score', 'super', 'majority']),
        'rank_handle_regex': Str,
        'new_rank_handle': Str % Choices(list(_rank_handles.keys())),
        'unclassified_label': Str
    },
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
                'taxonomy dataframe. "majority" finds the LCA consensus while '
                'giving preference to majority labels. ' + super_lca_desc,
        'rank_handle_regex': rank_handle_description + rank_handle_extra_note,
        'new_rank_handle': (
            'Specifies the set of rank handles to prepend to taxonomic labels '
            'at each rank. For example, "greengenes" will prepend 7-level '
            'greengenes-style rank handles. Note that merged taxonomies will '
            'only contain as many levels as there are handles if this '
            'parameter is used. So "greengenes" will trim all taxonomies to '
            'seven levels, even if longer annotations exist. Note that this '
            'parameter will prepend rank handles whether or not they already '
            'exist in the taxonomy, so should ALWAYS be used in conjunction '
            'with `rank_handle_regex` if rank handles exist in any of the '
            'inputs.'),
        'unclassified_label': 'Specifies what label should be used for '
                              'taxonomies that could not be resolved (when '
                              'LCA modes are used).'
    },
    name='Compare taxonomies and select the longest, highest scoring, or find '
         'the least common ancestor.',
    description='Compare taxonomy annotations and choose the best one. Can '
                'select the longest taxonomy annotation, the highest scoring, '
                'or the least common ancestor. Note: when a tie occurs, the '
                'last taxonomy added takes precedent.',
)


VSEARCH_PARAMS = {
    'threads': Int % Range(1, 256),
    'perc_identity': Float % Range(0, 1, inclusive_start=False,
                                   inclusive_end=True),
    'left_justify': Bool,
    'query_cov': Float % Range(0.0, 1.0, inclusive_end=True),
}

VSEARCH_PARAM_DESCRIPTIONS = {
    'threads': 'Number of computation threads to use (1 to 256). The '
               'number of threads should be lesser or equal to the number '
               'of available CPU cores.',
    'perc_identity': 'The percent identity at which clustering should be '
                     'performed. This parameter maps to vsearch\'s --id '
                     'parameter.',
    'left_justify': 'Reject match if the pairwise alignment begins with gaps',
    'query_cov': 'Reject match if query alignment coverage per high-scoring '
                 'pair is lower.',
}


plugin.methods.register_function(
    function=dereplicate,
    inputs={'sequences': FeatureData[Sequence],
            'taxa': FeatureData[Taxonomy]},
    parameters={
        'mode': Str % Choices(['uniq', 'lca', 'majority', 'super']),
        'threads': VSEARCH_PARAMS['threads'],
        'perc_identity': VSEARCH_PARAMS['perc_identity'],
        'derep_prefix': Bool,
        'rank_handles': Str % Choices(list(_rank_handles.keys()))},
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
                'pick the winner arbitrarily. ' + super_lca_desc,
        'threads': VSEARCH_PARAM_DESCRIPTIONS['threads'],
        'perc_identity': VSEARCH_PARAM_DESCRIPTIONS['perc_identity'],
        'derep_prefix': 'Merge sequences with identical prefixes. If a '
                        'sequence is identical to the prefix of two or more '
                        'longer sequences, it is clustered with the shortest '
                        'of them. If they are equally long, it is clustered '
                        'with the most abundant.',
        'rank_handles': (
            'Specifies the set of rank handles used to backfill '
            'missing ranks in the resulting dereplicated taxonomy. The '
            'default setting will backfill SILVA-style 7-level rank handles. '
            'Set to "none" to disable backfilling.')
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
    parameters={'labels': List[Str],
                'rank_handle_regex': Str},
    outputs=[('taxonomy_stats', Visualization)],
    input_descriptions={
        'taxonomies': 'One or more taxonomies to evaluate.'},
    parameter_descriptions={
        'labels': labels_description,
        'rank_handle_regex': rank_handle_description,
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


palettes = ['Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired',
            'Accent', 'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c',
            'viridis', 'plasma', 'inferno', 'magma', 'cividis', 'terrain',
            'rainbow', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']

plugin.visualizers.register_function(
    function=evaluate_seqs,
    inputs={'sequences': List[FeatureData[Sequence]]},
    parameters={'labels': List[Str],
                'kmer_lengths': List[Int % Range(1, None)],
                'subsample_kmers': Float % Range(0, 1, inclusive_start=False,
                                                 inclusive_end=True),
                'palette': Str % Choices(palettes)},
    input_descriptions={
        'sequences': 'One or more sets of sequences to evaluate.'},
    parameter_descriptions={
        'labels': labels_description,
        'kmer_lengths': 'Sequence kmer lengths to optionally use for entropy '
                        'calculation. Warning: kmer entropy calculations may '
                        'be time-consuming for large sequence sets.',
        'subsample_kmers': 'Optionally subsample sequences prior to kmer '
                           'entropy measurement. A fraction of the input '
                           'sequences will be randomly subsampled at the '
                           'specified value.',
        'palette': 'Color palette to use for plotting evaluation results.'
    },
    name='Compute summary statistics on sequence artifact(s).',
    description=(
        'Compute summary statistics on sequence artifact(s) and visualize. '
        'Summary statistics include the number of unique sequences, sequence '
        'entropy, kmer entropy, and sequence length distributions. This '
        'action is useful for both reference taxonomies and classification '
        'results.')
)


plugin.methods.register_function(
    function=cull_seqs,
    inputs={
        'sequences': FeatureData[Sequence | RNASequence]
    },
    parameters={
        'num_degenerates': Int % Range(1, None),
        'homopolymer_length': Int % Range(2, None),
        'n_jobs': Int % Range(1, None)
    },
    outputs=[('clean_sequences', FeatureData[Sequence])],
    input_descriptions={
        'sequences': 'DNA or RNA Sequences to be screened for removal based '
                     'on degenerate base and homopolymer screening criteria.'
    },
    parameter_descriptions={
        'num_degenerates': 'Sequences with N, or more, degenerate bases will '
                           'be removed.',
        'homopolymer_length': 'Sequences containing a homopolymer sequence of '
                              'length N, or greater, will be removed.',
        'n_jobs': 'Number of concurrent processes to use while processing '
                  'sequences. More is faster but typically should not be '
                  'higher than the number of available CPUs. Output sequence '
                  'order may change when using multiple jobs.'
    },
    output_descriptions={
        'clean_sequences': 'The resulting DNA sequences that pass degenerate '
                           'base and homopolymer screening criteria.'
    },
    name='Removes sequences that contain at least the specified number of '
         'degenerate bases and/or homopolymers of a given length.',
    description=(
        'Filter DNA or RNA sequences that contain ambiguous bases and '
        'homopolymers, and output filtered DNA sequences. Removes DNA '
        'sequences that have the specified number, or more, of IUPAC '
        'compliant degenerate bases. Remaining sequences are removed if they '
        'contain homopolymers equal to or longer than the specified length. '
        'If the input consists of RNA sequences, they are reverse transcribed '
        'to DNA before filtering.')
)


plugin.methods.register_function(
    function=degap_seqs,
    inputs={
        'aligned_sequences': FeatureData[AlignedSequence]
    },
    parameters={
        'min_length': Int % Range(1, None)
    },
    outputs=[('degapped_sequences', FeatureData[Sequence])],
    input_descriptions={
        'aligned_sequences': 'Aligned DNA Sequences to be degapped.'
    },
    parameter_descriptions={
        'min_length': 'Minimum length of sequence to be returned after '
                      'degapping.'},
    output_descriptions={
        'degapped_sequences': 'The resulting unaligned (degapped) DNA '
                              'sequences.'
    },
    name='Remove gaps from DNA sequence alignments.',
    description=('This method converts aligned DNA sequences to unaligned DNA '
                 'sequences by removing gaps ("-") and missing data (".") '
                 'characters from the sequences. Essentially, \'unaligning\' '
                 'the sequences.')
)


FILTER_PARAMS = {
    'global_min': Int % Range(1, None),
    'global_max': Int % Range(1, None)}

FILTER_PARAM_DESCRIPTIONS = {
    'global_min': 'The minimum length threshold for filtering all '
                  'sequences. Any sequence shorter than this length will '
                  'be removed.',
    'global_max': 'The maximum length threshold for filtering all '
                  'sequences. Any sequence longer than this length will '
                  'be removed.'}

FILTER_OUTPUT_DESCRIPTIONS = {
    'filtered_seqs': 'Sequences that pass the filtering thresholds.',
    'discarded_seqs': 'Sequences that fall outside the filtering thresholds.'}


plugin.methods.register_function(
    function=orient_seqs,
    inputs={'sequences': FeatureData[Sequence],
            'reference_sequences': FeatureData[Sequence]},
    parameters=VSEARCH_PARAMS,
    outputs=[('oriented_seqs', FeatureData[Sequence]),
             ('unmatched_seqs', FeatureData[Sequence])],
    input_descriptions={
        'sequences': 'Sequences to be oriented.',
        'reference_sequences': 'Reference sequences to orient against.'},
    parameter_descriptions=VSEARCH_PARAM_DESCRIPTIONS,
    output_descriptions={
        'oriented_seqs': 'Query sequences in same orientation as top matching '
                         'reference sequence.',
        'unmatched_seqs': 'Query sequences that fail to match at least one '
                          'reference sequence in either + or - orientation.'},
    name='Orient input sequences by comparison against reference.',
    description=(
        'Orient input sequences by comparison against a set of reference '
        'sequences using VSEARCH. This action can also be used to quickly '
        'filter out sequences that (do not) match a set of reference '
        'sequences in either orientation.'
    ),
    citations=[citations['rognes2016vsearch']]
)


plugin.methods.register_function(
    function=filter_seqs_length_by_taxon,
    inputs={'sequences': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy]},
    parameters={
        'labels': List[Str],
        'min_lens': List[Int % Range(1, None)],
        'max_lens': List[Int % Range(1, None)],
        **FILTER_PARAMS},
    outputs=[('filtered_seqs', FeatureData[Sequence]),
             ('discarded_seqs', FeatureData[Sequence])],
    input_descriptions={
        'sequences': 'Sequences to be filtered by length.',
        'taxonomy': 'Taxonomic classifications of sequences to be filtered.'},
    parameter_descriptions={
        'labels': 'One or more taxonomic labels to use for conditional '
                  'filtering. For example, use this option to set different '
                  'min/max filter settings for individual phyla. Must input '
                  'the same number of labels as min_lens and/or max_lens. If '
                  'a sequence matches multiple taxonomic labels, this method '
                  'will apply the most stringent threshold(s): the longest '
                  'minimum length and/or the shortest maximum length that is '
                  'associated with the matching labels.',
        'min_lens': 'Minimum length thresholds to use for filtering sequences '
                    'associated with each label. If any min_lens are '
                    'specified, must have the same number of min_lens as '
                    'labels. Sequences that contain this label in their '
                    'taxonomy will be removed if they are less than the '
                    'specified length.',
        'max_lens': 'Maximum length thresholds to use for filtering sequences '
                    'associated with each label. If any max_lens are '
                    'specified, must have the same number of max_lens as '
                    'labels. Sequences that contain this label in their '
                    'taxonomy will be removed if they are more than the '
                    'specified length.',
        **FILTER_PARAM_DESCRIPTIONS},
    output_descriptions=FILTER_OUTPUT_DESCRIPTIONS,
    name='Filter sequences by length and taxonomic group.',
    description=(
        'Filter sequences by length. Can filter both globally by minimum '
        'and/or maximum length, and set individual threshold for individual '
        'taxonomic groups (using the "labels" option). Note that filtering '
        'can be performed for multiple taxonomic groups simultaneously, and '
        'nested taxonomic filters can be applied (e.g., to apply a more '
        'stringent filter for a particular genus, but a less stringent filter '
        'for other members of the kingdom). For global length-based filtering '
        'without conditional taxonomic filtering, see filter_seqs_length.'),
)


plugin.methods.register_function(
    function=filter_seqs_length,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={
        **FILTER_PARAMS,
        'threads': VSEARCH_PARAMS['threads']},
    outputs=[('filtered_seqs', FeatureData[Sequence]),
             ('discarded_seqs', FeatureData[Sequence])],
    input_descriptions={
        'sequences': 'Sequences to be filtered by length.'},
    parameter_descriptions={
        **FILTER_PARAM_DESCRIPTIONS,
        'threads': VSEARCH_PARAM_DESCRIPTIONS['threads']},
    output_descriptions=FILTER_OUTPUT_DESCRIPTIONS,
    name='Filter sequences by length.',
    description=(
        'Filter sequences by length with VSEARCH. For a combination of global '
        'and conditional taxonomic filtering, see filter_seqs_length_by_taxon.'
    ),
    citations=[citations['rognes2016vsearch']]
)


INCLUDE_SPECIES_LABELS_DESCRIPTION = (
    'Include species rank labels in taxonomy output. Note: species-labels may '
    'not be reliable in all cases.')

RANK_PROPAGATE_DESCRIPTION = (
    'If a rank has no taxonomy associated with it, the taxonomy from the '
    'upper-level rank of that lineage, will be propagated downward. For '
    'example, if we are missing the genus label for '
    '\'f__Pasteurellaceae; g__\'then the \'f__\' rank will be propagated to '
    'become: \'f__Pasteurellaceae; g__Pasteurellaceae\'.'
)

RANK_DESCRIPTION = ('List of taxonomic ranks for building a taxonomy from the '
                    'SILVA Taxonomy database. Use \'include_species_labels\' '
                    'to append the organism name as the species label. '
                    "[default: '" +
                    "', '".join(DEFAULT_RANKS) + "']")

_SILVA_VERSIONS = ['128', '132', '138', '138.1']
_SILVA_TARGETS = ['SSURef_NR99', 'SSURef', 'LSURef_NR99', 'LSURef']

version_map, target_map, _ = TypeMap({
    (Str % Choices('128', '132'),
     Str % Choices('SSURef_NR99', 'SSURef', 'LSURef')): Visualization,
    (Str % Choices('138'),
     Str % Choices('SSURef_NR99', 'SSURef')): Visualization,
    (Str % Choices('138.1'),
     Str % Choices('SSURef_NR99', 'SSURef', 'LSURef_NR99',
                   'LSURef')): Visualization,
})


plugin.pipelines.register_function(
    function=get_silva_data,
    inputs={},
    parameters={
        'version': version_map,
        'target': target_map,
        'include_species_labels': Bool,
        'rank_propagation': Bool,
        'ranks': List[Str % Choices(ALLOWED_RANKS)],
        'download_sequences': Bool},
    outputs=[('silva_sequences', FeatureData[RNASequence]),
             ('silva_taxonomy', FeatureData[Taxonomy])],
    input_descriptions={},
    parameter_descriptions={
        'version': 'SILVA database version to download.',
        'target': 'Reference sequence target to download. SSURef = redundant '
                  'small subunit reference. LSURef = redundant large subunit '
                  'reference. SSURef_NR99 = non-redundant (clustered at 99% '
                  'similarity) small subunit reference.',
        'include_species_labels': INCLUDE_SPECIES_LABELS_DESCRIPTION,
        'rank_propagation': RANK_PROPAGATE_DESCRIPTION,
        'ranks': RANK_DESCRIPTION,
        'download_sequences': 'Toggle whether or not to download and import '
                              'the SILVA reference sequences associated with '
                              'the release. Skipping the sequences is useful '
                              'if you only want to download and parse the '
                              'taxonomy, e.g., a local copy of the sequences '
                              'already exists or for testing purposes. NOTE: '
                              'if this option is used, a `silva_sequences` '
                              'output is still created, but contains no '
                              'data.'},
    output_descriptions={
        'silva_sequences': 'SILVA reference sequences.',
        'silva_taxonomy': 'SILVA reference taxonomy.'},
    name='Download, parse, and import SILVA database.',
    description=(
        'Download, parse, and import SILVA database files, given a version '
        'number and reference target. Downloads data directly from SILVA, '
        'parses the taxonomy files, and outputs ready-to-use sequence and '
        'taxonomy artifacts. REQUIRES STABLE INTERNET CONNECTION. ' +
        SILVA_LICENSE_NOTE),
    citations=[citations['Pruesse2007'], citations['Quast2013']]
)


plugin.methods.register_function(
    function=parse_silva_taxonomy,
    inputs={
        'taxonomy_tree': Phylogeny[Rooted],
        'taxonomy_map': FeatureData[SILVATaxidMap],
        'taxonomy_ranks': FeatureData[SILVATaxonomy],
    },
    parameters={
        'include_species_labels': Bool,
        'rank_propagation': Bool,
        'ranks': List[Str % Choices(ALLOWED_RANKS)]
    },
    outputs=[('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'taxonomy_tree': 'SILVA hierarchical taxonomy tree. The SILVA '
                         'release filename typically takes the form '
                         'of: \'tax_slv_ssu_X.tre\', where \'X\' is '
                         'the SILVA version number.',
        'taxonomy_map': 'SILVA taxonomy map. This file contains a mapping '
                        'of the sequence accessions to the numeric taxonomy '
                        'identifiers and species label information. '
                        'The SILVA release filename is typically in the '
                        'form of: \'taxmap_slv_ssu_ref_X.txt\', or '
                        '\'taxmap_slv_ssu_ref_nr_X.txt\' where \'X\' is '
                        'the SILVA version number.',
        'taxonomy_ranks': 'SILVA taxonomy file. This file contains the '
                         'taxonomic rank information for each numeric '
                         'taxonomy identifier and the taxonomy. The SILVA '
                         ' filename typically takes the form of: '
                         '\'tax_slv_ssu_X.txt\', where \'X\' is the SILVA '
                         'version number.',
    },
    parameter_descriptions={
        'include_species_labels': INCLUDE_SPECIES_LABELS_DESCRIPTION,
        'rank_propagation': RANK_PROPAGATE_DESCRIPTION,
        'ranks': RANK_DESCRIPTION
    },
    output_descriptions={
        'taxonomy': 'The resulting fixed-rank formatted SILVA taxonomy.'
    },
    name='Generates a SILVA fixed-rank taxonomy.',
    description=(
        'Parses several files from the SILVA reference database to produce a '
        'GreenGenes-like fixed rank taxonomy that is 6 or 7 ranks deep, '
        'depending on whether or not `include_species_labels` is applied. '
        'The generated ranks (and the rank handles used to label these '
        'ranks in the resulting taxonomy) are: domain (d__), phylum (p__), '
        'class (c__), order (o__), family (f__), genus (g__), and species '
        '(s__). ' + SILVA_LICENSE_NOTE
    ),
    citations=[citations['Pruesse2007'],
               citations['Quast2013']]
)


rt_short_description = 'Reverse transcribe RNA to DNA sequences.'
rt_inputs, rt_outputs = TypeMap({AlignedRNASequence: AlignedSequence,
                                 RNASequence: Sequence})
plugin.methods.register_function(
    function=reverse_transcribe,
    inputs={'rna_sequences': FeatureData[rt_inputs]},
    parameters={},
    outputs=[('dna_sequences', FeatureData[rt_outputs])],
    input_descriptions={
        'rna_sequences': 'RNA Sequences to reverse transcribe to DNA.'},
    parameter_descriptions={},
    output_descriptions={
        'dna_sequences': 'Reverse-transcribed DNA sequences.'},
    name=rt_short_description,
    description=(rt_short_description + ' Accepts aligned or unaligned RNA '
                 'sequences as input.')
)


GET_NCBI_DATA_PARAMS = {
    'query': Str,
    'accession_ids': Metadata,
    'ranks': List[Str % Choices(_allowed_ranks)],
    'rank_propagation': Bool,
    'logging_level': Str % Choices([
        'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
    'n_jobs': Int % Range(1, None)
}
GET_NCBI_DATA_PARAM_DESCRIPTIONS_COMMON = {
    'ranks': 'List of taxonomic ranks for building a taxonomy from the '
             "NCBI Taxonomy database. [default: '" +
             "', '".join(_default_ranks) + "']",
    'rank_propagation': 'Propagate known ranks to missing ranks if true',
    'logging_level': 'Logging level, set to INFO for download progress or '
                        'DEBUG for copious verbosity',
    'n_jobs': 'Number of concurrent download connections. More is faster '
              'until you run out of bandwidth.'
}
GET_NCBI_DATA_PARAM_DESCRIPTIONS_DNA = {
        'query': 'Query on the NCBI Nucleotide database',
        'accession_ids': 'List of accession ids for sequences in the NCBI '
                         'Nucleotide database.',
        **GET_NCBI_DATA_PARAM_DESCRIPTIONS_COMMON
    }
GET_NCBI_DATA_PARAM_DESCRIPTIONS_PROTEIN = {
        'query': 'Query on the NCBI Protein database',
        'accession_ids': 'List of accession ids for sequences in the NCBI '
                         'Protein database.',
        **GET_NCBI_DATA_PARAM_DESCRIPTIONS_COMMON
    }
GET_NCBI_DATA_DISCLAIMER = (
    '\n\nPlease be aware of the NCBI Disclaimer and Copyright notice '
    '(https://www.ncbi.nlm.nih.gov/home/about/policies/), particularly '
    '"run retrieval scripts on weekends or between 9 pm and 5 am Eastern '
    'Time weekdays for any series of more than 100 requests". As a rough '
    'guide, if you are downloading more than 125,000 sequences, only run '
    'this method at those times.\n\nThe NCBI servers can be capricious '
    'but reward polite persistence. If the download fails and gives you '
    'a message that contains the words "Last exception was ReadTimeout", '
    'you should probably try again, maybe with more connections. '
    'If it fails for any other reason, please create an issue at '
    'https://github.com/bokulich-lab/RESCRIPt.'
)

plugin.methods.register_function(
    function=get_ncbi_data,
    inputs={},
    parameters=GET_NCBI_DATA_PARAMS,
    outputs=[('sequences', FeatureData[Sequence]),
             ('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={},
    parameter_descriptions=GET_NCBI_DATA_PARAM_DESCRIPTIONS_DNA,
    output_descriptions={
        'sequences': 'Sequences from the NCBI Nucleotide database',
        'taxonomy': 'Taxonomies from the NCBI Taxonomy database'},
    name='Download, parse, and import NCBI sequences and taxonomies',
    description=(
        'Download and import sequences from the NCBI Nucleotide database '
        'and download, parse, and import the corresponding taxonomies '
        'from the NCBI Taxonomy database.' + GET_NCBI_DATA_DISCLAIMER
    ),
    citations=[citations['ncbi2018database'], citations['benson2012genbank']]
)


plugin.methods.register_function(
    function=get_ncbi_data_protein,
    inputs={},
    parameters=GET_NCBI_DATA_PARAMS,
    outputs=[('sequences', FeatureData[ProteinSequence]),
             ('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={},
    parameter_descriptions=GET_NCBI_DATA_PARAM_DESCRIPTIONS_PROTEIN,
    output_descriptions={
        'sequences': 'Sequences from the NCBI Protein database',
        'taxonomy': 'Taxonomies from the NCBI Taxonomy database'},
    name='Download, parse, and import NCBI protein sequences and taxonomies',
    description=(
        'Download and import sequences from the NCBI Protein database '
        'and download, parse, and import the corresponding taxonomies '
        'from the NCBI Taxonomy database.' + GET_NCBI_DATA_DISCLAIMER
    ),
    citations=[citations['ncbi2018database'], citations['benson2012genbank']]
)


plugin.methods.register_function(
    function=filter_taxa,
    inputs={'taxonomy': FeatureData[Taxonomy]},
    parameters={
        'ids_to_keep': Metadata,
        'include': List[Str],
        'exclude': List[Str]},
    outputs=[('filtered_taxonomy', FeatureData[Taxonomy])],
    input_descriptions={'taxonomy': 'Taxonomy to filter.'},
    parameter_descriptions={
        'ids_to_keep': 'List of IDs to keep (as Metadata). Selecting these '
                       'IDs occurs after inclusion and exclusion filtering.',
        'include': 'List of search terms. Taxa containing one or more of '
                   'these terms will be retained. Inclusion filtering occurs '
                   'prior to exclusion filtering and selecting `ids_to_keep`.',
        'exclude': 'List of search terms. Taxa containing one or more of '
                   'these terms will be excluded. Exclusion filtering occurs '
                   'after inclusion filtering and prior to selecting '
                   '`ids_to_keep`.'},
    output_descriptions={
        'filtered_taxonomy': 'The filtered taxonomy.'},
    name='Filter taxonomy by list of IDs or search criteria.',
    description=('Filter taxonomy by list of IDs or search criteria.'),
)

plugin.pipelines.register_function(
    function=trim_alignment,
    inputs={'aligned_sequences': FeatureData[AlignedSequence], },
    parameters={
        'primer_fwd': Str,
        'primer_rev': Str,
        'position_start': Int % Range(1, None),
        'position_end': Int % Range(1, None)
    },
    outputs=[('trimmed_sequences', FeatureData[AlignedSequence]), ],
    input_descriptions={'aligned_sequences': 'Aligned DNA sequences.', },
    parameter_descriptions={
        'primer_fwd': 'Forward primer used to find the start position '
                      'for alignment trimming.',
        'primer_rev': 'Reverse primer used to find the end position '
                      'for alignment trimming.',
        'position_start': 'Position within the alignment where the trimming '
                          'will begin. If not provided, alignment will not'
                          'be trimmed at the beginning. If forward primer is'
                          'specified this parameter will be ignored.',
        'position_end': 'Position within the alignment where the trimming '
                        'will end. If not provided, alignment will not be '
                        'trimmed at the end. If reverse primer is specified '
                        'this parameter will be ignored.'
    },
    output_descriptions={
        'trimmed_sequences': 'Trimmed sequence alignment.', },
    name='Trim alignment based on provided primers or specific positions.',
    description=(
        "Trim an existing alignment based on provided primers or specific, p"
        "re-defined positions. Primers take precedence over the positions,"
        "i.e. if both are provided, positions will be ignored."
        "When using primers in combination with a DNA alignment, a new "
        "alignment will be generated to locate primer positions. "
        "Subsequently, start (5'-most) and end (3'-most) position from fwd "
        "and rev primer located within the new alignment is identified and "
        "used for slicing the original alignment."),
)

T = TypeMatch([AlignedSequence, Sequence])
plugin.methods.register_function(
    function=subsample_fasta,
    inputs={'sequences': FeatureData[T]},
    parameters={
        'subsample_size':
            Float % Range(0, 1, inclusive_start=False, inclusive_end=True),
        'random_seed': Int % Range(1, None)
    },
    outputs=[('sample_sequences', FeatureData[T])],
    input_descriptions={'sequences': 'Sequences to subsample from.'},
    parameter_descriptions={
        'subsample_size': 'Size of the random sample as a '
                          'fraction of the total count',
        'random_seed': 'Seed to be used for random sampling.'
    },
    output_descriptions={
        'sample_sequences': 'Sample of original sequences.', },
    name='Subsample an indicated number of sequences from a FASTA file.',
    description=(
        "Subsample a set of sequences (either plain or aligned DNA)"
        "based on a fraction of original sequences."),
)

# Registrations
plugin.register_semantic_types(SILVATaxonomy, SILVATaxidMap)
plugin.register_semantic_type_to_format(
    FeatureData[SILVATaxonomy],
    artifact_format=SILVATaxonomyDirectoryFormat)
plugin.register_semantic_type_to_format(
    FeatureData[SILVATaxidMap],
    artifact_format=SILVATaxidMapDirectoryFormat)
plugin.register_formats(SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat,
                        SILVATaxidMapFormat, SILVATaxidMapDirectoryFormat)
importlib.import_module('rescript.types._transformer')
