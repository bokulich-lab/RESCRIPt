# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import (Str, Plugin, Choices, List, Citations, Range, Int,
                           Float, Visualization, Bool, TypeMap)
from .merge import merge_taxa
from .dereplicate import dereplicate
from .evaluate import evaluate_taxonomy
from .screenseq import screen_sequences
from .degap import degap_sequences
from .parse_silva_taxonomy import parse_silva_taxonomy
from .get_data import get_silva_data
from .cross_validate import cross_validate, evaluate_classifications
from .filter_length import filter_seqs_length_by_taxon, filter_seqs_length
from .orient import orient_seqs
from q2_types.feature_data import (FeatureData, Taxonomy, Sequence,
                                   AlignedSequence)
from q2_types.tree import Phylogeny, Rooted
from q2_feature_classifier.classifier import (_parameter_descriptions,
                                              _classify_parameters)

import rescript
from rescript._utilities import _rank_handles
from rescript.types._format import (
    SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat, SILVATaxidMapFormat,
    SILVATaxidMapDirectoryFormat, RNAFASTAFormat, RNASequencesDirectoryFormat)
from rescript.types._type import SILVATaxonomy, SILVATaxidMap, RNASequence
from rescript.types.methods import reverse_transcribe


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
    'For example, "^[dkpcofgs]__" will recognize greengenes or silva rank '
    'handles. ')

rank_handle_extra_note = (
    'Note that rank_handles are removed but not replaced; use the '
    'new_rank_handle parameter to replace the rank handles.')

labels_description = (
    'List of labels to use for labeling taxonomic results in the resulting '
    'visualization. Taxonomies are labeled with labels in the order that each '
    'is input. If there are fewer labels than taxonomies (or no labels), '
    'unnamed taxonomies are labeled numerically in sequential order. Extra '
    'labels are ignored.')


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
        'mode': Str % Choices(['len', 'lca', 'score']),
        'rank_handle_regex': Str,
        'new_rank_handle': Str % Choices(list(_rank_handles.keys()))},
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
            'inputs.')},
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
        'mode': Str % Choices(['uniq', 'lca', 'majority']),
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
                'pick the winner arbitrarily.',
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

plugin.methods.register_function(
    function=screen_sequences,
    inputs={
        'sequences': FeatureData[Sequence | RNASequence]
        },
    parameters={
        'num_degenerates': Int % Range(1, None),
        'homopolymer_length': Int % Range(2, None)
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
    function=degap_sequences,
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

_SILVA_VERSIONS = ['128', '132', '138']
_SILVA_TARGETS = ['SSURef_NR99', 'SSURef', 'LSURef']

version_map, target_map, _ = TypeMap({
    (Str % Choices('128', '132'),
     Str % Choices('SSURef_NR99', 'SSURef', 'LSURef')): Visualization,
    (Str % Choices('138'), Str % Choices('SSURef_NR99', 'SSURef')):
        Visualization,
})


plugin.pipelines.register_function(
    function=get_silva_data,
    inputs={},
    parameters={
        'version': version_map,
        'target': target_map,
        'include_species_labels': Bool,
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
        'taxonomy artifacts. REQUIRES STABLE INTERNET CONNECTION.'),
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
        'include_species_labels': Bool
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
        '(s__).'
        ),
    citations=[citations['Pruesse2007'],
               citations['Quast2013']]
)


plugin.methods.register_function(
    function=reverse_transcribe,
    inputs={'rna_sequences': FeatureData[RNASequence]},
    parameters={},
    outputs=[('dna_sequences', FeatureData[Sequence])],
    input_descriptions={
        'rna_sequences': 'RNA Sequences to reverse transcribe to DNA.'},
    parameter_descriptions={},
    output_descriptions={
        'dna_sequences': 'Reverse-transcribed DNA sequences.'},
    name='Reverse transcribe RNA to DNA sequences.',
    description=('Reverse transcribe RNA to DNA sequences.')
)


# Registrations
plugin.register_semantic_types(SILVATaxonomy, SILVATaxidMap, RNASequence)
plugin.register_semantic_type_to_format(
    FeatureData[SILVATaxonomy],
    artifact_format=SILVATaxonomyDirectoryFormat)
plugin.register_semantic_type_to_format(
    FeatureData[SILVATaxidMap],
    artifact_format=SILVATaxidMapDirectoryFormat)
plugin.register_semantic_type_to_format(
    FeatureData[RNASequence],
    artifact_format=RNASequencesDirectoryFormat)
plugin.register_formats(SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat,
                        SILVATaxidMapFormat, SILVATaxidMapDirectoryFormat,
                        RNAFASTAFormat, RNASequencesDirectoryFormat)
importlib.import_module('rescript.types._transformer')
