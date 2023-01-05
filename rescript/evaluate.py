# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from multiprocessing import cpu_count

import pandas as pd
import numpy as np
import scipy
import warnings
import qiime2 as q2
import q2templates
from os.path import join

from joblib import Parallel, delayed
from q2_types.feature_data import DNAIterator
from sklearn.feature_extraction.text import HashingVectorizer
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from itertools import zip_longest
import pkg_resources

from ._utilities import _taxon_to_list

TEMPLATES = pkg_resources.resource_filename('rescript', 'assets')


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


def _process_labels(labels, taxonomies, description='taxonomies'):
    if labels is None:
        labels = []
    n_labels, n_taxonomies = len(labels), len(taxonomies)
    # warn if fewer labels than taxonomies
    if n_labels < n_taxonomies:
        _warn_uneven_length(description)
        labels = [name if name is not None else n for n, (name, t)
                  in enumerate(zip_longest(labels, taxonomies), 1)]
    # trim labels to number of Taxonomies
    elif n_labels > n_taxonomies:
        labels = labels[:n_taxonomies]
    return labels


def _warn_uneven_length(description):
    msg = ("The lists of input {0} and labels are different lengths. "
           "Additional {0} will be labeled numerically by their order in the "
           "inputs. Note that if these numbers match existing labels, those "
           "data will be grouped in the visualization.".format(description))
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


def evaluate_seqs(output_dir: str, sequences: DNAIterator, labels: list = None,
                  kmer_lengths: list = None, subsample_kmers: float = 1.0,
                  palette: str = 'viridis') -> None:
    labels = _process_labels(labels, sequences, description='sequences')
    results, lengths = _evaluate_seqs(
        sequences, labels, kmer_lengths, subsample_kmers)
    fig = _plot_eval_seqs(results, lengths, palette)
    _visualize(output_dir, results, fig)
    plt.close('all')


def _process_kmers(seqs, kmer_len, subsample_kmers):
    if subsample_kmers < 1:
        subsample_size = int(subsample_kmers * len(seqs))
        seqs = np.random.choice(seqs, subsample_size, replace=False)
    vectorizer = HashingVectorizer(
        alternate_sign=False,
        analyzer='char',
        ngram_range=[kmer_len, kmer_len])
    X = vectorizer.fit_transform(seqs)
    kmer_freq = X.sum(axis=0)
    return scipy.stats.entropy(kmer_freq, axis=1)[0]


def _evaluate_seqs(sequences, labels, kmer_lengths=None, subsample_kmers=1.0):
    if not kmer_lengths:
        kmer_lengths = []

    rownames = ['Length ' + n for n in ['min', '1%', '25%', 'median', '75%',
                                        '99%', 'max']]
    rownames += ['N uniques', 'Sequence Entropy']
    rownames += ['%smer Entropy' % k for k in kmer_lengths]

    lengths = []
    results = pd.DataFrame(index=rownames)
    for n, seqs in zip(labels, sequences):
        seqs = [str(s) for s in seqs]
        length_array = np.array([len(s) for s in seqs])
        lengths.append(length_array)
        len_quantiles = np.nanquantile(
            length_array, [0, 0.01, 0.25, 0.5, 0.75, 0.99, 1])
        uniqs = Counter(seqs).values()
        seq_entropy = scipy.stats.entropy(list(uniqs))
        res = list(len_quantiles) + [len(uniqs), seq_entropy]

        n_jobs = min(len(kmer_lengths), cpu_count())
        if n_jobs > 0:
            parallel = Parallel(n_jobs=n_jobs, backend='loky')
            kmer_freqs = parallel(
                delayed(_process_kmers)(
                    seqs, k, subsample_kmers) for k in kmer_lengths)
            res.extend(kmer_freqs)
        results[n] = res
    return results.round(2), lengths


def _plot_eval_seqs(results, lengths, palette):
    n_groups = len(lengths)
    cmap = sns.color_palette(palette, n_groups)

    # adjust figure width based on number of inputs — minimum 7 inches
    fig_width = 3 * n_groups + 3
    fig = plt.figure(constrained_layout=True, figsize=(fig_width, 6))
    gs1 = fig.add_gridspec(
        nrows=2, ncols=4, left=0.05, right=0.48, wspace=0.05)

    # length distribution histogram
    ax = fig.add_subplot(gs1[:-1, :])
    for dat, color in zip(lengths, cmap):
        sns.kdeplot(dat, shade=True, color=color, ax=ax, alpha=0.5)
    ax.set_title('Sequence Length Distribution')
    ax.set_ylabel('Proportion')
    ax.set_xlabel('Length (nt)')
    ax.set_xlim(results.loc['Length 1%'].min(),
                results.loc['Length 99%'].max())

    # barplot of seq counts
    ax = fig.add_subplot(gs1[-1, :-2])
    sns.barplot(data=results.loc[['N uniques']], palette=palette, ax=ax)
    ax.set_title('N Unique Sequences')
    ax.set_ylabel('Count')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')

    # lineplot/barplot of seq/kmer entropy
    ax = fig.add_subplot(gs1[-1, -2:])
    if len(results) > 9:
        data = results[8:]
        data.index = [i.split(' ')[0] for i in data.index]
        sns.lineplot(data=data, sort=False, palette=palette, ax=ax)
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        # plot barplot if kmer entropy is not calculated
        sns.barplot(data=results.loc[['Sequence Entropy']],
                    palette=palette, ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='right')
    ax.set_title('Entropy')
    ax.set_ylabel('H')

    return fig


def _visualize(output_dir, results, plot):

    pd.set_option('display.max_colwidth', -1)

    # save results
    results.to_csv(join(output_dir, 'evaluate_seqs_results.tsv'), sep='\t')
    results = q2templates.df_to_html(results, index=True)

    plot.savefig(join(output_dir, 'evaluate_seqs.png'), bbox_inches='tight')
    plot.savefig(join(output_dir, 'evaluate_seqs.pdf'), bbox_inches='tight')

    index = join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'title': 'Sequence Evaluation Results',
        'running_title': 'evaluate_seqs',
        'results': results,
    })
