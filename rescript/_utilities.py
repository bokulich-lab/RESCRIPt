# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import zip_longest, takewhile
from re import sub
import subprocess
import skbio
from collections import Counter
from q2_types.feature_data import DNAFASTAFormat, AlignedDNAFASTAFormat


_rank_handles = {
    'silva': [' d__', ' p__', ' c__', ' o__', ' f__', ' g__', ' s__'],
    'greengenes': ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'],
    'gtdb': ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'],
    'disable': None,
}


# modified version of _find_lca that prioritizes majority and supersets
def _find_super_lca(taxa, collapse_substrings=True):
    # collapse and count unique labels at each rank
    # yields list of ('labels', counts) sorted by most to least abundant
    taxa = [[t for t in r if t not in [None, '']] for r in zip_longest(*taxa)]
    if collapse_substrings:
        # find longest string in group of sub/superstrings, combine
        taxa = [[max({i for i in x if t in i}, key=len) for t in x]
                if x else '' for x in taxa]
    taxa_comparison = [Counter(t).most_common() for t in taxa]
    # return majority wherever a clear majority is found
    # terminate when no majority is found, that's your LCA
    # propagate empty ranks that maintain majority/consensus by inserting ''
    # to preserve empty ranks in the original taxonomies (e.g., when using a
    # rank handle unannotated levels will be removed, but as long as these pass
    # muster we want to preserve those labels)
    return [rank[0][0] if rank else '' for rank in takewhile(
        lambda x: len(x) < 2 or x[0][1] > x[1][1], taxa_comparison)]


# LCA majority is same as super majority without substring collapsing
def _find_lca_majority(taxa):
    return _find_super_lca(taxa, collapse_substrings=False)


def _find_lca(taxa):
    '''Find least common ancestor between two semicolon-delimited strings.'''
    # determine optimal zip mode. Normally, zip is best because trimming to
    # shortest is an inherent feature of LCA.
    # However if only one frame contains an assignment for feature x, we want
    # to just take that taxonomy. zip_longest will accomplish this while using
    # the same machinery... same deal if we want to take a superset taxonomy
    if '' in taxa:
        zip_it = zip_longest
    else:
        zip_it = zip
    # LCA ends where zipped taxonomy strings no longer converge to len == 1
    taxa_comparison = [set(rank) - {None} for rank in zip_it(*taxa)]
    return (rank.pop() for rank in takewhile(
        lambda x: len(x) == 1, taxa_comparison))


def _rank_length(t1, t2):
    '''Determine which semicolon-delimited string has more elements.'''
    len_first = len(set(t1['Taxon']) - {None, ''})
    len_second = len(set(t2['Taxon']) - {None, ''})
    return t1 if len_first > len_second else t2


def _find_top_score(t1, t2):
    return t1 if t1['Score'] > t2['Score'] else t2


def _taxon_to_list(taxon, rank_handle):
    '''Split taxonomy string into list of taxonomic labels'''
    if rank_handle != '':
        return [sub(rank_handle, '', t.strip()) for t in taxon.split(';')]
    else:
        return [t.strip() for t in taxon.split(';')]


def _majority(taxon):
    '''Find most common element in a list'''
    return max(set(taxon), key=taxon.count)


def run_command(cmd, verbose=True):
    print("Running external command line application. This may print "
          "messages to stdout and/or stderr.")
    print("The command being run is below. This command cannot "
          "be manually re-run as it will depend on temporary files that "
          "no longer exist.")
    print("\nCommand:", end=' ')
    print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def _read_fasta(path, constructor=skbio.DNA):
    return skbio.read(path, format='fasta', constructor=constructor)


def _read_dna_alignment_fasta(path):
    return skbio.TabularMSA.read(path, format='fasta', constructor=skbio.DNA)


def _rna_to_dna(iterator):
    ff = DNAFASTAFormat()
    generator = _rna_to_dna_iterator(iterator)
    skbio.io.write(iter(generator), format='fasta', into=str(ff))
    return ff


def _rna_to_dna_iterator(iterator):
    for seq in iterator:
        yield seq.reverse_transcribe()


def _dna_iterator_to_aligned_fasta(iterator):
    ff = AlignedDNAFASTAFormat()
    skbio.io.write(iter(iterator), format='fasta', into=str(ff))
    return ff
