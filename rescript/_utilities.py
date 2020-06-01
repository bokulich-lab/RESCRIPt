# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import zip_longest, takewhile
from re import sub
import subprocess
import skbio
from q2_types.feature_data import DNAFASTAFormat


def _find_lca_in_series(t1, t2):
    '''
    Find least common ancestor between two semicolon-delimited strings.
    Input consists of two series containing index "Taxon".
    Returns the LCA of those taxonomies as an iterator.
    '''
    t1 = t1['Taxon']
    t2 = t2['Taxon']
    return _find_lca([t1, t2])


def _find_lca(taxa):
    '''Find least common ancestor between two semicolon-delimited strings.'''
    # determine optimal zip mode. Normally, zip is best because trimming to
    # shortest is an inherent feature of LCA.
    # However if only one frame contains an assignment for feature x, we want
    # to just take that taxonomy. zip_longest will accomplish this while using
    # the same machinery...
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
    return t1 if t1['score'] > t2['score'] else t2


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


def _read_rna_fasta(path):
    return skbio.read(path, format='fasta', constructor=skbio.RNA)


def _read_dna_fasta(path):
    return skbio.read(path, format='fasta', constructor=skbio.DNA)


def _rna_to_dna(path):
    ff = DNAFASTAFormat()
    with ff.open() as outfasta:
        for seq in _read_rna_fasta(path):
            seq.reverse_transcribe().write(outfasta)
    return ff
