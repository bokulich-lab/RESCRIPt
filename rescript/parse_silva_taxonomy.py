# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.tree import TreeNode
import re
from q2_types.feature_data import Taxonomy
import rescript
from rescript.types._format import (
    SILVATaxonomyFormat, SILVATaxonomyDirectoryFormat, SILVATaxidMapFormat,
    SILVATaxidMapDirectoryFormat)
from rescript.types._type import SILVATaxonomy, SILVATaxidMap
import pandas as pd
from collections import OrderedDict

WHITESPACE_REGEX = re.compile(r'\s+')
ALLOWED_CHARS = set(''.join(['0123456789',
                             'abcdefghijklmnopqrstuvwxyz',
                             'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
                             '_-[]()/.\\']))

# Do not use "'major_clade': 'mc__'" as this appears at multiple levels, even
#      within the same lineage. Tough to disambiguate.
ALLOWED_RANKS = OrderedDict({'domain': 'd__', 'superkingdom': 'sk__',
                'kingdom': 'k__', 'subkingdom': 'ks__', 'superphylum': 'sp__',
                'phylum': 'p__', 'subphylum': 'ps__', 'infraphylum': 'pi__',
                'superclass': 'sc__', 'class': 'c__', 'subclass': 'cs__',
                'infraclass': 'ci__', 'superorder': 'so__', 'order': 'o__',
                'suborder': 'os__', 'superfamily': 'sf__', 'family': 'f__',
                'subfamily': 'fs__', 'genus': 'g__'})

SELECTED_RANKS = OrderedDict({'domain': 'd__', 'phylum': 'p__', 'class': 'c__',
                             'order': 'o__', 'family': 'f__', 'genus': 'g__'})


def _keep_allowed_chars(lin_name, allowed_chars = ALLOWED_CHARS,
                        whitespace_pattern = WHITESPACE_REGEX):
    # Only keep "safe" characters and replace whitespace with '_'
    new_lin_name = whitespace_pattern.sub("_", lin_name.strip())
    new_lin_name = ''.join(i for i in new_lin_name if i in allowed_chars)
    return new_lin_name


def _build_base_silva_taxonomy(tree, taxrank, allowed_ranks):
    # Return base silva taxonomy dataframe by traversing taxonomy tree
    # and forward filling the lower ranks with upper-level taxonomy, then
    # pre-pend the rank labels. May be able to forgo traversing taxonomy tree
    # in the future and parse the taxrank file directly.
    tid = {}
    for node in tree.postorder():
        if node.is_root():
            break
        rank_list = []

        # should probably add a try statement here
        # to catch the case for when a node in the tree is
        # not dound in the taxrank dataframe.
        rank, taxonomy = taxrank.loc[str(node.name), ['taxrank',
                                                      'taxid_taxonomy']]
        # only pull taxonomy from allowed ranks
        if rank in allowed_ranks.keys():
            rank_list.append((allowed_ranks[rank], taxonomy))
        for ancestor in node.ancestors():
            if ancestor.is_root():
                break
            else:
                arank, ataxonomy = taxrank.loc[str(ancestor.name),
                                              ['taxrank', 'taxid_taxonomy']]
                if arank in allowed_ranks.keys():
                    rank_list.append((allowed_ranks[arank], ataxonomy))
        tid[node.name.strip()] = dict(rank_list)
    silva_tax_id_df = pd.DataFrame.from_dict(tid, orient = 'index',
                                             columns=ALLOWED_RANKS.values())
    silva_tax_id_df.index.name = 'taxid'
    # forward fill missing taxonomy with upper level taxonomy
    silva_tax_id_df.ffill(axis=1, inplace=True)
    return silva_tax_id_df


def _prep_taxmap(taxmap):
    # Return silva taxonomy dataframe map indexed by accession.seqtart.seqstop
    # This is how the coreesponding FASTA file IDs are structured
    taxmap.index = taxmap.index + '.' + taxmap.start + '.' + taxmap.stop
    taxmap.index.name = 'Feature ID'
    taxmap = taxmap.loc[:,['organism_name','taxid']]
    # parse organism name (i.e. "species") by only returning the first two
    # words (many species labels are very long or contain sub-species info).
    taxmap.loc[:,'organism_name'] = taxmap.loc[:,'organism_name'].apply(
                lambda x: _keep_allowed_chars('_'.join(x.strip().split()[:2])))
    return taxmap


def _prep_taxranks(taxrank):
    # return only taxonomy rank information indexed by taxonomy ID
    taxrank = taxrank[[1, 2]].reset_index()
    taxrank.columns = ['taxid_taxonomy', 'taxid', 'taxrank']
    taxrank.set_index('taxid', inplace=True)
    # return only the last taxonomy from a given taxonomy path
    taxrank.loc[:, 'taxid_taxonomy'] = taxrank.loc[:,'taxid_taxonomy'].apply(
                            lambda x: _keep_allowed_chars(x.rsplit(';')[-2]))
    return taxrank


def _compile_taxonomy_output(updated_taxmap, include_species_labels=False,
                            selected_ranks=SELECTED_RANKS):
    sr = list(selected_ranks.values())
    # add rank prefixes (i.e. 'p__')
    updated_taxmap.loc[:, sr] = updated_taxmap.loc[:, sr].apply(
                                                    lambda x: x.name + x)
    if include_species_labels:
        updated_taxmap.loc[:, 'organism_name'] = updated_taxmap.loc[
                                :,'organism_name'].apply(lambda x: 's__' + x)
        sr.append('organism_name')
        taxonomy = updated_taxmap.loc[:, sr].agg('; '.join, axis=1)
        taxonomy.rename('Taxon', inplace=True)
        return taxonomy
    else:
        taxonomy = updated_taxmap.loc[:, sr].agg('; '.join, axis=1)
        taxonomy.rename('Taxon', inplace=True)
        return taxonomy

def parse_silva_taxonomy(taxonomy_tree: TreeNode,
                         taxonomy_map: pd.DataFrame,
                         taxonomy_ranks: pd.DataFrame,
                         include_species_labels: bool = False) -> pd.Series:
    taxrank = _prep_taxranks(taxonomy_ranks)
    silva_tax_id_df = _build_base_silva_taxonomy(taxonomy_tree, taxrank,
                                                 ALLOWED_RANKS)
    taxmap = _prep_taxmap(taxonomy_map)
    updated_taxmap = pd.merge(taxmap, silva_tax_id_df, left_on='taxid',
                              right_index=True)
    taxonomy = _compile_taxonomy_output(updated_taxmap,
                                        include_species_labels,
                                        SELECTED_RANKS)
    return taxonomy
