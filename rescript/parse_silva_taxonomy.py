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
ALLOWED_RANKS = OrderedDict({'domain': 'd__', 'phylum': 'p__', 'class': 'c__',
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
    # forward fill missing taxonomy with upper level taxonomy
    silva_tax_id_df.ffill(axis=1, inplace=True)
    # add rank prefixes (i.e. 'p__')
    silva_tax_id_df.loc[:,'full_tax_str'] = silva_tax_id_df.astype(
                        str).apply(lambda x: x.name + x).agg('; '.join, axis=1)
    return silva_tax_id_df


def _prep_taxmap(taxmap):
    # Return silva taxonomy dataframe map indexed by accession.seqtart.seqstop
    # This is how the coreesponding FASTA file IDs are structured
    taxmap.index = taxmap.index + '.' + taxmap.start + '.' + taxmap.stop
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
    if include_species_labels:
        updated_taxmap.loc[:,'full_tax_str_w_orgname'] = updated_taxmap.loc[
            :,'full_tax_str'] + '; s__' + updated_taxmap.loc[:,'organism_name']
        taxonomy_w_sp_label = updated_taxmap['full_tax_str_w_orgname']
        # rename indices for FeatureData[Taxonomy] format
        taxonomy_w_sp_label.rename('Taxon', inplace=True)
        taxonomy_w_sp_label.index.name = 'Feature ID'
        return taxonomy_w_sp_label
    else:
        taxonomy_wo_sp_label = updated_taxmap["full_tax_str"]
        # rename indices for FeatureData[Taxonomy] format
        taxonomy_wo_sp_label.rename('Taxon', inplace=True)
        taxonomy_wo_sp_label.index.name = 'Feature ID'
        return taxonomy_wo_sp_label
