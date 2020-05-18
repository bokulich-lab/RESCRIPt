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

WHITESPACE_REGEX = re.compile(r'\s+')
ALLOWED_CHARS = set('0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-[]()/.\\')
ALLOWED_RANKS = {'domain': 'd__', 'phylum': 'p__', 'class': 'c__',
                 'order': 'o__', 'family': 'f__', 'genus': 'g__'}


def _keep_allowed_chars(lin_name, allowed_chars=ALLOWED_CHARS,
                                 whitespace_pattern=WHITESPACE_REGEX):
    new_lin_name = whitespace_pattern.sub("_", lin_name.strip())
    new_lin_name = ''.join(i for i in new_lin_name if i in allowed_chars)
    return new_lin_name


def _build_base_silva_taxonomy(tree, taxrank, allowed_ranks):
    tid = {}
    for node in tree.postorder():
        if node.is_root():
            break
        rank_list = []
        rank, taxonomy = taxrank.loc[str(node.name), ['taxrank', 'taxid_taxonomy']]
        if rank in allowed_ranks.keys():
            rank_list.append((allowed_ranks[rank], taxonomy))
        for ancestor in node.ancestors():
            if ancestor.is_root():
                break
            else:
                arank, ataxonomy = taxrank.loc[str(ancestor.name), ['taxrank', 'taxid_taxonomy']]
                if arank in allowed_ranks.keys():
                    rank_list.append((allowed_ranks[arank], ataxonomy))
        tid[node.name.strip()] = dict(rank_list)
    return tid


def _prep_taxmap(taxmap):
    #taxmap = taxmap.view(pd.DataFrame)
    taxmap.index = taxmap.index + '.' + taxmap.start + '.' + taxmap.stop
    taxmap = taxmap.loc[:,['organism_name','taxid']]
    taxmap.loc[:,'organism_name'] = taxmap.loc[:,'organism_name'].apply(lambda x: _keep_allowed_chars('_'.join(x.strip().split()[:2])))
    return taxmap


def _prep_taxranks(taxrank):
    #taxrank = taxrank.view(pd.DataFrame)
    taxrank = taxrank[[1, 2]].reset_index()
    taxrank.columns = ['taxid_taxonomy', 'taxid', 'taxrank']
    taxrank = taxrank.set_index('taxid')
    taxrank.loc[:, 'taxid_taxonomy'] = taxrank.loc[:,'taxid_taxonomy'].apply(lambda x: _keep_allowed_chars(x.rsplit(';')[-2]))
    return taxrank


def parse_silva_taxonomy(taxonomy_tree: TreeNode,
                         taxonomy_map: pd.DataFrame,
                         taxonomy_ranks: pd.DataFrame,
                         include_species_labels: bool = False) -> pd.DataFrame:
    taxrank = _prep_taxranks(taxonomy_ranks)
    silva_tax_dict = _build_base_silva_taxonomy(taxonomy_tree, taxrank, ALLOWED_RANKS)
    silva_tax_id_df = pd.DataFrame.from_dict(silva_tax_dict, orient = 'index', columns=['d__', 'p__', 'c__', 'o__', 'f__', 'g__'])
    silva_tax_id_df.ffill(axis=1, inplace=True)
    silva_tax_id_df.loc[:,'full_tax_str'] = silva_tax_id_df.astype(str).apply(lambda x: x.name + x).agg('; '.join, axis=1)
    if include_species_labels:
        taxmap = _prep_taxmap(taxonomy_map)
        updated_taxmap = pd.merge(taxmap, silva_tax_id_df, left_on='taxid', right_index=True)
        updated_taxmap.loc[:,'full_tax_str_w_orgname'] = updated_taxmap.loc[:,'full_tax_str'] + '; s__' + updated_taxmap.loc[:,'organism_name']
        silva_tax_w_sp_label = updated_taxmap.loc[:,"full_tax_str_w_orgname"].to_frame()
        silva_tax_w_sp_label = pd.DataFrame(silva_tax_w_sp_label, index=['Feature ID'], columns = ['Taxonomy'])
        return silva_tax_w_sp_label
    else:
        silva_tax_wo_sp_label = silva_tax_id_df.loc[:,"full_tax_str"].to_frame()
        silva_tax_w_sp_label = pd.DataFrame(silva_tax_wo_sp_label, index=['Feature ID'], columns = ['Taxonomy'])
        return silva_tax_wo_sp_label
