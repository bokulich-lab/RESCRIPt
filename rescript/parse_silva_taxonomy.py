# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.tree import TreeNode
import re
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
                             'kingdom': 'k__', 'subkingdom': 'ks__',
                             'superphylum': 'sp__', 'phylum': 'p__',
                             'subphylum': 'ps__', 'infraphylum': 'pi__',
                             'superclass': 'sc__', 'class': 'c__',
                             'subclass': 'cs__', 'infraclass': 'ci__',
                             'superorder': 'so__', 'order': 'o__',
                             'suborder': 'os__', 'superfamily': 'sf__',
                             'family': 'f__', 'subfamily': 'fs__',
                             'genus': 'g__'})

SELECTED_RANKS = OrderedDict({'domain': 'd__', 'phylum': 'p__',
                              'class': 'c__',  'order': 'o__',
                              'family': 'f__', 'genus': 'g__'})


def _keep_allowed_chars(lin_name, allowed_chars=ALLOWED_CHARS,
                        whitespace_pattern=WHITESPACE_REGEX):
    # Only keep "safe" characters and replace whitespace with '_'
    new_lin_name = whitespace_pattern.sub("_", lin_name.strip())
    new_lin_name = ''.join(i for i in new_lin_name if i in allowed_chars)
    return new_lin_name


def _build_base_silva_taxonomy(tree, taxrank, allowed_ranks):
    # Return base silva taxonomy dataframe by traversing taxonomy tree
    # and forward filling the lower ranks with upper-level taxonomy, then
    # pre-pend the rank labels. May be able to forgo traversing taxonomy tree
    # in the future and parse the taxrank file directly.
    # tree : skbio.TreeNode
    # taxrank : output from _prep_taxranks
    # allowed_ranks : Ordered dict {'domain': 'd__', ...} of allowed ranks to
    # capture and parse
    # output : pandas DataFrame with 'taxid' as the index and the full taxonomy
    # path (i.e. domain-to-genus) as the values.
    tid = {}
    for node in tree.postorder():
        if node.is_root():
            break
        rank_list = []

        # should probably add a try statement here?
        # though the _validate_taxrank_taxtree func should run
        # before getting here
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
    silva_tax_id_df = pd.DataFrame.from_dict(tid, orient='index',
                                             columns=ALLOWED_RANKS.values())
    silva_tax_id_df.index.name = 'taxid'
    # forward fill missing taxonomy with upper level taxonomy
    silva_tax_id_df.ffill(axis=1, inplace=True)
    return silva_tax_id_df


def _get_clean_organism_name(name):
    # parse organism name (i.e. "species") by only returning the first two
    # words (many species labels are very long or contain sub-species info).
    return _keep_allowed_chars('_'.join(name.strip().split()[:2]))


def _prep_taxmap(taxmap):
    # Return silva taxonomy dataframe map indexed by accession.seqtart.seqstop
    # This is how the coreesponding FASTA file IDs are structured
    taxmap.index = taxmap.index + '.' + taxmap.start + '.' + taxmap.stop
    taxmap.index.name = 'Feature ID'
    taxmap = taxmap[['organism_name', 'taxid']]
    taxmap.loc[:, 'organism_name'] = taxmap['organism_name'].apply(
                                            _get_clean_organism_name)
    return taxmap


def _get_terminal_taxon(name):
    return _keep_allowed_chars(name.rsplit(';')[-2])


def _prep_taxranks(taxrank):
    # return only taxonomy rank information indexed by taxonomy ID
    taxrank = taxrank[[1, 2]].reset_index()
    taxrank.columns = ['taxid_taxonomy', 'taxid', 'taxrank']
    taxrank.set_index('taxid', inplace=True)
    # return only the last taxonomy from a given taxonomy path
    taxrank.loc[:, 'taxid_taxonomy'] = taxrank.loc[:, 'taxid_taxonomy'].apply(
                                            lambda x: _get_terminal_taxon(x))
    return taxrank


def _validate_taxrank_taxtree(prepped_taxrank, taxtree):
    # prepass to make sure there is 100 % agreement with the taxids
    # present within the taxranks and taxtree files.
    tree_taxids = {node.name for node in taxtree.postorder()}
    print(tree_taxids)
    # '1' is the root of the tree of life,
    # which has no tax rank info in any other files.
    tree_taxids.remove('1')
    ptrs = set(prepped_taxrank.index.unique())
    if tree_taxids != ptrs:
        diffs = tree_taxids.symmetric_difference(ptrs)
        d_str = ', '.join(diffs)
        raise ValueError("The taxids of the SILVA Taxonomy Tree file "
                         "and the Taxonomy Ranks file do not match! "
                         "The values that are missing in at least one or "
                         "the other files are: %s" % d_str)


def _compile_taxonomy_output(updated_taxmap, include_species_labels=False,
                             selected_ranks=SELECTED_RANKS):
    # Takes the updated_taxmap dataframe (i.e. the combined
    # output from _build_base_silva_taxonomy and _prep_taxmap) and only
    # returns the desired ranks (as a pandas Series) from selected_ranks, which
    # is an ordered dict. For now, SELECTED_RANKS is used. In the future, we
    # may provide the option to allow users to feed in their own dictionary
    # of selected_ranks. The user can opt to return the species labels.
    sr = list(selected_ranks.values())
    # add rank prefixes (i.e. 'p__')
    updated_taxmap.loc[:, sr] = updated_taxmap.loc[:, sr].apply(
                                                    lambda x: x.name + x)
    if include_species_labels:
        updated_taxmap.loc[:, 'organism_name'] = updated_taxmap.loc[
                           :, 'organism_name'].apply(lambda x: 's__' + x)
        sr.append('organism_name')
    taxonomy = updated_taxmap.loc[:, sr].agg('; '.join, axis=1)
    taxonomy.rename('Taxon', inplace=True)
    return taxonomy


def parse_silva_taxonomy(taxonomy_tree: TreeNode,
                         taxonomy_map: pd.DataFrame,
                         taxonomy_ranks: pd.DataFrame,
                         include_species_labels: bool = False) -> pd.Series:
    # Traverse the taxnomy hierarchy tree (taxonomy_tree) to obtain the
    # taxids. These will be used to look up the taxonomy and rank information
    # from the taxonomy_ranks file. Finally the taxonomy information is
    # mapped to each Accesioned sequence via the taxonomy_map. An option
    # to include the, potentially untrustworthy, species labels is provided.
    taxrank = _prep_taxranks(taxonomy_ranks)
    _validate_taxrank_taxtree(taxrank, taxonomy_tree)
    silva_tax_id_df = _build_base_silva_taxonomy(taxonomy_tree, taxrank,
                                                 ALLOWED_RANKS)
    taxmap = _prep_taxmap(taxonomy_map)
    updated_taxmap = pd.merge(taxmap, silva_tax_id_df, left_on='taxid',
                              right_index=True)
    taxonomy = _compile_taxonomy_output(updated_taxmap,
                                        include_species_labels,
                                        SELECTED_RANKS)
    return taxonomy
