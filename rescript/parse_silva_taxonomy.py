# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
import pandas as pd
from skbio.tree import TreeNode
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

DEFAULT_RANKS = ['domain', 'phylum', 'class',  'order', 'family', 'genus']


def _keep_allowed_chars(lin_name, allowed_chars=ALLOWED_CHARS,
                        whitespace_pattern=WHITESPACE_REGEX):
    # Only keep "safe" characters and replace whitespace with '_'
    new_lin_name = whitespace_pattern.sub("_", lin_name.strip())
    new_lin_name = ''.join(i for i in new_lin_name if i in allowed_chars)
    return new_lin_name


def _build_base_silva_taxonomy(tree, taxrank, allowed_ranks,
                               rank_propagation):
    # Return base silva taxonomy dataframe by traversing taxonomy tree
    # and forward filling the lower ranks with upper-level taxonomy, then
    # pre-pend the rank labels.
    # tree : skbio.TreeNode
    # taxrank : output from _prep_taxranks
    # allowed_ranks : Ordered dict {'domain': 'd__', ...} of allowed ranks to
    # capture and parse
    # output : pandas DataFrame with 'taxid' as the index and the full taxonomy
    # path (i.e. domain-to-genus) as the values.
    # rank_propagation : boolean value to determine if empty taxonomy
    # ranks should be filled from upper-level taxonomies.
    tid = {}
    for node in tree.postorder():
        if node.is_root():
            break
        rank_list = []
        # the _validate_taxrank_taxtree func should run
        # before getting here, so all looks up should pass
        rank, taxonomy = taxrank.loc[str(node.name), ['taxrank',
                                                      'taxid_taxonomy']]
        # only pull taxonomy from allowed ranks
        if rank in allowed_ranks:
            rank_list.append((allowed_ranks[rank], taxonomy))
        for ancestor in node.ancestors():
            if ancestor.is_root():
                break
            else:
                arank, ataxonomy = taxrank.loc[str(ancestor.name),
                                               ['taxrank', 'taxid_taxonomy']]
                if arank in allowed_ranks:
                    rank_list.append((allowed_ranks[arank], ataxonomy))
        tid[node.name.strip()] = dict(rank_list)
    silva_tax_id_df = pd.DataFrame.from_dict(tid, orient='index',
                                             columns=allowed_ranks.values())
    silva_tax_id_df.index.name = 'taxid'
    # forward fill missing taxonomy with upper level taxonomy
    if rank_propagation:
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
    taxmap = taxmap.loc[:, ['organism_name', 'taxid']]
    taxmap.loc[:, 'organism_name'] = taxmap.loc[:, 'organism_name'].apply(
        _get_clean_organism_name)
    return taxmap


def _get_terminal_taxon(name):
    # Return the last taxonomy of a taxonomy path. That is
    # return: 'Aenigmarchaeales'
    # from: 'Archaea;Aenigmarchaeota;Aenigmarchaeia;Aenigmarchaeales;'
    return _keep_allowed_chars(name.rsplit(';')[-2])


def _prep_taxranks(taxrank):
    # return only taxonomy rank information indexed by taxonomy ID
    taxrank = taxrank[[1, 2]].reset_index()
    taxrank.columns = ['taxid_taxonomy', 'taxid', 'taxrank']
    taxrank.set_index('taxid', inplace=True)
    # return only the last taxonomy from a given taxonomy path
    taxrank.loc[:, 'taxid_taxonomy'] = taxrank.loc[:, 'taxid_taxonomy'].apply(
        _get_terminal_taxon)
    return taxrank


def _validate_taxrank_taxmap_taxtree(prepped_taxrank, prepped_taxmap, taxtree):
    # Prepass to make sure there is 100 % agreement with the taxids
    # present within the taxranks and taxtree files. Then check that
    # the taxids in the tree are a superset of the taxmap, i.e. not necessarily
    # all taxids are represented by accessioned sequences.
    tree_taxids = {node.name for node in taxtree.postorder()}
    # '1' is the root of the tree of life,
    # which has no tax rank info in any other files.
    tree_taxids.remove('1')  # tree taxids
    ptrs = set(prepped_taxrank.index.unique())  # taxrand taxids
    ptms = set(prepped_taxmap.taxid.unique())  # taxmap taxids

    tree_rank_diffs = tree_taxids.symmetric_difference(ptrs)
    if any(tree_rank_diffs):
        raise ValueError("The taxids of the SILVA Taxonomy Tree file "
                         "and the Taxonomy Ranks file do not match! "
                         "Please check that you are using the same SILVA "
                         "release versions of your files! The taxids that "
                         "are missing in at least one or the other files "
                         "are: ", tree_rank_diffs)
    tree_map_diffs = ptms.difference(tree_taxids)
    if len(tree_map_diffs) > 0:
        raise ValueError("The SILVA Taxonomy Map file conains taxids not "
                         "present within the the SILVA Taxonomy Tree file! "
                         "Please check that you are using the same SILVA "
                         "release versions of your files! The missing "
                         "taxids are: ", tree_map_diffs)


def _compile_taxonomy_output(updated_taxmap, ranks,
                             include_species_labels=False):
    # Takes the updated_taxmap dataframe (i.e. the combined
    # output from _build_base_silva_taxonomy and _prep_taxmap) and only
    # returns the desired ranks (as a pandas Series).
    # Sort user ranks relative to allowed_ranks
    # this ensures that the taxonomic ranks supplied by the user are in order
    sorted_ranks = [p for r, p in ALLOWED_RANKS.items() if r in ranks]
    # Must replace NaNs with empty string prior to the following
    # operations, or we'll get a "expected str instance, float found" or
    # similar error.
    updated_taxmap.fillna('', inplace=True)
    # add rank prefixes (i.e. 'p__')
    updated_taxmap.loc[:, sorted_ranks] = \
        updated_taxmap.loc[:, sorted_ranks].apply(lambda x: x.name + x)
    if include_species_labels:
        updated_taxmap.loc[:, 'organism_name'] = updated_taxmap.loc[
            :, 'organism_name'].apply(lambda x: 's__' + x)
        sorted_ranks.append('organism_name')
    taxonomy = updated_taxmap.loc[:, sorted_ranks].agg('; '.join, axis=1)
    taxonomy.rename('Taxon', inplace=True)
    return taxonomy


def parse_silva_taxonomy(taxonomy_tree: TreeNode,
                         taxonomy_map: pd.DataFrame,
                         taxonomy_ranks: pd.DataFrame,
                         rank_propagation: bool = True,
                         ranks: list = None,
                         include_species_labels: bool = False) -> pd.Series:
    # Traverse the taxonomy hierarchy tree (taxonomy_tree) to obtain the
    # taxids. These will be used to look up the taxonomy and rank information
    # from the taxonomy_ranks file. Finally the taxonomy information is
    # mapped to each Accesioned sequence via the taxonomy_map. An option
    # to include the, potentially untrustworthy, species labels is provided.
    if ranks is None:
        ranks = DEFAULT_RANKS

    taxrank = _prep_taxranks(taxonomy_ranks)
    taxmap = _prep_taxmap(taxonomy_map)
    _validate_taxrank_taxmap_taxtree(taxrank, taxmap, taxonomy_tree)
    silva_tax_id_df = _build_base_silva_taxonomy(taxonomy_tree, taxrank,
                                                 ALLOWED_RANKS,
                                                 rank_propagation)
    updated_taxmap = pd.merge(taxmap, silva_tax_id_df, left_on='taxid',
                              right_index=True)
    taxonomy = _compile_taxonomy_output(updated_taxmap, ranks,
                                        include_species_labels)
    return taxonomy
