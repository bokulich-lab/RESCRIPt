# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pandas as pd
from ._utilities import (_rank_length, _taxon_to_list, _find_top_score,
                         _rank_handles, _find_lca, _find_super_lca,
                         _find_lca_majority)


MODE_ERROR_SCORE = (
    'mode "score" can only operate on dataframes with taxonomy classification '
    'scores in the second column. Use "qiime metadata tabulate" to visually '
    'confirm the structure of your data before using this command.')


def merge_taxa(data: pd.DataFrame,
               mode: str = 'len',
               rank_handle_regex: str = '^[dkpcofgs]__',
               new_rank_handle: str = None,
               unclassified_label: str = 'Unassigned') -> pd.DataFrame:
    # Convert taxonomies to list; optionally remove rank handle
    for d in data:
        # Convert taxonomies to list; optionally remove rank handle
        d['Taxon'] = d['Taxon'].apply(
            lambda x: _taxon_to_list(x, rank_handle=rank_handle_regex))

        # Capitalize column names for sorting consistency
        d.columns = [x.capitalize() for x in d.columns]

    # consensus and other dataset-specific data are meaningless after LCA
    # or majority so we will just drop them and apply functions across rows.
    if mode in ['lca', 'super', 'majority']:
        data = [d[['Taxon']] for d in data]
        data = pd.concat(data, axis=1, sort=True).fillna('')
        if mode == 'lca':
            func = _find_lca
        elif mode == 'super':
            func = _find_super_lca
        elif mode == 'majority':
            func = _find_lca_majority
        result = data.apply(lambda x: func([t for t in x if t != '']), axis=1)
        result = result.to_frame(name='Taxon')

    # len and score modes are computed pairwise to preserve other taxon info
    else:
        if mode == 'len':
            func, fill_value, overwrite = _rank_length, '', False
        # We want to copy scores to a uniformly labeled column so that we can
        # merge into a unified column, while still preserving original score
        # names (e.g., if merging confidence and consensus scores).
        if mode == 'score':
            func, fill_value, overwrite = _find_top_score, 0, True
            for d in data:
                try:
                    d['Score'] = pd.to_numeric(d.iloc[:, 1])
                # if single-column frame is encountered, raise error
                except IndexError:
                    raise IndexError(MODE_ERROR_SCORE)
        data = iter(data)
        result = next(data)
        for d in data:
            result = result.T.combine(
                d.T, func, fill_value=fill_value, overwrite=overwrite).T
            reordered_cols = ['Taxon'] + [x for x in result.columns
                                          if x != 'Taxon']
            result = result[reordered_cols]

    # Insert new rank handles if selected
    if new_rank_handle is not None:
        new_rank_handle = _rank_handles[new_rank_handle]
        result['Taxon'] = result['Taxon'].apply(
            lambda x: ';'.join([''.join(t) for t in zip(new_rank_handle, x)]))
    else:
        result['Taxon'] = result['Taxon'].apply(lambda x: ';'.join(x))

    # fill unassigned taxa, if any
    result['Taxon'].replace('', unclassified_label, inplace=True)

    # gotta please the type validator gods
    result.index.name = 'Feature ID'

    return result
