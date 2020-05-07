# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import pandas as pd
from ._utilities import (_find_lca_in_series, _rank_length, _taxon_to_list,
                         _find_top_score)


def merge_taxa(data: pd.DataFrame, mode: str = 'len', rank_handle: str = '',
               new_rank_handle: str = None) -> pd.DataFrame:
    # Convert taxonomies to list; optionally remove rank handle
    for d in data:
        d['Taxon'] = d['Taxon'].apply(
            lambda x: _taxon_to_list(x, rank_handle=rank_handle))

    # consensus and other dataset-specific data are meaningless after LCA
    # so we will just drop them.
    if mode == 'lca':
        data = [d[['Taxon']] for d in data]

    # We want to copy scores to a uniformly labeled column so that we can merge
    # merge into a unified column, while still preserving original score names
    # (e.g., if merging confidence and consensus scores).
    if mode == 'score':
        for d in data:
            try:
                d['score'] = pd.to_numeric(d.iloc[:, 1])
            # if single-column frame is encountered, raise error
            except IndexError:
                raise IndexError('mode "score" can only operate on dataframes '
                                 'with taxonomy classification scores in the '
                                 'second column. Use "qiime metadata '
                                 'tabulate" to visually confirm the structure '
                                 'of your data before using this command.')
    data = iter(data)
    result = next(data)
    for d in data:
        fill_value = ''
        overwrite = False
        # choose longest taxonomy
        if mode == 'len':
            func = _rank_length
        # or do LCA merging
        elif mode == 'lca':
            func = _find_lca_in_series
        # or choose taxonomy with highest confidence score
        elif mode == 'score':
            func = _find_top_score
            fill_value = 0
            overwrite = True
        # merge data using selected function
        result = result.T.combine(
            d.T, func, fill_value=fill_value, overwrite=overwrite).T

    # Insert new rank handles if selected
    if new_rank_handle is not None:
        result['Taxon'] = result['Taxon'].apply(
            lambda x: ';'.join([''.join(t) for t in zip(
                new_rank_handle.split(';'), x)]))
    else:
        result['Taxon'] = result['Taxon'].apply(lambda x: ';'.join(x))

    return result
