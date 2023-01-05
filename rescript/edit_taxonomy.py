# ----------------------------------------------------------------------------
# Copyright (c) 2019-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import MetadataColumn, Categorical
import warnings
import pandas as pd
import re
from qiime2.plugin import List


def make_search_replace_dict(search_strings, replacement_strings):
    ssl = len(search_strings)
    rsl = len(replacement_strings)

    if ssl == rsl:
        sr_dict = dict(zip(search_strings, replacement_strings))
        return sr_dict
    else:
        raise ValueError('Search list is not the same length as the '
                         'replacement list! The search and replacement '
                         'lists contain %i and %i items, respectively.'
                         % (ssl, rsl))


def edit_taxonomy(taxonomy: pd.Series,
                  replacement_map: MetadataColumn[Categorical] = None,
                  search_strings: List = None,
                  replacement_strings: List = None,
                  use_regex: bool = False,
                  ) -> pd.Series:

    rm_dict = dict()
    if replacement_map:
        rm_dict.update(replacement_map.to_series().to_dict())
        print('Processng strings from replacement map file.')
    if search_strings and replacement_strings:
        cl_dict = make_search_replace_dict(search_strings, replacement_strings)
        rm_dict.update(cl_dict)
        print('Processing strings from command line.')
    if len(rm_dict.keys()) == 0:
        raise ValueError('Either a replacement-map or both search-strings '
                         'and replacement-stings must be supplied.')

    # `use_regex` will basically toggle between using, or not, `re.escape`
    # prior to passing to `re.compile`. This will allow the user to
    # search and replace substrings w/o having to be knowledgable about
    # regular expresions. That is, passing through `re.escape` will behave as
    # a literal string match.
    new_rm_dict = {(re.compile(k) if use_regex else re.compile(
                            re.escape(k))): v for k, v in rm_dict.items()}

    # Note, becuase we are compiling a regex either with or without
    # `re.escape`, we must always set `regex=True` for `pd.replace`.
    updated_tax_series = taxonomy.replace(new_rm_dict, regex=True)

    # User sanity check. Though only appears with `--verbose`
    ranks_in = taxonomy.str.count(';')
    ranks_out = updated_tax_series.str.count(';')
    rank_diffs = ranks_in.compare(ranks_out)
    if not rank_diffs.empty:
        warnings.warn("Warning: the number of taxonomy ranks in the output "
                      "differs from the number of taxonomy ranks of the "
                      "input. Was this intented? The following lines are "
                      "different:\n %s" % rank_diffs, UserWarning)

    return updated_tax_series
