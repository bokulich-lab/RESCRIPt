# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import MetadataColumn, Categorical
import warnings
import pandas as pd
import re


def replace_taxonomy(taxonomy: pd.Series,
                     replacement_map: MetadataColumn[Categorical],
                     use_regex: bool = False,
                     ) -> pd.Series:

    rm_dict = replacement_map.to_series().to_dict()
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
    if not ranks_in.compare(ranks_out).empty:
        warnings.warn("Warning: the number of taxonomy ranks in the output "
                      "differs from the number of taxonomy ranks of the "
                      "input. Was this intented?")

    return updated_tax_series
