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


def replace_taxonomy(taxonomy: pd.Series,
                     replacement_map: MetadataColumn[Categorical],
                     regex: bool = True,
                     ) -> pd.Series:

    updated_tax_series = taxonomy.replace(replacement_map.to_series(
                                                    ).to_dict(), regex=regex)

    # User sanity check. Though only appears with `--verbose`
    ranks_in = taxonomy.str.count(';')
    ranks_out = updated_tax_series.str.count(';')
    if not ranks_in.compare(ranks_out).empty:
        warnings.warn("Warning: the number of taxonomy ranks in the output "
                      "differs from the number of taxonomy ranks of the "
                      "input. Was this intented?")

    return updated_tax_series
