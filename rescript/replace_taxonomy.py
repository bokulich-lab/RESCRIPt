# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import MetadataColumn, Categorical
import pandas as pd


def lineages_with_invalid_number_of_ranks(updated_tax_series,
                                          num_expected_ranks):
    invalid_lineages = [label for label, lineage in updated_tax_series.items()
                        if len(lineage.split(';')) != num_expected_ranks]
    return invalid_lineages


def replace_taxonomy(taxonomy: pd.Series,
                     replacement_map: MetadataColumn[Categorical],
                     regex: bool = True,
                     num_expected_ranks: int = 7,
                     ) -> pd.Series:

    updated_tax_series = taxonomy.replace(replacement_map.to_series(
                                                    ).to_dict(), regex=regex)

    irl = lineages_with_invalid_number_of_ranks(updated_tax_series,
                                                num_expected_ranks)
    if irl == []:
        return updated_tax_series
    else:
        raise ValueError("The following lineages do not contain "
                         "the expected number of ranks:", irl)
