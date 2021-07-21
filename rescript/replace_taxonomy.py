# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Metadata
import pandas as pd
import re


def make_substitutions_dict(taxonomy_replacement_map):
    if taxonomy_replacement_map.column_count < 1:
        raise ValueError("Metadata file should only contain one column of "
                         "taxonomy replacement strings, none found.")
    elif taxonomy_replacement_map.column_count > 1:
        raise ValueError("Metadata file should only contain one column of "
                         "taxonomy replacement strings, found %d columns."
                         % taxonomy_replacement_map.column_count)
    else:
        replacements = taxonomy_replacement_map.to_dataframe().iloc[:, 0]
        replacement_dict = replacements.to_dict()
        return replacement_dict


def make_regex(substitutions):
    strings = substitutions.keys()
    regex = re.compile('|'.join(map(re.escape, strings)))
    return regex


def replace_taxonomy(taxonomy: pd.Series,
                     taxonomy_replacement_map: Metadata
                     ) -> pd.Series:

    rd = make_substitutions_dict(taxonomy_replacement_map)
    mr = make_regex(rd)

    updated_tax_series = taxonomy.apply(lambda taxonomy:
                                        mr.sub(lambda match:
                                               rd[match.group(0)], taxonomy))

    return updated_tax_series
