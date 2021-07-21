# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import MetadataColumn, Categorical
import pandas as pd
import re


def make_substitutions_dict(taxonomy_replacement_map):
    replacement_dict = taxonomy_replacement_map.to_series().to_dict()
    return replacement_dict


def make_regex(substitutions):
    strings = substitutions.keys()
    regex = re.compile('|'.join(map(re.escape, strings)))
    return regex


def replace_taxonomy(taxonomy: pd.Series,
                     taxonomy_replacement_map: MetadataColumn[Categorical]
                     ) -> pd.Series:

    rd = make_substitutions_dict(taxonomy_replacement_map)
    mr = make_regex(rd)

    updated_tax_series = taxonomy.apply(lambda taxonomy:
                                        mr.sub(lambda match:
                                               rd[match.group(0)], taxonomy))

    return updated_tax_series
