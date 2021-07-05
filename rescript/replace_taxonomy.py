# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.feature_data import TSVTaxonomyFormat
from qiime2.plugin import Metadata
import re


def make_substitutions_dict(taxonomy_replacements):
    # Check that all items in first metadata column or series are unique and
    #     warn the user if not the case. Print repeated values.
    # Then convert tab-delim file into dictionary.
    return substitutions

def make_regex(substitutions):
    strings = substitutions.keys()
    regex = re.compile('|'.join(map(re.escape, strings)))
    return regex

def replace_taxonomy(taxonomy: TSVTaxonomyFormat,
                     taxonomy_replacements: TSVTaxonomyFormat) ->
                     TSVTaxonomyFormat:

    sd = make_substitutions_dict(taxonomy_replacements)
    strings = sd.keys()
    regex = re.compile('|'.join(map(re.escape, strings)))
    mr = make_regex(sd)

    updated_tax_series = tax_series.apply(lambda taxonomy:
                            mr.sub(lambda match: sd[match.group(0)], taxonomy))
    return updated_tax_series
