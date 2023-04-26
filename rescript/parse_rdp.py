# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
from collections import OrderedDict, defaultdict
from q2_types.feature_data import DNAIterator
from .parse_silva_taxonomy import _keep_allowed_chars


ALLOWED_RDP_RANKS = OrderedDict({'rootrank': 'root__', 'domain': 'd__',
                                 'phylum': 'p__', 'class': 'c__',
                                 'subclass': 'cs__', 'order': 'o__',
                                 'suborder': 'os__', 'family': 'f__',
                                 'genus': 'g__'})

DEFAULT_RDP_RANKS = ['domain', 'phylum', 'class', 'subclass', 'order',
                     'suborder', 'family', 'genus']


def parse_rdp_taxonomy(rdp_reference_sequences: DNAIterator,
                       rank_propagation: bool = True, ranks: list = None,
                       include_species_labels: bool = False) -> pd.Series:

    sp_label_dict = defaultdict()
    taxd = defaultdict()
    for seq in rdp_reference_sequences:
        # parse taxonomy
        seq_id = seq.metadata['id']
        sp_label, tax_data = seq.metadata['description'].split(
                                'Lineage=')
        new_sp_label = _keep_allowed_chars(
                                        "_".join(sp_label.strip().split()[:2]))
        sp_label_dict[seq_id] = new_sp_label

        known_lineage = tax_data.replace('"', '').rsplit(';')
        if tax_data.endswith(';'):
            known_lineage = known_lineage[:-1]

        taxd[seq_id] = dict(
                            [(x[1], x[0]) for x in zip(
                                known_lineage[0::2], known_lineage[1::2])])

    # make dataframe with ranks as column labels.
    df = pd.DataFrame.from_dict(taxd, orient='index',
                                columns=DEFAULT_RDP_RANKS)

    # propagate ranks (or not), fillna, and append prefixes
    if rank_propagation:
        df.ffill(axis=1, inplace=True)

    # rename columns with prefixes
    df.rename(columns=ALLOWED_RDP_RANKS, inplace=True)

    # In case we decide NOT to ffill. We'll need an empty string, so
    # that we can have empty prefixes (e.g. 'p__') appear in our
    # taxonomy file.Otherwise these will be empty cells in the
    # dataframe.
    df.fillna('', inplace=True)
    # use default ranks if no rank list supplied
    if ranks is None:
        ranks = DEFAULT_RDP_RANKS

    # Sort user ranks relative to allowed_ranks
    # this ensures that the taxonomic ranks supplied by the user
    # are in order
    sorted_ranks = [p for r, p in ALLOWED_RDP_RANKS.items()
                    if r in ranks]

    # If we wish to include species labels we need to add to the sorted_ranks
    # and the dataframe.
    if include_species_labels:
        sorted_ranks.append('s__')
        df.loc[:, 's__'] = df.index.map(sp_label_dict)

    # add prefixes to selected ranks
    df.loc[:, sorted_ranks] = \
        df.loc[:, sorted_ranks].apply(lambda x: x.name + x)

    # subset dataframe based on user-selected ranks
    rdp_taxonomy = df.loc[:, sorted_ranks].agg('; '.join, axis=1)
    rdp_taxonomy.index.name = 'Feature ID'
    rdp_taxonomy.rename('Taxon', inplace=True)

    return rdp_taxonomy
