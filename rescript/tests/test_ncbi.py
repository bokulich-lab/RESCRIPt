# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import time

from qiime2 import Metadata
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import rescript
from pandas import DataFrame
from q2_types.feature_data import DNAIterator


class TestNCBI(TestPluginBase):
    package = 'rescript.tests'

    def test_get_ncbi_data(self):
        df = DataFrame(index=['M59083.2', 'AJ234039.1'])
        df.index.name = 'id'
        md = Metadata(df)
        get_ncbi_data = rescript.methods.get_ncbi_data
        acc_seq, acc_tax = get_ncbi_data(accession_ids=md)
        time.sleep(1)  # to avoid a 429
        que_seq, que_tax = get_ncbi_data(query='M59083.2 OR AJ234039.1')

        acc_seq = {s.metadata['id']: str(s) for s in acc_seq.view(DNAIterator)}
        que_seq = {s.metadata['id']: str(s) for s in que_seq.view(DNAIterator)}
        self.assertEqual(acc_seq, que_seq)

        acc_tax = acc_tax.view(DataFrame).to_dict()
        que_tax = que_tax.view(DataFrame).to_dict()
        self.assertEqual(acc_tax, que_tax)

        self.assertEqual(len(acc_seq), 2)
        self.assertEqual(len(acc_tax['Taxon']), 2)
