#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2019--, QIIME 2 development team.
#
# Distributed under the terms of the Lesser GPL 3.0 licence.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import pandas as pd
import qiime2

from q2_types.feature_data import DNAFASTAFormat

from ._utilities import run_command


def dereplicate(sequences: DNAFASTAFormat, taxa: pd.DataFrame
                ) -> (pd.Series, pd.DataFrame):
    with tempfile.NamedTemporaryFile() as out_fasta:
        with tempfile.NamedTemporaryFile() as out_uc:
            # dereplicate sequences with vsearch
            cmd = ['vsearch',
                   '--derep_fulllength', str(sequences),
                   '--output', out_fasta.name,
                   '--uc', out_uc.name,
                   '--qmask', 'none',
                   '--xsize']
            run_command(cmd)
            out_uc.seek(0)

            # open dereplicated fasta and uc file
            uc = pd.read_csv(out_uc.name, sep='\t', header=None)
            derep_seqs = qiime2.Artifact.import_data(
                'FeatureData[Sequence]', out_fasta.name).view(pd.Series)

    # transform raw sequences to series for easy parsing
    sequences = qiime2.Artifact.import_data(
        'FeatureData[Sequence]', sequences).view(pd.Series)

    # dereplicate taxonomy and add back sequences with unique taxonomies
    derep_taxa, seqs_out = _dereplicate_taxa(taxa, sequences, derep_seqs, uc)

    return seqs_out, derep_taxa


def _dereplicate_taxa(taxa, raw_seqs, derep_seqs, uc):
    # grab hit and centroid IDs
    uc = uc[uc[0] == 'H'][[8, 9]]
    uc.columns = ['seqID', 'centroidID']

    # map to taxonomy labels
    uc['taxa'] = uc['seqID'].apply(lambda x: taxa.loc[x])
    uc['centroidtaxa'] = uc['centroidID'].apply(lambda x: taxa.loc[x])

    # filter out hits that do not match centroid taxonomy
    rereplicates = uc[uc['taxa'] != uc['centroidtaxa']]

    # grab associated seqs
    rereplicate_seqs = raw_seqs.loc[rereplicates['seqID']]

    # concatenate with dereplicated seqs
    seqs_out = pd.concat([derep_seqs, rereplicate_seqs])

    # generate taxa
    derep_taxa = taxa.reindex(seqs_out.index)
    derep_taxa.index.name = 'Feature ID'

    return derep_taxa, seqs_out
