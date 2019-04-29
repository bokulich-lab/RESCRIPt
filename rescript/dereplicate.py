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

from ._utilities import run_command, _find_lca, _majority


def dereplicate(sequences: DNAFASTAFormat, taxa: pd.DataFrame,
                mode: str = 'uniq') -> (pd.Series, pd.DataFrame):
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
    derep_taxa, seqs_out = _dereplicate_taxa(
        taxa, sequences, derep_seqs, uc, mode=mode)

    return seqs_out, derep_taxa


def _dereplicate_taxa(taxa, raw_seqs, derep_seqs, uc, mode):
    # grab hit and centroid IDs
    if mode == 'uniq':
        uc_types = ['H']
    else:
        # we want to keep the centroids if we perform LCA or majority
        uc_types = ['H', 'S']
        # centroid entries have no centroid ID; so list their own seq ID
        uc.loc[uc[9] == '*', 9] = uc[8]

    uc = uc[uc[0].isin(uc_types)][[8, 9]]
    uc.columns = ['seqID', 'centroidID']

    # map to taxonomy labels
    uc['Taxon'] = uc['seqID'].apply(lambda x: taxa.loc[x])
    uc['centroidtaxa'] = uc['centroidID'].apply(lambda x: taxa.loc[x])

    if mode == 'uniq':
        # filter out hits that do not match centroid taxonomy
        rereplicates = uc[uc['Taxon'] != uc['centroidtaxa']]
        # drop duplicates that share centroid ID and taxa assignment
        rereplicates = rereplicates.drop_duplicates(['centroidID', 'Taxon'])
        # grab associated seqs
        rereplicate_seqs = raw_seqs.loc[rereplicates['seqID']]
        # concatenate with dereplicated seqs
        seqs_out = pd.concat([derep_seqs, rereplicate_seqs])
        # generate list of dereplicated taxa
        derep_taxa = taxa.reindex(seqs_out.index)

    else:
        # group seqs that share centroids (this includes the centroid)
        derep_taxa = uc.groupby('centroidID')['Taxon'].apply(lambda x: list(x))
        # find LCA within each cluster
        if mode == 'lca':
            derep_taxa = derep_taxa.apply(lambda x: ';'.join(
                _find_lca([y.split(';') for y in x]))).to_frame()
            # find majority taxon within each cluster
        elif mode == 'majority':
            derep_taxa = derep_taxa.apply(lambda x: _majority(x)).to_frame()
        # LCA and majority do nothing with the seqs
        seqs_out = derep_seqs

    # gotta please the type validator gods
    derep_taxa.index.name = 'Feature ID'

    return derep_taxa, seqs_out
