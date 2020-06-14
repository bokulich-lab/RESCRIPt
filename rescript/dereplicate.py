# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import pandas as pd
import qiime2

from q2_types.feature_data import DNAFASTAFormat

from ._utilities import run_command, _find_lca, _majority, _rank_handles


def dereplicate(sequences: DNAFASTAFormat,
                taxa: pd.DataFrame,
                mode: str = 'uniq',
                perc_identity: float = 1.0,
                threads: int = 1,
                rank_handles: str = 'silva',
                derep_prefix: bool = False) -> (pd.Series, pd.DataFrame):
    with tempfile.NamedTemporaryFile() as out_fasta:
        with tempfile.NamedTemporaryFile() as out_uc:
            # dereplicate sequences with vsearch
            # note that multithreading is not supported, but we will leave the
            # threads argument in the command in case support is added one day.
            # (multithreading _is_ supported in _vsearch_cluster_size below)
            _vsearch_derep(str(sequences), out_fasta.name, out_uc.name,
                           str(threads), derep_prefix)
            out_uc.seek(0)

            uc = _parse_uc(out_uc.name)

            # optionally cluster seqs into OTUs
            if perc_identity < 1.0:
                clustered_seqs = DNAFASTAFormat()
                _vsearch_cluster_size(str(out_fasta.name), str(perc_identity),
                                      str(clustered_seqs), out_uc.name,
                                      str(threads))
                out_uc.seek(0)

                # re-map derep centroids to cluster centroids
                uc_clust = _parse_uc(out_uc.name).set_index('seqID')
                uc['centroidID'] = uc['centroidID'].apply(
                    lambda x: uc_clust.loc[x, 'centroidID'])
            else:
                clustered_seqs = out_fasta.name

            # open dereplicated fasta and uc file
            derep_seqs = qiime2.Artifact.import_data(
                'FeatureData[Sequence]', str(clustered_seqs)).view(pd.Series)

            # transform raw sequences to series for easy parsing
            sequences = qiime2.Artifact.import_data(
                'FeatureData[Sequence]', sequences).view(pd.Series)

            derep_taxa, seqs_out = _dereplicate_taxa(
                taxa, sequences, derep_seqs, uc, mode=mode)

            if rank_handles != 'disable':
                rank_handles = _rank_handles[rank_handles]
                derep_taxa.loc[:, 'Taxon'] = derep_taxa['Taxon'].apply(
                    _backfill_taxonomy, args=([rank_handles]))

    return seqs_out, derep_taxa


def _backfill_taxonomy(taxon, rank_handles):
    formatted_taxon = taxon.split(';')
    tax_len = len(formatted_taxon)
    if tax_len >= len(rank_handles):
        return taxon
    else:
        return ';'.join(formatted_taxon + rank_handles[tax_len:])


def _vsearch_derep(sequences_fp, out_fasta_fp, out_uc_fp, threads,
                   derep_prefix):
    cmd = ['vsearch',
           '--derep_fulllength', sequences_fp,
           '--output', out_fasta_fp,
           '--uc', out_uc_fp,
           '--qmask', 'none',
           '--xsize',
           '--threads', threads]
    if derep_prefix:
        cmd[1] = '--derep_prefix'
    run_command(cmd)


def _vsearch_cluster_size(sequences_fp, perc_identity, out_fasta_fp, out_uc_fp,
                          threads):
    cmd = ['vsearch',
           '--cluster_size', sequences_fp,
           '--id', perc_identity,
           '--centroids', out_fasta_fp,
           '--uc', out_uc_fp,
           '--qmask', 'none',
           '--xsize',
           '--threads', threads]
    run_command(cmd)


def _parse_uc(uc_fp):
    uc = pd.read_csv(uc_fp, sep='\t', header=None, dtype=object)
    # grab hit and centroid IDs
    uc = uc[uc[0].isin(['H', 'S'])][[8, 9]]
    # centroid entries have no centroid ID; so list their own seq ID
    uc.loc[uc[9] == '*', 9] = uc[8]
    uc.columns = ['seqID', 'centroidID']
    return uc


def _dereplicate_taxa(taxa, raw_seqs, derep_seqs, uc, mode):
    # we only want to grab hits for uniq mode
    if mode == 'uniq':
        uc = uc[uc['seqID'] != uc['centroidID']]
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
