# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import zipfile
import py7zr
from urllib.request import urlretrieve
import pandas as pd
from q2_types.feature_data import (TSVTaxonomyFormat, DNAFASTAFormat,
                                   DNAIterator)

RRNA_GENE_LIST = ['SSU', 'LSU', 'ITS', 'longread', 'all']


def _assemble_rrna_url(rrna_gene,
                       version='1.9.4',
                       ):

    base_url = ('https://sisu.ut.ee/wp-content/uploads/sites/643/'
                'General_EUK_{gene}_v{ver}.zip')

    base_format_dict = {'gene': rrna_gene,
                        'ver': version}

    rrna_url = base_url.format(**base_format_dict)
    return rrna_url


def _get_eukaryome_data_path(tmpdirname, url, basename):
    destination = os.path.join(tmpdirname, basename)
    try:
        print(' Retrieving from {0}'.format(url))
        urlretrieve(url, destination)
    except Exception as e:
        raise ValueError(
                    'Unable to retrieve the following file '
                    'from Eukaryome:\n{url}\n: {e}'.format(
                                      url=url, e=e))
    return destination


def _extract_zip_euk_seq_file(tmpdirname, uncompressed_base_fn, fasta_zfp):
    uncompressed_fn = uncompressed_base_fn + '.fasta'
    out_path = os.path.join(tmpdirname, uncompressed_fn)
    with zipfile.ZipFile(fasta_zfp, 'r') as zipf:
        zipf.extract(uncompressed_fn, tmpdirname)
        seqs = DNAFASTAFormat(out_path, mode="r").view(DNAIterator)
    return seqs


def _extract_7zip_euk_seq_file(tmpdirname, uncompressed_base_fn, fasta_zfp):
    uncompressed_7z_fn = uncompressed_base_fn + '.7z'
    uncompressed_fn = uncompressed_base_fn + '.fasta'
    compressed_out_path = os.path.join(tmpdirname,  uncompressed_7z_fn)
    uncompressed_out_path = os.path.join(tmpdirname, uncompressed_fn)

    with zipfile.ZipFile(fasta_zfp, 'r') as zipf:
        zipf.extract(uncompressed_7z_fn, tmpdirname)
        with py7zr.SevenZipFile(compressed_out_path, mode='r') as zip7f:
            zip7f.extract(path=tmpdirname, targets=[uncompressed_fn])
            seqs = DNAFASTAFormat(uncompressed_out_path,
                                  mode="r").view(DNAIterator)
    return seqs


def _make_fasta_str(seqid, seq_str):
    fasta_str = '>' + seqid + '\n' + seq_str + '\n'
    return fasta_str


def _make_taxonomy_df(tax_dict):
    parsed_taxonomy_df = pd.DataFrame.from_dict(tax_dict, orient='index')
    parsed_taxonomy_df.index.name = 'Feature ID'
    parsed_taxonomy_df.rename(columns={parsed_taxonomy_df.columns[0]:
                                       'Taxon'}, inplace=True)
    return parsed_taxonomy_df


def _process_eukaryome_data(euk_seq_data):
    proc_seqs = DNAFASTAFormat()
    proc_tax = TSVTaxonomyFormat()
    tax_dict = {}
    with proc_seqs.open() as out_fasta, proc_tax.open() as out_tax:
        for seq in euk_seq_data:
            seqid, tax = seq.metadata['id'].split(';', 1)
            seq_str = str(seq)
            fasta_str = _make_fasta_str(seqid, seq_str)
            out_fasta.write(fasta_str)

            tax_dict[seqid] = tax
        print('    Sequences processed.')
        print('    Processing taxonomy...')
        parsed_taxonomy_df = _make_taxonomy_df(tax_dict)
        parsed_taxonomy_df.to_csv(out_tax, sep='\t')
        print('    Taxonomy processed.')
        return proc_seqs, proc_tax


def _retrieve_data_from_eukaryome(rrna_url, gene):
    with tempfile.TemporaryDirectory() as tmpdirname:
        basename = os.path.basename(rrna_url)
        uncompressed_base_fn = basename.rsplit('.', 1)[0]

        print('Fetching \'{0}\' data ...'.format(gene))
        fasta_zfp = _get_eukaryome_data_path(tmpdirname, rrna_url,
                                             basename)

        print('  Processing \'{0}\' data ...'.format(gene))
        if gene == 'ITS':
            seqs = _extract_7zip_euk_seq_file(tmpdirname,
                                              uncompressed_base_fn,
                                              fasta_zfp)
        else:
            seqs = _extract_zip_euk_seq_file(tmpdirname,
                                             uncompressed_base_fn,
                                             fasta_zfp)

        print('    Processing sequences ...')
        proc_seqs, proc_tax = _process_eukaryome_data(seqs)
        return proc_seqs, proc_tax


def get_eukaryome_data(
    rrna_gene: list,
    version: str = '1.9.4',
        ) -> (DNAFASTAFormat, TSVTaxonomyFormat):

    seq_res = {}
    tax_res = {}

    if 'all' in rrna_gene:
        RRNA_GENE_LIST.pop()  # remove 'all' from choices.
        rrna_gene = RRNA_GENE_LIST

    for gene in rrna_gene:
        rrna_url = _assemble_rrna_url(
                                  rrna_gene=gene,
                                  version=version,
                                  )

        seqs, tax = _retrieve_data_from_eukaryome(rrna_url, gene)
        seq_res[gene + '_seqs'] = seqs
        tax_res[gene + '_taxa'] = tax

    print('\n Saving files...\n')
    return seq_res, tax_res
