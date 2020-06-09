# ----------------------------------------------------------------------------
# Copyright (c) 2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import re
import tempfile
import hashlib
import shutil
import gzip
import warnings

import qiime2
import pandas as pd
from skbio.tree import TreeNode
from collections import OrderedDict
from urllib.request import urlretrieve
from urllib.error import HTTPError
from .types._format import RNAFASTAFormat


WHITESPACE_REGEX = re.compile(r'\s+')
ALLOWED_CHARS = set(''.join(['0123456789',
                             'abcdefghijklmnopqrstuvwxyz',
                             'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
                             '_-[]()/.\\']))

# Do not use "'major_clade': 'mc__'" as this appears at multiple levels, even
#      within the same lineage. Tough to disambiguate.
ALLOWED_RANKS = OrderedDict({'domain': 'd__', 'superkingdom': 'sk__',
                             'kingdom': 'k__', 'subkingdom': 'ks__',
                             'superphylum': 'sp__', 'phylum': 'p__',
                             'subphylum': 'ps__', 'infraphylum': 'pi__',
                             'superclass': 'sc__', 'class': 'c__',
                             'subclass': 'cs__', 'infraclass': 'ci__',
                             'superorder': 'so__', 'order': 'o__',
                             'suborder': 'os__', 'superfamily': 'sf__',
                             'family': 'f__', 'subfamily': 'fs__',
                             'genus': 'g__'})

SELECTED_RANKS = OrderedDict({'domain': 'd__', 'phylum': 'p__',
                              'class': 'c__',  'order': 'o__',
                              'family': 'f__', 'genus': 'g__'})


def _keep_allowed_chars(lin_name, allowed_chars=ALLOWED_CHARS,
                        whitespace_pattern=WHITESPACE_REGEX):
    # Only keep "safe" characters and replace whitespace with '_'
    new_lin_name = whitespace_pattern.sub("_", lin_name.strip())
    new_lin_name = ''.join(i for i in new_lin_name if i in allowed_chars)
    return new_lin_name


def _build_base_silva_taxonomy(tree, taxrank, allowed_ranks):
    # Return base silva taxonomy dataframe by traversing taxonomy tree
    # and forward filling the lower ranks with upper-level taxonomy, then
    # pre-pend the rank labels.
    # tree : skbio.TreeNode
    # taxrank : output from _prep_taxranks
    # allowed_ranks : Ordered dict {'domain': 'd__', ...} of allowed ranks to
    # capture and parse
    # output : pandas DataFrame with 'taxid' as the index and the full taxonomy
    # path (i.e. domain-to-genus) as the values.
    tid = {}
    for node in tree.postorder():
        if node.is_root():
            break
        rank_list = []
        # the _validate_taxrank_taxtree func should run
        # before getting here, so all looks up should pass
        rank, taxonomy = taxrank.loc[str(node.name), ['taxrank',
                                                      'taxid_taxonomy']]
        # only pull taxonomy from allowed ranks
        if rank in allowed_ranks:
            rank_list.append((allowed_ranks[rank], taxonomy))
        for ancestor in node.ancestors():
            if ancestor.is_root():
                break
            else:
                arank, ataxonomy = taxrank.loc[str(ancestor.name),
                                               ['taxrank', 'taxid_taxonomy']]
                if arank in allowed_ranks:
                    rank_list.append((allowed_ranks[arank], ataxonomy))
        tid[node.name.strip()] = dict(rank_list)
    silva_tax_id_df = pd.DataFrame.from_dict(tid, orient='index',
                                             columns=allowed_ranks.values())
    silva_tax_id_df.index.name = 'taxid'
    # forward fill missing taxonomy with upper level taxonomy
    silva_tax_id_df.ffill(axis=1, inplace=True)
    return silva_tax_id_df


def _get_clean_organism_name(name):
    # parse organism name (i.e. "species") by only returning the first two
    # words (many species labels are very long or contain sub-species info).
    return _keep_allowed_chars('_'.join(name.strip().split()[:2]))


def _prep_taxmap(taxmap):
    # Return silva taxonomy dataframe map indexed by accession.seqtart.seqstop
    # This is how the coreesponding FASTA file IDs are structured
    taxmap.index = taxmap.index + '.' + taxmap.start + '.' + taxmap.stop
    taxmap.index.name = 'Feature ID'
    taxmap = taxmap[['organism_name', 'taxid']]
    taxmap.loc[:, 'organism_name'] = taxmap['organism_name'].apply(
                                            _get_clean_organism_name)
    return taxmap


def _get_terminal_taxon(name):
    # Return the last taxonomy of a taxonomy path. That is
    # return: 'Aenigmarchaeales'
    # from: 'Archaea;Aenigmarchaeota;Aenigmarchaeia;Aenigmarchaeales;'
    return _keep_allowed_chars(name.rsplit(';')[-2])


def _prep_taxranks(taxrank):
    # return only taxonomy rank information indexed by taxonomy ID
    taxrank = taxrank[[1, 2]].reset_index()
    taxrank.columns = ['taxid_taxonomy', 'taxid', 'taxrank']
    taxrank.set_index('taxid', inplace=True)
    # return only the last taxonomy from a given taxonomy path
    taxrank.loc[:, 'taxid_taxonomy'] = taxrank.loc[:, 'taxid_taxonomy'].apply(
                                             _get_terminal_taxon)
    return taxrank


def _validate_taxrank_taxmap_taxtree(prepped_taxrank, prepped_taxmap, taxtree):
    # Prepass to make sure there is 100 % agreement with the taxids
    # present within the taxranks and taxtree files. Then check that
    # the taxids in the tree are a superset of the taxmap, i.e. not necessarily
    # all taxids are represented by accessioned sequences.
    tree_taxids = {node.name for node in taxtree.postorder()}
    # '1' is the root of the tree of life,
    # which has no tax rank info in any other files.
    tree_taxids.remove('1')  # tree taxids
    ptrs = set(prepped_taxrank.index.unique())  # taxrand taxids
    ptms = set(prepped_taxmap.taxid.unique())  # taxmap taxids

    tree_rank_diffs = tree_taxids.symmetric_difference(ptrs)
    if any(tree_rank_diffs):
        raise ValueError("The taxids of the SILVA Taxonomy Tree file "
                         "and the Taxonomy Ranks file do not match! "
                         "Please check that you are using the same SILVA "
                         "release versions of your files! The taxids that "
                         "are missing in at least one or the other files "
                         "are: ", tree_rank_diffs)
    tree_map_diffs = ptms.difference(tree_taxids)
    if len(tree_map_diffs) > 0:
        raise ValueError("The SILVA Taxonomy Map file conains taxids not "
                         "present within the the SILVA Taxonomy Tree file! "
                         "Please check that you are using the same SILVA "
                         "release versions of your files! The missing "
                         "taxids are: ", tree_map_diffs)


def _compile_taxonomy_output(updated_taxmap, include_species_labels=False,
                             selected_ranks=SELECTED_RANKS):
    # Takes the updated_taxmap dataframe (i.e. the combined
    # output from _build_base_silva_taxonomy and _prep_taxmap) and only
    # returns the desired ranks (as a pandas Series) from selected_ranks, which
    # is an ordered dict. For now, SELECTED_RANKS is used. In the future, we
    # may provide the option to allow users to feed in their own dictionary
    # of selected_ranks. The user can opt to return the species labels.
    sr = list(selected_ranks.values())
    # add rank prefixes (i.e. 'p__')
    updated_taxmap.loc[:, sr] = updated_taxmap.loc[:, sr].apply(
                                                    lambda x: x.name + x)
    if include_species_labels:
        updated_taxmap.loc[:, 'organism_name'] = updated_taxmap.loc[
                           :, 'organism_name'].apply(lambda x: 's__' + x)
        sr.append('organism_name')
    taxonomy = updated_taxmap.loc[:, sr].agg('; '.join, axis=1)
    taxonomy.rename('Taxon', inplace=True)
    return taxonomy


def parse_silva_taxonomy(taxonomy_tree: TreeNode,
                         taxonomy_map: pd.DataFrame,
                         taxonomy_ranks: pd.DataFrame,
                         include_species_labels: bool = False) -> pd.Series:
    # Traverse the taxnomy hierarchy tree (taxonomy_tree) to obtain the
    # taxids. These will be used to look up the taxonomy and rank information
    # from the taxonomy_ranks file. Finally the taxonomy information is
    # mapped to each Accesioned sequence via the taxonomy_map. An option
    # to include the, potentially untrustworthy, species labels is provided.
    taxrank = _prep_taxranks(taxonomy_ranks)
    taxmap = _prep_taxmap(taxonomy_map)
    _validate_taxrank_taxmap_taxtree(taxrank, taxmap, taxonomy_tree)
    silva_tax_id_df = _build_base_silva_taxonomy(taxonomy_tree, taxrank,
                                                 ALLOWED_RANKS)
    updated_taxmap = pd.merge(taxmap, silva_tax_id_df, left_on='taxid',
                              right_index=True)
    taxonomy = _compile_taxonomy_output(updated_taxmap,
                                        include_species_labels,
                                        SELECTED_RANKS)
    return taxonomy


def get_silva_data(ctx,
                   version='138',
                   target='SSURef_NR99',
                   include_species_labels=False,
                   download_sequences=True):
    # download data from SILVA
    print('Downloading raw files may take some time... get some coffee.')
    queries = _assemble_silva_data_urls(version, target, download_sequences)
    results = _retrieve_data_from_silva(queries)
    # parse taxonomy
    parse_taxonomy = ctx.get_action('rescript', 'parse_silva_taxonomy')
    taxonomy, = parse_taxonomy(
        taxonomy_tree=results['taxonomy tree'],
        taxonomy_map=results['taxonomy map'],
        taxonomy_ranks=results['taxonomy ranks'],
        include_species_labels=include_species_labels)
    # if skipping sequences, need to output an empty sequence file.
    if not download_sequences:
        results['sequences'] = qiime2.Artifact.import_data(
            'FeatureData[RNASequence]', RNAFASTAFormat())
    return results['sequences'], taxonomy


def _assemble_silva_data_urls(version, target, download_sequences=True):
    '''Generate SILVA urls, given database version and reference target.'''
    # assemble target urls
    ref_map = {'SSURef_NR99': 'ssu_ref_nr',
               'SSURef_Nr99': 'ssu_ref_nr',
               'SSURef': 'ssu_ref',
               'LSURef': 'lsu_ref'}
    # handle silly inconsistencies in filenames between versions
    if target == 'SSURef_NR99' and int(version) < 138:
        target = 'SSURef_Nr99'
    insert = ref_map[target]
    # now compile URLs
    base_url = 'https://www.arb-silva.de/fileadmin/silva_databases/'\
               'release_{0}/Exports/'.format(version)
    base_url_seqs = base_url + 'SILVA_{0}_{1}_tax_silva.fasta.gz'.format(
        version, target)
    base_url_taxmap = '{0}taxonomy/taxmap_slv_{1}_{2}'.format(
        base_url, insert, version)
    # More SILVA release inconsistencies
    if target == 'SSURef' and version == '132':
        base_url_taxmap += '-corrected.txt.gz'
    else:
        base_url_taxmap += '.txt.gz'
    base_url_tax = '{0}taxonomy/tax_slv_{1}_{2}'.format(
        base_url, insert.split('_')[0], version)
    tree_url = base_url_tax + '.tre'
    tax_url = base_url_tax + '.txt'
    if version == '138':
        tree_url += '.gz'
        tax_url += '.gz'

    # download and validate silva files
    queries = [('sequences', base_url_seqs, 'FeatureData[RNASequence]'),
               ('taxonomy map', base_url_taxmap, 'FeatureData[SILVATaxidMap]'),
               ('taxonomy tree', tree_url, 'Phylogeny[Rooted]'),
               ('taxonomy ranks', tax_url, 'FeatureData[SILVATaxonomy]')]

    # optionally skip downloading sequences
    if not download_sequences:
        queries = queries[1:]

    return queries


def _retrieve_data_from_silva(queries):
    '''
    Download data from SILVA, given a list of queries.

    queries: list of tuples of (str, str, str)
        (name, urlpath, QIIME 2 artifact type)
    '''
    results = dict()
    with tempfile.TemporaryDirectory() as tmpdirname:
        for name, query, dtype in queries:
            print('retrieving {0} from: {1}'.format(name, query))
            # grab url
            destination = os.path.join(tmpdirname, os.path.basename(query))
            urlretrieve(query, destination)
            file_md5 = _get_md5(destination)
            # grab expected md5
            # NOTE: SILVA is missing md5s for some files, so we will just skip
            try:
                md5_destination = os.path.join(tmpdirname, 'md5')
                urlretrieve(query + '.md5', md5_destination)
                exp_md5 = _read_silva_md5(md5_destination)
                # validate md5 checksum
                _validate_md5(exp_md5, file_md5, query)
            # if we get an HTTPError just move on, the md5 file does not exist
            except HTTPError:
                msg = ("No md5 file was detected in the SILVA archive for the "
                       "following file. No action is required, but be aware "
                       "that md5-hash validation was not performed for this "
                       "file: " + query)
                warnings.warn(msg, UserWarning)
            # gunzip on demand (SILVA releases are inconsistently gzipped)
            try:
                unzipped_destination = os.path.splitext(destination)[0]
                _gzip_decompress(destination, unzipped_destination)
                destination = unzipped_destination
            except OSError:
                pass
            # import as artifacts
            results[name] = qiime2.Artifact.import_data(dtype, destination)
    return results


def _validate_md5(exp_md5, file_md5, filename):
    if not exp_md5 == file_md5:
        raise ValueError('md5 sums do not match. Manually verify md5 sums '
                         'before proceeding.\nTarget file: {0}\nExpected md5: '
                         '{1}\nObserved md5: {2}\n'.format(
                            filename, exp_md5, file_md5))


# This function is specific for reading the SILVA md5 record files, which are
# a single line txt file in the format "md5  filename"
def _read_silva_md5(file):
    with open(file, 'r') as _md5:
        _md5 = _md5.read().split(' ')[0]
    return _md5


def _get_md5(file, chunksize=8192):
    md5_hash = hashlib.md5()
    with open(file, "rb") as f:
        for chunk in iter(lambda: f.read(chunksize), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def _gzip_decompress(input_fp, output_fp):
    with gzip.open(input_fp, 'rt') as temp_in:
        with open(output_fp, 'w') as temp_out:
            shutil.copyfileobj(temp_in, temp_out)
