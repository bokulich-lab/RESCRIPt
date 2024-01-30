# RESCRIPt

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3891931.svg)](https://doi.org/10.5281/zenodo.3891931)
 ![lint-build-test](https://github.com/bokulich-lab/RESCRIPt/workflows/lint-build-test/badge.svg)
 [![DOI](https://img.shields.io/badge/DOI-10.1371/journal.pcbi.1009581-B31B1B)](http://dx.doi.org/10.1371/journal.pcbi.1009581)
 <p align="left"><img src="logo.png" height="150" /></p>

REference Sequence annotation and CuRatIon Pipeline

RESCRIPt is a python 3 package to support a variety of operations for managing and curating reference sequence databases, DNA/RNA sequence data, and taxonomic data. See citation information below for a full benchmark and description.

## Install from source

RESCRIPt will be installable as conda package in the near future. In the meantime, we provide two routes for source installation: a minimal RESCRIPt environment, or within an existing QIIME 2 environment:

### Option 1: Minimal RESCRIPt environment:
First create a conda environment and install relevant dependencies:

```
conda create -y -n rescript
conda activate rescript
conda install \
  -c conda-forge -c bioconda -c qiime2 \
  -c https://packages.qiime2.org/qiime2/2024.2/amplicon/passed/ \
  -c defaults qiime2 q2cli q2templates q2-types q2-longitudinal q2-feature-classifier \
  'q2-types-genomics>2023.2' "pandas>=0.25.3" xmltodict ncbi-datasets-pylib
```
Install source:

```
pip install git+https://github.com/bokulich-lab/RESCRIPt.git

qiime rescript --help
```

### Option 2: Install within QIIME 2 environments `2024.2` and later.
This approach assumes you have the latest version of QIIME 2 installed. If not, please refer to the [QIIME 2 installation instructions](https://docs.qiime2.org/2024.2/install/native/). *Make sure the qiime version within the 
http string matches the version of the active qiime environment.*

```
conda activate qiime2-amplicon-2024.2

conda install -c conda-forge -c bioconda -c qiime2 \
    -c https://packages.qiime2.org/qiime2/2024.2/amplicon/passed/ \
    -c defaults xmltodict 'q2-types-genomics>2023.2' ncbi-datasets-pylib
```
Install source:

```
pip install git+https://github.com/bokulich-lab/RESCRIPt.git

qiime rescript --help
```

### Prior versions of RESCRIPt.
For details on how to install prior versions of RESCRIPt, please see the `install-prior-versions.md` document.

### Read help documentation
To view a help menu for using rescript via the QIIME 2 CLI:
```
qiime dev refresh-cache
qiime --help
```

## Using RESCRIPt

These tutorials demonstrate some of the basic functionality of RESCRIPt, via the q2CLI (QIIME 2 command-line interface):
- [General Overview and working with SILVA data](https://forum.qiime2.org/t/rescript-sequence-reference-database-management-tutorial/15494)
- [Getting sequences and taxonomy with get-ncbi-data](https://forum.qiime2.org/t/using-rescript-to-compile-an-sequence-databases-and-taxonomy-classifiers-from-ncbi-genbank/15947)
- [Building a COI database with BOLD sequences](https://forum.qiime2.org/t/building-a-coi-database-from-bold-references/16129)
- [Building a COI database with NCBI sequences](https://forum.qiime2.org/t/building-a-coi-database-from-ncbi-references/16500)
- [Using RESCRIPt's 'extract-seq-segments' to extract reference sequences without PCR primer pairs](https://forum.qiime2.org/t/using-rescripts-extract-seq-segments-to-extract-reference-sequences-without-pcr-primer-pairs/23618)
- [How to train a GTDB SSU classifier using RESCRIPt](https://forum.qiime2.org/t/how-to-train-a-gtdb-ssu-classifier-using-rescript/25725)
- [Constructing an RDP classifier](https://forum.qiime2.org/t/importing-sequence-data-with-lower-case-nucleotide-characters-constructing-an-rdp-classifier-as-an-example/25158)
- [How to train a UNITE classifier using RESCRIPt](https://forum.qiime2.org/t/how-to-train-a-unite-classifier-using-rescript/28285)

Examples of visualizations produced by RESCRIPt actions can be found in this [Visualization Gallery](https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494#heading--seventeenth-header). Other code examples can be found [here](https://github.com/bokulich-lab/db-benchmarks-2020).

## Getting Help
Problem? Suggestion? Technical errors and user support requests can be filed on the [QIIME 2 Forum](https://forum.qiime2.org/).


## Citation

If you use RESCRIPt in your research, please cite the following:

Michael S Robeson II, Devon R O'Rourke, Benjamin D Kaehler, Michal Ziemski, Matthew R Dillon, Jeffrey T Foster, Nicholas A Bokulich. (2021) *RESCRIPt: Reproducible sequence taxonomy reference database management*. PLoS Computational Biology 17 (11): e1009581. doi: [10.1371/journal.pcbi.1009581](http://dx.doi.org/10.1371/journal.pcbi.1009581).


## License

RESCRIPt is released under a BSD-3-Clause license. See LICENSE for more details.

However, other resources accessible via RESCRIPt are released under different licenses, as detailed below.

**If using the SILVA database** (*e.g.*, with `get-silva-data`): Versions are released under different licenses. Refer to the [current SILVA release license information](https://www.arb-silva.de/silva-license-information/) for more details. [How to cite SILVA](https://www.arb-silva.de/contact/) in your work.

**If using NCBI Genbank data** (*e.g.*, with `get-ncbi-data`): See the [NCBI disclaimer and copyright notice](https://www.ncbi.nlm.nih.gov/home/about/policies/) for more details. [How to cite NCBI](https://support.nlm.nih.gov/knowledgebase/article/KA-03391/en-us) in your work.

**If using GTDB data** (*e.g.*, with `get-gtdb-data`): See the [GTDB "about" page](https://gtdb.ecogenomic.org/about) for more details. [How to cite GTDB](https://gtdb.ecogenomic.org/about) in your work.

**If using UNITE data** (*e.g.*, with `get-unite-data`): See the [UNITE citation page](https://unite.ut.ee/cite.php) for more details. [How to cite UNITE]((https://unite.ut.ee/cite.php)) in your work. 
