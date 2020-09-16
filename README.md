# RESCRIPt

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3891931.svg)](https://doi.org/10.5281/zenodo.3891931)
 ![lint-build-test](https://github.com/bokulich-lab/RESCRIPt/workflows/lint-build-test/badge.svg)

 <p align="left"><img src="logo.png" height="256" /></p>

REference Sequence annotation and CuRatIon Pipeline

**Note:** This is a beta release. Usage, and other details are forthcoming. See citation information below.

RESCRIPt is a python 3 package to support a variety of operations for managing and curating reference sequence databases, DNA/RNA sequence data, and taxonomic data.

## Install from source

RESCRIPt will be installable as conda package in the near future. In the meantime, we provide a source installation.

First create a conda environment and install relevant dependencies (you can skip this step and install RESCRIPt in a conda environment containing QIIME 2 2020.2+, which should contain the necessary dependencies/versions):

```
conda create -y -n rescript
conda activate rescript
conda install \
  -c conda-forge -c bioconda -c qiime2 -c defaults \
  qiime2 q2cli q2templates q2-types q2-longitudinal q2-feature-classifier "pandas>=0.25.3" xmltodict
```

Finally install from source:

```
pip install git+https://github.com/bokulich-lab/RESCRIPt.git
```

To view a help menu for using rescript via the QIIME 2 CLI:
```
qiime dev refresh-cache
qiime --help
```

## Using RESCRIPt

These tutorials demonstrate some of the basic functionality of RESCRIPt, via the q2CLI (QIIME 2 command-line interface):
[General Overview and working with SILVA data](https://forum.qiime2.org/t/rescript-sequence-reference-database-management-tutorial/15494)
[Getting sequences and taxonomy with get-ncbi-data](https://forum.qiime2.org/t/using-rescript-to-compile-an-sequence-databases-and-taxonomy-classifiers-from-ncbi-genbank/15947)


## Getting Help
Problem? Suggestion? Technical errors and user support requests can be filed on the [QIIME 2 Forum](https://forum.qiime2.org/).


## Citation

A proper software announcement is forthcoming. In the meantime, if you use RESCRIPt in your research, please cite the Zenodo record:

Nicholas Bokulich, Mike Robeson, Ben Kaehler, & Matthew Dillon. bokulich-lab/RESCRIPt. Zenodo. http://doi.org/10.5281/zenodo.3891931

## License

RESCRIPt is released under a BSD-3-Clause license. See LICENSE for more details.

However, other resources accessible via RESCRIPt are released under different licenses, as detailed below.

**The SILVA database** versions are released under different licenses. Refer to the [current SILVA release license information](https://www.arb-silva.de/silva-license-information/) for more details.

**If using NCBI Genbank data** (e.g., with `get-ncbi-data`): See the [NCBI disclaimer and copyright notice](https://www.ncbi.nlm.nih.gov/home/about/policies/)
