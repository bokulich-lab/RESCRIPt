# RESCRIPt

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3891932.svg)](https://doi.org/10.5281/zenodo.3891932)
 ![lint-build-test](https://github.com/bokulich-lab/RESCRIPt/workflows/lint-build-test/badge.svg)

REference Sequence annotation and CuRatIon Pipeline

**Note:** This is a beta release. Usage, and other details are forthcoming. See citation information below.

RESCRIPt is a python 3 package to support a variety of operations for managing and curating reference sequence databases, DNA/RNA sequence data, and taxonomic data.

## Install from source

RESCRIPt will be installable as conda package in the near future. In the meantime, we provide a source installation.

First create a conda environment and install relevant dependencies:

```
conda create -y -n rescript
conda activate rescript
conda install \
  -c conda-forge -c bioconda -c qiime2 -c defaults \
  qiime2 q2cli q2templates q2-types q2-longitudinal q2-feature-classifier
```

Finally install from source:

```
pip install git+https://github.com/bokulich-lab/RESCRIPt.git
```

## Citation

A proper software announcement is forthcoming. In the meantime, if you use RESCRIPt in your research, please cite the Zenodo record:

Nicholas Bokulich, Mike Robeson, & Matthew Dillon. (2020, June 12). bokulich-lab/RESCRIPt: 2020.6.0 (Version 2020.6.0). Zenodo. http://doi.org/10.5281/zenodo.3891932

## License

RESCRIPt is released under a BSD-3-Clause license. See LICENSE for more details.

However, other resources accessible via RESCRIPt are released under different licenses, as detailed below.

**The SILVA database** versions are released under different licenses. Refer to the [current SILVA release license information](https://www.arb-silva.de/silva-license-information/) for more details.
