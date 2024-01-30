# How to install prior versions of RESCRIPt.

## Install within the QIIME 2 `2023.9` release: 

For this version of QIIME 2, RESCRIPt is included within the [qiime2-shotgun-2023.9 distribution](https://docs.qiime2.org/2023.9/install/native/#qiime-2-shotgun-distribution). You'll be able to run your RESCRIPt commands within this environment, then switch back to [qiime2-amplicon-2023.9](https://docs.qiime2.org/2023.9/install/native/#qiime-2-amplicon-distribution) environment if needed. 

It is possible to install RESCRIPt within `qiime2-amplicon-2023.9`. YOu can try via the following commands:

```
conda activate qiime2-amplicon-2023.9

conda install -c conda-forge -c bioconda -c qiime2 \
    -c https://packages.qiime2.org/qiime2/2023.9/shotgun/passed/  \
    -c defaults   xmltodict 'q2-types-genomics>2023.5' ncbi-datasets-pylib
```

Install source:

```
pip install git+https://github.com/bokulich-lab/RESCRIPt.git

qiime rescript --help
```


## Install within the QIIME 2 releases `2023.2` - `2023.7`.
The typical apporoach will be something like the command below. *Make sure the qiime version within the http string matches the version of the active qiime environment.*

```
conda activate qiime2-2023.7

conda install -c conda-forge -c bioconda -c qiime2 \
    -c https://packages.qiime2.org/qiime2/2023.7/passed/ \
    -c defaults xmltodict 'q2-types-genomics>2023.2' ncbi-datasets-pylib
```
Install source:

```
pip install git+https://github.com/bokulich-lab/RESCRIPt.git

qiime rescript --help
```

## Install within prior releases of QIIME 2 `2022.8` and earlier.
Download any of the prior [Releases](https://github.com/bokulich-lab/RESCRIPt/releases) and consult that version's README file for appropriate Installation instructions.