# How to install other versions of RESCRIPt.

## Minimal RESCRIPt environment (`2023.9` and later)

To make use of the latest "minimal" RESCRIPt release, some components of QIIME 2 (releases `2023.9` and later) are required:

First create a conda environment and install relevant dependencies using either `conda` or `mamba`. 

*Note: update `{ENV_VERSION}` in the commands below to match the QIIME 2 release.*

```
conda create -y -n rescript
conda activate rescript

conda install \
  -c https://packages.qiime2.org/qiime2/{ENV_VERSION}/shotgun/passed/ \
  -c https://packages.qiime2.org/qiime2/{ENV_VERSION}/amplicon/passed/ \
  -c conda-forge -c bioconda -c qiime2 -c defaults \
  qiime2 q2cli q2templates q2-types q2-types-genomics q2-longitudinal q2-feature-classifier \
  "pandas>=0.25.3" xmltodict ncbi-datasets-pylib rescript

qiime rescript --help
```

## Install within the QIIME 2 `2023.9` release: 

For this version of QIIME 2, RESCRIPt is included within the [qiime2-shotgun-2023.9 distribution](https://docs.qiime2.org/2023.9/install/native/#qiime-2-shotgun-distribution). You'll be able to run your RESCRIPt commands within this environment, then switch back to [qiime2-amplicon-2023.9](https://docs.qiime2.org/2023.9/install/native/#qiime-2-amplicon-distribution) environment if needed. It is possible to install RESCRIPt within `qiime2-amplicon-2023.9`. You can try via the following commands:

```
conda activate qiime2-amplicon-2023.9

conda install -c conda-forge -c bioconda -c qiime2 \
    -c https://packages.qiime2.org/qiime2/2023.9/shotgun/passed/  \
    -c defaults   xmltodict 'q2-types-genomics>2023.5' ncbi-datasets-pylib rescript

qiime rescript --help
```

## Install within the QIIME 2 releases `2023.2` - `2023.7`.
Current releases of RESCRIPt should generally be compatable with these QIIME 2 versions. If not, please download the corresponding [RESCRIPt release](https://github.com/bokulich-lab/RESCRIPt/releases) and consult that version's README file for appropriate Installation instructions.

*Note: update `{ENV_VERSION}` in the commands below to match the QIIME 2 release.*

```
conda activate qiime2-{ENV_VERSION}

conda install -c conda-forge -c bioconda -c qiime2 \
    -c https://packages.qiime2.org/qiime2/{ENV_VERSION}/passed/ \
    -c defaults xmltodict 'q2-types-genomics>2023.2' ncbi-datasets-pylib
```
Install source:

```
pip install git+https://github.com/bokulich-lab/RESCRIPt.git

qiime rescript --help
```

## Install within prior releases of QIIME 2 `2022.8` and earlier.
Download any of the prior [Releases](https://github.com/bokulich-lab/RESCRIPt/releases) and consult that version's README file for appropriate Installation instructions.