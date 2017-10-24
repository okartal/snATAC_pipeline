#!/bin/bash
# Stop on error
set -e

## conda environment name

ENV_NAME=bds_scATAC
ENV_NAME_PY3=bds_scATAC_py3

## install wiggler or not

INSTALL_WIGGLER_AND_MCR=0
INSTALL_GEM=0
INSTALL_PEAKSEQ=0

## install packages from official channels (bioconda and r)

conda create -n ${ENV_NAME} --file requirements.txt -y -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer -c anaconda -c conda-forge
conda create -n ${ENV_NAME_PY3} --file requirements_py3.txt -y -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer 

### bash function definition

function add_to_activate {
  if [ ! -f $CONDA_INIT ]; then
    echo > $CONDA_INIT
  fi
  for i in "${CONTENTS[@]}"; do
    if [ $(grep "$i" "$CONDA_INIT" | wc -l ) == 0 ]; then
      echo $i >> "$CONDA_INIT"
    fi
  done
}

## install useful tools for BigDataScript

mkdir -p $HOME/.bds
cp --remove-destination ./utils/bds_scr ./utils/bds_scr_5min ./utils/kill_scr bds.config $HOME/.bds/
cp --remove-destination -rf ./utils/clusterGeneric/ $HOME/.bds/

## install additional packages

source activate $ENV_NAME
conda install -c conda-forge -c bioconda samtools bzip2 #samtools 1.6 to avoid lib error
source deactivate


echo == Installing dependencies has been successfully done. ==
