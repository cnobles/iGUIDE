#!/bin/bash
# Setup bash script for conda environment iGUIDEseq
conda config --prepend channels 'bushmanlab'
conda config --prepend channels 'r'
conda config --prepend channels 'bioconda'
conda create -n iGUIDEseq --file ./setup/requirements.txt --yes
sleep 1
source activate iGUIDEseq
cd scripts/
git clone https://github.com/cnobles/dualDemultiplexR.git
git clone https://github.com/cnobles/seqTrimR.git
git clone https://github.com/cnobles/seqConsolidateR.git
git clone https://github.com/cnobles/blatCoupleR.git
cd ../
sleep 1
Rscript ./setup/bioconductor_setup.R
