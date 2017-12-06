#!/bin/bash

# Setup bash script for conda environment iDSBseq
conda config --prepend channels 'bushmanlab'
conda config --prepend channels 'r'
conda config --prepend channels 'bioconda'
conda create -n idsb --file ./setup/requirements.txt --yes
sleep 1
source activate idsb

cd tools
git clone https://github.com/cnobles/dualDemultiplexR.git
git clone https://github.com/cnobles/seqTrimR.git
git clone https://github.com/cnobles/seqConsolidateR.git
git clone https://github.com/cnobles/blatCoupleR.git
cd ../

sleep 1
Rscript ./bin/bioconductor_setup.R
