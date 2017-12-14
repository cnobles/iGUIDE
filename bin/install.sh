#!/bin/bash

set -e

PREFIX=${HOME}/miniconda3

IDSBSEQ_ENV_NAME=${1-idsb}
OUTPUT=${2-/dev/stdout}

install_conda () {
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PREFIX} >> ${OUTPUT}
    export PATH=${PATH}:${PREFIX}/bin
    command -v conda > /dev/null 2>&1 || { echo "Conda still is not on the path, try installing manually"; exit 1; }
    rm Miniconda3-latest-Linux-x86_64.sh
}

command -v conda > /dev/null 2>&1 || { echo "Conda not installed, installing now ..."; install_conda; }

conda config --prepend channels 'bushmanlab'
conda config --prepend channels 'r'
conda config --prepend channels 'bioconda'

# Create enviroment if it does not exist
conda env list | grep -Fxq ${IDSBSEQ_ENV_NAME} || {
    conda create --name=${IDSBSEQ_ENV_NAME} --file=bin/requirements.txt --yes >> ${OUTPUT}
    source activate ${IDSB_ENV_NAME}
    cd tools
    git clone https://github.com/cnobles/dualDemultiplexR.git >> ${OUTPUT}
    git clone https://github.com/cnobles/seqTrimR.git >> ${OUTPUT}
    git clone https://github.com/cnobles/seqConsolidateR.git >> ${OUTPUT}
    git clone https://github.com/cnobles/blatCoupleR.git >> ${OUTPUT}
    cd ../
    echo "Working Dir:"
    pwd
    echo "Contents:"
    ls
    echo "Path:"
    echo ${PATH}
    echo "Home:"
    echo ${HOME}
    Rscript bin/setup_bioconductor.R >> ${OUTPUT}
    echo "iDSBseq successfully installed.";
}

echo "To get started, ensure ${PREFIX}/bin is in your path and run 'source activate ${IDSBSEQ_ENV_NAME}'"
echo "To ensure ${PREFIX}/bin is in your path each time you long in, append the following to your .bashrc or .bash_profile:"
echo "# Append miniconda3/bin to path"
echo "export PATH=~/miniconda3/bin:\$PATH"
