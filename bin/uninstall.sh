#!bin/bash
# Remove created environment for iDSBseq
source deactivate
conda env remove -n iguide --yes
conda config --remove channels 'bushmanlab'
conda config --remove channels 'bioconda'
conda config --remove channels 'r'

