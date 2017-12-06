#!bin/bash
# Remove created environment for iGUIDEseq
source deactivate
conda env remove -n iGUIDEseq --yes
conda config --remove channels 'bushmanlab'
conda config --remove channels 'bioconda'
conda config --remove channels 'r'

