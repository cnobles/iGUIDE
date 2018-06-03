#!bin/bash
# Remove created environment for iGUIDE
IGUIDE_ENV_NAME=${1-iguide} 

source deactivate
conda env remove -n ${IGUIDE_ENV_NAME} --yes
conda config --remove channels 'bushmanlab'
conda config --remove channels 'bioconda'
conda config --remove channels 'r'
conda config --remove channels 'conda-forge'
