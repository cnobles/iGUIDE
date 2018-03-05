#!bin/bash
# Remove created environment for iGUIDE
source deactivate
conda env remove -n iguide --yes
conda config --remove channels 'bushmanlab'
conda config --remove channels 'bioconda'
conda config --remove channels 'r'
conda config --remove channels 'conda-forge'
