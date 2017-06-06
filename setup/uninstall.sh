#!bin/bash
# Remove created environment for guidesnaker
source deactivate
conda env remove -n guidesnaker --yes
conda config --remove channels 'bushmanlab'
conda config --remove channels 'bioconda'
conda config --remove channels 'r'

