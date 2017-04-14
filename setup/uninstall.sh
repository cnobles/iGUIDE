#!bin/bash
# Remove created environment for guide-seq-pipe
source deactivate
conda env remove -n guide-seq-pipe --yes
conda config --remove channels 'bushmanlab'
conda config --remove channels 'bioconda'
conda config --remove channels 'r'

