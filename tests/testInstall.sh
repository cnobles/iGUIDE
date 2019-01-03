
IGUIDE_ENV_NAME=${1-iguide}

rm -rf analysis/simulation/
rm -rf tools/seq* tools/blat* tools/dual*

conda env list | grep -Fq ${IGUIDE_ENV_NAME} && {
  source deactivate
  conda env remove -y --name iguide
}

bash install.sh
source activate iguide
iguide setup
iguide run
