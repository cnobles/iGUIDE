


IGUIDE_ENV_NAME=${1-iguide}

conda env list | grep -Fq ${IGUIDE_ENV_NAME} || {
  bash install.sh
  exit
}

conda env remove --name iguide
rm -rf tools/seq* tools/blat* tools/dual*
bash install.sh
