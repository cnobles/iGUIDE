
IGUIDE_ENV_NAME=${1-iguide}

rm -rf analysis/simulation/

conda env list | grep -Fq ${IGUIDE_ENV_NAME} && {
  source deactivate
  conda env remove -y --name iguide
}

bash install.sh
source activate iguide
iguide setup configs/simulation.config.yml
iguide run configs/simulation.config.yml
