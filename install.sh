#!/usr/bin/env bash

__conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

source install_init.sh

function __err_report() {
    local error_code
    error_code=${?}
    error "Error in ${__file} in function ${1} on line ${2}"
    exit ${error_code}
}
trap '__err_report "${FUNCNAME:-.}" ${LINENO}' ERR  

# help mode
if [[ "${arg_h:?}" = "1" ]]; then
    # Help exists with code 1
    help "Help using ${0}:"
fi

# verbose mode
if [[ "${arg_v:?}" = "1" ]]; then
    LOG_LEVEL="7"
fi

# debug mode
if [[ "${arg_d:?}" = "1" ]]; then
    set -o xtrace
    LOG_LEVEL="7"
fi

function debug_capture () {
    debug "$(echo -e "$(${@})")"
}

function installation_error () {
    error "${1} failed!"
    if [[ "${arg_v:?}" != 1 && "${arg_d:?}" != 1 ]]; then
        error "Try re-running with -v or -d, or file an issue on Github."
    fi
    exit 1
}

# Set variables
__conda_path="${arg_c:-${HOME}/miniconda3}"
__iguide_dir="${arg_s:-$(readlink -f ${__dir})}"
__iguide_env="${arg_e:-iguide}"
__update_lib=false
__update_env=false
if [[ "${arg_u}" = "all" || "${arg_u}" = "env" ]]; then
    __update_lib=true
    __update_env=true
elif [[ "${arg_u}" = "lib" ]]; then
    __update_lib=true
fi

__old_path=$PATH
__output=${2-/dev/stdout}


function __test_conda() {
    command -v conda &> /dev/null && echo true || echo false
}

function __detect_conda_install() {
    local discovered=$(__test_conda)
    if [[ $discovered = true ]]; then
      	local conda_path="$(which conda)"
        echo ${conda_path%'/bin/conda'}
    fi
}    

function __test_env() {
    if [[ $(__test_conda) = true ]]; then
        $(conda env list \
        | cut -f1 -d' ' \
        | grep -Fxq $__iguide_env > /dev/null) && \
        echo true || echo false
    else
      	echo false
    fi
}

function __test_iguide() {
    if [[ $(__test_env) = true ]]; then
      	activate_iguide
      	command -v iguide &> /dev/null && echo true || echo false
      	deactivate_iguide
    else
      	echo false
    fi
}

function activate_iguide () {
    set +o nounset
    source activate $__iguide_env
    set -o nounset
}

function deactivate_iguide () {
    set +o nounset
    source deactivate
    set -o nounset
}

function activate_iguide () {
    set +o nounset
    source activate $__iguide_env
    set -o nounset
}

function install_conda () {
    local tmpdir=$(mktemp -d)
    debug "Downloading miniconda..."
    debug_capture wget -nv ${__conda_url} -O ${tmpdir}/miniconda.sh 2>&1
    debug "Installing miniconda..."
    debug_capture bash ${tmpdir}/miniconda.sh -b -p ${__conda_path} 2>&1
    if [[ $(__test_conda) != true ]]; then
      	installation_error "Environment creation"
    fi
    rm ${tmpdir}/miniconda.sh
}

function install_environment () {
    debug_capture conda env update --name=$__iguide_env \
			  --quiet --file build.v0.2.0.yml
			  #			  --quiet --file bin/requirements.yml
		source activate $__iguide_env
		Rscript bin/setup.R >> ${__output}
    git clone https://github.com/cnobles/dualDemultiplexR.git tools/dualDemultiplexR >> ${__output}
    git clone https://github.com/cnobles/seqTrimR.git tools/seqTrimR >> ${__output}
    git clone https://github.com/cnobles/seqFiltR.git tools/seqFiltR >> ${__output}
    git clone https://github.com/cnobles/seqConsolidateR.git tools/seqConsolidateR >> ${__output}
    git clone https://github.com/cnobles/blatCoupleR.git tools/blatCoupleR >> ${__output}
    echo -e "iGUIDE successfully installed.\n" ;
    if [[ $(__test_env) != true ]]; then
      	installation_error "Environment creation"
    fi
}

function install_env_vars () {
    activate_iguide
    echo -ne "#/bin/sh\nexport IGUIDE_DIR=${__iguide_dir}" > \
	      ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
	  mkdir -p ${CONDA_PREFIX}/etc/conda/deactivate.d/
    echo -ne "#/bin/sh\nunset IGUIDE_DIR" > \
	      ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh
}

function install_iguidelib () {
    activate_iguide
    debug_capture pip install --upgrade $__iguide_dir 2>&1
    if [[ $(__test_iguide) != true ]]; then
      	installation_error "Library installation"
    fi
}

info "Starting iguide installation..."
info "    Conda path:   ${__conda_path}"
info "    iguide src:  ${__iguide_dir}"
info "    iguide env:  '${__iguide_env}'"

debug "Components detected:"
__conda_installed=$(__test_conda)
debug "    Conda:       ${__conda_installed}"
__env_exists=$(__test_env)
debug "    Environment: ${__env_exists}"
__iguide_installed=$(__test_iguide)
debug "    Library:     ${__iguide_installed}"

__env_changed=false


# Install Conda if necessary
if [[ $__conda_installed = true ]]; then
    if [[ $(__detect_conda_install) != $__conda_path ]]; then
        warning "Found pre-existing Conda installation in $(__detect_conda_install)".
        warning "Ignoring specified Conda path in favor of existing Conda install."
        __conda_path=$(__detect_conda_install)
    fi
    info "Conda already installed."
else
    info "Installing Conda..."
    install_conda
    __env_changed=true
fi


# Create Conda environment for iGuide
if [[ $__env_exists = true && $__update_env = false ]]; then
    info "Specified environment already exists (use '--update env' to update)"
else
    info "Creating iguide environment..."
    install_environment
    __env_changed=true
    info "iguide environment created."
fi




# Install iguidelib into environment if changed or requested
if [[ $__env_changed = true ]]; then
    info "Environment installed/updated; (re)installing iguide library..."
    install_iguidelib
elif [[ $__iguide_installed = false ]]; then
    info "Installing iguide library..."
    install_iguidelib
elif [[ $__update_lib = true ]]; then
    info "Updating iguide library..."
    install_iguidelib
else
    info "iguide library already installed (use '--update lib' to update)"
fi

# Always update the env_vars.sh in the iguide environment
debug "Updating \$IGUIDE_DIR variable to point to ${__iguide_dir}"
install_env_vars

# Check if on pre-existing path
if [[ $__old_path != *"${__conda_path}/bin"* ]]; then
    warning "** Conda was not detected on your PATH. **"
    warning "This is normal if you haven't installed Conda before."
    warning "To add it to your path, run "
    warning "   'echo \"export PATH=\$PATH:${__conda_path}/bin\" >> ~/.bashrc'"
    warning "and close and re-open your terminal session to apply."
    warning "When finished, run 'source activate ${__iguide_env}' to begin."
else
    info "Done. Run 'source activate ${__iguide_env}' to begin."
fi

echo -e "To get started, ensure ${__conda_path}/bin is in your path and\n" \
  "run 'source activate ${__iguide_env}'\n\n" \
  "To ensure ${__conda_path}/bin is in your path each time you long in,\n" \
  "append the following to your .bashrc or .bash_profile:\n\n" \
  "# Append miniconda3/bin to path\n" \
  "export PATH='~/miniconda3/bin:${PATH}'\n"
