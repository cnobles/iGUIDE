#!/usr/bin/env bash

__conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

read -r -d '' __usage <<-'EOF'
  -e --environment  [arg] Environment to install to. Default: "iguide"
  -s --iguide_dir   [arg] Location of iguide source code. Default: this directory
  -c --conda  [arg]       Location of Conda installation. Default: ${PREFIX}
  -u --update [arg]       Update iguide [lib]rary, conda [env], or [all].
  -r --requirements       Install from requirements rather than build (slow).
  -t --test               After installation, run test to check functionality.
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script installs or upgrades iGUIDE, including Conda (if not installed).
 To upgrade, pass the '--upgrade all' option, then be sure to update your config
 files using 'iguide config update'.
EOF


# Load BASH3Boilerplate for command-line parsing and logging
source etc/b3bp.sh

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
__run_iguide_tests=false
__reqs_install=false
__update_lib=false
__update_env=false
__req_r_version="3.4.1"
__old_path=$PATH
__output=${2-/dev/stdout}

if [[ "${arg_t:?}" = "1" ]]; then
    __run_iguide_tests=true
fi

if [[ "${arg_r:?}" = "1" ]]; then
    __reqs_install=true
fi

if [[ "${arg_u}" = "all" || "${arg_u}" = "env" ]]; then
    __update_lib=true
    __update_env=true
elif [[ "${arg_u}" = "lib" ]]; then
    __update_lib=true
fi

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

function __test_r_version () {
    activate_iguide
    
    local sem_version=$(R --version | grep 'R version' | cut -d ' ' -f 3)
    
    if (( ${sem_version//./} >= ${__req_r_version//./} )); then
        echo true
    else
        echo false
    fi
    
    deactivate_iguide
}

function __test_r_packages () {
    activate_iguide
    
    $(Rscript ${__iguide_dir}/etc/check_for_required_packages.R \
        > /dev/null) && echo true || echo false
    
    deactivate_iguide
}

function __test_iguidelib() {
    if [[ $(__test_env) = true ]]; then
      	activate_iguide
      	command -v iguide &> /dev/null && echo true || echo false
      	deactivate_iguide
    else
      	echo false
    fi
}

function __test_iguide() {
    if [[ $(__test_env) = true ]]; then
      	$(bash ${__iguide_dir}/tests/test.sh ${__iguide_env} &> /dev/null) && \
      	    echo "pass" || echo "fail"
    else
      	echo "fail"
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
    if [[ $__reqs_install == "true" ]]; then
        local install_options="--quiet --file etc/requirements.yml"
    else
        local install_options="--quiet --file etc/build.v0.3.0.yml"
    fi
    
    debug_capture conda env update --name=$__iguide_env ${install_options} 2>&1
    
    if [[ $(__test_env) != true ]]; then
      	installation_error "Environment creation"
    fi
		
		if [[ $(__test_r_version) != true ]]; then
        installation_error "Insufficient R-program version"
    fi
    
    if [[ $(__test_r_packages) != true ]]; then
        installation_error "R-package installation"
    fi
    
    echo -e "iGUIDE installed.\n"
}

function install_env_vars () {
    activate_iguide
    
    echo -ne "#/bin/sh\nexport IGUIDE_DIR=${__iguide_dir}" > \
	      ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh
	  
	  mkdir -p ${CONDA_PREFIX}/etc/conda/deactivate.d/
    
    echo -ne "#/bin/sh\nunset IGUIDE_DIR" > \
	      ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh
	 
	  deactivate_iguide
}

function install_iguidelib () {
    activate_iguide
    
    debug_capture pip install --upgrade ${__iguide_dir}/tools/iguidelib/ 2>&1
    
    if [[ $(__test_iguidelib) != true ]]; then
      	installation_error "Library installation"
    fi
    
    deactivate_iguide
}

info "Starting iGUIDE installation..."
info "    Conda path:  ${__conda_path}"
info "    iGUIDE src:  ${__iguide_dir}"
info "    iGUIDE env:  '${__iguide_env}'"

debug "Components detected:"
__conda_installed=$(__test_conda)
debug "    Conda:         ${__conda_installed}"
__env_exists=$(__test_env)
debug "    Environment:   ${__env_exists}"
__iguidelib_installed=$(__test_iguidelib)
debug "    Library:       ${__iguidelib_installed}"

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


# Create Conda environment for iGUIDE
if [[ $__env_exists = true && $__update_env = false ]]; then
    info "Specified environment already exists (use '--update env' to update)"
else
    if [[ $__reqs_install = "true" ]]; then
        __build_source="etc/requirements.yml"
    else
        __build_source="etc/build.v0.3.0.yml"
    fi

    info "Creating iGUIDE environment..."
    info "    Building from: $__build_source"
    install_environment
    __env_changed=true
    info "$__iguide_env environment created."
fi


# Always update the env_vars.sh in the iGUIDE environment
debug "Updating \$IGUIDE_DIR variable to point to ${__iguide_dir}"
info "Setting environmental variables..."
install_env_vars


# Install iguidelib into environment if changed or requested
if [[ $__env_changed = true ]]; then
    info "Environment installed/updated; (re)installing iGUIDE library..."
    install_iguidelib
elif [[ $__iguidelib_installed = false ]]; then
    info "Installing iGUIDE library..."
    install_iguidelib
elif [[ $__update_lib = true ]]; then
    info "Updating iGUIDE library..."
    install_iguidelib
else
    info "iGUIDE library already installed (use '--update lib' to update)"
fi


# Run tests if requested
if [[ $__run_iguide_tests = true ]]; then
    info "Running iGUIDE tests...(this may take a few mins)"
    __iguide_tested=$(__test_iguide)
    
    if [[ $__iguide_tested = "pass" ]]; then
        info "    iGUIDE Tests:  ${__iguide_tested}"
    else
        warning "    iGUIDE Tests:  ${__iguide_tested}"
    fi
fi


# Check if on pre-existing path
if [[ $__old_path != *"${__conda_path}/bin"* ]]; then
    warning "** Conda was not detected on your PATH. **"
    warning "This is normal if you have not installed Conda before."
    warning "To add it to your path, run "
    warning "   'echo \"export PATH=\$PATH:${__conda_path}/bin\" >> ~/.bashrc'"
    warning "and close and re-open your terminal session to apply."
    warning "When finished, run 'source activate ${__iguide_env}' to begin."
else
    info "Done. Run 'source activate ${__iguide_env}' to begin."
fi


echo -e "To get started, ensure ${__conda_path}/bin is in your path and\n" \
  "run 'source activate ${__iguide_env}'\n\n"
  
if [[ $__old_path != *"${__conda_path}/bin"* ]]; then
  echo "To ensure ${__conda_path}/bin is in your path each time you login,\n" \
  "append the following to your .bashrc or .bash_profile:\n\n" \
  "# Append miniconda3/bin to path\n" \
  "export PATH='~/miniconda3/bin:${PATH}'\n"
fi
