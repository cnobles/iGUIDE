# This code could be inserted directly into install.sh
# It exists in this separate script because it kept screwing up syntax highlighting in install.sh

# It must exist in the same directory as install.sh for installation purposes

# This section is necessary to source etc/b3bp.sh

read -r -d '' __usage <<-'EOF'
  -e --environment  [arg] Environment to install to. Default: "iguide"
  -s --iguide_dir  [arg] Location of iguide source code. Default: this directory
  -c --conda  [arg]       Location of Conda installation. Default: ${PREFIX}
  -u --update [arg]       Update iguide [lib]rary, conda [env], or [all].
  -v --verbose            Show subcommand output
  -d --debug              Run in debug mode.
  -h --help               Display this message and exit.
EOF

read -r -d '' __helptext <<-'EOF'
 This script installs or upgrades iguide, including Conda (if not installed).
 To upgrade, pass the '--upgrade all' option, then be sure to update your config
 files using 'iguide config update'.
EOF

# Load BASH3Boilerplate for command-line parsing and logging
source etc/b3bp.sh

