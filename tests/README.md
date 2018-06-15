# Tests
This directory contains tests used by travis-ci and SnakeMake to verify that the iGUIDE software is functioning correctly.

To generate a test dataset, use the associated Rscript found in this directory (`simulate_incorp_data.R`). 

```
# Make sure you are within the appropriate environment:
source activate iguide

# For information about the script:
Rscript simulate_incorp_data.R -h

# To learn more about the parameters involved in the dataset generation:
# Please read the sim_config.yml
cat sim_config.yml

# To generate new data:
Rscript simulate_incorp_data.R sim_config.yml
# or
Rscript simulate_incorp_data.R sim_config.yml -o {output_folder_name} -s {numeric_seed}
```

An already generated dataset for testing has been deposited here and is used for the 'default' tests.

