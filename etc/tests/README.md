# Tests
This directory contains test data to verify that the iGUIDE software is functioning correctly.

To generate a test dataset, use the associated Rscript found in the construct_scripts directory (`simulate_incorp_data.R`). 

```
# Make sure you are within the appropriate environment:
source activate iguide
cd $IGUIDE_DIR

# For information about the script:
Rscript etc/tests/construct_scripts/simulate_incorp_data.R -h

# To learn more about the parameters involved in the dataset generation:
# Please read the sim_config.yml
cat etc/tests/construct_scripts/sim_config.yml

# To generate new data:
Rscript etc/tests/construct_scripts/simulate_incorp_data.R etc/tests/construct_scripts/sim_config.yml
# or
Rscript etc/tests/construct_scripts/simulate_incorp_data.R etc/tests/construct_scripts/sim_config.yml \
    -o {output_folder_name} -s {numeric_seed}
```

An already generated dataset for testing has been deposited here and is used for the 'default' tests.

