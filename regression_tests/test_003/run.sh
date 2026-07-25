#!/bin/bash
set -e

# Select the configuration files to run the test with.
# Multiple configuration files can be specified in the array below:
config_files=(
    input_config_A.nml
)

# Run the test for each configuration file specified in the array.
for config_file in "${config_files[@]}"; do
    echo
    echo "Running Test 003: $config_file"

    cd 01_reference
    ./run.sh "$config_file"

    cd ../02_run_test
    ./run.sh "$config_file"

    cd ../03_compare
    ./run.sh "$config_file"

    cd ..
done
