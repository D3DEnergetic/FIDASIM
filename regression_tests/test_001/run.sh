#!/bin/bash
set -e

# Run from regression_tests/test_001 after building FIDASIM and activating
# FIDASIM_env. FIDASIM_DIR must be set.

cd 02_run_test
./test_001
cd ..

python3 01_reference/produce_reference.py 01_reference/input_config.nml
python3 02_run_test/postprocess.py 02_run_test/input_config.nml
python3 03_compare/compare_mean_free_paths.py 03_compare/input_config.nml
