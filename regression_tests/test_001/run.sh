#!/bin/bash
set -e

# Run from regression_tests/test_001 after building FIDASIM and activating
# FIDASIM_env. FIDASIM_DIR must be set.

cd 02_run_test
./test_001
cd ..

python3 01_reference/python/produce_reference.py 01_reference/python/input_config.nml
python3 02_run_test/python/postprocess.py 02_run_test/python/input_config.nml
python3 03_compare/python/compare_mean_free_paths.py \
    03_compare/python/input_config.nml
