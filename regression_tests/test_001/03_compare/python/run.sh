#!/bin/bash
set -e

# Run from this directory after activating FIDASIM_env.
python3 compare_mean_free_paths.py input_config.nml
