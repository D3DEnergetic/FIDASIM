#!/bin/bash
set -e

config_file="${1:-input_config.nml}"
normalized_config="build/normalized_input_config.nml"

python3 normalize_config.py "${config_file}" "${normalized_config}"
./test_002 "${normalized_config}"
python3 plot_sampled_data.py "${config_file}"
