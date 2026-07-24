#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: ./run.sh <input_config.nml>"
    exit 1
fi

config_file="$1"
normalized_config="build/normalized_input_config.nml"

python3 normalize_config.py "${config_file}" "${normalized_config}"
./test_003 "${normalized_config}"
python3 plot_sampled_data.py "${config_file}"
