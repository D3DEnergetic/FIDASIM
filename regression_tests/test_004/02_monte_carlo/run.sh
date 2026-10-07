#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
  echo "Usage: ./run.sh <input_config.nml>" >&2
  exit 2
fi

config_file="$1"
normalized_config="build/normalized_input_config.nml"

if [ ! -x ./test_004 ]; then
  echo "Test 004 executable not found. Run 'make regression_tests' from the FIDASIM repository root first." >&2
  exit 2
fi

python3 normalize_config.py "${config_file}" "${normalized_config}"
./test_004 "${normalized_config}"

python3 postprocess_monte_carlo.py "${config_file}"
