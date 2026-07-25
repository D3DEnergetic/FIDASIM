#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
  echo "Usage: ./run.sh <input_config.nml>" >&2
  exit 2
fi

config_file="$1"
normalized_config="build/normalized_input_config.nml"

python3 normalize_config.py "${config_file}" "${normalized_config}"

# The Test 004 makefile is in the parent directory.
make -C ..

./test_004 "${normalized_config}"
