#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
  echo "Usage: ./run.sh <input_config.nml>" >&2
  exit 2
fi

python3 generate_deterministic.py "$1"
