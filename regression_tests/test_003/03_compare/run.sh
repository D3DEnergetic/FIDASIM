#!/bin/bash
set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: ./run.sh <input_config.nml>"
    exit 1
fi

python3 compare_distributions.py "$1"
