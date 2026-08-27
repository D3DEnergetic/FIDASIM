#!/bin/bash
set -e

# Run from this directory after activating FIDASIM_env.
python3 produce_reference.py input_config.nml
