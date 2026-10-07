#!/bin/bash
set -e

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_directory}"

if [ ! -x 02_monte_carlo/test_004 ]; then
  echo "Test 004 executable not found. Run 'make regression_tests' from the FIDASIM repository root first." >&2
  exit 2
fi

echo "Running Test 004 deterministic calculation"
(
  cd 01_deterministic
  ./run.sh ../input_config_A.nml
)

echo "Running Test 004 Monte Carlo calculation"
(
  cd 02_monte_carlo
  ./run.sh ../input_config_A.nml
)

echo "Running Test 004 comparison"
(
  cd 03_compare
  ./run.sh input_config_A.nml
)
