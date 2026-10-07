#!/bin/bash
set -e

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export FIDASIM_DIR="${FIDASIM_DIR:-$(cd "${script_directory}/.." && pwd)}"

tests=(
  test_001
  test_002
  test_003
  test_004
  test_005
  test_006
)

for test_name in "${tests[@]}"; do
  echo
  echo "Running ${test_name}"
  (
    cd "${script_directory}/${test_name}"
    ./run.sh
  )
done

echo
echo "All regression tests passed."
