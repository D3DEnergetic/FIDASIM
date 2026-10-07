#!/bin/bash
set -e

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${script_directory}"

artifact_file="${ARTIFACT_FILE:-fixtures/tiny_mesh.h5}"
behavior_file="${BEHAVIOR_FILE:-fixtures/mesh_behavior.nml}"
executable="build/test_006"

if [ ! -f "${artifact_file}" ]; then
  echo "Mesh artifact not found: ${artifact_file}" >&2
  exit 2
fi

if [ ! -f "${behavior_file}" ]; then
  echo "Mesh behavior file not found: ${behavior_file}" >&2
  exit 2
fi

if [ ! -x "${executable}" ]; then
  echo "Test 006 executable not found. Run 'make regression_tests' from the FIDASIM repository root first." >&2
  exit 2
fi

mkdir -p output_data
./"${executable}" "${artifact_file}" "${behavior_file}"
