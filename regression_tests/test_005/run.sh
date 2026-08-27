#!/bin/bash

# Run the prebuilt Test 005 executable and write its comparison report.
set -e

mkdir -p output_data
./build/test_005
