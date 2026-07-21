#!/bin/bash
set -e

cd 02_run_test
./run.sh

cd ../03_compare
./run.sh
