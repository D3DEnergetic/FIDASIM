#!/usr/bin/env python3
"""Generate the Test 002 single-location HDF5 reference files."""

import argparse
from pathlib import Path
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from reference_tools import generate_reference_files


def main():
    """Read the command line and run the Stage 1 extraction workflow."""
    parser = argparse.ArgumentParser(
        description="Extract compact single-location CQL3D F4D reference files."
    )
    parser.add_argument("config_path", help="Path to the namelist configuration")
    arguments = parser.parse_args()

    try:
        output_paths = generate_reference_files(config_path=arguments.config_path)
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None

    for output_path in output_paths:
        print(output_path)


if __name__ == "__main__":
    main()
