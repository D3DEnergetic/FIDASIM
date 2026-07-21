#!/usr/bin/env python3
"""Entry-point script for generating reference data for regression tests."""

import argparse
from pathlib import Path
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from reference_generator_tools import generate_outputs


def main():
    # Set up command-line argument parsing:
    parser = argparse.ArgumentParser(
        description=(
            "Export canonical energy-pitch HDF5 slices from a supported "
            "4D distribution."
        )
    )
    parser.add_argument("config_path", help="Path to the namelist config file")
    args = parser.parse_args()

    # Generate the reference outputs based on the provided configuration file:
    try:
        outputs = generate_outputs(
            config_path=args.config_path,
        )
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
