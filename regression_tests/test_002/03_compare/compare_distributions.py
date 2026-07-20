#!/usr/bin/env python3
"""Command-line entry point for the test_002 comparison workflow."""

import argparse
from pathlib import Path
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from comparison_tools import run_comparison


def main():
    # Read the path to the comparison configuration from the command line.
    parser = argparse.ArgumentParser(
        description="Compare reference and sampled energy-pitch distributions."
    )
    parser.add_argument("config_path", help="Path to the comparison namelist")
    arguments = parser.parse_args()

    # The package API performs the complete comparison workflow.
    try:
        run_comparison(config_path=arguments.config_path)
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None


if __name__ == "__main__":
    main()
