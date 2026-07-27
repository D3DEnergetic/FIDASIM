#!/usr/bin/env python3
"""Run the Test 002 CQL3D-to-FIDASIM conversion workflow."""

import argparse
from pathlib import Path
import sys

regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from conversion_tools import convert_reference_distributions_to_fidasim
from regression_test_tools import ConfigError


def main():
    """Read the command line and convert all configured reference files."""
    parser = argparse.ArgumentParser(description="Convert Test 002 references.")
    parser.add_argument("config_path", help="Path to the Stage 2 namelist")
    arguments = parser.parse_args()
    try:
        convert_reference_distributions_to_fidasim(
            config_path=arguments.config_path
        )
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None


if __name__ == "__main__":
    main()
