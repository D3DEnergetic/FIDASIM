#!/usr/bin/env python3
"""Command-line entry point for the Test 002 moment comparison."""

import argparse
from pathlib import Path
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from comparison_tools import run_comparison


def main():
    """Read the requested configuration and run the complete comparison."""
    parser = argparse.ArgumentParser(
        description="Compare CQL3D reference and converted FIDASIM moments."
    )
    parser.add_argument("config_path", help="Path to the comparison namelist")
    arguments = parser.parse_args()

    try:
        all_passed = run_comparison(config_path=arguments.config_path)
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None

    if not all_passed:
        raise SystemExit(
            "Test 002 comparison failed its configured moment tolerance."
        )


if __name__ == "__main__":
    main()
