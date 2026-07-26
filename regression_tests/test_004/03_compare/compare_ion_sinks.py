#!/usr/bin/env python3
"""Command-line entry point for Test 004 Stage 3."""

import argparse
from pathlib import Path
import sys

# Make the shared regression-test and Test 004 packages importable without
# installing either package into the active Python environment.
test_004_directory = Path(__file__).resolve().parents[1]
regression_tests_directory = test_004_directory.parent
sys.path.insert(0, str(test_004_directory))
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from comparison_tools import run_comparison


def main():
    """Run the configured comparison and return a regression-test exit code."""
    parser = argparse.ArgumentParser(
        description="Compare deterministic and Monte Carlo Test 004 sinks."
    )
    parser.add_argument(
        "config_path",
        help="Path to the Stage 3 comparison namelist",
    )
    arguments = parser.parse_args()

    try:
        comparison = run_comparison(arguments.config_path)
    except (ConfigError, OSError, ValueError) as error:
        raise SystemExit(f"Test 004 comparison error: {error}") from None

    if not comparison.all_passed:
        raise SystemExit(
            "Test 004 comparison failed its configured rate tolerances. "
            f"See {comparison.report_filename}."
        )


if __name__ == "__main__":
    main()
