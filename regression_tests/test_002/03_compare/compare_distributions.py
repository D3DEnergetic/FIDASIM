#!/usr/bin/env python3
"""Command-line entry point for the test_002 comparison workflow."""

import argparse

from comparison_tools import run_comparison


def main():
    # Read the path to the comparison configuration from the command line.
    parser = argparse.ArgumentParser(
        description="Compare reference and sampled energy-pitch distributions."
    )
    parser.add_argument("config_path", help="Path to the comparison namelist")
    arguments = parser.parse_args()

    # The package API performs the complete comparison workflow.
    run_comparison(arguments.config_path)


if __name__ == "__main__":
    main()
