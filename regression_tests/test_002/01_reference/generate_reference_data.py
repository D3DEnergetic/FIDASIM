#!/usr/bin/env python3
"""Entry-point script for generating reference data for regression tests."""

import argparse

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
    outputs = generate_outputs(args.config_path)
    for output in outputs:
        print(output)


if __name__ == "__main__":
    main()
