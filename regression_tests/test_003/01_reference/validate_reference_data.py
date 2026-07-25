#!/usr/bin/env python3
"""Validate the shared Test 002 distributions referenced by Test 003."""

import argparse
from pathlib import Path
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from reference_generator_tools import validate_references


def main():
    # Set up command-line argument parsing:
    parser = argparse.ArgumentParser(
        description="Validate the Test 002 distributions used by Test 003."
    )
    parser.add_argument("config_path", help="Path to the namelist config file")
    args = parser.parse_args()

    try:
        result = validate_references(
            config_path=args.config_path,
        )
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None
    print(f"Particle parameters: {result['particle']}")
    for index, case in enumerate(result["cases"], start=1):
        print(
            f"{index:03d}: {case['input_path'].name} "
            f"R={case['selected_r']:.6g} cm "
            f"Z={case['selected_z']:.6g} cm "
            f"shape={case['shape']} denf={case['density']:.6e} ions/cm^3"
        )
    print(f"Validated {len(result['cases'])} shared reference files.")


if __name__ == "__main__":
    main()
