#!/usr/bin/env python3
"""Load and display the canonical basic example configuration."""

import argparse
from pathlib import Path
from pprint import pprint
import sys


# This example is run directly from the source tree. Add the directory that
# contains namelist_config_tools to sys.path without requiring installation.
package_parent_directory = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(package_parent_directory))

from namelist_config_tools import ConfigError, print_config

from config import read_config


def main():
    """Read the requested namelist and display its canonical representation."""
    parser = argparse.ArgumentParser(
        description="Validate the generic configuration-tools example."
    )
    parser.add_argument("config_path", help="Path to the example namelist")
    arguments = parser.parse_args()

    try:
        config = read_config(config_filename=arguments.config_path)
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None

    print_config(config)
    print("\nCanonical values")
    pprint(config, sort_dicts=False)


if __name__ == "__main__":
    main()
