#!/usr/bin/env python3
"""Generate deterministic direct charge-exchange ion-sink references."""

from pathlib import Path
import sys

regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from reference_tools import run_reference_workflow


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: generate_reference.py <input_config.nml>")

    try:
        run_reference_workflow(config_filename=sys.argv[1])
    except (ConfigError, OSError, ValueError) as error:
        raise SystemExit(f"Reference calculation failed: {error}") from None


if __name__ == "__main__":
    main()
