#!/usr/bin/env python3
"""Generate deterministic direct charge-exchange ion-sink artifacts."""

from pathlib import Path
import sys

test_directory = Path(__file__).resolve().parents[1]
regression_tests_directory = test_directory.parent
sys.path.insert(0, str(test_directory))
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from deterministic_tools import run_deterministic_workflow


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: generate_deterministic.py <input_config.nml>")

    try:
        run_deterministic_workflow(config_filename=sys.argv[1])
    except (ConfigError, OSError, ValueError) as error:
        raise SystemExit(f"Deterministic calculation failed: {error}") from None


if __name__ == "__main__":
    main()
