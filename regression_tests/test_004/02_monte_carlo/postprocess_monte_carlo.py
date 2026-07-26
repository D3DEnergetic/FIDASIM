#!/usr/bin/env python3
"""Postprocess Test 004 production sink files."""

from pathlib import Path
import sys

test_directory = Path(__file__).resolve().parents[1]
regression_tests_directory = test_directory.parent
sys.path.insert(0, str(test_directory))
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import ConfigError
from monte_carlo_tools import postprocess_monte_carlo


def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: postprocess_monte_carlo.py <input_config.nml>"
        )
    try:
        postprocess_monte_carlo(sys.argv[1])
    except (ConfigError, OSError, ValueError) as error:
        raise SystemExit(f"Monte Carlo postprocessing error: {error}") from None


if __name__ == "__main__":
    main()
