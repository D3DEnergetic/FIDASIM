#!/usr/bin/env python3
"""Command-line entry point for the Test 001 scientific comparison.

This stage compares the independently generated ADF21 reference with the
post-processed FIDASIM result, persists its report and diagnostic figures, and
uses its exit status to expose the regression verdict to batch runners.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path


# The shared regression configuration helpers live three directories above
# this stage entry point. Add that stable repository location before importing
# the local configuration and workflow modules.
REGRESSION_TESTS_DIRECTORY = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REGRESSION_TESTS_DIRECTORY))

# Keep Matplotlib's writable cache with the ignored Python artifacts. This is
# set before importing the workflow because Matplotlib resolves it at import.
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parent / "output_data" / ".matplotlib"),
)

from regression_test_tools import ConfigError

from comparison_config import read_comparison_config
from comparison_workflow import run_comparison


def main(arguments: list[str]) -> int:
    """Run the configured comparison and return its regression status.

    Args:
        arguments: Command-line arguments including the program name.

    Returns:
        Zero when the regression passes, one when it fails or cannot run, and
        two when the command line is invalid.
    """

    if len(arguments) > 2:
        print("Usage: compare_mean_free_paths.py [input_config.nml]", file=sys.stderr)
        return 2

    if len(arguments) == 2:
        config_path = Path(arguments[1])
    else:
        config_path = Path(__file__).with_name("input_config.nml")

    try:
        configuration = read_comparison_config(config_path)
        regression_passed = run_comparison(configuration)
    except (ConfigError, OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1

    return 0 if regression_passed else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
