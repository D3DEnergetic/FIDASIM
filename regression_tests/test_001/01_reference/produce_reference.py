"""Command-line entry point for Test 001 reference generation.

The scientific work lives in :mod:`reference_workflow`; this file only checks
the command line, converts the namelist argument to a path, and returns a
conventional process status.
"""

from pathlib import Path
import sys

from reference_workflow import run_reference_workflow


def main(arguments: list[str]) -> int:
    """Run the reference workflow named on the command line.

    Args:
        arguments: Command-line arguments excluding the program name.

    Returns:
        Zero after a successful reference calculation.
    """
    if len(arguments) != 1:
        print("Usage: produce_reference.py <input_config.nml>", file=sys.stderr)
        return 2

    run_reference_workflow(Path(arguments[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
