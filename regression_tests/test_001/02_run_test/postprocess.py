#!/usr/bin/env python3
"""Command-line entry point for Test 001 FIDASIM post-processing.

The source ``test_001.h5`` is treated as a read-only regression artifact.  All
validation, mean-free-path calculation, metadata writing, and plotting are
implemented in :mod:`postprocess_support`.
"""

from __future__ import annotations

import sys
from pathlib import Path

from postprocess_support import run_postprocessing


def main(arguments: list[str]) -> int:
    """Run post-processing from one optional namelist configuration path.

    Args:
        arguments: Command-line arguments including the program name.

    Returns:
        Zero after successful processing, one after a processing error, or two
        when too many command-line arguments are supplied.
    """
    if len(arguments) > 2:
        print("Usage: postprocess.py [input_config.nml]", file=sys.stderr)
        return 2

    if len(arguments) == 2:
        config_path = Path(arguments[1])
    else:
        config_path = Path(__file__).with_name("input_config.nml")

    try:
        run_postprocessing(config_path)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
