#!/usr/bin/env python3
"""Validate Stage 2 settings and write a normalized Fortran namelist."""

from pathlib import Path
import re
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import (
    ConfigError,
    normalize_path,
    read_namelist,
    reject_unknown_fields,
    require_boolean,
    require_blocks,
    require_existing_file,
    require_fields,
    require_integer,
    require_string,
)

import f90nml


RUN_TEST_REQUIRED_FIELDS = [
    "input_distribution_config",
    "output_directory",
    "n_samples",
    "seed",
]

RUN_TEST_OPTIONAL_FIELDS = [
    "comment",
    "plot_data",
]

MAX_REFERENCE_FILES = 256


def read_config(config_filename):
    """Read and validate the configuration used by the Fortran sampler.

    Args:
        config_filename (str or Path): User-facing Stage 2 namelist file.

    Returns:
        dict: Canonical configuration containing a normalized ``run_test``
        block. All paths in this block are absolute.

    Raises:
        ConfigError: If required settings are missing or invalid.
    """
    config_path, blocks = read_namelist(
        config_path=config_filename,
    )
    require_blocks(
        blocks=blocks,
        required_blocks=["run_test"],
    )

    run_block = blocks["run_test"]
    require_fields(
        block=run_block,
        required_fields=RUN_TEST_REQUIRED_FIELDS,
        block_label="&run_test",
    )

    allowed_fields = []
    for field_name in RUN_TEST_REQUIRED_FIELDS:
        allowed_fields.append(field_name)
    for field_name in RUN_TEST_OPTIONAL_FIELDS:
        allowed_fields.append(field_name)

    reject_unknown_fields(
        block=run_block,
        allowed_fields=allowed_fields,
        block_label="&run_test",
    )

    n_samples = require_integer(
        value=run_block["n_samples"],
        field_label="n_samples",
    )
    if n_samples < 1:
        raise ConfigError("n_samples must be positive.")

    seed = require_integer(
        value=run_block["seed"],
        field_label="seed",
    )
    if seed < 1:
        raise ConfigError("seed must be positive.")

    plot_data = require_boolean(
        value=run_block.get("plot_data", False),
        field_label="plot_data",
    )

    distribution_config_value = require_string(
        value=run_block["input_distribution_config"],
        field_label="input_distribution_config",
    )
    distribution_config = normalize_path(
        value=distribution_config_value,
        config_path=config_path,
        field_label="input_distribution_config",
    )
    require_existing_file(
        path=distribution_config,
        field_label="input_distribution_config",
    )
    reference_files = _discover_reference_files(distribution_config)
    n_reference_files = len(reference_files)
    if n_reference_files > MAX_REFERENCE_FILES:
        raise ConfigError(
            f"No more than {MAX_REFERENCE_FILES} reference files are supported."
        )

    output_directory_value = require_string(
        value=run_block["output_directory"],
        field_label="output_directory",
    )
    output_directory = normalize_path(
        value=output_directory_value,
        config_path=config_path,
        field_label="output_directory",
    )

    return {
        "run_test": {
            "n_reference_files": n_reference_files,
            "reference_files": reference_files,
            "output_directory": str(output_directory),
            "n_samples": n_samples,
            "seed": seed,
            "plot_data": plot_data,
        }
    }


def _discover_reference_files(distribution_config):
    """Discover contiguous Test 002 Stage 2 outputs from its configuration."""
    distribution_path, blocks = read_namelist(config_path=distribution_config)
    require_blocks(blocks=blocks, required_blocks=["save_data_block"])
    save_block = blocks["save_data_block"]
    require_fields(
        block=save_block,
        required_fields=["output_filename"],
        block_label="Test 002 &save_data_block",
    )
    output_value = require_string(
        value=save_block["output_filename"],
        field_label="Test 002 output_filename",
    )
    output_base = normalize_path(
        value=output_value,
        config_path=distribution_path,
        field_label="Test 002 output_filename",
    )

    pattern = re.compile(
        rf"{re.escape(output_base.stem)}_(\d{{3}})"
        rf"{re.escape(output_base.suffix)}"
    )
    indexed_paths = []
    for path in output_base.parent.glob(
        f"{output_base.stem}_*{output_base.suffix}"
    ):
        match = pattern.fullmatch(path.name)
        if match is not None and path.is_file():
            indexed_paths.append((int(match.group(1)), path.resolve()))
    indexed_paths.sort(key=lambda item: item[0])

    indices = [index for index, _ in indexed_paths]
    if not indexed_paths:
        raise ConfigError(
            "No Test 002 Stage 2 outputs were found. Run "
            f"'{distribution_path.parent}/run.sh {distribution_path.name}' first."
        )
    if indices != list(range(1, len(indices) + 1)):
        raise ConfigError(
            "Test 002 output indices must be contiguous and start at 001."
        )
    return [str(path) for _, path in indexed_paths]


def write_config(config, output_filename):
    """Write the canonical sampler configuration as a Fortran namelist.

    Args:
        config (dict): Canonical configuration returned by ``read_config``.
        output_filename (str or Path): Destination for the generated namelist.

    Returns:
        Path: Absolute path of the generated namelist.
    """
    output_path = Path(output_filename).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    f90nml.write(config, output_path, force=True)
    return output_path


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: normalize_config.py <input_config.nml> "
            "<normalized_config.nml>"
        )

    try:
        config = read_config(
            config_filename=sys.argv[1],
        )
        output_path = write_config(
            config=config,
            output_filename=sys.argv[2],
        )
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None
    except OSError as error:
        raise SystemExit(f"Cannot write normalized configuration: {error}") from None

    print(f"Wrote normalized configuration: {output_path}")


if __name__ == "__main__":
    main()
