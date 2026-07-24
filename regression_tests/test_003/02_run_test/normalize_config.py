#!/usr/bin/env python3
"""Validate Stage 2 settings and write a normalized Fortran namelist."""

from pathlib import Path
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import (
    ConfigError,
    as_list,
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
    "n_reference_files",
    "reference_files",
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

    # Validate the scalar sampling controls.
    n_reference_files = require_integer(
        value=run_block["n_reference_files"],
        field_label="n_reference_files",
    )
    if n_reference_files < 1 or n_reference_files > MAX_REFERENCE_FILES:
        raise ConfigError(
            "n_reference_files must be from 1 through "
            f"{MAX_REFERENCE_FILES}."
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

    # Normalize and verify every reference path before starting Fortran.
    configured_reference_files = as_list(run_block["reference_files"])
    if len(configured_reference_files) != n_reference_files:
        raise ConfigError(
            "n_reference_files must equal the number of entries in "
            "reference_files."
        )

    reference_files = []
    for file_index, configured_path in enumerate(configured_reference_files):
        field_label = f"reference_files({file_index + 1})"
        path_value = require_string(
            value=configured_path,
            field_label=field_label,
        )
        reference_path = normalize_path(
            value=path_value,
            config_path=config_path,
            field_label=field_label,
        )
        require_existing_file(
            path=reference_path,
            field_label=field_label,
        )
        reference_files.append(str(reference_path))

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
