"""Read Stage 3 settings and discover the files produced by Stages 1 and 2."""

from pathlib import Path

import numpy as np

from regression_test_tools import (
    ConfigError,
    normalize_path,
    read_namelist,
    require_boolean,
    require_existing_file,
    require_string,
    validate_schema,
)


CONFIG_SCHEMA = {
    "compare": {
        "required": True,
        "required_fields": ["run_config", "output_directory", "generate_plot"],
        "optional_fields": [],
    },
}


def _indexed_path(base_path, case_index):
    """Insert a three-digit case number before a path's extension."""
    return base_path.with_name(
        f"{base_path.stem}_{case_index:03d}{base_path.suffix}"
    )


def read_config(config_filename):
    """Read and validate the Stage 3 namelist.

    Args:
        config_filename (str or Path): Stage 3 input configuration.

    Returns:
        dict: Canonical configuration organized under the ``compare`` block.
    """
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)
    compare_block = blocks["compare"]

    run_config_value = require_string(
        value=compare_block["run_config"], field_label="run_config"
    )
    run_config = normalize_path(
        value=run_config_value,
        config_path=config_path,
        field_label="run_config",
    )
    require_existing_file(path=run_config, field_label="run_config")

    output_value = require_string(
        value=compare_block["output_directory"], field_label="output_directory"
    )
    output_directory = normalize_path(
        value=output_value,
        config_path=config_path,
        field_label="output_directory",
    )
    generate_plot = require_boolean(
        value=compare_block["generate_plot"], field_label="generate_plot"
    )

    return {
        "compare": {
            "run_config": run_config,
            "output_directory": output_directory,
            "generate_plot": generate_plot,
        }
    }


def discover_file_pairs(run_config):
    """Recover corresponding Stage 1 and Stage 2 files from their namelists."""
    run_config_path, run_blocks = read_namelist(config_path=run_config)
    if "input" not in run_blocks or "save_data_block" not in run_blocks:
        raise ConfigError("The Stage 2 configuration is missing required blocks.")

    run_input = run_blocks["input"]
    run_save = run_blocks["save_data_block"]
    if "reference_config" not in run_input:
        raise ConfigError("The Stage 2 input block is missing reference_config.")
    if "output_filename" not in run_save:
        raise ConfigError("The Stage 2 save_data_block is missing output_filename.")

    reference_config = normalize_path(
        value=run_input["reference_config"],
        config_path=run_config_path,
        field_label="reference_config",
    )
    require_existing_file(path=reference_config, field_label="reference_config")
    converted_base = normalize_path(
        value=run_save["output_filename"],
        config_path=run_config_path,
        field_label="Stage 2 output_filename",
    )

    reference_config_path, reference_blocks = read_namelist(
        config_path=reference_config
    )
    if "input" not in reference_blocks or "save_data_block" not in reference_blocks:
        raise ConfigError("The Stage 1 configuration is missing required blocks.")

    reference_input = reference_blocks["input"]
    reference_save = reference_blocks["save_data_block"]
    if "r_locations" not in reference_input:
        raise ConfigError("The Stage 1 input block is missing r_locations.")
    if "output_filename" not in reference_save:
        raise ConfigError("The Stage 1 save_data_block is missing output_filename.")

    reference_base = normalize_path(
        value=reference_save["output_filename"],
        config_path=reference_config_path,
        field_label="Stage 1 output_filename",
    )
    number_of_cases = len(np.atleast_1d(reference_input["r_locations"]))

    file_pairs = []
    for case_index in range(1, number_of_cases + 1):
        reference_path = _indexed_path(reference_base, case_index)
        converted_path = _indexed_path(converted_base, case_index)
        require_existing_file(
            path=reference_path,
            field_label=f"reference file {case_index:03d}",
        )
        require_existing_file(
            path=converted_path,
            field_label=f"converted file {case_index:03d}",
        )
        file_pairs.append(
            {"reference": Path(reference_path), "converted": Path(converted_path)}
        )

    return file_pairs
