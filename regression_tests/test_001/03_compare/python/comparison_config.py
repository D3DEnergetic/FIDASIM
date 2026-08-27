"""Read and validate the Test 001 Stage 3 comparison namelist."""

import math
from collections.abc import Mapping
from pathlib import Path

from regression_test_tools import (
    ConfigError,
    normalize_path,
    read_namelist,
    require_boolean,
    require_existing_file,
    require_real,
    require_string,
    validate_schema,
)


CONFIG_SCHEMA = {
    "compare": {
        "required": True,
        "required_fields": [
            "reference_file",
            "postprocessed_file",
            "raw_test_file",
            "output_directory",
            "generate_plots",
            "mfp_relative_error_tolerance",
        ],
        "optional_fields": [],
    },
}


def _configured_path(
    block: Mapping[str, object],
    field_name: str,
    config_path: Path,
) -> Path:
    """Resolve one required path from the comparison block.

    Args:
        block: Mapping-like namelist block containing the requested field.
        field_name: Name of the field to read and resolve.
        config_path: Absolute path of the namelist that owns the field.

    Returns:
        The absolute path resolved relative to ``config_path``.

    Raises:
        ConfigError: If the field is not a nonempty string or path.
    """
    value = require_string(
        value=block[field_name],
        field_label=field_name,
    )
    return normalize_path(
        value=value,
        config_path=config_path,
        field_label=field_name,
    )


def read_comparison_config(
    config_filename: str | Path,
) -> dict[str, object]:
    """Read and validate the complete Stage 3 configuration.

    Every configured path is resolved relative to the namelist rather than
    the caller's working directory. Input artifacts must already exist. The
    output directory is allowed to be absent because the workflow creates it.

    Args:
        config_filename: Stage 3 Fortran namelist path.

    Returns:
        A flat dictionary containing absolute ``Path`` objects, the plotting
        flag, and the finite positive relative-error tolerance.

    Raises:
        ConfigError: If the namelist structure, a field value, or an input
            artifact is invalid.
    """
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    compare_block = blocks["compare"]
    reference_file = _configured_path(
        block=compare_block,
        field_name="reference_file",
        config_path=config_path,
    )
    postprocessed_file = _configured_path(
        block=compare_block,
        field_name="postprocessed_file",
        config_path=config_path,
    )
    raw_test_file = _configured_path(
        block=compare_block,
        field_name="raw_test_file",
        config_path=config_path,
    )
    output_directory = _configured_path(
        block=compare_block,
        field_name="output_directory",
        config_path=config_path,
    )

    require_existing_file(
        path=reference_file,
        field_label="reference_file",
    )
    require_existing_file(
        path=postprocessed_file,
        field_label="postprocessed_file",
    )
    require_existing_file(
        path=raw_test_file,
        field_label="raw_test_file",
    )

    generate_plots = require_boolean(
        value=compare_block["generate_plots"],
        field_label="generate_plots",
    )
    tolerance = require_real(
        value=compare_block["mfp_relative_error_tolerance"],
        field_label="mfp_relative_error_tolerance",
    )
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ConfigError(
            "mfp_relative_error_tolerance must be finite and greater than zero."
        )

    return {
        "reference_file": reference_file,
        "postprocessed_file": postprocessed_file,
        "raw_test_file": raw_test_file,
        "output_directory": output_directory,
        "generate_plots": generate_plots,
        "mfp_relative_error_tolerance": tolerance,
    }
