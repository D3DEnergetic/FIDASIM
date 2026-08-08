"""Read and validate the basic example configuration."""

import math

from namelist_config_tools import (
    ConfigError,
    as_list,
    normalize_path,
    normalize_string,
    read_namelist,
    require_blocks,
    require_boolean,
    require_choice,
    require_existing_file,
    require_integer,
    require_real,
    require_string,
    validate_schema,
)


CONFIG_SCHEMA = {
    "run": {
        "required": True,
        "required_fields": [
            "input_filename",
            "mode",
            "n_steps",
            "tolerance",
            "write_output",
        ],
        "optional_fields": ["comment", "labels"],
    },
    "output": {
        "required": False,
        "required_fields": ["output_directory"],
        "optional_fields": [],
    },
}

SUPPORTED_MODES = ["fast", "accurate"]


def _normalize_labels(value):
    """Normalize an optional scalar or sequence into nonempty label strings.

    A missing value becomes an empty list, while one label or a sequence of
    labels becomes a list whose entries are validated as nonempty strings.
    """
    labels = []
    for label_index, raw_label in enumerate(as_list(value), start=1):
        label = require_string(
            value=raw_label,
            field_label=f"labels entry {label_index}",
        )
        labels.append(label)
    return labels


def read_config(config_filename):
    """Read the example namelist and return its canonical configuration.

    Returns:
        dict: Validated settings with the following schema::

            {
                "run": {
                    "comment": str,
                    "config_path": Path,
                    "input_filename": Path,
                    "mode": str,
                    "n_steps": int,
                    "tolerance": float,
                    "labels": list[str],
                    "write_output": bool,
                },
                "output": {
                    "output_directory": Path or None,
                },
            }

    Raises:
        ConfigError: If the namelist does not satisfy the example contract.
    """
    # Parse the namelist first. ``config_path`` is the absolute path of the
    # file that was read, and ``blocks`` contains its namelist blocks and
    # variables as ordinary Python dictionaries.
    config_path, blocks = read_namelist(config_path=config_filename)

    # Check the overall structure against CONFIG_SCHEMA. This verifies that
    # all required blocks and variables are present and rejects block or
    # variable names that the application does not recognize.
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    # The schema guarantees that the required ``run`` block now exists. Read
    # its settings and convert each value to the type used by the application.
    run_block = blocks["run"]
    write_output = require_boolean(
        value=run_block["write_output"],
        field_label="write_output",
    )

    # The output block is structurally optional but becomes required when the
    # application is asked to write output.
    if write_output:
        require_blocks(blocks=blocks, required_blocks=["output"])

    # Namelist paths are strings. Convert the input path into an absolute Path
    # relative to the configuration file, then verify that its target exists.
    input_filename = normalize_path(
        value=require_string(
            value=run_block["input_filename"],
            field_label="input_filename",
        ),
        config_path=config_path,
        field_label="input_filename",
    )
    require_existing_file(
        path=input_filename,
        field_label="input_filename",
    )

    # Normalize the mode's spelling before checking that it is one of the
    # choices supported by this example application.
    mode = require_choice(
        value=normalize_string(run_block["mode"]),
        supported_values=SUPPORTED_MODES,
        field_label="mode",
    )

    # Type validation alone is not sufficient for every setting. After
    # requiring an integer, apply the application's positive-value rule.
    n_steps = require_integer(
        value=run_block["n_steps"],
        field_label="n_steps",
    )
    if n_steps < 1:
        raise ConfigError("n_steps must be positive.")

    # Likewise, require a real number and then enforce the application's
    # stronger rule that the tolerance must be finite and positive.
    tolerance = require_real(
        value=run_block["tolerance"],
        field_label="tolerance",
    )
    if not math.isfinite(tolerance) or tolerance <= 0.0:
        raise ConfigError("tolerance must be finite and greater than zero.")

    # Supply an application default for an optional scalar setting. When the
    # user provides it, validate it in exactly the same way as required text.
    comment = ""
    if "comment" in run_block:
        comment = require_string(
            value=run_block["comment"],
            field_label="comment",
        )

    # Normalize the optional labels so downstream code always receives a
    # list, whether the namelist supplied no label, one label, or many labels.
    labels = _normalize_labels(run_block.get("labels"))

    # The output block may be absent when output is disabled. Preserve that
    # state as None; otherwise normalize its directory relative to this file.
    output_directory = None
    if "output" in blocks:
        output_directory = normalize_path(
            value=require_string(
                value=blocks["output"]["output_directory"],
                field_label="output_directory",
            ),
            config_path=config_path,
            field_label="output_directory",
        )

    # Return one canonical representation. Application code can now consume
    # these values without knowing about namelist syntax or repeating checks.
    return {
        "run": {
            "comment": comment,
            "config_path": config_path,
            "input_filename": input_filename,
            "mode": mode,
            "n_steps": n_steps,
            "tolerance": tolerance,
            "labels": labels,
            "write_output": write_output,
        },
        "output": {
            "output_directory": output_directory,
        },
    }
