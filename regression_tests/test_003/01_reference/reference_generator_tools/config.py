"""Read and validate the Test 003 Stage 1 configuration."""

from regression_test_tools import (
    ConfigError,
    normalize_path,
    normalize_string,
    read_namelist,
    require_boolean,
    require_blocks,
    require_choice,
    require_existing_file,
    require_real,
    require_string,
    validate_schema,
)


CONFIG_SCHEMA = {
    "input": {
        "required": True,
        "required_fields": ["input_config"],
        "optional_fields": ["comment", "plot_data", "save_data"],
    },
    "plot_data_block": {
        "required": False,
        "required_fields": [],
        "optional_fields": [
            "scale",
            "fmin",
            "fmax",
            "enable_colorbar",
            "colormap",
            "emax",
        ],
    },
    "save_data_block": {
        "required": False,
        "required_fields": ["output_filename"],
        "optional_fields": [],
    },
}

SUPPORTED_SCALES = ["lin", "log"]
SUPPORTED_COLORMAPS = ["viridis", "viridis_r", "hot", "hot_r"]


def _normalize_plot_limit(value, field_label):
    """Return an automatic or numerical plotting limit."""
    if value is None:
        return None
    if isinstance(value, str):
        normalized_value = normalize_string(value=value)
        if normalized_value == "auto":
            return normalized_value
        raise ConfigError(f"{field_label} must be a real number or 'auto'.")
    return require_real(value=value, field_label=field_label)


def read_config(config_filename):
    """Read the Stage 1 namelist into a canonical block-structured dictionary.

    Args:
        config_filename (str or Path): Test 003 Stage 1 namelist file.

    Returns:
        dict: Validated ``input``, ``plot_data_block``, and
        ``save_data_block`` settings.

    Raises:
        ConfigError: If a required setting is absent or invalid.
    """
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    input_block = blocks["input"]
    plot_block = blocks.get("plot_data_block", {})
    save_block = blocks.get("save_data_block", {})

    input_value = require_string(
        value=input_block["input_config"], field_label="input_config"
    )
    input_config = normalize_path(
        value=input_value,
        config_path=config_path,
        field_label="input_config",
    )
    require_existing_file(path=input_config, field_label="input_config")

    plot_data = require_boolean(
        value=input_block.get("plot_data", False), field_label="plot_data"
    )
    save_data = require_boolean(
        value=input_block.get("save_data", False), field_label="save_data"
    )

    # Both HDF5 and PNG outputs use output_filename as their common basename.
    if save_data or plot_data:
        require_blocks(blocks=blocks, required_blocks=["save_data_block"])

    output_filename = None
    if "save_data_block" in blocks:
        output_value = require_string(
            value=save_block["output_filename"], field_label="output_filename"
        )
        output_filename = normalize_path(
            value=output_value,
            config_path=config_path,
            field_label="output_filename",
        )
        if output_filename.suffix.lower() != ".h5":
            raise ConfigError("output_filename must have a .h5 extension.")

    scale = normalize_string(value=plot_block.get("scale", "lin"))
    scale = require_choice(
        value=scale,
        supported_values=SUPPORTED_SCALES,
        field_label="scale",
    )
    colormap = normalize_string(value=plot_block.get("colormap", "viridis"))
    colormap = require_choice(
        value=colormap,
        supported_values=SUPPORTED_COLORMAPS,
        field_label="colormap",
    )
    enable_colorbar = require_boolean(
        value=plot_block.get("enable_colorbar", True),
        field_label="enable_colorbar",
    )
    emax = require_real(value=plot_block.get("emax", 150.0), field_label="emax")
    if emax <= 0.0:
        raise ConfigError("emax must be greater than zero.")

    return {
        "input": {
            "input_config": input_config,
            "plot_data": plot_data,
            "save_data": save_data,
        },
        "plot_data_block": {
            "scale": scale,
            "fmin": _normalize_plot_limit(plot_block.get("fmin"), "fmin"),
            "fmax": _normalize_plot_limit(plot_block.get("fmax"), "fmax"),
            "enable_colorbar": enable_colorbar,
            "colormap": colormap,
            "emax": emax,
        },
        "save_data_block": {"output_filename": output_filename},
    }
