"""Read and validate the Test 002 Stage 2 configuration."""

from regression_test_tools import (
    ConfigError,
    normalize_path,
    normalize_string,
    read_namelist,
    require_boolean,
    require_choice,
    require_existing_file,
    require_integer,
    require_real,
    require_string,
    validate_schema,
)


CONFIG_SCHEMA = {
    "input": {
        "required": True,
        "required_fields": ["reference_config", "nenergy", "npitch", "plot_data"],
        "optional_fields": ["comment"],
    },
    "plot_data_block": {
        "required": False,
        "required_fields": [
            "scale",
            "fmin",
            "fmax",
            "enable_colorbar",
            "colormap",
            "contour_levels",
            "emax",
        ],
        "optional_fields": [],
    },
    "save_data_block": {
        "required": True,
        "required_fields": ["output_filename"],
        "optional_fields": [],
    },
}

SUPPORTED_SCALES = ["lin", "log"]
SUPPORTED_COLORMAPS = ["viridis", "viridis_r", "hot", "hot_r"]


def _normalize_plot_limit(value, field_label):
    """Return either ``auto`` or a numerical plotting limit."""
    if isinstance(value, str):
        normalized_value = normalize_string(value=value)
        if normalized_value == "auto":
            return normalized_value
        raise ConfigError(f"{field_label} must be a real number or 'auto'.")

    return require_real(value=value, field_label=field_label)


def read_config(config_filename):
    """Read the Stage 2 namelist into a canonical block-structured dictionary."""
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    input_block = blocks["input"]
    plot_block = blocks.get("plot_data_block")
    save_block = blocks["save_data_block"]

    reference_config_value = require_string(
        value=input_block["reference_config"],
        field_label="reference_config",
    )
    reference_config = normalize_path(
        value=reference_config_value,
        config_path=config_path,
        field_label="reference_config",
    )
    require_existing_file(path=reference_config, field_label="reference_config")

    nenergy = require_integer(value=input_block["nenergy"], field_label="nenergy")
    npitch = require_integer(value=input_block["npitch"], field_label="npitch")
    if nenergy < 2 or npitch < 2:
        raise ConfigError("nenergy and npitch must both be at least 2.")

    plot_data = require_boolean(
        value=input_block["plot_data"], field_label="plot_data"
    )
    if plot_data and plot_block is None:
        raise ConfigError(
            "plot_data_block is required when input/plot_data is true."
        )

    plot_config = None
    if plot_block is not None:
        scale = normalize_string(value=plot_block["scale"])
        scale = require_choice(
            value=scale,
            supported_values=SUPPORTED_SCALES,
            field_label="scale",
        )
        colormap = normalize_string(value=plot_block["colormap"])
        colormap = require_choice(
            value=colormap,
            supported_values=SUPPORTED_COLORMAPS,
            field_label="colormap",
        )
        enable_colorbar = require_boolean(
            value=plot_block["enable_colorbar"],
            field_label="enable_colorbar",
        )
        contour_levels = require_integer(
            value=plot_block["contour_levels"],
            field_label="contour_levels",
        )
        if contour_levels < 2:
            raise ConfigError("contour_levels must be at least 2.")

        emax = require_real(value=plot_block["emax"], field_label="emax")
        if emax <= 0.0:
            raise ConfigError("emax must be greater than zero.")

        plot_config = {
            "scale": scale,
            "fmin": _normalize_plot_limit(plot_block["fmin"], "fmin"),
            "fmax": _normalize_plot_limit(plot_block["fmax"], "fmax"),
            "enable_colorbar": enable_colorbar,
            "colormap": colormap,
            "contour_levels": contour_levels,
            "emax": emax,
        }

    output_value = require_string(
        value=save_block["output_filename"],
        field_label="output_filename",
    )
    output_filename = normalize_path(
        value=output_value,
        config_path=config_path,
        field_label="output_filename",
    )
    if output_filename.suffix.lower() not in [".h5", ".hdf5"]:
        raise ConfigError("output_filename must have a .h5 or .hdf5 extension.")

    return {
        "input": {
            "reference_config": reference_config,
            "nenergy": nenergy,
            "npitch": npitch,
            "plot_data": plot_data,
        },
        "plot_data_block": plot_config,
        "save_data_block": {"output_filename": output_filename},
    }
