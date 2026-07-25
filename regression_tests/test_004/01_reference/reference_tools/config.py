"""Read and validate the Test 004 Stage 1 configuration."""

import math

from regression_test_tools import (
    ConfigError,
    normalize_path,
    normalize_string,
    read_namelist,
    require_boolean,
    require_blocks,
    require_choice,
    require_existing_file,
    require_integer,
    require_real,
    require_string,
    validate_schema,
)


CONFIG_SCHEMA = {
    "reference": {
        "required": True,
        "required_fields": [
            "input_distribution_config",
            "tables_filename",
            "n_gyro",
            "plot_data",
            "save_data",
        ],
        "optional_fields": ["comment"],
    },
    "neutrals": {
        "required": True,
        "required_fields": [
            "density",
            "energy",
            "injection_angle",
            "level_split_method",
        ],
        "optional_fields": ["level_decay"],
    },
    "plot_data_block": {
        "required": False,
        "required_fields": [],
        "optional_fields": [
            "scale",
            "emax",
            "fmin",
            "fmax",
            "enable_colorbar",
            "colormap",
            "contour_levels",
        ],
    },
    "save_data_block": {
        "required": False,
        "required_fields": ["output_directory"],
        "optional_fields": [],
    },
}

SUPPORTED_LEVEL_SPLITS = ["ground-only", "exponential"]
SUPPORTED_SCALES = ["lin", "log"]
SUPPORTED_COLORMAPS = ["viridis", "viridis_r", "hot", "hot_r"]


def _normalize_plot_limit(value, field_label):
    if value is None:
        return None
    if isinstance(value, str):
        normalized = normalize_string(value=value)
        if normalized == "auto":
            return normalized
        raise ConfigError(f"{field_label} must be a real number or 'auto'.")
    limit = require_real(value=value, field_label=field_label)
    if not math.isfinite(limit):
        raise ConfigError(f"{field_label} must be finite.")
    return limit


def read_config(config_filename):
    """Return a validated, canonical Stage 1 configuration."""
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    reference_block = blocks["reference"]
    neutral_block = blocks["neutrals"]
    plot_block = blocks.get("plot_data_block", {})
    save_block = blocks.get("save_data_block", {})

    comment = ""
    if "comment" in reference_block:
        comment = require_string(reference_block["comment"], "comment")

    distribution_config = normalize_path(
        value=require_string(
            reference_block["input_distribution_config"],
            "input_distribution_config",
        ),
        config_path=config_path,
        field_label="input_distribution_config",
    )
    require_existing_file(distribution_config, "input_distribution_config")

    tables_filename = normalize_path(
        value=require_string(
            reference_block["tables_filename"], "tables_filename"
        ),
        config_path=config_path,
        field_label="tables_filename",
    )
    require_existing_file(tables_filename, "tables_filename")

    n_gyro = require_integer(reference_block["n_gyro"], "n_gyro")
    if n_gyro < 1:
        raise ConfigError("n_gyro must be positive.")

    plot_data = require_boolean(reference_block["plot_data"], "plot_data")
    save_data = require_boolean(reference_block["save_data"], "save_data")
    if plot_data:
        require_blocks(blocks, ["plot_data_block"])
    if plot_data or save_data:
        require_blocks(blocks, ["save_data_block"])

    density = require_real(neutral_block["density"], "density")
    energy = require_real(neutral_block["energy"], "energy")
    injection_angle = require_real(
        neutral_block["injection_angle"], "injection_angle"
    )
    if not math.isfinite(density) or density <= 0.0:
        raise ConfigError("density must be finite and greater than zero.")
    if not math.isfinite(energy) or energy <= 0.0:
        raise ConfigError("energy must be finite and greater than zero.")
    if not math.isfinite(injection_angle):
        raise ConfigError("injection_angle must be finite.")

    level_split_method = require_choice(
        normalize_string(neutral_block["level_split_method"]),
        SUPPORTED_LEVEL_SPLITS,
        "level_split_method",
    )
    level_decay = None
    if level_split_method == "exponential":
        if "level_decay" not in neutral_block:
            raise ConfigError(
                "level_decay is required for exponential level splitting."
            )
        level_decay = require_real(neutral_block["level_decay"], "level_decay")
        if not math.isfinite(level_decay) or level_decay <= 0.0:
            raise ConfigError(
                "level_decay must be finite and greater than zero for "
                "exponential splitting."
            )
    elif "level_decay" in neutral_block:
        # Validate the supplied but inactive value without assigning it
        # physical meaning in the canonical configuration.
        inactive_decay = require_real(neutral_block["level_decay"], "level_decay")
        if not math.isfinite(inactive_decay):
            raise ConfigError("level_decay must be finite.")

    scale = require_choice(
        normalize_string(plot_block.get("scale", "lin")),
        SUPPORTED_SCALES,
        "scale",
    )
    colormap = require_choice(
        normalize_string(plot_block.get("colormap", "viridis")),
        SUPPORTED_COLORMAPS,
        "colormap",
    )
    enable_colorbar = require_boolean(
        plot_block.get("enable_colorbar", True), "enable_colorbar"
    )
    contour_levels = require_integer(
        plot_block.get("contour_levels", 100), "contour_levels"
    )
    if contour_levels < 2:
        raise ConfigError("contour_levels must be at least 2.")

    fmin = _normalize_plot_limit(plot_block.get("fmin"), "fmin")
    fmax = _normalize_plot_limit(plot_block.get("fmax"), "fmax")
    emax = _normalize_plot_limit(plot_block.get("emax"), "emax")
    if isinstance(fmin, float) and isinstance(fmax, float) and fmin >= fmax:
        raise ConfigError("fmin must be less than fmax.")
    if isinstance(emax, float) and emax <= 0.0:
        raise ConfigError("emax must be greater than zero.")

    output_directory = None
    if save_block:
        output_directory = normalize_path(
            value=require_string(
                save_block["output_directory"], "output_directory"
            ),
            config_path=config_path,
            field_label="output_directory",
        )

    return {
        "reference": {
            "comment": comment,
            "input_distribution_config": distribution_config,
            "tables_filename": tables_filename,
            "n_gyro": n_gyro,
            "plot_data": plot_data,
            "save_data": save_data,
        },
        "neutrals": {
            "density": density,
            "energy": energy,
            "injection_angle": injection_angle,
            "level_split_method": level_split_method,
            "level_decay": level_decay,
        },
        "plot_data_block": {
            "scale": scale,
            "emax": emax,
            "fmin": fmin,
            "fmax": fmax,
            "enable_colorbar": enable_colorbar,
            "colormap": colormap,
            "contour_levels": contour_levels,
        },
        "save_data_block": {"output_directory": output_directory},
    }
