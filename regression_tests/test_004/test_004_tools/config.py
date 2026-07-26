"""Read and validate the unified Test 004 configuration."""

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


SHARED_CONFIG_SCHEMA = {
    "test_case": {
        "required": True,
        "required_fields": [
            "input_distribution_config",
            "tables_filename",
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

DETERMINISTIC_FIELDS = ["n_gyro", "plot_data", "save_data"]
MONTE_CARLO_FIELDS = [
    "n_markers",
    "reservoir_size",
    "seed",
    "plot_data",
    "save_data",
]

DETERMINISTIC_BLOCK_SCHEMA = {
    "required": True,
    "required_fields": DETERMINISTIC_FIELDS,
    "optional_fields": ["comment"],
}
MONTE_CARLO_BLOCK_SCHEMA = {
    "required": True,
    "required_fields": MONTE_CARLO_FIELDS,
    "optional_fields": ["comment"],
}


def _optional_schema(block_schema):
    """Allow an unselected implementation block without requiring its fields."""
    return {
        "required": False,
        "required_fields": [],
        "optional_fields": (
            block_schema["required_fields"] + block_schema["optional_fields"]
        ),
    }


DETERMINISTIC_CONFIG_SCHEMA = {
    **SHARED_CONFIG_SCHEMA,
    "deterministic": DETERMINISTIC_BLOCK_SCHEMA,
    "monte_carlo": _optional_schema(MONTE_CARLO_BLOCK_SCHEMA),
}
MONTE_CARLO_CONFIG_SCHEMA = {
    **SHARED_CONFIG_SCHEMA,
    "deterministic": _optional_schema(DETERMINISTIC_BLOCK_SCHEMA),
    "monte_carlo": MONTE_CARLO_BLOCK_SCHEMA,
}
FULL_CONFIG_SCHEMA = {
    **SHARED_CONFIG_SCHEMA,
    "deterministic": DETERMINISTIC_BLOCK_SCHEMA,
    "monte_carlo": MONTE_CARLO_BLOCK_SCHEMA,
}

SUPPORTED_LEVEL_SPLITS = ["ground-only", "exponential"]
SUPPORTED_SCALES = ["lin", "log"]
SUPPORTED_COLORMAPS = ["viridis", "viridis_r", "hot", "hot_r"]
INT32_MAX = 2**31 - 1
INT64_MAX = 2**63 - 1


def _optional_comment(block, field_label):
    if "comment" not in block:
        return ""
    return require_string(block["comment"], field_label)


def _positive_integer(block, field_name, maximum=None):
    raw_value = block[field_name]
    if isinstance(raw_value, bool) or not isinstance(raw_value, (int, float)):
        raise ConfigError(f"{field_name} must be an integer value.")
    if isinstance(raw_value, float):
        if not math.isfinite(raw_value) or not raw_value.is_integer():
            raise ConfigError(f"{field_name} must be an integer value.")
    value = int(raw_value)
    if value < 1:
        raise ConfigError(f"{field_name} must be positive.")
    if maximum is not None and value > maximum:
        raise ConfigError(f"{field_name} must not exceed {maximum}.")
    return value


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


def _read_test_case(block, config_path):
    distribution_config = normalize_path(
        value=require_string(
            block["input_distribution_config"],
            "input_distribution_config",
        ),
        config_path=config_path,
        field_label="input_distribution_config",
    )
    require_existing_file(distribution_config, "input_distribution_config")

    tables_filename = normalize_path(
        value=require_string(block["tables_filename"], "tables_filename"),
        config_path=config_path,
        field_label="tables_filename",
    )
    require_existing_file(tables_filename, "tables_filename")

    return {
        "comment": _optional_comment(block, "test_case comment"),
        "config_path": config_path,
        "input_distribution_config": distribution_config,
        "tables_filename": tables_filename,
    }


def _read_neutrals(block):
    density = require_real(block["density"], "density")
    energy = require_real(block["energy"], "energy")
    injection_angle = require_real(
        block["injection_angle"],
        "injection_angle",
    )
    if not math.isfinite(density) or density <= 0.0:
        raise ConfigError("density must be finite and greater than zero.")
    if not math.isfinite(energy) or energy <= 0.0:
        raise ConfigError("energy must be finite and greater than zero.")
    if not math.isfinite(injection_angle):
        raise ConfigError("injection_angle must be finite.")

    level_split_method = require_choice(
        normalize_string(block["level_split_method"]),
        SUPPORTED_LEVEL_SPLITS,
        "level_split_method",
    )
    level_decay = None
    if level_split_method == "exponential":
        if "level_decay" not in block:
            raise ConfigError(
                "level_decay is required for exponential level splitting."
            )
        level_decay = require_real(block["level_decay"], "level_decay")
        if not math.isfinite(level_decay) or level_decay <= 0.0:
            raise ConfigError(
                "level_decay must be finite and greater than zero for "
                "exponential splitting."
            )
    elif "level_decay" in block:
        inactive_decay = require_real(block["level_decay"], "level_decay")
        if not math.isfinite(inactive_decay):
            raise ConfigError("level_decay must be finite.")

    return {
        "density": density,
        "energy": energy,
        "injection_angle": injection_angle,
        "level_split_method": level_split_method,
        "level_decay": level_decay,
    }


def _read_deterministic(block):
    return {
        "comment": _optional_comment(block, "deterministic comment"),
        "n_gyro": _positive_integer(block, "n_gyro"),
        "plot_data": require_boolean(block["plot_data"], "plot_data"),
        "save_data": require_boolean(block["save_data"], "save_data"),
    }


def _read_monte_carlo(block):
    plot_data = require_boolean(block["plot_data"], "plot_data")
    save_data = require_boolean(block["save_data"], "save_data")
    if plot_data and not save_data:
        raise ConfigError(
            "Monte Carlo plot_data requires save_data because the plot is "
            "generated from the sink HDF5 output."
        )

    return {
        "comment": _optional_comment(block, "monte_carlo comment"),
        "n_markers": _positive_integer(
            block,
            "n_markers",
            maximum=INT64_MAX,
        ),
        "reservoir_size": _positive_integer(
            block,
            "reservoir_size",
            maximum=INT32_MAX,
        ),
        "seed": _positive_integer(block, "seed", maximum=INT32_MAX),
        "plot_data": plot_data,
        "save_data": save_data,
    }


def _read_plot_data(block):
    scale = require_choice(
        normalize_string(block.get("scale", "lin")),
        SUPPORTED_SCALES,
        "scale",
    )
    colormap = require_choice(
        normalize_string(block.get("colormap", "viridis")),
        SUPPORTED_COLORMAPS,
        "colormap",
    )
    enable_colorbar = require_boolean(
        block.get("enable_colorbar", True),
        "enable_colorbar",
    )
    contour_levels = require_integer(
        block.get("contour_levels", 100),
        "contour_levels",
    )
    if contour_levels < 2:
        raise ConfigError("contour_levels must be at least 2.")

    fmin = _normalize_plot_limit(block.get("fmin"), "fmin")
    fmax = _normalize_plot_limit(block.get("fmax"), "fmax")
    emax = _normalize_plot_limit(block.get("emax"), "emax")
    if isinstance(fmin, float) and isinstance(fmax, float) and fmin >= fmax:
        raise ConfigError("fmin must be less than fmax.")
    if scale == "log":
        if isinstance(fmin, float) and fmin <= 0.0:
            raise ConfigError("fmin must be greater than zero for log scale.")
        if isinstance(fmax, float) and fmax <= 0.0:
            raise ConfigError("fmax must be greater than zero for log scale.")
    if isinstance(emax, float) and emax <= 0.0:
        raise ConfigError("emax must be greater than zero.")

    return {
        "scale": scale,
        "emax": emax,
        "fmin": fmin,
        "fmax": fmax,
        "enable_colorbar": enable_colorbar,
        "colormap": colormap,
        "contour_levels": contour_levels,
    }


def _read_blocks(config_filename, schema):
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=schema)
    return config_path, blocks


def _read_shared(config_path, blocks, plot_enabled, output_enabled):
    if plot_enabled:
        require_blocks(blocks, ["plot_data_block"])
    if output_enabled:
        require_blocks(blocks, ["save_data_block"])

    plot_data = _read_plot_data(blocks.get("plot_data_block", {}))

    output_directory = None
    if "save_data_block" in blocks:
        output_directory = normalize_path(
            value=require_string(
                blocks["save_data_block"]["output_directory"],
                "output_directory",
            ),
            config_path=config_path,
            field_label="output_directory",
        )

    shared = {
        "test_case": _read_test_case(blocks["test_case"], config_path),
        "neutrals": _read_neutrals(blocks["neutrals"]),
        "plot_data_block": plot_data,
        "save_data_block": {
            "output_directory": output_directory,
            "deterministic_directory": (
                None
                if output_directory is None
                else output_directory / "deterministic"
            ),
            "monte_carlo_directory": (
                None
                if output_directory is None
                else output_directory / "monte_carlo"
            ),
        },
    }
    return shared


def read_deterministic_config(config_filename):
    """Return shared settings and the validated deterministic block."""
    config_path, blocks = _read_blocks(
        config_filename,
        schema=DETERMINISTIC_CONFIG_SCHEMA,
    )
    deterministic = _read_deterministic(blocks["deterministic"])
    config = _read_shared(
        config_path=config_path,
        blocks=blocks,
        plot_enabled=deterministic["plot_data"],
        output_enabled=(
            deterministic["plot_data"] or deterministic["save_data"]
        ),
    )
    config["deterministic"] = deterministic
    return config


def read_monte_carlo_config(config_filename):
    """Return shared settings and the validated Monte Carlo block."""
    config_path, blocks = _read_blocks(
        config_filename,
        schema=MONTE_CARLO_CONFIG_SCHEMA,
    )
    monte_carlo = _read_monte_carlo(blocks["monte_carlo"])
    config = _read_shared(
        config_path=config_path,
        blocks=blocks,
        plot_enabled=monte_carlo["plot_data"],
        output_enabled=monte_carlo["plot_data"] or monte_carlo["save_data"],
    )
    config["monte_carlo"] = monte_carlo
    return config


def read_test_config(config_filename):
    """Return the fully validated configuration for one Test 004 case."""
    config_path, blocks = _read_blocks(
        config_filename,
        schema=FULL_CONFIG_SCHEMA,
    )
    deterministic = _read_deterministic(blocks["deterministic"])
    monte_carlo = _read_monte_carlo(blocks["monte_carlo"])
    config = _read_shared(
        config_path=config_path,
        blocks=blocks,
        plot_enabled=(
            deterministic["plot_data"] or monte_carlo["plot_data"]
        ),
        output_enabled=(
            deterministic["plot_data"]
            or deterministic["save_data"]
            or monte_carlo["plot_data"]
            or monte_carlo["save_data"]
        ),
    )
    config["deterministic"] = deterministic
    config["monte_carlo"] = monte_carlo
    return config
