"""Read and validate the Test 002 Stage 1 configuration."""

from regression_test_tools import (
    ConfigError,
    as_list,
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
        "required_fields": [
            "input_file_type",
            "input_filename",
            "species",
            "atomic_number",
            "mass_number",
            "charge_state",
            "r_locations",
            "z_locations",
            "save_data",
        ],
        "optional_fields": ["comment", "plot_data"],
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
            "contour_levels",
        ],
    },
    "save_data_block": {
        "required": True,
        "required_fields": ["output_filename"],
        "optional_fields": [],
    },
}

SUPPORTED_INPUT_FILE_TYPES = ["cql3d_f4d"]
SUPPORTED_SPECIES = ["h", "d", "t"]
SUPPORTED_SCALES = ["lin", "log"]
SUPPORTED_COLORMAPS = ["viridis", "viridis_r", "hot", "hot_r"]


def _normalize_locations(value, field_label):
    """Return a scalar or sequence of configured locations as floats."""
    locations = []
    for index, raw_value in enumerate(as_list(value=value), start=1):
        location = require_real(
            value=raw_value,
            field_label=f"{field_label} entry {index}",
        )
        locations.append(location)
    return locations


def _normalize_plot_limit(value, field_label):
    """Normalize an automatic or numerical plotting limit."""
    if value is None:
        return None
    if isinstance(value, str):
        normalized_value = normalize_string(value=value)
        if normalized_value == "auto":
            return normalized_value
        raise ConfigError(f"{field_label} must be a real number or 'auto'.")
    return require_real(value=value, field_label=field_label)


def _validate_particle_properties(input_config):
    """Validate relationships between particle identifiers."""
    atomic_number = input_config["atomic_number"]
    mass_number = input_config["mass_number"]
    charge_state = input_config["charge_state"]

    if atomic_number <= 0:
        raise ConfigError("atomic_number must be greater than zero.")
    if mass_number < atomic_number:
        raise ConfigError(
            "mass_number must be greater than or equal to atomic_number."
        )
    if charge_state < 0 or charge_state > atomic_number:
        raise ConfigError("charge_state must be between zero and atomic_number.")


def read_config(config_filename):
    """Read the canonical Stage 1 configuration.

    Args:
        config_filename (str or Path): Stage 1 Fortran namelist file.

    Returns:
        dict: Canonical configuration organized under ``input``,
        ``plot_data_block``, and ``save_data_block``.

    Raises:
        ConfigError: If the configuration structure or a value is invalid.
    """
    config_path, blocks = read_namelist(config_path=config_filename)
    validate_schema(blocks=blocks, schema=CONFIG_SCHEMA)

    input_block = blocks["input"]
    plot_block = blocks.get("plot_data_block", {})
    save_block = blocks.get("save_data_block", {})

    input_file_type = normalize_string(value=input_block["input_file_type"])
    input_file_type = require_choice(
        value=input_file_type,
        supported_values=SUPPORTED_INPUT_FILE_TYPES,
        field_label="input_file_type",
    )

    input_filename_value = require_string(
        value=input_block["input_filename"],
        field_label="input_filename",
    )
    input_filename = normalize_path(
        value=input_filename_value,
        config_path=config_path,
        field_label="input_filename",
    )
    require_existing_file(path=input_filename, field_label="input_filename")

    species = normalize_string(value=input_block["species"])
    species = require_choice(
        value=species,
        supported_values=SUPPORTED_SPECIES,
        field_label="species",
    )

    atomic_number = require_integer(
        value=input_block["atomic_number"], field_label="atomic_number"
    )
    mass_number = require_integer(
        value=input_block["mass_number"], field_label="mass_number"
    )
    charge_state = require_integer(
        value=input_block["charge_state"], field_label="charge_state"
    )
    r_locations = _normalize_locations(
        value=input_block["r_locations"], field_label="r_locations"
    )
    z_locations = _normalize_locations(
        value=input_block["z_locations"], field_label="z_locations"
    )
    if not r_locations or len(r_locations) != len(z_locations):
        raise ConfigError(
            "r_locations and z_locations must be nonempty and have equal length."
        )

    plot_data = require_boolean(
        value=input_block.get("plot_data", False), field_label="plot_data"
    )
    save_data = require_boolean(
        value=input_block["save_data"], field_label="save_data"
    )
    if not save_data:
        raise ConfigError("save_data must be true for the Stage 1 workflow.")

    scale = normalize_string(value=plot_block.get("scale", "lin"))
    scale = require_choice(
        value=scale, supported_values=SUPPORTED_SCALES, field_label="scale"
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
    fmin = _normalize_plot_limit(value=plot_block.get("fmin"), field_label="fmin")
    fmax = _normalize_plot_limit(value=plot_block.get("fmax"), field_label="fmax")
    contour_levels = require_integer(
        value=plot_block.get("contour_levels", 100),
        field_label="contour_levels",
    )
    if contour_levels < 2:
        raise ConfigError("contour_levels must be at least 2.")

    output_filename_value = require_string(
        value=save_block["output_filename"], field_label="output_filename"
    )
    output_filename = normalize_path(
        value=output_filename_value,
        config_path=config_path,
        field_label="output_filename",
    )
    if output_filename.suffix.lower() not in [".h5", ".hdf5"]:
        raise ConfigError("output_filename must have a .h5 or .hdf5 extension.")

    config = {
        "input": {
            "input_file_type": input_file_type,
            "input_filename": input_filename,
            "species": species,
            "atomic_number": atomic_number,
            "mass_number": mass_number,
            "charge_state": charge_state,
            "r_locations": r_locations,
            "z_locations": z_locations,
            "plot_data": plot_data,
            "save_data": save_data,
        },
        "plot_data_block": {
            "scale": scale,
            "fmin": fmin,
            "fmax": fmax,
            "enable_colorbar": enable_colorbar,
            "colormap": colormap,
            "contour_levels": contour_levels,
        },
        "save_data_block": {"output_filename": output_filename},
    }
    _validate_particle_properties(input_config=config["input"])
    return config
