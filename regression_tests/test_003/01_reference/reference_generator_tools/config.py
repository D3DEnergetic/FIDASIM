"""Configuration parsing for the reference data generator."""

from pathlib import Path

from regression_test_tools import (
    ConfigError,
    as_list,
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
    print_config
)


# Define the complete structure of the Stage 1 namelist in one place.
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
        ],
        "optional_fields": [
            "plot_data",
            "save_data",
        ],
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
        ],
    },
    "save_data_block": {
        "required": False,
        "required_fields": [
            "output_filename",
        ],
        "optional_fields": [],
    },
}


# These collections provide one place to discover and maintain accepted values.
SUPPORTED_SPECIES = [
    "h",
    "d",
    "t",
]
SUPPORTED_INPUT_FILE_TYPES = [
    "fidasim_h5",
    "cql3d_f4d",
]
SUPPORTED_PLOT_DATA_SCALES = [
    "lin",
    "log",
]
SUPPORTED_PLOT_DATA_COLORMAPS = [
    "viridis",
    "viridis_r",
    "hot",
    "hot_r",
]


def _normalize_locations(value, field_label):
    """Represent scalar or list locations as a list of floats.

    Args:
        value (object): Scalar or list value read from the namelist.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        list of float: Canonical spatial locations.

    Raises:
        ConfigError: If a location is not a real number.
    """
    raw_locations = as_list(
        value=value,
    )

    locations = []
    for location_index, raw_location in enumerate(raw_locations, start=1):
        location = require_real(
            value=raw_location,
            field_label=f"{field_label} entry {location_index}",
        )
        locations.append(location)

    return locations


def _normalize_plot_limit(value, field_label):
    """Normalize an automatic or numerical plot limit.

    Args:
        value (object): Plot limit read from the namelist. Accepted values are
        ``None``, the string ``"auto"``, or a real number.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        float, str, or None: Canonical plot limit.

    Raises:
        ConfigError: If the value is not numerical, ``None``, or ``"auto"``.
    """
    if value is None:
        return None

    if isinstance(value, str):
        normalized_value = normalize_string(
            value=value,
        )
        if normalized_value == "auto":
            return normalized_value
        raise ConfigError(f"{field_label} must be a real number or 'auto'.")

    return require_real(
        value=value,
        field_label=field_label,
    )


def _validate_particle_properties(input_config):
    """Validate relationships between the particle identifiers."""
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


def _validate_locations(input_config):
    """Require matching, nonempty radial and axial location lists."""
    r_locations = input_config["r_locations"]
    z_locations = input_config["z_locations"]

    if len(r_locations) != len(z_locations):
        raise ConfigError("r_locations and z_locations must have the same length.")
    if not r_locations:
        raise ConfigError("At least one (r, z) location must be provided.")


def read_config(config_filename):
    """Read and validate the Stage 1 reference-data configuration.

    Args:
        config_filename (str or Path): Path to the Stage 1 namelist file.

    Returns:
        dict: Canonical configuration organized under ``input``,
        ``plot_data_block``, and ``save_data_block`` keys that match the
        namelist block names.

    Raises:
        ConfigError: If the namelist structure or a configured value is
        invalid.
    """

    # Step 1: read the namelist and validate its Stage 1 block structure.
    config_path, blocks = read_namelist(
        config_path=config_filename,
    )
    validate_schema(
        blocks=blocks,
        schema=CONFIG_SCHEMA,
    )

    input_block = blocks["input"]
    plot_block = blocks.get("plot_data_block", {})
    save_block = blocks.get("save_data_block", {})

    # Step 2: normalize and validate the primary input settings.
    input_file_type = normalize_string(
        value=input_block["input_file_type"],
    )
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
    require_existing_file(
        path=input_filename,
        field_label="input_filename",
    )

    species = normalize_string(
        value=input_block["species"],
    )
    species = require_choice(
        value=species,
        supported_values=SUPPORTED_SPECIES,
        field_label="species",
    )

    atomic_number = require_integer(
        value=input_block["atomic_number"],
        field_label="atomic_number",
    )
    mass_number = require_integer(
        value=input_block["mass_number"],
        field_label="mass_number",
    )
    charge_state = require_integer(
        value=input_block["charge_state"],
        field_label="charge_state",
    )

    r_locations = _normalize_locations(
        value=input_block["r_locations"],
        field_label="r_locations",
    )
    z_locations = _normalize_locations(
        value=input_block["z_locations"],
        field_label="z_locations",
    )

    plot_data = require_boolean(
        value=input_block.get("plot_data", False),
        field_label="plot_data",
    )
    save_data = require_boolean(
        value=input_block.get("save_data", False),
        field_label="save_data",
    )

    # Step 3: require the save block only when HDF5 output is enabled. If the
    # optional block is present, validate_schema has already required its
    # output_filename field.
    if save_data:
        require_blocks(
            blocks=blocks,
            required_blocks=["save_data_block"],
        )

    # Step 4: normalize and validate the optional plotting settings.
    scale = normalize_string(
        value=plot_block.get("scale", "lin"),
    )
    scale = require_choice(
        value=scale,
        supported_values=SUPPORTED_PLOT_DATA_SCALES,
        field_label="scale",
    )

    colormap = normalize_string(
        value=plot_block.get("colormap", "viridis"),
    )
    colormap = require_choice(
        value=colormap,
        supported_values=SUPPORTED_PLOT_DATA_COLORMAPS,
        field_label="colormap",
    )

    enable_colorbar = require_boolean(
        value=plot_block.get("enable_colorbar", True),
        field_label="enable_colorbar",
    )
    fmin = _normalize_plot_limit(
        value=plot_block.get("fmin"),
        field_label="fmin",
    )
    fmax = _normalize_plot_limit(
        value=plot_block.get("fmax"),
        field_label="fmax",
    )

    # Step 5: normalize the save path when the optional block is present.
    output_filename = None
    if "save_data_block" in blocks:
        output_filename_value = require_string(
            value=save_block["output_filename"],
            field_label="output_filename",
        )
        output_filename = normalize_path(
            value=output_filename_value,
            config_path=config_path,
            field_label="output_filename",
        )
        if output_filename.suffix.lower() != ".h5":
            raise ConfigError("output_filename must have a .h5 extension.")

    # Step 6: assemble the canonical dictionary before applying relationships
    # between fields.
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
        },
        "save_data_block": {
            "output_filename": output_filename,
        },
    }

    _validate_particle_properties(
        input_config=config["input"],
    )
    _validate_locations(
        input_config=config["input"],
    )

    return config
