"""Configuration parsing for the reference data generator."""

from pathlib import Path

import f90nml


# These collections provide one place to discover and maintain accepted values.
SUPPORTED_SPECIES = ["h", "d", "t"]
SUPPORTED_INPUT_FILE_TYPES = [
    "fidasim_h5",
    "cql3d_f4d",
]
SUPPORTED_PLOT_DATA_SCALES = ["lin", "log"]
SUPPORTED_PLOT_DATA_COLORMAPS = [
    "viridis",
    "viridis_r",
    "hot",
    "hot_r",
]


def _normalize_string(value):
    """Return a canonical lowercase string while preserving invalid types."""
    if isinstance(value, str):
        return value.strip().lower()
    return value


def _normalize_path(value, config_dir):
    """Return an absolute path resolved from the namelist directory."""
    if not isinstance(value, str):
        return value

    stripped_value = value.strip()
    if not stripped_value:
        return stripped_value

    path = Path(stripped_value).expanduser()
    if not path.is_absolute():
        path = config_dir / path

    return str(path.resolve())


def _normalize_locations(value, label):
    """Represent scalar or array locations as a list of floats."""
    if isinstance(value, list):
        raw_locations = value
    elif value is None:
        raw_locations = []
    else:
        raw_locations = [value]

    locations = []
    for raw_location in raw_locations:
        try:
            location = float(raw_location)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{label} must contain only numeric values.") from exc
        locations.append(location)

    return locations


def _normalize_blocks(blocks, config_dir):
    """Normalize equivalent input representations into canonical values."""
    # Read each namelist block. Plotting and saving blocks are optional.
    input_block = blocks["input"]
    plot_block = blocks.get("plot_data_block", {})
    save_block = blocks.get("save_data_block", {})

    # Normalize the main input settings.
    input_config = {}
    input_config["input_file_type"] = _normalize_string(
        input_block.get("input_file_type")
    )
    input_config["input_filename"] = _normalize_path(
        input_block.get("input_filename"), config_dir
    )
    input_config["species"] = _normalize_string(input_block.get("species"))
    input_config["atomic_number"] = input_block.get("atomic_number")
    input_config["mass_number"] = input_block.get("mass_number")
    input_config["charge_state"] = input_block.get("charge_state")
    input_config["r_locations"] = _normalize_locations(
        input_block.get("r_locations"), "r_locations"
    )
    input_config["z_locations"] = _normalize_locations(
        input_block.get("z_locations"), "z_locations"
    )

    # Preserve logical values unchanged so validation can reject non-booleans.
    input_config["plot_data"] = input_block.get("plot_data", False)
    input_config["save_data"] = input_block.get("save_data", False)

    # Normalize optional plotting settings and provide their default values.
    plot_config = {}
    plot_config["scale"] = _normalize_string(plot_block.get("scale", "lin"))
    plot_config["fmin"] = plot_block.get("fmin")
    plot_config["fmax"] = plot_block.get("fmax")
    plot_config["enable_colorbar"] = plot_block.get("enable_colorbar", True)
    plot_config["colormap"] = _normalize_string(
        plot_block.get("colormap", "viridis")
    )

    # Normalize the optional output base path used when saving is enabled.
    save_config = {}
    save_config["output_filename"] = _normalize_path(
        save_block.get("output_filename"), config_dir
    )

    # Preserve the same three-section organization used by the namelist.
    return {
        "input": input_config,
        "plot_data": plot_config,
        "save_data": save_config,
    }


def _validate_required_fields(config_block, required_fields, block_name):
    """Raise if any required fields are missing from a config block."""
    missing_fields = []
    for field in required_fields:
        if config_block[field] is None:
            missing_fields.append(field)

    if missing_fields:
        missing_text = ", ".join(missing_fields)
        raise ValueError(f"Missing required field(s) in &{block_name}: {missing_text}")


def _validate_supported_value(value, label, supported_values):
    if not isinstance(value, str) or value not in supported_values:
        supported_text = ", ".join(supported_values)
        raise ValueError(
            f"Unsupported {label} '{value}'. Supported values: {supported_text}"
        )


def _validate_boolean(value, label):
    if not isinstance(value, bool):
        raise ValueError(f"{label} must be .true. or .false.")


def _validate_input_config(input_config):
    
    
    # Check that input block has required fields:
    required_fields = (
        "input_file_type",
        "input_filename",
        "species",
        "atomic_number",
        "mass_number",
        "charge_state",
    )
    _validate_required_fields(input_config, required_fields, "input")

    # Check that selector variables have supported values and are not empty:
    _validate_supported_value(
        value=input_config["input_file_type"],
        label="input file type",
        supported_values=SUPPORTED_INPUT_FILE_TYPES,
    )
    _validate_supported_value(
        value=input_config["species"],
        label="species",
        supported_values=SUPPORTED_SPECIES
    )
    _validate_boolean(input_config["plot_data"], "plot_data")
    _validate_boolean(input_config["save_data"], "save_data")

    # Check that input filename is specified:
    input_filename = input_config["input_filename"]
    if not isinstance(input_filename, str) or not input_filename:
        raise ValueError("input_filename must be a non-empty path string.")

    # Check that species values are integers:
    for field in ("atomic_number", "mass_number", "charge_state"):
        value = input_config[field]
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"{field} must be an integer.")

    # Check that species values are physically consistent:
    atomic_number = input_config["atomic_number"]
    mass_number = input_config["mass_number"]
    charge_state = input_config["charge_state"]

    if atomic_number <= 0:
        raise ValueError("atomic_number must be greater than zero.")
    if mass_number < atomic_number:
        raise ValueError("mass_number must be greater than or equal to atomic_number.")
    if charge_state < 0 or charge_state > atomic_number:
        raise ValueError("charge_state must be between zero and atomic_number.")

    # Check that r_locations and z_locations are consistent and non-empty:
    r_locations = input_config["r_locations"]
    z_locations = input_config["z_locations"]
    if len(r_locations) != len(z_locations):
        raise ValueError("r_locations and z_locations must have the same length.")
    if not r_locations:
        raise ValueError("At least one (r, z) location must be provided.")


def _validate_plot_config(plot_config):

    # Check that selector variables have supported values and are not empty:
    _validate_supported_value(
        value=plot_config["scale"],
        label="plot scale",
        supported_values=SUPPORTED_PLOT_DATA_SCALES
    )
    _validate_supported_value(
        value=plot_config["colormap"],
        label="plot colormap",
        supported_values=SUPPORTED_PLOT_DATA_COLORMAPS,
    )
    _validate_boolean(plot_config["enable_colorbar"], "enable_colorbar")


def _validate_save_config(save_config):

    # Check that output_filename is specified:
    output_filename = save_config["output_filename"]
    if not isinstance(output_filename, str) or not output_filename:
        raise ValueError(
            "output_filename must be a non-empty string when save_data=.true."
        )

    output_path = Path(output_filename)
    if output_path.suffix.lower() != ".h5":
        raise ValueError("output_filename must have a .h5 extension.")


def _validate_config(config):
    """Validate canonical values without modifying the configuration."""

    # Validate the input block first, since other blocks depend on it:
    input_config = config["input"]
    _validate_input_config(input_config)

    # Dependent blocks matter only when their corresponding feature is enabled.
    if input_config["plot_data"]:
        _validate_plot_config(config["plot_data"])
    if input_config["save_data"]:
        _validate_save_config(config["save_data"])


def parse_config(config_path):
    """Parse, normalize, validate, and return the canonical configuration."""
    config_path = Path(config_path)

    try:
        blocks = f90nml.read(config_path)
    except (OSError, ValueError) as exc:
        message = f"Unable to read namelist config '{config_path}': {exc}"
        raise ValueError(message) from exc

    if "input" not in blocks:
        raise ValueError("Missing &input block in config.")

    # Get absolute path of input configuration file:
    config_dir = config_path.parent.resolve()

    # Normalize input configuration data:
    normalized_config = _normalize_blocks(blocks, config_dir)

    # Validate normalized configuration data:
    _validate_config(normalized_config)

    return normalized_config
