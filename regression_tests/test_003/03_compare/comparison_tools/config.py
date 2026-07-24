"""Read comparison settings and identify reference-sampled file pairs."""

from dataclasses import dataclass
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
    require_existing_directory,
    require_existing_file,
    require_fields,
    require_integer,
    require_real,
    require_string,
    validate_schema,
)


# Define the complete structure of the Stage 3 namelist in one place.
CONFIG_SCHEMA = {
    "compare": {
        "required": True,
        "required_fields": [
            "sampling_config_file",
            "output_directory",
        ],
        "optional_fields": [
            "comment",
            "generate_plots",
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
            "emax",
        ],
    },
}

SUPPORTED_SCALES = [
    "lin",
    "log",
]
SUPPORTED_COLORMAPS = [
    "viridis",
    "viridis_r",
    "hot",
    "hot_r",
]


# Stage 3 consumes these fields from the Stage 2 configuration. Stage 2 owns
# the complete run_test schema, so Stage 3 does not reject its other fields.
SAMPLING_REQUIRED_BLOCKS = (
    "run_test",
)
SAMPLING_REQUIRED_FIELDS = (
    "n_reference_files",
    "reference_files",
    "output_directory",
)


@dataclass
class FilePair:
    reference: Path
    sampled: Path


def _normalize_plot_limit(value, field_label):
    """Normalize an automatic or numerical plot limit.

    Args:
        value (object): Plot limit read from the namelist. Accepted values are
        ``None``, the string ``"auto"``, or a real number.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        float, str, or None: ``None`` and ``"auto"`` retain their automatic
        meaning. Numerical values are returned as Python ``float`` objects.

    Raises:
        ConfigError: If a string other than ``"auto"`` is supplied or the
        value is not a real number.
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


def read_config(config_filename):
    """Read and validate the Stage 3 comparison configuration.

    Args:
        config_filename (str or Path): Path to the Stage 3 namelist file.

    Returns:
        dict: Canonical configuration organized under ``compare`` and
        ``plot_data_block`` keys that match the namelist block names.

    Raises:
        ConfigError: If a required block or field is missing, a block contains
        an unknown field, or a configured value is invalid.
    """

    # Step 1: read the namelist and validate its complete Stage 3 structure.
    config_path, blocks = read_namelist(
        config_path=config_filename,
    )
    validate_schema(
        blocks=blocks,
        schema=CONFIG_SCHEMA,
    )

    compare_block = blocks["compare"]
    plot_block = blocks.get("plot_data_block", {})

    # Step 2: validate and resolve the required paths.
    sampling_config_value = require_string(
        value=compare_block["sampling_config_file"],
        field_label="sampling_config_file",
    )
    output_directory_value = require_string(
        value=compare_block["output_directory"],
        field_label="output_directory",
    )

    sampling_config_file = normalize_path(
        value=sampling_config_value,
        config_path=config_path,
        field_label="sampling_config_file",
    )
    output_directory = normalize_path(
        value=output_directory_value,
        config_path=config_path,
        field_label="output_directory",
    )
    require_existing_file(
        path=sampling_config_file,
        field_label="sampling_config_file",
    )

    # Step 3: validate the optional comparison and plotting settings.
    generate_plots = require_boolean(
        value=compare_block.get("generate_plots", True),
        field_label="generate_plots",
    )
    enable_colorbar = require_boolean(
        value=plot_block.get("enable_colorbar", True),
        field_label="enable_colorbar",
    )

    scale = normalize_string(
        value=plot_block.get("scale", "lin"),
    )
    scale = require_choice(
        value=scale,
        supported_values=SUPPORTED_SCALES,
        field_label="scale",
    )

    colormap = normalize_string(
        value=plot_block.get("colormap", "viridis"),
    )
    colormap = require_choice(
        value=colormap,
        supported_values=SUPPORTED_COLORMAPS,
        field_label="colormap",
    )

    fmin = _normalize_plot_limit(
        value=plot_block.get("fmin"),
        field_label="fmin",
    )
    fmax = _normalize_plot_limit(
        value=plot_block.get("fmax"),
        field_label="fmax",
    )
    emax = require_real(
        value=plot_block.get("emax", 150.0),
        field_label="emax",
    )
    if emax <= 0.0:
        raise ConfigError("emax must be greater than zero.")

    # Step 4: preserve the namelist block structure in the returned canonical
    # configuration. Optional fields and blocks are populated with defaults.
    return {
        "compare": {
            "sampling_config_file": sampling_config_file,
            "output_directory": output_directory,
            "generate_plots": generate_plots,
        },
        "plot_data_block": {
            "scale": scale,
            "fmin": fmin,
            "fmax": fmax,
            "enable_colorbar": enable_colorbar,
            "colormap": colormap,
            "emax": emax,
        },
    }


def build_file_pairs(sampling_config_file):
    """Construct ordered reference-sampled pairs from a Stage 2 namelist.

    Args:
        sampling_config_file (str or Path): Path to the Stage 2 namelist used
        to generate the sampled distributions.

    Returns:
        list of FilePair: Reference and sampled paths in the order configured
        by Stage 2.

    Raises:
        ConfigError: If the required Stage 2 fields are invalid, directories
        or files are missing, or reference basenames are duplicated.
    """

    # Step 1: read the Stage 2 namelist and require only the fields consumed by
    # Stage 3. Other valid Stage 2 fields are deliberately ignored here.
    sampling_config_path, blocks = read_namelist(
        config_path=sampling_config_file,
    )
    require_blocks(
        blocks=blocks,
        required_blocks=SAMPLING_REQUIRED_BLOCKS,
    )

    run_block = blocks["run_test"]
    require_fields(
        block=run_block,
        required_fields=SAMPLING_REQUIRED_FIELDS,
        block_label="&run_test",
    )

    # Step 2: validate the file count, reference list, and sampled directory.
    number_of_files = require_integer(
        value=run_block["n_reference_files"],
        field_label="n_reference_files",
    )
    if number_of_files < 1:
        raise ConfigError("n_reference_files must be positive.")

    reference_values = as_list(
        value=run_block["reference_files"],
    )
    if len(reference_values) != number_of_files:
        raise ConfigError(
            "n_reference_files does not match the number of configured "
            "reference_files."
        )

    sampled_directory_value = require_string(
        value=run_block["output_directory"],
        field_label="Stage 2 output_directory",
    )
    sampled_directory = normalize_path(
        value=sampled_directory_value,
        config_path=sampling_config_path,
        field_label="Stage 2 output_directory",
    )
    require_existing_directory(
        path=sampled_directory,
        field_label="Stage 2 output_directory",
    )

    # Step 3: preserve the configured order while constructing each file pair.
    pairs = []
    used_basenames = set()

    for reference_index, reference_value in enumerate(reference_values, start=1):
        reference_label = f"reference_files entry {reference_index}"
        reference_value = require_string(
            value=reference_value,
            field_label=reference_label,
        )
        reference_file = normalize_path(
            value=reference_value,
            config_path=sampling_config_path,
            field_label=reference_label,
        )
        require_existing_file(
            path=reference_file,
            field_label=reference_label,
        )

        basename = reference_file.name
        if basename in used_basenames:
            raise ConfigError(f"Duplicate reference basename: {basename}")
        used_basenames.add(basename)

        sampled_file = sampled_directory / basename
        require_existing_file(
            path=sampled_file,
            field_label=f"sampled file for {basename}",
        )

        pair = FilePair(
            reference=reference_file,
            sampled=sampled_file,
        )
        pairs.append(pair)

    return pairs
