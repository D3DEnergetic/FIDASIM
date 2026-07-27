"""Read comparison settings and identify reference-sampled file pairs."""

from dataclasses import dataclass
import math
from pathlib import Path
import re

from regression_test_tools import (
    ConfigError,
    normalize_path,
    normalize_string,
    read_namelist,
    require_boolean,
    require_blocks,
    require_choice,
    require_existing_directory,
    require_existing_file,
    require_fields,
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
            "moment_relative_tolerance",
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
    "input_distribution_config",
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
    moment_relative_tolerance = require_real(
        value=compare_block["moment_relative_tolerance"],
        field_label="moment_relative_tolerance",
    )
    if (
        not math.isfinite(moment_relative_tolerance)
        or moment_relative_tolerance <= 0.0
    ):
        raise ConfigError(
            "moment_relative_tolerance must be finite and greater than zero."
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
            "moment_relative_tolerance": moment_relative_tolerance,
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

    # Step 2: discover the native Test 002 reference collection.
    distribution_config_value = require_string(
        value=run_block["input_distribution_config"],
        field_label="input_distribution_config",
    )
    distribution_config = normalize_path(
        value=distribution_config_value,
        config_path=sampling_config_path,
        field_label="input_distribution_config",
    )
    require_existing_file(
        path=distribution_config,
        field_label="input_distribution_config",
    )
    reference_files = _discover_reference_files(distribution_config)

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

    for reference_file in reference_files:
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


def _discover_reference_files(distribution_config):
    """Discover contiguous Test 002 outputs from a Stage 2 configuration."""
    distribution_path, blocks = read_namelist(config_path=distribution_config)
    require_blocks(blocks=blocks, required_blocks=["save_data_block"])
    save_block = blocks["save_data_block"]
    require_fields(
        block=save_block,
        required_fields=["output_filename"],
        block_label="Test 002 &save_data_block",
    )
    output_value = require_string(
        value=save_block["output_filename"],
        field_label="Test 002 output_filename",
    )
    output_base = normalize_path(
        value=output_value,
        config_path=distribution_path,
        field_label="Test 002 output_filename",
    )

    pattern = re.compile(
        rf"{re.escape(output_base.stem)}_(\d{{3}})"
        rf"{re.escape(output_base.suffix)}"
    )
    indexed_paths = []
    for path in output_base.parent.glob(
        f"{output_base.stem}_*{output_base.suffix}"
    ):
        match = pattern.fullmatch(path.name)
        if match is not None and path.is_file():
            indexed_paths.append((int(match.group(1)), path.resolve()))
    indexed_paths.sort(key=lambda item: item[0])

    indices = [index for index, _ in indexed_paths]
    if not indexed_paths:
        raise ConfigError("No Test 002 Stage 2 reference files were found.")
    if indices != list(range(1, len(indices) + 1)):
        raise ConfigError(
            "Test 002 output indices must be contiguous and start at 001."
        )
    return [path for _, path in indexed_paths]
