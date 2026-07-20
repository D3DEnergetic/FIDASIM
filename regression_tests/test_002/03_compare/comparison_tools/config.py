"""Read comparison settings and identify reference-sampled file pairs."""

from dataclasses import dataclass
from pathlib import Path

import f90nml


SUPPORTED_SCALES = {"lin", "log"}
SUPPORTED_COLORMAPS = {"viridis", "viridis_r", "hot", "hot_r"}


@dataclass
class ComparisonConfig:
    sampling_config_file: Path
    output_directory: Path
    generate_plots: bool
    scale: str
    fmin: object
    fmax: object
    enable_colorbar: bool
    colormap: str


@dataclass
class FilePair:
    reference: Path
    sampled: Path


def resolve_path(config_file, configured_path):
    """Resolve a configured path relative to its namelist file."""
    path = Path(configured_path).expanduser()
    if not path.is_absolute():
        path = config_file.parent / path
    return path.resolve()


def read_config(config_filename):
    """Read and validate the comparison and plotting settings."""
    config_file = Path(config_filename).resolve()
    blocks = f90nml.read(config_file)
    if "compare" not in blocks:
        raise ValueError(f"Missing &compare block in {config_file}.")

    compare_block = blocks["compare"]
    plot_block = blocks.get("plot_data_block", {})

    sampling_config_value = compare_block.get("sampling_config_file")
    output_directory_value = compare_block.get("output_directory")
    if not isinstance(sampling_config_value, str) or not sampling_config_value.strip():
        raise ValueError("sampling_config_file must be a non-empty path.")
    if not isinstance(output_directory_value, str) or not output_directory_value.strip():
        raise ValueError("output_directory must be a non-empty path.")

    generate_plots = compare_block.get("generate_plots", True)
    enable_colorbar = plot_block.get("enable_colorbar", True)
    if not isinstance(generate_plots, bool):
        raise ValueError("generate_plots must be .true. or .false.")
    if not isinstance(enable_colorbar, bool):
        raise ValueError("enable_colorbar must be .true. or .false.")

    scale = str(plot_block.get("scale", "lin")).strip().lower()
    colormap = str(plot_block.get("colormap", "viridis")).strip().lower()
    if scale not in SUPPORTED_SCALES:
        raise ValueError("scale must be 'lin' or 'log'.")
    if colormap not in SUPPORTED_COLORMAPS:
        choices = ", ".join(sorted(SUPPORTED_COLORMAPS))
        raise ValueError(f"colormap must be one of: {choices}.")

    return ComparisonConfig(
        sampling_config_file=resolve_path(config_file, sampling_config_value),
        output_directory=resolve_path(config_file, output_directory_value),
        generate_plots=generate_plots,
        scale=scale,
        fmin=plot_block.get("fmin"),
        fmax=plot_block.get("fmax"),
        enable_colorbar=enable_colorbar,
        colormap=colormap,
    )


def build_file_pairs(sampling_config_file):
    """Construct ordered reference-sampled pairs from the Stage 2 namelist."""
    if not sampling_config_file.is_file():
        raise FileNotFoundError(
            f"Sampling configuration not found: {sampling_config_file}"
        )

    blocks = f90nml.read(sampling_config_file)
    if "run_test" not in blocks:
        raise ValueError(f"Missing &run_test block in {sampling_config_file}.")
    run_block = blocks["run_test"]

    number_of_files = run_block.get("n_reference_files")
    reference_values = run_block.get("reference_files")
    sampled_directory_value = run_block.get("output_directory")

    if not isinstance(number_of_files, int) or number_of_files < 1:
        raise ValueError("n_reference_files in the sampling config must be positive.")
    if isinstance(reference_values, str):
        reference_values = [reference_values]
    if not isinstance(reference_values, list):
        raise ValueError("reference_files must be a path or an ordered path list.")
    if len(reference_values) != number_of_files:
        raise ValueError(
            "n_reference_files does not match the number of configured "
            "reference_files."
        )
    if not isinstance(sampled_directory_value, str) or not sampled_directory_value.strip():
        raise ValueError("The Stage 2 output_directory must be a non-empty path.")

    sampled_directory = resolve_path(sampling_config_file, sampled_directory_value)
    pairs = []
    used_basenames = set()

    # Preserve the explicit sampling order and ignore unrelated directory files.
    for reference_value in reference_values:
        if not isinstance(reference_value, str) or not reference_value.strip():
            raise ValueError("Every reference_files entry must be a non-empty path.")

        reference_file = resolve_path(sampling_config_file, reference_value)
        basename = reference_file.name
        if basename in used_basenames:
            raise ValueError(f"Duplicate reference basename: {basename}")
        used_basenames.add(basename)

        sampled_file = sampled_directory / basename
        if not reference_file.is_file():
            raise FileNotFoundError(f"Reference file not found: {reference_file}")
        if not sampled_file.is_file():
            raise FileNotFoundError(f"Sampled file not found: {sampled_file}")

        pair = FilePair(reference=reference_file, sampled=sampled_file)
        pairs.append(pair)

    return pairs
