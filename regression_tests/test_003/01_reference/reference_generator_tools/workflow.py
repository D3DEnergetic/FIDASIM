"""Validate the shared Test 002 distributions used by Test 003."""

from pathlib import Path
import re

import numpy as np

from regression_test_tools import (
    ConfigError,
    normalize_path,
    read_namelist,
    require_string,
)

from .config import read_config
from .readers import (
    load_fidasim_h5_distribution,
    load_fidasim_h5_species_parameters,
)


def _discover_test_002_outputs(run_config):
    """Return contiguous indexed outputs from a Test 002 Stage 2 namelist."""
    run_path, blocks = read_namelist(config_path=run_config)
    if "save_data_block" not in blocks:
        raise ConfigError(
            "Test 002 Stage 2 configuration is missing save_data_block."
        )
    save_block = blocks["save_data_block"]
    if "output_filename" not in save_block:
        raise ConfigError(
            "Test 002 save_data_block is missing output_filename."
        )

    output_value = require_string(
        value=save_block["output_filename"],
        field_label="Test 002 output_filename",
    )
    output_base = normalize_path(
        value=output_value,
        config_path=run_path,
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
        raise ConfigError(
            "No Test 002 Stage 2 outputs were found. Run "
            f"{run_path.parent}/run.sh {run_path.name} first."
        )
    if indices != list(range(1, len(indices) + 1)):
        raise ConfigError(
            "Test 002 output indices must be contiguous and start at 001."
        )
    return [path for _, path in indexed_paths]


def _validate_distribution(source_path):
    """Validate one native Test 002 energy-pitch distribution."""
    z, r, pitch, energy, f_array, density = load_fidasim_h5_distribution(
        source_path
    )
    if r.size != 1 or z.size != 1:
        raise ConfigError(
            f"Expected one R and one Z location in {source_path}."
        )

    # Test 002 writes the FIDASIM h5py schema as (z, r, pitch, energy).
    # Select its single spatial cell and normalize (pitch, energy) into the
    # Test 003 canonical (energy, pitch) order. This transpose reflects the
    # file schema; the Test 002 file itself was written by Python.
    values = np.transpose(f_array[0, 0, :, :])
    if energy.size < 2 or pitch.size < 2:
        raise ConfigError(f"Energy and pitch grids are too short in {source_path}.")
    if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(pitch)):
        raise ConfigError(f"A coordinate grid is non-finite in {source_path}.")
    if not np.all(np.isfinite(values)) or np.any(values < 0.0):
        raise ConfigError(f"The distribution is invalid in {source_path}.")
    if not np.isfinite(density[0, 0]) or density[0, 0] <= 0.0:
        raise ConfigError(f"denf must be finite and positive in {source_path}.")

    return {
        "input_path": Path(source_path),
        "selected_r": float(r[0]),
        "selected_z": float(z[0]),
        "shape": values.shape,
        "density": float(density[0, 0]),
    }


def validate_references(config_path):
    """Validate and summarize the shared Test 002 reference collection."""
    config = read_config(config_filename=config_path)
    distribution_config = config["reference"]["input_distribution_config"]
    source_paths = _discover_test_002_outputs(run_config=distribution_config)

    cases = []
    collection_particle = None
    for source_path in source_paths:
        cases.append(_validate_distribution(source_path=source_path))
        particle = load_fidasim_h5_species_parameters(input_path=source_path)
        if collection_particle is None:
            collection_particle = particle
        elif particle != collection_particle:
            raise ConfigError(
                f"{source_path}: species parameters differ from case 001."
            )

    return {
        "cases": cases,
        "particle": collection_particle,
        "input_distribution_config": distribution_config,
    }
