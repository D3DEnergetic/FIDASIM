"""Discover smooth Test 002 Stage 2 output distributions."""

from pathlib import Path
import re

import h5py
import numpy as np

from regression_test_tools import (
    ConfigError,
    normalize_path,
    read_namelist,
    require_string,
)


PARTICLE_DATASETS = (
    "species",
    "atomic_number",
    "mass_number",
    "charge_state",
    "A",
)


def _read_text_scalar(dataset, filename):
    """Read one scalar UTF-8 HDF5 dataset."""
    value = dataset[()]
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    value = str(value).strip().lower()
    if not value:
        raise ConfigError(f"{filename}: species must not be empty.")
    return value


def _read_scalar(dataset, filename):
    """Read one numerical HDF5 dataset containing exactly one value."""
    values = np.asarray(dataset[()]).reshape(-1)
    if values.size != 1:
        raise ConfigError(f"{filename}: {dataset.name} must be scalar.")
    return values[0]


def _read_case_metadata(filename):
    """Read the location and species parameters from one converted file."""
    required_datasets = (
        "energy",
        "pitch",
        "f",
        "denf",
        "r",
        "z",
        *PARTICLE_DATASETS,
    )
    try:
        with h5py.File(filename, mode="r") as h5file:
            missing = [
                dataset_name
                for dataset_name in required_datasets
                if dataset_name not in h5file
            ]
            if missing:
                raise ConfigError(
                    f"{filename} is missing datasets: {', '.join(missing)}"
                )

            return {
                "selected_r": float(_read_scalar(h5file["r"], filename)),
                "selected_z": float(_read_scalar(h5file["z"], filename)),
                "particle": {
                    "species": _read_text_scalar(
                        h5file["species"],
                        filename,
                    ),
                    "atomic_number": int(
                        _read_scalar(h5file["atomic_number"], filename)
                    ),
                    "mass_number": int(
                        _read_scalar(h5file["mass_number"], filename)
                    ),
                    "charge_state": int(
                        _read_scalar(h5file["charge_state"], filename)
                    ),
                    "A": float(_read_scalar(h5file["A"], filename)),
                },
            }
    except OSError as error:
        raise ConfigError(f"Could not read {filename}: {error}") from error


def _discover_indexed_paths(output_base, run_path):
    """Return matching three-digit output paths in contiguous index order."""
    filename_pattern = re.compile(
        rf"{re.escape(output_base.stem)}_(\d{{3}})"
        rf"{re.escape(output_base.suffix)}"
    )
    glob_pattern = f"{output_base.stem}_*{output_base.suffix}"

    indexed_paths = []
    for output_path in output_base.parent.glob(glob_pattern):
        match = filename_pattern.fullmatch(output_path.name)
        if match is not None and output_path.is_file():
            indexed_paths.append((int(match.group(1)), output_path.resolve()))
    indexed_paths.sort(key=lambda item: item[0])

    if not indexed_paths:
        command = f"cd {run_path.parent} && ./run.sh {run_path.name}"
        raise ConfigError(
            f"No indexed Test 002 Stage 2 outputs were found for {output_base}. "
            f"Generate them first from the Test 002 Stage 2 directory. "
            f"For example: {command}"
        )

    indices = [index for index, _ in indexed_paths]
    expected_indices = list(range(1, len(indexed_paths) + 1))
    if indices != expected_indices:
        raise ConfigError(
            "Test 002 Stage 2 output indices must be contiguous and start at 001."
        )
    return indexed_paths


def discover_distributions(run_config):
    """Return ordered paths, locations, and species parameters from Test 002.

    The returned dictionary has ``cases``, a list of dictionaries containing
    ``index``, ``input_path``, ``selected_r``, and ``selected_z``; ``particle``,
    the common ``species``, ``atomic_number``, ``mass_number``,
    ``charge_state``, and ``A`` values; and the resolved ``run_config`` path.
    """
    run_path, run_blocks = read_namelist(config_path=run_config)
    if "save_data_block" not in run_blocks:
        raise ConfigError(
            "Test 002 Stage 2 configuration is missing &save_data_block."
        )

    save_block = run_blocks["save_data_block"]
    if "output_filename" not in save_block:
        raise ConfigError(
            "Test 002 &save_data_block is missing output_filename."
        )

    output_filename = require_string(
        value=save_block["output_filename"],
        field_label="Test 002 output_filename",
    )
    output_base = normalize_path(
        value=output_filename,
        config_path=run_path,
        field_label="Test 002 output_filename",
    )

    cases = []
    collection_particle = None
    for index, input_path in _discover_indexed_paths(
        output_base=output_base,
        run_path=run_path,
    ):
        metadata = _read_case_metadata(filename=input_path)
        particle = metadata["particle"]
        if collection_particle is None:
            collection_particle = particle
        elif particle != collection_particle:
            raise ConfigError(
                f"{input_path}: species parameters differ from case 001."
            )

        cases.append(
            {
                "index": index,
                "input_path": input_path,
                "selected_r": metadata["selected_r"],
                "selected_z": metadata["selected_z"],
            }
        )

    return {
        "cases": cases,
        "particle": collection_particle,
        "run_config": run_path,
    }
