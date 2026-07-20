"""Read and validate the shared energy-pitch HDF5 data contract."""

from dataclasses import dataclass

import h5py
import numpy as np


REQUIRED_DATASETS = (
    "energy_grid",
    "pitch_grid",
    "f_array",
    "species",
    "atomic_number",
    "mass_number",
    "charge_state",
    "requested_r",
    "requested_z",
    "selected_r",
    "selected_z",
    "r_index",
    "z_index",
)

MATCHED_METADATA = (
    "species",
    "atomic_number",
    "mass_number",
    "charge_state",
    "requested_r",
    "requested_z",
    "selected_r",
    "selected_z",
    "r_index",
    "z_index",
)


@dataclass
class DistributionData:
    energy: np.ndarray
    pitch: np.ndarray
    values: np.ndarray
    metadata: dict
    units: str


def _decode_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _read_scalar(dataset):
    values = np.asarray(dataset[()]).reshape(-1)
    if values.size != 1:
        raise ValueError(f"{dataset.name} must contain exactly one value.")

    value = values[0]
    if isinstance(value, bytes) or isinstance(value, str):
        return _decode_text(value)
    return value.item() if hasattr(value, "item") else value


def _validate_uniform_grid(grid, grid_name, filename):
    """Require a finite, strictly monotonic, uniformly spaced grid."""
    if grid.ndim != 1 or grid.size < 2:
        raise ValueError(f"{filename}: {grid_name} must have at least two values.")
    if not np.all(np.isfinite(grid)):
        raise ValueError(f"{filename}: {grid_name} contains non-finite values.")

    differences = np.diff(grid)
    increasing = np.all(differences > 0.0)
    decreasing = np.all(differences < 0.0)
    if not increasing and not decreasing:
        raise ValueError(f"{filename}: {grid_name} must be strictly monotonic.")

    scale = max(1.0, float(np.max(np.abs(grid))), abs(float(differences[0])))
    tolerance = 100.0 * np.finfo(float).eps * scale
    if not np.allclose(differences, differences[0], rtol=0.0, atol=tolerance):
        raise ValueError(f"{filename}: {grid_name} must be uniformly spaced.")


def read_distribution(filename):
    """Read one distribution and enforce the shared Stage 1/Stage 2 schema."""
    with h5py.File(filename, "r") as h5file:
        missing = [name for name in REQUIRED_DATASETS if name not in h5file]
        if missing:
            names = ", ".join(missing)
            raise ValueError(f"{filename}: missing required datasets: {names}")

        energy = np.asarray(h5file["energy_grid"][:], dtype=float)
        pitch = np.asarray(h5file["pitch_grid"][:], dtype=float)
        values = np.asarray(h5file["f_array"][:], dtype=float)
        metadata = {
            name: _read_scalar(h5file[name]) for name in MATCHED_METADATA
        }
        units = _decode_text(h5file["f_array"].attrs.get("units", ""))

    _validate_uniform_grid(energy, "energy_grid", filename)
    _validate_uniform_grid(pitch, "pitch_grid", filename)

    expected_shape = (energy.size, pitch.size)
    if values.shape != expected_shape:
        raise ValueError(
            f"{filename}: f_array shape {values.shape} does not match "
            f"{expected_shape}."
        )
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{filename}: f_array contains non-finite values.")
    if np.any(values < 0.0):
        raise ValueError(f"{filename}: f_array contains negative values.")
    if float(np.sum(values)) <= 0.0:
        raise ValueError(f"{filename}: f_array must have a positive sum.")

    return DistributionData(
        energy=energy,
        pitch=pitch,
        values=values,
        metadata=metadata,
        units=units,
    )


def validate_pair(reference, sampled, pair):
    """Require identical grids, metadata, array shape, and distribution units."""
    if not np.array_equal(reference.energy, sampled.energy):
        raise ValueError(f"{pair.sampled}: energy_grid does not match the reference.")
    if not np.array_equal(reference.pitch, sampled.pitch):
        raise ValueError(f"{pair.sampled}: pitch_grid does not match the reference.")
    if reference.values.shape != sampled.values.shape:
        raise ValueError(f"{pair.sampled}: f_array shape does not match the reference.")
    if reference.units != sampled.units:
        raise ValueError(f"{pair.sampled}: f_array units do not match the reference.")

    for name in MATCHED_METADATA:
        reference_value = reference.metadata[name]
        sampled_value = sampled.metadata[name]
        if reference_value != sampled_value:
            raise ValueError(f"{pair.sampled}: {name} does not match the reference.")
