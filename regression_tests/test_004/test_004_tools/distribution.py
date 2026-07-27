"""Read and validate shared smooth Test 002 energy-pitch distributions."""

from dataclasses import dataclass

import h5py
import numpy as np

from regression_test_tools import ConfigError


@dataclass
class Distribution:
    energy: np.ndarray
    pitch: np.ndarray
    values: np.ndarray
    density: float
    atomic_mass: float
    species: str
    atomic_number: int
    mass_number: int
    charge_state: int
    selected_r: float
    selected_z: float


def _scalar(dataset):
    """Return the only numerical value stored in a scalar dataset."""
    values = np.asarray(dataset[...]).reshape(-1)
    if values.size != 1:
        raise ConfigError(f"{dataset.file.filename}: {dataset.name} must be scalar.")
    return float(values[0])


def read_distribution(filename):
    """Return one validated single-location Test 002 Stage 2 distribution."""
    with h5py.File(filename, "r") as h5file:
        required = [
            "energy",
            "pitch",
            "f",
            "denf",
            "species",
            "atomic_number",
            "mass_number",
            "charge_state",
            "A",
            "r",
            "z",
        ]
        missing = [name for name in required if name not in h5file]
        if missing:
            raise ConfigError(f"{filename} is missing: {', '.join(missing)}")

        energy = np.asarray(h5file["energy"][:], dtype=float)
        pitch = np.asarray(h5file["pitch"][:], dtype=float)
        raw = np.asarray(h5file["f"][:], dtype=float)
        density = _scalar(h5file["denf"])
        atomic_mass = _scalar(h5file["A"])
        atomic_number = int(_scalar(h5file["atomic_number"]))
        mass_number = int(_scalar(h5file["mass_number"]))
        charge_state = int(_scalar(h5file["charge_state"]))
        selected_r = _scalar(h5file["r"])
        selected_z = _scalar(h5file["z"])
        species = h5file["species"][()]

    if isinstance(species, bytes):
        species = species.decode("utf-8")
    species = str(species).strip().lower()

    if raw.shape != (1, 1, pitch.size, energy.size):
        raise ConfigError(
            f"{filename}: expected f shape (1,1,{pitch.size},{energy.size}), "
            f"got {raw.shape}."
        )

    # Test 002 writes the FIDASIM h5py schema as (z, r, pitch, energy).
    # Select the single spatial cell and normalize (pitch, energy) into the
    # Test 004 canonical calculation order (energy, pitch). The source file is
    # Python-written but deliberately uses the layout expected by Fortran.
    values = raw[0, 0, :, :].T

    arrays = [energy, pitch, values]
    if any(not np.all(np.isfinite(array)) for array in arrays):
        raise ConfigError(f"{filename}: distribution data must be finite.")
    if energy.size < 2 or pitch.size < 2:
        raise ConfigError(f"{filename}: energy and pitch need at least two points.")
    if np.any(values < 0.0) or np.sum(values) <= 0.0:
        raise ConfigError(f"{filename}: f must be nonnegative with positive sum.")
    if density <= 0.0 or not np.isfinite(density):
        raise ConfigError(f"{filename}: denf must be finite and positive.")
    if atomic_mass <= 0.0 or not np.isfinite(atomic_mass):
        raise ConfigError(f"{filename}: A must be finite and positive.")

    denergy = np.diff(energy)
    dpitch = np.diff(pitch)
    if np.any(denergy <= 0.0) or np.any(dpitch <= 0.0):
        raise ConfigError(f"{filename}: grids must be strictly increasing.")
    if not np.allclose(denergy, denergy[0]) or not np.allclose(
        dpitch, dpitch[0]
    ):
        raise ConfigError(f"{filename}: grids must be uniformly spaced.")

    return Distribution(
        energy=energy,
        pitch=pitch,
        values=values,
        density=density,
        atomic_mass=atomic_mass,
        species=species,
        atomic_number=atomic_number,
        mass_number=mass_number,
        charge_state=charge_state,
        selected_r=selected_r,
        selected_z=selected_z,
    )
