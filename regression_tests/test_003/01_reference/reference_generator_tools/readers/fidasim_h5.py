"""Reader for FIDASIM HDF5 distribution data."""

import numpy as np
import h5py


def load_fidasim_h5_species_parameters(input_path):
    """Load scalar particle metadata from a Test 002 Stage 2 output."""
    required_datasets = (
        "species",
        "atomic_number",
        "mass_number",
        "charge_state",
        "A",
    )
    with h5py.File(input_path, "r") as h5f:
        missing_datasets = [
            dataset_name
            for dataset_name in required_datasets
            if dataset_name not in h5f
        ]
        if missing_datasets:
            missing_names = ", ".join(missing_datasets)
            raise KeyError(f"Missing particle datasets: {missing_names}")

        species = h5f["species"][()]
        particle = {
            "atomic_number": int(h5f["atomic_number"][()]),
            "mass_number": int(h5f["mass_number"][()]),
            "charge_state": int(h5f["charge_state"][()]),
            "A": float(h5f["A"][()]),
        }

    if isinstance(species, bytes):
        species = species.decode("utf-8")
    particle["species"] = str(species).strip().lower()
    return particle


def load_fidasim_h5_distribution(input_path):
    """Load a FIDASIM file in its native h5py-visible array order.

    Args:
        input_path (Path): FIDASIM distribution HDF5 file to read.

    Returns:
        tuple: ``z, r, pitch, energy, f, denf``, where the one-dimensional
        arrays define the coordinate grids, ``f`` retains the FIDASIM file
        order ``(nz, nr, npitch, nenergy)``, and ``denf`` has shape
        ``(nz, nr)``. The caller selects the spatial cell and converts its
        distribution to canonical ``(nenergy, npitch)`` order.

    Raises:
        KeyError: If a required dataset, including ``denf``, is absent.
        ValueError: If the shape of ``f`` or ``denf`` does not match the
        coordinate grids.
    """
    with h5py.File(input_path, "r") as h5f:
        required = ["z", "r", "pitch", "energy", "f", "denf"]
        missing = [name for name in required if name not in h5f]
        if missing:
            raise KeyError(f"Missing datasets in input file: {missing}")

        z = np.asarray(h5f["z"], dtype=float)
        r = np.asarray(h5f["r"], dtype=float)
        pitch = np.asarray(h5f["pitch"], dtype=float)
        energy = np.asarray(h5f["energy"], dtype=float)
        f = np.asarray(h5f["f"], dtype=float)
        denf = np.asarray(h5f["denf"], dtype=float)

        expected_shape = (len(z), len(r), len(pitch), len(energy))
        if tuple(f.shape) != expected_shape:
            raise ValueError(
                f"Input f dataset has shape {f.shape}, expected {expected_shape}"
            )

        expected_denf_shape = (len(z), len(r))
        if tuple(denf.shape) != expected_denf_shape:
            raise ValueError(
                "Input denf dataset has shape "
                f"{denf.shape}, expected {expected_denf_shape}"
            )

    return z, r, pitch, energy, f, denf
