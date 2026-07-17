"""Reader for FIDASIM HDF5 distribution data."""

import numpy as np
import h5py


def load_fidasim_h5_distribution(input_path):
    """Load a FIDASIM file into the canonical energy-pitch representation."""
    with h5py.File(input_path, "r") as h5f:
        required = ["z", "r", "pitch", "energy", "f"]
        missing = [name for name in required if name not in h5f]
        if missing:
            raise KeyError(f"Missing datasets in input file: {missing}")

        z = np.asarray(h5f["z"], dtype=float)
        r = np.asarray(h5f["r"], dtype=float)
        pitch = np.asarray(h5f["pitch"], dtype=float)
        energy = np.asarray(h5f["energy"], dtype=float)
        f = np.asarray(h5f["f"], dtype=float)

        expected_shape = (len(z), len(r), len(pitch), len(energy))
        if tuple(f.shape) != expected_shape:
            raise ValueError(
                f"Input f dataset has shape {f.shape}, expected {expected_shape}"
            )

    return z, r, pitch, energy, f
