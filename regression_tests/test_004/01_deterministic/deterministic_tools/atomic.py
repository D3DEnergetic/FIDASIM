"""Independent reader and interpolator for FIDASIM H-H CX cross sections."""

from dataclasses import dataclass

import h5py
import numpy as np

from regression_test_tools import ConfigError


N_LEVELS = 6


@dataclass
class ChargeExchangeTable:
    energy: np.ndarray
    log_cross_section: np.ndarray
    minimum_log_cross_section: float


def read_charge_exchange_table(filename):
    """Read the six-level H-H CX table using its on-disk HDF5 schema."""
    with h5py.File(filename, "r") as h5file:
        group = "/cross/H_H"
        if group not in h5file:
            raise ConfigError(f"{filename}: missing {group}.")
        energy = np.asarray(h5file[f"{group}/energy"][:], dtype=float)
        raw = np.asarray(h5file[f"{group}/cx"][:], dtype=float)

    if raw.ndim != 3 or raw.shape[0] != energy.size:
        raise ConfigError(f"{filename}: unexpected /cross/H_H/cx shape {raw.shape}.")
    if raw.shape[1] < N_LEVELS or raw.shape[2] < N_LEVELS:
        raise ConfigError(f"{filename}: H-H CX table has fewer than six levels.")

    # The Fortran table writer stores cx(initial, final, energy). HDF5 exposes
    # those dimensions to h5py in reverse order as (energy, final, initial).
    # Swap the level axes here so the Python-facing table has the explicit
    # semantic order (energy, initial, final) used by the calculation.
    cross_section = raw[:, :N_LEVELS, :N_LEVELS].swapaxes(1, 2)
    positive = cross_section[cross_section > 0.0]
    if positive.size == 0:
        raise ConfigError(f"{filename}: H-H CX table has no positive entries.")
    minimum = float(np.min(positive))
    prepared = np.where(cross_section > 0.0, cross_section, 0.9 * minimum)

    return ChargeExchangeTable(
        energy=energy,
        log_cross_section=np.log10(prepared),
        minimum_log_cross_section=float(np.log10(minimum)),
    )


def interpolate_cross_sections(table, relative_energy):
    """Reproduce `bb_cx_rates` log interpolation and endpoint clamping."""
    query = np.asarray(relative_energy, dtype=float)
    if np.any(query <= 0.0) or not np.all(np.isfinite(query)):
        raise ValueError("Relative energies must be finite and positive.")

    log_grid = np.log10(table.energy)
    log_query = np.log10(query)
    upper = np.searchsorted(log_grid, log_query, side="right")
    upper = np.clip(upper, 1, log_grid.size - 1)
    lower = upper - 1

    below = log_query <= log_grid[0]
    above = log_query >= log_grid[-1]
    fraction = (log_query - log_grid[lower]) / (
        log_grid[upper] - log_grid[lower]
    )
    fraction = np.where(below, 0.0, fraction)
    fraction = np.where(above, 1.0, fraction)
    lower = np.where(above, log_grid.size - 2, lower)
    upper = np.where(above, log_grid.size - 1, upper)
    lower = np.where(below, 0, lower)
    upper = np.where(below, 1, upper)

    low_values = table.log_cross_section[lower]
    high_values = table.log_cross_section[upper]
    interpolated_log = low_values + fraction[..., None, None] * (
        high_values - low_values
    )
    return np.where(
        interpolated_log < table.minimum_log_cross_section,
        0.0,
        np.power(10.0, interpolated_log),
    )
