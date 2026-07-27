"""Remap transformed values onto uniform FIDASIM cell centers."""

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from regression_test_tools import ConfigError


def remap_to_uniform_grid(
    nonuniform_energy,
    nonuniform_pitch,
    nonuniform_distribution,
    energy_upper_edge,
    nenergy,
    npitch,
):
    """Interpolate transformed values onto uniform cell-centered grids.

    Args:
        nonuniform_energy (array-like): Increasing transformed energy values.
        nonuniform_pitch (array-like): Increasing transformed pitch values.
        nonuniform_distribution (array-like): Values with shape
            ``(source_npitch, source_nenergy)``.
        energy_upper_edge (float): Energy at the upper edge of the final
            source proper-velocity cell, in keV.
        nenergy (int): Number of uniform FIDASIM energy cells.
        npitch (int): Number of uniform FIDASIM pitch cells.

    Returns:
        dict: Uniform energy and pitch centers, spacings, and ``f_array`` with
        FIDASIM shape ``(nenergy, npitch)``.
    """
    nonuniform_energy = np.asarray(nonuniform_energy, dtype=float)
    nonuniform_pitch = np.asarray(nonuniform_pitch, dtype=float)
    nonuniform_distribution = np.asarray(nonuniform_distribution, dtype=float)

    if nonuniform_energy.ndim != 1 or nonuniform_energy.size < 2:
        raise ConfigError("The transformed energy grid must be one-dimensional.")
    if nonuniform_pitch.ndim != 1 or nonuniform_pitch.size < 2:
        raise ConfigError("The transformed pitch grid must be one-dimensional.")
    expected_shape = (nonuniform_pitch.size, nonuniform_energy.size)
    if nonuniform_distribution.shape != expected_shape:
        raise ConfigError(
            "The transformed distribution has shape "
            f"{nonuniform_distribution.shape}; expected {expected_shape}."
        )
    if nenergy < 2 or npitch < 2:
        raise ConfigError("nenergy and npitch must both be at least 2.")
    if not np.all(np.diff(nonuniform_energy) > 0.0):
        raise ConfigError("The transformed energy grid must be increasing.")
    if not np.all(np.diff(nonuniform_pitch) > 0.0):
        raise ConfigError("The transformed pitch grid must be increasing.")

    energy_upper_edge = float(energy_upper_edge)
    if energy_upper_edge <= nonuniform_energy[-1]:
        raise ConfigError(
            "energy_upper_edge must be greater than the final energy coordinate."
        )
    denergy = energy_upper_edge / nenergy
    dpitch = 2.0 / npitch

    energy = (np.arange(nenergy) + 0.5) * denergy
    pitch = -1.0 + (np.arange(npitch) + 0.5) * dpitch

    interpolator = RegularGridInterpolator(
        (nonuniform_pitch, nonuniform_energy),
        nonuniform_distribution,
        method="linear",
        bounds_error=False,
        fill_value=None,
    )

    # RegularGridInterpolator expects points as (pitch, energy). Build the
    # complete target mesh explicitly and restore its two-dimensional shape.
    pitch_2d, energy_2d = np.meshgrid(pitch, energy, indexing="ij")
    interpolation_points = np.column_stack(
        (pitch_2d.ravel(), energy_2d.ravel())
    )
    distribution_pitch_energy = interpolator(interpolation_points)
    distribution_pitch_energy = distribution_pitch_energy.reshape(npitch, nenergy)

    # Linear extrapolation is used only in the final half source cell. Clip any
    # small negative extrapolated values because a distribution cannot be negative.
    distribution_pitch_energy = np.maximum(distribution_pitch_energy, 0.0)

    # The interpolation is assembled as (pitch, energy). Expose the converted
    # distribution to the rest of Stage 2 in the canonical (energy, pitch)
    # calculation order.
    return {
        "energy": energy,
        "pitch": pitch,
        "denergy": denergy,
        "dpitch": dpitch,
        "f_array": distribution_pitch_energy.T,
    }
