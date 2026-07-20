"""Calculate marginals, physical moments, and comparison errors."""

from dataclasses import dataclass

import numpy as np


@dataclass
class PhysicalMoments:
    density: float
    parallel_temperature: float
    perpendicular_temperature: float


@dataclass
class ComparisonResult:
    basename: str
    reference_moments: PhysicalMoments
    sampled_moments: PhysicalMoments
    density_error: float
    parallel_temperature_error: float
    perpendicular_temperature_error: float


def calculate_marginals(distribution):
    """Integrate f(E, pitch) over pitch and energy, respectively."""
    denergy = abs(float(distribution.energy[1] - distribution.energy[0]))
    dpitch = abs(float(distribution.pitch[1] - distribution.pitch[0]))

    energy_marginal = np.sum(distribution.values, axis=1) * dpitch
    pitch_marginal = np.sum(distribution.values, axis=0) * denergy
    return energy_marginal, pitch_marginal


def calculate_physical_moments(distribution):
    """Calculate density, parallel temperature, and perpendicular temperature."""
    denergy = abs(float(distribution.energy[1] - distribution.energy[0]))
    dpitch = abs(float(distribution.pitch[1] - distribution.pitch[0]))
    cell_area = denergy * dpitch

    # Broadcast the one-dimensional coordinates over the 2D distribution.
    energy_2d = distribution.energy[:, np.newaxis]
    pitch_2d = distribution.pitch[np.newaxis, :]

    density = float(np.sum(distribution.values) * cell_area)
    parallel_energy_density = float(
        np.sum(energy_2d * pitch_2d**2 * distribution.values) * cell_area
    )
    perpendicular_energy_density = float(
        np.sum(energy_2d * (1.0 - pitch_2d**2) * distribution.values) * cell_area
    )

    parallel_temperature = 2.0 * parallel_energy_density / density
    perpendicular_temperature = perpendicular_energy_density / density

    return PhysicalMoments(
        density=density,
        parallel_temperature=parallel_temperature,
        perpendicular_temperature=perpendicular_temperature,
    )


def _relative_error(reference_value, sampled_value):
    if reference_value == 0.0:
        raise ValueError("Cannot calculate relative error from a zero reference value.")
    return abs(sampled_value - reference_value) / abs(reference_value)


def compare_moments(basename, reference_moments, sampled_moments):
    """Assemble the scalar comparison results for one file pair."""
    density_error = _relative_error(
        reference_moments.density,
        sampled_moments.density,
    )
    parallel_temperature_error = _relative_error(
        reference_moments.parallel_temperature,
        sampled_moments.parallel_temperature,
    )
    perpendicular_temperature_error = _relative_error(
        reference_moments.perpendicular_temperature,
        sampled_moments.perpendicular_temperature,
    )

    return ComparisonResult(
        basename=basename,
        reference_moments=reference_moments,
        sampled_moments=sampled_moments,
        density_error=density_error,
        parallel_temperature_error=parallel_temperature_error,
        perpendicular_temperature_error=perpendicular_temperature_error,
    )
