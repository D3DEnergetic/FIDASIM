"""Calculate physical moments from a single-location CQL3D distribution."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from regression_test_tools import ConfigError


# Physical constants and isotope masses match the values used by FIDASIM.
SPEED_OF_LIGHT_CM_PER_SECOND = 2.99792458e10
ATOMIC_MASS_ENERGY_KEV = 9.3149410242e5
SPECIES_MASS_AMU = {
    "h": 1.007276466879,
    "d": 2.013553212745,
    "t": 3.01550071632,
}


@dataclass
class PhysicalMoments:
    """Density and directional temperatures calculated from one distribution."""

    density: float
    parallel_temperature: float
    perpendicular_temperature: float


def write_physical_moments(output_path, moments):
    """Save calculated physical moments in a reference HDF5 file.

    Args:
        output_path (str or Path): Reference HDF5 file to update.
        moments (PhysicalMoments): Calculated density and temperatures.
    """
    output_path = Path(output_path)

    try:
        with h5py.File(output_path, mode="a") as h5file:
            moments_group = h5file.create_group("moments")

            density = moments_group.create_dataset("density", data=moments.density)
            density.attrs["units"] = "ions/cm^3"
            density.attrs["description"] = "Fast-ion number density"

            parallel_temperature = moments_group.create_dataset(
                "parallel_temperature",
                data=moments.parallel_temperature,
            )
            parallel_temperature.attrs["units"] = "keV"
            parallel_temperature.attrs["description"] = (
                "Parallel temperature from the relativistic pressure moment"
            )

            perpendicular_temperature = moments_group.create_dataset(
                "perpendicular_temperature",
                data=moments.perpendicular_temperature,
            )
            perpendicular_temperature.attrs["units"] = "keV"
            perpendicular_temperature.attrs["description"] = (
                "Perpendicular temperature from the relativistic pressure moment "
                "per perpendicular degree of freedom"
            )
    except OSError as error:
        raise ConfigError(
            f"Could not write physical moments to {output_path}: {error}"
        ) from error


def _read_moment_inputs(input_path):
    """Read and validate the arrays needed by the moment calculation."""
    required_datasets = [
        "u_norm",
        "u_bar",
        "theta",
        "u_bar_weight",
        "theta_weight",
        "f_u_theta",
    ]

    try:
        with h5py.File(input_path, mode="r") as h5file:
            missing_datasets = []
            for dataset_name in required_datasets:
                if dataset_name not in h5file:
                    missing_datasets.append(dataset_name)
            if missing_datasets:
                names = ", ".join(missing_datasets)
                raise ConfigError(
                    f"Reference file is missing moment datasets: {names}"
                )

            u_norm = float(h5file["u_norm"][()])
            u_bar = np.asarray(h5file["u_bar"][:], dtype=float)
            theta = np.asarray(h5file["theta"][:], dtype=float)

            # These are complete integration weights, not bare grid spacings.
            u_bar_weight = np.asarray(h5file["u_bar_weight"][:], dtype=float)
            theta_weight = np.asarray(h5file["theta_weight"][:], dtype=float)
            distribution = np.asarray(h5file["f_u_theta"][:], dtype=float)
    except OSError as error:
        raise ConfigError(f"Could not read reference moments from {input_path}: {error}") from error

    expected_shape = (
        theta.size,
        u_bar.size,
    )
    if distribution.shape != expected_shape:
        raise ConfigError(
            f"f4d in {input_path} has shape {distribution.shape}; "
            f"expected {expected_shape}."
        )
    if u_bar_weight.shape != u_bar.shape:
        raise ConfigError("u_bar_weight must have the same shape as u_bar.")
    if theta_weight.shape != theta.shape:
        raise ConfigError("theta_weight must have the same shape as theta.")
    if u_norm <= 0.0:
        raise ConfigError("u_norm must be greater than zero.")

    return {
        "u_norm": u_norm,
        "u_bar": u_bar,
        "theta": theta,
        "u_bar_weight": u_bar_weight,
        "theta_weight": theta_weight,
        "distribution": distribution,
    }


def calculate_cql3d_moments(input_path, species):
    """Calculate density, parallel temperature, and perpendicular temperature.

    The CQL3D arrays ``f4ddv`` and ``f4ddt`` already contain the normalized
    velocity-space cell weights ``u**2 du`` and ``2*pi*sin(theta)*dtheta``.
    Directional temperatures are obtained from the fully relativistic pressure
    moments ``p_parallel*v_parallel`` and ``p_perpendicular*v_perpendicular``.

    Args:
        input_path (str or Path): Single-location Test 002 HDF5 file.
        species (str): Supported isotope identifier: ``h``, ``d``, or ``t``.

    Returns:
        PhysicalMoments: Density in ions/cm^3 and temperatures in keV.

    Raises:
        ConfigError: If required data are invalid or density is not positive.
    """
    input_path = Path(input_path)
    arrays = _read_moment_inputs(input_path=input_path)

    u = arrays["u_bar"] * arrays["u_norm"]
    u_over_c = u / SPEED_OF_LIGHT_CM_PER_SECOND
    lorentz_factor = np.sqrt(1.0 + u_over_c**2)
    rest_mass_energy = SPECIES_MASS_AMU[species] * ATOMIC_MASS_ENERGY_KEV
    # The relativistic pressure moment contains p*v = m*u**2/gamma. Express
    # this factor in keV so the resulting directional temperatures are in keV.
    pressure_energy = (
        rest_mass_energy * u_over_c**2 / lorentz_factor
    )

    # Form the two-dimensional CQL3D velocity-space integration weights.
    cell_weight = (
        arrays["theta_weight"][:, np.newaxis]
        * arrays["u_bar_weight"][np.newaxis, :]
    )
    weighted_distribution = arrays["distribution"] * cell_weight

    density = float(np.sum(weighted_distribution))
    if density <= 0.0:
        raise ConfigError(
            f"Cannot calculate temperatures because density is not positive: {density}"
        )

    pitch_cosine = np.cos(arrays["theta"])
    parallel_pressure = float(
        np.sum(
            pressure_energy[np.newaxis, :]
            * pitch_cosine[:, np.newaxis] ** 2
            * weighted_distribution
        )
    )
    perpendicular_pressure = 0.5 * float(
        np.sum(
            pressure_energy[np.newaxis, :]
            * (1.0 - pitch_cosine[:, np.newaxis] ** 2)
            * weighted_distribution
        )
    )

    parallel_temperature = parallel_pressure / density
    perpendicular_temperature = perpendicular_pressure / density

    return PhysicalMoments(
        density=density,
        parallel_temperature=parallel_temperature,
        perpendicular_temperature=perpendicular_temperature,
    )
