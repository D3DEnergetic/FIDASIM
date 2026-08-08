"""Transform a CQL3D distribution to nonrelativistic energy-pitch coordinates."""

import numpy as np

from regression_test_tools import ConfigError


SPEED_OF_LIGHT_CM_PER_SECOND = 2.99792458e10
ERG_PER_KEV = 1.602176634e-9
ATOMIC_MASS_GRAMS = 1.66053906660e-24

# CQL3D may store the normalized velocity grid as either float32 or float64.
# For float32 grids, rounding causes small differences between spacings that
# are mathematically uniform, so the comparison must allow for that precision.
UNIFORM_SPACING_RELATIVE_TOLERANCE = 2.0e-5
UNIFORM_SPACING_ABSOLUTE_TOLERANCE = 1.0e-12


def _validate_inputs(f_u_theta, u_bar, theta, u_norm, mass_amu):
    """Check the arrays and physical values used by the transformation."""
    expected_shape = (theta.size, u_bar.size)
    if f_u_theta.shape != expected_shape:
        raise ConfigError(
            f"f_u_theta has shape {f_u_theta.shape}; expected {expected_shape}."
        )
    if u_bar.ndim != 1 or theta.ndim != 1:
        raise ConfigError("u_bar and theta must be one-dimensional arrays.")
    if u_bar.size < 2 or not np.all(np.diff(u_bar) > 0.0):
        raise ConfigError("u_bar must contain at least two increasing values.")
    du_bar = u_bar[1] - u_bar[0]
    if not np.allclose(
        np.diff(u_bar),
        du_bar,
        rtol=UNIFORM_SPACING_RELATIVE_TOLERANCE,
        atol=UNIFORM_SPACING_ABSOLUTE_TOLERANCE,
    ):
        raise ConfigError("u_bar must be uniformly spaced.")
    if u_norm <= 0.0:
        raise ConfigError("u_norm must be greater than zero.")
    if mass_amu <= 0.0:
        raise ConfigError("The isotope mass A must be greater than zero.")
    if not np.all(np.isfinite(f_u_theta)) or np.any(f_u_theta < 0.0):
        raise ConfigError("f_u_theta must contain finite, nonnegative values.")


def transform_to_nonrelativistic_energy_pitch(
    f_u_theta,
    u_bar,
    theta,
    u_norm,
    mass_amu,
):
    """Transform CQL3D values onto their corresponding nonuniform E-P points.

    This function applies the conventional nonrelativistic ``sqrt(E)``
    expression. It does not yet remap the result onto the uniform grids required
    by FIDASIM.

    Here ``epsilon`` is the energy conversion ``ERG_PER_KEV``. For
    ``E = m*u**2/(2*epsilon)`` and ``P = cos(theta)``, particle-number
    conservation gives

    ``F(E, P) dE dP = 2*pi*f_v(u, theta)*u**2 du dP``.

    The resulting phase-space Jacobian factor is

    ``2*pi*sqrt(2)*(epsilon/m)**(3/2)*sqrt(E)``.

    The factor includes integration over gyrophase and the velocity-to-energy
    coordinate Jacobian. It is multiplied by the physical velocity-space
    distribution at the point where ``f_energy_pitch`` is constructed below.

    Args:
        f_u_theta (array-like): CQL3D distribution with shape ``(ntheta, nu)``.
        u_bar (array-like): Normalized proper-velocity coordinates.
        theta (array-like): Pitch-angle coordinates in radians.
        u_norm (float): Proper-velocity normalization in cm/s.
        mass_amu (float): Physical isotope mass in atomic mass units.

    Returns:
        dict: Nonuniform energy grid, increasing pitch grid, transformed
        distribution, and relativistic comparison diagnostics.
    """
    f_u_theta = np.asarray(f_u_theta, dtype=float)
    u_bar = np.asarray(u_bar, dtype=float)
    theta = np.asarray(theta, dtype=float)
    u_norm = float(u_norm)
    mass_amu = float(mass_amu)

    _validate_inputs(
        f_u_theta=f_u_theta,
        u_bar=u_bar,
        theta=theta,
        u_norm=u_norm,
        mass_amu=mass_amu,
    )

    mass_grams = mass_amu * ATOMIC_MASS_GRAMS
    u = u_bar * u_norm

    # FIDASIM interprets energy using the nonrelativistic relation E=m*u^2/2.
    energy = 0.5 * mass_grams * u**2 / ERG_PER_KEV
    pitch = np.cos(theta)

    # The final source coordinate is a cell center. Transform the upper edge
    # half a source spacing beyond it for construction of the target grid.
    du_bar = u_bar[1] - u_bar[0]
    u_bar_upper_edge = u_bar[-1] + 0.5 * du_bar
    u_upper_edge = u_bar_upper_edge * u_norm
    energy_upper_edge = 0.5 * mass_grams * u_upper_edge**2 / ERG_PER_KEV

    # Undo CQL3D's u_norm^3 scaling to recover the physical velocity-space
    # distribution in ions/[cm^3*(cm/s)^3].
    f_velocity = f_u_theta / u_norm**3

    # Combine the velocity-to-energy Jacobian with the 2*pi gyrophase
    # integration. Multiplication by this factor converts f_velocity into the
    # FIDASIM distribution density per unit energy and pitch.
    phase_space_jacobian = (
        2.0
        * np.pi
        * np.sqrt(2.0)
        * (ERG_PER_KEV / mass_grams) ** 1.5
        * np.sqrt(energy)
    )
    f_energy_pitch = phase_space_jacobian[np.newaxis, :] * f_velocity

    # Increasing theta produces decreasing pitch. Reverse both pitch and the
    # corresponding distribution axis for the FIDASIM convention.
    if pitch[0] > pitch[-1]:
        pitch = pitch[::-1]
        f_energy_pitch = f_energy_pitch[::-1, :]

    # Quantify the approximation without using the relativistic result as the
    # FIDASIM output. The stable gamma-1 expression avoids cancellation.
    u_over_c = u / SPEED_OF_LIGHT_CM_PER_SECOND
    gamma = np.sqrt(1.0 + u_over_c**2)
    rest_mass_energy_kev = (
        mass_grams * SPEED_OF_LIGHT_CM_PER_SECOND**2 / ERG_PER_KEV
    )
    relativistic_energy = rest_mass_energy_kev * u_over_c**2 / (gamma + 1.0)
    relativistic_f = gamma[np.newaxis, :] * f_energy_pitch

    return {
        "energy": energy,
        "energy_upper_edge": energy_upper_edge,
        "pitch": pitch,
        "f_energy_pitch": f_energy_pitch,
        "relativistic_energy": relativistic_energy,
        "relativistic_f": relativistic_f,
    }
