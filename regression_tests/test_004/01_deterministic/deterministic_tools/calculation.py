"""Deterministic energy-pitch-gyrophase ion-sink calculation."""

from dataclasses import dataclass

import numpy as np

from .atomic import N_LEVELS, interpolate_cross_sections


MASS_U = 1.660539040e-27
ELEMENTARY_CHARGE = 1.60217733e-19
V2_TO_ENERGY_PER_AMU = MASS_U / (2.0 * ELEMENTARY_CHARGE * 1.0e3) * 1.0e-4


@dataclass
class DeterministicResult:
    sink_distribution: np.ndarray
    energy_marginal: np.ndarray
    pitch_marginal: np.ndarray
    kernel: np.ndarray
    total_rate: float
    neutral_velocity: np.ndarray
    neutral_energy: float
    neutral_pitch: float
    level_density: np.ndarray
    gyroangle: np.ndarray


def build_level_density(neutral_config):
    """Construct the six-level neutral density vector."""
    density = neutral_config["density"]
    method = neutral_config["level_split_method"]
    if method == "ground-only":
        fractions = np.zeros(N_LEVELS)
        fractions[0] = 1.0
    else:
        levels = np.arange(N_LEVELS, dtype=float)
        fractions = np.exp(-neutral_config["level_decay"] * levels)
        fractions /= np.sum(fractions)
    return density * fractions


def calculate_deterministic(distribution, table, neutral_config, n_gyro):
    """Calculate the deterministic reaction-weighted sink distribution."""
    angle = np.deg2rad(neutral_config["injection_angle"])
    neutral_speed = np.sqrt(
        neutral_config["energy"]
        / (V2_TO_ENERGY_PER_AMU * distribution.atomic_mass)
    )
    neutral_velocity = neutral_speed * np.array(
        [np.sin(angle), 0.0, np.cos(angle)]
    )
    neutral_pitch = float(neutral_velocity[2] / neutral_speed)
    level_density = build_level_density(neutral_config)

    gyroangle = 2.0 * np.pi * (
        np.arange(n_gyro, dtype=float) + 0.5
    ) / n_gyro
    energy = distribution.energy[:, None, None]
    pitch = distribution.pitch[None, :, None]
    gyro = gyroangle[None, None, :]

    ion_speed = np.sqrt(
        energy / (V2_TO_ENERGY_PER_AMU * distribution.atomic_mass)
    )
    perpendicular = np.sqrt(np.maximum(1.0 - pitch**2, 0.0))
    vx = ion_speed * perpendicular * np.cos(gyro)
    charge_sign = np.sign(distribution.charge_state)
    vy = -charge_sign * ion_speed * perpendicular * np.sin(gyro)
    vz = ion_speed * pitch * np.ones_like(gyro)

    relative_speed = np.sqrt(
        (vx - neutral_velocity[0]) ** 2
        + (vy - neutral_velocity[1]) ** 2
        + (vz - neutral_velocity[2]) ** 2
    )
    relative_energy = V2_TO_ENERGY_PER_AMU * relative_speed**2
    cross_section = interpolate_cross_sections(table, relative_energy)

    # The table reader presents the final two axes as (initial, final).
    state_rates = np.einsum(
        "...lm,l,...->...m",
        cross_section,
        level_density,
        relative_speed,
        optimize=True,
    )
    kernel = np.mean(np.sum(state_rates, axis=-1), axis=-1)

    denergy = float(distribution.energy[1] - distribution.energy[0])
    dpitch = float(distribution.pitch[1] - distribution.pitch[0])
    normalization = float(np.sum(distribution.values) * denergy * dpitch)
    probability_density = distribution.values / normalization
    sink_distribution = distribution.density * probability_density * kernel
    total_rate = float(np.sum(sink_distribution) * denergy * dpitch)
    energy_marginal = np.sum(sink_distribution, axis=1) * dpitch
    pitch_marginal = np.sum(sink_distribution, axis=0) * denergy

    return DeterministicResult(
        sink_distribution=sink_distribution,
        energy_marginal=energy_marginal,
        pitch_marginal=pitch_marginal,
        kernel=kernel,
        total_rate=total_rate,
        neutral_velocity=neutral_velocity,
        neutral_energy=neutral_config["energy"],
        neutral_pitch=neutral_pitch,
        level_density=level_density,
        gyroangle=gyroangle,
    )
