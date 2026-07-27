"""Read and validate deterministic and Monte Carlo ion-sink artifacts."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


DATASET_UNITS = {
    "energy": "keV",
    "pitch": "dimensionless",
    "f_array": "ions/(cm^3*keV*dP)",
    "denf": "ions/cm^3",
    "sink_distribution": "ions/(cm^3*s*keV*dP)",
    "energy_marginal": "ions/(cm^3*s*keV)",
    "pitch_marginal": "ions/(cm^3*s*dP)",
    "total_reaction_rate": "ions/(cm^3*s)",
    "neutral_velocity": "cm/s",
    "neutral_level_density": "neutrals/cm^3",
    "neutral_energy": "keV",
    "injection_angle": "degree",
    "atomic_number": "dimensionless",
    "mass_number": "dimensionless",
    "charge_state": "elementary charge",
    "A": "amu",
    "selected_r": "cm",
    "selected_z": "cm",
}

DETERMINISTIC_DATASET_UNITS = {
    "gyroangle": "rad",
}

MONTE_CARLO_DATASET_UNITS = {
    "sample_standard_deviation": "ions/(cm^3*s)",
    "standard_error": "ions/(cm^3*s)",
    "n_markers": "dimensionless",
    "reservoir_size": "dimensionless",
    "seed": "dimensionless",
    "neutral_density": "neutrals/cm^3",
}

SHARED_NUMERICAL_DATASETS = (
    "energy",
    "pitch",
    "f_array",
    "denf",
    "neutral_velocity",
    "neutral_level_density",
    "neutral_energy",
    "injection_angle",
    "atomic_number",
    "mass_number",
    "charge_state",
    "A",
    "selected_r",
    "selected_z",
)

PROVENANCE_ATTRIBUTES = (
    "source_distribution",
    "input_distribution_config",
    "atomic_tables_file",
    "test_case_config",
    "test_case_comment",
    "level_split_method",
)

INTEGRATION_RELATIVE_TOLERANCE = 1.0e-10
ATOMIC_MASS_UNIT = 1.660539040e-27
ELEMENTARY_CHARGE = 1.60217733e-19
V2_TO_ENERGY_PER_AMU = (
    ATOMIC_MASS_UNIT / (2.0 * ELEMENTARY_CHARGE * 1.0e3) * 1.0e-4
)
NUMBER_OF_ATOMIC_LEVELS = 6


@dataclass
class SinkData:
    """Validated data used to compare one ion-sink implementation."""

    filename: Path
    energy: np.ndarray
    pitch: np.ndarray
    sink_distribution: np.ndarray
    energy_marginal: np.ndarray
    pitch_marginal: np.ndarray
    total_rate: float
    standard_error: float
    neutral_energy: float
    neutral_velocity: np.ndarray
    shared_values: dict
    provenance: dict
    implementation_values: dict

    @property
    def neutral_pitch(self):
        """Return the injected neutral pitch on the ion velocity-space axes."""
        speed = float(np.linalg.norm(self.neutral_velocity))
        return float(self.neutral_velocity[2] / speed)


def _decode_text(value):
    """Return HDF5 byte strings as ordinary Python strings."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _read_scalar(dataset, filename):
    """Read exactly one scalar value from an HDF5 dataset."""
    values = np.asarray(dataset[()]).reshape(-1)
    if values.size != 1:
        raise ValueError(
            f"{filename}: {dataset.name} must contain exactly one value."
        )
    value = values[0]
    if isinstance(value, (bytes, str)):
        return _decode_text(value)
    return value.item() if hasattr(value, "item") else value


def _read_units(dataset):
    """Return a dataset's units attribute as text."""
    units = dataset.attrs.get("units", "")
    return _decode_text(units)


def _validate_grid(grid, name, filename):
    """Require a finite, uniformly increasing one-dimensional grid."""
    if grid.ndim != 1 or grid.size < 2:
        raise ValueError(
            f"{filename}: {name} must be a vector with at least two values."
        )
    if not np.all(np.isfinite(grid)):
        raise ValueError(f"{filename}: {name} contains nonfinite values.")

    spacing = np.diff(grid)
    if np.any(spacing <= 0.0):
        raise ValueError(f"{filename}: {name} must be strictly increasing.")
    scale = max(1.0, float(np.max(np.abs(grid))), float(spacing[0]))
    tolerance = 100.0 * np.finfo(float).eps * scale
    if not np.allclose(spacing, spacing[0], rtol=0.0, atol=tolerance):
        raise ValueError(f"{filename}: {name} must be uniformly spaced.")


def _is_finite_positive(value):
    """Return whether a scalar is finite and strictly positive."""
    return bool(np.isfinite(value) and value > 0.0)


def _validate_units(group, filename, expected_units_by_dataset):
    """Require the specified physical units from one HDF5 group."""
    for dataset_name, expected_units in expected_units_by_dataset.items():
        actual_units = _read_units(group[dataset_name])
        if actual_units != expected_units:
            raise ValueError(
                f"{filename}: {group[dataset_name].name} units are "
                f"'{actual_units}'; expected '{expected_units}'."
            )


def _validate_integrated_products(data):
    """Check that stored marginals and total rate reproduce the 2D sink."""
    denergy = float(data.energy[1] - data.energy[0])
    dpitch = float(data.pitch[1] - data.pitch[0])

    calculated_energy_marginal = (
        np.sum(data.sink_distribution, axis=1) * dpitch
    )
    calculated_pitch_marginal = (
        np.sum(data.sink_distribution, axis=0) * denergy
    )
    calculated_total_rate = float(
        np.sum(data.sink_distribution) * denergy * dpitch
    )

    if not np.allclose(
        data.energy_marginal,
        calculated_energy_marginal,
        rtol=INTEGRATION_RELATIVE_TOLERANCE,
        atol=0.0,
    ):
        raise ValueError(
            f"{data.filename}: stored energy_marginal does not reproduce "
            "sink_distribution."
        )
    if not np.allclose(
        data.pitch_marginal,
        calculated_pitch_marginal,
        rtol=INTEGRATION_RELATIVE_TOLERANCE,
        atol=0.0,
    ):
        raise ValueError(
            f"{data.filename}: stored pitch_marginal does not reproduce "
            "sink_distribution."
        )
    if not np.isclose(
        data.total_rate,
        calculated_total_rate,
        rtol=INTEGRATION_RELATIVE_TOLERANCE,
        atol=0.0,
    ):
        raise ValueError(
            f"{data.filename}: total_reaction_rate does not reproduce "
            "sink_distribution."
        )


def read_sink(filename, implementation):
    """Read one implementation's self-contained ion-sink data.

    Deterministic artifacts store their comparison datasets at the HDF5 root.
    Production Monte Carlo sink files retain the FIDASIM root schema and store
    the derived comparison contract under ``/test_004``. Both comparison
    contracts are written by Python in canonical ``(energy, pitch)`` order, so
    the two-dimensional datasets require no Fortran/HDF5 axis conversion here.
    """
    if implementation not in ("deterministic", "monte_carlo"):
        raise ValueError(f"Unknown implementation: {implementation}")

    with h5py.File(filename, "r") as h5file:
        if implementation == "deterministic":
            group = h5file
            implementation_units = DETERMINISTIC_DATASET_UNITS
        else:
            if "test_004" not in h5file:
                raise ValueError(f"{filename}: missing /test_004.")
            group = h5file["test_004"]
            implementation_units = MONTE_CARLO_DATASET_UNITS

        required = [*DATASET_UNITS, *implementation_units, "species"]
        missing = [name for name in required if name not in group]
        if missing:
            raise ValueError(
                f"{filename}: {group.name} is missing {', '.join(missing)}."
            )

        _validate_units(
            group,
            filename,
            {
                **DATASET_UNITS,
                **implementation_units,
            },
        )

        energy = np.asarray(group["energy"][:], dtype=float)
        pitch = np.asarray(group["pitch"][:], dtype=float)
        sink_distribution = np.asarray(
            group["sink_distribution"][:],
            dtype=float,
        )
        energy_marginal = np.asarray(
            group["energy_marginal"][:],
            dtype=float,
        )
        pitch_marginal = np.asarray(
            group["pitch_marginal"][:],
            dtype=float,
        )
        total_rate = float(_read_scalar(group["total_reaction_rate"], filename))
        standard_error = 0.0
        if implementation == "monte_carlo":
            standard_error = float(
                _read_scalar(group["standard_error"], filename)
            )

        shared_values = {
            name: np.asarray(group[name][()])
            for name in SHARED_NUMERICAL_DATASETS
        }
        shared_values["species"] = _read_scalar(group["species"], filename)

        missing_attributes = [
            name for name in PROVENANCE_ATTRIBUTES if name not in group.attrs
        ]
        if missing_attributes:
            raise ValueError(
                f"{filename}: {group.name} is missing attributes "
                f"{', '.join(missing_attributes)}."
            )
        provenance = {
            name: _decode_text(group.attrs[name])
            for name in PROVENANCE_ATTRIBUTES
        }

        if implementation == "deterministic":
            implementation_values = {
                "gyroangle": np.asarray(group["gyroangle"][:], dtype=float)
            }
        else:
            implementation_values = {
                name: _read_scalar(group[name], filename)
                for name in MONTE_CARLO_DATASET_UNITS
            }

    _validate_grid(energy, "energy", filename)
    _validate_grid(pitch, "pitch", filename)
    expected_shape = (energy.size, pitch.size)
    if sink_distribution.shape != expected_shape:
        raise ValueError(
            f"{filename}: sink_distribution has shape "
            f"{sink_distribution.shape}; expected {expected_shape}."
        )
    if energy_marginal.shape != energy.shape:
        raise ValueError(f"{filename}: energy_marginal has an invalid shape.")
    if pitch_marginal.shape != pitch.shape:
        raise ValueError(f"{filename}: pitch_marginal has an invalid shape.")

    for name, values in (
        ("sink_distribution", sink_distribution),
        ("energy_marginal", energy_marginal),
        ("pitch_marginal", pitch_marginal),
    ):
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{filename}: {name} contains nonfinite values.")
        if np.any(values < 0.0):
            raise ValueError(f"{filename}: {name} contains negative values.")
    if not _is_finite_positive(total_rate):
        raise ValueError(f"{filename}: total_reaction_rate must be positive.")
    if not np.isfinite(standard_error) or standard_error < 0.0:
        raise ValueError(f"{filename}: standard_error must be nonnegative.")

    if implementation == "deterministic":
        gyroangle = implementation_values["gyroangle"]
        if (
            gyroangle.ndim != 1
            or gyroangle.size < 1
            or not np.all(np.isfinite(gyroangle))
        ):
            raise ValueError(
                f"{filename}: gyroangle must be a nonempty finite vector."
            )
    else:
        n_markers = implementation_values["n_markers"]
        reservoir_size = implementation_values["reservoir_size"]
        seed = implementation_values["seed"]
        for name, value in (
            ("n_markers", n_markers),
            ("reservoir_size", reservoir_size),
            ("seed", seed),
        ):
            is_number = isinstance(
                value,
                (int, float, np.integer, np.floating),
            )
            if (
                isinstance(value, bool)
                or not is_number
                or not np.isfinite(value)
                or int(value) != value
                or value < 1
            ):
                raise ValueError(
                    f"{filename}: {name} must be a positive integer."
                )

        sample_standard_deviation = implementation_values[
            "sample_standard_deviation"
        ]
        if (
            not np.isfinite(sample_standard_deviation)
            or sample_standard_deviation < 0.0
        ):
            raise ValueError(
                f"{filename}: sample_standard_deviation must be nonnegative."
            )
        expected_standard_error = sample_standard_deviation / np.sqrt(
            n_markers
        )
        if not np.isclose(
            standard_error,
            expected_standard_error,
            rtol=INTEGRATION_RELATIVE_TOLERANCE,
            atol=0.0,
        ):
            raise ValueError(
                f"{filename}: standard_error does not equal "
                "sample_standard_deviation/sqrt(n_markers)."
            )
        neutral_density = implementation_values["neutral_density"]
        if not _is_finite_positive(neutral_density):
            raise ValueError(
                f"{filename}: neutral_density must be finite and positive."
            )

    neutral_velocity = np.asarray(shared_values["neutral_velocity"], dtype=float)
    if neutral_velocity.shape != (3,) or not np.all(
        np.isfinite(neutral_velocity)
    ):
        raise ValueError(f"{filename}: neutral_velocity must have shape (3,).")
    if float(np.linalg.norm(neutral_velocity)) <= 0.0:
        raise ValueError(f"{filename}: neutral_velocity must be nonzero.")

    f_array = np.asarray(shared_values["f_array"], dtype=float)
    if f_array.shape != expected_shape:
        raise ValueError(f"{filename}: f_array has an invalid shape.")
    if not np.all(np.isfinite(f_array)) or np.any(f_array < 0.0):
        raise ValueError(
            f"{filename}: f_array must contain finite, nonnegative values."
        )
    if not _is_finite_positive(float(shared_values["denf"])):
        raise ValueError(f"{filename}: denf must be finite and positive.")

    neutral_level_density = np.asarray(
        shared_values["neutral_level_density"],
        dtype=float,
    )
    if neutral_level_density.shape != (6,):
        raise ValueError(
            f"{filename}: neutral_level_density must have shape (6,)."
        )
    if (
        not np.all(np.isfinite(neutral_level_density))
        or np.any(neutral_level_density < 0.0)
        or np.sum(neutral_level_density) <= 0.0
    ):
        raise ValueError(
            f"{filename}: neutral_level_density must be finite, "
            "nonnegative, and have a positive sum."
        )

    data = SinkData(
        filename=Path(filename),
        energy=energy,
        pitch=pitch,
        sink_distribution=sink_distribution,
        energy_marginal=energy_marginal,
        pitch_marginal=pitch_marginal,
        total_rate=total_rate,
        standard_error=standard_error,
        neutral_energy=float(shared_values["neutral_energy"]),
        neutral_velocity=neutral_velocity,
        shared_values=shared_values,
        provenance=provenance,
        implementation_values=implementation_values,
    )
    _validate_integrated_products(data)
    return data


def validate_pair(deterministic, monte_carlo, pair):
    """Require identical physical inputs, grids, and provenance."""
    for name in (*SHARED_NUMERICAL_DATASETS, "species"):
        deterministic_value = deterministic.shared_values[name]
        monte_carlo_value = monte_carlo.shared_values[name]
        if not np.array_equal(deterministic_value, monte_carlo_value):
            raise ValueError(
                f"{pair.basename}: shared dataset {name} differs between "
                "the deterministic and Monte Carlo artifacts."
            )

    for name in PROVENANCE_ATTRIBUTES:
        if deterministic.provenance[name] != monte_carlo.provenance[name]:
            raise ValueError(
                f"{pair.basename}: provenance attribute {name} differs "
                "between the deterministic and Monte Carlo artifacts."
            )


def _expected_neutral_values(neutral_config, atomic_mass):
    """Construct current neutral values independently of either artifact."""
    if neutral_config["level_split_method"] == "ground-only":
        fractions = np.zeros(NUMBER_OF_ATOMIC_LEVELS)
        fractions[0] = 1.0
    else:
        levels = np.arange(NUMBER_OF_ATOMIC_LEVELS, dtype=float)
        fractions = np.exp(-neutral_config["level_decay"] * levels)
        fractions /= np.sum(fractions)

    speed = np.sqrt(
        neutral_config["energy"]
        / (V2_TO_ENERGY_PER_AMU * atomic_mass)
    )
    angle = np.deg2rad(neutral_config["injection_angle"])
    velocity = speed * np.array(
        [np.sin(angle), 0.0, np.cos(angle)]
    )
    return {
        "neutral_energy": neutral_config["energy"],
        "injection_angle": neutral_config["injection_angle"],
        "neutral_velocity": velocity,
        "neutral_level_density": neutral_config["density"] * fractions,
    }


def validate_current_inputs(
    deterministic,
    monte_carlo,
    pair,
    source_distribution,
    test_config,
):
    """Require both artifacts to represent the currently selected inputs.

    Pair validation alone only proves that two artifacts agree with each
    other. These checks prevent an old but mutually consistent pair from being
    compared after the unified configuration or source files have changed.
    """
    for data in (deterministic, monte_carlo):
        for name, expected_value in pair.expected_provenance.items():
            if data.provenance[name] != expected_value:
                raise ValueError(
                    f"{data.filename}: provenance attribute {name} does not "
                    "match the current Test 004 inputs."
                )

    expected_source_values = {
        "energy": source_distribution.energy,
        "pitch": source_distribution.pitch,
        "f_array": source_distribution.values,
        "denf": source_distribution.density,
        "species": source_distribution.species,
        "atomic_number": source_distribution.atomic_number,
        "mass_number": source_distribution.mass_number,
        "charge_state": source_distribution.charge_state,
        "A": source_distribution.atomic_mass,
        "selected_r": source_distribution.selected_r,
        "selected_z": source_distribution.selected_z,
    }
    for name, expected_value in expected_source_values.items():
        if not np.array_equal(
            deterministic.shared_values[name],
            expected_value,
        ):
            raise ValueError(
                f"{pair.basename}: shared dataset {name} does not match "
                "the current source distribution."
            )

    expected_neutral_values = _expected_neutral_values(
        neutral_config=test_config["neutrals"],
        atomic_mass=source_distribution.atomic_mass,
    )
    for name, expected_value in expected_neutral_values.items():
        if not np.array_equal(
            deterministic.shared_values[name],
            expected_value,
        ):
            raise ValueError(
                f"{pair.basename}: shared dataset {name} does not match "
                "the current neutral configuration."
            )

    n_gyro = test_config["deterministic"]["n_gyro"]
    expected_gyroangle = 2.0 * np.pi * (
        np.arange(n_gyro, dtype=float) + 0.5
    ) / n_gyro
    if not np.array_equal(
        deterministic.implementation_values["gyroangle"],
        expected_gyroangle,
    ):
        raise ValueError(
            f"{pair.deterministic}: gyroangle does not match the current "
            "n_gyro setting."
        )

    monte_carlo_config = test_config["monte_carlo"]
    for name in ("n_markers", "reservoir_size", "seed"):
        if monte_carlo.implementation_values[name] != monte_carlo_config[name]:
            raise ValueError(
                f"{pair.monte_carlo}: {name} does not match the current "
                "Monte Carlo configuration."
            )
    if (
        monte_carlo.implementation_values["neutral_density"]
        != test_config["neutrals"]["density"]
    ):
        raise ValueError(
            f"{pair.monte_carlo}: neutral_density does not match the current "
            "neutral configuration."
        )
