"""Write self-describing deterministic ion-sink reference files."""

from pathlib import Path

import h5py
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]


def _repository_relative(path):
    """Return stable repository-relative provenance when possible."""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _dataset(group, name, value, units, description):
    dataset = group.create_dataset(name, data=value)
    dataset.attrs["units"] = units
    dataset.attrs["description"] = description
    return dataset


def _text_dataset(group, name, value, description):
    """Create one scalar UTF-8 metadata dataset."""
    string_type = h5py.string_dtype(encoding="utf-8")
    dataset = group.create_dataset(name, data=value, dtype=string_type)
    dataset.attrs["description"] = description
    return dataset


def write_reference(
    filename,
    distribution,
    result,
    config,
    case,
    provenance,
):
    """Write one trusted Stage 1 HDF5 result."""
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5file:
        h5file.attrs["description"] = (
            "Deterministic direct charge-exchange ion-sink reference"
        )
        h5file.attrs["source_distribution"] = _repository_relative(
            case["input_path"]
        )
        h5file.attrs["input_distribution_config"] = _repository_relative(
            provenance["run_config"]
        )
        h5file.attrs["atomic_tables_file"] = _repository_relative(
            config["reference"]["tables_filename"]
        )
        h5file.attrs["level_split_method"] = config["neutrals"][
            "level_split_method"
        ]
        h5file.attrs["comment"] = config["reference"]["comment"]

        _dataset(h5file, "energy", distribution.energy, "keV", "Ion energy grid")
        _dataset(
            h5file, "pitch", distribution.pitch, "dimensionless", "Ion pitch grid"
        )
        _dataset(
            h5file,
            "f_array",
            distribution.values,
            "ions/(cm^3*keV*dP)",
            "Smooth Test 002 energy-pitch distribution",
        )
        _dataset(
            h5file,
            "denf",
            distribution.density,
            "ions/cm^3",
            "Authoritative ion density from Test 002",
        )
        _dataset(
            h5file,
            "sink_distribution",
            result.sink_distribution,
            "ions/(cm^3*s*keV*dP)",
            "Reaction-weighted energy-pitch ion-sink distribution",
        )
        _dataset(
            h5file,
            "cx_kernel",
            result.kernel,
            "1/s",
            "Gyrophase-averaged charge-exchange rate kernel",
        )
        _dataset(
            h5file,
            "energy_marginal",
            result.energy_marginal,
            "ions/(cm^3*s*keV)",
            "Sink distribution integrated over pitch",
        )
        _dataset(
            h5file,
            "pitch_marginal",
            result.pitch_marginal,
            "ions/(cm^3*s*dP)",
            "Sink distribution integrated over energy",
        )
        _dataset(
            h5file,
            "total_reaction_rate",
            result.total_rate,
            "ions/(cm^3*s)",
            "Total volumetric direct charge-exchange ion-sink rate",
        )
        _dataset(
            h5file,
            "gyroangle",
            result.gyroangle,
            "rad",
            "Deterministic midpoint gyrophase grid",
        )
        _dataset(
            h5file,
            "neutral_velocity",
            result.neutral_velocity,
            "cm/s",
            "Neutral velocity in Cartesian x-y-z coordinates",
        )
        _dataset(
            h5file,
            "neutral_level_density",
            result.level_density,
            "neutrals/cm^3",
            "Neutral density in atomic levels 1 through 6",
        )
        _dataset(
            h5file,
            "neutral_energy",
            config["neutrals"]["energy"],
            "keV",
            "Total neutral kinetic energy",
        )
        _dataset(
            h5file,
            "injection_angle",
            config["neutrals"]["injection_angle"],
            "degree",
            "Signed angle from +z toward +x",
        )
        _text_dataset(
            h5file,
            "species",
            distribution.species,
            "Canonical fast-ion species identifier",
        )
        _dataset(
            h5file,
            "atomic_number",
            distribution.atomic_number,
            "dimensionless",
            "Number of protons in the ion nucleus",
        )
        _dataset(
            h5file,
            "mass_number",
            distribution.mass_number,
            "dimensionless",
            "Integer isotope mass number",
        )
        _dataset(
            h5file,
            "charge_state",
            distribution.charge_state,
            "elementary charge",
            "Ion charge state",
        )
        _dataset(
            h5file,
            "A",
            distribution.atomic_mass,
            "amu",
            "Physical isotope mass",
        )
        _dataset(
            h5file,
            "selected_r",
            distribution.selected_r,
            "cm",
            "Selected source radial location",
        )
        _dataset(
            h5file,
            "selected_z",
            distribution.selected_z,
            "cm",
            "Selected source axial location",
        )
    return path
