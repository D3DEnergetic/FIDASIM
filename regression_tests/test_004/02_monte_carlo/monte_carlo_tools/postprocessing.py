"""Convert production sink particles into deterministic-compatible products."""

from pathlib import Path

import h5py
import numpy as np

from regression_test_tools import ConfigError
from test_004_tools import (
    discover_distributions,
    read_distribution,
    read_monte_carlo_config,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
MAXIMUM_RELATIVE_RATE_ERROR = 1.0e-12
ATOMIC_MASS_UNIT = 1.660539040e-27
ELEMENTARY_CHARGE = 1.60217733e-19
V2_TO_ENERGY_PER_AMU = (
    ATOMIC_MASS_UNIT / (2.0 * ELEMENTARY_CHARGE * 1.0e3) * 1.0e-4
)


def _repository_relative(path):
    """Return stable repository-relative provenance when possible."""
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(REPOSITORY_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _dataset(group, name, value, units, description):
    """Create one documented numerical dataset."""
    dataset = group.create_dataset(name, data=value)
    dataset.attrs["units"] = units
    dataset.attrs["description"] = description
    return dataset


def _text_dataset(group, name, value, description):
    """Create one documented scalar UTF-8 dataset."""
    string_type = h5py.string_dtype(encoding="utf-8")
    dataset = group.create_dataset(name, data=value, dtype=string_type)
    dataset.attrs["description"] = description
    return dataset


def _grid_edges(centers):
    """Return bin edges whose interior bins are centered on the input grid."""
    spacing = np.diff(centers)
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * spacing[0]
    edges[-1] = centers[-1] + 0.5 * spacing[-1]
    return edges


def _cell_volume(h5file, filename):
    """Return the uniform Cartesian beam-cell volume stored in a sink file."""
    widths = []
    for coordinate in ("x", "y", "z"):
        dataset_name = f"grid/{coordinate}"
        if dataset_name not in h5file:
            raise ConfigError(f"{filename}: missing /{dataset_name}.")
        centers = np.asarray(h5file[dataset_name][:], dtype=float)
        if centers.size < 2 or not np.all(np.isfinite(centers)):
            raise ConfigError(
                f"{filename}: /{dataset_name} must contain a finite grid."
            )
        spacing = np.diff(centers)
        if np.any(spacing <= 0.0) or not np.allclose(spacing, spacing[0]):
            raise ConfigError(
                f"{filename}: /{dataset_name} must be uniformly increasing."
            )
        widths.append(float(spacing[0]))
    return float(np.prod(widths))


def _read_sink_particles(h5file, filename, expected_markers):
    """Read and validate the production datasets used by postprocessing."""
    required = ["n_sink", "dens", "energy", "pitch", "weight", "ind"]
    missing = [name for name in required if name not in h5file]
    if missing:
        raise ConfigError(f"{filename} is missing: {', '.join(missing)}")

    number_of_particles = int(np.asarray(h5file["n_sink"][()]).reshape(-1)[0])
    if number_of_particles != expected_markers:
        raise ConfigError(
            f"{filename}: n_sink={number_of_particles}, expected "
            f"{expected_markers}."
        )

    energy = np.asarray(h5file["energy"][:], dtype=float)
    pitch = np.asarray(h5file["pitch"][:], dtype=float)
    weight = np.asarray(h5file["weight"][:], dtype=float)

    # write_sink_profile uses the Fortran layouts ind(component, particle) and
    # dens(component, x, y, z). h5py exposes these as (particle, component) and
    # (z, y, x, component). The particle indices are already in the desired
    # Python order; density is only reduced over every axis, so no reorder is
    # required for either dataset.
    indices = np.asarray(h5file["ind"][:], dtype=int)
    density = np.asarray(h5file["dens"][:], dtype=float)

    for name, values in (
        ("energy", energy),
        ("pitch", pitch),
        ("weight", weight),
        ("dens", density),
    ):
        if not np.all(np.isfinite(values)):
            raise ConfigError(f"{filename}: /{name} contains nonfinite values.")
    if energy.shape != (number_of_particles,):
        raise ConfigError(f"{filename}: /energy has an incorrect shape.")
    if pitch.shape != energy.shape or weight.shape != energy.shape:
        raise ConfigError(f"{filename}: particle arrays have inconsistent shapes.")
    if indices.shape != (number_of_particles, 3):
        raise ConfigError(f"{filename}: /ind must have shape (n_sink,3).")
    if np.any(indices != np.array([2, 2, 2])):
        raise ConfigError(f"{filename}: sink particles are not all in cell [2,2,2].")
    if np.any(weight <= 0.0):
        raise ConfigError(f"{filename}: /weight must be strictly positive.")
    if np.count_nonzero(density) != 1:
        raise ConfigError(
            f"{filename}: /dens must contain exactly one nonzero cell."
        )

    return {
        "energy": energy,
        "pitch": pitch,
        "weight": weight,
        "density": density,
    }


def _calculate_products(distribution, particles, cell_volume, filename):
    """Return weighted energy-pitch products and total-rate uncertainty."""
    energy_edges = _grid_edges(distribution.energy)
    pitch_edges = _grid_edges(distribution.pitch)
    weighted_bins, _, _ = np.histogram2d(
        particles["energy"],
        particles["pitch"],
        bins=(energy_edges, pitch_edges),
        weights=particles["weight"],
    )

    particle_weight = float(np.sum(particles["weight"]))
    binned_weight = float(np.sum(weighted_bins))
    relative_binning_error = abs(binned_weight - particle_weight) / particle_weight
    if relative_binning_error > MAXIMUM_RELATIVE_RATE_ERROR:
        raise ConfigError(
            f"{filename}: energy-pitch bins do not contain all particle weight."
        )

    energy_width = np.diff(energy_edges)
    pitch_width = np.diff(pitch_edges)
    bin_area = energy_width[:, None] * pitch_width[None, :]
    sink_distribution = weighted_bins / (cell_volume * bin_area)
    energy_marginal = np.sum(
        sink_distribution * pitch_width[None, :],
        axis=1,
    )
    pitch_marginal = np.sum(
        sink_distribution * energy_width[:, None],
        axis=0,
    )
    total_rate = particle_weight / cell_volume

    marker_rates = (
        particles["weight"] * particles["weight"].size / cell_volume
    )
    if marker_rates.size > 1:
        sample_standard_deviation = float(np.std(marker_rates, ddof=1))
        standard_error = sample_standard_deviation / np.sqrt(marker_rates.size)
    else:
        sample_standard_deviation = 0.0
        standard_error = 0.0

    profile_rate = float(np.sum(particles["density"]))
    relative_profile_error = abs(profile_rate - total_rate) / total_rate
    if relative_profile_error > MAXIMUM_RELATIVE_RATE_ERROR:
        raise ConfigError(
            f"{filename}: /dens and particle weights give different rates."
        )

    return {
        "sink_distribution": sink_distribution,
        "energy_marginal": energy_marginal,
        "pitch_marginal": pitch_marginal,
        "total_rate": total_rate,
        "sample_standard_deviation": sample_standard_deviation,
        "standard_error": standard_error,
    }


def _neutral_parameters(config, distribution):
    """Return the configured neutral velocity and six-level density."""
    neutral = config["neutrals"]
    speed = np.sqrt(
        neutral["energy"] / (V2_TO_ENERGY_PER_AMU * distribution.atomic_mass)
    )
    angle = np.deg2rad(neutral["injection_angle"])
    velocity = np.array(
        [speed * np.sin(angle), 0.0, speed * np.cos(angle)]
    )

    if neutral["level_split_method"] == "ground-only":
        fractions = np.zeros(6)
        fractions[0] = 1.0
    else:
        levels = np.arange(6, dtype=float)
        fractions = np.exp(-neutral["level_decay"] * levels)
        fractions /= np.sum(fractions)
    return velocity, neutral["density"] * fractions


def _write_products(
    h5file,
    distribution,
    result,
    config,
    case,
    provenance,
    cell_volume,
):
    """Replace the derived /test_004 group in one production sink file."""
    neutral_velocity, neutral_level_density = _neutral_parameters(
        config,
        distribution,
    )
    if "test_004" in h5file:
        del h5file["test_004"]
    group = h5file.create_group("test_004")
    group.attrs["description"] = (
        "Test 004 Monte Carlo energy-pitch ion-sink products"
    )
    group.attrs["source_distribution"] = _repository_relative(case["input_path"])
    group.attrs["input_distribution_config"] = _repository_relative(
        provenance["run_config"]
    )
    group.attrs["atomic_tables_file"] = _repository_relative(
        config["test_case"]["tables_filename"]
    )
    group.attrs["test_case_config"] = _repository_relative(
        config["test_case"]["config_path"]
    )
    group.attrs["test_case_comment"] = config["test_case"]["comment"]
    group.attrs["implementation_comment"] = config["monte_carlo"]["comment"]
    group.attrs["level_split_method"] = config["neutrals"]["level_split_method"]

    _dataset(group, "energy", distribution.energy, "keV", "Ion energy grid")
    _dataset(
        group, "pitch", distribution.pitch, "dimensionless", "Ion pitch grid"
    )
    _dataset(
        group,
        "f_array",
        distribution.values,
        "ions/(cm^3*keV*dP)",
        "Smooth Test 002 energy-pitch distribution",
    )
    _dataset(
        group,
        "denf",
        distribution.density,
        "ions/cm^3",
        "Authoritative ion density from Test 002",
    )
    _dataset(
        group,
        "sink_distribution",
        result["sink_distribution"],
        "ions/(cm^3*s*keV*dP)",
        "Weighted Monte Carlo energy-pitch ion-sink distribution",
    )
    _dataset(
        group,
        "energy_marginal",
        result["energy_marginal"],
        "ions/(cm^3*s*keV)",
        "Sink distribution integrated over pitch",
    )
    _dataset(
        group,
        "pitch_marginal",
        result["pitch_marginal"],
        "ions/(cm^3*s*dP)",
        "Sink distribution integrated over energy",
    )
    _dataset(
        group,
        "total_reaction_rate",
        result["total_rate"],
        "ions/(cm^3*s)",
        "Total volumetric Monte Carlo ion-sink rate",
    )
    _dataset(
        group,
        "sample_standard_deviation",
        result["sample_standard_deviation"],
        "ions/(cm^3*s)",
        "Sample standard deviation of the unscaled marker rate",
    )
    _dataset(
        group,
        "standard_error",
        result["standard_error"],
        "ions/(cm^3*s)",
        "Monte Carlo standard error of the total reaction-rate estimate",
    )
    _dataset(
        group,
        "cell_volume",
        cell_volume,
        "cm^3",
        "Central Cartesian beam-cell volume",
    )
    _dataset(
        group,
        "n_markers",
        config["monte_carlo"]["n_markers"],
        "dimensionless",
        "Number of sampled ion markers",
    )
    _dataset(
        group,
        "reservoir_size",
        config["monte_carlo"]["reservoir_size"],
        "dimensionless",
        "Number of markers in the central-cell neutral reservoir",
    )
    _dataset(
        group,
        "seed",
        config["monte_carlo"]["seed"],
        "dimensionless",
        "Serial random-number seed",
    )
    _dataset(
        group,
        "neutral_density",
        config["neutrals"]["density"],
        "neutrals/cm^3",
        "Configured total neutral density",
    )
    _dataset(
        group,
        "neutral_velocity",
        neutral_velocity,
        "cm/s",
        "Neutral velocity in Cartesian x-y-z coordinates",
    )
    _dataset(
        group,
        "neutral_level_density",
        neutral_level_density,
        "neutrals/cm^3",
        "Neutral density in atomic levels 1 through 6",
    )
    _dataset(
        group,
        "neutral_energy",
        config["neutrals"]["energy"],
        "keV",
        "Configured neutral kinetic energy",
    )
    _dataset(
        group,
        "injection_angle",
        config["neutrals"]["injection_angle"],
        "degree",
        "Signed neutral injection angle from +z toward +x",
    )
    _text_dataset(
        group,
        "species",
        distribution.species,
        "Canonical fast-ion species identifier",
    )
    _dataset(
        group,
        "atomic_number",
        distribution.atomic_number,
        "dimensionless",
        "Number of protons in the ion nucleus",
    )
    _dataset(
        group,
        "mass_number",
        distribution.mass_number,
        "dimensionless",
        "Integer isotope mass number",
    )
    _dataset(
        group,
        "charge_state",
        distribution.charge_state,
        "elementary charge",
        "Ion charge state",
    )
    _dataset(group, "A", distribution.atomic_mass, "amu", "Physical isotope mass")
    _dataset(
        group,
        "selected_r",
        distribution.selected_r,
        "cm",
        "Selected source radial location",
    )
    _dataset(
        group,
        "selected_z",
        distribution.selected_z,
        "cm",
        "Selected source axial location",
    )


def postprocess_monte_carlo(config_filename):
    """Append derived products and optionally plot every production sink file."""
    config = read_monte_carlo_config(config_filename)
    if not config["monte_carlo"]["save_data"]:
        print("Monte Carlo save_data is disabled; no sink files to postprocess.")
        return []

    provenance = discover_distributions(
        config["test_case"]["input_distribution_config"]
    )
    output_directory = config["save_data_block"]["monte_carlo_directory"]
    case_files = [
        (
            case,
            output_directory
            / f"{Path(case['input_path']).stem}_ion_sink.h5",
        )
        for case in provenance["cases"]
    ]
    missing_files = [
        str(filename) for _, filename in case_files if not filename.is_file()
    ]
    if missing_files:
        raise ConfigError(
            "Missing production sink files: " + ", ".join(missing_files)
        )

    results = []

    for case, filename in case_files:
        distribution = read_distribution(case["input_path"])

        with h5py.File(filename, "r+") as h5file:
            cell_volume = _cell_volume(h5file, filename)
            particles = _read_sink_particles(
                h5file,
                filename,
                expected_markers=config["monte_carlo"]["n_markers"],
            )
            result = _calculate_products(
                distribution,
                particles,
                cell_volume,
                filename,
            )
            _write_products(
                h5file,
                distribution,
                result,
                config,
                case,
                provenance,
                cell_volume,
            )

        print(
            f"{case['index']:03d}: {filename.name} "
            f"R={result['total_rate']:.10e} "
            f"+/- {result['standard_error']:.3e} ions/(cm^3*s)"
        )
        plot_filename = None
        if config["monte_carlo"]["plot_data"]:
            # Delay importing Matplotlib until plotting has been requested.
            from .plotting import plot_monte_carlo

            plot_filename = filename.with_suffix(".png")
            plot_monte_carlo(
                filename,
                plot_filename,
                config["plot_data_block"],
            )
            print(f"     plot: {plot_filename.name}")

        results.append(
            {
                "case": case,
                "filename": filename,
                "plot_filename": plot_filename,
                "result": result,
            }
        )
    return results
