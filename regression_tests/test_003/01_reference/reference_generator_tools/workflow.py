"""Core workflow for generating per-location reference HDF5 datasets."""

from pathlib import Path

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from regression_test_tools import (
    ConfigError,
    as_list,
    normalize_path,
    normalize_string,
    read_namelist,
    require_choice,
    require_existing_file,
    require_integer,
    require_real,
    require_string,
)

from .config import read_config
from .readers import load_fidasim_h5_distribution


def _build_output_path(output_filename, index):
    """Insert a numeric index before the configured HDF5 extension."""
    output_path = Path(output_filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    indexed_name = f"{output_path.stem}_{index:03d}{output_path.suffix}"
    return output_path.parent / indexed_name


def _read_locations(value, field_label):
    """Return a scalar or sequence of configured locations as floats."""
    locations = []
    for entry_index, raw_value in enumerate(as_list(value=value), start=1):
        location = require_real(
            value=raw_value,
            field_label=f"{field_label} entry {entry_index}",
        )
        locations.append(location)
    return locations


def _discover_test_002_outputs(run_config):
    """Discover Test 002 outputs, locations, and particle metadata."""
    run_config_path, run_blocks = read_namelist(config_path=run_config)
    if "input" not in run_blocks or "save_data_block" not in run_blocks:
        raise ConfigError("The Test 002 Stage 2 configuration is missing blocks.")

    run_input = run_blocks["input"]
    run_save = run_blocks["save_data_block"]
    if "reference_config" not in run_input:
        raise ConfigError("The Test 002 input block is missing reference_config.")
    if "output_filename" not in run_save:
        raise ConfigError(
            "The Test 002 save_data_block is missing output_filename."
        )

    reference_config = normalize_path(
        value=run_input["reference_config"],
        config_path=run_config_path,
        field_label="Test 002 reference_config",
    )
    require_existing_file(
        path=reference_config,
        field_label="Test 002 reference_config",
    )
    output_base = normalize_path(
        value=run_save["output_filename"],
        config_path=run_config_path,
        field_label="Test 002 output_filename",
    )

    reference_config_path, reference_blocks = read_namelist(
        config_path=reference_config
    )
    if "input" not in reference_blocks:
        raise ConfigError("The Test 002 Stage 1 configuration has no input block.")
    reference_input = reference_blocks["input"]
    required_fields = [
        "species",
        "atomic_number",
        "mass_number",
        "charge_state",
        "r_locations",
        "z_locations",
    ]
    for field_name in required_fields:
        if field_name not in reference_input:
            raise ConfigError(
                f"The Test 002 Stage 1 input is missing {field_name}."
            )

    species_value = require_string(
        value=reference_input["species"], field_label="Test 002 species"
    )
    species = normalize_string(value=species_value)
    species = require_choice(
        value=species,
        supported_values=["h", "d", "t"],
        field_label="Test 002 species",
    )
    particle = {
        "species": species,
        "atomic_number": require_integer(
            value=reference_input["atomic_number"],
            field_label="Test 002 atomic_number",
        ),
        "mass_number": require_integer(
            value=reference_input["mass_number"],
            field_label="Test 002 mass_number",
        ),
        "charge_state": require_integer(
            value=reference_input["charge_state"],
            field_label="Test 002 charge_state",
        ),
    }
    if particle["atomic_number"] <= 0:
        raise ConfigError("Test 002 atomic_number must be greater than zero.")
    if particle["mass_number"] < particle["atomic_number"]:
        raise ConfigError(
            "Test 002 mass_number must be at least atomic_number."
        )
    if not 0 <= particle["charge_state"] <= particle["atomic_number"]:
        raise ConfigError(
            "Test 002 charge_state must be between zero and atomic_number."
        )

    requested_r = _read_locations(
        value=reference_input["r_locations"],
        field_label="Test 002 r_locations",
    )
    requested_z = _read_locations(
        value=reference_input["z_locations"],
        field_label="Test 002 z_locations",
    )
    if not requested_r or len(requested_r) != len(requested_z):
        raise ConfigError(
            "Test 002 r_locations and z_locations must be nonempty and equal length."
        )

    number_of_cases = len(requested_r)
    output_paths = []
    for case_index in range(1, number_of_cases + 1):
        output_path = _build_output_path(output_base, case_index)
        require_existing_file(
            path=output_path,
            field_label=f"Test 002 output {case_index:03d}",
        )
        output_paths.append(output_path)

    return {
        "output_paths": output_paths,
        "requested_r": requested_r,
        "requested_z": requested_z,
        "particle": particle,
        "run_config": run_config_path,
        "reference_config": reference_config_path,
    }


def _read_source_species(source_path):
    """Read and normalize the species attribute from a Test 002 output."""
    with h5py.File(source_path, mode="r") as h5file:
        if "species" not in h5file.attrs:
            raise ConfigError(f"Missing species attribute in {source_path}.")
        species = h5file.attrs["species"]
    if isinstance(species, bytes):
        species = species.decode("utf-8")
    return str(species).strip().lower()


def _validate_source_data(source_path, energy, pitch, distribution, density):
    """Validate the numerical data required by the Test 003 sampler."""
    if energy.size < 2 or pitch.size < 2:
        raise ConfigError(f"Energy and pitch grids are too short in {source_path}.")
    if distribution.shape != (energy.size, pitch.size):
        raise ConfigError(
            f"Distribution shape {distribution.shape} does not match "
            f"({energy.size}, {pitch.size}) in {source_path}."
        )
    if not np.all(np.isfinite(energy)) or not np.all(np.isfinite(pitch)):
        raise ConfigError(f"A coordinate grid is non-finite in {source_path}.")
    if not np.all(np.isfinite(distribution)):
        raise ConfigError(f"The distribution is non-finite in {source_path}.")
    if np.any(distribution < 0.0):
        raise ConfigError(f"The distribution is negative in {source_path}.")
    if not np.isfinite(density) or density <= 0.0:
        raise ConfigError(f"denf must be finite and positive in {source_path}.")


def _resolve_plot_limits(values, fmin, fmax):
    values = np.asarray(values, dtype=float)

    if fmin in (None, "", "auto", "AUTO"):
        fmin = np.nanmin(values)
    if fmax in (None, "", "auto", "AUTO"):
        fmax = np.nanmax(values)

    return fmin, fmax


def _plot_distribution(
    data,
    pitch,
    energy,
    selected_r,
    selected_z,
    output_path,
    config,
):
    plot_values = np.array(data, dtype=float)
    scale = config["scale"].lower()

    if scale == "log":
        positive = np.where(plot_values > 0.0, plot_values, np.nan)
        plot_values = np.log10(positive)

    vmin, vmax = _resolve_plot_limits(plot_values, config["fmin"], config["fmax"])

    fig, ax = plt.subplots(figsize=(6, 4))
    im = ax.imshow(
        plot_values,
        origin="lower",
        aspect="auto",
        extent=[pitch[0], pitch[-1], energy[0], energy[-1]],
        vmin=vmin,
        vmax=vmax,
        cmap=config["colormap"],
    )

    ax.set_xlabel("pitch")
    ax.set_ylabel("energy")
    ax.set_ylim(0.0, config["emax"])
    title = f"f(E, pitch) at R = {selected_r:.2f} cm, Z = {selected_z:.2f} cm"
    ax.set_title(title)
    if config["enable_colorbar"]:
        fig.colorbar(im, ax=ax, label="f")

    fig.tight_layout()
    fig.savefig(output_path.with_suffix(".png"))
    plt.close(fig)


def _write_output_file(
    output_path,
    energy,
    pitch,
    distribution,
    requested_r,
    requested_z,
    selected_r,
    selected_z,
    r_index,
    z_index,
    particle,
    input_path,
    source_density,
):
    """Write one selected distribution and its metadata to HDF5."""
    with h5py.File(output_path, "w") as h5f:
        h5f.create_dataset("energy_grid", data=energy)
        h5f.create_dataset("pitch_grid", data=pitch)
        h5f.create_dataset("f_array", data=distribution)
        h5f.create_dataset("denf", data=np.array([source_density], dtype=float))

        string_type = h5py.string_dtype(encoding="utf-8")
        h5f.create_dataset(
            "species", data=particle["species"], dtype=string_type
        )
        h5f.create_dataset("atomic_number", data=particle["atomic_number"])
        h5f.create_dataset("mass_number", data=particle["mass_number"])
        h5f.create_dataset("charge_state", data=particle["charge_state"])

        h5f.create_dataset("requested_r", data=np.array([requested_r], dtype=float))
        h5f.create_dataset("requested_z", data=np.array([requested_z], dtype=float))
        h5f.create_dataset("selected_r", data=np.array([selected_r], dtype=float))
        h5f.create_dataset("selected_z", data=np.array([selected_z], dtype=float))
        h5f.create_dataset("r_index", data=np.array([r_index], dtype=int))
        h5f.create_dataset("z_index", data=np.array([z_index], dtype=int))

        # Store units beside each dataset so the file remains self-describing.
        h5f["energy_grid"].attrs["units"] = "keV"
        h5f["pitch_grid"].attrs["units"] = "dimensionless"
        h5f["f_array"].attrs["units"] = "ions/(cm^3*keV*dP)"
        h5f["denf"].attrs["units"] = "ions/cm^3"
        h5f["species"].attrs["units"] = "n/a"
        h5f["atomic_number"].attrs["units"] = "dimensionless"
        h5f["mass_number"].attrs["units"] = "dimensionless"
        h5f["charge_state"].attrs["units"] = "elementary charge"

        position_datasets = (
            "requested_r",
            "requested_z",
            "selected_r",
            "selected_z",
        )
        for dataset_name in position_datasets:
            h5f[dataset_name].attrs["units"] = "cm"

        h5f["r_index"].attrs["units"] = "dimensionless"
        h5f["z_index"].attrs["units"] = "dimensionless"

        # Descriptions make the purpose of each dataset available in the file.
        h5f["energy_grid"].attrs["description"] = "Fast-ion energy grid"
        h5f["pitch_grid"].attrs["description"] = "Fast-ion pitch grid"
        h5f["f_array"].attrs["description"] = "Selected fast-ion distribution"
        h5f["denf"].attrs["description"] = (
            "Fast-ion density at the selected spatial grid point"
        )
        h5f["species"].attrs["description"] = "Normalized particle species label"
        h5f["atomic_number"].attrs["description"] = "Number of protons"
        h5f["mass_number"].attrs["description"] = (
            "Total number of protons and neutrons"
        )
        h5f["charge_state"].attrs["description"] = "Particle charge state"
        h5f["requested_r"].attrs["description"] = "Requested radial position"
        h5f["requested_z"].attrs["description"] = "Requested axial position"
        h5f["selected_r"].attrs["description"] = "Selected radial grid value"
        h5f["selected_z"].attrs["description"] = "Selected axial grid value"
        h5f["r_index"].attrs["description"] = "Selected radial grid index"
        h5f["z_index"].attrs["description"] = "Selected axial grid index"

        # Root attributes contain provenance for the output file as a whole.
        h5f.attrs["data_source_type"] = "test_002_outputs"
        h5f.attrs["data_source_name"] = str(input_path)
        h5f.attrs["description"] = "Extracted energy-pitch distribution slice"


def generate_outputs(config_path):
    """Generate compact sampler fixtures from one Test 002 output collection."""

    # Read the input configuration file:
    config_path = Path(config_path)
    config = read_config(
        config_filename=config_path,
    )

    # Extract configuration blocks:
    input_config = config["input"]
    plot_config = config["plot_data_block"]
    save_config = config["save_data_block"]

    output_filename = save_config["output_filename"]
    source = _discover_test_002_outputs(run_config=input_config["input_config"])
    generated_paths = []

    for output_index, source_path in enumerate(source["output_paths"], start=1):
        z, r, pitch, energy, f, denf = load_fidasim_h5_distribution(source_path)
        if r.size != 1 or z.size != 1:
            raise ConfigError(
                f"Expected one R and one Z location in {source_path}."
            )

        selected_r = float(r[0])
        selected_z = float(z[0])
        distribution = np.transpose(f[0, 0, :, :])
        source_density = float(denf[0, 0])
        _validate_source_data(
            source_path=source_path,
            energy=energy,
            pitch=pitch,
            distribution=distribution,
            density=source_density,
        )
        source_species = _read_source_species(source_path=source_path)
        if source_species != source["particle"]["species"]:
            raise ConfigError(
                f"Species in {source_path} is '{source_species}', expected "
                f"'{source['particle']['species']}'."
            )

        output_path = None
        if input_config["save_data"] or input_config["plot_data"]:
            output_path = _build_output_path(output_filename, output_index)

        if input_config["save_data"]:
            _write_output_file(
                output_path=output_path,
                energy=energy,
                pitch=pitch,
                distribution=distribution,
                requested_r=source["requested_r"][output_index - 1],
                requested_z=source["requested_z"][output_index - 1],
                selected_r=selected_r,
                selected_z=selected_z,
                r_index=0,
                z_index=0,
                particle=source["particle"],
                input_path=source_path,
                source_density=source_density,
            )
            generated_paths.append(output_path)

        if input_config["plot_data"]:
            _plot_distribution(
                data=distribution,
                pitch=pitch,
                energy=energy,
                selected_r=selected_r,
                selected_z=selected_z,
                output_path=output_path,
                config=plot_config,
            )

    return generated_paths
