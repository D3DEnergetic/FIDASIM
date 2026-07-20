"""Core workflow for generating per-location reference HDF5 datasets."""

from pathlib import Path

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from regression_test_tools import print_config

from .config import read_config
from .readers import load_input_distribution


def select_nearest_index(values, target):
    values = np.asarray(values)
    return int(np.argmin(np.abs(values - target)))


def _build_output_path(output_filename, index):
    """Insert a numeric index before the configured HDF5 extension."""
    output_path = Path(output_filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    indexed_name = f"{output_path.stem}_{index:03d}{output_path.suffix}"
    return output_path.parent / indexed_name


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
    input_config,
    input_path,
):
    """Write one selected distribution and its metadata to HDF5."""
    with h5py.File(output_path, "w") as h5f:
        h5f.create_dataset("energy_grid", data=energy)
        h5f.create_dataset("pitch_grid", data=pitch)
        h5f.create_dataset("f_array", data=distribution)

        string_type = h5py.string_dtype(encoding="utf-8")
        h5f.create_dataset(
            "species", data=input_config["species"], dtype=string_type
        )
        h5f.create_dataset("atomic_number", data=input_config["atomic_number"])
        h5f.create_dataset("mass_number", data=input_config["mass_number"])
        h5f.create_dataset("charge_state", data=input_config["charge_state"])

        h5f.create_dataset("requested_r", data=np.array([requested_r], dtype=float))
        h5f.create_dataset("requested_z", data=np.array([requested_z], dtype=float))
        h5f.create_dataset("selected_r", data=np.array([selected_r], dtype=float))
        h5f.create_dataset("selected_z", data=np.array([selected_z], dtype=float))
        h5f.create_dataset("r_index", data=np.array([r_index], dtype=int))
        h5f.create_dataset("z_index", data=np.array([z_index], dtype=int))

        # Store units beside each dataset so the file remains self-describing.
        h5f["energy_grid"].attrs["units"] = "keV"
        h5f["pitch_grid"].attrs["units"] = "dimensionless"
        h5f["f_array"].attrs["units"] = "fast-ions/(dE*dP*cm^3)"
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
        h5f.attrs["data_source_type"] = input_config["input_file_type"]
        h5f.attrs["data_source_name"] = str(input_path)
        h5f.attrs["description"] = "Extracted energy-pitch distribution slice"


def generate_outputs(config_path):

    # Read the input configuration file:
    config_path = Path(config_path)
    config = read_config(
        config_filename=config_path,
    )

    # Extract configuration blocks:
    input_config = config["input"]
    plot_config = config["plot_data_block"]
    save_config = config["save_data_block"]

    input_path = Path(input_config["input_filename"])
    output_filename = save_config["output_filename"]

    # The reader dispatcher uses the configured input type to select the
    # appropriate reader. That reader returns the common representation used
    # by the format-independent workflow below.
    z, r, pitch, energy, f = load_input_distribution(config)

    r_locations = input_config["r_locations"]
    z_locations = input_config["z_locations"]
    output_paths = []

    # Process corresponding R and Z entries as one requested sample location:
    num_locations = len(r_locations)
    for ii in range(num_locations):
        requested_r = r_locations[ii]
        requested_z = z_locations[ii]
        output_index = ii + 1

        z_index = select_nearest_index(z, requested_z)
        r_index = select_nearest_index(r, requested_r)

        selected_z = float(z[z_index])
        selected_r = float(r[r_index])

        f_slice = f[z_index, r_index, :, :]
        f_energy_pitch = np.transpose(f_slice)

        # Use the same indexed base name for this location's HDF5 and PNG files.
        if input_config["save_data"] or input_config["plot_data"]:
            output_path = _build_output_path(output_filename, output_index)

        if input_config["save_data"]:

            # Append the output path to the list of generated outputs for reporting back to the caller:
            output_paths.append(output_path)

            # Write the selected distribution and its metadata to an HDF5 file:
            _write_output_file(
                output_path,
                energy,
                pitch,
                f_energy_pitch,
                requested_r,
                requested_z,
                selected_r,
                selected_z,
                r_index,
                z_index,
                input_config,
                input_path,
            )

        # If requested, generate a PNG plot of the selected distribution:
        if input_config["plot_data"]:
            _plot_distribution(
                f_energy_pitch,
                pitch,
                energy,
                selected_r,
                selected_z,
                output_path,
                plot_config,
            )

    return output_paths
