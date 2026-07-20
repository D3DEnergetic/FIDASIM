#!/usr/bin/env python3
"""Plot the sampled two-dimensional distributions produced by test_002."""

from pathlib import Path
import sys

import f90nml
import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SUPPORTED_SCALES = {"lin", "log"}
SUPPORTED_COLORMAPS = {"viridis", "viridis_r", "hot", "hot_r"}


def read_config(config_filename):

    # Read input configuration file:
    blocks = f90nml.read(config_filename)

    # Check that the &run_test block exists:
    if "run_test" not in blocks:
        raise KeyError("Missing run_test block in configuration file.")
    
    # Get blocks:
    run_config = blocks["run_test"]
    plot_config = blocks.get("plot_data_block", {})

    # Check if plotting is enabled:
    enabled = run_config.get("plot_data", False)
    if not isinstance(enabled, bool):
        raise ValueError("plot_data must be .true. or .false.")

    # Assemble configuration dictionary:
    config = {
        "plot_enabled": enabled,
        "output_directory": Path(run_config["output_directory"]),
        "scale": str(plot_config.get("scale", "lin")).strip().lower(),
        "fmin": plot_config.get("fmin"),
        "fmax": plot_config.get("fmax"),
        "enable_colorbar": plot_config.get("enable_colorbar", True),
        "colormap": str(plot_config.get("colormap", "viridis")).strip().lower(),
    }

    # Validate against supported values:
    if config["scale"] not in SUPPORTED_SCALES:
        raise ValueError("scale must be 'lin' or 'log'.")
    if config["colormap"] not in SUPPORTED_COLORMAPS:
        choices = ", ".join(sorted(SUPPORTED_COLORMAPS))
        raise ValueError(f"colormap must be one of: {choices}.")
    if not isinstance(config["enable_colorbar"], bool):
        raise ValueError("enable_colorbar must be .true. or .false.")

    return config

def resolve_plot_limits(values, configured_minimum, configured_maximum):
    if configured_minimum in (None, "", "auto", "AUTO"):
        configured_minimum = np.nanmin(values)
    if configured_maximum in (None, "", "auto", "AUTO"):
        configured_maximum = np.nanmax(values)
    return float(configured_minimum), float(configured_maximum)

def plot_sampled_file(filename, config):
    
    # Open the HDF5 file and read the data:
    with h5py.File(filename, "r") as h5file:

        # Extract data and normalize them to numpy arrays:
        energy = np.asarray(h5file["energy_grid"][:], dtype=float)
        pitch = np.asarray(h5file["pitch_grid"][:], dtype=float)
        distribution = np.asarray(h5file["f_array"][:], dtype=float)

        # Extract metadata for the plot title:
        selected_r = h5file["selected_r"][0]
        selected_z = h5file["selected_z"][0]
        units = h5file["f_array"].attrs.get("units", "f")

    # Validate that the distribution shape matches the energy-pitch grid:
    expected_shape = (energy.size, pitch.size)
    if distribution.shape != expected_shape:
        raise ValueError(
            f"{filename}: f_array shape {distribution.shape} does not match "
            f"energy-pitch shape {expected_shape}."
        )

    # Prepare the data for plotting:
    plot_values = distribution.copy()
    colorbar_label = str(units)
    if config["scale"] == "log":
        plot_values = np.log10(np.where(plot_values > 0.0, plot_values, np.nan))
        colorbar_label = f"log10({colorbar_label})"

    # Determine the color scale limits for the plot:
    vmin, vmax = resolve_plot_limits(
        plot_values, config["fmin"], config["fmax"]
    )
    
    # Create the plot:
    figure, axes = plt.subplots(figsize=(6, 4))
    image = axes.imshow(
        plot_values,
        origin="lower",
        aspect="auto",
        extent=[pitch[0], pitch[-1], energy[0], energy[-1]],
        vmin=vmin,
        vmax=vmax,
        cmap=config["colormap"],
    )
    axes.set_xlabel("pitch")
    axes.set_ylabel("energy [keV]")
    axes.set_title(
        f"Sampled f(E, pitch) at R = {selected_r:.2f} cm, "
        f"Z = {selected_z:.2f} cm"
    )
    if config["enable_colorbar"]:
        figure.colorbar(image, ax=axes, label=colorbar_label)
    figure.tight_layout()

    # Save the plot to a PNG file:
    figure.savefig(filename.with_suffix(".png"))
    plt.close(figure)

    # Print a message indicating that the plot was saved:
    print(f"Wrote plot: {filename.with_suffix('.png')}")

def main():
    # Check that the user provided a configuration file as a command-line argument:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: plot_sampled_data.py <input_config.nml>")

    # Check that the configuration file exists:
    input_arg = sys.argv[1]
    if not Path(input_arg).is_file():
        raise FileNotFoundError(f"Configuration file not found: {input_arg}")
    
    # Read the configuration file and check if plotting is enabled:
    config = read_config(input_arg)
    if not config["plot_enabled"]:
        print("Plotting disabled by plot_data = .false.")
        return

    # Check that the output files exist:
    output_files = sorted(config["output_directory"].glob("*.h5"))
    if not output_files:
        raise FileNotFoundError(
            f"No sampled HDF5 files found in {config['output_directory']}."
        )
    
    # Plot each sampled HDF5 file:
    for filename in output_files:
        plot_sampled_file(filename, config)


if __name__ == "__main__":
    main()
