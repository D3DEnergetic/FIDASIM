#!/usr/bin/env python3
"""Plot the sampled two-dimensional distributions produced by test_002."""

from pathlib import Path
import sys

# Make the shared regression-test tools importable without installing a package.
regression_tests_directory = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(regression_tests_directory))

from regression_test_tools import (
    ConfigError,
    normalize_path,
    normalize_string,
    print_config,
    read_namelist,
    reject_unknown_fields,
    require_boolean,
    require_blocks,
    require_choice,
    require_existing_directory,
    require_fields,
    require_real,
    require_string,
)

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# The Python plotter consumes only these fields from the Fortran-owned block.
RUN_TEST_CONSUMED_FIELDS = [
    "output_directory",
    "plot_data",
]

# The plotter owns the complete structure of plot_data_block.
PLOT_DATA_SCHEMA = {
    "required_fields": [],
    "optional_fields": [
        "scale",
        "fmin",
        "fmax",
        "enable_colorbar",
        "colormap",
    ],
}

SUPPORTED_SCALES = [
    "lin",
    "log",
]
SUPPORTED_COLORMAPS = [
    "viridis",
    "viridis_r",
    "hot",
    "hot_r",
]


def _normalize_plot_limit(value, field_label):
    """Normalize an automatic or numerical plot limit.

    Args:
        value (object): Plot limit read from the namelist. Accepted values are
        ``None``, the string ``"auto"``, or a real number.
        field_label (str): Human-readable field label used in error messages.

    Returns:
        float, str, or None: Canonical plot limit.

    Raises:
        ConfigError: If the value is not numerical, ``None``, or ``"auto"``.
    """
    if value is None:
        return None

    if isinstance(value, str):
        normalized_value = normalize_string(
            value=value,
        )
        if normalized_value == "auto":
            return normalized_value
        raise ConfigError(f"{field_label} must be a real number or 'auto'.")

    return require_real(
        value=value,
        field_label=field_label,
    )


def read_config(config_filename):
    """Read the canonical configuration consumed by the Stage 2 plotter.

    Args:
        config_filename (str or Path): Path to the Stage 2 namelist file.

    Returns:
        dict: Canonical configuration organized under ``run_test`` and
        ``plot_data_block`` keys that match the namelist block names.

    Raises:
        ConfigError: If a consumed field or plotting setting is invalid.
    """

    # Step 1: read the namelist and require the Fortran-owned run_test block.
    config_path, blocks = read_namelist(
        config_path=config_filename,
    )
    require_blocks(
        blocks=blocks,
        required_blocks=["run_test"],
    )

    run_block = blocks["run_test"]
    plot_block = blocks.get("plot_data_block", {})

    # Require only the run_test fields consumed by this Python plotter. Other
    # run_test fields belong to the Fortran sampler and are accepted unchanged.
    require_fields(
        block=run_block,
        required_fields=RUN_TEST_CONSUMED_FIELDS,
        block_label="&run_test",
    )

    # Validate every field in the plotter-owned plot_data_block.
    plot_required_fields = PLOT_DATA_SCHEMA["required_fields"]
    plot_optional_fields = PLOT_DATA_SCHEMA["optional_fields"]
    plot_allowed_fields = []
    for field_name in plot_required_fields:
        plot_allowed_fields.append(field_name)
    for field_name in plot_optional_fields:
        plot_allowed_fields.append(field_name)

    require_fields(
        block=plot_block,
        required_fields=plot_required_fields,
        block_label="&plot_data_block",
    )
    reject_unknown_fields(
        block=plot_block,
        allowed_fields=plot_allowed_fields,
        block_label="&plot_data_block",
    )

    # Step 2: normalize the consumed run_test settings.
    plot_data = require_boolean(
        value=run_block["plot_data"],
        field_label="plot_data",
    )
    output_directory_value = require_string(
        value=run_block["output_directory"],
        field_label="output_directory",
    )
    output_directory = normalize_path(
        value=output_directory_value,
        config_path=config_path,
        field_label="output_directory",
    )

    if plot_data:
        require_existing_directory(
            path=output_directory,
            field_label="output_directory",
        )

    # Step 3: normalize and validate the plotting settings.
    scale = normalize_string(
        value=plot_block.get("scale", "lin"),
    )
    scale = require_choice(
        value=scale,
        supported_values=SUPPORTED_SCALES,
        field_label="scale",
    )

    colormap = normalize_string(
        value=plot_block.get("colormap", "viridis"),
    )
    colormap = require_choice(
        value=colormap,
        supported_values=SUPPORTED_COLORMAPS,
        field_label="colormap",
    )

    enable_colorbar = require_boolean(
        value=plot_block.get("enable_colorbar", True),
        field_label="enable_colorbar",
    )
    fmin = _normalize_plot_limit(
        value=plot_block.get("fmin"),
        field_label="fmin",
    )
    fmax = _normalize_plot_limit(
        value=plot_block.get("fmax"),
        field_label="fmax",
    )

    # Step 4: return canonical blocks that preserve the namelist structure.
    return {
        "run_test": {
            "output_directory": output_directory,
            "plot_data": plot_data,
        },
        "plot_data_block": {
            "scale": scale,
            "fmin": fmin,
            "fmax": fmax,
            "enable_colorbar": enable_colorbar,
            "colormap": colormap,
        },
    }


def resolve_plot_limits(values, configured_minimum, configured_maximum):
    """Return explicit plotting limits, using the data for automatic limits."""
    if configured_minimum in (None, "auto"):
        configured_minimum = np.nanmin(values)
    if configured_maximum in (None, "auto"):
        configured_maximum = np.nanmax(values)
    return float(configured_minimum), float(configured_maximum)


def plot_sampled_file(filename, plot_config):
    """Read and plot one sampled energy-pitch distribution.

    Args:
        filename (Path): Sampled HDF5 file to plot.
        plot_config (dict): Canonical ``plot_data_block`` settings.

    Returns:
        None.
    """

    # Open the HDF5 file and read the data.
    with h5py.File(filename, "r") as h5file:
        energy = np.asarray(h5file["energy_grid"][:], dtype=float)
        pitch = np.asarray(h5file["pitch_grid"][:], dtype=float)
        distribution = np.asarray(h5file["f_array"][:], dtype=float)

        # Each sampled file represents one selected spatial location.
        selected_r = h5file["selected_r"][0]
        selected_z = h5file["selected_z"][0]
        units = h5file["f_array"].attrs.get("units", "f")

    # Validate that the distribution shape matches the energy-pitch grid.
    expected_shape = (energy.size, pitch.size)
    if distribution.shape != expected_shape:
        raise ValueError(
            f"{filename}: f_array shape {distribution.shape} does not match "
            f"energy-pitch shape {expected_shape}."
        )

    # Prepare the data for plotting.
    plot_values = distribution.copy()
    colorbar_label = str(units)
    if plot_config["scale"] == "log":
        positive_values = np.where(plot_values > 0.0, plot_values, np.nan)
        plot_values = np.log10(positive_values)
        colorbar_label = f"log10({colorbar_label})"

    # Determine the color scale limits for the plot.
    vmin, vmax = resolve_plot_limits(
        values=plot_values,
        configured_minimum=plot_config["fmin"],
        configured_maximum=plot_config["fmax"],
    )

    # Create the plot.
    figure, axes = plt.subplots(figsize=(6, 4))
    image = axes.imshow(
        plot_values,
        origin="lower",
        aspect="auto",
        extent=[pitch[0], pitch[-1], energy[0], energy[-1]],
        vmin=vmin,
        vmax=vmax,
        cmap=plot_config["colormap"],
    )
    axes.set_xlabel("pitch")
    axes.set_ylabel("energy [keV]")
    axes.set_title(
        f"Sampled f(E, pitch) at R = {selected_r:.2f} cm, "
        f"Z = {selected_z:.2f} cm"
    )
    if plot_config["enable_colorbar"]:
        figure.colorbar(image, ax=axes, label=colorbar_label)
    figure.tight_layout()

    # Save the plot to a PNG file.
    output_filename = filename.with_suffix(".png")
    figure.savefig(output_filename)
    plt.close(figure)
    print(f"Wrote plot: {output_filename}")


def main():
    # Check that the user provided a configuration file as a command-line argument.
    if len(sys.argv) != 2:
        raise SystemExit("Usage: plot_sampled_data.py <input_config.nml>")

    # Read the canonical configuration and report configuration errors clearly.
    try:
        config = read_config(
            config_filename=sys.argv[1],
        )
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None

    run_config = config["run_test"]
    plot_config = config["plot_data_block"]

    if not run_config["plot_data"]:
        print("Plotting disabled by plot_data = .false.")
        return

    # Check that sampled HDF5 files are available.
    output_directory = run_config["output_directory"]
    output_files = sorted(output_directory.glob("*.h5"))
    if not output_files:
        raise FileNotFoundError(
            f"No sampled HDF5 files found in {output_directory}."
        )

    # Plot each sampled HDF5 file.
    for filename in output_files:
        plot_sampled_file(
            filename=filename,
            plot_config=plot_config,
        )


if __name__ == "__main__":
    main()
