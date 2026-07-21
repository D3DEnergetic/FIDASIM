"""Create diagnostic plots of converted FIDASIM distributions."""

from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from regression_test_tools import ConfigError


def _resolve_plot_limits(values, configured_minimum, configured_maximum):
    """Replace automatic limits with the minimum and maximum finite values."""
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        raise ConfigError("The distribution contains no finite values to plot.")

    if configured_minimum == "auto":
        plot_minimum = float(np.min(finite_values))
    else:
        plot_minimum = configured_minimum

    if configured_maximum == "auto":
        plot_maximum = float(np.max(finite_values))
    else:
        plot_maximum = configured_maximum

    if plot_minimum >= plot_maximum:
        raise ConfigError("The plot minimum must be smaller than the plot maximum.")

    return plot_minimum, plot_maximum


def plot_converted_distribution(input_path, plot_config):
    """Plot one saved FIDASIM distribution in energy-pitch coordinates.

    Args:
        input_path (str or Path): Converted single-location HDF5 file.
        plot_config (dict): Validated settings from ``plot_data_block``.

    Returns:
        Path: Path of the generated PNG file.
    """
    input_path = Path(input_path)

    with h5py.File(input_path, mode="r") as h5file:
        energy = np.asarray(h5file["energy"][:], dtype=float)
        pitch = np.asarray(h5file["pitch"][:], dtype=float)
        distribution = np.asarray(h5file["f"][0, 0, :, :], dtype=float).T
        selected_r = float(h5file["r"][0])
        selected_z = float(h5file["z"][0])
        density = float(h5file["moments/density"][()])
        parallel_temperature = float(
            h5file["moments/parallel_temperature"][()]
        )
        perpendicular_temperature = float(
            h5file["moments/perpendicular_temperature"][()]
        )

    expected_shape = (energy.size, pitch.size)
    if distribution.shape != expected_shape:
        raise ConfigError(
            f"f has shape {distribution.shape}; expected {expected_shape}."
        )

    plot_values = np.array(distribution, dtype=float)
    colorbar_label = r"$F(E,P)$ [ions/(cm$^3$ keV $dP$)]"
    if plot_config["scale"] == "log":
        positive_values = np.where(plot_values > 0.0, plot_values, np.nan)
        plot_values = np.log10(positive_values)
        colorbar_label = r"$\log_{10} F(E,P)$"

    plot_minimum, plot_maximum = _resolve_plot_limits(
        values=plot_values,
        configured_minimum=plot_config["fmin"],
        configured_maximum=plot_config["fmax"],
    )
    contour_values = np.linspace(
        plot_minimum,
        plot_maximum,
        plot_config["contour_levels"],
    )

    figure, axes = plt.subplots(figsize=(7, 5.5))
    contour = axes.contourf(
        energy,
        pitch,
        plot_values.T,
        levels=contour_values,
        cmap=plot_config["colormap"],
        extend="both",
    )
    axes.set_xlabel("Energy, E [keV]")
    axes.set_ylabel("Pitch, P")
    axes.set_xlim(0.0, plot_config["emax"])
    axes.set_title(
        f"FIDASIM distribution at R = {selected_r:.2f} cm, "
        f"Z = {selected_z:.2f} cm"
    )

    moment_text = (
        f"n = {density:.2e} ions/cm³\n"
        f"T∥ = {parallel_temperature:.2e} keV\n"
        f"T⊥ = {perpendicular_temperature:.2e} keV"
    )
    axes.text(
        0.98,
        0.98,
        moment_text,
        transform=axes.transAxes,
        horizontalalignment="right",
        verticalalignment="top",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )

    if plot_config["enable_colorbar"]:
        figure.colorbar(contour, ax=axes, label=colorbar_label)

    figure.tight_layout()
    output_path = input_path.with_suffix(".png")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path
