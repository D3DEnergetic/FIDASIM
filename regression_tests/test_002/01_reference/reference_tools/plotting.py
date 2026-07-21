"""Create diagnostic plots of Test 002 reference distributions."""

from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from regression_test_tools import ConfigError


def _resolve_plot_limits(values, configured_minimum, configured_maximum):
    """Replace automatic plot limits with the finite data limits."""
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        raise ConfigError("The distribution contains no finite values to plot.")

    if configured_minimum is None or configured_minimum == "auto":
        plot_minimum = float(np.min(finite_values))
    else:
        plot_minimum = configured_minimum

    if configured_maximum is None or configured_maximum == "auto":
        plot_maximum = float(np.max(finite_values))
    else:
        plot_maximum = configured_maximum

    if plot_minimum >= plot_maximum:
        raise ConfigError("The plot minimum must be smaller than the plot maximum.")

    return plot_minimum, plot_maximum


def plot_reference_distribution(input_path, plot_config):
    """Plot one reference distribution in normalized parallel-perpendicular space.

    Args:
        input_path (str or Path): Reference HDF5 file containing the distribution.
        plot_config (dict): Canonical settings from ``plot_data_block``.

    Returns:
        Path: Path of the generated PNG file.
    """
    input_path = Path(input_path)

    with h5py.File(input_path, mode="r") as h5file:
        u_bar = np.asarray(h5file["u_bar"][:], dtype=float)
        theta = np.asarray(h5file["theta"][:], dtype=float)
        distribution = np.asarray(h5file["f_u_theta"][:], dtype=float)
        selected_r = float(h5file["selected_r"][()])
        selected_z = float(h5file["selected_z"][()])
        density = float(h5file["moments/density"][()])
        parallel_temperature = float(
            h5file["moments/parallel_temperature"][()]
        )
        perpendicular_temperature = float(
            h5file["moments/perpendicular_temperature"][()]
        )

    expected_shape = (theta.size, u_bar.size)
    if distribution.shape != expected_shape:
        raise ConfigError(
            f"f_u_theta has shape {distribution.shape}; expected {expected_shape}."
        )

    # Transform every (u_bar, theta) point into parallel-perpendicular space.
    u_bar_2d = u_bar[np.newaxis, :]
    theta_2d = theta[:, np.newaxis]
    u_bar_parallel = u_bar_2d * np.cos(theta_2d)
    u_bar_perpendicular = u_bar_2d * np.sin(theta_2d)

    plot_values = np.array(distribution, dtype=float)
    colorbar_label = r"$f(\bar{u}_{\parallel},\bar{u}_{\perp})$"
    if plot_config["scale"] == "log":
        positive_values = np.where(plot_values > 0.0, plot_values, np.nan)
        plot_values = np.log10(positive_values)
        colorbar_label = r"$\log_{10} f(\bar{u}_{\parallel},\bar{u}_{\perp})$"

    plot_minimum, plot_maximum = _resolve_plot_limits(
        values=plot_values,
        configured_minimum=plot_config["fmin"],
        configured_maximum=plot_config["fmax"],
    )
    contour_levels = np.linspace(
        plot_minimum,
        plot_maximum,
        plot_config["contour_levels"],
    )

    figure, axes = plt.subplots(figsize=(7, 5.5))
    contour = axes.contourf(
        u_bar_parallel,
        u_bar_perpendicular,
        plot_values,
        levels=contour_levels,
        cmap=plot_config["colormap"],
        extend="both",
    )

    axes.set_xlabel(r"Normalized parallel proper velocity, $\bar{u}_{\parallel}$")
    axes.set_ylabel(
        r"Normalized perpendicular proper velocity, $\bar{u}_{\perp}$"
    )
    axes.set_aspect("equal", adjustable="box")
    axes.set_title(
        f"CQL3D distribution at R = {selected_r:.2f} cm, "
        f"Z = {selected_z:.2f} cm"
    )

    moment_text = (
        f"n = {density:.2e} ions/cm³\n"
        f"T∥ = {parallel_temperature:.2e} keV\n"
        f"T⊥ = {perpendicular_temperature:.2e} keV"
    )
    axes.text(
        0.02,
        0.98,
        moment_text,
        transform=axes.transAxes,
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
