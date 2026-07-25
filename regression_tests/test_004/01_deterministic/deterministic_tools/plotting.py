"""Plot deterministic reaction-weighted ion-sink distributions."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
import numpy as np


class _TwoDecimalScalarFormatter(ScalarFormatter):
    """Keep scientific-notation colorbar mantissas at two decimal places."""

    def _set_format(self):
        self.format = "%.2f"


def _plot_limits(values, plot_config):
    plotted = np.asarray(values, dtype=float)
    if plot_config["scale"] == "log":
        plotted = np.where(plotted > 0.0, np.log10(plotted), np.nan)
    finite = plotted[np.isfinite(plotted)]
    if finite.size == 0:
        raise ValueError("No finite values are available for plotting.")
    minimum = plot_config["fmin"]
    maximum = plot_config["fmax"]
    if minimum in (None, "auto"):
        minimum = float(np.min(finite))
    elif plot_config["scale"] == "log":
        minimum = float(np.log10(minimum))
    if maximum in (None, "auto"):
        maximum = float(np.max(finite))
    elif plot_config["scale"] == "log":
        maximum = float(np.log10(maximum))
    if maximum <= minimum:
        raise ValueError("Plot fmax must be greater than fmin.")
    return plotted, minimum, maximum


def plot_deterministic(filename, distribution, result, plot_config):
    """Write the principal 2D reaction-weighted sink plot."""
    maximum_energy = plot_config["emax"]
    if maximum_energy in (None, "auto"):
        maximum_energy = float(distribution.energy[-1])
    if maximum_energy <= distribution.energy[0]:
        raise ValueError(
            "Plot emax must be greater than the minimum distribution energy."
        )

    values, minimum, maximum = _plot_limits(
        result.sink_distribution, plot_config
    )
    levels = np.linspace(minimum, maximum, plot_config["contour_levels"])

    figure, axis = plt.subplots(figsize=(8.0, 5.5), constrained_layout=True)
    contour = axis.contourf(
        distribution.pitch,
        distribution.energy,
        values,
        levels=levels,
        cmap=plot_config["colormap"],
        extend="both",
    )
    axis.set_xlabel("Pitch")
    axis.set_ylabel("Energy [keV]")
    axis.set_ylim(float(distribution.energy[0]), maximum_energy)
    axis.scatter(
        result.neutral_pitch,
        result.neutral_energy,
        s=90,
        marker="o",
        facecolor="limegreen",
        edgecolor="black",
        linewidth=1.0,
        label="Injected neutral",
        zorder=5,
    )
    axis.legend(loc="upper right", framealpha=0.9)
    axis.set_title(
        "Deterministic ion sink\n"
        f"R = {result.total_rate:.2e} ions cm$^{{-3}}$ s$^{{-1}}$"
    )
    if plot_config["enable_colorbar"]:
        label = r"$S(E,p)$ [ions cm$^{-3}$ s$^{-1}$ keV$^{-1}$]"
        tick_format = _TwoDecimalScalarFormatter()
        tick_format.set_powerlimits((0, 0))
        if plot_config["scale"] == "log":
            label = r"$\log_{10} S(E,p)$"
            tick_format = FormatStrFormatter("%.2f")
        figure.colorbar(
            contour,
            ax=axis,
            label=label,
            format=tick_format,
        )

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path
