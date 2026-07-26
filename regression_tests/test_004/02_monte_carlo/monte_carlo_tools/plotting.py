"""Plot postprocessed Monte Carlo ion-sink distributions."""

from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
import numpy as np


class _TwoDecimalScalarFormatter(ScalarFormatter):
    """Keep scientific-notation colorbar mantissas at two decimal places."""

    def _set_format(self):
        self.format = "%.2f"


def _read_plot_data(filename):
    """Read the self-contained Test 004 datasets needed by the plotter."""
    with h5py.File(filename, "r") as h5file:
        if "test_004" not in h5file:
            raise ValueError(f"{filename}: missing /test_004.")
        group = h5file["test_004"]
        required = [
            "energy",
            "pitch",
            "sink_distribution",
            "total_reaction_rate",
            "standard_error",
            "neutral_energy",
            "neutral_velocity",
        ]
        missing = [name for name in required if name not in group]
        if missing:
            raise ValueError(
                f"{filename}: /test_004 is missing {', '.join(missing)}."
            )

        data = {
            "energy": np.asarray(group["energy"][:], dtype=float),
            "pitch": np.asarray(group["pitch"][:], dtype=float),
            "sink_distribution": np.asarray(
                group["sink_distribution"][:],
                dtype=float,
            ),
            "total_rate": float(group["total_reaction_rate"][()]),
            "standard_error": float(group["standard_error"][()]),
            "neutral_energy": float(group["neutral_energy"][()]),
            "neutral_velocity": np.asarray(
                group["neutral_velocity"][:],
                dtype=float,
            ),
        }

    expected_shape = (data["energy"].size, data["pitch"].size)
    if data["sink_distribution"].shape != expected_shape:
        raise ValueError(
            f"{filename}: /test_004/sink_distribution has shape "
            f"{data['sink_distribution'].shape}; expected {expected_shape}."
        )
    neutral_speed = float(np.linalg.norm(data["neutral_velocity"]))
    if neutral_speed <= 0.0:
        raise ValueError(f"{filename}: neutral speed must be positive.")
    data["neutral_pitch"] = data["neutral_velocity"][2] / neutral_speed
    return data


def _plot_limits(values, plot_config):
    """Apply the requested scale and return finite contour limits."""
    plotted = np.asarray(values, dtype=float)
    if plot_config["scale"] == "log":
        positive = plotted > 0.0
        logarithm = np.full(plotted.shape, np.nan)
        logarithm[positive] = np.log10(plotted[positive])
        plotted = logarithm

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


def plot_monte_carlo(filename, output_filename, plot_config):
    """Write the Monte Carlo energy-pitch ion-sink plot."""
    data = _read_plot_data(filename)
    maximum_energy = plot_config["emax"]
    if maximum_energy in (None, "auto"):
        maximum_energy = float(data["energy"][-1])
    if maximum_energy <= data["energy"][0]:
        raise ValueError(
            "Plot emax must be greater than the minimum distribution energy."
        )

    values, minimum, maximum = _plot_limits(
        data["sink_distribution"],
        plot_config,
    )
    levels = np.linspace(minimum, maximum, plot_config["contour_levels"])

    figure, axis = plt.subplots(figsize=(8.0, 5.5), constrained_layout=True)
    contour = axis.contourf(
        data["pitch"],
        data["energy"],
        values,
        levels=levels,
        cmap=plot_config["colormap"],
        extend="both",
    )
    axis.set_xlabel("Pitch")
    axis.set_ylabel("Energy [keV]")
    axis.set_ylim(float(data["energy"][0]), maximum_energy)
    axis.scatter(
        data["neutral_pitch"],
        data["neutral_energy"],
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
        "Monte Carlo ion sink\n"
        f"R = {data['total_rate']:.2e} "
        f"+/- {data['standard_error']:.2e}"
        " ions cm$^{-3}$ s$^{-1}$"
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

    path = Path(output_filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path
