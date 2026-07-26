"""Create Test 004 ion-sink comparison figures."""

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


def _maximum_energy(energy, configured_maximum):
    """Return the requested upper energy extent after validating it."""
    if configured_maximum in (None, "auto"):
        return float(energy[-1])
    if configured_maximum <= energy[0]:
        raise ValueError(
            "Plot emax must be greater than the minimum distribution energy."
        )
    return float(configured_maximum)


def _transform_distribution(values, scale):
    """Return linear values or log10 values with empty bins masked."""
    plotted = np.asarray(values, dtype=float)
    if scale == "lin":
        return plotted.copy()

    transformed = np.full(plotted.shape, np.nan)
    positive = plotted > 0.0
    transformed[positive] = np.log10(plotted[positive])
    return transformed


def _shared_plot_limits(deterministic, monte_carlo, plot_config):
    """Return one color range shared by both two-dimensional panels."""
    deterministic_values = _transform_distribution(
        deterministic,
        plot_config["scale"],
    )
    monte_carlo_values = _transform_distribution(
        monte_carlo,
        plot_config["scale"],
    )
    combined = np.concatenate(
        (deterministic_values.ravel(), monte_carlo_values.ravel())
    )
    finite = combined[np.isfinite(combined)]
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

    return (
        deterministic_values,
        monte_carlo_values,
        float(minimum),
        float(maximum),
    )


def _mark_neutral(axis, data):
    """Mark the injected neutral on one energy-pitch panel."""
    axis.scatter(
        data.neutral_pitch,
        data.neutral_energy,
        s=75,
        marker="o",
        facecolor="limegreen",
        edgecolor="black",
        linewidth=1.0,
        label="Injected neutral",
        zorder=5,
    )


def _output_path(output_directory, filename):
    """Create the requested output directory and return one child path."""
    directory = Path(output_directory)
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename


def plot_distributions(
    pair,
    deterministic,
    monte_carlo,
    result,
    output_directory,
    plot_config,
):
    """Plot deterministic and Monte Carlo 2D sinks on one color scale."""
    (
        deterministic_values,
        monte_carlo_values,
        minimum,
        maximum,
    ) = _shared_plot_limits(
        deterministic.sink_distribution,
        monte_carlo.sink_distribution,
        plot_config,
    )
    levels = np.linspace(
        minimum,
        maximum,
        plot_config["contour_levels"],
    )
    maximum_energy = _maximum_energy(
        deterministic.energy,
        plot_config["emax"],
    )

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12.0, 5.2),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    contour = axes[0].contourf(
        deterministic.pitch,
        deterministic.energy,
        deterministic_values,
        levels=levels,
        cmap=plot_config["colormap"],
        extend="both",
    )
    axes[1].contourf(
        monte_carlo.pitch,
        monte_carlo.energy,
        monte_carlo_values,
        levels=levels,
        cmap=plot_config["colormap"],
        extend="both",
    )

    axes[0].set_title(
        "Deterministic\n"
        f"R = {result.deterministic_rate:.2e} ions cm$^{{-3}}$ s$^{{-1}}$"
    )
    axes[1].set_title(
        "Monte Carlo\n"
        f"R = {result.monte_carlo_rate:.2e} "
        f"+/- {result.monte_carlo_standard_error:.2e} "
        "ions cm$^{-3}$ s$^{-1}$"
    )
    for axis, data in zip(axes, (deterministic, monte_carlo)):
        axis.set_xlabel("Pitch")
        axis.set_ylim(float(deterministic.energy[0]), maximum_energy)
        _mark_neutral(axis, data)
    axes[0].set_ylabel("Energy [keV]")
    axes[1].legend(loc="upper right", framealpha=0.9)

    status = "PASS" if result.passed else "FAIL"
    figure.suptitle(
        f"{pair.basename}: ion-sink distributions [{status}]\n"
        f"relative rate difference = "
        f"{100.0 * result.signed_relative_rate_difference:.2f}%, "
        f"standardized difference = "
        f"{result.standardized_rate_difference:.2f}"
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
            ax=axes,
            label=label,
            format=tick_format,
            shrink=0.9,
        )

    output_filename = _output_path(
        output_directory,
        f"{pair.stem}_distributions.png",
    )
    figure.savefig(output_filename, dpi=180)
    plt.close(figure)
    print(f"  wrote plot: {output_filename.name}")
    return output_filename


def plot_marginals(
    pair,
    deterministic,
    monte_carlo,
    result,
    output_directory,
    plot_config,
):
    """Overlay deterministic and Monte Carlo energy and pitch marginals."""
    maximum_energy = _maximum_energy(
        deterministic.energy,
        plot_config["emax"],
    )
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(11.5, 4.8),
        constrained_layout=True,
    )

    axes[0].plot(
        deterministic.energy,
        deterministic.energy_marginal,
        linewidth=2.0,
        label="Deterministic",
    )
    axes[0].plot(
        monte_carlo.energy,
        monte_carlo.energy_marginal,
        "--",
        linewidth=1.7,
        label="Monte Carlo",
    )
    axes[0].set_xlim(float(deterministic.energy[0]), maximum_energy)
    axes[0].set_xlabel("Energy [keV]")
    axes[0].set_ylabel(
        r"$S_E(E)$ [ions cm$^{-3}$ s$^{-1}$ keV$^{-1}$]"
    )
    axes[0].set_title(
        "Pitch-integrated sink\n"
        f"normalized L1 difference = "
        f"{result.energy_marginal_l1_difference:.2e}"
    )

    axes[1].plot(
        deterministic.pitch,
        deterministic.pitch_marginal,
        linewidth=2.0,
        label="Deterministic",
    )
    axes[1].plot(
        monte_carlo.pitch,
        monte_carlo.pitch_marginal,
        "--",
        linewidth=1.7,
        label="Monte Carlo",
    )
    axes[1].set_xlim(
        float(deterministic.pitch[0]),
        float(deterministic.pitch[-1]),
    )
    axes[1].set_xlabel("Pitch")
    axes[1].set_ylabel(
        r"$S_P(p)$ [ions cm$^{-3}$ s$^{-1}$]"
    )
    axes[1].set_title(
        "Energy-integrated sink\n"
        f"normalized L1 difference = "
        f"{result.pitch_marginal_l1_difference:.2e}"
    )

    for axis in axes:
        axis.grid(True, alpha=0.3)
        if plot_config["scale"] == "log":
            axis.set_yscale("log")
    axes[0].legend()
    axes[1].legend()

    figure.suptitle(
        f"{pair.basename}: deterministic and Monte Carlo marginals"
    )
    output_filename = _output_path(
        output_directory,
        f"{pair.stem}_marginals.png",
    )
    figure.savefig(output_filename, dpi=180)
    plt.close(figure)
    print(f"  wrote plot: {output_filename.name}")
    return output_filename


def plot_rate_summary(results, output_directory, relative_tolerance):
    """Plot rates and relative differences across the complete case set."""
    case_indices = np.array([result.case_index for result in results])
    deterministic_rates = np.array(
        [result.deterministic_rate for result in results]
    )
    monte_carlo_rates = np.array(
        [result.monte_carlo_rate for result in results]
    )
    standard_errors = np.array(
        [result.monte_carlo_standard_error for result in results]
    )
    relative_differences = 100.0 * np.array(
        [result.signed_relative_rate_difference for result in results]
    )

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(8.5, 7.0),
        sharex=True,
        constrained_layout=True,
    )
    axes[0].plot(
        case_indices,
        deterministic_rates,
        "o-",
        linewidth=1.8,
        label="Deterministic",
    )
    axes[0].errorbar(
        case_indices,
        monte_carlo_rates,
        yerr=standard_errors,
        fmt="s--",
        linewidth=1.5,
        capsize=3,
        label="Monte Carlo",
    )
    axes[0].set_ylabel(r"Rate [ions cm$^{-3}$ s$^{-1}$]")
    axes[0].set_title("Total ion-sink reaction rates")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    relative_limit_percent = 100.0 * relative_tolerance
    axes[1].axhspan(
        -relative_limit_percent,
        relative_limit_percent,
        color="tab:green",
        alpha=0.15,
        label="Relative tolerance",
    )
    axes[1].axhline(0.0, color="black", linewidth=1.0)
    axes[1].plot(case_indices, relative_differences, "o-", color="tab:red")
    axes[1].set_xlabel("Distribution case index")
    axes[1].set_ylabel(r"$(R_{\rm MC}-R_{\rm det})/R_{\rm det}$ [%]")
    axes[1].set_xticks(case_indices)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    output_filename = _output_path(
        output_directory,
        "reaction_rate_summary.png",
    )
    figure.savefig(output_filename, dpi=180)
    plt.close(figure)
    print(f"Wrote plot: {output_filename.name}")
    return output_filename
