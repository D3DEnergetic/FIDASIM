"""Create marginal and two-dimensional distribution comparison figures."""

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .analysis import calculate_marginals


def _resolve_plot_limits(reference_values, sampled_values, fmin, fmax):
    combined = np.concatenate(
        (reference_values.ravel(), sampled_values.ravel())
    )
    finite_values = combined[np.isfinite(combined)]
    if finite_values.size == 0:
        raise ValueError("No finite values are available for plotting.")

    if fmin in (None, "", "auto", "AUTO"):
        fmin = float(np.min(finite_values))
    if fmax in (None, "", "auto", "AUTO"):
        fmax = float(np.max(finite_values))

    return float(fmin), float(fmax)


def _transform_distribution_for_plot(values, scale):
    if scale == "log":
        positive_values = np.where(values > 0.0, values, np.nan)
        return np.log10(positive_values)
    return values.copy()


def plot_marginals(
    pair,
    reference,
    sampled,
    result,
    output_directory,
    plot_config,
):
    """Overlay reference and sampled f_E and f_P on a two-panel figure.

    Args:
        pair (FilePair): Paths identifying the current reference-sampled pair.
        reference (DistributionData): Validated reference distribution.
        sampled (DistributionData): Validated sampled distribution.
        result (ComparisonResult): Scalar comparison errors for the pair.
        output_directory (Path): Directory in which to write the PNG file.
        plot_config (dict): Canonical ``plot_data_block`` settings.

    Returns:
        None.
    """
    reference_energy, reference_pitch = calculate_marginals(
        distribution=reference,
    )
    sampled_energy, sampled_pitch = calculate_marginals(
        distribution=sampled,
    )

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    # First panel: integrate each distribution over pitch and compare f_E(E).
    axes[0].plot(reference.energy, reference_energy, label="Reference", linewidth=2)
    axes[0].plot(
        sampled.energy,
        sampled_energy,
        "--",
        label="Sampled",
        linewidth=2,
    )
    axes[0].set_xlabel("Energy [keV]")
    axes[0].set_xlim(0.0, plot_config["emax"])
    axes[0].set_ylabel(r"$f_E(E)$")
    axes[0].set_title("Pitch-integrated distribution")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    # Second panel: integrate each distribution over energy and compare f_P(P).
    axes[1].plot(reference.pitch, reference_pitch, label="Reference", linewidth=2)
    axes[1].plot(
        sampled.pitch,
        sampled_pitch,
        "--",
        label="Sampled",
        linewidth=2,
    )
    axes[1].set_xlim(reference.pitch[0], reference.pitch[-1])
    axes[1].set_xlabel("Pitch")
    axes[1].set_ylabel(r"$f_P(P)$")
    axes[1].set_title("Energy-integrated distribution")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    if plot_config["scale"] == "log":
        axes[0].set_yscale("log")
        axes[1].set_yscale("log")

    selected_r = reference.metadata["selected_r"]
    selected_z = reference.metadata["selected_z"]
    figure.suptitle(
        f"{pair.reference.name}: R = {selected_r:.2f} cm, "
        f"Z = {selected_z:.2f} cm"
    )

    metrics_text = (
        f"Status: {'PASS' if result.passed else 'FAIL'} | "
        f"Relative errors: density = {result.density_error:.2e}, "
        f"T_parallel = {result.parallel_temperature_error:.2e}, "
        f"T_perpendicular = {result.perpendicular_temperature_error:.2e}"
    )
    figure.text(0.5, 0.02, metrics_text, ha="center")
    figure.tight_layout(rect=(0.0, 0.08, 1.0, 0.93))

    output_filename = output_directory / (
        f"{pair.reference.stem}_marginals.png"
    )
    figure.savefig(output_filename)
    plt.close(figure)
    print(f"Wrote plot: {output_filename}")


def plot_moment_relative_errors(
    results,
    relative_tolerance,
    output_directory,
):
    """Plot all moment errors together with the acceptance tolerance."""
    case_indices = [result.case_index for result in results]
    density_errors = [result.density_error for result in results]
    parallel_errors = [
        result.parallel_temperature_error for result in results
    ]
    perpendicular_errors = [
        result.perpendicular_temperature_error for result in results
    ]

    figure, axes = plt.subplots(figsize=(7.5, 5.0))
    axes.plot(case_indices, density_errors, "o-", label="Density")
    axes.plot(
        case_indices,
        parallel_errors,
        "s-",
        label=r"$T_\parallel$",
    )
    axes.plot(
        case_indices,
        perpendicular_errors,
        "^-",
        label=r"$T_\perp$",
    )
    axes.axhline(
        relative_tolerance,
        color="black",
        linestyle="--",
        linewidth=1.5,
        label="Acceptance tolerance",
    )
    axes.set_xlabel("Case index")
    axes.set_xticks(case_indices)
    axes.set_ylabel("Absolute relative error")
    axes.set_title("Moment agreement after distribution sampling")
    axes.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes.grid(True, alpha=0.3)
    axes.legend()
    figure.tight_layout()

    output_filename = output_directory / "moment_relative_errors.png"
    figure.savefig(output_filename, dpi=150)
    plt.close(figure)
    print(f"Wrote plot: {output_filename}")


def plot_distributions(
    pair,
    reference,
    sampled,
    output_directory,
    plot_config,
):
    """Plot reference and sampled f(E, P) with one shared color scale.

    Args:
        pair (FilePair): Paths identifying the current reference-sampled pair.
        reference (DistributionData): Validated reference distribution.
        sampled (DistributionData): Validated sampled distribution.
        output_directory (Path): Directory in which to write the PNG file.
        plot_config (dict): Canonical ``plot_data_block`` settings.

    Returns:
        None.
    """
    reference_values = _transform_distribution_for_plot(
        values=reference.values,
        scale=plot_config["scale"],
    )
    sampled_values = _transform_distribution_for_plot(
        values=sampled.values,
        scale=plot_config["scale"],
    )
    vmin, vmax = _resolve_plot_limits(
        reference_values=reference_values,
        sampled_values=sampled_values,
        fmin=plot_config["fmin"],
        fmax=plot_config["fmax"],
    )

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(11, 4.5),
        sharex=True,
        sharey=True,
        layout="constrained",
    )
    extent = [
        reference.pitch[0],
        reference.pitch[-1],
        reference.energy[0],
        reference.energy[-1],
    ]

    # Both panels use the same vmin, vmax, and colormap for direct comparison.
    reference_image = axes[0].imshow(
        reference_values,
        origin="lower",
        aspect="auto",
        extent=extent,
        vmin=vmin,
        vmax=vmax,
        cmap=plot_config["colormap"],
    )
    axes[1].imshow(
        sampled_values,
        origin="lower",
        aspect="auto",
        extent=extent,
        vmin=vmin,
        vmax=vmax,
        cmap=plot_config["colormap"],
    )

    axes[0].set_title("Reference f(E, P)")
    axes[1].set_title("Sampled f(E, P)")
    for axis in axes:
        axis.set_xlabel("Pitch")
        axis.set_ylim(0.0, plot_config["emax"])
    axes[0].set_ylabel("Energy [keV]")

    selected_r = reference.metadata["selected_r"]
    selected_z = reference.metadata["selected_z"]
    figure.suptitle(
        f"{pair.reference.name}: R = {selected_r:.2f} cm, "
        f"Z = {selected_z:.2f} cm"
    )

    if plot_config["enable_colorbar"]:
        colorbar_label = reference.units
        if plot_config["scale"] == "log":
            colorbar_label = f"log10({colorbar_label})"
        figure.colorbar(
            reference_image,
            ax=axes,
            label=colorbar_label,
            shrink=0.9,
            pad=0.03,
        )

    output_filename = output_directory / (
        f"{pair.reference.stem}_distributions.png"
    )
    figure.savefig(output_filename)
    plt.close(figure)
    print(f"Wrote plot: {output_filename}")
