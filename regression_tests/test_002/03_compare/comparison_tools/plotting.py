"""Plot relative moment differences across the selected locations."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_moment_values(results, output_directory):
    """Compare reference and converted moment values versus case index."""
    case_indices = [result["case_index"] for result in results]

    moment_panels = [
        ("density", "Density [ions/cm³]"),
        ("parallel_temperature", r"$T_\parallel$ [keV]"),
        ("perpendicular_temperature", r"$T_\perp$ [keV]"),
    ]

    figure, axes = plt.subplots(3, 1, figsize=(8.0, 9.0), sharex=True)
    for axis, (moment_name, y_label) in zip(axes, moment_panels):
        reference_values = [
            result["reference"]["moments"][moment_name] for result in results
        ]
        converted_values = [
            result["converted"]["moments"][moment_name] for result in results
        ]

        axis.plot(
            case_indices,
            reference_values,
            "o-",
            linewidth=2,
            label="Reference",
        )
        axis.plot(
            case_indices,
            converted_values,
            "s--",
            linewidth=2,
            label="Converted",
        )
        axis.set_ylabel(y_label)
        axis.grid(True, alpha=0.3)
        axis.legend()

    axes[-1].set_xlabel("Case index")
    axes[-1].set_xticks(case_indices)
    figure.suptitle("Physical moments before and after coordinate conversion")
    figure.tight_layout()

    output_path = output_directory / "moment_values.png"
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    print(f"Wrote plot: {output_path}")
    return output_path


def plot_relative_differences(results, output_directory):
    """Plot density and temperature relative differences versus case index."""
    case_indices = [result["case_index"] for result in results]
    density_errors = [result["errors"]["density"] for result in results]
    parallel_errors = [
        result["errors"]["parallel_temperature"] for result in results
    ]
    perpendicular_errors = [
        result["errors"]["perpendicular_temperature"] for result in results
    ]

    figure, axes = plt.subplots(figsize=(7.5, 5.0))
    axes.plot(case_indices, density_errors, "o-", label="Density")
    axes.plot(case_indices, parallel_errors, "s-", label=r"$T_\parallel$")
    axes.plot(case_indices, perpendicular_errors, "^-", label=r"$T_\perp$")
    axes.set_xlabel("Case index")
    axes.set_xticks(case_indices)
    axes.set_ylabel("Absolute relative difference")
    axes.set_title("Moment preservation after coordinate conversion")
    axes.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes.grid(True, alpha=0.3)
    axes.legend()
    figure.tight_layout()

    output_path = output_directory / "moment_relative_differences.png"
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    print(f"Wrote plot: {output_path}")
    return output_path
