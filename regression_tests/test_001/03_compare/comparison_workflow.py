"""Load, compare, report, and plot the Test 001 mean-free-path results."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib
import numpy as np
from numpy.typing import NDArray

matplotlib.use("Agg")
import matplotlib.pyplot as plt


FloatArray = NDArray[np.float64]
TEMPERATURE_GRID_RELATIVE_TOLERANCE = 1.0e-7


def read_dataset(source: h5py.File, name: str) -> FloatArray:
    """Read one required HDF5 dataset as a double-precision array.

    Args:
        source: Open HDF5 file.
        name: Absolute dataset path.

    Returns:
        Dataset values represented as ``numpy.float64``.

    Raises:
        ValueError: If the named object is absent or is not a dataset.
    """

    if name not in source:
        raise ValueError(f"HDF5 file is missing required dataset {name}")
    dataset = source[name]
    if not isinstance(dataset, h5py.Dataset):
        raise ValueError(f"HDF5 object {name} is not a dataset")
    return np.asarray(dataset[...], dtype=np.float64)


def read_reference(path: Path) -> dict[str, FloatArray]:
    """Read the Python ADF21 reference in semantic axis order.

    Args:
        path: Reference HDF5 file produced by Stage 1.

    Returns:
        Density, energy, temperature, beam mass, and mean free path. Mean free
        path is ordered ``(density, energy, temperature)``.
    """

    with h5py.File(path, mode="r") as source:
        density = read_dataset(source, "/ne").reshape(-1)
        energy = read_dataset(source, "/eb").reshape(-1)
        temperature = read_dataset(source, "/T").reshape(-1)
        beam_mass = read_dataset(source, "/ab").reshape(-1)
        raw_mfp = read_dataset(source, "/mfp")

    # The established file layout stores the semantic (density, energy,
    # temperature) array in raw (temperature, energy, density) order.
    expected_raw_shape = (temperature.size, energy.size, density.size)
    if raw_mfp.shape != expected_raw_shape:
        raise ValueError(
            f"Reference /mfp has h5py shape {raw_mfp.shape}; "
            f"expected {expected_raw_shape} for [temperature, energy, density]"
        )

    return {
        "density": density,
        "energy": energy,
        "temperature": temperature,
        "beam_mass": beam_mass,
        "mfp": raw_mfp.transpose(2, 1, 0).copy(),
    }


def read_postprocessed(path: Path) -> dict[str, FloatArray]:
    """Read the Python-derived FIDASIM mean free paths in semantic order.

    Args:
        path: Postprocessed HDF5 file produced by Stage 2.

    Returns:
        Density, energy, temperature, beam mass, and mean free path. Mean free
        path is ordered ``(density, energy, temperature)``.
    """

    with h5py.File(path, mode="r") as source:
        density = read_dataset(source, "/dene").reshape(-1)
        energy = read_dataset(source, "/eb").reshape(-1)
        temperature = read_dataset(source, "/temperature").reshape(-1)
        beam_mass = read_dataset(source, "/beam_mass").reshape(-1)
        raw_mfp = read_dataset(source, "/mfp")

    # Stage 2 preserves the established HDF5 axis order. Make the conversion
    # to semantic order visible rather than relying on shape alone.
    expected_raw_shape = (temperature.size, energy.size, density.size)
    if raw_mfp.shape != expected_raw_shape:
        raise ValueError(
            f"Postprocessed /mfp has h5py shape {raw_mfp.shape}; "
            f"expected {expected_raw_shape} for [temperature, energy, density]"
        )

    return {
        "density": density,
        "energy": energy,
        "temperature": temperature,
        "beam_mass": beam_mass,
        "mfp": raw_mfp.transpose(2, 1, 0).copy(),
    }


def read_raw_profiles(path: Path) -> dict[str, FloatArray]:
    """Read the raw FIDASIM coordinates and profiles without modifying them.

    Args:
        path: Existing raw ``test_001.h5`` file.

    Returns:
        Raw grids plus spatial coordinates in logical ``(x, density)`` order
        and total neutral density in logical
        ``(x, density, energy, temperature)`` order.
    """

    with h5py.File(path, mode="r") as source:
        electron_density = read_dataset(source, "/dene").reshape(-1)
        energy = read_dataset(source, "/eb").reshape(-1)
        temperature = read_dataset(source, "/temperature").reshape(-1)
        beam_mass = read_dataset(source, "/beam_mass").reshape(-1)
        raw_x = read_dataset(source, "/x")
        raw_neutral_density = read_dataset(source, "/denn")

    if raw_x.ndim != 2:
        raise ValueError(f"Raw /x must have rank 2; found {raw_x.shape}")
    if raw_neutral_density.ndim != 4:
        raise ValueError(
            f"Raw /denn must have rank 4; found {raw_neutral_density.shape}"
        )

    # In h5py, /x is [density, x] and /denn is
    # [temperature, energy, density, x]. These are the reversals of the
    # semantic arrays [x, density] and [x, density, energy, temperature].
    expected_raw_x_shape = (electron_density.size, raw_x.shape[1])
    if raw_x.shape != expected_raw_x_shape:
        raise ValueError(
            f"Raw /x has h5py shape {raw_x.shape}; expected "
            f"{expected_raw_x_shape} for [density, x]"
        )
    expected_raw_density_shape = (
        temperature.size,
        energy.size,
        electron_density.size,
        raw_x.shape[1],
    )
    if raw_neutral_density.shape != expected_raw_density_shape:
        raise ValueError(
            f"Raw /denn has h5py shape {raw_neutral_density.shape}; expected "
            f"{expected_raw_density_shape} for "
            "[temperature, energy, density, x]"
        )

    return {
        "electron_density": electron_density,
        "energy": energy,
        "temperature": temperature,
        "beam_mass": beam_mass,
        "x": raw_x.transpose(1, 0).copy(),
        "neutral_density": raw_neutral_density.transpose(3, 2, 1, 0).copy(),
    }


def require_finite_positive(name: str, values: FloatArray) -> None:
    """Require an array to contain only finite positive values.

    Args:
        name: Human-readable array label.
        values: Numerical values to validate.

    Raises:
        ValueError: If any element is nonfinite, zero, or negative.
    """

    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} contains nonfinite values")
    if not np.all(values > 0.0):
        raise ValueError(f"{name} contains zero or negative values")


def validate_inputs(
    reference: dict[str, FloatArray],
    test: dict[str, FloatArray],
    raw_profiles: dict[str, FloatArray],
) -> None:
    """Validate grids, shapes, masses, and profile values before subtraction.

    Args:
        reference: Stage 1 values.
        test: Stage 2 values.
        raw_profiles: Raw FIDASIM spatial coordinates and total density.

    Raises:
        ValueError: If any schema or numerical invariant is invalid.
    """

    for dataset_name in ("density", "energy", "temperature", "beam_mass", "mfp"):
        require_finite_positive(f"reference {dataset_name}", reference[dataset_name])
        require_finite_positive(f"test {dataset_name}", test[dataset_name])

    for dataset_name in (
        "electron_density",
        "energy",
        "temperature",
        "beam_mass",
        "neutral_density",
    ):
        require_finite_positive(
            f"raw {dataset_name}",
            raw_profiles[dataset_name],
        )

    if not np.array_equal(reference["density"], test["density"]):
        raise ValueError("Reference and FIDASIM electron-density grids differ")
    if not np.array_equal(reference["energy"], test["energy"]):
        raise ValueError("Reference and FIDASIM beam-energy grids differ")
    if not np.array_equal(reference["beam_mass"], test["beam_mass"]):
        raise ValueError("Reference and FIDASIM beam masses differ")
    if not np.array_equal(reference["density"], raw_profiles["electron_density"]):
        raise ValueError("Reference and raw electron-density grids differ")
    if not np.array_equal(reference["energy"], raw_profiles["energy"]):
        raise ValueError("Reference and raw beam-energy grids differ")
    if not np.array_equal(reference["beam_mass"], raw_profiles["beam_mass"]):
        raise ValueError("Reference and raw beam masses differ")
    if not np.allclose(
        reference["temperature"],
        test["temperature"],
        rtol=TEMPERATURE_GRID_RELATIVE_TOLERANCE,
        atol=0.0,
    ):
        maximum_difference = np.max(
            np.abs(reference["temperature"] - test["temperature"])
        )
        raise ValueError(
            "Reference and FIDASIM temperature grids differ beyond "
            f"rtol={TEMPERATURE_GRID_RELATIVE_TOLERANCE:.1e}; maximum "
            f"difference is {maximum_difference:.6e} eV"
        )
    if not np.allclose(
        reference["temperature"],
        raw_profiles["temperature"],
        rtol=TEMPERATURE_GRID_RELATIVE_TOLERANCE,
        atol=0.0,
    ):
        maximum_difference = np.max(
            np.abs(reference["temperature"] - raw_profiles["temperature"])
        )
        raise ValueError(
            "Reference and raw temperature grids differ beyond "
            f"rtol={TEMPERATURE_GRID_RELATIVE_TOLERANCE:.1e}; maximum "
            f"difference is {maximum_difference:.6e} eV"
        )

    expected_mfp_shape = (
        reference["density"].size,
        reference["energy"].size,
        reference["temperature"].size,
    )
    if reference["mfp"].shape != expected_mfp_shape:
        raise ValueError(
            f"Reference mfp has shape {reference['mfp'].shape}; "
            f"expected {expected_mfp_shape}"
        )
    if test["mfp"].shape != expected_mfp_shape:
        raise ValueError(
            f"FIDASIM mfp has shape {test['mfp'].shape}; "
            f"expected {expected_mfp_shape}"
        )

    x = raw_profiles["x"]
    total_density = raw_profiles["neutral_density"]
    if x.shape[0] < 2:
        raise ValueError("Raw x must contain at least two spatial samples")
    if x.shape[1] != expected_mfp_shape[0]:
        raise ValueError(
            f"Raw x has shape {x.shape}; expected density extent "
            f"{expected_mfp_shape[0]}"
        )
    expected_profile_shape = (x.shape[0],) + expected_mfp_shape
    if total_density.shape != expected_profile_shape:
        raise ValueError(
            f"Raw total density has shape {total_density.shape}; "
            f"expected {expected_profile_shape}"
        )
    if not np.all(np.isfinite(x)):
        raise ValueError("Raw x contains nonfinite values")
    if not np.all(x >= 0.0):
        raise ValueError("Raw x contains negative coordinates")
    for density_index in range(x.shape[1]):
        if not np.all(np.diff(x[:, density_index]) > 0.0):
            raise ValueError(
                "Raw x must increase strictly for every electron-density case"
            )

    if reference["energy"].size > 6:
        raise ValueError("Legacy comparison plots support at most six energies")


def calculate_statistics(
    reference_mfp: FloatArray,
    test_mfp: FloatArray,
) -> dict[str, object]:
    """Calculate the legacy absolute-relative-error summary.

    Args:
        reference_mfp: ADF21 mean free paths in cm.
        test_mfp: FIDASIM mean free paths in cm.

    Returns:
        Relative-error array, maximum, mean, sample standard deviation, and
        semantic index of the maximum error.
    """

    relative_error = np.abs(test_mfp - reference_mfp) / np.abs(reference_mfp)
    if not np.all(np.isfinite(relative_error)):
        raise ValueError("Mean-free-path relative error contains nonfinite values")

    # Preserve the established column-major tie-breaking order explicitly.
    maximum_flat_index = int(np.argmax(relative_error.ravel(order="F")))
    maximum_index = np.unravel_index(
        maximum_flat_index,
        relative_error.shape,
        order="F",
    )
    return {
        "relative_error": relative_error,
        "maximum": float(np.max(relative_error)),
        "mean": float(np.mean(relative_error)),
        "sample_standard_deviation": float(np.std(relative_error, ddof=1)),
        "maximum_index": maximum_index,
    }


def format_report(statistics: dict[str, object], tolerance: float) -> str:
    """Format the persistent report using the established text contract.

    Args:
        statistics: Output from :func:`calculate_statistics`.
        tolerance: Maximum accepted absolute relative error.

    Returns:
        Formatted text report ending with a blank line.
    """

    maximum_percent = statistics["maximum"] * 100.0
    mean_percent = statistics["mean"] * 100.0
    standard_deviation_percent = statistics["sample_standard_deviation"] * 100.0
    tolerance_percent = tolerance * 100.0
    regression_passed = bool(
        np.isfinite(statistics["maximum"])
        and statistics["maximum"] <= tolerance
    )
    status = "PASS" if regression_passed else "FAIL"

    return (
        "==============================================\n"
        f" OVERALL REGRESSION STATUS: {status}\n"
        f" Maximum MFP relative error [%]: {maximum_percent:.2f}\n"
        f" Acceptance tolerance [%]       : {tolerance_percent:.2f}\n"
        f" Acceptance check: {maximum_percent:.2f}% <= {tolerance_percent:.2f}%\n"
        "==============================================\n"
        "\n"
        "----------------------------------------------\n"
        " Mean Free Path Regression Test Summary\n"
        "----------------------------------------------\n"
        f"  Max relative error [%]   : {maximum_percent:.2f}\n"
        f"  Mean relative error [%]  : {mean_percent:.2f}\n"
        f"  Std. relative error [%]  : {standard_deviation_percent:.2f}\n"
        "----------------------------------------------\n"
        "\n"
    )


def save_figure_pair(figure: plt.Figure, figure_directory: Path, stem: str) -> None:
    """Save one comparison figure as 600-DPI PNG and vector PDF.

    Args:
        figure: Completed Matplotlib figure.
        figure_directory: Existing destination directory.
        stem: Legacy basename without an extension.
    """

    figure.savefig(
        figure_directory / f"{stem}.png",
        dpi=600,
        bbox_inches="tight",
    )
    figure.savefig(
        figure_directory / f"{stem}.pdf",
        dpi=600,
        bbox_inches="tight",
    )


def plot_mfp_comparisons(
    figure_directory: Path,
    reference: dict[str, FloatArray],
    test: dict[str, FloatArray],
) -> None:
    """Create the three legacy reference-versus-FIDASIM MFP figures.

    Args:
        figure_directory: Destination for all figure files.
        reference: Stage 1 grids and mean free paths.
        test: Stage 2 grids and mean free paths.
    """

    colors = ("black", "red", "green", "blue", "magenta", "cyan")
    for density_index, electron_density in enumerate(reference["density"]):
        figure, axes = plt.subplots(figsize=(5.7, 7.5), facecolor="white")
        maximum_value = 0.0

        for energy_index, beam_energy in enumerate(reference["energy"]):
            reference_values = reference["mfp"][density_index, energy_index, :]
            test_values = test["mfp"][density_index, energy_index, :]
            maximum_value = max(
                maximum_value,
                float(np.max(reference_values)),
                float(np.max(test_values)),
            )
            label = (
                rf"(REF) $n_e$: {electron_density:.1e} [cm$^{{-3}}$], "
                rf"$E_b$: {beam_energy * 1.0e-3:.1e} [keV]"
            )
            axes.plot(
                reference["temperature"] * 1.0e-3,
                reference_values,
                color=colors[energy_index],
                linewidth=3.0,
                label=label,
            )
            axes.plot(
                test["temperature"] * 1.0e-3,
                test_values,
                color=colors[energy_index],
                linewidth=2.0,
                linestyle="--",
            )

        axes.set_title(r"Mean free path $\lambda_s$ [cm]", fontsize=16)
        axes.set_xlabel("T [keV]", fontsize=16)
        axes.set_ylabel(r"$\lambda_s$ [cm]", fontsize=16)
        axes.set_ylim(0.0, 1.2 * maximum_value)
        axes.tick_params(labelsize=14)
        axes.grid(True)
        axes.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.16),
            fontsize=14,
        )
        figure.subplots_adjust(bottom=0.36)
        try:
            save_figure_pair(
                figure,
                figure_directory,
                f"comparison__MFP_{density_index + 1}",
            )
        finally:
            plt.close(figure)


def plot_relative_errors(
    figure_directory: Path,
    reference: dict[str, FloatArray],
    relative_error: FloatArray,
    tolerance: float,
) -> None:
    """Create the three legacy relative-error figures.

    Args:
        figure_directory: Destination for all figure files.
        reference: Stage 1 density, energy, and temperature grids.
        relative_error: Absolute relative error in semantic MFP order.
        tolerance: Regression tolerance used as the fixed vertical limit.
    """

    colors = ("black", "red", "green", "blue", "magenta", "cyan")
    for density_index, electron_density in enumerate(reference["density"]):
        figure, axes = plt.subplots(figsize=(5.6, 4.2), facecolor="white")
        for energy_index, beam_energy in enumerate(reference["energy"]):
            label = (
                rf"$n_e$: {electron_density:.1e} [cm$^{{-3}}$], "
                rf"$E_b$: {beam_energy * 1.0e-3:.1e} [keV]"
            )
            axes.plot(
                reference["temperature"] * 1.0e-3,
                relative_error[density_index, energy_index, :] * 100.0,
                color=colors[energy_index],
                linewidth=3.0,
                label=label,
            )

        axes.set_title(r"Relative error [$\%$]", fontsize=16)
        axes.set_xlabel("T [keV]", fontsize=16)
        axes.set_ylabel(r"rel. err. [$\%$]", fontsize=16)
        axes.set_ylim(0.0, tolerance * 100.0)
        axes.tick_params(labelsize=14)
        axes.grid(True)
        axes.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.27),
            fontsize=14,
        )
        figure.subplots_adjust(bottom=0.15)
        try:
            save_figure_pair(
                figure,
                figure_directory,
                f"comparison__relative_error_{density_index + 1}",
            )
        finally:
            plt.close(figure)


def plot_maximum_error_profile(
    figure_directory: Path,
    reference: dict[str, FloatArray],
    raw_profiles: dict[str, FloatArray],
    maximum_index: tuple[int, ...],
) -> None:
    """Create the normalized profile at the maximum MFP error condition.

    Args:
        figure_directory: Destination for all figure files.
        reference: Stage 1 grids and MFP values.
        raw_profiles: Raw FIDASIM x coordinates and total neutral density.
        maximum_index: ``(density, energy, temperature)`` error index.
    """

    density_index, energy_index, temperature_index = maximum_index
    x = raw_profiles["x"][:, density_index]
    profile = raw_profiles["neutral_density"][
        :, density_index, energy_index, temperature_index
    ]
    normalized_test = profile / np.max(profile)
    reference_mfp = reference["mfp"][
        density_index,
        energy_index,
        temperature_index,
    ]
    normalized_reference = np.exp(-x / reference_mfp)

    electron_density = reference["density"][density_index]
    density_exponent = int(np.floor(np.log10(electron_density)))
    density_mantissa = electron_density / (10.0**density_exponent)
    temperature_kev = reference["temperature"][temperature_index] * 1.0e-3
    energy_kev = reference["energy"][energy_index] * 1.0e-3
    title = (
        rf"$n_e = {density_mantissa:.2f}\times 10^{{{density_exponent}}}\,"
        rf"\mathrm{{cm}}^{{-3}},\ T = {temperature_kev:.3g}\,\mathrm{{keV}},\ "
        rf"E_b = {energy_kev:.3g}\,\mathrm{{keV/amu}}$"
    )

    figure, axes = plt.subplots(figsize=(5.6, 4.2), facecolor="white")
    axes.plot(x, normalized_reference, color="black", linewidth=3.0, label="REF")
    axes.plot(
        x,
        normalized_test,
        color="black",
        linewidth=3.0,
        linestyle="--",
        label="TEST",
    )
    axes.set_xlim(0.0, 3.0 * reference_mfp)
    axes.set_yscale("linear")
    axes.set_xlabel("x [cm]", fontsize=16)
    axes.set_ylabel(r"$n_n/n_{n0}$", fontsize=24)
    axes.set_title(title, fontsize=16)
    axes.tick_params(labelsize=14)
    axes.legend(fontsize=14)
    figure.tight_layout()
    try:
        save_figure_pair(figure, figure_directory, "max_relative_error_profiles")
    finally:
        plt.close(figure)


def write_outputs(
    output_directory: Path,
    report_text: str,
    generate_plots: bool,
    reference: dict[str, FloatArray],
    test: dict[str, FloatArray],
    raw_profiles: dict[str, FloatArray],
    statistics: dict[str, object],
    tolerance: float,
) -> None:
    """Persist the report and optional figures before returning test status.

    Args:
        output_directory: Root directory for Stage 3 generated artifacts.
        report_text: Formatted persistent comparison report.
        generate_plots: Whether the seven figure pairs should be created.
        reference: Stage 1 data.
        test: Stage 2 data.
        raw_profiles: Raw FIDASIM coordinates and density profiles.
        statistics: Relative-error statistics.
        tolerance: Maximum accepted relative error.
    """

    output_directory.mkdir(parents=True, exist_ok=True)
    report_path = output_directory / "mfp_regression_report.txt"
    report_path.write_text(report_text, encoding="utf-8")

    if not generate_plots:
        return

    figure_directory = output_directory / "figures"
    figure_directory.mkdir(parents=True, exist_ok=True)
    plot_mfp_comparisons(figure_directory, reference, test)
    plot_relative_errors(
        figure_directory,
        reference,
        statistics["relative_error"],
        tolerance,
    )
    plot_maximum_error_profile(
        figure_directory,
        reference,
        raw_profiles,
        statistics["maximum_index"],
    )


def run_comparison(configuration: dict[str, object]) -> bool:
    """Run the complete Stage 3 comparison and write all diagnostics.

    Args:
        configuration: Validated result from ``read_comparison_config``.

    Returns:
        ``True`` when the maximum relative error is finite and no greater than
        the configured tolerance; otherwise ``False`` after outputs are saved.
    """

    reference = read_reference(configuration["reference_file"])
    test = read_postprocessed(configuration["postprocessed_file"])
    raw_profiles = read_raw_profiles(configuration["raw_test_file"])
    validate_inputs(reference, test, raw_profiles)

    statistics = calculate_statistics(reference["mfp"], test["mfp"])
    tolerance = float(configuration["mfp_relative_error_tolerance"])
    report_text = format_report(statistics, tolerance)
    print(report_text, end="")
    write_outputs(
        output_directory=configuration["output_directory"],
        report_text=report_text,
        generate_plots=bool(configuration["generate_plots"]),
        reference=reference,
        test=test,
        raw_profiles=raw_profiles,
        statistics=statistics,
        tolerance=tolerance,
    )

    regression_passed = bool(
        np.isfinite(statistics["maximum"])
        and statistics["maximum"] <= tolerance
    )
    return regression_passed
