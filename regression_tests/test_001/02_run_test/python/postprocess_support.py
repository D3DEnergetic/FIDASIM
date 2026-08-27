"""Explicit numerical and plotting operations for Test 001 post-processing.

The source HDF5 file was written through the Fortran HDF5 interface.  h5py
therefore exposes its multidimensional arrays with the reverse of the logical
Fortran dimension order.  This module converts those arrays immediately after
reading them, validates the logical arrays, and applies the established
mean-free-path estimator without introducing a fitted estimator.
"""

from __future__ import annotations

import os
from datetime import datetime, timezone
from pathlib import Path

# Matplotlib requires a writable configuration directory even with the Agg
# backend.  Keep its cache beside the derived, ignored outputs when the caller
# has not selected another location.
if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = str(
        Path(__file__).resolve().parent / "output_data" / ".matplotlib_cache"
    )

import f90nml
import h5py
import matplotlib
import numpy as np
from numpy.typing import NDArray

matplotlib.use("Agg")
import matplotlib.pyplot as plt


FloatArray = NDArray[np.float64]

GRID_UNIFORMITY_RELATIVE_TOLERANCE = 1.0e-12
LEVEL_SUM_RELATIVE_TOLERANCE = 1.0e-12


def resolve_configured_path(value: object, config_directory: Path, field: str) -> Path:
    """Resolve one configured path relative to the configuration file.

    Args:
        value: Value read by ``f90nml`` for the path field.
        config_directory: Directory containing the namelist file.
        field: Configuration field name used in error messages.

    Returns:
        An absolute, normalized path.

    Raises:
        ValueError: If the configured value is not a non-empty string.
    """
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"&postprocess field '{field}' must be a path string")

    configured_path = Path(value).expanduser()
    if not configured_path.is_absolute():
        configured_path = config_directory / configured_path

    return configured_path.resolve()


def read_configuration(config_path: Path) -> dict[str, Path]:
    """Read and validate the Stage 2 namelist configuration.

    Args:
        config_path: Namelist file supplied on the command line.

    Returns:
        Absolute paths for the configuration, source HDF5 file, output HDF5
        file, and figure directory.

    Raises:
        OSError: If the namelist cannot be opened.
        ValueError: If its block or fields do not match the expected schema.
    """
    resolved_config_path = config_path.expanduser().resolve()
    if not resolved_config_path.is_file():
        raise ValueError(f"Configuration file does not exist: {resolved_config_path}")

    try:
        namelist = f90nml.read(resolved_config_path)
    except (OSError, ValueError) as error:
        raise ValueError(
            f"Unable to read configuration '{resolved_config_path}': {error}"
        ) from error

    block_names = set(namelist.keys())
    if block_names != {"postprocess"}:
        raise ValueError(
            "Configuration must contain exactly one &postprocess block; "
            f"found {sorted(block_names)}"
        )

    block = namelist["postprocess"]
    required_fields = {"source_h5", "output_h5", "figures_dir"}
    actual_fields = set(block.keys())
    missing_fields = sorted(required_fields - actual_fields)
    unknown_fields = sorted(actual_fields - required_fields)

    if missing_fields:
        raise ValueError(
            "Missing &postprocess field(s): " + ", ".join(missing_fields)
        )
    if unknown_fields:
        raise ValueError(
            "Unknown &postprocess field(s): " + ", ".join(unknown_fields)
        )

    config_directory = resolved_config_path.parent
    source_h5 = resolve_configured_path(
        block["source_h5"], config_directory, "source_h5"
    )
    output_h5 = resolve_configured_path(
        block["output_h5"], config_directory, "output_h5"
    )
    figures_dir = resolve_configured_path(
        block["figures_dir"], config_directory, "figures_dir"
    )

    if source_h5 == output_h5:
        raise ValueError("source_h5 and output_h5 must identify different files")

    return {
        "config_path": resolved_config_path,
        "source_h5": source_h5,
        "output_h5": output_h5,
        "figures_dir": figures_dir,
    }


def read_numeric_dataset(h5file: h5py.File, name: str) -> FloatArray:
    """Read one required HDF5 dataset as a double-precision NumPy array.

    Args:
        h5file: Source HDF5 file opened read-only.
        name: Absolute HDF5 dataset name.

    Returns:
        Dataset values converted to ``numpy.float64``.

    Raises:
        ValueError: If the object is missing, is not a dataset, or cannot be
        represented as double-precision numeric data.
    """
    if name not in h5file:
        raise ValueError(f"Source HDF5 file is missing required dataset {name}")

    dataset = h5file[name]
    if not isinstance(dataset, h5py.Dataset):
        raise ValueError(f"Source HDF5 object {name} is not a dataset")

    try:
        values = np.asarray(dataset[...], dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"Source HDF5 dataset {name} is not numeric: {error}"
        ) from error

    return values


def read_source_h5(source_path: Path) -> dict[str, FloatArray]:
    """Read the raw FIDASIM data without opening the source for modification.

    Args:
        source_path: Existing ``test_001.h5`` path.

    Returns:
        Raw h5py-visible arrays required for the mean-free-path calculation.

    Raises:
        OSError: If HDF5 cannot open or read the source.
        ValueError: If a required dataset is absent or nonnumeric.
    """
    if not source_path.is_file():
        raise ValueError(f"Source HDF5 file does not exist: {source_path}")

    raw_data: dict[str, FloatArray] = {}
    required_names = (
        "/temperature",
        "/eb",
        "/dene",
        "/beam_mass",
        "/x",
        "/denn",
        "/denn_n",
    )

    with h5py.File(source_path, mode="r") as h5file:
        for name in required_names:
            raw_data[name] = read_numeric_dataset(h5file, name)

    return raw_data


def transpose_source_arrays(raw_data: dict[str, FloatArray]) -> dict[str, FloatArray]:
    """Convert raw h5py axis order to the documented logical array order.

    Args:
        raw_data: Arrays returned by :func:`read_source_h5`.

    Returns:
        Logical arrays with ``x(nx, nden)``,
        ``denn(nx, nden, neb, nt)``, and
        ``denn_n(nlev, nx, nden, neb, nt)``.

    Raises:
        ValueError: If a raw multidimensional array has the wrong rank.
    """
    raw_x = raw_data["/x"]
    raw_denn = raw_data["/denn"]
    raw_denn_n = raw_data["/denn_n"]

    if raw_x.ndim != 2:
        raise ValueError(f"Raw /x must have rank 2; found shape {raw_x.shape}")
    if raw_denn.ndim != 4:
        raise ValueError(f"Raw /denn must have rank 4; found shape {raw_denn.shape}")
    if raw_denn_n.ndim != 5:
        raise ValueError(
            f"Raw /denn_n must have rank 5; found shape {raw_denn_n.shape}"
        )

    logical_data = {
        "temperature": raw_data["/temperature"].copy(),
        "eb": raw_data["/eb"].copy(),
        "dene": raw_data["/dene"].copy(),
        "beam_mass": raw_data["/beam_mass"].copy(),
        # h5py (nden, nx) -> logical (nx, nden)
        "x": raw_x.transpose(1, 0).copy(),
        # h5py (nt, neb, nden, nx) -> logical (nx, nden, neb, nt)
        "denn": raw_denn.transpose(3, 2, 1, 0).copy(),
        # h5py (nt, neb, nden, nx, nlev) -> logical
        # (nlev, nx, nden, neb, nt)
        "denn_n": raw_denn_n.transpose(4, 3, 2, 1, 0).copy(),
    }

    return logical_data


def require_finite_positive(name: str, values: FloatArray) -> None:
    """Require every array value to be finite and strictly positive.

    Args:
        name: Human-readable dataset name.
        values: Values to validate.

    Raises:
        ValueError: If an element is nonfinite or nonpositive.
    """
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} contains nonfinite values")
    if not np.all(values > 0.0):
        raise ValueError(f"{name} contains zero or negative values")


def validate_logical_data(logical_data: dict[str, FloatArray]) -> float:
    """Validate logical shapes, grids, profiles, and level-density sums.

    Args:
        logical_data: Transposed source arrays.

    Returns:
        Maximum pointwise relative mismatch between ``sum(denn_n)`` and
        ``denn``.

    Raises:
        ValueError: If any shape or numerical invariant is invalid.
    """
    temperature = logical_data["temperature"]
    eb = logical_data["eb"]
    dene = logical_data["dene"]
    beam_mass = logical_data["beam_mass"]
    x = logical_data["x"]
    denn = logical_data["denn"]
    denn_n = logical_data["denn_n"]

    if temperature.ndim != 1:
        raise ValueError(
            f"/temperature must be one-dimensional; found {temperature.shape}"
        )
    if eb.ndim != 1:
        raise ValueError(f"/eb must be one-dimensional; found {eb.shape}")
    if dene.ndim != 1:
        raise ValueError(f"/dene must be one-dimensional; found {dene.shape}")
    if beam_mass.shape != (1,):
        raise ValueError(f"/beam_mass must have shape (1,); found {beam_mass.shape}")

    nt = temperature.size
    neb = eb.size
    nden = dene.size
    if x.ndim != 2:
        raise ValueError(f"Logical x must have rank 2; found shape {x.shape}")
    nx = x.shape[0]

    expected_x_shape = (nx, nden)
    expected_denn_shape = (nx, nden, neb, nt)
    if x.shape != expected_x_shape:
        raise ValueError(f"Logical x has shape {x.shape}; expected {expected_x_shape}")
    if denn.shape != expected_denn_shape:
        raise ValueError(
            f"Logical denn has shape {denn.shape}; expected {expected_denn_shape}"
        )

    if denn_n.ndim != 5:
        raise ValueError(
            f"Logical denn_n must have rank 5; found shape {denn_n.shape}"
        )
    nlev = denn_n.shape[0]
    expected_denn_n_shape = (nlev, nx, nden, neb, nt)
    if denn_n.shape != expected_denn_n_shape:
        raise ValueError(
            "Logical denn_n has shape "
            f"{denn_n.shape}; expected {expected_denn_n_shape}"
        )

    if nx < 2:
        raise ValueError(
            f"Legacy estimator requires at least two x points; found nx={nx}"
        )

    require_finite_positive("/temperature", temperature)
    require_finite_positive("/eb", eb)
    require_finite_positive("/dene", dene)
    require_finite_positive("/beam_mass", beam_mass)
    require_finite_positive("logical denn", denn)
    require_finite_positive("logical denn_n", denn_n)

    if not np.all(np.isfinite(x)):
        raise ValueError("logical x contains nonfinite values")

    for density_index in range(nden):
        spacings = np.diff(x[:, density_index])
        if not np.all(spacings > 0.0):
            raise ValueError(f"x[:, {density_index}] is not strictly increasing")
        if not np.allclose(
            spacings,
            spacings[0],
            rtol=GRID_UNIFORMITY_RELATIVE_TOLERANCE,
            atol=0.0,
        ):
            maximum_deviation = float(np.max(np.abs(spacings - spacings[0])))
            raise ValueError(
                f"x[:, {density_index}] is not uniform; maximum spacing "
                f"deviation is {maximum_deviation:.6e} cm"
            )

    summed_levels = np.sum(denn_n, axis=0)
    relative_level_error = np.abs(summed_levels - denn) / denn
    maximum_relative_level_error = float(np.max(relative_level_error))
    if not np.allclose(
        summed_levels,
        denn,
        rtol=LEVEL_SUM_RELATIVE_TOLERANCE,
        atol=0.0,
    ):
        raise ValueError(
            "sum(denn_n, axis=level) does not reproduce denn within "
            f"rtol={LEVEL_SUM_RELATIVE_TOLERANCE:.1e}; maximum relative "
            f"error is {maximum_relative_level_error:.6e}"
        )

    return maximum_relative_level_error


def compute_legacy_mfp(x: FloatArray, denn: FloatArray) -> FloatArray:
    """Apply the established Test 001 mean-free-path estimator exactly.

    For each density, energy, and temperature profile, the procedure divides
    by that profile's maximum, takes the natural logarithm, forms ``nx - 1``
    local forward-difference slopes using the first uniform grid spacing,
    negates their reciprocals, and takes the ordinary arithmetic mean.  The
    current Test 001 source has ``nx = 100`` and therefore uses 99 slopes.  No
    regression fit or downstream asymptotic estimator is applied.

    Args:
        x: Spatial coordinates in logical ``(nx, nden)`` order.
        denn: Total neutral density in logical ``(nx, nden, neb, nt)`` order.

    Returns:
        Mean free path in logical ``(nden, neb, nt)`` order, in centimetres.

    Raises:
        ValueError: If a local reciprocal slope is nonfinite or nonpositive.
    """
    nden = denn.shape[1]
    neb = denn.shape[2]
    nt = denn.shape[3]
    expected_local_slope_count = x.shape[0] - 1
    mfp = np.empty((nden, neb, nt), dtype=np.float64)

    for density_index in range(nden):
        dx = x[1, density_index] - x[0, density_index]

        for energy_index in range(neb):
            for temperature_index in range(nt):
                profile = denn[:, density_index, energy_index, temperature_index]
                normalized_profile = profile / np.max(profile)
                natural_log_profile = np.log(normalized_profile)
                local_log_slopes = np.diff(natural_log_profile) / dx
                local_reciprocal_slopes = -1.0 / local_log_slopes

                if local_reciprocal_slopes.size != expected_local_slope_count:
                    raise ValueError(
                        "Legacy estimator did not produce exactly "
                        f"{expected_local_slope_count} local slopes"
                    )
                if not np.all(np.isfinite(local_reciprocal_slopes)):
                    raise ValueError(
                        "Legacy estimator produced a nonfinite local reciprocal "
                        f"slope at indices ({density_index}, {energy_index}, "
                        f"{temperature_index})"
                    )
                if not np.all(local_reciprocal_slopes > 0.0):
                    raise ValueError(
                        "Legacy estimator produced a nonpositive local reciprocal "
                        f"slope at indices ({density_index}, {energy_index}, "
                        f"{temperature_index})"
                    )

                mfp[density_index, energy_index, temperature_index] = np.mean(
                    local_reciprocal_slopes
                )

    require_finite_positive("computed mfp", mfp)
    return mfp


def set_vector_attributes(
    dataset: h5py.Dataset, description: str, units: str, dimension: str
) -> None:
    """Attach consistent metadata to a one-dimensional output dataset.

    Args:
        dataset: Newly created HDF5 dataset.
        description: Human-readable physical meaning.
        units: Physical units.
        dimension: Logical dimension label.
    """
    dataset.attrs["description"] = description
    dataset.attrs["units"] = units
    dataset.attrs["dimensions"] = dimension
    dataset.attrs["h5py_axis_order"] = dimension


def write_postprocessed_h5(
    output_path: Path,
    source_path: Path,
    config_path: Path,
    logical_data: dict[str, FloatArray],
    computed_mfp: FloatArray,
) -> None:
    """Write a separate, self-contained derived HDF5 artifact.

    Args:
        output_path: Destination that is distinct from the raw source file.
        source_path: Read-only raw FIDASIM input used for provenance.
        config_path: Namelist used for this run.
        logical_data: Validated grid arrays and beam mass.
        computed_mfp: Logical ``(nden, neb, nt)`` mean free path.
    Raises:
        OSError: If the output directory or HDF5 file cannot be written.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    relative_source = os.path.relpath(source_path, output_path.parent)
    relative_config = os.path.relpath(config_path, output_path.parent)
    created_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    local_slope_count = logical_data["x"].shape[0] - 1

    # Store the array in (nt, neb, nden) order while numerical work uses the
    # documented semantic (nden, neb, nt) order.
    mfp_h5py_order = computed_mfp.transpose(2, 1, 0).copy()

    with h5py.File(output_path, mode="w") as h5file:
        h5file.attrs["artifact_kind"] = "FIDASIM test_001 derived mean free path"
        h5file.attrs["created_by"] = (
            "regression_tests/test_001/02_run_test/python/postprocess.py"
        )
        h5file.attrs["created_utc"] = created_utc
        h5file.attrs["configuration_file"] = relative_config
        h5file.attrs["source_h5"] = relative_source
        h5file.attrs["source_h5_open_mode"] = "read-only"
        h5file.attrs["estimator"] = (
            "normalize each profile by its maximum; natural logarithm; "
            f"{local_slope_count} forward-difference slopes; negative "
            "reciprocal; ordinary mean"
        )
        h5file.attrs["estimator_local_slope_count"] = local_slope_count
        h5file.attrs["downstream_fit_applied"] = False

        temperature_dataset = h5file.create_dataset(
            "/temperature", data=logical_data["temperature"]
        )
        set_vector_attributes(
            temperature_dataset,
            "plasma temperature grid",
            "eV",
            "[temperature]",
        )

        energy_dataset = h5file.create_dataset("/eb", data=logical_data["eb"])
        set_vector_attributes(
            energy_dataset,
            "beam energy grid",
            "eV/amu",
            "[eb]",
        )

        density_dataset = h5file.create_dataset("/dene", data=logical_data["dene"])
        set_vector_attributes(
            density_dataset,
            "electron density grid",
            "cm^-3",
            "[dene]",
        )

        beam_mass_dataset = h5file.create_dataset(
            "/beam_mass", data=logical_data["beam_mass"]
        )
        set_vector_attributes(
            beam_mass_dataset,
            "beam atomic mass",
            "amu",
            "[1]",
        )

        mfp_dataset = h5file.create_dataset("/mfp", data=mfp_h5py_order)
        mfp_dataset.attrs["description"] = (
            "Mean free path from the legacy ordinary mean of "
            f"{local_slope_count} local negative reciprocal logarithmic slopes"
        )
        mfp_dataset.attrs["units"] = "cm"
        mfp_dataset.attrs["dimensions"] = "[dene, eb, temperature]"
        mfp_dataset.attrs["h5py_axis_order"] = "[temperature, eb, dene]"


def save_figure_pair(figure: plt.Figure, figures_dir: Path, stem: str) -> None:
    """Save one Matplotlib figure in both PNG and PDF formats.

    Args:
        figure: Completed figure.
        figures_dir: Output directory.
        stem: Filename without an extension.
    """
    png_path = figures_dir / f"{stem}.png"
    pdf_path = figures_dir / f"{stem}.pdf"
    figure.savefig(png_path, dpi=600, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")


def plot_mfp_figures(
    figures_dir: Path,
    temperature: FloatArray,
    eb: FloatArray,
    dene: FloatArray,
    computed_mfp: FloatArray,
) -> None:
    """Reproduce the three legacy FIDASIM mean-free-path figures.

    Args:
        figures_dir: Destination for PNG and PDF files.
        temperature: Temperature grid in eV.
        eb: Beam-energy grid in eV/amu.
        dene: Electron-density grid in cm^-3.
        computed_mfp: Logical ``(nden, neb, nt)`` mean free path in cm.
    """
    colors = ("black", "red", "green", "blue", "magenta", "cyan")

    for density_index in range(dene.size):
        figure, axes = plt.subplots(figsize=(10.0, 7.0), facecolor="white")

        for energy_index in range(eb.size):
            label = (
                rf"$n_e$: {dene[density_index]:.1e} cm$^{{-3}}$, "
                rf"$E_b$: {eb[energy_index] * 1.0e-3:.1e} keV/amu"
            )
            axes.plot(
                temperature * 1.0e-3,
                computed_mfp[density_index, energy_index, :],
                color=colors[energy_index],
                linewidth=3.0,
                label=label,
            )

        maximum_mfp = float(np.max(computed_mfp[density_index, :, :]))
        axes.set_ylim(0.0, 1.2 * maximum_mfp)
        axes.set_xlabel("T [keV]", fontsize=16)
        axes.set_ylabel(r"$\lambda_s$ [cm]", fontsize=16)
        axes.set_title(r"Mean free path $\lambda_s$ [cm]", fontsize=16)
        axes.tick_params(labelsize=14)
        axes.grid(True)
        axes.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.16),
            fontsize=11,
            frameon=True,
            ncol=1,
        )
        figure.subplots_adjust(bottom=0.34)

        stem = f"FIDASIM_MFP_{density_index + 1}"
        save_figure_pair(figure, figures_dir, stem)
        plt.close(figure)


def plot_example_profile(
    figures_dir: Path,
    x: FloatArray,
    denn: FloatArray,
    temperature: FloatArray,
    eb: FloatArray,
    dene: FloatArray,
) -> None:
    """Reproduce the first-profile legacy log10 neutral-density figure.

    Args:
        figures_dir: Destination for PNG and PDF files.
        x: Logical spatial grid in cm.
        denn: Logical total neutral-density profiles.
        temperature: Temperature grid in eV.
        eb: Beam-energy grid in eV/amu.
        dene: Electron-density grid in cm^-3.
    """
    density_index = 0
    energy_index = 0
    temperature_index = 0

    profile = denn[:, density_index, energy_index, temperature_index]
    normalized_profile = profile / np.max(profile)
    log10_profile = np.log10(normalized_profile)

    density_value = float(dene[density_index])
    density_exponent = int(np.floor(np.log10(density_value)))
    density_mantissa = density_value / (10.0**density_exponent)
    temperature_kev = float(temperature[temperature_index] * 1.0e-3)
    energy_kev_per_amu = float(eb[energy_index] * 1.0e-3)

    title = (
        rf"$n_e = {density_mantissa:.2f}\times 10^{{{density_exponent}}}\,"
        rf"\mathrm{{cm}}^{{-3}},\ T = {temperature_kev:.3g}\,"
        rf"\mathrm{{keV}},\ E_b = {energy_kev_per_amu:.3g}\,"
        rf"\mathrm{{keV/amu}}$"
    )

    figure, axes = plt.subplots(figsize=(10.0, 7.0), facecolor="white")
    axes.plot(
        x[:, density_index],
        log10_profile,
        color="tab:blue",
        linewidth=3.0,
    )
    axes.set_xlabel("x [cm]", fontsize=16)
    axes.set_ylabel(r"$\log_{10}(n_n/n_{n0})$", fontsize=18)
    axes.set_title(title, fontsize=16)
    axes.tick_params(labelsize=14)
    axes.grid(True)
    figure.tight_layout()

    save_figure_pair(figure, figures_dir, "Example_neutral_density_profile")
    plt.close(figure)


def run_postprocessing(config_path: Path) -> None:
    """Execute the complete read-only-source Stage 2 workflow.

    Args:
        config_path: Path to ``input_config.nml``.

    Raises:
        OSError: If an input cannot be read or an output cannot be written.
        ValueError: If configuration, dimensions, or invariants fail.
    """
    paths = read_configuration(config_path)
    source_path = paths["source_h5"]
    output_path = paths["output_h5"]
    figures_dir = paths["figures_dir"]

    raw_data = read_source_h5(source_path)
    logical_data = transpose_source_arrays(raw_data)
    maximum_level_sum_error = validate_logical_data(logical_data)
    computed_mfp = compute_legacy_mfp(logical_data["x"], logical_data["denn"])

    write_postprocessed_h5(
        output_path=output_path,
        source_path=source_path,
        config_path=paths["config_path"],
        logical_data=logical_data,
        computed_mfp=computed_mfp,
    )

    figures_dir.mkdir(parents=True, exist_ok=True)
    plot_mfp_figures(
        figures_dir=figures_dir,
        temperature=logical_data["temperature"],
        eb=logical_data["eb"],
        dene=logical_data["dene"],
        computed_mfp=computed_mfp,
    )
    plot_example_profile(
        figures_dir=figures_dir,
        x=logical_data["x"],
        denn=logical_data["denn"],
        temperature=logical_data["temperature"],
        eb=logical_data["eb"],
        dene=logical_data["dene"],
    )

    print("PASS: Test 001 Python FIDASIM post-processing completed")
    print(f"  configuration: {paths['config_path']}")
    print(f"  read-only source: {source_path}")
    print(f"  logical x shape: {logical_data['x'].shape}")
    print(f"  logical denn shape: {logical_data['denn'].shape}")
    print(f"  logical denn_n shape: {logical_data['denn_n'].shape}")
    print(
        "  maximum relative level-sum error: "
        f"{maximum_level_sum_error:.6e}"
    )
    print(f"  computed logical mfp shape: {computed_mfp.shape}")
    print(f"  local reciprocal slopes per profile: {logical_data['x'].shape[0] - 1}")
    print(f"  written /mfp h5py shape: {computed_mfp.transpose(2, 1, 0).shape}")

    print(f"  derived HDF5: {output_path}")
    print(f"  figures: {figures_dir} (4 PNG and 4 PDF)")
