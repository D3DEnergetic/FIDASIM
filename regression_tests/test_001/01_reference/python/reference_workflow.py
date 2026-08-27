"""Build, write, and plot the Test 001 ADF21 mean-free-path reference."""

from datetime import datetime
from pathlib import Path

import f90nml
import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from adf21 import ADF21Data, interpolate_adf21, read_adf21


PROTON_MASS_KG = 1.6726e-27
ELEMENTARY_CHARGE_C = 1.6020e-19
CENTIMETRES_PER_METRE = 1.0e2


ReferenceConfig = dict[str, Path]
ReferenceGrid = dict[str, object]


def read_reference_config(config_path: Path) -> ReferenceConfig:
    """Read and resolve the Python reference-stage path configuration.

    Args:
        config_path: Fortran namelist containing one ``&reference`` block.

    Returns:
        Dictionary of absolute paths resolved relative to the configuration
        file rather than the process working directory.

    Raises:
        FileNotFoundError: If the configuration or an input file is absent.
        ValueError: If the required block or a required field is absent.
    """
    resolved_config_path = config_path.expanduser().resolve()
    if not resolved_config_path.is_file():
        raise FileNotFoundError(
            f"Reference configuration does not exist: '{resolved_config_path}'."
        )

    namelist = f90nml.read(resolved_config_path)
    if "reference" not in namelist:
        raise ValueError("Reference configuration requires an &reference block.")

    reference_block = namelist["reference"]
    required_fields = (
        "grid_input_file",
        "adf21_file",
        "output_file",
        "figure_directory",
    )
    for field_name in required_fields:
        if field_name not in reference_block:
            raise ValueError(
                f"Reference configuration is missing '{field_name}'."
            )

    config_directory = resolved_config_path.parent
    resolved_paths: ReferenceConfig = {"config_file": resolved_config_path}
    for field_name in required_fields:
        configured_path = Path(str(reference_block[field_name])).expanduser()
        if not configured_path.is_absolute():
            configured_path = config_directory / configured_path
        resolved_paths[field_name] = configured_path.resolve()

    for input_field in ("grid_input_file", "adf21_file"):
        if not resolved_paths[input_field].is_file():
            raise FileNotFoundError(
                f"Configured {input_field} does not exist: "
                f"'{resolved_paths[input_field]}'."
            )

    return resolved_paths


def read_reference_grid(grid_input_file: Path) -> ReferenceGrid:
    """Read the legacy MATLAB grid definition from its shared namelist.

    Args:
        grid_input_file: Existing Test 001 ``input.nml`` path.

    Returns:
        Dictionary containing density, energy, and temperature grids together
        with the output beam-mass metadata value.

    Raises:
        ValueError: If required values are absent or grid sizes are invalid.
    """
    namelist = f90nml.read(grid_input_file)
    if "run_input" not in namelist:
        raise ValueError("Grid input requires a &run_input block.")

    run_input = namelist["run_input"]
    required_fields = (
        "nden",
        "dene",
        "tmin",
        "tmax",
        "nt",
        "emin",
        "emax",
        "neb",
        "beam_mass",
    )
    for field_name in required_fields:
        if field_name not in run_input:
            raise ValueError(f"Grid input is missing '{field_name}'.")

    number_of_densities = int(run_input["nden"])
    number_of_temperatures = int(run_input["nt"])
    number_of_energies = int(run_input["neb"])
    if min(number_of_densities, number_of_temperatures, number_of_energies) < 1:
        raise ValueError("Grid dimensions nden, nt, and neb must be positive.")

    density = np.atleast_1d(
        np.asarray(run_input["dene"], dtype=np.float64)
    )
    if density.size != number_of_densities:
        raise ValueError(
            f"nden is {number_of_densities}, but dene contains "
            f"{density.size} values."
        )

    temperature = np.linspace(
        float(run_input["tmin"]),
        float(run_input["tmax"]),
        number_of_temperatures,
        dtype=np.float64,
    )
    energy = np.linspace(
        float(run_input["emin"]),
        float(run_input["emax"]),
        number_of_energies,
        dtype=np.float64,
    )

    return {
        "density": density,
        "temperature": temperature,
        "energy": energy,
        "beam_mass": float(run_input["beam_mass"]),
    }


def validate_reference_grid(grid: ReferenceGrid, data: ADF21Data) -> None:
    """Check that every requested condition lies within the ADF21 table.

    Args:
        grid: Grid returned by :func:`read_reference_grid`.
        data: Parsed ADF21 table returned by :func:`read_adf21`.

    Returns:
        None.

    Raises:
        ValueError: If a requested density, energy, or temperature is outside
        the corresponding ADF21 interval.
    """
    density = np.asarray(grid["density"], dtype=np.float64)
    energy = np.asarray(grid["energy"], dtype=np.float64)
    temperature = np.asarray(grid["temperature"], dtype=np.float64)
    table_density = np.asarray(data["density"], dtype=np.float64)
    table_energy = np.asarray(data["energy"], dtype=np.float64)
    table_temperature = np.asarray(data["temperature"], dtype=np.float64)

    checks = (
        (density.min(), table_density.min(), "minimum density", ">="),
        (density.max(), table_density.max(), "maximum density", "<="),
        (temperature.min(), table_temperature.min(), "minimum temperature", ">="),
        (temperature.max(), table_temperature.max(), "maximum temperature", "<="),
        (energy.min(), table_energy.min(), "minimum beam energy", ">="),
        (energy.max(), table_energy.max(), "maximum beam energy", "<="),
    )

    for requested, table_limit, label, relation in checks:
        if relation == ">=" and requested < table_limit:
            raise ValueError(
                f"Requested {label} must be {relation} {table_limit:.2e}."
            )
        if relation == "<=" and requested > table_limit:
            raise ValueError(
                f"Requested {label} must be {relation} {table_limit:.2e}."
            )


def calculate_mean_free_path(
    grid: ReferenceGrid,
    data: ADF21Data,
) -> np.ndarray:
    """Calculate the MATLAB-equivalent logical mean-free-path array.

    Args:
        grid: Density, energy, and temperature values for the reference table.
        data: Parsed ADF21 stopping-coefficient data.

    Returns:
        ``float64`` array in logical ``(nden, neb, nt)`` order, measured in cm.
    """
    density = np.asarray(grid["density"], dtype=np.float64)
    energy = np.asarray(grid["energy"], dtype=np.float64)
    temperature = np.asarray(grid["temperature"], dtype=np.float64)
    mean_free_path = np.empty(
        (density.size, energy.size, temperature.size), dtype=np.float64
    )

    beam_mass_kg = PROTON_MASS_KG * float(data["beam_mass_amu"])
    for density_index, electron_density in enumerate(density):
        for energy_index, beam_energy in enumerate(energy):
            for temperature_index, plasma_temperature in enumerate(temperature):
                stopping_coefficient = interpolate_adf21(
                    data,
                    float(beam_energy),
                    float(electron_density),
                    float(plasma_temperature),
                )
                beam_velocity = CENTIMETRES_PER_METRE * np.sqrt(
                    2.0
                    * ELEMENTARY_CHARGE_C
                    * beam_energy
                    / beam_mass_kg
                )
                mean_free_path[
                    density_index,
                    energy_index,
                    temperature_index,
                ] = beam_velocity / (electron_density * stopping_coefficient)

    return mean_free_path


def _write_matlab_text_attribute(
    target: h5py.AttributeManager,
    name: str,
    value: str,
) -> None:
    """Write a one-element string attribute like MATLAB string metadata.

    Args:
        target: HDF5 attribute manager receiving the value.
        name: Attribute name.
        value: Text stored as a one-element UTF-8 vector.

    Returns:
        None.
    """
    string_type = h5py.string_dtype(encoding="utf-8")
    target.create(name, np.asarray([value], dtype=string_type))


def _create_numeric_dataset(
    output_file: h5py.File,
    name: str,
    values: np.ndarray,
    description: str,
    units: str,
) -> h5py.Dataset:
    """Create one numeric dataset with legacy description and units metadata.

    Args:
        output_file: Open output HDF5 file.
        name: Root dataset name.
        values: Numeric data in h5py-visible file order.
        description: Legacy MATLAB dataset description.
        units: Legacy MATLAB units string.

    Returns:
        The created HDF5 dataset.
    """
    dataset = output_file.create_dataset(name, data=values, dtype=np.float64)
    _write_matlab_text_attribute(dataset.attrs, "description", description)
    _write_matlab_text_attribute(dataset.attrs, "units", units)
    return dataset


def write_reference_h5(
    output_path: Path,
    adf21_filename: str,
    grid: ReferenceGrid,
    mean_free_path: np.ndarray,
) -> None:
    """Write the MATLAB-compatible reference HDF5 schema with h5py.

    Args:
        output_path: Destination for the Python-owned reference artifact.
        adf21_filename: Basename recorded in root provenance metadata.
        grid: Reference grid and output beam-mass metadata.
        mean_free_path: Logical ``(nden, neb, nt)`` array in centimetres.

    Returns:
        None.
    """
    density = np.asarray(grid["density"], dtype=np.float64)
    energy = np.asarray(grid["energy"], dtype=np.float64)
    temperature = np.asarray(grid["temperature"], dtype=np.float64)

    # MATLAB's HDF5 interface reverses array dimensions at the file boundary.
    # h5py must therefore receive (T, eb, ne) to reproduce the existing raw
    # file shape while the documented logical order remains (ne, eb, T).
    raw_mean_free_path = np.transpose(mean_free_path, axes=(2, 1, 0))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, mode="w") as output_file:
        _create_numeric_dataset(
            output_file,
            "ne",
            density.reshape(1, density.size),
            "electron density vector",
            "[cm^-3]",
        )
        _create_numeric_dataset(
            output_file,
            "eb",
            energy.reshape(1, energy.size),
            "beam energy vector in eV/amu",
            "[eV/amu]",
        )
        _create_numeric_dataset(
            output_file,
            "ab",
            np.asarray([grid["beam_mass"]], dtype=np.float64),
            "beam mass in amu",
            "amu",
        )
        _create_numeric_dataset(
            output_file,
            "T",
            temperature.reshape(1, temperature.size),
            "plasma temperature in eV (Te = Ti)",
            "[eV]",
        )
        mean_free_path_dataset = _create_numeric_dataset(
            output_file,
            "mfp",
            raw_mean_free_path,
            "mean free path",
            "[cm]",
        )
        _write_matlab_text_attribute(
            mean_free_path_dataset.attrs,
            "dimension_order",
            "(ne, eb, T) = (nden, neb, nt)",
        )

        _write_matlab_text_attribute(
            output_file.attrs,
            "description",
            "Contains mean free paths derived from an ADF21 dataset.",
        )
        output_file.attrs["creation_date"] = datetime.now().strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        output_file.attrs["adf21_file"] = adf21_filename


def save_reference_figures(
    figure_directory: Path,
    figure_id: str,
    grid: ReferenceGrid,
    mean_free_path: np.ndarray,
) -> list[Path]:
    """Save MATLAB-equivalent temperature scans as PNG and PDF figures.

    Args:
        figure_directory: Directory receiving all generated figures.
        figure_id: Filename prefix derived from the ADF21 filename.
        grid: Density, beam-energy, and temperature vectors.
        mean_free_path: Logical ``(nden, neb, nt)`` array in centimetres.

    Returns:
        Paths of every PNG and PDF file written by the function.
    """
    density = np.asarray(grid["density"], dtype=np.float64)
    energy = np.asarray(grid["energy"], dtype=np.float64)
    temperature = np.asarray(grid["temperature"], dtype=np.float64)
    line_colors = ("k", "r", "g", "b", "m", "c")
    if energy.size > len(line_colors):
        raise ValueError(
            f"The legacy plot style supports at most {len(line_colors)} energies."
        )

    figure_directory.mkdir(parents=True, exist_ok=True)
    saved_paths = []
    for density_index, electron_density in enumerate(density):
        figure, axes = plt.subplots(figsize=(10.0, 7.0), facecolor="white")
        maximum_value = 0.0
        for energy_index, beam_energy in enumerate(energy):
            plot_values = mean_free_path[density_index, energy_index, :]
            maximum_value = max(maximum_value, float(plot_values.max()))
            legend_label = (
                rf"$n_e$: {electron_density:.1e} [cm$^{{-3}}$], "
                rf"$E_b$: {beam_energy * 1.0e-3:.1e} [keV]"
            )
            axes.plot(
                temperature * 1.0e-3,
                plot_values,
                color=line_colors[energy_index],
                linewidth=3.0,
                label=legend_label,
            )

        axes.grid(visible=True)
        axes.tick_params(axis="both", labelsize=14)
        axes.set_title(r"Mean free path $\lambda_s$ [cm]", fontsize=16)
        axes.set_ylabel(r"$\lambda_s$ [cm]", fontsize=16)
        axes.set_xlabel("T [keV]", fontsize=16)
        axes.set_ylim(0.0, 1.2 * maximum_value)
        axes.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.16),
            fontsize=14,
        )
        figure.tight_layout()

        figure_name = f"{figure_id}_MFP_{density_index + 1}"
        for extension in ("png", "pdf"):
            figure_path = figure_directory / f"{figure_name}.{extension}"
            figure.savefig(figure_path, dpi=600, bbox_inches="tight")
            saved_paths.append(figure_path)
        plt.close(figure)

    return saved_paths


def run_reference_workflow(config_path: Path) -> None:
    """Run the complete Python ADF21 reference-generation stage.

    Args:
        config_path: Path to the reference-stage Fortran namelist.

    Returns:
        None.
    """
    config = read_reference_config(config_path)
    grid = read_reference_grid(config["grid_input_file"])
    data = read_adf21(config["adf21_file"])
    validate_reference_grid(grid, data)
    mean_free_path = calculate_mean_free_path(grid, data)

    adf21_path = config["adf21_file"]
    output_path = config["output_file"]
    figure_directory = config["figure_directory"]
    write_reference_h5(
        output_path,
        adf21_path.name,
        grid,
        mean_free_path,
    )
    figure_id = adf21_path.stem.replace("#", "_")
    figure_paths = save_reference_figures(
        figure_directory,
        figure_id,
        grid,
        mean_free_path,
    )

    print(f"Wrote reference data: {output_path}")
    print(f"Wrote {len(figure_paths)} figures under: {figure_directory}")
