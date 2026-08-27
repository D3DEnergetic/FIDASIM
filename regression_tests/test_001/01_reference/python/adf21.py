"""Read and interpolate the ADAS ADF21 table used by Test 001.

ADF21 stores an energy-density scan at a reference temperature and a separate
temperature scan at a reference energy and density.  Test 001 combines those
two scans using the established Test 001 factorization.

This module intentionally implements only the small ADF21 interface needed by
the regression test.  It has no dependency on the compiled ADAS Python
library, which keeps the test runnable in the standard ``FIDASIM_env``.
"""

from pathlib import Path
from typing import TextIO

import numpy as np
from scipy.interpolate import interpn


ADF21Data = dict[str, object]


def _next_line(file: TextIO, path: Path) -> str:
    """Read one record and normalize Fortran ``D`` exponent markers.

    Args:
        file: Open ADF21 text stream.
        path: Source path included in an unexpected-EOF error.

    Returns:
        The next record with signed ``D``/``d`` exponents converted to the
        ``E`` form accepted by Python's :class:`float`.

    Raises:
        ValueError: If the file ends before the expected record is present.
    """
    line = file.readline()
    if not line:
        raise ValueError(f"Unexpected end-of-file in '{path}'.")
    return (
        line.replace("D+", "E+")
        .replace("D-", "E-")
        .replace("d+", "E+")
        .replace("d-", "E-")
    )


def _keyword(line: str, name: str, default: str | None = None) -> str:
    """Extract the first whitespace-delimited value following ``NAME=``.

    Args:
        line: ADF21 header or dimensions record.
        name: Keyword name without the equals sign.
        default: Value returned when an optional keyword is absent.  Leaving
            this as ``None`` makes the keyword mandatory.

    Returns:
        The keyword value as text.  Conversion to a numeric type remains the
        caller's responsibility.

    Raises:
        ValueError: If a required keyword is absent.
    """
    marker = f"{name}="
    if marker not in line:
        if default is not None:
            return default
        raise ValueError(f"Missing {marker} in ADF21 record: {line.rstrip()}")
    return line.split(marker, 1)[1].split()[0]


def _read_values(file: TextIO, path: Path, count: int) -> np.ndarray:
    """Read a fixed number of numeric values spanning any number of records.

    ADF21 wraps long vectors according to its text format, so a logical vector
    cannot be assumed to occupy one physical line.

    Args:
        file: Open ADF21 text stream positioned at numeric data.
        path: Source path used for unexpected-EOF reporting.
        count: Number of values in the logical vector.

    Returns:
        One-dimensional ``float64`` array containing exactly ``count`` values.
    """
    values: list[float] = []
    while len(values) < count:
        values.extend(float(value) for value in _next_line(file, path).split())
    return np.asarray(values[:count], dtype=np.float64)


def _beam_species(path: Path) -> tuple[str, float]:
    """Derive the beam symbol and atomic mass from an ADAS filename.

    Supported names include forms such as ``bms93#h_h1.dat`` and
    ``bms10#he_ne10.dat``.  The portion before ``#`` identifies the dataset;
    the portion between ``#`` and ``_`` identifies the beam.

    Args:
        path: Path to the ADF21 file.

    Returns:
        Display symbol and atomic mass in atomic mass units.

    Raises:
        ValueError: If the filename identifies an unsupported beam species.
    """
    beam = path.stem.split("#")[-1].split("_", 1)[0].lower()
    try:
        return {
            "h": ("H", 1.00784),
            "d": ("D", 2.01410),
            "t": ("T", 3.01605),
            "he": ("He", 4.00260),
        }[beam]
    except KeyError as error:
        raise ValueError(f"Unknown beam species in '{path.name}'.") from error


def read_adf21(adf21_path: Path) -> ADF21Data:
    """Read one ADAS ADF21 beam-stopping table.

    The returned arrays use semantic Python ordering.  In particular,
    ``energy_density_coefficient`` has shape ``(n_energy, n_density)`` to
    match the physical description of the table.

    Args:
        adf21_path: Existing ADF21 data file.

    Returns:
        Dictionary containing provenance, beam and target information,
        reference conditions, scan grids, and stopping coefficients.

    Raises:
        FileNotFoundError: If ``adf21_path`` does not identify a file.
        ValueError: If the filename or a required ADF21 record is invalid.
    """
    path = adf21_path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"ADF21 file does not exist: '{path}'.")

    beam_symbol, beam_mass_amu = _beam_species(path)

    with path.open(encoding="utf-8") as file:
        header = _next_line(file, path)
        target_charge = int(header.split()[0])
        reference_coefficient = float(
            _keyword(header, "SCREF" if "SCREF=" in header else "SVREF")
        )
        target_species = _keyword(header, "SPEC", default="")

        # First block: k(E, n_e) tabulated at the reference temperature.
        _next_line(file, path)  # separator
        dimensions = _next_line(file, path)
        number_of_energies = int(dimensions[:5])
        number_of_densities = int(dimensions[5:10])
        reference_temperature = float(_keyword(dimensions, "TREF"))

        _next_line(file, path)  # labels
        energy = _read_values(file, path, number_of_energies)
        density = _read_values(file, path, number_of_densities)
        _next_line(file, path)  # labels

        energy_density_coefficient = np.column_stack(
            [
                _read_values(file, path, number_of_energies)
                for _ in range(number_of_densities)
            ]
        )

        # Second block: k(T) tabulated at the reference energy and density.
        _next_line(file, path)  # separator
        temperature_dimensions = _next_line(file, path)
        number_of_temperatures = int(temperature_dimensions[:5])
        temperature_reference_energy = float(
            _keyword(temperature_dimensions, "EREF")
        )
        temperature_reference_density = float(
            _keyword(temperature_dimensions, "NREF")
        )

        _next_line(file, path)  # labels
        temperature = _read_values(file, path, number_of_temperatures)
        _next_line(file, path)  # labels
        temperature_coefficient = _read_values(file, path, number_of_temperatures)

    return {
        "path": path,
        "beam_symbol": beam_symbol,
        "beam_mass_amu": beam_mass_amu,
        "target_charge": target_charge,
        "target_species": target_species,
        "reference_coefficient": reference_coefficient,
        "reference_temperature": reference_temperature,
        "energy": energy,
        "density": density,
        "energy_density_coefficient": energy_density_coefficient,
        "temperature_reference_energy": temperature_reference_energy,
        "temperature_reference_density": temperature_reference_density,
        "temperature": temperature,
        "temperature_coefficient": temperature_coefficient,
    }


def interpolate_adf21(
    data: ADF21Data,
    beam_energy: float,
    electron_density: float,
    temperature: float,
) -> float:
    """Interpolate the factorized ADF21 stopping coefficient.

    The legacy workflow defines

    ``k(E, n_e, T) = k_T(T) / k_ref * k_EN(E, n_e)``.

    Linear interpolation is used for both scans.  The caller validates that
    all requested conditions are inside the table bounds before reaching this
    function; ``interpn`` nevertheless retains its own bounds check.

    Args:
        data: Parsed table returned by :func:`read_adf21`.
        beam_energy: Beam energy in eV/amu.
        electron_density: Electron density in ``cm^-3``.
        temperature: Equal electron and ion temperature in eV.

    Returns:
        Effective stopping coefficient in ``cm^3/s``.
    """
    # The temperature correction is a one-dimensional scan through the
    # reference energy-density condition.
    coefficient_at_temperature = np.interp(
        temperature,
        np.asarray(data["temperature"]),
        np.asarray(data["temperature_coefficient"]),
    )
    # interpn expects the values to follow (density, energy), hence the
    # transpose from the reader's semantic (energy, density) layout.
    coefficient_at_energy_density = interpn(
        (np.asarray(data["density"]), np.asarray(data["energy"])),
        np.asarray(data["energy_density_coefficient"]).T,
        np.array([electron_density, beam_energy]),
        method="linear",
        bounds_error=True,
    ).item()

    # Scale the energy-density surface by the normalized temperature scan.
    return float(
        coefficient_at_temperature
        / float(data["reference_coefficient"])
        * coefficient_at_energy_density
    )
