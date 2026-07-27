"""Read and compare the physical moments stored in one file pair."""

import h5py

from regression_test_tools import ConfigError


MOMENT_NAMES = (
    "density",
    "parallel_temperature",
    "perpendicular_temperature",
)


def _read_species(h5file, filename):
    """Read and normalize the scalar species dataset."""
    if "species" not in h5file:
        raise ConfigError(f"{filename} is missing the species dataset.")
    species = h5file["species"][()]
    if isinstance(species, bytes):
        species = species.decode("utf-8")
    return str(species).strip().lower()


def _read_moments(filename, file_kind):
    """Read location, species, and moments from one Stage 1 or Stage 2 file."""
    try:
        with h5py.File(filename, mode="r") as h5file:
            if file_kind == "reference":
                selected_r = float(h5file["selected_r"][()])
                selected_z = float(h5file["selected_z"][()])
            else:
                selected_r = float(h5file["r"][0])
                selected_z = float(h5file["z"][0])

            moments = {}
            for moment_name in MOMENT_NAMES:
                moments[moment_name] = float(
                    h5file[f"moments/{moment_name}"][()]
                )

            return {
                "path": filename,
                "species": _read_species(
                    h5file=h5file,
                    filename=filename,
                ),
                "selected_r": selected_r,
                "selected_z": selected_z,
                "moments": moments,
            }
    except (KeyError, OSError, TypeError, ValueError) as error:
        raise ConfigError(f"Could not read moments from {filename}: {error}") from error


def _relative_difference(reference_value, converted_value, quantity):
    """Calculate an absolute relative difference using the reference value."""
    if reference_value == 0.0:
        raise ConfigError(
            f"Cannot calculate the relative difference for zero reference {quantity}."
        )
    return abs(converted_value - reference_value) / abs(reference_value)


def compare_file_pair(pair, case_index, relative_tolerance):
    """Validate one file pair and compare its three physical moments."""
    reference = _read_moments(pair["reference"], file_kind="reference")
    converted = _read_moments(pair["converted"], file_kind="converted")

    if reference["species"] != converted["species"]:
        raise ConfigError(f"Species differ for comparison case {case_index:03d}.")

    location_tolerance = 1.0e-10
    r_difference = abs(reference["selected_r"] - converted["selected_r"])
    z_difference = abs(reference["selected_z"] - converted["selected_z"])
    if r_difference > location_tolerance or z_difference > location_tolerance:
        raise ConfigError(f"Locations differ for comparison case {case_index:03d}.")

    errors = {}
    moment_passed = {}
    for moment_name in MOMENT_NAMES:
        errors[moment_name] = _relative_difference(
            reference_value=reference["moments"][moment_name],
            converted_value=converted["moments"][moment_name],
            quantity=moment_name,
        )
        moment_passed[moment_name] = (
            errors[moment_name] <= relative_tolerance
        )

    return {
        "case_index": case_index,
        "reference": reference,
        "converted": converted,
        "errors": errors,
        "moment_passed": moment_passed,
        "passed": all(moment_passed.values()),
    }
