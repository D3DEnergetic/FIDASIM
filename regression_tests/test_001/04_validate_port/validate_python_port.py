#!/usr/bin/env python3
"""Perform the temporary MATLAB-to-Python numerical parity check.

This script exists only while Test 001 is being migrated.  It compares the
numerical arrays produced by the legacy MATLAB workflow with those produced by
Python.  It does not validate metadata, reports, figures, or provenance, and it
is not part of the normal three-stage Python workflow.

Once the Python implementation is accepted and the MATLAB implementation is
retired, this script and its configuration can be removed.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import f90nml
import h5py
import numpy as np


DATASET_NAMES = ("density", "energy", "temperature", "beam_mass", "mfp")
REQUIRED_CONFIG_FIELDS = (
    "matlab_reference_file",
    "python_reference_file",
    "matlab_test_file",
    "python_postprocessed_file",
    "relative_tolerance",
    "absolute_tolerance",
)


def parse_arguments() -> argparse.Namespace:
    """Read the parity namelist path from the command line."""
    parser = argparse.ArgumentParser(
        description="Compare MATLAB and Python Test 001 numerical arrays."
    )
    parser.add_argument("config", type=Path, help="Parity namelist path")
    return parser.parse_args()


def read_config(config_path: Path) -> dict[str, object]:
    """Read paths and tolerances from the small migration-test namelist.

    Relative paths are resolved from the directory containing the namelist so
    the check behaves the same regardless of the caller's working directory.

    Args:
        config_path: Path to a namelist containing one ``&parity`` block.

    Returns:
        Resolved MATLAB/Python artifact paths and numerical tolerances.

    Raises:
        ValueError: If the namelist is missing fields or names a missing file.
    """
    path = config_path.expanduser().resolve()
    blocks = f90nml.read(path)
    if set(blocks) != {"parity"}:
        raise ValueError("Configuration must contain exactly one &parity block.")

    block = blocks["parity"]
    missing = [field for field in REQUIRED_CONFIG_FIELDS if field not in block]
    if missing:
        raise ValueError(f"Missing &parity fields: {', '.join(missing)}")

    config: dict[str, object] = {}
    for field in REQUIRED_CONFIG_FIELDS[:4]:
        value = Path(str(block[field])).expanduser()
        if not value.is_absolute():
            value = path.parent / value
        value = value.resolve()
        if not value.is_file():
            raise ValueError(f"{field} does not exist: {value}")
        config[field] = value

    for field in REQUIRED_CONFIG_FIELDS[4:]:
        value = float(block[field])
        if not np.isfinite(value) or value < 0.0:
            raise ValueError(f"{field} must be finite and non-negative.")
        config[field] = value

    return config


def read_reference(path: Path) -> dict[str, np.ndarray]:
    """Read a MATLAB- or Python-generated ADF21 reference artifact.

    h5py exposes the three-dimensional ``mfp`` array in the reverse of the
    MATLAB logical order, so it is transposed to ``(density, energy,
    temperature)`` immediately after reading.
    """
    with h5py.File(path, "r") as source:
        return {
            "density": np.asarray(source["ne"], dtype=float).reshape(-1),
            "energy": np.asarray(source["eb"], dtype=float).reshape(-1),
            "temperature": np.asarray(source["T"], dtype=float).reshape(-1),
            "beam_mass": np.asarray(source["ab"], dtype=float).reshape(-1),
            "mfp": np.asarray(source["mfp"], dtype=float).transpose(2, 1, 0),
        }


def read_test_result(path: Path) -> dict[str, np.ndarray]:
    """Read a MATLAB- or Python-generated FIDASIM MFP artifact."""
    with h5py.File(path, "r") as source:
        if "mfp" not in source:
            raise ValueError(f"Parity artifact has no /mfp dataset: {path}")
        return {
            "density": np.asarray(source["dene"], dtype=float).reshape(-1),
            "energy": np.asarray(source["eb"], dtype=float).reshape(-1),
            "temperature": np.asarray(source["temperature"], dtype=float).reshape(-1),
            "beam_mass": np.asarray(source["beam_mass"], dtype=float).reshape(-1),
            "mfp": np.asarray(source["mfp"], dtype=float).transpose(2, 1, 0),
        }


def compare_artifacts(
    label: str,
    matlab: dict[str, np.ndarray],
    python: dict[str, np.ndarray],
    relative_tolerance: float,
    absolute_tolerance: float,
) -> bool:
    """Compare every numerical array and print a compact review summary.

    Args:
        label: Human-readable artifact name printed with each result.
        matlab: Arrays read from the legacy MATLAB artifact.
        python: Arrays read from the Python artifact.
        relative_tolerance: Relative tolerance passed to ``numpy.allclose``.
        absolute_tolerance: Absolute tolerance passed to ``numpy.allclose``.

    Returns:
        ``True`` only when all shapes and values agree within tolerance.
    """
    passed = True
    for name in DATASET_NAMES:
        expected = matlab[name]
        actual = python[name]
        shape_matches = expected.shape == actual.shape
        values_match = shape_matches and np.allclose(
            actual,
            expected,
            rtol=relative_tolerance,
            atol=absolute_tolerance,
            equal_nan=False,
        )

        if shape_matches:
            difference = np.abs(actual - expected)
            scale = np.maximum(np.abs(expected), np.finfo(float).tiny)
            maximum_absolute = float(np.max(difference))
            maximum_relative = float(np.max(difference / scale))
            detail = f"max abs={maximum_absolute:.6e}, max rel={maximum_relative:.6e}"
        else:
            detail = f"shapes differ: MATLAB {expected.shape}, Python {actual.shape}"

        status = "PASS" if values_match else "FAIL"
        print(f"{status}: {label} /{name}: {detail}")
        passed = passed and bool(values_match)

    return passed


def main() -> int:
    """Run the two numerical comparisons and return their combined status."""
    try:
        config = read_config(parse_arguments().config)
        reference_passed = compare_artifacts(
            "reference",
            read_reference(config["matlab_reference_file"]),
            read_reference(config["python_reference_file"]),
            config["relative_tolerance"],
            config["absolute_tolerance"],
        )
        postprocessing_passed = compare_artifacts(
            "postprocessing",
            read_test_result(config["matlab_test_file"]),
            read_test_result(config["python_postprocessed_file"]),
            config["relative_tolerance"],
            config["absolute_tolerance"],
        )
    except (OSError, TypeError, ValueError) as error:
        print(f"PYTHON PORT PARITY: FAIL ({error})")
        return 1

    passed = reference_passed and postprocessing_passed
    print(f"PYTHON PORT PARITY: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
