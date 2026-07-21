"""Coordinate conversion of Stage 1 references into FIDASIM HDF5 files."""

from pathlib import Path

import h5py
import numpy as np

from regression_test_tools import ConfigError, normalize_path, read_namelist

from .config import read_config
from .plotting import plot_converted_distribution
from .remapping import remap_to_uniform_grid
from .reporting import write_moment_report
from .transformation import (
    ERG_PER_KEV,
    SPECIES_MASS_AMU,
    transform_to_nonrelativistic_energy_pitch,
)


def _indexed_path(base_path, case_index):
    """Insert a three-digit case number before the filename extension."""
    return base_path.with_name(
        f"{base_path.stem}_{case_index:03d}{base_path.suffix}"
    )


def _discover_reference_files(reference_config):
    """Use the Stage 1 configuration to reconstruct its indexed output paths."""
    config_path, blocks = read_namelist(config_path=reference_config)
    if "input" not in blocks or "save_data_block" not in blocks:
        raise ConfigError("The Stage 1 configuration is missing required blocks.")

    input_block = blocks["input"]
    save_block = blocks["save_data_block"]
    if "r_locations" not in input_block or "output_filename" not in save_block:
        raise ConfigError("The Stage 1 configuration lacks location or output data.")

    output_base = normalize_path(
        value=save_block["output_filename"],
        config_path=config_path,
        field_label="Stage 1 output_filename",
    )
    number_of_cases = len(np.atleast_1d(input_block["r_locations"]))

    reference_paths = []
    for case_index in range(1, number_of_cases + 1):
        reference_path = _indexed_path(output_base, case_index)
        if not reference_path.is_file():
            raise ConfigError(f"Reference file does not exist: {reference_path}")
        reference_paths.append(reference_path)
    return reference_paths


def _calculate_output_moments(energy, pitch, f_array, mass_grams):
    """Calculate density and relativistic pressure moments on the E-P grid."""
    denergy = energy[1] - energy[0]
    dpitch = pitch[1] - pitch[0]
    cell_area = denergy * dpitch
    density = float(np.sum(f_array) * cell_area)

    energy_2d = energy[:, np.newaxis]
    pitch_2d = pitch[np.newaxis, :]
    u_squared = 2.0 * ERG_PER_KEV * energy_2d / mass_grams
    speed_of_light = 2.99792458e10
    gamma = np.sqrt(1.0 + u_squared / speed_of_light**2)
    pressure_energy = 2.0 * energy_2d / gamma

    parallel_pressure = float(
        np.sum(pressure_energy * pitch_2d**2 * f_array) * cell_area
    )
    perpendicular_pressure = 0.5 * float(
        np.sum(pressure_energy * (1.0 - pitch_2d**2) * f_array) * cell_area
    )
    return density, parallel_pressure / density, perpendicular_pressure / density


def _write_output(output_path, reference_path, reference, uniform, moments):
    """Write one single-location distribution using the FIDASIM HDF5 schema."""
    energy = uniform["energy"]
    pitch = uniform["pitch"]
    f_array = uniform["f_array"]
    density, parallel_temperature, perpendicular_temperature = moments

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, mode="w") as h5file:
        h5file.create_dataset("nenergy", data=energy.size)
        h5file.create_dataset("npitch", data=pitch.size)
        h5file.create_dataset("nr", data=1)
        h5file.create_dataset("nz", data=1)
        h5file.create_dataset("energy", data=energy)
        h5file.create_dataset("pitch", data=pitch)
        h5file.create_dataset("r", data=np.array([reference["selected_r"]]))
        h5file.create_dataset("z", data=np.array([reference["selected_z"]]))
        h5file.create_dataset("denf", data=np.array([[density]]))
        h5file.create_dataset("f", data=f_array.T[np.newaxis, np.newaxis, :, :])
        h5file.create_dataset("A", data=SPECIES_MASS_AMU[reference["species"]])

        moments_group = h5file.create_group("moments")
        moments_group.create_dataset("density", data=density)
        moments_group.create_dataset("parallel_temperature", data=parallel_temperature)
        moments_group.create_dataset(
            "perpendicular_temperature", data=perpendicular_temperature
        )

        h5file["energy"].attrs["units"] = "keV"
        h5file["pitch"].attrs["units"] = "dimensionless"
        h5file["r"].attrs["units"] = "cm"
        h5file["z"].attrs["units"] = "cm"
        h5file["denf"].attrs["units"] = "ions/cm^3"
        h5file["f"].attrs["units"] = "ions/(cm^3*keV*dP)"
        moments_group["density"].attrs["units"] = "ions/cm^3"
        moments_group["parallel_temperature"].attrs["units"] = "keV"
        moments_group["perpendicular_temperature"].attrs["units"] = "keV"
        h5file.attrs["source_reference"] = str(reference_path.resolve())
        h5file.attrs["species"] = reference["species"]


def run_conversion(config_path):
    """Convert every Stage 1 reference file and return generated output paths."""
    config = read_config(config_filename=config_path)
    input_config = config["input"]
    output_base = Path(config["save_data_block"]["output_filename"])
    references = _discover_reference_files(input_config["reference_config"])

    output_paths = []
    report_results = []
    for case_index, reference_path in enumerate(references, start=1):
        with h5py.File(reference_path, mode="r") as h5file:
            species = str(h5file.attrs["species"]).lower()
            reference = {
                "species": species,
                "selected_r": float(h5file["selected_r"][()]),
                "selected_z": float(h5file["selected_z"][()]),
            }
            reference_moments = (
                float(h5file["moments/density"][()]),
                float(h5file["moments/parallel_temperature"][()]),
                float(h5file["moments/perpendicular_temperature"][()]),
            )
            transformed = transform_to_nonrelativistic_energy_pitch(
                f_u_theta=h5file["f_u_theta"][:],
                u_bar=h5file["u_bar"][:],
                theta=h5file["theta"][:],
                u_norm=float(h5file["u_norm"][()]),
                species=species,
            )

        uniform = remap_to_uniform_grid(
            nonuniform_energy=transformed["energy"],
            nonuniform_pitch=transformed["pitch"],
            nonuniform_distribution=transformed["f_energy_pitch"],
            energy_upper_edge=transformed["energy_upper_edge"],
            nenergy=input_config["nenergy"],
            npitch=input_config["npitch"],
        )
        mass_grams = SPECIES_MASS_AMU[species] * 1.66053906660e-24
        moments = _calculate_output_moments(
            energy=uniform["energy"],
            pitch=uniform["pitch"],
            f_array=uniform["f_array"],
            mass_grams=mass_grams,
        )
        output_path = _indexed_path(output_base, case_index)
        _write_output(output_path, reference_path, reference, uniform, moments)
        print(f"Wrote converted file: {output_path}")
        output_paths.append(output_path)

        if input_config["plot_data"]:
            plot_path = plot_converted_distribution(
                input_path=output_path,
                plot_config=config["plot_data_block"],
            )
            print(f"Wrote plot: {plot_path}")

        report_results.append(
            {
                "reference_path": reference_path,
                "output_path": output_path,
                "selected_r": reference["selected_r"],
                "selected_z": reference["selected_z"],
                "reference_moments": reference_moments,
                "converted_moments": moments,
            }
        )

    report_path = output_base.parent / "conversion_moments.txt"
    write_moment_report(results=report_results, output_path=report_path)
    print(f"Wrote report: {report_path}")
    return output_paths
