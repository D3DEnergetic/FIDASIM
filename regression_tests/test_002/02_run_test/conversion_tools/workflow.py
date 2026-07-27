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
    ATOMIC_MASS_GRAMS,
    ERG_PER_KEV,
    transform_to_nonrelativistic_energy_pitch,
)


SPECIES_DATASETS = (
    "species",
    "atomic_number",
    "mass_number",
    "charge_state",
    "A",
)


def _read_species_metadata(h5file, label):
    """Read and normalize species metadata from an open HDF5 file.

    Args:
        h5file (h5py.File or h5py.Group): Open object containing the required
            scalar species datasets.
        label (str or Path): Diagnostic name used to identify the source in
            validation errors.

    Returns:
        dict: Species values and their original HDF5 dataset attributes with
        the following schema::

            {
                "species": str,
                "atomic_number": int,
                "mass_number": int,
                "charge_state": int,
                "A": float,
                "attributes": {
                    "species": dict,
                    "atomic_number": dict,
                    "mass_number": dict,
                    "charge_state": dict,
                    "A": dict,
                },
            }

        ``species`` is decoded, stripped, and converted to lowercase. Each
        nested attribute dictionary contains the attributes copied from its
        corresponding source dataset.

    Raises:
        ConfigError: If any required species dataset is missing.
    """
    missing_datasets = []
    for dataset_name in SPECIES_DATASETS:
        if dataset_name not in h5file:
            missing_datasets.append(dataset_name)

    if missing_datasets:
        missing_names = ", ".join(missing_datasets)
        raise ConfigError(
            f"{label} is missing species datasets: {missing_names}"
        )

    # Depending on the HDF5 string type, h5py may return raw bytes instead of
    # a Python string, for example b"D" instead of "D". Decode bytes so
    # downstream metadata always uses str.
    species = h5file["species"][()]
    if isinstance(species, bytes):
        species = species.decode("utf-8")

    metadata = {
        "species": str(species).strip().lower(),
        "atomic_number": int(h5file["atomic_number"][()]),
        "mass_number": int(h5file["mass_number"][()]),
        "charge_state": int(h5file["charge_state"][()]),
        "A": float(h5file["A"][()]),
    }
    metadata["attributes"] = {}
    for dataset_name in SPECIES_DATASETS:
        source_attributes = dict(h5file[dataset_name].attrs)
        metadata["attributes"][dataset_name] = source_attributes

    return metadata


def _write_species_metadata(h5file, species_metadata):
    """Copy scalar species metadata into one Stage 2 output file."""
    string_type = h5py.string_dtype(encoding="utf-8")
    for dataset_name in SPECIES_DATASETS:
        if dataset_name == "species":
            output_dataset = h5file.create_dataset(
                dataset_name,
                data=species_metadata[dataset_name],
                dtype=string_type,
            )
        else:
            output_dataset = h5file.create_dataset(
                dataset_name,
                data=species_metadata[dataset_name],
            )
        source_attributes = species_metadata["attributes"][dataset_name]
        for attribute_name, value in source_attributes.items():
            output_dataset.attrs[attribute_name] = value


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


def _calculate_output_moments(
    energy,
    pitch,
    f_pitch_energy,
    mass_grams,
):
    """Calculate density and relativistic pressure moments on the E-P grid."""
    denergy = energy[1] - energy[0]
    dpitch = pitch[1] - pitch[0]
    cell_area = denergy * dpitch
    density = float(np.sum(f_pitch_energy) * cell_area)

    energy_2d = energy[np.newaxis, :]
    pitch_2d = pitch[:, np.newaxis]
    u_squared = 2.0 * ERG_PER_KEV * energy_2d / mass_grams
    speed_of_light = 2.99792458e10
    gamma = np.sqrt(1.0 + u_squared / speed_of_light**2)
    pressure_energy = 2.0 * energy_2d / gamma

    parallel_pressure = float(
        np.sum(
            pressure_energy * pitch_2d**2 * f_pitch_energy
        )
        * cell_area
    )
    perpendicular_pressure = 0.5 * float(
        np.sum(
            pressure_energy
            * (1.0 - pitch_2d**2)
            * f_pitch_energy
        )
        * cell_area
    )
    return density, parallel_pressure / density, perpendicular_pressure / density


def _write_fidasim_distribution(
    output_path,
    reference_path,
    distribution_metadata,
    uniform_distribution,
    moments,
):
    """Write one energy-pitch distribution using the FIDASIM HDF5 schema.

    Args:
        output_path (Path): Destination HDF5 file.
        reference_path (Path): Stage 1 file from which the distribution came.
        distribution_metadata (dict): Species parameters and selected spatial
            location associated with the distribution.
        uniform_distribution (dict): Distribution sampled on the uniform
            energy-pitch grid. It contains ``energy``, ``pitch``, and
            ``f_pitch_energy``.
        moments (tuple): Density, parallel temperature, and perpendicular
            temperature calculated from the uniform distribution.
    """
    energy = uniform_distribution["energy"]
    pitch = uniform_distribution["pitch"]
    f_pitch_energy = uniform_distribution["f_pitch_energy"]
    density, parallel_temperature, perpendicular_temperature = moments

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_path, mode="w") as h5file:
        h5file.create_dataset("nenergy", data=energy.size)
        h5file.create_dataset("npitch", data=pitch.size)
        h5file.create_dataset("nr", data=1)
        h5file.create_dataset("nz", data=1)
        h5file.create_dataset("energy", data=energy)
        h5file.create_dataset("pitch", data=pitch)
        h5file.create_dataset(
            "r",
            data=np.array([distribution_metadata["selected_r"]]),
        )
        h5file.create_dataset(
            "z",
            data=np.array([distribution_metadata["selected_z"]]),
        )
        h5file.create_dataset("denf", data=np.array([[density]]))

        # The distribution represents one spatial location. Add singleton z
        # and r axes to the existing (pitch, energy) array. The stored dataset
        # therefore has the h5py-visible FIDASIM shape
        # (z, r, pitch, energy). The Fortran HDF5 interface reverses these
        # file dimensions when reading them and presents the FIDASIM code with
        # its logical shape (energy, pitch, r, z).
        stored_distribution = f_pitch_energy[
            np.newaxis,
            np.newaxis,
            :,
            :,
        ]
        h5file.create_dataset("f", data=stored_distribution)
        _write_species_metadata(
            h5file=h5file,
            species_metadata=distribution_metadata,
        )

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


def convert_reference_distributions_to_fidasim(config_path):
    """Convert the complete Stage 1 collection into FIDASIM distributions.

    The Stage 2 configuration points to the Stage 1 configuration that created
    the indexed single-location CQL3D reference files. This function discovers
    that collection and, for each reference:

    1. Reads the ``f(u, theta)`` distribution, location, species metadata,
       and reference moments.
    2. Transforms the distribution into nonrelativistic energy-pitch
       coordinates.
    3. Remaps it onto the configured uniform FIDASIM energy-pitch grid.
    4. Recalculates the density and directional temperature moments.
    5. Writes a self-describing HDF5 file using the native FIDASIM
       distribution schema and optionally creates its diagnostic plot.

    After all cases have been converted, the function writes one collection
    report comparing the Stage 1 and Stage 2 moments.

    Args:
        config_path (str or Path): Stage 2 namelist configuration file.

    Raises:
        ConfigError: If the configuration, reference collection, or required
            HDF5 data does not satisfy the Stage 2 input contract.
    """
    config = read_config(config_filename=config_path)
    input_config = config["input"]
    output_base = Path(config["save_data_block"]["output_filename"])
    references = _discover_reference_files(input_config["reference_config"])

    report_results = []
    for case_index, reference_path in enumerate(references, start=1):
        with h5py.File(reference_path, mode="r") as h5file:

            # Read the species parameters associated with this distribution.
            distribution_metadata = _read_species_metadata(
                h5file=h5file,
                label=reference_path,
            )

            selected_r = float(h5file["selected_r"][()])
            selected_z = float(h5file["selected_z"][()])

            # Combine the species parameters and spatial location into the
            # metadata describing this single-location distribution.
            distribution_metadata["selected_r"] = selected_r
            distribution_metadata["selected_z"] = selected_z

            reference_moments = (
                float(h5file["moments/density"][()]),
                float(h5file["moments/parallel_temperature"][()]),
                float(h5file["moments/perpendicular_temperature"][()]),
            )

            f_u_theta = h5file["f_u_theta"][:]
            u_bar = h5file["u_bar"][:]
            theta = h5file["theta"][:]
            u_norm = float(h5file["u_norm"][()])

        # u_bar and theta from the CQL3D dataset are uniform grids.
        # They can be convered into non-uniform E and P grids

        # Convert the CQL3D f(u, theta) distribution into F(E, P) at the
        # corresponding nonuniform energy-pitch coordinates. The result
        # retains (pitch, energy) axis order; interpolation onto the
        # uniform grid occurs next.
        transformed = transform_to_nonrelativistic_energy_pitch(
            f_u_theta=f_u_theta,
            u_bar=u_bar,
            theta=theta,
            u_norm=u_norm,
            mass_amu=distribution_metadata["A"],
        )

        # Interpolate F(E(u), P(theta)) from the nonuniform energy-pitch
        # coordinates produced by transforming the source u and theta grids
        # onto the configured uniform energy and pitch cell centers. The
        # returned distribution retains (pitch, energy) axis order.
        uniform_distribution = remap_to_uniform_grid(
            nonuniform_energy=transformed["energy"],
            nonuniform_pitch=transformed["pitch"],
            nonuniform_distribution=transformed["f_energy_pitch"],
            energy_upper_edge=transformed["energy_upper_edge"],
            nenergy=input_config["nenergy"],
            npitch=input_config["npitch"],
        )
        mass_grams = distribution_metadata["A"] * ATOMIC_MASS_GRAMS
        moments = _calculate_output_moments(
            energy=uniform_distribution["energy"],
            pitch=uniform_distribution["pitch"],
            f_pitch_energy=uniform_distribution["f_pitch_energy"],
            mass_grams=mass_grams,
        )
        output_path = _indexed_path(output_base, case_index)
        _write_fidasim_distribution(
            output_path=output_path,
            reference_path=reference_path,
            distribution_metadata=distribution_metadata,
            uniform_distribution=uniform_distribution,
            moments=moments,
        )
        print(f"Wrote converted file: {output_path}")

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
                "selected_r": distribution_metadata["selected_r"],
                "selected_z": distribution_metadata["selected_z"],
                "reference_moments": reference_moments,
                "converted_moments": moments,
            }
        )

    report_path = output_base.parent / "conversion_moments.txt"
    write_moment_report(results=report_results, output_path=report_path)
    print(f"Wrote report: {report_path}")
