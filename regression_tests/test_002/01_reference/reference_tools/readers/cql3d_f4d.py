"""Extract single-location fixtures from CQL3D F4D NetCDF files."""

from pathlib import Path

import h5py
import numpy as np
from netCDF4 import Dataset

from regression_test_tools import ConfigError

from ..particle_metadata import write_particle_metadata


R_DIMENSION = "dim_nr_f4d"
Z_DIMENSION = "dim_nz_f4d"


def _validate_schema(source):
    """Check that a NetCDF dataset contains the CQL3D data used here."""
    required_variables = [
        "mnemonic",
        "version",
        "enorm",
        "k",
        "vnorm",
        "f4dr",
        "f4dz",
        "f4dv",
        "f4dt",
        "f4ddv",
        "f4ddt",
        "f4d",
    ]

    for variable_name in required_variables:
        if variable_name not in source.variables:
            raise ConfigError(
                f"CQL3D F4D file is missing required variable: {variable_name}"
            )

    expected_f4d_dimensions = (
        "dim_nt_f4d",
        "dim_nv_f4d",
        Z_DIMENSION,
        R_DIMENSION,
    )
    if source.variables["f4d"].dimensions != expected_f4d_dimensions:
        raise ConfigError(
            "CQL3D variable f4d must have dimensions "
            "(dim_nt_f4d, dim_nv_f4d, dim_nz_f4d, dim_nr_f4d)."
        )


def _nearest_grid_index(grid, requested_location):
    """Return the index and value of the grid point nearest a location."""
    index = int(np.argmin(np.abs(grid - requested_location)))
    return index, float(grid[index])


def _read_text(source_variable):
    """Convert a fixed-width NetCDF character array to a Python string."""
    characters = source_variable[:]
    if np.ma.isMaskedArray(characters):
        characters = characters.filled(b" ")
    return characters.tobytes().decode("utf-8").rstrip("\x00 ")


def _write_reference_data(source, destination, r_index, z_index):
    """Write the compact HDF5 data required by this regression test."""
    string_type = h5py.string_dtype(encoding="utf-8")

    destination.create_dataset(
        "mnemonic", data=_read_text(source.variables["mnemonic"]), dtype=string_type
    )
    destination.create_dataset(
        "version", data=_read_text(source.variables["version"]), dtype=string_type
    )
    destination.create_dataset("enorm", data=source.variables["enorm"].getValue())
    destination.create_dataset("k", data=source.variables["k"].getValue())

    destination.create_dataset("u_norm", data=source.variables["vnorm"].getValue())
    destination.create_dataset("selected_r", data=source.variables["f4dr"][r_index])
    destination.create_dataset("selected_z", data=source.variables["f4dz"][z_index])
    destination.create_dataset("u_bar", data=source.variables["f4dv"][:])
    destination.create_dataset("theta", data=source.variables["f4dt"][:])
    destination.create_dataset("u_bar_weight", data=source.variables["f4ddv"][:])
    destination.create_dataset("theta_weight", data=source.variables["f4ddt"][:])
    destination.create_dataset(
        "f_u_theta",
        data=source.variables["f4d"][:, :, z_index, r_index],
    )

    destination["enorm"].attrs["units"] = "keV"
    destination["u_norm"].attrs["units"] = "cm/s"
    destination["selected_r"].attrs["units"] = "cm"
    destination["selected_z"].attrs["units"] = "cm"
    destination["theta"].attrs["units"] = "radians"
    destination["u_bar_weight"].attrs["description"] = "u_bar^2 du_bar"
    destination["theta_weight"].attrs["description"] = "2*pi*sin(theta) dtheta"
    destination["f_u_theta"].attrs["axis_order"] = "theta, u_bar"
    destination["f_u_theta"].attrs["units"] = (
        "ions*u_norm**3/(cm**3*(cm/sec)**3)"
    )


def extract_cql3d_f4d_location(
    source_path,
    output_path,
    requested_r,
    requested_z,
    particle,
):
    """Extract the nearest CQL3D spatial cell into a compact HDF5 file.

    Args:
        source_path (str or Path): Parent CQL3D F4D NetCDF file.
        output_path (str or Path): Single-location HDF5 file to create.
        requested_r (float): Requested major radius in centimetres.
        requested_z (float): Requested axial location in centimetres.
        particle (dict): Species, atomic number, mass number, and charge state.

    Returns:
        dict: Requested locations, selected indices, and actual grid locations.

    Raises:
        ConfigError: If the source schema is invalid or the file cannot be read.
    """
    source_path = Path(source_path)
    output_path = Path(output_path)

    try:
        with Dataset(source_path, mode="r") as source:
            
            # Check that the source file has the expected CQL3D F4D schema before extracting:
            _validate_schema(source=source)

            r_grid = np.asarray(source.variables["f4dr"][:], dtype=float)
            z_grid = np.asarray(source.variables["f4dz"][:], dtype=float)
            r_index, selected_r = _nearest_grid_index(
                grid=r_grid,
                requested_location=requested_r,
            )
            z_index, selected_z = _nearest_grid_index(
                grid=z_grid,
                requested_location=requested_z,
            )

            output_path.parent.mkdir(parents=True, exist_ok=True)
            with h5py.File(output_path, mode="w") as destination:
                _write_reference_data(
                    source=source,
                    destination=destination,
                    r_index=r_index,
                    z_index=z_index,
                )

                # Add regression-test provenance without changing the CQL3D variables.
                destination.attrs["source_file"] = str(source_path.resolve())
                destination.attrs["requested_r_cm"] = requested_r
                destination.attrs["requested_z_cm"] = requested_z
                destination.attrs["selected_r_index"] = r_index
                destination.attrs["selected_z_index"] = z_index
                write_particle_metadata(
                    h5file=destination,
                    particle=particle,
                )
    except (OSError, RuntimeError) as error:
        raise ConfigError(f"Could not extract CQL3D F4D data: {error}") from error

    return {
        "requested_r": requested_r,
        "requested_z": requested_z,
        "r_index": r_index,
        "z_index": z_index,
        "selected_r": selected_r,
        "selected_z": selected_z,
    }
