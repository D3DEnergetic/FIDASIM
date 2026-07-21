"""Coordinate Stage 1 reference-file extraction."""

from pathlib import Path

from regression_test_tools import ConfigError

from .config import read_config
from .moments import calculate_cql3d_moments, write_physical_moments
from .plotting import plot_reference_distribution
from .readers import extract_cql3d_f4d_location
from .reporting import write_moment_report


def _indexed_output_path(base_path, output_index):
    """Create the output filename for one selected spatial location.

    The configuration provides one base filename for all generated files. This
    function inserts a three-digit case number before the file extension so
    that every selected location receives a unique filename. For example,
    ``cql3d_f4d.h5`` and index ``1`` produce ``cql3d_f4d_001.h5``.

    Args:
        base_path (Path): Configured output path without a case number.
        output_index (int): One-based case number for the selected location.

    Returns:
        Path: Output path containing the three-digit case number.
    """
    return base_path.with_name(
        f"{base_path.stem}_{output_index:03d}{base_path.suffix}"
    )


def generate_reference_files(config_path):
    """Generate one compact CQL3D reference fixture per requested location.

    Args:
        config_path (str or Path): Stage 1 namelist configuration file.

    Returns:
        list[Path]: Paths of the generated single-location HDF5 files.

    Raises:
        ConfigError: If saving is disabled or the input format is unsupported.
    """
    config = read_config(config_filename=config_path)
    input_config = config["input"]
    output_base = config["save_data_block"]["output_filename"]

    if not input_config["save_data"]:
        raise ConfigError(
            "save_data must be true while Stage 1 only performs file extraction."
        )

    particle = {
        "species": input_config["species"],
        "atomic_number": input_config["atomic_number"],
        "mass_number": input_config["mass_number"],
        "charge_state": input_config["charge_state"],
    }

    generated_paths = []
    moment_results = []
    number_of_locations = len(input_config["r_locations"])
    for location_index in range(number_of_locations):
        requested_r = input_config["r_locations"][location_index]
        requested_z = input_config["z_locations"][location_index]
        output_index = location_index + 1

        output_path = _indexed_output_path(
            base_path=Path(output_base),
            output_index=output_index,
        )

        # Select the nearest spatial cell and write a compact HDF5 file containing the
        # requested location and the CQL3D variables.
        if input_config["input_file_type"] == "cql3d_f4d":
            selection = extract_cql3d_f4d_location(
                source_path=input_config["input_filename"],
                output_path=output_path,
                requested_r=requested_r,
                requested_z=requested_z,
                particle=particle,
            )
        else:
            raise ConfigError(
                f"Unsupported input_file_type: {input_config['input_file_type']}"
            )

        # Print the requested and selected locations for each generated file:
        print(
            f"Case {output_index:03d}: requested "
            f"(R, Z) = ({selection['requested_r']:.6g}, "
            f"{selection['requested_z']:.6g}) cm; selected "
            f"({selection['selected_r']:.6g}, "
            f"{selection['selected_z']:.6g}) cm"
        )
        generated_paths.append(output_path)

        # Calculate the physical moments for the selected location and write them to the HDF5 file: 
        moments = calculate_cql3d_moments(
            input_path=output_path,
            species=input_config["species"],
        )
        write_physical_moments(
            output_path=output_path,
            moments=moments,
        )

        # Plot the selected distribution after its moments have been stored:
        if input_config["plot_data"]:
            plot_path = plot_reference_distribution(
                input_path=output_path,
                plot_config=config["plot_data_block"],
            )
            print(f"Wrote plot: {plot_path}")

        # Append the results to the moment report for all selected locations:
        moment_results.append(
            {
                "path": output_path,
                "selected_r": selection["selected_r"],
                "selected_z": selection["selected_z"],
                "moments": moments,
            }
        )

    # Write a report of the physical moments for all selected locations:
    report_path = output_base.parent / "reference_moments.txt"
    write_moment_report(results=moment_results, output_path=report_path)
    print(f"Wrote report: {report_path}")

    return generated_paths
