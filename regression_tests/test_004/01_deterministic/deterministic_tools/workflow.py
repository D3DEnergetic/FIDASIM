"""Orchestrate the Test 004 deterministic calculation."""

from pathlib import Path

from regression_test_tools import ConfigError, print_config

from .atomic import read_charge_exchange_table
from test_004_tools import discover_distributions, read_deterministic_config

from .calculation import calculate_deterministic
from .data import read_distribution
from .output import write_deterministic
from .plotting import plot_deterministic


def run_deterministic_workflow(config_filename):
    """Process every smooth converted distribution selected through Test 002."""
    config = read_deterministic_config(config_filename)
    print_config(config)

    input_distribution_config = config["test_case"][
        "input_distribution_config"
    ]

    # Get the ordered Test 002 cases and their shared species parameters.
    provenance = discover_distributions(input_distribution_config)

    table = read_charge_exchange_table(config["test_case"]["tables_filename"])
    output_directory = config["save_data_block"]["deterministic_directory"]

    results = []
    for case in provenance["cases"]:
        distribution = read_distribution(case["input_path"])
        particle = provenance["particle"]
        distribution_particle = {
            "species": distribution.species,
            "atomic_number": distribution.atomic_number,
            "mass_number": distribution.mass_number,
            "charge_state": distribution.charge_state,
            "A": distribution.atomic_mass,
        }
        if distribution_particle != particle:
            raise ConfigError(
                f"{case['input_path']}: species parameters differ from "
                "the discovered collection metadata."
            )
        if (
            distribution.selected_r != case["selected_r"]
            or distribution.selected_z != case["selected_z"]
        ):
            raise ConfigError(
                f"{case['input_path']}: selected coordinates differ from "
                "the discovered case metadata."
            )

        result = calculate_deterministic(
            distribution=distribution,
            table=table,
            neutral_config=config["neutrals"],
            n_gyro=config["deterministic"]["n_gyro"],
        )
        stem = f"{Path(case['input_path']).stem}_ion_sink"
        hdf5_path = None
        plot_path = None
        if config["deterministic"]["save_data"]:
            hdf5_path = write_deterministic(
                filename=output_directory / f"{stem}.h5",
                distribution=distribution,
                result=result,
                config=config,
                case=case,
                provenance=provenance,
            )
        if config["deterministic"]["plot_data"]:
            plot_path = plot_deterministic(
                filename=output_directory / f"{stem}.png",
                distribution=distribution,
                result=result,
                plot_config=config["plot_data_block"],
            )
        print(
            f"{case['index']:03d}: {case['input_path'].name} "
            f"R={result.total_rate:.10e} ions/(cm^3*s)"
        )
        results.append(
            {
                "case": case,
                "result": result,
                "hdf5_path": hdf5_path,
                "plot_path": plot_path,
            }
        )
    return results
