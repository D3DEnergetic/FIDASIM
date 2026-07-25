#!/usr/bin/env python3
"""Normalize the Test 004 Monte Carlo configuration for Fortran."""

from pathlib import Path
import sys

test_directory = Path(__file__).resolve().parents[1]
regression_tests_directory = test_directory.parent
sys.path.insert(0, str(test_directory))
sys.path.insert(0, str(regression_tests_directory))

import f90nml

from regression_test_tools import ConfigError
from test_004_tools import discover_distributions, read_monte_carlo_config

# The maximum number of input distributions that can be processed in one Monte Carlo run:
MAX_CASES = 256

# From the src, the maximum length of a FIDASIM character field is 200 characters:
FIDASIM_CHARLIM = 200


def _require_fidasim_length(value, label):
    """Require a string to fit in FIDASIM's fixed Fortran character fields.

    FIDASIM declares production values such as ``inputs%runid``,
    ``inputs%result_dir``, ``inputs%tables_file``, and the filename assembled
    by ``write_sink_profile`` as ``character(charlim)``, where ``charlim`` is
    200. Rejecting longer normalized values here prevents silent truncation
    when Python strings cross the Fortran interface.
    """
    if len(value) > FIDASIM_CHARLIM:
        raise ConfigError(
            f"{label} exceeds FIDASIM's {FIDASIM_CHARLIM}-character "
            "limit."
        )


def _runid(input_path):
    """Return the FIDASIM run ID that produces the shared sink basename."""
    return f"{Path(input_path).stem}_ion"


def normalize_and_flatten_config(config_filename):
    """Normalize Test 004 inputs and flatten them into one Fortran block."""

    # Read the Monte Carlo configuration:
    config = read_monte_carlo_config(config_filename)

    # Discover the ordered Test 002 distribution cases and their species parameters:
    provenance = discover_distributions(
        config["test_case"]["input_distribution_config"]
    )

    # Make sure the number of cases is less than MAX_CASES:
    cases = provenance["cases"]
    n_cases = len(cases)
    if n_cases > MAX_CASES:
        raise ConfigError(
            f"No more than {MAX_CASES} input distributions are supported."
        )

    # Get paths to the input distribution files and their run IDs::
    distribution_files = []
    runids = []
    for case in cases:
        input_path = case["input_path"]
        distribution_file = str(input_path)
        runid = _runid(input_path)

        distribution_files.append(distribution_file)
        runids.append(runid)

    # Verify that every run ID is unique:
    if len(set(runids)) != n_cases:
        raise ConfigError(
            "Every input distribution must produce a unique run ID."
        )

    # Verify that every run ID and the tables filename fit in FIDASIM's character limit:
    for runid in runids:
        _require_fidasim_length(
            value=runid,
            label="runid",
        )

    monte_carlo = config["monte_carlo"]
    output_directory = config["save_data_block"]["monte_carlo_directory"]
    if monte_carlo["save_data"]:
        output_directory_value = str(output_directory)
    else:
        output_directory_value = ""

    tables_filename = str(config["test_case"]["tables_filename"])
    _require_fidasim_length(
        value=tables_filename,
        label="tables_filename",
    )
    if monte_carlo["save_data"]:
        _require_fidasim_length(
            value=output_directory_value,
            label="output_directory",
        )
        for runid in runids:
            sink_filename = str(
                Path(output_directory_value) / f"{runid}_sink.h5"
            )
            _require_fidasim_length(
                value=sink_filename,
                label="sink filename",
            )

    level_decay = config["neutrals"]["level_decay"]
    if level_decay is None:
        level_decay = 0.0

    return {
        "run_test": {
            "n_cases": n_cases,
            "distribution_files": distribution_files,
            "runids": runids,
            "tables_filename": tables_filename,
            "test_config": str(config["test_case"]["config_path"]),
            "input_distribution_config": str(provenance["run_config"]),
            "output_directory": output_directory_value,
            "case_comment": config["test_case"]["comment"],
            "implementation_comment": monte_carlo["comment"],
            "n_markers": monte_carlo["n_markers"],
            "reservoir_size": monte_carlo["reservoir_size"],
            "seed": monte_carlo["seed"],
            "save_data": monte_carlo["save_data"],
            "neutral_density": config["neutrals"]["density"],
            "neutral_energy": config["neutrals"]["energy"],
            "injection_angle": config["neutrals"]["injection_angle"],
            "level_split_method": config["neutrals"][
                "level_split_method"
            ],
            "level_decay": level_decay,
        }
    }


def _require_safe_output_path(config, output_path):
    """Prevent the generated namelist from overwriting an input asset.

    ``test_config`` is the unified, user-facing Test 004 configuration passed
    to this normalizer, such as ``regression_tests/test_004/input_config_A.nml``.
    ``input_distribution_config`` is the Test 002 Stage 2 configuration used
    to discover the smooth distributions.

    The protected inputs are:

    - the unified Test 004 configuration;
    - the Test 002 Stage 2 distribution configuration;
    - the FIDASIM atomic-tables HDF5 file;
    - every discovered Test 002 distribution HDF5 file.

    """
    run_test = config["run_test"]

    test_004_config_path = Path(run_test["test_config"]).resolve()
    test_002_config_path = Path(
        run_test["input_distribution_config"]
    ).resolve()
    atomic_tables_path = Path(run_test["tables_filename"]).resolve()

    protected_paths = {
        test_004_config_path,
        test_002_config_path,
        atomic_tables_path,
    }

    distribution_files = run_test["distribution_files"]
    for distribution_file in distribution_files:
        distribution_path = Path(distribution_file).resolve()
        protected_paths.add(distribution_path)

    if output_path in protected_paths:
        raise ConfigError(
            "The normalized configuration must not overwrite an input file."
        )


def write_normalized_config(config, output_filename):
    """Write the flattened configuration."""

    # Get the output path:
    output_path = Path(output_filename).expanduser().resolve()

    # Protect against overwriting any configured input:
    _require_safe_output_path(config=config, output_path=output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create the output directory if saving is enabled:
    run_test = config["run_test"]
    if run_test["save_data"]:
        Path(run_test["output_directory"]).mkdir(parents=True, exist_ok=True)

    # Write the normalized configuration to the output path:
    f90nml.write(config, output_path, force=True)

    return output_path


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: normalize_config.py <input_config.nml> "
            "<normalized_config.nml>"
        )

    try:
        config = normalize_and_flatten_config(
            config_filename=sys.argv[1]
        )
        output_path = write_normalized_config(
            config=config,
            output_filename=sys.argv[2],
        )
    except ConfigError as error:
        raise SystemExit(f"Configuration error: {error}") from None
    except OSError as error:
        raise SystemExit(
            f"Cannot write normalized configuration: {error}"
        ) from None

    print(f"Wrote normalized configuration: {output_path}")


if __name__ == "__main__":
    main()
