"""Coordinate the complete Test 004 ion-sink comparison workflow."""

from dataclasses import dataclass
from pathlib import Path

from regression_test_tools import print_config
from test_004_tools import read_distribution

from .analysis import compare_sinks
from .config import build_file_pairs, read_comparison_config
from .data import read_sink, validate_current_inputs, validate_pair
from .reporting import write_report


@dataclass
class ComparisonRun:
    """Outputs and overall status returned by one comparison workflow."""

    results: list
    report_filename: Path
    all_passed: bool


def run_comparison(config_path):
    """Compare every deterministic–Monte Carlo pair selected by Stage 3."""
    # Step 1: validate the Stage 3 interface and derive exact upstream files.
    config = read_comparison_config(config_path)
    print_config(config)
    compare_config = config["compare"]
    plot_config = config["plot_data_block"]
    inputs = build_file_pairs(
        compare_config["test_config"]
    )

    # Step 2: read, validate, and analyze the complete collection before
    # writing anything. A malformed late case therefore cannot leave a
    # partially updated report or plot collection.
    comparison_records = []
    for pair in inputs.file_pairs:
        print(f"Comparing case {pair.case_index:03d}: {pair.basename}")
        source_distribution = read_distribution(pair.source_distribution)
        deterministic = read_sink(
            pair.deterministic,
            implementation="deterministic",
        )
        monte_carlo = read_sink(
            pair.monte_carlo,
            implementation="monte_carlo",
        )
        validate_pair(deterministic, monte_carlo, pair)
        validate_current_inputs(
            deterministic=deterministic,
            monte_carlo=monte_carlo,
            pair=pair,
            source_distribution=source_distribution,
            test_config=inputs.test_config,
        )
        result = compare_sinks(
            pair=pair,
            deterministic=deterministic,
            monte_carlo=monte_carlo,
            relative_tolerance=compare_config[
                "rate_relative_tolerance"
            ],
            sigma_tolerance=compare_config["rate_sigma_tolerance"],
        )
        comparison_records.append(
            {
                "pair": pair,
                "deterministic": deterministic,
                "monte_carlo": monte_carlo,
                "result": result,
            }
        )
        print(
            f"  relative rate difference: "
            f"{100.0 * result.signed_relative_rate_difference:.2f}%"
        )
        print(
            f"  standardized difference: "
            f"{result.standardized_rate_difference:.2f}"
        )
        print(f"  status: {'PASS' if result.passed else 'FAIL'}")

    # Step 3: generate diagnostic plots only after every pair has passed
    # structural validation.
    inputs.output_directory.mkdir(parents=True, exist_ok=True)
    if compare_config["generate_plots"]:
        # Import Matplotlib only when the user has requested figures.
        from .plotting import (
            plot_distributions,
            plot_marginals,
            plot_rate_summary,
        )

        for record in comparison_records:
            plot_marginals(
                pair=record["pair"],
                deterministic=record["deterministic"],
                monte_carlo=record["monte_carlo"],
                result=record["result"],
                output_directory=inputs.output_directory,
                plot_config=plot_config,
            )
            plot_distributions(
                pair=record["pair"],
                deterministic=record["deterministic"],
                monte_carlo=record["monte_carlo"],
                result=record["result"],
                output_directory=inputs.output_directory,
                plot_config=plot_config,
            )

        plot_rate_summary(
            results=[
                record["result"] for record in comparison_records
            ],
            output_directory=inputs.output_directory,
            relative_tolerance=compare_config[
                "rate_relative_tolerance"
            ],
        )

    # Step 4: always write one text report containing the numerical results.
    results = [record["result"] for record in comparison_records]
    report_filename = write_report(
        results=results,
        comparison_config=compare_config,
        output_directory=inputs.output_directory,
    )
    return ComparisonRun(
        results=results,
        report_filename=report_filename,
        all_passed=all(result.passed for result in results),
    )
