"""Coordinate the complete reference-sampled comparison workflow."""

from .analysis import calculate_physical_moments, compare_moments
from .config import build_file_pairs, read_config
from .data import read_distribution, validate_pair
from .plotting import (
    plot_distributions,
    plot_marginals,
    plot_moment_relative_errors,
)
from .reporting import write_report


def run_comparison(config_path):
    """Run every comparison and return the overall regression status."""

    # Step 1: read configuration and recover the exact Stage 2 file list.
    config = read_config(
        config_filename=config_path,
    )
    compare_config = config["compare"]
    plot_config = config["plot_data_block"]

    file_pairs = build_file_pairs(
        sampling_config_file=compare_config["sampling_config_file"],
    )
    output_directory = compare_config["output_directory"]
    output_directory.mkdir(parents=True, exist_ok=True)

    results = []
    for case_index, pair in enumerate(file_pairs, start=1):
        print(
            f"Comparing case {case_index:03d}: "
            f"{pair.reference.name}"
        )

        # Step 2: read and validate the shared HDF5 data contract.
        reference = read_distribution(
            filename=pair.reference,
        )
        sampled = read_distribution(
            filename=pair.sampled,
        )
        validate_pair(
            reference=reference,
            sampled=sampled,
            pair=pair,
        )

        # Step 3: calculate the physical quantities used for comparison.
        reference_moments = calculate_physical_moments(
            distribution=reference,
        )
        sampled_moments = calculate_physical_moments(
            distribution=sampled,
        )
        result = compare_moments(
            case_index=case_index,
            basename=pair.reference.name,
            reference_moments=reference_moments,
            sampled_moments=sampled_moments,
            relative_tolerance=compare_config[
                "moment_relative_tolerance"
            ],
        )
        results.append(result)
        print(f"  status: {'PASS' if result.passed else 'FAIL'}")

        # Step 4: generate the two requested diagnostic figures.
        if compare_config["generate_plots"]:
            plot_marginals(
                pair=pair,
                reference=reference,
                sampled=sampled,
                result=result,
                output_directory=output_directory,
                plot_config=plot_config,
            )
            plot_distributions(
                pair=pair,
                reference=reference,
                sampled=sampled,
                output_directory=output_directory,
                plot_config=plot_config,
            )

    # Step 5: generate the collection-level acceptance plot and write one
    # report after every pair has been processed.
    if compare_config["generate_plots"]:
        plot_moment_relative_errors(
            results=results,
            relative_tolerance=compare_config[
                "moment_relative_tolerance"
            ],
            output_directory=output_directory,
        )

    write_report(
        results=results,
        relative_tolerance=compare_config["moment_relative_tolerance"],
        output_directory=output_directory,
    )

    all_passed = all(result.passed for result in results)
    print()
    print(f"Test 003 comparison status: {'PASS' if all_passed else 'FAIL'}")
    return all_passed
