"""Coordinate the complete Test 002 moment-comparison workflow."""

from .comparison import compare_file_pair
from .config import discover_file_pairs, read_config
from .plotting import plot_moment_values, plot_relative_differences
from .reporting import write_report


def run_comparison(config_path):
    """Compare all configured file pairs and return the overall pass status."""
    config = read_config(config_filename=config_path)
    compare_config = config["compare"]
    file_pairs = discover_file_pairs(run_config=compare_config["run_config"])

    results = []
    for case_index, pair in enumerate(file_pairs, start=1):
        print(f"Comparing case {case_index:03d}: {pair['reference'].name}")
        result = compare_file_pair(
            pair=pair,
            case_index=case_index,
            relative_tolerance=compare_config[
                "moment_relative_tolerance"
            ],
        )
        results.append(result)
        print(f"  status: {'PASS' if result['passed'] else 'FAIL'}")

    output_directory = compare_config["output_directory"]
    output_directory.mkdir(parents=True, exist_ok=True)
    write_report(
        results=results,
        relative_tolerance=compare_config["moment_relative_tolerance"],
        output_directory=output_directory,
    )

    if compare_config["generate_plot"]:
        plot_moment_values(
            results=results,
            output_directory=output_directory,
        )
        plot_relative_differences(
            results=results,
            relative_tolerance=compare_config[
                "moment_relative_tolerance"
            ],
            output_directory=output_directory,
        )

    all_passed = all(result["passed"] for result in results)
    print()
    print(f"Test 002 comparison status: {'PASS' if all_passed else 'FAIL'}")
    return all_passed
