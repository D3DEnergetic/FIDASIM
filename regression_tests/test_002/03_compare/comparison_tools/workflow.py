"""Coordinate the complete Test 002 moment-comparison workflow."""

from .comparison import compare_file_pair
from .config import discover_file_pairs, read_config
from .plotting import plot_moment_values, plot_relative_differences
from .reporting import write_report


def run_comparison(config_path):
    """Compare every reference-converted pair configured by Stage 2."""
    config = read_config(config_filename=config_path)
    compare_config = config["compare"]
    file_pairs = discover_file_pairs(run_config=compare_config["run_config"])

    results = []
    for case_index, pair in enumerate(file_pairs, start=1):
        print(f"Comparing case {case_index:03d}: {pair['reference'].name}")
        result = compare_file_pair(pair=pair, case_index=case_index)
        results.append(result)

    output_directory = compare_config["output_directory"]
    output_directory.mkdir(parents=True, exist_ok=True)
    write_report(results=results, output_directory=output_directory)

    if compare_config["generate_plot"]:
        plot_moment_values(
            results=results,
            output_directory=output_directory,
        )
        plot_relative_differences(
            results=results,
            output_directory=output_directory,
        )

    return results
