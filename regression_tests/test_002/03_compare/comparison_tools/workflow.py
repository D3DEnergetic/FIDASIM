"""Coordinate the complete reference-sampled comparison workflow."""

from .analysis import calculate_physical_moments, compare_moments
from .config import build_file_pairs, read_config
from .data import read_distribution, validate_pair
from .plotting import plot_distributions, plot_marginals
from .reporting import write_report


def run_comparison(config_path):
    """Run every comparison configured by the Stage 3 namelist."""

    # Step 1: read configuration and recover the exact Stage 2 file list.
    config = read_config(config_path)
    file_pairs = build_file_pairs(config.sampling_config_file)
    config.output_directory.mkdir(parents=True, exist_ok=True)

    results = []
    for pair in file_pairs:
        print(f"Comparing {pair.reference.name}")

        # Step 2: read and validate the shared HDF5 data contract.
        reference = read_distribution(pair.reference)
        sampled = read_distribution(pair.sampled)
        validate_pair(reference, sampled, pair)

        # Step 3: calculate the physical quantities used for comparison.
        reference_moments = calculate_physical_moments(reference)
        sampled_moments = calculate_physical_moments(sampled)
        result = compare_moments(
            pair.reference.name,
            reference_moments,
            sampled_moments,
        )
        results.append(result)

        # Step 4: generate the two requested diagnostic figures.
        if config.generate_plots:
            plot_marginals(pair, reference, sampled, result, config)
            plot_distributions(pair, reference, sampled, config)

    # Step 5: write one report after every pair has been processed.
    write_report(results, config.output_directory)
