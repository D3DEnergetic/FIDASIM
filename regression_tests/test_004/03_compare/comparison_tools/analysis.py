"""Calculate scalar and marginal Test 004 comparison metrics."""

from dataclasses import dataclass

import numpy as np


@dataclass
class ComparisonResult:
    """Comparison metrics for one deterministic–Monte Carlo pair."""

    case_index: int
    basename: str
    deterministic_rate: float
    monte_carlo_rate: float
    monte_carlo_standard_error: float
    signed_rate_difference: float
    absolute_rate_difference: float
    signed_relative_rate_difference: float
    absolute_relative_rate_difference: float
    standardized_rate_difference: float
    energy_marginal_l1_difference: float
    pitch_marginal_l1_difference: float
    relative_rate_passed: bool
    sigma_rate_passed: bool
    passed: bool


def _normalized_l1_difference(reference, candidate, spacing, total_rate):
    """Return the integrated absolute difference normalized by reference rate."""
    integrated_difference = float(
        np.sum(np.abs(candidate - reference)) * spacing
    )
    return integrated_difference / total_rate


def compare_sinks(
    pair,
    deterministic,
    monte_carlo,
    relative_tolerance,
    sigma_tolerance,
):
    """Calculate rate acceptance metrics and diagnostic marginal distances."""
    signed_difference = (
        monte_carlo.total_rate - deterministic.total_rate
    )
    absolute_difference = abs(signed_difference)
    signed_relative_difference = (
        signed_difference / deterministic.total_rate
    )
    absolute_relative_difference = abs(signed_relative_difference)

    if monte_carlo.standard_error > 0.0:
        standardized_difference = (
            signed_difference / monte_carlo.standard_error
        )
    elif signed_difference == 0.0:
        standardized_difference = 0.0
    else:
        standardized_difference = float(
            np.copysign(np.inf, signed_difference)
        )

    denergy = float(deterministic.energy[1] - deterministic.energy[0])
    dpitch = float(deterministic.pitch[1] - deterministic.pitch[0])
    energy_l1_difference = _normalized_l1_difference(
        reference=deterministic.energy_marginal,
        candidate=monte_carlo.energy_marginal,
        spacing=denergy,
        total_rate=deterministic.total_rate,
    )
    pitch_l1_difference = _normalized_l1_difference(
        reference=deterministic.pitch_marginal,
        candidate=monte_carlo.pitch_marginal,
        spacing=dpitch,
        total_rate=deterministic.total_rate,
    )

    relative_rate_passed = (
        absolute_relative_difference <= relative_tolerance
    )
    sigma_rate_passed = abs(standardized_difference) <= sigma_tolerance

    return ComparisonResult(
        case_index=pair.case_index,
        basename=pair.basename,
        deterministic_rate=deterministic.total_rate,
        monte_carlo_rate=monte_carlo.total_rate,
        monte_carlo_standard_error=monte_carlo.standard_error,
        signed_rate_difference=signed_difference,
        absolute_rate_difference=absolute_difference,
        signed_relative_rate_difference=signed_relative_difference,
        absolute_relative_rate_difference=absolute_relative_difference,
        standardized_rate_difference=standardized_difference,
        energy_marginal_l1_difference=energy_l1_difference,
        pitch_marginal_l1_difference=pitch_l1_difference,
        relative_rate_passed=relative_rate_passed,
        sigma_rate_passed=sigma_rate_passed,
        passed=relative_rate_passed and sigma_rate_passed,
    )
