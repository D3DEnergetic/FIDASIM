"""Write a readable Test 004 ion-sink comparison report."""

from datetime import datetime
from pathlib import Path


def _status(passed):
    """Return a compact status label."""
    return "PASS" if passed else "FAIL"


def _format_result(result):
    """Format every metric for one matched result pair."""
    return [
        f"Case {result.case_index:03d}: {result.basename}",
        f"  Status:                              {_status(result.passed)}",
        "  Total reaction rate [ions/(cm^3*s)]:",
        f"    Deterministic:                     "
        f"{result.deterministic_rate:.10e}",
        f"    Monte Carlo:                       "
        f"{result.monte_carlo_rate:.10e}",
        f"    Monte Carlo standard error:        "
        f"{result.monte_carlo_standard_error:.10e}",
        f"    Signed difference (MC-det):        "
        f"{result.signed_rate_difference:.10e}",
        f"    Absolute difference:               "
        f"{result.absolute_rate_difference:.10e}",
        f"    Signed relative difference:        "
        f"{result.signed_relative_rate_difference:.10e}",
        f"    Absolute relative difference:      "
        f"{result.absolute_relative_rate_difference:.10e} "
        f"[{_status(result.relative_rate_passed)}]",
        f"    Standardized difference:           "
        f"{result.standardized_rate_difference:.10e} "
        f"[{_status(result.sigma_rate_passed)}]",
        "  Diagnostic marginal differences:",
        f"    Energy normalized L1 difference:   "
        f"{result.energy_marginal_l1_difference:.10e}",
        f"    Pitch normalized L1 difference:    "
        f"{result.pitch_marginal_l1_difference:.10e}",
        "",
    ]


def write_report(results, comparison_config, output_directory):
    """Write per-case metrics, acceptance criteria, and overall status."""
    output_directory = Path(output_directory)
    timestamp = datetime.now().astimezone().strftime(
        "%Y-%m-%d %H:%M:%S %Z"
    )
    all_passed = all(result.passed for result in results)

    lines = [
        "=" * 46,
        f"OVERALL REGRESSION STATUS: {_status(all_passed)}",
        "=" * 46,
        "",
        "test_004 ion-sink comparison",
        "=" * 32,
        f"Computed: {timestamp}",
        "",
        f"Comment: {comparison_config['comment']}",
        f"Unified Test 004 config: {comparison_config['test_config']}",
        "",
        "Rate acceptance criteria",
        "------------------------",
        "A case passes only when both conditions below are satisfied:",
        "  abs(R_MC - R_det) / R_det <= "
        f"{comparison_config['rate_relative_tolerance']:.6g}",
        "  abs(R_MC - R_det) / SE_MC <= "
        f"{comparison_config['rate_sigma_tolerance']:.6g}",
        "",
        "The marginal L1 differences are diagnostic and do not determine "
        "pass/fail.",
        "",
    ]

    for result in results:
        lines.extend(_format_result(result))

    maximum_relative = max(
        result.absolute_relative_rate_difference for result in results
    )
    maximum_sigma = max(
        abs(result.standardized_rate_difference) for result in results
    )
    maximum_energy_l1 = max(
        result.energy_marginal_l1_difference for result in results
    )
    maximum_pitch_l1 = max(
        result.pitch_marginal_l1_difference for result in results
    )
    failed_cases = [
        f"{result.case_index:03d}" for result in results if not result.passed
    ]

    lines.extend(
        [
            "Overall summary",
            "---------------",
            f"Status:                                {_status(all_passed)}",
            f"Compared file pairs:                   {len(results)}",
            f"Failed case indices:                   "
            f"{', '.join(failed_cases) if failed_cases else 'none'}",
            f"Maximum absolute relative rate error:  "
            f"{maximum_relative:.10e}",
            f"Maximum absolute standardized error:   "
            f"{maximum_sigma:.10e}",
            f"Maximum energy marginal L1 difference: "
            f"{maximum_energy_l1:.10e}",
            f"Maximum pitch marginal L1 difference:  "
            f"{maximum_pitch_l1:.10e}",
            "",
        ]
    )

    report_filename = output_directory / "comparison_report.txt"
    report_filename.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote report: {report_filename}")
    return report_filename
