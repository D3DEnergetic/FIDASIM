"""Format and write the text comparison report."""

from datetime import datetime


def _status(passed):
    """Return a compact pass/fail label."""
    return "PASS" if passed else "FAIL"


def _format_result(result):
    """Format one reference-sampled comparison as a readable section."""
    reference = result.reference_moments
    sampled = result.sampled_moments

    return [
        f"Case {result.case_index:03d}: {result.basename}",
        f"  Status:          {_status(result.passed)}",
        "  Density:",
        f"    Reference:       {reference.density:.2e}",
        f"    Sampled:         {sampled.density:.2e}",
        f"    Relative error:  {result.density_error:.2e} "
        f"[{_status(result.density_passed)}]",
        "  Parallel temperature [keV]:",
        f"    Reference:       {reference.parallel_temperature:.2e}",
        f"    Sampled:         {sampled.parallel_temperature:.2e}",
        f"    Relative error:  {result.parallel_temperature_error:.2e} "
        f"[{_status(result.parallel_temperature_passed)}]",
        "  Perpendicular temperature [keV]:",
        f"    Reference:       {reference.perpendicular_temperature:.2e}",
        f"    Sampled:         {sampled.perpendicular_temperature:.2e}",
        f"    Relative error:  {result.perpendicular_temperature_error:.2e} "
        f"[{_status(result.perpendicular_temperature_passed)}]",
        "",
    ]


def write_report(results, relative_tolerance, output_directory):
    """Write results, acceptance criteria, and the overall pass status."""
    computation_time = datetime.now().astimezone()
    formatted_time = computation_time.strftime("%Y-%m-%d %H:%M:%S %Z")

    all_passed = all(result.passed for result in results)
    failed_cases = [
        f"{result.case_index:03d}" for result in results if not result.passed
    ]
    maximum_density_error = max(result.density_error for result in results)
    maximum_parallel_error = max(
        result.parallel_temperature_error for result in results
    )
    maximum_perpendicular_error = max(
        result.perpendicular_temperature_error for result in results
    )

    report_lines = [
        "=" * 46,
        f"OVERALL REGRESSION STATUS: {_status(all_passed)}",
        "=" * 46,
        "",
        "test_003 distribution comparison",
        "=" * 32,
        f"Computed: {formatted_time}",
        "",
        "REGRESSION TEST SUMMARY",
        "=======================",
        f"Compared file pairs: {len(results)}",
        "Failed case indices: "
        f"{', '.join(failed_cases) if failed_cases else 'none'}",
        f"Moment relative tolerance: {relative_tolerance:.6g}",
        f"Maximum density relative error: {maximum_density_error:.2e}",
        "Maximum T_parallel relative error: "
        f"{maximum_parallel_error:.2e}",
        "Maximum T_perpendicular relative error: "
        f"{maximum_perpendicular_error:.2e}",
        "",
        "Acceptance criterion",
        "--------------------",
        "Every moment in every case must satisfy:",
        f"  absolute relative error <= {relative_tolerance:.6g}",
        "",
    ]

    # Add one readable section for every configured reference-sampled pair.
    for result in results:
        result_lines = _format_result(result)
        report_lines.extend(result_lines)

    report_lines.extend(
        [
            "Overall summary",
            "---------------",
            f"Status:                                {_status(all_passed)}",
            f"Compared file pairs:                   {len(results)}",
            "Failed case indices:                   "
            f"{', '.join(failed_cases) if failed_cases else 'none'}",
            f"Maximum density relative error:        {maximum_density_error:.2e}",
            "Maximum T_parallel relative error:     "
            f"{maximum_parallel_error:.2e}",
            "Maximum T_perpendicular relative error: "
            f"{maximum_perpendicular_error:.2e}",
            "",
        ]
    )

    report_filename = output_directory / "comparison_report.txt"
    report_filename.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"Wrote report: {report_filename}")
