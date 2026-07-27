"""Write the Test 002 physical-moment comparison report."""

from datetime import datetime

from .comparison import MOMENT_NAMES


MOMENT_LABELS = {
    "density": "Density [ions/cm^3]",
    "parallel_temperature": "Parallel temperature [keV]",
    "perpendicular_temperature": "Perpendicular temperature [keV]",
}


def _status(passed):
    """Return a compact pass/fail label."""
    return "PASS" if passed else "FAIL"


def write_report(results, relative_tolerance, output_directory):
    """Write comparison metrics, acceptance criteria, and overall status."""
    computation_time = datetime.now().astimezone()
    formatted_time = computation_time.strftime("%Y-%m-%d %H:%M:%S %Z")

    all_passed = all(result["passed"] for result in results)
    failed_cases = [
        f"{result['case_index']:03d}"
        for result in results
        if not result["passed"]
    ]
    maximum_errors = {}
    for moment_name in MOMENT_NAMES:
        maximum_errors[moment_name] = max(
            result["errors"][moment_name] for result in results
        )

    report_lines = [
        "=" * 46,
        f"OVERALL REGRESSION STATUS: {_status(all_passed)}",
        "=" * 46,
        "",
        "test_002 CQL3D-to-FIDASIM moment comparison",
        "=" * 46,
        f"Computed: {formatted_time}",
        "",
        "REGRESSION TEST SUMMARY",
        "=======================",
        f"Compared file pairs: {len(results)}",
        "Failed case indices: "
        f"{', '.join(failed_cases) if failed_cases else 'none'}",
        f"Moment relative tolerance: {relative_tolerance:.6g}",
        "Maximum density relative difference: "
        f"{maximum_errors['density']:.2e}",
        "Maximum T_parallel relative difference: "
        f"{maximum_errors['parallel_temperature']:.2e}",
        "Maximum T_perpendicular relative difference: "
        f"{maximum_errors['perpendicular_temperature']:.2e}",
        "",
        "Acceptance criterion",
        "--------------------",
        "Every moment in every case must satisfy:",
        "  absolute relative difference <= "
        f"{relative_tolerance:.6g}",
        "",
    ]

    for result in results:
        reference = result["reference"]
        converted = result["converted"]
        report_lines.extend(
            [
                f"Case {result['case_index']:03d}",
                f"  Status: {_status(result['passed'])}",
                f"  Reference file: {reference['path']}",
                f"  Converted file: {converted['path']}",
                "  Selected location [cm]:",
                f"    R:                          {reference['selected_r']:.6g}",
                f"    Z:                          {reference['selected_z']:.6g}",
            ]
        )

        for moment_name in MOMENT_NAMES:
            report_lines.extend(
                [
                    f"  {MOMENT_LABELS[moment_name]}:",
                    "    Reference:                  "
                    f"{reference['moments'][moment_name]:.2e}",
                    "    Converted:                  "
                    f"{converted['moments'][moment_name]:.2e}",
                    "    Relative difference:        "
                    f"{result['errors'][moment_name]:.2e} "
                    f"[{_status(result['moment_passed'][moment_name])}]",
                ]
            )
        report_lines.append("")

    report_lines.extend(["Overall summary", "---------------"])
    report_lines.append(f"Status: {_status(all_passed)}")
    report_lines.append(f"Compared file pairs: {len(results)}")
    report_lines.append(
        "Failed case indices: "
        f"{', '.join(failed_cases) if failed_cases else 'none'}"
    )
    for moment_name in MOMENT_NAMES:
        report_lines.append(
            f"Maximum {MOMENT_LABELS[moment_name]} relative difference: "
            f"{maximum_errors[moment_name]:.2e}"
        )
    report_lines.append("")

    output_path = output_directory / "moment_comparison.txt"
    output_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"Wrote report: {output_path}")
