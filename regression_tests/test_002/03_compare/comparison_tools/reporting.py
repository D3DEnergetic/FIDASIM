"""Write the Test 002 physical-moment comparison report."""

from datetime import datetime

from .comparison import MOMENT_NAMES


MOMENT_LABELS = {
    "density": "Density [ions/cm^3]",
    "parallel_temperature": "Parallel temperature [keV]",
    "perpendicular_temperature": "Perpendicular temperature [keV]",
}


def write_report(results, output_directory):
    """Write every file-pair comparison and the maximum observed errors."""
    computation_time = datetime.now().astimezone()
    formatted_time = computation_time.strftime("%Y-%m-%d %H:%M:%S %Z")
    report_lines = [
        "test_002 CQL3D-to-FIDASIM moment comparison",
        "=" * 46,
        f"Computed: {formatted_time}",
        "",
    ]

    for result in results:
        reference = result["reference"]
        converted = result["converted"]
        report_lines.extend(
            [
                f"Case {result['case_index']:03d}",
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
                    f"{result['errors'][moment_name]:.2e}",
                ]
            )
        report_lines.append("")

    report_lines.extend(["Overall summary", "---------------"])
    report_lines.append(f"Compared file pairs: {len(results)}")
    for moment_name in MOMENT_NAMES:
        maximum_error = max(result["errors"][moment_name] for result in results)
        report_lines.append(
            f"Maximum {MOMENT_LABELS[moment_name]} relative difference: "
            f"{maximum_error:.2e}"
        )
    report_lines.append("")

    output_path = output_directory / "moment_comparison.txt"
    output_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"Wrote report: {output_path}")
    return output_path
