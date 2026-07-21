"""Write the Test 002 Stage 1 reference-moment report."""

from datetime import datetime


def write_moment_report(results, output_path):
    """Write selected locations and physical moments to a readable text file.

    Args:
        results (list[dict]): File path, selected location, and moments per case.
        output_path (Path): Text report path.

    Returns:
        Path: Path of the report that was written.
    """
    computation_time = datetime.now().astimezone()
    formatted_time = computation_time.strftime("%Y-%m-%d %H:%M:%S %Z")

    report_lines = [
        "test_002 CQL3D reference moments",
        "=" * 34,
        f"Computed: {formatted_time}",
        "",
    ]

    for case_index, result in enumerate(results, start=1):
        moments = result["moments"]
        report_lines.extend(
            [
                f"Case {case_index:03d}: {result['path'].name}",
                "  Selected location [cm]:",
                f"    R:                          {result['selected_r']:.6g}",
                f"    Z:                          {result['selected_z']:.6g}",
                "  Density [ions/cm^3]:",
                f"    Reference:                  {moments.density:.2e}",
                "  Parallel temperature [keV]:",
                f"    Reference:                  {moments.parallel_temperature:.2e}",
                "  Perpendicular temperature [keV]:",
                f"    Reference:                  {moments.perpendicular_temperature:.2e}",
                "",
            ]
        )

    output_path.write_text("\n".join(report_lines), encoding="utf-8")
    return output_path
