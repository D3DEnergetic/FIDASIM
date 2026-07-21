"""Write the Test 002 Stage 2 conversion-moment report."""

from datetime import datetime


def _relative_error(reference_value, converted_value):
    """Return the absolute relative difference between two scalar values."""
    if reference_value == 0.0:
        if converted_value == 0.0:
            return 0.0
        return float("inf")

    return abs(converted_value - reference_value) / abs(reference_value)


def write_moment_report(results, output_path):
    """Write source and converted moments for every selected location.

    Args:
        results (list[dict]): Reference path, output path, location, and moment
            values for each converted case.
        output_path (Path): Text report path.

    Returns:
        Path: Path of the report that was written.
    """
    computation_time = datetime.now().astimezone()
    formatted_time = computation_time.strftime("%Y-%m-%d %H:%M:%S %Z")

    report_lines = [
        "test_002 CQL3D-to-FIDASIM conversion moments",
        "=" * 47,
        f"Computed: {formatted_time}",
        "",
    ]

    for case_index, result in enumerate(results, start=1):
        reference = result["reference_moments"]
        converted = result["converted_moments"]

        density_error = _relative_error(reference[0], converted[0])
        parallel_error = _relative_error(reference[1], converted[1])
        perpendicular_error = _relative_error(reference[2], converted[2])

        report_lines.extend(
            [
                f"Case {case_index:03d}: {result['output_path'].name}",
                f"  Reference file: {result['reference_path']}",
                "  Selected location [cm]:",
                f"    R:                          {result['selected_r']:.6g}",
                f"    Z:                          {result['selected_z']:.6g}",
                "  Density [ions/cm^3]:",
                f"    Reference:                  {reference[0]:.2e}",
                f"    Converted:                  {converted[0]:.2e}",
                f"    Relative difference:        {density_error:.2e}",
                "  Parallel temperature [keV]:",
                f"    Reference:                  {reference[1]:.2e}",
                f"    Converted:                  {converted[1]:.2e}",
                f"    Relative difference:        {parallel_error:.2e}",
                "  Perpendicular temperature [keV]:",
                f"    Reference:                  {reference[2]:.2e}",
                f"    Converted:                  {converted[2]:.2e}",
                f"    Relative difference:        {perpendicular_error:.2e}",
                "",
            ]
        )

    output_path.write_text("\n".join(report_lines), encoding="utf-8")
    return output_path
