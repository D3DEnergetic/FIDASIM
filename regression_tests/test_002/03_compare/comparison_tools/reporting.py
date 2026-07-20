"""Format and write the text comparison report."""


def _format_result(result):
    """Format one reference-sampled comparison as a readable section."""
    reference = result.reference_moments
    sampled = result.sampled_moments

    return [
        f"File: {result.basename}",
        "  Density:",
        f"    Reference:       {reference.density:.2e}",
        f"    Sampled:         {sampled.density:.2e}",
        f"    Relative error:  {result.density_error:.2e}",
        "  Parallel temperature [keV]:",
        f"    Reference:       {reference.parallel_temperature:.2e}",
        f"    Sampled:         {sampled.parallel_temperature:.2e}",
        f"    Relative error:  {result.parallel_temperature_error:.2e}",
        "  Perpendicular temperature [keV]:",
        f"    Reference:       {reference.perpendicular_temperature:.2e}",
        f"    Sampled:         {sampled.perpendicular_temperature:.2e}",
        f"    Relative error:  {result.perpendicular_temperature_error:.2e}",
        "",
    ]


def write_report(results, output_directory):
    """Write per-file moments and an overall maximum-error summary."""
    report_lines = ["test_002 distribution comparison", "=" * 32, ""]

    # Add one readable section for every configured reference-sampled pair.
    for result in results:
        result_lines = _format_result(result)
        report_lines.extend(result_lines)

    # Finish with the largest error observed for each physical quantity.
    maximum_density_error = max(result.density_error for result in results)
    maximum_parallel_error = max(
        result.parallel_temperature_error for result in results
    )
    maximum_perpendicular_error = max(
        result.perpendicular_temperature_error for result in results
    )

    report_lines.extend(
        [
            "Overall summary",
            "---------------",
            f"Compared file pairs:                   {len(results)}",
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
