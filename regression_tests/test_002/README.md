[Regression tests](../README.md) / Test 002

**Navigation:** [Previous: Test 001](../test_001/README.md) | [Up: Regression tests](../README.md) | [Next: Reference data](01_reference/README.md)

---

# Test 002: CQL3D-to-FIDASIM distribution conversion

`test_002` validates the conversion of a CQL3D distribution expressed on
normalized proper-velocity and pitch-angle coordinates, `f(u_bar, theta)`, into
the energy-pitch density, `F(E, pitch)`, required by FIDASIM.

The workflow will establish the relativistic coordinate transformation, the
phase-space Jacobian, the output units, and preservation of density, parallel
temperature, and perpendicular temperature.

## Workflow

| Stage | Purpose | Principal outputs |
| --- | --- | --- |
| [`01_reference`](01_reference/README.md) | Extract trusted single-location CQL3D distributions and calculate their velocity-space moments. | Compact HDF5 reference files, contour plots, and a moment report. |
| [`02_run_test`](02_run_test/README.md) | Convert each reference into the FIDASIM energy-pitch representation and calculate its moments. | FIDASIM HDF5 files, contour plots, and a moment report. |
| [`03_compare`](03_compare/README.md) | Compare the moments before and after conversion. | Summary report and relative-difference plot. |

Together, the three stages isolate the coordinate conversion and quantify the
moment error introduced by remapping onto the uniform FIDASIM grid.

---

**Navigation:** [Previous: Test 001](../test_001/README.md) | [Up: Regression tests](../README.md) | [Next: Reference data](01_reference/README.md)
