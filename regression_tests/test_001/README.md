[Regression tests](../README.md) / Test 001

**Navigation:** Previous: — | [Up: Regression tests](../README.md) | [Next: Test 002](../test_002/README.md)

---

# Test 001: mean-free-path regression

`test_001` is the mean-free-path regression workflow.

## FIDASIM code under test

| Module | Procedure | Description |
| --- | --- | --- |
| `libfida` | `colrad` | Evolves the level-resolved neutral states over each spatial step and returns the neutral densities used to calculate the mean free path. |

The regression comparison is therefore focused on the neutral attenuation and
mean free path produced by `colrad`. Atomic-table loading, construction of the
input `LocalProfiles` state, grid generation, HDF5 writing, and the procedures
called internally by `colrad` support the test or its implementation; they are
not independently tested by this regression workflow.

## Workflow

The files are arranged in three stages:

```text
test_001/
├── 01_reference/    # trusted reference calculation and data
├── 02_run_test/     # FIDASIM test program and generated result
└── 03_compare/      # comparison script, figures, and report
```

Build the regression programs from the FIDASIM repository root:

```bash
cd "$FIDASIM_DIR"
make regression_tests
```

The stage directories contain the legacy MATLAB and Fortran workflow used by
this test.

Stage 3 passes when the maximum relative error in the mean free path is no
greater than `0.12` (12%). Its terminal output and optional text report begin
with the overall `PASS` or `FAIL` status, measured maximum, and tolerance.
Diagnostic figures and report output are produced before MATLAB raises an
error for a failed comparison, so a batch run returns a nonzero exit status
without discarding the useful failure details.

---

**Navigation:** Previous: — | [Up: Regression tests](../README.md) | [Next: Test 002](../test_002/README.md)
