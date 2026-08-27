# FIDASIM regression tests

This directory contains reproducible workflows for checking selected FIDASIM
calculations against trusted reference data or embedded expected results. A
test may generate reference data, run a purpose-built executable linked to
FIDASIM modules, and compare its outputs in separate stages. Self-checking
tests may instead evaluate all expected results in one executable.

## Available tests

| Test | Purpose | Documentation |
| --- | --- | --- |
| `test_001` | Python-driven beam attenuation and mean-free-path regression workflow. | [Test 001 portal](test_001/README.md) |
| `test_002` | Conversion of CQL3D velocity-space distributions into the FIDASIM energy-pitch representation. | [Test 002 portal](test_002/README.md) |
| `test_003` | Two-dimensional inverse-transform sampling of fast-ion energy-pitch distributions. | [Test 003 portal](test_003/README.md) |
| `test_004` | Charge-exchange ion-sink rate and energy-pitch distribution. | [Test 004 portal](test_004/README.md) |
| `test_005` | Ray-triangle intersection and mesh collision dispatch. | [Test 005 portal](test_005/README.md) |

## Shared Python tools

[`regression_test_tools`](regression_test_tools/namelist_config_tools/README.md)
exposes the common namelist configuration utilities used across the
workflows. The reusable implementation is kept in its nested
`namelist_config_tools` package, together with a guide and runnable example.

## Architecture and learning guides

The regression workflows cross Python, compiled Fortran, and HDF5 boundaries.
Understanding those boundaries helps explain why the tests normalize inputs,
use purpose-built executables, preserve native outputs, and perform comparison
and visualization in a separate Python layer.

The [regression-test learning guides](docs/README.md) provide that shared
background:

- [Regression-test architecture](docs/regression_test_architecture.md)
  explains the input, compiled calculation, and output layers.
- [HDF5 array ordering](docs/hdf5_array_ordering.md) explains how Fortran and
  Python expose multidimensional datasets and how to preserve semantic axes.
- [Build system](docs/build_system.md) explains the environment scripts,
  Makefile hierarchy, HDF5 selection, and known limitations.

These pages document reusable concepts. Test-specific physics, data contracts,
commands, and pass criteria remain in each test portal.

## Build quick start

Configure and source the appropriate platform environment file in each new
build shell, then build from the FIDASIM repository root:

```bash
cd "$FIDASIM_DIR"
# Inspect and configure the selected file before sourcing it.
source env/<platform>.sh
make
```

The environment file must be sourced so its compiler and library settings
remain available to Make. See the [environment setup](../env/Readme.md) and
the detailed [build-system guide](docs/build_system.md) for HDF5 selection and
platform-specific considerations.

After a complete build has established the selected HDF5 installation and
FIDASIM objects, rebuild the regression targets with:

```bash
make regression_tests
```

The current target compiles Tests 001, 003, and 004 without running their
separate scientific comparison workflows. Test 005 is self-checking, so the
same target executes it, verifies its persistent report, and prints a summary
that distinguishes these verification levels. Test 002 consists of Python
workflows and has no compiled target in the regression-test Makefile.

Run instructions and dependencies specific to each test are documented in its
portal or stage README.
