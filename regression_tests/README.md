# FIDASIM regression tests

This directory contains reproducible workflows for checking selected FIDASIM
calculations against trusted reference data. A regression test may generate or
store its reference data, run the implementation under test, and compare the
result in separate stages.

## Available tests

| Test | Purpose | Documentation |
| --- | --- | --- |
| `test_001` | Beam attenuation regression workflow. | [Test 001 portal](test_001/README.md) |
| `test_002` | Conversion of CQL3D velocity-space distributions into the FIDASIM energy-pitch representation. | [Test 002 portal](test_002/README.md) |
| `test_003` | Two-dimensional inverse-transform sampling of fast-ion energy-pitch distributions. | [Test 003 portal](test_003/README.md) |
| `test_004` | Charge-exchange ion-sink rate and energy-pitch distribution. | [Test 004 portal](test_004/README.md) |

## Build

Build the regression programs through FIDASIM's top-level makefile so the
configured compiler, module paths, HDF5 libraries, and FIDASIM objects are
available:

```bash
cd "$FIDASIM_DIR"
make regression_tests
```

Run instructions and dependencies specific to each test are documented in its
portal or stage README.
