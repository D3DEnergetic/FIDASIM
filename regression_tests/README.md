# FIDASIM regression tests

This directory contains reproducible workflows for checking selected FIDASIM
calculations against trusted reference data. A regression test may generate or
store its reference data, run the implementation under test, and compare the
result in separate stages.

## Available tests

| Test | Purpose | Documentation |
| --- | --- | --- |
| `test_001` | Mean-free-path regression workflow. | [Test 001 portal](test_001/README.md) |
| `test_002` | Two-dimensional inverse-transform sampling of fast-ion energy-pitch distributions. | [Test 002 portal](test_002/README.md) |

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
