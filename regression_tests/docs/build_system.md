[Regression tests](../README.md) / [Learning guides](README.md) / Build system

**Navigation:** [Previous: HDF5 array ordering](hdf5_array_ordering.md) | [Up: Learning guides](README.md) | Next: —

---

# Regression-test build system

FIDASIM uses a hierarchy of Makefiles controlled by the repository-level
Makefile. Platform environment scripts keep compiler and library settings
explicit, reproducible, and portable instead of placing them in a user's shell
startup files.

The normal commands remain in the [regression-test portal](../README.md). This
guide explains what happens underneath them.

## Environment selection

Configure and source the appropriate platform environment file in each new
build shell:

```bash
cd "$FIDASIM_DIR"
# Inspect and configure the selected file before sourcing it.
source env/<platform>.sh
make
```

The script must be sourced rather than executed so its exported variables
remain available in the current shell when Make starts.

The environment file selects the compilers and one of two HDF5 build paths:

| `USE_SYSTEM_HDF5` | HDF5 path |
| --- | --- |
| `0` | Build and use the bundled HDF5 installation under `deps/hdf5`. |
| `1` | Use an existing HDF5 installation described by `HDF5_INCLUDE`, `HDF5_LIB`, and, when needed, `HDF5_EXTRA_LIBS`. |

The bundled HDF5 version is a legacy dependency and can fail to build with
newer compilers. Configure `USE_SYSTEM_HDF5=1` and valid system HDF5 paths on
affected platforms. See the [environment setup](../../env/Readme.md) and
[historical build notes](../../env/Historial_note.md) for platform-specific
configuration.

Sourcing the environment file exports the selected values before Make starts.
Assignments using `?=` supply defaults only when a value has not already been
defined, so they preserve the environment selection. The repository-level
Makefile then constructs and exports the compiler, module-path, and linker
settings consumed by the Makefiles below it.

## HDF5-dependent build graph

```text
configure and source env/<platform>.sh
                    |
       repository Makefile: fidasim
                    |
       +------------+------------+
       |                         |
USE_SYSTEM_HDF5=0        USE_SYSTEM_HDF5=1
       |                         |
add the deps target        omit the deps target
       |                         |
build and use              use HDF5_INCLUDE,
deps/hdf5                  HDF5_LIB, and
                           HDF5_EXTRA_LIBS
```

After applying that conditional selection, the repository-level Makefile
requests these targets:

```text
repository target: fidasim
├── deps                         [bundled-HDF5 builds only]
├── src
├── tables                       [depends on src]
├── python
└── regression_tests             [depends on src]
    ├── test_001
    ├── test_003
    ├── test_004
    ├── test_005
    └── test_006
```

The branches are sibling targets requested by the repository-level Makefile;
one sibling does not invoke the next. The explicit dependency labels mean that
`tables` and `regression_tests` wait for `src`. The `src` Makefile builds the
FIDASIM objects and executable. The repository-level `python` target
independently creates the configured Python executable link.

Test 002 consists of Python workflows and therefore has no compiled target in
the regression-test Makefile.

## Rebuilding the regression programs

After the complete build has established the selected HDF5 infrastructure and
FIDASIM objects, rebuild the regression programs from the repository root:

```bash
make regression_tests
```

This target invokes the `src` submake first, so changed FIDASIM sources are
rebuilt before the regression programs. It does not invoke `deps`, and
therefore assumes that the selected HDF5 installation is already available.

The target compiles Tests 001, 003, 004, 005, and 006. It does not execute any
regression test. Activate the Python environment required by the workflows,
then run the complete suite separately with:

```bash
conda activate FIDASIM_env
./regression_tests/run.sh
```

The execution runner includes the Python-only Test 002 and orders the tests so
that its generated distributions are available to Tests 003 and 004.

## Known limitations

- In a bundled-HDF5 build, `deps` and `src` are separate prerequisites of the
  top-level target; `src` does not explicitly depend on `deps`. A parallel
  clean build can therefore start compiling FIDASIM before the bundled HDF5
  module files and libraries have been created.
Until the bundled dependency ordering is corrected, use the full non-parallel
repository build when creating the bundled HDF5 installation from a clean
checkout.

---

**Navigation:** [Previous: HDF5 array ordering](hdf5_array_ordering.md) | [Up: Learning guides](README.md) | Next: —
