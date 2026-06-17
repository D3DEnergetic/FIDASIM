# Environment Setup Scripts

This directory contains platform-specific environment setup scripts used when building FIDASIM.

Before compiling FIDASIM, source the script corresponding to your platform:

```bash
source env/<platform>.sh
make
```

For example:

```bash
source env/perlmutter.sh
make
```

or

```bash
source env/workstation.sh
make
```

## Purpose

FIDASIM requires an HDF5 library to read and write its input and output files.

Depending on the platform, HDF5 may be:

* provided by the operating system or HPC environment,
* installed locally by the user,
* or built from the bundled dependency distributed with FIDASIM.

The environment setup scripts provide a convenient mechanism for configuring the build environment and selecting the desired HDF5 installation.

In addition to HDF5 configuration, these scripts may define:

* compiler variables (`FC`, `CC`, `CXX`),
* module loads,
* include and library search paths,
* platform-specific build options.

## Available Scripts

### `perlmutter.sh`

Environment configuration for NERSC Perlmutter.

This script:

* loads the required Perlmutter modules,
* selects the GNU programming environment,
* configures the system-provided HDF5 installation,
* defines compiler variables,
* removes stack size limits commonly required by scientific Fortran codes.

### `workstation.sh`

Environment configuration for a local Linux workstation.

This script:

* configures the local compiler environment,
* defines the location of a locally installed HDF5 library,
* sets compiler variables and library paths as required by the local system.

## Common Variables

The following variables may be defined by the environment scripts:

| Variable          | Description                                                         |
| ----------------- | ------------------------------------------------------------------- |
| `FC`              | Fortran compiler                                                    |
| `CC`              | C compiler                                                          |
| `CXX`             | C++ compiler                                                        |
| `USE_SYSTEM_HDF5` | Selects whether FIDASIM uses the bundled HDF5 dependency (`0`) or an existing HDF5 installation (`1`) |
| `HDF5_DIR`        | Root directory of the system HDF5 installation                                                        |
| `HDF5_INCLUDE`    | Include directory used when `USE_SYSTEM_HDF5=1`                                                       |
| `HDF5_LIB`        | Library directory used when `USE_SYSTEM_HDF5=1`                                                       |


The exact values depend on the target platform.

## Notes

Platform-specific implementation details are documented within the individual environment scripts.

Historical information regarding platform-specific build issues, workarounds, and installation notes can be found in the documentation file `Historical_note.md`.
