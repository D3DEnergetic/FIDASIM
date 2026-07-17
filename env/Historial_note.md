## Historical Note: FIDASIM Build on NERSC Perlmutter (June 2026)

In June 2026, FIDASIM was installed and tested on NERSC Perlmutter.

At the time, the standard FIDASIM build process relied on the bundled HDF5 dependency distributed with the source code:

```text
deps/hdf5-1.8.16.tar.gz
```

The Perlmutter software environment used during testing was:

| Component               | Version            |
| ----------------------- | ------------------ |
| Programming Environment | PrgEnv-gnu/8.6.0   |
| GNU Compiler Module     | gcc-native/14      |
| GNU Fortran             | 14.3.0             |
| MPI                     | cray-mpich/9.0.1   |
| System HDF5 Module      | cray-hdf5/1.12.2.9 |

When FIDASIM attempted to build the bundled HDF5 dependency, the build failed during the HDF5 compilation/install process.

The relevant error reported in `deps/hdf5_build.log` was:

```text
../src/H5Epublic.h:174:57: note: expected 'H5E_auto2_t'
but argument is of type 'herr_t (**)(hid_t, void *)'

make[3]: *** [Makefile:1667: testframe.lo] Error 1
make[2]: *** [Makefile:576: install-recursive] Error 1
```

This indicated that the bundled HDF5 version (1.8.16) could not be built successfully in the tested Perlmutter environment.

### Resolution

Rather than modifying or upgrading the bundled HDF5 build system, support was added for building FIDASIM against a system-provided HDF5 installation.

Perlmutter already provides a maintained HDF5 installation through the module system:

```bash
module load cray-hdf5
```

The build system was modified to support two modes:

1. Use the bundled HDF5 dependency (default behaviour).
2. Use an existing HDF5 installation provided by the host system.

The new build option introduced was:

```bash
export USE_SYSTEM_HDF5=1
```

When enabled, the FIDASIM build system uses the HDF5 installation specified by:

```bash
export HDF5_DIR=<path-to-hdf5>
```

instead of building the bundled HDF5 1.8.16 dependency.

### Outcome

Using the Perlmutter-provided HDF5 installation (`cray-hdf5/1.12.2.9`) successfully compiled and linked FIDASIM.

As of June 2026, the recommended build procedure on Perlmutter is:

```bash
source env/perlmutter.sh
make
```

where `env/perlmutter.sh` configures FIDASIM to use the system-provided HDF5 library.
