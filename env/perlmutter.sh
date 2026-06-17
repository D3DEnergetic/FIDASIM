#!/usr/bin/env bash

#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# FIDASIM build environment for <platform>.
#
# This script configures the compiler environment and HDF5 selection for this
# platform. FIDASIM can be built in either of two HDF5 modes:
#
#   USE_SYSTEM_HDF5=0
#       Build and use the bundled HDF5 dependency distributed with FIDASIM.
#
#   USE_SYSTEM_HDF5=1
#       Use an existing HDF5 installation provided by the system, module
#       environment, package manager, or user.
#
# Platform scripts may choose the mode that is most appropriate for the target
# machine. Users can edit this file or override variables at the command line
# when testing a different configuration.
#
# Usage:
#
#     source env/<platform>.sh
#     make
# -----------------------------------------------------------------------------

# Select the GNU compiler environment on Perlmutter.
module load PrgEnv-gnu

# Load the Perlmutter-provided HDF5 module. This defines HDF5_DIR.
module load cray-hdf5

export FC=gfortran
export CC=gcc
export CXX=g++

# -----------------------------------------------------------------------------
# HDF5 configuration
#
# FIDASIM supports two HDF5 modes:
#
#   USE_SYSTEM_HDF5=0
#       Use the bundled HDF5 dependency distributed with FIDASIM.
#
#   USE_SYSTEM_HDF5=1
#       Use an existing HDF5 installation provided by the system.
#
# On Perlmutter, the recommended option is USE_SYSTEM_HDF5=1.
# Users may modify this setting if they wish to test the bundled
# HDF5 build path.
# -----------------------------------------------------------------------------

export USE_SYSTEM_HDF5=1

# -----------------------------------------------------------------------------
# System HDF5 location
#
# The following variables are only used when USE_SYSTEM_HDF5=1.
#
# The cray-hdf5 module defines HDF5_DIR, which points to the root
# directory of the HDF5 installation provided by Perlmutter.
#
# FIDASIM uses:
#
#   HDF5_INCLUDE : location of HDF5 header/module files
#   HDF5_LIB     : location of HDF5 libraries
#
# These variables are ignored when USE_SYSTEM_HDF5=0.
# -----------------------------------------------------------------------------

export HDF5_INCLUDE=$(HDF5_DIR)/include
export HDF5_LIB=$(HDF5_DIR)/lib

ulimit -s unlimited
