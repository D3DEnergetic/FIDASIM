#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# FIDASIM build environment for a local Linux workstation.
#
# This script configures the compiler environment and HDF5 selection for this
# platform. FIDASIM can be built in either of two HDF5 modes:
#
#   USE_SYSTEM_HDF5=0
#       Build and use the bundled HDF5 dependency distributed with FIDASIM.
#
#   USE_SYSTEM_HDF5=1
#       Use an existing HDF5 installation provided by the operating system,
#       package manager, or a user-installed HDF5 build.
#
# Users may edit this file to select the desired HDF5 mode and configure
# the location of a system HDF5 installation.
#
# Usage:
#
#     source env/workstation.sh
#     make
# -----------------------------------------------------------------------------

# GNU compilers.
export FC=gfortran
export CC=gcc
export CXX=g++

# -----------------------------------------------------------------------------
# HDF5 configuration
#
# FIDASIM supports two HDF5 modes:
#
#   USE_SYSTEM_HDF5=0 (RECOMMENDED default)
#       Use the bundled HDF5 dependency distributed with FIDASIM.
#
#   USE_SYSTEM_HDF5=1
#       Use an existing HDF5 installation already available on the system.
#
# Users may change this variable as required for their installation.
# -----------------------------------------------------------------------------

export USE_SYSTEM_HDF5=1

# -----------------------------------------------------------------------------
# System HDF5 location
#
# The following variables are only used when USE_SYSTEM_HDF5=1.
#
# If HDF5 is not already installed, common installation methods are:
#
# Ubuntu/Debian:
#
#     sudo apt install libhdf5-dev
#
# Conda:
#
#     conda install hdf5
#
# After installation, locate the HDF5 Fortran wrapper:
#
#     which h5fc
#
# Then inspect the compiler and linker flags used by HDF5:
#
#     h5fc -show
#
# The output will contain flags similar to:
#
#     -I<include_path>
#     -L<library_path>
#
# These correspond to:
#
#     HDF5_INCLUDE=<include_path>
#     HDF5_LIB=<library_path>
#
# Example:
#
#     h5fc -show
#
# may report:
#
#     -I/usr/include/hdf5/serial
#     -L/usr/lib/x86_64-linux-gnu/hdf5/serial
#
# giving:
#
#     export HDF5_INCLUDE=/usr/include/hdf5/serial
#     export HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial
#
# The values below are examples and should be updated to match the
# HDF5 installation on your system.
# -----------------------------------------------------------------------------

export HDF5_INCLUDE=/usr/include/hdf5/serial
export HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial

# -----------------------------------------------------------------------------
# Additional HDF5 link libraries
#
# Some HDF5 installations require additional libraries beyond the core
# HDF5 libraries. These dependencies can be determined from:
#
#     h5fc -show
#
# Any additional "-l..." flags appearing after the HDF5 libraries may be
# added here.
#
# Example (Ubuntu 22.04):
#
#     h5fc -show
#
# reports:
#
#     -lcrypto -lcurl -lpthread -lsz -lz -ldl -lm
#
# Since FIDASIM already links against -lz and -ldl, we only need:
#
#     -lcrypto -lcurl -lpthread -lsz -lm
#
# Leave this variable empty if no additional libraries are required.
# -----------------------------------------------------------------------------

export HDF5_EXTRA_LIBS="-lcrypto -lcurl -lpthread -lsz -lm"

# Remove stack size limits
ulimit -s unlimited
