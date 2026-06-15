#!/usr/bin/env bash

module load PrgEnv-gnu
module load cray-hdf5

export FC=gfortran
export CC=gcc
export CXX=g++

export USE_SYSTEM_HDF5=1
export HDF5_INCLUDE="$HDF5_DIR/include"
export HDF5_LIB="$HDF5_DIR/lib"

ulimit -s unlimited
