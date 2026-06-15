#!/usr/bin/env bash

export FC=gfortran
export CC=gcc
export CXX=g++

export USE_SYSTEM_HDF5=1

# Edit these for your system
export HDF5_DIR=/usr
export HDF5_INCLUDE=/usr/include/hdf5/serial
export HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial
