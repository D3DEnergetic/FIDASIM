SHELL = /bin/sh
DEBUG = n
EX = fidasim

SRC_DIR = $(FIDASIM_DIR)/src
HDF5_LIB = $(FIDASIM_DIR)/deps/hdf5/lib
HDF5_INCLUDE = $(FIDASIM_DIR)/deps/hdf5/include
HDF5_FLAGS = -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz -ldl
export LD_LIBRARY_PATH=$(HDF5_LIB):$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$(HDF5_INCLUDE):$LD_LIBRARY_PATH

VERSION := $(shell [ -e $(FIDASIM_DIR)/VERSION ] && cat $(FIDASIM_DIR)/VERSION)
ifneq ($(VERSION),)
	UFLAGS = -D_VERSION=\"$(VERSION)\" 
endif

BUILD := $(shell command -v git >/dev/null 2>&1 && \
	[ -d $(FIDASIM_DIR)/.git ] && \
	git --git-dir=$(FIDASIM_DIR)/.git describe --tags --always)

ifneq ($(BUILD),)
	UFLAGS = -D_VERSION=\"$(BUILD)\"
endif

ifeq ($(FIDASIM_COMPILER),gfortran)
	LFLAGS = $(HDF5_FLAGS) -lm
	CFLAGS = -Ofast -fopenmp -g -fbacktrace -cpp -D_OMP $(UFLAGS)
	DEBUG_CFLAGS = -O0 -g -cpp -fbacktrace -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow
endif

ifeq ($(FIDASIM_COMPILER),ifort)
	LFLAGS = $(HDF5_FLAGS) -limf -lm
	CFLAGS = -O2 -fpp -g -traceback -openmp -D_OMP $(UFLAGS)
	DEBUG_CFLAGS = -O0 -fpp -g -traceback -debug all -check all -check bounds -fpe:0 -warn
endif

ifeq ($(DEBUG),y)
	CFLAGS = $(DEBUG_CFLAGS) $(UFLAGS)
	EX = fidasim_debug
endif

$(FIDASIM_DIR)/$(EX): fidasim.o eigensystem.o parallel_rng.o hdf5_extra.o
	$(FIDASIM_COMPILER) $(CFLAGS) $^ -o $@ -L$(HDF5_LIB) $(LFLAGS) -L$(SRC_DIR)

fidasim.o: fidasim.f90 eigensystem.mod parallel_rng.mod hdf5_extra.mod
	$(FIDASIM_COMPILER) $(CFLAGS) -c -I$(HDF5_INCLUDE) -I$(SRC_DIR) $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(FIDASIM_COMPILER) $(CFLAGS) -c -I$(SRC_DIR) $<

parallel_rng.mod parallel_rng.o: parallel_rng.f90
	$(FIDASIM_COMPILER) $(CFLAGS) -c $<

hdf5_extra.mod hdf5_extra.o: hdf5_extra.f90
	$(FIDASIM_COMPILER) $(CFLAGS) -c $< -I$(HDF5_INCLUDE) -I$(SRC_DIR)

clean:
	-rm -f *.mod *.o
