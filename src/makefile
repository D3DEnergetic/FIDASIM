SHELL = /bin/bash
DEBUG = n
EX = fidasim

VERSION := $(shell [ -e $(FIDASIM_DIR)/VERSION ] && cat $(FIDASIM_DIR)/VERSION)
ifneq ($(VERSION),)
	UFLAGS = -D_VERSION=\"$(VERSION)\" 
endif

BUILD := $(shell command -v git >/dev/null 2>&1 && \
	[ -d $(FIDASIM_DIR)/.git ] && \
	git describe --tags --always --dirty)

ifneq ($(BUILD),)
	UFLAGS = -D_VERSION=\"$(BUILD)\"
endif

ifneq ($(findstring gfortran, $(FC)),)
	LFLAGS = $(HDF5_FLAGS) -L$(SRC_DIR) -lm 
	CFLAGS = -Ofast -fopenmp -g -fbacktrace -cpp -D_OMP $(UFLAGS) #-pg #gprof profiling
	DEBUG_CFLAGS = -O0 -g -cpp -fbacktrace -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow
endif

ifneq ($(findstring ifort, $(FC)),)
	LFLAGS = $(HDF5_FLAGS) -L$(SRC_DIR) -limf -lm
	CFLAGS = -O2 -fpp -g -traceback -openmp -D_OMP $(UFLAGS) #-p #gprof profiling
	DEBUG_CFLAGS = -O0 -fpp -g -traceback -debug all -check all -check bounds -fpe:0 -warn
endif

ifeq ($(DEBUG),y)
	CFLAGS = $(DEBUG_CFLAGS) $(UFLAGS)
	EX = fidasim_debug
endif

$(FIDASIM_DIR)/$(EX): fidasim.o eigensystem.o parallel_rng.o hdf5_extra.o
	$(FC) $(CFLAGS) $^ -o $@ $(LFLAGS)

fidasim.o: fidasim.f90 eigensystem.mod parallel_rng.mod hdf5_extra.mod
	$(FC) $(CFLAGS) -c -I$(HDF5_INCLUDE) -I$(SRC_DIR) $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(FC) $(CFLAGS) -c -I$(SRC_DIR) $<

parallel_rng.mod parallel_rng.o: parallel_rng.f90
	$(FC) $(CFLAGS) -c $<

hdf5_extra.mod hdf5_extra.o: hdf5_extra.f90
	$(FC) $(CFLAGS) -c $< -I$(HDF5_INCLUDE) -I$(SRC_DIR)

clean:
	-rm -f *.mod *.o
