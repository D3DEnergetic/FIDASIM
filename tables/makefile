SHELL = /bin/sh

SRC_DIR = $(FIDASIM_DIR)/src
TABLES_DIR = $(FIDASIM_DIR)/tables

HDF5_FLAGS = -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz

ifneq ($(findstring gfortran, $(FC)),)
	LFLAGS = $(HDF5_FLAGS) -lm
	CFLAGS = -Ofast -fopenmp -g -fbacktrace -cpp -D_OMP
endif

ifneq ($(findstring ifort, $(FC)),)
	LFLAGS = $(HDF5_FLAGS) -limf -lm
	CFLAGS = -O2 -fpp -g -traceback -openmp -D_OMP
endif

all: generate_tables table_settings

generate_tables: atomic_tables.o $(SRC_DIR)/hdf5_extra.o
	$(FC) $(CFLAGS) $^ -o $@ -L$(HDF5_LIB) -L$(SRC_DIR) $(LFLAGS)

atomic_tables.o: atomic_tables.f90 $(SRC_DIR)/hdf5_extra.mod
	$(FC) $(CFLAGS) -c $< -I$(HDF5_INCLUDE) -I$(SRC_DIR)

.PHONY: table_settings
table_settings: table_settings.dat generate_tables

table_settings.dat:
	$(TABLES_DIR)/generate_tables > table_settings.dat

.PHONY: atomic_tables
atomic_tables:
	$(TABLES_DIR)/generate_tables $(TABLES_DIR)/table_settings.dat $(NTHREADS)

clean:
	-rm -f *.mod *.o *.dat generate_tables
