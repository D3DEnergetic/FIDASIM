SHELL = /bin/sh

SRC_DIR = $(FIDASIM_DIR)/src
TABLES_DIR = $(FIDASIM_DIR)/tables

HDF5_FLAGS = -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz

define TABLES_DEFAULT = 
&general_settings
n_max = 12,    !Number of initial atomic energy levels
m_max = 12,    !Number of final atomic energy levels
tables_file = "$(FIDASIM_DIR)/tables/atomic_tables.h5"
/
!Hydrogen-Hydrogen Cross Sections
&H_H_cross
nenergy = 200, !Number of energy values
emin = 1.0E-3, !Minimum energy [keV]
emax = 8.0E2   !Maximum energy [keV]
/
!Hydrogen-Electron Cross Sections
&H_e_cross
nenergy = 200, !Number of energy values
emin = 1.0E-3, !Minimum energy [keV]
emax = 8.0E2   !Maximum energy [keV]
/
!Hydrogen-Impurity Cross Sections
&H_Aq_cross
q = 6,         !Impurity charge: Boron: 5, Carbon: 6, ...
nenergy = 200, !Number of energy values
emin = 1.0E-3, !Minimum energy [keV]
emax = 8.0E2   !Maximum energy [keV]
/
!Hydrogen-Hydrogen Reaction Rates
&H_H_rates
nenergy = 100, !Number of energy values
emin = 1.0E-3, !Minimum energy [keV]
emax = 4.0E2,  !Maximum energy [keV]
ntemp = 100,   !Number of temperature values
tmin = 1.0E-3, !Minimum ion temperature [keV]
tmax = 2.0E1   !Maximum ion temperature [keV]
/
!Hydrogen-Electron Reaction Rates
&H_e_rates
nenergy = 100, !Number of energy values
emin = 1.0E-3, !Minimum energy [keV]
emax = 4.0E2,  !Maximum energy [keV]
ntemp = 100,   !Number of temperature values
tmin = 1.0E-3, !Minimum electron temperature [keV]
tmax = 2.0E1   !Maximum electron temperature [keV]
/
!Hydrogen-Impurity Reaction Rates
&H_Aq_rates
q = 6,         !Impurity charge: Boron: 5, Carbon: 6, ...
mass = 12.011, !Impurity mass [amu]
nenergy = 100, !Number of energy values
emin = 1.0E-3, !Minimum energy [keV]
emax = 4.0E2,  !Maximum energy [keV]
ntemp = 100,   !Number of temperature values
tmin = 1.0E-3, !Minimum ion temperature [keV]
tmax = 2.0E1   !Maximum ion temperature [keV]
/
endef
export TABLES_DEFAULT

ifeq ($(FIDASIM_COMPILER),gfortran)
	LFLAGS = $(HDF5_FLAGS) -lm
	CFLAGS = -Ofast -fopenmp -g -fbacktrace -cpp -D_OMP
endif

ifeq ($(FIDASIM_COMPILER),ifort)
	LFLAGS = $(HDF5_FLAGS) -limf -lm
	CFLAGS = -O2 -fpp -g -traceback -openmp -D_OMP
endif

all: generate_tables table_settings

generate_tables: atomic_tables.o $(SRC_DIR)/hdf5_extra.o
	$(FIDASIM_COMPILER) $(CFLAGS) $^ -o $@ -L$(HDF5_LIB) -L$(SRC_DIR) $(LFLAGS)

atomic_tables.o: atomic_tables.f90 $(SRC_DIR)/hdf5_extra.mod
	$(FIDASIM_COMPILER) $(CFLAGS) -c $< -I$(HDF5_INCLUDE) -I$(SRC_DIR)

.PHONY: table_settings
table_settings: table_settings.dat

table_settings.dat:
	@echo "$$TABLES_DEFAULT" > table_settings.dat

.PHONY: atomic_tables
atomic_tables:
	$(TABLES_DIR)/generate_tables $(TABLES_DIR)/table_settings.dat $(NTHREADS)

clean:
	-rm -f *.mod *.o *.dat generate_tables
