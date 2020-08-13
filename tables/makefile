SHELL = /bin/bash

EX = $(TABLES_DIR)/generate_tables$(EXEEXT)

TBDEPS=$(SRC_DIR)/hdf5_utils.o $(SRC_DIR)/utilities.o
TBDEPMODS=$(SRC_DIR)/hdf5_utils.mod $(SRC_DIR)/utilities.mod
ifeq ($(USE_MPI),y)
  TBDEPS    := $(SRC_DIR)/mpi_utils.o $(TBDEPS)
  TBDEPMODS := $(SRC_DIR)/mpi_utils.mod $(TBDEPMODS)
endif

all: generate_tables$(EXEEXT) table_settings

generate_tables$(EXEEXT): atomic_tables.o $(TBDEPS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

atomic_tables.o: atomic_tables.f90 $(TBDEPMODS)
	$(MPI_FC) $(C_FLAGS) -c $< $(I_FLAGS)

.PHONY: table_settings
table_settings: table_settings.dat generate_tables$(EXEEXT)

table_settings.dat:
	$(EX) > table_settings.dat

clean:
	-rm -f *.mod *.o *.dat generate_tables$(EXEEXT)
