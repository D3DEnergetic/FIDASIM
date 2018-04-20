SHELL = /bin/bash

EX = $(TABLES_DIR)/generate_tables

TBDEPS=$(SRC_DIR)/hdf5_extra.o
TBDEPMODS=$(SRC_DIR)/hdf5_extra.mod
ifeq ($(USE_MPI),y)
  TBDEPS    := $(SRC_DIR)/fidampi.o $(TBDEPS)
  TBDEPMODS := $(SRC_DIR)/fidampi.mod $(TBDEPMODS)
endif

all: generate_tables table_settings

generate_tables: atomic_tables.o $(TBDEPS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

atomic_tables.o: atomic_tables.f90 $(TBDEPMODS)
	$(MPI_FC) $(C_FLAGS) -c $< $(I_FLAGS)

.PHONY: table_settings
table_settings: table_settings.dat generate_tables

table_settings.dat:
	$(EX) > table_settings.dat

clean:
	-rm -f *.mod *.o *.dat generate_tables
