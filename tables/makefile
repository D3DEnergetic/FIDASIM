SHELL = /bin/bash

EX = $(TABLES_DIR)/generate_tables

ifeq ($(USE_MPI),y)
	EX := cafrun -np 1 $(EX)
endif

all: generate_tables table_settings

generate_tables: atomic_tables.o $(SRC_DIR)/hdf5_extra.o
	$(MPI_FC) $(CFLAGS) $^ -o $@ $(LFLAGS)

atomic_tables.o: atomic_tables.f90 $(SRC_DIR)/hdf5_extra.mod
	$(MPI_FC) $(CFLAGS) -c $< $(IFLAGS)

.PHONY: table_settings
table_settings: table_settings.dat generate_tables

table_settings.dat:
	$(EX) > table_settings.dat

.PHONY: atomic_tables
atomic_tables:
	$(TABLES_DIR)/generate_tables $(TABLES_DIR)/table_settings.dat $(NTHREADS)

clean:
	-rm -f *.mod *.o *.dat generate_tables
