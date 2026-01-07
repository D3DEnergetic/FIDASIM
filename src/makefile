SHELL := /bin/bash

FIDADEPS=fidasim.o libfida.o eigensystem.o utilities.o hdf5_utils.o

ifeq ($(USE_MPI),y)
   FIDADEPS := mpi_utils.o $(FIDADEPS)
endif

# Only require mpi_utils.mod in MPI builds;
# In non-MPI builds the USE mpi_utils lines are removed by preprocessing, so no module or object is needed
MPIUTIL_MOD_DEP :=
ifeq ($(USE_MPI),y)
  MPIUTIL_MOD_DEP := mpi_utils.mod
endif

$(FIDASIM_DIR)/fidasim$(EXEEXT): $(FIDADEPS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

mpi_utils.mod mpi_utils.o: mpi_utils.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

fidasim.o: fidasim.f90 libfida.mod utilities.mod hdf5_utils.mod $(MPIUTIL_MOD_DEP)
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

libfida.mod libfida.o: libfida.f90 hdf5_utils.mod eigensystem.mod utilities.mod $(MPIUTIL_MOD_DEP)
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

utilities.mod utilities.o: utilities.f90
	$(MPI_FC) $(C_FLAGS) -c $<

hdf5_utils.mod hdf5_utils.o: hdf5_utils.f90
	$(MPI_FC) $(C_FLAGS) -c $< $(I_FLAGS)

.PHONY: clean
clean:
	-rm -f *.mod *.o
