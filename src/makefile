SHELL = /bin/bash

FIDADEPS=fidasim.o libfida.o eigensystem.o utilities.o hdf5_utils.o

ifeq ($(USE_MPI),y)
   FIDADEPS := mpi_utils.o $(FIDADEPS)
endif


$(FIDASIM_DIR)/fidasim$(EXEEXT): $(FIDADEPS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

mpi_utils.o: mpi_utils.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

fidasim.o: fidasim.f90 libfida.mod eigensystem.mod utilities.mod hdf5_utils.mod mpi_utils.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

libfida.mod libfida.o: libfida.f90 eigensystem.mod utilities.mod hdf5_utils.mod mpi_utils.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

utilities.mod utilities.o: utilities.f90
	$(MPI_FC) $(C_FLAGS) -c $<

hdf5_utils.mod hdf5_utils.o: hdf5_utils.f90
	$(MPI_FC) $(C_FLAGS) -c $< $(I_FLAGS)

clean:
	-rm -f *.mod *.o
