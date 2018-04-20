SHELL = /bin/bash

FIDADEPS=fidasim.o eigensystem.o utilities.o hdf5_extra.o

ifeq ($(USE_MPI),y)
   FIDADEPS := fidampi.o $(FIDADEPS)
endif


$(FIDASIM_DIR)/fidasim: $(FIDADEPS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

fidampi.o: fidampi.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

fidasim.o: fidasim.f90 eigensystem.mod utilities.mod hdf5_extra.mod fidampi.f90 fidampi.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(MPI_FC) $(C_FLAGS) -c $(I_FLAGS) $<

utilities.mod utilities.o: utilities.f90
	$(MPI_FC) $(C_FLAGS) -c $<

hdf5_extra.mod hdf5_extra.o: hdf5_extra.f90
	$(MPI_FC) $(C_FLAGS) -c $< $(I_FLAGS)

clean:
	-rm -f *.mod *.o
