SHELL = /bin/bash

FIDADEPS=fidasim.o eigensystem.o utilities.o hdf5_extra.o

ifeq ($(USE_MPI),y)
   FIDADEPS := fidampi.o $(FIDADEPS)
endif


$(FIDASIM_DIR)/fidasim: $(FIDADEPS)
	$(MPI_FC) $(CFLAGS) $^ -o $@ $(LFLAGS)

fidampi.o: fidampi.f90
	$(MPI_FC) $(CFLAGS) -c $(IFLAGS) $< 

fidasim.o: fidasim.f90 eigensystem.mod utilities.mod hdf5_extra.mod fidampi.f90 fidampi.f90
	$(MPI_FC) $(CFLAGS) -c $(IFLAGS) $< 

eigensystem.mod eigensystem.o: eigensystem.f90
	$(MPI_FC) $(CFLAGS) -c $(IFLAGS) $<

utilities.mod utilities.o: utilities.f90
	$(MPI_FC) $(CFLAGS) -c $<

hdf5_extra.mod hdf5_extra.o: hdf5_extra.f90
	$(MPI_FC) $(CFLAGS) -c $< $(IFLAGS)

clean:
	-rm -f *.mod *.o
