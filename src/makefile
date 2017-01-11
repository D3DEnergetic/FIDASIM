SHELL = /bin/bash

$(FIDASIM_DIR)/fidasim: fidasim.o eigensystem.o utilities.o hdf5_extra.o
	$(FC) $(CFLAGS) $^ -o $@ $(LFLAGS)

fidasim.o: fidasim.f90 eigensystem.mod utilities.mod hdf5_extra.mod
	$(FC) $(CFLAGS) -c $(IFLAGS)  $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(FC) $(CFLAGS) -c $(IFLAGS) $<

utilities.mod utilities.o: utilities.f90
	$(FC) $(CFLAGS) -c $<

hdf5_extra.mod hdf5_extra.o: hdf5_extra.f90
	$(FC) $(CFLAGS) -c $< $(IFLAGS)

clean:
	-rm -f *.mod *.o
