SHELL = /bin/bash

FC = gfortran
CC = gcc

$(FIDASIM_DIR)/fidasim: fidasim.o eigensystem.o utilities.o hdf5_extra.o libfida.o
	$(FC) $(CFLAGS) $^ -o $@ $(LFLAGS)

fidasim.o: fidasim.f90 eigensystem.mod utilities.mod hdf5_extra.mod libfida.mod
	$(FC) $(CFLAGS) -c $(IFLAGS)  $<

libfida.mod libfida.o: libfida.f90 eigensystem.o utilities.o hdf5_extra.o
	$(FC) $(CFLAGS) -c $(IFLAGS) $<

fidanet.mod fidanet.o: fidanet.f90 libfida.o
	$(FC) $(CFLAGS) -c $(IFLAGS) $(LFLAGS) $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(FC) $(CFLAGS) -c $(IFLAGS) $<

utilities.mod utilities.o: utilities.f90
	$(FC) $(CFLAGS) -c $<

hdf5_extra.mod hdf5_extra.o: hdf5_extra.f90
	$(FC) $(CFLAGS) -c $< $(IFLAGS)

_fidanet.so:
	f90wrap -m fidanet -k kind_map fidanet.f90
	f2py-f90wrap --fcompiler=$(FC) --f90flags="$(CFLAGS)" $(IFLAGS) -L../deps/hdf5/lib/ -lhdf5 -lhdf5_fortran -lhdf5hl_fortran -lz -ldl -c -m _fidanet f90wrap*.f90 *.o

./fidatest: fidatest.o fidanet.o eigensystem.o utilities.o hdf5_extra.o libfida.o
	$(FC) $(CFLAGS) $^ -o $@ $(LFLAGS)

fidatest.o: fidatest.f90 eigensystem.mod utilities.mod hdf5_extra.mod libfida.mod fidanet.mod
	$(FC) $(CFLAGS) -c $(IFLAGS)  $<

.PHONY: fidanet
fidanet: fidanet.o _fidanet.so ./fidatest

clean:
	-rm -f *.mod *.o fidanet.py _fidanet.so