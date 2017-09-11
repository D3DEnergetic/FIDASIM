SHELL = /bin/bash

$(FIDASIM_DIR)/fidasim: fidasim.o eigensystem.o utilities.o hdf5_extra.o _libfida*.so
	$(FC) $(CFLAGS) $^ -o $@ $(LFLAGS)

fidasim.o: fidasim.f90 eigensystem.mod utilities.mod hdf5_extra.mod
	$(FC) $(CFLAGS) -c $(IFLAGS)  $<

eigensystem.mod eigensystem.o: eigensystem.f90
	$(FC) $(CFLAGS) -c $(IFLAGS) $<

utilities.mod utilities.o: utilities.f90
	$(FC) $(CFLAGS) -c $<

hdf5_extra.mod hdf5_extra.o: hdf5_extra.f90
	$(FC) $(CFLAGS) -c $< $(IFLAGS)

_libfida*.so:
	f90wrap -m libfida_core -k kind_map fidasim.f90 utilities.f90
	f2py-f90wrap -m _libfida_core --fcompiler=$(FC) --f90flags="$(CFLAGS)" $(IFLAGS) -lgomp -L../deps/hdf5/lib/ -lhdf5 -lhdf5_fortran -lhdf5hl_fortran -lz -ldl -c f90wrap*.f90 *.o

clean:
	-rm -f *.mod *.o f90wrap*.f90 libfida.py *.so
