SHELL = /bin/bash

DEPS_DIR = $(FIDASIM_DIR)/deps

all: hdf5 openblas

hdf5: 
	@printf "\nBuilding HDF5...\n"
	tar -zxvf $(DEPS_DIR)/hdf5-1.8.16.tar.gz > /dev/null 2>&1
	cd $(DEPS_DIR)/hdf5-1.8.16; ./configure --prefix=$(DEPS_DIR)/hdf5 --enable-fortran --enable-cxx > /dev/null 2>&1
	cd $(DEPS_DIR)/hdf5-1.8.16; make > /dev/null 2>&1; make install > /dev/null 2>&1
	cp $(DEPS_DIR)/hdf5-1.8.16/COPYING $(DEPS_DIR)/hdf5
	cp $(DEPS_DIR)/hdf5-1.8.16/README.txt $(DEPS_DIR)/hdf5
	-rm -rf hdf5-1.8.16
	@printf "HDF5 build complete\n\n"

openblas:
	@printf "Building OpenBLAS...\n"
	tar -zxvf $(DEPS_DIR)/OpenBLAS-0.2.18.tar.gz > /dev/null 2>&1
	cd $(DEPS_DIR)/OpenBLAS-0.2.18; make DYNAMIC_ARCH=1 USE_OPENMP=1
	cd $(DEPS_DIR)/OpenBLAS-0.2.18; make install PREFIX=$(DEPS_DIR)/openblas
	-rm -rf OpenBLAS-0.2.18
	@printf "OpenBLAS build complete\n\n"

clean:
	-rm -rf hdf5 hdf5-1.8.16 openblas OpenBLAS-0.2.18
