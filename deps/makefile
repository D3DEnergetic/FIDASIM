SHELL = /bin/sh

DEPS_DIR = $(FIDASIM_DIR)/deps

ifeq ($(FIDASIM_COMPILER),gfortran)
	export CC=gcc
	export FC=gfortran
endif

ifeq ($(FIDASIM_COMPILER),ifort)
	export CC=icc
	export FC=ifort
	export CXX=icpc
endif

hdf5: 
	tar -zxvf $(DEPS_DIR)/hdf5-1.8.16.tar.gz
	cd $(DEPS_DIR)/hdf5-1.8.16; ./configure --prefix=$(DEPS_DIR)/hdf5 --enable-fortran --enable-cxx
	cd $(DEPS_DIR)/hdf5-1.8.16; make; make install
	cp $(DEPS_DIR)/hdf5-1.8.16/COPYING $(DEPS_DIR)/hdf5
	cp $(DEPS_DIR)/hdf5-1.8.16/README.txt $(DEPS_DIR)/hdf5
	-rm -rf hdf5-1.8.16

clean:
	-rm -rf hdf5 
