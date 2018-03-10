SHELL = /bin/bash

DEPS_DIR = $(FIDASIM_DIR)/deps
HDF5_LOG = $(DEPS_DIR)/hdf5_build.log
OPENBLAS_LOG = $(DEPS_DIR)/openblas_build.log
OPENCOARRAYS_LOG = $(DEPS_DIR)/opencoarrays_build.log

# Needed for OpenMPI/MPICH
export OMPI_FC = $(FC)
export OMPI_CC = $(CC)
export OMPI_CXX = $(CXX)

TARGETS = hdf5 openblas
ifeq ($(USE_MPI),y)
    TARGETS := $(TARGETS) opencoarrays
endif

all: $(TARGETS)

hdf5:
	@printf "\nBuilding HDF5...\n"
	tar -zxvf $(DEPS_DIR)/hdf5-1.8.16.tar.gz >> /dev/null 2>&1
	cd $(DEPS_DIR)/hdf5-1.8.16; ./configure --prefix=$(DEPS_DIR)/hdf5 --enable-fortran --enable-cxx >> $(HDF5_LOG) 2>&1
	cd $(DEPS_DIR)/hdf5-1.8.16; make >> $(HDF5_LOG) 2>&1; make install >> $(HDF5_LOG) 2>&1
	cp $(DEPS_DIR)/hdf5-1.8.16/COPYING $(DEPS_DIR)/hdf5
	cp $(DEPS_DIR)/hdf5-1.8.16/README.txt $(DEPS_DIR)/hdf5
	-rm -rf hdf5-1.8.16
	@printf "HDF5 build complete\n\n"

openblas:
	@printf "Building OpenBLAS...\n"
	tar -zxvf $(DEPS_DIR)/OpenBLAS-0.2.20.tar.gz > /dev/null 2>&1
	cd $(DEPS_DIR)/OpenBLAS-0.2.20; make USE_THREAD=0 >> $(OPENBLAS_LOG) 2>&1
	cd $(DEPS_DIR)/OpenBLAS-0.2.20; make install PREFIX=$(DEPS_DIR)/openblas >> $(OPENBLAS_LOG) 2>&1
	-rm -rf OpenBLAS-0.2.20
	@printf "OpenBLAS build complete\n\n"

opencoarrays:
	@printf "\nBuilding OpenCoarrays...\n"
	tar -zvxf $(DEPS_DIR)/OpenCoarrays-1.9.3.tar.gz >> /dev/null 2>&1
	cd $(DEPS_DIR)/OpenCoarrays-1.9.3; ./install.sh --install-prefix=$(DEPS_DIR)/opencoarrays \
	    -y --with-fortran $(FC) --with-c $(CC) --with-cxx $(CXX) >> $(OPENCOARRAYS_LOG) 2>&1
	cp $(DEPS_DIR)/OpenCoarrays-1.9.3/LICENSE $(DEPS_DIR)/opencoarrays
	cp $(DEPS_DIR)/OpenCoarrays-1.9.3/GETTING_STARTED.md $(DEPS_DIR)/opencoarrays
	-rm -rf OpenCoarrays-1.9.3
	@printf "OpenCoarrays build is complete\n\n"

clean: clean_hdf5 clean_openblas clean_opencoarrays

clean_hdf5:
	-rm -rf hdf5 hdf5-1.8.16 $(HDF5_LOG)

clean_openblas:
	-rm -rf openblas OpenBLAS-0.2.20 $(OPENBLAS_LOG)

clean_opencoarrays:
	-rm -rf opencoarrays OpenCoarrays-1.9.3 $(OPENCOARRAYS_LOG)
