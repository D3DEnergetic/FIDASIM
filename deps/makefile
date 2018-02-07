SHELL = /bin/bash

DEPS_DIR = $(FIDASIM_DIR)/deps
HDF5_LOG = $(DEPS_DIR)/hdf5_build.log
OPENCOARRAYS_LOG = $(DEPS_DIR)/opencoarrays_build.log

# Needed for OpenMPI/MPICH
export OMPI_FC = $(FC)
export OMPI_CC = $(CC)
export OMPI_CXX = $(CXX)

all: hdf5 opencoarrays

hdf5:
	@printf "\nBuilding HDF5...\n"
	tar -zxvf $(DEPS_DIR)/hdf5-1.8.16.tar.gz >> $(HDF5_LOG) 2>&1
	cd $(DEPS_DIR)/hdf5-1.8.16; ./configure --prefix=$(DEPS_DIR)/hdf5 --enable-fortran --enable-cxx >> $(HDF5_LOG) 2>&1
	cd $(DEPS_DIR)/hdf5-1.8.16; make >> $(HDF5_LOG) 2>&1; make install >> $(HDF5_LOG) 2>&1
	cp $(DEPS_DIR)/hdf5-1.8.16/COPYING $(DEPS_DIR)/hdf5
	cp $(DEPS_DIR)/hdf5-1.8.16/README.txt $(DEPS_DIR)/hdf5
	-rm -rf hdf5-1.8.16
	@printf "HDF5 build complete\n\n"

opencoarrays:
	@printf "\nBuilding OpenCoarrays...\n"
	tar -zvxf $(DEPS_DIR)/OpenCoarrays-1.9.3.tar.gz >> $(OPENCOARRAYS_LOG) 2>&1
	cd $(DEPS_DIR)/OpenCoarrays-1.9.3; ./install.sh --install-prefix=$(DEPS_DIR)/opencoarrays \
	    -y --with-fortran $(FC) --with-c $(CC) --with-cxx $(CXX) >> $(OPENCOARRAYS_LOG) 2>&1
	cp $(DEPS_DIR)/OpenCoarrays-1.9.3/LICENSE $(DEPS_DIR)/opencoarrays
	cp $(DEPS_DIR)/OpenCoarrays-1.9.3/GETTING_STARTED.md $(DEPS_DIR)/opencoarrays
	-rm -rf OpenCoarrays-1.9.3
	@printf "OpenCoarrays build is complete\n\n"

clean: clean_hdf5 clean_opencoarrays

clean_hdf5:
	-rm -rf hdf5 hdf5-1.8.16 $(HDF5_LOG)

clean_opencoarrays:
	-rm -rf opencoarrays OpenCoarrays-1.9.3 $(OPENCOARRAYS_LOG)
