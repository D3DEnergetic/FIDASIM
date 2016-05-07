SHELL = /bin/sh

SUPPORTED_FC = gfortran ifort
SUPPORTED_CC = gcc icc
SUPPORTED_CXX = g++ icpc

HAS_FC := $(strip $(foreach SC, $(SUPPORTED_FC), $(findstring $(SC), $(FC))))
ifeq ($(HAS_FC),)
    $(error Fortran compiler $(FC) is not supported. Set FC to gfortran or ifort)
endif

HAS_CC := $(strip $(foreach SC, $(SUPPORTED_CC), $(findstring $(SC), $(CC))))
ifeq ($(HAS_CC),)
    $(error C compiler $(CC) is not supported. Set CC to gcc or icc)
endif

HAS_CXX := $(strip $(foreach SC, $(SUPPORTED_CXX), $(findstring $(SC), $(CXX))))
ifeq ($(HAS_CXX),)
    $(error C++ compiler $(CXX) is not supported. Set CXX to g++ or icpc)
endif

# directories
FIDASIM_DIR := $(shell pwd)
SRC_DIR = $(FIDASIM_DIR)/src
DEPS_DIR = $(FIDASIM_DIR)/deps
TABLES_DIR = $(FIDASIM_DIR)/tables
LIB_DIR = $(FIDASIM_DIR)/lib
DOCS_DIR = $(FIDASIM_DIR)/docs

# HDF5 variables
HDF5_LIB = $(DEPS_DIR)/hdf5/lib
HDF5_INCLUDE = $(DEPS_DIR)/hdf5/include
HDF5_FLAGS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz -ldl -Wl,-rpath,$(HDF5_LIB)

# atomic table variables
OUTPUT_DIR = $(TABLES_DIR)
NTHREADS = 1000 

# FORD documentation variables
FORD_FLAGS = -d $(SRC_DIR) -d $(TABLES_DIR) -d $(LIB_DIR) -p $(DOCS_DIR)/user-guide -o $(DOCS_DIR)/html

export FIDASIM_DIR
export SRC_DIR
export DEPS_DIR
export TABLES_DIR
export HDF5_LIB
export HDF5_INCLUDE
export HDF5_FLAGS
export OUTPUT_DIR
export NTHREADS

fidasim: deps src tables

debug: clean
debug: fidasim_debug

fidasim_debug: deps
	@cd $(SRC_DIR); make DEBUG=y

.PHONY: deps
deps:
	@cd $(DEPS_DIR); make

.PHONY: src
src:
	@cd $(SRC_DIR); make

.PHONY: tables
tables: src
	@cd $(TABLES_DIR); make

.PHONY: atomic_tables
atomic_tables:
	@cd $(TABLES_DIR); make atomic_tables

.PHONY: docs
docs:
	ford $(FORD_FLAGS) $(DOCS_DIR)/fidasim.md
	@echo "Checking for broken links..."
	linkchecker --check-extern $(DOCS_DIR)/html/index.html

clean_all: clean clean_deps clean_docs

clean: clean_src clean_tables
	-rm -f *.mod *.o fidasim fidasim_debug

clean_src:
	@cd $(SRC_DIR); make clean

clean_deps:
	@cd $(DEPS_DIR); make clean

clean_tables:
	@cd $(TABLES_DIR); make clean

clean_docs:
	-rm -f $(DOCS_DIR)/html
