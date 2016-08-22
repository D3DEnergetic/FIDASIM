SHELL = /bin/bash

#Default compiling options
USE_OPENMP = y
PROFILE = n
DEBUG = n

# Directories
FIDASIM_DIR := $(shell pwd)
SRC_DIR = $(FIDASIM_DIR)/src
DEPS_DIR = $(FIDASIM_DIR)/deps
TABLES_DIR = $(FIDASIM_DIR)/tables
LIB_DIR = $(FIDASIM_DIR)/lib
DOCS_DIR = $(FIDASIM_DIR)/docs

#Compilers
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

# Compiler Flags
# User defined Flags
VERSION := $(shell [ -e $(FIDASIM_DIR)/VERSION ] && cat $(FIDASIM_DIR)/VERSION)
ifneq ($(VERSION),)
	UFLAGS := -D_VERSION=\"$(VERSION)\" 
endif

BUILD := $(shell command -v git >/dev/null 2>&1 && \
	[ -d $(FIDASIM_DIR)/.git ] && \
	git describe --tags --always --dirty)

ifneq ($(BUILD),)
	UFLAGS := -D_VERSION=\"$(BUILD)\"
endif

# HDF5 variables
HDF5_LIB = $(DEPS_DIR)/hdf5/lib
HDF5_INCLUDE = $(DEPS_DIR)/hdf5/include
HDF5_FLAGS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz -ldl -Wl,-rpath,$(HDF5_LIB)

ifneq ($(findstring gfortran, $(FC)),)
	LFLAGS = -lm
	CFLAGS = -Ofast -g -fbacktrace -cpp
	DEBUG_CFLAGS = -O0 -g -cpp -fbacktrace -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -D_DEBUG
	OPENMP_FLAGS = -fopenmp -D_OMP
	PROF_FLAGS = -pg -D_PROF
endif

ifneq ($(findstring ifort, $(FC)),)
	LFLAGS = -limf -lm
	CFLAGS = -O2 -g -traceback -fpp
	DEBUG_CFLAGS = -O0 -g -fpp -traceback -debug all -check all -check bounds -fpe0 -warn -D_DEBUG
	OPENMP_FLAGS = -openmp -D_OMP
	PROF_FLAGS = -p -D_PROF
endif

ifeq ($(PROFILE),n)
ifeq ($(USE_OPENMP),y)
	CFLAGS := $(CFLAGS) $(UFLAGS) $(OPENMP_FLAGS)
else
	CFLAGS := $(CFLAGS) $(UFLAGS)
endif
else
	CFLAGS := $(CFLAGS) $(UFLAGS) $(PROF_FLAGS)
endif

ifeq ($(DEBUG),y)
	CFLAGS := $(DEBUG_CFLAGS) $(UFLAGS)
endif

LFLAGS := $(LFLAGS) $(HDF5_FLAGS) -L$(SRC_DIR)
IFLAGS := -I$(HDF5_INCLUDE) -I$(SRC_DIR)

# atomic table variables
NTHREADS = 1000 

# FORD documentation variables
FORD_FLAGS = -d $(SRC_DIR) -d $(TABLES_DIR) -d $(LIB_DIR) -p $(DOCS_DIR)/user-guide -o $(DOCS_DIR)/html
CHECK_LINKS = y

export FIDASIM_DIR
export SRC_DIR
export DEPS_DIR
export TABLES_DIR
export CFLAGS
export LFLAGS
export IFLAGS
export NTHREADS

fidasim: deps src tables

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
	@ if [[ $(CHECK_LINKS) == [yY]* ]]; then \
		echo "Checking for broken links..."; \
		linkchecker --ignore-url=/master $(DOCS_DIR)/html/index.html ; \
	  fi

clean_all: clean clean_deps clean_docs

clean: clean_src clean_tables
	-rm -f *.mod *.o fidasim

clean_src:
	@cd $(SRC_DIR); make clean

clean_deps:
	@cd $(DEPS_DIR); make clean

clean_tables:
	@cd $(TABLES_DIR); make clean

clean_docs:
	-rm -rf $(DOCS_DIR)/html
