SHELL = /bin/bash

#Default compiling options
USE_OPENMP = y
USE_MPI = n
ARCH = n
PROFILE = n
DEBUG = n

# Directories
FIDASIM_DIR := $(shell pwd)
SRC_DIR = $(FIDASIM_DIR)/src
DEPS_DIR = $(FIDASIM_DIR)/deps
TABLES_DIR = $(FIDASIM_DIR)/tables
LIB_DIR = $(FIDASIM_DIR)/lib
DOCS_DIR = $(FIDASIM_DIR)/docs

PYTHON_EXEC = $(shell which python3 python 2> /dev/null | head -1)
ifeq ($(PYTHON_EXEC),)
    $(error Python3 executable was not found)
endif

#Operating Systems
OS := $(shell uname)

#Compilers
SUPPORTED_FC = gfortran pgf90 ifort
SUPPORTED_CC = gcc pgcc
SUPPORTED_CXX = g++ pgc++

HAS_FC := $(strip $(foreach SC, $(SUPPORTED_FC), $(findstring $(SC), $(FC))))
ifeq ($(HAS_FC),)
    $(error Fortran compiler $(FC) is not supported. Set FC to gfortran)
endif

HAS_CC := $(strip $(foreach SC, $(SUPPORTED_CC), $(findstring $(SC), $(CC))))
ifeq ($(HAS_CC),)
    $(error C compiler $(CC) is not supported. Set CC to gcc)
endif

HAS_CXX := $(strip $(foreach SC, $(SUPPORTED_CXX), $(findstring $(SC), $(CXX))))
ifeq ($(HAS_CXX),)
    $(error C++ compiler $(CXX) is not supported. Set CXX to g++)
endif

# Compiler Flags
# User defined Flags
VERSION := $(shell [ -e $(FIDASIM_DIR)/VERSION ] && cat $(FIDASIM_DIR)/VERSION)
ifneq ($(VERSION),)
	U_FLAGS := -D_VERSION=\"$(VERSION)\"
endif

BUILD := $(shell command -v git >/dev/null 2>&1 && \
	[ -d $(FIDASIM_DIR)/.git ] && \
	git describe --tags --always --dirty)

ifneq ($(BUILD),)
	U_FLAGS := -D_VERSION=\"$(BUILD)\"
endif

# HDF5 variables
HDF5_LIB = $(DEPS_DIR)/hdf5/lib
HDF5_INCLUDE = $(DEPS_DIR)/hdf5/include
ifeq ($(OS),Linux)
	HDF5_FLAGS = -L$(HDF5_LIB) -Wl,-Bstatic -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -Wl,-Bdynamic -lz -ldl
endif
ifneq ($(findstring MINGW, $(OS)),)
	HDF5_FLAGS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lm -lz -lws2_32
	U_FLAGS=-static
	EXEEXT=.exe
	CFLAGS=-DH5_HAVE_WIN32_API
	export CFLAGS
	export EXEEXT
endif
ifeq ($(OS),Darwin)
	HDF5_FLAGS = -L/usr/lib -L$(HDF5_LIB) -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz -ldl -Wl,-rpath,$(HDF5_LIB)
endif

ifneq ($(findstring gfortran, $(FC)),)
        L_FLAGS = -lm
        COMMON_CFLAGS = -Ofast -g -fbacktrace -cpp -Wfatal-errors
        DEBUG_CFLAGS = -O0 -g -cpp -fbacktrace -fcheck=all -Wall -ffpe-trap=invalid,zero,overflow -D_DEBUG
        OPENMP_FLAGS = -fopenmp -D_OMP
        MPI_FLAGS = -D_MPI
        PROF_FLAGS = -pg -D_PROF
ifneq ($(ARCH),n)
        COMMON_CFLAGS := $(COMMON_CFLAGS) -march=$(ARCH)
endif
endif
ifneq ($(findstring pgf90, $(FC)),)
        L_FLAGS = -lm -llapack -lblas
        COMMON_CFLAGS = -O3 -Mpreprocess -D_DEF_INTR -D_USE_BLAS
        DEBUG_CFLAGS = -O0 -g -Mpreprocess -traceback -D_DEF_INTR -D_USE_BLAS -D_DEBUG
        OPENMP_FLAGS = -mp -D_OMP
        MPI_FLAGS = -D_MPI
        PROF_FLAGS = -pg -D_PROF
ifneq ($(ARCH),n)
        COMMON_CFLAGS := $(COMMON_CFLAGS) -tp=$(ARCH)
endif
endif
ifneq ($(findstring ifort, $(FC)),)
        L_FLAGS = -lm -mkl
        COMMON_CFLAGS = -Ofast -g -traceback -fpp -D_USE_BLAS
        DEBUG_CFLAGS = -O0 -g -fpp -traceback -check all -D_USE_BLAS -D_DEBUG
        OPENMP_FLAGS = -qopenmp -D_OMP
        MPI_FLAGS = -D_MPI
        PROF_FLAGS = -p -D_PROF
ifneq ($(ARCH),n)
        COMMON_CFLAGS := $(COMMON_CFLAGS) -x$(ARCH)
endif
endif


MPI_FC = $(FC)
C_FLAGS = $(COMMON_CFLAGS) $(U_FLAGS)

ifeq ($(DEBUG),y)
	USE_OPENMP = n
	USE_MPI = n
	PROFILE = n
	MPI_FC = $(FC)
	C_FLAGS = $(DEBUG_CFLAGS) $(UFLAGS)
endif

ifeq ($(PROFILE),y)
	USE_OPENMP = n
	USE_MPI = n
	MPI_FC = $(FC)
	C_FLAGS = $(COMMON_CFLAGS) $(U_FLAGS) $(PROF_FLAGS)
endif

ifeq ($(USE_MPI),y)
	USE_OPENMP = n
	MPI_FC = mpif90
	C_FLAGS = $(COMMON_CFLAGS) $(U_FLAGS) $(MPI_FLAGS)
	# >>>>>> JFCM 2024_04_16:
	# To allow debugging with MPI enabled:
	ifeq ($(DEBUG),y)
	 C_FLAGS = $(DEBUG_CFLAGS) $(U_FLAGS) $(MPI_FLAGS)
	endif
	# <<<<<< JFCM 2024_04_16:
endif

ifeq ($(USE_OPENMP),y)
	MPI_FC = $(FC)
	C_FLAGS = $(COMMON_CFLAGS) $(U_FLAGS) $(OPENMP_FLAGS)
endif

L_FLAGS := $(L_FLAGS) $(HDF5_FLAGS) -L$(SRC_DIR)
I_FLAGS := -I$(HDF5_INCLUDE) -I$(SRC_DIR)

# atomic table variables
NTHREADS = 1000

# FORD documentation variables
FORD_FLAGS = -d $(SRC_DIR) -d $(TABLES_DIR) -d $(LIB_DIR)/idl -d $(LIB_DIR)/python/fidasim -p $(DOCS_DIR)/user-guide -o $(DOCS_DIR)/html
CHECK_LINKS = y

export FIDASIM_DIR
export SRC_DIR
export DEPS_DIR
export TABLES_DIR
export MPI_FC
export USE_MPI
export C_FLAGS
export L_FLAGS
export I_FLAGS
export NTHREADS

fidasim: deps src tables python

.PHONY: deps
deps:
	@cd $(DEPS_DIR); make

.PHONY: src
src:
	@cd $(SRC_DIR); make

.PHONY: tables
tables: src
	@cd $(TABLES_DIR); make

.PHONY: docs
docs:
	ford $(FORD_FLAGS) $(DOCS_DIR)/fidasim.md
	@ if [[ $(CHECK_LINKS) == [yY]* ]]; then \
		echo "Checking for broken links..."; \
		linkchecker --ignore-url=/master $(DOCS_DIR)/html/index.html ; \
	  fi

.PHONY: python
python:
	@echo "Linking python executable: $(PYTHON_EXEC) --> $(DEPS_DIR)/python"
	@ln -sf $(PYTHON_EXEC) $(DEPS_DIR)/python

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

help:
	@echo ""
	@echo "Makefile Targets:"
	@echo "    fidasim"
	@echo "        Builds fidasim executable (Default Target)"
	@echo "    docs"
	@echo "        Builds FIDASIM documentation website in $(DOCS_DIR)/html"
	@echo "    clean"
	@echo "        Deletes object and module files in the src and tables directories and deletes the fidasim executable"
	@echo "    clean_all"
	@echo "        Cleans src, tables, deps, and docs directories"
	@echo "    clean_src"
	@echo "        Cleans src directory"
	@echo "    clean_tables"
	@echo "        Cleans the tables directory"
	@echo "    clean_deps"
	@echo "        Cleans deps directory including HDF5 build"
	@echo "    clean_docs"
	@echo "        Deletes FIDASIM documentation website directory"
	@echo ""
	@echo "Compiling Options:"
	@echo "    USE_MPI = (y|n)"
	@echo "        Turns on MPI parallelization and turns off OpenMP"
	@echo "        Default: n"
	@echo "    USE_OPENMP = (y|n)"
	@echo "        Turns on OpenMP"
	@echo "        Default: y"
	@echo "    ARCH = (<TARGET>|n)"
	@echo "        Compile code optimized for <TARGET> e.g. ARCH=haswell/CORE-AVX2/native"
	@echo "        Options differ depending on compiler"
	@echo "        Generated code may not run on different architectures"
	@echo "        Default: n"
	@echo "    DEBUG = (y|n)"
	@echo "        Turns off all optimizations and OpenMP"
	@echo "        Turns on bounds checking and floating point exceptions"
	@echo "        Default: n"
	@echo "    PROFILE = (y|n)"
	@echo "        Turns on gprof profiling and turns off OpenMP"
	@echo "        Default: n"
	@echo "    CHECK_LINKS = (y|n)"
	@echo "        Checks whether documentation website hyperlinks are active using linkchecker"
	@echo "        Default: y"
