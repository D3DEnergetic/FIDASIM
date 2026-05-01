SHELL = /bin/bash

TEST_NAME := test_001
TEST_SRC  := 02_run_test
FIDA_SRC  := $(FIDASIM_DIR)/src

# Test sources/outputs live inside regression_tests/test_001/src/
TEST_F90  := $(TEST_SRC)/$(TEST_NAME).f90
TEST_OBJ  := $(TEST_SRC)/$(TEST_NAME).o
TEST_EXE  := $(TEST_SRC)/$(TEST_NAME)$(EXEEXT)

# Modules:
MOD_DIR := $(TEST_SRC)

# Objects already built by $(FIDASIM_DIR)/src/makefile
FIDAOBJS = $(FIDA_SRC)/libfida.o \
           $(FIDA_SRC)/eigensystem.o \
           $(FIDA_SRC)/utilities.o \
           $(FIDA_SRC)/hdf5_utils.o

ifeq ($(USE_MPI),y)
  FIDAOBJS := $(FIDA_SRC)/mpi_utils.o $(FIDAOBJS)
endif

.PHONY: all clean
all: $(TEST_EXE)

# Link regression executable with libfida and deps
$(TEST_EXE): $(TEST_OBJ) $(FIDAOBJS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

# Compile regression test; needs .mod files from $(FIDASIM_DIR)/src
$(TEST_OBJ): $(TEST_F90)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) -I$(MOD_DIR) -J$(MOD_DIR) $< -o $@

clean:
	-rm -f $(TEST_SRC)/*.o $(TEST_SRC)/*.mod $(TEST_EXE)
