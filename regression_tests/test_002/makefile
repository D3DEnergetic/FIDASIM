SHELL = /bin/bash

TEST_NAME := test_002
TEST_DIR  := 02_run_test
TEST_SRC  := $(TEST_DIR)/src
FIDA_SRC  := $(FIDASIM_DIR)/src

TEST_F90 := $(TEST_SRC)/$(TEST_NAME).f90
TEST_OBJ := $(TEST_SRC)/$(TEST_NAME).o
HDF5_F90 := $(TEST_SRC)/$(TEST_NAME)_hdf5.f90
HDF5_OBJ := $(TEST_SRC)/$(TEST_NAME)_hdf5.o
SAMPLING_F90 := $(TEST_SRC)/$(TEST_NAME)_sampling.f90
SAMPLING_OBJ := $(TEST_SRC)/$(TEST_NAME)_sampling.o
TEST_EXE := $(TEST_DIR)/$(TEST_NAME)$(EXEEXT)
MOD_DIR  := $(TEST_SRC)

FIDAOBJS := $(FIDA_SRC)/utilities.o

.PHONY: all clean
all: $(TEST_EXE)

$(TEST_EXE): $(HDF5_OBJ) $(SAMPLING_OBJ) $(TEST_OBJ) $(FIDAOBJS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

$(HDF5_OBJ): $(HDF5_F90)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) \
		-I$(MOD_DIR) -J$(MOD_DIR) $< -o $@

$(SAMPLING_OBJ): $(SAMPLING_F90) $(HDF5_OBJ) $(FIDAOBJS)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) \
		-I$(MOD_DIR) -J$(MOD_DIR) $< -o $@

$(TEST_OBJ): $(TEST_F90) $(HDF5_OBJ) $(SAMPLING_OBJ)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) \
		-I$(MOD_DIR) -J$(MOD_DIR) $< -o $@

clean:
	-rm -f $(TEST_SRC)/*.o $(TEST_SRC)/*.mod $(TEST_EXE)
