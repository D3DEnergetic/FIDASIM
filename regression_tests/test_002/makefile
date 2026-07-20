SHELL = /bin/bash

TEST_NAME  := test_002
TEST_DIR   := 02_run_test
SRC_DIR    := $(TEST_DIR)/src
MODULE_DIR := $(SRC_DIR)/modules
BUILD_DIR  := $(TEST_DIR)/build
FIDA_SRC   := $(FIDASIM_DIR)/src

TEST_F90     := $(SRC_DIR)/$(TEST_NAME).f90
CONFIG_F90   := $(MODULE_DIR)/$(TEST_NAME)_config.f90
HDF5_F90     := $(MODULE_DIR)/$(TEST_NAME)_hdf5.f90
SAMPLING_F90 := $(MODULE_DIR)/$(TEST_NAME)_sampling.f90

TEST_OBJ     := $(BUILD_DIR)/$(TEST_NAME).o
CONFIG_OBJ   := $(BUILD_DIR)/$(TEST_NAME)_config.o
HDF5_OBJ     := $(BUILD_DIR)/$(TEST_NAME)_hdf5.o
SAMPLING_OBJ := $(BUILD_DIR)/$(TEST_NAME)_sampling.o
TEST_EXE     := $(TEST_DIR)/$(TEST_NAME)$(EXEEXT)

FIDAOBJS := $(FIDA_SRC)/utilities.o

.PHONY: all clean
all: $(TEST_EXE)

$(TEST_EXE): $(CONFIG_OBJ) $(HDF5_OBJ) $(SAMPLING_OBJ) $(TEST_OBJ) $(FIDAOBJS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

$(BUILD_DIR):
	mkdir -p $@

$(CONFIG_OBJ): $(CONFIG_F90) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) \
		-I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(HDF5_OBJ): $(HDF5_F90) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) \
		-I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(SAMPLING_OBJ): $(SAMPLING_F90) $(HDF5_OBJ) $(FIDAOBJS) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) \
		-I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(TEST_OBJ): $(TEST_F90) $(CONFIG_OBJ) $(HDF5_OBJ) $(SAMPLING_OBJ) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) \
		-I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

clean:
	-rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod $(TEST_EXE)
