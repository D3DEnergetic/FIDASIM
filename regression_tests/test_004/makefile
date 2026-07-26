SHELL = /bin/bash

TEST_NAME  := test_004
TEST_DIR   := 02_monte_carlo
SRC_DIR    := $(TEST_DIR)/src
MODULE_DIR := $(SRC_DIR)/modules
BUILD_DIR  := $(TEST_DIR)/build
FIDASIM_DIR ?= $(abspath ../..)
FIDA_SRC   := $(FIDASIM_DIR)/src

MPI_FC  ?= gfortran
C_FLAGS ?= -O0 -g -fbacktrace -fcheck=all -Wall -Wextra -std=f2008 -fopenmp
I_FLAGS ?= -I$(FIDASIM_DIR)/deps/hdf5/include -I$(FIDA_SRC)
L_FLAGS ?= -L$(FIDASIM_DIR)/deps/hdf5/lib -Wl,-Bstatic \
	-lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 \
	-Wl,-Bdynamic -lz -ldl -lm -fopenmp

TYPES_F90  := $(MODULE_DIR)/$(TEST_NAME)_types.f90
CONFIG_F90 := $(MODULE_DIR)/$(TEST_NAME)_config.f90
HDF5_UTILS_F90 := $(MODULE_DIR)/$(TEST_NAME)_hdf5_utils.f90
HDF5_F90   := $(MODULE_DIR)/$(TEST_NAME)_hdf5.f90
NEUTRAL_F90 := $(MODULE_DIR)/$(TEST_NAME)_neutral.f90
SETUP_F90 := $(MODULE_DIR)/$(TEST_NAME)_setup.f90
SAMPLING_F90 := $(MODULE_DIR)/$(TEST_NAME)_sampling.f90
TEST_F90   := $(SRC_DIR)/$(TEST_NAME).f90

TYPES_OBJ  := $(BUILD_DIR)/$(TEST_NAME)_types.o
CONFIG_OBJ := $(BUILD_DIR)/$(TEST_NAME)_config.o
HDF5_UTILS_OBJ := $(BUILD_DIR)/$(TEST_NAME)_hdf5_utils.o
HDF5_OBJ   := $(BUILD_DIR)/$(TEST_NAME)_hdf5.o
NEUTRAL_OBJ := $(BUILD_DIR)/$(TEST_NAME)_neutral.o
SETUP_OBJ := $(BUILD_DIR)/$(TEST_NAME)_setup.o
SAMPLING_OBJ := $(BUILD_DIR)/$(TEST_NAME)_sampling.o
TEST_OBJ   := $(BUILD_DIR)/$(TEST_NAME).o
TEST_EXE   := $(TEST_DIR)/$(TEST_NAME)$(EXEEXT)

FIDA_OBJS := $(FIDA_SRC)/libfida.o $(FIDA_SRC)/eigensystem.o \
	$(FIDA_SRC)/utilities.o $(FIDA_SRC)/hdf5_utils.o

.PHONY: all clean
all: $(TEST_EXE)

$(TEST_EXE): $(TYPES_OBJ) $(CONFIG_OBJ) $(HDF5_UTILS_OBJ) $(HDF5_OBJ) \
		$(NEUTRAL_OBJ) \
		$(SETUP_OBJ) $(SAMPLING_OBJ) $(TEST_OBJ) $(FIDA_OBJS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

$(BUILD_DIR):
	mkdir -p $@

$(TYPES_OBJ): $(TYPES_F90) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(CONFIG_OBJ): $(CONFIG_F90) $(TYPES_OBJ) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(HDF5_UTILS_OBJ): $(HDF5_UTILS_F90) $(FIDA_OBJS) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(HDF5_OBJ): $(HDF5_F90) $(TYPES_OBJ) $(HDF5_UTILS_OBJ) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(NEUTRAL_OBJ): $(NEUTRAL_F90) $(TYPES_OBJ) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(SETUP_OBJ): $(SETUP_F90) $(TYPES_OBJ) $(FIDA_OBJS) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(FIDA_SRC) \
		-I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(SAMPLING_OBJ): $(SAMPLING_F90) $(TYPES_OBJ) $(FIDA_OBJS) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(FIDA_SRC) \
		-I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

$(TEST_OBJ): $(TEST_F90) $(TYPES_OBJ) $(CONFIG_OBJ) $(HDF5_UTILS_OBJ) \
		$(HDF5_OBJ) \
		$(NEUTRAL_OBJ) $(SETUP_OBJ) $(SAMPLING_OBJ) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) $(I_FLAGS) -c -I$(BUILD_DIR) -J$(BUILD_DIR) $< -o $@

clean:
	-rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod $(TEST_EXE)
