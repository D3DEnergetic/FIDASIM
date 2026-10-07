SHELL = /bin/bash

TEST_NAME := test_006
TEST_DIR := src
BUILD_DIR := build
REPORT_DIR := output_data
REPORT_FILE := $(REPORT_DIR)/comparison_report.txt
FIDASIM_DIR ?= $(abspath ../..)
FIDA_SRC := $(FIDASIM_DIR)/src
TEST_F90 := $(TEST_DIR)/$(TEST_NAME).f90
TEST_OBJ := $(BUILD_DIR)/$(TEST_NAME).o
TEST_EXE := $(BUILD_DIR)/$(TEST_NAME)$(EXEEXT)
FIDA_OBJS := $(FIDA_SRC)/libfida.o $(FIDA_SRC)/eigensystem.o \
	$(FIDA_SRC)/utilities.o $(FIDA_SRC)/hdf5_utils.o

.PHONY: all clean
all: $(TEST_EXE)

$(TEST_EXE): $(TEST_OBJ) $(FIDA_OBJS)
	$(MPI_FC) $(C_FLAGS) $^ -o $@ $(L_FLAGS)

$(TEST_OBJ): $(TEST_F90) | $(BUILD_DIR)
	$(MPI_FC) $(C_FLAGS) -c -I$(FIDA_SRC) $(I_FLAGS) $< -o $@

$(BUILD_DIR):
	mkdir -p $@

clean:
	-rm -f $(BUILD_DIR)/*.o $(TEST_EXE) $(REPORT_FILE)
