cmake_minimum_required(VERSION 3.12)
project(FIDASIM_Tests Fortran)

# Find pFUnit
set(PFUNIT_DIR "/home/gwin/FIDASIM/deps/pFUnit/install/PFUNIT-4.12")
find_package(PFUNIT REQUIRED HINTS ${PFUNIT_DIR})

# Add pFUnit module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Set up includes and libraries
include_directories(
    ${PFUNIT_DIR}/include
    ${CMAKE_SOURCE_DIR}/../../src
    ${CMAKE_SOURCE_DIR}/../../deps/hdf5/include
)

link_directories(
    ${CMAKE_SOURCE_DIR}/../../src
    ${CMAKE_SOURCE_DIR}/../../deps/hdf5/lib
)

# Enable testing
enable_testing()

# Add test sources
set(TEST_SOURCES
    test_utilities.pf
    test_eigensystem.pf
    test_fidasim.pf
)

# Create test executable for utilities
add_pfunit_ctest(test_utilities
    TEST_SOURCES test_utilities.pf
    LINK_LIBRARIES utilities eigensystem hdf5_fortran hdf5
)

# Create test executable for eigensystem
add_pfunit_ctest(test_eigensystem
    TEST_SOURCES test_eigensystem.pf
    LINK_LIBRARIES eigensystem
)

# Create test executable for fidasim
add_pfunit_ctest(test_fidasim
    TEST_SOURCES test_fidasim.pf
    LINK_LIBRARIES fidasim utilities eigensystem hdf5_utils hdf5_fortran hdf5hl_fortran hdf5_hl hdf5
)
