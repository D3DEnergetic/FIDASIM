# FIDASIM Unit Tests

This directory contains unit tests for FIDASIM using the pFUnit testing framework.

## Prerequisites

- pFUnit 4.12 (already installed in `/home/gwin/FIDASIM/deps/pFUnit/`)
- CMake 3.12 or higher
- gfortran compiler
- HDF5 libraries

## Test Coverage

### test_utilities.pf
Tests for the utilities module including:
- `ind2sub` and `sub2ind` - Array indexing conversions
- `cumsum` - Cumulative sum function
- Random number generator initialization and generation
- Sparse array creation and value retrieval

### test_eigensystem.pf
Tests for the eigensystem module including:
- `RSWAP` - Value swapping
- `balance` - Matrix balancing
- Matrix inversion
- Matrix multiplication

### test_fidasim.pf
Tests for the main fidasim module including:
- `approx_eq`, `approx_ge`, `approx_le` - Approximate comparison functions
- `cross_product` - Vector cross product
- Coordinate transformations (xyz ↔ uvw, xyz ↔ cylindrical)
- Grid boundary checking (`in_grid`)
- Grid indexing (`get_indices`)
- Basis transformations (Tait-Bryan angles, line basis, plane basis)
- Physical constants validation

## Running Tests

### Using Make

From the FIDASIM root directory:
```bash
make test
```

This will:
1. Build the FIDASIM source files if needed
2. Configure the test build with CMake
3. Compile the test executables
4. Run all tests using CTest

### Manual Testing

From this directory:
```bash
mkdir -p build
cd build
cmake ..
make
ctest --output-on-failure
```

### Running Individual Tests

To run a specific test suite:
```bash
cd build
ctest -R test_utilities --output-on-failure
ctest -R test_eigensystem --output-on-failure
ctest -R test_fidasim --output-on-failure
```

### Verbose Output

For detailed test output:
```bash
ctest -V
```

## Adding New Tests

1. Create a new `.pf` file with your tests using pFUnit macros:
   - `@test` - Marks a test subroutine
   - `@assertEqual` - Assert equality with optional tolerance
   - `@assertTrue` / `@assertFalse` - Boolean assertions

2. Update `CMakeLists.txt` to include your new test file

3. Update `testSuites.inc` if needed

4. Example test structure:
```fortran
@test
subroutine test_my_function()
    use my_module
    use pfunit_mod
    implicit none
    
    real(8) :: result, expected
    
    result = my_function(1.0d0)
    expected = 2.0d0
    
    @assertEqual(expected, result, tolerance=1.0d-10)
    
end subroutine test_my_function
```

## Cleaning Test Artifacts

To remove all test build files:
```bash
make clean_tests
```

Or manually:
```bash
rm -rf build/
```
