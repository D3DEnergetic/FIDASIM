# FIDASIM Python Unit Tests

This directory contains unit tests for the FIDASIM Python utilities and preprocessing modules.

## Prerequisites

The tests require Python 3.6+ and the following packages:
- pytest (testing framework)
- numpy
- scipy
- h5py
- matplotlib

Install test requirements:
```bash
pip install -r requirements.txt
```

Or using the FIDASIM Python:
```bash
$FIDASIM_DIR/deps/python -m pip install -r requirements.txt
```

## Test Structure

### test_utils.py
Comprehensive tests for the `fidasim.utils` module:

- **TestCOCOS**: Tests for COCOS (COordinate COnventions) class
  - Default initialization
  - Various COCOS indices (1, 2, 5, 11, etc.)
  - Sign conventions

- **TestGetFidasimDir**: Tests for FIDASIM directory detection

- **TestGetVersion**: Tests for version retrieval
  - Git-based version
  - Fallback to VERSION file

- **TestCoordinateTransforms**: Tests for coordinate transformations
  - Tait-Bryan angles (tb_zyx)
  - uvw ↔ xyz transformations
  - Identity transformations
  - Inverse relationship verification

- **TestLineBasis**: Tests for line basis construction
  - Orthonormal basis from line direction
  - Various orientations

- **TestRZGrid**: Tests for R-Z grid creation
  - Basic 2D grids
  - 3D grids with phi dimension
  - Grid spacing and ranges

- **TestBeamGrid**: Tests for neutral beam grid creation
  - Grid orientation and basis
  - Grid dimensions and spacing

- **TestColoredOutput**: Tests for colored terminal output
  - info, warn, error, success functions
  - ANSI color codes

- **TestWriteData**: Tests for HDF5 data writing
  - Scalar, array, and matrix data
  - Attributes (descriptions, units)

- **TestReadNCDF**: Tests for NetCDF file reading (with mocking)

- **TestAABBIntersect**: Tests for axis-aligned bounding box ray intersection
  - Hit detection
  - Miss detection

### test_preprocessing.py
Tests for the `fidasim.preprocessing` module:

- **TestCheckDictSchema**: Tests for dictionary validation
  - Valid dictionaries
  - Missing keys
  - Wrong types
  - Wrong shapes
  - Extra keys (warnings)
  - NaN/Inf detection

- **TestCheckInputs**: Tests for input parameter validation
  - Complete valid inputs
  - Missing required fields
  - Path validation

- **TestCheckPlasma**: Tests for plasma profile validation
  - Valid plasma dictionaries
  - Shape consistency checks
  - Required fields

- **TestCheckFields**: Tests for electromagnetic field validation
  - Valid field dictionaries
  - Missing fields
  - Shape validation

- **TestCheckDistribution**: Tests for distribution function validation
  - Type 1: Fast-ion distribution function
  - Type 2: Guiding center Monte Carlo
  - Required arrays and dimensions

- **TestCheckNBI**: Tests for neutral beam injection validation
  - Beam geometry parameters
  - Aperture specifications
  - Array dimensions

- **TestCheckSpec**: Tests for spectroscopy chord validation
  - Channel configurations
  - Required fields

- **TestPrefidaIntegration**: Integration tests for the main prefida function
  - Minimal valid configuration
  - File I/O mocking
  - End-to-end validation

## Running Tests

### Run all Python tests:
```bash
make test-python
```

### Run specific test file:
```bash
cd test/python_tests
python -m pytest test_utils.py -v
```

### Run specific test class:
```bash
python -m pytest test_utils.py::TestCOCOS -v
```

### Run specific test method:
```bash
python -m pytest test_utils.py::TestCOCOS::test_cocos_defaults -v
```

### Run with coverage report:
```bash
python -m pytest --cov=fidasim --cov-report=html
```
This creates an HTML coverage report in `htmlcov/index.html`

### Run with different verbosity levels:
```bash
# Quiet mode
python -m pytest -q

# Normal mode
python -m pytest

# Verbose mode
python -m pytest -v

# Very verbose mode
python -m pytest -vv
```

### Run tests matching a pattern:
```bash
python -m pytest -k "cocos" -v  # Run all tests with "cocos" in the name
```

## Test Output

Tests use pytest's assertion introspection for clear failure messages. When a test fails, pytest shows:
- The assertion that failed
- The actual vs expected values
- The full traceback (can be controlled with --tb option)

Example output:
```
test_utils.py::TestCOCOS::test_cocos_defaults PASSED
test_utils.py::TestCOCOS::test_cocos_index_1 PASSED
test_preprocessing.py::TestCheckDictSchema::test_valid_dict PASSED
```

## Writing New Tests

1. Create test functions/methods starting with `test_`
2. Use pytest fixtures for setup/teardown
3. Use `pytest.raises` for exception testing
4. Use `unittest.mock` for mocking external dependencies
5. Use `numpy.testing.assert_*` for numerical comparisons

Example:
```python
def test_my_function():
    """Test description"""
    # Arrange
    input_data = create_test_data()

    # Act
    result = my_function(input_data)

    # Assert
    assert result.shape == (10, 20)
    np.testing.assert_array_almost_equal(result, expected)
```

## Continuous Integration

These tests can be integrated into CI/CD pipelines:
```yaml
# Example GitHub Actions workflow
- name: Run Python tests
  run: |
    pip install -r test/python_tests/requirements.txt
    make test-python
```
## Troubleshooting
### Import Errors
If you get import errors, ensure:
1. FIDASIM_DIR environment variable is set
2. Python path includes lib/python directory
3. Dependencies are installed
### Test Discovery Issues
pytest automatically discovers tests in files matching:
- `test_*.py` or `*_test.py`
- Test functions/methods starting with `test_`

### Mocking Issues
Some tests use mocking to avoid file I/O. If mocks fail:
1. Check mock patch paths match actual import paths
2. Ensure mocked methods match actual signatures

## Clean Up

Remove test artifacts:
```bash
make clean_tests
```

Or manually:
```bash
rm -rf __pycache__ .pytest_cache
find . -name "*.pyc" -delete
```
