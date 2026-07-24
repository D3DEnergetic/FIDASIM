[Regression tests](../../README.md) / [Test 003](../README.md) / Reference data

**Navigation:** Previous: — | [Up: Test 003](../README.md) | [Next: Run the sampler](../02_run_test/README.md)

---

# Generate sampling references from Test 002

Stage 1 adapts the correctly transformed FIDASIM distributions produced by
Test 002 into compact two-dimensional fixtures for the Test 003 Fortran
sampler. It does not transform, interpolate, or renormalize the distributions.

The selected Test 002 Stage 2 configuration is the single scientific input.
Through that file, this stage discovers the ordered output files and inherits
the requested locations and particle metadata from Test 002 Stage 1.

## Regenerate the trusted references

First run the selected Test 002 Stage 2 workflow so its indexed
`fidasim_f4d_*.h5` files exist. Then run:

```bash
cd regression_tests/test_003/01_reference
./run.sh input_config_A.nml
```

The regenerated HDF5 files and diagnostic PNGs are trusted artifacts and are
committed. Ordinary Test 003 runs begin at Stage 2 and use these committed
fixtures; they do not rerun Test 002 or regenerate Stage 1.

## Input configuration

```fortran
&input
  comment = 'Dataset A: 60 keV neutral beam at 45 degrees'
  input_config = '../../test_002/02_run_test/input_config_A.nml'
  plot_data = .true.
  save_data = .true.
/

&plot_data_block
  scale = 'lin'
  fmin = 'auto'
  fmax = 'auto'
  enable_colorbar = .true.
  colormap = 'hot_r'
  emax = 150.0
/

&save_data_block
  output_filename = 'output_data/dataset_A/f_array.h5'
/
```

Relative paths are resolved from the directory containing the configuration
file. To use a different collection of CQL3D conditions, generate it through
Test 002 and change only `input_config` to the corresponding Test 002 Stage 2
configuration.

### `&input` schema

| Variable | Type | Required | Meaning |
| --- | --- | --- | --- |
| `comment` | String | No | Human-readable description of the dataset collection. |
| `input_config` | String | Yes | Test 002 Stage 2 configuration that identifies the source output collection. |
| `plot_data` | Logical | No | Generate diagnostic PNGs; defaults to `.false.`. |
| `save_data` | Logical | No | Generate compact HDF5 fixtures; defaults to `.false.`. |

### `&plot_data_block` schema

| Variable | Type | Required | Meaning |
| --- | --- | --- | --- |
| `scale` | String | No | `lin` or `log`; defaults to `lin`. |
| `fmin` | Real or `auto` | No | Lower color limit; defaults to the data minimum. |
| `fmax` | Real or `auto` | No | Upper color limit; defaults to the data maximum. |
| `enable_colorbar` | Logical | No | Enable the colorbar; defaults to `.true.`. |
| `colormap` | String | No | `viridis`, `viridis_r`, `hot`, or `hot_r`. |
| `emax` | Real | No | Maximum displayed energy in keV; defaults to `150.0` and must be positive. |

### `&save_data_block` schema

`output_filename` is required whenever HDF5 or PNG output is enabled. It must
have a `.h5` extension. A three-digit case index is inserted before the
extension.

## Data discovery and validation

The adapter follows `reference_config` in the selected Test 002 Stage 2
configuration. The Test 002 Stage 1 configuration supplies:

- case count and ordering;
- requested R and Z locations;
- species, atomic number, mass number, and charge state.

For every indexed Test 002 Stage 2 output, the adapter requires:

- one R location and one Z location;
- one-dimensional energy and pitch grids with at least two points;
- `f` with shape `(1, 1, npitch, nenergy)`;
- finite, nonnegative distribution values;
- a finite, positive `denf`;
- a species attribute matching the Test 002 Stage 1 configuration.

The single spatial slice is converted only by array selection and transpose:

```python
f_array = f[0, 0, :, :].T
```

Thus, `f_array` has shape `(nenergy, npitch)`. No values are rescaled.

## Generated HDF5 schema

| Dataset | Dimensions | Units | Description |
| --- | --- | --- | --- |
| `energy_grid` | `(nenergy)` | `keV` | Cell-centred energy grid. |
| `pitch_grid` | `(npitch)` | Dimensionless | Cell-centred pitch grid. |
| `f_array` | `(nenergy, npitch)` | `ions/(cm^3*keV*dP)` | Distribution consumed by the sampler. |
| `denf` | `(1)` | `ions/cm^3` | Density supplied by Test 002. |
| `species` | Scalar | n/a | Species inherited from Test 002. |
| `atomic_number` | Scalar | Dimensionless | Proton count inherited from Test 002. |
| `mass_number` | Scalar | Dimensionless | Nucleon count inherited from Test 002. |
| `charge_state` | Scalar | Elementary charge | Charge state inherited from Test 002. |
| `requested_r`, `requested_z` | `(1)` | `cm` | Locations requested by Test 002 Stage 1. |
| `selected_r`, `selected_z` | `(1)` | `cm` | Locations stored in the Test 002 output. |
| `r_index`, `z_index` | `(1)` | Dimensionless | Indices of the singleton spatial axes. |

Root attributes identify the Test 002 output used for each fixture. Every
dataset also contains `units` and `description` attributes.

## Code organization

- `config.py` validates this stage's input, plotting, and output controls.
- `workflow.py` follows the Test 002 configurations, validates each source,
  and writes the compact fixtures and plots.
- `readers/fidasim_h5.py` reads the internal FIDASIM HDF5 schema.

Test 002 owns the CQL3D-to-FIDASIM transformation. New source formats or
conversion methods belong there rather than in this adapter.

---

**Navigation:** Previous: — | [Up: Test 003](../README.md) | [Next: Run the sampler](../02_run_test/README.md)
