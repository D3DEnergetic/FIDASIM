[Regression tests](../../README.md) / [Test 003](../README.md) / Reference data

**Navigation:** Previous: — | [Up: Test 003](../README.md) | [Next: Run the sampler](../02_run_test/README.md)

---

# 4D distribution reference-data generator

This tool ingests a supported four-dimensional distribution and extracts a
two-dimensional `f(E, pitch)` slice at each requested `(r, z)` location. It
writes every slice to a separate HDF5 file and can create an associated PNG
plot.

The generated HDF5 files provide trusted reference data for subsequent
regression-test work. Source formats may use different velocity-space
coordinates or dimension ordering; their readers convert them to the common
energy-pitch output representation.

## Quick start

The tool requires Python with these packages:

- `f90nml`
- `numpy`
- `h5py`
- `matplotlib`

For example, create and activate a Conda environment:

```bash
conda create -n reference_data_env -c conda-forge \
  python numpy h5py matplotlib f90nml
conda activate reference_data_env
```

Run the tool from its directory:

```bash
cd regression_tests/test_003/01_reference
python generate_reference_data.py input_config.nml
```

Set `input_filename` in [input_config.nml](input_config.nml) to the relative or
absolute path of the source distribution.

## Input configuration

The input configuration file uses standard Fortran namelist syntax:

- Blocks begin with `&name` and end with `/`.
- Comments begin with `!`.
- Strings use single or double quotes.
- Logical values use `.true.` or `.false.`.
- Arrays use comma-separated values.

```fortran
&input
  input_file_type = 'fidasim_h5'
  input_filename = 'input_data/distribution.h5'
  species = 'D'
  atomic_number = 1
  mass_number = 2
  charge_state = 1
  r_locations = 0.0, 0.1
  z_locations = 0.0, 10.0
  plot_data = .true.
  save_data = .true.
/

&plot_data_block
  scale = 'lin'
  fmin =
  fmax =
  enable_colorbar = .true.
  colormap = 'viridis'
/

&save_data_block
  output_filename = 'output_data/f_array.h5'
/
```

Supported string choices are case-insensitive.

### Path resolution

**Every relative path is resolved relative to the directory containing the
input configuration file.** It is not resolved relative to the terminal's
current directory. This gives each relative path a stable meaning.

For example, if the configuration is:

```text
/project/reference/input_config.nml
```

then:

```fortran
input_filename = 'input_data/distribution.h5'
output_filename = 'output_data/f_array.h5'
```

resolves to:

```text
/project/reference/input_data/distribution.h5
/project/reference/output_data/f_array.h5
```

Absolute paths remain absolute, and `~` is expanded to the user home
directory.

### Input block

| Variable | Type | Required | Units | Allowed values and behavior |
| --- | --- | --- | --- | --- |
| `input_file_type` | string | Yes | n/a | `fidasim_h5` or `cql3d_f4d`. See supported formats below. |
| `input_filename` | string | Yes | n/a | Absolute path or path relative to the input configuration file. |
| `species` | string | Yes | n/a | `H`, `D`, or `T`; stored in lowercase. |
| `atomic_number` | integer | Yes | n/a | Positive proton count, $Z$. |
| `mass_number` | integer | Yes | n/a | Proton-plus-neutron count, $A$, greater than or equal to `atomic_number`. |
| `charge_state` | integer | Yes | elementary charge | Integer from zero through `atomic_number`. |
| `r_locations` | real array | Yes | `cm` | Radial locations; must match the length of `z_locations`. |
| `z_locations` | real array | Yes | `cm` | Axial locations; must match the length of `r_locations`. |
| `plot_data` | logical | No | n/a | Defaults to `.false.`; enables PNG output. |
| `save_data` | logical | No | n/a | Defaults to `.false.`; enables HDF5 output. |

At least one `(r, z)` pair is required. The nearest spatial grid point is used;
the tool does not interpolate.

### Plot-data block

This block is validated when `plot_data = .true.`.

| Variable | Type | Required | Allowed values and behavior |
| --- | --- | --- | --- |
| `scale` | string | No | `lin` or `log`; defaults to `lin`. Log scale plots $\log_{10}(f)$ for positive values. |
| `fmin` | real | No | Lower plot limit. Empty or `auto` uses the data minimum. |
| `fmax` | real | No | Upper plot limit. Empty or `auto` uses the data maximum. |
| `enable_colorbar` | logical | No | Defaults to `.true.`. |
| `colormap` | string | No | `viridis`, `viridis_r`, `hot`, or `hot_r`; defaults to `viridis`. `_r` reverses the color order. |

Plot titles contain `f(E, pitch)` and the selected R and Z grid values in cm.

### Save-data block

This block is validated when `save_data = .true.`.

| Variable | Type | Required | Allowed values and behavior |
| --- | --- | --- | --- |
| `output_filename` | string | Yes | Path with a required `.h5` extension. May be absolute or relative to the input configuration file. |

The output directory is created automatically. A numeric suffix is inserted
before the extension:

```text
output_data/f_array_001.h5
output_data/f_array_002.h5
```

## Supported input formats

### `fidasim_h5`

This reader is implemented. It requires the following datasets:

| Dataset | Type | Dimensions | Units | Description |
| --- | --- | --- | --- | --- |
| `z` | Numeric, read as Float64 | `(nz)` | `cm` | Axial grid. |
| `r` | Numeric, read as Float64 | `(nr)` | `cm` | Radial grid. |
| `pitch` | Numeric, read as Float64 | `(npitch)` | Dimensionless | Pitch grid, $v_\parallel/v$. |
| `energy` | Numeric, read as Float64 | `(nenergy)` | `keV` | Fast-ion energy grid. |
| `f` | Numeric, read as Float64 | `(nz, nr, npitch, nenergy)` | `fast-ions/(dE*dP*cm^3)` | Distribution on the four grids. |
| `denf` | Numeric, read as Float64 | `(nz, nr)` | `cm^-3` | Fast-ion density on the spatial grid. |

The reader verifies that the shapes of `f` and `denf` match the corresponding
grid lengths. Additional datasets such as `r2d`, `z2d`, and `time` are
ignored.

The required input datasets do not reliably identify the particle species,
atomic number, mass number, or charge state. Supply these values in the
`&input` block; the generator records them in every output file.

### `cql3d_f4d`

This input type is recognized, but its reader has not been implemented.
Selecting it raises:

```text
NotImplementedError: The input distribution reader for 'cql3d_f4d' has not been implemented.
```

Its schema will be documented when the reader is implemented.

## Generated output files

For each requested location, the tool selects the nearest spatial indices,
extracts the corresponding distribution, and writes it with dimensions
`(nenergy, npitch)`.

### Dataset schema

| Dataset | HDF5 type | Dimensions | Units | Description |
| --- | --- | --- | --- | --- |
| `energy_grid` | Float64 | `(nenergy)` | `keV` | Energy coordinates of `f_array`. |
| `pitch_grid` | Float64 | `(npitch)` | Dimensionless | Pitch coordinates of `f_array`. |
| `f_array` | Float64 | `(nenergy, npitch)` | `fast-ions/(dE*dP*cm^3)` | Distribution at the selected spatial grid point. |
| `denf` | Float64 | `(1)` | `cm^-3` | Fast-ion density supplied at the selected spatial grid point. |
| `species` | UTF-8 string | Scalar | n/a | Normalized species label. |
| `atomic_number` | Integer | Scalar | Dimensionless | Number of protons, $Z$. |
| `mass_number` | Integer | Scalar | Dimensionless | Number of protons and neutrons, $A$. |
| `charge_state` | Integer | Scalar | Elementary charge | Particle charge state. |
| `requested_r` | Float64 | `(1)` | `cm` | Requested radial position. |
| `requested_z` | Float64 | `(1)` | `cm` | Requested axial position. |
| `selected_r` | Float64 | `(1)` | `cm` | Selected radial grid value. |
| `selected_z` | Float64 | `(1)` | `cm` | Selected axial grid value. |
| `r_index` | Integer | `(1)` | Dimensionless | Zero-based radial grid index. |
| `z_index` | Integer | `(1)` | Dimensionless | Zero-based axial grid index. |

Every dataset has `units` and `description` string attributes. Root attributes
contain file-level provenance:

| Root attribute | Type | Description |
| --- | --- | --- |
| `data_source_type` | String | Input distribution type. |
| `data_source_name` | String | Absolute source-file path. |
| `description` | String | Description of the generated file. |

### Reading an output file

```python
import h5py

with h5py.File("f_array_001.h5", "r") as h5f:
    distribution = h5f["f_array"][:]
    energy_units = h5f["energy_grid"].attrs["units"]
    species = h5f["species"].asstr()[()]
    mass_number = h5f["mass_number"][()]
    source_file = h5f.attrs["data_source_name"]
```

## Developer notes

### Reader interface

Each reader may interpret its source velocity coordinates and dimension order
differently, but it must return:

```python
z, r, pitch, energy, f, denf
```

with this canonical representation:

| Value | Dimensions | Representation |
| --- | --- | --- |
| `z` | `(nz)` | Axial grid in `cm` |
| `r` | `(nr)` | Radial grid in `cm` |
| `pitch` | `(npitch)` | Dimensionless pitch grid |
| `energy` | `(nenergy)` | Energy grid in `keV` |
| `f` | `(nz, nr, npitch, nenergy)` | Distribution on the canonical grids |
| `denf` | `(nz, nr)` | Fast-ion density on the spatial grid |

The reader owns any coordinate conversion, unit conversion, and array
reordering needed by its source format.

### Code organization

```text
reference_generator_tools/
├── __init__.py
├── config.py
├── workflow.py
└── readers/
    ├── __init__.py
    ├── dispatcher.py
    └── fidasim_h5.py
```

- `config.py` handles the input configuration.
- `workflow.py` selects slices, writes HDF5 files, and creates plots.
- `readers/dispatcher.py` selects the reader for `input_file_type`.
- `readers/fidasim_h5.py` implements the `fidasim_h5` reader.

### Adding another reader

1. Add a reader module that returns the canonical reader interface.
2. Add its name to `SUPPORTED_INPUT_FILE_TYPES` in `config.py`.
3. Add its selection branch to `load_input_distribution()` in
   `readers/dispatcher.py`.
4. Document the source schema under “Supported input formats.”

---

**Navigation:** Previous: — | [Up: Test 003](../README.md) | [Next: Run the sampler](../02_run_test/README.md)
