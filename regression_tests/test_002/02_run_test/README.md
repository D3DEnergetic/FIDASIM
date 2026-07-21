[Regression tests](../../README.md) / [Test 002](../README.md) / Run sampler

**Navigation:** [Previous: Reference data](../01_reference/README.md) | [Up: Test 002](../README.md) | [Next: Compare results](../03_compare/README.md)

---

# Energy-pitch distribution sampler

This stage applies FIDASIM's two-dimensional inverse-transform sampler to the
energy-pitch distributions produced by `01_reference`. It uses `randind`,
`randu`, and FIDASIM's serial random-number infrastructure, bins the sampled
energy and pitch values back onto the reference grids, and writes one
reconstructed HDF5 file per reference file.

## Quick start

Activate a Python environment containing the plotting dependencies listed
below, then build from the FIDASIM repository root and run from `02_run_test`:

```bash
cd "$FIDASIM_DIR"
make regression_tests

cd regression_tests/test_002/02_run_test
./run.sh
```

`run.sh` is the supported wrapper around the Fortran executable. Run it from
`02_run_test`; it validates and normalizes the user configuration before
starting Fortran, then runs the optional plotting stage. In outline, it
executes:

```bash
python3 normalize_config.py input_config.nml build/normalized_input_config.nml
./test_002 build/normalized_input_config.nml
python3 plot_sampled_data.py input_config.nml
```

To use another configuration file, pass its path to the wrapper:

```bash
./run.sh path/to/input_config.nml
```

## Dependencies

Build the test through FIDASIM's top-level makefile. It exports the configured
Fortran compiler and flags (`MPI_FC`, `C_FLAGS`, `I_FLAGS`, and `L_FLAGS`) to
the regression-test makefiles. These settings supply the FIDASIM module path,
the selected bundled or system HDF5 include and library paths, and the required
linker libraries. The test also links the already-built
`FIDASIM/src/utilities.o` object to use `rng`, `rng_init`, `randind`, and
`randu`.

Plotting requires Python with:

- `f90nml`
- `h5py`
- `numpy`
- `matplotlib`

Activate an environment containing these packages before running `run.sh`.

## Directory layout

```text
02_run_test/
├── input_config.nml
├── normalize_config.py
├── plot_sampled_data.py
├── run.sh
├── build/                        # generated object and module files
├── src/
│   ├── test_002.f90
│   └── modules/
│       ├── test_002_config.f90
│       ├── test_002_hdf5.f90
│       └── test_002_sampling.f90
└── output_data/                 # generated HDF5 files and PNG plots
```

The executable is built as `02_run_test/test_002`. The main program is kept at
the top of `src/`, reusable test modules are under `src/modules/`, and compiler
artifacts are isolated in `build/`.

## Input configuration

The user configuration uses Fortran namelist syntax. The Python normalization
step interprets relative paths from the directory containing this file and
writes absolute paths to `build/normalized_input_config.nml`. The generated
file contains only the `run_test` block consumed by Fortran.

```fortran
&run_test
  n_reference_files = 4
  reference_files =
    '../01_reference/output_data/f_array_001.h5',
    '../01_reference/output_data/f_array_002.h5',
    '../01_reference/output_data/f_array_003.h5',
    '../01_reference/output_data/f_array_004.h5'
  output_directory = 'output_data'
  n_samples = 1000000
  seed = 12345
  plot_data = .true.
/

&plot_data_block
  scale = 'lin'
  fmin =
  fmax =
  enable_colorbar = .true.
  colormap = 'viridis'
/
```

### `run_test` schema

| Variable | Type | Required | Default | Allowed values and behavior |
| --- | --- | --- | --- | --- |
| `n_reference_files` | Integer (Int32) | Yes | None | Number of reference files. Must be from 1 through 256. |
| `reference_files` | String array | Yes | None | Ordered list of reference HDF5 paths. The first `n_reference_files` entries must all be nonempty. |
| `output_directory` | String | Yes | None | Nonempty directory path for sampled HDF5 files and plots. It is created when needed. |
| `n_samples` | Integer (Int64) | Yes | None | Number of samples drawn independently from each reference distribution. Must be positive. |
| `seed` | Integer (Int32) | Yes | None | Deterministic serial RNG seed. Must be positive. |
| `plot_data` | Logical | No | `.false.` | Enables or disables the Python plotting stage. |

The serial RNG stream is reinitialized with the same `seed` before sampling
each reference file. A file's results therefore do not depend on its position
in the list or on which other files are processed.

### `plot_data_block` schema

| Variable | Type | Required | Default | Allowed values and behavior |
| --- | --- | --- | --- | --- |
| `scale` | String | No | `'lin'` | `'lin'` plots `f`; `'log'` plots `log10(f)` for positive values. Choices are case-insensitive. |
| `fmin` | Real or string | No | Automatic | Lower color limit. Leave empty or use `'auto'` to use the minimum of the plotted values. |
| `fmax` | Real or string | No | Automatic | Upper color limit. Leave empty or use `'auto'` to use the maximum of the plotted values. |
| `enable_colorbar` | Logical | No | `.true.` | Enables or disables the color bar. |
| `colormap` | String | No | `'viridis'` | Must be `'viridis'`, `'viridis_r'`, `'hot'`, or `'hot_r'`. Choices are case-insensitive. |

This block is consumed by `plot_sampled_data.py` and affects the workflow only
when `plot_data = .true.` in the `run_test` block.

The plotter processes every `.h5` file in `output_directory` and writes a PNG
with the same basename.

### Path resolution

The Python wrapper resolves relative paths from the directory containing the
user configuration file. It validates that every reference file exists and
passes a generated namelist containing absolute paths to Fortran. The plotter
applies the same path rule when it reads the original configuration. Absolute
paths are accepted unchanged.

## HDF5 data contract

The reference input and sampled output use the same energy-pitch distribution
schema documented under
[Generated output files](../01_reference/README.md#generated-output-files) in
the `01_reference` README. Grid and particle metadata, datatypes, units, and
descriptions are copied from the reference file. The sampled file replaces
the numerical contents and description of `f_array`. It also replaces `denf`
with the density calculated from the sampled energy-pitch distribution:

```text
denf = sum(f_array) * abs(denergy * dpitch)
```

Stage 3 does not compare this value yet. The density-preservation check will
be enabled after the CQL3D-to-FIDASIM distribution convention is established
by the dedicated conversion regression test.

Before sampling, this stage additionally verifies that:

- One-dimensional `energy_grid` and `pitch_grid` datasets with at least two
  finite values are present.
- `f_array` has Python/HDF5 shape `(nenergy, npitch)`.
- The distribution values are finite and nonnegative, with a positive sum.
- Both grids are strictly monotonic and uniformly spaced. Increasing and
  decreasing grids are supported.

The HDF5 Fortran interface reports the two-dimensional file shape in reversed
order. The reader checks this explicitly and stores the Fortran array as
`f_array(energy_index, pitch_index)`.

The sampled file adds root-level provenance:

- `description`
- `data_source_type = sampled_reconstruction`
- `data_source_name` and `reference_file`
- `n_samples`
- `rng_seed`

The bundled HDF5 1.8 Fortran interface cannot write an Int64 attribute, so
`n_samples` is stored as its full decimal string. `rng_seed` is stored as a
numeric Int32 attribute.

## Sampling and normalization

The method is a two-dimensional discrete inverse-transform sampler. `randind`
constructs the cumulative distribution from the flattened 2D `f_array`, draws
a weighted energy-pitch bin, and converts the selected flat index back into its
two array indices. `randu` then places the sample uniformly within that bin.

For each sample, the program performs the same operations used by FIDASIM's
nonthermal energy-pitch sampler:

```fortran
call randind(f_array, ep_ind)
call randu(randomu3)

energy = energy_grid(ep_ind(1,1)) + denergy * (randomu3(1) - 0.5)
pitch  = pitch_grid(ep_ind(2,1))   + dpitch  * (randomu3(2) - 0.5)
```

The sampled values are assigned to the reference bins. Histogram counts are
stored as Int64, every sample must fall inside the grid domain, and the final
count sum must equal `n_samples`.

For uniform bins, the sampled distribution is reconstructed in the reference
units using:

```text
f_sampled(i,j) = counts(i,j) / n_samples * sum(f_reference)
```

Equivalently, this preserves
`sum(f) * abs(denergy * dpitch)`, so the sampled and reference total integrals
match apart from floating-point rounding.

## Output files

Each sampled file uses the reference basename:

```text
01_reference/output_data/f_array_001.h5
02_run_test/output_data/f_array_001.h5
02_run_test/output_data/f_array_001.png
```

The HDF5 file follows the shared data contract above. The PNG is generated only
when `plot_data = .true.`.

---

**Navigation:** [Previous: Reference data](../01_reference/README.md) | [Up: Test 002](../README.md) | [Next: Compare results](../03_compare/README.md)
