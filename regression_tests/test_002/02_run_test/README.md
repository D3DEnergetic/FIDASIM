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

Running from `02_run_test` is required because the namelist paths are relative
to the working directory. `run.sh` executes:

```bash
./test_002 input_config.nml
python3 plot_sampled_data.py input_config.nml
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
├── plot_sampled_data.py
├── run.sh
├── src/
│   ├── test_002.f90
│   ├── test_002_hdf5.f90
│   └── test_002_sampling.f90
└── output_data/                 # generated HDF5 files and PNG plots
```

The executable is built as `02_run_test/test_002`. Object and module files are
kept under `src/`.

## Input configuration

The configuration uses Fortran namelist syntax. Paths are interpreted relative
to the directory from which the program is run; the supported workflow runs it
from `02_run_test`.

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

### Run-test fields

| Field | Meaning |
| --- | --- |
| `n_reference_files` | Positive number of reference files, up to 256. |
| `reference_files` | Explicit ordered list of reference HDF5 paths. The first `n_reference_files` entries are used. |
| `output_directory` | Directory for sampled HDF5 files and plots. It is created when needed. |
| `n_samples` | Positive Int64 sample count used independently for every reference distribution. |
| `seed` | Positive deterministic Int32 RNG seed. |
| `plot_data` | Enables or disables the Python plotting stage. |

The serial RNG stream is reinitialized with the same `seed` before sampling
each reference file. A file's results therefore do not depend on its position
in the list or on which other files are processed.

### Plot-data fields

| Field | Meaning |
| --- | --- |
| `scale` | `lin` plots `f`; `log` plots `log10(f)` for positive values. |
| `fmin` | Lower color limit. Leave empty or use `auto` for the data minimum. |
| `fmax` | Upper color limit. Leave empty or use `auto` for the data maximum. |
| `enable_colorbar` | Enables or disables the color bar. |
| `colormap` | `viridis`, `viridis_r`, `hot`, or `hot_r`. |

The plotter processes every `.h5` file in `output_directory` and writes a PNG
with the same basename.

## HDF5 data contract

The reference input and sampled output use the same energy-pitch distribution
schema documented under
[Generated output files](../01_reference/README.md#generated-output-files) in
the `01_reference` README. Grid and particle metadata, datatypes, units, and
descriptions are copied from the reference file. The sampled file replaces
only the numerical contents and description of `f_array`.

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
