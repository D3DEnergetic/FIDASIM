[Regression tests](../../README.md) / [Test 003](../README.md) / Run sampler

**Navigation:** [Previous: Reference data](../01_reference/README.md) | [Up: Test 003](../README.md) | [Next: Compare results](../03_compare/README.md)

---

# Energy-pitch distribution sampler

This stage samples the native Test 002 Stage 2 energy-pitch distributions
using FIDASIM's `randind`, `randu`, and serial random-number infrastructure.
It reads the singleton spatial slice from each Test 002 `f` dataset directly;
no copied Test 003 reference files are used.

## Run

Generate the Test 002 Stage 2 outputs, build the regression executables, and
run:

```bash
cd regression_tests/test_003/02_run_test
./run.sh input_config_A.nml
```

The Python normalization step discovers the contiguous indexed Test 002 files
and writes their absolute paths into the private namelist consumed by Fortran.

## Input configuration

```fortran
&run_test
  comment = 'Dataset A: 60 keV neutral beam at 45 degrees'
  input_distribution_config = '../../test_002/02_run_test/input_config_A.nml'
  output_directory = 'output_data/dataset_A'
  n_samples = 1000000
  seed = 12345
  plot_data = .true.
/
```

### `run_test` schema

| Variable | Type | Required | Description |
| --- | --- | --- | --- |
| `comment` | String | No | Human-readable collection description. |
| `input_distribution_config` | String | Yes | Test 002 Stage 2 configuration whose indexed outputs are sampled. |
| `output_directory` | String | Yes | Directory for sampled HDF5 files and plots. |
| `n_samples` | Integer (Int64) | Yes | Positive number of samples drawn from each distribution. |
| `seed` | Integer (Int32) | Yes | Positive deterministic serial RNG seed. |
| `plot_data` | Logical | No | Generate sampled-distribution plots; default `.false.`. |

Relative paths are resolved from this configuration. The indexed Test 002
filenames must begin at `001` and be contiguous.

The existing `plot_data_block` accepts `scale`, `fmin`, `fmax`,
`enable_colorbar`, `colormap`, and `emax`.

## HDF5 contracts

Reference inputs use the native Test 002 schema:

- `energy`, `pitch`;
- `f(nz=1,nr=1,npitch,nenergy)` and `denf`;
- `r`, `z`;
- `species`, `atomic_number`, `mass_number`, `charge_state`, and `A`.

Sampled outputs use a compact rank-two schema:

- `energy_grid(nenergy)`, `pitch_grid(npitch)`;
- `f_array(nenergy,npitch)` and reconstructed `denf`;
- `selected_r`, `selected_z`;
- the five copied species-parameter datasets.

The sampled filename matches its Test 002 reference basename.

---

**Navigation:** [Previous: Reference data](../01_reference/README.md) | [Up: Test 003](../README.md) | [Next: Compare results](../03_compare/README.md)
