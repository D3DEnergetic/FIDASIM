[Regression tests](../README.md) / Test 003

**Navigation:** [Previous: Test 002](../test_002/README.md) | [Up: Regression tests](../README.md) | [Next: Reference data](01_reference/README.md)

---

# Test 003: energy-pitch distribution sampling

`test_003` checks FIDASIM's statistical sampling of a two-dimensional fast-ion
distribution, `f(E, pitch)`, where `pitch = v_parallel/v`. The test starts from
trusted distributions, samples them with the same random-selection primitives
used by FIDASIM, reconstructs them on their original grids, and compares
physically meaningful results.

The test targets the sampling primitive itself: the two-dimensional discrete
inverse-transform selection performed by `randind`, the uniform within-bin
sampling performed with `randu`, and FIDASIM's serial RNG infrastructure. It
does not exercise the complete guiding-centre or plasma integration path.

## FIDASIM code under test

| Module | Public procedure | Resolved implementation | Description |
| --- | --- | --- | --- |
| `utilities` | `rng_init` | `rng_init` | Initializes the deterministic serial RNG stream before each reference distribution. |
| `utilities` | `randind` | `randind_w_2` | Performs weighted two-dimensional inverse-transform selection of an energy-pitch bin from `f_array`. |
| `utilities` | `randu` | `randu_arr` | Supplies the uniform deviates used to place the sampled energy and pitch within the selected bin. |

The cumulative-distribution and RNG-state helpers called by these procedures
are implementation details; they are not independently tested by this
regression workflow.

The test reproduces the energy-pitch sampling statements used by
`libfida::mc_sample_ion_f4d_gc`, but deliberately does not call that complete
routine. HDF5 reading, validation, histogramming, normalization, and plotting
belong to the regression-test harness rather than the FIDASIM sampling code
under test.

## Quick start

The complete runner requires a Conda environment containing these Python
packages:

- `f90nml`
- `h5py`
- `numpy`
- `matplotlib`

Activate that environment before building and running the test. For example,
if these dependencies are installed in `FIDASIM_env`:

```bash
conda activate FIDASIM_env

cd "$FIDASIM_DIR"
make regression_tests

cd regression_tests/test_003
./run.sh
```

If a different Conda environment contains the required packages, activate
that environment instead.

The top-level runner executes reference validation, sampling, and comparison
for each configuration listed in its short `config_files` array. Each listed
filename must exist in all three stages.

## Dataset collections

The same letter identifies a collection across Test 002 and Test 003. Test 003
Stage 1 points to the matching Test 002 Stage 2 configuration and writes its
compact sampler references under the corresponding output directory.

| Collection | Output directory | Description |
| --- | --- | --- |
| `A` | `dataset_A` | 60 keV neutral beam at 45°, inherited from Test 002 collection A. |

The descriptive `comment` field inside each configuration makes its physical
meaning visible without encoding those details in the filename.

## Workflow

| Stage | Purpose | Main outputs |
| --- | --- | --- |
| [`01_reference`](01_reference/README.md) | Validate the shared Test 002 Stage 2 distribution collection. | Validation summary; no copied data. |
| [`02_run_test`](02_run_test/README.md) | Sample each reference distribution and reconstruct it on the same grid. | Sampled HDF5 files and optional PNG plots. |
| [`03_compare`](03_compare/README.md) | Compare reference and sampled marginals, density, parallel temperature, and perpendicular temperature. | Comparison figures and a text report. |

```text
Correct Test 002 FIDASIM distributions
            │
            ▼
      01_reference
    Trusted f(E, pitch)
            │
            ▼
       02_run_test
  Sampled reconstructions
            │
            ▼
       03_compare
 Marginals, moments, report
```

## Data contract

Test 002 Stage 2 owns the native reference schema. Test 003 Stage 2 reads that
schema directly and writes compact sampled reconstructions. Stage 3 reads both
layouts and compares their common physical data and metadata.

Test 003 depends on the generated Test 002 Stage 2 outputs. Its Stage 1
validator makes this dependency explicit and checks the collection before use.

## Directory structure

```text
test_003/
├── README.md                       # this portal
├── makefile
├── run.sh                          # run sampling and comparison for every listed collection
├── 01_reference/
│   ├── README.md
│   ├── input_config_A.nml
│   ├── validate_reference_data.py
│   └── reference_generator_tools/
├── 02_run_test/
│   ├── README.md
│   ├── input_config_A.nml
│   ├── run.sh
│   ├── plot_sampled_data.py
│   ├── src/
│   └── output_data/dataset_A/
└── 03_compare/
    ├── README.md
    ├── input_config_A.nml
    ├── run.sh
    ├── compare_distributions.py
    ├── comparison_tools/
    └── output_data/dataset_A/
```

The tree omits individual generated output files and compiler artifacts for
clarity.

## Navigate the test

1. [01_reference](01_reference/README.md): Generate the reference distributions
2. [02_run_test](02_run_test/README.md): Run the inverse-transform sampler
3. [03_compare](03_compare/README.md): Compare the reconstructed distributions

---

**Navigation:** [Previous: Test 002](../test_002/README.md) | [Up: Regression tests](../README.md) | [Next: Reference data](01_reference/README.md)
