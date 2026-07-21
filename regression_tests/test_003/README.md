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

The top-level runner uses the committed reference HDF5 files. It runs the
sampling workflow in `02_run_test` and then the comparison workflow in
`03_compare`. It does not regenerate the trusted Stage 1 reference data.

## Workflow

| Stage | Purpose | Main outputs |
| --- | --- | --- |
| [`01_reference`](01_reference/README.md) | Extract trusted `f(E, pitch)` slices from a supported source distribution. | Reference HDF5 files and optional PNG plots. |
| [`02_run_test`](02_run_test/README.md) | Sample each reference distribution and reconstruct it on the same grid. | Sampled HDF5 files and optional PNG plots. |
| [`03_compare`](03_compare/README.md) | Compare reference and sampled marginals, density, parallel temperature, and perpendicular temperature. | Comparison figures and a text report. |

```text
Source fast-ion distribution
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

Stage 1 defines the common energy-pitch HDF5 schema under
[Generated output files](01_reference/README.md#generated-output-files).
Stage 2 preserves that dataset schema and metadata while replacing `f_array`
with the sampled reconstruction and adding sampling provenance. Stage 3 uses
the Stage 2 input configuration as the ordered source of file pairs.

## Directory structure

```text
test_003/
├── README.md                       # this portal
├── makefile
├── run.sh                          # run sampling and comparison
├── 01_reference/
│   ├── README.md
│   ├── input_config.nml
│   ├── generate_reference_data.py
│   ├── reference_generator_tools/
│   └── output_data/
├── 02_run_test/
│   ├── README.md
│   ├── input_config.nml
│   ├── run.sh
│   ├── plot_sampled_data.py
│   ├── src/
│   └── output_data/
└── 03_compare/
    ├── README.md
    ├── input_config.nml
    ├── run.sh
    ├── compare_distributions.py
    ├── comparison_tools/
    └── output_data/
```

The tree omits individual generated output files and compiler artifacts for
clarity.

## Navigate the test

1. [01_reference](01_reference/README.md): Generate the reference distributions
2. [02_run_test](02_run_test/README.md): Run the inverse-transform sampler
3. [03_compare](03_compare/README.md): Compare the reconstructed distributions

---

**Navigation:** [Previous: Test 002](../test_002/README.md) | [Up: Regression tests](../README.md) | [Next: Reference data](01_reference/README.md)
