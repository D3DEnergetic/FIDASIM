[Regression tests](../../README.md) / [Test 003](../README.md) / Compare results

**Navigation:** [Previous: Run the sampler](../02_run_test/README.md) | [Up: Test 003](../README.md) | Next: —

---

# Compare reference and sampled distributions

This stage compares every reference distribution sampled by `02_run_test`
with its reconstructed output. It avoids pointwise relative errors, which are
poorly behaved where the reference distribution approaches zero. Instead, it
compares energy and pitch marginals and three integrated physical quantities.

## Quick start

Run the comparison from this directory in a Python environment containing the
dependencies listed below:

```bash
cd regression_tests/test_003/03_compare
./run.sh input_config_A.nml
```

`run.sh` executes:

```bash
./run.sh input_config_A.nml
```

## Dependencies

- `f90nml`
- `h5py`
- `numpy`
- `matplotlib`

## Code organization

```text
03_compare/
├── compare_distributions.py       # thin command-line entry point
├── input_config_A.nml
├── run.sh
├── comparison_tools/
│   ├── __init__.py                # public run_comparison API
│   ├── config.py                  # configuration and file pairing
│   ├── data.py                    # HDF5 reading and validation
│   ├── analysis.py                # marginals, moments, and errors
│   ├── plotting.py                # comparison figures
│   ├── reporting.py               # text report
│   └── workflow.py                # complete comparison sequence
└── output_data/                   # generated figures and report
```

The command-line script only reads the configuration argument and calls
`comparison_tools.run_comparison`. The five main workflow steps remain visible
in `workflow.py`; the supporting modules hide the implementation details for
each step.

## Input configuration

```fortran
&compare
  comment = 'Dataset A: 60 keV neutral beam at 45 degrees'
  sampling_config_file = '../02_run_test/input_config_A.nml'
  output_directory = 'output_data/dataset_A'
  generate_plots = .true.
/

&plot_data_block
  scale = 'lin'
  fmin =
  fmax =
  enable_colorbar = .true.
  colormap = 'viridis'
  emax = 150.0
/
```

Paths are resolved relative to the configuration file that contains them.
`sampling_config_file` points to the exact Stage 2 configuration used for the
sampling run. The comparison follows its `input_distribution_config`,
discovers the contiguous native Test 002 references, and pairs each with the
sampled file of the same basename.

### `compare` schema

| Variable | Type | Required | Default | Allowed values and behavior |
| --- | --- | --- | --- | --- |
| `comment` | String | No | None | Human-readable description of the dataset collection. |
| `sampling_config_file` | String | Yes | None | Stage 2 configuration used for sampling. Stage 3 discovers the same Test 002 references and sampled-output directory from it. |
| `output_directory` | String | Yes | None | Nonempty directory path for the comparison report and figures. It is created when needed. |
| `generate_plots` | Logical | No | `.true.` | Enables or disables both comparison figures for every file pair. The text report is always written. |

### `plot_data_block` schema

| Variable | Type | Required | Default | Allowed values and behavior |
| --- | --- | --- | --- | --- |
| `scale` | String | No | `'lin'` | `'lin'` plots the original values; `'log'` uses logarithmic marginal y-axes and plots `log10(f)` for positive 2D distribution values. Choices are case-insensitive. |
| `fmin` | Real or string | No | Automatic | Lower limit of the shared 2D color scale. Leave empty or use `'auto'` to use the minimum across both distributions. It does not set a marginal-axis limit. |
| `fmax` | Real or string | No | Automatic | Upper limit of the shared 2D color scale. Leave empty or use `'auto'` to use the maximum across both distributions. It does not set a marginal-axis limit. |
| `enable_colorbar` | Logical | No | `.true.` | Enables or disables the shared color bar on the 2D comparison figure. |
| `colormap` | String | No | `'viridis'` | Must be `'viridis'`, `'viridis_r'`, `'hot'`, or `'hot_r'`. Choices are case-insensitive. |
| `emax` | Real | No | `150.0` | Maximum displayed energy in keV for the marginal and two-dimensional comparison figures. Must be positive. |

The plotting choices intentionally match Stage 2. This block affects the
workflow only when `generate_plots = .true.` in the `compare` block.

### Path resolution

Stage 3 resolves each relative path against the configuration file containing
it. Path meanings therefore do not depend on the terminal's working directory.

## Validation

Before calculating results, the script checks that each pair has:

- the native Test 002 reference schema and compact Test 003 sampled schema;
- matching energy and pitch grids and `f_array` dimensions;
- matching `f_array` units;
- matching species parameters and selected-position metadata;
- finite, nonnegative distribution values with a positive integral; and
- strictly monotonic, uniformly spaced energy and pitch grids.

## Quantities compared

For uniform grid spacings `dE` and `dP`, the script calculates:

```text
f_E(E) = sum over pitch of f(E,P) * abs(dP)
f_P(P) = sum over energy of f(E,P) * abs(dE)

density = sum f(E,P) * abs(dE*dP)

T_parallel = 2/density * sum E*P^2*f(E,P) * abs(dE*dP)
T_perpendicular = 1/density * sum E*(1-P^2)*f(E,P) * abs(dE*dP)
```

It reports absolute relative errors for density, `T_parallel`, and
`T_perpendicular`. No acceptance threshold or pass/fail decision is imposed;
the report provides diagnostic values for review.

## Output files

For each configured pair, two figures are written:

- `<basename>_marginals.png` overlays reference and sampled `f_E(E)` in the
  first panel and reference and sampled `f_P(P)` in the second.
- `<basename>_distributions.png` places reference and sampled `f(E,P)` side by
  side using the same colormap and exactly the same color limits.

`comparison_report.txt` records the reference value, sampled value, and
relative error for each physical quantity, followed by the maximum relative
error across all file pairs.

Generated files are written under `output_data/` and are excluded from Git.

---

**Navigation:** [Previous: Run the sampler](../02_run_test/README.md) | [Up: Test 003](../README.md) | Next: —
