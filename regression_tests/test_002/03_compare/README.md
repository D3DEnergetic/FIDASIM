[Regression tests](../../README.md) / [Test 002](../README.md) / Compare results

**Navigation:** [Previous: Run the sampler](../02_run_test/README.md) | [Up: Test 002](../README.md) | Next: —

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
cd regression_tests/test_002/03_compare
./run.sh
```

`run.sh` executes:

```bash
python3 compare_distributions.py input_config.nml
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
├── input_config.nml
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
  sampling_config_file = '../02_run_test/input_config.nml'
  output_directory = 'output_data'
  generate_plots = .true.
/

&plot_data_block
  scale = 'lin'
  fmin =
  fmax =
  enable_colorbar = .true.
  colormap = 'viridis'
/
```

Paths are resolved relative to the configuration file that contains them.
`sampling_config_file` points to the exact Stage 2 configuration used for the
sampling run. The comparison reads its ordered `reference_files` list and
`output_directory`, then pairs every reference file with the sampled file of
the same basename. It does not scan either directory for unrelated HDF5 files.

### `compare` schema

| Variable | Type | Required | Default | Allowed values and behavior |
| --- | --- | --- | --- | --- |
| `sampling_config_file` | String | Yes | None | Nonempty path to the Stage 2 configuration used for sampling. Stage 3 reads its ordered reference-file list and sampled-output directory. |
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

The plotting choices intentionally match Stage 2. This block affects the
workflow only when `generate_plots = .true.` in the `compare` block.

### Path resolution

Stage 3 resolves a relative `sampling_config_file` or `output_directory`
against the directory containing its own configuration file. When it reads the
Stage 2 configuration, it likewise resolves Stage 2 reference and output paths
against the directory containing that configuration. The meanings of these
paths therefore do not depend on the terminal's current working directory.

## Validation

Before calculating results, the script checks that each pair has:

- the shared HDF5 schema defined by Stage 1;
- matching energy and pitch grids and `f_array` dimensions;
- matching `f_array` units;
- matching species, ion, selected-position, requested-position, and grid-index
  metadata;
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

**Navigation:** [Previous: Run the sampler](../02_run_test/README.md) | [Up: Test 002](../README.md) | Next: —
