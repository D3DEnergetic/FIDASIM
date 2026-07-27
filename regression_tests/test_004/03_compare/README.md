[Regression tests](../../README.md) / [Test 004](../README.md) / Comparison

**Navigation:** [Previous: Monte Carlo calculation](../02_monte_carlo/README.md) | [Up: Test 004](../README.md) | Next: —

---

# Compare deterministic and Monte Carlo ion sinks

This stage pairs the deterministic and Monte Carlo ion-sink artifacts,
validates that they represent identical physical inputs, and compares their
total reaction rates and energy-pitch distributions.

The total-rate metrics provide the regression decision. The two-dimensional
and marginal comparisons remain diagnostic because pointwise relative errors
are poorly behaved in empty or sparsely sampled Monte Carlo bins.

## Running the comparison

Run Stage 1 and Stage 2 with the shared Test 004 configuration before running:

```bash
conda activate FIDASIM_env
cd regression_tests/test_004/03_compare
./run.sh input_config_A.nml
```

The script writes the requested plots and the text report before returning a
nonzero exit status when one or more cases fail the configured rate criteria.

## Input configuration

```fortran
&compare
  comment = 'Dataset A: deterministic and Monte Carlo ion-sink comparison'
  test_config = '../input_config_A.nml'
  rate_relative_tolerance = 0.02
  rate_sigma_tolerance = 3.0
  generate_plots = .true.
/

&plot_data_block
  scale = 'lin'
  emax = 150
  fmin = 'auto'
  fmax = 'auto'
  enable_colorbar = .true.
  colormap = 'hot_r'
  contour_levels = 30
/
```

All relative paths are resolved from the comparison configuration file.
`test_config` points to the unified configuration used to run both upstream
implementations. Stage 3 derives the ordered source collection and the
`deterministic`, `monte_carlo`, and `comparison` output directories from that
single file.

### `compare` schema

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `comment` | String scalar | — | No | Empty | Human-readable description of the comparison collection. |
| `test_config` | String scalar | — | Yes | None | Unified Test 004 configuration used by both calculation stages. |
| `rate_relative_tolerance` | Real scalar | dimensionless | Yes | None | Positive maximum absolute relative difference in total reaction rate. |
| `rate_sigma_tolerance` | Real scalar | standard errors | Yes | None | Positive maximum absolute rate difference measured in Monte Carlo standard errors. |
| `generate_plots` | Logical scalar | — | No | `.true.` | Enables the per-case distribution and marginal figures and the collection rate summary. |

### `plot_data_block` schema

The block is optional. Defaults are applied when plotting is enabled and the
block is absent.

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `scale` | String scalar | — | No | `lin` | `lin` or `log`; log mode masks zero-valued 2D bins and uses logarithmic marginal axes. |
| `emax` | Real or `auto` | `keV` | No | Automatic | Maximum displayed energy. |
| `fmin` | Real or `auto` | Sink units | No | Automatic | Shared lower color limit for each deterministic–Monte Carlo pair. |
| `fmax` | Real or `auto` | Sink units | No | Automatic | Shared upper color limit for each deterministic–Monte Carlo pair. |
| `enable_colorbar` | Logical scalar | — | No | `.true.` | Enables the shared two-dimensional colorbar. |
| `colormap` | String scalar | — | No | `viridis` | `viridis`, `viridis_r`, `hot`, or `hot_r`. |
| `contour_levels` | Integer scalar | — | No | `30` | Number of filled contour levels; must be at least 2. |

Explicit `fmin` and `fmax` values are always expressed in physical sink
units. The plotting layer performs the base-10 transformation internally
when `scale='log'`.

## Pairing and validation

For each smooth distribution selected through the unified Test 004
configuration, Stage 3 constructs the common basename
`fidasim_f4d_NNN_ion_sink.h5` beneath both implementation directories. It
requires the complete expected collection before processing any case.

The deterministic comparison contract is stored at the HDF5 root. The
production Monte Carlo schema remains at its root, while its self-contained
comparison contract is under `/test_004`. Both layouts are normalized to one
in-memory object.

Every pair must have identical:

- energy and pitch grids;
- smooth input distribution and authoritative ion density;
- isotope, charge, and selected spatial position;
- neutral energy, velocity, and six-level density;
- source-distribution, atomic-table, and configuration provenance; and
- dataset units.

Agreement between two files is not sufficient on its own because both could
be stale. Stage 3 therefore also checks their provenance, source distribution,
neutral parameters, gyrophase grid, and Monte Carlo controls against the
currently selected unified configuration and Test 002 files.

Each file is checked internally: the stored energy marginal, pitch marginal,
and total rate must reproduce the two-dimensional sink integral. The Monte
Carlo standard error must also reproduce
`sample_standard_deviation/sqrt(n_markers)`.

## Metrics and acceptance

For deterministic rate $R_{\rm det}$, Monte Carlo rate $R_{\rm MC}$, and
Monte Carlo standard error $\mathrm{SE}_{\rm MC}$, Stage 3 reports

$$
\delta_R
=
\frac{R_{\rm MC}-R_{\rm det}}{R_{\rm det}},
$$

and

$$
z_R
=
\frac{R_{\rm MC}-R_{\rm det}}{\mathrm{SE}_{\rm MC}}.
$$

A case passes only when both conditions hold:

$$
|\delta_R|
\leq
\mathtt{rate\_relative\_tolerance},
$$

$$
|z_R|
\leq
\mathtt{rate\_sigma\_tolerance}.
$$

For either energy or pitch marginal, Stage 3 also reports the integrated
absolute difference normalized by the deterministic total rate. For example,

$$
D_E
=
\frac{
\int
\left|S_{E,\rm MC}(E)-S_{E,\rm det}(E)\right|\ dE
}{
R_{\rm det}
}.
$$

These marginal distances describe shape and normalization disagreement but
do not currently determine pass/fail.

## Outputs

Generated files are written beneath
`save_data_block/output_directory/comparison`:

- `<stem>_marginals.png` overlays the deterministic and Monte Carlo energy
  and pitch marginals;
- `<stem>_distributions.png` shows both 2D sinks with a shared color scale and
  the injected-neutral marker;
- `reaction_rate_summary.png` compares rates and relative differences across
  all cases; and
- `comparison_report.txt` begins with a prominent overall `PASS` or `FAIL`
  banner, then records all criteria, metrics, and per-case details.

## Code organization

```text
03_compare/
├── compare_ion_sinks.py
├── input_config_A.nml
├── run.sh
└── comparison_tools/
    ├── __init__.py
    ├── config.py
    ├── data.py
    ├── analysis.py
    ├── plotting.py
    ├── reporting.py
    └── workflow.py
```

The command-line program is deliberately thin. `workflow.py` exposes the
ordered comparison steps, while the supporting modules isolate configuration,
HDF5 contracts, numerical metrics, plotting, and reporting.

---

**Navigation:** [Previous: Monte Carlo calculation](../02_monte_carlo/README.md) | [Up: Test 004](../README.md) | Next: —
