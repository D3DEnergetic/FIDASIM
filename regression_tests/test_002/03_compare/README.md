[Regression tests](../../README.md) / [Test 002](../README.md) / Compare the moments

**Navigation:** [Previous: Run the conversion](../02_run_test/README.md) | [Up: Test 002](../README.md) | Next: —

---

# Compare moments before and after conversion

Stage 3 determines whether converting the CQL3D distribution from
$f(\bar{u},\theta)$ to the FIDASIM representation $F(E,P)$ preserves density,
parallel temperature, and perpendicular temperature.

The workflow reads the physical moments directly from the Stage 1 and Stage 2
HDF5 files. The text reports from those stages are intended for inspection and
are not used as numerical inputs.

## Run the comparison

The Stage 1 and Stage 2 workflows must have been run first. From this directory,
run:

```bash
./run.sh
```

This is equivalent to:

```bash
python3 compare_moments.py input_config.nml
```

## Input configuration

```fortran
&compare
  run_config = '../02_run_test/input_config.nml'
  output_directory = 'output_data'
  generate_plot = .true.
/
```

### `&compare` schema

| Field | Type | Meaning |
| --- | --- | --- |
| `run_config` | string | Path to the Stage 2 configuration. Relative paths are resolved from the Stage 3 configuration file. |
| `output_directory` | string | Directory for the comparison report and plot. Relative paths are resolved from the Stage 3 configuration file. |
| `generate_plot` | logical | Generate the relative-difference summary plot when `.true.`. |

The Stage 2 configuration identifies both its converted output basename and
the Stage 1 configuration. Stage 3 follows these links to construct every
reference-converted file pair; the file list is therefore not duplicated in
the comparison configuration.

## Comparison

For each quantity $q$, the absolute relative difference is

\[
\epsilon_q=
\left|
\frac{q_{\mathrm{converted}}-q_{\mathrm{reference}}}
{q_{\mathrm{reference}}}
\right|.
\]

Before calculating the differences, each pair is checked for matching species
and selected $(R,Z)$ location. Stage 3 does not compare the distribution arrays
point by point because the source and converted arrays use different coordinate
systems and grids.

## Generated outputs

| File | Description |
| --- | --- |
| `output_data/moment_comparison.txt` | Reference value, converted value, and relative difference for all three moments and every case, followed by the maximum observed differences. |
| `output_data/moment_values.png` | Reference and converted density, $T_\parallel$, and $T_\perp$ versus case index. |
| `output_data/moment_relative_differences.png` | Density, $T_\parallel$, and $T_\perp$ relative differences versus case index. |

Case index is used instead of spatial position because the selected reference
locations may be nonconsecutive or may not follow a one-dimensional spatial
scan.

No pass/fail tolerance is imposed at this stage. The report records the
observed error introduced by coordinate transformation and interpolation.

---

**Navigation:** [Previous: Run the conversion](../02_run_test/README.md) | [Up: Test 002](../README.md) | Next: —
