[Regression tests](../README.md) / Test 004

**Navigation:** [Previous: Test 003](../test_003/README.md) | [Up: Regression tests](../README.md) | [Next: Deterministic calculation](01_deterministic/README.md)

---

# Test 004: Charge Exchange Ion Sink Verification Test

`test_004` validates the charge-exchange ion-sink calculation performed by
FIDASIM. This calculation appears in both `calculate_dcx_process` and
`calculate_halo_process`. The test takes the steps in
`calculate_dcx_process` as the representative workflow and exercises them in
an intentionally stripped-down calculation.

The test isolates ion sampling, charge-exchange reactivity, and ion-sink
storage. Neutral tracking, attenuation, photon production, and subsequent
neutral and ion-birth storage are outside its scope.

## Workflow

| Stage | Purpose |
| --- | --- |
| [`01_deterministic`](01_deterministic/README.md) | Integrate the smooth distributions deterministically over energy, pitch, and gyrophase. |
| [`02_monte_carlo`](02_monte_carlo/README.md) | Run the stripped-down Monte Carlo ion-sink calculation using FIDASIM procedures. |
| [`03_compare`](03_compare/README.md) | Compare the deterministic and Monte Carlo distributions, marginals, and total rates. |

## FIDASIM code under test

| Module | Procedure | Role |
| --- | --- | --- |
| `libfida` | `mc_sample_ion_f4d_gc` | Samples ion position and velocity from the configured distribution. |
| `libfida` | `get_total_cx_rate` | Evaluates the reservoir-averaged, level-resolved CX rate. |
| `libfida` | `store_sinks` | Accumulates the cell-resolved sink-rate density. |
| `libfida` | `store_sink_particle` | Records each weighted ion-sink particle. |
| `libfida` | `write_sink_profile` | Writes the normal FIDASIM sink HDF5 result. |

## Problem statement

The energy-pitch-resolved ion-sink distribution is

$$
S(E,p)=
n_i\ g(E,p)
\sum_{l,m=1}^{6}n_l
\left\langle K_{m\leftarrow l}\right\rangle_\phi(E,p),
$$

and the total volumetric ion-sink reaction rate is

$$
R=\int S(E,p)\ dE\ dp.
$$

Here $n_i$ is the ion density, $g(E,p)$ is the normalized ion distribution,
and $\left\langle K_{m\leftarrow l}\right\rangle_\phi(E,p)$ is the
state-resolved, gyrophase-averaged charge-exchange reactivity:

$$
\left\langle K_{m\leftarrow l}\right\rangle_\phi(E,p)=
\left\langle
\sigma_{m\leftarrow l}(\varepsilon_{\rm rel})v_{\rm rel}
\right\rangle_{\phi}.
$$

For any gyrophase-dependent quantity $q(\phi)$, the gyro-averaging operator
is

$$
\left\langle q\right\rangle_\phi
=\frac{1}{2\pi}\int_0^{2\pi}q(\phi)\ d\phi.
$$

The index $l$ identifies the initial atomic energy level of the target
neutral and $n_l$ is its number density. The index $m$ identifies the
atomic energy level populated by charge exchange.
$\sigma_{m\leftarrow l}$ is directional: it is the cross section for

$$
\mathrm{H}^{+}+\mathrm{H}(l)
\longrightarrow
\mathrm{H}(m)+\mathrm{H}^{+}.
$$

For ion mass $M_i$ and neutral mass $M_n$, the physical centre-of-mass
collision energy is

$$
E_{\rm rel}=\frac{1}{2}\mu v_{\rm rel}^{2},
\quad
\mu=\frac{M_iM_n}{M_i+M_n}.
$$

FIDASIM's H-H cross-section tables use the relative collision energy per
reduced-mass atomic mass unit:

$$
\varepsilon_{\rm rel}
=\frac{E_{\rm rel}}{\mu/m_u}
=\frac{1}{2}m_u v_{\rm rel}^{2},
$$

expressed in `keV/amu`, where $m_u$ is the atomic mass constant. The relative
speed is

$$
v_{\rm rel}
=\left|\boldsymbol v_i(E,p,\phi)-\boldsymbol v_n\right|.
$$

Finally, the normalized ion distribution is

$$
g(E,p)=\frac{f(E,p)}
{\int f(E,p)\ dE\ dp},
\qquad
\int g(E,p)\ dE\ dp=1.
$$

The deterministic and Monte Carlo implementations evaluate these same
quantities for matching ion and neutral isotopes.

## Data dependency

Test 002 Stage 1 contains the committed CQL3D single-slice distributions.
Test 002 Stage 2 converts them into generated, smooth FIDASIM
energy-pitch distributions. Those converted files are the common inputs to
both Test 004 implementations:

> The **Test 2, Stage 2** smooth distribution files are the **input data files**
> for the **Test 4** workflow:

The workflow is illustrated below:

```text
Test 002 Stage 1 CQL3D slices (git-tracked reference files)
                |
                v
Test 002 Stage 2 smooth F(E,p)
                |
        +-------+-------+
        |               |
        v               v
01_deterministic   02_monte_carlo
        |               |
        +-------+-------+
                |
                v
           03_compare
```

The two calculation implementations (`01_deterministic`, `02_monte_carlo`)
are independent peers and may run in either order.

Their HDF5 files and plots are generated test artifacts and are
not committed. The comparison stage is the only stage that requires both
output collections.

## Shared input configuration

Both calculation implementations share
[`input_config_A.nml`](input_config_A.nml). The `test_case` and `neutrals`
blocks define common physical inputs;
`deterministic` and `monte_carlo` contain implementation-specific controls.
The shared `plot_data_block` gives the two principal sink plots consistent
styling, while `save_data_block/output_directory` defines their common output
root.

All relative paths are resolved from the directory containing
`input_config_A.nml`. Selector values are case-insensitive.

### `test_case` schema

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `comment` | String scalar | — | No | Empty | Human-readable description of the common physical test case. |
| `input_distribution_config` | String scalar | — | Yes | None | Path to the Test 002 Stage 2 configuration that identifies the smooth energy-pitch distributions. |
| `tables_filename` | String scalar | — | Yes | None | Path to the FIDASIM atomic-tables HDF5 file used by both implementations. |

### `neutrals` schema

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `density` | Real scalar | `cm^-3` | Yes | None | Total neutral density; must be positive. |
| `energy` | Real scalar | `keV` | Yes | None | Total kinetic energy of each neutral; must be positive. |
| `injection_angle` | Real scalar | degrees | Yes | None | Signed angle from `+z` toward `+x`; the velocity is proportional to `[sin(theta),0,cos(theta)]`. |
| `level_split_method` | String scalar | — | Yes | None | `ground-only` or `exponential`. |
| `level_decay` | Real scalar | — | Conditional | None | Required and positive for `exponential`; ignored for `ground-only`. |

`ground-only` assigns the total density to level 1. `exponential` assigns
normalized fractions proportional to
$\exp[-\mathtt{level\_decay}(l-1)]$ over FIDASIM's six atomic levels. The
neutral isotope is taken from each source distribution and must match its ion
isotope.

### `plot_data_block` schema

This block is required when an enabled implementation requests plotting.

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `scale` | String scalar | — | No | `lin` | `lin` or `log`. |
| `emax` | Real or `auto` | `keV` | No | Automatic | Maximum displayed energy; it must be positive and greater than the distribution's minimum energy. This affects only the plot. |
| `fmin` | Real or `auto` | Sink units | No | Automatic | Lower plotted color limit. `auto` selects it from the data; an explicit value must be positive for `log` scale. |
| `fmax` | Real or `auto` | Sink units | No | Automatic | Upper plotted color limit. `auto` selects it from the data; an explicit value must be positive for `log` scale. |
| `enable_colorbar` | Logical scalar | — | No | `.true.` | Enables the plot color bar. |
| `colormap` | String scalar | — | No | `viridis` | `viridis`, `viridis_r`, `hot`, or `hot_r`. |
| `contour_levels` | Integer scalar | — | No | `100` | Number of filled contour levels; must be at least 2. |

Explicit `fmin` and `fmax` values always use physical sink units. For a log
plot, the plotting layer converts them to base-10 logarithms internally.

### `save_data_block` schema

This block is required when an enabled implementation requests HDF5 or plot
output.

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `output_directory` | String scalar | — | Yes | None | Common output root. Deterministic and Monte Carlo files are written beneath `deterministic/` and `monte_carlo/`, respectively. |

The implementation-specific schemas are documented with the
[deterministic interface](01_deterministic/README.md#input-configuration) and
[Monte Carlo interface](02_monte_carlo/README.md#input-configuration).
The comparison stage owns a
[separate configuration and schema](03_compare/README.md#input-configuration)
under `03_compare`.

## How Test 004 finds and validates its input distributions

The variable `input_distribution_config` in the `test_case` block of the
input configuration file identifies a Test 002 Stage 2 configuration rather than
an individual HDF5 distribution. Both Test 004 implementations use that
configuration to construct the same ordered input collection:

1. Read the Test 002 Stage 2 configuration selected by
   `input_distribution_config`.
2. Read its `save_data_block/output_filename`, which defines the converted
   output basename.
3. Find the corresponding three-digit indexed files, such as
   `fidasim_f4d_001.h5`, `fidasim_f4d_002.h5`, and so on.
4. Require the indices to begin at `001` and remain contiguous.
5. Read `energy`, `pitch`, `f`, `denf`, spatial coordinates, and species
   metadata directly from each HDF5 file.
6. Require the species parameters to be consistent across the complete
   collection.

No Test 002 Stage 1 configuration lookup is required because each converted
HDF5 file contains the necessary location and species metadata.

If the converted files have not yet been generated, run:

```bash
cd regression_tests/test_002/02_run_test
./run.sh input_config_A.nml
```

---

**Navigation:** [Previous: Test 003](../test_003/README.md) | [Up: Regression tests](../README.md) | [Next: Deterministic calculation](01_deterministic/README.md)
