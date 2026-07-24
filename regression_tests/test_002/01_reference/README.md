[Regression tests](../../README.md) / [Test 002](../README.md) / Reference data

**Navigation:** Previous: — | [Up: Test 002](../README.md) | [Next: Run the conversion](../02_run_test/README.md)

---

# CQL3D reference distributions

This stage will extract selected `(R,Z)` locations from a large CQL3D F4D
NetCDF file. Each output is a compact HDF5 file containing only the coordinates,
integration weights, distribution, and metadata needed by this regression
test. The selected distribution is stored directly as a two-dimensional array.

## Requirements

Run this workflow from the `FIDASIM_env` Conda environment defined in
`regression_tests/environment.yml`. To create it from the repository root:

```bash
cd regression_tests
conda env create -f environment.yml
conda activate FIDASIM_env
```

The source CQL3D F4D NetCDF file must be available locally at the path given by
`input_filename`.

## Run Stage 1

From the Stage 1 directory:

```bash
cd regression_tests/test_002/01_reference
./run.sh input_config_A.nml
```

`run.sh` passes the selected configuration to `generate_reference_data.py`.
Paths are resolved relative to that configuration file.

## Configuration

The configuration contains three namelist blocks:

- `input` identifies the source format, particle, and requested locations.
- `plot_data_block` controls diagnostic contours.
- `save_data_block` controls the indexed HDF5 output names.

The first supported `input_file_type` is `cql3d_f4d`. Additional source
formats can be added through format-specific readers without changing the
user-facing workflow.

### `input` block

| Field | Required | Description |
|---|---|---|
| `comment` | No | Human-readable description of the dataset collection |
| `input_file_type` | Yes | Input format; currently `cql3d_f4d` |
| `input_filename` | Yes | CQL3D F4D NetCDF file |
| `species` | Yes | Fast-ion species: `H`, `D`, or `T` |
| `atomic_number` | Yes | Number of protons; positive integer |
| `mass_number` | Yes | Number of protons and neutrons; not smaller than `atomic_number` |
| `charge_state` | Yes | Charge state from zero through `atomic_number` |
| `r_locations` | Yes | Requested radial locations in cm |
| `z_locations` | Yes | Requested axial locations in cm; one per radial location |
| `plot_data` | No | Generate diagnostic PNG files; default `.false.` |
| `save_data` | Yes | Must be `.true.` so Stage 1 can store its reference fixtures and moments |

### `save_data_block`

This block is required by the Stage 1 workflow.

| Field | Required | Description |
|---|---|---|
| `output_filename` | Yes | Base `.h5` or `.hdf5` path; a three-digit case number is inserted before the extension |

For example, `output_data/dataset_A/cql3d_f4d.h5` produces
`output_data/dataset_A/cql3d_f4d_001.h5`,
`cql3d_f4d_002.h5`, and so on in the same dataset directory.

## Generated output files

For every configured location, Stage 1 writes:

- One indexed HDF5 reference file containing the distribution and moments.
- One same-basename PNG when `plot_data = .true.`.

After all locations are processed, `reference_moments.txt` is written beside
the indexed files and summarizes their density, parallel temperature, and
perpendicular temperature.

## Reference HDF5 schema

Each selected location produces one indexed HDF5 file with the following
datasets:

| Dataset | Shape | Description |
|---|---:|---|
| `mnemonic` | scalar | CQL3D run identifier |
| `version` | scalar | CQL3D version identifier |
| `enorm` | scalar | CQL3D energy normalization in keV |
| `k` | scalar | CQL3D species number |
| `u_norm` | scalar | Proper-velocity normalization in cm/s |
| `selected_r` | scalar | Selected radial coordinate in cm |
| `selected_z` | scalar | Selected axial coordinate in cm |
| `u_bar` | `(nu,)` | Normalized proper-velocity grid |
| `theta` | `(ntheta,)` | Pitch-angle grid in radians |
| `u_bar_weight` | `(nu,)` | Integration weights $\bar{u}^2d\bar{u}$ |
| `theta_weight` | `(ntheta,)` | Integration weights $2\pi\sin\theta\,d\theta$ |
| `f_u_theta` | `(ntheta, nu)` | Distribution at the selected location, in `ions*u_norm**3/(cm**3*(cm/sec)**3)` |
| `moments/density` | scalar | Fast-ion density in ions/cm³ |
| `moments/parallel_temperature` | scalar | Relativistic parallel pressure-equivalent temperature in keV |
| `moments/perpendicular_temperature` | scalar | Relativistic perpendicular pressure-equivalent temperature in keV |

The file attributes record the source path, requested location, selected source
indices, species, atomic number, mass number, and charge state.

The HDF5 moment datasets retain full numerical precision and are the
authoritative inputs for later comparisons. `reference_moments.txt` presents
the same values in a concise human-readable form and is not parsed by the
comparison workflow.

## Diagnostic plots

When `plot_data = .true.`, Stage 1 produces one PNG beside each reference HDF5
file. The stored coordinates are transformed according to

\[
\bar{u}_{\parallel}=\bar{u}\cos\theta,
\qquad
\bar{u}_{\perp}=\bar{u}\sin\theta.
\]

The horizontal axis is the signed normalized parallel proper velocity, while
the vertical axis is the nonnegative normalized perpendicular proper velocity.
Each figure is annotated with the selected location and the stored density,
parallel temperature, and perpendicular temperature.

The `plot_data_block` accepts the following settings:

| Field | Required | Default | Description |
|---|---|---|---|
| `scale` | No | `lin` | Plot using `lin` or base-10 `log` values |
| `fmin` | No | `auto` | Lower color limit |
| `fmax` | No | `auto` | Upper color limit |
| `enable_colorbar` | No | `.true.` | Show or hide the colorbar |
| `colormap` | No | `viridis` | Matplotlib colormap |
| `contour_levels` | No | `100` | Number of filled contour levels; must be at least 2 |

## Analytical formulation

### Notation

| Symbol | Meaning |
|---|---|
| $\bar{u}$ | Normalized proper velocity stored in `f4dv` |
| $u_{\mathrm{norm}}$ | Proper-velocity normalization stored in `vnorm` |
| $u=\bar{u}u_{\mathrm{norm}}$ | Dimensional proper velocity |
| $v$ | Ordinary particle velocity |
| $p=mu$ | Relativistic momentum magnitude |
| $\theta$ | Pitch angle stored in `f4dt` |
| $\gamma$ | Lorentz factor |
| $f(\bar{u},\theta)$ | Scaled CQL3D distribution |

### Coordinates and density

Let $f(\bar{u},\theta)$ be the scaled CQL3D gyrotropic distribution expressed
in normalized proper velocity $\bar{u}$ and pitch angle $\theta$. Integration
over gyrophase gives the normalized velocity-space element

\[
d^3\bar{u}
=2\pi \bar{u}^2\sin(\theta)\,d\bar{u}\,d\theta.
\]

The number density is therefore

\[
n=\int f(\bar{\boldsymbol{u}})\,d^3\bar{u}
=2\pi\int_0^\infty\int_0^\pi
f(\bar{u},\theta)\bar{u}^2\sin(\theta)\,
d\theta\,d\bar{u}.
\]

### Relativistic momentum and velocity

CQL3D stores the normalized proper velocity $\bar{u}$ in `f4dv`. The
dimensional proper velocity $u$, particle momentum, Lorentz factor, and
ordinary velocity are

\[
u=\bar{u}u_{\mathrm{norm}},
\qquad
p=mu,
\]

\[
\gamma=\sqrt{1+\left(\frac{u}{c}\right)^2},
\qquad
v=\frac{u}{\gamma}.
\]

Here, $m$ is the ion mass and $c$ is the speed of light.

### Relativistic pressure moments

Relativistic directional temperatures are defined from the pressure tensor,
whose kinetic integrand is $p_i v_j$. For a gyrotropic distribution, the
parallel and perpendicular pressures are

\[
P_{\parallel}
=\int p_{\parallel}v_{\parallel}
f(\bar{\boldsymbol{u}})\,d^3\bar{u}
=\int \frac{m u^2}{\gamma}\cos^2(\theta)
f(\bar{\boldsymbol{u}})\,d^3\bar{u}
\]

and

\[
P_{\perp}
=\frac{1}{2}\int p_{\perp}v_{\perp}
f(\bar{\boldsymbol{u}})\,d^3\bar{u}
=\frac{1}{2}\int \frac{m u^2}{\gamma}\sin^2(\theta)
f(\bar{\boldsymbol{u}})\,d^3\bar{u}.
\]

The parallel direction contains one degree of freedom. The perpendicular plane
contains two, and the factor $1/2$ defines $P_{\perp}$ per perpendicular degree
of freedom. The pressure-equivalent directional temperatures are

\[
T_{\parallel}=\frac{P_{\parallel}}{n},
\qquad
T_{\perp}=\frac{P_{\perp}}{n}.
\]

Temperatures are expressed directly in energy units, so the Boltzmann constant
is absorbed into these definitions. For a general non-Maxwellian distribution,
these are pressure-equivalent moments rather than equilibrium thermodynamic
temperatures.

### Nonrelativistic limit

When $u\ll c$, $\gamma\rightarrow1$ and

\[
E=\frac{p^2}{2m}=\frac{1}{2}mu^2.
\]

The relativistic pressure definitions then reduce to the familiar expressions

\[
T_{\parallel}
=\frac{2}{n}\int E\cos^2(\theta)
f(\bar{\boldsymbol{u}})\,d^3\bar{u}
\]

and

\[
T_{\perp}
=\frac{1}{n}\int E\sin^2(\theta)
f(\bar{\boldsymbol{u}})\,d^3\bar{u}.
\]

The different factors arise because $T_{\parallel}$ describes one degree of
freedom, whereas $T_{\perp}$ describes the common temperature of two
perpendicular degrees of freedom.

## Numerical implementation

Let $f_{ij}=f(\bar{u}_j,\theta_i)$. CQL3D directly supplies the cell weights

\[
\Delta V_{\bar{u},j}=\bar{u}_j^2\,\Delta\bar{u}_j
\]

and

\[
\Delta V_{\theta,i}=2\pi\sin(\theta_i)\,\Delta\theta_i
\]

in `f4ddv` and `f4ddt`, respectively. The density is evaluated as

\[
n=\sum_i\sum_j
f_{ij}\,\Delta V_{\bar{u},j}\,\Delta V_{\theta,i}.
\]

Define the relativistic pressure-energy factor

\[
u_j=\bar{u}_ju_{\mathrm{norm}},
\qquad
\gamma_j=\sqrt{1+\left(\frac{u_j}{c}\right)^2},
\qquad
Q_j=\frac{m u_j^2}{\gamma_j}.
\]

The discrete directional temperatures are

\[
T_{\parallel}
=\frac{1}{n}\sum_i\sum_j
Q_j\cos^2(\theta_i)
f_{ij}\,\Delta V_{\bar{u},j}\,\Delta V_{\theta,i}
\]

and

\[
T_{\perp}
=\frac{1}{2n}\sum_i\sum_j
Q_j\sin^2(\theta_i)
f_{ij}\,\Delta V_{\bar{u},j}\,\Delta V_{\theta,i}.
\]

The implementation evaluates $Q_j$ in `keV`, so density is reported in
`ions/cm^3` and both directional temperatures are reported in `keV`.

---

**Navigation:** Previous: — | [Up: Test 002](../README.md) | [Next: Run the conversion](../02_run_test/README.md)
