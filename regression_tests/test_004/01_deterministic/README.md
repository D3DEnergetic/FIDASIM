[Regression tests](../../README.md) / [Test 004](../README.md) / Deterministic calculation

**Navigation:** Previous: — | [Up: Test 004](../README.md) | [Next: Monte Carlo calculation](../02_monte_carlo/README.md)

---

# Deterministic ion-sink calculation

The deterministic implementation reads every smooth FIDASIM energy-pitch
distribution identified by a Test 002 Stage 2 configuration. It calculates
the ion-sink distribution and total reaction rate by deterministic quadrature
over energy, pitch, and gyrophase.

The resulting HDF5 files and principal plots are generated artifacts under the
common Test 004 output directory. They are not committed to the repository.

## Velocity construction

Assume that the magnetic field is directed along the Cartesian `+z` direction,

$$
\boldsymbol B=B\widehat{\boldsymbol z},
\qquad B>0.
$$

The distribution pitch is therefore

$$
p=\frac{v_\parallel}{v}
=\frac{v_z}{v}.
$$

For an ion with kinetic energy $E$, mass $M_i$, and charge $q_i$, define

$$
s_q=\frac{q_i}{|q_i|}.
$$

The ion speed and its parallel and perpendicular components are

$$
v=\sqrt{\frac{2E}{M_i}},
\qquad
v_\parallel=pv,
\qquad
v_\perp=v\sqrt{1-p^2}.
$$

The gyrophase $\phi\in[0,2\pi)$ is treated as a positive geometric angle.
The direction of gyromotion is included explicitly through the ion charge
sign. The ion velocity is therefore

$$
\boldsymbol v_i(E,p,\phi)=
\begin{bmatrix}
v_\perp\cos\phi\\
-s_qv_\perp\sin\phi\\
v_\parallel
\end{bmatrix}.
$$

For a positive ion, $s_q=+1$, so an ion initially moving along `+x` turns
toward `-y`. For a negative ion, $s_q=-1$, and the direction of gyromotion
is reversed.

### General neutral velocity

First consider an arbitrary neutral velocity with all three Cartesian
components:

$$
\boldsymbol v_n=
\begin{bmatrix}
v_{n,x}\\
v_{n,y}\\
v_{n,z}
\end{bmatrix}.
$$

The relative velocity vector is

$$
\boldsymbol v_{\rm rel}
=\boldsymbol v_i-\boldsymbol v_n,
$$

or explicitly,

$$
\boldsymbol v_{\rm rel}
=
\begin{bmatrix}
v_\perp\cos\phi-v_{n,x}\\
-s_qv_\perp\sin\phi-v_{n,y}\\
v_\parallel-v_{n,z}
\end{bmatrix}.
$$

The corresponding relative speed is

$$
v_{\rm rel}=
\sqrt{
\left(v_\perp\cos\phi-v_{n,x}\right)^2+
\left(-s_qv_\perp\sin\phi-v_{n,y}\right)^2+
\left(v_\parallel-v_{n,z}\right)^2
}.
$$

Expanding the squared relative speed gives

$$
\begin{aligned}
v_{\rm rel}^2
={}&
v_\perp^2\cos^2\phi
-2v_\perp v_{n,x}\cos\phi
+v_{n,x}^2\\
&{}+
v_\perp^2\sin^2\phi
+2s_qv_\perp v_{n,y}\sin\phi
+v_{n,y}^2\\
&{}+
v_\parallel^2
-2v_\parallel v_{n,z}
+v_{n,z}^2.
\end{aligned}
$$

Using

$$
\cos^2\phi+\sin^2\phi=1,
\qquad
v_\perp^2+v_\parallel^2=v^2,
$$

and

$$
v_n^2=v_{n,x}^2+v_{n,y}^2+v_{n,z}^2,
$$

the result becomes

$$
\boxed{
v_{\rm rel}^2=
v^2+v_n^2
-2v_\perp v_{n,x}\cos\phi
+2s_qv_\perp v_{n,y}\sin\phi
-2v_\parallel v_{n,z}.
}
$$

The term

$$
2s_qv_\perp v_{n,y}\sin\phi
$$

shows explicitly where the ion charge sign enters the relative speed at a
particular gyrophase.

### Neutral injection in the `x-z` plane

For the present geometry, the neutral is injected at a signed angle $\theta$
measured from `+z` toward `+x`. For neutral kinetic energy $E_n$ and mass
$M_n$,

$$
v_n=\sqrt{\frac{2E_n}{M_n}},
$$

and the neutral velocity is

$$
\boldsymbol v_n=
\begin{bmatrix}
v_n\sin\theta\\
0\\
v_n\cos\theta
\end{bmatrix}.
$$

Thus,

$$
v_{n,x}=v_n\sin\theta,
\qquad
v_{n,y}=0,
\qquad
v_{n,z}=v_n\cos\theta.
$$

Substitution into the general result gives

$$
v_{\rm rel}^2=
v^2+v_n^2
-2v_\perp v_n\sin\theta\cos\phi
-2v_\parallel v_n\cos\theta.
$$

Therefore,

$$
\boxed{
v_{\rm rel}=
\sqrt{
v^2+v_n^2
-2v_n\left(
v_\perp\sin\theta\cos\phi+
v_\parallel\cos\theta
\right)
}.
}
$$

Equivalently, before expansion,

$$
v_{\rm rel}=
\sqrt{
\left(v_\perp\cos\phi-v_n\sin\theta\right)^2+
v_\perp^2\sin^2\phi+
\left(v_\parallel-v_n\cos\theta\right)^2
}.
$$

Because the neutral has no `y` velocity component, $v_{n,y}=0$, the
explicitly charge-dependent term vanishes. Consequently, the relative speed
at a given sampled gyrophase is independent of the charge sign for this
particular neutral injection geometry.

The charge sign should nevertheless be retained in the ion velocity vector so
that the velocity construction represents the physically correct direction
of gyromotion and remains valid for more general neutral velocity vectors.

The implementation converts `keV` and isotope mass in `amu` to `cm/s` before
forming these vectors. The deterministic gyro average evaluates this
expression on the configured midpoint gyrophase grid.

## Atomic-rate evaluation

The atomic table is written through the Fortran HDF5 interface in the logical
order `cx(initial level, final level, relative energy)`. Consequently, h5py
exposes `/cross/H_H/cx` in the reversed storage order
`(relative energy, final level, initial level)`. The deterministic reader
swaps the two level axes once and presents the calculation with the canonical
order `(relative energy, initial level, final level)`, matching FIDASIM's
internal `bb_cx_rates` matrix.

For each $E,p,\phi$, the deterministic calculation follows `bb_cx_rates`: it
linearly interpolates $\log_{10}\sigma_{m\leftarrow l}$ on the table's
uniformly spaced $\log_{10}\varepsilon_{\rm rel}$ grid, clamps an out-of-range
energy to the nearest endpoint, restores the cross section in `cm^2`, and
evaluates

$$
r_m(E,p,\phi)
=v_{\rm rel}\sum_{l=1}^{6}\sigma_{m\leftarrow l}
(\varepsilon_{\rm rel})n_l.
$$

Table entries below the smallest positive tabulated cross section are restored
as zero, matching `bb_cx_rates`. The scalar kernel used by the ion sink is

$$
K(E,p)=\frac{1}{N_\phi}
\sum_{j=1}^{N_\phi}\sum_{m=1}^{6}r_m(E,p,\phi_j),
\qquad
\phi_j=\frac{2\pi(j-\tfrac12)}{N_\phi}.
$$

This is the deterministic counterpart of `get_total_cx_rate` for one neutral
type. In the Monte Carlo calculation, identical reservoir-marker velocities
make the reservoir-weighted average equal to this single-velocity rate,
independently of `reservoir_size`.

## Input configuration

The deterministic implementation reads the unified
[`../input_config_A.nml`](../input_config_A.nml). The common `test_case` and
`neutrals` blocks are shared unchanged with the Monte Carlo implementation.
The `deterministic` block contains controls specific to this calculation.
`plot_data_block` supplies common styling for both implementation plots, and
`save_data_block` defines their common output root.

The [top-level Test 004 interface](../README.md#shared-input-configuration)
documents every shared block, including units and path resolution.

### `deterministic` schema

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `comment` | String scalar | — | No | Empty | Human-readable description of the deterministic implementation. |
| `n_gyro` | Integer scalar | — | Yes | None | Number of midpoint gyrophase points; must be positive. |
| `plot_data` | Logical scalar | — | Yes | None | Enables diagnostic plots. |
| `save_data` | Logical scalar | — | Yes | None | Enables deterministic HDF5 output. |

The deterministic implementation contains no random sampling, so its block
has no `seed` field. Random-number control belongs to `monte_carlo`.

One indexed HDF5 file and, when enabled, one same-basename PNG are generated
for each converted Test 002 distribution. Their basenames follow
`fidasim_f4d_NNN_ion_sink`.

The principal plot marks the injected neutral with a filled green circle at

$$
(p_n,E_n)=\left(\frac{v_{n,z}}{|\boldsymbol v_n|},E_n\right)
=\left(\cos\theta,E_n\right).
$$

This places the neutral velocity on the same energy-pitch axes as the
reaction-weighted ion-sink distribution.

## Running the calculation

The top-level README explains
[how Test 004 finds and validates its shared input distributions](../README.md#how-test-004-finds-and-validates-its-input-distributions).
After those converted distributions are available, run:

```bash
cd regression_tests/test_004/01_deterministic
./run.sh ../input_config_A.nml
```

## Output schema

Each indexed deterministic HDF5 file contains the following datasets:

| Dataset | Shape | Units | Description |
| --- | --- | --- | --- |
| `energy` | `(n_energy)` | `keV` | Ion energy grid. |
| `pitch` | `(n_pitch)` | dimensionless | Ion pitch grid. |
| `f_array` | `(n_energy,n_pitch)` | `ions/(cm^3*keV*dP)` | Smooth Test 002 distribution. |
| `denf` | scalar | `ions/cm^3` | Authoritative ion density from Test 002. |
| `sink_distribution` | `(n_energy,n_pitch)` | `ions/(cm^3*s*keV*dP)` | Reaction-weighted ion-sink distribution. |
| `cx_kernel` | `(n_energy,n_pitch)` | `1/s` | Gyrophase-averaged CX rate kernel. |
| `energy_marginal` | `(n_energy)` | `ions/(cm^3*s*keV)` | Sink distribution integrated over pitch. |
| `pitch_marginal` | `(n_pitch)` | `ions/(cm^3*s*dP)` | Sink distribution integrated over energy. |
| `total_reaction_rate` | scalar | `ions/(cm^3*s)` | Total volumetric ion-sink rate. |
| `gyroangle` | `(n_gyro)` | `rad` | Midpoint gyrophase grid. |
| `neutral_velocity` | `(3)` | `cm/s` | Neutral Cartesian velocity `[vx,vy,vz]`. |
| `neutral_level_density` | `(6)` | `neutrals/cm^3` | Neutral densities for levels 1 through 6. |
| `neutral_energy` | scalar | `keV` | Total neutral kinetic energy. |
| `injection_angle` | scalar | degrees | Signed angle from `+z` toward `+x`. |
| `species` | scalar string | — | Canonical isotope identifier. |
| `atomic_number` | scalar | dimensionless | Nuclear proton number. |
| `mass_number` | scalar | dimensionless | Integer isotope mass number. |
| `charge_state` | scalar | elementary charge | Ion charge state. |
| `A` | scalar | `amu` | Physical isotope mass. |
| `selected_r` | scalar | `cm` | Source radial coordinate. |
| `selected_z` | scalar | `cm` | Source axial coordinate. |

Every dataset carries a description attribute and, except for `species`, a
units attribute. Root attributes record repository-relative source
distribution, Test 002 configuration, and atomic-table paths together with
the unified Test 004 configuration, neutral-level selector, common test-case
comment, and deterministic implementation comment. The legacy `comment`
attribute is retained as an alias for the implementation comment.

The HDF5 files and principal PNG plots are regenerated test artifacts. The
committed source of the distribution data remains Test 002 Stage 1.

---

**Navigation:** Previous: — | [Up: Test 004](../README.md) | [Next: Monte Carlo calculation](../02_monte_carlo/README.md)
