[Regression tests](../README.md) / Test 004

**Navigation:** [Previous: Test 003](../test_003/README.md) | [Up: Regression tests](../README.md) | [Next: Reference calculation](01_reference/README.md)

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

## Problem statement

The energy-pitch-resolved ion-sink distribution is

$$
S(E,p)=
n_i\,g(E,p)
\sum_{l,m=1}^{6}n_l
\left\langle K_{m\leftarrow l}\right\rangle_\phi(E,p),
$$

and the total volumetric ion-sink reaction rate is

$$
R=\int S(E,p)\,dE\,dp.
$$

Here \(n_i\) is the ion density, \(g(E,p)\) is the normalized ion distribution,
and \(\left\langle K_{m\leftarrow l}\right\rangle_\phi(E,p)\) is the
state-resolved, gyrophase-averaged charge-exchange reactivity:

$$
\left\langle K_{m\leftarrow l}\right\rangle_\phi(E,p)=
\left\langle
\sigma_{m\leftarrow l}(\varepsilon_{\rm rel})v_{\rm rel}
\right\rangle_{\phi}.
$$

For any gyrophase-dependent quantity \(q(\phi)\), the gyro-averaging operator
is

$$
\left\langle q\right\rangle_\phi
=\frac{1}{2\pi}\int_0^{2\pi}q(\phi)\,d\phi.
$$

The index \(l\) identifies the initial atomic energy level of the target
neutral and \(n_l\) is its number density. The index \(m\) identifies the
atomic energy level populated by charge exchange.
\(\sigma_{m\leftarrow l}\) is directional: it is the cross section for

$$
\mathrm{H}^{+}+\mathrm{H}(l)
\longrightarrow
\mathrm{H}(m)+\mathrm{H}^{+}.
$$

For ion mass \(M_i\) and neutral mass \(M_n\), the physical centre-of-mass
collision energy is

$$
E_{\rm rel}=\frac{1}{2}\mu v_{\rm rel}^{2},
\qquad
\mu=\frac{M_iM_n}{M_i+M_n}.
$$

FIDASIM's H-H cross-section tables use the relative collision energy per
reduced-mass atomic mass unit:

$$
\varepsilon_{\rm rel}
=\frac{E_{\rm rel}}{\mu/m_u}
=\frac{1}{2}m_u v_{\rm rel}^{2},
$$

expressed in `keV/amu`, where \(m_u\) is the atomic mass constant. The relative
speed is

$$
v_{\rm rel}
=\left|\boldsymbol v_i(E,p,\phi)-\boldsymbol v_n\right|.
$$

Finally, the normalized ion distribution is

$$
g(E,p)=\frac{f(E,p)}
{\int f(E,p)\,dE\,dp},
\qquad
\int g(E,p)\,dE\,dp=1.
$$

The deterministic and Monte Carlo implementations evaluate these same
quantities for matching ion and neutral isotopes.

## FIDASIM code under test

| Module | Procedure | Role |
| --- | --- | --- |
| `libfida` | `mc_sample_ion_f4d_gc` | Samples ion position and velocity from the configured distribution. |
| `libfida` | `get_total_cx_rate` | Evaluates the reservoir-averaged, level-resolved CX rate. |
| `libfida` | `store_sinks` | Accumulates the cell-resolved sink-rate density. |
| `libfida` | `store_sink_particle` | Records each weighted ion-sink particle. |
| `libfida` | `write_sink_profile` | Writes the normal FIDASIM sink HDF5 result. |

## Workflow

| Stage | Purpose |
| --- | --- |
| [`01_reference`](01_reference/README.md) | Integrate the smooth distributions deterministically over energy, pitch, and gyrophase and store trusted reference data. |
| [`02_run_test`](02_run_test/README.md) | Run the stripped-down Monte Carlo ion-sink calculation using FIDASIM procedures. |
| [`03_compare`](03_compare/README.md) | Compare the reference and Monte Carlo distributions, marginals, and total rates. |

---

**Navigation:** [Previous: Test 003](../test_003/README.md) | [Up: Regression tests](../README.md) | [Next: Reference calculation](01_reference/README.md)
