[Regression tests](../../README.md) / [Test 002](../README.md) / Run the conversion

**Navigation:** [Previous: Reference data](../01_reference/README.md) | [Up: Test 002](../README.md) | [Next: Compare the moments](../03_compare/README.md)

---

# Convert CQL3D distributions to FIDASIM coordinates

Stage 2 converts each reference distribution from normalized proper velocity
and pitch angle, $f_{\mathrm{CQL}}(\bar{u},\theta)$, to the FIDASIM
energy-pitch density, $F(E,P)$. The transformation must preserve the number
of particles represented by each phase-space region.

## Relativistic coordinate transformation

### 1. Source density element

The CQL3D distribution is expressed in normalized proper velocity
$\bar{u}$, pitch angle $\theta$, and gyrophase $\phi$. Its normalized
proper-velocity volume element is

\[
d^3\bar{u}
=\bar{u}^2\sin\theta\,
d\bar{u}\,d\theta\,d\phi.
\]

The distribution is gyrotropic, so integration over $0\leq\phi<2\pi$
gives the particle-density element

\[
dn
=2\pi f_{\mathrm{CQL}}(\bar{u},\theta)
\bar{u}^2\sin\theta\,
d\bar{u}\,d\theta.
\]

### 2. Target density element

FIDASIM uses kinetic energy $E$ and pitch $P$. Define $F(E,P)$ by

\[
dn=F(E,P)\,dE\,dP.
\]

Equating the source and target descriptions gives

\[
F(E,P)\,dE\,dP
=2\pi f_{\mathrm{CQL}}(\bar{u},\theta)
\bar{u}^2\sin\theta\,
d\bar{u}\,d\theta.
\]

### 3. Energy and pitch coordinates

The dimensional proper velocity is

\[
u=\bar{u}u_{\mathrm{norm}}.
\]

Define the rest-mass energy in keV and the normalized proper-velocity scale as

\[
M=mc^2,
\qquad
\beta_{\mathrm{norm}}=\frac{u_{\mathrm{norm}}}{c}.
\]

The Lorentz factor can then be written as

\[
\gamma
=\sqrt{1+\left(\frac{u}{c}\right)^2}
=\sqrt{1+\beta_{\mathrm{norm}}^2\bar{u}^2}.
\]

The target coordinates are

\[
E(\bar{u})=M(\gamma-1)
\]

and

\[
P(\theta)=\cos\theta.
\]

### 4. Energy derivative

Differentiate the relativistic energy with respect to $\bar{u}$:

\[
\frac{dE}{d\bar{u}}
=M\frac{d\gamma}{d\bar{u}}.
\]

Since

\[
\gamma=(1+\beta_{\mathrm{norm}}^2\bar{u}^2)^{1/2},
\]

its derivative is

\[
\frac{d\gamma}{d\bar{u}}
=\frac{\beta_{\mathrm{norm}}^2\bar{u}}{\gamma}.
\]

Therefore,

\[
\boxed{
\frac{dE}{d\bar{u}}
=M\frac{\beta_{\mathrm{norm}}^2\bar{u}}{\gamma}
}.
\]

### 5. Pitch derivative and sign

Differentiating $P=\cos\theta$ gives

\[
\boxed{
\frac{dP}{d\theta}=-\sin\theta
}.
\]

The negative sign records the reversal of coordinate direction:

\[
\theta:0\rightarrow\pi
\qquad\Longleftrightarrow\qquad
P:1\rightarrow-1.
\]

For example,

\[
\int_0^\pi g(\theta)\sin\theta\,d\theta
=\int_1^{-1}g(P)(-dP)
=\int_{-1}^{1}g(P)\,dP.
\]

Thus, the negative derivative reverses the integration limits; it does not
make the transformed distribution negative.

### 6. Jacobian determinant

The two-dimensional Jacobian matrix is

\[
\frac{\partial(E,P)}{\partial(\bar{u},\theta)}
=
\begin{pmatrix}
\dfrac{M\beta_{\mathrm{norm}}^2\bar{u}}{\gamma} & 0
\\[6pt]
0 & -\sin\theta
\end{pmatrix}.
\]

Its determinant is

\[
\det\left(
\frac{\partial(E,P)}{\partial(\bar{u},\theta)}
\right)
=-
\frac{M\beta_{\mathrm{norm}}^2\bar{u}}{\gamma}
\sin\theta.
\]

A phase-space cell has positive area, so the differential volume uses the
absolute value of this determinant:

\[
J
=\left|
\det\left(
\frac{\partial(E,P)}{\partial(\bar{u},\theta)}
\right)
\right|
=\frac{M\beta_{\mathrm{norm}}^2\bar{u}}{\gamma}
\sin\theta.
\]

Consequently,

\[
dE\,dP=J\,d\bar{u}\,d\theta.
\]

### 7. Transformed distribution

Substitute the Jacobian into particle-number conservation:

\[
F(E,P)
\frac{M\beta_{\mathrm{norm}}^2\bar{u}}{\gamma}
\sin\theta
=2\pi f_{\mathrm{CQL}}(\bar{u},\theta)
\bar{u}^2\sin\theta.
\]

Cancel $\sin\theta$ and one power of $\bar{u}$:

\[
\boxed{
F(E,P)
=\frac{2\pi\gamma\bar{u}}
{M\beta_{\mathrm{norm}}^2}
f_{\mathrm{CQL}}(\bar{u},\theta)
}.
\]

The expression is evaluated at corresponding coordinates

\[
E=M(\gamma-1),
\qquad
P=\cos\theta.
\]

Because $M$ is expressed in keV and pitch is dimensionless, the resulting
FIDASIM distribution has units

\[
[F]
=\frac{\text{ions}}
{\mathrm{cm}^3\,\mathrm{keV}\,dP}.
\]

## Nonrelativistic approximation

The nonrelativistic limit provides an independent check of the relativistic
derivation. When $u\ll c$, define

\[
x=\beta_{\mathrm{norm}}^2\bar{u}^2\ll1.
\]

Expanding the Lorentz factor gives

\[
\gamma=\sqrt{1+x}\simeq1+\frac{x}{2}.
\]

The relativistic kinetic energy therefore reduces to

\[
E=M(\gamma-1)
\simeq\frac{1}{2}M\beta_{\mathrm{norm}}^2\bar{u}^2.
\]

Thus, the nonrelativistic energy is

\[
\boxed{
E_{\mathrm{NR}}(\bar{u})
=\frac{1}{2}M\beta_{\mathrm{norm}}^2\bar{u}^2
}
\]

and its derivative is

\[
\boxed{
\frac{dE_{\mathrm{NR}}}{d\bar{u}}
=M\beta_{\mathrm{norm}}^2\bar{u}
}.
\]

The pitch transformation is unchanged, so the nonrelativistic Jacobian
magnitude is

\[
J_{\mathrm{NR}}
=M\beta_{\mathrm{norm}}^2\bar{u}\sin\theta.
\]

Substitution into the same particle-number conservation equation gives

\[
\boxed{
F_{\mathrm{NR}}(E_{\mathrm{NR}},P)
=\frac{2\pi\bar{u}}
{M\beta_{\mathrm{norm}}^2}
f_{\mathrm{CQL}}(\bar{u},\theta)
}.
\]

The relativistic and nonrelativistic results are related by

\[
F(E,P)=\gamma F_{\mathrm{NR}}(E_{\mathrm{NR}},P).
\]

As $u/c\rightarrow0$, $\gamma\rightarrow1$, so both the energy mapping and
the transformed distribution approach their nonrelativistic forms.

### Connection to the conventional $\sqrt{E}$ density integral

Proper velocity is related to ordinary velocity by $u=\gamma v$. In the
nonrelativistic limit, $\gamma\simeq1$, so

\[
u\simeq v.
\]

The physical CQL3D proper-velocity distribution can therefore be treated as
the ordinary velocity-space distribution $f_v(\mathbf{v})$ in this limit. Let
$\varepsilon$ convert the numerical energy coordinate to physical energy. For
example, $\varepsilon=e$ for energy in eV and $\varepsilon=1000e$ for energy
in keV. Then

\[
v=\sqrt{\frac{2\varepsilon E}{m}},
\qquad
v^2dv
=\sqrt{2}\left(\frac{\varepsilon}{m}\right)^{3/2}
\sqrt{E}\,dE.
\]

After integration over gyrophase, the density is therefore

\[
\boxed{
n
=2\pi\sqrt{2}
\left(\frac{\varepsilon}{m}\right)^{3/2}
\int_0^\infty\int_{-1}^{1}
f_v(E,P)\sqrt{E}\,dP\,dE
}.
\]

CQL3D stores the scaled distribution

\[
f_{\mathrm{CQL}}=u_{\mathrm{norm}}^3f_v.
\]

Using

\[
M\beta_{\mathrm{norm}}^2
=\frac{mu_{\mathrm{norm}}^2}{\varepsilon},
\]

the $u_{\mathrm{norm}}^3$ factors cancel from the CQL3D form of the
nonrelativistic integral, recovering the boxed conventional expression above.
The two formulations are therefore equivalent; they use differently normalized
distribution functions.

## Numerical implementation

The production conversion follows the nonrelativistic convention used by
FIDASIM. First, the stored CQL3D distribution is converted back to the physical
velocity-space distribution:

\[
f_v=\frac{f_{\mathrm{CQL}}}{u_{\mathrm{norm}}^3}.
\]

### Unit convention

The implementation uses cgs mechanical units together with keV because these
are the native conventions of the CQL3D and FIDASIM data:

| Quantity | Numerical unit |
|---|---|
| $u_{\mathrm{norm}}$ and $u$ | cm/s |
| Ion mass $m$ | g |
| Energy $E$ | keV |
| Density $n$ | ions/cm³ |

The ion mass is calculated from the isotope mass in atomic mass units using

\[
m=A m_{\mathrm{u}},
\qquad
m_{\mathrm{u}}=1.66053906660\times10^{-24}\ \mathrm{g}.
\]

The symbol $\varepsilon$ is only a unit-conversion factor; it is not an
additional physical term. One electron-volt is defined as

\[
1\ \mathrm{eV}
=1.602176634\times10^{-19}\ \mathrm{J}.
\]

Therefore,

\[
1\ \mathrm{keV}
=10^3\ \mathrm{eV}
=1.602176634\times10^{-16}\ \mathrm{J}.
\]

Since

\[
1\ \mathrm{J}=10^7\ \mathrm{erg},
\]

one keV expressed in erg is

\[
1\ \mathrm{keV}
=1.602176634\times10^{-9}\ \mathrm{erg}.
\]

We collect this conversion into

\[
\varepsilon
=1.602176634\times10^{-9}\ \mathrm{erg/keV},
\]

which is named `ERG_PER_KEV` in the code. Multiplying a value in keV by
$\varepsilon$ converts it to erg; dividing a value in erg by $\varepsilon$
converts it to keV.

The nonrelativistic kinetic energy calculated from mass in grams and velocity
in cm/s is

\[
K_{\mathrm{erg}}=\frac{1}{2}mu^2.
\]

This result is in erg because

\[
1\ \mathrm{erg}=1\ \mathrm{g\,cm^2/s^2},
\]

so the numerical energy in keV is

\[
E_{\mathrm{keV}}
=\frac{K_{\mathrm{erg}}}{\varepsilon}
=\frac{mu^2}{2\varepsilon}.
\]

For example, dividing $1.602176634\times10^{-9}$ erg by
$1.602176634\times10^{-9}$ erg/keV gives exactly 1 keV.

The CQL3D distribution and the physical velocity-space distribution have units

\[
[f_{\mathrm{CQL}}]
=\frac{\mathrm{ions}\,u_{\mathrm{norm}}^3}
{\mathrm{cm}^3(\mathrm{cm/s})^3},
\qquad
[f_v]
=\frac{\mathrm{ions}}
{\mathrm{cm}^3(\mathrm{cm/s})^3}.
\]

In the conventional transformation,

\[
\left(\frac{\varepsilon}{m}\right)^{3/2}\sqrt{E}
\]

has units $(\mathrm{cm/s})^3/\mathrm{keV}$. Multiplication by $f_v$
therefore gives

\[
[F]
=\frac{\mathrm{ions}}
{\mathrm{cm}^3\,\mathrm{keV}\,dP},
\]

which is the FIDASIM energy-pitch distribution convention. Pitch $P$ is
dimensionless; $dP$ is retained in the unit label to identify the coordinate
density.

### Transformation evaluated by the code

Using this unit convention, the code evaluates

\[
E=\frac{mu^2}{2\varepsilon}
\]

and

\[
F(E,P)
=2\pi\sqrt{2}
\left(\frac{\varepsilon}{m}\right)^{3/2}
\sqrt{E}\,f_v.
\]

Pitch is calculated as $P=\cos\theta$. Because increasing $\theta$ produces
decreasing $P$, the pitch array and the corresponding distribution axis are
reversed before output so pitch increases from $-1$ to $+1$.

The implementation also calculates the corresponding relativistic energy and
distribution values as diagnostics, but these are not used as the FIDASIM
output.

The analytical transformation produces nonuniform energy and pitch coordinates
from the uniform CQL3D $\bar{u}$ and $\theta$ grids. FIDASIM sampling uses
constant energy and pitch spacings, so Stage 2 must subsequently remap the
transformed distribution onto uniform $E$ and $P$ grids while preserving
particle density as accurately as the grid resolution allows.

### Uniform cell-centered grids

FIDASIM samples uniformly within a selected cell using the stored coordinate
as its center. Energy centers must therefore be offset from zero so the first
cell cannot produce a negative sampled energy.

The CQL3D proper-velocity coordinates are uniformly spaced. The final
coordinate is the center of the final source cell, so its upper edge is

\[
\bar{u}_{\max,\mathrm{edge}}
=\bar{u}_{N-1}+\frac{\Delta\bar{u}}{2}.
\]

The upper energy boundary is obtained by applying the same nonrelativistic
energy transformation to this edge:

\[
E_{\max}
=\frac{m
\left(\bar{u}_{\max,\mathrm{edge}}u_{\mathrm{norm}}\right)^2}
{2\varepsilon}.
\]

For $N_E$ target energy cells,

\[
\Delta E=\frac{E_{\max}}{N_E},
\qquad
E_i=\left(i+\frac{1}{2}\right)\Delta E.
\]

The first energy cell then covers $[0,\Delta E]$. For $N_P$ target pitch
cells,

\[
\Delta P=\frac{2}{N_P},
\qquad
P_j=-1+\left(j+\frac{1}{2}\right)\Delta P.
\]

These pitch cells exactly cover $[-1,1]$. The nonuniform transformed
distribution is linearly interpolated onto the target cell centers. The output
array is stored in FIDASIM order, `f_array(nenergy, npitch)`.

---

**Navigation:** [Previous: Reference data](../01_reference/README.md) | [Up: Test 002](../README.md) | [Next: Compare the moments](../03_compare/README.md)
