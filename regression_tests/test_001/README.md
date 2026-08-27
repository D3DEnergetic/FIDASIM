[Regression tests](../README.md) / Test 001

**Navigation:** Previous: — | [Up: Regression tests](../README.md) |
[Next: Test 002](../test_002/README.md)

---

# Test 001: mean-free-path regression

Test 001 checks the neutral-beam stopping mean free path produced by
`libfida::colrad` against an independent ADAS ADF21 stopping-coefficient
reference. The Python implementation provides the reference, post-processing,
and comparison functionality and is the canonical high-level implementation
of the regression.

## Implementation status

- **Python workflow:** Canonical implementation for reference generation,
  FIDASIM post-processing, scientific comparison, and reporting.
- **FIDASIM calculation:** The top-level runner executes the compiled Test 001,
  then passes the newly generated `test_001.h5` through the Python workflow.
- **Legacy MATLAB sources:** Archived under each stage's `legacy/` directory;
  they are retained for historical comparison and are not run by `run.sh`.

The canonical workflow preserves the established scientific definition and
acceptance criterion. The Fortran calculation remains unchanged.

## Objective, scope, and test boundary

The regression asks whether the FIDASIM attenuation result agrees with the
ADF21-derived reference to within the established 12% tolerance.

The FIDASIM procedure under test is:

| Module | Procedure | Responsibility exercised by this test |
| --- | --- | --- |
| `libfida` | `colrad` | Evolves six level-resolved neutral populations through each spatial step and returns the densities used to measure attenuation. |

The test covers the resulting neutral attenuation and mean free path across
the configured density, energy, and temperature grid. Atomic-table loading,
construction of the `LocalProfiles` state, grid generation, HDF5 writing, and
procedures called internally by `colrad` support the test; they are not
independent regression targets.

The complete workflow starts with the compiled Test 001 calculation and uses
Python for all reference generation, post-processing, comparison, and reports.

## Scientific definition

### Notation and units

| Symbol | Implementation value | Units | Meaning |
| --- | --- | --- | --- |
| $E_b$ | `eb` | `eV/amu` | Beam energy per atomic mass unit. |
| $n_e$ | `dene` or `ne` | `cm^-3` | Electron density. |
| $T$ | `temperature` or `T` | `eV` | Equal electron and ion temperature, $T_e=T_i$. |
| $k_s$ | stopping coefficient | `cm^3/s` | Effective ADF21 beam-stopping rate coefficient. |
| $v_b$ | beam velocity | `cm/s` | Neutral-beam speed used by the reference. |
| $n_n(x)$ | `denn` | `cm^-3` | Total neutral density along the ray. |
| $n_{n,l}(x)$ | `denn_n` | `cm^-3` | Neutral density in atomic level $l$. |
| $\lambda_s$ | `mfp` | `cm` | Stopping mean free path. |

### ADF21 reference

The reference uses the factorized ADF21 coefficient

$$
k_s(E_b,n_e,T)
=
\frac{k_T(T)}{k_{\mathrm{ref}}}
k_{EN}(E_b,n_e),
$$

where $k_T$ is linearly interpolated on the temperature table and $k_{EN}$ is
linearly interpolated on the energy-density table. Requested grid points must
remain inside the ADF21 table bounds.

The reference obtains the beam mass $A_b$ from the ADF21 filename,
sets $M_b=m_p A_b$, and evaluates

$$
v_b
=
10^2\sqrt{\frac{2eE_b}{M_b}},
$$

with $m_p=1.6726\times10^{-27}\ \mathrm{kg}$,
$e=1.6020\times10^{-19}\ \mathrm{C}$, and the factor $10^2$ converting
metres per second to centimetres per second. The reference mean free path is

$$
\lambda_s^{\mathrm{ADF21}}
=
\frac{v_b}{n_e k_s}.
$$

The Python reference stage uses these established constants, interpolations,
and operations without revising the physics definition.

### FIDASIM profile estimator

For one positive neutral-density profile $n_n(x_i)$ on a uniform grid, first
normalize the profile:

$$
z_i
=
\frac{n_n(x_i)}{\max_j n_n(x_j)}.
$$

Using the first grid spacing $\Delta x=x_1-x_0$, form every local natural-log
slope and its negative reciprocal:

$$
s_i
=
\frac{\ln z_{i+1}-\ln z_i}{\Delta x},
\qquad
\ell_i=-\frac{1}{s_i},
\qquad
i=0,\ldots,n_x-2.
$$

The reported FIDASIM mean free path is the ordinary arithmetic mean

$$
\lambda_s^{\mathrm{FIDASIM}}
=
\frac{1}{n_x-1}
\sum_{i=0}^{n_x-2}\ell_i.
$$

The canonical input has $n_x=100$, so each profile contributes exactly 99
local reciprocal slopes. The algorithm uses `log`, not `log10`, and an
ordinary mean, not a global log-linear fit, filtered fit, weighted mean, or
asymptotic fit.

Before applying the estimator, the Python postprocessor requires:

- expected array ranks and compatible grid extents;
- finite, strictly positive `denn` and `denn_n` values;
- finite, strictly increasing, uniform `x` columns; and
- `sum(denn_n, axis=level)` to reproduce `denn` within relative tolerance
  `1e-12`.

### Scientific comparison and acceptance tolerance

For every density, energy, and temperature point, Stage 3 calculates

$$
\epsilon
=
\frac{
\left|\lambda_s^{\mathrm{FIDASIM}}
-\lambda_s^{\mathrm{ADF21}}\right|
}{
\left|\lambda_s^{\mathrm{ADF21}}\right|
}.
$$

With $N=n_{den}n_{eb}n_t$, the report contains

$$
\epsilon_{\max}=\max_i\epsilon_i,
\qquad
\overline{\epsilon}=\frac{1}{N}\sum_i\epsilon_i,
$$

and the sample standard deviation

$$
s_{\epsilon}
=
\sqrt{
\frac{1}{N-1}
\sum_i(\epsilon_i-\overline{\epsilon})^2
}.
$$

The test passes when $\epsilon_{\max}\leq0.12$. This 12% criterion measures
scientific agreement between FIDASIM and the independent ADF21 reference.

## Canonical test conditions

[`01_reference/input.nml`](01_reference/input.nml) supplies the shared grid to
the reference and compiled Test 001 calculation.

| Input | Current value | Units and interpretation |
| --- | --- | --- |
| `nden` | `3` | Number of electron-density cases. |
| `dene` | `1e12, 1e13, 1e14` | `cm^-3`; explicit density vector. |
| `Lx` | `8e3, 8e2, 6e1` | `cm`; ray length for each density case. |
| `nx` | `100` | Spatial points per density case. |
| `tmin`, `tmax`, `nt` | `1e3`, `19e3`, `10` | `eV`; inclusive linear temperature grid. |
| `emin`, `emax`, `neb` | `5e3`, `1e5`, `5` | `eV/amu`; inclusive linear energy grid. |
| `beam_mass` | `1.0` | `amu`; FIDASIM and output metadata mass. |
| `ion_mass` | `1` | `amu`; thermal-ion mass used by FIDASIM. |

The ADF21 reference is
[`bms93#h_h1.dat`](01_reference/adas/bms93%23h_h1.dat). The Python reference
parser recognizes the beam species from the filename and the current table
describes hydrogen on hydrogen.

## Python workflow

```text
input.nml + bms93#h_h1.dat       compiled Test 001
              |                         |
              v                         v
      Python Stage 1               test_001.h5
        reference.h5                    |
              |                         v
              |                  Python Stage 2
              |                  postprocessed.h5
              |                         |
              +------------+------------+
                           |
                           v
                    Python Stage 3
                  figures and report
```

| Stage | Implementation | Output |
| --- | --- | --- |
| 1. ADF21 reference | [`produce_reference.py`](01_reference/produce_reference.py) and support modules | `01_reference/output_data/reference.h5` and three figure pairs |
| 2. FIDASIM result and post-processing | [`test_001.f90`](02_run_test/test_001.f90) produces raw profiles; [`postprocess.py`](02_run_test/postprocess.py) derives mean free paths | `02_run_test/output_data/postprocessed.h5` and four figure pairs |
| 3. Scientific comparison | [`compare_mean_free_paths.py`](03_compare/compare_mean_free_paths.py) | `03_compare/output_data/mfp_regression_report.txt` and seven figure pairs |

All Python-stage paths are configured independently. Relative paths are
resolved from the directory containing the corresponding namelist, not from
the process working directory.

## Prerequisites and fresh-clone setup

### Python environment

The workflow uses Python 3.12 with `numpy`, `scipy`, `h5py`, `matplotlib`, and
`f90nml`. The repository environment is defined in
[`regression_tests/environment.yml`](../environment.yml).

Create it once if it does not already exist, then activate it in each new
shell:

```bash
cd "$FIDASIM_DIR/regression_tests"
conda env create -f environment.yml
conda activate FIDASIM_env
```

The run scripts call `python3` directly and therefore require
`FIDASIM_env` to be active. The required Python packages are:

- `numpy`
- `scipy`
- `h5py`
- `matplotlib`
- `f90nml`

### Generated artifacts

HDF5 and generated-output files are intentionally ignored by Git. The complete
runner overwrites `02_run_test/test_001.h5` with the output of the current
compiled Test 001 calculation. The shared build procedure is documented in the
[regression-test build guide](../docs/build_system.md).

## Run the complete workflow

After building FIDASIM and activating the Python environment, run the compiled
Test 001 calculation and all three Python stages from the Test 001 directory:

```bash
conda activate FIDASIM_env
cd "$FIDASIM_DIR/regression_tests/test_001"
./run.sh
```

[`run.sh`](run.sh) uses simple paths relative to the Test 001 directory. It
runs the regression executable from `02_run_test`, and then runs the three
Python stages in dependency order. It returns the Stage 3 scientific regression
status. Execution, reference, or post-processing errors are immediate failures.

### Run or debug one stage

From `regression_tests/test_001`, run an individual Python entry point with its
configuration:

```bash
python3 01_reference/produce_reference.py 01_reference/input_config.nml
python3 02_run_test/postprocess.py 02_run_test/input_config.nml
python3 03_compare/compare_mean_free_paths.py 03_compare/input_config.nml
```

Stages 1–3 may be rerun independently once their input artifacts exist.

### Legacy MATLAB sources

The former MATLAB scripts are preserved under `01_reference/legacy/`,
`02_run_test/legacy/`, and `03_compare/legacy/`. They retain the source as it
existed at the point of migration and are not maintained as runnable parallel
workflows. The Python stages and top-level `run.sh` are authoritative.

## Configuration interface

The checked-in namelists are complete, valid examples. Each stage validates
required values before performing numerical work. Every relative path is
resolved from the configuration file that contains it.

### Stage 1: `&reference`

Configuration:
[`01_reference/input_config.nml`](01_reference/input_config.nml)

| Field | Type | Required | Meaning and validation |
| --- | --- | --- | --- |
| `grid_input_file` | Path string | Yes | Existing shared `&run_input` namelist. |
| `adf21_file` | Path string | Yes | Existing ADF21 stopping-coefficient table. |
| `output_file` | Path string | Yes | Python-owned reference HDF5 destination; its parent is created. |
| `figure_directory` | Path string | Yes | Destination for reference PNG/PDF pairs; it is created. |

The stage requires positive grid sizes, requires exactly `nden` density
values, and rejects requested density, energy, or temperature points outside
the ADF21 table.

### Stage 2: `&postprocess`

Configuration:
[`02_run_test/input_config.nml`](02_run_test/input_config.nml)

| Field | Type | Required | Meaning and validation |
| --- | --- | --- | --- |
| `source_h5` | Path string | Yes | Existing raw `test_001.h5`, opened read-only. |
| `output_h5` | Path string | Yes | Separate derived HDF5 destination; it must not resolve to `source_h5`. |
| `figures_dir` | Path string | Yes | Destination for four PNG/PDF figure pairs. |

The namelist must contain exactly one `&postprocess` block and no unknown
fields. The raw source is opened read-only, and the derived output must have a
different path.

### Stage 3: `&compare`

Configuration:
[`03_compare/input_config.nml`](03_compare/input_config.nml)

| Field | Type | Required | Meaning and validation |
| --- | --- | --- | --- |
| `reference_file` | Path string | Yes | Existing Stage 1 Python reference HDF5. |
| `postprocessed_file` | Path string | Yes | Existing Stage 2 derived HDF5. |
| `raw_test_file` | Path string | Yes | Existing raw FIDASIM file used only for maximum-error profile diagnostics. |
| `output_directory` | Path string | Yes | Destination for the report and optional `figures/` directory. |
| `generate_plots` | Logical scalar | Yes | `.true.` writes seven figure pairs; `.false.` writes only the report. |
| `mfp_relative_error_tolerance` | Positive real | Yes | Scientific maximum-error tolerance; the checked-in value is `0.12`. |

Density, energy, and beam-mass grids must match exactly. Temperature grids
must agree within relative tolerance `1e-7`, and both mean-free-path arrays
must have semantic shape `(nden, neb, nt)`.

## HDF5 and array-order contract

The Fortran-written logical array dimensions appear in the opposite order in
the HDF5 file as viewed through h5py. Python transposes immediately after
reading so all numerical work uses semantic axes. See the shared
[HDF5 array-ordering guide](../docs/hdf5_array_ordering.md) for the general
rule.

| Quantity | Semantic order used by Python | Raw h5py order | Canonical shapes |
| --- | --- | --- | --- |
| `x` | `(nx, nden)` | `(nden, nx)` | `(100, 3)` logical; `(3, 100)` raw |
| `denn` | `(nx, nden, neb, nt)` | `(nt, neb, nden, nx)` | `(100, 3, 5, 10)` logical; `(10, 5, 3, 100)` raw |
| `denn_n` | `(nlev, nx, nden, neb, nt)` | `(nt, neb, nden, nx, nlev)` | `(6, 100, 3, 5, 10)` logical; `(10, 5, 3, 100, 6)` raw |
| `mfp` | `(nden, neb, nt)` | `(nt, neb, nden)` | `(3, 5, 10)` logical; `(10, 5, 3)` raw |

The explicit h5py-to-semantic transposes are:

```python
x = raw_x.transpose(1, 0)
denn = raw_denn.transpose(3, 2, 1, 0)
denn_n = raw_denn_n.transpose(4, 3, 2, 1, 0)
mfp = raw_mfp.transpose(2, 1, 0)
```

Python writes `/mfp` in raw `(nt, neb, nden)` order and records the semantic
axis order `(nden, neb, nt)` in the dataset metadata.

### Python HDF5 artifacts

| Artifact | Required datasets | Contract |
| --- | --- | --- |
| Stage 1 `reference.h5` | `/ne`, `/eb`, `/T`, `/ab`, `/mfp` | ADF21 reference. Vectors retain the established row-vector file layout; `/mfp` has raw shape `(nt, neb, nden)`. |
| Stage 2 `postprocessed.h5` | `/dene`, `/eb`, `/temperature`, `/beam_mass`, `/mfp` | Self-contained derived mean-free-path result. `/mfp` has raw shape `(nt, neb, nden)` and records both semantic and h5py axis metadata. |

Stage 2 also records its estimator, slope count, creation time, configured
source path, and read-only source mode as root attributes. It never appends
data to the raw HDF5 file.

## Generated outputs and figures

Generated HDF5, PNG, PDF, report, and `output_data/` artifacts are ignored by
Git. They remain local review evidence and are regenerated by their owning
stage.

| Stage | Output | Contents |
| --- | --- | --- |
| 1 | `01_reference/output_data/reference.h5` | ADF21 grids, beam mass, and reference MFP. |
| 1 | `01_reference/output_data/figures/bms93_h_h1_MFP_{1..3}.{png,pdf}` | Reference MFP versus temperature, one density per figure. |
| 2 | `02_run_test/output_data/postprocessed.h5` | FIDASIM grids, beam mass, and derived MFP. |
| 2 | `02_run_test/output_data/figures/FIDASIM_MFP_{1..3}.{png,pdf}` | FIDASIM MFP versus temperature, one density per figure. |
| 2 | `02_run_test/output_data/figures/Example_neutral_density_profile.{png,pdf}` | `log10` of the first normalized FIDASIM neutral-density profile. |
| 3 | `03_compare/output_data/mfp_regression_report.txt` | PASS/FAIL, maximum, mean, sample standard deviation, and 12% acceptance check. |
| 3 | `03_compare/output_data/figures/comparison__MFP_{1..3}.{png,pdf}` | Solid ADF21 reference and dashed FIDASIM MFP curves. |
| 3 | `03_compare/output_data/figures/comparison__relative_error_{1..3}.{png,pdf}` | Relative error versus temperature with the vertical range fixed to the tolerance. |
| 3 | `03_compare/output_data/figures/max_relative_error_profiles.{png,pdf}` | ADF21 exponential and normalized raw FIDASIM profile at the maximum-error condition. |

Each scientific figure is written as a 600-DPI PNG and a vector PDF.

## Pass, failure, and reporting behavior

- Stage 1 fails on missing files, malformed ADF21 records, inconsistent grid
  sizes, or requests outside the table domain.
- Stage 2 prints `PASS` only after shape, positivity, grid, level-sum,
  estimator, and output checks complete.
- Stage 3 always writes the text report and, when enabled, all diagnostic
  figures before returning its scientific PASS/FAIL status.
- Stage 3 returns 0 when the finite maximum relative error is no greater than
  the configured tolerance and 1 otherwise.

Do not infer success only from the presence of figures. Review the Stage 3
report for scientific status.

## Debugging guide

### Import or interpreter failure

Confirm that `FIDASIM_env` is active and provides the required packages:

```bash
conda activate FIDASIM_env
command -v python3
python3 -c "import f90nml, h5py, matplotlib, numpy, scipy"
```

If an import fails, recreate or update `FIDASIM_env` from
[`regression_tests/environment.yml`](../environment.yml).

### Missing executable

If `02_run_test/test_001` is missing, build FIDASIM and its regression targets
before running the workflow. The normal runner generates `test_001.h5`; it
does not require a preserved raw result.

### Configuration path failure

Read the path from the namelist that reports the error and resolve it relative
to that namelist's directory. Running the same stage from another working
directory does not change configured-path semantics.

### Shape or axis failure

Compare the h5py-visible shape with the
[array-order table](#hdf5-and-array-order-contract). Do not remove or guess a
transpose to make two shapes superficially match; first identify each
physical axis.

### Stage 3 scientific failure

Open `mfp_regression_report.txt`, then inspect the relative-error figures and
`max_relative_error_profiles`. The report and plots persist even though the
stage returns nonzero. Do not relax the 12% threshold merely to make the test
pass.

## Limitations and deferred work

- The regression covers the canonical hydrogen-on-hydrogen grid produced by
  the current Test 001 calculation; it is not a parameter-space survey.
- The ADF21 reference assumes $T_e=T_i$ and linear interpolation inside the
  table domain.
- The legacy profile estimator averages local reciprocal slopes directly. It
  does not identify an entrance transient, reject noisy tail points, or fit an
  asymptotic attenuation interval.
An improved transient detector and asymptotic attenuation fit are explicitly
deferred to a separate scientific-method milestone. That change requires its
own theory, configuration, validation cases, tolerance review, and reference
baseline and must remain separate from the language port.

## Directory map

```text
test_001/
├── README.md
├── makefile
├── run.sh
├── 01_reference/
│   ├── input.nml
│   ├── input_config.nml
│   ├── adas/
│   ├── legacy/                       # archived MATLAB sources
│   ├── produce_reference.py
│   ├── adf21.py
│   ├── reference_workflow.py
│   └── output_data/                  # ignored generated outputs
├── 02_run_test/
│   ├── test_001.f90
│   ├── test_001.h5                   # ignored generated raw result
│   ├── legacy/                       # archived MATLAB source
│   ├── input_config.nml
│   ├── postprocess.py
│   ├── postprocess_support.py
│   └── output_data/                  # ignored generated outputs
└── 03_compare/
    ├── legacy/                       # archived MATLAB source
    ├── input_config.nml
    ├── compare_mean_free_paths.py
    ├── comparison_config.py
    ├── comparison_workflow.py
    └── output_data/                  # ignored report and figures
```

---

**Navigation:** Previous: — | [Up: Regression tests](../README.md) |
[Next: Test 002](../test_002/README.md)
