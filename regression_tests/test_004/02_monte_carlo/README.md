[Regression tests](../../README.md) / [Test 004](../README.md) / Monte Carlo test

**Navigation:** [Previous: Deterministic calculation](../01_deterministic/README.md) | [Up: Test 004](../README.md) | [Next: Comparison](../03_compare/README.md)

---

# Monte Carlo ion-sink test

The Monte Carlo implementation samples the smooth Test 002 energy-pitch
distributions and executes the ion-sink portion of `calculate_dcx_process`
using FIDASIM source procedures. Neutral tracking, attenuation, and ion-birth
processing are excluded.

The Monte Carlo and deterministic implementations read the same unified Test
004 configuration and independently consume the same generated Test 002 Stage
2 distributions. Neither implementation requires the other implementation's
outputs.

## Running the workflow

The Monte Carlo workflow uses Python wrappers around the Fortran calculation:

1. Before Fortran runs, Python reads and validates the user-facing Test 004
   configuration, discovers the Test 002 distributions, and writes a
   flattened temporary namelist for Fortran.
2. Fortran reads that normalized namelist and performs the Monte Carlo
   calculation.
3. After Fortran runs, Python reads the output files, constructs the
   energy-pitch sink products, and produces the requested plots.

At the present development stage, the Python input wrapper, Fortran
configuration and HDF5 readers, neutral-parameter construction, and artificial
FIDASIM grid/plasma/field/FBM test setup are implemented. The configured
atomic tables are loaded once before the distribution cases are processed.
For each case, the type-1 neutral density and reservoir are populated only in
the central beam cell and then released. The program prints and verifies a
summary of each stage. One ion is sampled through `mc_sample_ion_f4d_gc` to
show the production sampling interface. The complete configured marker loop
then samples the nonthermal distribution, calculates type-1 CX rates, and
stores both the central-cell sink density and production sink particles. The
particle count, spatial support, and agreement between the density and summed
particle weights are verified. When saving is enabled, the production
`write_sink_profile` routine writes and releases each case's sink data. The
Python postprocessor then validates the complete output collection, appends
the weighted energy-pitch products under `/test_004`, and, when requested,
writes one PNG ion-sink plot per case.

After reading the normalized configuration, `configure_fidasim` translates
the shared Test 004 settings into the case-independent FIDASIM controls. This
explicit boundary selects the configured atomic tables, full-distribution
charge-exchange sampling without FLR displacement, and disables unrelated
beam-stopping and nuclear calculations. Atomic tables are then loaded before
the case loop. Isotope-dependent state remains part of each case setup.

Run the current workflow from the FIDASIM repository root:

```bash
conda activate FIDASIM_env
cd regression_tests/test_004/02_monte_carlo
./run.sh ../input_config_A.nml
```

## Input configuration

The implementation-specific portion of
[`../input_config_A.nml`](../input_config_A.nml) is:

```fortran
&monte_carlo
  comment = 'Monte Carlo ion sampling and sink storage'
  n_markers = 5e5
  reservoir_size = 50
  seed = 12345
  plot_data = .true.
  save_data = .true.
/
```

The shared `test_case`, `neutrals`, `plot_data_block`, and `save_data_block`
settings are documented in the
[top-level Test 004 interface](../README.md#shared-input-configuration).
All relative paths are resolved from the unified configuration file.

### `monte_carlo` schema

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `comment` | String scalar | — | No | Empty | Human-readable description of the Monte Carlo implementation. |
| `n_markers` | Integer-valued scalar | — | Yes | None | Number of ion markers sampled independently for each distribution; range `1` to `2^63-1`. Integral scientific notation such as `1e5` is accepted and normalized to an integer. |
| `reservoir_size` | Integer-valued scalar | — | Yes | None | Number of representative neutral particles in the central-cell reservoir; range `1` to `2^31-1`. |
| `seed` | Integer-valued scalar | — | Yes | None | Serial RNG seed in the range `1` to `2^31-1`; reset for every distribution. |
| `plot_data` | Logical scalar | — | Yes | None | Enables the weighted Monte Carlo energy-pitch sink plot; requires `save_data=.true.`. |
| `save_data` | Logical scalar | — | Yes | None | Enables the production sink HDF5 output. |

Monte Carlo files are written beneath
`save_data_block/output_directory/monte_carlo`.

## Fixed numerical test setup

The following values define the regression problem and are not configurable:

- serial execution;
- one positive ion species, with isotope and charge inherited from each
  distribution;
- a `3x3x3` Cartesian beam grid with 1 cm cells and bounds
  `[-1.5,1.5] cm` on every axis;
- center test cell `[2,2,2]` in Fortran indexing;
- identity beam-grid transformation;
- `flr=0`;
- an axisymmetric interpolation grid covering `R=0...2.5 cm` and
  `Z=-2...2 cm`, with `dR=dZ=0.5 cm` and `nphi=1`;
- a uniform plasma, mask, magnetic field, and spatially replicated FBM over
  the interpolation grid;
- `B=[0,0,1] T`, zero electric field, and zero plasma rotation;
- ion density equal to the authoritative `denf` stored in the current Test 002
  distribution;
- only the type-1 population `neut%full`, populated in the center cell.

The equilibrium mask is initialized before `make_beam_grid`, and
`make_beam_grid` is called before `init_neutral_population`; this ordering is
required because reservoir storage is allocated only for beam cells marked as
plasma.

The center population is built by calling `update_neutrals` exactly
`reservoir_size` times. Every call supplies the common neutral velocity and
the six-level density vector divided by `reservoir_size`. The accumulated
center-cell density is therefore the configured level-density vector, while
the reservoir contains identical velocity markers whose weights sum to the
configured total neutral density. No other neutral type or beam cell is
populated.

## Per-case workflow

For every smooth distribution discovered through
`test_case/input_distribution_config`:

1. Read its energy grid, pitch grid, distribution, `denf`, isotope, and charge.
2. Construct the beam grid, interpolation grid, uniform equilibrium, fields,
   and spatially replicated FBM.
3. Set the production controls required by the stripped-down path:
   `n_thermal=1`, the isotope mass, `inputs%flr=0`,
   `inputs%non_thermal_cx_sampling=1`, and `reservoir_size`. The shared atomic
   tables were loaded once before entering the case loop.
4. Initialize `neut%full`, then populate its center-cell density and reservoir
   through `update_neutrals` as described above.
5. Allocate `sink%part(n_markers)` and
   `sink%dens(1,3,3,3)`, set `sink%cnt=1`, and zero the density.
6. Reset the serial FIDASIM RNG to the configured seed. This occurs after
   reservoir construction so its internal random draws cannot alter the ion
   sample sequence.
7. For each marker:
   - call `mc_sample_ion_f4d_gc` in cell `[2,2,2]`;
   - require a nonthermal sample;
   - require `ind_p=[2,2,2]` and `rp=rgc` because `flr=0`;
   - call `get_total_cx_rate` for `[nbif_type]`;
   - calculate the per-marker volumetric sink weight
     `denf*sum(rates)/n_markers`;
   - call `store_sinks`;
   - call `store_sink_particle`, passing the same volumetric weight and the
     sampled `denf4d` values expected by its production interface.
8. Require `sink%cnt-1=n_markers`, save the center-cell sink-rate density, and
   verify all other sink-density cells are zero.
9. If `save_data=.true.`, call `write_sink_profile`, then reopen the sink file
   and append the Test 004 metadata described below. Otherwise, reset
   `sink%cnt` and explicitly deallocate `sink%part` and `sink%dens` without
   calling the writer.
10. Call `free_neutral_population(neut%full)` and release the remaining
    case-specific equilibrium, FBM, interpolation-grid, and beam-grid
    allocations before loading the
    next distribution.

`write_sink_profile` resets `sink%cnt` to one and deallocates both
`sink%part` and `sink%dens`. The next case must allocate them again. The Monte
Carlo implementation must not manually deallocate those two arrays after the
writer returns. Explicit cleanup is used only on the no-save path, where the
writer is never called.

The workflow does not call `track_to_wall`, COLRAD, `store_neutrals`,
`store_births`, `store_birth_particle`, or
`write_particle_tracks_to_file`.

## Output contract

Set `inputs%result_dir` to the configured output directory and set
`inputs%runid` to the source distribution stem followed by `_ion`. For
example, the run ID for `fidasim_f4d_001.h5` is
`fidasim_f4d_001_ion`. Because the production writer appends `_sink.h5`, it
then creates:

```text
fidasim_f4d_001_ion_sink.h5
fidasim_f4d_002_ion_sink.h5
...
```

The deterministic and Monte Carlo directories therefore contain identical
basenames in the same case order, allowing the two output collections to be
paired without a separate mapping.

When `monte_carlo/plot_data=.true.`, each HDF5 file is accompanied by a PNG
with the same stem, for example `fidasim_f4d_001_ion_sink.png`. The plot reads
only the self-contained `/test_004` group. It shows the weighted ion-sink
distribution in energy-pitch space and marks the injected neutral with a
green circle. Zero-valued bins are omitted on a logarithmic scale.

The standard FIDASIM sink datasets remain authoritative:

| Dataset | Units | Description |
| --- | --- | --- |
| `/n_sink` | — | Number of stored ion-sink particles; it must equal `n_markers`. |
| `/dens` | `ions/(s*cm^3)` | Fortran logical shape `(species,x,y,z)`; only species 1, cell `[2,2,2]` is nonzero. |
| `/energy` | `keV` | Sampled ion energies. |
| `/pitch` | dimensionless | Sampled ion pitches relative to the magnetic field. |
| `/weight` | `ions/s` | Particle sink weights, equal to the per-marker volumetric rate times the 1 `cm^3` test-cell volume. |
| `/vi` | `cm/s` | Ion velocities written in cylindrical components. |
| `/ri`, `/ri_gc` | `cm`, radians | Particle and guiding-center positions in `(R,Z,Phi)`. |
| `/ind` | — | Beam-grid indices; every particle must originate from `[2,2,2]`. |
| `/atomic_mass` | `amu` | Ion mass associated with every sink particle. |
| `/grid/*` | mixed | Artificial beam-grid definition written by `write_beam_grid`. |

After Fortran finishes, the Python postprocessor checks that every expected
production file exists before modifying any of them. It then validates the
particle count, central-cell indices, positive weights, spatial support, and
agreement between `/dens` and the summed particle weights. It replaces the
derived `/test_004` group on each run, leaving all production datasets intact.

The group contains the same `energy`, `pitch`, `f_array`, `denf`,
`sink_distribution`, `energy_marginal`, `pitch_marginal`, and
`total_reaction_rate` names and units used by the deterministic artifacts. It
also stores `sample_standard_deviation`, `standard_error`, the cell volume,
Monte Carlo controls, neutral parameters, isotope metadata, selected source
location, and input provenance.

This makes every Monte Carlo file self-contained for plotting and comparison.

## Internal design

The Monte Carlo implementation uses a Python normalization layer followed by
one serial Fortran executable.

### Source tree

The input normalizer, Fortran calculation, production writer, Python
postprocessor, and plotting layer are implemented:

```text
02_monte_carlo/
├── normalize_config.py
├── postprocess_monte_carlo.py
├── run.sh
├── build/
├── monte_carlo_tools/
│   ├── __init__.py
│   ├── postprocessing.py
│   └── plotting.py
└── src/
    ├── test_004.f90
    └── modules/
        ├── test_004_types.f90
        ├── test_004_config.f90
        ├── test_004_hdf5_utils.f90
        ├── test_004_hdf5.f90
        ├── test_004_neutral.f90
        ├── test_004_setup.f90
        ├── test_004_sampling.f90
        └── test_004_sink.f90
```

The Test 004 makefile builds the implemented modules and links the FIDASIM and
HDF5 objects required by the artificial test setup.

### Python normalization

`normalize_config.py` validates and normalizes the unified Test 004
configuration, then flattens the required values into one Fortran-facing
block. It uses `test_004_tools.read_monte_carlo_config` and reuses the shared
Test 002
[input-distribution discovery](../README.md#how-test-004-finds-and-validates-its-input-distributions).
It does not open deterministic outputs or require the deterministic
implementation to have run.

Run the normalizer independently with:

```bash
cd regression_tests/test_004/02_monte_carlo
python3 normalize_config.py \
  ../input_config_A.nml \
  build/normalized_input_config.nml
```

The generated file is the only configuration that the Fortran executable will
read and contains one flattened `run_test` block:

| Field | Fortran type/size | Source |
| --- | --- | --- |
| `n_cases` | `integer(Int32)` scalar | Number of discovered distributions. |
| `distribution_files` | Character array `(n_cases)` | Test 002 Stage 2 outputs. |
| `runids` | Character array `(n_cases)` | Source stems followed by `_ion`; `write_sink_profile` appends `_sink.h5`. |
| `tables_filename` | Character scalar | Shared `test_case/tables_filename`. |
| `test_config` | Character scalar | Unified Test 004 configuration, resolved absolutely. |
| `input_distribution_config` | Character scalar | Shared Test 002 Stage 2 configuration. |
| `output_directory` | Character scalar | Derived `monte_carlo` output directory; empty when saving is disabled. |
| `case_comment` | Character scalar | Shared test-case comment. |
| `implementation_comment` | Character scalar | Monte Carlo comment. |
| `n_markers` | `integer(Int64)` scalar | `monte_carlo/n_markers`. |
| `reservoir_size` | `integer(Int32)` scalar | `monte_carlo/reservoir_size`. |
| `seed` | `integer(Int32)` scalar | `monte_carlo/seed`. |
| `save_data` | Logical scalar | `monte_carlo/save_data`. |
| `neutral_density` | `real(Float64)` scalar | Shared neutral density in `cm^-3`. |
| `neutral_energy` | `real(Float64)` scalar | Shared neutral energy in `keV`. |
| `injection_angle` | `real(Float64)` scalar | Shared signed injection angle in degrees. |
| `level_split_method` | Character scalar | Canonical shared selector. |
| `level_decay` | `real(Float64)` scalar | Shared value; zero when inactive. |

All normalized paths are absolute because the file is a transient
Python-to-Fortran interface rather than committed provenance. At most 256
distribution cases are supported. When `save_data=.true.`, normalization also
creates the Monte Carlo output directory required by `write_sink_profile`.
`plot_data` remains Python-facing and is read from the unified configuration
by the postprocessing layer. Production-facing table paths, output paths, run
IDs, and complete sink filenames are limited to 200 characters to match
FIDASIM's fixed character fields.

### Fortran data types

`test_004_types` defines `MonteCarloConfig`, which owns every normalized
scalar and exact-size arrays for the distribution paths and run IDs, and
`DistributionCase`, which owns one loaded smooth distribution and its
location and numerical isotope metadata, and `NeutralParameters`, which owns the
six-level density vector, mass, speed, and Cartesian velocity. It does not use
`libfida` or HDF5. The sink-summary type will be introduced with the module
that needs it.

`test_004_config` provides:

```fortran
call read_config(filename, config)
call print_config(config)
```

Python performs the semantic validation before writing the normalized file.
The Fortran reader only checks the 256-case staging-array limit and verifies
that each case has a distribution path and run ID before allocating the
exact-size arrays. The print routine then exposes every value crossing the
Python-to-Fortran boundary.

Build and exercise this boundary from the Test 004 directory with:

```bash
make
02_monte_carlo/test_004 \
  02_monte_carlo/build/normalized_input_config.nml
```

The current executable reads and prints the normalized configuration, then
loads every discovered Test 002 HDF5 file and prints its grid, distribution,
density, location, and numerical isotope summary.

`test_004_hdf5_utils` contains the test-local, schema-independent HDF5
operations used to read real vectors, one-value arrays, real and integer
scalars, and array dimensions. It also provides common HDF5 status handling.
Its diagnostics describe only the failed HDF5 operation, dataset, and file;
they do not refer to Test 002 or Test 004.

`test_004_hdf5` provides the Test 002 distribution-specific interface:

```fortran
call read_distribution(filename, distribution)
call print_distribution(distribution)
call release_distribution(distribution)
```

It reads the energy and pitch grids, the single spatial slice of the smooth
distribution, `denf`, selected `R` and `Z`, atomic number, mass number, charge
state, and physical mass `A`. The human-readable `species` dataset is retained
in the upstream file for provenance but is not required or read by the Fortran
calculation. The future output writer will add the corresponding Test 004
metadata to the sink files.

`test_004_neutral` provides:

```fortran
call build_neutral_parameters(config, distribution, neutral)
call print_neutral_parameters(neutral)
```

It uses the distribution's physical isotope mass and the configured total
neutral energy to calculate the speed in `cm/s`. The signed injection angle is
then applied as
`[vx,vy,vz] = speed*[sin(angle),0,cos(angle)]`. The configured total density is
split across exactly six atomic levels using the same `ground-only` or
normalized exponential rule as the deterministic implementation.

`test_004_setup` constructs the artificial FIDASIM state required by each
distribution case. It provides:

```fortran
call initialize_test_setup(distribution)
call print_test_setup()
call teardown_test_setup()
```

`initialize_test_setup` initializes the production controls and species mass,
the `6x9x1` interpolation grid, uniform equilibrium, `3x3x3` beam grid, and
the spatially replicated FBM. Its cylindrical interpolation domain,
`R=[0,2.5] cm` and `Z=[-2,2] cm`, encloses the Cartesian beam-grid volume.
The equilibrium mask is established before `make_beam_grid` classifies its 27
cells. The source distribution is deliberately replicated over the complete
interpolation grid to represent a spatially uniform ion population.

`test_004_sink` provides:

```fortran
call calculate_ion_sink(config, distribution)
call print_ion_sink(config)
call finalize_ion_sink(config)
```

These routines own sink allocation, the production marker loop, runtime
assertions, conditional output, and sink-array cleanup. The production writer
owns cleanup when saving is enabled; otherwise `finalize_ion_sink` releases
the arrays directly.

### Main-program lifecycle

`test_004.f90` performs the following orchestration:

```text
read normalized configuration
for each case
    read distribution
    construct neutral parameters
    initialize test setup
    initialize the type-1 neutral population
    reset the RNG and calculate the ion sink
    verify and optionally write the production sink file
    release the neutral population
    teardown test setup
    release distribution
end for
```

After the executable returns, `postprocess_monte_carlo.py` validates the
complete production output collection and appends `/test_004`.

The launcher performs:

```text
normalize unified configuration
build the executable when needed
run the executable with the normalized configuration
plot weighted sink histograms when monte_carlo/plot_data is enabled
```

`plot_monte_carlo.py` bins `/energy`, `/pitch`, and `/weight` on the original
source grid and uses the shared `plot_data_block`. It overlays the same filled
green neutral marker as the deterministic plot at
`(pitch,energy)=(cos(injection_angle),neutral_energy)`.
Comparison-specific plots remain the responsibility of the comparison stage
and will use a separate configuration.

---

**Navigation:** [Previous: Deterministic calculation](../01_deterministic/README.md) | [Up: Test 004](../README.md) | [Next: Comparison](../03_compare/README.md)
