[Regression tests](../../README.md) / [Test 004](../README.md) / Monte Carlo test

**Navigation:** [Previous: Reference calculation](../01_reference/README.md) | [Up: Test 004](../README.md) | [Next: Comparison](../03_compare/README.md)

---

# Monte Carlo ion-sink test

Stage 2 samples the smooth Test 002 energy-pitch distributions and executes the
ion-sink portion of `calculate_dcx_process` using FIDASIM source procedures.
Neutral tracking, attenuation, and ion-birth processing are excluded.

The Stage 2 configuration points to the matching Test 004 Stage 1
configuration. Stage 2 inherits the upstream distribution collection, atomic
tables, neutral physical parameters, and reference output location through
that file. This provides one source of truth for every physical condition
shared by the deterministic and Monte Carlo calculations.

## Input configuration

```fortran
&run_test
  comment = 'Dataset A: Monte Carlo ion-sink calculation'
  reference_config = '../01_reference/input_config_A.nml'
  n_markers = 1000000
  reservoir_size = 50
  seed = 12345
/

&save_data_block
  output_directory = 'output_data/dataset_A'
/
```

Relative paths are resolved from the configuration file containing them.

### `run_test` schema

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `comment` | String scalar | — | No | None | Human-readable description of the Monte Carlo collection. |
| `reference_config` | String scalar | — | Yes | None | Existing Test 004 Stage 1 configuration defining the matching deterministic case. |
| `n_markers` | Integer scalar | — | Yes | None | Number of ion markers sampled independently for each distribution; must be positive. |
| `reservoir_size` | Integer scalar | — | Yes | None | Number of representative neutral particles stored in the central-cell reservoir; must be positive. |
| `seed` | Integer scalar | — | Yes | None | Positive serial RNG seed. The RNG is reset to this seed for every distribution. |

The referenced Stage 1 configuration must be valid and must enable reference
data output. Stage 2 obtains from it:

- the Test 002 Stage 2 configuration and ordered smooth distributions;
- the standard FIDASIM atomic-tables file;
- neutral density, energy, injection angle, and level-split model;
- the directory containing the matching deterministic outputs.

### `save_data_block` schema

| Variable | Type/size | Units | Required | Default | Allowed values and description |
| --- | --- | --- | --- | --- | --- |
| `output_directory` | String scalar | — | Yes | None | Destination for Monte Carlo sink HDF5 files; created when needed. |

Output is mandatory for this regression stage.

## Fixed numerical fixture

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

The central neutral density is assigned once as a six-level vector.
The reservoir contains `reservoir_size` identical neutral velocities with
equal weights `density/reservoir_size`.

## Per-case workflow

For every smooth distribution discovered through `reference_config`:

1. Read its energy grid, pitch grid, distribution, `denf`, isotope, and charge.
2. Construct the beam grid, interpolation grid, uniform equilibrium, fields,
   and spatially replicated FBM.
3. Initialize the center-cell type-1 neutral density and reservoir.
4. Allocate `sink%part(n_markers)` and
   `sink%dens(1,3,3,3)`, set `sink%cnt=1`, and zero the density.
5. Reset the serial FIDASIM RNG to the configured seed.
6. For each marker:
   - call `mc_sample_ion_f4d_gc` in cell `[2,2,2]`;
   - require a nonthermal sample;
   - require `ind_p=[2,2,2]` and `rp=rgc` because `flr=0`;
   - call `get_total_cx_rate` for `[nbif_type]`;
   - calculate the per-marker volumetric sink weight
     `denf*sum(rates)/n_markers`;
   - call `store_sinks`;
   - call `store_sink_particle`.
7. Save the total center-cell sink rate before writing.
8. Call `write_sink_profile`.
9. Reopen the sink file and append the Test 004 metadata described below.
10. Release the remaining case-specific fixture allocations before loading the
    next distribution.

`write_sink_profile` resets `sink%cnt` and deallocates both `sink%part` and
`sink%dens`. The next case must allocate them again. Stage 2 must not manually
deallocate those two arrays after the writer returns.

The workflow does not call `track_to_wall`, COLRAD, `store_neutrals`,
`store_births`, `store_birth_particle`, or
`write_particle_tracks_to_file`.

## Output contract

Set `inputs%result_dir` to the configured output directory and set
`inputs%runid` to the source distribution stem. The production writer then
creates correspondingly named files:

```text
fidasim_f4d_001_sink.h5
fidasim_f4d_002_sink.h5
...
```

The source basename and case ordering match Stage 1, allowing deterministic
and Monte Carlo files to be paired without a separate mapping.

The standard FIDASIM sink datasets remain authoritative:

| Dataset | Description |
| --- | --- |
| `/n_sink` | Number of stored ion-sink particles. |
| `/dens` | Species- and beam-cell-resolved volumetric sink rate. |
| `/energy`, `/pitch`, `/weight` | Particle coordinates and reaction weights used to construct the MC sink histogram. |
| `/vi`, `/ri`, `/ri_gc`, `/ind` | Stored particle velocity, position, guiding-centre position, and beam-cell index. |
| `/atomic_mass` | Ion mass associated with every sink particle. |
| `/grid/*` | Artificial beam-grid definition. |

After `write_sink_profile` closes the file, a test-specific writer reopens it
and appends `/test_004` metadata:

- source energy grid, pitch grid, smooth distribution, and `denf`;
- source distribution and configuration provenance;
- isotope and charge;
- neutral velocity and six-level density vector;
- neutral energy, injection angle, and level-split settings;
- `n_markers`, `reservoir_size`, and seed;
- saved scalar center-cell and total sink rates;
- fixed fixture and coordinate-convention metadata.

This preserves the production sink format while making every Monte Carlo file
self-contained for the future Python plotting and comparison workflow.

---

**Navigation:** [Previous: Reference calculation](../01_reference/README.md) | [Up: Test 004](../README.md) | [Next: Comparison](../03_compare/README.md)
