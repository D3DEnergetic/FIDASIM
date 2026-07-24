[Regression tests](../README.md) / Test 002

**Navigation:** [Previous: Test 001](../test_001/README.md) | [Up: Regression tests](../README.md) | [Next: Reference data](01_reference/README.md)

---

# Test 002: CQL3D-to-FIDASIM distribution conversion

`test_002` validates the conversion of a CQL3D distribution expressed on
normalized proper-velocity and pitch-angle coordinates, `f(u_bar, theta)`, into
the energy-pitch density, `F(E, pitch)`, required by FIDASIM.

The workflow will establish the relativistic coordinate transformation, the
phase-space Jacobian, the output units, and preservation of density, parallel
temperature, and perpendicular temperature.

## Workflow

| Stage | Purpose | Principal outputs |
| --- | --- | --- |
| [`01_reference`](01_reference/README.md) | Extract trusted single-location CQL3D distributions and calculate their velocity-space moments. | Compact HDF5 reference files, contour plots, and a moment report. |
| [`02_run_test`](02_run_test/README.md) | Convert each reference into the FIDASIM energy-pitch representation and calculate its moments. | FIDASIM HDF5 files, contour plots, and a moment report. |
| [`03_compare`](03_compare/README.md) | Compare the moments before and after conversion. | Summary report and relative-difference plot. |

Together, the three stages isolate the coordinate conversion and quantify the
moment error introduced by remapping onto the uniform FIDASIM grid.

## Run all dataset collections

From the `test_002` directory, run:

```bash
./run.sh
```

The launcher runs the configuration filenames listed in its `config_files`
array through `02_run_test` and `03_compare`, in that order. Add or remove
entries in that short list to control the collections included in a regression
run. Each listed filename must exist in both stages.

Stage 1 reference generation is intentionally excluded. A developer must run
the appropriate `01_reference/run.sh input_config_<letter>.nml` command
manually when the committed reference data needs to be regenerated.

## Dataset collections

Each collection uses a letter consistently in its configuration filenames and
output directories across all three stages. The descriptive `comment` field in
each configuration records the physical meaning of that collection.

| Collection | Output directory | Beam energy | Injection angle | Source description |
| --- | --- | ---: | ---: | --- |
| `A` | `dataset_A` | 60 keV | 45° | Current CQL3D WHAM expander F4D case, sampled at seven axial locations from 0 through 60 cm. |

The configuration chain for this collection is:

```text
01_reference/input_config_A.nml
                         ↓
02_run_test/input_config_A.nml
                         ↓
03_compare/input_config_A.nml
```

New conditions should use the next letter in the same three configuration
filenames and in each stage's `output_data/dataset_<letter>/` directory. Add
the configuration filename to `run.sh` and document its meaning in this table.

---

**Navigation:** [Previous: Test 001](../test_001/README.md) | [Up: Regression tests](../README.md) | [Next: Reference data](01_reference/README.md)
