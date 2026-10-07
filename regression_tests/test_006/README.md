# Test 006: surface-mesh setup and tracking

Test 006 verifies FIDASIM's surface-mesh loading, setup, and particle-tracking path using its maintained test inputs. It has no runtime dependency on Surface Mesh Tools or any other project.

## What the test checks

The test calls the production `setup_vessel_mesh` routine and verifies that
FIDASIM:

1. loads the HDF5 vertices, triangles, and named-region catalogue;
2. matches each behavior block to a mesh region by its exact name;
3. loads the deterministic absorbing-wall parameters;
4. reconstructs grid-cell triangle candidates for the active region;
5. excludes the inactive region from those candidates; and
6. tracks one particle into the active wall and absorbs it.

Test 005 separately checks detailed ray-triangle cases, candidate ordering,
and an active transparent opening. Test 006 checks file loading and the
complete setup-to-tracking path.

## Maintained inputs

The `fixtures/` directory contains:

- `tiny_mesh.h5`, the FIDASIM-owned version-0 mesh artifact;
- `mesh_behavior.nml`, the run-specific behavior assigned to its named
  regions.

The HDF5 fixture was generated from the small tagged Gmsh source used during
2026_surface_mesh_tools project development. The generated artifact is now a maintained
FIDASIM regression input. Regeneration must preserve the values documented
below or update this README and the test together.

## Fixture structure

This section explains the contents of the maintained Test 006 fixture. It is
not the complete schema for every valid FIDASIM mesh artifact. The complete
artifact schema belongs in the IPS-FIDASIM reference documentation.

The fixture contains five top-level HDF5 groups:

| Group | Purpose |
|---|---|
| `/metadata` | Identifies the artifact format and records its conventions and provenance |
| `/mesh` | Stores vertices, triangles, triangle-to-region assignments, and source-placement provenance |
| `/regions` | Maps internal region IDs to stable physical-region names |
| `/grid` | Records the FIDASIM beam grid for which the triangle index was constructed |
| `/voxel_index` | Stores the compact mapping from grid cells to candidate triangles |

The metadata values used by this fixture are:

| Attribute | Value | Meaning |
|---|---|---|
| `schema_name` | `fidasim_mesh_artifact` | Identifies the artifact format |
| `schema_version` | `0` | Selects the implemented version-0 contract |
| `length_units` | `cm` | States the coordinate unit |
| `vertex_frame` | `beam_grid_local` | States the coordinate frame used by `/mesh/vertices` |
| `voxel_order` | `i_fastest_then_j_then_k` | Defines the order used to flatten the three-dimensional grid |
| `triangle_id_base` | `1` | Triangle and vertex identifiers are one based |
| `offset_base` | `0` | Compact-index offsets start at zero |

The runtime geometry and region datasets are:

| Dataset | Contents |
|---|---|
| `/mesh/vertices` | Vertex coordinates in centimetres |
| `/mesh/triangles` | One-based vertex IDs for each triangle |
| `/mesh/triangle_region_id` | One-based region ID for each triangle |
| `/regions/id` | Region catalogue IDs |
| `/regions/name` | Region catalogue names |
| `/voxel_index/cell_offsets` | Compact cell-to-triangle index offsets |
| `/voxel_index/triangle_ids` | One-based triangle IDs in the compact index |

The region catalogue maps each region ID to exactly one name. Triangle region
IDs must reference that catalogue, and every catalogue region must have a
corresponding behavior block in `mesh_behavior.nml`. The artifact version and
dataset layout are part of the maintained regression-test contract.

## Mesh geometry

Coordinates are centimetres in the beam-grid-local frame. Vertex and triangle
identifiers are one based.

| Vertex ID | X | Y | Z |
|---:|---:|---:|---:|
| 1 | 0 | 0 | 0 |
| 2 | 1 | 0 | 0 |
| 3 | 1 | 1 | 0 |
| 4 | 0 | 1 | 0 |

The vertices form a 1 cm square in the plane `Z = 0`.

| Triangle ID | Vertex IDs | Region ID | Region name |
|---:|---|---:|---|
| 1 | 1, 2, 3 | 1 | `wall` |
| 2 | 1, 3, 4 | 2 | `opening` |

The HDF5 relationship is:

```text
/mesh/triangle_region_id = [1, 2]
/regions/id               = [1, 2]
/regions/name             = ["wall", "opening"]
```

`triangle_region_id` states which named mesh region contains each triangle.
The behavior file uses the region names rather than relying on catalogue or
block order.

## Source frame

`/mesh/source_frame` records how the original mesh coordinate system was
placed in the laboratory frame before the vertices were converted into
beam-grid-local coordinates:

```text
original mesh coordinates
    -> source-frame placement
laboratory coordinates
    -> beam-grid transformation
stored beam-grid-local vertices
```

The stored `/mesh/vertices` are the final coordinates used by FIDASIM. The
current FIDASIM reader does not use `/mesh/source_frame` during the simulation;
the group is retained as provenance explaining how those vertices were
produced.

For this fixture, the source origin and all three rotation angles are zero.
The source basis and inverse basis are identity matrices. The original mesh
axes, laboratory axes, and beam-grid axes are therefore aligned.

## Behavior model used by the test

| Region | Active? | Behavior | Purpose in this test |
|---|---:|---|---|
| `wall` | Yes | Wall with absorption probability 1 | Deterministic wall collision and absorption |
| `opening` | No | Opening | Verify that an inactive region is excluded from runtime candidates |

An inactive region is ignored by particle tracking: its triangles are not
placed in FIDASIM's runtime grid-cell candidate lists. This differs from an
active opening. An active opening remains part of the geometric index and is
recognized as a transparent surface. Test 005 covers that active-opening
case.

Every artifact region has an explicit behavior block. A missing, duplicate,
or unknown region name is an input error; omission does not mean inactive.

## Grid and voxel index

The voxel index depends on both the raw surface mesh and the selected FIDASIM
beam grid. The production order is therefore:

1. define the beam grid that the FIDASIM calculation will use;
2. provide that same grid and the raw tagged mesh to the artifact producer;
3. construct the cell-to-triangle index and store the grid with it; and
4. load the artifact in FIDASIM and confirm that its stored grid matches the
   active simulation grid.

If the stored and active grids differ, FIDASIM rejects the artifact. It must
not use a triangle index constructed for a different grid.

The fixture uses an identity-aligned `2 × 2 × 2` grid:

| Axis | Minimum | Maximum | Cells |
|---|---:|---:|---:|
| X | -0.5 cm | 1.5 cm | 2 |
| Y | -0.5 cm | 1.5 cm | 2 |
| Z | -0.5 cm | 1.5 cm | 2 |

The grid origin and rotation angles are zero. Its basis and inverse basis are
the identity matrix.

The stored compact candidate index is:

```text
/voxel_index/cell_offsets = [0, 2, 4, 6, 8, 8, 8, 8, 8]
/voxel_index/triangle_ids = [1, 2, 1, 2, 1, 2, 1, 2]
```

Both triangles touch the first four grid cells at `Z = 0`. After behavior is
loaded, triangle 2 is removed because its region is inactive. Triangle 1
remains available in those cells.

## Deterministic particle track

The particle begins at:

```text
position = [0.25, 0.10, 0.25] cm
velocity = [0, 0, -1] cm/s
```

The starting point is above the interior of triangle 1 rather than the
diagonal shared by the two triangles. It travels downward to triangle 1. The
wall absorption probability is one, so the expected result is one
positive-duration track segment followed by absorption.

## What this fixture does not prove

The two triangles form a flat open square, not a closed vessel. This test does
not check:

- watertightness;
- manifold edges;
- consistent winding across a closed surface;
- globally inward- or outward-facing normals; or
- the physical validity of a complete vessel mesh.

Those properties require a separate closed-mesh fixture and focused producer
and consumer tests.

## Running the test

After sourcing the appropriate platform environment, build the maintained
regression executables from the FIDASIM repository root:

```bash
make regression_tests
```

Then run Test 006:

```bash
cd regression_tests/test_006
./run.sh
```

The runner checks the maintained inputs, creates the output directory, and
executes the prebuilt test. The test writes
`output_data/comparison_report.txt` and returns a nonzero status if any check
fails.

The root Makefile supplies the compiler and HDF5 settings required by the test
Makefile. The maintained fixture paths are the defaults. The runner accepts
`ARTIFACT_FILE` and `BEHAVIOR_FILE` environment overrides for diagnostic work.

## Command-line interface

The executable accepts two required positional arguments, in this order:

```text
test_006 <mesh_artifact.h5> <mesh_behavior.nml>
```

| Argument | Meaning |
|---|---|
| 1 | Path to the version-0 HDF5 mesh artifact |
| 2 | Path to the Fortran namelist assigning behavior to the artifact's named regions |

The executable reads these paths with Fortran's `get_command_argument` and
then loads the referenced files. It does not read the Makefile or the runner's
environment variables. The Makefile builds the executable; `run.sh` selects
the paths and passes them as command-line arguments. The behavior settings
remain in the namelist file.

To invoke the executable directly from this test directory:

```bash
mkdir -p output_data
./build/test_006 fixtures/tiny_mesh.h5 fixtures/mesh_behavior.nml
```

Relative input paths and `output_data/comparison_report.txt` are resolved
against the executable's current working directory. Direct callers must
create `output_data` before running. Quote paths containing spaces. Each
input path must fit in the executable's 512-character argument buffer.
Missing, empty, or unreadable/truncated arguments cause an error stop.
Additional arguments are ignored; there are no named options or `--help`
handling.

The runner changes to this test directory, checks that both input files and
the executable exist, creates `output_data`, and invokes the executable. Its
defaults are `fixtures/tiny_mesh.h5` and `fixtures/mesh_behavior.nml`. Override
them for a diagnostic run as follows:

```bash
ARTIFACT_FILE=/absolute/path/tiny_mesh.h5 \
BEHAVIOR_FILE=/absolute/path/mesh_behavior.nml ./run.sh
```

Relative override paths are resolved against this test directory, even when
the runner is launched from elsewhere. Unset or empty overrides use the
defaults. Alternate inputs must preserve the grid, geometry, named regions,
and deterministic behavior expected by this test, as documented above; this
executable is a fixture-specific regression check rather than a general mesh
validation tool.

On a completed run, the executable prints each `PASS` or `FAIL` check and
replaces `output_data/comparison_report.txt`. Exit status zero indicates that
all checks passed; a nonzero status indicates an input, setup, output, or
assertion failure. The runner exits with status 2 when an input file or the
executable is missing, and otherwise propagates the executable's status.
Failures before report generation can leave an older report in place, so
check the exit status before interpreting the report as the current result.
