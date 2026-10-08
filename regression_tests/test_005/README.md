[Regression tests](../README.md) / Test 005

**Navigation:** [Previous: Test 004](../test_004/README.md) | [Up: Regression tests](../README.md) | Next: —

---

# Test 005: surface geometry and mesh collision dispatch

`test_005` is a deterministic, self-checking regression test for surface
geometry helpers and the mesh collision capability. It checks analytic
wall-source angle sampling, isolates the ray-triangle intersection primitive,
and verifies the policy that selects a physical wall collision from candidate
triangles in one mesh voxel.

The expected answers are embedded in the Fortran test program. The test does
not require an external reference dataset or a separate scientific comparison
script.

## Scope and terminology

The test uses the following terms deliberately:

- A **candidate triangle** is a broad-phase possibility associated with a
  voxel. It is not yet a confirmed hit.
- A **geometric intersection** is the mathematical result returned by the
  ray-triangle primitive.
- A **physical collision event** is an accepted intersection after mesh-region
  behavior has been applied. Openings are transparent; walls may collide.
- `s_coll` is the physical collision time in seconds. Travelled distance is
  ray speed multiplied by `s_coll`.
- The **primitive solver** evaluates one ray and one triangle. The **mesh
  dispatcher** evaluates candidates, applies region behavior, and retains the
  nearest accepted wall.

Implemented coverage includes full-circle and directed-angle sampling,
ray-triangle geometry, boundary and tolerance behavior, scale robustness,
opening transparency, and nearest-wall selection.
The test does not cover artifact loading, voxel traversal, full particle
tracking, or validation that a vessel mesh is watertight and consistently
oriented.

## FIDASIM code under test

| Module | Procedure | Responsibility |
| --- | --- | --- |
| `libfida` | `sample_angle_between_dir` | Maps a unit random coordinate onto a full, clockwise, or counter-clockwise angular interval. |
| `libfida` | `solve_ray_triangle_intersection` | Computes the geometric intersection between one directed ray and one triangle. |
| `libfida` | `find_ray_surface_intersection` | Dispatches mesh candidates, rejects openings, and returns the nearest accepted wall collision. |

Supporting mesh allocations, synthetic voxel data, result accumulation, and
report writing belong to the test harness. They are not independent production
capabilities tested here.

## Concept-to-code mapping

| Concept | Implementation name | Units or meaning |
| --- | --- | --- |
| Ray origin | `ray%origin` | Cartesian position, `cm` |
| Ray velocity | `ray%v` | Cartesian velocity, `cm/s` |
| Collision time | `collision%event(1)%s_coll` | `s` |
| Collision point | `collision%event(1)%p_coll` | Cartesian position, `cm` |
| Triangle vertices | `vessel%mesh%vertices` | Cartesian position, `cm` |
| Triangle connectivity | `vessel%mesh%triangles` | Three vertex IDs per triangle |
| Triangle-to-region map | `vessel%mesh%triangle_region_id` | Region ID for each triangle |
| Region policy | `behavior_type` | `wall` or `opening` in this test |

## Test fixtures

### Ray-triangle primitive

The reusable fixture is a right triangle in the XY plane with vertices
`(0,0,0)`, `(2L,0,0)`, and `(0,2L,0)` centimetres. The scale factor `L` is set
to `1`, `10^-6`, and `10^6` to exercise ordinary, very small, and very large
triangles without changing the underlying geometry.

### Mesh dispatch

The dispatch fixture contains two parallel triangles with identical XY
coordinates. The near triangle lies at `z = 0 cm`; the far triangle lies at
`z = -1 cm`. A ray starts one centimetre above the near triangle and travels
straight downward at one centimetre per second.

One synthetic voxel contains both triangle IDs in far-then-near order. This
reversed candidate order proves that the dispatcher selects the nearest wall
by collision time rather than accepting the first candidate in the list.

## Protected behavior

| Group | Named checks | Purpose |
| --- | --- | --- |
| Angular sampling | `full-circle angle sampling`, `full-circle capitalization`, `directed angle sampling unchanged` | Confirm uniform full-circle quantiles, accepted capitalization, and preservation of the clockwise and counter-clockwise paths. |
| Ordinary hits | `interior hit`, `back-side hit` | Confirm interior intersections from both triangle faces. |
| Boundaries and tolerance | `edge hit`, `barycentric-tolerance hit`, `vertex hit`, `zero-time hit` | Exercise exact boundaries, the accepted near-edge tolerance, and a ray beginning on the triangle. |
| Time scaling | `velocity-scaled hit` | Confirm that doubling ray speed halves collision time while preserving the collision point. |
| Misses | `outside miss`, `parallel miss`, `behind miss` | Reject intersections outside the triangle, parallel rays, and intersections behind the ray origin. |
| Scale robustness | `small triangle interior hit`, `large triangle interior hit` | Check the same interior geometry at the two scale extremes. |
| Mesh dispatch | `opening exposes farther wall`, `candidate order preserves nearest wall`, `all openings are transparent` | Check opening transparency, nearest-wall selection, and the no-wall result. |

For an expected primitive hit, the test checks the event count, collision time,
point, normal, surface ID, and region ID. Expected misses require zero events.
Dispatch hits check the selected surface ID and time.

## Workflow

```text
Embedded geometry and expected values
                  |
                  v
        Test 005 Fortran executable
                  |
                  v
       Record every named check result
                  |
                  v
       output_data/comparison_report.txt
                  |
                  v
       PASS, or nonzero failure status
```

The program does not stop at the first failed assertion. Each helper records
one named result, later cases continue to run, and the complete report is
written before the program returns a nonzero status. This preserves diagnostic
evidence while still allowing an automated runner to detect failure.

## Directory structure

```text
test_005/
├── README.md                    # this portal
├── makefile                     # compilation and cleanup
├── run.sh                       # execute the prebuilt self-checking test
├── src/
│   └── test_005.f90             # fixtures, assertions, and report writer
├── build/                       # generated objects and executable
└── output_data/                 # generated comparison report
```

The `build` and `output_data` contents are generated and ignored by Git.

## Build and run

Configure the supported compiler and HDF5 environment, then build the
regression targets from the FIDASIM repository root:

```bash
cd "$FIDASIM_DIR"
source env/<platform>.sh
make regression_tests
```

The regression target compiles Test 005 but does not execute it.

To rerun only the already-built Test 005 executable:

```bash
cd "$FIDASIM_DIR/regression_tests/test_005"
./run.sh
```

The runner creates `output_data` when needed and executes
`build/test_005`. Rebuild through the repository-level makefile after changing
the test or FIDASIM source.

Remove Test 005 compiler artifacts and its generated report with the Test 005
makefile after the build environment has been configured:

```bash
make clean
```

## Comparison report

The persistent report is written to
`output_data/comparison_report.txt`. It contains:

- the overall regression status;
- the number of named, passed, and failed checks;
- one status and diagnostic detail for each named check; and
- a repeated overall summary at the end.

The executable writes the report before returning. It returns a nonzero status
when any named check fails, so `run.sh` directly provides the regression
verdict to both users and the complete-suite runner.

## Debugging one failed check

1. Run `./run.sh` from the Test 005 directory.
2. Find the failed case label in the terminal output or comparison report.
3. Read its expected-versus-actual detail.
4. Find the call with the same `label` keyword in `src/test_005.f90`.
5. Trace that case through `expect_hit`, `expect_miss`, or the dispatch helper
   and then into the relevant `libfida` procedure.

---

**Navigation:** [Previous: Test 004](../test_004/README.md) | [Up: Regression tests](../README.md) | Next: —
