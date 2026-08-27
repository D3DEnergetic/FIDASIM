program test_005
  !+ Verify the mesh ray-triangle primitive and mesh collision dispatch.
  !+
  !+ The primitive cases cover ordinary interior and boundary hits, misses,
  !+ flight-time scaling, and robustness for very small and very large
  !+ triangles. The dispatch cases then verify transparent openings and
  !+ nearest-wall selection within one synthetic voxel. This test does not
  !+ exercise voxel traversal or end-to-end particle tracking.
  use iso_fortran_env,only:real64
  use libfida,only:ray_type,ray_collision_type,vessel, &
    solve_ray_triangle_intersection,find_ray_surface_intersection
  implicit none
  type(ray_type)::ray
  type(ray_collision_type)::collision

  ! Build one reusable triangle fixture. set_triangle reconfigures the same
  ! three vertices for the unit, small, and large scale test phases below.
  allocate(vessel%mesh%vertices(3,3))
  allocate(vessel%mesh%triangles(3,1))
  allocate(vessel%mesh%triangle_region_id(1))
  vessel%mesh%triangles(:,1)=[1,2,3]
  vessel%mesh%triangle_region_id(1)=1

  ! Unit-scale primitive cases: exercise ordinary hits, triangle boundaries,
  ! the barycentric tolerance, time scaling, and the principal miss modes.
  call set_triangle(scale=1.d0)
  call expect_hit(origin=[0.5d0,0.5d0,1.d0], &
    velocity=[0.d0,0.d0,-1.d0],expected_time=1.d0, &
    expected_point=[0.5d0,0.5d0,0.d0],label='interior hit')
  call expect_hit(origin=[0.5d0,0.5d0,-1.d0], &
    velocity=[0.d0,0.d0,1.d0],expected_time=1.d0, &
    expected_point=[0.5d0,0.5d0,0.d0],label='back-side hit')
  call expect_hit(origin=[1.d0,1.d0,1.d0], &
    velocity=[0.d0,0.d0,-1.d0],expected_time=1.d0, &
    expected_point=[1.d0,1.d0,0.d0],label='edge hit')
  call expect_hit(origin=[1.d0,1.d0+1.d-12,1.d0], &
    velocity=[0.d0,0.d0,-1.d0],expected_time=1.d0, &
    expected_point=[1.d0,1.d0+1.d-12,0.d0], &
    label='barycentric-tolerance hit')
  call expect_hit(origin=[0.d0,0.d0,1.d0], &
    velocity=[0.d0,0.d0,-1.d0],expected_time=1.d0, &
    expected_point=[0.d0,0.d0,0.d0],label='vertex hit')
  call expect_hit(origin=[0.5d0,0.5d0,0.d0], &
    velocity=[0.d0,0.d0,1.d0],expected_time=0.d0, &
    expected_point=[0.5d0,0.5d0,0.d0],label='zero-time hit')
  call expect_hit(origin=[0.5d0,0.5d0,1.d0], &
    velocity=[0.d0,0.d0,-2.d0],expected_time=0.5d0, &
    expected_point=[0.5d0,0.5d0,0.d0],label='velocity-scaled hit')
  call expect_miss(origin=[1.5d0,1.5d0,1.d0], &
    velocity=[0.d0,0.d0,-1.d0],label='outside miss')
  call expect_miss(origin=[0.5d0,0.5d0,1.d0], &
    velocity=[1.d0,0.d0,0.d0],label='parallel miss')
  call expect_miss(origin=[0.5d0,0.5d0,1.d0], &
    velocity=[0.d0,0.d0,1.d0],label='behind miss')

  ! Small-scale robustness: reconfigure the fixture as a 10^-6 scale
  ! triangle and verify that an ordinary interior hit remains stable.
  call set_triangle(scale=1.d-6)
  call expect_hit(origin=[0.5d-6,0.5d-6,1.d0], &
    velocity=[0.d0,0.d0,-1.d0],expected_time=1.d0, &
    expected_point=[0.5d-6,0.5d-6,0.d0], &
    label='small triangle interior hit')

  ! Large-scale robustness: repeat the same interior-hit geometry with a
  ! 10^6 scale triangle to exercise the opposite numerical extreme.
  call set_triangle(scale=1.d6)
  call expect_hit(origin=[0.5d6,0.5d6,1.d0], &
    velocity=[0.d0,0.d0,-1.d0],expected_time=1.d0, &
    expected_point=[0.5d6,0.5d6,0.d0], &
    label='large triangle interior hit')

  ! Geometry-aware dispatch: move beyond the isolated primitive and verify
  ! behavior filtering and nearest-hit selection within a synthetic voxel.
  call test_mesh_dispatch()

  print *,'test_005: all ray-triangle and mesh-dispatch checks passed'

contains

  subroutine set_triangle(scale)
    !+ Configure the reusable fixture as a right triangle in the XY plane.
    !+
    !+ scale: Dimensionless factor applied to the two-centimetre X and Y edge
    !+ lengths. It supports the unit, very small, and very large test phases.
    real(real64),intent(in)::scale
    vessel%mesh%vertices(:,1)=[0.d0,0.d0,0.d0]
    vessel%mesh%vertices(:,2)=[2.d0*scale,0.d0,0.d0]
    vessel%mesh%vertices(:,3)=[0.d0,2.d0*scale,0.d0]
  end subroutine set_triangle

  subroutine expect_hit(origin,velocity,expected_time,expected_point,label)
    !+ Assert one ray-triangle hit and validate its complete collision record.
    !+
    !+ origin: Ray origin [cm].
    !+ velocity: Physical ray velocity [cm/s].
    !+ expected_time: Expected flight time to the triangle [s].
    !+ expected_point: Expected collision point [cm].
    !+ label: Human-readable test-case name included in failure messages.
    !+
    !+ The helper checks event count, time, point, normal, triangle ID, and
    !+ region ID using the shared unit-triangle fixture.
    real(real64),intent(in)::origin(3),velocity(3),expected_time,expected_point(3)
    character(len=*),intent(in)::label
    real(real64),parameter::tol=1.d-11
    ray%origin=origin
    ray%v=velocity
    call solve_ray_triangle_intersection(ray,1,collision)
    if(collision%n_events/=1) error stop label//': expected a hit'
    if(abs(collision%event(1)%s_coll-expected_time)>tol) &
      error stop label//': incorrect time'
    if(maxval(abs(collision%event(1)%p_coll-expected_point))>tol) &
      error stop label//': incorrect point'
    if(maxval(abs(collision%event(1)%normal-[0.d0,0.d0,1.d0]))>tol) &
      error stop label//': incorrect normal'
    if(collision%event(1)%surface_id/=1 .or. &
       collision%event(1)%region_id/=1) error stop label//': incorrect IDs'
  end subroutine expect_hit

  subroutine expect_miss(origin,velocity,label)
    !+ Assert that a ray produces no intersection with the fixture triangle.
    !+
    !+ origin: Ray origin [cm].
    !+ velocity: Physical ray velocity [cm/s].
    !+ label: Human-readable test-case name included in failure messages.
    real(real64),intent(in)::origin(3),velocity(3)
    character(len=*),intent(in)::label
    ray%origin=origin
    ray%v=velocity
    call solve_ray_triangle_intersection(ray,1,collision)
    if(collision%n_events/=0) error stop label//': expected a miss'
  end subroutine expect_miss

  subroutine test_mesh_dispatch()
    !+ Verify mesh behavior filtering and nearest-wall selection.
    !+
    !+ One synthetic voxel contains two congruent triangles parallel to the XY
    !+ plane. Triangle 1 is nearer the downward ray at z=0; triangle 2 is
    !+ farther away at z=-1. Their candidate IDs are deliberately stored in
    !+ the reverse order [2,1]. The cases verify that a nearer opening exposes
    !+ the wall behind it, candidate order cannot override collision distance,
    !+ and two openings produce no wall collision.

    ! Replace the primitive fixture with two parallel triangles and a one-cell
    ! candidate map. The ray begins one centimetre above the nearer triangle
    ! and travels perpendicularly downward at one centimetre per second.
    deallocate(vessel%mesh%vertices,vessel%mesh%triangles, &
      vessel%mesh%triangle_region_id)
    allocate(vessel%mesh%vertices(3,6),vessel%mesh%triangles(3,2), &
      vessel%mesh%triangle_region_id(2),vessel%mesh%region(2))
    allocate(vessel%map(1,1,1))
    allocate(vessel%map(1,1,1)%surface_id,source=[2,1])
    vessel%map(1,1,1)%has_surfaces=.true.
    vessel%geometry_type='mesh'

    vessel%mesh%vertices(:,1)=[0.d0,0.d0,0.d0]
    vessel%mesh%vertices(:,2)=[2.d0,0.d0,0.d0]
    vessel%mesh%vertices(:,3)=[0.d0,2.d0,0.d0]
    vessel%mesh%vertices(:,4)=[0.d0,0.d0,-1.d0]
    vessel%mesh%vertices(:,5)=[2.d0,0.d0,-1.d0]
    vessel%mesh%vertices(:,6)=[0.d0,2.d0,-1.d0]
    vessel%mesh%triangles(:,1)=[1,2,3]
    vessel%mesh%triangles(:,2)=[4,5,6]
    vessel%mesh%triangle_region_id=[1,2]
    vessel%mesh%region(1)%behavior_type='opening'
    vessel%mesh%region(2)%behavior_type='wall'
    ray%origin=[0.5d0,0.5d0,1.d0]
    ray%v=[0.d0,0.d0,-1.d0]

    ! Case 1: triangle 1 is an opening, so the ray passes through it and hits
    ! triangle 2, the farther wall, at two seconds.
    call find_ray_surface_intersection(ray,[1,1,1],collision)
    if(collision%n_events/=1 .or. collision%event(1)%surface_id/=2 .or. &
       abs(collision%event(1)%s_coll-2.d0)>1.d-12) &
      error stop 'mesh dispatch: opening did not expose farther wall'

    ! Case 2: both triangles are walls. Triangle 1 must win at one second even
    ! though triangle 2 appears first in the candidate list.
    vessel%mesh%region(1)%behavior_type='wall'
    call find_ray_surface_intersection(ray,[1,1,1],collision)
    if(collision%n_events/=1 .or. collision%event(1)%surface_id/=1 .or. &
       abs(collision%event(1)%s_coll-1.d0)>1.d-12) &
      error stop 'mesh dispatch: did not select nearest wall'

    ! Case 3: both triangles are openings, so no wall collision is returned.
    vessel%mesh%region(1)%behavior_type='opening'
    vessel%mesh%region(2)%behavior_type='opening'
    call find_ray_surface_intersection(ray,[1,1,1],collision)
    if(collision%n_events/=0) error stop 'mesh dispatch: opening was not transparent'
  end subroutine test_mesh_dispatch

end program test_005
