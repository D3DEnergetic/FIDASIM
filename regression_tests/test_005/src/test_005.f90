program test_005
  !+ Verify surface geometry helpers and mesh collision dispatch.
  !+
  !+ The primitive cases cover ordinary interior and boundary hits, misses,
  !+ flight-time scaling, and robustness for very small and very large
  !+ triangles. The dispatch cases then verify transparent openings and
  !+ nearest-wall selection within one synthetic voxel. This test does not
  !+ exercise voxel traversal or end-to-end particle tracking.
  use iso_fortran_env,only:real64
  use libfida,only:pi,ray_type,ray_collision_type,vessel, &
    solve_ray_triangle_intersection,find_ray_surface_intersection, &
    sample_angle_between_dir
  implicit none
  integer,parameter::max_results=32
  character(len=*),parameter::report_path= &
    'output_data/comparison_report.txt'
  type::test_result_type
    character(len=64)::label=''
    logical::passed=.false.
    character(len=256)::detail=''
  end type test_result_type
  type(ray_type)::ray
  type(ray_collision_type)::collision
  type(test_result_type)::results(max_results)
  integer::number_of_results=0
  logical::all_passed

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

  ! Exercise the angle sampler used to place analytic wall-source particles.
  ! In particular, protect the explicit full-circle branch so a 360-degree
  ! source cannot collapse to a zero-width directed interval.
  call test_directed_angle_sampling()

  ! Persist every result before returning a failing process status. This keeps
  ! the diagnostic evidence available to both people and automated runners.
  call write_comparison_report(report_path,all_passed)
  if(.not.all_passed) error stop 'test_005: one or more checks failed'
  print *,'test_005: all ray-triangle and mesh-dispatch checks passed'

contains

  subroutine test_directed_angle_sampling()
    !+ Verify full-circle, clockwise, and counter-clockwise angle sampling.
    !+
    !+ The full-circle case uses equal input and output bounds deliberately:
    !+ ``full`` must ignore those bounds and map the unit random coordinate
    !+ uniformly over [0, 2*pi). Four deterministic quantiles establish the
    !+ complete mapping. The two directed checks protect the existing paths.
    real(real64),parameter::tol=1.d-12
    real(real64),parameter::unit_samples(4)=[0.d0,0.25d0,0.5d0,0.75d0]
    real(real64),parameter::expected_angles(4)=[0.d0,0.5d0*pi,pi,1.5d0*pi]
    real(real64)::actual_angles(4),ccw_angle,cw_angle
    character(len=256)::detail
    integer::sample_index

    do sample_index=1,size(unit_samples)
      actual_angles(sample_index)=sample_angle_between_dir( &
        1.234d0,1.234d0,'full',unit_samples(sample_index))
    enddo
    write(detail,'(a,es12.4)') 'maximum angular error [rad]: ', &
      maxval(abs(actual_angles-expected_angles))
    call record_result('full-circle angle sampling', &
      maxval(abs(actual_angles-expected_angles))<tol,detail)

    call record_result('full-circle capitalization', &
      abs(sample_angle_between_dir(0.d0,0.d0,'FULL',0.125d0) &
      -0.25d0*pi)<tol, &
      'uppercase FULL maps the unit sample over the complete circle')

    ccw_angle=sample_angle_between_dir(1.5d0*pi,0.5d0*pi,'ccw',0.5d0)
    cw_angle=sample_angle_between_dir(1.5d0*pi,0.5d0*pi,'cw',0.5d0)
    write(detail,'(a,es12.4,a,es12.4)') 'ccw midpoint: ',ccw_angle, &
      '; cw midpoint: ',cw_angle
    call record_result('directed angle sampling unchanged', &
      abs(ccw_angle)<tol .and. abs(cw_angle-pi)<tol,detail)
  end subroutine test_directed_angle_sampling

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
    !+ Evaluate one ray-triangle hit and record its complete collision result.
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
    character(len=256)::detail

    ray%origin=origin
    ray%v=velocity
    call solve_ray_triangle_intersection(ray,1,collision)

    if(collision%n_events/=1) then
      write(detail,'(a,i0)') 'expected one collision event; actual events: ', &
        collision%n_events
      call record_result(label,.false.,detail)
      return
    endif
    if(abs(collision%event(1)%s_coll-expected_time)>tol) then
      write(detail,'(a,es12.4,a,es12.4)') 'expected time: ',expected_time, &
        '; actual time: ',collision%event(1)%s_coll
      call record_result(label,.false.,detail)
      return
    endif
    if(maxval(abs(collision%event(1)%p_coll-expected_point))>tol) then
      write(detail,'(a,3(es12.4,1x),a,3(es12.4,1x))') &
        'expected point: ',expected_point,'; actual point: ', &
        collision%event(1)%p_coll
      call record_result(label,.false.,detail)
      return
    endif
    if(maxval(abs(collision%event(1)%normal-[0.d0,0.d0,1.d0]))>tol) then
      write(detail,'(a,3(es12.4,1x),a,3(es12.4,1x))') &
        'expected normal: ',[0.d0,0.d0,1.d0],'; actual normal: ', &
        collision%event(1)%normal
      call record_result(label,.false.,detail)
      return
    endif
    if(collision%event(1)%surface_id/=1 .or. &
       collision%event(1)%region_id/=1) then
      write(detail,'(a,i0,a,i0,a,i0,a,i0)') &
        'expected surface/region IDs: ',1,'/',1,'; actual IDs: ', &
        collision%event(1)%surface_id,'/',collision%event(1)%region_id
      call record_result(label,.false.,detail)
      return
    endif
    call record_result(label,.true., &
      'event count, time, point, normal, and IDs match')
  end subroutine expect_hit

  subroutine expect_miss(origin,velocity,label)
    !+ Evaluate one expected miss and record whether no collision was returned.
    !+
    !+ origin: Ray origin [cm].
    !+ velocity: Physical ray velocity [cm/s].
    !+ label: Human-readable test-case name included in failure messages.
    real(real64),intent(in)::origin(3),velocity(3)
    character(len=*),intent(in)::label
    character(len=256)::detail

    ray%origin=origin
    ray%v=velocity
    call solve_ray_triangle_intersection(ray,1,collision)
    if(collision%n_events==0) then
      call record_result(label,.true.,'no collision event returned')
    else
      write(detail,'(a,i0)') 'expected zero collision events; actual events: ', &
        collision%n_events
      call record_result(label,.false.,detail)
    endif
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
    call expect_dispatch_hit( &
      label='opening exposes farther wall', &
      expected_surface_id=2,expected_time=2.d0)

    ! Case 2: both triangles are walls. Triangle 1 must win at one second even
    ! though triangle 2 appears first in the candidate list.
    vessel%mesh%region(1)%behavior_type='wall'
    call find_ray_surface_intersection(ray,[1,1,1],collision)
    call expect_dispatch_hit( &
      label='candidate order preserves nearest wall', &
      expected_surface_id=1,expected_time=1.d0)

    ! Case 3: both triangles are openings, so no wall collision is returned.
    vessel%mesh%region(1)%behavior_type='opening'
    vessel%mesh%region(2)%behavior_type='opening'
    call find_ray_surface_intersection(ray,[1,1,1],collision)
    call expect_no_dispatch_hit(label='all openings are transparent')
  end subroutine test_mesh_dispatch

  subroutine expect_dispatch_hit(label,expected_surface_id,expected_time)
    !+ Record whether the current dispatch result is the expected nearest hit.
    character(len=*),intent(in)::label
    integer,intent(in)::expected_surface_id
    real(real64),intent(in)::expected_time
    real(real64),parameter::tol=1.d-12
    character(len=256)::detail

    if(collision%n_events/=1) then
      write(detail,'(a,i0)') 'expected one collision event; actual events: ', &
        collision%n_events
      call record_result(label,.false.,detail)
      return
    endif
    if(collision%event(1)%surface_id/=expected_surface_id) then
      write(detail,'(a,i0,a,i0)') 'expected surface ID: ',expected_surface_id, &
        '; actual surface ID: ',collision%event(1)%surface_id
      call record_result(label,.false.,detail)
      return
    endif
    if(abs(collision%event(1)%s_coll-expected_time)>tol) then
      write(detail,'(a,es12.4,a,es12.4)') 'expected time: ',expected_time, &
        '; actual time: ',collision%event(1)%s_coll
      call record_result(label,.false.,detail)
      return
    endif
    call record_result(label,.true., &
      'nearest accepted wall has the expected surface ID and time')
  end subroutine expect_dispatch_hit

  subroutine expect_no_dispatch_hit(label)
    !+ Record whether the current dispatch result contains no wall collision.
    character(len=*),intent(in)::label
    character(len=256)::detail

    if(collision%n_events==0) then
      call record_result(label,.true.,'no wall collision returned')
    else
      write(detail,'(a,i0)') 'expected zero collision events; actual events: ', &
        collision%n_events
      call record_result(label,.false.,detail)
    endif
  end subroutine expect_no_dispatch_hit

  subroutine record_result(label,passed,detail)
    !+ Append one complete named result to the bounded result collection.
    character(len=*),intent(in)::label,detail
    logical,intent(in)::passed
    character(len=4)::status

    if(number_of_results>=max_results) &
      error stop 'test_005: result capacity exceeded'
    number_of_results=number_of_results+1
    results(number_of_results)%label=label
    results(number_of_results)%passed=passed
    results(number_of_results)%detail=detail
    if(passed) then
      status='PASS'
    else
      status='FAIL'
    endif
    write(*,'(2x,a,2x,a)') status,trim(label)
  end subroutine record_result

  subroutine write_comparison_report(path,passed)
    !+ Write all case results, summary counts, and the overall test status.
    !+
    !+ path: Destination report filename. Its parent directory must exist.
    !+ passed: True only when every recorded case passed.
    character(len=*),intent(in)::path
    logical,intent(out)::passed
    integer::report_unit,ios,ii,number_passed,number_failed
    character(len=4)::status

    number_passed=count(results(1:number_of_results)%passed)
    number_failed=number_of_results-number_passed
    passed=number_failed==0
    if(passed) then
      status='PASS'
    else
      status='FAIL'
    endif

    open(newunit=report_unit,file=path,status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'test_005: could not open comparison report'
    write(report_unit,'(a)') '=============================================='
    write(report_unit,'(a,a)') 'OVERALL REGRESSION STATUS: ',status
    write(report_unit,'(a)') '=============================================='
    write(report_unit,'(a)') ''
    write(report_unit,'(a)') 'Test 005 ray-triangle and mesh-dispatch report'
    write(report_unit,'(a)') '=============================================='
    write(report_unit,'(a,i0)') 'Named checks: ',number_of_results
    write(report_unit,'(a,i0)') 'Passed:       ',number_passed
    write(report_unit,'(a,i0)') 'Failed:       ',number_failed
    write(report_unit,'(a)') ''

    do ii=1,number_of_results
      if(results(ii)%passed) then
        status='PASS'
      else
        status='FAIL'
      endif
      write(report_unit,'(a,i0,a,a)') 'Case ',ii,': ',trim(results(ii)%label)
      write(report_unit,'(a,a)') '  Status: ',status
      write(report_unit,'(a,a)') '  Detail: ',trim(results(ii)%detail)
      write(report_unit,'(a)') ''
    enddo

    if(passed) then
      status='PASS'
    else
      status='FAIL'
    endif
    write(report_unit,'(a)') 'Overall summary'
    write(report_unit,'(a)') '---------------'
    write(report_unit,'(a,a)') 'Status: ',status
    write(report_unit,'(a,i0)') 'Named checks: ',number_of_results
    write(report_unit,'(a,i0)') 'Passed:       ',number_passed
    write(report_unit,'(a,i0)') 'Failed:       ',number_failed
    close(report_unit)
    write(*,'(a,a)') 'Wrote report: ',trim(path)
  end subroutine write_comparison_report

end program test_005
