program test_005
  ! Verify the mesh ray-triangle primitive against the reviewed Python cases.
  use iso_fortran_env,only:real64,int32
  use libfida,only:ray_type,ray_collision_type,vessel, &
    solve_ray_triangle_intersection
  implicit none
  type(ray_type)::ray
  type(ray_collision_type)::collision

  allocate(vessel%mesh%vertices(3,3))
  allocate(vessel%mesh%triangles(3,1))
  allocate(vessel%mesh%triangle_region_id(1))
  vessel%mesh%triangles(:,1)=[1_int32,2_int32,3_int32]
  vessel%mesh%triangle_region_id(1)=1

  call set_triangle(1.d0)
  call expect_hit([0.5d0,0.5d0,1.d0],[0.d0,0.d0,-1.d0],1.d0, &
    [0.5d0,0.5d0,0.d0],'interior hit')
  call expect_hit([0.5d0,0.5d0,-1.d0],[0.d0,0.d0,1.d0],1.d0, &
    [0.5d0,0.5d0,0.d0],'back-side hit')
  call expect_hit([1.d0,1.d0,1.d0],[0.d0,0.d0,-1.d0],1.d0, &
    [1.d0,1.d0,0.d0],'edge hit')
  call expect_hit([1.d0,1.d0+1.d-12,1.d0],[0.d0,0.d0,-1.d0],1.d0, &
    [1.d0,1.d0+1.d-12,0.d0],'barycentric-tolerance hit')
  call expect_hit([0.d0,0.d0,1.d0],[0.d0,0.d0,-1.d0],1.d0, &
    [0.d0,0.d0,0.d0],'vertex hit')
  call expect_hit([0.5d0,0.5d0,0.d0],[0.d0,0.d0,1.d0],0.d0, &
    [0.5d0,0.5d0,0.d0],'zero-time hit')
  call expect_hit([0.5d0,0.5d0,1.d0],[0.d0,0.d0,-2.d0],0.5d0, &
    [0.5d0,0.5d0,0.d0],'velocity-scaled hit')
  call expect_miss([1.5d0,1.5d0,1.d0],[0.d0,0.d0,-1.d0],'outside miss')
  call expect_miss([0.5d0,0.5d0,1.d0],[1.d0,0.d0,0.d0],'parallel miss')
  call expect_miss([0.5d0,0.5d0,1.d0],[0.d0,0.d0,1.d0],'behind miss')

  call set_triangle(1.d-6)
  call expect_hit([0.5d-6,0.5d-6,1.d0],[0.d0,0.d0,-1.d0],1.d0, &
    [0.5d-6,0.5d-6,0.d0],'small triangle')
  call set_triangle(1.d6)
  call expect_hit([0.5d6,0.5d6,1.d0],[0.d0,0.d0,-1.d0],1.d0, &
    [0.5d6,0.5d6,0.d0],'large triangle')

  print *,'test_005: all ray-triangle checks passed'

contains

  subroutine set_triangle(scale)
    real(real64),intent(in)::scale
    vessel%mesh%vertices(:,1)=[0.d0,0.d0,0.d0]
    vessel%mesh%vertices(:,2)=[2.d0*scale,0.d0,0.d0]
    vessel%mesh%vertices(:,3)=[0.d0,2.d0*scale,0.d0]
  end subroutine set_triangle

  subroutine expect_hit(origin,velocity,expected_time,expected_point,label)
    real(real64),intent(in)::origin(3),velocity(3),expected_time,expected_point(3)
    character(len=*),intent(in)::label
    real(real64),parameter::tol=1.d-11
    ray%origin=origin
    ray%v=velocity
    call solve_ray_triangle_intersection(ray,1_int32,collision)
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
    real(real64),intent(in)::origin(3),velocity(3)
    character(len=*),intent(in)::label
    ray%origin=origin
    ray%v=velocity
    call solve_ray_triangle_intersection(ray,1_int32,collision)
    if(collision%n_events/=0) error stop label//': expected a miss'
  end subroutine expect_miss

end program test_005
