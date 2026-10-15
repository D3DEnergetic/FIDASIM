program test_006
  !+ Verify FIDASIM's complete surface-mesh setup and tracking path.
  !+
  !+ The test owns a small HDF5 mesh artifact and a run-specific behavior
  !+ namelist under ``fixtures/``. It loads both through the production setup
  !+ routine, inspects the reconstructed runtime state, and tracks one
  !+ particle into a deterministic absorbing wall.

  use iso_fortran_env,only:int32,real64
  use libfida,only:ParticleTrack,apply_wall_interaction,beam_grid, &
    collision_event_type,cross_product,inputs,setup_vessel_mesh,track_to_wall, &
    track_exit_grid,track_exit_path_limit,vessel,wall_type
  use utilities,only:rng,rng_init
  implicit none

  integer,parameter::max_results=32
  integer,parameter::number_of_samples=100000
  character(len=*),parameter::report_path='output_data/comparison_report.txt'

  type::test_result_type
    character(len=72)::label=''
    logical::passed=.false.
    character(len=256)::detail=''
  end type test_result_type

  type(test_result_type)::results(max_results)
  type(ParticleTrack),allocatable::tracks(:)
  character(len=512)::artifact_path,behavior_path,detail
  integer::argument_status,i,j,k,reference_index,exit_reason,periodic_wraps
  integer::number_of_results=0,active_references,inactive_references,sample_index
  integer(int32)::ntrack
  logical::absorbed,all_passed
  real(real64)::incident_speed,maximum_speed_error,minimum_mu,mean_mu
  real(real64)::mean_azimuth_x,mean_azimuth_y,mu,azimuth_x,azimuth_y
  real(real64)::normal(3),tangent_x(3),tangent_y(3),velocity(3),direction(3)
  type(collision_event_type)::event
  type(wall_type)::elastic_wall

  ! Read the command-line fixture paths supplied by run.sh. Normal regression
  ! runs use the maintained files under fixtures/, while explicit arguments
  ! allow diagnostic runs against another compatible artifact or behavior.
  call get_command_argument(1,artifact_path,status=argument_status)
  if(argument_status/=0 .or. len_trim(artifact_path)==0) &
    error stop 'test_006: artifact path argument is required'
  call get_command_argument(2,behavior_path,status=argument_status)
  if(argument_status/=0 .or. len_trim(behavior_path)==0) &
    error stop 'test_006: behavior path argument is required'

  ! Recreate the exact beam grid recorded in the HDF5 fixture. Production
  ! setup must reject an artifact whose stored grid differs from this active
  ! FIDASIM grid because its cell-to-triangle index would no longer be valid.
  call initialize_test_grid()
  vessel%geometry_type='mesh'
  vessel%geometry_file=trim(artifact_path)
  vessel%behavior_file=trim(behavior_path)
  allocate(vessel%map(beam_grid%nx,beam_grid%ny,beam_grid%nz))
  inputs%verbose=0
  inputs%enable_wall_interactions=.true.
  ! Wall interaction samples the FIDASIM random stream before selecting an
  ! outcome, even when absorption probability is one. This test runs serially,
  ! so one deterministic stream supplies the required runtime state.
  allocate(rng(1))
  call rng_init(rng(1),12345_int32)

  ! Exercise artifact loading, name-based behavior resolution, and native
  ! voxel reconstruction through the production setup coordinator.
  call setup_vessel_mesh()

  call record_result('artifact geometry loaded', &
    allocated(vessel%mesh%vertices) .and. &
    allocated(vessel%mesh%triangles) .and. &
    size(vessel%mesh%vertices,2)==4 .and. &
    size(vessel%mesh%triangles,2)==2, &
    'expected four vertices and two triangles')

  call record_result('artifact region catalogue loaded', &
    allocated(vessel%mesh%region) .and. size(vessel%mesh%region)==2, &
    'expected wall and opening regions')

  if(allocated(vessel%mesh%region) .and. size(vessel%mesh%region)==2) then
    call record_result('behavior resolved by physical name', &
      trim(vessel%mesh%region(1)%name)=='wall' .and. &
      trim(vessel%mesh%region(1)%behavior_type)=='wall' .and. &
      vessel%mesh%region(1)%is_active .and. &
      trim(vessel%mesh%region(2)%name)=='opening' .and. &
      trim(vessel%mesh%region(2)%behavior_type)=='opening' .and. &
      .not.vessel%mesh%region(2)%is_active, &
      'expected active wall followed by inactive opening')
    call record_result('absorbing wall parameters loaded', &
      abs(vessel%mesh%region(1)%wall%p_absorb-1.d0)<1.d-15 .and. &
      abs(vessel%mesh%region(1)%wall%p_specular)<1.d-15 .and. &
      abs(vessel%mesh%region(1)%wall%p_elastic_diffuse)<1.d-15, &
      'expected absorption one and omitted reflection probabilities zero')
  else
    call record_result('behavior resolved by physical name',.false., &
      'region catalogue is unavailable')
    call record_result('absorbing wall parameters loaded',.false., &
      'region catalogue is unavailable')
  endif

  ! Triangle 1 belongs to the active wall; triangle 2 belongs to the inactive
  ! opening. Count their reconstructed references to prove that setup retained
  ! active geometry and removed inactive geometry from every grid cell.
  active_references=0
  inactive_references=0
  do k=1,beam_grid%nz
    do j=1,beam_grid%ny
      do i=1,beam_grid%nx
        if(.not.allocated(vessel%map(i,j,k)%surface_id)) cycle
        do reference_index=1,size(vessel%map(i,j,k)%surface_id)
          select case(vessel%map(i,j,k)%surface_id(reference_index))
          case(1)
            active_references=active_references+1
          case(2)
            inactive_references=inactive_references+1
          end select
        enddo
      enddo
    enddo
  enddo
  write(detail,'(a,i0)') 'active triangle references: ',active_references
  call record_result('active triangle reconstructed', &
    active_references>0,detail)
  write(detail,'(a,i0)') 'inactive triangle references: ',inactive_references
  call record_result('inactive triangle filtered', &
    inactive_references==0,detail)

  ! Exercise voxel traversal, exact mesh intersection, behavior lookup,
  ! collision padding, and wall physics in one deterministic track. The start
  ! point lies above the interior of triangle 1 rather than the diagonal shared
  ! by the fixture's two triangles.
  allocate(tracks(beam_grid%ntrack))
  call track_to_wall( &
    [0.25d0,0.10d0,0.25d0], &
    [0.d0,0.d0,-1.d0], &
    tracks,ntrack,absorbed)
  write(detail,'(a,l1,a,i0)') 'absorbed: ',absorbed, &
    '; track segments: ',ntrack
  call record_result('complete track reaches absorbing wall', &
    absorbed .and. ntrack==1,detail)

  if(ntrack>=1) then
    write(detail,'(a,3(i0,1x),a,es12.4)') 'cell: ',tracks(1)%ind, &
      '; time: ',tracks(1)%time
    call record_result('track segment geometry is correct', &
      all(tracks(1)%ind==[1,1,1]) .and. &
      tracks(1)%time>0.d0 .and. tracks(1)%time<0.25d0,detail)
  else
    call record_result('track segment geometry is correct',.false., &
      'no track segment was returned')
  endif

  ! Periodicity changes only grid-face behavior. A physical wall intersected
  ! before that face must therefore retain priority.
  beam_grid%periodic_axes=[.false.,.false.,.true.]
  call track_to_wall( &
    [0.25d0,0.10d0,0.25d0], &
    [0.d0,0.d0,-1.d0], &
    tracks,ntrack,absorbed,exit_reason,periodic_wraps)
  call record_result('physical wall precedes periodic boundary', &
    absorbed .and. ntrack==1 .and. periodic_wraps==0, &
    'expected the absorbing wall before the periodic Z face')

  ! With wall interactions disabled, the unchanged default exits at the
  ! nonperiodic Z face and records no wraps.
  inputs%enable_wall_interactions=.false.
  beam_grid%periodic_axes=.false.
  call track_to_wall( &
    [0.25d0,0.10d0,0.25d0], &
    [0.d0,0.d0,1.d0], &
    tracks,ntrack,absorbed,exit_reason,periodic_wraps)
  call record_result('default grid boundary remains nonperiodic', &
    .not.absorbed .and. exit_reason==track_exit_grid .and. &
    periodic_wraps==0 .and. ntrack==2, &
    'expected two Z-cell segments followed by a grid exit')

  ! A Z-periodic ray wraps between the paired faces until the configured
  ! segment capacity is exhausted.
  beam_grid%periodic_axes=[.false.,.false.,.true.]
  beam_grid%ntrack=5
  call track_to_wall( &
    [0.25d0,0.10d0,0.25d0], &
    [0.d0,0.d0,1.d0], &
    tracks,ntrack,absorbed,exit_reason,periodic_wraps)
  call record_result('periodic Z boundary wraps the track', &
    .not.absorbed .and. exit_reason==track_exit_path_limit .and. &
    periodic_wraps==2 .and. ntrack==5 .and. &
    all([tracks(1)%ind(3),tracks(2)%ind(3),tracks(3)%ind(3), &
      tracks(4)%ind(3),tracks(5)%ind(3)]==[1,2,1,2,1]), &
    'expected two wraps and termination at the five-segment path limit')

  ! Exercise pure elastic-diffuse dispatch and its shared cosine sampler with
  ! an arbitrary normal. This direct native regression is the isolated planar
  ! validation: it avoids a separate simulation fixture while checking the
  ! physical sampling contract over a large deterministic ensemble.
  call rng_init(rng(1),24680_int32)
  normal=[1.d0,2.d0,3.d0]/sqrt(14.d0)
  tangent_x=[-2.d0,1.d0,0.d0]/sqrt(5.d0)
  tangent_y=cross_product(normal,tangent_x)
  incident_speed=7.5d0
  event%normal=normal
  elastic_wall%temp=1.d-3
  elastic_wall%p_absorb=0.d0
  elastic_wall%p_specular=0.d0
  elastic_wall%p_elastic_diffuse=1.d0
  elastic_wall%p_thermal=0.d0
  maximum_speed_error=0.d0
  minimum_mu=1.d0
  mean_mu=0.d0
  mean_azimuth_x=0.d0
  mean_azimuth_y=0.d0

  do sample_index=1,number_of_samples
    velocity=-incident_speed*normal
    call apply_wall_interaction(velocity,event,elastic_wall,absorbed)
    if(absorbed) error stop 'test_006: pure elastic-diffuse wall absorbed a particle'
    maximum_speed_error=max(maximum_speed_error,abs(norm2(velocity)-incident_speed))
    direction=velocity/norm2(velocity)
    mu=dot_product(direction,normal)
    minimum_mu=min(minimum_mu,mu)
    mean_mu=mean_mu+mu
    azimuth_x=dot_product(direction,tangent_x)
    azimuth_y=dot_product(direction,tangent_y)
    mean_azimuth_x=mean_azimuth_x+azimuth_x
    mean_azimuth_y=mean_azimuth_y+azimuth_y
  enddo
  mean_mu=mean_mu/number_of_samples
  mean_azimuth_x=mean_azimuth_x/number_of_samples
  mean_azimuth_y=mean_azimuth_y/number_of_samples

  write(detail,'(a,es12.4)') 'maximum speed error: ',maximum_speed_error
  call record_result('elastic diffuse preserves speed', &
    maximum_speed_error<1.d-12,detail)
  write(detail,'(a,es12.4)') 'minimum direction cosine: ',minimum_mu
  call record_result('elastic diffuse points inward',minimum_mu>=0.d0,detail)
  write(detail,'(a,f10.6)') 'mean direction cosine: ',mean_mu
  call record_result('elastic diffuse has cosine-law mean', &
    abs(mean_mu-2.d0/3.d0)<5.d-3,detail)
  write(detail,'(a,2(es12.4,1x))') 'mean tangent components: ', &
    mean_azimuth_x,mean_azimuth_y
  call record_result('elastic diffuse has unbiased azimuth', &
    abs(mean_azimuth_x)<5.d-3 .and. abs(mean_azimuth_y)<5.d-3,detail)

  call write_comparison_report(all_passed)
  if(.not.all_passed) error stop 'test_006: one or more checks failed'
  print *,'test_006: self-contained mesh acceptance passed'

contains

  subroutine initialize_test_grid()
    !+ Construct the two-by-two-by-two grid encoded in the test artifact.
    !+
    !+ The routine initializes only the BeamGrid fields consumed by artifact
    !+ validation, voxel reconstruction, and track_to_wall. Coordinates use
    !+ centimetres in the beam-grid-local frame.
    beam_grid%nx=2
    beam_grid%ny=2
    beam_grid%nz=2
    beam_grid%dims=[2,2,2]
    beam_grid%ngrid=8
    beam_grid%ntrack=300
    beam_grid%xmin=-0.5d0
    beam_grid%xmax=1.5d0
    beam_grid%ymin=-0.5d0
    beam_grid%ymax=1.5d0
    beam_grid%zmin=-0.5d0
    beam_grid%zmax=1.5d0
    beam_grid%origin=[0.d0,0.d0,0.d0]
    beam_grid%alpha=0.d0
    beam_grid%beta=0.d0
    beam_grid%gamma=0.d0
    beam_grid%basis=0.d0
    beam_grid%basis(1,1)=1.d0
    beam_grid%basis(2,2)=1.d0
    beam_grid%basis(3,3)=1.d0
    beam_grid%inv_basis=beam_grid%basis
    beam_grid%dr=[1.d0,1.d0,1.d0]
    beam_grid%ds=1.d0
    beam_grid%lwh=[2.d0,2.d0,2.d0]
    beam_grid%periodic_axes=.false.
    allocate(beam_grid%xc(2),beam_grid%yc(2),beam_grid%zc(2))
    beam_grid%xc=[0.d0,1.d0]
    beam_grid%yc=[0.d0,1.d0]
    beam_grid%zc=[0.d0,1.d0]
  end subroutine initialize_test_grid

  subroutine record_result(label,passed,description)
    !+ Store one named acceptance assertion for terminal and file reporting.
    character(len=*),intent(in)::label,description
    logical,intent(in)::passed

    if(number_of_results>=max_results) &
      error stop 'test_006: result capacity exceeded'
    number_of_results=number_of_results+1
    results(number_of_results)%label=label
    results(number_of_results)%passed=passed
    results(number_of_results)%detail=description
    if(passed) then
      write(*,'(a,a)') 'PASS: ',trim(label)
    else
      write(*,'(a,a,a)') 'FAIL: ',trim(label),' - '//trim(description)
    endif
  end subroutine record_result

  subroutine write_comparison_report(passed)
    !+ Write the persistent Test 006 result summary and every named check.
    logical,intent(out)::passed
    integer::report_unit,result_index,number_passed

    passed=all(results(1:number_of_results)%passed)
    number_passed=count(results(1:number_of_results)%passed)
    open(newunit=report_unit,file=report_path,status='replace', &
      action='write')
    if(passed) then
      write(report_unit,'(a)') 'OVERALL REGRESSION STATUS: PASS'
    else
      write(report_unit,'(a)') 'OVERALL REGRESSION STATUS: FAIL'
    endif
    write(report_unit,'(a,i0)') 'Checks: ',number_of_results
    write(report_unit,'(a,i0)') 'Passed: ',number_passed
    write(report_unit,'(a,i0)') 'Failed: ',number_of_results-number_passed
    write(report_unit,'(a)') ''
    do result_index=1,number_of_results
      if(results(result_index)%passed) then
        write(report_unit,'(a,a)') 'PASS: ', &
          trim(results(result_index)%label)
      else
        write(report_unit,'(a,a)') 'FAIL: ', &
          trim(results(result_index)%label)
      endif
      write(report_unit,'(a,a)') '  ', &
        trim(results(result_index)%detail)
    enddo
    close(report_unit)
  end subroutine write_comparison_report

end program test_006
