program test_006
  !+ Verify FIDASIM's complete surface-mesh setup and tracking path.
  !+
  !+ The test owns a small HDF5 mesh artifact and a run-specific behavior
  !+ namelist under ``fixtures/``. It loads both through the production setup
  !+ routine, inspects the reconstructed runtime state, and tracks one
  !+ particle into a deterministic absorbing wall.

  use iso_fortran_env,only:int32,real64
  use libfida,only:ParticleTrack,beam_grid,inputs,setup_vessel_mesh, &
    track_to_wall,track_exit_grid,track_exit_path_limit,vessel
  use utilities,only:rng,rng_init
  implicit none

  integer,parameter::max_results=24
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
  integer::number_of_results=0,active_references,inactive_references
  integer(int32)::ntrack
  logical::absorbed,all_passed

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
      abs(vessel%mesh%region(1)%wall%p_specular)<1.d-15, &
      'expected deterministic absorption probability one')
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
