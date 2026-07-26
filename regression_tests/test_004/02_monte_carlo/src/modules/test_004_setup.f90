module test_004_setup
  use iso_fortran_env, &
    only: &
      Int32 => int32, &
      Float64 => real64, &
      output_unit
  use libfida, &
    only: &
      beam_grid, &
      inter_grid, &
      equil, &
      fbm, &
      inputs, &
      tables, &
      impurity_charge, &
      reservoir_size, &
      n_thermal, &
      thermal_mass, &
      read_tables, &
      make_beam_grid, &
      LocalEMFields, &
      LocalProfiles, &
      get_fields, &
      get_plasma, &
      get_distribution
  use test_004_types, &
    only: &
      MonteCarloConfig, &
      DistributionCase
  use utilities, &
    only: &
      rng, &
      rng_init
  implicit none
  private

  integer, parameter :: beam_cells_per_axis = 3
  integer, parameter :: radial_points = 6
  integer, parameter :: axial_points = 9
  integer, parameter :: toroidal_points = 1
  integer(Int32), parameter :: center_cell_index(3) = [2, 2, 2]
  real(Float64), parameter :: beam_half_extent_cm = 1.5_Float64
  real(Float64), parameter :: interpolation_spacing_cm = 0.5_Float64
  real(Float64), parameter :: interpolation_z_min_cm = -2.0_Float64
  real(Float64), parameter :: magnetic_field_z_tesla = 1.0_Float64
  real(Float64), parameter :: two_pi = &
    2.0_Float64 * acos(-1.0_Float64)
  integer, parameter :: carbon_charge_state = 6

  public :: &
    configure_fidasim, &
    initialize_serial_rng, &
    reset_serial_rng, &
    initialize_atomic_tables, &
    print_atomic_table_setup, &
    initialize_test_setup, &
    print_test_setup, &
    teardown_test_setup

contains

  subroutine configure_fidasim(config)
    !+ Translate the Test 004 configuration into FIDASIM global input configuration.
    type(MonteCarloConfig), intent(in) :: config
      !+ Validated, normalized configuration shared by all distribution cases.

    ! Set the FIDASIM input configuration to match the Test 004 requirements.
    inputs%tables_file = trim(config%tables_filename)
    inputs%full_f = 1 ! Distribution is a full distirbution not a correction to a thermal distribution.
    inputs%non_thermal_beam_stopping = 0 ! Non-thermal beam stopping is not used in this test.
    inputs%non_thermal_cx_sampling = 1 ! Non-thermal CX sampling is used in this test.
    inputs%flr = 0 ! Non-thermal FLR is not used in this test, approximates a very high B field case.
    inputs%dist_type = 1 ! Distribution is grid based not monte-carlo sampled.
    inputs%calc_neutron = 0 ! Neutron calculation is not used in this test.
    inputs%calc_cfpd = 0 ! CFPD calculation is not used in this test.
    inputs%verbose = 0 ! Verbose output is not used in this test.
    inputs%seed = config%seed
    inputs%reservoir_size = config%reservoir_size
    n_thermal = 1 ! Only one ion species is used in this test.
    reservoir_size = config%reservoir_size

    ! read_tables always loads one impurity-transition table. Carbon is used
    ! only to satisfy that production interface; the test impurity density is
    ! zero, so this table does not contribute to the ion-sink calculation.
    impurity_charge = carbon_charge_state
  end subroutine configure_fidasim

  subroutine initialize_serial_rng()
    !+ Initialize the single FIDASIM random-number stream from the test seed.

    if (.not. allocated(rng)) then
      allocate(rng(1))
    else if (size(rng) /= 1) then
      deallocate(rng)
      allocate(rng(1))
    end if

    call rng_init(rng(1), inputs%seed)
  end subroutine initialize_serial_rng

  subroutine reset_serial_rng()
    !+ Restore the configured seed before sampling ions for one case.
    call rng_init(rng(1), inputs%seed)
  end subroutine reset_serial_rng

  subroutine initialize_atomic_tables()
    !+ Load the production atomic data selected by configure_fidasim.
    call read_tables()
  end subroutine initialize_atomic_tables

  subroutine print_atomic_table_setup()
    !+ Print a concise summary confirming which atomic data were loaded.

    write(output_unit, '(a)') 'Atomic-table setup:'
    write(output_unit, '(a,a)') '  file: ', trim(inputs%tables_file)
    write(output_unit, '(a,i0)') '  atomic levels: ', &
      tables%H_H_cx_cross%n_max
    write(output_unit, '(a,i0)') '  H-H CX cross-section energies: ', &
      tables%H_H_cx_cross%nenergy
    write(output_unit, '(a,i0)') '  H-H CX rate energies: ', &
      tables%H_H_cx_rate%nenergy
    write(output_unit, '(a,i0)') '  H-H CX rate temperatures: ', &
      tables%H_H_cx_rate%ntemp
  end subroutine print_atomic_table_setup

  subroutine initialize_test_setup(distribution)
    !+ Populate the minimal FIDASIM global state used by the ion-sink test.
    type(DistributionCase), intent(in) :: distribution
      !+ Smooth Test 002 distribution and isotope metadata for the current case.

    ! Use the current distribution isotope for the configured ion species.
    thermal_mass = 0.0_Float64
    thermal_mass(1) = distribution%atomic_mass

    call setup_interpolation_grid()
    call setup_uniform_equilibrium(distribution)

    ! make_beam_grid queries the equilibrium mask while classifying its cells.
    call setup_beam_grid()
    call setup_uniform_fbm(distribution)
  end subroutine initialize_test_setup

  subroutine setup_interpolation_grid()
    !+ Construct a uniform axisymmetric interpolation grid around the beam grid.
    integer :: radial_index, axial_index

    ! R=[0,2.5] cm encloses the largest beam-grid radius,
    ! sqrt(1.5^2 + 1.5^2) cm, and Z=[-2,2] cm encloses its vertical extent.
    inter_grid%nr = radial_points
    inter_grid%nz = axial_points
    inter_grid%nphi = toroidal_points
    inter_grid%dims = [radial_points, axial_points, toroidal_points]
    inter_grid%dr = interpolation_spacing_cm
    inter_grid%dz = interpolation_spacing_cm
    ! One toroidal point represents the complete axisymmetric interval.
    inter_grid%dphi = two_pi
    inter_grid%da = inter_grid%dr * inter_grid%dz
    inter_grid%dv = inter_grid%da * inter_grid%dphi
    inter_grid%ngrid = product(inter_grid%dims)

    allocate(inter_grid%r(radial_points))
    allocate(inter_grid%z(axial_points))
    allocate(inter_grid%phi(toroidal_points))

    ! Set up the interpolation grid coordinates:
    do radial_index = 1, radial_points
      inter_grid%r(radial_index) = &
        real(radial_index - 1, Float64) * interpolation_spacing_cm
    end do
    do axial_index = 1, axial_points
      inter_grid%z(axial_index) = interpolation_z_min_cm + &
        real(axial_index - 1, Float64) * interpolation_spacing_cm
    end do
    inter_grid%phi = 0.0_Float64

  end subroutine setup_interpolation_grid

  subroutine setup_uniform_equilibrium(distribution)
    !+ Define uniform plasma profiles and fields over the interpolation grid.
    type(DistributionCase), intent(in) :: distribution
      !+ Distribution whose validated density defines the uniform profiles.

    allocate(equil%fields(radial_points, axial_points, toroidal_points))
    allocate(equil%plasma(radial_points, axial_points, toroidal_points))
    allocate(equil%mask(radial_points, axial_points, toroidal_points))

    ! The identity beam-grid basis makes cylindrical +Z the Cartesian lab +z.
    equil%mask = 1.0_Float64
    equil%fields%br = 0.0_Float64
    equil%fields%bt = 0.0_Float64
    equil%fields%bz = magnetic_field_z_tesla
    equil%fields%er = 0.0_Float64
    equil%fields%et = 0.0_Float64
    equil%fields%ez = 0.0_Float64

    ! Use the validated Test 002 density uniformly at every spatial point.
    ! Other profile components and field derivatives retain their zero defaults.
    equil%plasma%dene = distribution%denf
    equil%plasma%deni(1) = distribution%denf
    equil%plasma%denf = distribution%denf
    equil%plasma%zeff = 1.0_Float64
    equil%plasma%vr = 0.0_Float64
    equil%plasma%vt = 0.0_Float64
    equil%plasma%vz = 0.0_Float64
  end subroutine setup_uniform_equilibrium

  subroutine setup_beam_grid()
    !+ Construct the 3x3x3 Cartesian beam grid centered on the lab origin.

    beam_grid%nx = beam_cells_per_axis
    beam_grid%ny = beam_cells_per_axis
    beam_grid%nz = beam_cells_per_axis
    beam_grid%xmin = -beam_half_extent_cm
    beam_grid%xmax = beam_half_extent_cm
    beam_grid%ymin = -beam_half_extent_cm
    beam_grid%ymax = beam_half_extent_cm
    beam_grid%zmin = -beam_half_extent_cm
    beam_grid%zmax = beam_half_extent_cm
    beam_grid%alpha = 0.0_Float64
    beam_grid%beta = 0.0_Float64
    beam_grid%gamma = 0.0_Float64
    beam_grid%origin = 0.0_Float64

    call make_beam_grid()
  end subroutine setup_beam_grid

  subroutine setup_uniform_fbm(distribution)
    !+ Replicate one smooth F(E,p) distribution over the interpolation grid.
    type(DistributionCase), intent(in) :: distribution
      !+ Smooth single-location distribution replicated at every spatial point.
    integer :: radial_index, axial_index

    ! Test 002 guarantees non-singleton, uniformly spaced energy and pitch grids.
    fbm%A = distribution%atomic_mass
    fbm%nenergy = size(distribution%energy)
    fbm%npitch = size(distribution%pitch)
    fbm%nr = radial_points
    fbm%nz = axial_points
    fbm%nphi = toroidal_points
    fbm%dE = abs(distribution%energy(2) - distribution%energy(1))
    fbm%dp = abs(distribution%pitch(2) - distribution%pitch(1))
    fbm%dr = inter_grid%dr
    fbm%dz = inter_grid%dz
    fbm%dphi = inter_grid%dphi
    fbm%emin = minval(distribution%energy)
    fbm%emax = maxval(distribution%energy)
    fbm%e_range = fbm%emax - fbm%emin
    fbm%pmin = minval(distribution%pitch)
    fbm%pmax = maxval(distribution%pitch)
    fbm%p_range = fbm%pmax - fbm%pmin
    fbm%rmin = minval(inter_grid%r)
    fbm%rmax = maxval(inter_grid%r)
    fbm%r_range = fbm%rmax - fbm%rmin
    fbm%zmin = minval(inter_grid%z)
    fbm%zmax = maxval(inter_grid%z)
    fbm%z_range = fbm%zmax - fbm%zmin
    fbm%phimin = minval(inter_grid%phi)
    fbm%phimax = maxval(inter_grid%phi)
    fbm%phi_range = fbm%phimax - fbm%phimin

    allocate(fbm%energy(fbm%nenergy), fbm%pitch(fbm%npitch))
    allocate(fbm%r(fbm%nr), fbm%z(fbm%nz), fbm%phi(fbm%nphi))
    allocate(fbm%denf(fbm%nr, fbm%nz, fbm%nphi))
    allocate(fbm%f(fbm%nenergy, fbm%npitch, fbm%nr, fbm%nz, fbm%nphi))

    fbm%energy = distribution%energy
    fbm%pitch = distribution%pitch
    fbm%r = inter_grid%r
    fbm%z = inter_grid%z
    fbm%phi = inter_grid%phi
    fbm%denf = distribution%denf

    ! The source distribution's selected R and Z identify its provenance but do
    ! not constrain this uniform test setup. Copy the same velocity-space
    ! distribution to every R-Z point of the axisymmetric grid.
    do axial_index = 1, fbm%nz
      do radial_index = 1, fbm%nr
        fbm%f(:,:,radial_index,axial_index,1) = distribution%f
      end do
    end do
  end subroutine setup_uniform_fbm

  subroutine print_test_setup()
    !+ Print center-cell values obtained through FIDASIM interpolation routines.
    type(LocalEMFields) :: local_fields
    type(LocalProfiles) :: local_plasma
    real(Float64), allocatable :: local_distribution(:,:)
    real(Float64) :: local_denf

    allocate(local_distribution(fbm%nenergy, fbm%npitch))

    ! Get the test setup's center-cell fields, profiles, and distribution:
    call get_fields(local_fields, ind=center_cell_index)
    call get_plasma(local_plasma, ind=center_cell_index)
    call get_distribution( &
      local_distribution, local_denf, ind=center_cell_index)

    ! Print the center-cell values and other test-setup metadata:
    write(output_unit, '(a,3(1x,i0))') '  beam-grid dimensions:', &
      beam_grid%dims
    write(output_unit, '(a,3(1x,es14.6))') '  beam-grid spacing [cm]:', &
      beam_grid%dr
    write(output_unit, '(a,es14.6)') '  beam-cell volume [cm^3]: ', &
      beam_grid%dv
    write(output_unit, '(a,i0)') '  beam cells in plasma: ', &
      count(beam_grid%in_plasma == 1)
    write(output_unit, '(a,3(1x,i0))') '  interpolation-grid dimensions:', &
      inter_grid%dims
    write(output_unit, '(a,2(1x,es14.6))') '  interpolation R range [cm]:', &
      minval(inter_grid%r), maxval(inter_grid%r)
    write(output_unit, '(a,2(1x,es14.6))') '  interpolation Z range [cm]:', &
      minval(inter_grid%z), maxval(inter_grid%z)
    write(output_unit, '(a,3(1x,es14.6))') '  center B [T]:', &
      local_fields%br, local_fields%bt, local_fields%bz
    write(output_unit, '(a,es14.6)') '  center species-1 density [cm^-3]: ', &
      local_plasma%deni(1)
    write(output_unit, '(a,es14.6)') '  center FBM density [cm^-3]: ', &
      local_denf
    write(output_unit, '(a,es14.6)') '  center FBM maximum: ', &
      maxval(local_distribution)
  end subroutine print_test_setup

  subroutine teardown_test_setup()
    if (allocated(fbm%energy)) deallocate(fbm%energy)
    if (allocated(fbm%pitch)) deallocate(fbm%pitch)
    if (allocated(fbm%r)) deallocate(fbm%r)
    if (allocated(fbm%z)) deallocate(fbm%z)
    if (allocated(fbm%phi)) deallocate(fbm%phi)
    if (allocated(fbm%denf)) deallocate(fbm%denf)
    if (allocated(fbm%f)) deallocate(fbm%f)

    if (allocated(beam_grid%xc)) deallocate(beam_grid%xc)
    if (allocated(beam_grid%yc)) deallocate(beam_grid%yc)
    if (allocated(beam_grid%zc)) deallocate(beam_grid%zc)
    if (allocated(beam_grid%xe)) deallocate(beam_grid%xe)
    if (allocated(beam_grid%ye)) deallocate(beam_grid%ye)
    if (allocated(beam_grid%ze)) deallocate(beam_grid%ze)
    if (allocated(beam_grid%in_plasma)) deallocate(beam_grid%in_plasma)

    if (allocated(equil%fields)) deallocate(equil%fields)
    if (allocated(equil%plasma)) deallocate(equil%plasma)
    if (allocated(equil%mask)) deallocate(equil%mask)

    if (allocated(inter_grid%r)) deallocate(inter_grid%r)
    if (allocated(inter_grid%z)) deallocate(inter_grid%z)
    if (allocated(inter_grid%phi)) deallocate(inter_grid%phi)
  end subroutine teardown_test_setup

end module test_004_setup
