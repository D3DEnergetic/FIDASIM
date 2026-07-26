module test_004_neutral
  use iso_fortran_env, &
    only: &
      Int32 => int32, &
      Float64 => real64, &
      output_unit
  use libfida, &
    only: &
      neut, &
      init_neutral_population, &
      update_neutrals, &
      free_neutral_population
  use test_004_types, only: MonteCarloConfig, DistributionCase, &
    NeutralParameters, number_of_atomic_levels
  implicit none
  private

  real(Float64), parameter :: atomic_mass_unit = 1.660539040e-27_Float64
  real(Float64), parameter :: elementary_charge = 1.60217733e-19_Float64
  real(Float64), parameter :: v2_to_energy_per_amu = atomic_mass_unit / &
    (2.0_Float64 * elementary_charge * 1.0e3_Float64) * 1.0e-4_Float64
  integer(Int32), parameter :: center_cell_index(3) = [2, 2, 2]

  public :: &
    build_neutral_parameters, &
    print_neutral_parameters, &
    initialize_neutral_population, &
    print_neutral_population, &
    release_neutral_population

contains

  subroutine build_neutral_parameters(config, distribution, neutral)
    type(MonteCarloConfig), intent(in) :: config
    type(DistributionCase), intent(in) :: distribution
    type(NeutralParameters), intent(out) :: neutral

    real(Float64) :: angle_radians
    real(Float64) :: level_fractions(number_of_atomic_levels)
    integer :: level_index

    neutral%atomic_mass = distribution%atomic_mass
    neutral%speed = sqrt(config%neutral_energy / &
      (v2_to_energy_per_amu * neutral%atomic_mass))

    angle_radians = config%injection_angle * &
      acos(-1.0_Float64) / 180.0_Float64
    neutral%velocity(1) = neutral%speed * sin(angle_radians)
    neutral%velocity(2) = 0.0_Float64
    neutral%velocity(3) = neutral%speed * cos(angle_radians)

    level_fractions = 0.0_Float64
    if (trim(config%level_split_method) == 'ground-only') then
      level_fractions(1) = 1.0_Float64
    else
      do level_index = 1, number_of_atomic_levels
        level_fractions(level_index) = exp( &
          -config%level_decay * real(level_index - 1, Float64))
      end do
      level_fractions = level_fractions / sum(level_fractions)
    end if

    neutral%level_density = config%neutral_density * level_fractions
  end subroutine build_neutral_parameters

  subroutine print_neutral_parameters(neutral)
    type(NeutralParameters), intent(in) :: neutral
    integer :: level_index

    write(output_unit, '(a,es14.6)') '  neutral mass [amu]: ', &
      neutral%atomic_mass
    write(output_unit, '(a,es14.6)') '  neutral speed [cm/s]: ', neutral%speed
    write(output_unit, '(a,3(1x,es14.6))') '  neutral velocity [cm/s]:', &
      neutral%velocity
    do level_index = 1, number_of_atomic_levels
      write(output_unit, '(a,i0,a,es14.6)') '  level_density(', &
        level_index, ') [cm^-3]: ', neutral%level_density(level_index)
    end do
    write(output_unit, '(a,es14.6)') '  total neutral density [cm^-3]: ', &
      sum(neutral%level_density)
  end subroutine print_neutral_parameters

  subroutine initialize_neutral_population(config, neutral)
    !+ Populate the central-cell type-1 neutral density and reservoir.
    type(MonteCarloConfig), intent(in) :: config
      !+ Configuration defining the number of neutral reservoir markers.
    type(NeutralParameters), intent(in) :: neutral
      !+ Case-specific neutral velocity and six-level density.
    real(Float64) :: density_per_marker(number_of_atomic_levels)
    integer :: marker_index

    call init_neutral_population(neut%full)

    density_per_marker = neutral%level_density / &
      real(config%reservoir_size, Float64)

    ! Supplying exactly reservoir_size particles fills every allocated entry.
    do marker_index = 1, config%reservoir_size
      call update_neutrals( &
        neut%full, center_cell_index, neutral%velocity, density_per_marker)
    end do
  end subroutine initialize_neutral_population

  subroutine print_neutral_population(config, neutral)
    !+ Verify and print the populated central-cell neutral reservoir.
    type(MonteCarloConfig), intent(in) :: config
      !+ Configuration defining the expected reservoir size.
    type(NeutralParameters), intent(in) :: neutral
      !+ Neutral values expected in every stored marker.
    integer :: marker_index, stored_markers
    real(Float64) :: center_density(number_of_atomic_levels)
    real(Float64) :: center_total, off_center_total, stored_weight
    real(Float64) :: maximum_velocity_difference
    real(Float64) :: tolerance
    real(Float64), allocatable :: off_center_density(:,:,:,:)

    stored_markers = min( &
      neut%full%res(2,2,2)%n, neut%full%res(2,2,2)%k)
    center_density = neut%full%dens(:,2,2,2)
    center_total = sum(center_density)
    allocate(off_center_density, source=neut%full%dens)
    off_center_density(:,2,2,2) = 0.0_Float64
    off_center_total = sum(abs(off_center_density))
    stored_weight = sum( &
      neut%full%res(2,2,2)%R(1:stored_markers)%w)
    maximum_velocity_difference = 0.0_Float64
    do marker_index = 1, stored_markers
      maximum_velocity_difference = max( &
        maximum_velocity_difference, &
        maxval(abs( &
          neut%full%res(2,2,2)%R(marker_index)%v - neutral%velocity)))
    end do
    tolerance = 100.0_Float64 * epsilon(1.0_Float64) * &
      max(1.0_Float64, sum(neutral%level_density))

    if (stored_markers /= config%reservoir_size) then
      error stop 'The central-cell neutral reservoir is not full'
    end if
    if (maxval(abs(center_density - neutral%level_density)) > tolerance) then
      error stop 'The central-cell neutral level densities are incorrect'
    end if
    if (abs(off_center_total) > tolerance) then
      error stop 'Neutral density was populated outside the central cell'
    end if
    if (abs(stored_weight - center_total) > tolerance) then
      error stop 'Neutral reservoir weights do not sum to the density'
    end if
    if (maximum_velocity_difference > 0.0_Float64) then
      error stop 'Neutral reservoir velocities are inconsistent'
    end if

    write(output_unit, '(a)') '  type-1 neutral population:'
    write(output_unit, '(a,i0)') '    stored central-cell markers: ', &
      stored_markers
    write(output_unit, '(a,es14.6)') &
      '    central-cell density [cm^-3]: ', center_total
    write(output_unit, '(a,es14.6)') &
      '    summed reservoir weight [cm^-3]: ', stored_weight
    write(output_unit, '(a,es14.6)') &
      '    density outside central cell [cm^-3]: ', off_center_total
  end subroutine print_neutral_population

  subroutine release_neutral_population()
    !+ Release the case-specific type-1 neutral population.
    call free_neutral_population(neut%full)
  end subroutine release_neutral_population

end module test_004_neutral
