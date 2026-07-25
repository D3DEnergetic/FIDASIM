module test_004_neutral
  use iso_fortran_env, only: Float64 => real64, output_unit
  use test_004_types, only: MonteCarloConfig, DistributionCase, &
    NeutralParameters, number_of_atomic_levels
  implicit none
  private

  real(Float64), parameter :: atomic_mass_unit = 1.660539040e-27_Float64
  real(Float64), parameter :: elementary_charge = 1.60217733e-19_Float64
  real(Float64), parameter :: v2_to_energy_per_amu = atomic_mass_unit / &
    (2.0_Float64 * elementary_charge * 1.0e3_Float64) * 1.0e-4_Float64

  public :: build_neutral_parameters, print_neutral_parameters

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

end module test_004_neutral
