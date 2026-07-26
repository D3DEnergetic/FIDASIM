module test_004_sampling
  use iso_fortran_env, &
    only: &
      Int32 => int32, &
      Float64 => real64, &
      output_unit
  use libfida, &
    only: &
      LocalEMFields, &
      mc_sample_ion_f4d_gc, &
      v2_to_E_per_amu
  use test_004_types, &
    only: &
      DistributionCase
  implicit none
  private

  integer(Int32), parameter :: center_cell_index(3) = [2, 2, 2]
  integer, parameter :: ion_species_index = 1
  real(Float64), parameter :: maximum_relative_density_error = 1.0e-12_Float64
  real(Float64), parameter :: maximum_position_difference_cm = 1.0e-12_Float64

  type, public :: IonSample
    real(Float64) :: guiding_center_position(3) = 0.0_Float64
    real(Float64) :: particle_position(3) = 0.0_Float64
    real(Float64) :: velocity(3) = 0.0_Float64
    integer(Int32) :: particle_cell_index(3) = 0
    real(Float64) :: ion_density = 0.0_Float64
    real(Float64) :: distribution_value = 0.0_Float64
    real(Float64) :: energy = 0.0_Float64
    real(Float64) :: pitch = 0.0_Float64
    logical :: nonthermal_distribution_used = .false.
    type(LocalEMFields) :: fields
  end type IonSample

  public :: sample_ion, validate_ion_sample, print_ion_sample

contains

  subroutine sample_ion(distribution, sample)
    !+ Sample one central-cell ion through the production FIDASIM routine.
    type(DistributionCase), intent(in) :: distribution
      !+ Current distribution, used here for its isotope mass.
    type(IonSample), intent(out) :: sample
      !+ Sampled ion position, velocity, fields, density, energy, and pitch.

    call mc_sample_ion_f4d_gc( &
      center_cell_index, &
      ion_species_index, &
      sample%guiding_center_position, &
      sample%ion_density, &
      sample%fields, &
      sample%particle_position, &
      sample%velocity, &
      sample%particle_cell_index, &
      sample%distribution_value, &
      sample%nonthermal_distribution_used)

    sample%energy = dot_product(sample%velocity, sample%velocity) * &
      v2_to_E_per_amu * distribution%atomic_mass
    sample%pitch = dot_product( &
      sample%fields%b_norm, sample%velocity / norm2(sample%velocity))
  end subroutine sample_ion

  subroutine validate_ion_sample(distribution, sample)
    !+ Require the guarantees expected from the Test 004 sampling setup.
    type(DistributionCase), intent(in) :: distribution
      !+ Distribution defining the expected density and velocity-space bounds.
    type(IonSample), intent(in) :: sample
      !+ Sample returned by sample_ion.
    real(Float64) :: relative_density_error, position_difference
    real(Float64) :: energy_half_spacing, pitch_half_spacing

    relative_density_error = abs( &
      sample%ion_density - distribution%denf) / distribution%denf
    position_difference = maxval(abs( &
      sample%particle_position - sample%guiding_center_position))

    ! The sampler draws continuously within the selected energy-pitch bin.
    ! The valid sampling domain therefore extends half a grid spacing beyond
    ! the first and last bin centers.
    energy_half_spacing = 0.5_Float64 * &
      abs(distribution%energy(2) - distribution%energy(1))
    pitch_half_spacing = 0.5_Float64 * &
      abs(distribution%pitch(2) - distribution%pitch(1))

    if (.not. sample%nonthermal_distribution_used) then
      error stop 'Ion sampling did not use the nonthermal distribution'
    end if
    if (norm2(sample%velocity) <= 0.0_Float64) then
      error stop 'Ion sampling returned a zero velocity'
    end if
    if (any(sample%particle_cell_index /= center_cell_index)) then
      error stop 'The sampled ion is outside the central beam-grid cell'
    end if
    if (position_difference > maximum_position_difference_cm) then
      error stop 'Particle and guiding-center positions differ with flr=0'
    end if
    if (relative_density_error > maximum_relative_density_error) then
      error stop 'The sampled ion density does not match denf'
    end if
    if (sample%energy < minval(distribution%energy) - energy_half_spacing .or. &
        sample%energy > maxval(distribution%energy) + energy_half_spacing) then
      error stop 'The sampled ion energy is outside the distribution grid'
    end if
    if (sample%pitch < minval(distribution%pitch) - pitch_half_spacing .or. &
        sample%pitch > maxval(distribution%pitch) + pitch_half_spacing) then
      error stop 'The sampled ion pitch is outside the distribution grid'
    end if
  end subroutine validate_ion_sample

  subroutine print_ion_sample(sample)
    !+ Print the validated sample used to exercise the production interface.
    type(IonSample), intent(in) :: sample

    write(output_unit, '(a)') '  ion-sampling interface check:'
    write(output_unit, '(a,es14.6)') '    energy [keV]: ', sample%energy
    write(output_unit, '(a,es14.6)') '    pitch: ', sample%pitch
    write(output_unit, '(a,es14.6)') &
      '    sampled ion density [cm^-3]: ', sample%ion_density
    write(output_unit, '(a,3(1x,es14.6))') &
      '    guiding-center position [cm]:', sample%guiding_center_position
    write(output_unit, '(a,3(1x,es14.6))') &
      '    velocity [cm/s]:', sample%velocity
    write(output_unit, '(a,3(1x,i0))') &
      '    particle cell index:', sample%particle_cell_index
  end subroutine print_ion_sample

end module test_004_sampling
