module test_004_sink
  use iso_fortran_env, &
    only: &
      Int32 => int32, &
      Int64 => int64, &
      Float64 => real64, &
      output_unit
  use libfida, &
    only: &
      beam_grid, &
      sink, &
      nbif_type, &
      get_total_cx_rate, &
      store_sinks, &
      store_sink_particle
  use test_004_types, &
    only: &
      MonteCarloConfig, &
      DistributionCase, &
      number_of_atomic_levels
  use test_004_sampling, &
    only: &
      IonSample, &
      sample_ion, &
      validate_ion_sample, &
      print_ion_sample
  implicit none
  private

  integer(Int32), parameter :: ion_species_index = 1
  integer(Int32), parameter :: neutral_types(1) = [nbif_type]
  real(Float64), parameter :: maximum_relative_rate_error = 1.0e-12_Float64

  public :: calculate_ion_sink, print_ion_sink, release_ion_sink

contains

  subroutine calculate_ion_sink(config, distribution)
    !+ Sample ions and store their type-1 charge-exchange sink contributions.
    type(MonteCarloConfig), intent(in) :: config
      !+ Configuration defining the number of Monte Carlo markers.
    type(DistributionCase), intent(in) :: distribution
      !+ Distribution sampled by the production FIDASIM routine.
    type(IonSample) :: sample
    real(Float64) :: rates(number_of_atomic_levels)
    real(Float64) :: sink_rate_per_marker
    real(Float64) :: distribution_value_per_marker
    integer(Int64) :: marker_index

    allocate(sink%part(config%n_markers))
    allocate(sink%dens( &
      1, beam_grid%nx, beam_grid%ny, beam_grid%nz))
    sink%cnt = 1
    sink%dens = 0.0_Float64

    do marker_index = 1_Int64, config%n_markers
      call sample_ion(distribution, sample)
      call validate_ion_sample(distribution, sample)

      call get_total_cx_rate( &
        sample%particle_cell_index, &
        sample%particle_position, &
        sample%velocity, &
        neutral_types, &
        rates)

      if (sum(rates) <= 0.0_Float64) then
        error stop 'The sampled ion has a non-positive charge-exchange rate'
      end if

      sink_rate_per_marker = sample%ion_density * sum(rates) / &
        real(config%n_markers, Float64)
      distribution_value_per_marker = sample%distribution_value / &
        real(config%n_markers, Float64)

      call store_sinks( &
        sample%particle_cell_index, ion_species_index, sink_rate_per_marker)
      call store_sink_particle( &
        sample%particle_cell_index, &
        sample%particle_position, &
        sample%velocity, &
        ion_species_index, &
        sink_rate_per_marker, &
        sample%distribution_value, &
        distribution_value_per_marker, &
        sample%fields)

      if (marker_index == 1_Int64) call print_ion_sample(sample)
    end do
  end subroutine calculate_ion_sink

  subroutine print_ion_sink(config)
    !+ Verify and print the completed central-cell ion sink.
    type(MonteCarloConfig), intent(in) :: config
      !+ Configuration defining the expected number of stored particles.
    integer(Int64) :: stored_particles
    real(Float64) :: central_rate, particle_rate, off_center_rate
    real(Float64) :: relative_rate_error
    real(Float64), allocatable :: off_center_density(:,:,:,:)

    stored_particles = int(sink%cnt - 1, Int64)
    if (stored_particles /= config%n_markers) then
      error stop 'The number of stored ion-sink particles is incorrect'
    end if

    central_rate = sink%dens(1,2,2,2)
    if (central_rate <= 0.0_Float64) then
      error stop 'The central-cell ion-sink rate is not positive'
    end if

    particle_rate = sum(sink%part(1:stored_particles)%weight) / beam_grid%dv

    allocate(off_center_density, source=sink%dens)
    off_center_density(1,2,2,2) = 0.0_Float64
    off_center_rate = sum(abs(off_center_density))
    relative_rate_error = abs(particle_rate - central_rate) / central_rate

    if (off_center_rate > 0.0_Float64) then
      error stop 'Ion-sink density was stored outside the central cell'
    end if
    if (relative_rate_error > maximum_relative_rate_error) then
      error stop 'Particle weights do not reproduce the ion-sink density'
    end if

    write(output_unit, '(a)') '  ion-sink calculation:'
    write(output_unit, '(a,i0)') '    stored particles: ', stored_particles
    write(output_unit, '(a,es14.6)') &
      '    central-cell rate [ions/(s cm^3)]: ', central_rate
    write(output_unit, '(a,es14.6)') &
      '    particle-weight rate [ions/(s cm^3)]: ', particle_rate
    write(output_unit, '(a,es14.6)') &
      '    relative rate difference: ', relative_rate_error
  end subroutine print_ion_sink

  subroutine release_ion_sink()
    !+ Release sink storage when no production HDF5 writer is called.
    sink%cnt = 1
    if (allocated(sink%part)) deallocate(sink%part)
    if (allocated(sink%dens)) deallocate(sink%dens)
  end subroutine release_ion_sink

end module test_004_sink
