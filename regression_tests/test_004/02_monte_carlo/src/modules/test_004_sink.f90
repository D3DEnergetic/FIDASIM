module test_004_sink
  !+ Test how FIDASIM records ions removed by charge exchange.
  !+
  !+ The test uses a small artificial grid and arranges for every sampled ion to
  !+ be removed in the middle cell. Each sampled ion has a weight and kinetic
  !+ energy. FIDASIM should add its weighted particle-removal rate and removed
  !+ kinetic power to that middle cell, record one contribution, and save the
  !+ particle in the sink list.
  !+
  !+ After all ions have been sampled, the test adds the saved particle weights
  !+ and energies independently. Those sums must reproduce the values stored on
  !+ the grid. The test also checks that no removal was recorded in another
  !+ cell, that the squared-contribution tallies were filled, that the normal
  !+ sink HDF5 file can be written, and that the temporary sink arrays are freed.
  use iso_fortran_env, &
    only: &
      Int32 => int32, &
      Int64 => int64, &
      Float64 => real64, &
      output_unit
  use libfida, &
    only: &
      beam_grid, &
      inputs, &
      sink, &
      nbif_type, &
      get_total_cx_rate, &
      store_sinks, &
      store_sink_particle, &
      write_sink_profile
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
  real(Float64), parameter :: joules_per_kev = 1.60217733e-16_Float64

  public :: calculate_ion_sink, print_ion_sink, finalize_ion_sink

contains

  subroutine calculate_ion_sink(config, distribution)
    !+ Sample ions and accumulate their charge-exchange sink contributions.
    !+
    !+ For each requested marker, sample the prescribed non-thermal ion
    !+ distribution with the production FIDASIM sampler, validate the sampled
    !+ phase-space point, and evaluate the six-level charge-exchange rate against
    !+ the central-cell neutral reservoir. Convert that rate to one marker's
    !+ particle-rate-density contribution, then call the production storage
    !+ routines with the sampled energy and particle metadata.
    !+
    !+ The production routines accumulate particle-rate density, kinetic-power
    !+ density, contribution count, and squared rate and power contributions in
    !+ the central grid cell. They also retain the weighted sink particle used by
    !+ the independent checks in [[print_ion_sink]].
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
    allocate(sink%power_dens( &
      1, beam_grid%nx, beam_grid%ny, beam_grid%nz))
    allocate(sink%contribution_count( &
      1, beam_grid%nx, beam_grid%ny, beam_grid%nz))
    allocate(sink%rate_squared_sum( &
      1, beam_grid%nx, beam_grid%ny, beam_grid%nz))
    allocate(sink%power_squared_sum( &
      1, beam_grid%nx, beam_grid%ny, beam_grid%nz))
    sink%cnt = 1
    sink%dens = 0.0_Float64
    sink%power_dens = 0.0_Float64
    sink%contribution_count = 0.0_Float64
    sink%rate_squared_sum = 0.0_Float64
    sink%power_squared_sum = 0.0_Float64

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
        sample%particle_cell_index, ion_species_index, sink_rate_per_marker, &
        sample%energy)
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
    !+ Verify and summarize the completed central-cell ion sink.
    !+
    !+ This routine checks that one sink particle was stored for every requested
    !+ marker and that all particle-rate density is confined to the central cell.
    !+ It reconstructs the central-cell rate from the stored particle weights and
    !+ independently reconstructs the kinetic-power density from particle weights
    !+ and energies. Both reconstructed values must agree with the gridded
    !+ production moments within [[maximum_relative_rate_error]].
    !+
    !+ Because this stripped-down test makes exactly one sink contribution per
    !+ marker, the central contribution count must equal `config%n_markers`. The
    !+ rate-squared and power-squared sums must also be positive. These last checks
    !+ establish that the sampling moments are populated; they do not validate a
    !+ derived uncertainty estimator.
    type(MonteCarloConfig), intent(in) :: config
      !+ Configuration defining the expected number of stored particles.
    integer(Int64) :: stored_particles
    real(Float64) :: central_rate, particle_rate, off_center_rate
    real(Float64) :: central_power, particle_power
    real(Float64) :: relative_rate_error, relative_power_error
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
    central_power = sink%power_dens(1,2,2,2)
    particle_power = sum( &
      sink%part(1:stored_particles)%weight * &
      sink%part(1:stored_particles)%energy) * joules_per_kev / beam_grid%dv

    allocate(off_center_density, source=sink%dens)
    off_center_density(1,2,2,2) = 0.0_Float64
    off_center_rate = sum(abs(off_center_density))
    relative_rate_error = abs(particle_rate - central_rate) / central_rate
    relative_power_error = abs(particle_power - central_power) / central_power

    if (off_center_rate > 0.0_Float64) then
      error stop 'Ion-sink density was stored outside the central cell'
    end if
    if (relative_rate_error > maximum_relative_rate_error) then
      error stop 'Particle weights do not reproduce the ion-sink density'
    end if
    if (central_power <= 0.0_Float64 .or. &
        relative_power_error > maximum_relative_rate_error) then
      error stop 'Particle energies do not reproduce the ion-sink power density'
    end if
    if (abs(sink%contribution_count(1,2,2,2) - &
        real(config%n_markers, Float64)) > 0.5_Float64) then
      error stop 'The ion-sink contribution count is incorrect'
    end if
    if (sink%rate_squared_sum(1,2,2,2) <= 0.0_Float64 .or. &
        sink%power_squared_sum(1,2,2,2) <= 0.0_Float64) then
      error stop 'The ion-sink squared sums are not positive'
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

  subroutine finalize_ion_sink(config)
    !+ Complete production output handling and release all sink storage.
    !+
    !+ When output is enabled, call [[write_sink_profile]], require the normal
    !+ FIDASIM sink HDF5 file to exist, and verify that the writer reset the sink
    !+ counter and deallocated the particle, rate, power, count, and squared-sum
    !+ arrays. When output is disabled, perform the same cleanup locally. This
    !+ checks the production writer lifecycle in addition to the numerical sink
    !+ accumulation.
    type(MonteCarloConfig), intent(in) :: config
      !+ Configuration selecting whether the production output is saved.
    character(len=4096) :: output_filename
    logical :: output_exists

    if (config%save_data) then
      output_filename = trim(inputs%result_dir)//'/'// &
        trim(inputs%runid)//'_sink.h5'

      ! The production writer resets sink%cnt and releases all sink storage.
      call write_sink_profile()

      inquire(file=trim(output_filename), exist=output_exists)
      if (.not. output_exists) then
        error stop 'The production ion-sink file was not created'
      end if
      if (sink%cnt /= 1 .or. allocated(sink%part) .or. &
          allocated(sink%dens) .or. allocated(sink%power_dens) .or. &
          allocated(sink%contribution_count) .or. &
          allocated(sink%rate_squared_sum) .or. &
          allocated(sink%power_squared_sum)) then
        error stop 'The production sink writer did not release sink storage'
      end if

      write(output_unit, '(a,a)') '    production sink file: ', &
        trim(output_filename)
    else
      sink%cnt = 1
      if (allocated(sink%part)) deallocate(sink%part)
      if (allocated(sink%dens)) deallocate(sink%dens)
      if (allocated(sink%power_dens)) deallocate(sink%power_dens)
      if (allocated(sink%contribution_count)) deallocate(sink%contribution_count)
      if (allocated(sink%rate_squared_sum)) deallocate(sink%rate_squared_sum)
      if (allocated(sink%power_squared_sum)) deallocate(sink%power_squared_sum)
    end if
  end subroutine finalize_ion_sink

end module test_004_sink
