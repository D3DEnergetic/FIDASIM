module test_002_sampling
  use iso_fortran_env, only: Int32 => int32, Int64 => int64, &
    Float64 => real64, output_unit
  use utilities, only: rng, rng_init, randind, randu
  use test_002_hdf5, only: DistributionFunction2D
  implicit none
  private

  public :: initialize_serial_rng, sample_distribution, print_sampling_summary

contains

  subroutine initialize_serial_rng(seed)
    integer(Int32), intent(in) :: seed

    if (.not. allocated(rng)) then
      allocate(rng(1))
    else if (size(rng) /= 1) then
      deallocate(rng)
      allocate(rng(1))
    end if
    call rng_init(rng(1), seed)
  end subroutine initialize_serial_rng

  subroutine sample_distribution(distribution, number_of_samples, counts, &
      sampled_f_array)
    type(DistributionFunction2D), intent(in) :: distribution
    integer(Int64), intent(in) :: number_of_samples
    integer(Int64), allocatable, intent(out) :: counts(:,:)
    real(Float64), allocatable, intent(out) :: sampled_f_array(:,:)

    integer, dimension(2,1) :: ep_ind
    real(Float64), dimension(3) :: randomu3
    real(Float64) :: sampled_energy, sampled_pitch
    integer(Int64) :: i
    integer :: energy_bin, pitch_bin, nenergy, npitch

    nenergy = size(distribution%energy_grid)
    npitch = size(distribution%pitch_grid)
    allocate(counts(nenergy, npitch), source=0_Int64)

    do i = 1, number_of_samples
      call randind(distribution%f_array, ep_ind)
      call randu(randomu3)

      sampled_energy = distribution%energy_grid(ep_ind(1,1)) &
        + distribution%denergy * (randomu3(1) - 0.5_Float64)
      sampled_pitch = distribution%pitch_grid(ep_ind(2,1)) &
        + distribution%dpitch * (randomu3(2) - 0.5_Float64)

      energy_bin = histogram_bin(sampled_energy, distribution%energy_grid(1), &
        distribution%denergy, nenergy)
      pitch_bin = histogram_bin(sampled_pitch, distribution%pitch_grid(1), &
        distribution%dpitch, npitch)
      if (energy_bin < 1 .or. energy_bin > nenergy .or. &
          pitch_bin < 1 .or. pitch_bin > npitch) then
        write(*, '(a,2(1x,es14.6))') 'Sample outside histogram domain:', &
          sampled_energy, sampled_pitch
        error stop 'Unable to assign sampled energy and pitch to histogram bins'
      end if
      counts(energy_bin, pitch_bin) = counts(energy_bin, pitch_bin) + 1_Int64
    end do

    if (sum(counts) /= number_of_samples) then
      error stop 'Histogram count total does not equal n_samples'
    end if

    allocate(sampled_f_array(nenergy, npitch))
    sampled_f_array = real(counts, Float64) / real(number_of_samples, Float64) &
      * sum(distribution%f_array)
  end subroutine sample_distribution

  integer function histogram_bin(value, first_center, spacing, number_of_bins)
    real(Float64), intent(in) :: value, first_center, spacing
    integer, intent(in) :: number_of_bins
    real(Float64) :: fractional_bin

    fractional_bin = (value - (first_center - 0.5_Float64 * spacing)) / spacing
    if (fractional_bin < 0.0_Float64 .or. &
        fractional_bin >= real(number_of_bins, Float64)) then
      histogram_bin = 0
    else
      histogram_bin = floor(fractional_bin) + 1
    end if
  end function histogram_bin

  subroutine print_sampling_summary(distribution, number_of_samples, counts, &
      sampled_f_array)
    type(DistributionFunction2D), intent(in) :: distribution
    integer(Int64), intent(in) :: number_of_samples
    integer(Int64), intent(in) :: counts(:,:)
    real(Float64), intent(in) :: sampled_f_array(:,:)

    real(Float64) :: reference_integral, sampled_integral, bin_area

    bin_area = abs(distribution%denergy * distribution%dpitch)
    reference_integral = sum(distribution%f_array) * bin_area
    sampled_integral = sum(sampled_f_array) * bin_area

    write(output_unit, '(/,a)') '  sampling complete'
    write(output_unit, '(a,i0)') '  requested samples: ', number_of_samples
    write(output_unit, '(a,i0)') '  histogram counts: ', sum(counts)
    write(output_unit, '(a,es14.6)') '  reference integral:', reference_integral
    write(output_unit, '(a,es14.6)') '  sampled integral:  ', sampled_integral
    write(output_unit, '(a,es14.6)') '  sampled minimum:   ', minval(sampled_f_array)
    write(output_unit, '(a,es14.6)') '  sampled maximum:   ', maxval(sampled_f_array)
  end subroutine print_sampling_summary

end module test_002_sampling
