program test_003
  use iso_fortran_env, only: Int64 => int64, Float64 => real64
  use test_003_config, only: SamplerConfig, read_config, print_config
  use test_003_hdf5, only: DistributionFunction2D, &
    read_reference_distribution, print_reference_distribution, &
    write_sampled_distribution
  use test_003_sampling, only: initialize_serial_rng, sample_distribution, &
    print_sampling_summary
  implicit none

  type(SamplerConfig) :: config
  type(DistributionFunction2D) :: distribution
  integer(Int64), allocatable :: counts(:,:)
  real(Float64), allocatable :: sampled_f_array(:,:)
  character(len=4096) :: config_filename
  integer :: i

  if (command_argument_count() /= 1) then
    write(*, '(a)') 'Usage: test_003 <input_config.nml>'
    error stop 'Expected exactly one command-line argument'
  end if

  call get_command_argument(1, config_filename)
  call read_config(trim(config_filename), config)
  call print_config(config)

  do i = 1, config%n_reference_files
    call read_reference_distribution(trim(config%reference_files(i)), distribution)
    call print_reference_distribution(trim(config%reference_files(i)), distribution)
    call initialize_serial_rng(config%seed)
    call sample_distribution(distribution, config%n_samples, counts, &
      sampled_f_array)
    call print_sampling_summary(distribution, config%n_samples, counts, &
      sampled_f_array)
    call write_sampled_distribution(trim(config%reference_files(i)), &
      trim(config%output_directory), sampled_f_array, &
      distribution%f_array_dimensions, distribution%denergy, &
      distribution%dpitch, config%n_samples, config%seed)
  end do
end program test_003
