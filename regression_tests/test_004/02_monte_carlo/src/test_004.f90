program test_004
  use test_004_types
  use test_004_config, &
    only: &
      read_config, &
      print_config
  use test_004_hdf5, &
    only: &
      read_distribution, &
      print_distribution, &
      release_distribution
  use test_004_neutral, &
    only: &
      build_neutral_parameters, &
      print_neutral_parameters
  use test_004_setup, &
    only: &
      initialize_test_setup, &
      print_test_setup, &
      teardown_test_setup
  implicit none

  type(MonteCarloConfig) :: config
  type(DistributionCase) :: distribution
  type(NeutralParameters) :: neutral
  character(len=4096) :: config_filename
  integer :: case_index

  if (command_argument_count() /= 1) then
    write(*, '(a)') 'Usage: test_004 <normalized_input_config.nml>'
    error stop 'Expected exactly one command-line argument'
  end if

  call get_command_argument(1, config_filename)
  call read_config(trim(config_filename), config)
  call print_config(config)

  ! Loop over each distribution case specified in the configuration file:
  do case_index = 1, config%n_cases
    call read_distribution( &
      trim(config%distribution_files(case_index)), distribution)
    call print_distribution(distribution)
    call build_neutral_parameters(config, distribution, neutral)
    call print_neutral_parameters(neutral)
    call initialize_test_setup(distribution)
    call print_test_setup()
    call teardown_test_setup()
    call release_distribution(distribution)
  end do

end program test_004
