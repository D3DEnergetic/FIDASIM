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
      print_neutral_parameters, &
      initialize_neutral_population, &
      print_neutral_population, &
      release_neutral_population
  use test_004_setup, &
    only: &
      configure_fidasim, &
      initialize_serial_rng, &
      reset_serial_rng, &
      initialize_atomic_tables, &
      print_atomic_table_setup, &
      initialize_test_setup, &
      print_test_setup, &
      teardown_test_setup
  use test_004_sampling, &
    only: &
      IonSample, &
      sample_ion, &
      validate_ion_sample, &
      print_ion_sample
  implicit none

  type(MonteCarloConfig) :: config
  type(DistributionCase) :: distribution
  type(NeutralParameters) :: neutral
  type(IonSample) :: ion_sample
  character(len=4096) :: config_filename
  integer :: case_index

  if (command_argument_count() /= 1) then
    write(*, '(a)') 'Usage: test_004 <normalized_input_config.nml>'
    error stop 'Expected exactly one command-line argument'
  end if

  call get_command_argument(1, config_filename)
  call read_config(trim(config_filename), config)
  call print_config(config)
  call configure_fidasim(config)
  call initialize_atomic_tables()
  call print_atomic_table_setup()
  call initialize_serial_rng()

  ! Loop over each distribution case specified in the configuration file:
  do case_index = 1, config%n_cases
    call read_distribution( &
      trim(config%distribution_files(case_index)), distribution)
    call print_distribution(distribution)
    call build_neutral_parameters(config, distribution, neutral)
    call print_neutral_parameters(neutral)
    call initialize_test_setup(distribution)
    call print_test_setup()
    call initialize_neutral_population(config, neutral)
    call print_neutral_population(config, neutral)
    call reset_serial_rng()
    call sample_ion(distribution, ion_sample)
    call validate_ion_sample(distribution, ion_sample)
    call print_ion_sample(ion_sample)
    call release_neutral_population()
    call teardown_test_setup()
    call release_distribution(distribution)
  end do

end program test_004
