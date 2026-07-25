module test_004_config
  use iso_fortran_env, only: Int32 => int32, Int64 => int64, &
    Float64 => real64, output_unit
  use test_004_types, only: MonteCarloConfig, path_length, &
    fidasim_string_length, comment_length, selector_length
  implicit none
  private

  ! This must remain synchronized with MAX_CASES in normalize_config.py.
  integer, parameter :: max_cases = 256

  public :: read_config, print_config

contains

  subroutine read_config(filename, config)
    character(len=*), intent(in) :: filename
    type(MonteCarloConfig), intent(out) :: config

    integer(Int32) :: n_cases, reservoir_size, seed
    integer(Int64) :: n_markers
    logical :: save_data
    real(Float64) :: neutral_density, neutral_energy
    real(Float64) :: injection_angle, level_decay
    character(len=path_length) :: distribution_files(max_cases)
    character(len=fidasim_string_length) :: runids(max_cases)
    character(len=fidasim_string_length) :: tables_filename
    character(len=path_length) :: test_config
    character(len=path_length) :: input_distribution_config
    character(len=fidasim_string_length) :: output_directory
    character(len=comment_length) :: case_comment, implementation_comment
    character(len=selector_length) :: level_split_method
    character(len=512) :: iomsg
    integer :: unit, ios, i

    namelist /run_test/ n_cases, distribution_files, runids, &
      tables_filename, test_config, input_distribution_config, &
      output_directory, case_comment, implementation_comment, n_markers, &
      reservoir_size, seed, save_data, neutral_density, neutral_energy, &
      injection_angle, level_split_method, level_decay

    n_cases = 0
    distribution_files = ''
    runids = ''
    tables_filename = ''
    test_config = ''
    input_distribution_config = ''
    output_directory = ''
    case_comment = ''
    implementation_comment = ''
    n_markers = 0_Int64
    reservoir_size = 0
    seed = 0
    save_data = .false.
    neutral_density = 0.0_Float64
    neutral_energy = 0.0_Float64
    injection_angle = 0.0_Float64
    level_split_method = ''
    level_decay = 0.0_Float64

    open(newunit=unit, file=trim(filename), status='old', action='read', &
      iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
      write(output_unit, '(a)') 'Cannot open configuration file: '//trim(filename)
      write(output_unit, '(a)') trim(iomsg)
      error stop 'Configuration file open failed'
    end if

    read(unit, nml=run_test, iostat=ios, iomsg=iomsg)
    close(unit)
    if (ios /= 0) then
      write(output_unit, '(a)') 'Cannot read /run_test/ namelist: '//trim(iomsg)
      error stop 'Configuration namelist read failed'
    end if

    ! Python has already validated the normalized values. These checks only
    ! protect the fixed-size staging arrays used by the Fortran namelist.
    if (n_cases < 1) error stop 'n_cases must be positive'
    if (n_cases > max_cases) then
      error stop 'n_cases exceeds the supported maximum of 256'
    end if

    do i = 1, n_cases
      if (len_trim(distribution_files(i)) == 0) then
        write(output_unit, '(a,i0)') 'Missing distribution_files entry ', i
        error stop 'Every configured case must have a distribution path'
      end if
      if (len_trim(runids(i)) == 0) then
        write(output_unit, '(a,i0)') 'Missing runids entry ', i
        error stop 'Every configured case must have a run ID'
      end if
    end do

    config%n_cases = n_cases
    config%tables_filename = trim(tables_filename)
    config%test_config = trim(test_config)
    config%input_distribution_config = trim(input_distribution_config)
    config%output_directory = trim(output_directory)
    config%case_comment = trim(case_comment)
    config%implementation_comment = trim(implementation_comment)
    config%n_markers = n_markers
    config%reservoir_size = reservoir_size
    config%seed = seed
    config%save_data = save_data
    config%neutral_density = neutral_density
    config%neutral_energy = neutral_energy
    config%injection_angle = injection_angle
    config%level_split_method = trim(level_split_method)
    config%level_decay = level_decay

    allocate(config%distribution_files(n_cases))
    allocate(config%runids(n_cases))
    do i = 1, n_cases
      config%distribution_files(i) = trim(distribution_files(i))
      config%runids(i) = trim(runids(i))
    end do
  end subroutine read_config

  subroutine print_config(config)
    type(MonteCarloConfig), intent(in) :: config
    integer :: i

    write(output_unit, '(a)') 'test_004 Monte Carlo configuration'
    write(output_unit, '(a,i0)') '  n_cases: ', config%n_cases
    do i = 1, config%n_cases
      write(output_unit, '(a,i0,a)') '  case ', i, ':'
      write(output_unit, '(a,a)') '    distribution_file: ', &
        trim(config%distribution_files(i))
      write(output_unit, '(a,a)') '    runid: ', trim(config%runids(i))
    end do
    write(output_unit, '(a,a)') '  tables_filename: ', &
      trim(config%tables_filename)
    write(output_unit, '(a,a)') '  test_config: ', trim(config%test_config)
    write(output_unit, '(a,a)') '  input_distribution_config: ', &
      trim(config%input_distribution_config)
    write(output_unit, '(a,a)') '  output_directory: ', &
      trim(config%output_directory)
    write(output_unit, '(a,a)') '  case_comment: ', trim(config%case_comment)
    write(output_unit, '(a,a)') '  implementation_comment: ', &
      trim(config%implementation_comment)
    write(output_unit, '(a,i0)') '  n_markers: ', config%n_markers
    write(output_unit, '(a,i0)') '  reservoir_size: ', config%reservoir_size
    write(output_unit, '(a,i0)') '  seed: ', config%seed
    write(output_unit, '(a,l1)') '  save_data: ', config%save_data
    write(output_unit, '(a,es14.6)') '  neutral_density [cm^-3]: ', &
      config%neutral_density
    write(output_unit, '(a,es14.6)') '  neutral_energy [keV]: ', &
      config%neutral_energy
    write(output_unit, '(a,es14.6)') '  injection_angle [degrees]: ', &
      config%injection_angle
    write(output_unit, '(a,a)') '  level_split_method: ', &
      trim(config%level_split_method)
    write(output_unit, '(a,es14.6)') '  level_decay: ', config%level_decay
  end subroutine print_config

end module test_004_config
