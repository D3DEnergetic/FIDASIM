module test_002_config
  use iso_fortran_env, only: Int32 => int32, Int64 => int64, output_unit
  implicit none
  private

  integer, parameter :: max_reference_files = 256
  integer, parameter :: path_length = 4096

  type, public :: SamplerConfig
    integer(Int32) :: n_reference_files = 0
    character(len=path_length), allocatable :: reference_files(:)
    character(len=path_length) :: output_directory = ''
    integer(Int64) :: n_samples = 0_Int64
    integer(Int32) :: seed = 0
    logical :: plot_data = .false.
  end type SamplerConfig

  public :: read_config, print_config

contains

  subroutine read_config(filename, config)
    character(len=*), intent(in) :: filename
    type(SamplerConfig), intent(out) :: config

    integer(Int32) :: n_reference_files, seed
    integer(Int64) :: n_samples
    logical :: plot_data
    character(len=path_length) :: reference_files(max_reference_files)
    character(len=path_length) :: output_directory
    character(len=512) :: iomsg
    integer :: unit, ios, i
    namelist /run_test/ n_reference_files, reference_files, &
      output_directory, n_samples, seed, plot_data

    n_reference_files = 0
    reference_files = ''
    output_directory = ''
    n_samples = 0_Int64
    seed = 0
    plot_data = .false.

    open(newunit=unit, file=trim(filename), status='old', action='read', &
      iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
      write(*, '(a)') 'Cannot open configuration file: '//trim(filename)
      error stop trim(iomsg)
    end if

    read(unit, nml=run_test, iostat=ios, iomsg=iomsg)
    close(unit)
    if (ios /= 0) error stop 'Cannot read /run_test/ namelist: '//trim(iomsg)

    if (n_reference_files < 1) error stop 'n_reference_files must be positive'
    if (n_reference_files > max_reference_files) then
      error stop 'n_reference_files exceeds the supported maximum of 256'
    end if
    if (n_samples < 1_Int64) error stop 'n_samples must be positive'
    if (seed < 1) error stop 'seed must be positive'
    if (len_trim(output_directory) == 0) error stop 'output_directory must not be empty'

    do i = 1, n_reference_files
      if (len_trim(reference_files(i)) == 0) then
        write(*, '(a,i0)') 'Missing reference_files entry ', i
        error stop 'Every configured reference file must have a path'
      end if
    end do

    config%n_reference_files = n_reference_files
    config%n_samples = n_samples
    config%seed = seed
    config%plot_data = plot_data
    config%output_directory = trim(output_directory)
    allocate(config%reference_files(n_reference_files))
    do i = 1, n_reference_files
      config%reference_files(i) = trim(reference_files(i))
    end do
  end subroutine read_config

  subroutine print_config(config)
    type(SamplerConfig), intent(in) :: config
    integer :: i

    write(output_unit, '(a)') 'test_002 configuration'
    write(output_unit, '(a,i0)') '  n_reference_files: ', config%n_reference_files
    do i = 1, config%n_reference_files
      write(output_unit, '(a,i0,a,a)') '  reference_files(', i, '): ', &
        trim(config%reference_files(i))
    end do
    write(output_unit, '(a,a)') '  output_directory: ', trim(config%output_directory)
    write(output_unit, '(a,i0)') '  n_samples: ', config%n_samples
    write(output_unit, '(a,i0)') '  seed: ', config%seed
    write(output_unit, '(a,l1)') '  plot_data: ', config%plot_data
  end subroutine print_config

end module test_002_config

program test_002
  use iso_fortran_env, only: Int64 => int64, Float64 => real64
  use test_002_config, only: SamplerConfig, read_config, print_config
  use test_002_hdf5, only: DistributionFunction2D, &
    read_reference_distribution, print_reference_distribution, &
    write_sampled_distribution
  use test_002_sampling, only: initialize_serial_rng, sample_distribution, &
    print_sampling_summary
  implicit none

  type(SamplerConfig) :: config
  type(DistributionFunction2D) :: distribution
  integer(Int64), allocatable :: counts(:,:)
  real(Float64), allocatable :: sampled_f_array(:,:)
  character(len=4096) :: config_filename
  integer :: i

  if (command_argument_count() /= 1) then
    write(*, '(a)') 'Usage: test_002 <input_config.nml>'
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
      distribution%f_array_dimensions, config%n_samples, config%seed)
  end do
end program test_002
