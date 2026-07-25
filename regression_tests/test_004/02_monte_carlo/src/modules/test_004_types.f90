module test_004_types
  use iso_fortran_env, only: Int32 => int32, Int64 => int64, Float64 => real64
  implicit none
  private

  integer, parameter, public :: path_length = 4096
  integer, parameter, public :: fidasim_string_length = 200
  integer, parameter, public :: comment_length = 1024
  integer, parameter, public :: selector_length = 64
  integer, parameter, public :: species_length = 16
  integer, parameter, public :: number_of_atomic_levels = 6

  type, public :: MonteCarloConfig
    integer(Int32) :: n_cases = 0
    character(len=path_length), allocatable :: distribution_files(:)
    character(len=fidasim_string_length), allocatable :: runids(:)

    character(len=fidasim_string_length) :: tables_filename = ''
    character(len=path_length) :: test_config = ''
    character(len=path_length) :: input_distribution_config = ''
    character(len=fidasim_string_length) :: output_directory = ''
    character(len=comment_length) :: case_comment = ''
    character(len=comment_length) :: implementation_comment = ''

    integer(Int64) :: n_markers = 0_Int64
    integer(Int32) :: reservoir_size = 0
    integer(Int32) :: seed = 0
    logical :: save_data = .false.

    real(Float64) :: neutral_density = 0.0_Float64
    real(Float64) :: neutral_energy = 0.0_Float64
    real(Float64) :: injection_angle = 0.0_Float64
    character(len=selector_length) :: level_split_method = ''
    real(Float64) :: level_decay = 0.0_Float64
  end type MonteCarloConfig

  type, public :: DistributionCase
    character(len=path_length) :: filename = ''
    real(Float64), allocatable :: energy(:)
    real(Float64), allocatable :: pitch(:)
    real(Float64), allocatable :: f(:,:)
    real(Float64) :: denf = 0.0_Float64
    real(Float64) :: selected_r = 0.0_Float64
    real(Float64) :: selected_z = 0.0_Float64
    character(len=species_length) :: species = ''
    integer(Int32) :: atomic_number = 0
    integer(Int32) :: mass_number = 0
    integer(Int32) :: charge_state = 0
    real(Float64) :: atomic_mass = 0.0_Float64
  end type DistributionCase

  type, public :: NeutralParameters
    real(Float64) :: atomic_mass = 0.0_Float64
    real(Float64) :: speed = 0.0_Float64
    real(Float64) :: velocity(3) = 0.0_Float64
    real(Float64) :: level_density(number_of_atomic_levels) = 0.0_Float64
  end type NeutralParameters

end module test_004_types
