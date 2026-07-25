module test_004_hdf5
  use iso_fortran_env, only: Int32 => int32, Float64 => real64, output_unit
  use iso_c_binding, only: c_ptr, c_null_ptr, c_loc, c_f_pointer, &
    c_size_t, c_int, c_int64_t, c_char, c_associated
  use hdf5
  use h5lt
  use test_004_types, only: DistributionCase
  implicit none
  private

  public :: read_distribution, print_distribution, release_distribution

  interface
    function h5dread_c(dataset_id, memory_type_id, memory_space_id, &
        file_space_id, transfer_id, buffer) bind(C, name='H5Dread') &
        result(error)
      import :: c_ptr, c_int, c_int64_t
      integer(c_int64_t), value :: dataset_id, memory_type_id
      integer(c_int64_t), value :: memory_space_id, file_space_id, transfer_id
      type(c_ptr), value :: buffer
      integer(c_int) :: error
    end function h5dread_c

    function string_length_c(string_pointer) bind(C, name='strlen') &
        result(string_length)
      import :: c_ptr, c_size_t
      type(c_ptr), value :: string_pointer
      integer(c_size_t) :: string_length
    end function string_length_c

    function h5free_memory_c(memory_pointer) bind(C, name='H5free_memory') &
        result(error)
      import :: c_ptr, c_int
      type(c_ptr), value :: memory_pointer
      integer(c_int) :: error
    end function h5free_memory_c
  end interface

contains

  subroutine read_distribution(filename, distribution)
    character(len=*), intent(in) :: filename
    type(DistributionCase), intent(out) :: distribution

    integer(HID_T) :: file_id
    integer :: error

    distribution%filename = trim(filename)

    call h5open_f(error)
    call check_hdf5(error, 'initializing HDF5', filename)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
    call check_hdf5(error, 'opening distribution file', filename)

    call read_vector(file_id, 'energy', distribution%energy, filename)
    call read_vector(file_id, 'pitch', distribution%pitch, filename)
    call read_distribution_array(file_id, distribution, filename)
    call read_single_value(file_id, 'denf', distribution%denf, filename)
    call read_single_value(file_id, 'r', distribution%selected_r, filename)
    call read_single_value(file_id, 'z', distribution%selected_z, filename)
    call read_string(file_id, 'species', distribution%species, filename)
    call read_integer(file_id, 'atomic_number', &
      distribution%atomic_number, filename)
    call read_integer(file_id, 'mass_number', &
      distribution%mass_number, filename)
    call read_integer(file_id, 'charge_state', &
      distribution%charge_state, filename)
    call read_single_value(file_id, 'A', distribution%atomic_mass, filename)

    call h5fclose_f(file_id, error)
    call check_hdf5(error, 'closing distribution file', filename)
    call h5close_f(error)
    call check_hdf5(error, 'closing HDF5', filename)
  end subroutine read_distribution

  subroutine read_vector(file_id, dataset_name, values, filename)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    real(Float64), allocatable, intent(out) :: values(:)

    integer(HSIZE_T), allocatable :: dimensions(:)
    integer :: error

    call read_dimensions(file_id, dataset_name, dimensions, filename)
    if (size(dimensions) /= 1) then
      write(output_unit, '(a)') trim(dataset_name)//' must be one-dimensional'
      error stop 'Unexpected Test 002 dataset shape'
    end if

    allocate(values(int(dimensions(1))))
    call h5ltread_dataset_double_f(file_id, trim(dataset_name), values, &
      dimensions, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
  end subroutine read_vector

  subroutine read_distribution_array(file_id, distribution, filename)
    integer(HID_T), intent(in) :: file_id
    type(DistributionCase), intent(inout) :: distribution
    character(len=*), intent(in) :: filename

    integer(HID_T) :: dataset_id, file_space_id, memory_space_id
    integer(HSIZE_T), allocatable :: dimensions(:)
    integer(HSIZE_T) :: memory_dimensions(2)
    integer(HSIZE_T) :: start(4), count(4)
    integer :: error, nenergy, npitch

    nenergy = size(distribution%energy)
    npitch = size(distribution%pitch)
    call read_dimensions(file_id, 'f', dimensions, filename)

    if (size(dimensions) /= 4) then
      error stop 'The Test 002 f dataset must be four-dimensional'
    end if
    if (dimensions(1) /= nenergy .or. dimensions(2) /= npitch .or. &
        dimensions(3) /= 1 .or. dimensions(4) /= 1) then
      write(output_unit, '(a,4(1x,i0))') 'Unexpected f dimensions:', dimensions
      error stop 'The Test 002 f dataset has an unexpected shape'
    end if

    allocate(distribution%f(nenergy, npitch))
    memory_dimensions = [int(nenergy, HSIZE_T), int(npitch, HSIZE_T)]
    start = [0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T]
    count = dimensions

    call h5dopen_f(file_id, 'f', dataset_id, error)
    call check_hdf5(error, 'opening f', filename)
    call h5dget_space_f(dataset_id, file_space_id, error)
    call check_hdf5(error, 'getting f dataspace', filename)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, start, &
      count, error)
    call check_hdf5(error, 'selecting f spatial slice', filename)
    call h5screate_simple_f(2, memory_dimensions, memory_space_id, error)
    call check_hdf5(error, 'creating f memory space', filename)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, distribution%f, &
      memory_dimensions, error, memory_space_id, file_space_id)
    call check_hdf5(error, 'reading f spatial slice', filename)

    call h5sclose_f(memory_space_id, error)
    call check_hdf5(error, 'closing f memory space', filename)
    call h5sclose_f(file_space_id, error)
    call check_hdf5(error, 'closing f dataspace', filename)
    call h5dclose_f(dataset_id, error)
    call check_hdf5(error, 'closing f', filename)
  end subroutine read_distribution_array

  subroutine read_single_value(file_id, dataset_name, value, filename)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    real(Float64), intent(out) :: value

    integer(HSIZE_T), allocatable :: dimensions(:)
    real(Float64), allocatable :: values(:)
    integer :: error, number_of_values

    call read_dimensions(file_id, dataset_name, dimensions, filename)
    number_of_values = int(product(dimensions))
    if (number_of_values /= 1) then
      write(output_unit, '(a)') trim(dataset_name)//' must contain one value'
      error stop 'Unexpected Test 002 dataset size'
    end if

    allocate(values(number_of_values))
    call h5ltread_dataset_double_f(file_id, trim(dataset_name), values, &
      dimensions, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
    value = values(1)
  end subroutine read_single_value

  subroutine read_integer(file_id, dataset_name, value, filename)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    integer(Int32), intent(out) :: value

    integer(HSIZE_T), parameter :: dimensions(1) = [1_HSIZE_T]
    integer :: error
    integer :: values(1)

    call h5ltread_dataset_int_f(file_id, trim(dataset_name), values, &
      dimensions, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
    value = int(values(1), Int32)
  end subroutine read_integer

  subroutine read_string(file_id, dataset_name, value, filename)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    character(len=*), intent(out) :: value

    integer(HID_T) :: dataset_id, datatype_id
    type(c_ptr), target :: string_pointer
    character(kind=c_char), pointer :: characters(:)
    integer(c_size_t) :: string_length
    integer :: error, character_index, characters_to_copy

    value = ''
    string_pointer = c_null_ptr
    call h5dopen_f(file_id, trim(dataset_name), dataset_id, error)
    call check_hdf5(error, 'opening '//trim(dataset_name), filename)
    call h5dget_type_f(dataset_id, datatype_id, error)
    call check_hdf5(error, 'getting '//trim(dataset_name)//' type', filename)

    error = h5dread_c(int(dataset_id, c_int64_t), &
      int(datatype_id, c_int64_t), 0_c_int64_t, 0_c_int64_t, &
      0_c_int64_t, c_loc(string_pointer))
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
    if (.not. c_associated(string_pointer)) then
      error stop 'The species dataset contains a null string'
    end if

    string_length = string_length_c(string_pointer)
    characters_to_copy = min(int(string_length), len(value))
    call c_f_pointer(string_pointer, characters, [int(string_length)])
    do character_index = 1, characters_to_copy
      value(character_index:character_index) = characters(character_index)
    end do

    error = h5free_memory_c(string_pointer)
    call check_hdf5(error, 'releasing '//trim(dataset_name)//' string', &
      filename)
    call h5tclose_f(datatype_id, error)
    call check_hdf5(error, 'closing '//trim(dataset_name)//' type', filename)
    call h5dclose_f(dataset_id, error)
    call check_hdf5(error, 'closing '//trim(dataset_name), filename)
  end subroutine read_string

  subroutine read_dimensions(file_id, dataset_name, dimensions, filename)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    integer(HSIZE_T), allocatable, intent(out) :: dimensions(:)

    integer(HID_T) :: dataset_id, dataspace_id
    integer(HSIZE_T), allocatable :: maximum_dimensions(:)
    integer :: error, rank

    call h5dopen_f(file_id, trim(dataset_name), dataset_id, error)
    call check_hdf5(error, 'opening '//trim(dataset_name), filename)
    call h5dget_space_f(dataset_id, dataspace_id, error)
    call check_hdf5(error, 'getting '//trim(dataset_name)//' dataspace', &
      filename)
    call h5sget_simple_extent_ndims_f(dataspace_id, rank, error)
    call check_hdf5(error, 'getting '//trim(dataset_name)//' rank', filename)

    if (rank == 0) then
      allocate(dimensions(1), maximum_dimensions(1))
      dimensions = 1_HSIZE_T
      maximum_dimensions = 1_HSIZE_T
    else
      allocate(dimensions(rank), maximum_dimensions(rank))
      call h5sget_simple_extent_dims_f(dataspace_id, dimensions, &
        maximum_dimensions, error)
      call check_hdf5(error, 'getting '//trim(dataset_name)//' dimensions', &
        filename)
    end if

    call h5sclose_f(dataspace_id, error)
    call check_hdf5(error, 'closing '//trim(dataset_name)//' dataspace', &
      filename)
    call h5dclose_f(dataset_id, error)
    call check_hdf5(error, 'closing '//trim(dataset_name), filename)
  end subroutine read_dimensions

  subroutine print_distribution(distribution)
    type(DistributionCase), intent(in) :: distribution

    write(output_unit, '(/,a)') trim(distribution%filename)
    write(output_unit, '(a,i0)') '  energy points: ', size(distribution%energy)
    write(output_unit, '(a,2(1x,es14.6))') '  energy range [keV]:', &
      minval(distribution%energy), maxval(distribution%energy)
    write(output_unit, '(a,i0)') '  pitch points: ', size(distribution%pitch)
    write(output_unit, '(a,2(1x,es14.6))') '  pitch range:', &
      minval(distribution%pitch), maxval(distribution%pitch)
    write(output_unit, '(a,2(1x,i0))') '  Fortran f shape:', &
      shape(distribution%f)
    write(output_unit, '(a,es14.6)') '  f minimum: ', minval(distribution%f)
    write(output_unit, '(a,es14.6)') '  f maximum: ', maxval(distribution%f)
    write(output_unit, '(a,es14.6)') '  denf [ions/cm^3]: ', distribution%denf
    write(output_unit, '(a,2(1x,es14.6))') '  location R,Z [cm]:', &
      distribution%selected_r, distribution%selected_z
    write(output_unit, '(a,a)') '  species: ', trim(distribution%species)
    write(output_unit, '(a,i0)') '  atomic number: ', &
      distribution%atomic_number
    write(output_unit, '(a,i0)') '  mass number: ', distribution%mass_number
    write(output_unit, '(a,i0)') '  charge state: ', distribution%charge_state
    write(output_unit, '(a,es14.6)') '  atomic mass [amu]: ', &
      distribution%atomic_mass
  end subroutine print_distribution

  subroutine release_distribution(distribution)
    type(DistributionCase), intent(inout) :: distribution

    if (allocated(distribution%energy)) deallocate(distribution%energy)
    if (allocated(distribution%pitch)) deallocate(distribution%pitch)
    if (allocated(distribution%f)) deallocate(distribution%f)
  end subroutine release_distribution

  subroutine check_hdf5(error, operation, filename)
    integer, intent(in) :: error
    character(len=*), intent(in) :: operation, filename

    if (error < 0) then
      write(output_unit, '(a)') 'HDF5 error while '//trim(operation)
      write(output_unit, '(a)') '  file: '//trim(filename)
      error stop 'Unable to read Test 002 distribution'
    end if
  end subroutine check_hdf5

end module test_004_hdf5
