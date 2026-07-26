module test_004_hdf5_utils
  use iso_fortran_env, only: Int32 => int32, Float64 => real64, output_unit
  use hdf5
  use h5lt
  use hdf5_utils, only: h5ltread_dataset_double_scalar_f
  implicit none
  private

  public :: read_real_vector
  public :: read_single_real_value
  public :: read_scalar_real
  public :: read_scalar_integer
  public :: read_array_dimensions
  public :: check_hdf5

contains

  subroutine read_real_vector(file_id, dataset_name, values, filename)
    !+ Read an allocatable rank-one real dataset.
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    real(Float64), allocatable, intent(out) :: values(:)

    integer(HSIZE_T), allocatable :: dimensions(:)
    integer :: error

    call read_array_dimensions( &
      file_id, dataset_name, dimensions, filename)
    if (size(dimensions) /= 1) then
      write(output_unit, '(a)') trim(dataset_name)//' must be one-dimensional'
      error stop 'Unexpected HDF5 dataset shape'
    end if

    allocate(values(int(dimensions(1))))
    call h5ltread_dataset_double_f( &
      file_id, trim(dataset_name), values, dimensions, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
  end subroutine read_real_vector

  subroutine read_single_real_value( &
      file_id, dataset_name, value, filename)
    !+ Read the only real value contained in a non-scalar dataset.
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    real(Float64), intent(out) :: value

    integer(HSIZE_T), allocatable :: dimensions(:)
    real(Float64), allocatable :: values(:)
    integer :: error, number_of_values

    call read_array_dimensions( &
      file_id, dataset_name, dimensions, filename)
    number_of_values = int(product(dimensions))
    if (number_of_values /= 1) then
      write(output_unit, '(a)') &
        trim(dataset_name)//' must contain exactly one value'
      error stop 'Unexpected HDF5 dataset size'
    end if

    allocate(values(number_of_values))
    call h5ltread_dataset_double_f( &
      file_id, trim(dataset_name), values, dimensions, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
    value = values(1)
  end subroutine read_single_real_value

  subroutine read_scalar_real(file_id, dataset_name, value, filename)
    !+ Read a scalar real dataset.
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    real(Float64), intent(out) :: value

    integer :: error

    call h5ltread_dataset_double_scalar_f( &
      file_id, trim(dataset_name), value, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
  end subroutine read_scalar_real

  subroutine read_scalar_integer(file_id, dataset_name, value, filename)
    !+ Read a scalar integer dataset into a 32-bit integer.
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    integer(Int32), intent(out) :: value

    integer(HSIZE_T), parameter :: dimensions(1) = [1_HSIZE_T]
    integer :: error
    integer :: values(1)

    call h5ltread_dataset_int_f( &
      file_id, trim(dataset_name), values, dimensions, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
    value = int(values(1), Int32)
  end subroutine read_scalar_integer

  subroutine read_array_dimensions( &
      file_id, dataset_name, dimensions, filename)
    !+ Read and return the dimensions of a non-scalar HDF5 dataset.
    integer(HID_T), intent(in) :: file_id
      !+ Identifier of the open HDF5 file.
    character(len=*), intent(in) :: dataset_name, filename
      !+ Dataset name and source filename used for error reporting.
    integer(HSIZE_T), allocatable, intent(out) :: dimensions(:)
      !+ Dataset extent allocated with one entry per array dimension.

    integer(HID_T) :: dataset_id, dataspace_id
    integer(HSIZE_T), allocatable :: maximum_dimensions(:)
    integer :: error, rank

    call h5dopen_f(file_id, trim(dataset_name), dataset_id, error)
    call check_hdf5(error, 'opening '//trim(dataset_name), filename)
    call h5dget_space_f(dataset_id, dataspace_id, error)
    call check_hdf5( &
      error, 'getting '//trim(dataset_name)//' dataspace', filename)
    call h5sget_simple_extent_ndims_f(dataspace_id, rank, error)
    call check_hdf5( &
      error, 'getting '//trim(dataset_name)//' rank', filename)

    if (rank < 1) then
      write(output_unit, '(a)') trim(dataset_name)//' must be an array'
      error stop 'Expected an HDF5 array dataset'
    end if

    allocate(dimensions(rank), maximum_dimensions(rank))
    call h5sget_simple_extent_dims_f( &
      dataspace_id, dimensions, maximum_dimensions, error)
    call check_hdf5( &
      error, 'getting '//trim(dataset_name)//' dimensions', filename)

    call h5sclose_f(dataspace_id, error)
    call check_hdf5( &
      error, 'closing '//trim(dataset_name)//' dataspace', filename)
    call h5dclose_f(dataset_id, error)
    call check_hdf5(error, 'closing '//trim(dataset_name), filename)
  end subroutine read_array_dimensions

  subroutine check_hdf5(error, operation, filename)
    !+ Stop with file and operation context when an HDF5 call fails.
    integer, intent(in) :: error
      !+ HDF5 status code; a negative value indicates failure.
    character(len=*), intent(in) :: operation, filename
      !+ Description of the attempted operation and the affected file.

    if (error < 0) then
      write(output_unit, '(a)') 'HDF5 error while '//trim(operation)
      write(output_unit, '(a)') '  file: '//trim(filename)
      error stop 'HDF5 operation failed'
    end if
  end subroutine check_hdf5

end module test_004_hdf5_utils
