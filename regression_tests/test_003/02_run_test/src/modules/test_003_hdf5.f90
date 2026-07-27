module test_003_hdf5
  use iso_fortran_env, only: Int32 => int32, Int64 => int64, &
    Float64 => real64, output_unit
  use ieee_arithmetic, only: ieee_is_finite
  use hdf5
  use h5lt
  implicit none
  private

  type, public :: DistributionFunction2D
    real(Float64), allocatable :: energy_grid(:)
    real(Float64), allocatable :: pitch_grid(:)
    real(Float64), allocatable :: f_array(:,:)
    integer(HSIZE_T), allocatable :: f_array_dimensions(:)
    real(Float64) :: denergy = 0.0_Float64
    real(Float64) :: dpitch = 0.0_Float64
  end type DistributionFunction2D

  public :: read_reference_distribution, print_reference_distribution
  public :: write_sampled_distribution

contains

  subroutine read_reference_distribution(filename, distribution)
    character(len=*), intent(in) :: filename
    type(DistributionFunction2D), intent(out) :: distribution

    integer(HID_T) :: file_id
    integer :: error

    call h5open_f(error)
    call check_hdf5(error, 'initializing HDF5', filename)
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
    call check_hdf5(error, 'opening reference file', filename)

    call read_grid(file_id, 'energy', distribution%energy_grid, filename)
    call read_grid(file_id, 'pitch', distribution%pitch_grid, filename)
    call read_dataset_dimensions(file_id, 'f', &
      distribution%f_array_dimensions, filename)
    call read_f_array(file_id, distribution, filename)
    call validate_reference_data(distribution, filename)

    call h5fclose_f(file_id, error)
    call check_hdf5(error, 'closing reference file', filename)
    call h5close_f(error)
    call check_hdf5(error, 'closing HDF5', filename)
  end subroutine read_reference_distribution

  subroutine read_grid(file_id, dataset_name, grid, filename)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    real(Float64), allocatable, intent(out) :: grid(:)
    integer(HSIZE_T), allocatable :: dimensions(:)
    integer :: error

    call read_dataset_dimensions(file_id, dataset_name, dimensions, filename)
    if (size(dimensions) /= 1) then
      error stop trim(dataset_name)//' must be one-dimensional in '//trim(filename)
    end if
    allocate(grid(int(dimensions(1))))
    call h5ltread_dataset_double_f(file_id, trim(dataset_name), grid, &
      dimensions, error)
    call check_hdf5(error, 'reading '//trim(dataset_name), filename)
  end subroutine read_grid

  subroutine read_f_array(file_id, distribution, filename)
    integer(HID_T), intent(in) :: file_id
    type(DistributionFunction2D), intent(inout) :: distribution
    character(len=*), intent(in) :: filename
    integer(HID_T) :: dataset_id, file_space_id, memory_space_id
    integer(HSIZE_T) :: memory_dimensions(2)
    integer(HSIZE_T) :: start(4), count(4)
    integer :: error, nenergy, npitch

    nenergy = size(distribution%energy_grid)
    npitch = size(distribution%pitch_grid)
    if (size(distribution%f_array_dimensions) /= 4) then
      error stop 'f must be four-dimensional in '//trim(filename)
    end if

    ! Test 002 writes the h5py-visible order (z, r, pitch, energy). The
    ! Fortran HDF5 interface presents those dimensions here in reverse order
    ! as (energy, pitch, r, z).
    if (distribution%f_array_dimensions(1) /= nenergy .or. &
        distribution%f_array_dimensions(2) /= npitch .or. &
        distribution%f_array_dimensions(3) /= 1 .or. &
        distribution%f_array_dimensions(4) /= 1) then
      write(*, '(a)') 'f dimensions do not match the expected Test 002 layout'
      write(*, '(a,4(1x,i0))') '  Fortran HDF5 dimensions:', &
        distribution%f_array_dimensions
      write(*, '(a,4(1x,i0))') '  Expected:', nenergy, npitch, 1, 1
      error stop trim(filename)
    end if

    allocate(distribution%f_array(nenergy, npitch))
    memory_dimensions = [int(nenergy, HSIZE_T), int(npitch, HSIZE_T)]
    start = [0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T, 0_HSIZE_T]
    count = distribution%f_array_dimensions

    call h5dopen_f(file_id, 'f', dataset_id, error)
    call check_hdf5(error, 'opening f', filename)
    call h5dget_space_f(dataset_id, file_space_id, error)
    call check_hdf5(error, 'getting f dataspace', filename)
    call h5sselect_hyperslab_f(file_space_id, H5S_SELECT_SET_F, start, &
      count, error)
    call check_hdf5(error, 'selecting f spatial slice', filename)
    call h5screate_simple_f(2, memory_dimensions, memory_space_id, error)
    call check_hdf5(error, 'creating f memory space', filename)
    call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, distribution%f_array, &
      memory_dimensions, error, memory_space_id, file_space_id)
    call check_hdf5(error, 'reading f spatial slice', filename)
    call h5sclose_f(memory_space_id, error)
    call check_hdf5(error, 'closing f memory space', filename)
    call h5sclose_f(file_space_id, error)
    call check_hdf5(error, 'closing f dataspace', filename)
    call h5dclose_f(dataset_id, error)
    call check_hdf5(error, 'closing f', filename)

    ! Reuse this array as the Fortran dimension argument for the sampled
    ! output. Supplying (pitch, energy) makes the HDF5 file visible to h5py as
    ! the desired (energy, pitch) shape.
    distribution%f_array_dimensions = [int(npitch, HSIZE_T), &
      int(nenergy, HSIZE_T)]
  end subroutine read_f_array

  subroutine validate_reference_data(distribution, filename)
    type(DistributionFunction2D), intent(inout) :: distribution
    character(len=*), intent(in) :: filename
    real(Float64) :: total

    if (size(distribution%energy_grid) < 2) then
      error stop 'energy_grid must contain at least two points in '//trim(filename)
    end if
    if (size(distribution%pitch_grid) < 2) then
      error stop 'pitch_grid must contain at least two points in '//trim(filename)
    end if
    if (.not. all(ieee_is_finite(distribution%energy_grid))) then
      error stop 'energy_grid contains a non-finite value in '//trim(filename)
    end if
    if (.not. all(ieee_is_finite(distribution%pitch_grid))) then
      error stop 'pitch_grid contains a non-finite value in '//trim(filename)
    end if

    distribution%denergy = distribution%energy_grid(2) - distribution%energy_grid(1)
    distribution%dpitch = distribution%pitch_grid(2) - distribution%pitch_grid(1)
    call validate_strict_monotonicity(distribution%energy_grid, &
      distribution%denergy, 'energy_grid', filename)
    call validate_strict_monotonicity(distribution%pitch_grid, &
      distribution%dpitch, 'pitch_grid', filename)
    call validate_uniform_spacing(distribution%energy_grid, &
      distribution%denergy, 'energy_grid', filename)
    call validate_uniform_spacing(distribution%pitch_grid, &
      distribution%dpitch, 'pitch_grid', filename)

    if (.not. all(ieee_is_finite(distribution%f_array))) then
      error stop 'f_array contains a non-finite value in '//trim(filename)
    end if
    if (any(distribution%f_array < 0.0_Float64)) then
      error stop 'f_array contains a negative value in '//trim(filename)
    end if
    total = sum(distribution%f_array)
    if (.not. ieee_is_finite(total) .or. total <= 0.0_Float64) then
      error stop 'f_array must have a finite, positive sum in '//trim(filename)
    end if
  end subroutine validate_reference_data

  subroutine validate_strict_monotonicity(grid, spacing, grid_name, filename)
    real(Float64), intent(in) :: grid(:), spacing
    character(len=*), intent(in) :: grid_name, filename

    if (spacing > 0.0_Float64) then
      if (any(grid(2:) <= grid(:size(grid)-1))) then
        error stop trim(grid_name)//' must be strictly monotonic in '//trim(filename)
      end if
    else if (spacing < 0.0_Float64) then
      if (any(grid(2:) >= grid(:size(grid)-1))) then
        error stop trim(grid_name)//' must be strictly monotonic in '//trim(filename)
      end if
    else
      error stop trim(grid_name)//' must be strictly monotonic in '//trim(filename)
    end if
  end subroutine validate_strict_monotonicity

  subroutine validate_uniform_spacing(grid, spacing, grid_name, filename)
    real(Float64), intent(in) :: grid(:), spacing
    character(len=*), intent(in) :: grid_name, filename
    real(Float64) :: tolerance, scale

    scale = max(1.0_Float64, maxval(abs(grid)), abs(spacing))
    tolerance = 100.0_Float64 * epsilon(1.0_Float64) * scale
    if (any(abs((grid(2:) - grid(:size(grid)-1)) - spacing) > tolerance)) then
      error stop trim(grid_name)//' must be uniformly spaced in '//trim(filename)
    end if
  end subroutine validate_uniform_spacing

  subroutine read_dataset_dimensions(file_id, dataset_name, dimensions, filename)
    !+ Return an HDF5 dataset's extents in Fortran interface order.
    !+
    !+ For multidimensional data, this is the reverse of the raw file order
    !+ exposed by C-based readers such as h5py.
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: dataset_name, filename
    integer(HSIZE_T), allocatable, intent(out) :: dimensions(:)
    integer(HID_T) :: dataset_id, dataspace_id
    integer(HSIZE_T), allocatable :: maximum_dimensions(:)
    integer :: error, rank

    call h5dopen_f(file_id, trim(dataset_name), dataset_id, error)
    call check_hdf5(error, 'opening '//trim(dataset_name), filename)
    call h5dget_space_f(dataset_id, dataspace_id, error)
    call check_hdf5(error, 'getting dataspace for '//trim(dataset_name), filename)
    call h5sget_simple_extent_ndims_f(dataspace_id, rank, error)
    call check_hdf5(error, 'getting rank for '//trim(dataset_name), filename)
    allocate(dimensions(rank), maximum_dimensions(rank))
    call h5sget_simple_extent_dims_f(dataspace_id, dimensions, &
      maximum_dimensions, error)
    call check_hdf5(error, 'getting dimensions for '//trim(dataset_name), filename)
    call h5sclose_f(dataspace_id, error)
    call check_hdf5(error, 'closing dataspace for '//trim(dataset_name), filename)
    call h5dclose_f(dataset_id, error)
    call check_hdf5(error, 'closing '//trim(dataset_name), filename)
  end subroutine read_dataset_dimensions

  subroutine print_reference_distribution(filename, distribution)
    character(len=*), intent(in) :: filename
    type(DistributionFunction2D), intent(in) :: distribution
    integer :: i

    write(output_unit, '(/,a)') trim(filename)
    write(output_unit, '(a,i0)') '  energy points: ', size(distribution%energy_grid)
    write(output_unit, '(a,i0)') '  pitch points:  ', size(distribution%pitch_grid)
    write(output_unit, '(a)', advance='no') '  f_array dimensions:'
    do i = 1, size(distribution%f_array_dimensions)
      write(output_unit, '(1x,i0)', advance='no') distribution%f_array_dimensions(i)
    end do
    write(output_unit, *)
    write(output_unit, '(a,2(1x,i0))') '  Fortran f_array shape:', &
      shape(distribution%f_array)
    write(output_unit, '(a,es14.6)') '  minimum value: ', minval(distribution%f_array)
    write(output_unit, '(a,es14.6)') '  maximum value: ', maxval(distribution%f_array)
    write(output_unit, '(a,es14.6)') '  sum:           ', sum(distribution%f_array)
    write(output_unit, '(a,es14.6)') '  energy spacing:', distribution%denergy
    write(output_unit, '(a,es14.6)') '  pitch spacing: ', distribution%dpitch
    write(output_unit, '(a)') '  validation: passed'
  end subroutine print_reference_distribution

  subroutine write_sampled_distribution(reference_filename, output_directory, &
      sampled_f_array, f_array_dimensions, denergy, dpitch, &
      number_of_samples, seed)
    character(len=*), intent(in) :: reference_filename, output_directory
    real(Float64), intent(in) :: sampled_f_array(:,:)
    integer(HSIZE_T), intent(in) :: f_array_dimensions(:)
    real(Float64), intent(in) :: denergy, dpitch
    integer(Int64), intent(in) :: number_of_samples
    integer(Int32), intent(in) :: seed

    character(len=4096) :: output_filename
    character(len=32) :: number_of_samples_string
    character(len=32), parameter :: metadata_names(5) = [character(len=32) :: &
      'species', 'atomic_number', 'mass_number', 'charge_state', 'A']
    real(Float64), allocatable :: output_f_array(:,:)
    integer(HID_T) :: reference_id, output_id, f_array_id, f_array_space_id
    integer(HSIZE_T), parameter :: scalar_dimensions(1) = [1_HSIZE_T]
    real(Float64) :: sampled_density
    integer :: error, i, exit_status

    call execute_command_line('mkdir -p "'//trim(output_directory)//'"', &
      exitstat=exit_status)
    if (exit_status /= 0) error stop 'Unable to create output directory'
    output_filename = trim(output_directory)//'/'//basename(reference_filename)

    call h5open_f(error)
    call check_hdf5(error, 'initializing HDF5', output_filename)
    call h5fopen_f(trim(reference_filename), H5F_ACC_RDONLY_F, &
      reference_id, error)
    call check_hdf5(error, 'opening reference file', reference_filename)
    call h5fcreate_f(trim(output_filename), H5F_ACC_TRUNC_F, output_id, error)
    call check_hdf5(error, 'creating sampled file', output_filename)

    call h5ocopy_f(reference_id, 'energy', output_id, 'energy_grid', error)
    call check_hdf5(error, 'copying energy grid', output_filename)
    call h5ocopy_f(reference_id, 'pitch', output_id, 'pitch_grid', error)
    call check_hdf5(error, 'copying pitch grid', output_filename)
    call h5ocopy_f(reference_id, 'r', output_id, 'selected_r', error)
    call check_hdf5(error, 'copying selected R', output_filename)
    call h5ocopy_f(reference_id, 'z', output_id, 'selected_z', error)
    call check_hdf5(error, 'copying selected Z', output_filename)
    do i = 1, size(metadata_names)
      call h5ocopy_f(reference_id, trim(metadata_names(i)), output_id, &
        trim(metadata_names(i)), error)
      call check_hdf5(error, 'copying '//trim(metadata_names(i)), output_filename)
    end do

    call h5screate_simple_f(2, f_array_dimensions, f_array_space_id, error)
    call check_hdf5(error, 'creating f_array dataspace', output_filename)
    call h5dcreate_f(output_id, 'f_array', H5T_NATIVE_DOUBLE, &
      f_array_space_id, f_array_id, error)
    call check_hdf5(error, 'creating sampled f_array', output_filename)

    ! sampled_f_array is indexed as (energy, pitch). The output dataspace uses
    ! the reversed Fortran dimensions (pitch, energy), so transpose the buffer
    ! to preserve f_array[energy,pitch] when h5py reads the resulting file.
    output_f_array = transpose(sampled_f_array)
    call h5dwrite_f(f_array_id, H5T_NATIVE_DOUBLE, output_f_array, &
      f_array_dimensions, error)
    call check_hdf5(error, 'writing sampled f_array', output_filename)
    call h5dclose_f(f_array_id, error)
    call check_hdf5(error, 'closing sampled f_array', output_filename)
    call h5sclose_f(f_array_space_id, error)
    call check_hdf5(error, 'closing f_array dataspace', output_filename)
    call h5ltset_attribute_string_f(output_id, 'f_array', 'units', &
      'ions/(cm^3*keV*dP)', error)
    call check_hdf5(error, 'setting f_array units', output_filename)
    call h5ltset_attribute_string_f(output_id, 'f_array', 'description', &
      'Sampled reconstruction of the fast-ion distribution', error)
    call check_hdf5(error, 'setting f_array description', output_filename)

    ! Store the density represented by the reconstructed energy-pitch array.
    sampled_density = sum(sampled_f_array) * abs(denergy * dpitch)
    call h5ltmake_dataset_double_f(output_id, 'denf', 1, scalar_dimensions, &
      [sampled_density], error)
    call check_hdf5(error, 'creating sampled denf', output_filename)
    call h5ltset_attribute_string_f(output_id, 'denf', 'units', &
      'ions/cm^3', error)
    call check_hdf5(error, 'setting denf units', output_filename)
    call h5ltset_attribute_string_f(output_id, 'denf', 'description', &
      'Density calculated from the sampled energy-pitch distribution', error)
    call check_hdf5(error, 'setting denf description', output_filename)

    call h5ltset_attribute_string_f(output_id, '/', 'description', &
      'Sampled energy-pitch distribution reconstruction', error)
    call check_hdf5(error, 'setting file description', output_filename)
    call h5ltset_attribute_string_f(output_id, '/', 'data_source_type', &
      'sampled_reconstruction', error)
    call check_hdf5(error, 'setting data source type', output_filename)
    call h5ltset_attribute_string_f(output_id, '/', 'data_source_name', &
      trim(reference_filename), error)
    call check_hdf5(error, 'setting data source name', output_filename)
    call h5ltset_attribute_string_f(output_id, '/', 'reference_file', &
      trim(reference_filename), error)
    call check_hdf5(error, 'setting reference file', output_filename)
    write(number_of_samples_string, '(i0)') number_of_samples
    call h5ltset_attribute_string_f(output_id, '/', 'n_samples', &
      trim(number_of_samples_string), error)
    call check_hdf5(error, 'setting n_samples', output_filename)
    call h5ltset_attribute_int_f(output_id, '/', 'rng_seed', [seed], &
      int(1, HSIZE_T), error)
    call check_hdf5(error, 'setting rng_seed', output_filename)

    call h5fclose_f(output_id, error)
    call check_hdf5(error, 'closing sampled file', output_filename)
    call h5fclose_f(reference_id, error)
    call check_hdf5(error, 'closing reference file', reference_filename)
    call h5close_f(error)
    call check_hdf5(error, 'closing HDF5', output_filename)

    write(output_unit, '(a,a)') '  wrote sampled file: ', trim(output_filename)
  end subroutine write_sampled_distribution

  function basename(path) result(name)
    character(len=*), intent(in) :: path
    character(len=4096) :: name
    integer :: slash

    slash = scan(trim(path), '/', back=.true.)
    name = path(slash+1:len_trim(path))
  end function basename

  subroutine check_hdf5(error, operation, filename)
    integer, intent(in) :: error
    character(len=*), intent(in) :: operation, filename

    if (error < 0) then
      write(*, '(a)') 'HDF5 error while '//trim(operation)//': '//trim(filename)
      error stop 'Unable to inspect reference HDF5 file'
    end if
  end subroutine check_hdf5

end module test_003_hdf5
