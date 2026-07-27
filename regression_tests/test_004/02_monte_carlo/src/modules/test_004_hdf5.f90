module test_004_hdf5
  use iso_fortran_env, only: Float64 => real64, output_unit
  use hdf5
  use h5lt
  use test_004_types, only: DistributionCase
  use test_004_hdf5_utils, only: &
    read_real_vector, &
    read_single_real_value, &
    read_scalar_real, &
    read_scalar_integer, &
    read_array_dimensions, &
    check_hdf5
  implicit none
  private

  public :: read_distribution, print_distribution, release_distribution

contains

  subroutine read_distribution(filename, distribution)
    !+ Read one Test 002 smooth energy-pitch distribution and its numerical
    !+ isotope metadata from an HDF5 file.
    character(len=*), intent(in) :: filename
      !+ Path to the Test 002 Stage 2 HDF5 file.
    type(DistributionCase), intent(out) :: distribution
      !+ Distribution data and metadata populated from the file.

    integer(HID_T) :: file_id
    integer :: error

    distribution%filename = trim(filename)

    call h5open_f(error)
    call check_hdf5(error, 'initializing HDF5', filename)

    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
    call check_hdf5(error, 'opening distribution file', filename)

    ! Read the distribution data:
    call read_real_vector(file_id, 'energy', distribution%energy, filename)
    call read_real_vector(file_id, 'pitch', distribution%pitch, filename)
    call read_distribution_array(file_id, distribution, filename)
    call read_single_real_value( &
      file_id, 'denf', distribution%denf, filename)
    call read_single_real_value( &
      file_id, 'r', distribution%selected_r, filename)
    call read_single_real_value( &
      file_id, 'z', distribution%selected_z, filename)

    ! Read the distribution metadata:
    call read_scalar_integer(file_id, 'atomic_number', &
      distribution%atomic_number, filename)
    call read_scalar_integer(file_id, 'mass_number', &
      distribution%mass_number, filename)
    call read_scalar_integer(file_id, 'charge_state', &
      distribution%charge_state, filename)
    call read_scalar_real(file_id, 'A', distribution%atomic_mass, filename)

    ! Close the HDF5 file and clean up:
    call h5fclose_f(file_id, error)
    call check_hdf5(error, 'closing distribution file', filename)
    call h5close_f(error)
    call check_hdf5(error, 'closing HDF5', filename)
  end subroutine read_distribution

  subroutine read_distribution_array(file_id, distribution, filename)
    !+ Read the single-location Test 002 distribution into a rank-two array.
    integer(HID_T), intent(in) :: file_id
      !+ Identifier of the open Test 002 Stage 2 HDF5 file.
    type(DistributionCase), intent(inout) :: distribution
      !+ Distribution object to be updated with the rank-two F(E,p) array. Its
      !+ loaded energy and pitch grids define the expected array dimensions.
    character(len=*), intent(in) :: filename
      !+ Source filename used for error reporting.

    integer(HSIZE_T), allocatable :: dimensions(:)
    real(Float64), allocatable :: stored_distribution(:,:,:,:)
    integer :: error, nenergy, npitch

    ! Test 002 writes the h5py-visible order (z, r, pitch, energy). The
    ! Fortran HDF5 interface presents those dimensions here in reverse order
    ! as (energy, pitch, r, z).
    nenergy = size(distribution%energy)
    npitch = size(distribution%pitch)
    call read_array_dimensions(file_id, 'f', dimensions, filename)
    if (size(dimensions) /= 4) then
      error stop 'The Test 002 f dataset must be four-dimensional'
    end if
    if (dimensions(1) /= nenergy .or. dimensions(2) /= npitch .or. &
        dimensions(3) /= 1 .or. dimensions(4) /= 1) then
      write(output_unit, '(a,4(1x,i0))') 'Unexpected f dimensions:', dimensions
      error stop 'The Test 002 f dataset has an unexpected shape'
    end if

    ! Read the complete stored array, then remove the singleton spatial
    ! dimensions used by the general FIDASIM distribution schema.
    allocate(stored_distribution(nenergy, npitch, 1, 1))
    call h5ltread_dataset_double_f( &
      file_id, 'f', stored_distribution, dimensions, error)
    call check_hdf5(error, 'reading f', filename)

    ! Populate the distribution object with the rank-two F(E,p) array:
    allocate(distribution%f(nenergy, npitch))
    distribution%f = stored_distribution(:,:,1,1)

  end subroutine read_distribution_array

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
    write(output_unit, '(a,i0)') '  atomic number: ', &
      distribution%atomic_number
    write(output_unit, '(a,i0)') '  mass number: ', distribution%mass_number
    write(output_unit, '(a,i0)') '  charge state: ', distribution%charge_state
    write(output_unit, '(a,es14.6)') '  atomic mass [amu]: ', &
      distribution%atomic_mass
  end subroutine print_distribution

  subroutine release_distribution(distribution)
    !+ Release arrays owned by a loaded distribution after processing a case.
    type(DistributionCase), intent(inout) :: distribution
      !+ Distribution object whose energy, pitch, and F(E,p) arrays are freed.

    if (allocated(distribution%energy)) deallocate(distribution%energy)
    if (allocated(distribution%pitch)) deallocate(distribution%pitch)
    if (allocated(distribution%f)) deallocate(distribution%f)
  end subroutine release_distribution

end module test_004_hdf5
