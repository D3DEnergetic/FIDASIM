!+This file contains HDF5 helper routines for writing compressed data files
MODULE hdf5_utils
    !+ A library for writing compressed HDF5 files

USE H5LT
USE HDF5

IMPLICIT NONE

public :: h5ltmake_compressed_dataset_double_f
public :: h5ltmake_compressed_dataset_double_f_1
public :: h5ltmake_compressed_dataset_double_f_2
public :: h5ltmake_compressed_dataset_double_f_3
public :: h5ltmake_compressed_dataset_double_f_4
public :: h5ltmake_compressed_dataset_double_f_5
public :: h5ltmake_compressed_dataset_double_f_6
public :: h5ltmake_compressed_dataset_double_f_7
public :: h5ltmake_compressed_dataset_int_f
public :: h5ltmake_compressed_dataset_int_f_1
public :: h5ltmake_compressed_dataset_int_f_2
public :: h5ltmake_compressed_dataset_int_f_3
public :: h5ltmake_compressed_dataset_int_f_4
public :: h5ltmake_compressed_dataset_int_f_5
public :: h5ltmake_compressed_dataset_int_f_6
public :: h5ltmake_compressed_dataset_int_f_7
public :: h5ltread_dataset_int_scalar_f
public :: h5ltread_dataset_double_scalar_f
public :: check_compression_availability

integer, parameter, private   :: long      = 4 !bytes = 32 bits (-2,147,483,648 to 2,147,483,647)
integer, parameter, private   :: long64    = 8 !bytes = 64 bits (-9,223,372,036,854,775,808 to 9,223,372,036,854,775,807)
integer, parameter, private   :: float     = 4 !bytes = 32 bits (1.2E-38 to 3.4E+38) at 6 decimal places
integer, parameter, private   :: double    = 8 !bytes = 64 bits (2.3E-308 to 1.7E+308) at 15 decimal places
logical, private :: compress_data = .True.

interface h5ltmake_compressed_dataset_double_f
    !+ Write a compressed datasets of 64-bit floats
    module procedure h5ltmake_compressed_dataset_double_f_1
    module procedure h5ltmake_compressed_dataset_double_f_2
    module procedure h5ltmake_compressed_dataset_double_f_3
    module procedure h5ltmake_compressed_dataset_double_f_4
    module procedure h5ltmake_compressed_dataset_double_f_5
    module procedure h5ltmake_compressed_dataset_double_f_6
    module procedure h5ltmake_compressed_dataset_double_f_7
end interface

interface h5ltmake_compressed_dataset_int_f
    !+ Write a compressed dataset of 32-bit integers
    module procedure h5ltmake_compressed_dataset_int_f_1
    module procedure h5ltmake_compressed_dataset_int_f_2
    module procedure h5ltmake_compressed_dataset_int_f_3
    module procedure h5ltmake_compressed_dataset_int_f_4
    module procedure h5ltmake_compressed_dataset_int_f_5
    module procedure h5ltmake_compressed_dataset_int_f_6
    module procedure h5ltmake_compressed_dataset_int_f_7
end interface

contains

subroutine check_compression_availability
    !+ Checks whether dataset compression is available

    IMPLICIT NONE

    logical :: shuffle_avail, gzip_avail
    integer :: gzip_info, shuf_info, filter_info_both
    integer :: error

    call h5open_f(error)

    filter_info_both = ior(H5Z_FILTER_ENCODE_ENABLED_F,H5Z_FILTER_DECODE_ENABLED_F)

    !! Check for GZIP filter
    call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, gzip_avail, error)
    call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, gzip_info, error)
    if(.not.gzip_avail) then
        print*,'HDF5: gzip filter is not available'
        compress_data = .False.
    endif
    if (filter_info_both .ne. gzip_info) then
        print*,'HDF5: gzip filter is not available for encoding and decoding'
        compress_data = .False.
    endif

    !! Check for SHUFFLE filter
    call h5zfilter_avail_f(H5Z_FILTER_SHUFFLE_F, shuffle_avail, error)
    call h5zget_filter_info_f(H5Z_FILTER_SHUFFLE_F, shuf_info, error)
    if(.not.shuffle_avail) then
        print*,'HDF5: shuffle filter is not available'
        compress_data = .False.
    endif
    if (filter_info_both .ne. shuf_info) then
        print*,'HDF5: shuffle filter is not available for encoding and decoding'
        compress_data = .False.
    endif

    if(.not.compress_data) then
        print*,'HDF5: Compression is not available. Proceeding without compression.'
    endif

    call h5close_f(error)

end subroutine check_compression_availability

subroutine h5ltread_dataset_int_scalar_f(loc_id, dset_name, x, error)
    !+ Write a scalar 32-bit integer

    IMPLICIT NONE

    integer(HID_T), intent(in)   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in) :: dset_name
        !+ Name of the dataset to create
    integer, intent(inout)       :: x
        !+ Data to be written to the dataset
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HSIZE_T), dimension(1) :: dims(1) = 1
    integer, dimension(1) :: dummy

    call h5ltread_dataset_int_f(loc_id, dset_name, dummy, dims, error)
    x = dummy(1)

end subroutine h5ltread_dataset_int_scalar_f

subroutine h5ltread_dataset_double_scalar_f(loc_id, dset_name, x, error)
    !+ Write a scalar 64-bit float
    IMPLICIT NONE

    integer(HID_T), intent(in)   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in) :: dset_name
        !+ Name of the dataset to create
    real(double), intent(inout)  :: x
        !+ Data to be written to the dataset
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HSIZE_T), dimension(1) :: dims(1) = 1
    real(double), dimension(1) :: dummy

    call h5ltread_dataset_double_f(loc_id, dset_name, dummy, dims, error)
    x = dummy(1)

end subroutine h5ltread_dataset_double_scalar_f

!Compressed Doubles
subroutine h5ltmake_compressed_dataset_double_f_1(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 64-bit float dataset of dimension 1

    IMPLICIT NONE

    integer(hid_t), intent(in)                   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                 :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                          :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)   :: dims
        !+ Array of the size of each dimension
    real(double), dimension(dims(1)), intent(in) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                         :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id
    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_1

subroutine h5ltmake_compressed_dataset_double_f_2(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 64-bit float dataset of dimension 2

    IMPLICIT NONE

    integer(hid_t), intent(in)                           :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                         :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                  :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)           :: dims
        !+ Array of the size of each dimension
    real(double), dimension(dims(1),dims(2)), intent(in) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                 :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_2

subroutine h5ltmake_compressed_dataset_double_f_3(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 64-bit float dataset of dimension 3

    IMPLICIT NONE

    integer(hid_t), intent(in)                                   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                                 :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                          :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)                   :: dims
        !+ Array of the size of each dimension
    real(double), dimension(dims(1),dims(2),dims(3)), intent(in) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                         :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_3

subroutine h5ltmake_compressed_dataset_double_f_4(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 64-bit float dataset of dimension 4

    IMPLICIT NONE

    integer(hid_t), intent(in)                     :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                   :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                            :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)     :: dims
        !+ Array of the size of each dimension
    real(double), intent(in), &
        dimension(dims(1),dims(2),dims(3),dims(4)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                           :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_4

subroutine h5ltmake_compressed_dataset_double_f_5(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 64-bit float dataset of dimension 5

    IMPLICIT NONE

    integer(hid_t), intent(in)                           :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                         :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                  :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)           :: dims
        !+ Array of the size of each dimension
    real(double), intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                 :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_5

subroutine h5ltmake_compressed_dataset_double_f_6(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 64-bit float dataset of dimension 6

    IMPLICIT NONE

    integer(hid_t), intent(in)                                   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                                 :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                          :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)                   :: dims
        !+ Array of the size of each dimension
    real(double), intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                         :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_6

subroutine h5ltmake_compressed_dataset_double_f_7(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 64-bit float dataset of dimension 7

    IMPLICIT NONE

    integer(hid_t), intent(in)                                           :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                                         :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                                  :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)                           :: dims
        !+ Array of the size of each dimension
    real(double), intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                                 :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_7

!Compressed Integers
subroutine h5ltmake_compressed_dataset_int_f_1(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 32-bit integer dataset of dimension 1

    IMPLICIT NONE

    integer(hid_t), intent(in)                 :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)               :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                        :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in) :: dims
        !+ Array of the size of each dimension
    integer, dimension(dims(1)), intent(in)    :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                       :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_1

subroutine h5ltmake_compressed_dataset_int_f_2(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 32-bit integer dataset of dimension 2

    IMPLICIT NONE

    integer(hid_t), intent(in)                      :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                    :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                             :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)      :: dims
        !+ Array of the size of each dimension
    integer, dimension(dims(1),dims(2)), intent(in) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                            :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_2

subroutine h5ltmake_compressed_dataset_int_f_3(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 32-bit integer dataset of dimension 3

    IMPLICIT NONE

    integer(hid_t), intent(in)                              :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                            :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                     :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)              :: dims
        !+ Array of the size of each dimension
    integer, dimension(dims(1),dims(2),dims(3)), intent(in) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                    :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_3

subroutine h5ltmake_compressed_dataset_int_f_4(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 32-bit integer dataset of dimension 4

    IMPLICIT NONE

    integer(hid_t), intent(in)                                      :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                                    :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                             :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)                      :: dims
        !+ Array of the size of each dimension
    integer, dimension(dims(1),dims(2),dims(3),dims(4)), intent(in) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                            :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id
    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif
end subroutine h5ltmake_compressed_dataset_int_f_4

subroutine h5ltmake_compressed_dataset_int_f_5(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 32-bit integer dataset of dimension 5

    IMPLICIT NONE

    integer(hid_t), intent(in)                           :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                         :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                  :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)           :: dims
        !+ Array of the size of each dimension
    integer, intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                 :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif
end subroutine h5ltmake_compressed_dataset_int_f_5

subroutine h5ltmake_compressed_dataset_int_f_6(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 32-bit integer dataset of dimension 6

    IMPLICIT NONE

    integer(hid_t), intent(in)                                   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                                 :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                          :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)                   :: dims
        !+ Array of the size of each dimension
    integer, intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                         :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_6

subroutine h5ltmake_compressed_dataset_int_f_7(loc_id,&
           dset_name, rank, dims, buf, error )
    !+ Write a compressed 32-bit integer dataset of dimension 7

    IMPLICIT NONE

    integer(hid_t), intent(in)                                           :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in)                                         :: dset_name
        !+ Name of the dataset to create
    integer, intent(in)                                                  :: rank
        !+ Number of dimensions of dataspace
    integer(HSIZE_T), dimension(*), intent(in)                           :: dims
        !+ Array of the size of each dimension
    integer, intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                                 :: error
        !+ HDF5 error code

    integer(HID_T) :: did, sid, plist_id

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        call h5pset_shuffle_f(plist_id, error)
        call h5pset_deflate_f(plist_id, 9, error)
        call h5pset_chunk_f(plist_id, rank, dims, error)

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_7

END MODULE hdf5_utils
