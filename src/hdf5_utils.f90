!+This file contains HDF5 helper routines for writing compressed data files
MODULE hdf5_utils
    !+ A library for writing compressed HDF5 files

USE H5LT
USE HDF5
USE H5DS

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
public :: h5_set_dimension_name
public :: h5_make_dimension_scale
public :: h5_attach_dimension_scale

integer, parameter, private   :: Int32   = 4 !bytes = 32 bits (-2,147,483,648 to 2,147,483,647)
integer, parameter, private   :: Int64   = 8 !bytes = 64 bits (-9,223,372,036,854,775,808 to 9,223,372,036,854,775,807)
integer, parameter, private   :: Float32 = 4 !bytes = 32 bits (1.2E-38 to 3.4E+38) at 6 decimal places
integer, parameter, private   :: Float64 = 8 !bytes = 64 bits (2.3E-308 to 1.7E+308) at 15 decimal places
logical, private :: compress_data = .True.
logical, private :: default_compress = .True.
integer, parameter, private   :: default_compression_level = 4

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
    real(Float64), intent(inout)  :: x
        !+ Data to be written to the dataset
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HSIZE_T), dimension(1) :: dims(1) = 1
    real(Float64), dimension(1) :: dummy

    call h5ltread_dataset_double_f(loc_id, dset_name, dummy, dims, error)
    x = dummy(1)

end subroutine h5ltread_dataset_double_scalar_f

subroutine chunk_size(elsize, dims, cdims)
    integer, intent(in)                          :: elsize
        !+ Size of elements in bytes
    integer(HSIZE_T), dimension(*), intent(in)   :: dims
        !+ Dimensions of dataset
    integer(HSIZE_T), dimension(:), intent(out)  :: cdims
        !+ Maximum allowed chunk size/dims

    real, parameter :: max_bytes = 4*1e9 !GigaBytes

    integer :: d
    real(Float64) :: nbytes


    d = size(cdims)
    cdims(1:d) = dims(1:d)
    nbytes = elsize*product(1.d0*cdims)
    do while ((nbytes.gt.max_bytes).and.(d.gt.0))
        cdims(d) = max(floor(cdims(d)*max_bytes/nbytes,Int32),1)
        nbytes = elsize*product(1.d0*cdims)
        d = d - 1
    enddo

end subroutine chunk_size

!Compressed Doubles
subroutine h5ltmake_compressed_dataset_double_f_1(loc_id,&
           dset_name, rank, dims, buf, error, compress, level )
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
    real(Float64), dimension(dims(1)), intent(in) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                         :: error
        !+ HDF5 error code
    logical, intent(in), optional                :: compress
        !+ Flag to compress
    integer, intent(in), optional                :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(1)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Float64, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_1

subroutine h5ltmake_compressed_dataset_double_f_2(loc_id,&
           dset_name, rank, dims, buf, error, compress, level )
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
    real(Float64), dimension(dims(1),dims(2)), intent(in):: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                 :: error
        !+ HDF5 error code
    logical, intent(in), optional                        :: compress
        !+ Flag to compress
    integer, intent(in), optional                        :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(2)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Float64, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_2

subroutine h5ltmake_compressed_dataset_double_f_3(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    real(Float64), dimension(dims(1),dims(2),dims(3)), intent(in):: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                         :: error
        !+ HDF5 error code
    logical, intent(in), optional                                :: compress
        !+ Flag to compress
    integer, intent(in), optional                                :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(3)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Float64, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_3

subroutine h5ltmake_compressed_dataset_double_f_4(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    real(Float64), intent(in), &
        dimension(dims(1),dims(2),dims(3),dims(4)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                           :: error
        !+ HDF5 error code
    logical, intent(in), optional                  :: compress
        !+ Flag to compress
    integer, intent(in), optional                  :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(4)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Float64, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_4

subroutine h5ltmake_compressed_dataset_double_f_5(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    real(Float64), intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                 :: error
        !+ HDF5 error code
    logical, intent(in), optional                        :: compress
        !+ Flag to compress
    integer, intent(in), optional                        :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(5)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Float64, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_5

subroutine h5ltmake_compressed_dataset_double_f_6(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    real(Float64), intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                         :: error
        !+ HDF5 error code
    logical, intent(in), optional                                :: compress
        !+ Flag to compress
    integer, intent(in), optional                                :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(6)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Float64, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_DOUBLE, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_DOUBLE, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_double_f_6

subroutine h5ltmake_compressed_dataset_double_f_7(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    real(Float64), intent(in), &
      dimension(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)) :: buf
        !+ Buffer with data to be written to the dataset
    integer, intent(out)                                                 :: error
        !+ HDF5 error code
    logical, intent(in), optional                                        :: compress
        !+ Flag to compress
    integer, intent(in), optional                                        :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(7)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Float64, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

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
           dset_name, rank, dims, buf, error, compress, level  )
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
    logical, intent(in), optional              :: compress
        !+ Flag to compress
    integer, intent(in), optional              :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(1)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Int32, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_1

subroutine h5ltmake_compressed_dataset_int_f_2(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    logical, intent(in), optional                   :: compress
        !+ Flag to compress
    integer, intent(in), optional                   :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(2)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Int32, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_2

subroutine h5ltmake_compressed_dataset_int_f_3(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    logical, intent(in), optional                           :: compress
        !+ Flag to compress
    integer, intent(in), optional                           :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(3)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Int32, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_3

subroutine h5ltmake_compressed_dataset_int_f_4(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    logical, intent(in), optional                                   :: compress
        !+ Flag to compress
    integer, intent(in), optional                                   :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(4)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Int32, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif
end subroutine h5ltmake_compressed_dataset_int_f_4

subroutine h5ltmake_compressed_dataset_int_f_5(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    logical, intent(in), optional                        :: compress
        !+ Flag to compress
    integer, intent(in), optional                        :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(5)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Int32, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif
end subroutine h5ltmake_compressed_dataset_int_f_5

subroutine h5ltmake_compressed_dataset_int_f_6(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    logical, intent(in), optional                                :: compress
        !+ Flag to compress
    integer, intent(in), optional                                :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(6)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Int32, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_6

subroutine h5ltmake_compressed_dataset_int_f_7(loc_id,&
           dset_name, rank, dims, buf, error, compress, level  )
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
    logical, intent(in), optional                                        :: compress
        !+ Flag to compress
    integer, intent(in), optional                                        :: level
        !+ Compression level

    integer(HID_T) :: did, sid, plist_id
    integer(HSIZE_T) :: cdims(7)

    logical :: do_chunk
    integer :: c

    do_chunk=default_compress
    if(present(compress)) then
        do_chunk = compress
    endif

    c = default_compression_level
    if(present(level)) then
        c = level
    endif

    if(.not.compress_data) then
        call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims, buf, error)
    else
        call h5screate_simple_f(rank, dims, sid, error)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

        if(do_chunk) then
            call h5pset_shuffle_f(plist_id, error)
            call h5pset_deflate_f(plist_id, c, error)

            call chunk_size(Int32, dims, cdims)
            call h5pset_chunk_f(plist_id, rank, cdims, error)
        endif

        call h5dcreate_f(loc_id, dset_name, H5T_NATIVE_INTEGER, sid, &
             did, error, dcpl_id=plist_id)

        call h5dwrite_f(did, H5T_NATIVE_INTEGER, buf, dims, error)

        call h5sclose_f(sid, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(did, error)
    endif

end subroutine h5ltmake_compressed_dataset_int_f_7

subroutine h5_set_dimension_name(loc_id, dset_name, dim_idx, dim_name, error)
    !+ Set a label/name for a specific dimension of a dataset
    !+ This uses the HDF5 Dimension Scale API (H5DSset_label)

    IMPLICIT NONE

    integer(HID_T), intent(in)   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in) :: dset_name
        !+ Name of the dataset
    integer, intent(in)          :: dim_idx
        !+ Dimension index (1-based Fortran indexing)
    character(len=*), intent(in) :: dim_name
        !+ Name/label for the dimension
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HID_T) :: dset_id

    ! Open the dataset
    call h5dopen_f(loc_id, dset_name, dset_id, error)
    if (error .ne. 0) return

    ! Set the dimension label (HDF5 Fortran API uses 1-based indexing)
    call h5dsset_label_f(dset_id, dim_idx, dim_name, error)

    ! Close the dataset
    call h5dclose_f(dset_id, error)

end subroutine h5_set_dimension_name

subroutine h5_make_dimension_scale(loc_id, dset_name, scale_name, error)
    !+ Convert a dataset into a dimension scale
    !+ This marks a 1D coordinate array as a dimension scale that can be attached to other datasets

    IMPLICIT NONE

    integer(HID_T), intent(in)   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in) :: dset_name
        !+ Name of the dataset to convert to a dimension scale
    character(len=*), intent(in) :: scale_name
        !+ Name for this dimension scale (can be same as dset_name)
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HID_T) :: dset_id

    ! Open the dataset
    call h5dopen_f(loc_id, dset_name, dset_id, error)
    if (error .ne. 0) return

    ! Make it a dimension scale
    call h5dsset_scale_f(dset_id, error, scale_name)

    ! Close the dataset
    call h5dclose_f(dset_id, error)

end subroutine h5_make_dimension_scale

subroutine h5_attach_dimension_scale(loc_id, dset_name, scale_name, dim_idx, error)
    !+ Attach a dimension scale to a specific dimension of a dataset
    !+ This creates the link between a data array and its coordinate array

    IMPLICIT NONE

    integer(HID_T), intent(in)   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in) :: dset_name
        !+ Name of the dataset to attach scale to
    character(len=*), intent(in) :: scale_name
        !+ Name of the dimension scale dataset
    integer, intent(in)          :: dim_idx
        !+ Dimension index (1-based Fortran indexing)
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HID_T) :: dset_id, scale_id

    ! Open the dataset
    call h5dopen_f(loc_id, dset_name, dset_id, error)
    if (error .ne. 0) return

    ! Open the dimension scale dataset
    call h5dopen_f(loc_id, scale_name, scale_id, error)
    if (error .ne. 0) then
        call h5dclose_f(dset_id, error)
        return
    endif

    ! Attach the scale to the dimension (HDF5 Fortran API uses 1-based indexing)
    call h5dsattach_scale_f(dset_id, scale_id, dim_idx, error)

    ! Close both datasets
    call h5dclose_f(scale_id, error)
    call h5dclose_f(dset_id, error)

end subroutine h5_attach_dimension_scale

subroutine h5_create_soft_link(loc_id, target_path, link_name, error)
    !+ Create a soft (symbolic) link to an existing dataset
    !+ This allows coordinates to appear in multiple groups without data duplication

    use hdf5
    IMPLICIT NONE

    integer(HID_T), intent(in)   :: loc_id
        !+ HDF5 file or group identifier
    character(len=*), intent(in) :: target_path
        !+ Path to the target dataset (e.g., "/grid/x")
    character(len=*), intent(in) :: link_name
        !+ Name of the link to create (e.g., "x")
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HID_T) :: lcpl_id, lapl_id

    ! Create link creation property list
    call h5pcreate_f(H5P_LINK_CREATE_F, lcpl_id, error)
    if (error .ne. 0) return

    ! Create link access property list
    call h5pcreate_f(H5P_LINK_ACCESS_F, lapl_id, error)
    if (error .ne. 0) then
        call h5pclose_f(lcpl_id, error)
        return
    endif

    ! Create the soft link
    call h5lcreate_soft_f(target_path, loc_id, link_name, error, lcpl_id, lapl_id)

    ! Close property lists
    call h5pclose_f(lapl_id, error)
    call h5pclose_f(lcpl_id, error)

end subroutine h5_create_soft_link

subroutine h5_create_coordinate_links(group_id, grid_path, error)
    !+ Create soft links to grid coordinates within a group
    !+ This makes the group self-contained for xarray compatibility

    IMPLICIT NONE

    integer(HID_T), intent(in)   :: group_id
        !+ HDF5 group identifier where links will be created
    character(len=*), intent(in) :: grid_path
        !+ Path to grid group (e.g., "/grid")
    integer, intent(out)         :: error
        !+ HDF5 error code

    ! Create soft links for standard grid coordinates
    call h5_create_soft_link(group_id, trim(grid_path)//"/x", "x", error)
    if (error .ne. 0) return

    call h5_create_soft_link(group_id, trim(grid_path)//"/y", "y", error)
    if (error .ne. 0) return

    call h5_create_soft_link(group_id, trim(grid_path)//"/z", "z", error)

end subroutine h5_create_coordinate_links

subroutine h5_attach_local_scales(group_id, dset_name, dim_names, ndims, error)
    !+ Attach dimension scales that are local to the group
    !+ Assumes coordinate arrays exist as soft links in the same group

    use H5DS
    IMPLICIT NONE

    integer(HID_T), intent(in)   :: group_id
        !+ HDF5 group identifier
    character(len=*), intent(in) :: dset_name
        !+ Name of the dataset (relative to group)
    character(len=*), intent(in) :: dim_names(:)
        !+ Names of dimensions to attach
    integer, intent(in)          :: ndims
        !+ Number of dimensions
    integer, intent(out)         :: error
        !+ HDF5 error code

    integer(HID_T) :: dset_id, scale_id
    integer :: i
    character(len=128) :: coord_name

    ! Open the dataset
    call h5dopen_f(group_id, dset_name, dset_id, error)
    if (error .ne. 0) return

    ! Attach each dimension scale
    do i = 1, ndims
        coord_name = trim(dim_names(i))

        ! Skip if dimension name is empty or not a coordinate
        if (len_trim(coord_name) == 0) cycle
        if (coord_name == 'atomic_level' .or. &
            coord_name == 'beam_component' .or. &
            coord_name == 'particle' .or. &
            coord_name == 'channel') cycle

        ! Try to open and attach the coordinate
        call h5dopen_f(group_id, coord_name, scale_id, error)
        if (error == 0) then
            call h5dsattach_scale_f(dset_id, scale_id, i, error)
            call h5dclose_f(scale_id, error)
        endif
    enddo

    ! Close the dataset
    call h5dclose_f(dset_id, error)
    error = 0  ! Reset error since some scales might not exist

end subroutine h5_attach_local_scales

END MODULE hdf5_utils
