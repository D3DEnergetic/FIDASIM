module mpi_utils

  integer, parameter, private   :: Float64 = 8
  real, parameter, private      :: max_bytes = 2*1e9 !!GigaBytes
  integer, private :: numranks, myrank

  interface parallel_sum
     module procedure parallel_sum_d0, parallel_sum_d1, parallel_sum_d2, &
                      parallel_sum_d3, parallel_sum_d4, parallel_sum_d5,&
                      parallel_sum_i0, parallel_sum_i1, parallel_sum_i2
  end interface

contains

  subroutine init_mpi()
    use mpi
    implicit none

    integer :: provided, ierr

#ifdef _OMP
    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
#else
    call MPI_INIT(ierr)
#endif

   numranks = 1
   myrank = 0
   if (ierr/=0) then
     write(*,*) "MPI initialization failed, assuming single MPI process"
   else
     call MPI_COMM_SIZE(MPI_COMM_WORLD,numranks,ierr)
     if (numranks>1) then
       call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
     endif
   endif

  end subroutine

  subroutine cleanup_mpi()
    use mpi
    implicit none

    integer :: ierr

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call MPI_FINALIZE(ierr)
  end subroutine

  recursive function my_rank() result (n)
    n = myrank
  end function

  recursive function num_ranks() result (n)
    n = numranks
  end function

  recursive subroutine parallel_sum_d0(A)
    use mpi
    implicit none

    real(Float64), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = 1

    if (numranks>1) then
       call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_d1(A)
    use mpi
    use iso_c_binding
    implicit none

    real(Float64), dimension(:), intent(inout) :: A

    integer :: sizeA,h,ierr
    integer(C_SIZE_T) :: nbytes

    sizeA = size(A,1)

    if (numranks>1) then
        nbytes = Float64*sizeA
        if (nbytes.gt.max_bytes) then
            h = int(sizeA/2)
            call parallel_sum_d1(A(1:h))
            call parallel_sum_d1(A((h+1):sizeA))
        else
            call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
        endif
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_d2(A)
    use mpi
    use iso_c_binding
    implicit none

    real(Float64), dimension(:,:), target, intent(inout) :: A

    integer :: sizeA,h,ierr
    integer(C_SIZE_T) :: nbytes
    real(Float64), pointer :: A_ptr(:)

    sizeA = size(A,1)*size(A,2)

    if (numranks>1) then
        nbytes = Float64*sizeA
        if (nbytes.gt.max_bytes) then
            h = int(sizeA/2)
            call c_f_pointer(c_loc(A), A_ptr, [sizeA])
            call parallel_sum_d1(A_ptr(1:h))
            call parallel_sum_d1(A_ptr((h+1):sizeA))
        else
            call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
        endif
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_d3(A)
    use mpi
    use iso_c_binding
    implicit none

    real(Float64), dimension(:,:,:), target, intent(inout) :: A

    integer :: sizeA,h,ierr
    integer(C_SIZE_T) :: nbytes
    real(Float64), pointer :: A_ptr(:)

    sizeA = size(A,1)*size(A,2)*size(A,3)

    if (numranks>1) then
        nbytes = Float64*sizeA
        if (nbytes.gt.max_bytes) then
            h = int(sizeA/2)
            call c_f_pointer(c_loc(A), A_ptr, [sizeA])
            call parallel_sum_d1(A_ptr(1:h))
            call parallel_sum_d1(A_ptr((h+1):sizeA))
        else
            call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
        endif
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_d4(A)
    use mpi
    use iso_c_binding
    implicit none

    real(Float64), dimension(:,:,:,:), target, intent(inout) :: A

    integer :: sizeA,h,ierr
    integer(C_SIZE_T) :: nbytes
    real(Float64), pointer :: A_ptr(:)

    sizeA = size(A,1)*size(A,2)*size(A,3)*size(A,4)

    if (numranks>1) then
        nbytes = Float64*sizeA
        if (nbytes.gt.max_bytes) then
            h = int(sizeA/2)
            call c_f_pointer(c_loc(A), A_ptr, [sizeA])
            call parallel_sum_d1(A_ptr(1:h))
            call parallel_sum_d1(A_ptr((h+1):sizeA))
        else
            call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
        endif
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_d5(A)
    use mpi
    use iso_c_binding
    implicit none

    real(Float64), dimension(:,:,:,:,:), target, intent(inout) :: A

    integer :: sizeA,h,ierr
    integer(C_SIZE_T) :: nbytes
    real(Float64), pointer :: A_ptr(:)

    sizeA = size(A,1)*size(A,2)*size(A,3)*size(A,4)*size(A,5)

    if (numranks>1) then
        nbytes = Float64*sizeA
        if (nbytes.gt.max_bytes) then
            h = int(sizeA/2)
            call c_f_pointer(c_loc(A), A_ptr, [sizeA])
            call parallel_sum_d1(A_ptr(1:h))
            call parallel_sum_d1(A_ptr((h+1):sizeA))
        else
            call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
        endif
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_i0(A)
    use mpi
    implicit none

    integer, intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = 1

    if (numranks>1) then
       call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_INTEGER,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_i1(A)
    use mpi
    use iso_c_binding
    implicit none

    integer, dimension(:), intent(inout) :: A

    integer :: sizeA,h,ierr
    integer(C_SIZE_T) :: nbytes

    sizeA = size(A,1)

    if (numranks>1) then
        nbytes = kind(1)*sizeA
        if (nbytes.gt.max_bytes) then
            h = int(sizeA/2)
            call parallel_sum_i1(A(1:h))
            call parallel_sum_i1(A((h+1):sizeA))
        else
            call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_INTEGER,MPI_Sum,MPI_COMM_WORLD,ierr)
        endif
    endif ! else nothing to do

  end subroutine

  recursive subroutine parallel_sum_i2(A)
    use mpi
    use iso_c_binding
    implicit none

    integer, dimension(:,:), target, intent(inout) :: A

    integer :: sizeA,h,ierr
    integer(C_SIZE_T) :: nbytes
    integer, pointer :: A_ptr(:)

    sizeA = size(A,1)*size(A,2)

    if (numranks>1) then
        nbytes = kind(1)*sizeA
        if (nbytes.gt.max_bytes) then
            h = int(sizeA/2)
            call c_f_pointer(c_loc(A), A_ptr, [sizeA])
            call parallel_sum_i1(A_ptr(1:h))
            call parallel_sum_i1(A_ptr((h+1):sizeA))
        else
            call MPI_Allreduce(MPI_IN_PLACE,A,sizeA,MPI_INTEGER,MPI_Sum,MPI_COMM_WORLD,ierr)
        endif
    endif ! else nothing to do

  end subroutine

end module
