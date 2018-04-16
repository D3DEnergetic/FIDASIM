module fidampi

  integer, parameter, private   :: Float64 = 8
  integer, private :: num_ranks, my_rank

  interface fidampi_sum
     module procedure fidampi_sum_d0, fidampi_sum_d1, fidampi_sum_d2, &
                      fidampi_sum_d3, fidampi_sum_d4, fidampi_sum_d5,&
                      fidampi_sum_i0, fidampi_sum_i1, fidampi_sum_i2
  end interface

contains

  subroutine fidampi_init()
    use mpi
    implicit none

    integer :: provided, ierr

#ifdef _OMP
    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
#else
    call MPI_INIT(ierr)
#endif

   num_ranks = 1
   my_rank = 0
   if (ierr/=0) then
     write(*,*) "MPI initialization failed, assuming single MPI process"
   else
     call MPI_COMM_SIZE(MPI_COMM_WORLD,num_ranks,ierr)
     if (num_ranks>1) then
       call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
     endif
   endif

  end subroutine


  recursive function fidampi_my_rank() result (n)
    n = my_rank
  end function

  recursive function fidampi_num_ranks() result (n)
    n = num_rank
  end function

  recursive subroutine fidampi_sum_d0(A)
    use mpi
    implicit none

    real(Float64), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = 1

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_d1(A)
    use mpi
    implicit none  

    real(Float64), dimension(:), intent(inout) :: A
    
    integer :: sizeA,ierr

    sizeA = size(A,1)

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_d2(A)
    use mpi
    implicit none

    real(Float64), dimension(:,:), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = size(A,1)*size(A,2)

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_d3(A)
    use mpi
    implicit none

    real(Float64), dimension(:,:,:), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = size(A,1)*size(A,2)*size(A,3)

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_d4(A)
    use mpi
    implicit none

    real(Float64), dimension(:,:,:,:), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = size(A,1)*size(A,2)*size(A,3)*size(A,4)

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_d5(A)
    use mpi
    implicit none

    real(Float64), dimension(:,:,:,:,:), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = size(A,1)*size(A,2)*size(A,3)*size(A,4)*size(A,5)

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_DOUBLE,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_i0(A)
    use mpi
    implicit none

    integer, intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = 1

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_INTEGER,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_i1(A)
    use mpi
    implicit none

    integer, dimension(:), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = size(A,1)

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_INTEGER,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

  recursive subroutine fidampi_sum_i2(A)
    use mpi
    implicit none

    integer, dimension(:,:), intent(inout) :: A

    integer :: sizeA,ierr

    sizeA = size(A,1)*size(A,2)

    if (num_ranks>1) then
       call MPI_Allreduce(A,A,sizeA,MPI_INTEGER,MPI_Sum,MPI_COMM_WORLD,ierr)
    endif ! else nothing to do

  end subroutine

end module
