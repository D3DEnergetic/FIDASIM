!+This file contains routines for parallel random number generation and a basic
!+sparse array implementation
module utilities
    !+ Utilities for parallel random number generation and sparse arrays
#ifdef _OMP
  use omp_lib
#endif

implicit none
private
public :: ind2sub, sub2ind, time, cumsum
public :: rng_type, rng_init, rng_seed, get_rng, rng, randind_cdf
public :: rng_uniform, rng_normal, randu, randn, randind
public :: SparseArray, get_value, sparse
public :: deriv
#ifdef _DEF_INTR
public :: norm2
#endif

integer, parameter :: Int32 = 4
integer, parameter :: Int64 = kind(int8(1))
integer, parameter :: Float32 = kind(1.e0)
integer, parameter :: Float64 = kind(1.d0)

integer(Int32), parameter :: IA = 16807
integer(Int32), parameter :: IM = 2147483647
integer(Int32), parameter :: IQ = 127773
integer(Int32), parameter :: IR = 2836

integer, parameter :: ns = 2

type :: rng_type
    !+ Random Number Generator Derived Type
    integer(Int32) :: seed
    integer(Int32), dimension(ns) :: state
end type rng_type

type(rng_type), dimension(:), allocatable :: rng

type SparseArray
    integer :: nnz = 0
        !+ Number of non-zero elements
    integer :: nd = 0
        !+ Number of dimensions
    integer, dimension(:), allocatable :: dims
        !+ Dimensions of array
    integer, dimension(:), allocatable :: inds
        !+ Linear index of non-zero elements
    real(Float64), dimension(:), allocatable :: vals
        !+ Array values
end type SparseArray

interface randu
    module procedure randu_arr
    module procedure randu_r_arr
end interface

interface randn
    module procedure randn_arr
    module procedure randn_r_arr
end interface

interface randind_cdf
    !+ Procedure for generating a random array index/subscripts
    module procedure randind_r_cdf_1
    module procedure randind_cdf_1
end interface

interface randind
    !+ Procedure for generating a random array index/subscripts
    module procedure randind_n
    module procedure randind_w_1
    module procedure randind_w_2
    module procedure randind_w_3
    module procedure randind_r_n
    module procedure randind_r_w_1
    module procedure randind_r_w_2
    module procedure randind_r_w_3
end interface

interface search_sorted_first
    !+ Function for searching a sorted array
    module procedure search_sorted_first_integer
    module procedure search_sorted_first_float64
end interface

interface sparse
    !+ Creates a sparse array from a dense array
    module procedure sparse_1
    module procedure sparse_2
    module procedure sparse_3
    module procedure sparse_4
end interface

interface deriv
    !+ Procedure for finding derivatives from an array
    module procedure deriv_1d
    module procedure deriv_2d
end interface

contains

!============================================================================
!---------------------------Array Indexing Routines--------------------------
!============================================================================
subroutine ind2sub(dims,ind,subs)
    !+ Calculate the subscripts `subs` into an array with dimensions `dims`
    !+ given the corresponding linear index `ind`
    integer, dimension(:), intent(in)  :: dims
        !+ Dimensions of array
    integer, intent(in)                :: ind
        !+ Linear index
    integer, dimension(:), intent(out) :: subs
        !+ Subscripts corresponding to the linear index

    integer :: i, ndims, ind1, ind2

    ind1=ind
    ndims = size(dims)
    do i=1,ndims-1
        ind2 = (ind1-1)/dims(i) + 1
        subs(i) = ind1 - dims(i)*(ind2-1)
        ind1 = ind2
    enddo
    subs(ndims) = ind1

end subroutine ind2sub

function sub2ind(dims, subs) result (ind)
    !+ Calculates the linear index of an array with dimensions `dims` and
    !+ subcripts `subs`
    integer, dimension(:), intent(in) :: dims
        !+ Dimension of Array
    integer, dimension(:), intent(in) :: subs
        !+ Subscripts to convert
    integer :: ind
        !+ Linear index

    integer :: k, l, p

    ind = subs(1)
    do k=2,size(dims)
        p = dims(1)
        do l=2, k-1
            p = p*dims(l)
        enddo
        ind = ind + p*(subs(k)-1)
    enddo

end function sub2ind

!============================================================================
!-------------------------------Search Routines------------------------------
!============================================================================
function search_sorted_first_integer(x, v) result(hi)
    !+ Returns the index of the first value in `x` greater or equal to `v`.
    !+ Returns `length(x)+1` if `v` is greater then all values in `x`.
    integer, dimension(:), intent(in) :: x
        !+ Monotonically increasing array
    integer, intent(in)               :: v
        !+ Value to search

    integer :: hi, lo, m

    lo = 0
    hi = size(x)+1

    do while (lo.lt.(hi-1))
        m = rshift(lo+hi,1)
        if(x(m).lt.v) then
            lo = m
        else
            hi = m
        endif
    enddo

end function search_sorted_first_integer

function search_sorted_first_float64(x, v) result(hi)
    !+ Returns the index of the first value in `x` greater or equal to `v`.
    !+ Returns `length(x)+1` if `v` is greater then all values in `x`.
    real(Float64), dimension(:), intent(in) :: x
        !+ Monotonically increasing array
    real(Float64), intent(in)               :: v
        !+ Value to search

    integer :: hi, lo, m

    lo = 0
    hi = size(x)+1

    do while (lo.lt.(hi-1))
        m = rshift(lo+hi,1)
        if(x(m).lt.v) then
            lo = m
        else
            hi = m
        endif
    enddo

end function search_sorted_first_float64

!============================================================================
!-----------------------Parallel Random Number Routines----------------------
!============================================================================
function rng_seed() result (seed)
    !+ Generates random 32-bit integer seed from `/dev/urandom`
    integer(Int32) :: seed
        !+ Seed value

    open(89, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old")
    read(89) seed
    close(89)
    seed = abs(seed)

end function rng_seed

subroutine rng_init(self, seed)
#ifdef _MPI
    use mpi_utils
#endif

    !+ Procedure to initialize a random number generator with a seed.
    !+ If seed is negative then random seed is used
    type(rng_type), intent(inout) :: self
        !+ Random Number Generator
    integer(Int32), intent(in)    :: seed
        !+ Initial Seed Value

    integer(Int32) :: s

    if(seed.lt.0) then
        s = rng_seed()
    else
#ifdef _MPI
        s = seed + my_rank()
#else
        s = seed
#endif
    endif


    self%seed = s
    self%state(1) = ieor(777755555,abs(s))
    self%state(2) = ior(ieor(888889999,abs(s)),1)

end subroutine rng_init

function get_rng() result(r)
    type(rng_type) :: r
    integer :: thread_id

#ifdef _OMP
    thread_id = OMP_get_thread_num() + 1
#else
    thread_id = 1
#endif

    r = rng(thread_id)

end function get_rng

subroutine update_rng(r)
    type(rng_type), intent(in) :: r
    integer :: thread_id

#ifdef _OMP
    thread_id = OMP_get_thread_num() + 1
#else
    thread_id = 1
#endif

    rng(thread_id) = r

end subroutine update_rng

function rng_uniform(self) result(u)
    !+ Generate a uniformally-distributed random number in the range [0,1)
    type(rng_type), intent(inout) :: self
        !+ Random Number Generator
    real(Float64)                 :: u
        !+ Uniform random deviate

    integer(Int32) :: ix,iy,k
    real(Float64) :: am 

    ix = self%state(1)
    iy = self%state(2)

    ix = ieor(ix,ishft(ix,13))
    ix = ieor(ix,ishft(ix,-17))
    ix = ieor(ix,ishft(ix,5))
    k=iy/IQ
    iy=IA*(iy - k*IQ) - IR*k
    if(iy.lt.0) iy = iy + IM
    self%state(1) = ix
    self%state(2) = iy

    am = nearest(1.0,-1.0)/IM
    u = am*ior(iand(IM,ieor(ix,iy)),1)

end function rng_uniform

function rng_normal(self) result(n)
    !+ Generate a normally-distributed random number with mean 0 and standard deviation 1
    type(rng_type), intent(inout) :: self
        !+ Random Number Generator
    real(Float64)                 :: n
        !+ Normal random deviate

    real(Float64), parameter :: s = 0.449871d0
    real(Float64), parameter :: t = 0.386595d0
    real(Float64), parameter :: a = 0.196000d0
    real(Float64), parameter :: b = 0.254720d0
    real(Float64), parameter :: r1 = 0.27597d0
    real(Float64), parameter :: r2 = 0.27846d0
    real(Float64) :: u, v, x, y, q

    do
        u = rng_uniform(self)
        v = rng_uniform(self)
        v = 1.7156d0 * (v - 0.5d0)

        x = u - s
        y = abs(v) + t
        q = x**2 + y*(a*y - b*x)

        if (q.lt.r1) exit
        if (q.gt.r2) cycle
        if((v**2).lt.(-4.0*log(u)*u**2)) exit
    enddo

    n = v/u

end function rng_normal

subroutine randu_r_arr(r, randomu)
    !+ Generate an array of uniformally-distributed random deviates
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:), intent(out) :: randomu
        !+ Array of uniform random deviates

    integer :: i

    randomu = 0.d0
    do i=1,size(randomu)
        randomu(i) = rng_uniform(r)
    enddo

end subroutine randu_r_arr

subroutine randu_arr(randomu)
    !+ Generate an array of uniformally-distributed random deviates
    real(Float64), dimension(:), intent(out) :: randomu
        !+ Array of uniform random deviates

    type(rng_type) ::r

    r = get_rng()
    call randu_r_arr(r, randomu)
    call update_rng(r)

end subroutine randu_arr

subroutine randn_r_arr(r, randomn)
    !+ Generate an array of normally-distributed random deviates
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:), intent(out) :: randomn
        !+ Array of normal random deviates

    integer :: i

    randomn = 0.d0
    do i=1,size(randomn)
        randomn(i) = rng_normal(r)
    enddo

end subroutine randn_r_arr

subroutine randn_arr(randomn)
    !+ Generate an array of normally-distributed random deviates
    real(Float64), dimension(:), intent(out) :: randomn
        !+ Array of normal random deviates

    type(rng_type) :: r

    r = get_rng()
    call randn_r_arr(r,randomn)
    call update_rng(r)

end subroutine randn_arr

subroutine randind_r_n(r, n, randomi)
    !+ Generate a array of uniformally-distributed random integers in the range [1, n]
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    integer, intent(in)                :: n
        !+ Largest possible value
    integer, dimension(:), intent(out) :: randomi
        !+ Array of uniform deviates

    integer :: i
    real(Float64), dimension(1) :: randomu

    randomi = 0
    do i=1,size(randomi)
        call randu_r_arr(r, randomu)
        randomi(i) = ceiling(randomu(1)*n)
    enddo

end subroutine randind_r_n

subroutine randind_n(n, randomi)
    !+ Generate a array of uniformally-distributed random integers in the range [1, n]
    integer, intent(in)                :: n
        !+ Largest possible value
    integer, dimension(:), intent(out) :: randomi
        !+ Array of uniform deviates

    type(rng_type) :: r

    r = get_rng()
    call randind_r_n(r, n, randomi)
    call update_rng(r)

end subroutine randind_n

subroutine randind_r_cdf_1(r, cdf, randomi)
    !+ Generate an array of random indices of an 1D array distributed according to `cdf`
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:), intent(in) :: cdf
        !+ 1D array of index weights
    integer, dimension(:), intent(out)      :: randomi
        !+ Random indices

    integer :: i, n
    real(Float64) :: cdf_val
    real(Float64), dimension(1) :: randomu

    n = size(cdf)

    randomi = 0
    do i=1,size(randomi)
        call randu_r_arr(r,randomu)
        cdf_val = randomu(1)*cdf(n)
        randomi(i) = min(search_sorted_first(cdf,cdf_val),n)
    enddo

end subroutine randind_r_cdf_1

subroutine randind_cdf_1(cdf, randomi)
    !+ Generate an array of random indices of an 1D array distributed according to `cdf`
    real(Float64), dimension(:), intent(in) :: cdf
        !+ 1D array of index weights
    integer, dimension(:), intent(out)      :: randomi
        !+ Random indices

    type(rng_type) :: r

    r = get_rng()
    call randind_r_cdf_1(r, cdf, randomi)
    call update_rng(r)

end subroutine randind_cdf_1

subroutine cumsum(x, cs)
    !+ Calculate cumulative sum
    real(Float64), dimension(:), intent(in)  :: x
        !+ Array to sum
    real(Float64), dimension(:), intent(out) :: cs
        !+ Cumulative sum of `x`

    integer :: i, n
    real(Float64) :: cdf_val, t

    n = size(x)
    t = 0.d0
    do i=1, n
        cs(i) = t + x(i)
        t = cs(i)
    enddo

end subroutine cumsum

subroutine randind_r_w_1(r, w, randomi)
    !+ Generate an array of random indices of an 1D array distributed according to `w`
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:), intent(in) :: w
        !+ 1D array of index weights
    integer, dimension(:), intent(out)      :: randomi
        !+ Random indices

    real(Float64), dimension(size(w)) :: cdf

    call cumsum(w, cdf)
    call randind_r_cdf_1(r, cdf, randomi)

end subroutine randind_r_w_1

subroutine randind_w_1(w, randomi)
    !+ Generate an array of random indices of an 1D array distributed according to `w`
    real(Float64), dimension(:), intent(in) :: w
        !+ 1D array of index weights
    integer, dimension(:), intent(out)      :: randomi
        !+ Random indices

    type(rng_type) :: r

    r = get_rng()
    call randind_r_w_1(r, w, randomi)
    call update_rng(r)

end subroutine randind_w_1

subroutine randind_r_w_2(r, w, randomi)
    !+ Generate an array of random subscripts of an 2D array distributed according to `w`
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:,:), intent(in) :: w
        !+ 2D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    integer :: i,nw
    integer, dimension(2) :: subs
    integer, dimension(size(randomi,2)) :: randi

    randomi = 0
    nw = size(w)

    call randind_r_w_1(r, reshape(w,[nw]),randi)
    do i=1,size(randomi,2)
        call ind2sub(shape(w),randi(i),subs)
        randomi(:,i) = subs
    enddo

end subroutine randind_r_w_2

subroutine randind_w_2(w, randomi)
    !+ Generate an array of random subscripts of an 2D array distributed according to `w`
    real(Float64), dimension(:,:), intent(in) :: w
        !+ 2D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    type(rng_type) :: r

    r = get_rng()
    call randind_r_w_2(r, w, randomi)
    call update_rng(r)

end subroutine randind_w_2

subroutine randind_r_w_3(r, w, randomi)
    !+ Generate an array of random subscripts of an 3D array distributed according to `w`
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:,:,:), intent(in) :: w
        !+ 3D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    integer :: i,nw
    integer, dimension(3) :: subs
    integer, dimension(size(randomi,2)) :: randi

    randomi = 0
    nw = size(w)

    call randind_r_w_1(r,reshape(w,[nw]),randi)
    do i=1,size(randomi,2)
        call ind2sub(shape(w),randi(i),subs)
        randomi(:,i) = subs
    enddo

end subroutine randind_r_w_3

subroutine randind_w_3(w, randomi)
    !+ Generate an array of random subscripts of an 3D array distributed according to `w`
    real(Float64), dimension(:,:,:), intent(in) :: w
        !+ 3D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    type(rng_type) :: r

    r = get_rng()
    call randind_r_w_3(r, w, randomi)
    call update_rng(r)

end subroutine randind_w_3

!============================================================================
!------------------------------Sparse Routines-------------------------------
!============================================================================
subroutine sparse_1(A,SA)
    !+ Routine to create a 1D sparse array from a 1D dense array
    real(Float64), dimension(:), intent(in) :: A
        !+ Dense Array
    type(SparseArray), intent(out) :: SA
        !+ Sparse Array

    integer :: n,i,c

    SA%nd = 1
    allocate(SA%dims(SA%nd))
    SA%dims = shape(A)

    SA%nnz = count(A.ne.0.d0)
    if(SA%nnz.eq.0) return
    allocate(SA%vals(SA%nnz),SA%inds(SA%nnz))

    n = size(A)
    c = 1
    do i=1,n
        if(A(i).ne.0.d0)then
            SA%inds(c) = i
            SA%vals(c) = A(i)
            c = c + 1
        endif
        if(c.gt.SA%nnz) exit
    enddo

end subroutine sparse_1

subroutine sparse_2(A,SA)
    !+ Routine to create a 2D sparse array from a 2D dense array
    real(Float64), dimension(:,:), intent(in) :: A
        !+ Dense Array
    type(SparseArray), intent(out) :: SA
        !+ Sparse Array

    integer :: subs(2)
    integer :: n,i,c

    SA%nd = 2
    allocate(SA%dims(SA%nd))
    SA%dims = shape(A)

    SA%nnz = count(A.ne.0.d0)
    if(SA%nnz.eq.0) return
    allocate(SA%vals(SA%nnz),SA%inds(SA%nnz))

    n = size(A)
    c = 1
    do i=1,n
        call ind2sub(SA%dims,i,subs)
        if(A(subs(1),subs(2)).ne.0.d0)then
            SA%inds(c) = i
            SA%vals(c) = A(subs(1),subs(2))
            c = c + 1
        endif
        if(c.gt.SA%nnz) exit
    enddo

end subroutine sparse_2

subroutine sparse_3(A,SA)
    !+ Routine to create a 3D sparse array from a 3D dense array
    real(Float64), dimension(:,:,:), intent(in) :: A
        !+ Dense Array
    type(SparseArray), intent(out) :: SA
        !+ Sparse Array

    integer :: subs(3)
    integer :: n,i,c

    SA%nd = 3
    allocate(SA%dims(SA%nd))
    SA%dims = shape(A)

    SA%nnz = count(A.ne.0.d0)
    if(SA%nnz.eq.0) return
    allocate(SA%vals(SA%nnz),SA%inds(SA%nnz))

    n = size(A)
    c = 1
    do i=1,n
        call ind2sub(SA%dims,i,subs)
        if(A(subs(1),subs(2),subs(3)).ne.0.d0)then
            SA%inds(c) = i
            SA%vals(c) = A(subs(1),subs(2),subs(3))
            c = c + 1
        endif
        if(c.gt.SA%nnz) exit
    enddo

end subroutine sparse_3

subroutine sparse_4(A,SA)
    !+ Routine to create a 4D sparse array from a 4D dense array
    real(Float64), dimension(:,:,:,:), intent(in) :: A
        !+ Dense Array
    type(SparseArray), intent(out) :: SA
        !+ Sparse Array

    integer :: subs(4)
    integer :: n,i,c

    SA%nd = 4
    allocate(SA%dims(SA%nd))
    SA%dims = shape(A)

    SA%nnz = count(A.ne.0.d0)
    if(SA%nnz.eq.0) return
    allocate(SA%vals(SA%nnz),SA%inds(SA%nnz))

    n = size(A)
    c = 1
    do i=1,n
        call ind2sub(SA%dims,i,subs)
        if(A(subs(1),subs(2),subs(3),subs(4)).ne.0.d0)then
            SA%inds(c) = i
            SA%vals(c) = A(subs(1),subs(2),subs(3),subs(4))
            c = c + 1
        endif
        if(c.gt.SA%nnz) exit
    enddo

end subroutine sparse_4

function get_value(SA, subs) result (val)
    !+ Gets value of sparse array `SA` at the subscripts `subs`
    type(SparseArray), intent(in)     :: SA
        !+ Sparse Array
    integer, dimension(:), intent(in) :: subs
        !+ Subscripts of Sparse Array
    real(Float64) :: val
        !+ Value of `SA` at `subs`

    integer :: ind, cind

    val = 0.d0
    if(SA%nnz.eq.0) return

    ind = sub2ind(SA%dims, subs)
    cind = search_sorted_first(SA%inds, ind)
    if(ind.eq.SA%inds(cind))then
        val = SA%vals(cind)
    endif

end function get_value

!============================================================================
!--------------------------------Deriv Routines------------------------------
!============================================================================

subroutine deriv_1d(x,y,yp)
    !+ Uses 3 point lagrangian method to calculate the derivative of an array
    real(Float64), dimension(:),intent(in)  :: x
        !+ X Values
    real(Float64), dimension(:),intent(in)  :: y
        !+ Y Values
    real(Float64), dimension(:),intent(out) :: yp
        !+ Derivative of Y w.r.t. X

    integer :: i,n
        !! temporary values for loops
    real(Float64) :: p1,p2,p3
        !! intermeadiate values for 3 point lagrangian

    n = size(x)-1
    do i = 2,n
        p1 = x(i-1)
        p2 = x(i)
        p3 = x(i+1)
        yp(i) = (y(i-1)*(p2-p3)/((p1-p2)*(p1-p3))) + &
                (y(i)*((1/(p2-p3))-(1/(p1-p2)))) -   &
                (y(i+1)*(p1-p2)/((p1-p3)*(p2-p3)))
    enddo

    yp(1) = (y(1)*((x(1)-x(2))+(x(1)-x(3)))/((x(1)-x(2))*(x(1)-x(3)))) - &
            (y(2)*(x(1)-x(3))/((x(1)-x(2))*(x(2)-x(3)))) +               &
            (y(3)*(x(1)-x(2))/((x(1)-x(3))*(x(2)-x(3))))
    yp(n+1) = -(y(n-1)*(x(n)-x(n+1))/((x(n-1)-x(n))*(x(n-1)-x(n+1)))) +   &
              (y(n)*(x(n-1)-x(n+1))/((x(n-1)-x(n))*(x(n)-x(n+1)))) -     &
              (y(n+1)*((x(n-1)-x(n+1))+(x(n)-x(n+1)))/((x(n-1)-x(n+1)) * &
              (x(n)-x(n+1))))

end subroutine deriv_1d

subroutine deriv_2d(x,y,z,zxp,zyp)
    !+ Uses 3 point lagrangian method to calculate the partial derivative
    !+ of an array Z w.r.t X and Y
    real(Float64), dimension(:), intent(in)  :: x
        !+ X Values
    real(Float64), dimension(:), intent(in)  :: y
        !+ Y Values
    real(Float64), dimension(:,:), intent(in)  :: z
        !+ Z Values
    real(Float64), dimension(:,:), intent(out) :: zxp
        !+ Derivative of Z w.r.t. X
    real(Float64), dimension(:,:), intent(out) :: zyp
        !+ Derivative of Z w.r.t. Y

    integer :: i,n
        !! temporary values for loops

    n = size(y)
    do i = 1,n
        call deriv_1d(x,z(:,i),zxp(:,i))
    enddo

    n = size(x)
    do i = 1,n
        call deriv_1d(y,z(i,:),zyp(i,:))
    enddo

end subroutine deriv_2d

!============================================================================
!------------------------------ Misc. Routines ------------------------------
!============================================================================
function time(time_start) result (time_str)
    !+ Returns time string
    integer, dimension(8), intent(in), optional :: time_start
        !+ Optional start time
    character(30) :: time_str
        !+ Time string

    integer :: ts(8), ta(8), hour, minu, sec

    ts = 0
    if(present(time_start)) then
        ts = time_start
    endif

    call date_and_time(values=ta)
    hour = ta(5) - ts(5)
    minu = ta(6) - ts(6)
    sec  = ta(7) - ts(7)
    if (minu.lt.0.) then
        minu = minu + 60
        hour = hour - 1
    endif
    if (sec.lt.0.) then
        sec  = sec + 60
        minu = minu - 1
    endif

    if(present(time_start)) then
        write(time_str,'(I2,":",I2.2,":",I2.2," --- elapsed:",I2,":",I2.2,":",I2.2)') &
            ta(5),ta(6),ta(7), hour,minu,sec
    else
        write(time_str,'(I2,":",I2.2,":",I2.2)') &
            ta(5),ta(6),ta(7)
    endif

end function time

#ifdef _DEF_INTR
! define missing intrinsics

function norm2( in ) result ( res )
  implicit none
  real(Float64),dimension(:) :: in
  real(Float64) :: res
  res = sqrt(sum( in(:)**2 ))
end function norm2

#endif

end module utilities
