!+This file contains routines for parallel random number generation and a basic
!+sparse array implementation
module utilities
    !+ Utilities for parallel random number generation and sparse arrays
#ifdef _OMP
  use omp_lib
#endif

use iso_c_binding
implicit none
private
public :: ind2sub, sub2ind, time_string, cumsum
public :: rng_type, rng_init, rng_seed, get_rng, rng, randind_cdf
public :: rng_uniform, rng_normal, randu, randn, randind
public :: SparseArray, get_value, sparse
public :: deriv
public :: InterpolCoeffs1D, InterpolCoeffs2D, InterpolCoeffs3D
public :: interpol, interpol_coeff
#ifdef _DEF_INTR
public :: norm2
#endif

!============================================================================
!------------------------------- Parameters ---------------------------------
!============================================================================

integer, parameter :: Int32 = 4
integer, parameter :: Int64 = kind(int8(1))
integer, parameter :: Float32 = kind(1.e0)
integer, parameter :: Float64 = kind(1.d0)

integer(Int32), parameter :: IA = 16807
integer(Int32), parameter :: IM = 2147483647
integer(Int32), parameter :: IQ = 127773
integer(Int32), parameter :: IR = 2836

integer, parameter :: ns = 2

!============================================================================
!---------------------------------- Types -----------------------------------
!============================================================================

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

type InterpolCoeffs1D
    !+ Linear Interpolation Coefficients and indices
    integer :: i = 0
        !+ Index of position right before `xout`
    real(Float64) :: b1 = 0.d0
        !+ Coefficient for y(i) term
    real(Float64) :: b2 = 0.d0
        !+ Coefficient for y(i+1) term
end type InterpolCoeffs1D

type InterpolCoeffs2D
    !+ 2D Linear Interpolation Coefficients and indices
    integer :: i = 0
        !+ Index of abscissa before `xout`
    integer :: j = 0
        !+ Index of ordinate before `yout`
    real(Float64) :: b11 = 0.d0
        !+ Coefficient for z(i,j) term
    real(Float64) :: b12 = 0.d0
        !+ Coefficient for z(i,j+1) term
    real(Float64) :: b21 = 0.d0
        !+ Coefficient for z(i+1,j) term
    real(Float64) :: b22 = 0.d0
        !+ Coefficient for z(i+1,j+1) term
end type InterpolCoeffs2D

type InterpolCoeffs3D
    !+ 3D Cylindrical Interpolation Coefficients and indices
    integer :: i = 0
        !+ Index of R before `rout`
    integer :: j = 0
        !+ Index of Z before `zout`
    integer :: k = 0
        !+ Index of Phi before `phiout`
    real(Float64) :: b111 = 0.d0
        !+ Coefficient for z(i,j,k) term
    real(Float64) :: b121 = 0.d0
        !+ Coefficient for z(i,j+1,k) term
    real(Float64) :: b112 = 0.d0
        !+ Coefficient for z(i,j,k+1) term
    real(Float64) :: b122 = 0.d0
        !+ Coefficient for z(i,j+1,k+1) term
    real(Float64) :: b211 = 0.d0
        !+ Coefficient for z(i+1,j,k) term
    real(Float64) :: b212 = 0.d0
        !+ Coefficient for z(i+1,j,k+1) term
    real(Float64) :: b221 = 0.d0
        !+ Coefficient for z(i+1,j+1,k) term
    real(Float64) :: b222 = 0.d0
        !+ Coefficient for z(i+1,j+1,k+1) term
end type InterpolCoeffs3D

!============================================================================
!-------------------------------- Interfaces --------------------------------
!============================================================================

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
    module procedure randind_w_4
    module procedure randind_w_5
    module procedure randind_r_n
    module procedure randind_r_w_1
    module procedure randind_r_w_2
    module procedure randind_r_w_3
    module procedure randind_r_w_4
    module procedure randind_r_w_5
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
    module procedure deriv_3d
end interface

interface interpol_coeff
    !+ Calculates interpolation coefficients
    module procedure interpol1D_coeff, interpol1D_coeff_arr
    module procedure interpol2D_coeff, interpol2D_coeff_arr
    module procedure cyl_interpol3D_coeff, cyl_interpol3D_coeff_arr
end interface

interface interpol
    !+ Performs linear/bilinear/cylindrical interpolation
    module procedure interpol1D_arr
    module procedure interpol2D_arr, interpol2D_2D_arr
    module procedure interpol3D_arr, interpol3D_2D_arr
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
    real(Float64), dimension(:,:), target, intent(in) :: w
        !+ 2D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    integer :: i,nw
    integer, dimension(2) :: subs
    integer, dimension(size(randomi,2)) :: randi
    real(Float64), pointer :: w_ptr(:)

    randomi = 0
    nw = size(w)

    call c_f_pointer(c_loc(w), w_ptr, [nw])
    call randind_r_w_1(r, w_ptr,randi)
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
    real(Float64), dimension(:,:,:), target, intent(in) :: w
        !+ 3D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    integer :: i,nw
    integer, dimension(3) :: subs
    integer, dimension(size(randomi,2)) :: randi
    real(Float64), pointer :: w_ptr(:)

    randomi = 0
    nw = size(w)

    call c_f_pointer(c_loc(w), w_ptr, [nw])
    call randind_r_w_1(r, w_ptr, randi)
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

subroutine randind_r_w_4(r, w, randomi)
    !+ Generate an array of random subscripts of an 4D array distributed according to `w`
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:,:,:,:), target, intent(in) :: w
        !+ 4D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    integer :: i,nw
    integer, dimension(4) :: subs
    integer, dimension(size(randomi,2)) :: randi
    real(Float64), pointer :: w_ptr(:)

    randomi = 0
    nw = size(w)

    call c_f_pointer(c_loc(w), w_ptr, [nw])
    call randind_r_w_1(r, w_ptr, randi)
    do i=1,size(randomi,2)
        call ind2sub(shape(w),randi(i),subs)
        randomi(:,i) = subs
    enddo

end subroutine randind_r_w_4

subroutine randind_w_4(w, randomi)
    !+ Generate an array of random subscripts of an 4D array distributed according to `w`
    real(Float64), dimension(:,:,:,:), intent(in) :: w
        !+ 4D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    type(rng_type) :: r

    r = get_rng()
    call randind_r_w_4(r, w, randomi)
    call update_rng(r)

end subroutine randind_w_4

subroutine randind_r_w_5(r, w, randomi)
    !+ Generate an array of random subscripts of an 5D array distributed according to `w`
    type(rng_type), intent(inout) :: r
        !+ Random Number Generator
    real(Float64), dimension(:,:,:,:,:), target, intent(in) :: w
        !+ 5D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    integer :: i,nw
    integer, dimension(5) :: subs
    integer, dimension(size(randomi,2)) :: randi
    real(Float64), pointer :: w_ptr(:)

    randomi = 0
    nw = size(w)

    call c_f_pointer(c_loc(w), w_ptr, [nw])
    call randind_r_w_1(r, w_ptr, randi)
    do i=1,size(randomi,2)
        call ind2sub(shape(w),randi(i),subs)
        randomi(:,i) = subs
    enddo

end subroutine randind_r_w_5

subroutine randind_w_5(w, randomi)
    !+ Generate an array of random subscripts of an 5D array distributed according to `w`
    real(Float64), dimension(:,:,:,:,:), intent(in) :: w
        !+ 5D array of subscript weights
    integer, dimension(:,:), intent(out)      :: randomi
        !+ A 2D (ndim, :) array of random subscripts

    type(rng_type) :: r

    r = get_rng()
    call randind_r_w_5(r, w, randomi)
    call update_rng(r)

end subroutine randind_w_5

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

subroutine deriv_3d(r,z,phi,f,frp,fzp,fphip)
    !+ Uses 3 point lagrangian method to calculate the partial derivative
    !+ of an array F w.r.t R, Z and Phi
    real(Float64), dimension(:), intent(in)  :: r
        !+ R Values
    real(Float64), dimension(:), intent(in)  :: z
        !+ Z Values
    real(Float64), dimension(:), intent(in)  :: phi
        !+ Phi Values
    real(Float64), dimension(:,:,:), intent(in)  :: f
        !+ F Values
    real(Float64), dimension(:,:,:), intent(out) :: frp
        !+ Derivative of F w.r.t. R
    real(Float64), dimension(:,:,:), intent(out) :: fzp
        !+ Derivative of F w.r.t. Z
    real(Float64), dimension(:,:,:), intent(out) :: fphip
        !+ Derivative of F w.r.t. Phi

    integer :: i,n
        !! temporary values for loops
    if (size(phi) .gt. 2) then

        n = size(phi)
        do i = 1,n
            call deriv_2d(r,z,f(:,:,i),frp(:,:,i),fzp(:,:,i))
        enddo

        n = size(z)
        do i = 1,n
            call deriv_2d(r,phi,f(:,i,:),frp(:,i,:),fphip(:,i,:))
        enddo

        n = size(r)
        do i = 1,n
            call deriv_2d(z,phi,f(i,:,:),fzp(i,:,:),fphip(i,:,:))
        enddo
    else
        fphip = 0.0d0
        call deriv_2d(r,z,f(:,:,1),frp(:,:,1),fzp(:,:,1))

    endif

end subroutine deriv_3d

!============================================================================
!---------------------------Interpolation Routines---------------------------
!============================================================================
subroutine interpol1D_coeff(xmin,dx,nx,xout,c,err)
    !+ Linear interpolation coefficients and index for a 1D grid y(x)
    real(Float64), intent(in)           :: xmin
        !+ Minimum abscissa value
    real(Float64), intent(in)           :: dx
        !+ Absissa spacing
    integer, intent(in)                 :: nx
        !+ Number of abscissa
    real(Float64), intent(in)           :: xout
        !+ Abscissa value to interpolate
    type(InterpolCoeffs1D), intent(out) :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional      :: err
        !+ Error code

    real(Float64) :: x1, xp, b1, b2
    integer :: i, err_status

    err_status = 1
    xp = max(xout,xmin)
    i = floor((xp - xmin)/dx)+1

    if ((i.gt.0).and.(i.le.(nx-1))) then
        x1 = xmin + (i-1)*dx

        b2 = (xp - x1)/dx
        b1 = (1.0 - b2)

        c%i = i
        c%b1 = b1
        c%b2 = b2
        err_status = 0
    endif

    if(present(err)) err = err_status

end subroutine interpol1D_coeff

subroutine interpol1D_coeff_arr(x,xout,c,err)
    !+ Linear interpolation coefficients and index for a 1D grid y(x)
    real(Float64), dimension(:), intent(in) :: x
        !+ Abscissa values
    real(Float64), intent(in)               :: xout
        !+ Abscissa value to interpolate
    type(InterpolCoeffs1D), intent(out)     :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional          :: err
        !+ Error code

    real(Float64) :: xmin, dx
    integer :: sx,err_status

    err_status = 1
    sx = size(x)
    xmin = x(1)
    dx = abs(x(2)-x(1))

    call interpol1D_coeff(xmin, dx, sx, xout, c, err_status)

    if(present(err)) err = err_status

end subroutine interpol1D_coeff_arr

subroutine interpol2D_coeff(xmin,dx,nx,ymin,dy,ny,xout,yout,c,err)
    !+ Bilinear interpolation coefficients and indicies for a 2D grid z(x,y)
    real(Float64), intent(in)           :: xmin
        !+ Minimum abscissa
    real(Float64), intent(in)           :: dx
        !+ Abscissa spacing
    integer, intent(in)                 :: nx
        !+ Number of abscissa
    real(Float64), intent(in)           :: ymin
        !+ Minimum ordinate
    real(Float64), intent(in)           :: dy
        !+ Ordinate spacing
    integer, intent(in)                 :: ny
        !+ Number of ordinates points
    real(Float64), intent(in)           :: xout
        !+ Abscissa value to interpolate
    real(Float64), intent(in)           :: yout
        !+ Ordinate value to interpolate
    type(InterpolCoeffs2D), intent(out) :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional      :: err
        !+ Error code

    real(Float64) :: x1, x2, y1, y2, xp, yp
    integer :: i, j, err_status

    err_status = 1
    xp = max(xout,xmin)
    yp = max(yout,ymin)
    i = floor((xp-xmin)/dx)+1
    j = floor((yp-ymin)/dy)+1

    if (((i.gt.0).and.(i.le.(nx-1))).and.((j.gt.0).and.(j.le.(ny-1)))) then
        x1 = xmin + (i-1)*dx
        x2 = x1 + dx
        y1 = ymin + (j-1)*dy
        y2 = y1 + dy

        c%b11 = ((x2 - xp) * (y2 - yp))/(dx*dy)
        c%b21 = ((xp - x1) * (y2 - yp))/(dx*dy)
        c%b12 = ((x2 - xp) * (yp - y1))/(dx*dy)
        c%b22 = ((xp - x1) * (yp - y1))/(dx*dy)
        c%i = i
        c%j = j
        err_status = 0
    endif

    if(present(err)) err = err_status

end subroutine interpol2D_coeff

subroutine interpol2D_coeff_arr(x,y,xout,yout,c,err)
    !!Bilinear interpolation coefficients and indicies for a 2D grid z(x,y)
    real(Float64), dimension(:), intent(in) :: x
        !+ Abscissa values
    real(Float64), dimension(:), intent(in) :: y
        !+ Ordinate values
    real(Float64), intent(in)               :: xout
        !+ Abscissa value to interpolate
    real(Float64), intent(in)               :: yout
        !+ Ordinate value to interpolate
    type(InterpolCoeffs2D), intent(out)     :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional          :: err
        !+ Error code

    real(Float64) :: xmin, ymin, dx, dy
    integer :: sx, sy, err_status

    err_status = 1
    sx = size(x)
    sy = size(y)
    xmin = x(1)
    ymin = y(1)
    dx = abs(x(2)-x(1))
    dy = abs(y(2)-y(1))

    call interpol2D_coeff(xmin, dx, sx, ymin, dy, sy, xout, yout, c, err_status)

    if(present(err)) err = err_status

end subroutine interpol2D_coeff_arr

subroutine cyl_interpol3D_coeff(rmin,dr,nr,zmin,dz,nz,phimin,dphi,nphi,rout,zout,phiout,c,err)
    !+ Cylindrical interpolation coefficients and indicies for a 3D grid
    real(Float64), intent(in)           :: rmin
        !+ Minimum R
    real(Float64), intent(in)           :: dr
        !+ R spacing
    integer, intent(in)                 :: nr
        !+ Number of R points
    real(Float64), intent(in)           :: zmin
        !+ Minimum Z
    real(Float64), intent(in)           :: dz
        !+ Z spacing
    integer, intent(in)                 :: nz
        !+ Number of Z points
    real(Float64), intent(in)           :: phimin
        !+ Minimum phi
    real(Float64), intent(in)           :: dphi
        !+ Phi spacing
    integer, intent(in)                 :: nphi
        !+ Number of phi points
    real(Float64), intent(in)           :: rout
        !+ R value to interpolate
    real(Float64), intent(in)           :: zout
        !+ Z value to interpolate
    real(Float64), intent(in)           :: phiout
        !+ Phi value to interpolate
    type(InterpolCoeffs3D), intent(out) :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional      :: err
        !+ Error code

    type(InterpolCoeffs2D) :: b
    real(Float64) :: r1, r2, phi1, phi2, z1, z2, rp, phip, zp, dV
    real(Float64) :: phi
    integer :: i, j, k, err_status

    err_status = 1

    rp = max(rout,rmin)
    zp = max(zout,zmin)
    phip = max(phiout,phimin)
    i = floor((rp-rmin)/dr)+1
    j = floor((zp-zmin)/dz)+1
    k = floor((phip-phimin)/dphi)+1

    if (nphi .eq. 1) then
        if (((i.gt.0).and.(i.le.(nr-1))).and.((j.gt.0).and.(j.le.(nz-1)))) then
            call interpol2D_coeff(rmin, dr, nr, zmin, dz, nz, rout, zout, b, err_status)
            c%b111 = b%b11
            c%b121 = b%b12
            c%b221 = b%b22
            c%b211 = b%b21
            c%b212 = 0
            c%b222 = 0
            c%b122 = 0
            c%b112 = 0
            c%i = b%i
            c%j = b%j
            c%k = 1
            err_status = 0
        endif
    else
        if ((((i.gt.0).and.(i.le.(nr-1))).and.((j.gt.0).and.(j.le.(nz-1)))).and.((k.gt.0).and.(k.le.(nphi-1)))) then
            r1 = rmin + (i-1)*dr
            r2 = r1 + dr
            z1 = zmin + (j-1)*dz
            z2 = z1 + dz
            phi1 = phimin + (k-1)*dphi
            phi2 = phi1 + dphi
            dV = ((r2**2 - r1**2) * (phi2 - phi1) * (z2 - z1))

            !! Both volume elements have a factor of 1/2 that cancels out
            c%b111 = ((r2**2 - rp**2) * (phi2 - phip) * (z2 - zp)) / dV
            c%b121 = ((r2**2 - rp**2) * (phi2 - phip) * (zp - z1)) / dV
            c%b221 = ((rp**2 - r1**2) * (phi2 - phip) * (zp - z1)) / dV
            c%b211 = ((rp**2 - r1**2) * (phi2 - phip) * (z2 - zp)) / dV
            c%b212 = ((rp**2 - r1**2) * (phip - phi1) * (z2 - zp)) / dV
            c%b222 = ((rp**2 - r1**2) * (phip - phi1) * (zp - z1)) / dV
            c%b122 = ((r2**2 - rp**2) * (phip - phi1) * (zp - z1)) / dV
            c%b112 = ((r2**2 - rp**2) * (phip - phi1) * (z2 - zp)) / dV
            c%i = i
            c%j = j
            c%k = k
            err_status = 0
        endif
    endif

    if(present(err)) err = err_status

end subroutine cyl_interpol3D_coeff

subroutine cyl_interpol3D_coeff_arr(r,z,phi,rout,zout,phiout,c,err)
    !+ Cylindrical interpolation coefficients and indicies for a 3D grid
    real(Float64), dimension(:), intent(in) :: r
        !+ R values
    real(Float64), dimension(:), intent(in) :: z
        !+ Z values
    real(Float64), dimension(:), intent(in) :: phi
        !+ Phi values
    real(Float64), intent(in)               :: rout
        !+ R value to interpolate
    real(Float64), intent(in)               :: zout
        !+ Z value to interpolate
    real(Float64), intent(in)               :: phiout
        !+ Phi value to interpolate
    type(InterpolCoeffs3D), intent(out)     :: c
        !+ Interpolation Coefficients
    integer, intent(out), optional          :: err
        !+ Error code

    type(InterpolCoeffs2D) :: b
    real(Float64) :: rmin, phimin, zmin, dr, dphi, dz
    integer :: sr, sphi, sz, err_status

    err_status = 1
    sr = size(r)
    sphi = size(phi)
    sz = size(z)

    rmin = r(1)
    zmin = z(1)
    dr = abs(r(2)-r(1))
    dz = abs(z(2)-z(1))

    if (sphi .eq. 1) then
        call interpol2D_coeff(rmin, dr, sr, zmin, dz, sz, rout, zout, b, err_status)
        c%b111 = b%b11
        c%b121 = b%b12
        c%b221 = b%b22
        c%b211 = b%b21
        c%b212 = 0
        c%b222 = 0
        c%b122 = 0
        c%b112 = 0
        c%i = b%i
        c%j = b%j
        c%k = 1
    else
        phimin = phi(1)
        dphi = abs(phi(2)-phi(1))
        call cyl_interpol3D_coeff(rmin, dr, sr, zmin, dz, sz, phimin, dphi, sphi, rout, zout, phiout, c, err_status)
    endif

    if(present(err)) err = err_status

end subroutine cyl_interpol3D_coeff_arr

subroutine interpol1D_arr(x, y, xout, yout, err, coeffs)
    !+ Performs linear interpolation on a uniform 1D grid y(x)
    real(Float64), dimension(:), intent(in)      :: x
        !+ The abscissa values of `y`
    real(Float64), dimension(:), intent(in)      :: y
        !+ Values at abscissa values `x`: y(x)
    real(Float64), intent(in)                    :: xout
        !+ Abscissa value to interpolate
    real(Float64), intent(out)                   :: yout
        !+ Interpolant: y(xout)
    integer, intent(out), optional               :: err
        !+ Error code
    type(InterpolCoeffs1D), intent(in), optional :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    type(InterpolCoeffs1D) :: c
    integer :: i, err_status

    err_status = 1
    if(present(coeffs)) then
        c = coeffs
        err_status = 0
    else
        call interpol_coeff(x,xout,c,err_status)
    endif

    if(err_status.eq.0) then
        i = c%i
        yout = c%b1*y(i) + c%b2*y(i+1)
    else
        yout = 0.d0
    endif

    if(present(err)) err = err_status

end subroutine interpol1D_arr

subroutine interpol2D_arr(x, y, z, xout, yout, zout, err, coeffs)
    !+ Performs bilinear interpolation on a 2D grid z(x,y)
    real(Float64), dimension(:), intent(in)   :: x
        !+ The abscissa values of `z`
    real(Float64), dimension(:), intent(in)   :: y
        !+ The ordinate values of `z`
    real(Float64), dimension(:,:), intent(in) :: z
        !+ Values at the abscissa/ordinates: z(x,y)
    real(Float64), intent(in)                 :: xout
        !+ The abscissa value to interpolate
    real(Float64), intent(in)                 :: yout
        !+ The ordinate value to interpolate
    real(Float64), intent(out)                :: zout
        !+ Interpolant: z(xout,yout)
    integer, intent(out), optional            :: err
        !+ Error code
    type(InterpolCoeffs2D), intent(in), optional :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    type(InterpolCoeffs2D) :: c
    integer :: i, j, err_status

    err_status = 1
    if(present(coeffs)) then
        c = coeffs
        err_status = 0
    else
        call interpol_coeff(x,y,xout,yout,c,err_status)
    endif

    if(err_status.eq.0) then
        i = c%i
        j = c%j
        zout = c%b11*z(i,j) + c%b12*z(i,j+1) + c%b21*z(i+1,j) + c%b22*z(i+1,j+1)
    else
        zout = 0.d0
    endif

    if(present(err)) err = err_status

end subroutine interpol2D_arr

subroutine interpol2D_2D_arr(x, y, z, xout, yout, zout, err, coeffs)
    !+ Performs bilinear interpolation on a 2D grid of 2D arrays z(:,:,x,y)
    real(Float64), dimension(:), intent(in)       :: x
        !+ The abscissa values of `z`
    real(Float64), dimension(:), intent(in)       :: y
        !+ The ordinate values of `z`
    real(Float64), dimension(:,:,:,:), intent(in) :: z
        !+ Values at the abscissa/ordinates: z(:,:,x,y)
    real(Float64), intent(in)                     :: xout
        !+ The abscissa value to interpolate
    real(Float64), intent(in)                     :: yout
        !+ The ordinate value to interpolate
    real(Float64), dimension(:,:), intent(out)    :: zout
        !+ Interpolant: z(:,:,xout,yout)
    integer, intent(out), optional                :: err
        !+ Error code
    type(InterpolCoeffs2D), intent(in), optional :: coeffs
        !+ Precomputed Linear Interpolation Coefficients

    type(InterpolCoeffs2D) :: c
    integer :: i, j, err_status

    err_status = 1
    if(present(coeffs)) then
        c = coeffs
        err_status = 0
    else
        call interpol_coeff(x,y,xout,yout,c,err_status)
    endif

    if(err_status.eq.0) then
        i = c%i
        j = c%j
        zout = c%b11*z(:,:,i,j) + c%b12*z(:,:,i,j+1) + c%b21*z(:,:,i+1,j) + c%b22*z(:,:,i+1,j+1)
    else
        zout = 0.0
    endif

    if(present(err)) err = err_status

end subroutine interpol2D_2D_arr

subroutine interpol3D_arr(r, z, phi, d, rout, zout, phiout, dout, err, coeffs)
    !+ Performs cylindrical interpolation on a 3D grid f(r,z,phi)
    real(Float64), dimension(:), intent(in) :: r
        !+ R values
    real(Float64), dimension(:), intent(in) :: z
        !+ Z values
    real(Float64), dimension(:), intent(in) :: phi
        !+ Phi values
    real(Float64), dimension(:,:,:), intent(in) :: d
        !+ Values at r,z,phi: d(r,z,phi)
    real(Float64), intent(in)               :: rout
        !+ R value to interpolate
    real(Float64), intent(in)               :: zout
        !+ Z value to interpolate
    real(Float64), intent(in)               :: phiout
        !+ Phi value to interpolate
    real(Float64), intent(out)                :: dout
        !+ Interpolant: d(rout,zout,phiout)
    integer, intent(out), optional          :: err
        !+ Error code
    type(InterpolCoeffs3D), intent(in), optional :: coeffs
        !+ Precomputed Interpolation Coefficients

    type(InterpolCoeffs3D) :: b
    integer :: i, j, k, k2, err_status
    integer :: nphi

    err_status = 1

    nphi = size(phi)

    if(present(coeffs)) then
        b = coeffs
        if(nphi .eq. 1) then
            b%b212 = 0
            b%b222 = 0
            b%b122 = 0
            b%b112 = 0
            b%k = 1
        endif
        err_status = 0
    else
        call interpol_coeff(r,z,phi,rout,zout,phiout,b,err_status)
    endif


    if(err_status.eq.0) then
        i = b%i
        j = b%j
        k = b%k
        if(nphi .eq. 1) then
            k2 = min(k+1,nphi)
        else
            k2 = k+1
        endif

        dout = b%b111*d(i,j,k)    + b%b121*d(i,j+1,k) +   &
               b%b112*d(i,j,k2)   + b%b122*d(i,j+1,k2) +  &
               b%b211*d(i+1,j,k)  + b%b221*d(i+1,j+1,k) + &
               b%b212*d(i+1,j,k2) + b%b222*d(i+1,j+1,k2)
    else
        dout = 0.d0
    endif

    if(present(err)) err = err_status

end subroutine interpol3D_arr

subroutine interpol3D_2D_arr(r, z, phi, f, rout, zout, phiout, fout, err, coeffs)
    !+ Performs cylindrical interpolation on a 3D grid of 2D arrays
    !+ f(:,:,r,z,phi)
    real(Float64), dimension(:), intent(in) :: r
        !+ R values
    real(Float64), dimension(:), intent(in) :: z
        !+ Z values
    real(Float64), dimension(:), intent(in) :: phi
        !+ Phi values
    real(Float64), dimension(:,:,:,:,:), intent(in) :: f
        !+ Values at r,z,phi: f(:,:,r,z,phi)
    real(Float64), intent(in)               :: rout
        !+ R value to interpolate
    real(Float64), intent(in)               :: zout
        !+ Z value to interpolate
    real(Float64), intent(in)               :: phiout
        !+ Phi value to interpolate
    real(Float64), dimension(:,:), intent(out)    :: fout
        !+ Interpolant: f(:,:,rout,zout,phiout)
    integer, intent(out), optional          :: err
        !+ Error code
    type(InterpolCoeffs3D), intent(in), optional :: coeffs
        !+ Precomputed Interpolation Coefficients

    type(InterpolCoeffs3D) :: b
    integer :: i, j, k, k2, err_status
    integer :: nphi

    err_status = 1

    nphi = size(phi)

    if(present(coeffs)) then
        b = coeffs
        if(nphi .eq. 1) then
            b%b212 = 0
            b%b222 = 0
            b%b122 = 0
            b%b112 = 0
            b%k = 1
        endif
        err_status = 0
    else
        call interpol_coeff(r,z,phi,rout,zout,phiout,b,err_status)
    endif

    if(err_status.eq.0) then
        i = b%i
        j = b%j
        k = b%k
        if(nphi .eq. 1) then
            k2 = min(k+1,nphi)
        else
            k2 = k+1
        endif
        fout = b%b111*f(:,:,i,j,k)    + b%b121*f(:,:,i,j+1,k) +   &
               b%b112*f(:,:,i,j,k2)   + b%b122*f(:,:,i,j+1,k2) +  &
               b%b211*f(:,:,i+1,j,k)  + b%b221*f(:,:,i+1,j+1,k) + &
               b%b212*f(:,:,i+1,j,k2) + b%b222*f(:,:,i+1,j+1,k2)
    else
        fout = 0.0
    endif

    if(present(err)) err = err_status

end subroutine interpol3D_2D_arr

!============================================================================
!------------------------------ Misc. Routines ------------------------------
!============================================================================
function time_string(time_start) result (time_str)
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

end function time_string

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
