!+ Specialized 6x6 matrix routines optimized for performance
!+ These routines replace general-purpose linear algebra operations
!+ for the fixed-size 6x6 matrices used in colrad calculations
module matrix6x6
    use iso_fortran_env, only: real64
    implicit none

    integer, parameter :: Float64 = real64
    integer, parameter :: N6 = 6
    real(Float64), parameter :: ZERO = 0.0d0
    real(Float64), parameter :: ONE = 1.0d0
    real(Float64), parameter :: TWO = 2.0d0
    real(Float64), parameter :: HALF = 0.5d0
    real(Float64), parameter :: TOL = 1.0d-14
    real(Float64), parameter :: SMALL = 1.0d-20

    private
    public :: eigen6x6, linsolve6x6, matmul6x6_vec, matmul6x6
    public :: matrix_is_symmetric

contains

    !========================================================================
    ! Helper routines
    !========================================================================

    pure function matrix_is_symmetric(A) result(is_sym)
        !+ Check if a 6x6 matrix is symmetric
        real(Float64), dimension(6,6), intent(in) :: A
        logical :: is_sym
        integer :: i, j
        real(Float64) :: max_diff

        max_diff = ZERO
        do i = 1, 6
            do j = i+1, 6
                max_diff = max(max_diff, abs(A(i,j) - A(j,i)))
            enddo
        enddo

        is_sym = (max_diff < TOL)
    end function matrix_is_symmetric

    pure subroutine swap_rows6(A, b, i, j)
        !+ Swap rows i and j in matrix A and vector b
        real(Float64), dimension(6,6), intent(inout) :: A
        real(Float64), dimension(6), intent(inout) :: b
        integer, intent(in) :: i, j
        real(Float64) :: temp
        integer :: k

        if (i == j) return

        do k = 1, 6
            temp = A(i,k)
            A(i,k) = A(j,k)
            A(j,k) = temp
        enddo

        temp = b(i)
        b(i) = b(j)
        b(j) = temp
    end subroutine swap_rows6

    pure subroutine jacobi_rotation(A, V, p, q, c, s)
        !+ Apply Jacobi rotation to matrix A and accumulate in V
        real(Float64), dimension(6,6), intent(inout) :: A, V
        integer, intent(in) :: p, q
        real(Float64), intent(in) :: c, s
        integer :: i
        real(Float64) :: temp1, temp2

        ! Update rows p and q of A
        do i = 1, 6
            if (i /= p .and. i /= q) then
                temp1 = c * A(p,i) - s * A(q,i)
                temp2 = s * A(p,i) + c * A(q,i)
                A(p,i) = temp1
                A(q,i) = temp2
            endif
        enddo

        ! Update columns p and q of A
        do i = 1, 6
            if (i /= p .and. i /= q) then
                temp1 = c * A(i,p) - s * A(i,q)
                temp2 = s * A(i,p) + c * A(i,q)
                A(i,p) = temp1
                A(i,q) = temp2
            endif
        enddo

        ! Update diagonal elements
        temp1 = c * c * A(p,p) + s * s * A(q,q) - TWO * s * c * A(p,q)
        temp2 = s * s * A(p,p) + c * c * A(q,q) + TWO * s * c * A(p,q)
        A(p,p) = temp1
        A(q,q) = temp2
        A(p,q) = ZERO
        A(q,p) = ZERO

        ! Accumulate eigenvectors
        do i = 1, 6
            temp1 = c * V(i,p) - s * V(i,q)
            temp2 = s * V(i,p) + c * V(i,q)
            V(i,p) = temp1
            V(i,q) = temp2
        enddo
    end subroutine jacobi_rotation

    pure subroutine max_offdiag(A, p, q, max_val)
        !+ Find maximum off-diagonal element and its indices
        real(Float64), dimension(6,6), intent(in) :: A
        integer, intent(out) :: p, q
        real(Float64), intent(out) :: max_val
        integer :: i, j

        max_val = ZERO
        p = 1
        q = 2

        do i = 1, 5
            do j = i+1, 6
                if (abs(A(i,j)) > max_val) then
                    max_val = abs(A(i,j))
                    p = i
                    q = j
                endif
            enddo
        enddo
    end subroutine max_offdiag

    !========================================================================
    ! Eigenvalue decomposition for symmetric matrices (Jacobi method)
    !========================================================================

    pure subroutine eigen6x6_symmetric(A, eigvec, eigval)
        !+ Specialized eigenvalue decomposition for 6x6 symmetric matrices
        !+ Uses Jacobi method optimized for small matrices
        real(Float64), dimension(6,6), intent(in) :: A
        real(Float64), dimension(6,6), intent(out) :: eigvec
        real(Float64), dimension(6), intent(out) :: eigval

        real(Float64), dimension(6,6) :: work
        real(Float64) :: c, s, tau, t, theta, max_val
        integer :: p, q, i, iter
        integer, parameter :: max_iter = 50

        ! Initialize
        work = A
        eigvec = ZERO
        do i = 1, 6
            eigvec(i,i) = ONE
        enddo

        ! Jacobi iterations
        do iter = 1, max_iter
            ! Find largest off-diagonal element
            call max_offdiag(work, p, q, max_val)

            ! Check convergence
            if (max_val < TOL) exit

            ! Calculate rotation angle
            if (abs(work(p,q)) > SMALL) then
                theta = (work(q,q) - work(p,p)) / (TWO * work(p,q))
                t = sign(ONE, theta) / (abs(theta) + sqrt(ONE + theta**2))
                c = ONE / sqrt(ONE + t**2)
                s = t * c

                ! Apply rotation
                call jacobi_rotation(work, eigvec, p, q, c, s)
            endif
        enddo

        ! Extract eigenvalues from diagonal
        do i = 1, 6
            eigval(i) = work(i,i)
        enddo
    end subroutine eigen6x6_symmetric

    !========================================================================
    ! Eigenvalue decomposition for general matrices (uses existing eigen)
    !========================================================================

    subroutine eigen6x6_general(A, eigvec, eigval)
        !+ Eigenvalue decomposition for general 6x6 matrices
        !+ Falls back to general eigen routine from eigensystem module
        use eigensystem, only: eigen
        real(Float64), dimension(6,6), intent(in) :: A
        real(Float64), dimension(6,6), intent(out) :: eigvec
        real(Float64), dimension(6), intent(out) :: eigval

        ! Use general eigen routine for non-symmetric matrices
        call eigen(6, A, eigvec, eigval)
    end subroutine eigen6x6_general

    !========================================================================
    ! Main eigenvalue decomposition interface
    !========================================================================

    subroutine eigen6x6(A, eigvec, eigval)
        !+ Eigenvalue decomposition for 6x6 matrices
        !+ Automatically selects symmetric or general algorithm
        real(Float64), dimension(6,6), intent(in) :: A
        real(Float64), dimension(6,6), intent(out) :: eigvec
        real(Float64), dimension(6), intent(out) :: eigval

        ! Check if matrix is symmetric
        if (matrix_is_symmetric(A)) then
            ! Use optimized symmetric algorithm
            call eigen6x6_symmetric(A, eigvec, eigval)
        else
            ! Use general algorithm
            call eigen6x6_general(A, eigvec, eigval)
        endif
    end subroutine eigen6x6

    !========================================================================
    ! Linear solver using Gaussian elimination with partial pivoting
    !========================================================================

    pure subroutine linsolve6x6(A, b, x)
        !+ Specialized linear solver for 6x6 system Ax = b
        !+ Uses Gaussian elimination with partial pivoting
        real(Float64), dimension(6,6), intent(in) :: A
        real(Float64), dimension(6), intent(in) :: b
        real(Float64), dimension(6), intent(out) :: x

        real(Float64), dimension(6,6) :: LU
        real(Float64), dimension(6) :: work
        real(Float64) :: factor, pivot
        integer :: i, j, k, imax

        ! Copy to preserve input
        LU = A
        work = b

        ! Forward elimination with partial pivoting
        ! Column 1
        imax = 1
        do i = 2, 6
            if (abs(LU(i,1)) > abs(LU(imax,1))) imax = i
        enddo
        call swap_rows6(LU, work, 1, imax)

        if (abs(LU(1,1)) < SMALL) then
            x = ZERO
            return
        endif

        factor = ONE / LU(1,1)
        do i = 2, 6
            LU(i,1) = LU(i,1) * factor
            do j = 2, 6
                LU(i,j) = LU(i,j) - LU(i,1) * LU(1,j)
            enddo
            work(i) = work(i) - LU(i,1) * work(1)
        enddo

        ! Column 2
        imax = 2
        do i = 3, 6
            if (abs(LU(i,2)) > abs(LU(imax,2))) imax = i
        enddo
        call swap_rows6(LU, work, 2, imax)

        if (abs(LU(2,2)) < SMALL) then
            x = ZERO
            return
        endif

        factor = ONE / LU(2,2)
        do i = 3, 6
            LU(i,2) = LU(i,2) * factor
            do j = 3, 6
                LU(i,j) = LU(i,j) - LU(i,2) * LU(2,j)
            enddo
            work(i) = work(i) - LU(i,2) * work(2)
        enddo

        ! Column 3
        imax = 3
        do i = 4, 6
            if (abs(LU(i,3)) > abs(LU(imax,3))) imax = i
        enddo
        call swap_rows6(LU, work, 3, imax)

        if (abs(LU(3,3)) < SMALL) then
            x = ZERO
            return
        endif

        factor = ONE / LU(3,3)
        do i = 4, 6
            LU(i,3) = LU(i,3) * factor
            do j = 4, 6
                LU(i,j) = LU(i,j) - LU(i,3) * LU(3,j)
            enddo
            work(i) = work(i) - LU(i,3) * work(3)
        enddo

        ! Column 4
        imax = 4
        do i = 5, 6
            if (abs(LU(i,4)) > abs(LU(imax,4))) imax = i
        enddo
        call swap_rows6(LU, work, 4, imax)

        if (abs(LU(4,4)) < SMALL) then
            x = ZERO
            return
        endif

        factor = ONE / LU(4,4)
        do i = 5, 6
            LU(i,4) = LU(i,4) * factor
            do j = 5, 6
                LU(i,j) = LU(i,j) - LU(i,4) * LU(4,j)
            enddo
            work(i) = work(i) - LU(i,4) * work(4)
        enddo

        ! Column 5
        imax = 5
        if (abs(LU(6,5)) > abs(LU(5,5))) imax = 6
        call swap_rows6(LU, work, 5, imax)

        if (abs(LU(5,5)) < SMALL) then
            x = ZERO
            return
        endif

        factor = ONE / LU(5,5)
        LU(6,5) = LU(6,5) * factor
        LU(6,6) = LU(6,6) - LU(6,5) * LU(5,6)
        work(6) = work(6) - LU(6,5) * work(5)

        ! Back substitution (fully unrolled)
        if (abs(LU(6,6)) < SMALL) then
            x = ZERO
            return
        endif

        x(6) = work(6) / LU(6,6)
        x(5) = (work(5) - LU(5,6)*x(6)) / LU(5,5)
        x(4) = (work(4) - LU(4,5)*x(5) - LU(4,6)*x(6)) / LU(4,4)
        x(3) = (work(3) - LU(3,4)*x(4) - LU(3,5)*x(5) - LU(3,6)*x(6)) / LU(3,3)
        x(2) = (work(2) - LU(2,3)*x(3) - LU(2,4)*x(4) - LU(2,5)*x(5) - LU(2,6)*x(6)) / LU(2,2)
        x(1) = (work(1) - LU(1,2)*x(2) - LU(1,3)*x(3) - LU(1,4)*x(4) - LU(1,5)*x(5) - LU(1,6)*x(6)) / LU(1,1)
    end subroutine linsolve6x6

    !========================================================================
    ! Optimized matrix-vector multiplication
    !========================================================================

    pure function matmul6x6_vec(A, x) result(y)
        !+ Optimized 6x6 matrix-vector multiplication
        !+ Fully unrolled for maximum performance
        real(Float64), dimension(6,6), intent(in) :: A
        real(Float64), dimension(6), intent(in) :: x
        real(Float64), dimension(6) :: y

        ! Fully unrolled for maximum performance
        y(1) = A(1,1)*x(1) + A(1,2)*x(2) + A(1,3)*x(3) + A(1,4)*x(4) + A(1,5)*x(5) + A(1,6)*x(6)
        y(2) = A(2,1)*x(1) + A(2,2)*x(2) + A(2,3)*x(3) + A(2,4)*x(4) + A(2,5)*x(5) + A(2,6)*x(6)
        y(3) = A(3,1)*x(1) + A(3,2)*x(2) + A(3,3)*x(3) + A(3,4)*x(4) + A(3,5)*x(5) + A(3,6)*x(6)
        y(4) = A(4,1)*x(1) + A(4,2)*x(2) + A(4,3)*x(3) + A(4,4)*x(4) + A(4,5)*x(5) + A(4,6)*x(6)
        y(5) = A(5,1)*x(1) + A(5,2)*x(2) + A(5,3)*x(3) + A(5,4)*x(4) + A(5,5)*x(5) + A(5,6)*x(6)
        y(6) = A(6,1)*x(1) + A(6,2)*x(2) + A(6,3)*x(3) + A(6,4)*x(4) + A(6,5)*x(5) + A(6,6)*x(6)
    end function matmul6x6_vec

    !========================================================================
    ! Optimized matrix-matrix multiplication
    !========================================================================

    pure function matmul6x6(A, B) result(C)
        !+ Optimized 6x6 matrix-matrix multiplication
        real(Float64), dimension(6,6), intent(in) :: A, B
        real(Float64), dimension(6,6) :: C
        integer :: i, j

        ! Unrolled inner loops for better performance
        do i = 1, 6
            do j = 1, 6
                C(i,j) = A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j) + &
                         A(i,4)*B(4,j) + A(i,5)*B(5,j) + A(i,6)*B(6,j)
            enddo
        enddo
    end function matmul6x6

end module matrix6x6
