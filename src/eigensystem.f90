!+ This file contains all the linear algebra routines FIDASIM uses
module eigensystem
  !+ A basic libary for calculating matrix eigen-decompositions and inverses
  implicit none
  !!Definition for the kind of the variables:
  integer , parameter   :: long      = kind(int(1))
  integer , parameter   :: float     = kind(1.e0)
  integer , parameter   :: double    = kind(1.d0)
  !! eigenvalue decomposition values
  real(double),parameter:: ONE=1.d0,TWO=2.d0,ZERO=0.d0
  real(double),parameter:: XMACH_EPS=2.22d-16
  integer , parameter   :: MAXIT=50
contains
  ! first subroutines for eigenvalue decomposition
  subroutine RSWAP(a,b)
    !+ Swaps values `a` and `b`
    real(double) :: a,b, t
    t=a; a=b; b=t
  end subroutine RSWAP
  subroutine balance(n,  mat, scal, low, high )
    !+Balances the matrix so that the rows with zero entries
    !+off the diagonal are isolated and the remaining columns and rows
    !+are resized to have one norm close to 1.
    integer, intent(in)  :: n
        !+ Dimension of `mat`
    real(double)         :: mat(0:n-1,0:n-1)
        !+ `n`x`n` scaled matrix
    real(double)         :: scal(0:n-1)
        !+ Contains isolated eigenvalue in the positions 0-`low` and `high`-`n`-1
        !+ its other components contain the scaling factors for transforming `mat`
    integer, intent(out) :: high
    integer, intent(out) :: low


    integer, parameter   :: basis = 2
    real(double)         :: b2, r, c, f, g, s
    integer              :: m, k, i, j, iter
    scal=0.d0
    b2 = basis * basis
    m = 0
    k = n - 1
    iter=1
    do while(iter==1)
       iter = 0
       do j = k, 0, -1
          r = ZERO
          do i = 0, k
             if (i.ne.j)  r = r + DABS(mat(j,i))
          enddo
          if (r == ZERO) then
             scal(k) = j
             if (j.ne.k) then
                do i = 0, k
                   call RSWAP(mat(i,j), mat(i,k))
                enddo
                do i = m, n-1
                   call RSWAP(mat(j,i), mat(k,i))
                enddo
             endif
             k=k-1
             iter = 1
          endif
       enddo !j loop
    enddo !while iter=1
    iter=1
    do while (iter==1)
       iter = 0
       do j = m, k
          c = ZERO
          do i = m, k
             if (i.ne.j)  c = c + DABS(mat(i,j))
          enddo
          if (c == ZERO) then
             scal(m) = j
             if (j.ne.m) then
                do i = 0, k
                   call RSWAP(mat(i,j), mat(i,m))
                enddo
                do i = m, n-1
                   call RSWAP(mat(j,i), mat(m,i))
                enddo
             endif
             m = m + 1
             iter = 1
          endif
       enddo !j loop
    enddo !while iter=1
    low = m
    high = k
    do i = m, k
       scal(i) = ONE
    enddo
    iter=1
    do while (iter==1)
       iter = 0
       do i = m, k
          c=ZERO; r=ZERO
          do j = m, k
             if (j.ne.i) then
                c = c + DABS(mat(j,i))
                r = r + DABS(mat(i,j))
             endif
          enddo
          g = r / basis
          f = ONE
          s = c + r
          do while (c < g)
             f = f * basis
             c = c * b2
          enddo
          g = r * basis
          do while (c >= g)
             f = f / basis
             c = c / b2
          enddo
          if ((c + r) / f < 0.95 * s) then
             g = ONE / f
             scal(i) = scal(i) * f
             iter = 1
             do j = m, n-1
                mat(i,j) = mat(i,j) * g
             enddo
             do j = 0, k
                mat(j,i) = mat(j,i) * f
             enddo
          endif
       enddo !i loop
    enddo !while iter=1
    return
  end subroutine balance
  subroutine balback(n, low, high, scal, eivec )
    !+  Reverses the balancing of balance for the eigenvectors
    integer,intent(in)         :: n
      !+ Dimension of matrix
    integer,intent(in)         :: low
      !+ First nonzero row
    integer,intent(in)         :: high
      !+ Last nonzero row
    real(double), intent(in)   ::  scal(0:n-1)
      !+ Scaling data from balance
    real(double), intent(inout):: eivec(0:n-1,0:n-1)
      !+ Input: n x n matrix of eigenvectors, as computed in qr2
      !+ Output: Non-normalized eigenvectors of the original matrix

    real(double) :: s
    integer      :: i,j,k
    do i = low, high
       s = scal(i)
       do j = 0, n-1
      eivec(i,j) = eivec(i,j) * s
       enddo
    enddo
    do i = low-1, 0, -1
       k = Int(scal(i))
       if (k.ne.i) then
          do j = 0, n-1
             call RSWAP(eivec(i,j), eivec(k,j))
          enddo
       endif
    enddo
    do i = high + 1, n-1
       k = Int(scal(i))
       if (k.ne.i) then
          do j = 0, n-1
             call RSWAP(eivec(i,j), eivec(k,j))
          enddo
       endif
    enddo
    return
  end subroutine balback
  subroutine elmhes(n, low, high, mat, perm )
    !+Transforms the matrix `mat` to upper Hessenberg form.
    integer,intent(in)   :: n
      !+Dimension of `mat`
    integer,intent(in)   :: low
      !+First nonzero row
    integer,intent(in)   :: high
      !+Last nonzero row
    real(double), intent(inout):: mat(0:n-1,0:n-1)
      !+Input: `n`x`n` matrix
      !+Output: Upper Hessenberg matrix; additional information on the tranformation
      !+is stored in the lower triangle
    integer,intent(out)  :: perm(0:n-1)
      !+Permutation vector for elmtrans

    integer              :: i, j, m
    real(double)         ::  x, y
    do m = low + 1, high-1
       i = m
       x = ZERO
       do j = m, high
          if (DABS(mat(j,m-1)) > DABS (x)) then
             x = mat(j,m-1)
             i = j
          endif
       enddo
       perm(m) = i
       if (i.ne.m) then
          do j = m - 1, n-1
             call RSWAP(mat(i,j), mat(m,j))
          enddo
          do j = 0, high
             call RSWAP(mat(j,i), mat(j,m))
          enddo
       endif
       if (x.ne.ZERO) then
          do i = m + 1, high
             y = mat(i,m-1)
             if (y.ne.ZERO) then
                y = y / x
                mat(i,m-1) = y
                do j = m, n-1
                   mat(i,j) = mat(i,j) - y * mat(m,j)
                enddo
                do j = 0, high
                   mat(j,m) = mat(j,m) + y * mat(j,i)
                enddo
             endif
          enddo !i loop
       endif !x <> ZERO
    enddo !m loop
  end subroutine elmhes
  Subroutine elmtrans(n, low, high, mat, perm, h  )
    !+  Elmtrans copies the Hessenberg matrix stored in `mat` to `h`
    integer,intent(in)         :: n
      !+ Dimension of mat
    integer,intent(in)         :: low
      !+ First nonzero row
    integer,intent(in)         :: high
      !+ Last nonzero row
    real(double), intent(in)   :: mat(0:n-1,0:n-1)
      !+ `n`x`n` input matrix
    integer,intent(in)         :: perm(0:n-1)
      !+ Permutation data from elmhes
    real(double),intent(out)   :: h(0:n-1,0:n-1)
      !+ Hessenberg matrix
    integer                    :: i, j, k
    do i = 0, n-1
       do k = 0, n-1
      h(i,k) = ZERO
       enddo
       h(i,i) = ONE
    enddo
    do i = high - 1, low+1, -1
       j = perm(i)
       do k = i + 1, high
      h(k,i) = mat(k,i-1)
       enddo
       if (i.ne.j) then
          do k = i, high
             h(i,k) = h(j,k)
             h(j,k) = ZERO
          enddo
          h(j,i) = ONE
       endif
    enddo
  end subroutine elmtrans
  subroutine Comdiv(ar, ai, br, bi, cr, ci, rc )
    !+ Performs complex division `c` = `a` / `b`
    real(double) ::  ar
      !+ Real part of numerator
    real(double) ::  ai
      !+ Imaginary part of numerator
    real(double) ::  br
      !+ Real part of denominator
    real(double) ::  bi
      !+ Imaginary part of denominator
    real(double) ::  cr
      !+ Real part of quotient
    real(double) ::  ci
      !+ Imaginary part of quotient
    integer      :: rc
      !+ return code

    real(double) :: tmp
    if (br == ZERO.AND.bi == ZERO) then
       rc = 1
       return
    endif
    if (dabs(br) > dabs(bi)) then
       tmp = bi / br
       br  = tmp * bi + br
       cr  = (ar + tmp * ai) / br
       ci  = (ai - tmp * ar) / br
    else
       tmp = br / bi
       bi  = tmp * br + bi
       cr  = (tmp * ar + ai) / bi
       ci  = (tmp * ai - ar) / bi
    endif
    rc = 0
  end subroutine Comdiv !Comdiv
  function comabs(ar,ai)
    !+ Calculates absolute value of a complex number `a`
    real(double) :: ar
      !+ Real part of `a`
    real(double) :: ai
      !+ Imaginary part of `a`
    real(double) :: comabs
      !+ Absolute value of `a`

    if (ar == ZERO.and.ai == ZERO) then
       Comabs = ZERO
       return
    endif
    ar = DABS(ar)
    ai = DABS(ai)
    if (ai > ar) then                                  !Switch  ai and ar
       call RSWAP(ai, ar)
    endif
    if (ai == ZERO) then
       Comabs = ar
    else
       Comabs = ar * DSQRT(ONE + ai / ar * ai / ar)
    endif
  end function comabs
  subroutine  hqrvec(n,     & !Dimension of matrix .......
       low,   & !first nonzero row .........
       high,  & !last nonzero row ..........
       h,     & !upper Hessenberg matrix ...
       wr,    & !Real parts of evalues .....
       wi,    & !Imaginary parts of evalues
       eivec, & !Eigenvectors ..............
       rc  )   !return code ...............
    !+Computes the eigenvectors for the eigenvalues found in hqr2
    !+
    !+###Input parameters
    !+   n     :   int n;  ( n > 0 )
    !+         :   Dimension of  mat and eivec, number of eigenvalues.
    !+
    !+   low   :   int low;
    !+
    !+   high  :   int high; see  balance
    !+
    !+   h     :   n x n upper Hessenberg matrix
    !+
    !+   wr    :   vector of size n;
    !+         :   Real parts of the n eigenvalues.
    !+
    !+   wi    :   vector of size n; Imaginary parts of the n eigenvalues.
    !+
    !+###Output parameter:
    !+   eivec :  n x n matrix, whose columns are the eigenvectors
    integer,intent(in)         :: n
    integer,intent(in)         :: high, low
    real(double), intent(in)   :: wr(0:n-1),wi(0:n-1)
    real(double), intent(out)  :: eivec(0:n-1,0:n-1)
    real(double)  :: h(0:n-1,0:n-1)
    integer :: rc
    integer :: i, j, m, k, na, l
    integer :: code, en
    real(double)  :: p, q, r, s, t, w, x, y, z, ra, sa, vr, vi, norm, temp

    r=ZERO; s=ZERO; z=ZERO; norm=ZERO
    do i = 0, n-1                               !find norm of h
       do j = i, n-1
          norm = norm + DABS(h(i,j))
       enddo
    enddo
    if (norm == ZERO) then
       rc = 1                                    !zero matrix
       return
    endif
    do en = n-1, 0, -1                          !transform back
       p = wr(en)
       q = wi(en)
       na = en - 1
       if (q == ZERO) then
          m = en
          h(en,en) = ONE
          do i = na, 0, -1
             w = h(i,i) - p
             r = h(i,en)
             do j = m, na
                r = r + h(i,j) * h(j,en)
             enddo
             if (wi(i) < ZERO) then
                z = w
                s = r
             else
                m = i
                if (wi(i) == ZERO) then
                   if (w.ne.ZERO) then
                      temp = w
                   else
                      temp=XMACH_EPS * norm
                   endif
                   h(i,en) = -r/temp
                else
                   !Solve the linear system:
                   !| w   x |  | h[i][en]   |   | -r |
                   !|       |  |            | = |    |
                   !| y   z |  | h[i+1][en] |   | -s |
                   x = h(i,i+1)
                   y = h(i+1,i)
                   q = (wr(i) - p)**2 + wi(i)**2
                   h(i,en) = (x * s - z * r) / q
                   t = h(i,en)
                   if (DABS(x) > DABS(z)) then
                      temp = (-r -w * t) / x
                   else
                      temp = (-s -y * t) / z
                   endif
                   h(i+1,en) = temp
                endif
             endif !wi[i] < 0
          enddo !i loop
       else if (q < ZERO) then
          m = na
          if (DABS(h(en,na)) > DABS(h(na,en))) then
             h(na,na) = - (h(en,en) - p) / h(en,na)
             h(na,en) = - q / h(en,na)
          else
             call Comdiv(-h(na,en),0.d0, h(na,na)-p, q, h(na,na), h(na,en),code)
          endif
          h(en,na) = ONE
          h(en,en) = ZERO
          do i = na - 1, 0, -1
             w = h(i,i) - p
             ra = h(i,en)
             sa = ZERO
             do j = m, na
                ra = ra + h(i,j) * h(j,na)
                sa = sa + h(i,j) * h(j,en)
             enddo
             if (wi(i) < ZERO) then
                z = w
                r = ra
                s = sa
             else
                m = i
                if (wi(i) == ZERO) then
                   call Comdiv(-ra, -sa, w, q, h(i,na), h(i,en),code)
                else
            !  solve complex linear system:
                   !| w+i*q     x | | h[i][na] + i*h[i][en]  |   | -ra+i*sa |
                   !|             | |                        | = |          |
            !|   y    z+i*q| | h[i+1][na]+i*h[i+1][en]|   | -r+i*s   |
                   x = h(i,i+1)
                   y = h(i+1,i)
                   vr = (wr(i) - p)**2 + wi(i)**2 - q*q
                   vi = TWO * q * (wr(i) - p)
                   if (vr == ZERO.AND.vi == ZERO) then
                      vr = XMACH_EPS * norm * (DABS(w) + DABS(q)  &
                           + DABS(x) + DABS(y) + DABS(z))
                   endif

                   call Comdiv (x*r-z*ra+q*sa,x*s-z*sa-q*ra &
                        ,vr,vi,h(i,na),h(i,en),code)
                   if (DABS(x) > DABS(z) + DABS(q)) then
                      h(i+1,na) = (-ra - w * h(i,na) + q * h(i,en)) / x
                      h(i+1,en) = (-sa - w * h(i,en) - q * h(i,na)) / x
                   else
                      call Comdiv (-r - y * h(i,na), -s - y * h(i,en) &
                           , z, q, h(i+1,na), h(i+1,en),code)
                   endif
                endif !wi[i] = 0
             endif !wi[i] < 0
          enddo !i loop
       endif !else if q < 0
    enddo !en loop
    do i = 0, n-1                        !Eigenvectors for the evalues for
       if (i < low.or.i > high) then      !rows < low  and rows  > high
          do k = i + 1, n-1
             eivec(i,k) = h(i,k)
          enddo
       endif
    enddo
    j = n-1
    do while (j>=low)
       if(j<=high)then
          m =j
       else
          j = high
       endif
       if(j < 0) exit
       if (wi(j) < ZERO) then
          l=j-1
          do i = low, high
             y=ZERO; z=ZERO
             do k = low, m
                y = y + eivec(i,k) * h(k,l)
                z = z + eivec(i,k) * h(k,j)
             enddo
             eivec(i,l) = y
             eivec(i,j) = z
          enddo
       else
          if (wi(j) == ZERO) then
             do i = low, high
                z = ZERO
                do k = low, m
                   z = z + eivec(i,k) * h(k,j)
                enddo
                eivec(i,j) = z
             enddo
          endif
       endif
       j = j - 1
    enddo !j loop
    rc = 0
  end subroutine hqrvec
  subroutine hqr2(n,    &  !Dimension of matrix .........
       low,  &  !first nonzero row ...........
       high, &  !last nonzero row ............
       h,    &  !Hessenberg matrix ...........
       wr,   &  !Real parts of eigenvalues ...
       wi,   &  !Imaginary parts of evalues ..
       eivec,&  !Matrix of eigenvectors ......
       cnt,  &  !Iteration counter ...........
       rc   )   !return code .................
    !+Computes the eigenvalues and (if vec = True) the eigenvectors
    !+of an  n * n upper Hessenberg matrix.
    !+
    !+###Input parameters
    !+   n     :  integer;  ( n > 0 )
    !+            Dimension of  h and eivec,
    !+            length of the real parts vector  wr and of the
    !+            imaginary parts vector  wi of the eigenvalues.
    !+
    !+   low   :  integer;
    !+
    !+   high  :  integer;  see balance
    !+
    !+   h     :  n x n matrix;
    !+            upper Hessenberg matrix as output of Elmhes
    !+            (destroyed in the process).
    !+###Output parameters
    !+   eivec :  n x n matrix;  (only if vec = 1)
    !+            Matrix, which for vec = 1 contains the
    !+            eigenvectors as follows:
    !+            For real eigebvalues the corresponding column
    !+            contains the corresponding eigenvactor, while for
    !+            complex eigenvalues the corresponding column contains
    !+            the real part of the eigenvactor with its imaginary
    !+            part is stored in the subsequent column of eivec.
    !+            The eigenvactor for the complex conjugate eigenvactor
    !+            is given by the complex conjugate eigenvactor.
    !+
    !+   wr    :  vector of size n;
    !+            Real part of the n eigenvalues.
    !+
    !+   wi    :  vector of size n;
    !+            Imaginary parts of the eigenvalues
    !+
    !+   cnt   :  Integer vector of size n;
    !+            vector of iterations used for each eigenvalue.
    !+            For a complex conjugate eigenvalue pair the second
    !+            entry is negative.
    integer,intent(in)          :: n
    integer,intent(in)          :: high, low
    real(double) ,intent(out)   :: h(0:n-1,0:n-1)
    real(double), intent(out)   :: wr(0:n-1),wi(0:n-1)
    real(double), intent(out)   :: eivec(0:n-1,0:n-1)
    integer,intent(out)         :: rc
    integer,intent(out)         :: cnt(0:n-1)
    integer :: en
    integer :: i, j, na, iter, l, ll, m, k
    real(double)  :: p, q, r, s, t, w, x, y, z
    real(double)  :: r_p,r_r,r_s,r_x,r_z

    p=ZERO; q=ZERO; r=ZERO
    do i = 0, n-1
       if (i < low.or.i > high) then
          wr(i) = h(i,i)
          wi(i) = ZERO
          cnt(i) = 0
       endif
    enddo
    en = high
    t = ZERO
    do while (en >= low)
       iter = 0
       na = en - 1
       do while(1<2)
          ll=999
          do l = en, low+1, -1                      !search for small
             !subdiagonal element
             if(DABS(h(l,l-1))<=XMACH_EPS*(DABS(h(l-1,l-1))+DABS(h(l,l))))then
                ll=l;      !save current index
                goto 10    !exit l loop
             endif
          enddo
10        if(ll.ne.999)then
             l=ll
      else
             l=0          !restore l
          endif
          x = h(en,en)
          if (l == en) then                         !found one evalue
             wr(en) = x + t
             h(en,en) = x + t
             wi(en) = ZERO
             cnt(en) = iter
             en = en - 1
             goto 15      !exit from loop while(True)
          endif
          y = h(na,na)
          w = h(en,na) * h(na,en)
          if (l == na) then                         !found two evalues
             p = (y - x) * 0.5d0
             q = p * p + w
             z = DSQRT(DABS(q))
             x = x + t
             h(en,en) = x + t
             h(na,na) = y + t
             cnt(en) = -iter
             cnt(na) = iter
             if (q >= ZERO) then                     !real eigenvalues
                if (p<ZERO) then
                   z=p-z
                else
                   z=p+z
                endif
                r_z = 1.0d0/z
                wr(na) = x + z
                wr(en) = x - w * r_z
                s = w - w * r_z
                wi(na) = ZERO
                wi(en) = ZERO
                x = h(en,na)
                r = DSQRT (x * x + z * z)
                r_r = 1.0d0/r
                p = x * r_r
                q = z * r_r
                do j = na, n-1
                   z = h(na,j)
                   h(na,j) = q * z + p * h(en,j)
                   h(en,j) = q * h(en,j) - p * z
                enddo
                do i = 0, en
                   z = h(i,na)
                   h(i,na) = q * z + p * h(i,en)
                   h(i,en) = q * h(i,en) - p * z
                enddo
                do i = low, high
                   z = eivec(i,na)
                   eivec(i,na) = q * z + p * eivec(i,en)
                   eivec(i,en) = q * eivec(i,en) - p * z
                enddo
             else                                  !pair of complex
                wr(na) = x + p
                wr(en) = x + p
                wi(na) =   z
                wi(en) = - z
             endif !if q>=ZERO
             en = en - 2
             goto 15                               !exit while(1<2)
          endif !if l = na
          if (iter >= MAXIT) then
             cnt(en) = MAXIT + 1
             rc = en
             write(*,*) ' stop at iter >= MAXIT.'
             return
          endif
          if (iter.ne.0.and.MOD(iter,10) == 0) then
             t = t + x
             do i = low, en
                h(i,i) = h(i,i) - x
             enddo
             s = DABS(h(en,na)) + DABS(h(na,en-2))
             x = 0.75d0 * s; y = x
             w = -0.4375d0 * s * s
          endif
          iter = iter + 1
          do m = en - 2, l, -1
             z = h(m,m)
             r = x - z
             s = y - z
             p = ( r * s - w ) / h(m+1,m) + h(m,m+1)
             q = h(m + 1,m + 1) - z - r - s
             r = h(m + 2,m + 1)
             s = DABS(p) + DABS(q) + DABS (r)
             r_s = 1.0d0/s
             p = p * r_s
             q = q * r_s
             r = r * r_s
             if (m == l)  goto 12
             if (DABS(h(m,m-1)) * (DABS(q) + DABS(r)) <= XMACH_EPS * DABS(p) &
                  * (DABS(h(m-1,m-1)) + DABS(z) + DABS(h(m+1,m+1)))) then
                goto 12                !exit m loop
             endif
          enddo
12        do i = m + 2, en
             h(i,i-2) = ZERO
          enddo
          do i = m + 3, en
             h(i,i-3) = ZERO
          enddo
          do k = m, na
             if(k.ne.m)then!double QR step, for rows l to en and columns m to en
                p = h(k,k-1)
                q = h(k+1,k-1)
                if (k.ne.na) then
                   r = h(k+2,k-1)
                else
                   r = ZERO
                endif
                x = DABS(p) + DABS(q) + DABS(r)
                if (x == ZERO) goto 30                  !next k
                r_x = 1.0d0/x
                p = p * r_x
                q = q * r_x
                r = r * r_x
             endif
             s = DSQRT(p * p + q * q + r * r)
             if (p < ZERO) s = -s
             r_s = 1.0d0/s
             if (k.ne.m) then
                h(k,k-1) = -s * x
             else if (l.ne.m) then
                h(k,k-1) = -h(k,k-1)
             endif
             p = p + s
             r_p = 1.0d0/p
             x = p * r_s
             y = q * r_s
             z = r * r_s
             q = q * r_p
             r = r * r_p
             do j = k, n-1                          !modify rows
                p = h(k,j) + q * h(k+1,j)
                if (k.ne.na) then
                   p = p + r * h(k+2,j)
                   h(k+2,j) = h(k+2,j) - p * z
                endif
                h(k+1,j) = h(k+1,j) - p * y
                h(k,j)   = h(k,j) - p * x
             enddo
             if (k+3 < en) then
                j=k+3
             else
                j=en
             endif
             do i = 0, j                            !modify columns
                p = x * h(i,k) + y * h(i,k+1)
                if (k.ne.na) then
                   p = p + z * h(i,k+2)
                   h(i,k+2) = h(i,k+2) - p * r
                endif
                h(i,k+1) = h(i,k+1) - p * q
                h(i,k)   = h(i,k) - p
             enddo
             do i = low, high
                p = x * eivec(i,k) + y * eivec(i,k+1)
                if (k.ne.na) then
                   p = p + z * eivec(i,k+2)
                   eivec(i,k+2) = eivec(i,k+2) - p * r
                endif
                eivec(i,k+1) = eivec(i,k+1) - p * q
                eivec(i,k)   = eivec(i,k) - p
             enddo
30           continue
          enddo !k loop
       enddo !while(1<2)
15  continue
    enddo !while en >= low                         All evalues found
    !transform evectors back
    call hqrvec (n, low, high, h, wr, wi, eivec,rc)
  end subroutine hqr2

  subroutine eigen (n, matrix, eigvec, eigval)
    !+The subroutine eigen  determines all eigenvalues and (if desired)
    !+all eigenvectors of a real square  n * n  matrix via the QR method
    !+in the version of Martin, Parlett, Peters, Reinsch and Wilkinson.
    !+
    !+###Literature
    !+1. Peters, Wilkinson: Eigenvectors of real and complex
    !+   matrices by LR and QR triangularisations,
    !+   Num. Math. 16, p.184-204, (1970); [PETE70]; contribution
    !+   II/15, p. 372 - 395 in [WILK71].
    !+2. Martin, Wilkinson: Similarity reductions of a general
    !+   matrix to Hessenberg form, Num. Math. 12, p. 349-368,(1968)
    !+   [MART 68]; contribution II,13, p. 339 - 358 in [WILK71].
    !+3. Parlett, Reinsch: Balancing a matrix for calculations of
    !+   eigenvalues and eigenvectors, Num. Math. 13, p. 293-304,
    !+   (1969); [PARL69]; contribution II/11, p.315 - 326 in
    !+   [WILK71].
    !+
    !+###Input parameters
    !+   n     :  integer; ( n > 0 )
    !+            size of matrix, number of eigenvalues
    !+
    !+   mat   :  n x n matrix;
    !+            input matrix
    !+
    !+###Output parameters
    !+   eivec :  n x n matrix;     (only if vec = 1)
    !+            matrix, if  vec = 1  that holds the eigenvectors
    !+            thus :
    !+            If the jth eigenvalue of the matrix is real then the
    !+            jth column is the corresponding real eigenvector;
    !+            if the jth eigenvalue is complex then the jth column
    !+            of eivec contains the real part of the eigenvector
    !+            while its imaginary part is in column j+1.
    !+            (the j+1st eigenvector is the complex conjugate
    !+            vector.)
    !+
    !+   valre :  vector of size n;
    !+            Real parts of the eigenvalues.
    !+
    !+   valim :  vector of size n;
    !+            Imaginary parts of the eigenvalues
    !+
    !+   cnt   :  Integer vector of size n;
    !+            vector containing the number of iterations for each
    !+            eigenvalue. (for a complex conjugate pair the second
    !+            entry is negative).
    integer      ,intent(in)                 :: n ! nlevels
    real(double) ,intent(in),dimension(n,n)  :: matrix
    real(double) ,intent(out),dimension(n,n) :: eigvec
    real(double) ,intent(out),dimension(n)   :: eigval
    real(double)     :: mat(0:n-1,0:n-1)
    real(double)     :: eivec(0:n-1,0:n-1)
    real(double)     :: valre(0:n-1) !real parts of eigenvalues
    real(double)     :: valim(0:n-1) !imaginary parts of eigenvalues
    integer          :: rc             !return code
    integer          :: cnt(0:n-1)   !Iteration counter
    integer          :: high, low
    real(double)     :: d(0:n-1), scale(0:n-1)
    integer          :: perm(0:n-1)

    cnt=0 ; d=0.d0
    mat(0:n-1,0:n-1)=matrix(1:n,1:n)
    !balance mat for nearly
    call balance(n, mat, scale, low, high)      !equal row and column
    !reduce mat to upper
    call elmhes(n, low, high, mat, perm)        !reduce mat to upper
    !Hessenberg form
    call elmtrans(n, low, high, mat, perm, eivec)
    !QR algorithm for eigenvalues and eigenvectors
    call hqr2(n, low, high, mat, valre, valim, eivec, cnt,rc)
    !reverse balancing to determine eigenvectors
    call balback(n, low, high, scale, eivec)
    if (rc.ne.0) then
      print*, 'matrix = '
      print*, matrix
      stop 'problem in eigen!'
    endif
    eigval(1:n)=valre(0:n-1)
    eigvec(1:n,1:n)=eivec(0:n-1,0:n-1)
  end subroutine eigen

  function outerprod(a,b)
    !+ Calculates outer product
    real(double), dimension(:), intent(IN)   :: a,b
    real(double), dimension(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  end function outerprod
  subroutine swap(a,b)
    !+Swap arrays `a` and `b`
    real(double), dimension(:), intent(INOUT) :: a,b
    real(double), dimension(size(a))          :: dum
    dum=a
    a=b
    b=dum
  end subroutine swap
  subroutine ludcmp(a,indx,d)
    !+Calculates LU decomposition
    real(double), dimension(:,:),intent(INOUT):: a
    integer,dimension(:),  intent(OUT)        :: indx
    real(double),                intent(OUT)  :: d
    real(double), dimension(size(a,1))        :: vv
    integer,dimension(1)                      :: imaxloc
    integer :: j,n,imax
    n=size(indx)
    d=1.0
    vv=maxval(abs(a),dim=2)
    if(any(vv.eq.0.))stop 'singular matrix in ludcmp'
    vv=1.d0/vv
    do j=1,n
       imaxloc=maxloc(vv(j:n)*abs(a(j:n,j)))
       imax=(j-1)+imaxloc(1)
       if (j /= imax) then
          call swap(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if (a(j,j) == 0.0) a(j,j)=1.0d-20
       a(j+1:n,j)=a(j+1:n,j)/a(j,j)
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    enddo
  end subroutine ludcmp
  subroutine lubksb(a,indx,b)
    !+ Does LU back substitution
    real(double), dimension(:,:),intent(IN)   :: a
    integer,dimension(:),  intent(IN)         :: indx
    real(double), dimension(:),  intent(INOUT):: b
    integer       :: i,n,ii,ll
    real(double)  :: summ
    n=size(indx)
    ii=0
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i
       endif
       b(i)=summ
    enddo
    do i=n,1,-1
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    enddo
  end subroutine lubksb
  subroutine matinv(a, b)
    !+ Matrix inversion with LU-decomposition
    !====================================================
    real(double), dimension(:,:), intent(IN)             :: a
    real(double), dimension(:,:), intent(OUT)            :: b
    real(double), dimension(size(a,dim=1),size(a,dim=2)) :: ah, y
    integer                                              :: i, N
    integer, dimension(size(a,dim=1))                    :: indx
    real(double)                                         :: d
    N = size(a,dim=1)
    if (N /= size(a,dim=2)) stop 'SUB matinv: ludcmp matrix must be square!'
    ah = a
    y  = 0.
    do i = 1, N
       y(i,i) = 1.d0
    enddo
    call ludcmp(ah,indx,d)
    do i = 1, N
       call lubksb(ah,indx,y(:,i))
    enddo
    b = y
  end subroutine matinv


  subroutine linsolve(a,b,x)
   !+ Solve linear equations A * X = B
   real(double), dimension(:,:), intent(IN)  :: a ! assume a square
   real(double), dimension(:), intent(IN)    :: b ! assume size(b) == size(a,1)
   real(double), dimension(:), intent(OUT)   :: x ! assume size(x) == size(b)

#ifdef _USE_BLAS
   real(double), dimension(size(b),size(b)) :: lu
   integer, dimension(size(b)) :: ipiv
   integer :: n,info

   n = size(b)

   ! first factorize a
   lu(:,:) = a(:,:)
   call DGETRF(n,n,lu,n,ipiv,info)
   if (info /= 0) stop 'sub linsolve: DGETRF failed!'

   x(:) = b(:)
   call DGETRS('N',n,1,lu,n,ipiv,x,n,info)
   if (info /= 0) stop 'sub linsolve: DGETRS failed!'

#else
   real(double), dimension(size(b),size(b)) :: a_inv
   call matinv(a, a_inv)
   x = matmul(a_inv, b)!coeffs determined from states at t=0
#endif

  end subroutine linsolve

end module eigensystem
