!  ====================================================================
!
!     Library  WuL :  Rewrite NAG library fucnctions
!     Module   lnag_local  :      Subroutines for interactive programs
!
!  ====================================================================

!     1.0.  General Functions
!              FuVal, FuText, FuTxt, FuInv, FuHelp, FuAsk
!
!     2.0.  Functions for solid state physics
!              u2Debye, DebyeInt

!  ====================================================================
!     1.0.  General Functions
!  ====================================================================

        subroutine M01DAF_local(RV, M1, M2, ORDER, IRANK, IFAIL)
!-------------------------------------------------------------
! Computes UNIQUE stable ranks of a real array section RV(M1:M2)
!
! ORDER: 'A' = ascending, 'D' = descending
! IRANK(i) contains a unique rank in 1..(M2-M1+1), i=M1..M2
!
! Ranking rule (stable, unique):
!   ASC:  rank(i) = 1 + count{ j : RV(j) < RV(i)
!                               OR (RV(j)==RV(i) AND j<i) }
!   DESC: rank(i) = 1 + count{ j : RV(j) > RV(i)
!                               OR (RV(j)==RV(i) AND j<i) }
!
! This guarantees all ranks are unique and preserves original
! order among equal values.
!-------------------------------------------------------------
        implicit none
        integer, intent(in)              :: M1, M2
        real(8), intent(in)              :: RV(M1:M2)
        character(len=1), intent(in)     :: ORDER
        integer, intent(out)             :: IRANK(M1:M2)
        integer, intent(out)             :: IFAIL
        integer                          :: i, j

        ! Validate ORDER
        if (ORDER /= 'A' .and. ORDER /= 'D' .and.
     *      ORDER /= 'a' .and. ORDER /= 'd') then
            print *, 'Invalid ORDER = ', ORDER
            IFAIL = 1
            return
        end if
        IFAIL = 0

! Compute unique stable ranks
        do i = M1, M2
            IRANK(i) = 1
            do j = M1, M2
            if (ORDER == 'A' .or. ORDER == 'a') then
                if ( (RV(j) < RV(i)) .or.
     *               ((RV(j) == RV(i)) .and. (j < i)) ) then
                IRANK(i) = IRANK(i) + 1
                end if
            else
                if ( (RV(j) > RV(i)) .or.
     *               ((RV(j) == RV(i)) .and. (j < i)) ) then
                IRANK(i) = IRANK(i) + 1
                end if
            end if
            end do
        end do
        end subroutine M01DAF_local

        subroutine M01DEF_local(RM, IDM, M1, M2, N1, N2, ORDER,
     *                        IRANK, IFAIL)
          !--------------------------------------------------------------------
          ! Stable ranking of rows I=M1..M2 using lexicographic keys from columns
          ! J=N1..N2 (primary key N1, then N1+1, ..., up to N2).
          !
          ! ORDER:
          !   'A' = ascending  (smaller values first)
          !   'D' = descending (larger values first)
          !
          ! Stability:
          !   If two rows are equal in all key columns, their relative order in the
          !   input (increasing I) is preserved in the output ranking.
          !
          ! Output:
          !   IRANK(i, N1) = unique rank (1..NROW) of row i among rows M1..M2.
          !                 (Other IRANK entries are not modified.)
          !
          ! Notes:
          !   - Ranks are unique even when rows tie (stable tie-breaking by input order).
          !   - Only RM(:,N1..N2) is used as the key.
          !--------------------------------------------------------------------
          implicit none

          integer,          intent(in)  :: IDM, M1, M2, N1, N2
          real(8),          intent(in)  :: RM(1:IDM, 1:N2)
          character(len=1), intent(in)  :: ORDER
          integer,          intent(out) :: IRANK(1:IDM, 1:N2)
          integer,          intent(out) :: IFAIL

          integer :: i, k, pos, nrow, keyrow
          integer, allocatable :: idx(:)

          ! Validate ORDER
          if (ORDER /= 'A' .and. ORDER /= 'D') then
            IFAIL = 1
            return
          end if
          IFAIL = 0

          nrow = M2 - M1 + 1
          if (nrow <= 0) then
            ! Empty row range: nothing to do, but not an error
            return
          end if

          if (N2 < N1) then
            ! Empty key range: rank by input order (stable, unique)
            do i = M1, M2
              IRANK(i, N1) = i - M1 + 1
            end do
            return
          end if

          allocate(idx(nrow))

          ! Initialize ordered list of row indices with the first row
          idx(1) = M1
          pos = 1

          ! Insertion-sort remaining rows into idx(1:pos) using stable lexicographic compare
          do i = M1 + 1, M2
            pos = pos + 1
            idx(pos) = i

            keyrow = i
            k = pos

            if (ORDER == 'A') then
              ! Ascending, STABLE: shift while keyrow is STRICTLY less than previous row
              do while (k > 1)
                if (row_less_lex(RM, IDM, keyrow, idx(k-1), N1, N2))
     *           then
                  idx(k) = idx(k-1)
                  k = k - 1
                else
                  exit
                end if
              end do
            else
              ! Descending, STABLE: shift while keyrow is STRICTLY greater than previous row
              do while (k > 1)
                if (row_greater_lex(RM, IDM, keyrow, idx(k-1), N1, N2))
     *           then
                  idx(k) = idx(k-1)
                  k = k - 1
                else
                  exit
                end if
              end do
            end if

            idx(k) = keyrow
          end do

          ! Convert sorted row index list -> unique ranks
          do k = 1, nrow
            IRANK(idx(k), N1) = k
          end do

          deallocate(idx)

        contains

          logical function row_less_lex(RMloc, IDMloc, ia, ib, c1, c2)
            ! TRUE if row ia is lexicographically STRICTLY less than row ib
            implicit none
            integer, intent(in) :: IDMloc, ia, ib, c1, c2
            real(8), intent(in) :: RMloc(1:IDMloc, 1:c2)
            integer :: j
            do j = c1, c2
              if (RMloc(ia,j) < RMloc(ib,j)) then
                row_less_lex = .true.
                return
              else if (RMloc(ia,j) > RMloc(ib,j)) then
                row_less_lex = .false.
                return
              end if
            end do
            ! All key columns equal => NOT strictly less (ensures stability)
            row_less_lex = .false.
          end function row_less_lex

          logical function row_greater_lex(RMloc, IDMloc, ia, ib,
     *                                     c1, c2)
            ! TRUE if row ia is lexicographically STRICTLY greater than row ib
            implicit none
            integer, intent(in) :: IDMloc, ia, ib, c1, c2
            real(8), intent(in) :: RMloc(1:IDMloc, 1:c2)
            integer :: j
            do j = c1, c2
              if (RMloc(ia,j) > RMloc(ib,j)) then
                row_greater_lex = .true.
                return
              else if (RMloc(ia,j) < RMloc(ib,j)) then
                row_greater_lex = .false.
                return
              end if
            end do
            ! All key columns equal => NOT strictly greater (ensures stability)
            row_greater_lex = .false.
          end function row_greater_lex

        end subroutine M01DEF_local

        subroutine M01EAF_local(RV, M1, M2, IRANK, IFAIL)
          !-------------------------------------------------------------
          ! Reorder RV(M1:M2) in-place according to permutation IRANK(M1:M2)
          ! using the "scatter" convention:
          !
          !   RV_after(IRANK(i)) = RV_before(i)
          !
          ! i.e. IRANK(i) is the DESTINATION position for the element
          ! originally at position i.
          !
          ! Requirements:
          !   IRANK(M1:M2) must be a permutation of [M1..M2].
          !
          ! On error:
          !   IFAIL = 1 if IRANK is not a valid permutation.
          !-------------------------------------------------------------
          implicit none

          integer, intent(in)    :: M1, M2
          real(8), intent(inout) :: RV(M1:M2)
          integer, intent(in)    :: IRANK(M1:M2)
          integer, intent(out)   :: IFAIL

          integer :: i, n, j
          real(8), allocatable :: temp(:)
          logical, allocatable :: seen(:)

          IFAIL = 0
          n = M2 - M1 + 1
          if (n <= 0) return

          ! Validate permutation
          allocate(seen(n))
          seen = .false.

          do i = M1, M2
            if (IRANK(i) < M1 .or. IRANK(i) > M2) then
              IFAIL = 1
              deallocate(seen)
              return
            end if
            j = IRANK(i) - M1 + 1
            if (seen(j)) then
              IFAIL = 1
              deallocate(seen)
              return
            end if
            seen(j) = .true.
          end do
          deallocate(seen)

          ! Apply permutation (scatter)
          allocate(temp(n))
          do i = M1, M2
            temp(IRANK(i) - M1 + 1) = RV(i)
          end do
          do i = M1, M2
            RV(i) = temp(i - M1 + 1)
          end do
          deallocate(temp)
        end subroutine M01EAF_local

        subroutine M01EBF_local(IV, M1, M2, IRANK, IFAIL)
          !-------------------------------------------------------------
          ! Reorder IV(M1:M2) in-place according to permutation IRANK(M1:M2)
          ! using the "scatter" convention:
          !
          !   IV_after(IRANK(i)) = RV_before(i)
          !
          ! i.e. IRANK(i) is the DESTINATION position for the element
          ! originally at position i.
          !
          ! Requirements:
          !   IRANK(M1:M2) must be a permutation of [M1..M2].
          !
          ! On error:
          !   IFAIL = 1 if IRANK is not a valid permutation.
          !-------------------------------------------------------------
          implicit none

          integer, intent(in)    :: M1, M2
          integer, intent(inout) :: IV(M1:M2)
          integer, intent(in)    :: IRANK(M1:M2)
          integer, intent(out)   :: IFAIL

          integer :: i, n, j
          integer, allocatable :: temp(:)
          logical, allocatable :: seen(:)

          IFAIL = 0
          n = M2 - M1 + 1
          if (n <= 0) return

          ! Validate permutation
          allocate(seen(n))
          seen = .false.

          do i = M1, M2
            if (IRANK(i) < M1 .or. IRANK(i) > M2) then
              IFAIL = 1
              deallocate(seen)
              return
            end if
            j = IRANK(i) - M1 + 1
            if (seen(j)) then
              IFAIL = 1
              deallocate(seen)
              return
            end if
            seen(j) = .true.
          end do
          deallocate(seen)

          ! Apply permutation (scatter)
          allocate(temp(n))
          do i = M1, M2
            temp(IRANK(i) - M1 + 1) = IV(i)
          end do
          do i = M1, M2
            IV(i) = temp(i - M1 + 1)
          end do
          deallocate(temp)
        end subroutine M01EBF_local


      subroutine M01ZAF_local(IPERM, M1, M2, IFAIL)
        !------------------------------------------------------------
        ! Invert a permutation in place over indices M1..M2.
        !
        ! Input:
        !   IPERM(M1:M2) is a permutation of [M1..M2]
        !
        ! Output:
        !   IPERM becomes its inverse over M1..M2, i.e.
        !     IPERM(old_IPERM(i)) = i
        !
        ! IFAIL:
        !   0 = success
        !   1 = IPERM is not a valid permutation of [M1..M2]
        !------------------------------------------------------------
        implicit none

        integer, intent(in)    :: M1, M2
        integer, intent(inout) :: IPERM(M1:M2)
        integer, intent(out)   :: IFAIL

        integer :: i
        integer, allocatable :: temp(:)
        logical, allocatable :: seen(:)

        IFAIL = 0
        if (M2 < M1) return

        !-----------------------------
        ! Validate permutation property
        !-----------------------------
        allocate(seen(M1:M2))
        seen = .false.

        do i = M1, M2
          if (IPERM(i) < M1 .or. IPERM(i) > M2) then
            IFAIL = 1
            deallocate(seen)
            return
          end if
          if (seen(IPERM(i))) then
            IFAIL = 1
            deallocate(seen)
            return
          end if
          seen(IPERM(i)) = .true.
        end do

        deallocate(seen)

        !-----------------------------
        ! Compute inverse using a temp
        !-----------------------------
        allocate(temp(M1:M2))
        temp = IPERM

        do i = M1, M2
          IPERM(temp(i)) = i
        end do

        deallocate(temp)
      end subroutine M01ZAF_local


      subroutine D01GAF_local(X, Y, N, ANS, ER, IFAIL)
        !--------------------------------------------------------
        ! Integrate tabulated data (X(i), Y(i)) for i=1..N using
        ! the Gill & Miller (1972) algorithm for unequally spaced
        ! data based on local cubic interpolation.
        !
        ! Input:
        !   X(1:N) : strictly monotone abscissae (increasing or decreasing),
        !            with all X(i) distinct
        !   Y(1:N) : ordinates corresponding to X(1:N)
        !   N      : number of points (must be >= 4)
        !
        ! Output:
        !   ANS    : integral estimate over [X(1), X(N)]
        !            (NAG-style: ANS = I + ER, where I is the basic GM estimate)
        !   ER     : integrated truncation-error estimate (indicator, not a bound)
        !
        ! IFAIL:
        !   0 = success
        !   1 = N < 4
        !   2 = X not strictly monotone
        !   3 = repeated X values encountered
        !--------------------------------------------------------
        implicit none
        integer, intent(in)  :: N
        real(8), intent(in)  :: X(1:N), Y(1:N)
        real(8), intent(out) :: ANS, ER
        integer, intent(out) :: IFAIL

        integer :: i
        real(8) :: e, h1, h2, h3, h4, r1, r2, r3, r4, d1, d2, d3,
     *             c, s, intP
        logical :: inc, dec

        ANS   = 0.0d0
        ER    = 0.0d0
        IFAIL = 0

        if (N < 4) then
          IFAIL = 1
          return
        end if

        ! Monotonicity: strictly increasing OR strictly decreasing
        inc = .true.
        dec = .true.
        do i = 1, N-1
          if (X(i+1) <= X(i)) inc = .false.
          if (X(i+1) >= X(i)) dec = .false.
          if (X(i+1) == X(i)) then
            IFAIL = 3   ! distinct points required
            return
          end if
        end do
        if (.not.(inc .or. dec)) then
          IFAIL = 2
          return
        end if

        ! Integrate over full range: ia=1, ib=N (NAG D01GAF integrates whole specified range)
        e    = 0.0d0
        intP = 0.0d0
        s    = 0.0d0

        ! For full-range case, Gill-Miller reduces to:
        !   j = 3, k = N-1 (1-based) with special handling inside loop.
        do i = 3, N-1

          if (i == 3) then
            ! Build initial divided differences around points 1..4 (forward start)
            h2 = X(2) - X(1)
            d3 = (Y(2) - Y(1)) / h2

            h3 = X(3) - X(2)
            d1 = (Y(3) - Y(2)) / h3

            h1 = h2 + h3
            d2 = (d1 - d3) / h1

            h4 = X(4) - X(3)
            r1 = (Y(4) - Y(3)) / h4
            r2 = (r1 - d1) / (h4 + h3)

            h1 = h1 + h4
            r3 = (r2 - d2) / h1

            ! ia == 1 => include first interval contribution using forward formula
            intP = h2 * ( Y(1) + h2 * ( d3/2.0d0 - h2 *
     *                  ( d2/6.0d0 - (h2 + 2.0d0*h3)*r3/12.0d0 ) ) )

            s = -(h2**3) *
     * ( h2*(3.0d0*h2 + 5.0d0*h4) + 10.0d0*h3*h1 ) / 60.0d0

          else
            ! Advance divided differences (centered interior)
            h4 = X(i+1) - X(i)
            r1 = (Y(i+1) - Y(i)) / h4

            r4 = h4 + h3
            r2 = (r1 - d1) / r4

            r4 = r4 + h2
            r3 = (r2 - d2) / r4

            r4 = r4 + h1
            r4 = (r3 - d3) / r4
          end if

          ! Include interval (i-1, i) for i=2..N (here i runs 3..N-1 so covers 2..N-2),
          ! plus end-handling later.
          ! Condition in original: if (i <= ib && i > ia)
          if (i <= N .and. i > 1) then
            intP = intP + h3 *
     * ( (Y(i) + Y(i-1))/2.0d0 - (h3**2) * (d2 + r2 +
     * (h2 - h4)*r3) / 12.0d0 )

            c = (h3**3) *
     * ( 2.0d0*(h3**2) + 5.0d0*( h3*(h4 + h2) +
     *   2.0d0*h4*h2 ) ) / 120.0d0
            e = e + (c + s) * r4

            if (i == 3) then
              s = 2.0d0*c + s
            else
              s = c
            end if
          else
            e = e + r4*s
          end if

          if (i == N-1) then
            ! Backward end formula adds the last interval contribution (N-1, N)
            intP = intP + h4 *
     * ( Y(N) - h4 * ( r1/2.0d0 + h4 * ( r2/6.0d0 +
     * (2.0d0*h3 + h4)*r3/12.0d0 ) ) )

            e = e - (h4**3) * r4 * ( h4*(3.0d0*h4 + 5.0d0*h2) +
     * 10.0d0*h3*(h2 + h3 + h4) ) / 60.0d0

            ! ib >= N-1 is true for full range, so add final term
            e = e + s*r4
          else
            ! Shift state for next step
            h1 = h2
            h2 = h3
            h3 = h4
            d1 = r1
            d2 = r2
            d3 = r3
          end if

        end do

        ER  = e
        ANS = intP + ER
      end subroutine D01GAF_local

      subroutine E01BAF_local(N, X, Y, LAMDA, C, LCK, WORK, WRK, IFAIL)
      !-----------------------------------------------------------------
      ! E01BAF_local: Cubic spline interpolant in B-spline form (NAG-like)
      !
      ! Input:
      !   N        : number of data points (must be >= 4)
      !   X(1:N)   : strictly increasing abscissae
      !   Y(1:N)   : ordinates
      !   LCK      : declared size of LAMDA and C
      !   WRK      : declared size of WORK
      !
      ! Output:
      !   LAMDA(1:N+4) : knots λ(1..N+4) with
      !       λ(1:4)     = X(1)
      !       λ(5:N)     = X(3),...,X(N-2)
      !       λ(N+1:N+4) = X(N)
      !   C(1:N)       : B-spline coefficients c1..cN in
      !       s(x) = sum_{i=1..N} C(i)*N_i(x)   (degree 3)
      !       (elements C(N+1..LCK) are not used)
      !
      ! WORK:
      !   Uses band matrix storage and vectors; requires WRK >= 10*N + 50
      !
      ! IFAIL:
      !   0 = success
      !   1 = invalid N or sizes (N<4, LCK<N+4, WRK too small)
      !   2 = X not strictly increasing
      !   3 = singular/near-singular system detected
      !-----------------------------------------------------------------
        implicit none
        integer, intent(in)    :: N, LCK, WRK
        real(8), intent(in)    :: X(1:N), Y(1:N)
        real(8), intent(out)   :: LAMDA(1:LCK)
        real(8), intent(out)   :: C(1:LCK)
        real(8), intent(inout) :: WORK(1:WRK)
        integer, intent(out)   :: IFAIL

        integer, parameter :: p = 3              ! degree
        integer, parameter :: k = 4              ! order = p+1
        integer, parameter :: kl = 3, ku = 3     ! band (safe for this construction)
        integer            :: i, j, span, j0, r
        integer            :: ldab, offAB, offRHS
        real(8)            :: u, piv, mult, tol
        real(8)            :: Nvals(0:p)

        ! AB has (kl+ku+1) rows, N cols (LAPACK-style band storage)
        ldab  = kl + ku + 1
        offAB = 1
        offRHS = offAB + ldab*N

        tol = 1.0d-14

        IFAIL = 0

        !---------------------------
        ! Validate sizes/inputs
        !---------------------------
        if (N < 4) then
          IFAIL = 1
          return
        end if
        if (LCK < N+4) then
          IFAIL = 1
          return
        end if
        if (WRK < (ldab*N + N + 50)) then
          IFAIL = 1
          return
        end if
        do i = 2, N
          if (X(i) <= X(i-1)) then
            IFAIL = 2
            return
          end if
        end do

        !---------------------------
        ! Build knot vector LAMDA(1..N+4) exactly as NAG E01BAF describes
        ! λ(1:4)=x1, λ(5..N)=x3..x_{N-2}, λ(N+1:N+4)=xN
        !---------------------------
        LAMDA(1) = X(1)
        LAMDA(2) = X(1)
        LAMDA(3) = X(1)
        LAMDA(4) = X(1)

        if (N > 4) then
          ! Fill λ(5..N) with X(3..N-2). This loop runs only if N>=6.
          do i = 5, N
            LAMDA(i) = X(i-2)
          end do
        end if

        LAMDA(N+1) = X(N)
        LAMDA(N+2) = X(N)
        LAMDA(N+3) = X(N)
        LAMDA(N+4) = X(N)

        ! Clear unused part of C
        do i = 1, LCK
          C(i) = 0.0d0
        end do

        !---------------------------
        ! Build band matrix A(i,j)=B-spline basis value N_j(X(i))
        ! Store in WORK as AB(ldab, N) with:
        !   A(i,j) -> AB(ku+1+i-j, j)
        !---------------------------
        do j = 1, ldab*N
          WORK(offAB + j - 1) = 0.0d0
        end do

        do i = 1, N
          u    = X(i)
          span = find_span(N, p, LAMDA, u)        ! span in [p+1 .. N] (1-based)
          call basis_funs(span, u, p, LAMDA, Nvals)

          j0 = span - p                           ! first nonzero basis index (1-based)
          do r = 0, p
            j = j0 + r
            if (j >= 1 .and. j <= N) then
              call set_band(kl, ku, N, i, j, Nvals(r), WORK(offAB))
            end if
          end do
        end do

        ! RHS = Y
        do i = 1, N
          WORK(offRHS + i - 1) = Y(i)
        end do

        !---------------------------
        ! Solve band system (no pivot): AB * C = RHS
        ! Forward elimination
        !---------------------------
        do j = 1, N-1
          piv = get_band(kl, ku, N, j, j, WORK(offAB))
          if (abs(piv) <= tol) then
            IFAIL = 3
            return
          end if

          do i = j+1, min(N, j+kl)
            mult = get_band(kl, ku, N, i, j, WORK(offAB)) / piv
            call set_band(kl, ku, N, i, j, mult, WORK(offAB))

            ! Update row i in columns j+1..min(N, j+ku)
            do r = j+1, min(N, j+ku)
              call set_band(kl, ku, N, i, r,
     *      get_band(kl, ku, N, i, r, WORK(offAB)) -
     * mult*get_band(kl, ku, N, j, r, WORK(offAB)), WORK(offAB))
            end do

            ! Update RHS
            WORK(offRHS + i - 1) = WORK(offRHS + i - 1) -
     * mult*WORK(offRHS + j - 1)
          end do
        end do

        piv = get_band(kl, ku, N, N, N, WORK(offAB))
        if (abs(piv) <= tol) then
          IFAIL = 3
          return
        end if

        ! Back substitution
        C(N) = WORK(offRHS + N - 1) /
     * get_band(kl, ku, N, N, N, WORK(offAB))
        do i = N-1, 1, -1
          piv = get_band(kl, ku, N, i, i, WORK(offAB))
          if (abs(piv) <= tol) then
            IFAIL = 3
            return
          end if
          WORK(offRHS + i - 1) = WORK(offRHS + i - 1) -
     * sum_upper(i, N, kl, ku, WORK(offAB), C)
          C(i) = WORK(offRHS + i - 1) / piv
        end do

      contains

        integer function find_span(ncoef, deg, t, u)
          implicit none
          integer, intent(in) :: ncoef, deg
          real(8), intent(in) :: t(1:ncoef+deg+1)
          real(8), intent(in) :: u
          integer :: low, high, mid, nlast

          ! ncoef = N, deg = 3, knots length = N+4 = ncoef+deg+1
          nlast = ncoef
          if (u >= t(nlast+1)) then
            find_span = nlast
            return
          end if

          low  = deg + 1
          high = nlast + 1
          mid  = (low + high)/2

          do while (u < t(mid) .or. u >= t(mid+1))
            if (u < t(mid)) then
              high = mid
            else
              low = mid
            end if
            mid = (low + high)/2
          end do
          find_span = mid
        end function find_span

        subroutine basis_funs(span, u, deg, t, Nvals)
          implicit none
          integer, intent(in) :: span, deg
          real(8), intent(in) :: u
          real(8), intent(in) :: t(:)
          real(8), intent(out):: Nvals(0:deg)
          real(8) :: left(1:deg), right(1:deg), saved, temp
          integer :: j, r

          Nvals = 0.0d0
          Nvals(0) = 1.0d0

          do j = 1, deg
            left(j)  = u - t(span+1-j)
            right(j) = t(span+j) - u
            saved = 0.0d0
            do r = 0, j-1
              temp = Nvals(r) / (right(r+1) + left(j-r))
              Nvals(r) = saved + right(r+1)*temp
              saved = left(j-r)*temp
            end do
            Nvals(j) = saved
          end do
        end subroutine basis_funs

        subroutine set_band(kl, ku, n, i, j, val, ab)
          implicit none
          integer, intent(in) :: kl, ku, n, i, j
          real(8), intent(in) :: val
          real(8), intent(inout) :: ab(1:(kl+ku+1)*n)
          integer :: row, idx, ldab
          ldab = kl + ku + 1
          row  = ku + 1 + i - j
          if (row < 1 .or. row > ldab) return
          idx = (j-1)*ldab + row
          ab(idx) = val
        end subroutine set_band

        real(8) function get_band(kl, ku, n, i, j, ab)
          implicit none
          integer, intent(in) :: kl, ku, n, i, j
          real(8), intent(in) :: ab(1:(kl+ku+1)*n)
          integer :: row, idx, ldab
          ldab = kl + ku + 1
          row  = ku + 1 + i - j
          if (row < 1 .or. row > ldab) then
            get_band = 0.0d0
            return
          end if
          idx = (j-1)*ldab + row
          get_band = ab(idx)
        end function get_band

        real(8) function sum_upper(i, n, kl, ku, ab, cvec)
          implicit none
          integer, intent(in) :: i, n, kl, ku
          real(8), intent(in) :: ab(1:(kl+ku+1)*n)
          real(8), intent(in) :: cvec(1:LCK)
          integer :: j
          real(8) :: s
          s = 0.0d0
          do j = i+1, min(n, i+ku)
            s = s + get_band(kl, ku, n, i, j, ab) * cvec(j)
          end do
          sum_upper = s
        end function sum_upper

      end subroutine E01BAF_local

      subroutine E02BBF_local(N, LAMDA, C, XVAL, SVAL, IFAIL)
      !------------------------------------------------------------
      ! Evaluate a cubic spline in B-spline form (as produced by E01BAF)
      ! at a single point XVAL using de Boor's algorithm.
      !
      ! Input:
      !   N            : number of B-spline coefficients (must be >= 4)
      !   LAMDA(1:N+4) : knot vector (nondecreasing), cubic => order 4
      !   C(1:N)       : B-spline coefficients
      !   XVAL         : evaluation point
      !
      ! Output:
      !   SVAL         : spline value S(XVAL)
      !
      ! IFAIL:
      !   0 = success
      !   1 = invalid N
      !   2 = invalid knot vector (not nondecreasing)
      !   3 = XVAL outside supported knot span
      !------------------------------------------------------------
        implicit none
        integer, intent(in)  :: N
        real(8), intent(in)  :: LAMDA(1:N+4), C(1:N), XVAL
        real(8), intent(out) :: SVAL
        integer, intent(out) :: IFAIL

        integer, parameter :: p = 3          ! degree
        integer, parameter :: k = 4          ! order = p+1
        integer :: i, j, r, s, span, j0
        real(8) :: d(0:p), alpha, denom
        real(8) :: u

        IFAIL = 0
        SVAL  = 0.0d0

        if (N < 4) then
          IFAIL = 1
          return
        end if

        ! Check knots nondecreasing
        do i = 1, N+3
          if (LAMDA(i+1) < LAMDA(i)) then
            IFAIL = 2
            return
          end if
        end do

        u = XVAL

        ! Supported domain (open clamped knot vector): [LAMDA(4), LAMDA(N+1)]
        if (u < LAMDA(k) .or. u > LAMDA(N+1)) then
          IFAIL = 3
          return
        end if

        ! Find span: span in [k .. N] such that LAMDA(span) <= u < LAMDA(span+1),
        ! with special case u == LAMDA(N+1) -> span = N
        span = find_span(N, p, LAMDA, u)

        ! Local coefficients d(0..p) = C(j0..j0+p), where j0 = span-p
        j0 = span - p
        do r = 0, p
          d(r) = C(j0 + r)
        end do

        ! de Boor recursion
        do r = 1, p
          do j = p, r, -1
            denom = LAMDA(j0 + j + k - r) - LAMDA(j0 + j)
            if (denom == 0.0d0) then
              alpha = 0.0d0
            else
              alpha = (u - LAMDA(j0 + j)) / denom
            end if
            d(j) = (1.0d0 - alpha)*d(j-1) + alpha*d(j)
          end do
        end do

        SVAL = d(p)

      contains

        integer function find_span(ncoef, deg, t, u)
          implicit none
          integer, intent(in) :: ncoef, deg
          real(8), intent(in) :: t(1:ncoef+deg+1)
          real(8), intent(in) :: u
          integer :: low, high, mid, nlast

          nlast = ncoef
          if (u >= t(nlast+1)) then
            find_span = nlast
            return
          end if

          low  = deg + 1        ! = k
          high = nlast + 1
          mid  = (low + high)/2

          do while (u < t(mid) .or. u >= t(mid+1))
            if (u < t(mid)) then
              high = mid
            else
              low = mid
            end if
            mid = (low + high)/2
          end do
          find_span = mid
        end function find_span

      end subroutine E02BBF_local


      subroutine CALCJAC_C(fcn, M, N, X, FJAC)
        !--------------------------------------------------------------------
        ! Compute Jacobian J(M,N) of a vector function FCN at point X(1:N)
        ! using centered finite differences:
        !
        !   J(i,j) ≈ [ F_i(x + h e_j) - F_i(x - h e_j) ] / (2h)
        !
        ! Inputs:
        !   fcn   : user routine  subroutine fcn(m,n,x,f,iflag)
        !   M     : number of function components
        !   N     : number of variables
        !   X(N)  : evaluation point
        !
        ! Output:
        !   FJAC(M,N) : Jacobian matrix
        !
        ! Notes / fixes vs original:
        !   - Uses an explicit interface for FCN (type safety)
        !   - Uses array assignment xh = x instead of manual copy loop
        !   - Step size uses a standard scale: h = sqrt(eps)*max(1,|xj|)
        !--------------------------------------------------------------------
        implicit none

        integer, intent(in)  :: M, N
        real(8), intent(in)  :: X(1:N)
        real(8), intent(out) :: FJAC(1:M,1:N)

        interface
          subroutine fcn(m, n, x, f, iflag)
            implicit none
            integer, intent(in)  :: m, n
            real(8), intent(in)  :: x(1:n)
            real(8), intent(out) :: f(1:m)
            integer, intent(inout) :: iflag
          end subroutine fcn
        end interface

        integer :: i, j, iflag
        real(8) :: xh(1:N), fh(1:M), fl(1:M)
        real(8) :: h, eps
        ! Baseline scale for central differences
        eps = sqrt(epsilon(1.0d0))

        do j = 1, N
          ! Scaled step for variable j
          h = eps * max(1.0d0, abs(X(j)))

          ! x + h e_j
          xh    = X
          xh(j) = X(j) + h
          iflag = 0
          call fcn(M, N, xh, fh, iflag)
          if (iflag /= 0) then
            FJAC(:,j) = 0.0d0
            cycle
          end if

          ! x - h e_j
          xh(j) = X(j) - h
          iflag = 0
          call fcn(M, N, xh, fl, iflag)
          if (iflag /= 0) then
            FJAC(:,j) = 0.0d0
            cycle
          end if

          ! Central difference column j
          do i = 1, M
            FJAC(i,j) = (fh(i) - fl(i)) / (2.0d0*h)
          end do
        end do

      end subroutine CALCJAC_C


      subroutine E04YCF_local(JOB, M, N, FVEC, FJAC, CJ, IWORK,
     *                        RWORK, IFAIL)
        !---------------------------------------------------------------------
        ! Calculate diagonal of covariance matrix for least-squares coefficients.
        !
        ! JOB = 0:
        !   JTJ = FJAC^T * FJAC
        !   inv(JTJ) via LINPACK DGECO/DGEDI (LU-based)
        !   SIG2 = sum(FVEC^2)/(M-N)
        !   CJ(i) = SIG2 * inv(JTJ)(i,i)
        !
        ! Outputs:
        !   CJ(1:N) = variances of coefficients (diagonal only)
        !   IFAIL = 0 success, >0 error
        !---------------------------------------------------------------------
        implicit none

        !-------------------------
        ! Dummy arguments (styled like your M01DAF_local)
        !-------------------------
        integer,          intent(in)    :: JOB, M, N
        real(8),          intent(in)    :: FVEC(1:M)
        real(8),          intent(in)    :: FJAC(1:M,1:N)
        real(8),          intent(out)   :: CJ(1:N)
        integer,          intent(in)    :: IWORK(1:N)   ! kept for compatibility; not used
        real(8),          intent(inout) :: RWORK(1:N)   ! LINPACK workspace
        integer,          intent(out)   :: IFAIL

        !-------------------------
        ! Locals
        !-------------------------
        real(8)  :: JTJ(1:N,1:N)
        integer  :: IPTV(1:N)
        real(8)  :: RCOND
        real(8)  :: DET(2)
        real(8)  :: SIG2, SUM
        integer  :: I, J, K
        integer  :: IJOB

        IJOB  = 11      ! DGEDI: 01 inverse, 10 determinant, 11 both
        IFAIL = 0

        ! Basic checks
        if (JOB /= 0) then
          print *, 'E04YCF_local: JOB != 0 not implemented'
          IFAIL = 2
          return
        end if

        if (M <= N) then
          print *, 'E04YCF_local: Need M > N for SIG2 = SSE/(M-N)'
          IFAIL = 3
          return
        end if

        !------------------------------------------------------------
        ! 1) Form JTJ = FJAC^T * FJAC (symmetric)
        !------------------------------------------------------------
        do I = 1, N
          do J = I, N
            SUM = 0.0D0
            do K = 1, M
              SUM = SUM + FJAC(K,I) * FJAC(K,J)
            end do
            JTJ(I,J) = SUM
            JTJ(J,I) = SUM
          end do
        end do

        !------------------------------------------------------------
        ! 2) Invert JTJ using LINPACK: DGECO + DGEDI
        !------------------------------------------------------------
        call DGECO(JTJ, N, N, IPTV, RCOND, RWORK)

        ! Singularity test (classic LINPACK)
        if (1.0D0 + RCOND == 1.0D0) then
          print *, 'E04YCF_local: JTJ numerically singular (RCOND ~ 0)'
          IFAIL = 4
          return
        end if

        call DGEDI(JTJ, N, N, IPTV, DET, RWORK, IJOB)

        !------------------------------------------------------------
        ! 3) Residual variance SIG2 = sum(FVEC^2)/(M-N)
        !------------------------------------------------------------
        SIG2 = 0.0D0
        do I = 1, M
          SIG2 = SIG2 + FVEC(I)*FVEC(I)
        end do
        SIG2 = SIG2 / dble(M - N)

        !------------------------------------------------------------
        ! 4) Output diagonal: CJ(i) = SIG2 * inv(JTJ)(i,i)
        !------------------------------------------------------------
        do I = 1, N
          CJ(I) = SIG2 * JTJ(I,I)
        end do

      end subroutine E04YCF_local

      subroutine E04YCF_weighed_local(JOB, M, N, FVEC, FJAC, CJ, IWORK,
     *                        RWORK, IFAIL)
        !---------------------------------------------------------------------
        ! Calculate diagonal of covariance matrix for least-squares coefficients.
        !
        ! JOB = 0:
        !   JTJ = FJAC^T * FJAC
        !   inv(JTJ) via LINPACK DGECO/DGEDI (LU-based)
        !   SIG2 = sum(FVEC^2)/(M-N)
        !   CJ(i) = SIG2 * inv(JTJ)(i,i)
        !
        ! Outputs:
        !   CJ(1:N) = variances of coefficients (diagonal only)
        !   IFAIL = 0 success, >0 error
        !---------------------------------------------------------------------
        implicit none

        !-------------------------
        ! Dummy arguments (styled like your M01DAF_local)
        !-------------------------
        integer,          intent(in)    :: JOB, M, N
        real(8),          intent(in)    :: FVEC(1:M)
        real(8),          intent(in)    :: FJAC(1:M,1:N)
        real(8),          intent(out)   :: CJ(1:N)
        integer,          intent(in)    :: IWORK(1:N)   ! kept for compatibility; not used
        real(8),          intent(inout) :: RWORK(1:N)   ! LINPACK workspace
        integer,          intent(out)   :: IFAIL

        !-------------------------
        ! Locals
        !-------------------------
        real(8)  :: JTJ(1:N,1:N)
        integer  :: IPTV(1:N)
        real(8)  :: RCOND
        real(8)  :: DET(2)
        real(8)  :: SIG2, SUM
        integer  :: I, J, K
        integer  :: IJOB

        IJOB  = 11      ! DGEDI: 01 inverse, 10 determinant, 11 both
        IFAIL = 0

        ! Basic checks
        if (JOB /= 0) then
          print *, 'E04YCF_local: JOB != 0 not implemented'
          IFAIL = 2
          return
        end if

        if (M <= N) then
          print *, 'E04YCF_local: Need M > N for SIG2 = SSE/(M-N)'
          IFAIL = 3
          return
        end if

        !------------------------------------------------------------
        ! 1) Form JTJ = FJAC^T * FJAC (symmetric)
        !------------------------------------------------------------
        do I = 1, N
          do J = I, N
            SUM = 0.0D0
            do K = 1, M
              SUM = SUM + FJAC(K,I) * FJAC(K,J)
            end do
            JTJ(I,J) = SUM
            JTJ(J,I) = SUM
          end do
        end do

        !------------------------------------------------------------
        ! 2) Invert JTJ using LINPACK: DGECO + DGEDI
        !------------------------------------------------------------
        call DGECO(JTJ, N, N, IPTV, RCOND, RWORK)

        ! Singularity test (classic LINPACK)
        if (1.0D0 + RCOND == 1.0D0) then
          print *, 'E04YCF_local: JTJ numerically singular (RCOND ~ 0)'
          IFAIL = 4
          return
        end if

        call DGEDI(JTJ, N, N, IPTV, DET, RWORK, IJOB)

        !------------------------------------------------------------
        ! 3) Residual variance SIG2 = sum(FVEC^2)/(M-N)
        !------------------------------------------------------------
        SIG2 = 0.0D0
        do I = 1, M
          SIG2 = SIG2 + FVEC(I)*FVEC(I)
        end do
        SIG2 = SIG2 / dble(M - N)
        SIG2 = 1.0D0

        !------------------------------------------------------------
        ! 4) Output diagonal: CJ(i) = SIG2 * inv(JTJ)(i,i)
        !------------------------------------------------------------
        do I = 1, N
          CJ(I) = SIG2 * JTJ(I,I)
        end do

      end subroutine E04YCF_weighed_local


      integer function JMOD(I, J)
        implicit none
        integer, intent(in) :: I, J

        JMOD = mod(I, J)
      end function JMOD


      subroutine LINTERP(N, X, Y, xval, yval)
        !------------------------------------------------------------
        ! Linear interpolation on tabulated data (X(i),Y(i)), i=1..N.
        !
        ! Behavior:
        !   - If xval is outside [X(1), X(N)], returns yval = 0.0
        !   - Otherwise finds i such that X(i) <= xval < X(i+1) and interpolates
        !
        ! Assumes X is strictly increasing.
        !------------------------------------------------------------
        implicit none

        integer, intent(in)  :: N
        real(8), intent(in)  :: X(1:N), Y(1:N)
        real(8), intent(in)  :: xval
        real(8), intent(out) :: yval

        integer :: i, lo, hi, mid
        real(8) :: t

        ! Handle degenerate input
        if (N < 2) then
          yval = 0.0d0
          return
        end if

        ! Out-of-bounds handling (match your original behavior)
        if (xval <= X(1) .or. xval >= X(N)) then
          yval = 0.0d0
          return
        end if

        ! Binary search for interval: find i with X(i) <= xval < X(i+1)
        lo = 1
        hi = N
        do while (hi - lo > 1)
          mid = (lo + hi) / 2
          if (xval >= X(mid)) then
            lo = mid
          else
            hi = mid
          end if
        end do
        i = lo

        ! Linear interpolation
        t = (xval - X(i)) / (X(i+1) - X(i))
        yval = (1.0d0 - t)*Y(i) + t*Y(i+1)

      end subroutine LINTERP


      SUBROUTINE CONVOLVE_SAME(A, B, N, C)
      INTEGER N, I, J, K, FULL
      REAL*8 A(N), B(N), C(N)
      PARAMETER (MAXN = 1000)
      REAL*8 TMP(2*MAXN-1)
      print *, 'CONVOLVOLUTION INCORRECT!!!'
!-- Full convolution: length is 2*N-1
      DO I = 1, 2*N-1
         TMP(I) = 0.0
      ENDDO

      DO I = 1, N
         DO J = 1, N
            TMP(I+J-1) = TMP(I+J-1) + A(I) * B(J)
         ENDDO
      ENDDO

!-- Extract centered "same" part
      FULL = 2*N-1
      DO I = 1, N
         K = I + N - 1   ! Centered index in TMP
         C(I) = TMP(K)
      ENDDO

      RETURN
      END !CONVOLVE_SAME

      subroutine CONVOLVE(A, NA, B, NB, C, NC, IFAIL)
        !------------------------------------------------------------
        ! Linear convolution of two real vectors:
        !   C(k) = sum_{j} A(j) * B(k-j+1)
        !
        ! Inputs:
        !   A(1:NA), B(1:NB)
        !
        ! Outputs:
        !   C(1:NC) where NC = NA+NB-1
        !   IFAIL = 0 success
        !           1 invalid sizes (NA<1 or NB<1)
        !------------------------------------------------------------
        implicit none

        integer, intent(in)  :: NA, NB
        real(8), intent(in)  :: A(1:NA), B(1:NB)
        integer, intent(out) :: NC, IFAIL
        real(8), intent(out) :: C(1:NA+NB-1)

        integer :: i, j, kmin, kmax

        IFAIL = 0

        if (NA < 1 .or. NB < 1) then
          NC = 0
          IFAIL = 1
          return
        end if

        NC = NA + NB - 1
        C  = 0.0d0

        ! For each output index i, only a restricted range of j contributes.
        do i = 1, NC
          kmin = max(1, i - NB + 1)   ! smallest j such that (i-j+1) <= NB
          kmax = min(NA, i)           ! largest  j such that (i-j+1) >= 1
          do j = kmin, kmax
            C(i) = C(i) + A(j) * B(i - j + 1)
          end do
        end do

      end subroutine CONVOLVE


      subroutine CLEAN_STRING(S, NKEPT)
  !-------------------------------------------------------------
  ! Remove non-printable ASCII from S (keep codes 32..126 only).
  !
  ! On return:
  !   S     = cleaned string (same declared length, padded with blanks)
  !   NKEPT = number of characters kept (logical length)
  !-------------------------------------------------------------
      implicit none

      character(len=*), intent(inout) :: S
      integer,          intent(out)   :: NKEPT

      integer :: i, j, c, l

      l = len(S)
      j = 0

      do i = 1, l
          c = iachar(S(i:i))                 ! ASCII code of this character
          if (c >= 32 .and. c <= 126) then   ! keep printable ASCII only
               j = j + 1
               S(j:j) = S(i:i)
          end if
      end do

  ! Pad remainder with blanks
      if (j < l) S(j+1:l) = ' '

      NKEPT = j
      end subroutine CLEAN_STRING

      subroutine REMOVE_SPACES(OS, IS)
        !----------------------------------------------------------
        ! Remove all space characters (' ') from input string IS and
        ! write the result to OS. Remaining characters in OS are
        ! padded with blanks.
        !----------------------------------------------------------
        implicit none

        character(len=*), intent(out) :: OS
        character(len=*), intent(in)  :: IS

        integer :: i, j, l1, l2

        l1 = len(IS)
        l2 = len(OS)
        j  = 0

        ! Clear output first (ensures fully padded even if OS shorter)
        OS = ' '

        ! Copy only non-space characters, up to OS length
        do i = 1, l1
          if (IS(i:i) /= ' ') then
            if (j < l2) then
              j = j + 1
              OS(j:j) = IS(i:i)
            else
              exit
            end if
          end if
        end do

      end subroutine REMOVE_SPACES


      double precision function G05DDF_local(A, B)
      !-------------------------------------------------------------
      ! Replacement for NAG G05DDF:
      !   Returns a pseudo-random DOUBLE PRECISION number from
      !   Normal(mean=A, std=B) using:
      !     - SLATEC FNLIB RAND(R) for U(0,1)
      !     - Box–Muller transform for N(0,1)
      !   Caches the second variate for efficiency.
      !
      ! Notes:
      !   - Requires SLATEC rand.f linked in.
      !   - B < 0 is treated as an error; this function will STOP.
      !   - Seed SLATEC RAND once before use (see seed helper below).
      !-------------------------------------------------------------
            implicit none
            double precision, intent(in) :: A, B
            real, external :: rand

            double precision :: u1, u2, r, theta, z
            double precision, parameter :: pi =
     *  3.1415926535897932384626433832795d0

            logical :: has_spare
            double precision :: spare
            save has_spare, spare
            data has_spare /.false./

            integer :: tries

            if (B .lt. 0.0d0) then
                print *, "ERROR: G05DDF_local: B < 0"
                stop 1
            end if

            if (has_spare) then
                has_spare = .false.
                G05DDF_local = A + B * spare
                return
            end if

            tries = 0
   10       continue
            tries = tries + 1
            if (tries .gt. 100000) then
                print *,
     * "ERROR: G05DDF_local: RAND returned U<=0 too often"
                stop 2
            end if

            u1 = dble( rand(0.0) )
            if (u1 .le. 0.0d0) goto 10
            u2 = dble( rand(0.0) )

            r     = dsqrt( -2.0d0 * dlog(u1) )
            theta = 2.0d0 * pi * u2

            z     = r * dcos(theta)
            spare = r * dsin(theta)
            has_spare = .true.

            G05DDF_local = A + B * z
            return
            end function G05DDF_local
