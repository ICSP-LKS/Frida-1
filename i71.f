C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i71   :     Fourier transforms
C
C  ====================================================================

C      Contents :
C         TraFourier, TraFilon

C      Aenderungsverzeichnis :
C         JWu   nov91 : Ida5.182 -> Ida7
C         JWu   sep91 : slow FT (non equidistant)
C         JWu   jul91 : new FFT routine
C         JWu   jun91 : -> Modul Ida5; Trafo w<->tof
C         JWu   nov90 : FT begonnen
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE TraFFTsingle (nJlist, JList, Fehler)
C  --------------------------------------------------------------------
            ! JWu completely new jul91. Remained always a weak point.
            ! Small repairs 19feb92. Simplifications for non-Pairs 11jul95.
            ! The mathematics in two separate subroutines 9nov96
            ! Pairs and Non-Pairs separately 10nov96

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*)  Fehler
      INTEGER        JList(*)
      PARAMETER    (MC2=2*MC)
      DIMENSION     DD(MC2), DD1(MC2), DD2(MC2), YY(MC2), Z1(MZ)
      CHARACTER*80  aus
      DATA          tol /1.d-1/

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop over files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL FourierUnits (j, jout, qT2F, fac, Fehler)
         IF (qT2F) THEN
            CALL OlfComAdd (jout, 'F', 'F-trafo H(t)->R(w)', Fehler)
         ELSE
            CALL OlfComAdd (jout, 'F', 'F-trafo H(t)->R(w)', Fehler)
            ENDIF

         DO K = 1, nK

            ncut = 0
 10         CONTINUE ! loop FFT -> modification of nC
            CALL TakeTwoSpectra
     * (j, 0, K, K, nZ, Z1, nC, X1, Y1, Y2, DD1, DD2, tol, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            nC = nC - ncut

            CALL CheckScale (nC, X1, 1.d-1, qEDi, dx)
            IF (.not.qEDi .or. dx.eq.0.) THEN
               Fehler = ' x-scales not equidistant'
               RETURN
               ENDIF

            IF (qT2F) THEN
               CALL TraFourierT2F (MC, nC, nCt, ncut, 0, 2, dx, fac,
     *                             X1, Y1, Y2, DD, DD1, DD2, Fehler)
cdeb               type *, ' T2F cut : ', ncut
            ELSE
               iOdd = 0
               CALL TraFourierF2T (MC, nC, nCt, ncut, iOdd, 2, dx, fac,
     *                           X1, YY, Y1, Y2, DD, DD1, DD2, Fehler)
cdeb               type *, ' F2T cut : ', ncut
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            IF (ncut.ne.0) GOTO 10

            CALL OlfPutSpe (jout, K, nZ, Z1, nCt, X1, Y1, DD1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO
            CALL OlfClos (jout, nK, Fehler)

         ENDDO ! lj

      END ! TraFFTsingle

C  --------------------------------------------------------------------
      SUBROUTINE TraFFTpair (nJlist, JList, Fehler)
C  --------------------------------------------------------------------
            ! JWu completely new jul91. Remained always a weak point.
            ! Small repairs 19feb92. Simplifications for non-Pairs 11jul95.
            ! The mathematics in two separate subroutines 9nov96
            ! Dialogue simplified (uncomplete pairs disallowed) 10nov96

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*)  Fehler
      INTEGER        JList(*)
      PARAMETER    (MC2=2*MC)
      DIMENSION     DD(MC2), DD1(MC2), DD2(MC2), YY(MC2), Z1(MZ)
      CHARACTER*80  aus

      DATA          tol /1.d-1/

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Direction of transform :
      qH2R = qAskDi (' Real->Hermitian(0) or Hermitian->Real(1)',
     *               intq(qH2R))
         ! ' Hermitian data are complex(3)'
         iIn = 3
            iOut = 0
         ! ' Real data are full(0)'
         iIn = 0
         iOut = 3

C  Loop over files :
      DO lj = 1, nJList

         ! one file from list, the other must be asked for :
         j1 = JList(lj)
         IF (qH2R) THEN
            j2 = iAsk (' Imaginary part for real file '//cl3(j1))
         ELSE
            j2 = 0
            ENDIF

c96/7         CALL ParJoin (j1, j2, 1, ..)

         CALL OlfHeadDup (j1, .false., j3, nK, Kout3, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! coordinate names and units (recopied from TraFilon 11jul95) :
         CALL FourierUnits (j1, j3, qT2F, fac, Fehler)
         IF (qH2R .neqv. qT2F) THEN
            CALL Gong (2)
            Print *, 'WARNING/ qH2R <> qT2F'
            ENDIF

         IF (qH2R) THEN
            CALL OlfComAdd (j3, 'F',  'F-trafo H->R', Fehler)
         ELSE
            CALL OlfComAdd (j3, 'F',  'F-trafo R->H', Fehler)
            ENDIF

         !  separate documentation for files 1, 2 :
         IF (.not.qH2R) THEN
            CALL OlfHeadDup (j3, .false., j4, nK, Kout4, Fehler)
            IF (Fehler.ne.'&ff') RETURN
c96/7       append _r, _i !!
         ELSE
            j4 = 0
            ENDIF

         DO K = 1, nK

            ncut = 0
 10         CONTINUE ! loop FFT -> modification of nC
            CALL TakeTwoSpectra
     * (j1, j2, K, K, nZ, Z1, nC, X1, Y1, Y2, DD1, DD2, tol, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            nC = nC - ncut

            CALL CheckScale (nC, X1, 1.d-1, qEDi, dx)
            IF (.not.qEDi .or. dx.eq.0.) THEN
               Fehler = ' x-scales not equidistant'
               RETURN
               ENDIF

            IF (qH2R) THEN
               CALL TraFourierT2F (MC, nC, nCt, ncut, 0, iOut, dx, fac,
     *                             X1, Y1, Y2, DD, DD1, DD2, Fehler)
            ELSE
               iOdd = 0
               CALL TraFourierF2T (MC, nC, nCt, ncut, iOdd, iIn, dx,
     *                      fac, X1, YY, Y1, Y2, DD, DD1, DD2, Fehler)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            IF (ncut.ne.0) GOTO 10

            CALL OlfCopZ   (j, j3, K, K, Fehler)
            CALL OlfPutXYD (j3, K, nCt, X1, Y1, DD1, Fehler)
            IF (j4.ne.0) THEN
               CALL OlfCopZ   (j, j4, K, K, Fehler)
               CALL OlfPutXYD (j4, K, nCt, X1, Y2, DD2, Fehler)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (j3, nK, Fehler)
         IF (j4.ne.0) CALL OlfClos (j4, nK, Fehler)

         ENDDO ! lj

      END ! TraFFTpair

      SUBROUTINE TraFourierT2F (MC, nin, nout, ncut, iOdd, iOut, dx,
     *                          fac, X1, Y1, Y2, DD, DD1, DD2, Fehler)
C     -------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler
      INTEGER       MC, iOdd, iOut, nin, nout, N, L, i,
     *              nsoll1, nsoll2, ncut, nsoll
      REAL*8        xofft, dxt, twopi, fac, dx, dsum,
     *              X1(*), Y1(*), Y2(*), DD(*), DD1(*), DD2(*)
      DATA          twopi /6.2831853/

      N = nin
      L = 2*N - 2 + iOdd
      nout = L
      IF (L.gt.MC) THEN
         Fehler = ' Wo kommen soviele Kanaele her ?'
         RETURN
         ENDIF
      ! new x-scale :
      xofft = 0 ! provis___
      dxt   = twopi / (fac * dx * L)
      DO i = 1, N
         X1(i) = xofft + (i-1) * dxt ! positive half
         ENDDO
      DO i = N+1, L
         X1(i) = xofft + (i-1-L) * dxt ! continue on negative side
         ENDDO
      ! y-data into hermitian sequence :
      DO i = 1, L-N
         Y1 (L+1-i) = Y2 (1+i)  ! here no minus sign
         ENDDO
      ! y-Trafo :
      CALL TraFT (Y1, L, 1/(L*dxt), 'FBF', nsoll, Fehler)
      IF (nsoll.ne.0) THEN
         ncut = (L-nsoll) / 2
cdeb         Print *, ' L nsoll N ', L, nsoll, N
         RETURN
      ELSE
         ncut = 0
         ENDIF
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'data transform/ ')
         RETURN
         ENDIF
      ! error transform (JWu 29/30jul91) :
      ! transform squared errors :
      DO i = 1, N
         DD1(i) = DD1(i)**2
         DD2(i) = DD2(i)**2
         ENDDO
      DO i = N+1, L
         DD1(i) = 0.
         DD2(i) = 0.
         ENDDO
      dsum = 0.
      DO i = 2, (L+1)/2
         dsum = dsum + DD1(i) + DD2(i)
         ENDDO
      ! FT as of real data :
      CALL TraFT (DD1, L, 1.d0, 'FBF', nsoll1, Fehler)
      CALL TraFT (DD2, L, 1.d0, 'FBF', nsoll2, Fehler)
      IF (nsoll1.ne.0 .or. nsoll2.ne.0) THEN
         Fehler = ' nsoll<>0 after error transform'
         RETURN
         ENDIF
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'error transform/ ')
         RETURN
         ENDIF
      ! extract transform with double period :
      CALL rCopy (DD,  1, L+1-N, 1, DD1, 1,       2)
      CALL rCopy (DD,  L+2-N, L, 1, DD1, L+3-2*N, 2)
      CALL rCopy (DD1, 1, L,     1, DD,  1,       1)
      CALL rCopy (DD,  1, L+1-N, 1, DD2, 1,       2)
      CALL rCopy (DD,  L+2-N, L, 1, DD2, L+3-2*N, 2)
      CALL rCopy (DD2, 1, L,     1, DD,  1,       1)
      DO i = 1, L
         DD1(i) = dsqrt0 (2*dsum + DD1(i) + DD2(i)) / (L*dxt)
         ENDDO

      IF     (iOut.eq.0) THEN
         ! full wave : everything is prepared in X(..) :
         CALL SortChannels (X1, Y1, DD1, L, .true.)
      ELSEIF (iOut.eq.2) THEN
         ! only real half wave :
         L = N ! throw away the rest (simplest solution, 26oct93)
      ELSE
         Fehler = 'split output/ unfertich'
         RETURN
         ENDIF

      END ! TraFourierT2F

      SUBROUTINE TraFourierF2T (MC, nin, nout, ncut, iOdd, iIn, dx,
     *                     fac, X1, YY, Y1, Y2, DD, DD1, DD2, Fehler)
C     ------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler
      INTEGER       MC, iOdd, iIn, nin, nout, N, L, ncut, nsoll, i
      REAL*8        xofft, dxt, twopi, fac, dx, dsum,
     *              X1(*), YY(*), Y1(*), Y2(*), DD(*), DD1(*), DD2(*)
      DATA          twopi /6.2831853/

      IF (iIn.ne.0) THEN
      ! real half wave -> hermitian component (corrected 27aug92) :
         N   = nin
         L   = 2*N - 2 + iOdd
         DO i = 2, L-N+1  ! set Y(L)..Y(N+1) := Y1(2)..Y2(N-1+iodd)
            YY (L+2-i) = Y1(i) - Y2(i)
            DD1(L+2-i) = dsqrt0 ( DD1(i)**2 + DD2(i)**2 )
            ENDDO
         YY(N) = Y1(N)
         DD(N) = DD1(N)
         DO i = 1, L-N+1  ! set Y(1)..Y(N-1+iodd)
            YY (i) = Y1(i) + Y2(i)
            DD1(i) = dsqrt0 ( DD1(i)**2 + DD2(i)**2 )
            ENDDO
      ELSE
      ! real full wave -> hermitian (JWu 28jul91) :
         L   = nin
         N   = L/2+1
         iOdd= mod(L,2)
         CALL rCopy (YY, 1, L, 1, Y1, 1, 1)
         ENDIF
      nout= N
      ! new x-scale :
      xofft = 0.
      dxt   = twopi / (fac * dx * L)
      DO i = 1, N
         X1(i) = xofft + (i-1) * dxt
         ENDDO
      ! y-Trafo :
      CALL TraFT (YY, L, dx, 'FAF', nsoll, Fehler)
      IF (nsoll.ne.0) THEN ! try it again
         IF (iIn.ne.0) THEN
            ncut = (L - nsoll)/2
         ELSE
            ncut = L - nsoll ! cured endless loop 11jul94
            ENDIF
cdeb      type *, 'N L nsoll ncut ', N, L, nsoll, ncut
         RETURN
      ELSE
         ncut = 0
         ENDIF
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'data transform/ ')
         RETURN
         ENDIF
      ! decompose the hermitian sequence :
      DO i = 0, L/2
         Y1(1+i) = YY(1+i)
         ENDDO
      Y2(1)   = 0.
      Y2(N) = 0.
      DO i = 1, (L-1)/2
         Y2 (1+i) = - YY (L+1-i)    ! here a minus sign
         ENDDO
      ! error transform (JWu 29/30jul91) :
      ! transform squared errors :
      dsum = 0.
      DO i = 1, L
         DD1(i) = DD1(i)**2
         dsum  = dsum + DD1(i)
         ENDDO
      ! copy D1 into D2 with twice the period :
      DO i = 2, L, 2
         DD2(i) = 0.
         ENDDO
      DD2(1) = DD1(1)
      DO i = 3, L, 2
         DD2(i) = DD1((i+1)/2) + DD1(L+2-(i+1)/2)
         ENDDO
      IF (.not.qintr(iOdd)) DD2(L-1) = DD1(N) ! overwrite
      CALL TraFT (DD1, L, 1.d0, 'FAF', nsoll, Fehler)
         IF (nsoll.ne.0) THEN
            Fehler = ' nsoll<>0 after error transform'
            RETURN
            ENDIF
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'error transform/ ')
         RETURN
         ENDIF
      ! errors of Re and Im :
      DO i = 1, N
         DD1(i) = dsqrt0 ((dsum+DD2(i))/2) * dx
         DD2(i) = dsqrt0 ((dsum-DD2(i))/2) * dx
         ENDDO

      END ! TraFourierF2T

C  --------------------------------------------------------------------
      SUBROUTINE TraFT (Y, n, fact, NAG, nsoll, Fehler)
C  --------------------------------------------------------------------
            ! JWu 29jul91
         ! Call to a NAG Fourier transform routine;
         ! remove the stupid factor of sqrt(n), and
         ! multiply the result by fact.
         ! If n is too odd, return a better nsoll.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      PARAMETER    (MC2=2*MC)
      DIMENSION     Y(MC2), Work(MC2+15), Trig(MC2)    !Artem Work(MC2) -> Work(MC2+15) for rfftf()
      CHARACTER     NAG*(*), Fehler*(*), cFirstCall*1, cl4*4

      DATA          nOld /-1/

      nsoll = 0 ! normal status

      IF (n.le.0) CALL Absturz ('TraFT', 'n.le.0')
      IF (n.eq.nOld) THEN
         cFirstCall = 's'
      ELSE
         cFirstCall = 'i'
         ENDIF

      ifail = 1 ! silent exit

c      DO i=1,n
c         Print *, '>> FAF in  ', i, Y(i)
c         ENDDO
      IF     (NAG.eq.'FAF') THEN
         CALL RFFTI (n, Work)    !Artem: add for rfftf()
         CALL RFFTF (n, Y, Work) !Artem: Replace with rfftf function from slatec/fishfft/rfftf.f:  CALL C06FAF (Y, n, Work, ifail)
         IF (ifail.eq.1) GOTO 91
      ELSEIF (NAG.eq.'FBF') THEN
         CALL CFFTI (n, Work)
         CALL CFFTF (n, Y, Work)  !Artem: Replace with cfftf function from slatec/fishfft/rfftf.f:  CALL C06FBF (Y, n, Work, ifail)
         IF (ifail.eq.1) GOTO 91
      ELSE
          CALL Absturz ('TraFT', 'parameter NAG o.o.r.')
          ENDIF
c      DO i=1,n
c         Print *, '>> FAF out ', i, Y(i)
c         ENDDO

      IF (ifail.ne.0) THEN
         Fehler = ' Failure in Fourier transform C06'//NAG//
     *      ' : ifail = '//cl4(ifail)
         RETURN
         ENDIF

C  Correct for prefactor :
      DO i = 1, n
         Y(i) = dsqrt(dble(n)) * fact * Y(i)
         ENDDO

      RETURN ! normal exit

C  Try to find a better number of spectra (=nsoll since 27aug92) :
 91   CONTINUE
      DO nsoll = n-2, n/2, -2
         ifail = 1 ! silent exit
         CALL RFFTI (nsoll, Work)    !Artem: add for rfftf()
         CALL RFFTF (nsoll, Trig, Work) !Artem: Replace with rfftf function from slatec/fishfft/rfftf.f:  CALL C06FAF (Trig, nsoll, Work, ifail)
C         CALL C06FAF (Trig, nsoll, Work, ifail) ! Trig as dummy; don't change Y
         IF (ifail.eq.0) RETURN ! found good nsoll
         ENDDO
      Fehler = ' FFT/ cannot find a good number of channels'

      END ! TraFT

C  --------------------------------------------------------------------
      SUBROUTINE SetGridFT (n1, X1, n2, X2, qBack, fac, Fehler)
C  --------------------------------------------------------------------
            ! JWu 12sep91

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      DIMENSION     X1(MC), X2(MC)
      CHARACTER     Fehler*(*)

C  Analyse X1 :
      IF (irSorted(X1,n1).lt.2) THEN
         Fehler = ' x scale is not sorted'
         RETURN
         ENDIF
      IF (X1(2).le.0.) THEN
         Fehler = ' x scale contains non-positive points'
         RETURN
         ENDIF
      dx11  = X1(2) - X1(1)
      dx1mi = dx11            ! minimum
      dx1ma = dx11            ! maximum
      dx1av = dx11            ! average absolute step width
      dx1ra = dx11/X1(2)      ! average relative step width
      DO i = 1, n1-1
         dx = X1(i+1) - X1(i)
         dx1mi = dmin1 (dx1mi, dx)
         dx1ma = dmax1 (dx1ma, dx)
         dx1av = dx1av + dx
         dx1ra = dx1ra + dx/X1(i+1)
         ENDDO
      dx1av = dx1av / n1
      dx1ra = dx1ra / n1

      dk1mi = fac / dx1mi
      dk1ma = fac / dx1ma
      dk1av = fac / dx1av

      Print *, ' reciprocal of original x-scale characteristics :'
      Print '(2(a,g10.2,2x))',
     *    '    2pi/dx_min = ', dk1mi, ' 2pi/dx_max = ', dk1ma
      Print '(2(a,g10.2,2x))',
     *    '    2pi/ <dx>  = ', dk1av, ' <dx/x>     = ', dx1ra

C  Ask for new scale :
      IF ( x2ma.le.0.)  x2ma = dk1mi
c????      IF (dx2mi.le.0.) dx2mi = dx2ma/30  !  FEHLER : dx2ma not assigned
      Print *, ' new scale :'
      x2ma  = rAskD (' Maximum x', x2ma)
      delta = rAskD (' Delta (crossover linear - logarithmic)', delta)
      IF (x2ma.le.0. .or. delta.le.0.) THEN
         Fehler = ' '
         RETURN
         ENDIF
      n2    = iAskDMu (' Number of points', n2, 3, MC)

C  Set the scale :
         ! log (x+1) - scale : see Protokollbuch C2,31.
      alpha = dlog (x2ma/delta + 1) / (n2 - 1)
      X2(1) = 0.
      DO i = 2, n2
         X2(i) = delta * (dexp(alpha*(i-1)) - 1)
         ENDDO

      END ! SetGridFT

C  --------------------------------------------------------------------
      SUBROUTINE TraSlowFT ( n1, X1, Y1, n2, X2, Y2,
     *                       qBack, fac, Fehler)
C  --------------------------------------------------------------------
            ! JWu 12sep91
         ! Fourier transform with Y1 approximated by a polygon.
         ! See fax by M.Kiebel(26jul91) and Protokollbuch C2,30.
         ! ---- out of use ----

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      DIMENSION     X1(MC), Y1(MC), X2(MC), Y2(MC)
      CHARACTER     Fehler*(*)

      DATA          twopi /6.2831853d0/

      DO ii = 1, n2
         t  = X2(ii)
         th = X2(ii) / fac

         IF (t.eq.0.) THEN
            s = Y1(n1) * X1(n1) - Y1(1) * X1(1)
            DO i = 1, n1-1
               s = s + Y1(i)*X1(i+1) - Y1(i+1)*X1(i)
               ENDDO

         ELSE
            s = Y1(n1) * dsin(X1(n1)*th) - Y1(1)*dsin(X1(1)*th)
            DO i = 1, n1-1
               s = s + (Y1(i+1)-Y1(i)) *
     *          (dcos(X1(i+1)*th)-dcos(X1(i)*th)) /
     *          ( (X1(i+1)-X1(i)) * th )
               ENDDO
            s = s * 2 / th
            ENDIF

         IF (.not.qBack) THEN
            Y2(ii) = s
         ELSE
            Y2(ii) = s / (twopi * fac) ! siehe Protokollbuch C2,30
            ENDIF

         CALL Counter (ii, 20, n2, ' ')

         ENDDO

      END ! TraSlowFT

C  --------------------------------------------------------------------
      SUBROUTINE TraFilon (nJlist, JList, qOv, qCos, Fehler)
C  --------------------------------------------------------------------
            ! JWu 27-apr93. The algorithm itself by courtesy of M.Fuchs
            ! qCos : cosine or sine transform ? 24jan95

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      PARAMETER    (MLB=128)

      DIMENSION     I1B(MLB)

      CHARACTER     doc*40, aus*80

      COMMON / Filon / fico(11), twopi

      DATA          iGriX /4/

C  Vorbelegung einiger haeufig benoetigter Zahlen :
      fico(1) = 2.d0/45.
      fico(2) = -2.d0/315.
      fico(3) = 2.d0/4725.
      fico(4) = 2.d0/3.
      fico(5) = 2.d0/15.
      fico(6) = -4.d0/105.
      fico(7) = 2.d0/567.
      fico(8) = 4.d0/3.
      fico(9) = -2.d0/15.
      fico(10)= 1.d0/210.
      fico(11)= -1.d0/11340.
      twopi   = 8 * datan(1.d0)

C  Set output grid :
      aus =
     * ' Choose output grid: lin(1) 1/2-log(2) log(3) '//
     * 'lin-blocks(4) lin-FFT(5)'
 12   iGriX = iAskDMu (aus, iGriX, 0, 5)
      IF (iGriX.eq.0) THEN
         Fehler = ' '
         RETURN
      ELSEIF (iGriX.eq.5) THEN
      ELSE
         CALL SetGridReg (iGriX, sg1, sgn, nsg, X1, MC, n1,doc,Fehler)
         IF (Fehler.ne.'&ff') THEN
            CALL FehlerGong (Fehler,1)
            GOTO 12
            ENDIF
         ENDIF

C  Loop over files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL FourierUnits (j, jout, qT2F, fac, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (qCos) THEN
            CALL OlfComAdd (j, 'F', 'Filon-Cos-FT', Fehler)
         ELSE
            CALL OlfComAdd (j, 'F', 'Filon-Sin-FT', Fehler)
            ENDIF

C  Loop spectra :
         q0warn = .false.
         DO K = 1, nK
            CALL Counter (K, 1, nK, '... working hard')

            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

C  Check input grid :
            IF (n.lt.3) THEN
               aus = ' # channels <= 3'
               GOTO 990
               ENDIF
            IF (irSorted(X, n).ne.2) THEN
               aus = ' input grid not in ascending order'
               GOTO 990
               ENDIF
            IF (X(1).lt.0.) THEN
               aus = ' input grid starts with negative value'
               GOTO 990
               ENDIF
            IF (.not.qEqTol(X(1), 0.d0, X(2)/100) .and. .not.q0warn)
     *         THEN
               IF (.not.qAsk (
     *   ' input grid doesn''t start with 0 -- continue ?')) THEN
                  aus = ' agreed on exit'
                  GOTO 990
                  ENDIF
               q0warn = .true.
               ENDIF
C  Find blocking of input data :
            ! block 1 starts at X(1) :
            nB = 1
            I1B(1) = 1
            step = X(2) - X(1)
            DO i = 3, n
               IF (.not.qEqTol(X(i),X(i-1)+step,step/100)) THEN
                  ! new block found
                  IF (i.le.I1B(nB)+2) THEN ! block has less than 3 elements
                     DO ii = max0(i-4,1), i
                        Print '(a,i3,2x,g12.6)', ' i X : ', ii, X(ii)
                        ENDDO
                     aus = ' block too short'
                     GOTO 990
                     ENDIF
                  nB = nB + 1
                  IF (nB.gt.MLB) THEN
                     Fehler = ' too many blocks'
                     RETURN
                     ENDIF
                  I1B(nB) = i-1
                  step = X(i) - X(i-1)
                  ENDIF
               ENDDO
            I1B(nB+1) = n ! end of last block

C  Automatic choice of grid (as for FFT) :
            IF (iGriX.eq.5) THEN
               n1 = n
               step = X(n) - X(1)
               step = twopi / fac / step
               DO i = 1, n1
                  X1(i) = (i-1)*step
                  ENDDO
               ENDIF

C  Perform the Fourier integration :
            CALL IntgrFilon (qCos, MLB, nB, I1B, X, Y, qT2F,
     *                       fac, n1, X1, Y1)
            CALL rSet (D1, 1, n1, 1, 0.d0)

C  Save result :
            CALL OlfCopZ   (j, jout, K, K , Fehler)
            CALL OlfPutXYD (jout, K, n1, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO
C  End loop spectra.

         CALL OlfClos (jout, nK, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO
C  End loop files.

      RETURN
 990  CONTINUE
      CALL Compose3 (Fehler, 'file '//cl3(j), ', spectrum '//
     *               cl3(K), ': '//aus)
      RETURN

      END ! TraFilon

C  --------------------------------------------------------------------
      SUBROUTINE FourierUnits (jin, jout, qT2F, fac, Fehler)
C  --------------------------------------------------------------------
            ! From TraFilon. As a subroutine 11jul95. Use UnitSI 11jun96
         ! Input  :  jin, jout :    file names
         ! Output :  qT2F :        time -> frequency or backwards ?
         !           fac :          conversion factor

      IMPLICIT NONE
      CHARACTER*(*) Fehler
      CHARACTER*40  coX, unX, coY, unY, coi, uni, aux
      REAL*8        fac, faci, faco, twopi, hbar
      INTEGER       jin, jout, intq
      LOGICAL       qT2F, qAskDi
      DATA          twopi /6.2831853/, hbar /0.658218/, unX / ' ' /

      CALL OlfCnuG (jin, 'x', coi, uni, Fehler)
      CALL OlfCnuG (jin, 'y', coY, unY, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  which direction ?
      IF     (coi(1:1).eq.'t') THEN
         qT2F = .true.
      ELSEIF (coi(1:1).eq.'w' .or. coi(1:1).eq.'f') THEN
         qT2F = .false.
      ELSE
         Print *, ' input x coordinate is ', coi
         qT2F = qAskDi (
     *    ' Treat as frequency/energy(0) or as time(1)', intq(qT2F))
         ENDIF

C  which is the input unit ?
      IF (qT2F) THEN
         CALL UnitSI ('time', uni, 1, faci)
      ELSE
         CALL UnitSI ('frequency', uni, 2, faci)
         ENDIF
      IF (faci.le.0) THEN
         Fehler = ' invalid unit'
         RETURN
         ENDIF

C  Output x coordinate and unit :
      IF (qT2F) THEN
         CALL ReplaceT (coX, 't', 'f')
      ELSE
         coX = 't'
         ENDIF
      CALL FrageCD (' Output x coordinate', coX, coX)

      CALL FrageCD (' Output x unit', unX, unX)

      IF (qT2F) THEN
         CALL UnitSI (coX, unX, 2, faco)
      ELSE
         CALL UnitSI (coX, unX, 1, faco)
         ENDIF
      IF (faco.le.0) THEN
         Fehler = ' invalid unit'
         RETURN
         ENDIF

C  Output y unit :
      CALL ReplaceT (coY, coi, coX)
      CALL FrageCD (' Output y coordinate ?', coY, coY)

C  Output y unit :
      CALL Compose2 (aux, uni, '-1')
      IF     (unY.eq.aux) THEN ! [input-y] = [input-x]^-1
         unY = ' ' ! output-y is dimensionless
      ELSEIF (unY.eq.' ') THEN ! input-y is dimensionless
         CALL Compose2 (unY, unX, '-1')
      ELSE
         CALL Append (unY, '*'//uni)
         ENDIF
      CALL FrageCD (' And its unit', unY, unY)

C  conversion factor :
      fac = faci * faco
      IF (fac.le.0) THEN
         Fehler = 'F-Units/ conversion factor <= 0'
         RETURN
         ENDIF

      CALL OlfCnuP (jout, 'x', coX, unX, Fehler)
      CALL OlfCnuP (jout, 'y', coY, unY, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! FourierUnits

C  --------------------------------------------------------------------
      SUBROUTINE IntgrFilon (qCos, MLB, nB, I1B, X, Y, qT2F,
     *                       fac, m, om, YT)
C  --------------------------------------------------------------------

         ! Import :
         !    qCos     cos or sin ?
         !    MLB, nB   Zahl von Bloecken (max / actual)
         !    I1B(nB)  Numbers of first channels of blocks
         !    X        input channels
         !    Y        input data
         !    m        Zahl zu berechnenden Frequenzen
         !    om(m)    diese Frequenzen
         !    qT2F    time -> freq/energy (else  <- )
         !    fac      unit of t*f
         ! Export :
         !    YT       transform of Y

      IMPLICIT NONE

      REAL*8     fmin, eps
      COMPLEX*16 ci,cr,c0
      LOGICAL    qCos

      INCLUDE 'i_dim.f'
      PARAMETER (cr = (1.d0,0.d0), ci=(0.d0,1.d0), c0=(0.d0,0.d0) )
         !    complex numbers 1, i, 0
      PARAMETER (fmin=1.d-10 , eps=1.d-16 )
         !    fmin     Abbruchkriterium, falls f(t) < fmin     ( f(t) <=1 ).

      LOGICAL    qT2F
      INTEGER    MLB, iB,j,in,m,nB,I1B(MLB),i,n2
      REAL*8     X(MC),f(MC),Y(MC),YT(MC),t0,f0, om(MC), fac,
     *           h,alpha,beta,gamma,fico, twopi
      COMPLEX*16 cf(MC),cr0,crut,c2rut,c1,c2,c2nrut,cdum

      COMMON / Filon / fico(11), twopi

C  Transform will be summed in cf(j) :
      DO j = 1 , m
         cf(j) = c0
         ENDDO

C  Loop over blocks :
      DO iB = 1 , nB
         ! transcribe block :
         n2 = I1B(iB+1) - I1B(iB)
         f0 = Y(I1B(iB))
         DO i = 1, n2
            f(i) = Y(I1B(iB)+i)
            ENDDO
         h  = X(I1B(iB)+1) - X(I1B(iB))
         t0 = X(I1B(iB))

         if ( abs(f(1)) .le. fmin ) goto 9 !
c  Beginn der Filon Integration eines einzelnen Teilintegrals a bis b
c  fuer all frequencies om(j)
         DO j = 1 , m
c  Berechnung der Koeffizienten alpha,beta,gamma (siehe Abramowitz)
            CALL FilonCoeff (h*om(j), alpha, beta, gamma)
c  komplexe Berechnung von exp(i*omega*h), um alle noetigen cos und sin durch
c  Multiplikation zu erhalten
            cr0 = om(j) * t0 * ( ci - eps * cr ) * fac
            cr0 = exp( cr0 ) ! ??????????????????????????????????????????
            crut = h * om(j) * ( ci - eps * cr ) * fac
            crut = exp( crut ) ! ????????????????????????????????????????
            c2rut = crut * crut
            c2nrut =  n2 * h * om(j) * ( ci - eps * cr ) * fac
            c2nrut = cdexp( c2nrut )

c  Beginn der Aufsummationen
            c2 = .5 *  f(n2) * cr
            DO in = n2-2 , 2 , -2
               c2 = c2 * c2rut + f(in) * cr
               ENDDO
            c2 = c2 * c2rut + .5 *  f0 * cr
            c1 = f(n2-1) * cr
            DO in = n2 - 3 , 1 , - 2
               c1 = c1 * c2rut + f(in) * cr
               ENDDO
            c1 = c1 * crut
            cdum = gamma * c1 + beta * c2
     *           - alpha * ci * ( f(n2) * c2nrut - f0 * cr )
c  ein einzelnes Teilintegral fuer eine Frequenz ist mit Filon nun berechnet

c die Teilintegrale werden in cf(j) aufsummiert
            cf(j) = cf(j) + h * cdum * cr0
            ENDDO
         ENDDO
    9 CONTINUE

C  The complex cf(j) contains cosine and sine transforms :
      DO j = 1 , m
         IF (qCos) THEN
            IF (qT2F) THEN
               YT(j) = 2 * real (cf(j)) / twopi * fac
               ! factor 2 because result is interpreted as int_{-\infty}^\infty
            ELSE
               YT(j) = 2 * real (cf(j))
               ENDIF
         ELSE
            IF (qT2F) THEN
               YT(j) = 2 * dimag (cf(j)) / twopi * fac
            ELSE
               YT(j) =-2 * dimag (cf(j))
               ENDIF
            ENDIF
         ENDDO

      END ! IntgrFilon

      SUBROUTINE FilonCoeff (dum, alp, bet, gam)
C     ------------------------------------------
         ! calculate the coeffiecients alpha, betha, gamma

      REAL*8 fico, twopi, dum, dumq, alp, bet, gam, ds, dc
      COMMON / Filon / fico(11), twopi

      IF (dum.gt.0.01) THEN ! Abramowitz 25.4.53 : small theta expansion
         ds = sin(dum)
         dc = cos(dum)
         alp = ((dc-2*ds/dum)*ds/dum+1)/dum
         bet = 2*(1+dc*dc-2*ds*dc/dum)/dum/dum
         gam = 4*(ds/dum-dc)/dum/dum
      ELSE                  ! Abramowitz 25.4.52 : general formulae
         dumq = dum * dum
         alp =  ((fico( 3)*dumq+fico( 2))*dumq+fico(1))*dumq*dum
         bet = (((fico( 7)*dumq+fico( 6))*dumq+fico(5))*dumq+fico(4))
         gam = (((fico(11)*dumq+fico(10))*dumq+fico(9))*dumq+fico(8))
         ENDIF

      END ! FilonCoeff
