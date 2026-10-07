C  ====================================================================
C
C     Library  WuL :  General FORTRAN Library
C     Module   L6  :      Subroutines for interactive programs
C
C  ====================================================================

C     6.5.  General Functions
C              FuVal, FuText, FuTxt, FuInv, FuHelp, FuAsk
C
C     6.7.  Functions for solid state physics
C              u2Debye, DebyeInt

C  ====================================================================
C     6.5.  General Functions
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE FuVal (i, y, dy, x, dx, z, dz)
C     -----------------------------------------
            ! JWu 1990.
         ! Calculate y, dy from x, dx, z, dz according to function no. i

      IMPLICIT REAL*8 (a-h,o-p,r-z)

      IF (dx.lt.0.) Print *, 'Warning from FuVal : dx<0'
      IF (dz.lt.0.) Print *, 'Warning from FuVal : dz<0'

C  One argument (i=0,..,49) :
      IF     (i.eq. 0) THEN
         y    = x
         dy   = dx
      ELSEIF (i.eq. 1) THEN
         y    = dln0(x)
         dy   = dquot0(dx,x)
      ELSEIF (i.eq. 2) THEN
         y    = dlg0(x)
         dy   = dquot0 (dx, x*2.3026) ! / ln 10
      ELSEIF (i.eq. 3) THEN
         y    = dexp1(x)
         dy   = y * dx
      ELSEIF (i.eq. 4) THEN
         y    = x**2
         dy   = 2*dabs(x)*dx
      ELSEIF (i.eq. 5) THEN
         y    = dsqrt0(x)
         dy   = 0.5*dquot0(dx,y)
      ELSEIF (i.eq. 6) THEN
         y    = dexp1(x*2.302585)
         dy   = 2.3*y*dx
      ELSEIF (i.eq. 7) THEN
         y    = -x
         dy   = dx
      ELSEIF (i.eq. 8) THEN
         y    = dquot0 (1.d0, x)
         dy   = dx * y**2
      ELSEIF (i.eq. 9) THEN
         y    = dabs (x)
         dy   = dx
      ELSEIF (i.eq.11) THEN
         y    = dsin (x)
         dy   = dabs (dcos (x) * dx)
      ELSEIF (i.eq.12) THEN
         IF (x.lt.-1 .or. x.gt.1) GOTO 90
         y    = dasin (x)
         dy   = 0.
      ELSEIF (i.eq.13) THEN
         y    = dcos (x)
         dy   = dabs (dsin (x) * dx)
      ELSEIF (i.eq.14) THEN
         IF (x.lt.-1 .or. x.gt.1) GOTO 90
         y    = dacos (x)
         dy   = 0.
      ELSEIF (i.eq.16) THEN
         IF (x.lt.-1 .or. x.gt.1) GOTO 90
         y    = datan (x)
         dy   = 0.
      ELSEIF (i.eq.21) THEN
         y    = dsind (x)
         dy   = 0
      ELSEIF (i.eq.22) THEN
         IF (x.lt.-1 .or. x.gt.1) GOTO 90
         y    = dasind (x)
         dy   = 0.
      ELSEIF (i.eq.23) THEN
         y    = dcosd (x)
         dy   = 0
      ELSEIF (i.eq.24) THEN
         IF (x.lt.-1 .or. x.gt.1) GOTO 90
         y    = dacosd (x)
         dy   = 0.
      ELSEIF (i.eq.31) THEN
         y    = dsinh (x)
         dy   = 0
      ELSEIF (i.eq.32) THEN
         y    = dcosh (x)
         dy   = 0
      ELSEIF (i.eq.41) THEN
         y    = dgamma1 (x)
         dy   = 0.
C  Two arguments (i=50,..,99) :
      ELSEIF (i.eq.50) THEN
         y    = z
         dy   = dz
      ELSEIF (i.eq.51) THEN
         y    = x + z
         dy   = dsqrt0 (dx**2 + dz**2)
      ELSEIF (i.eq.52) THEN
         y    = x - z
         dy   = dsqrt0 (dx**2 + dz**2)
      ELSEIF (i.eq.53) THEN
         y    = z - x
         dy   = dsqrt0 (dx**2 + dz**2)
      ELSEIF (i.eq.54) THEN
         y    = x * z
         dy   = dsqrt0 ( z**2 * dx**2 + x**2 * dz**2 )
      ELSEIF (i.eq.55) THEN
         y    = dquot0 (x, z)
         dy   = dquot0 ( dsqrt0 ( z**2 * dx**2 + x**2 * dz**2 ), z**2 )
      ELSEIF (i.eq.56) THEN
         y    = dquot0 (z, x)
         dy   = dquot0 ( dsqrt0 ( z**2 * dx**2 + x**2 * dz**2 ), x**2 )
      ELSEIF (i.eq.57) THEN
         y    = dpow0 (x, z)
         dy   = dsqrt0 ( (z*dpow0(x,z-1)*dx)**2 + (dln0(x)*y*dz)**2 )
      ELSEIF (i.eq.58) THEN ! the error calculation is ugly
         IF     (x.lt.z) THEN
            y    = x
            dy   = dx
         ELSEIF (x.gt.z) THEN
            y    = z
            dy   = dz
         ELSE
            y    = (x+z)/2
            dy   = dsqrt(x**2+z**2)/2
            ENDIF
      ELSEIF (i.eq.59) THEN
         IF     (x.gt.z) THEN
            y    = x
            dy   = dx
         ELSEIF (x.lt.z) THEN
            y    = z
            dy   = dz
         ELSE
            y    = (x+z)/2
            dy   = dsqrt(x**2+z**2)/2
            ENDIF
      ELSEIF (i.eq.60) THEN
         y    = dmod0 (x, z)
         dy   = dx
      ELSE
         Print *, ' function no. i = ', i
         CALL Absturz ('FuVal', 'Call to undefined function')
         ENDIF

      RETURN ! regular end

 90   CONTINUE
      y  = 0
      dy = 0

      END ! FuVal

      BLOCK DATA FuText
C     -----------------
            ! JWu 4sep93

      PARAMETER (LF=99)
      CHARACTER FText*20, FSigl*5
      COMMON   /LibFu/ FText(0:LF), FSigl(0:LF)

      DATA   FText / 'p1',
     *               'ln p1', 'lg p1', 'exp p1', 'p1^2', 'p1^1/2',
     *               '10^p1', '-p1',   '1/p1',  '|p1|', ' ',
     *               'sin p1', 'asin p1', 'cos p1', 'acos p1', ' ',
     *               'atan p1', 4*' ',
     *               'sind p1', 'asind p1', 'cosd p1', 'acosd p1', ' ',
     *               5*' ',
     *               'sinh p1', 'cosh p1', ' ', ' ', ' ',
     *               5*' ',
     *               'Gamma(p1)', ' ', ' ', ' ', ' ',
     *               4*' ',
     *               'p2',
     *               'p1+p2', 'p1-p2', 'p2-p1', 'p1*p2', 'p1/p2',
     *               'p2/p1', 'p1^p2', 'min(p1,p2)', 'max(p1,p2)',
     *               'p1|p2', 39*' ' /

      DATA   FSigl / '=',
     *               'ln', 'lg', 'e^', '^2', '^1/2',
     *               '10^', '0-', '1/', '||', ' ',
     *               'sin', 'asin', 'cos', 'acos', ' ',
     *               'atan', 4*' ',
     *               'sind', 'asind', 'cosd', 'acosd', ' ',
     *               5*' ',
     *               'sinh', 'cosh', ' ', ' ', ' ',
     *               5*' ',
     *               'Gamma', ' ', ' ', ' ', ' ',
     *               4*' ',
     *               '~',
     *               '+', '-', '~-', '*', '/',
     *               '~/', '^', 'min', 'max', 'mod',
     *               39*' ' /

      END ! FuText

      SUBROUTINE FuTxt (i, text, t1, t2)
C     ----------------------------------
            ! separated from FuVal : JWu 29oct92. Renewed 4sep93.
         ! Return description of function no. i
         ! as text('p1','p2') = y(x,z)

      IMPLICIT LOGICAL (q)

      PARAMETER (LF=99)
      CHARACTER text*(*), FText*20, FSigl*5, t1*(*), t2*(*)
      COMMON   /LibFu/ FText(0:LF), FSigl(0:LF)

      IF (qiinside(i,0,LF)) THEN
         IF (FText(i).ne.' ') THEN
            text = FText(i)
            CALL ReplaceT (text, 'p1', t1)
            CALL ReplaceT (text, 'p2', t2)
         ELSE
            text = '&undefined'
            ENDIF
      ELSE
         text = '&undefined'
         ENDIF

      END ! FuTxt

      SUBROUTINE FuInv (i, x, dx, y, dy, z, dz)
C     -----------------------------------------
         ! Calculate x, dx from y, dy, z, dz according to
         ! the inverted function no. i

      IMPLICIT REAL*8 (a-h,o-p,r-z)

      IF (dx.lt.0.) Print *, 'Warning from FuVal : dx<0'
      IF (dz.lt.0.) Print *, 'Warning from FuVal : dz<0'

C  One argument (i=0,..,49) :
      IF     (i.eq. 0) THEN
         x    = y
         dx   = dy
         ! inverse of y = 'p1'
      ELSEIF (i.eq. 1) THEN
         x    = dexp(y)
         dx   = x * dy
         ! inverse of y = 'ln p1'
      ELSEIF (i.eq. 2) THEN
         x    = 10**y
         dx   = 2.3026 * x * dy
         ! inverse of y = 'lg p1'
      ELSEIF (i.eq. 3) THEN
         x    = dln0 (y)
         dx   = dquot0 (dy,y)
         ! inverse of y = 'exp p1'
      ELSEIF (i.eq. 4) THEN
         x    = dsqrt0 (y)
         dx   = dquot0 (dy, 2*x)
         ! inverse of y = 'p1^2'
      ELSEIF (i.eq. 5) THEN
         x    = y**2
         dx   = 2 * y * dy
         ! inverse of y = 'p1^1/2'
      ELSEIF (i.eq. 6) THEN
         x    = dln0(y) / 2.3026
         dx   = dquot0(dy,y) / 2.3026
         ! inverse of y = '10^x'
      ELSEIF (i.eq. 7) THEN
         x    = -y
         dx   = dy
         ! inverse of y = '-x'
      ELSEIF (i.eq. 8) THEN
         x    = dquot0 (1.d0, y)
         dx   = dy * x**2
         ! inverse of y = '1/x'
      ELSEIF (i.eq. 9) THEN
         x    = y
         dx   = dy
         ! inverse of y = '|x|' is unknown
      ELSEIF (i.eq.11) THEN
         x    = dasin (y)
         dx   = 0
      ELSEIF (i.eq.12) THEN
         x    = dsin (x)
         dx   = 0.
      ELSEIF (i.eq.13) THEN
         x    = dacos (y)
         dx   = 0
      ELSEIF (i.eq.14) THEN
         x    = dcos (y)
         dx   = 0.
      ELSEIF (i.eq.16) THEN
         x    = dtan (y)
         dx   = 0
      ELSEIF (i.eq.21) THEN
         x    = dasind (y)
         dx   = 0
      ELSEIF (i.eq.22) THEN
         x    = dsind (x)
         dx   = 0.
      ELSEIF (i.eq.23) THEN
         x    = dacosd (y)
         dx   = 0
      ELSEIF (i.eq.24) THEN
         x    = dcosd (y)
         dx   = 0.
      ELSEIF (i.eq.31) THEN
         x    = y
         dx   = dy
         ! inverse of y = sinh(x) is unknown
      ELSEIF (i.eq.32) THEN
         x    = y
         dx   = dy
         ! inverse of y = sinh(x) is unknown
      ELSEIF (i.eq.41) THEN
         x    = y
         dx   = dy
         ! inverse of y = Gamma(x) is unknown

C  Two arguments (i=50,..,99) :
      ELSEIF (i.eq.50) THEN
         x    = z
         dx   = dz
         ! inverse of y = 'p1^p2'
      ELSEIF (i.eq.51) THEN
         x    = y - z
         dx   = dsqrt0 (dy**2 + dz**2)
         ! inverse of y = 'p1+p2'
      ELSEIF (i.eq.52) THEN
         x    = y + z
         dx   = dsqrt0 (dy**2 + dz**2)
         ! inverse of y = 'p1-p2'
      ELSEIF (i.eq.53) THEN
         x    = z - y
         dx   = dsqrt0 (dy**2 + dz**2)
         ! inverse of y = 'p2-p1'
      ELSEIF (i.eq.54) THEN
         x    = dquot0 (y, z)
         dx   = dquot0 ( dsqrt0 ( z**2 * dy**2 + y**2 * dz**2 ), z**2 )
         ! inverse of y = 'p1*p2'
      ELSEIF (i.eq.55) THEN
         x    = y * z
         dx   = dsqrt0 ( z**2 * dy**2 + y**2 * dz**2 )
         ! inverse of y = 'p1/p2'
      ELSEIF (i.eq.56) THEN
         x    = dquot0 (z, y)
         dx   = dquot0 ( dsqrt ( z**2 * dy**2 + y**2 * dz**2 ), y**2 )
         ! inverse of y = 'p2/p1'
      ELSEIF (i.eq.57) THEN
         x    = dpow0 (dquot0(1.d0,z), y)
         dx   = 0. ! keine Lust
         ! inverse of y = 'p1^p2'
      ELSEIF (i.eq.58) THEN
         x    = 0.
         dx   = 0.
         Print *, 'WARNING/ Inverse of function no. 58 undefined'
      ELSEIF (i.eq.59) THEN
         x    = 0.
         dx   = 0.
         Print *, 'WARNING/ Inverse of function no. 59 undefined'
      ELSEIF (i.eq.60) THEN
         x    = 0.
         dx   = 0.
         Print *, 'WARNING/ Inverse of function no. 60 undefined'
      ELSE
         CALL Absturz ('FuInv', 'iFu o.o.r.')
         ENDIF

      END ! FuInv

      SUBROUTINE FuHelp (t1, t2)
C     --------------------------
         ! list of available functions

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      PARAMETER (LF=99)
      CHARACTER FText*20, FSigl*5, t1*(*), t2*(*), text*40
      COMMON   /LibFu/ FText(0:LF), FSigl(0:LF)

      Print *, ' Functions of one argument :'
      DO i = 0, 49
         CALL FuTxt (i, text, t1, t2)
         IF (text.ne.'&undefined') Print '(a,i2,a,a5,2a)',
     *        '  ',i,'  (', FSigl(i), ')   lhs = ', text
         ENDDO

      Print *, ' Binary Operations :'
      DO i = 50, 99
         CALL FuTxt (i, text, t1, t2)
         IF (text.ne.'&undefined') Print '(a,i2,a,a5,2a)',
     *        '  ',i,'  (', FSigl(i), ')   lhs = ', text
         ENDDO
      Print *

      END ! FuHelp

      SUBROUTINE FuAsk (que, inp, ifu, r2a)
C     -------------------------------------
            ! JWu 21jul93, 4sep93
         ! Ask for a general function.
         ! If (inp<>' ') then decode only the answer.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      PARAMETER (LF=99)
      CHARACTER FText*20, FSigl*5
      CHARACTER  que*(*), inp*(*), ein*80, hilf*20

      COMMON   /LibFu/ FText(0:LF), FSigl(0:LF)

      ein = inp

 1    CONTINUE

      IF (ein.eq.' ') THEN
         IF (qiinside(ifu,0,LF)) THEN
            CALL FrageHD (que, ein, FSigl(ifu))
         ELSE
            CALL FrageH  (que, ein)
            ENDIF
         ENDIF

      IF     (ein.eq.' ') THEN
         GOTO 9 ! not allowed
      ELSEIF (ein(1:1).eq.'?' .or. ein(1:1).eq.'h') THEN
         Print *, 'INPUT HELP/'
         Print *, 
     * '   required input : the number or the symbol of a function'
         Print *, '   here is a list of all allowed functions :'
         CALL FuHelp ('A', 'B')
         ein = ' '
         GOTO 1
         ENDIF

      hilf = ein
      CALL Fi1I (hilf, ifu)
      IF (hilf.eq.'#') THEN ! integer given
         IF (FText(ifu).ne.' ') GOTO 2 ! everything ok.
         GOTO 9 ! invalid answer
         ENDIF
      ! text given :
      DO i = 0, LF
         IF (ein.eq.FSigl(i)) THEN
            ifu = i
            GOTO 2 ! Sigl is valid
            ENDIF
         ENDDO
      GOTO 9 ! Sigl not reckognized

C  Part 2 : Determine 2nd argument :
 2    CONTINUE
      IF (ifu.ge.50) THEN
         r2a = rAskD (' 2nd argument', r2a)
         ENDIF
      RETURN

 9    CONTINUE
      CALL Gong (2)
      ein = ' '
      GOTO 1

      END ! FuAsk

C  ====================================================================
C     6.7.  Functions for solid state physics
C  ====================================================================

      REAL*8 FUNCTION u2Debye (T, TD, rMass, eps)
C     -------------------------------------------
            ! JWu 26sep91 for SQW, 22oct91 for WuLib.
         ! Calculate <u_x^2> from the temperature T[K], the Debye
         ! temperature TD[K] and the atomic mass rMass[a.m.u.].
         ! Required precision is eps.

      IMPLICIT REAL*8 (a-h,o-p,r-z)

      IF (rMass.le.0. .or. TD.le.0. .or. eps.le.0. .or. eps.ge.1.) THEN
         u2Debye = -.777777 ! warning
         RETURN
         ENDIF

      IF (T.lt.0.) THEN
         f = 0.  ! no DWF correction wanted
      ELSEIF (T.lt.1.d-3*TD) THEN
         f = .25 ! + exponentially vanishing contribution
      ELSE
         xD = TD / T
         ifail = -1 ! soft exit
c         f = D01AHF (0., xD, 1.d-5, nPts, relErr, DebyeInt, 0, ifail)
         ! Primitiv-Integration :
         M  = idint (1/ eps)
         f  = 0.
         dx = xD / (M+1)
         DO i = 1, M
            x = (i-.5) * dx
            f = f + DebyeInt(x)*dx
            ENDDO
         ! end integration
         f = f / ( 2 * xD**2 )
         ENDIF
      u2Debye = 145.532 * f / ( rMass * TD )  ! 145.5 = 3 hbar^2 / amu / k_B
      END ! u2Debye

      REAL*8 FUNCTION DebyeInt (x)
C     ----------------------------
         ! integrand for u2Debye
      REAL*8 x

      IF (x.lt.1.d-8) THEN
         DebyeInt = 2.
      ELSE
         DebyeInt = x / dtanh(x/2)
         ENDIF

      END ! DebyeInt
