C  ====================================================================
C
C      Library  IDA   :  Inelastic data treatment
C      Modul    i32   :     general functions and symbolic calculation
C
C  ====================================================================

C     Contents :
C        1.  General function IdaFu :
C               IdaFuVal, IdaFuTxt, IdaFuCoord, IdaFuAsk, IdaFuNArg
C        2.  Symbolic calculations :
C               MultUnits, DepowUnit, EnpowUnit

C  ====================================================================
C  i32 / 1 :   General function IdaFu
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE IdaFuVal (i, y, dy, x, dx, z, dz)
C     --------------------------------------------
         ! Calculate y, dy from x, dx, z, dz according to function no. i
         ! For i= 0,..,79, the function is given by FuVal (in Wul6),
         ! for i=80,..,99, the function is defined here.

      IMPLICIT REAL*8 (a-h,o-p,r-z)

      IF     (i.lt.0) THEN
         CALL Absturz ('IdaFuVal', 'function no. < 0')
      ELSEIF (i.le.79) THEN
         CALL FuVal (i, y, dy, x, dx, z, dz)
      ELSEIF (i.eq.81) THEN
         y    = yQ_of_0 (z, x)
         dy   = 0.
      ELSEIF (i.eq.82) THEN
         y    = Tau_of_W (x, z)
         dy   = 0.
      ELSEIF (i.eq.83) THEN
         y    = W_of_Tau (x, z)
         dy   = 0.
      ELSEIF (i.eq.84) THEN
         y    = DKi_of_W (x, z)
         dy   = 0.
      ELSEIF (i.eq.85) THEN
         y    = Th2_of_Q0 (z, x)
         dy   = 0.
      ELSEIF (i.eq.86) THEN
         CALL GoetzeExp (x, y, dum) ! lambda -> a
         dy   = 0.
      ELSEIF (i.eq.87) THEN
         CALL GoetzeExp (x, dum, y) ! lambda -> b
         dy   = 0.
      ELSE
         CALL Absturz ('IdaFuVal', 'function no. > max')
         ENDIF

      END ! IdaFuVal

      SUBROUTINE IdaFuTxt (i, text, t1, t2)
C     -------------------------------------
         ! Return the text describing function no. i,
         ! 'p1', 'p2' will be replaced by t1, t2.
         ! For i= 0,..,79, the function is given by FuVal (in Wul6),
         ! for i=80,..,99, the function is defined here.

      CHARACTER*(*)    text, t1, t2

      IF     (i.lt.0) THEN
         CALL Absturz ('IdaFuTxt', 'function no. < 0')
      ELSEIF (i.le.79) THEN
         CALL FuTxt (i, text, t1, t2)
         RETURN
      ELSEIF (i.eq.81) THEN
         text = 'q(2th=p1;E0=p2)'
      ELSEIF (i.eq.82) THEN
         text = 'tof(w=p1;E0=p2)'
      ELSEIF (i.eq.83) THEN
         text = 'w(tof=p1;E0=p2)'
      ELSEIF (i.eq.84) THEN
         text = 'Dki(w=p1;E0=p2)'
      ELSEIF (i.eq.85) THEN
         text = '2th(q=p1;E0=p2)'
      ELSEIF (i.eq.86) THEN
         text = 'a_mct(lambda=p1)'
      ELSEIF (i.eq.87) THEN
         text = 'b_mct(lambda=p1)'
      ELSE
         text = '&undefined'
         ENDIF
      CALL ReplaceT (text, 'p1', t1)
      CALL ReplaceT (text, 'p2', t2)

      END ! IdaFuTxt

      SUBROUTINE IdaFuCoord (i, LHS, co, un,
     *                       coIn1, unIn1, coIn2, unIn2, Fehler)
C     ----------------------------------------------------------
            ! JWu 1990/91. Rewritten (simplified) 30jun92.
         ! For function no. i, find new coordinate name and unit.
         ! Fi=unctions 101-120 are for OprIntegral.
         ! If necessary, the user is consulted.

      IMPLICIT REAL*8  (a-h,o-p,r-z) ! needed for dummies
      IMPLICIT LOGICAL (q)
      CHARACTER*(*)    LHS, co, coIn1, coIn2, un, unIn1, unIn2, Fehler
      CHARACTER*40     coE, unE, coG, unG, h1, h2
      DATA             iOptN /3/, iOptU /3/, coE /' '/, unE /' '/

      IF     (i.eq. 0) THEN
         coG = coIn1
         unG = unIn1
      ELSEIF (i.eq. 1) THEN
         coG = 'ln '//coIn1
         unG = ' '
      ELSEIF (i.eq. 2) THEN
         coG = 'lg '//coIn1
         unG = ' '
      ELSEIF (i.eq. 3) THEN
         coG = 'exp '//coIn1
         unG = ' '
      ELSEIF (i.eq. 4) THEN
         CALL Compose2 (coG, coIn1, '^2')
         IF (unIn1.ne.' ') CALL Compose2 (unG, unIn1, '^2')
      ELSEIF (i.eq. 5) THEN
         CALL Compose2 (coG, coIn1, '^1/2')
         IF (unIn1.ne.' ') CALL Compose2 (unG, unIn1, '^1/2')
      ELSEIF (i.eq. 6) THEN
         coG = '10^'//coIn1
         unG = ' '
      ELSEIF (i.eq. 7) THEN
         coG = '-'//coIn1
         unG = unIn1
      ELSEIF (i.eq. 8) THEN
         coG = '1/'//coIn1
         IF (unIn1.ne.' ') unG = '1/'//unIn1
      ELSEIF (i.eq. 9) THEN
         CALL Compose2 (coG, '|'//coIn1, '|')
         unG = unIn1
      ELSEIF (i.eq.11) THEN
         coG = 'sin '//coIn1
         unG = ' '
      ELSEIF (i.eq.12) THEN
         coG = 'asin '//coIn1
         unG = ' '
      ELSEIF (i.eq.13) THEN
         coG = 'cos '//coIn1
         unG = ' '
      ELSEIF (i.eq.14) THEN
         coG = 'acos '//coIn1
         unG = ' '
      ELSEIF (i.eq.16) THEN
         coG = 'acos '//coIn1
         unG = ' '
      ELSEIF (i.eq.21) THEN
         coG = 'sind '//coIn1
         unG = ' '
      ELSEIF (i.eq.22) THEN
         coG = 'asind '//coIn1
         unG = ' '
      ELSEIF (i.eq.23) THEN
         coG = 'cosd '//coIn1
         unG = ' '
      ELSEIF (i.eq.24) THEN
         coG = 'acosd '//coIn1
         unG = ' '
      ELSEIF (i.eq.31) THEN
         coG = 'sinh '//coIn1
         unG = ' '
      ELSEIF (i.eq.32) THEN
         coG = 'cosh '//coIn1
         unG = ' '
      ELSEIF (i.eq.41) THEN
         coG = 'Gamma '//coIn1
         unG = ' '
      ELSEIF (i.eq.50) THEN
         coG = coIn2
         unG = unIn2
      ELSEIF (i.eq.51) THEN
         IF     (coIn1.eq.coIn2) THEN
            coG = coIn1
         ELSEIF (coIn2.eq.' ') THEN
            coG = coIn1
         ELSE
            CALL Compose3 (coG, coIn1, '+', coIn2)
            ENDIF
         unG = unIn1
      ELSEIF (i.eq.52) THEN
         IF     (coIn1.eq.coIn2) THEN
            coG = coIn1
         ELSEIF (coIn2.eq.' ') THEN
            coG = coIn1
         ELSE
            CALL Compose3 (coG, coIn1, '-', coIn2)
            ENDIF
         unG = unIn1
      ELSEIF (i.eq.53) THEN
         IF (coIn1.eq.coIn2) THEN
            coG = coIn1
         ELSE
            CALL Compose3 (coG, coIn2, '-', coIn1)
            ENDIF
         unG = unIn1
      ELSEIF (i.eq.54) THEN
         IF     (coIn2.eq.' ') THEN
            coG = coIn1
         ELSE
            CALL Compose3 (coG, coIn1, '*', coIn2)
            ENDIF
         CALL Compose2 (h1, unIn1, '-1')
         CALL Compose2 (h2, unIn2, '-1')
         IF     (h1.eq.unIn2 .or. h2.eq.unIn1) THEN
            unG = ' '
         ELSEIF (unIn2.eq.' ') THEN
            unG = unIn1
         ELSEIF (unIn1.eq.' ') THEN
            unG = unIn2
         ELSE
            CALL Compose3 (unG, unIn1, '*', unIn2)
            ENDIF

      ELSEIF (i.eq.55) THEN
         IF     (coIn2.eq.' ') THEN
            coG = coIn1
         ELSEIF (coIn1.eq.' ') THEN
            coG = '1/'//coIn2
         ELSEIF (coIn1.eq.coIn2) THEN
            coG = coIn1//'/idem'
         ELSE
            CALL Compose3 (coG, coIn1, '/', coIn2)
            ENDIF
         IF     (unIn2.eq.' ') THEN
            unG = unIn1
         ELSEIF (unIn1.eq.' ') THEN
            unG = '1/'//unIn2
         ELSEIF (unIn1.eq.unIn2) THEN
            unG = ' '
         ELSE
            CALL Compose3 (unG, unIn1, '/', unIn2)
            ENDIF
      ELSEIF (i.eq.56) THEN
         IF     (coIn2.eq.' ') THEN
            coG = '1/'//coIn1
         ELSEIF (coIn1.eq.' ') THEN
            coG = coIn2
         ELSEIF (coIn1.eq.coIn2) THEN
            coG = ' '
         ELSE
            CALL Compose3 (coG, coIn2, '/', coIn1)
            ENDIF
         IF     (unIn2.eq.' ') THEN
            unG = '1/'//unIn1
         ELSEIF (unIn1.eq.' ') THEN
            unG = unIn2
         ELSEIF (unIn1.eq.unIn2) THEN
            unG = ' '
         ELSE
            CALL Compose3 (unG, unIn2, '/', unIn1)
            ENDIF
      ELSEIF (i.eq.57) THEN
         CALL Compose3 (coG, coIn1, '^', coIn2)
         CALL Compose3 (unG, unIn1, '^', coIn2)
      ELSEIF (i.eq.58) THEN
         IF (coIn1.eq.coIn2) THEN
            coG = coIn1
         ELSE
            CALL Compose4 (coG, 'min_', coIn1, ',', coIn2)
            ENDIF
         unG = unIn1
      ELSEIF (i.eq.59) THEN
         IF (coIn1.eq.coIn2) THEN
            coG = coIn1
         ELSE
            CALL Compose4 (coG, 'max_', coIn1, ',', coIn2)
            ENDIF
         unG = unIn1
      ELSEIF (i.eq.60) THEN
         CALL Compose3 (coG, coIn1, ' mod ', coIn2)
         CALL Compose3 (unG, unIn1, '/', unIn2)
      ELSEIF (i.eq.81) THEN
         coG = 'q'
         unG = 'A-1'
      ELSEIF (i.eq.82) THEN
         coG = 'tof'
         unG = 'msec/m'
      ELSEIF (i.eq.83) THEN
         coG = 'w'
         unG = 'meV'
      ELSEIF (i.eq.84) THEN
         coG = 'Dki'
         unG = 'A-1'
      ELSEIF (i.eq.85) THEN
         coG = '2th'
         unG = ' '
      ELSEIF (i.eq.86) THEN
         coG = 'a_mct'
         unG = ' '
      ELSEIF (i.eq.87) THEN
         coG = 'b_mct'
         unG = ' '
      ELSE
         Fehler = ' idaFuCoord/ i o.o.r.'
         RETURN
         ENDIF

      co = coG
      un = unG

      END ! IdaFuCoord

      INTEGER FUNCTION IdaFuAsk (que, inp, t1, t2, ifd, narmax)
C     ---------------------------------------------------------
            ! JWu 4/6sep93
         ! Get a function number, default is ifd, max # arguments narmax.
         ! Result must be tested for -1 (escape)

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)
      PARAMETER (LF=99)
      CHARACTER que*(*), inp*(*), t1*(*), t2*(*), text*40, 
     *          ein*40, hilf*40
      CHARACTER FText*20, FSigl*5, IFSigl(0:LF)*5

      COMMON   /LibFu/ FText(0:LF), FSigl(0:LF) ! initialized in L6

      SAVE     q1stCall, IFSigl

      DATA     q1stCall / .true. /

C  Initializations :
      IF (q1stCall) THEN
         DO i=0,LF
            IFSigl(i) = FSigl(i)
            ENDDO
         IFSigl(81) = 'q'
         IFSigl(82) = 'tof_w'
         IFSigl(83) = 'w_tof'
         IFSigl(84) = 'Dki'
         IFSigl(85) = '2th'
         IFSigl(86) = 'a_mct'
         IFSigl(87) = 'b_mct'
         q1stCall = .false.
         ENDIF

      ein = inp
      ifu = ifd ! default

C  Loop : ask until answer is valid :
 1    CONTINUE

      IF (ein.eq.' ') THEN
         IF (qiinside(ifu,0,LF)) THEN
            CALL FrageHD (que, ein, IFSigl(ifu))
         ELSE
            CALL FrageH  (que, ein)
            ENDIF
         ENDIF

      IF     (ein.eq.' ') THEN
         GOTO 9 ! not allowed
      ELSEIF (ein.eq.'^]' .or. ein.eq.CHAR(27)) THEN
         IdaFuAsk = -1
         RETURN ! escape
      ELSEIF (ein(1:1).eq.'?' .or. ein(1:1).eq.'h') THEN
         Print *, 'INPUT HELP/'
         Print *, 
     * '   required input : the number or the symbol of a function'
         Print *, '   use ''^]'' to escape'
         Print *, '   here is a list of all allowed functions :'
         CALL FuHelp (t1, t2)
         Print *, ' Ida-defined functions :'
         DO i = 80, 99
            CALL IdaFuTxt (i, text, t1, t2)
            IF (text.ne.'&undefined') Print '(a,i2,a,a5,2a)',
     *        '  ',i,'(', IFSigl(i), ') :  lhs = ', text
            ENDDO
         Print *
         ein = ' '
         GOTO 1
         ENDIF

      hilf = ein
      CALL Fi1I (hilf, ifu)
      IF (hilf.eq.'#') THEN ! integer given
         IF (IFSigl(ifu).ne.' ') GOTO 2 ! everything ok.
         GOTO 9 ! invalid answer
         ENDIF
      ! text given :
      DO i = 0, LF
         IF (ein.eq.IFSigl(i)) THEN
            ifu = i
            GOTO 2 ! Sigl is valid
            ENDIF
         ENDDO
      GOTO 9 ! Sigl not reckognized

C  Valid answer :
 2    CONTINUE
      IdaFuAsk = ifu
      RETURN

 9    CONTINUE
      CALL Gong (2)
      ein = ' '
      GOTO 1

      END ! IdaFuAsk

      INTEGER FUNCTION IdaFuNArg (ifu)
C     --------------------------------
            ! JWu 4sep93
         ! Returns number of arguments of function ifu.

      CHARACTER*40 text

      CALL IdaFuTxt (ifu, text, 'p1', 'p2')

      IF (text.eq.'&undefined') THEN
         IdaFuNArg = 0
      ELSEIF (ifu.lt.50) THEN
         IdaFuNArg = 1
      ELSE
         IdaFuNArg = 2
         ENDIF

      END ! IdaFuNArg

C  ====================================================================
C  i32 / 2 :   Symbolic Calculation
C  ====================================================================

      SUBROUTINE MultUnits (out, in1, jp1, in2, jp2)
C     ----------------------------------------------
            ! JWu 19may95
         ! propose out = in1^jp1 * in2^jp2

      CHARACTER*(*)   out, in1, in2
      CHARACTER*40    au1, au2

      au1 = in1
      au2 = in2

      CALL DepowUnit (au1, ja1)
      CALL DepowUnit (au2, ja2)

      ja1 = ja1*jp1
      ja2 = ja2*jp2

      IF (au1.eq.au2) THEN
         CALL EnpowUnit (au1, ja1+ja2)
         out = au1
      ! ELSE Mult..
      ELSEIF (ja1.eq.0) THEN
         CALL EnpowUnit (au2, ja2)
         out = au2
      ELSEIF (ja2.eq.0) THEN
         CALL EnpowUnit (au1, ja1)
         out = au1
      ELSE
         CALL EnpowUnit (au1, ja1)
         CALL EnpowUnit (au2, ja2)
         CALL Compose2 (out, au1, '*'//au2)
         ENDIF

      END ! MultUnits

      SUBROUTINE DepowUnit (io, jp)
C     -----------------------------
            ! JWu 19may95
         ! decode io as io^jp

      CHARACTER io*(*), aux*40
      INTEGER   IAux(7)

      jp = 1 ! default, for case of early return

      aux = io
      CALL FindN (aux, 7, niaux, IAux)
      IF (niaux.le.0) RETURN

      jaux = lenU(aux)
      IF (jaux.le.2) RETURN
      IF (aux(jaux:jaux).ne.'#') RETURN
      IF     (aux(jaux-1:jaux-1).eq.'^' .or. 
     *        aux(jaux-1:jaux-1).eq.'+') THEN
         jp =  IAux(niaux)
         io = aux(1:jaux-2)
      ELSEIF (aux(jaux-1:jaux-1).eq.'-') THEN
         jp = -IAux(niaux)
         io = aux(1:jaux-2)
         ENDIF
      IF (jaux.le.3) RETURN
      IF (aux(jaux-2:jaux-1).eq.'^-') THEN
         jp = -IAux(niaux)
         io = aux(1:jaux-3)
         ENDIF

      END ! DepowUnit

      SUBROUTINE EnpowUnit (io, jp)
C     -----------------------------
            ! JWu 19may95 immer noch vorm Fr"uhst"uck
         ! replace io by io^jp

      CHARACTER io*(*), cl3*3

      IF     (jp.eq.0) THEN
         io = ' '
      ELSEIF (jp.eq.1) THEN
         ! io = io
      ELSEIF (jp.gt.1) THEN
         CALL Append (io, '^'//cl3(jp))
      ELSEIF (jp.lt.0) THEN
         CALL Append (io, '-'//cl3(-jp))
         ENDIF

      END ! DepowUnit
