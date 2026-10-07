C  ====================================================================
C
C     Library  WuLib :  General FORTRAN Library
C     Module   WuL1  :  Elementary Operations
C
C  ====================================================================

C     J.Wuttke - 1989ff.

C     Outline :
C        This module is used by all my programs. It contains
C        type conversion functions, elementary calculations,
C        and the worst case subroutine Absturz1.
C        No other error treatment than Absturz1 is foreseen,
C        however, most of the subroutines return a defined value
C        if the input is out of range.

C     Contents :
C
C        1.1.   Type conversion :
C        1.1.1.    int  -> char :
C                     ch1, cr2..cr8, cl2..cl8, cv2..cv8,
C                        (auxiliary : convICr/l/v)
C        1.1.2.    real -> char :
C                     convRCs (ausser Betrieb), NiceNum, Enc..
C        1.1.4.    char -> int :
C                     ichar1
C        1.1.5.    bool -> int  :
C                     intq
C        1.1.6.    int  -> bool :
C                     qintr
C        1.1.7.   date & time -> string
C                     Datum, Zeit
C
C        1.2.   Integer functions :
C        1.2.1.    int -> int :
C                     iInside
C        1.2.2.    real-> int :
C                     iPint
C        1.2.3.    on the decimal representation :
C                     jdigit, ndigits
C        1.2.4.    on the binary representation :
C                     SetBit, iGetBit
C
C        1.3.   Real functions :
C        1.3.1.    returning 0 on error :
C                     dpow0, dpowii, dquot0, dln0, dlg0, dsqrt0, dmod0, dLorentz
C        1.3.2.    returning other well defined values :
C                     dquot1, dexp1, dexp2
C        1.3.3.    of complex argument (Re,Im) :
C                     dphase
C        1.3.4.    like min, max :
C                     dinside
C        1.3.5.    rounding :
C                     rund, RoundLin, RoundLog
C        1.3.6.    interpolation, integration :
C                     LinIntPol, -Array, rIntegralXY
C        1.3.7.    with error propagation :
C                     DEmult, DEquot0
C        1.3.8.    mathematical and physical constants :
C                     R8Const
C
C
C        1.4.   Idem for REAL*4 .
C                     s..
C
C        1.5.   Logical functions :
C                     qEq(Eps|Tol), q(i|r)(Nin|In)(Open|Clos|Range)
C
C        1.7.   Operations on vectors :
C        1.7.1.    Integer :
C                     iSum, iMinMax
C        1.7.2.    Real*8 :
C                     rCopy, rSet, rSum, Sum2, rMinMax, rStepAt,
C                     rCycling, irPos(|Opt), rSort, irSorted, CheckScale
C        1.7.4.    Logical :
C                     iqSum, qCopy, qSet
C
C        1.9.  Control :
C                     Absturz1 (wo, warum)

C     History :
C        JWu 1989   Module A_
C        JWu oct90  VMS version WuL1
C        JWu 1990   vector operations
C        JWu jul91  compact lay-out

C  ====================================================================
C        1.1.  Typumwandlung :
C  ====================================================================

C  --------------------------------------------------------------------
C        1.1.1.   int  -> char : ch, cr, cl, cv, (convICl), (convICv)
C  --------------------------------------------------------------------
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      CHARACTER*1 FUNCTION ch1 (i)
C     ----------------------------
         IF (0.le.i .and. i.le.9) THEN
            ch1 = char(48+i)
         ELSE
            ch1 = '*'
            ENDIF
         END

C     Rechtsbuendig :
C     ---------------

      CHARACTER*2 FUNCTION cr2 (i)
         CHARACTER*2 stri
         CALL convICr (i, stri)
         cr2 = stri
         END
      CHARACTER*3 FUNCTION cr3 (i)
         CHARACTER*3 stri
         CALL convICr (i, stri)
         cr3 = stri
         END
      CHARACTER*4 FUNCTION cr4 (i)
         CHARACTER*4 stri
         CALL convICr (i, stri)
         cr4 = stri
         END
      CHARACTER*5 FUNCTION cr5 (i)
         CHARACTER*5 stri
         CALL convICr (i, stri)
         cr5 = stri
         END
      CHARACTER*6 FUNCTION cr6 (i)
         CHARACTER*6 stri
         CALL convICr (i, stri)
         cr6 = stri
         END
      CHARACTER*7 FUNCTION cr7 (i)
         CHARACTER*7 stri
         CALL convICr (i, stri)
         cr7 = stri
         END
      CHARACTER*8 FUNCTION cr8 (i)
         CHARACTER*8 stri
         CALL convICr (i, stri)
         cr8 = stri
         END


C     Linksbuendig :
C     --------------

      CHARACTER*2 FUNCTION cl2 (i)
         CALL convICl (i, cl2)
         END
      CHARACTER*3 FUNCTION cl3 (i)
         CALL convICl (i, cl3)
         END
      CHARACTER*4 FUNCTION cl4 (i)
         CALL convICl (i, cl4)
         END
      CHARACTER*5 FUNCTION cl5 (i)
         CALL convICl (i, cl5)
         END
      CHARACTER*6 FUNCTION cl6 (i)
         CALL convICl (i, cl6)
         END
      CHARACTER*7 FUNCTION cl7 (i)
         CALL convICl (i, cl7)
         END
      CHARACTER*8 FUNCTION cl8 (i)
         CALL convICl (i, cl8)
         END


C     Mit Nullen gefuellt :
C     ---------------------

      CHARACTER*2 FUNCTION cv2 (i)
         CALL convICv (i, cv2)
         END
      CHARACTER*3 FUNCTION cv3 (i)
         CALL convICv (i, cv3)
         END
      CHARACTER*4 FUNCTION cv4 (i)
         CALL convICv (i, cv4)
         END
      CHARACTER*5 FUNCTION cv5 (i)
         CALL convICv (i, cv5)
         END
      CHARACTER*6 FUNCTION cv6 (i)
         CALL convICv (i, cv6)
         END
      CHARACTER*7 FUNCTION cv7 (i)
         CALL convICv (i, cv7)
         END
      CHARACTER*8 FUNCTION cv8 (i)
         CALL convICv (i, cv8)
         END


C     auxiliary routines convICr/l/v :
C     --------------------------------

      SUBROUTINE convICr (iin, stri)
C     ------------------------------
            ! JWu 10mar93 to prevent internal write
         ! Schreibt i rechtsbuendig nach stri

      CHARACTER  stri*(*), ch1*1

      ls  = len(stri)
      stri = ' '
      IF (iin.eq.0) THEN
         stri(ls:ls) = '0'
         RETURN
         ENDIF
      la = iabs (ndigits(iin))
      IF (la.gt.ls .or. (iin.lt.0.and.la+1.gt.ls)) THEN ! overflow
         DO  j = 1, ls
            stri (j:j) = '*'
            ENDDO
         RETURN
         ENDIF
      IF (iin.lt.0) stri(ls-la:ls-la) = '-'
      ia = abs(iin)
      DO j = 1, la
         stri(ls+1-j:ls+1-j) = ch1( jdigit(ia,j) )
         ENDDO

      END ! convICr

      SUBROUTINE convICl ( i, c)
C     --------------------------
            ! rev 16.10.90
         ! Schreibt i linksbuendig nach c

      CHARACTER  c *(*), ch1*1

      ls  = len (c)
      IF (i.eq.0) THEN
         c = '0'
            RETURN
         ENDIF
      la = iabs (ndigits(i))
      IF (i.gt.0) THEN
         lf = la            ! last digit
         l0 = 1             ! first digit
         c  = ' '
      ELSE  !(i.lt.0)
         lf = la + 1
         l0 = 2
         c  = '-'
         ENDIF
      IF (lf.gt.ls) THEN
         DO  8 j = 1, ls
            c (j:j) = '*'
  8         CONTINUE
         RETURN
         ENDIF
      ia   = abs (i)
      DO 10 j = l0,lf
         c (j:j) = ch1 ( jdigit (ia, la+l0-j) )
 10      CONTINUE

      END ! convICl

      SUBROUTINE convICv ( i, c)
C     --------------------------
C        schreibt i nach c, rechtsbuendig und mit Nullen aufgefuellt

         CHARACTER c *(*), form *24, ch1 *1
         l = len (c)
         form = '(i' // ch1(l) // ')'
         write (c,form,err=91) i
         DO 1 j = 1,l
            IF (c(j:j).eq.' ') THEN
               c(j:j) = '0'
            ELSEIF (c(j:j).eq.'-') THEN
                        c(j:j) = '0'
                        c(1:1) = '-'
            ELSE
               RETURN
               ENDIF
 1          CONTINUE
         CALL Absturz1 ('convICv', 'unexpected exit from DO loop')
 91      DO 911 j = 1,l
            c(j:j) = '*'
 911        CONTINUE
         END ! convICv

C  --------------------------------------------------------------------
C        1.1.2.   real -> char : convRCs
C  --------------------------------------------------------------------

      SUBROUTINE NiceNum (X, Stri, icNum)
C     -----------------------------------
            ! H.P.Schildberg =< 1987.
            ! Revisions JWu 14apr91, 13feb92, 13oct92.
                  ! WARNUNG : benoetigt z.Z. WuL2
      ! Value X is written into the string Stri which must
      ! have at least mNum elements. The number of non-blank
      ! characters of Stri will be given by icNum.

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER       (mNum=13)
      CHARACTER*(*)    Stri       ! as such 13feb92
      CHARACTER        c*1

      icNum = 0
      lNum  = len(Stri)
      IF (lNum.lt.mNum) CALL Absturz1 ('NiceNum', 'String too short')

      x10 = dlg0(dabs(X))
      IF (-8.99.lt.x10 .and. x10.lt.8.99) THEN
         write (Stri, '(1pg12.5e1)') X
      ELSE
         write (Stri, '(1pg12.5e2)') X
         ENDIF

cdeb      print *, 'NiceNum / X = ', X
cdeb      print *, 'NiceNum / S = ', Stri

      IF (Stri(2:3).eq.' .') Stri(2:2) = '0' ! f"ur Andreas, f"ur HP
      IF (Stri(3:4).eq.' .') Stri(3:3) = '0'

C  Suppress leading blanks :
      i = 0
10    CONTINUE
         i = i + 1
         IF (Stri(i:i).ne.' ') GOTO 19  ! first non-blank character
         IF (i.eq.mNum) GOTO 99         ! only blanks in String
         GOTO 10
19    CONTINUE
      i = i-1                           ! number of leading blanks
      IF (i.gt.0) THEN                  ! suppress them
         CALL DelVonBis (Stri, 1, i)
         ENDIF
      nNum = mNum - i                   ! remaining characters

C  Position of 'E' and '.' ?
      iE  = nNum
      iDP = nNum
      DO i = 1,nNum
         IF ((Stri(i:i).eq.'E') .or. (Stri(i:i).eq.'e')) THEN
            iE = i
            Stri(i:i) = 'e'             ! replace 'E' by 'e'
            ENDIF
         IF (Stri(i:i).eq.'.') iDP = i
         ENDDO

C  Delete terminating zeros :
      iTZ = iE
60    CONTINUE
         iTZ = iTZ - 1
         IF (iTZ.lt.1) THEN
            CALL Absturz1 ('NiceNum', ' Endless loop 60')
            ENDIF
         c = Stri(iTZ:iTZ)
         IF (c.eq.'0') GOTO 60
         IF (c.eq.' ') GOTO 60
      ! now Stri(iTZ:iTZ) contains DP or significant digit.
      IF (Stri(iTZ:iTZ).eq.'.') THEN
         iTZ = iTZ - 1     ! delete the DP
c         iTZ = iTZ + 1     ! leave one '0' after DP
c         IF (Stri(iTZ:iTZ).eq.' ') Stri(iTZ:iTZ)='0'
c            ! If necessary, insert this '0'. Don't care about
c            ! shifting the exponent part : the case of an
c            ! exponent following a DP will never happen.
         ENDIF
      IF (iE.eq.nNum) THEN
         icNum = iTZ
         IF (icNum.lt.lNum) Stri(icNum+1:lNum) = ' '
         GOTO 99                   ! Finished, if no exponent
         ENDIF
      ! Shift the exponent :
      IF (iTZ.lt.iE-1) THEN
         CALL DelVonBis (Stri, iTZ+1, iE-1) ! flushleft
         iE = iTZ + 1
         ENDIF

C  Now make the exponent nicer :
      icNum = iE + 3
      IF     (icNum.gt.nNum) THEN
         icNum = nNum
      ELSEIF (Stri(icNum:icNum).eq.' ') THEN
         icNum = icNum - 1
         ENDIF

      IF     (icNum.eq.iE+3) THEN
         IF (Stri(iE+2:iE+3).eq.'00')  icNum = iE - 1  ! discard 'E+00'
      ELSEIF (icNum.eq.iE+2) THEN
         IF (Stri(iE+2:iE+2).eq.'0')   icNum = iE - 1  ! discard 'E+0'
         ENDIF
      IF (icNum.ge.iE+2) THEN
         IF (Stri(iE+1:iE+1).eq.'+') THEN
            CALL DelVonBis (Stri, iE+1, iE+1)          ! discard '+'
            icNum = icNum - 1
            ENDIF
         ENDIF
      IF (icNum.ge.iE+2) THEN
         IF (Stri(iE+1:iE+1).eq.'0') THEN
            CALL DelVonBis (Stri, iE+1, iE+1)          ! discard '0'
            icNum = icNum - 1
            ENDIF
         ENDIF

      DO i = icNum + 1, lNum
         Stri(i:i) = ' ' ! 13feb92
         ENDDO

 99   CONTINUE

      END ! NiceNum

      SUBROUTINE EncC (aus, ia, ni, stri)
C     -----------------------------------
            ! 15apr94
         ! put stri flushright into aus(ia:..)
      CHARACTER     aus*(*), stri*(*)

      ls = min0 (ni, len(stri))
      aus (ia+ni-ls:ia+ni-1) = stri
      ia = ia + ni

      END ! EncC

      SUBROUTINE EncI (aus, ia, ni, ival)
C     -----------------------------------
            ! 14apr94
         ! encode ival into aus(ia:..)
      CHARACTER     aus*(*), form*12, stri*12, cv2*2

      stri = ' '
      ivmax =   10**(ni)
      ivmin = -(10**(ni))
      IF (ival.eq.0) THEN
         stri(ni:ni) = '0'
      ELSEIF (ival.le.ivmin .or. ival.ge.ivmax) THEN
         stri(ni:ni) = '*'
      ELSE
         form = '(i'//cv2(ni)//')'
         write (stri, form, err=9) ival
         ENDIF
      GOTO 10
 9    CONTINUE
      stri(ni:ni) = '%'
 10   CONTINUE

      IF (ni+ia.gt.len(aus)) CALL Absturz ('EncI', 'string too long')
      aus (ia:ia+ni-1) = stri
      ia = ia + ni

      END ! EncI

      SUBROUTINE Enc2I (aus, ia, ni, ival)
C     ------------------------------------
            ! 14apr94
         ! encode ival into aus(ia:..)
         ! 2 means INTEGER*2
      CHARACTER     aus*(*), form*12, stri*12, cv2*2
      INTEGER*2     ival

      stri = ' '
      ivmax =   10**(ni)
      ivmin = -(10**(ni))
      IF (ival.eq.0) THEN
         stri(ni:ni) = '0'
      ELSEIF (ival.lt.ivmin .or. ival.gt.ivmax) THEN
         stri(ni:ni) = '*'
      ELSE
         form = '(i'//cv2(ni)//')'
         write (stri, form, err=9) ival
         ENDIF
      GOTO 10
 9    CONTINUE
      stri(ni:ni) = '%'
 10   CONTINUE

      IF (ni+ia.gt.len(aus)) CALL Absturz ('Enc2I', 'string too long')
      aus (ia:ia+ni-1) = stri
      ia = ia + ni

      END ! Enc2I

      SUBROUTINE EncF (aus, ia, ni, nk, val)
C     --------------------------------------
            ! 14apr94
         ! encode val into aus(ia:..)
      CHARACTER     aus*(*), form*12, stri*22, cv2*2
      REAL*8        val, vmax, vmin, dpow0

      IF (nk.ge.ni) CALL Absturz ('EncF','Bad format')
      stri = ' '
      vmax =  dpow0 (10.d0, dfloat(ni-nk-1))
      vmin = -dpow0 (10.d0, dfloat(ni-nk-2))
      IF (val.eq.0.0) THEN
         stri(ni-nk-1:ni-nk) = '0'
      ELSEIF (val.le.vmin .or. val.ge.vmax) THEN
         stri(ni-nk:ni-nk) = '*'
      ELSE
         form = '(f'//cv2(ni)//'.'//cv2(nk)//')'
         write (stri, form, err=9) val
         ENDIF
      GOTO 10
 9    CONTINUE
      stri(ni-nk:ni-nk) = '%'
 10   CONTINUE

      IF (ni+ia.gt.len(aus)) CALL Absturz ('EncF', 'string too long')
      aus (ia:ia+ni-1) = stri
      ia = ia + ni

      END ! EncF

      SUBROUTINE Enc4F (aus, ia, ni, nk, val)
C     ---------------------------------------
            ! 14apr94
         ! encode val into aus(ia:..)
         ! 4 means REAL*4
      CHARACTER     aus*(*), form*12, stri*22, cv2*2
      REAL*4        val
      REAL*8        vmax, vmin, dpow0

      IF (nk.ge.ni) CALL Absturz ('Enc4F','Bad format')
      stri = ' '
      vmax =  dpow0 (10.d0, dfloat(ni-nk-1))
      vmin = -dpow0 (10.d0, dfloat(ni-nk-2))
      IF (val.eq.0.0) THEN
         stri(ni-nk-1:ni-nk-1) = '0'
      ELSEIF (val.le.vmin .or. val.ge.vmax) THEN
         stri(ni-nk:ni-nk) = '*'
      ELSE
         form = '(f'//cv2(ni)//'.'//cv2(nk)//')'
         write (stri, form, err=9) val
         ENDIF
      GOTO 10
 9    CONTINUE
      stri(ni-nk:ni-nk) = '%'
 10   CONTINUE

      IF (ni+ia.gt.len(aus)) CALL Absturz ('Enc4F', 'string too long')
      aus (ia:ia+ni-1) = stri
      ia = ia + ni

      END ! Enc4F

C  --------------------------------------------------------------------
C        1.1.4.   char -> int : ichar1
C  --------------------------------------------------------------------

      INTEGER FUNCTION ichar1 (c)
         ! convert '0' .. '9' to 0 .. 9 ; on error, return 1000+..
      CHARACTER c*(*)

      IF (len(c).ne.1) THEN
         ichar1 = 1000+len(c)
      ELSE
         i = ichar(c) - ichar('0')
         IF (i.lt.0 .or. i.gt.9) THEN
            ichar1 = 2000+i
         ELSE
            ichar1 = i
            ENDIF
         ENDIF
      END

C  --------------------------------------------------------------------
C        1.1.5.   bool -> int  : intq
C  --------------------------------------------------------------------

      INTEGER FUNCTION intq (q) ! JWu90
         LOGICAL q
         IF (q) THEN
            intq = 1  !  on
         ELSE
            intq = 0  !  off
            ENDIF
         END

C  --------------------------------------------------------------------
C        1.1.6.   int -> bool    : qintr
C  --------------------------------------------------------------------

      LOGICAL FUNCTION qintr (i) ! JWu90
         qintr = (i.ne.0) ! C-convention
         END

C  --------------------------------------------------------------------
C        1.1.7.   date & time -> string
C  --------------------------------------------------------------------

      CHARACTER*3 FUNCTION Month (im)
C     -------------------------------
      CHARACTER MDAT(12)*3
      DATA   MDAT / 'Jan', 'Feb', 'Mrz', 'Apr',
     *              'Mai', 'Jun', 'Jul', 'Aug',
     *              'Sep', 'Okt', 'Nov', 'Dez' /
      IF (im.lt.1 .or. im.gt.12) THEN
         Month = '***'
      ELSE
         Month = MDAT(im)
         ENDIF
      END ! Month

      CHARACTER*10 FUNCTION Datum (ilen)
C     ----------------------------------
      CHARACTER Month*3, cr2*2, cv2*2

      CALL SysDat (id, im, iy)
      IF (ilen.le.7) THEN
         Datum = cr2(id) // Month(im) // cv2(mod(iy,100))
      ELSE
         Datum = cr2(id) // '.' // Month(im) // ' ' // cv2(mod(iy,100))
         ENDIF

      END ! Datum

      CHARACTER*8 FUNCTION Zeit (ilen)
C     --------------------------------
      CHARACTER sh*2, sm*2, ss*2
      CALL SysTim (sh, sm, ss)
      IF (ilen.le.7) THEN
         Zeit = sh // ':' // sm
      ELSE
         Zeit = sh // ':' // sm // ':' //ss
         ENDIF
      END ! Zeit

C  ====================================================================
C        1.2.   Integer functions :
C  ====================================================================

C  --------------------------------------------------------------------
C        1.2.1.    int -> int
C  --------------------------------------------------------------------

      INTEGER FUNCTION iinside (i, imi, ima)
         ! returns i, as far as imi=<i=<ima
         iinside = max0 ( min0(i,ima), imi )
         END

C  --------------------------------------------------------------------
C        1.2.2.    real -> int
C  --------------------------------------------------------------------

      INTEGER FUNCTION iPint (r)
         ! preceding integer
         REAL*8 r
         ii = idInt(r)
         IF (r.ge.0.) THEN
            iPint = ii
         ELSE
            IF (dfloat(ii).eq.r) THEN
               iPint = ii
            ELSE
               iPint = ii - 1
               ENDIF
            ENDIF
         END ! iPint

C  --------------------------------------------------------------------
C        1.2.3.    on the decimal representation : jdigit, ndigits
C  --------------------------------------------------------------------

      INTEGER FUNCTION jdigit ( n, j ) ! j-letzte Stelle von n. JWu89
         CHARACTER muell *80
         IF (j.le.0 .or. j.ge.19) THEN
            write (muell,'(a,i7,a)') 'j = ',j,'; zulaessig: 1..19'
            CALL Absturz1 ('jdigit',muell)
            ENDIF
         nloc   = n
         kv     = 10 ** (j-1)
         kvv    = 10 ** (j)
         nloc   = nloc - kvv * ( nloc / kvv )
         jdigit = nloc / kv
         END

      INTEGER FUNCTION ndigits ( n )   ! # Stellen von n; -# f"ur n<0. JWu89
         IF     (n.gt.0) THEN
            ndigits = idint (dlog10(n+1.d-12)) + 1
         ELSEIF (n.eq.0) THEN
            ndigits = 0
         ELSE
            ndigits = - idint (dlog10(-n+1.d-12)) - 1
            ENDIF
         END ! ndigits

C  --------------------------------------------------------------------
C        1.2.4.    on the binary representation : SetBit, iGetBit
C  --------------------------------------------------------------------
              ! JWu 25jan95
         ! the low-order bit is position 0.

      SUBROUTINE BitSet (iLong, nPos, iVal)

      IF     (iVal.eq.0) THEN
         iLong = ibclr (iLong, nPos)
      ELSEIF (iVal.eq.1) THEN
         iLong = ibset (iLong, nPos)
      ELSE
         CALL Absturz1 ('SetBit', 'Illegal bit value')
         ENDIF

      END ! BitSet

      INTEGER FUNCTION iBitGet (iLong, nPos)

      iBitGet = ibits (iLong, nPos, 1)

      END ! iBitGet

      LOGICAL FUNCTION qBitGet (iLong, nPos)

      qBitGet = (iBitGet (iLong, nPos) .ne. 0)

      END ! qBitGet

C  ====================================================================
C        1.3.  Reelle Rechnungen :
C  ====================================================================

            !     Note : the allowed REAL*8 data range is
            !            machine dependent; in VAX FORTRAN,
            !            there are even two different compiler
            !            options. In the default mode (D_FLOATING,
            !            which seems also to be used by NAG)
            !            reals may vary from .29d-38 to 1.7d38

C  --------------------------------------------------------------------
C        1.3.1.   Reelle Funktionen : dpow0, dquot0, dln0, dsqrt0
C  --------------------------------------------------------------------
            !     berechnen einfache mathematische Ausdruecke, und
            !     geben ersatzweise 0. zurueck, wenn der Ausdruck
            !     undefiniert ist (z.B. Division durch 0.)

      REAL*8 FUNCTION dpowii (x,jn)
C     ----------------------------  ! x**c, for Kohlrausch-FT

      INTEGER x,jn, j
      REAL*8 p
      p=1
      DO j=1, jn
         p=p*x
         ENDDO
      dpowii = p

      END ! dpowii

      REAL*8 FUNCTION dquot0 (x,y)          ! x/y
         REAL*8 x,y
         IF (y.eq.0.) THEN
            dquot0 = 0.
         ELSE
            dquot0 = x / y
            ENDIF
         END ! dquot0

      REAl*8 FUNCTION dasin0 (x) ! arc sin (RAD)
          REAL*8 x
          IF ((x.lt.-1.).or.(x.gt.1.)) THEN
             dasin0 = 0.
          ELSE
             dasin0 = dasin(x)
             ENDIF
          END ! dasin0

      REAL*8 FUNCTION dln0 (x)             ! ln x
         REAL*8 x
         IF (x.le.0.) THEN
            dln0 = 0.
         ELSE
            dln0 = dlog (x)
            ENDIF
         END ! dln0

      REAL*8 FUNCTION dlg0 (x)             ! log[10] x
         REAL*8 x
         IF (x.le.0.) THEN
            dlg0 = 0.
         ELSE
            dlg0 = dlog10 (x)
            ENDIF
         END ! dlg0

      REAL*8 FUNCTION dsqrt0 (x)           ! x**1/2
         REAL*8 x
         IF (x.le.0.) THEN
            dsqrt0 = 0.
         ELSE
            dsqrt0 = sqrt (x)
            ENDIF
         END ! dsqrt0

      REAL*8 FUNCTION dmod0 (x,y)  ! x modulo y (1aug91).
         REAL*8 x, y
         IF (y.eq.0.) THEN
            dmod0 = 0.
         ELSE
            dmod0 = dmod (x, y)
            IF (dmod0.lt.0.) dmod0 = dmod0 + dabs(y)
            ENDIF
         END ! dmod0

      REAL*8 FUNCTION dLorentz (x, x1)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         dLorentz = dquot0 (dabs(x1), x**2 + x1**2)
         END ! dLorentz

C  --------------------------------------------------------------------
C        1.3.2.    returning other well defined values
C  --------------------------------------------------------------------

      REAL*8 FUNCTION dquot1 (x,y)         ! O(y)/y
         ! to be used, if x/y->1 for y->0, i.e. x = y + O(y^2)
         REAL*8 x,y
         IF (dabs(y).lt.1.d-14) THEN
            dquot1 = 1.
         ELSE
            dquot1 = x / y
            ENDIF
         END ! dquot1

      REAL*8 FUNCTION dexp1 (x)
         ! To prevent overflow, JWu 16jul91
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         PARAMETER (expmax=40.) ! actually, the maximum is ca. 88
         IF     (x.lt.-expmax) THEN
            dexp1 = 0.
         ELSEIF (x.gt. expmax) THEN
            dexp1 = dexp (expmax)
         ELSE
            dexp1 = dexp(x)
            ENDIF
         END ! dexp1

      REAL*8 FUNCTION dexp2 (x)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         PARAMETER (expmax=20.) ! for even more security 31dec92
         IF     (x.lt.-expmax) THEN
            dexp2 = 0.
         ELSEIF (x.gt. expmax) THEN
            dexp2 = dexp(expmax)
         ELSE
            dexp2 = dexp(x)
            ENDIF
         END ! dexp2

      REAL*8 FUNCTION dgamma1 (x) ! the Gamma function 15dec95
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         IF (x.lt.1d-198 .or. x.gt.120) THEN
            dgamma1 = 1.d199
            RETURN
            ENDIF
         dgamma1 = gamma(x)
         END ! dgamma1

      REAL*8 FUNCTION dnquot (x,n)          ! x/n
         REAL*8 x
         IF (n.eq.0) THEN
            dnquot = 0.
         ELSE
            dnquot = x / n
            ENDIF
         END ! dnquot

C  --------------------------------------------------------------------
C        1.3.4.    like min, max
C  --------------------------------------------------------------------

      REAL*8 FUNCTION dinside (x, xmi, xma)
         ! returns x, as far as xmi=<x=<xma
         IMPLICIT REAL*8 (x)
         dinside = dmax1 ( dmin1(x,xma), xmi )
         END

C  --------------------------------------------------------------------
C        1.3.5.   Rounding operations
C  --------------------------------------------------------------------

      REAL*8 FUNCTION roundN (x, j)         ! auf j Stellen runden
C     -----------------------------
         ! JWu 1989/90 ? Renovation jan92.
      IMPLICIT REAL*8 (a-h,o-p,r-z)

      IF (x.eq.0. .or. j.lt.1) THEN
         roundN = 0.
      ELSEIF (j.gt.8) THEN ! maxint=2.1E9 (18sep91)
         roundN = x
      ELSE
         xa = dabs(x)
         IF (xa.ge.1.) THEN
            n = j - idint(dlog10(xa)) - 1
         ELSE
            n = j - idint(dlog10(xa))
            ENDIF
         xr =  dfloat( idnint(xa * 10.d0**n) ) / 10.d0**n
         roundN = dsign (xr, x)
         ENDIF

      END ! roundN

      COMPLEX*16 FUNCTION roundZN (z, j) ! Re and Im auf j Stellen runden
C     ----------------------------------
         ! JWu 2025
      COMPLEX*16 z
      REAL*8     roundN
      INTEGER    j

      roundZN  = CMPLX(roundN(dble(z), j), roundN(DIMAG(z), j))

      END ! roundZN


      REAL*8 FUNCTION roundUD (x, j, qup) ! round up / down
C     -----------------------------------
         ! JWu 7sep98
      IMPLICIT NONE
      INTEGER          j, n
      REAL*8           x, xa, xr
      LOGICAL          qup

      IF (x.eq.0. .or. j.lt.1) THEN
         roundUD = 0.
      ELSEIF (j.gt.8) THEN ! maxint=2.1E9 (18sep91)
         roundUD = x
      ELSE
         xa = dabs(x)
         IF (xa.ge.1.) THEN
            n = j - idint(dlog10(xa)) - 1
         ELSE
            n = j - idint(dlog10(xa))
            ENDIF
         IF (qup) THEN ! aufrunden
            xr = dfloat( idnint(xa * 10.d0**n) + 1 ) / 10.d0**n
         ELSE          ! abrunden
            xr = dfloat( idnint(xa * 10.d0**n) - 1 ) / 10.d0**n
            ENDIF
         roundUD = dsign (xr, x)
         ENDIF

      END ! roundUD

      SUBROUTINE RoundLin (xmi, xma, relIn, idigIn)
C     ---------------------------------------------
         ! extends the limits [xmi..xma] by a proportion rel,
         ! and rounds the result to idig digits.

      IMPLICIT NONE
      REAL*8        xmi, xmid, xmir, xma, xmad, xmar, rel, dx, e, r,
     *              relIn, h, roundN, roundUD
      INTEGER       idigIn, idig

      idig = idigIn
      rel  = relIn
      IF (rel.lt.0.) CALL Absturz1 ('RoundLin', 'rel o.o.r.')
      IF (idig.lt.0) CALL Absturz1 ('RoundLin', 'idig o.o.r.')

      IF (xmi.eq.xma) THEN
         xma = roundUD (xma, idig, .true.)
         xmi = roundUD (xmi, idig, .false.)
         RETURN
         ENDIF

      IF (xmi.gt.xma) THEN ! sort limits
         h   = xmi
         xmi = xma
         xma =   h
         ENDIF

      dx = xma - xmi
      e  = 1.
      r  = 2.

 13   CONTINUE
         xmid = xmi - rel*dx
         xmad = xma + rel*dx
         IF (dsign(e,xmi-r*dx).ne.dsign(e,xmi)) xmid = 0.  ! 5/12aug91
         IF (dsign(e,xma+r*dx).ne.dsign(e,xma)) xmad = 0.
         IF (idig.gt.0) THEN
            xmir = roundN (xmid,idig)
            xmar = roundN (xmad,idig)
            ENDIF
         IF ((xmir.ge.xmi .or. xmar.le.xma) .and. idig.lt.10) THEN
            rel  = rel / 10
            idig = idig + 1
            GOTO 13
            ENDIF

      IF (xmir.eq.xmar) THEN
         xmir = - 1.5 * dabs(xma)
         xmir = roundN (xmir, 3)
         xmar = - xmir
         ENDIF

      xmi = xmir
      xma = xmar

      END ! RoundLin

      SUBROUTINE RoundLog (xmi, xma, relIn, idigIn)
C     ---------------------------------------------
         ! extends the limits [xmi..xma] by a proportion rel,
         ! and rounds the result to idig digits.

      IMPLICIT REAL *8 (a-h,o-p,r-z)

      idig = idigIn
      rel  = relIn
      IF (rel.lt.0.)  CALL Absturz1 ('RoundLog', 'rel o.o.r.')
      IF (idig.lt.0)  CALL Absturz1 ('RoundLog', 'idig o.o.r.')
 7    CONTINUE
      IF (xmi.lt.0.0) THEN
         Print *, 'SEVERE/ PROGRAM ERROR/ RoundLog/ min<0 :', xmi
         CALL Gong(19)
         xmi = -xmi
         ENDIF
      IF (xmi.gt.xma) THEN
         Print *, 'SCHOENHEITSFEHLER/ RoundLog/ max<min'
         CALL Gong (3)
         h   = xmi
         xmi = xma
         xma =   h
         GOTO 7
         ENDIF

      dx  = xma / xmi
      xmir= xmi / dx**rel
      xmar= xma * dx**rel
      IF (idig.gt.0) THEN
         xmir = roundN (xmir,idig)
         xmar = roundN (xmar,idig)
         ENDIF
      xmi = xmir
      xma = xmar
      IF (xmi.eq.xma) THEN
         xma = dabs(xma) * sqrt(10d0)
         xma = roundN (xma, 2)
         xmi = xma / 10
         ENDIF

      IF (xma.le.xmi .or. xmi.le.0) THEN
         Print *, 'SEVERE/ PROGRAM ERROR/ RoundLog/ must overwrite:'
         Print *, ' xmi = ', xmi, ', xma = ', xma
         CALL Gong (19)
         xmi = .1
         xma = 10
         ENDIF

      END ! RoundLog

C  --------------------------------------------------------------------
C  1.3.6. Interpolation
C  --------------------------------------------------------------------

      SUBROUTINE LinIntPol (x, xl, xh, yl, yh, dl, dh, y, d)
C     ------------------------------------------------------
         ! JWu 31mai91
         ! linear interpolation for xl < x < xh
         ! input : two data points with error (xl,yl,dl), (xh,yh,dh),
         !         new x
         ! output: y,d for the new x

      IMPLICIT REAL *8 (a-h,o-p,r-z)

      IF (xl.eq.xh) THEN
         ! simply the average of yl and yh :
         y = .5 * (yl + yh)
         d = .5 * dsqrt (dl**2 + dh**2)
      ELSE
         ! normal case
         al = (xh-x) / (xh-xl)
         ah = (x-xl) / (xh-xl)
c        y  =        al    * yl    + ah    * yh ! ill-conditioned
         y  = yl + ah * (yh-yl) ! 6sep93 This form should guarantee that
                                !        y in (yl:yh) if x in (xl:xh).
         d  = dsqrt (al**2 * dl**2 + ah**2 * dh**2)
         ENDIF

      END ! LinIntpol

      SUBROUTINE LinIntPol1 (x, xl, xh, yl, yh, y)
C     --------------------------------------------
            ! JWu 21mai95
         ! Accelerated version of LinIntPol : no error calculation

      IMPLICIT REAL *8 (a-h,o-p,r-z)

      IF (xl.eq.xh) THEN
         ! simply the average of yl and yh :
         y = .5 * (yl + yh)
      ELSE
         ! normal case
         al = (xh-x) / (xh-xl)
         ah = (x-xl) / (xh-xl)
c        y  =        al    * yl    + ah    * yh ! ill-conditioned
         y  = yl + ah * (yh-yl) ! 6sep93 This form should guarantee that
                                !        y in (yl:yh) if x in (xl:xh).
         ENDIF

      END ! LinIntpol1

      SUBROUTINE LinIntPolArray (X, Y, n, xi, yi, ovh, Fehler)
C     --------------------------------------------------------
            ! JWu 19feb92.
         ! Interpolation within array (X,Y,1..n).
         ! Extrapolation is forbidden, but
         ! an overhang ovh*dx can be allowed.

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER        Fehler*(*), aux*20
      DIMENSION        X(n), Y(n)

C  Checks :
      IF (Fehler.ne.'&ff') THEN
         RETURN
      ELSEIF (ovh.lt.0.) THEN
         CALL Absturz1 ('LinIntPolArray', 'negative overhang')
      ELSEIF (n.lt.2) THEN
         Fehler = ' LinIntPolArray/ n<2'
         RETURN
      ELSEIF (irSorted(X,n).ne.2) THEN
         Fehler = ' LinIntPolArray/ at present, data must be sorted'
         RETURN
         ENDIF

C  Find neighbours :
      il = irPosOpt (X, n, xi, 'r', il)
      IF     (il.lt.1) THEN
         IF (xi.ge.X(1)-ovh*(X(2)-X(1))) THEN
            il = 1
         ELSE
            CALL NiceNum (X(1)-ovh*(X(2)-X(1)), aux, ia)
            Fehler = ' LinIntPolArray/ x below '//aux(1:ia)
            RETURN
            ENDIF
      ELSEIF (il.gt.n-1) THEN
         IF (xi.le.X(n)+ovh*(X(n)-X(n-1))) THEN
            il = n-1
         ELSE
            CALL NiceNum (X(n)+ovh*(X(n)-X(n-1)), aux, ia)
            Fehler = ' LinIntPolArray/ x above '//aux(1:ia)
            RETURN
            ENDIF
         ENDIF

C  Interpolate :
      CALL LinIntPol (xi, X(il), X(il+1), Y(il), Y(il+1),
     *                dl, dh, yi, di)

      END ! LinIntPolArray

      REAL*8 FUNCTION rIntegralXY (X, Y, n)
C     -------------------------------------
            ! JWu 3feb95
         ! Estimate the integral I dx Y(x)

      IMPLICIT NONE
      REAL*8  X(*), Y(*), sum
      INTEGER n, i

      sum = 0
      DO i = 2, n-1
         sum = sum + (X(i)-X(i-1)) * (Y(i)+Y(i-1))/2
         ENDDO
      rIntegralXY = sum

      END ! rIntegralXY

C  --------------------------------------------------------------------
C  1.3.7. With error propagation :
C  --------------------------------------------------------------------
         ! This section opened 5aug91.
         ! Input :  pairs y1,d1; y2,d2;...
         ! Output : pair  y0,d0
         ! WARNING : the subroutines have to be written such that
         !           they may be called with an input argument
         !           equal to the output argument, e.g. DEquot (y,d,y,d,z,0)

      SUBROUTINE DEquot0 (y0,d0, y1,d1, y2,d2)
            ! y0 = y1/y2. - 5aug91
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         y = dquot0 (y1,y2)
         d = dquot0 (1.d0,y2) * dsqrt (d1**2 + y**2 * d2**2)
         y0= y
         d0= d
         END ! DEquot0

C  --------------------------------------------------------------------
C  L1.3.8.  Mathematical and physical constants
C  --------------------------------------------------------------------

      BLOCK DATA R8Const
C     ------------------
      IMPLICIT REAL*8 (a-h,o-p,r-z)

      COMMON / MathConst / twopi
      COMMON / PhysConst / val_hbar, val_kB

      DATA   twopi    / 6.28318530718  /
      DATA   val_hbar / 0.658218   /,  ! meV*psec
     *       val_kB   / 0.08617    /   ! meV/K

      END ! R8Const


C  --------------------------------------------------------------------
C  1.5.   Logical Functions :
C  --------------------------------------------------------------------

      LOGICAL FUNCTION qEqEps (x,y)
C     -----------------------------
         ! ( x .eq. y ) with precision eps
            ! qrEq JWu 3.11.90, renamed 6sep91
      REAL *8   eps, x, y, d, s
      eps = 1.d-11
      d = dAbs (x-y)
      s = dAbs(x) + dAbs(y)
      qEqEps = (d.le.s*eps)
      END ! qEqEps

      LOGICAL FUNCTION qEqTol (x,y, tol)
C     ----------------------------------
         ! ( x .eq. y ) with precision tol
            ! JWu 6sep91
      REAL *8   tol, x, y, d, s
      d = dAbs (x-y)
      s = dAbs(x) + dAbs(y)
      qEqTol = (d.le.s*tol)
      END ! qEqTol

      LOGICAL FUNCTION qioutside (i, imi, ima)
C     ---------------------------------------- ! JWu15nov90
         qioutside = (i.lt.imi .or. ima.lt.i)
         END ! qioutside

      LOGICAL FUNCTION qiinside (i, imi, ima)
C     ---------------------------------------- ! JWu19nov90
         qiinside = (imi.le.i .and. i.le.ima)
         END ! qiinside

      LOGICAL FUNCTION qroutside (r, rmi, rma)
C     ----------------------------------------  ! JWu15nov90
         REAL*8 r, rmi, rma
         qroutside = r.lt.rmi .or. rma.lt.r
         END ! qroutside

      LOGICAL FUNCTION qrInClos (r, rmi, rma)
C     ---------------------------------------
         ! r in closed interval [rmi,rma]
            ! qrinside 11feb91, negative direction 31aug93, qrInClos 21may94
      REAL*8 r, rmi, rma
      qrInClos = (rmi.le.r .and. r.le.rma) .or.
     *           (rma.le.r .and. r.le.rmi)
      END ! qrInClos

      LOGICAL FUNCTION qrInOpen (r, rmi, rma)
C     --------------------------------------- ! JWu 21may94
         ! r in open interval ]rmi,rma[
      REAL*8 r, rmi, rma
      qrInOpen = (rmi.lt.r .and. r.lt.rma) .or.
     *           (rma.lt.r .and. r.lt.rmi)
      END ! qrInOpen

      LOGICAL FUNCTION qrinside (r, rmi, rma)
C     --------------------------------------- ! JWu11feb91, l"auft aus
      REAL*8 r, rmi, rma
      qrinside = (rmi.le.r  .and. r.le.rma) .or.
     *           (rma.le.r  .and. r.le.rmi)
      END ! qrinside

      LOGICAL FUNCTION qrInRange (r, rmi, rma)
C     ---------------------------------------- ! JWu13sep91
         ! Convention : rmi=rma=0. means full range
      REAL*8 r, rmi, rma
      qrInRange = (r.ge.rmi .and. r.le.rma) .or.
     *            (rmi.eq.0. .and. rma.eq.0.)
      END                       ! qrInRange

C  ====================================================================
C        1.7.  Vector Operations
C                 MinMax, Cycling, Decompose
C  ====================================================================

C  --------------------------------------------------------------------
C        1.7.1.   INTEGER (i..)
C  --------------------------------------------------------------------

      SUBROUTINE iCopy (II, ivon, ibis, incr, JJ, jvon, jncr)
C     -------------------------------------------------------
         ! II(i) := JJ(j)  (JWu 28jan92)
      INTEGER II(*), JJ(*)
      DO k = 0, (ibis-ivon)
         II(ivon+k*incr) = JJ(jvon+k*jncr)
         ENDDO
      END ! iCopy

      INTEGER FUNCTION iSum (I, jvon, jbis, jncr)
C     ------------------------------------------- ! JWu13dec90
         INTEGER I(*)
         iS = 0
         DO j = jvon, jbis, jncr
            iS = iS + I(j)
            ENDDO
         iSum = iS
         END ! iSum

      SUBROUTINE iMinMax (I, n, imi, ima, jmi, jma)
C     ---------------------------------------------
         INTEGER    I(*)
         imi = I(1)
         jmi = 1
         ima = I(1)
         jma = 1
         DO j = 2, n
            IF (I(j).lt.imi) THEN
               imi = I(j)
               jmi = j
               ENDIF
            IF (I(j).gt.ima) THEN
               ima = I(j)
               jma = j
               ENDIF
            ENDDO
         END ! iMinMax

C  --------------------------------------------------------------------
C        1.7.2.   Real *8  (r..)
C  --------------------------------------------------------------------

      SUBROUTINE rCopy (Y, ivon, ibis, incr, X, jvon, jncr)
C     -----------------------------------------------------
            ! JWu15nov90. Correct 13mar92.
         ! Y := X
         REAL *8 X(*), Y(*)
         IF (incr.eq.0) CALL Absturz1 ('rCopy', 'incr=0')
         DO k = 0, (ibis-ivon)/incr
            Y(ivon+k*incr) = X(jvon+k*jncr)
            ENDDO
         END ! rCopy

      REAL*8 FUNCTION rSum  (X, ivon, ibis, incr)    ! = sum_i X(i)
C     -------------------------------------------
            ! JWu 21jan91
         REAL*8 rs, X(*)
         rs = 0.
         DO i = ivon, ibis, incr
            rs = rs + X(i)
            ENDDO
         rSum = rs
         END ! rSum

      SUBROUTINE rSet (Y, ivon, ibis, incr, x)
C     ----------------------------------------
            ! JWu15nov90. Correct 10mar92.
         ! Y := x
         REAL *8 Y(*), x
         DO i = ivon, ibis, incr
            Y(i) = x
            ENDDO
         END ! rSet

      SUBROUTINE rMinMax (X, n, xmi, xma, jmi, jma)
C     ---------------------------------------------
         REAL *8    X(*), xmi, xma
         xmi = X(1)
         jmi = 1
         xma = X(1)
         jma = 1
         DO j = 2, n
            IF (X(j).lt.xmi) THEN
               xmi = X(j)
               jmi = j
               ENDIF
            IF (X(j).gt.xma) THEN
               xma = X(j)
               jma = j
               ENDIF
            ENDDO
         END ! rMinMax

      REAL*8 FUNCTION rStepAt (X, n, i)
C     ---------------------------------
            ! X-stepwidth at X(i). - JWu 6aug91
         REAL *8    X(*)
         IF     (n.le.1) THEN
            rStepAt = 0.
         ELSEIF (i.eq.1) THEN
            rStepAt = X(2) - X(1)
         ELSEIF (i.eq.n) THEN
            rStepAt = X(n) - X(n-1)
         ELSE
            rStepAt = (X(i+1) - X(i-1)) / 2
            ENDIF
         END ! rStepAt

      INTEGER FUNCTION irPos (X, n, val, dir)
C     ---------------------------------------
         ! HPS =< 1987, renewed JWu nov90, 16feb91
      ! This function returns the position of value val in array X.
      ! X contains n elements and has to be given in ascending order.
      ! dir determines whether the position to be returned is that
      ! of the 'r'ight or the 'l'eft or the 'n'earest neighbour :
      ! dir = 'r' :                     X(irPos) <= val < X(irPos+1)
      !       'l' : X(irPos-1) < val <= X(irPos)
      !       'n' :            | val -  X(irPos) | -> min !

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      REAL *8           X(*)
      CHARACTER*1       dir

      IF     (dir.eq.'l') THEN
         irPos = 0
10       CONTINUE
            irPos = irPos + 1
            IF (irPos.gt.n) RETURN
            IF (X(irPos).lt.val) GOTO 10

      ELSEIF (dir.eq.'r') THEN
         irPos = n+1
30       CONTINUE
            irPos = irPos - 1
            IF (irPos.LT.1) RETURN
            IF (X(irPos).gt.val) GOTO 30

      ELSEIF (dir.eq.'n') THEN
         xp = X(1)
         IF (val.le.xp) THEN
            irPos = 1
            RETURN
            ENDIF
         DO i = 2, n
            xm = xp
            xp = X(i)
            IF (xm.gt.xp) CALL Absturz1 ('irPos', 'unsorted data')
            IF (xm.le.val .and. val .le.xp) THEN
               IF ((val-xm).lt.(xp-val)) THEN
                  irPos = i-1
               ELSE
                  irPos = i
                  ENDIF
               RETURN
               ENDIF
            ENDDO
         irPos = n

      ELSE
         CALL Absturz1 ('irPos', 'dir o.o.r.')
         ENDIF

      END ! irPos

      INTEGER FUNCTION irPosOpt (X, n, val, dir, iGuessIn)
C     ----------------------------------------------------
         ! HPS =< 1987, renewed JWu nov90, 16feb91
         ! optimized algorithm 3sep91
      ! This function returns the position of value val in array X.
      ! X contains n elements and has to be given in ascending order.
      !    !!! the ascending ORDER of X will NOT be CHECKED !!!
      ! dir determines whether the position to be returned is that
      ! of the 'r'ight or the 'l'eft or the 'n'earest neighbour :
      ! dir = 'r' :                     X(irPos) <= val < X(irPos+1)
      !       'l' : X(irPos-1) < val <= X(irPos)
      !       'n' :            | val -  X(irPos) | -> min !
      ! iGuess contains an initial guess for irPos.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      REAL *8           X(*)
      CHARACTER*1       dir

      iGuess = iinside ( iGuessIn, 1, n)

      IF     (dir.eq.'l') THEN
         irPosOpt = 0
10       CONTINUE
            irPosOpt = irPosOpt + 1
            IF (irPosOpt.gt.n) RETURN
            IF (X(irPosOpt).lt.val) GOTO 10

      ELSEIF (dir.eq.'r') THEN
         IF (val.ge.X(iGuess)) THEN
            irPosOpt = iGuess - 1
31          CONTINUE
            irPosOpt = irPosOpt + 1
            IF (irPosOpt.ge.n) RETURN
            IF (val.ge.X(irPosOpt+1)) GOTO 31
         ELSE
            irPosOpt = iGuess
32          CONTINUE
            irPosOpt = irPosOpt - 1
            IF (irPosOpt.lt.1) RETURN
            IF (val.lt.X(irPosOpt)) GOTO 32
            ENDIF

      ELSEIF (dir.eq.'n') THEN ! 6sep91 -  even simplified by the opimization
         IF (val.lt.X(iGuess)) THEN   ! down search
            DO i = iGuess, 2, -1
               IF (X(i)-val.le.val-X(i-1)) THEN
                     ! note that X(i)-val may also be negative
                  irPosOpt = i
                  RETURN
                  ENDIF
               ENDDO
            irPosOpt = 1
         ELSE                         ! up search
            DO i = iGuess, n-1
               IF (val-X(i).lt.X(i+1)-val) THEN
                  irPosOpt = i
                  RETURN
                  ENDIF
               ENDDO
            irPosOpt = n
            ENDIF

      ELSE
         CALL Absturz1 ('irPosOpt', 'dir o.o.r.')
         ENDIF

      END ! irPosOpt

      INTEGER FUNCTION irSorted (X, n)
C     --------------------------------
            ! JWu 6jul91
         ! check the order of X(1,..,n) :
         ! (-2) x1>x2>.., (-1)x1>=x2>=.., (0) unsorted, (1) x1<=x2<=..,
         !  (2) x1<x2<.., (3) only one value.

      REAL*8 X(*)

      IF (n.le.1) THEN
         irSorted = 3
         RETURN
         ENDIF
      IF (X(1).gt.X(n)) THEN ! descending order
         irSorted = -2
         DO i = 2, n
            IF (X(i).eq.X(i-1)) THEN
               irSorted = -1
            ELSEIF (X(i).gt.X(i-1)) THEN
               irSorted = 0
               RETURN
               ENDIF
            ENDDO
      ELSE                   ! ascending order
         irSorted = 2
         DO i = 2, n
            IF (X(i).eq.X(i-1)) THEN
               irSorted = 1
            ELSEIF (X(i).lt.X(i-1)) THEN
               irSorted = 0
               RETURN
               ENDIF
            ENDDO
         ENDIF

      END ! irSorted

      SUBROUTINE CheckScale (n, X, tolerance, qEqui, delX)
C     ----------------------------------------------------
      ! check whether the X values are equidistant and determine
      ! the step width delX

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)
      DIMENSION         X(*)

      IF (n.le.1) THEN
         qEqui = .false.
         RETURN
         ENDIF
      delX = ( X(n) - X(1) ) / ( n - 1 )
      qEqui = .true.
      DO i = 1, n-1
         delXi = X(i+1) - X(i)
         IF ( dAbs(delXi - delX) .gt. tolerance*dAbs(delX) ) THEN
            qEqui = .false.
            RETURN
            ENDIF
         ENDDO

      END ! CheckScale

C  --------------------------------------------------------------------
C        1.7.4.   Logical vectors
C  --------------------------------------------------------------------

      INTEGER FUNCTION iqSum (Q, n) ! number of true Q(i). JWu 4jul91.
         LOGICAL Q(*)
         iqSum = 0
         DO i = 1, n
            IF (Q(i)) iqSum = iqSum + 1
            ENDDO
         END ! iqSum

      SUBROUTINE qCopy (Qi, ivon, ibis, incr, Qj, jvon, jncr) ! JWu 2apr92
         LOGICAL Qj(*), Qi(*)
         IF (incr.eq.0) CALL Absturz1 ('qCopy', 'incr=0')
         DO k = 0, (ibis-ivon)/incr
            Qi(ivon+k*incr) = Qj(jvon+k*jncr)
            ENDDO
         END ! qCopy

      SUBROUTINE qSet (Q, ivon, ibis, incr, qVal) ! JWu 8nov91
         LOGICAL Q(*), qVal
         DO i = ivon, ibis, incr
            Q(i) = qVal
            ENDDO
         END ! qSet

C  ====================================================================
C        1.9.   Control : Absturz1
C  ====================================================================

         SUBROUTINE Absturz1 (wo, warum)
C        ------------------------------                                    oct90
         CHARACTER*(*) wo, warum
         CHARACTER     OS*4
         COMMON /Impltn / OS

C  Overlay out, full scroll, goto last line :
         Print '(1x,5a1,$)',
     *      Char(29),Char(27),Char(92),Char(51),Char(24)
         Print '(1x,2(a1,a))',
     *      Char(27), '[0;25r', Char(27), '[25;1H'

C  Error message :
         Print *, 'CRASH in : ', wo
         Print *, 'MESSAGE  : ', warum
         Print *, 'Now provocing another error in order'//
     *           ' to obtain a system traceback :'
C  Crash :
         IF (OS.eq.'VMS') THEN
            Print *, 'Now provocing another error in order'//
     *              ' to obtain a system traceback :'
            read (wo, '(i3)') wo ! a guaranteed error
         ELSEIF (OS.eq.'X11') THEN
            Print *, ' Hit RETURN to stop'
            Read (*,*)
            STOP
         ELSE
            STOP
            ENDIF

         END ! Absturz1
