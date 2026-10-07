C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i70   :     rescaling operations, convolutions pp.
C
C  ====================================================================

C     Contents :
C        1.  Rescaling:
C               TraUnits, TraAchseX
C        2.  Convolution
C               TraSymm, TraDouble, TraRepr, TraDeconv, TraConv
C        3.  Monte-Carlo Convolution
C               MC_Conv

C     Aenderungsverzeichnis :
C     JWu   nov91 : Ida5.182 -> Ida7
C     JWu   sep91 : slow FT (non equidistant)
C     JWu   jul91 : new FFT routine
C     JWu   jun91 : -> Modul Ida5; Trafo w<->tof
C     JWu   nov90 : FT begonnen
C  16.02.2026 Artem Panchenko: Corrected several line breaks
C  ====================================================================
C  i70 / 1 :   rescaling
C  ====================================================================

C  --------------------------------------------------------------------
      SUBROUTINE TraUnits (nJlist, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 10jun96.
         ! Nonlinear transform of x axis.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)
      CHARACTER*20   coX, coY, unX, unY, coX2, coY2, unX2, unY2

      DATA          iTra /1/

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop over files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (j, 'x', coX, unX, Fehler)
         CALL OlfCnuG (j, 'y', coY, unY, Fehler)
         coX2 = coX
         coY2 = coY

         ! Identify dimension of x :
 12      CONTINUE
         IF     (coX(1:1).eq.'t') THEN
            idX = 1
         ELSEIF (coX(1:1).eq.'w' .or. coX(1:1).eq.'f') THEN
            idX = 2
         ELSE
            IF (coX.ne.' ') Print *, ' input x coordinate is ', coX
            Print *, ' Dimension is'
            Print *, '    (1) time            (2) frequency/energy'
            Print *, '    (0) anything else'
            idIn = iAskDMu (' Option', idX, -1, 2)
            IF (idIn.eq.-1) THEN
               Fehler = ' '
               RETURN
               ENDIF
            idX = idIn
            ENDIF

         ! Identify unit of x :
         CALL UnitSI ('input x coordinate', unX, idX, siX)
         IF (siX.le.0) THEN
            coX = ' '
            GOTO 12 ! retry (maybe, coX was wrongly identified)
            ENDIF

         ! Identify unit of y :
         ! ... should lead to guess : ipYguess
         ipYguess = -1

         ! Ask for dimension of y :
         Print *, ' input y coordinate is '//coY
         ipY = iAskD (' Scales with which power of x', ipYguess)

         ! Ask for new unit of x, calculate conversion factor :
         IF (idX.ne.idXalt) unX2 = ' '
         CALL FrageCD (' New unit of x', unX2, unX2)
         CALL UnitSI ('output x coordinate', unX2, idX, siX2)
         IF (siX2.le.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         facX = siX / siX2

         ! New unit of y :
         facY = facX**ipY
         CALL ReplaceT (unY, unX, unX2)
         CALL FrageCD (' New unit of y', unY2, unY)

         CALL OlfCnuP (jout, 'x', coX2, unX2, Fehler)
         CALL OlfCnuP (jout, 'y', coY2, unY2, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! Now the conversion :
         DO K = 1, nK

            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            DO i = 1, n
               X(i) = X(i) * facX
               Y(i) = Y(i) * facY
               D(i) = D(i) * facY
               ENDDO ! i

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         CALL OlfClos (jout, nK, Fehler)

         idXalt = idX
         ENDDO ! files

      END ! TraUnits

C  --------------------------------------------------------------------
      SUBROUTINE UnitSI (coi, uni, idi, sii)
C  --------------------------------------------------------------------
            ! JWu 11jun96.
         ! Identify unit (uni) of a coordinate (coi) of given dimension (idi),
         ! express the result (sii) in SI units.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      CHARACTER*(*)  coi, uni

      DATA          twopi /6.2831853/, hbar /0.658218/

      IF     (idi.eq.1) THEN ! time
         IF     (uni.eq.'sec') THEN
            sii = 1
         ELSEIF (uni.eq.'msec') THEN
            sii = 1.d-3
         ELSEIF (uni.eq.'usec') THEN
            sii = 1.d-6
         ELSEIF (uni.eq.'nsec') THEN
            sii = 1.d-9
         ELSEIF (uni.eq.'psec') THEN
            sii = 1.d-12
         ELSE
            CALL Say3 (' unit of '//coi, ' is "'//uni, '"')
            s1 = rAskD (' How many sec is that', s1)
            sii = s1
            ENDIF

      ELSEIF (idi.eq.2) THEN ! frequency/energy
         IF     (uni.eq.'Hz') THEN
            sii = 1     * twopi
         ELSEIF (uni.eq.'kHz') THEN
            sii = 1.d3  * twopi
         ELSEIF (uni.eq.'MHz') THEN
            sii = 1.d6  * twopi
         ELSEIF (uni.eq.'GHz') THEN
            sii = 1.d9  * twopi
         ELSEIF (uni.eq.'THz') THEN
            sii = 1.d12 * twopi
         ELSEIF (uni.eq.'sec-1') THEN
            sii = 1.d0
         ELSEIF (uni.eq.'msec-1') THEN
            sii = 1.d3
         ELSEIF (uni.eq.'usec-1') THEN
            sii = 1.d6
         ELSEIF (uni.eq.'nsec-1') THEN
            sii = 1.d9
         ELSEIF (uni.eq.'psec-1') THEN
            sii = 1.d12
         ELSEIF (uni.eq.'eV') THEN
            sii = 1.d15 / hbar
         ELSEIF (uni.eq.'meV') THEN
            sii = 1.d12 / hbar
         ELSEIF (uni.eq.'ueV') THEN
            sii = 1.d9  / hbar
         ELSE
            CALL Say3 (' unit of '//coi, ' is "'//uni, '"')
            s2 = rAskD (' How many sec-1 is that', s2)
            sii = s2
            ENDIF

      ELSEIF (idi.eq.0) THEN ! "Apfel oder Birnen
            CALL Say3 (' unit of '//coi, ' is "'//uni, '"')
            s0 = rAskD (' How many basic units is that', s0)
            sii = s0

      ELSE
         CALL Absturz ('UnitSI', 'dimension idi unknown')

         ENDIF

      END ! UnitSI

C  --------------------------------------------------------------------
      SUBROUTINE UnitConv (uni, uno, fac, Fehler)
C  --------------------------------------------------------------------

            ! JWu 20feb97
         ! determine converse factor fac from units uni to uno
         ! such that data(uno) = data(uni) * fac

      IMPLICIT NONE
      CHARACTER*(*)  uni, uno, Fehler
      REAL*8         fac, twopi, hbar
      INTEGER        MUC, ii, io
      PARAMETER     (MUC=8)
      CHARACTER*8    Units(MUC)
      INTEGER        Dimns(MUC)
      REAL*8         Facts(MUC)
      LOGICAL        FirstCall

      DATA           twopi /6.2831853/, hbar /0.658218/,
     *               FirstCall /.true./
      DATA           Units / 'neV', 'ueV', 'meV', 'eV',
     *                       'MHz', 'GHz', 'THz', 'msec/m' /
      DATA           Dimns / 2, 2, 2, 2,
     *                       2, 2, 2, 2    /

      IF (FirstCall) THEN
         Facts( 1) = 1.d6  / hbar
         Facts( 2) = 1.d9  / hbar
         Facts( 3) = 1.d12 / hbar
         Facts( 4) = 1.d15 / hbar
         Facts( 5) = 1.d6  * twopi
         Facts( 6) = 1.d9  * twopi
         Facts( 7) = 1.d12 * twopi
         Facts( 8) = 1.d12 / hbar
         ENDIF

      DO ii = 1, MUC
         IF (uni.eq.Units(ii)) GOTO 11
         ENDDO
      CALL Compose2 (Fehler,'UnitConv/ input unit '//uni,' not known')
      RETURN
 11   CONTINUE

      DO io = 1, MUC
         IF (uno.eq.Units(io)) GOTO 12
         ENDDO
      CALL Compose2 (Fehler, 'UnitConv/ output unit '//
     *               uno, ' not known')
      RETURN
 12   CONTINUE

      IF (Dimns(ii).ne.Dimns(io)) THEN
         CALL Compose2 (Fehler,
     *        'unitConv/ different dimensions: '//uni, ' <-> '//uno)
         RETURN
         ENDIF

      fac = Facts(ii) / Facts(io)

      END ! UnitConv

C  --------------------------------------------------------------------
      SUBROUTINE TraAchseX (nJlist, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 20jun91.
         ! Transforms the scale of the x-axis (e.g. from energy to
         ! time-of-flight). The y-data are transformed accordingly.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'

      CHARACTER*(*)  Fehler
      CHARACTER*40   coX, coY, unX, unY
      INTEGER        JList(*)

      DATA          iTra /1/

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Choose mode :
      Print *, ' transforms :'
      Print *,
     * '    (1)  w   -> tof   (-1)  tof -> w   (S(q,w)<->tof-counts)'

      iTra = iAskDMu (' Option', iTra, -1, 1)
      IF (iTra.eq.0) RETURN

C  Loop over files :
      DO lj = 1, nJList
      j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (j, 'x', coX, unX, Fehler)
         CALL OlfCnuG (j, 'y', coY, unY, Fehler)

C  Documentation :
         IF     (iTra.eq. 1) THEN
            IF (coY.ne.'S(q,w)' .or. unY.ne.'meV-1') THEN
               Fehler = ' call with S(q,w) [meV-1]'
               RETURN
               ENDIF
            CALL OlfCnuP (jout, 'x', 'tof', 'msec/m', Fehler)
            CALL OlfCnuP (jout, 'y', 'dsig/dQdw', ' ', Fehler)
            CALL OlfComAdd (jout, 'T', 'trafo w->tof', Fehler)
         ELSEIF (iTra.eq.-1) THEN
            IF (coY.ne.'dsig/dQdw' .or. unY.ne.' ') THEN
               Fehler = ' call with dsig/dQdw[ ]'
               RETURN
               ENDIF
            CALL OlfCnuP (jout, 'x', 'w', 'meV', Fehler)
            CALL OlfCnuP (jout, 'y', 'S(q,w)', 'meV-1', Fehler)
            CALL OlfComAdd (jout, 'T', 'trafo tof->w', Fehler)
         ELSE
            Fehler = 'unexpected y-coordinate '
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

         IF (Fehler.ne.'&ff') RETURN

C  Information needed for conversion :
         iTrabs = iabs(iTra)
         IF (iTrabs.eq.1) THEN
            E0 = rOlfGG (j, 'E0', 'meV', Fehler)
            IF (E0.le.0.) THEN
               Fehler = 'E0 <= 0.'
               RETURN
               ENDIF
            halfmass = 5.2271  ! 0.5*mass_of_neutron [meV,millisec,m]
            fact     = halfmass * dsqrt (halfmass/E0)
               ! phase space factor ki/kf included 5jul91
            ENDIF

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            DO i = 1, n

               IF     (iTra.eq. 1) THEN  ! w -> tof
                  X(i) = Tau_of_W (X(i), E0)     !  t^2 = (m/2)/(E0+w)
                  dwdt = 2*fact/X(i)**4          ! |dw/dt| = 2*(m/2)/t^3; ki/kf
                  Y(i) = Y(i) * dwdt
                  D(i) = D(i) * dwdt

               ELSEIF (iTra.eq.-1) THEN  ! tof -> w
                  IF (X(i).le.0.) THEN
                     Fehler = 'time-of-flight <= 0.'
                     RETURN
                     ENDIF
                  dwdt = 2*fact/X(i)**4
                  X(i) = W_of_Tau (X(i), E0)
                  Y(i) = Y(i) / dwdt
                  D(i) = D(i) / dwdt

               ELSE
                  CALL Absturz ('TraAchseX', 'Option o.o.r.')
                  ENDIF

               ENDDO ! i

            CALL OlfPutXYD (jout, Kout, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (jout, Kout, Fehler)

         ENDDO ! lj

      END ! TraAchseX

C  ====================================================================
C  i70 / 2 :   convolutions
C  ====================================================================

C  --------------------------------------------------------------------
      SUBROUTINE TraSymm (nJList, JList, Fehler)
C  --------------------------------------------------------------------
            ! new version 19nov92
         ! Divide a spectrum into symmetric and antisymmetric part.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      DATA          q1 /.true./, q2 /.false./, tol /1.d-5/

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Common Questionary :
      q1 = qAskD (' Take symmetric part', intq(q1))
      q2 = qAskD (' Take antisymmetric part', intq(q2))
      qV = qAskD (' Set also x<0 side', intq(qV))
      tol= rAskDMu (' Tolerance for x-scale', tol, 1.d-14, 1.d0)

C  Loop over files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, .false., j1, nK, K1, Fehler)
         CALL OlfHeadDup (j, .false., j2, nK, K2, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (j1, '|', 'symmetric part', Fehler)
         CALL OlfComAdd (j2, '|', 'antisymmetric part', Fehler)
c96/7         CALL Append (tPar2(1), '_u')

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ! Find elastic channel :
            IF (n.le.3) THEN
               Fehler = 'spectrum too short'
               RETURN
               ENDIF
            IF (irSorted(X,n).ne.2) THEN
               Fehler = 'spectrum not sorted'
               RETURN
               ENDIF
            i0 = irPosOpt(X, n, 0.d0, 'n', n/2)
            IF (i0.le.2 .or. i0.ge.n-1) THEN
               Print '(a,2i5)', ' elastic channel, total channels :',
     *               i0,n
               Fehler = 'not enough channels in spectrum '//cl3(K)
               RETURN
               ENDIF
            dx = (X(i0+1) - X(i0-1))/2  ! typical step
            IF (dabs(X(i0)).gt.tol*dx) THEN
               Fehler = 'no entry x=0 in spectrum '//cl3(K)
               RETURN
               ENDIF

            ! Division in symmetric and antisymmetric parts :
            ! number of output channels :
               ! i0 - 1 negativ channels
               ! n - i0 positiv channels
            ng = max0 (i0, n-i0+1) ! even channels
            nu = min0 (i0, n-i0+1) ! odd  channels
            ! zero channel :
            X1 (1) = 0.
            Y1 (1) = Y(i0)
            D1 (1) = D(i0)
            X2 (1) = 0.
            Y2 (1) = 0.
            D2 (1) = 0.
            ! area covered by both even and odd channels :
            DO i = 2, nu
               ip = i0-1+i
               in = i0+1-i
               IF (X(ip)+X(in).gt.tol*X(ip)) THEN
                  Fehler = 'grid not symmetric'
                  RETURN
                  ENDIF
               X1 (i) = X(ip)
               Y1 (i) = ( Y(ip) + Y(in) ) / 2
               D1 (i) = dsqrt ( D(ip)**2 + D(in)**2 ) / 2
               X2 (i) = X(ip)
               Y2 (i) = ( Y(ip) - Y(in) ) / 2
               D2 (i) = D1(i)
               ENDDO
            ! rest of even channels :
            IF (i0 .gt. n-i0+1) THEN
               isi = -1 ! using the negativ side
            ELSE
               isi = +1
               ENDIF
            DO i = nu+1, ng
               ig = i0 + isi*(i-1)
               X1 (i) = isi * X(ig)
               Y1 (i) = Y(ig)
               D1 (i) = D(ig)
               ENDDO

            ! Double the spectra ? (6feb92)
            IF (qV) THEN
               Fehler = 'doubling spectra out of use'
               IF (0.eq.0) RETURN
               nC1old = nC1
               nC2old = nC2
               nC1 = 2*nC1 - 1
               nC2 = 2*nC2 - 1
               IF (nC1.gt.MC .or. nC2.gt.MC) THEN
                  Fehler = ' two-sided spectra became too long'
                  RETURN
                  ENDIF
               DO i = nC1, nC1old, -1 ! shift the x>0 side
                  X1(i) =   X1(i-nC1old+1)
                  Y1(i) =   Y1(i-nC1old+1)
                  D1(i) =   D1(i-nC1old+1)
                  ENDDO
               DO i = 1, nC1old-1     ! fill the x<0 side
                  X1(i) = - X1(nC1+1-i)
                  Y1(i) =   Y1(nC1+1-i)
                  D1(i) =   D1(nC1+1-i)
                  ENDDO
               DO i = nC2, nC2old, -1
                  X2(i) =   X2(i-nC2old+1)
                  Y2(i) =   Y2(i-nC2old+1)
                  D2(i) =   D2(i-nC2old+1)
                  ENDDO
               DO i = 1, nC2old-1
                  X2(i) = - X2(nC2+1-i)
                  Y2(i) = - Y2(nC2+1-i) ! antisymmetric filling
                  D2(i) =   D2(nC2+1-i)
                  ENDDO
               ENDIF

            ! Save the results :
            IF (q1) THEN
               CALL OlfCopZ   (j, j1, K, K, Fehler)
               CALL OlfPutXYD (j1, K, ng, X1, Y1, D1, Fehler)
               ENDIF
            IF (q2) THEN
               CALL OlfCopZ   (j, j2, K, K, Fehler)
               CALL OlfPutXYD (j2, K, nu, X2, Y2, D2, Fehler)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! K

         IF (q1) CALL OlfClos (j1, K, Fehler)
         IF (q2) CALL OlfClos (j2, K, Fehler)

         ENDDO ! lj

      END ! TraSymm

C  --------------------------------------------------------------------
      SUBROUTINE TraDouble (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! add negative part -X(n)..0 to a spectrum 0=X(1)..X(n)
             ! JWu Grenoble 28feb96

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      CHARACTER      LisDoc*40, aus*80
      DIMENSION      qList(MC)

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (j, 't', 'spectrum doubled (mirror at x=0)',
     *                   Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (irSorted(X,n).lt.2) THEN
               Fehler = ' x scale is not in ascending order'
               RETURN
               ENDIF
            IF (X(1).ne.0.) THEN
               Fehler = ' x scale does not start with 0'
               RETURN
               ENDIF

            nn = 2*n - 1
            IF (nn.gt.MC) THEN
               Fehler = ' too many points on x-scale'
               RETURN
               ENDIF

            DO i = n, 1, -1     ! shift old data
               X(n-1+i) = X(i)
               Y(n-1+i) = Y(i)
               D(n-1+i) = D(i)
               ENDDO
            DO i = 1, n-1       ! insert new data
               X(i) = -X(nn+1-i)
               Y(i) =  Y(nn+1-i)
               D(i) =  D(nn+1-i)
               ENDDO

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, nn, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (jout, nK, Fehler)

C  End loop files :
         ENDDO

      END ! TraDouble

C  ====================================================================
C  i70 / 3 :   Monte-Carlo convolution
C  ====================================================================

C  --------------------------------------------------------------------
      SUBROUTINE MC_Conv (nJList, JList, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'l_def.f'
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INTEGER       nJList, JList(*), KrK(MK)
      INTEGER       count, rate, maxc !Artem add for rand()
      REAL          rseed !Artem add for rand()
      CHARACTER*(*) Fehler
      CHARACTER*40  docG
      real, external :: rand !Artem add for rand()

      DATA  jr /-1/, nmc /-1/, jmc /0/

      IF (nJlist.eq.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      nmc = iAskDMu ('Number of Monte-Carlo runs', nmc, 0, 2**15)
      IF (nmc.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF
      jmc = iAskD ('Initialize random generator with (0=clock)', jmc)
      IF (jmc.eq.0) THEN
         call system_clock(count, rate, maxc) !Artem: add for rand()
         rseed = ( real(mod(count, 4194304)) + 0.5 ) !Artem: add for rand()
     *           / 4194304.0
         if (rseed <= 0.0) rseed = 0.5/4194304.0 !Artem: add for rand()
         ran = rand (rseed) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CCF()
      ELSE
         rseed = ( real(mod(abs(jmc), 4194304)) + 0.5 ) / 4194304.0 !Artem: add for rand()
         ran = rand (rseed) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CBF (jmc)
         ENDIF

      Print *, ' choose output grid :'
      CALL SetGridChoice (.true., docG, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO lj = 1, nJList
         jin = JList(lj)

         CALL OlfHeadDup (jin, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         jr = iAskDMu ('Resolution from file', jr, 0, MF)
         CALL GetK2K (jin, jr, nK, KrK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, '*', 'MC convolution', Fehler)

         CALL SetGridJ (jin, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         DO K = 1, nK

            CALL SetGridK (K, n1, Fehler) ! -> X1(1:n1)
            IF (Fehler.ne.'&ff') RETURN
            DO i = 1, n1
               Y1(i) = 0
               ENDDO

            ! scattering law :
            CALL OlfGetXY (jin, K, n, X, Y, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ! resolution function :
            CALL OlfGetXY (jr, KrK(K), n2, X2, Y2, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            CALL CoordBins (MC, n2, X2, X3, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            CALL CoordBins (MC, n, X, X4, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            CALL CoordBins (MC, n1, X1, D3, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ! probability distribution = normalized integral
            sum = 0
            DO i = 1, n2
               IF (Y2(i).lt.0) THEN
                  Fehler = ' negative entries in the resolution'
                  RETURN
                  ENDIF
               sum = sum + Y2(i) * (X3(i+1)-X3(i))
               Y3(i) = sum
               ENDDO
            IF (sum.le.0) THEN
               Fehler = 'vanishing resolution'
               RETURN
               ENDIF
            DO i = 1, n2
               Y3(i) = Y3(i) / sum
               ENDDO

            sum = 0
            DO i = 1, n
               IF (Y(i).lt.0) THEN
                  Fehler = ' negative entries in the scattering law'
                  RETURN
                  ENDIF
               sum = sum + Y(i) * (X4(i+1)-X4(i))
               Y4(i) = sum
               ENDDO
            IF (sum.le.0) THEN
               Fehler = 'vanishing scattering law'
               RETURN
               ENDIF
            DO i = 1, n
               Y4(i) = Y4(i) / sum
               ENDDO

            ! The Monte-Carlo loop:
            DO imc = 1, nmc
               ! random selection from resolution distribution:
               ran  = dble(rand(0.0)) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
               iin  = irPos (Y3, n2, ran, 'l')
               ran  = dble(rand(0.0)) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
               xin  = X3(iin)*ran + X(iin+1)*(ran-1)

               ! random selection from scattering law:
               ran  = dble(rand(0.0)) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
               isc  = irPos (Y4, n, ran, 'l')
               ran  = dble(rand(0.0)) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
               xsc  = X4(isc)*ran + X4(isc+1)*(ran-1)

               iout = irPos (D3, n1, xin+xsc, 'l')
               Y1(iout) = Y1(iout) + 1
               ENDDO ! Monte-Carlo

            DO i = 1, n1
               D1(i) = dsqrt(Y1(i))
               Y1(i) = Y1(i) / nmc
               D1(i) = D1(i) / nmc
               ENDDO

            CALL OlfCopZ (jin, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n1, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO ! K

         CALL OlfClos (jout, nK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! files

      END ! MC_Conv
