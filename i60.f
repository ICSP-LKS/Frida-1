C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i60   :     curves and fits
C
C  ====================================================================

C     Contents :
C        1.  Functions :
C               ExtTabXZ
C        2.  Special functions :
C               Goetze_, HavNeg, Voigt
C        3.  Convolution and Par-File :
C               jCuConvAsk, CuConvPrep, CuConvVal, CuFPFPrep, CuNice
C        4.  Fit :
C               CuFitExe, CuFitFunction, CuFitMonit
C        5.  Interface :
C               CuSetPar, CuSetAux, CuSetFit, CuGetPar, CuCreate, CuFitCall

C     Aenderungsverzeichnis :
C        JWu  3feb95 : Softening for convolution
C        JWu 24jan95 : function from external table
C        JWu  9dec92 : Convolution
C        JWu  6jun91 : Parameterzugriff neugestaltet; freie und feste Parameter
C        JWu  7mar91 : Plots und einfache Fits funktionieren
C        JWu 23feb91 : Speicherung von Curves als Files beschlossen
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  ====================================================================
C  i60 / 1 :   Special functions
C  ====================================================================

      SUBROUTINE ExtTabXZ (Table, XE, iiE, ifE, zE,
     *                     xoff, xmul, yoff, ymul, YE, Fehler)
C     --------------------------------------------------------
         ! interpolate function from externally supplied xy-table
            ! JWu 24/25jan95

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'

      CHARACTER*(*) Table, Fehler
      CHARACTER*80  Path, fileT, TableOld
      CHARACTER*40  DirNum

      DIMENSION   X1(MC), Y1(MC), X2(MC), Y2(MC), XE(*), YE(*), ZT(MK)

      IF (Fehler.ne.'&ff') RETURN

C  Search Table :
      CALL ExeML ('\p dir-num', DirNum)
      CALL Compose2 (Path, DirNum, Table)
      nF = MemBlockInq ('nF')
      IF (qiinside(jT,1,nF)) THEN ! try the same jT as before
         CALL tOlfG (jT, 'fil', fileT, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (fileT.eq.Table .or. ! PROVISORISCHHHHHHHHHHHHH ******************
     *     fileT(1:lenU(fileT)).eq.Path(1:lenU(Path))) GOTO 11 ! found
         ENDIF
      DO jT = 1, nF ! try all files
         CALL tOlfG (jT, 'fil', fileT, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (fileT.eq.Table) GOTO 11 ! found
         ENDDO

C  Load new table :
      Print *, ' .. loading numeric table'
      CALL FileLoad ('&int '//Path, Fehler) ! load as protected file
      IF (Fehler.ne.'&ff') THEN
         Print *, Fehler
         Fehler = 'Cannot load numeric table'
         RETURN
         ENDIF
      jT = nF + 1

      ! trivial check :
      CALL tOlfG (jT, 'fil', fileT, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      IF (fileT.ne.Table .and.
     *    fileT(1:lenU(fileT)).ne.Path(1:lenU(Path))) THEN
                             ! letzteres PROVISORISCH
         Print *, 'load table :', Table
         Print *, 'load path  :', Path
         Print *, 'internal   :', fileT
         Fehler = 'PROGR ERR/ num-tab file-name inconsistent'
         RETURN
         ENDIF

C  Nontrivial checks :
      CALL OlfGet1ZofK (jT, 1, nK, ZT, Fehler)
      DO K = 1, nK
         CALL OlfGetXY (jT, K, n1, X1, Y1, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (n1.lt.2) THEN
            Fehler = 'Numeric table contains spectrum with n<2'
            RETURN
         ELSEIF (irSorted(X1,n1).ne.2) THEN
            Fehler = 'Numeric tables must be sorted'
            RETURN
            ENDIF
         ENDDO

 11   CONTINUE ! so oder so - Tabelle geladen.

C  Interpolate in z :
      K = irPosOpt (ZT, nK, zE, 'r', K)
      IF (K.lt.1)  K = 1
      IF (K.ge.nK) K = nK-1

C  Load Y(X) for z(K) and z(K+1) :
      IF (Table.eq.TableOld .and. K.eq.Kold) THEN
         ! reuse old X1.., Y1.. % Hoffnung auf Beschleunigung 19mai95
      ELSE
         CALL OlfGetXY (jT, K, n1, X1, Y1, Fehler)
         CALL OlfGetXY (jT, K+1, n2, X2, Y2, Fehler)
         IF (Fehler.ne.'&ff') RETURN ! extrapolation -> 0.
         ENDIF

C  Interpolate in x :
      DO iE = iiE, ifE
         xarg = ( XE(iE) - xoff ) * xmul ! logischer waere div statt mul
                                         ! FK diese Aussage ist quatsch!!
                                         !
                                         ! external KWW table is wrong!
                                         ! this leads to a FT problem
                                         ! KWW tau values shown are a
                                         ! factor e[meV]/hbar too large in
                                         ! S(q,w). The value used for the
                                         ! fit is by the inverse of this factor
                                         ! lower. Hence instead of tau
                                         ! tau*e/hbar shall be used.(quick fix
                                         ! only for kww_sqw routines.
         i1 = irPosOpt (X1, n1, xarg, 'r', i1)
         IF     (i1.lt.1) THEN
            ye1 = Y1(1)
         ELSEIF (i1.gt.n1-1) THEN
            ye1 = Y1(n1)
         ELSE
            CALL LinIntPol1 (xarg, X1(i1), X1(i1+1),
     *                       Y1(i1), Y1(i1+1), ye1)
            ENDIF
         i2 = irPosOpt (X2, n2, xarg, 'r', i1)
         IF     (i2.lt.1) THEN
            ye2 = Y2(1)
         ELSEIF (i1.gt.n2-1) THEN
            ye2 = Y2(n2)
         ELSE
            CALL LinIntPol1 (xarg, X2(i2), X2(i2+1),
     *                       Y2(i2), Y2(i2+1), ye2)
            ENDIF

         YE(iE) = yoff + ymul * (ye1 + (zE-ZT(K)) /
     *                 (ZT(K+1)-ZT(K)) * (ye2-ye1))

         ENDDO

      TableOld = Table
      Kold     = K

      END ! ExtTabXZ

C  ====================================================================
C  i6.3.  Convolution
C  ====================================================================

      INTEGER FUNCTION jCuConvAsk (qCurv, qConv)
C     ------------------------------------------
            ! JWu 9dec92
         ! a piece of dialogue for CuFit, IdaPlot, OrgGrid
      IMPLICIT LOGICAL (q)

      IF (qCurv .and. qConv) THEN
         j2 = iAskD (' Convolute with file (0=no convolution)', j2)
         jCuConvAsk = j2
      ELSE
         jCuConvAsk = 0
         ENDIF
      END ! jCuConvAsk

      SUBROUTINE CuConvPrep (X1, n1, j2, K2, Fehler)
C     ----------------------------------------------
            ! JWu 9dec92
         ! prepare convolution :
         ! get spectrum of conv-file, check grid, transfer to COMMON

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      PARAMETER        (MC2=2*MC)
      DIMENSION         Ydummy(MC2), Wdummy(MC2)

      DIMENSION         X1(*)
      CHARACTER         Fehler*(*), cl6*6

      DATA          twopi /6.2831853/

C  Communication between CuConvPrep and CuConvVal :
      COMMON / CuFuCo / qConv, qEq, X3(MC), n3, n2, X2(MC), Y2(MC)

      DATA        tol / 1.d-5 /

      IF (n1.le.0 .or. n1.gt.MC)
     *    CALL Absturz ('CuConvPrep', 'n1 o.o.r.')

      qConv = (j2.ne.0) ! convolute or not ?

C  If no convolution :
      IF (.not.qConv) THEN
         ! copy x-scale to COMMON :
         n3 = n1
         DO i = 1, n1
            X3(i) = X1(i)
            ENDDO
         RETURN
         ENDIF

C  Examine data-file :
      CALL CheckScale (n1, X1, tol, qEq1, dX1)
      IF (qEq1 .and. dX1.le.0.) THEN
         Fehler = '1st file not in ascending order'
         GOTO 900
         ENDIF

C  Examine conv-file :
      CALL OlfGetXY  (j2, K2, n2, X2, Y2, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      IF (irSorted(X2,n2).ne.2) THEN
         Fehler = 'conv-file not sorted'
         GOTO 900
         ENDIF
      IF (X2(1).gt.0.d0) THEN
         Fehler = 'conv-file doesn''t contain 0.'
         GOTO 900
         ENDIF
      i02 = irPosOpt(X2, n2, 0.d0, 'n', n2/2)
      IF (i02.le.2 .or. i02.ge.n2-1) THEN
         Fehler = 'not enough positiv or negativ channels in conv-file'
         GOTO 900
         ENDIF
      sx2 = (X2(i02+1) - X2(i02-1))/2  ! typical step
 212  CONTINUE
      IF (dabs(X2(i02)).gt.tol*sx2) THEN
         CALL Gong (9)
         Print *, ' grid of conv-file doesn''t contain an entry 0.'
         rin = rAskDMu (' Tolerance (0=quit)', tol, 0.d0, 1.d0)
         IF (rin.le.1.d-12) GOTO 901
         tol = rin
         GOTO 212
         ENDIF

      CALL CheckScale (n2, X2, tol, qEq2, dX2)

C  Normalization :
      IF (qEq2) THEN
         y2sum = rSum (Y2, 1, n2, 1)
         IF (y2sum.le.0.) THEN
            Fehler = 'conv-file data are zero or negative'
            GOTO 900
            ENDIF
         DO i = 1, n2
            Y2(i) = Y2(i) / y2sum
            ENDDO
      ELSE
         Print *, ' WARNING/ Conv File not equidist/ not tested'
         CALL Gong (2)
         y2sum = rIntegralXY (X2, Y2, n2)
         DO i = 1, n2
            Y2(i) = Y2(i) / y2sum
            ENDDO
         ENDIF

C  Both files together :
      qEq = qEq1 .and. qEq2 .and. qEqEps(dX1,dX2)
      IF (qEq) THEN
         n3 = n1 + n2 - 1
         dX3 = (dX1 + dX2) / 2
         IF (n3.gt.MC) THEN
            Fehler = 'spectra too long'
            GOTO 900
            ENDIF
         x3i = X1(1) - (n2-i02) * dX3
         DO i = 1, n3
            X3(i) = x3i + (i-1) * dX3
            ENDDO

      ELSE
         n3 = n1
         DO i = 1, n1
            X3(i) = X1(i)
            ENDDO
         ENDIF

      RETURN

C  Errors :
 900  CONTINUE
      CALL Insert (Fehler, 1, 'Convolution impossible/ ')
      RETURN
 901  CONTINUE
      Fehler = ' '
      RETURN

      END ! CuConvPrep

      SUBROUTINE CuConvVal (ifc, P, Y, Fehler)
C     ----------------------------------------
            ! JWu 10dec92
         ! Convolute CuVal with conv-file.

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      DIMENSION         P(*), Y(*), X1(MC), Y1(MC)
      CHARACTER         Fehler*(*)

      COMMON / CuFuCo / qConv, qEq, X3(MC), n3, n2, X2(MC), Y2(MC)

      IF (qConv) THEN
         IF (qEq) THEN
            CALL CuFuVal (ifc, P, X3, Y1, n3, Fehler)
            DO i = 1, n3-n2+1
               ! convolution : see E2,33
               Y(i) = 0.
               DO ii = 1, n2
                  Y(i) = Y(i) + Y2(ii) * Y1(i-ii+n2)
                  ENDDO
               ENDDO

         ELSE ! nonequidistant convolution 11dec92.
            DO i = 1, n3
               Y(i) = 0.
               ENDDO
            DO ii = 1, n2
               DO i = 1, n3
                  X1(i) = X3(i) - X2(ii)
                  ENDDO
               CALL CuFuVal (ifc, P, X1, Y1, n3, Fehler)
               DO i = 1, n3
                  Y(i) = Y(i) + Y1(i)*Y2(ii)
                  ENDDO
               ENDDO

            ENDIF

      ELSE
         CALL CuFuVal (ifc, P, X3, Y, n3, Fehler)

         ENDIF ! qConv or not

      END ! CuConvVal

      SUBROUTINE CuFPFPrep (jFPF, KFPF, Fehler)
C     -----------------------------------------
         ! prepare par-file for use in CuFuVal

      IMPLICIT REAL*8 (a-h,o-p,r-z)

      INCLUDE 'i_dim.f'

      CHARACTER Fehler*(*)

      ! used only here and in CuFuVal :
      COMMON / FitParFil / XFPF(MC), YFPF(MC), nFPF

      IF (jFPF.eq.0) RETURN

      CALL OlfGetXY (jFPF, KFPF, nFPF, XFPF, YFPF, Fehler)
        IF (Fehler.ne.'&ff') RETURN

      ! no checks whatsoever

      END ! CuFPFPrep

      SUBROUTINE CuNice (ifu, Y, j2, K2K, j3, K3K,
     *                   xmi, xma, igfx, gfpx, liloX,
     *                   ymi, yma, igfy, gfpy, liloY,
     *                   nX, relY,
     *                   MC, nC1, X1, Y1, D1, Fehler)
C     -----------------------------------------------
            ! JWu 25aug93
         ! evaluate curve on a nice x-grid
            ! verdiente, allgemeiner programmiert
            ! zu werden (z.B. f"ur Integration).

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      CHARACTER Fehler*(*)
      DIMENSION Y(*), X1(*), Y1(*), D1(*)

      ! Set maximum y step, as given by input parameter relY :
      IF (liloY.eq.1) THEN
         IF (yma.lt.ymi .or. ymi.lt.0.) THEN
            Fehler = 'CuNice/ log y-range is bad'
            RETURN
            ENDIF
            dy = (yma / ymi) ** relY
      ELSE
         IF (yma.lt.ymi) THEN
            Fehler = 'CuNice/ lin y-range is bad'
            RETURN
            ENDIF
         dy = (yma-ymi) * relY
         ENDIF

      ! Set minimal x step :
      IF (liloX.eq.1) THEN
         IF (xma.lt.xmi .or. xmi.lt.0.) THEN
            Fehler = 'CuNice/ log x-range is bad'
            RETURN
            ENDIF
            dx = (xma / xmi) ** relY / 1024
      ELSE
         IF (xma.lt.xmi) THEN
            Fehler = 'CuNice/ lin x-range is bad'
            RETURN
            ENDIF
         dx = (xma-xmi) * relY / 1024
         ENDIF

      ! # interpol points :
      IF (nX.le.2) THEN
         Fehler = 'Calculate curve/ '//
     *            'less than 2 interpolation points demanded'
         RETURN
         ENDIF
      nC1 = nX

      ! Set regular x-grid :
      IF (liloX.eq.0) THEN
         xstep = (xma-xmi)/(nC1-1)
         DO i = 1, nC1
            x1i = xmi + (i-1)*xstep
            CALL FuInv (igfx, X1(i), dummy, x1i, 0.d0, gfpx, 0.d0)
            ENDDO
      ELSE
         IF (xmi.le.0. .or. xma.le.xmi) THEN
            Fehler = 'Negative plotrange <-> log mode'
            RETURN
            ENDIF
         xstep = (xma/xmi) ** (1./(nC1-1))
         DO i = 1, nC1
            x1i = xmi * xstep ** (i-1)
            CALL FuInv (igfx, X1(i), dummy, x1i, 0.d0, gfpx, 0.d0)
            ENDDO
         ENDIF

C  Loop : nice and nicer curve.
      iInsert = 0
 70   CONTINUE

C  Calculate Y1 on grid X1 :
         CALL CuFPFPrep (j3, K3K, Fehler)
         CALL CuConvPrep (X1, nC1, j2, K2K, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL CuConvVal (ifu, Y, Y1, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Check stepwidth in Y :
         i  = 1 ! scan Y1
         nn = 0 ! new x-grid -> D1
C  Subloop : check step from i to i+1 :
 75      CONTINUE
C  Copy X1 to D1 :
            nn = nn + 1
            D1(nn) = X1(i)

            IF (( (liloY.eq.0 .and. dabs(       Y1(i+1)-Y1(i) ).gt.dy)
     *    .or. (liloY.eq.1 .and. dabs(dquot0(Y1(i+1),Y1(i))).gt.dy))
     *    .and. ( (liloX.eq.0 .and. dabs(       X1(i+1)-X1(i) ).ge.dx)
     *    .or. (liloX.eq.1 .and. dabs(dquot0(X1(i+1),X1(i))).ge.dx))
     *          ) THEN
C  Step in Y1 is too large (and step in X1 not too small) :
               IF (nn+5.gt.MC) THEN ! Fehlermeldung -> Abbruch ??
                  Print *, ' .. dense oscillations ?'
                  GOTO 91
                  ENDIF
C  Insert new point in D1 :
               dx1 = (X1(i+1) - X1(i)) ! old x-step
               IF (i.gt.1) THEN ! insert in preceeding interval
                  D1(nn+1) = D1(nn)
                  D1(nn)   = D1(nn) - dx1/3
                  nn       = nn+1
                  ENDIF
               IF (nn+2.gt.MC) GOTO 91
               D1(nn+1) = D1(nn) +  dx1/3
               D1(nn+2) = D1(nn) +2*dx1/3
               nn       = nn+2
               IF (i.lt.nC1-1) THEN ! insert in following interval
                  IF (nn+2.gt.MC) GOTO 91
                  D1(nn+1) = X1(i+1)
                  D1(nn+2) = X1(i+1) + dx1/3
                  nn       = nn+2
                  i        = i+1 ! don't check next interval
                  ENDIF
               ENDIF ! insertion
            IF (i.lt.nC1-1) THEN
               i = i + 1
               GOTO 75
               ENDIF
            IF (nn+1.gt.MC) GOTO 91
            D1(nn+1) = X1(nC1)
            nn       = nn+1
C  End subloop.

C  If some points have been inserted, go back to recalculate Y1 :
         IF (nn.gt.nC1 .and. nn.le.MC) THEN
            nC1 = nn
            CALL rCopy (X1, 1, nC1, 1, D1, 1, 1)
            IF (iInsert.lt.10) THEN
               Print *, ' .. improving grid for curve'
               iInsert = iInsert + 1
               GOTO 70
               ENDIF
            Print *, ' .. escape from improvements'
            ENDIF
C  End loop.
 91   CONTINUE
      CALL rSet (D1, 1, nC1, 1, 0.d0)

      END ! CuNice

C  ====================================================================
C  i6.4.  Curve/ Fit
C  ====================================================================

      BLOCK DATA CuFitPreset
C     ----------------------

      IMPLICIT NONE
      REAL*8          tolFit, stpFit, etaFit1, etaFitN
      INTEGER         mclFit, imoFit, jmoFit, idiFit
      LOGICAL         qRPFit

      COMMON / FitSet / tolFit, stpFit, etaFit1, etaFitN, mclFit,
     *                  qRPFit, imoFit, jmoFit, idiFit

      DATA     tolFit / 1.d-3 /, stpFit / 1.d3 /,
     *         etaFit1/ 0.d0  /, etaFitN/ .5d0  /, mclFit /  100  /,
     *         qRPFit / .true. /, imoFit / -1 /,
     *         jmoFit / 1 /, idiFit / 0 /

      END ! CuFitPreset

      SUBROUTINE CuFitExe (jd, jc, K, jconv, Kconv,
     *                     jFPF, KFPF, SetupLine, Fehler)
C     -----------------------------------------------------------------
            ! JWu mar91
         ! fit data (jd,K) by curve (jc,K) convoluted with (jconv,Kconv)
         ! eventually further data are taken from (jFPF,KFPF) (jul91/aug93).

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  SetupLine, Fehler

C  Data file and curve file :
      CHARACTER*40   File, h1, h2, Un
      DIMENSION      X1(MC), Y1(MC), D1(MC), D(MC)

C  Arrays and matrices for use in fitroutine :
      DIMENSION      Delta(MC), FitPar(MFP), Fjac(MC,MFP),
     *               Sing(MP), Fsing(MFP, MFP), CoVar(MFP)

C  Work space and subroutines for E04FCF :
      PARAMETER     (MIwork=1, MPwork=8,
     *            MRwork=6*MPwork+MC*MPwork+2*MC+MPwork*(MPwork-1)/2)
      DIMENSION      Iwork(MIwork), Rwork(MRwork)
      EXTERNAL       CuFitFunction, CuFitMonit, CuFitFunction_MP
      REAL*8         diag(MFP), qtf(MFP), wa1(MFP),
     *               wa2(MFP), wa3(MFP), wa4(MC)
      REAL*8         ftol, xtol, gtol, epsfcn, factor
      INTEGER        ipvt(MFP), mode, info !Artem: Add additional parameters for lmdif

C  Communication :
      CHARACTER*80   aus, FitFehler

C  Communication with CuFitFunction, CuFitMonit :
      COMMON / CuFiCo / XFi(MC), YFi(MC), Weight(MC), ifc, qFiLoY,
     *          jFitMon, nCuP, qFix(MFP), FixPar(MFP), ScaPar(MFP)
      COMMON / CuFiFe / FitFehler

C  and with CuSetFit :
      COMMON / FitSet / tolFit, stpFit, etaFit1, etaFitN, mclFit,
     *                  qRPFit, imoFit, jmoFit, idiFit

C  Load data :
      CALL OlfGetXYD  (jd, K, nC, XFi, YFi, D, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Load curve :
      CALL OlfGetXYD  (jc, K, nC1, X1, Y1, D1, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      nK1    = iOlfG (jc, '#spectra', Fehler)
      ifc    = iOlfG (jc, 'fu#', Fehler)
      nCuP   = iOlfG (jc, '#fit-par', Fehler)
      iFixed = iOlfG (jc, '@fixed', Fehler)
      qFiLoY = qOlfG (jc, '?weight-log-y', Fehler)
      qFiErY = qOlfG (jc, '?weight-err-y', Fehler)
      qFiStX = qOlfG (jc, '?weight-stp-x', Fehler)
      fitRi  = rOlfG (jc, 'fit-i', Un, Fehler)
      fitRf  = rOlfG (jc, 'fit-f', Un, Fehler)
      CALL     tOlfG (jc, 'fil', file, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      ! write setting into one line :
      CALL NiceNum (fitRi, h1, ih1)
      CALL NiceNum (fitRf, h2, ih2)
      CALL Compose2 (SetupLine, 'fit to '//File,
     *                  ' over '//h1(1:ih1)//'..'//h2(1:ih2) )
      IF (qFiErY) CALL Insert (SetupLine, 1, 'E') ! weight dy
      IF (qFiStX) CALL Insert (SetupLine, 1, 'S') ! weight dx
      IF (qFiLoY) CALL Insert (SetupLine, 1, 'L') ! log y

      ! check compatibility :
c      IF (nK1.ne.1 .and. z1.ne.z) THEN ! -Z-MODIF-Z-
c         IF (.not.qAsk('Data and curve have different z - continue ?')) THEN
c             Fehler = ' '
c           RETURN
c            ENDIF
c         ENDIF

C  Limits for fit-range (rewritten 30dec92 for unsorted spectra 24may98) :
      IF (fitRi.ne.0. .or. fitRf.ne.0.) THEN
         nCFit = 0
         DO i = 1, nC
            IF (XFi(i).ge.fitRi .and. XFi(i).le.fitRf) THEN
               nCFit = nCFit + 1
               XFi(nCFit) = XFi(i)
               YFi(nCFit) = YFi(i)
               D  (nCFit) = D  (i)
               ENDIF
            ENDDO
         IF (nCFit.le.0) THEN
            CALL NiceNum (fitRf, h2, ih2)
            Fehler = ' No points in fitrange '//h1(1:ih1)//
     *               ' .. '//h2(1:ih2)
            RETURN
         ELSEIF (nCFit.le.3) THEN
            Print *, ' WARNING/ fitting only '//ch1(nCFit)//
     *               ' points'
            ENDIF

      ELSE
         nCFit = nC
         ENDIF

C  Prepare convolution :
      CALL CuConvPrep (XFi, nCFit, jconv, Kconv, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Prepare par-file :
      CALL CuFPFPrep (jFPF, KFPF, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Set weight :
      IF (qFiStX) THEN
         ! weighting with x-stepwidth :
         Weight (1) = (XFi(2)-XFi(1)) / 2
         DO i = 2, nCFit-1
            Weight (i) = (XFi(i+1)-XFi(i-1)) / 2
            ENDDO
         Weight(nCFit) = (XFi(nCFit)-XFi(nCFit-1)) / 2
         DO i = 1, nCFit
            IF (Weight(i).lt.0.) THEN
               Fehler = ' weighting/ data not sorted'
               RETURN
               ENDIF
            Weight(i) = dsqrt(Weight(i)) ! residuals are Weight*(y-f)
            ENDDO
      ELSE
         DO i = 1, nCFit
            Weight (i) = 1.
            ENDDO
         ENDIF

      IF (qFiErY) THEN
         ! Weighting with error bars :
         nD0 = 0
         DO i = 1, nCFit
            IF     (D(i).lt.1.d-10) THEN
               Weight(i) = 0.
               nD0 = nD0 + 1
            ELSEIF (D(i).eq.0.) THEN
               Print *, ' There are entries with error<0.'
               Fehler = ' Weighting with y-error not possible'
               RETURN
            ELSE
               Weight(i) = Weight(i) / D(i) ! corrected 3jul91
               ENDIF
            ENDDO
         IF (nD0.ge.nCFit/3) THEN
            CALL Compose2 (aus, ' WARNING/ there are '//cl6(nD0),
     *         ' channels with error = 0')
            Print *, aus
         ENDIF
      ELSE
         nD0 = 0
         ENDIF

C  Initialize fitparameter :
      IF (nCuP.le.0 .or. nCuP.gt.MFP) THEN
         Fehler = ' nCuP o.o.r. : = '//cl3(nCuP)
         RETURN
         ENDIF
      nFitP = 0
      DO i = 1, nCuP
         qFix(i) = qBitGet (iFixed, i)
         IF (.not.qFix(i)) THEN
            nFitP = nFitP + 1
            IF     (.not.qRPFit) THEN
               ScaPar(nFitP) = 1.
               FitPar(nFitP) = Y1(i)
            ELSEIF (dabs(Y1(i)).ge.1.d-8) THEN
               ! reduce parameters to O(1) 30dec92 :
               ScaPar(nFitP) = Y1(i)
               FitPar(nFitP) = 1.
            ELSE
               ScaPar(nFitP) = 1.d-4
               FitPar(nFitP) = Y1(i) / ScaPar(nFitP)
               ENDIF
         ELSE
            FixPar(i) = Y1(i)
               ! Nicht sehr logisch : FitPar hat lueckenlose eigene Nummern,
               ! waehrend FixPar genau wie Y1 numeriert wird.
               ! Fuer FitPar ist die dichtere Packung fuer NAG erforderlich.
            ENDIF
         ENDDO
      IF (nFitP.eq.0) THEN
         Fehler = ' All parameters are fixed'
         RETURN
         ENDIF
      nFree = nCFit - nFitP - nD0
      IF (nFree.le.0) THEN
         Fehler = ' There are not more datapoints than free parameters'
         RETURN
         ENDIF

C  Prepare work space :
      nIwork = 1
      IF (nFitP.gt.1) THEN
         nRwork = nFitP*(6+nFitP/2) + nCFit*(nFitP+2)
      ELSE
         nRwork = 7+3*nCFit
         ENDIF
      IF (nRwork.gt.MRwork) THEN
         Print *, 'nFitPP nCFit nRwork MRwork', nFitP,
     *            nCFit, nRwork, MRwork
         Fehler = 'Preparing for E04FCF : nRwork too large'
         RETURN
         ENDIF

C  Fits :
      MaxCall = max0 (MaxCall, mclFit*nFree)
      nCall = MaxCall
      jFitMon = jmoFit ! COMMON -> COMMON
      IF (nFitP.eq.1) THEN
         etaFit = etaFit1
      ELSE
         etaFit = etaFitN
         ENDIF
 10   CONTINUE
      ifail = 0
      ftol = 0.d0
	  gtol = 0.d0
	  epsfcn = 0.d0
	  mode = 1
	  factor = 100.d0
	  info = 0 !Artem: Add additional parameters for lmdif
      IF (imoFit.gt.1) Print *
C      CALL E04FCF (nCFit, nFitP, CuFitFunction, CuFitMonit,
C     *             imoFit, MaxCall, etaFit, tolFit, stpFit, FitPar,
C     *             Sum, Delta, Fjac, MC, Sing, Fsing, MFP,
C     *             nIter, nCalls, Iwork, nIwork, Rwork, nRwork, ifail)
      call lmdif(CuFitFunction_MP, nCFit, nFitP, FitPar, Delta, ftol,
     *          tolFit, gtol, MaxCall,
     *          epsfcn, diag, mode, factor, imoFit, info, nIter, Fjac,
     *          nCFit, ipvt, qtf, wa1, wa2, wa3, wa4) !Artem: Replace with lmdif from minpack/lmdif.f
      Sum = ENORM(nCFit, Delta) !Artem: Calculation of squared sum
      Sum = Sum*Sum
      IF     (ifail.eq.0) THEN
         h1 = 'success'
         MaxCall = max0 (MaxCall*86/100, MaxCall-20, mclFit*nFree)
      ELSEIF (ifail.eq.1) THEN
         Fehler = cr3(K)//' : very bad error in E04FCF: ifail = 1'
         RETURN
      ELSEIF (ifail.eq.2) THEN
         Print *, cr3(K)//' : fit not yet converged'
         aus = ' How many more calls'
         ia = iAskD (aus, mclFit*nFree)
         IF (ia.gt.0) THEN
            nCall = ia
            MaxCall = MaxCall + ia
            GOTO 10
         ELSE
            h1 = 'tired'
            ENDIF
      ELSEIF (ifail.eq.3) THEN
         h1 = 'no conv'
      ELSEIF (ifail.eq.4) THEN
         h1 = 'trapped'
      ELSEIF (ifail.eq.-1685) THEN
         Fehler = 'fit stopped/ '//FitFehler
         RETURN
      ELSE
         Fehler = 'unforeseen ifail ='//cr3(ifail)
         RETURN
         ENDIF

C  Show result (31dec92) :
      IF (nFitP.gt.4) THEN
         nFPshow = min0 (nFitP, 5)
         Print '(i3,a,g16.7,a,a8,5(g11.3))',
     *      K, '(', Sum/nFree, ') ', h1,
     *      (ScaPar(ifp)*FitPar(ifp), ifp=1,nFPshow) !Artem: g9.2 -> g16.7 increase output precision
      ELSE
         Print '(i3,a,g16.7,a,a8,3(g14.6))',
     *      K, '(', Sum/nFree, ') ', h1,
     *      (ScaPar(ifp)*FitPar(ifp), ifp=1,nFitP)  !Artem: g9.2 -> g16.7 increase output precision
         ENDIF

C  Estimate errors (20nov91, see C2,33) :
      CALL CALCJAC_C(CuFitFunction_MP, nCFit, nFitP, FitPar, Fjac) !Artem: calculation of Jacobian Fjac(M,N) of a vector function CuFitFunction_MP at point FitPar(1:N)
      job   = 0 ! calculate diagonal elements
      ifail = 1 ! silent exit

C      CALL E04YCF (job, nCFit, nFitP, Sum, Sing, FSing, MFP,
C     *             CoVar, Rwork, ifail)

      IF (qFiErY) THEN
         call E04YCF_weighed_local(job, nCFit, nFitP, Delta, Fjac,
     *     CoVar, ipvt, wa1, ifail)  !Artem: Replace with a self-made subroutine from lnag_local.f
      ELSE
        call E04YCF_local(job, nCFit, nFitP, Delta, Fjac, CoVar, ipvt,
     * wa1, ifail)  !Artem: Replace with a self-made subroutine from lnag_local.f
      ENDIF
      IF     (ifail.eq.1) THEN
         Fehler = ' Bad parameters in E04YCF'
         RETURN
      ELSEIF (ifail.ge.2) THEN
         IF (nFitP+2-ifail.eq.1) THEN
            Print *,' One singular value of the Jacobian is degenerate'
         ELSE
            CALL Compose2 (aus, cl3(nFitP+2-ifail),
     *         ' singular values of the Jacobian are degenerate')
            Print *, aus
            ENDIF
         ENDIF

      IF (idiFit.eq.1) THEN
         Print *, '  >  Singular values of the Jacobian :'
         Print '(a,7g10.2)', '  >> ', (Sing(i), i=1,nFitP)
         ENDIF

C  Save results (curve) :
      nFitP = 0
      DO i = 1, nCuP
         IF (.not.qFix(i)) THEN
            nFitP = nFitP + 1
            Y1(i) = ScaPar(nFitP) * FitPar(nFitP)
C s. busch changed 2010-JAN-12
            ! D1(i) = dabs(ScaPar(nFitP)) * dsqrt (2*Sum/nFree *CoVar(nFitP))
            D1(i) = dabs(ScaPar(nFitP)) * dsqrt (CoVar(nFitP))
         ELSE
            ! of course FixP is unchanged and has not to be returned
            D1(i) = 0.
            ENDIF
         ENDDO
      X1(1) = Sum / nFree  !  since 20nov91

      CALL OlfPutXYD (jc, K, nC1, X1, Y1, D1, Fehler)
         ! The file jc may become a mixture of fitted
         ! curves and curves that are given by guess only.
         ! However, the state of each spectrum can be seen
         ! from X1(1) = Sum
      IF (Fehler.ne.'&ff') RETURN

      END ! CuFitExe

      SUBROUTINE CuFitFunction (iflag, nC, nFitP, FitPar, Delta,
     *                          Iwork, nIwork, Rwork, nRwork)
C     ----------------------------------------------------------
         ! called by E04FCF in CuFit
         ! calculates the residuals Delta(x) = Curve(x) - Dat(x)

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      DIMENSION    FitPar(nFitP), Delta(nC), Y2(MC),
     *             Iwork(nIwork), Rwork(nRwork), CuPar(MFP)
      CHARACTER    FitFehler*80

      COMMON / CuFiCo / XFi(MC), YFi(MC), Weight(MC), ifc, qFiLoY,
     *          jFitMon, nCuP, qFix(MFP), FixPar(MFP), ScaPar(MFP)
      COMMON / CuFiFe / FitFehler

C  Transfer fit_parameter -> curve_parameter :
      iFitP = 0
      DO i = 1, nCuP
         IF (qFix(i)) THEN
            CuPar(i) = FixPar(i)
         ELSE
            iFitP = iFitP + 1
            CuPar(i) = ScaPar(iFitP) * FitPar(iFitP)
            ENDIF
         ENDDO
      IF (iFitP.ne.nFitP) CALL Absturz ('CuFitFunction',
     *   'nFitP not reproducable')

C  Calculate curve :
      FitFehler = '&ff'
      CALL CuConvVal (ifc, CuPar, Y2, FitFehler)
      IF (FitFehler.ne.'&ff') THEN
         iflag = -1685
         RETURN
         ENDIF

C  Calculate the residuals :
      IF (qFiLoY) THEN ! logarithmic 20mar93
         DO i = 1, nC
            IF     (qrinside(Y2(i),-1.d-15,1.d-15)) THEN
               Delta(i) = 1.d15
            ELSEIF (Y2(i).lt.0.) THEN
               Delta(i) = 1.d18
            ELSEIF (YFi(i).lt.1.d-15) THEN
               Delta(i) = 0.
            ELSE
               ! fit only positive data
               Delta(i) = dlog10(YFi(i) / Y2(i)) * Weight(i)
               ENDIF
            ENDDO
      ELSE
         DO i = 1, nC
            Delta(i) = ( Y2(i) - YFi(i) ) * Weight(i)
            Delta(i) = dinside (Delta(i), -1.d15, 1.d15) ! for security 31dec92
            ENDDO
         ENDIF

C  Monit :
      IF (jFitMon.ge.3) THEN
         sum = rSum (Delta, 1, nC, 1)
         Print '(a,g9.2,a,10(g15.8))',
     *   ' >>>>     (', sum/nFitP, ')', (CuPar(j), j=1, nCuP)
         ENDIF
      IF (jFitMon.ge.4) THEN
         Print '(a,5a15)', ' >>>>>   ', 'i', 'X', 'Y', 'Y2', 'W', 'Res'
         DO i = 1, nC
      Print '(i5,5(g15.8))', i, XFi(i), YFi(i), Y2(i),
     *      Weight(i), Delta(i)
            ENDDO
         ENDIF

      END ! CuFitFunction

      SUBROUTINE CuFitFunction_MP (nC, nFitP, FitPar, Delta, iflag) !Artem: chenge interface
C     ----------------------------------------------------------
         ! called by E04FCF in CuFit
         ! calculates the residuals Delta(x) = Curve(x) - Dat(x)

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      DIMENSION    FitPar(nFitP), Delta(nC), Y2(MC), CuPar(MFP)
      CHARACTER    FitFehler*80

      COMMON / CuFiCo / XFi(MC), YFi(MC), Weight(MC), ifc, qFiLoY,
     *          jFitMon, nCuP, qFix(MFP), FixPar(MFP), ScaPar(MFP)
      COMMON / CuFiFe / FitFehler

C  Transfer fit_parameter -> curve_parameter :
      iFitP = 0
      DO i = 1, nCuP
         IF (qFix(i)) THEN
            CuPar(i) = FixPar(i)
         ELSE
            iFitP = iFitP + 1
            CuPar(i) = ScaPar(iFitP) * FitPar(iFitP)
            ENDIF
         ENDDO
      IF (iFitP.ne.nFitP) CALL Absturz ('CuFitFunction_MP',
     *   'nFitP not reproducable')

C  Calculate curve :
      FitFehler = '&ff'
      CALL CuConvVal (ifc, CuPar, Y2, FitFehler)
      IF (FitFehler.ne.'&ff') THEN
         iflag = -1685
         RETURN
         ENDIF

C  Calculate the residuals :
      IF (qFiLoY) THEN ! logarithmic 20mar93
         DO i = 1, nC
            IF     (qrinside(Y2(i),-1.d-15,1.d-15)) THEN
               Delta(i) = 1.d15
            ELSEIF (Y2(i).lt.0.) THEN
               Delta(i) = 1.d18
            ELSEIF (YFi(i).lt.1.d-15) THEN
               Delta(i) = 0.
            ELSE
               ! fit only positive data
               Delta(i) = dlog10(YFi(i) / Y2(i)) * Weight(i)
               ENDIF
            ENDDO
      ELSE
         DO i = 1, nC
            Delta(i) = ( Y2(i) - YFi(i) ) * Weight(i)
            Delta(i) = dinside (Delta(i), -1.d15, 1.d15) ! for security 31dec92
            ENDDO
         ENDIF

C  Monit :
      IF (jFitMon.ge.3) THEN
         sum = rSum (Delta, 1, nC, 1)
         Print '(a,g9.2,a,10(g15.8))',
     *   ' >>>>     (', sum/nFitP, ')', (CuPar(j), j=1, nCuP)
         ENDIF
      IF (jFitMon.ge.4) THEN
         Print '(a,5a15)', ' >>>>>   ', 'i', 'X', 'Y', 'Y2', 'W', 'Res'
         DO i = 1, nC
      Print '(i5,5(g15.8))', i, XFi(i), YFi(i), Y2(i),
     *      Weight(i), Delta(i)
            ENDDO
         ENDIF

      END ! CuFitFunction_MP

      SUBROUTINE CuFitMonit ( nC, nFitP, FitPar, Delta, Fjac, nFjac,
     *                        Sing, iGrade, nIter, nCalls,
     *                        Iwork, nIwork, Rwork, nRwork)
C     -----------------------------------------------------------
         ! called by E04FCF in CuFit
         ! to monit the fit process

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE     'i_dim.f'
      DIMENSION    FitPar(nFitP), Delta(nC), Fjac(nFjac,nFitP),
     *             Sing(nFitP), Iwork(nIwork), Rwork(nRwork)

      COMMON / CuFiCo / XFi(MC), YFi(MC), Weight(MC), ifc, qFiLoY,
     *         jFitMon, nCuP, qFix(MFP), FixPar(MFP), ScaPar(MFP)

      sum = rSum (Delta, 1, nC, 1)
      Print '(a,i3,a,i3,a,g9.2,a,10(g10.2))',
     *   ' >>> ', nIter, '/', nCalls, '(', sum/nFitP, ')',
     *   (FitPar(j), j=1, nFitP)
      IF (jFitMon.gt.1) THEN
         Print '(a7,6a10)', 'i', 'X', 'Y', 'R', 'dR/dP1','dR/dP2','..'
         DO i = 1, nC
            Print '(i7,13(g10.2))',
     *            i, XFi(i), YFi(i), Delta(i), (Fjac(i,j), j=1,nFitP)
            ENDDO
         Print *
         ENDIF

      END ! CuFitMonit

C  ====================================================================
C  i6.5.  Curve/ Interface
C  ====================================================================

      SUBROUTINE EditCPar (nJList, JList, qNewIn, Fehler)
C     ---------------------------------------------------
            ! JWu CuSetPar 6jun91; recopied from EditRPar 4aug99
         ! List and change curve parameters

      IMPLICIT NONE
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      CHARACTER*(*) Fehler
      CHARACTER     cfix*1, h1*16, h2*16
      CHARACTER*40  file, Co2, Un2, Co(MFP), Un(MFP), h3, ein, ein1
      CHARACTER*80  ffo, ffo2, aus
      LOGICAL       qCurve, qNew, qNewIn, qPuni(MFP)
      INTEGER       nJList, JList(*), lj, j, j1, ifu, ifu2, nP, nP2,
     *              iP, niPList, iPList(MFP), ipl, K, iFixed,
     *              iFix(MFP), iZ, j2, K2, nK2, nC2, nK, KK1(MK)
      REAL*8        rval, z, ZZ1(MK), tol

      qNew = qNewIn

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Learn from file 1 :
      j1 = JList(1)

      qCurve = qOlfGdef (j1, '?cu', 0, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF (.not.qCurve) THEN
         Fehler = ' File is not a curve'
         RETURN
         ENDIF

      ifu = iOlfG (j1, 'fu#', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      nP = iOlfG (j1, '#fit-par', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL tOlfG (j1, 'fit-formula', ffo, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      iFixed = iOlfG (j1, '@fixed', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO iP = 1, nP
         CALL OlfCnuG (j1, 'p'//cl2(iP), Co(iP), Un(iP), Fehler)
         iFix(iP)   = iBitGet (iFixed, iP)
         ENDDO ! iP

C  Compare with other files :
      DO lj = 2, nJList
         j = JList(lj)

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (.not.qCurve) THEN
            CALL Compose2 (Fehler, ' File '//cl3(j), ' is not a curve')
            RETURN
            ENDIF

         ifu2 = iOlfG (JList(1), 'fu#', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (ifu2.ne.ifu) THEN
            CALL Compose2 (Fehler, ' File '//cl3(j),
     *           ' is another function than file '//cl3(j1))
            RETURN
            ENDIF

         nP2 = iOlfG (j, '#fit-par', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (nP2.ne.nP) THEN
            Fehler = 'PROGRAM ERROR/ nP<> although ifu=='
            RETURN
            ENDIF

         CALL tOlfG (j, 'fit-formula', ffo2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (ffo2.ne.ffo) THEN
            Fehler = 'PROGRAM ERROR/ formulae <>'
            RETURN
            ENDIF

         iFixed = iOlfG (j, '@fixed', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         DO iP = 1, nP
            CALL OlfCnuG (j, 'p'//cl2(iP), Co2, Un2, Fehler)
            IF (Co2.ne.Co(iP) .or. Un2.ne.Un(iP)) THEN
               CALL Compose2 (Fehler, ' File '//cl3(j),
     *           ' has other parameter(s) than file '//cl3(j1))
               ENDIF
            IF (iBitGet(iFixed,iP).ne.iFix(iP)) iFix(iP) = -1 ! varies
            ENDDO

         ENDDO ! weitere files

C  Another loop over all files, prepare for modifications :
      DO lj = 1, nJList
         j = JList(lj)

         NofJ(lj) = iOlfG (j, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfOpen (j, 1, K, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO

      Print *, 'Modify curve(s): '//aus(1:lenU(aus))
      Print *

C  Preset ?
      IF (qNew) THEN
         DO iP = 1, nP
            X(iP) = 0
            Y(iP) = 0
            ENDDO
         DO lj = 1, nJList
            DO K = 1, NofJ(lj)
               CALL OlfPutXY0 (JList(lj), K, nP, X, Y, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO
            ENDDO
         ENDIF

C  Main loop (code efficient, not execution time efficient) :
 10   CONTINUE

C  Compare values :
      CALL OlfGetY (j1, 1, nP2, Y1, Fehler) ! will become min
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfGetY (j1, 1, nP2, Y2, Fehler) ! will become max
      IF (Fehler.ne.'&ff') RETURN
      DO lj = 1, nJList
         DO K = 1, NofJ(lj)
            CALL OlfGetY (JList(lj), K, nP2, Y, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (nP2.ne.nP) THEN
               Fehler = 'PROGRAM ERROR / nP<>'
               RETURN
               ENDIF
            DO iP = 1, nP
               Y1(iP) = dmin1 (Y1(iP), Y(iP))
               Y2(iP) = dmax1 (Y2(iP), Y(iP))
               ENDDO
            ENDDO
         ENDDO

C  List :
      DO iP = 1, nP
         IF (Y1(iP).ne.Y2(iP)) THEN
            write (h1, '(g16.8)') Y1(iP)
            write (h2, '(g16.8)') Y2(iP)
            h3 = 'varies from '//h1//' to '//h2
            qPuni(iP) = .false.
         ELSE
            write (h3, '(g16.8)') Y1(iP)
            qPuni(iP) = .true.
            ENDIF

         IF     (iFix(iP).eq.0) THEN ! free
            cfix = ' '
         ELSEIF (iFix(iP).eq.1) THEN ! fixed
            cfix = '!'
         ELSE                        ! varies
            cfix = '?'
            ENDIF

         Print '(i3,3x,a16,2x,a12,1x,a1,1x,a)', iP, Co(iP),
     *         Un(iP), cfix, h3
         ENDDO
      Print *

C  Menu :
 30      CONTINUE
         IF (qNew) THEN
            ein1 = '*'
            qNew = .false.
            Print *, ' Set initial values :'
         ELSE
            CALL FrageCD ('  Modify parameters (<list>,-,a,p,f)',
     *                    ein1, '-')
            ENDIF
         Print *

         IF (ein1.eq.'?') THEN
            Print *, '    <list>  modify existing parameters'
            Print *, '    -       quit'
            Print *, '    a       add to plot'
            Print *, '    p       plot'
            Print *, '    f       fit'
            ein1 = ' '
         ELSEIF (ein1.eq.'a') THEN ! add
            CALL IdaPlot ('a', nJList, JList, ' ', Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (ein1.eq.'p') THEN ! plot
            CALL IdaPlot ('p', nJList, JList, ' ', Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (ein1.eq.'f') THEN ! fit
            CALL CuFitCall (nJList, JList, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (ein1.eq.'-') THEN
            RETURN
         ELSEIF (ein1.eq.' ') THEN ! re-display
         ELSE ! list given or bad input
 304        CONTINUE
            CALL DecJList (ein1, MFP, niPList, iPList, 1, nP, Fehler)
            ein1 = ' '
            IF (Fehler.ne.'&ff') THEN
               CALL FehlerGong (Fehler, 1)
               Print *
               GOTO 10
               ENDIF
            DO ipl = 1, niPList
               iP = iPList(ipl)

C  Inner loop: modifications on one parameter
 31            CONTINUE
C  Submenu :
               CALL FrageC (' '//Co(iP)(1:lenU(Co(iP)))//
     *             ' (val,f,d,z,i,s,n,u,l) [none]', ein)
               CALL Fi1R (ein, rval)

               IF (ein.eq.'#') THEN ! global set
                  DO lj = 1, nJList
                     DO K = 1, NofJ(lj)
                        CALL OlfGetXY (JList(lj), K, nP2, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        Y(iP) = rval
                        CALL OlfPutXY0 (JList(lj), K, nP, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        ENDDO
                     ENDDO
               ELSEIF (ein.eq.'f') THEN
                  DO lj = 1, nJList
                     rval = rAsk ('   Value for file '//cl4(JList(lj)))
                     DO K = 1, NofJ(lj)
                        CALL OlfGetXY (JList(lj), K, nP2, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        Y(iP) = rval
                        CALL OlfPutXY0 (JList(lj), K, nP, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        ENDDO
                     ENDDO
               ELSEIF (ein.eq.'d') THEN
                  rval = 0
                  DO lj = 1, nJList
                     DO K = 1, NofJ(lj)
                        CALL OlfGetXY (JList(lj), K, nP2, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        CALL Compose2 (aus,
     *                       '   Value for file '//cl3(JList(lj)),
     *                       ' spectrum '//cl3(K))
                        Y(iP) = rAskD (aus, rval)
                        CALL OlfPutXY0 (JList(lj), K, nP, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        ENDDO
                     ENDDO
               ELSEIF (ein.eq.'z') THEN
                  iZ = iAskDMu ('Which z', iZ, 1, MZ)
                  DO lj = 1, nJList
                     DO K = 1, NofJ(lj)
                        CALL OlfGet1Z (JList(lj), K, iZ, z, Fehler)
                        CALL OlfGetXY (JList(lj), K, nP2, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        Y(iP) = z
                        CALL OlfPutXY0 (JList(lj), K, nP, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        ENDDO
                     ENDDO
               ELSEIF (ein.eq.'i') THEN
                  DO lj = 1, nJList
                     j = JList(lj)
                     Print *, ' operating on file ', j
                     j2 = iAskD (' Values from integral file', j2)
                     nK2 = iOlfG (j2, '#spectra', Fehler)
                     IF (Fehler.ne.'&ff') THEN
                        CALL FehlerGong (Fehler,1)
                        GOTO 31
                        ENDIF
                     IF (nK2.gt.1) THEN
                        aus = ' Take which of the '//cl3(nK2)
                        CALL Append (aus, ' spectra in the 2nd file')
                        K2 = iAskDMu (aus, mod(K2,nK2)+1, 0, nK2)
                        IF (K2.eq.0) RETURN
                     ELSE
                        K2 = 1
                        ENDIF
                     CALL OlfGetXY (j2, K2, nC2, X2, Y2, Fehler)
                     IF (iabs(irSorted(X2,nC2)).lt.2)
     *                  Fehler = ' 2nd file contains double entries'
                     IF (Fehler.ne.'&ff') THEN
                        CALL FehlerGong (Fehler,1)
                        GOTO 31
                        ENDIF
                     CALL OlfGet1ZofK (j, 1, nK, ZZ1, Fehler)
                     IF (Fehler.ne.'&ff') THEN
                        CALL FehlerGong (Fehler,1)
                        GOTO 31
                        ENDIF
                     tol = 1.d-8
                     CALL GetIndex (' searching z',
     *                    ZZ1, nK, X2, nC2, KK1, tol, Fehler)
                     IF (Fehler.ne.'&ff') THEN
                        CALL FehlerGong (Fehler,1)
                        GOTO 31
                        ENDIF
                     DO K = 1, nK
                        CALL OlfGetXY (JList(lj), K, nP2, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        Y(iP) = Y2(KK1(K))
                        CALL OlfPutXY0 (JList(lj), K, nP, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        ENDDO
                     ENDDO ! lj
c                     qModif = .true.
               ELSEIF (ein.eq.'a') THEN
                  DO lj = 1, nJList
                     DO K = 1, NofJ(lj)
                        CALL OlfGetXY (JList(lj), K, nP2, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        Y(iP) = dabs (Y(iP))
                        CALL OlfPutXY0 (JList(lj), K, nP, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        ENDDO
                     ENDDO
               ELSEIF (ein.eq.'s') THEN
                  IF (iFix(iP).eq.0) THEN
                     iFix(iP) = 1
                  ELSE
                     iFix(iP) = 0
                     ENDIF
                  DO lj = 1, nJList
                     iFixed = iOlfG (JList(lj), '@fixed', Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     CALL BitSet (iFixed, iP, iFix(iP))
                     CALL iOlfP (JList(lj), '@fixed', iFixed, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO
               ELSEIF (ein.eq.'n') THEN
                  CALL FrageCD ('   Coordinate', Co(iP), Co(iP))
                  CALL FrageCD ('   And its unit', Un(iP), Un(iP))
                  DO lj = 1, nJList
                     CALL OlfCnuP (JList(lj),
     *                    'p'//cl2(iP), Co(iP), Un(iP), Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO
               ELSEIF (ein.eq.'u') THEN
                  CALL FrageCD ('   Unit of coordinate '//Co(iP),
     *                 Un(iP), Un(iP))
                  DO lj = 1, nJList
                     CALL OlfCnuP (JList(lj),
     *                    'p'//cl2(iP), Co(iP), Un(iP), Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO
               ELSEIF (ein.eq.'l') THEN
                  Print '(a)', '  #F #S     z1     '//Co(iP)
                  DO lj = 1, nJList
                     DO K = 1, NofJ(lj)
                        CALL OlfGet1Z (JList(lj), K, 1, z, Fehler)
                        CALL OlfGetXY (JList(lj), K, nP2, X, Y, Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                        Print '(2i3,2g12.4)', j, K, z, Y(iP)
                        ENDDO
                     ENDDO
                  Print *
               ELSEIF (ein.eq.'h' .or. ein.eq.'?') THEN
                  Print *,
     * '    <real value>  global parameter value'
                  Print *,
     * '    f             different value for each file'
                  Print *,
     * '    d             different value for each spectrum'
                  Print *, '    z             from z'
                  Print *,'    i             from y(z) (integral file)'
                  Print *, '    a             reset p:=|p|'
                  Print *, '    s             status (free/fixed)'
                  Print *,
     * '    sf            status (free/fixed) per file'
                  Print *, '    n             coordinate name and unit'
                  Print *, '    u             coordinate unit'
                  Print *, '    l             list current values'
                  Print *, '    <RETURN>      no modification'
                  GOTO 31
               ELSEIF (ein.eq.' ') THEN
                  DO lj = 1, nJList
                     CALL tOlfP (JList(lj), 'fit-range&weight',
     *                    'parameters as given', Fehler)
                     CALL OlfClos (JList(lj), NofJ(lj), Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO
                  RETURN
                  ! regular exit
               ELSE
                  CALL Gong (3)
                  GOTO 31
                  ENDIF
               ENDDO ! iP
            ENDIF
            Print *
            GOTO 10

      END ! EditCPar

      SUBROUTINE CuSetPar (nJList, JList, qNew, Fehler)
C     -------------------------------------------------
         ! JWu 6jun91.
         ! ask for the functional parameters (name/value/status)
         ! if qNew then default values are set

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'

      CHARACTER*(*) Fehler
      INTEGER       JList(*), KK1(MK)
      REAL*8        X(MC), Y(MC), ZZ1(MK), X2(MC), Y2(MC), PVal(MK,MFP)
      LOGICAL       qConst(MFP), qLisP(MFP), qLisK(MK), qFix(MFP)

      CHARACTER     cval*12, cfix*1, line(MFP)*40
      CHARACTER*40  Co(MFP), Un(MFP)
      CHARACTER*80  aus, ein, cLisP2, cLisP3, cLisP4, cLisK1

      DATA          cLisP2 /' '/, cLisP3 /' '/, cLisP4 /'*'/,
     *              cLisK1 /'*'/, qSame / .true./, iMod4/1/

C  Which curves ?
      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop curves :
      DO lj = 1, nJList
         j = JList(lj)

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (.not.qCurve) THEN
            Fehler = ' File is not a curve'
            RETURN
            ENDIF

         CALL tOlfG (j, 'fit-formula', aus, Fehler)
         Print *, ' Curve '//cl3(j)//': '//aus
         Print *

         nK = iOlfG (j, '#spectra', Fehler)
         nP = iOlfG (j, '#fit-par', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! preset of parameters: from curve, from previous curve, or zero :

         IF (lj.gt.1)
     * qSame = qAskD (' The same setting as before', intq(qSame))

         IF     (lj.gt.1 .and. qSame) THEN
            ! let everything unchanged
            ! but ... setting of names not yet included
            IF (nK.gt.nKold) THEN
               CALL Gong (1)
               Print *, 'WARNING/ parameters for new K set to old K=1'
               DO K = nKold+1, nK
                  DO i = 1, nP
                     PVal(K,i) = PVal(1,i)
                     ENDDO
                  ENDDO
               ENDIF

         ELSE
            DO i = 1, nP
               CALL OlfCnuG (j, 'p'//cl2(i), Co(i), Un(i), Fehler)
               ENDDO
            IF (qNew) THEN
               ! set default values :
               DO i = 1, nP
                  DO K = 1, nK
                     PVal(K,i) = 0.
                     ENDDO
                  ENDDO
               iFixed = 0
            ELSE
               ! read the values :
               iFixed = iOlfG (j, '@fixed', Fehler)
               DO K = 1, nK
                  CALL OlfGetXY (j, K, nC, X, Y, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  DO i = 1, nP
                     PVal(K,i) = Y(i)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF

         DO i = 1, nP
            qFix(i) = qintr (iBitGet (iFixed, i))
            ENDDO

         ! Create table of parameters :
         qFirstSet = qNew .and. (lj.eq.1 .or. .not. qSame)
 2       CONTINUE
         DO i = 1, nP
            val = PVal(1,i)
            qConst(i) = .true.
            cfix      = ' '
            IF (qFix(i)) cfix = '!'
            IF (nK.gt.1) THEN
               ! check whether the value is the same for all spectra
               DO K = 2, nK
                  IF (PVal(K,i).ne.val)  qConst(i) = .false.
                  ENDDO
               ENDIF
            IF (qConst(i)) THEN
               write (cval,'(g12.5)') val
            ELSE
               cval = 'par = par(K)'
               ENDIF

            line(i) = cr2(i)//' '//Co(i)(1:10)//' '//Un(i)(28:35)//
     *                cfix//' '//cval
            ENDDO

         ! Show the table :
         aus = ' # name       unit     s value        '
         IF (nP.le.8) THEN
            Print *, aus
            DO i = 1, nP
               Print *, line(i)
               ENDDO
         ELSE
            Print *, aus(1:38), aus(1:38)
            DO i = 1, nP/2
               Print *, line(i)(1:38), line(nP-nP/2+i)(1:38)
               ENDDO
            IF (nP-nP/2.ne.nP/2) Print *, line(nP-nP/2)(1:38)
            ENDIF
         Print *

C  Menu (renewed 6aug93) :
 21      CONTINUE
         IF (qFirstSet) THEN
            CALL qSet (qLisP, 1, nP, 1, .true.) ! initialize ALL parameters
            Print *, ' Set initial values :'
            qFirstSet = .false.
         ELSE
            cLisP2 = '-'
            CALL GetNList (' Modify parameters', cLisP2, qLisP, nP)
            ENDIF

         IF (iqSum(qLisP,nP).eq.0) THEN
            ! save full setup :
            CALL OlfOpen (j, 1, Kout, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            DO i = 1, nP
               CALL BitSet (iFixed, i, intq(qFix(i)))
               CALL OlfCnuP (j, 'p'//cl2(i), Co(i), Un(i), Fehler)
               ENDDO
            CALL iOlfP (j, '@fixed', iFixed, Fehler)
            CALL tOlfP (j, 'fit-range&weight',
     *                  'parameters as given', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            DO K = 1, nK
               CALL OlfGetXY (j, K, nC, X, Y, Fehler) ! to get X
               DO i = 1, nP
                  Y(i) = PVal(K,i)
                  ENDDO
               CALL OlfPutXY0 (j, K, nP, X, Y, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO ! K
            CALL OlfClos (j, nK, Fehler)
            GOTO 9
            ENDIF

C  Loop over parameters :
         qModif = .false.
         DO i = 1, nP
         IF (qLisP(i)) THEN

 110        CONTINUE
            CALL Compose2 (aus, ' Enter '//Co(i),
     *                     ' (l,<value>,d,f,z,a,s,n,u,?)')
            CALL FrageH (aus, ein) ! this way 20mar93

            CALL Fi1R (ein, val)

            IF     (ein.eq.'?') THEN
                Print *, ' INPUT HELP/'
                Print *,
     * '    enter either a real value or one of the following :'
                Print *,'    l   list values without changing anything'
                Print *,
     * '    d   enter different values for each spectrum'
                Print *, '    f   take values from an integral file'
                Print *, '    z   take coordinate z'
                Print *, '    a   reset  p := |p|'
                Print *, '    s   toggle the status (fixed/free)'
                Print *, '    n   rename the parameter'
                Print *, '    u   rename its unit'
                GOTO 110
            ELSEIF (ein.eq.'#') THEN
                     DO K = 1, nK
                        PVal (K,i) = val
                        ENDDO
               qModif = .true.
            ELSEIF (ein.eq.'d') THEN
                     IF (nK.gt.2) THEN
                        aus = ' Change parameters of spectra'
                        CALL GetNList (aus, cLisK1, qLisK, nK)
                     ELSE
                        CALL qSet (qLisK, 1, nK, 1, .true.)
                        ENDIF
                     DO K = 1, nK
                        IF (qLisK(K)) THEN
                           aus = ' Value for spectrum '//cl6(K)
                           PVal(K,i) = rAskD (aus, PVal(K,i))
                           ENDIF
                        ENDDO
               qModif = .true.
            ELSEIF (ein.eq.'l') THEN ! 8aug93
               Print '(a)', '  #     z      '//Co(i)
               DO K = 1, nK
                  CALL OlfGet1Z (j, K, 1, z, Fehler)
                  Print '(i3,2g12.4)', K, z, PVal(K,i)
                  ENDDO
               Print *
               GOTO 110
            ELSEIF (ein.eq.'f') THEN
               j2 = iAskD (' Values from integral file', j2)
               nK2 = iOlfG (j2, '#spectra', Fehler)
               IF (Fehler.ne.'&ff') THEN
                  CALL FehlerGong (Fehler,1)
                  GOTO 110
                  ENDIF
               IF (nK2.gt.1) THEN
                  aus = ' Take which of the '//cl3(nK2)
                  CALL Append (aus, ' spectra in the 2nd file')
                  K2 = iAskDMu (aus, mod(K2,nK2)+1, 0, nK2)
                  IF (K2.eq.0) RETURN
               ELSE
                  K2 = 1
                  ENDIF
               CALL OlfGetXY (j2, K2, nC2, X2, Y2, Fehler)
               IF (iabs(irSorted(X2,nC2)).lt.2)
     *            Fehler = ' 2nd file contains double entries'
               IF (Fehler.ne.'&ff') THEN
                  CALL FehlerGong (Fehler,1)
                  GOTO 110
                  ENDIF
               CALL OlfGet1ZofK (j, 1, nK, ZZ1, Fehler)
               IF (Fehler.ne.'&ff') THEN
                  CALL FehlerGong (Fehler,1)
                  GOTO 110
                  ENDIF
               tol = 1.d-8
               CALL GetIndex (' searching z',
     *              ZZ1, nK, X2, nC2, KK1, tol, Fehler)
               IF (Fehler.ne.'&ff') THEN
                  CALL FehlerGong (Fehler,1)
                  GOTO 110
                  ENDIF
               DO K = 1, nK
                  PVal(K,i) = Y2(KK1(K))
                  ENDDO
               qModif = .true.
            ELSEIF (ein.eq.'z') THEN
               DO K = 1, nK
                  CALL OlfGet1Z (j, K, 1, z, Fehler)
                  PVal(K,i) = z
                  ENDDO
            ELSEIF (ein.eq.'a') THEN
                     DO K = 1, nK
                        PVal(K,i) = dabs ( PVal(K,i) )
                        ENDDO
               qModif = .true.
            ELSEIF (ein.eq.'s') THEN ! toggle status
               qFix(i) = .not.qFix(i)
               qModif = .true.
               GOTO 110
            ELSEIF (ein.eq.'n') THEN ! Modify names :
               ! should be done rueckwaerts : erst p10, dann p1 ersetzen
               CALL FrageTD (' Name for parameter '//
     *                       cl2(i), Co(i), Co(i))
c96/7               CALL ReplaceT (formula, Co, ein)
               qModif = .true.
               GOTO 110
            ELSEIF (ein.eq.'u') THEN ! Modify units :
               CALL FrageCD (' Unit of parameter '//Co(i),Un(i),Un(i))
               qModif = .true.
               GOTO 110
            ELSEIF (ein.eq.' ') THEN
               ! do nothing
            ELSE
               CALL Gong (3)
               GOTO 110
               ENDIF

         ENDIF
         ENDDO
C  End loop over parameters.

         IF (qModif) THEN ! something has been changed
            GOTO 2  ! show modified table
         ELSE
            GOTO 21 ! directly to menu (8aug93)
            ENDIF
C  End loop modifications.

 9       CONTINUE
         nKold = nK
         ENDDO
C  End loop curves.

      END ! CuSetPar

      SUBROUTINE CuRedef (nJList, JList, Fehler)
C     ------------------------------------------
         ! JWu 6jun97.
         ! call `cnn': change function no. of existing curve
         ! ask for the functional parameters (name/value/status)
         ! if qNew then default values are set

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'

      CHARACTER*(*) Fehler
      INTEGER       JList(*), KK1(MK)
      REAL*8        X(MC), Y(MC), ZZ1(MK), X2(MC), Y2(MC), PVal(MK,MFP)
      LOGICAL       qConst(MFP), qLisP(MFP), qLisK(MK), qFix(MFP)

      CHARACTER     cval*12, cfix*1, line(MFP)*40
      CHARACTER*40  Co(MFP), Un(MFP)
      CHARACTER*80  aus, ein, cLisP2, cLisP3, cLisP4, cLisK1

      DATA          cLisP2 /' '/, cLisP3 /' '/, cLisP4 /'*'/,
     *              cLisK1 /'*'/, qSame / .true./, iMod4/1/

C  Which curves ?
      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop curves :
      DO lj = 1, nJList
         j = JList(lj)

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (.not.qCurve) THEN
            Fehler = ' File is not a curve'
            RETURN
            ENDIF

         CALL tOlfG (j, 'fit-formula', aus, Fehler)
         Print *, ' Curve '//cl3(j)//': '//aus
            ! save full setup :

         CALL OlfOpen (j, 1, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         nP = iOlfG (j, '#fit-par', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         DO i = 1, nP
            CALL BitSet (iFixed, i, intq(qFix(i)))
            CALL OlfCnuP (j, 'p'//cl2(i), Co(i), Un(i), Fehler)
            ENDDO
         CALL iOlfP (j, '@fixed', iFixed, Fehler)
         CALL tOlfP (j, 'fit-range&weight', 'parameters as given',
     *               Fehler)
         CALL OlfComAdd (j, ' ', 'ACHTUNG: Funktion umdefiniert',
     *                   Fehler)
         CALL OlfClos (j, nK, Fehler)

         ENDDO
C  End loop curves.

      END ! CuRedef

      SUBROUTINE CuSetAux (nJList, JList, Fehler)
C     -------------------------------------------
         ! JWu 6jun91
         ! ask for all auxiliary parameters of a curve

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*) Fehler
      INTEGER       JList(*)
      CHARACTER     aus*80, line*80, un*40
      DATA          qSame /.true./

C  Which curves ?
      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop curves :
      DO lj = 1, nJList
         j = JList(lj)

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (.not.qCurve) THEN
            Fehler = ' File is not a curve'
            RETURN
            ENDIF
         CALL OlfOpen (j, 1, Kout, Fehler)
         nK = iOlfG (j, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (lj.gt.1) THEN
            CALL Compose2 (aus, ' Modify curve '//cl2(j),
     *                     ' : same answers as before')
            qSame = qAskD (aus, intq(qSame))
            ENDIF

         IF (lj.eq.1 .or. .not. qSame) THEN
            ri = rOlfG (j, 'fit-i', un, Fehler)
            rf = rOlfG (j, 'fit-f', un, Fehler)
            CALL rAskRgeFull (' Range of fit', ri, rf, ri, rf)
            CALL rOlfP (j, 'fit-i', ' ', ri, Fehler)
            CALL rOlfP (j, 'fit-f', ' ', rf, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            qwsx = qOlfGdef (j, '?weight-stp-x', 0, Fehler)
            qwey = qOlfGdef (j, '?weight-err-y', 0, Fehler)
            qwly = qOlfGdef (j, '?weight-log-y', 0, Fehler)
            qcon = qOlfGdef (j, '?conv',         0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            qwsx = qAskD (' Weighting with x-stepwidth', intq(qwsx))
            qwey = qAskD (' Weighting with y-error', intq(qwey))
            qwly = qAskD (' Logarithmic weighting of ', intq(qwly))
            qcon = qAskD (' Convolute with resolution', intq(qcon))
            CALL iOlfP (j, '?weight-stp-x', intq(qwsx), Fehler)
            CALL iOlfP (j, '?weight-err-y', intq(qwey), Fehler)
            CALL iOlfP (j, '?weight-log-y', intq(qwly), Fehler)
            CALL iOlfP (j, '?conv',         intq(qcon), Fehler)

         ELSE
            CALL iOlfCopy (jold, j, '?weight-stp-x', Fehler)
            CALL iOlfCopy (jold, j, '?weight-err-y', Fehler)
            CALL iOlfCopy (jold, j, '?weight-log-y', Fehler)
            CALL iOlfCopy (jold, j, '?conv',         Fehler)
            CALL rOlfCopy (jold, j, 'fit-i',         Fehler)
            CALL rOlfCopy (jold, j, 'fit-f',         Fehler)
            ENDIF

         CALL OlfClos (j, nK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         jold = j
         ENDDO

      END ! CuSetAux

      SUBROUTINE CuSetFit (Fehler)
C     ----------------------------
            ! JWu 25feb92
         ! setup for fitroutine (via COMMON /FitSet/)

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      CHARACTER*(*)     Fehler

      COMMON / FitSet / tolFit, stpFit, etaFit1, etaFitN,
     *                  mclFit, qRPFit,
     *                  imoFit, jmoFit, idiFit

      Print *, ' setup for fitroutine : '

      tolFit = rAskDMu (' Tolerance', tolFit, 0.d0, 1.d3)
      stpFit = rAskDMu (' Maximal variation', stpFit, 1.d-3, 1.d9)
      etaFitN= rAskDMu (
     *  ' Eta (small: more differentiations - large: more trials)',
     *    etaFitN, 0.d0, 1.d0)
      etaFit1= rAskDMu (' dito (for case of only one fitparameter)',
     *    etaFit1, 0.d0, 1.d0)
      qRPFit = qAskD   (' Reduce parameters to 1', intq(qRPFit))
      mclFit = iAskDMu (' Function calls', mclFit, 1, 10000)
      imoFit = iAskDMu (' Frequency of monitor calls',imoFit,-1,1000)
      IF (imoFit.gt.0)
     *jmoFit = iAskDMu (' Monitor what', jmoFit, 1, 4)
      idiFit = iAskDMu (' Display matrices', idiFit, 0, 1)

      END ! CuSetFit

      SUBROUTINE CuGetPar (nJList, JList, Fehler)
C     -------------------------------------------
         ! JWu 8jul91. Aufgeblaeht 27jan92 (2 loops).
         ! curve -> normal data file, with one spectrum,
         ! containing the results for one parameter.

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       MFP1
      PARAMETER    (MFP1=MFP+1)
      CHARACTER*(*) Fehler
      CHARACTER*40  Un, Co, file
      INTEGER       nJList, JList(*), nFJdef, j, j1, lj, nK, K, Kout,
     *              ifc, ifc2, iP, iiP, nP, nC, iOlfG, ih1, ih2
      REAL*8        X(MC), Y(MC), D(MC), X1(MC), Y1(MC), D1(MC),
     *              rOlfG, fri, frf
      LOGICAL       qLisP(MFP1), qCurve, qOlfGdef
      CHARACTER     aus*80, cLisP*80, h1*16, h2*16, hc1*1, hc2*1, hc3*1

      DATA          cLisP /' '/, nFJdef /4/

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  List parameters from first of given curves :
      DO lj = 1, nJList
         j = JList(lj)
         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurve) GOTO 11
         ENDDO
      Fehler = 'no curve given'
      RETURN
 11   CONTINUE

      nP = iOlfG (j, '#fit-par', Fehler)
      ifc= iOlfG (j, 'fu#', Fehler)
      ! construct list :
      aus = ' parameters are '
      DO iP = 1, nP
         CALL OlfCnuG (j, 'p'//cl2(iP), Co, Un, Fehler)
         CALL Append (aus, ', '//Co)
         ENDDO
      Print *, aus

      CALL Compose2 (aus, ' Extract which parameters ('//
     *               cl6(MFP1), '=sum)')
      CALL GetNList (aus, cLisP, qLisP, MFP1)

C  Outer loop : parameters :
      DO iP = 1, MFP1
      IF (qLisP(iP)) THEN

C  Loop over curves - exceptionally, this is an inner loop :
         DO lj = 1, nJList
            j = JList(lj)
            qCurve = qOlfGdef (j, '?cu', 0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (.not.qCurve) GOTO 59 ! ' File is not a curve'

            ! test ifc=ifc, then we are sure that nP=nP :
            ifc2 = iOlfG (j, 'fu#', Fehler)
            IF (ifc2.ne.ifc) THEN
               Fehler = ' not the same function'
               RETURN
               ENDIF

            ! open output data-file :
            CALL OlfHeadDup (j, .false., j1, nK, Kout, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            CALL iOlfP (j1, '?cu', 0, Fehler)
c96/7            CALL iOlfDel (j1, 'fu#', Fehler)
c96/7            CALL iOlfDel (j1, '?conv', Fehler)
c96/7            CALL iOlfDel (j1, 'fit-i', Fehler)
c96/7            ...
            CALL iOlfP (j1, 'plot-sy#', 0, Fehler)

            IF (iP.eq.MFP1) THEN
               CALL OlfCnuP (j1, 'y', 'Sum', ' ', Fehler)
            ELSE
               CALL OlfCnuG (j,  'p'//cl2(iP), Co, Un, Fehler)
               IF (Un.eq.'?') THEN
                  CALL FrageC  (' Unit of '//Co, Un)
                  ! overwrite the curve file :
                  CALL OlfOpen (j, 1, Kout, Fehler)
                  CALL OlfCnuP (j, 'p'//cl2(iP), Co, Un, Fehler)
                  CALL OlfClos (j, 0, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               CALL OlfCnuP (j1, 'y', Co, Un, Fehler)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            DO iiP = 1, nP
               CALL OlfCnuP (j1, 'p'//cl2(iiP), ' ', ' ', Fehler) ! 26feb99
               ENDDO

            CALL tOlfG (j, 'fil', file, Fehler)
            CALL Append (file, Co)
            CALL tOlfP (j1, 'fil', file, Fehler)

C  Documentation (26feb99) :
            hc1 = '-'
            hc2 = '-'
            hc3 = '-'
            IF (iOlfG(j, '?weight-stp-x', Fehler).ne.0) hc1 = 'x'
            IF (iOlfG(j, '?weight-err-y', Fehler).ne.0) hc2 = 'd'
            IF (iOlfG(j, '?weight-log-y', Fehler).ne.0) hc3 = 'l'
            IF (Fehler.ne.'&ff') RETURN
            CALL Compose3 (aus, 'p'//cl2(iP), ' from fit #'//cl3(ifc),
     *           ' weighted '//hc1//hc2//hc3)

            fri = rOlfG (j, 'fit-i', Un, Fehler)
            frf = rOlfG (j, 'fit-f', Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (fri.le.frf) THEN
               CALL Append (aus, ' over full range')
            ELSE
               CALL NiceNum (fri, h1, ih1)
               CALL NiceNum (frf, h2, ih2)
               CALL Append (aus, ' from '//h1(1:ih1)//
     *          ' to '//h2(1:ih2)//' '//Un)
               ENDIF

            CALL OlfComAdd (j1, '=', aus, Fehler)

C  Now extract the result :
            DO K = 1, nK
               CALL OlfGetXYD (j, K, nC, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               IF (iP.eq.MFP1) THEN
                  Y1(K) = X(1)
                  D1(K) = 0.
               ELSE
                  Y1(K) = Y(iP)
                  D1(K) = D(iP)
                  ENDIF
               ENDDO ! K

            CALL TensorSaveInt (j, j1, nK, Kout, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL OlfClos (j1, Kout, Fehler)
            IF (Fehler.ne.'&ff') RETURN

 59         CONTINUE
            ENDDO ! lj

         ENDIF
         ENDDO ! iP

      END ! CuGetPar

      SUBROUTINE CuCreate (nJList, JList, qOv, Fehler)
C     ------------------------------------------------
            ! Early 1991. As subroutine 28jan92.
            ! qOv 9nov96 shall allow for call by FileMake

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      CHARACTER*40   aux, Co, Un, CoX, UnX, CoY, UnY
      DIMENSION      Z0(MK), NofK0(MK), X1(MC), Y1(MC), D1(MC), IAux(2)
      CHARACTER*80   aus, Name, Formula, ParDef, ParUni

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Which function ?
 61   ifc = iAskDu (' Function no. [help]', 0)
      IF (ifc.eq.0) THEN ! help
         ishow = 0
         i = 0
 62      CONTINUE
            i = i+1
            CALL CuFuText (i, Name, Formula, ParDef, ParUni, nP)
            IF (Name.eq.'&lastentry') GOTO 61 ! regular exit from loop
            IF (Name.ne.'&undefined') THEN
               Print '(1x,i3,2x,a26,a47)', i, Name, Formula
               ishow = ishow + 1
               IF (mod(ishow,23).eq.0) THEN
                  ifc = iAskDu (' Function no. [continue help]', 0)
                  IF     (ifc.gt.0) THEN
                     GOTO 69
                  ELSEIF (ifc.lt.0) THEN
                     Fehler = ' '
                     RETURN
                     ENDIF
                  ENDIF
               ENDIF
            GOTO 62 ! endless loop
         ENDIF
 69      CONTINUE
      CALL CuFuText (ifc, Name, Formula, ParDef, ParUni, nP)
      IF (Name.eq.'&undefined') THEN
         CALL Gong (5)
         GOTO 61
         ENDIF

C  Create curves :
      DO lj = 1, nJList

         ! get data file :
         jd = JList(lj)

         qCurve = qOlfGdef (jd, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurve) THEN
            Fehler = ' referred file is a curve, not a data file'
            RETURN
            ENDIF

         CALL OlfHeadDup (jd, qOv, jc, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL iOlfP (jc, '?cu',           1,   Fehler)
         CALL iOlfP (jc, 'fu#',           ifc, Fehler)
         CALL iOlfP (jc, '#fit-par',      nP,  Fehler)
         CALL iOlfP (jc, '@fixed',        0,   Fehler)
         CALL iOlfP (jc, 'fit-dat-file#', jd,  Fehler)
         CALL iOlfP (jc, 'plot-#pts',     100, Fehler)
         CALL iOlfP (jc, 'plot-sy#',      1,   Fehler)

         CALL rOlfP (jc, 'fit-i',   ' ', 0.d0, Fehler)
         CALL rOlfP (jc, 'fit-f',   ' ', 0.d0, Fehler)
         CALL rOlfP (jc, 'plot-i',  ' ', 0.d0, Fehler)
         CALL rOlfP (jc, 'plot-f',  ' ', 0.d0, Fehler)

         ! file name and title :
         CALL tOlfG (jd, 'fil', aus, Fehler)
         CALL Append (aus, '-')
         CALL tOlfP (jc, 'fil', aus, Fehler)

         CALL tOlfG (jd, 'tit', aus, Fehler)
         CALL Insert (aus, 1, 'fit to ')
         CALL tOlfP (jc, 'tit', aus, Fehler)

         CALL OlfCnuG (jd, 'x', CoX, UnX, Fehler)
         CALL OlfCnuG (jd, 'y', CoY, UnY, Fehler)

         ! build comment :
         CALL CuFuText (ifc, Name, Formula, ParDef, ParUni, nP)

         jp = 0
         IF (ParDef(1:1).eq.'&') THEN
            CALL TakeVorDel (ParDef, aux, ';')
            IF     (aux.eq.'&j') THEN
               jp = iAskMu (' Parameter file', 1, MF)
            ELSE
               Fehler = 'unknown &-expression in ParDef : '//aux
               RETURN
               ENDIF
            ENDIF
         CALL iOlfP (jc, 'fit-par-file#', jp,  Fehler)

         ! coordonate names and units :
         DO i = 1, nP
            CALL TakeVorDel (ParDef, Co, ';')
            CALL ReplaceT (Formula, 'p'//cl2(i), Co)
            IF (ParUni.ne.' ') THEN ! unit - since 5jan93
               CALL TakeVorDel (ParUni, aux, ';')
               CALL FindI (aux, 2, niaux, IAux)
               IF (niaux.lt.2) THEN
                  Un = '??'
               ELSEIF (aux.eq.'#,#') THEN
                  CALL MultUnits (Un, UnX, IAux(1), UnY, IAux(2))
               ELSE
                  Un = '???'
                  ENDIF
            ELSE
               Un = '?'
               ENDIF
            CALL OlfCnuP (jc, 'p'//cl2(i), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         CALL tOlfP (jc, 'fit-name', 'curve : '//Name, Fehler)
         CALL tOlfP (jc, 'fit-formula', Formula, Fehler)
         CALL tOlfP (jc, 'fit-range&weight', 'parameters as given',
     *               Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComLinDelAll (jc, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfComAdd (jc, ' ', Formula, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         DO K = 1, nK
            CALL OlfCopZ   (jd, jc, K, K, Fehler)
            CALL OlfPutXYD (jc, K, nP, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO ! K
         CALL OlfClos (jc, nK, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         JList(lj) = jc ! replace file by curve as default
         ENDDO ! lj

C  Now the subroutines for modifying the parameters :
      CALL CuSetPar (nJList, JList, .true., Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL CuSetAux (nJList, JList, Fehler)

      END ! CuCreate

      SUBROUTINE CuFitCall (nJList, JList, Fehler)
C     --------------------------------------------
            ! Early 1991. As subroutine 28jan92.
         ! Interface between Ida and CuFitExe

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*), K2K(MK), K3K(MK), KList(MK)
      CHARACTER*80   aus, cLisK, Line

      DATA           cLisK /'*'/

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      DO lj = 1, nJList
         jc = JList(lj)

         CALL OlfOpen (jc, 1, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! get curve :
         qCurve = qOlfGdef (jc, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (.not.qCurve) THEN
            Fehler = ' file '//cl2(jc)//' is no curve'
            RETURN
            ENDIF
         nK     = iOlfG (jc, '#spectra', Fehler)
         ifc    = iOlfG (jc, 'fu#', Fehler)
         jddef  = iOlfG (jc, 'fit-dat-file#', Fehler)
         jpdef  = iOlfG (jc, 'fit-par-file#', Fehler)
         qConv  = qOlfG (jc, '?conv', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! get data :
         CALL Compose2 (aus, ' Fit curve '//cl2(jc), ' to data file')
         jd = iAskD (aus, jddef)
         IF (jd.eq.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         qCurv1 = qOlfGdef (jd, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurv1) THEN
            Fehler = ' file '//cl2(jd)//' is a curve'
            RETURN
            ENDIF
         ! check compatibility :
         nK1    = iOlfG (jd, '#spectra', Fehler)
         IF (nK1.ne.nK) THEN
            Fehler =
     * ' Data and curve file inconsistent - # spectra different'
            RETURN
            ENDIF
         ! data file seems o.k.
         CALL iOlfP (jc, 'fit-dat-file#', jd, Fehler)

         ! get conv-file :
         jconv = jCuConvAsk (qCurve, qConv)
         IF (jconv.ne.0) THEN
            CALL GetK2K (jc, jconv, nK, K2K, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDIF

         ! get par-file :
         IF (jpdef.ne.0) THEN
            jp = iAskDMu ('Parameter file', jpdef, 1, MF)
            CALL GetK2K (jc, jp, nK, K3K, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL iOlfP (jc, 'fit-par-file#', jp, Fehler)
         ELSE
            jp = 0
            ENDIF

         ! Fit which spectra ?
         IF (nK.gt.1) THEN
            CALL GetJList(' Fit which spectra',
     *                    cLisK, MK, nKList, KList, 1, nK)
            IF (nKList.lt.1) THEN
               Fehler = ' '
               RETURN
               ENDIF
         ELSE
            nKList   = 1
            KList(1) = 1
            ENDIF

         ! Loop K - fit spectra - call CuFitExe :
         DO lK = 1, nKList
            K = KList(lK)
            CALL CuFitExe (jd, jc, K, jconv, K2K(K), jp,
     *                     K3K(K), Line, Fehler)
            IF (Fehler.ne.'&ff') THEN
               Print *, ' Error during fit of spectrum '//cl6(K)
               RETURN
               ENDIF
            ENDDO ! lK

         CALL tOlfP (jc, 'fit-range&weight', Line, Fehler)

         ! close curve file :
         CALL OlfClos (jc, nK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! lj

      END ! CuFitCall
