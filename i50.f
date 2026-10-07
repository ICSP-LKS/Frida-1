C  ====================================================================
C
C      Library  IDA   :  Inelastic data treatment
C      Modul    i50   :     functional operations on data
C
C  ====================================================================

C     Contents :
C        Functional operations  Opr...
C            ...Integral, Differential, Pointwise, Tensorial

C     Record of development:
C     JWu   nov91 : operations from Ida4(200) to Ida5
C     JWu   mai91 : splitting of modules
C     JWu   feb91 : improvement, implementation in IDA
C     JWu   jan91 : main components
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE OprIntegral (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! any integral property, i.e. extraction of one value from a
         ! whole spectrum, e.g. the integral, or a value y(x_0),..
         ! JWu 7-9feb91

C     Input :  file j :    spectra  x(i,K) -> y(i,K) with z=z(K)
C     Output : file jout : spectrum z(K)   -> zi(K),
C        where zi is a functional of y(x(i,K))

      IMPLICIT NONE

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      INCLUDE 'i_wrk.f'

      CHARACTER*(*) Fehler
      INTEGER       iShift(MK), IzK(MK), nJList, JList(*), nJListO,
     *              JListO(MF), in, iMod, modIntNor, modIntDif,
     *              modIntWgt, modSav, MiMaMo, nMiMaGru, Knorm,
     *              j, j2, jout, lj, K, KK, K2, Kout, nK, iZ, iZ1,
     *              iZ2, nZ, nZ1, n, n2, na, i, i2, ii, if, iX0, iX,
     *              imax, nrel, nCom, ih1, ih2, ih3, ih4, ixdummy,
     *              nK1dummy
      REAL*8        Z(MZ), facIntNor, divIntNor, wIntL, wIntU,
     *              xi, xf, rmi, rma, rX0, y1norm, ymax,
     *              xdummy, ydummy, sdummy, ddummy, yZaehler, yNenner,
     *              yMatch, dwInt, ya, da
      LOGICAL       qOv, qFullRange, qSameRange, qCom2, qComN, qFJ,
     *              qSortmi
      CHARACTER*40  comA, h, h1, h2, h3, h4,
     *              coX, unX, coY, unY, coZ, unZ, CoG, UnG
      CHARACTER     aus*80, cZ*2

C  Defaults for dialogue :
      DATA          modIntDif /2/, modSav /1/, MiMaMo /1/,
     *              nMiMaGru /1/, Knorm /1/, iZ2 /1/, iZ /1/,
     *              facIntNor /1.d0/, divIntNor /1.d0/

C  Input :
      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Choice of the functional :
      Print *, ' creating a new file with one spectrum y''(x''),'
      Print *, ' the ordinate will be x''=z.'
      Print *, ' options for the abszissa y'' are :'
      Print *, '    (1) simply z              (2) one value x'
      Print *, '    (3) one value y(x)        (4) integral dx y(x)'
      Print *, '    (5) value y_max           (6) value y_min'
      Print *, '    (7) position x of y_max   (8) position x of y_min'
      Print *, '    (9) match #K to #K-1     (10) quality of match'
      Print *, '   (11) match #K to #K''      (12) quality of match'
      Print *, '   (13) half width           (14) integral-de-luxe'
      Print *, '   (15) sum                  (16) average'

      in = iAskDMu (' Option', iMod, 0, 16)
      IF     (in.eq. 0) THEN
         Fehler = ' '
         RETURN
         ENDIF
      iMod = in

C  Enter parameter :
      IF (iMod.eq. 1) THEN  ! just z
         iZ2 = iAskDMu (' Which z', iZ2, 1, MZ)
      ELSEIF (iMod.eq. 2) THEN  ! one x
         iX0 = iAskMu (' Channel no.', 1, MC)
      ELSEIF (iMod.eq. 3) THEN  ! one y
         iX0 = iAskDuMu(' Channel no. [by value]', 0, 0, MC)
         IF (iX0.eq.0) THEN
            aus = ' x-Coordinate ' ! ' in units of '//tPar(8)
            rX0 = rAskD (aus, rX0)
            ENDIF
      ELSEIF (iMod.eq. 4 .or. iMod.eq.14) THEN  ! int dx y(x)
         IF (iMod.eq.14) THEN
 111        modIntDif = iAskD (' 2-point or 4-point differences ?',
     *                         modIntDif)
            IF (modIntDif.ne.2 .and. modIntDif.ne.4) THEN
               CALL Gong (3)
               GOTO 111
               ENDIF
            ENDIF
 112        CONTINUE
         CALL rAskRgeFull (' Range for integration',
     *            wIntL, wIntU, wIntL, wIntU)
         IF (iMod.eq.14) THEN
            modIntWgt = iAskD (' Weight with power of x', modIntWgt)
            modIntNor = iAskDMu (
     *   ' Normalisation: none (0) x-range (1) multiply (2) divide (3)',
     *                           modIntNor, 0, 3)
            IF     (modIntNor.eq.0) THEN
               facIntNor = 1
            ELSEIF (modIntNor.eq.1) THEN
               IF (wIntL.eq.0 .and. wIntU.eq.0) THEN
                  Print *,
     *             'this choice requires finite integration range'
                  CALL Gong (3)
                  GOTO 112
                  ENDIF
               IF (modIntWgt.eq.-1) THEN
                  facIntNor = 1 / dlog (wIntU / wIntL)
               ELSE
                  facIntNor = (modIntWgt + 1) /
     *                 ( wIntU**(modIntWgt+1) - wIntL**(modIntWgt+1) )
                  ENDIF
            ELSEIF (modIntNor.eq.2) THEN
               facIntNor = rAskD (' Multiply by correction factor',
     *                            facIntNor)
            ELSEIF (modIntNor.eq.3) THEN
               divIntNor = rAskDMu (' Divide by correction factor',
     *                              divIntNor, 1.d-80, 1.d80)
               facIntNor = 1 / divIntNor
               ENDIF
         ELSE
            modIntWgt = 0
            modIntDif = 2
            facIntNor = 1
            ENDIF
      ELSEIF (qiinside(iMod,5,8)) THEN  ! options for IdaMinMax
         Print *, ' minimax options :'  ! these options added JWu 6jul91
         Print *, '    (1) simply take min/max'
         Print *, '    (2) average over the n smallest/largest values'
         Print *, '    (3) barycenter of the peak'
         MiMaMo = iAskDMu (' Option', MiMaMo, 0, 3)
         IF     (MiMaMo.eq.0) THEN
            RETURN
         ELSEIF (MiMaMo.eq.2) THEN
            nMiMaGru = iAskDMu (' Enter n', nMiMaGru, 1, MC)
            ENDIF
      ELSEIF (iMod.eq. 9 .or. iMod.eq.11) THEN ! match
c         IF (nK.gt.2) THEN
            qSameRange = qAskD (' The same x-range for all matches ?',
     *                          intq(qSameRange)) ! 6may96
c         ELSE
c            qSameRange = .true.
c            ENDIF
         IF (qSameRange) THEN
            CALL rAskRgeFull (' Which x-range for average',
     *                        rmi, rma, rmi, rma)
            qFullRange = (rmi.eq.0. .and. rma.eq.0.)
            ENDIF
         IF (iMod.eq.9) THEN
            Knorm = iAskDMu (' Normalize to weight of spectrum',
     *                       Knorm, 1, MK)
            ENDIF
      ELSEIF (iMod.eq.10 .or. iMod.eq.12) THEN ! test match
         qSameRange = .true. ! anders noch nicht m"oglich
         CALL rAskRgeFull (' Which x-range for average',
     *                     rmi, rma, rmi, rma)
         qFullRange = (rmi.eq.0. .and. rma.eq.0.)
         Print *, ' calculating sum of (y(K)-y(K+1))^2 / <y>^2'
      ELSEIF (iMod.eq.15) THEN  ! Sum y(x)
         CALL rAskRgeFull (' Range for sum',
     *            wIntL, wIntU, wIntL, wIntU)
      ELSEIF (iMod.eq.16) THEN  ! <y(x)>
         CALL rAskRgeFull (' Range for average',
     *            wIntL, wIntU, wIntL, wIntU)
         ENDIF

      modSav = iAskDMu (
     * ' Save result as file (1), as z or r (2), or not at all (3)',
     * modSav, 0, 3)
      IF (modSav.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop in j :
      DO lj = 1, nJList

         j = JList(lj)
         nK = iOlfG (j, '#spectra', Fehler)

         ! Second file needed ?
         IF (iMod.eq.11 .or. iMod.eq.12) THEN
            j2 = iAskDMu ('Second file', j2, 0, MF)
            IF (j2.eq.0) THEN
               Fehler = ' '
               RETURN
               ENDIF
            CALL GetK2K (j, j2, nK1dummy, IzK, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDIF

         CALL OlfCnuG (j, 'x', coX, unX, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuG (j, 'y', coY, unY, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! guess y coord and unit, compose comment :
         IF     (iMod.eq. 1) THEN
            CALL OlfCnuG (j, 'z'//cl2(iZ2), CoG, UnG, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            comA = 'z'//cl2(iZ2)
         ELSEIF (iMod.eq. 2) THEN
            CoG = coX
            UnG = unX
            CALL Compose2 (comA, 'x(ch'//cl3(iX0), ')')
         ELSEIF (iMod.eq. 3) THEN
            CoG = coY
            UnG = unY
            IF (iX0.eq.0) THEN ! by x-value
               CALL NiceNum (rX0, h1, ih1)
               comA = 'y('//h1(1:ih1)//')'
            ELSE ! by channel number
               CALL Compose2 (comA, 'y(ch'//cl3(iX0), ')')
               ENDIF
         ELSEIF (iMod.eq. 4 .or. iMod.eq.14) THEN
            CALL Compose2 (CoG, 'I d'//coX, ' '//coY)
            CALL Compose2 (h, unX, '-1')   ! e.g. 'meV-1*meV'
            IF (h.eq.unY) THEN
               UnG = ' '
            ELSE
               CALL Compose2 (h, unY, '-1')
               IF (h.eq.unX) THEN
                  UnG = ' '
               ELSE
                  CALL Compose2 (UnG, unX, '*'//unY)
                  ENDIF
               ENDIF
            ! comA contains integration limits :
            CALL NiceNum (wIntL, h1, ih1)
            CALL NiceNum (wIntU, h2, ih2)
            CALL NiceNum (dwInt, h3, ih3)
            IF (facIntNor.ne.1.d0) THEN
               CALL NiceNum (facIntNor, h4, ih4)
               h4 = h4(1:ih4)//' * int '
               ih4 = ih4 + 6
            ELSE
               h4 = 'int '
               ih4 = 4
               ENDIF
            comA = h4(1:ih4)//h1(1:ih1)//'..'//h2(1:ih2)//
     *             '+-'//h3(1:ih3)
         ELSEIF (iMod.eq. 5) THEN
            CALL Compose2 (CoG, coY, '_max')
            UnG = unY
            comA = 'y_max'
         ELSEIF (iMod.eq. 6) THEN
            CALL Compose2 (CoG, coY, '_min')
            UnG = unY
            comA = 'y_min'
         ELSEIF (iMod.eq. 7) THEN
            CALL Compose2 (CoG, coX, '[max]')
            UnG = unX
            comA = 'x_max'
         ELSEIF (iMod.eq. 8) THEN
            CALL Compose2 (CoG, coX, '[min]')
            UnG = unX
            comA = 'x_min'
         ELSEIF (iMod.eq. 9) THEN
            CALL Compose2 (CoG, 'match('//coZ, ')')
            UnG = ' '
            comA = '< y(x;K)/y(x;1) >'
            IF (qSameRange) THEN
               IF (qFullRange) THEN
                  CALL Append (comA, '_all x')
               ELSE
                  CALL NiceNum (rmi, h1, ih1)
                  CALL NiceNum (rma, h2, ih2)
                  CALL Append (comA, ' x='//h1(1:ih1)//'..'//h2(1:ih2))
                  ENDIF
            ELSE
               CALL Append (comA, ' indiv. ranges')
               ENDIF
         ELSEIF (iMod.eq.10) THEN
            CALL Compose2 (CoG, 'quality[match]('//coZ, ')')
            UnG = ' '
            comA = '< (y(K)-y(K-1))^2 / y^p >'
            IF (qSameRange) THEN
               IF (qFullRange) THEN
                  CALL Append (comA, '_all x')
               ELSE
                  CALL NiceNum (rmi, h1, ih1)
                  CALL NiceNum (rma, h2, ih2)
                  CALL Append (comA, ' x='//h1(1:ih1)//'..'//h2(1:ih2))
                  ENDIF
            ELSE
               CALL Append (comA, ' indiv. ranges')
               ENDIF
         ELSEIF (iMod.eq.11) THEN
            CALL Compose2 (CoG, 'Match('//coZ, ')')
            UnG = ' '
            comA = '< y(x;K)/y''(x;K) >'
            IF (qFullRange) THEN
               CALL Append (comA, '_all x')
            ELSE
               CALL NiceNum (rmi, h1, ih1)
               CALL NiceNum (rma, h2, ih2)
               CALL Append (comA, ' x='//h1(1:ih1)//'..'//h2(1:ih2))
               ENDIF
         ELSEIF (iMod.eq.13) THEN
            CALL Compose3 (CoG, 'fwhm[', coY, ']')
            UnG = unX
            comA = 'fwhm'
         ELSEIF (iMod.eq.15) THEN
            CALL Compose2 (CoG, 'Sum ', coY)
            UnG = unY
            comA = 'sum'
         ELSEIF (iMod.eq.16) THEN
            CALL Compose3 (CoG, '<', coY, '>')
            UnG = unY
            comA = 'average'
            ENDIF

C  Checks and preliminary calculations :
         IF (iMod.eq.9) THEN
            ! Checks :
            IF (Knorm.gt.nK) THEN
               Fehler = ' cannot normalize, there are only '//cl2(nK)
               CALL Append (Fehler, ' spectra')
               RETURN
               ENDIF
            ! Common x - scale of input spectra :
            CALL CommonScale (
     *         j, .true., qComN, nCom, qCom2, n2, X2, iShift, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (.not.qCom2) THEN
               Fehler =
     * ' no common scale - the rel. weight is not defined'
               RETURN
               ENDIF
            ENDIF

C  Loop spectra :
         KK = 0 ! not every input spectrum must yield an output point {2feb00}
         DO K = 1, nK
            KK = KK + 1
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF     (iMod.eq. 1) THEN
               CALL OlfGet1Z (j, K, iZ2, Y1(KK), Fehler)
               IF (Fehler.ne.'&ff') RETURN
            ELSEIF (iMod.eq. 2) THEN
               IF (.not.qiinside (iX0, 1, n)) THEN
                  Print *,
     * 'channel number out of range in spectrum '//cl4(K)
                  KK = KK - 1
                  GOTO 199
                  ENDIF
               Y1(KK) = X(iX0)
               D1(KK) = 0.
            ELSEIF (iMod.eq. 3) THEN
               IF (iX0.eq.0) THEN
                  IF (rX0.lt.X(1) .or. rX0.gt.X(n)) THEN
                     Print *,
     * 'channel not in x-range in spectrum '//cl4(K)
                     KK = KK - 1
                     GOTO 199
                     ENDIF
                  iX = irPosOpt (X, n, rX0, 'n', iX)
               ELSE
                  IF (.not.qiinside (iX0, 1, n)) THEN
                     Print *,
     * 'channel number out of range in spectrum '//cl4(K)
                     KK = KK -1
                     GOTO 199
                     ENDIF
                  iX = iX0
                  ENDIF
               Y1(KK) = Y(iX)
               D1(KK) = D(iX)
            ELSEIF (iMod.eq. 4 .or. iMod.eq.14) THEN
               IF (modIntWgt.ne.0) THEN
                  DO i = 1, n
                     IF (X(i).eq.0) THEN
                        Y(i) = 0
                     ELSE
                        Y(i) = X(i)**modIntWgt * Y(i)
                        ENDIF
                     ENDDO
                  ENDIF
               CALL Integrate (X, Y, D, n, modIntDif, wIntL, wIntU,
     *            Y1(KK), D1(KK), Fehler)
               Y1(KK) = Y1(KK) * facIntNor
               D1(KK) = D1(KK) * facIntNor
               IF (Fehler.ne.'&ff') THEN
                  CALL Say2 (Fehler, ' in spectrum '//cl4(K))
                  Fehler = '&ff'
                  KK = KK - 1
                  GOTO 199
                  ENDIF
            ELSEIF (iMod.eq. 5) THEN
               CALL IdaMinMax (.false., MiMaMo, nMiMaGru, X, Y, D, n,
     *             ixdummy, xdummy, sdummy, Y1(KK), D1(KK), Fehler)
            ELSEIF (iMod.eq. 6) THEN
               CALL IdaMinMax (.true., MiMaMo, nMiMaGru, X, Y, D, n,
     *             ixdummy, xdummy, sdummy, Y1(KK), D1(KK), Fehler)
            ELSEIF (iMod.eq. 7) THEN
               CALL IdaMinMax (.false., MiMaMo, nMiMaGru, X, Y, D, n,
     *             ixdummy, Y1(KK), D1(KK), ydummy, ddummy, Fehler)
            ELSEIF (iMod.eq. 8) THEN
               CALL IdaMinMax (.true., MiMaMo, nMiMaGru, X, Y, D, n,
     *             ixdummy, Y1(KK), D1(KK), ydummy, ddummy, Fehler)
            ELSEIF (qiinside(iMod, 9, 12)) THEN
               ! relative weight : 11feb91
               ! 24feb/11mar (K/K-1) - 10jul91 changed back to (K-1/K)
               IF (iMod.le.10 .and. K.eq.1) THEN
                  IF     (iMod.eq. 9) THEN
                     Y1(KK) = 1.
                     D1(KK) = 0.
                  ELSEIF (iMod.eq.10) THEN
                     Y1(KK) = 0.
                     D1(KK) = 0.
                     ENDIF
               ELSE
                  IF (iMod.le.10) THEN
                     j2 = j
                     K2 = K-1
                  ELSE
                     K2 = IzK(K)
                     ENDIF
                  CALL OlfGetXY (j2, K2, n2, X2, Y2, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  nrel = 0
                  yZaehler = 0.  ! for iMod  9
                  yNenner  = 0.  !      "
                  yMatch   = 0.  ! for iMod 10
                  ! ask for range ?
                  IF (.not.qSameRange) THEN
                     CALL Compose2 (aus,
     *                    ' x-Range for match of spectra '//cl4(K),
     *                    ' and '//cl4(K2))
                     CALL rAskRgeFull (aus, rmi, rma, rmi, rma)
                     qFullRange = (rmi.eq.0. .and. rma.eq.0.)
                     ENDIF
                  DO i = 1, n
                        ! noch nich auf iShift umgestellt
                     IF (qrinside(X(i), X2(1), X2(n2)).and.
     *                   qrInRange(X(i),rmi, rma)           ) THEN
                        i2 = irPosOpt (X2, n2, X(i), 'n', i2+1)
                        nrel = nrel + 1
                           ! new Ansatz 9/10jul91 :
                           ! least squares (Y-(yZaehler/yNenner)*Y2)**2
                        yZaehler = yZaehler + Y2(i2)*Y (i)
                        yNenner  = yNenner  + Y2(i2)*Y2(i2)
                        yMatch   = yMatch   + dquot0(
     *                       (Y(i)-Y2(i2))   **2, ((Y(i)+Y2(i2))/2)**2)
                        ENDIF
                     ENDDO
                  IF (nrel.eq.0) THEN
                     Print *, 'no overlap of spectra'//cr2(K-1)//'f.'
                     KK = KK - 1
                     GOTO 199
                     ENDIF
                  IF     (iMod.eq. 9) THEN
                     Y1(KK) = Y1(KK-1) * dquot0 (yZaehler, yNenner)
                     D1(KK) = 0. ! Fehlerrechnung unklar
                  ELSEIF (iMod.eq.11) THEN
                     Y1(KK) = dquot0 (yZaehler, yNenner)
                     D1(KK) = 0.
                  ELSEIF (iMod.eq.10 .or. iMod.eq.12) THEN
                     Y1(KK) = yMatch / nrel
                     D1(KK) = 0. ! Fehlerrechnung erst recht unklar
                     ENDIF
                  ENDIF
            ELSEIF (iMod.eq.13) THEN ! fwhm 15dez97
               D1(KK) = 0
               Y1(KK) = 0
               CALL IdaMinMax (.false., 1, 1, X, Y, D, n,
     *             imax, xdummy, sdummy, ymax, ddummy, Fehler)
               DO ii = imax-1, 1, -1
                  IF (Y(ii).lt.ymax/2) GOTO 131
                  ENDDO
               GOTO 139
 131           CONTINUE
               DO if = imax+1, n
                  IF (Y(if).lt.ymax/2) GOTO 133
                  ENDDO
               GOTO 139
 133           CONTINUE
               CALL LinIntPol1 (ymax/2, Y(ii+1), Y(ii), X(ii+1),
     *                          X(ii), xi)
               CALL LinIntPol1 (ymax/2, Y(if-1), Y(if), X(if-1),
     *                          X(if), xf)
               Y1(KK) = xf - xi
 139           CONTINUE
            ELSEIF (iMod.eq.15) THEN ! sum
               na = 0
               ya = 0
               da = 0
               DO i = 1, n
                  IF (qrInRange(X(i),wIntL,wIntU)) THEN
                     na = na + 1
                     ya = ya + Y(i)
                     da = da + D(i)**2
                     ENDIF
                  ENDDO
               IF (na.le.0) THEN
                  Print *, 'sum over empty range in spectrum '//cl4(K)
                  GOTO 199
                  ENDIF
               Y1(KK) = ya
               D1(KK) = dsqrt0(da)
            ELSEIF (iMod.eq.16) THEN ! average
               na = 0
               ya = 0
               da = 0
               DO i = 1, n
                  IF (qrInRange(X(i),wIntL,wIntU)) THEN
                     na = na + 1
                     ya = ya + Y(i)
                     da = da + D(i)**2
                     ENDIF
                  ENDDO
               IF (na.le.0) THEN
                  Print *,
     * 'average over empty range in spectrum '//cl4(K)
                  GOTO 199
                  ENDIF
               Y1(KK) = ya / na
               D1(KK) = dsqrt0(da/na)
            ELSE
               CALL Absturz ('OprInt', 'iMod oor')
               ENDIF ! iMod
            IF (Fehler.ne.'&ff') RETURN
 199        CONTINUE
            ENDDO

         IF (KK.le.0) THEN
            Fehler = 'no points on output'
            RETURN
            ENDIF

C  Normalize :
         IF (iMod.eq.9) THEN
            y1norm = Y1(Knorm)
            IF (y1norm.eq.0.) THEN
               Fehler = ' Normalization channel has result 0.0'
               RETURN
               ENDIF
            DO K = 1, KK
               Y1(K) = Y1(K)/y1norm
               ENDDO
            ENDIF

C  Save :
         IF     (modSav.eq.1) THEN ! the traditional way: save as file
            CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
             IF (Fehler.ne.'&ff') RETURN
            CALL iOlfP (jout, 'plot-sy#', 0, Fehler)
            CALL OlfCnuP (jout, 'y', CoG, UnG, Fehler)
            CALL OlfComAdd (jout, '_'//cl2(iMod), 'int/ y='//
     *                      comA, Fehler)
             IF (Fehler.ne.'&ff') RETURN
            CALL TensorSaveInt (j, jout, KK, Kout, Y1, D1, Fehler)
             IF (Fehler.ne.'&ff') RETURN
            CALL OlfClos (jout, Kout, Fehler)
             IF (Fehler.ne.'&ff') RETURN
         ELSEIF (modSav.eq.2) THEN ! save as z/r (14sep98)
            IF (nK.eq.1) THEN
               CALL rOlfP (j, CoG, UnG, Y1(1), Fehler)
                IF (Fehler.ne.'&ff') RETURN
            ELSE ! nK>1
               CALL OlfOpen (j, 1, K, Fehler)
                IF (Fehler.ne.'&ff') RETURN
               nZ = iOlfG (j, '#Z', Fehler)
               CALL OlfCnuP (j, 'z+', CoG, UnG, Fehler)
                IF (Fehler.ne.'&ff') RETURN
               DO K = 1, KK
                  CALL OlfGetZ (j, K, nZ, Z, Fehler)
                   IF (Fehler.ne.'&ff') RETURN
                  Z(nZ+1) = Y1(K)
                  CALL OlfPutZ (j, K, nZ+1, Z, Fehler)
                   IF (Fehler.ne.'&ff') RETURN
                  ENDDO
               CALL OlfClos (j, KK, Fehler)
               ENDIF
         ELSEIF (modSav.eq.3) THEN ! display
            IF (nK.eq.1) THEN
               CALL NiceNum (Y1(1), h1, ih1)
               CALL Say3 (' result : '//CoG, ' = '//h1, ' '//UnG)
            ELSE
               CALL Say3 (' result: K -> '//CoG, ' ('//UnG, ')')
               DO K = 1, KK
                  Print '(2x,i3,g14.5)', K, Y1(K)
                  ENDDO
               ENDIF
            ENDIF

         ENDDO ! lj (loop files)

      END ! OprIntegral

C  --------------------------------------------------------------------
      SUBROUTINE OprDifferential (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! y -> dy/dx u.a.
         ! JWu 17apr91. Integration 19feb92.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*) Fehler
      DIMENSION     JList(*)
      CHARACTER*40  coX, coY, unX, unY, CoG, CoA, UnG, UnA, h1

C  Files :
      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Operation modes :
      Print *, ' Mappings y(x) -> '
      Print *, '    (1) delta x           (2) delta y'
      Print *, '    (3) dy / dx           (4) I_a^x dx'' y(x'')'
      Print *, '    (5) y*f(x)            (6) y/f(x)'
      iD = iAskDMu (' Option', iD, 0, 6)
      IF (iD.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Further input :
      IF (iD.eq.4) THEN
         xLow = rAskD (' Lower integration boundary', xLow)
      ELSEIF (iD.eq.5 .or. iD.eq.6) THEN
         iFx = iAskDMu (' Function f(x) no.', iFx, 0, 49)
         ENDIF

C  Loop Files :
      DO lj = 1, nJList
      j = JList(lj)

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Documentation :
         IF     (iD.eq.1) THEN
            CALL OlfComAdd (jout, 'd', 'y : = delta x', Fehler)
            CoG = 'delta '//coX
            UnG = unX
         ELSEIF (iD.eq.2) THEN
            CALL OlfComAdd (jout, 'd', 'y : = delta y', Fehler)
            CoG = 'delta '//coY
            UnG = unY
         ELSEIF (iD.eq.3) THEN
            CALL OlfComAdd (jout, 'd', 'y : = dy/dx', Fehler)
            CALL Compose2 (CoG,  'd'//coY, '/d'//coX)
            CALL Compose2 (UnG, unY, '/'//unX)
         ELSEIF (iD.eq.4) THEN
            CALL NiceNum (x0, h1, ih1)
            CALL OlfComAdd (jout, 'i',
     *         'y := Int '//h1(1:ih1)//'..x of y', Fehler)
            CALL Compose3 (CoG, 'I{'//coY, '}('//coX, ')')
            CALL Compose2 (h1, unX,'-1')
            IF (unY.eq.h1) THEN
               UnG = ' '
            ELSE
               CALL Compose2 (UnG, unY, '*'//unX)
               ENDIF
         ELSEIF (iD.eq.5) THEN
            CALL IdaFuTxt (iFx, h1, 'x', '&exclu')
            IF (h1.eq.'&undefined') THEN
               Fehler = ' chosen function f(x) is undefined'
               RETURN
               ENDIF
            CALL OlfComAdd (jout, 'f', 'y := y * '//h1, Fehler)
            CALL IdaFuTxt (iFx, h1, coX, '&exclu')
            CALL Compose2 (CoG, coY, '*'//h1)
            CALL IdaFuTxt (iFx, h1, unX, '&exclu')
            CALL Compose2 (UnG, unY, '*'//h1)
         ELSEIF (iD.eq.6) THEN
            CALL IdaFuTxt (iFx, h1, 'x', '&exclu')
            IF (h1.eq.'&undefined') THEN
               Fehler = ' chosen function f(x) is undefined'
               RETURN
               ENDIF
            CALL OlfComAdd (jout, 'f', 'y := y / '//h1, Fehler)
            CALL IdaFuTxt (iFx, h1, coX, '&exclu')
            CALL Compose2 (CoG, coY, '/'//h1)
            CALL IdaFuTxt (iFx, h1, unX, '&exclu')
            CALL Compose2 (UnG, unY, '/'//h1)
            ENDIF

         CALL OlfCnuP (jout, 'y', CoG, UnG, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (iD.le.4) THEN
               IF (irSorted(X,n).ne.2) THEN
                  Fehler = ' spectrum'//cl3(K)//' is not sorted'
                  RETURN
                  ENDIF
               ENDIF

C  Operate :
            IF (iD.eq.1 .or. iD.eq.2 .or. iD.eq.3) THEN
               n = n - 1
               DO i = 1, n
                  X1(i) = (X(i+1) + X(i)) / 2
                  ENDDO ! i
               IF     (iD.eq.1) THEN
                  DO i = 1, n
                     Y1(i) = X(i+1) - X(i)
                     D1(i) = 0
                     ENDDO ! i
               ELSEIF (iD.eq.2) THEN
                  DO i = 1, n
                     Y1(i) = Y(i+1) - Y(i)
                     D1(i) = dsqrt (Y(i+1)**2 + Y(i)**2)
                     ENDDO ! i
               ELSEIF (iD.eq.3) THEN
                  DO i = 1, n
                     Y1(i) = (Y(i+1) - Y(i)) / (X(i+1) - X(i))
                     D1(i) = dsqrt (Y(i+1)**2 + Y(i)**2) /
     *                             (X(i+1) - X(i))
                     ENDDO ! i
                  ENDIF
               CALL rCopy (X, 1, n, 1, X1, 1, 1)
            ELSEIF (iD.eq.3) THEN
               ! Quotient of differences :
               n = n - 1
               DO i = 1, n
                  X1(i) = (X(i+1) + X(i)) / 2
                  Y1(i) = (Y(i+1) - Y(i)) / (X(i+1) - X(i))
                  D1(i) = dsqrt (Y(i+1)**2 + Y(i)**2) / (X(i+1) - X(i))
                  ENDDO ! i
               CALL rCopy (X, 1, n, 1, X1, 1, 1)
            ELSEIF (iD.eq.4) THEN
               ! integrate :
               Y1(1) = 0
               D1(1) = 0
               DO i = 2, n
                  Y1(i) = Y1(i-1) + (Y(i-1)+Y(i))/2 * (X(i)-X(i-1))
                  D1(i) = 0. ! unfertich
                  ENDDO
               ! subtract value at lower boundary :
               CALL LinIntPolArray (X, Y1, n, xLow, yLow, .5d0, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               DO i = 1, n
                  Y1(i) = Y1(i) - yLow
                  ENDDO
            ELSEIF (iD.eq.5) THEN
               DO i = 1, n
                  CALL IdaFuVal (iFx, fx, dfx, X(i), 0.d0, 0.d0, 0.d0)
                  Y1(i) = Y(i) * fx
                  D1(i) = D(i) * fx
                  ENDDO
            ELSEIF (iD.eq.6) THEN
               DO i = 1, n
                  CALL IdaFuVal (iFx, fx, dfx, X(i), 0.d0, 0.d0, 0.d0)
                  Y1(i) = dquot0 (Y(i), fx)
                  D1(i) = dquot0 (D(i), fx)
                  ENDDO
               ENDIF

            CALL OlfCopZ   (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n, X, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! K

         CALL OlfClos (jout, nK, Fehler)
         ENDDO
C  End loop files.

      END ! OprDifferential

C  --------------------------------------------------------------------
      SUBROUTINE OprPointwise  (LHS, nJList, JList, qSR, qOv,
     *                          Object, Fehler)
C  --------------------------------------------------------------------
         ! FK Dec 2002 option for q added: takes directly E0 from file
         ! JWu 1990/91. Dialogue revised 21jun91
         ! any pointwise function of x, y, or z, entered by hand.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER    (MRL=10)
      CHARACTER*(*) Fehler, Object, LHS
      DIMENSION     Ix2(MC), IZofKorC(MC), qLiSR(MC), Zz1(MK),
     *              IZXofZ(MZ), JList(*), J2List(MF), JListO(MF),
     *              RL(MRL), Z(MZ)
      CHARACTER*80  aus, zahl
      CHARACTER*40  textRHS, Co1arg, Un1arg, Co2arg, Un2arg,
     *              doc2arg, doc, Co2glo, word, ein, docSR, file2,
     *              Co2arg3, Co2arg4
      CHARACTER     typ2arg*4

      DATA          iFu /-1/, Co2glo /'#'/, iInpModY /1/, iInpModD /2/,
     *              r2arg /0.d0/, d2arg /0.d0/, qSortX /.true./

      IF (qErrEntry('OprPointwise', Fehler)) RETURN
      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Check available space in memory (11mar93) :
      IF (.not.qOv) THEN
         newE = 0
         DO lj = 1, nJList
            CALL OlfGetNofK (JList(lj), nK, NofK, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            newE = newE + iSum (NofK, 1, NK, 1)
            ENDDO
         IF (qMemFill(nJList, newE, Fehler)) RETURN
         ENDIF

C  Initial values :
      tolX = 1.d-6
      tolZ = 1.d-3
      qSortXAsked = .false.

C  operate on which coordinate (left-hand-side) :
      iLHS = iLabOfCnu (LHS)
      IF (iLHS.lt.-1) THEN
         Fehler = 'invalid LHS '//cl2(iLHS)
         RETURN
         ENDIF

C  Common questionary for all files :

C  Function :
      CALL TakeVorDel (Object, word, ' ')
      in = IdaFuAsk (' Function', word, LHS, '#2', iFu, 2)
      IF (in.lt.0) THEN
         Fehler = ' '
         RETURN
         ENDIF
      iFu = in
      CALL IdaFuTxt (iFu, textRHS, LHS, '#2')
      Print *, ' function = '//textRHS

C  Arguments of function :
      IF (iFu.lt.50) THEN
         !  Only one argument :
         i2arg = 0
      ELSE ! renewed dialogue 6sep93.
         ein = Object
 461     CONTINUE
         IF (ein.eq.' ') THEN
            CALL FrageHD ('2nd argument', ein, typ2arg)
            IF (ein.eq.'^[' .or. ein.eq.CHAR(27)) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF
         typ2arg = ein
         IF     (ein.eq.'ec') THEN
            i2arg =  1  ! fix
         ELSEIF (ein.eq.'ef') THEN
            i2arg =  2  ! per file
         ELSEIF (ein.eq.'r') THEN
            i2arg =  3
         ELSEIF (ein.eq.'r2y') THEN
            i2arg =  4  ! y' of rPar
         ELSEIF (ein.eq.'if') THEN
            i2arg =  5  ! y'(1-point) of no-z
         ELSEIF (ein.eq.'es') THEN
            i2arg =  6  ! enter per spectrum
         ELSEIF (ein(1:1).eq.'z') THEN
            i2arg =  7  ! z#
            IF (ein.eq.'z') THEN
               i2arg2 = 1
               typ2arg = 'z1'
            ELSE
               CALL Fi1N (ein, i2arg2)
               IF (ein.ne.'z#') THEN
                  Fehler = ' invalid z-par'
                  RETURN
                  ENDIF
               ENDIF
         ELSEIF (ein.eq.'i') THEN ! a more systematic name would be 'z2y' ...
            i2arg =  8  ! y' of z
         ELSEIF (ein.eq.'i1') THEN
            i2arg =  9  ! y' of z, old version
         ELSEIF (ein.eq.'ep') THEN
            i2arg = 11  ! enter per point
         ELSEIF (ein.eq.'y2') THEN
            i2arg = 13  ! y' of x
         ELSEIF (ein.eq.'d2') THEN
            i2arg = 14  ! dy of x
         ELSEIF (ein.eq.'x') THEN
            i2arg = 15
         ELSEIF (ein.eq.'y') THEN
            i2arg = 16
         ELSEIF (ein.eq.'d') THEN
            i2arg = 17
         ELSEIF (ein.eq.'n') THEN
            i2arg = 18
         ELSEIF (ein.eq.'?' .or. ein.eq.'h') THEN
            Print *, 'INPUT HELP/'
            Print *, '   required input : the kind of 2nd argument'
            Print *, '   here is a list of allowed answers :'
            Print *,
     * '    ec  : one value for all files, enter once'
            Print *,
     * '    ef  : one value per file, enter per file'
            Print *,
     * '    if  : one value per file, take from one-point-file'
            Print *,
     * '    r   : one value per file, take a real parameter'
            Print *,
     * '    r2y : one value per file, take y''(real-par) from 2nd file'
            Print *,
     * '    es  : one value per spectrum, enter per spectrum'
            Print *,
     * '    z<n>: one value per spectrum, take z<n>'
            Print *,
     * '    i   : one value per spectrum, take y''(z) from 2nd file'//
     *                             ' (integral file)'
            Print *,
     * '    i1  : one value per spectrum, take y''(z) from 2nd file'//
     *                             ' (integral file, old 1dim version)'
            Print *,
     * '    ep  : one value per point, enter per point'
            Print *,
     * '    x   : one value per point, take x'
            Print *,
     * '    y   : one value per point, take y'
            Print *,
     * '    d   : one value per point, take d'
            Print *,
     * '    n   : one value per point, take number of point'
            Print *,
     * '    y2  : one value per point, take y''(x) from 2nd file'
            Print *,
     * '    d2  : one value per point, take the error dy''(x)'//
     * ' from 2nd file'
            Print *,
     * '    ^[  : escape'
            ein = ' ' ! damit danach wieder gefragt wird
            typ2arg = ' '
            GOTO 461
         ELSE
            CALL Gong (4)
            Print *, ' Type ''?'' for help'
            ein = ' '
            typ2arg = ' '
            GOTO 461
            ENDIF

         ENDIF

C  Reorganise description of 2nd argument :
      qJdep = (i2arg.ge. 2)
      qKdep = (i2arg.ge. 6)
      qXdep = (i2arg.ge.11)
      q2file = (i2arg.eq.5 .or. i2arg.eq.8 .or. i2arg.eq.9.or.
     *         qiinside(i2arg,12,14))
      q2same = qiinside (i2arg, 15, 18)

      IF (qiinside(i2arg,12,14)) i2col = i2arg-11 ! 2nd column = x,y, or d
      IF (qiinside(i2arg,15,17)) i2col = i2arg-14 ! i2col      = 1,2, or 3
      IF (i2arg.eq.18) i2col = 0
      IF (i2arg.eq.19) CALL Absturz('OprPointwise', 'iarg=19')
      IF (qiinside(i2arg,12,18)) i2arg=19         ! pointwise from spectrum

C  Enter global 2nd argument :
! FK changed dec02
      IF (.not.qJdep) THEN
           IF ((i2arg.eq.1).and.(iFu.eq.81)) THEN
            r2glo = rOlfGG (j,'E0', 'meV', Fehler)
            IF (Fehler(1:4).eq.'&pnf')
     *          Fehler = 'please define r-parameter "E0"'
            IF (LHS.eq.'y') d2glo = rAskD (' Its error', d2glo)
            CALL NiceNum (r2glo, aus, ih2)
            doc2arg = aus
            CALL FrageCD  (' Its name (#=value)', Co2glo, Co2glo)
            CALL FrageCD  (' And its unit', Un2arg, Un2arg)
            IF (Co2glo.ne.'#') THEN
              Co2arg = Co2glo
              CALL Append (doc2arg, ' = '//Co2arg)
              IF (Un2arg.ne.' ') CALL Append(doc2arg, Un2arg)
            ELSE
              Co2arg = aus
              Un2arg = ' '
              ENDIF
            ELSEIF((i2arg.eq.1).and.(iFu.ne.81)) THEN
             r2glo = rAskD  (' The argument', r2glo)
             IF (LHS.eq.'y') d2glo = rAskD (' Its error', d2glo)
             CALL NiceNum (r2glo, aus, ih2)
             doc2arg = aus
             CALL FrageCD  (' Its name (#=value)', Co2glo, Co2glo)
             CALL FrageCD  (' And its unit', Un2arg, Un2arg)
             IF (Co2glo.ne.'#') THEN
               Co2arg = Co2glo
               CALL Append (doc2arg, ' = '//Co2arg)
               IF (Un2arg.ne.' ') CALL Append(doc2arg, Un2arg)
             ELSE
               Co2arg = aus
               Un2arg = ' '
               ENDIF
            ENDIF
C  Prepare 2nd argument per file :
      ELSE
         IF     (i2arg.eq. 2) THEN
            CALL FrageC (' Name of the 2nd argument', Co2arg)
            IF (Co2arg.ne.' ') THEN
               CALL FrageCD  (' And its unit', Un2arg, Un2arg)
            ELSE
               Un2arg = ' '
               ENDIF
         ELSEIF (i2arg.eq. 3) THEN
            CALL FrageCD (' 2nd argument is r-parameter',
     *                    Co2arg3, Co2arg3)
            Co2arg = Co2arg3
            Un2arg = '[rPar]' ! c96/7
            doc2arg  = Co2arg
            aus = '  -> argument = '//Co2arg
            CALL Append (aus, Un2arg)
            Print *, aus
         ELSEIF (i2arg.eq. 4) THEN ! JWu 15oct91
            CALL FrageCD (
     *        ' 2nd argument y''(r) where r is r-parameter',
     *        Co2arg4, Co2arg4)
            j2 = iAskD (' 2nd argument from file no.', j2)
            nK2 = iOlfG (j2, '#spectra', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (nK2.gt.1) THEN
               aus = ' Take which of the '//cl3(nK2)
               CALL Append (aus, ' spectra in the 2nd file')
               K2 = iAskDMu (aus, mod(K2,nK2)+1, 0, nK2)
               IF (K2.eq.0) RETURN
            ELSE
               K2 = 1
               ENDIF
            CALL OlfGetXYD (j2, K2, n2, X2, Y2, D2, Fehler)
            CALL OlfCnuG (j2, 'y', Co2arg, Un2arg, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (iabs(irSorted(X2,n2)).lt.2) THEN
               Fehler = ' 2nd file is not sorted'
               RETURN
               ENDIF
            CALL tOlfG (j2, 'fil', file2, Fehler)
            CALL Compose2 (doc2arg, 'y='//Co2arg4, ' from '//file2)
            ENDIF
         ENDIF

C  Files for 2nd argument :
      IF (q2file) THEN
         aus = ' 2nd Argument'
         CALL GetJ2J (aus, JList, nJList, J2List, iKmod, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF

C  Loop over files :
      nJListO = 0
      DO lj = 1, nJList
      j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (j, LHS, Co1arg, Un1arg, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF     (q2file) THEN
            j2 = J2List(lj)
         ELSEIF (q2same) THEN
            j2 = j
            ENDIF

         IF (.not.qJdep) THEN
            r2arg = r2glo
            IF (LHS.eq.'y') d2arg = d2glo
         ELSE
            IF (i2arg.eq. 2) THEN
               r2arg = rAskD  (' The 2nd argument', r2arg)
               IF (LHS.eq.'y') d2arg = rAskD (' Its error', d2arg)
               CALL NiceNum (r2arg, aus, ih2)
               doc2arg = aus(1:ih2)
               IF (Co2arg.ne.aus) THEN
                  CALL Append (doc2arg, ' = '//Co2arg)
                  CALL Append (doc2arg, Un2arg)
               ELSE
                  Un2arg = ' '
                  ENDIF
            ELSEIF (i2arg.eq. 3) THEN
               r2arg = rOlfG (j, Co2arg3, Un2arg, Fehler)
               d2arg = 0.
            ELSEIF (i2arg.eq. 4) THEN
               CALL OlfGet1ZofK (j, 1, nK, Zz1, Fehler)
               Zz1(1) = rOlfG (j, Co2arg4, Un2arg, Fehler)
               CALL GetIndex (aus, Zz1, 1, X2, n2, IZofKorC,
     *                        tolZ, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               r2arg = Y2(IZofKorC(1))
               d2arg = D2(IZofKorC(1))
            ELSEIF (i2arg.eq. 5) THEN ! ip {2feb00}
               IF (nK.ne.1) THEN
                  Fehler = '2nd arg "if" requires 1 spectrum / file'
                  RETURN
                  ENDIF
               nK2 = iOlfG (j2, '#spectra', Fehler)
               IF (nK2.ne.1) THEN
                  Fehler = 'one-point-file has several specctra'
                  RETURN
                  ENDIF
               CALL OlfGetXYD (j2, 1, n2, X2, Y2, D2, Fehler)
               IF (n2.ne.1) THEN
                  Fehler = '2nd file contains not just one point'
                  RETURN
                  ENDIF
               r2arg = Y2(1)
               d2arg = D2(1)
               ! doc :
               CALL OlfCnuG (j2, 'y', Co2arg, Un2arg, Fehler)
               CALL tOlfG (j2, 'fil', file2, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL Compose2 (doc2arg, 'y='//Co2arg, ' from '//file2)
            ELSEIF (i2arg.eq. 6) THEN
               CALL FrageC (' Name of the argument', Co2arg)
               CALL FrageCD  (' And its unit', Un2arg, Un2arg)
               doc2arg = Co2arg
            ELSEIF (i2arg.eq. 7) THEN
               CALL OlfCnuG (j, typ2arg, Co2arg, Un2arg, Fehler)
               doc2arg  = Co2arg
            ELSEIF (i2arg.eq. 8) THEN
               CALL SetZXofZ (j, j2, IZXofZ, iz2ofx, nZ2, Fehler)
               IF (iz2ofx.le.0) THEN ! nicht weiter zurueckverfolgt, 7mai98
                  Fehler = 'SetZXofZ returns invalid iz2ofx='//
     *                     cl2(iz2ofx)
                  RETURN
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN
               nK2 = iOlfG (j2, '#spectra', Fehler) ! necessary for case nZ2=0
               IF (Fehler.ne.'&ff') RETURN
               DO iZ2 = 1, nZ2
                  CALL OlfGet1ZofK (j2, iZ2, nK2, ZZofK(1,iZ2), Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDDO
               ! doc :
               CALL OlfCnuG (j2, 'y', Co2arg, Un2arg, Fehler)
               CALL tOlfG (j2, 'fil', file2, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL Compose2 (doc2arg, 'y='//Co2arg, ' from '//file2)
            ELSEIF (i2arg.eq. 9) THEN
               nK2 = iOlfG (j2, '#spectra', Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (nK2.gt.1) THEN
                  aus = ' Take which of the '//cl3(nK2)
                  CALL Append (aus, ' spectra in the 2nd file')
                  K2 = iAskDMu (aus, mod(K2,nK2)+1, 0, nK2)
                  IF (K2.eq.0) RETURN
               ELSE
                  K2 = 1
                  ENDIF
               CALL OlfGetXYD (j2, K2, n2, X2, Y2, D2, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (iabs(irSorted(X2,n2)).lt.2) THEN
                  Fehler = ' 2nd file is not sorted'
                  RETURN
                  ENDIF
               CALL OlfCnuG (j2, 'y', Co2arg, Un2arg, Fehler)
               CALL tOlfG (j2, 'fil', file2, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL Compose2 (doc2arg, 'y='//Co2arg, ' from '//file2)
               aus = ' getting z'
               CALL OlfGet1ZofK (j, 1, nK, Zz1, Fehler)
               CALL GetIndex (aus, Zz1, nK, X2, n2, IZofKorC,
     *                        tolZ, Fehler)
               IF (Fehler.ne.'&ff') RETURN
            ELSEIF (i2arg.eq.11) THEN
               CALL FrageC (' Name of the argument', Co2arg)
               CALL FrageCD  (' And its unit', Un2arg, Un2arg)
                  ! option choice introduced 10jan95 :
               Print *, ' Input mode for y :'
               Print *, '   (1) enter individually'
               Print *, '   (2) y from multicolumn'
               Print *, '   (3) y-d from multicolumn'
               iOpt = iAskDMu (' Option', iInpModY, 0, 3)
               IF (iOpt.eq.0) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               iInpModY = iOpt
               IF     (iInpModY.eq.2) THEN
                  IF (iMulColY.le.0) iMulColY = 1
                  iMulColY = iAskDMu (' Read y from column', iMulColY,
     *                                0, MRL)
                  IF (iMulColY.le.0) THEN
                     Fehler = ' '
                     RETURN
                     ENDIF
               ELSEIF (iInpModY.eq.3) THEN
                  IF (iMulColY.le.0) iMulColY = 1
                  IF (iMulColD.lt.0) iMulColD = 0
                  CALL i2FrageD (' Read y,d from columns',
     *               iMulColY, iMulColD, iMulColY, iMulColD)
                  IF     (iMulColY.le.0 .or. iMulColD.lt.0) THEN
                     Fehler = ' '
                     RETURN
                  ELSEIF (iMulColY.gt.MRL .or. iMulColD.gt.MRL) THEN
                     Fehler = ' max. number of columns exceeded'
                     RETURN
                     ENDIF
                  ENDIF
               IF (iInpModY.ne.3) THEN
                  Print *, ' Input mode for d :'
                  Print *, '   (1) enter individually'
                  Print *, '   (2) common value 0.0'
                  iOpt = iAskDMu (' Option', iInpModD, 0, 2)
                  IF (iOpt.eq.0) THEN
                     Fehler = ' '
                     RETURN
                     ENDIF
                  iInpModD = iOpt
                  ENDIF

            ELSEIF (i2arg.eq.19) THEN
               qOneOne = .true.  ! option not 1:1 ausser Betrieb
               IF (qOneOne) THEN
                  ! get the z pairs : (20jan92)
                  CALL GetK2K (j, j2, nKdummy, IZofKorC, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               ! get the name of the 2nd argument :
               IF    (i2col.eq.1) THEN
                  CALL OlfCnuG (j2, 'x', Co2arg, Un2arg, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  doc2arg  = 'x'
               ELSEIF (i2col.eq.2) THEN
                  CALL OlfCnuG (j2, 'y', Co2arg, Un2arg, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  doc2arg  = 'y'
               ELSEIF (i2col.eq.3) THEN
                  CALL OlfCnuG (j2, 'y', Co2arg, Un2arg, Fehler)
                  Co2arg = 'd'//Co2arg ! dy : option added 2jul91
                  doc2arg  = 'dy'
               ELSEIF (i2col.eq.0) THEN
                  doc2arg  = 'n'
                  Co2arg = 'no. of pt.'
                  Un2arg = ' '
               ELSE
                  Print *, ' >> i2c n2a : ', i2col, Co2arg
                  CALL Absturz('OprPointwise', 'i2col o.o.r.')
                  ENDIF
               CALL Append (doc2arg, '='//Co2arg)
               IF     (q2same) THEN
                  CALL Append (doc2arg, ' from same file')
               ELSE
                  CALL tOlfG (j2, 'fil', file2, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  CALL Compose2 (doc2arg, doc2arg, ' from '//file2)
                  ENDIF
            ELSE
               Fehler = ' option not yet implemented'
               RETURN
               ENDIF
            ENDIF ! 2nd argument j-dependent

C  Subrange ?
         IF (qSR) CALL SelectChSet ('Operate on ', j, lj,
     *                              docSR, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  - Co-ordinate name and unit :
         CALL IdaFuCoord (iFu, LHS, Co1arg, Un1arg,
     *      Co1arg, Un1arg, Co2arg, Un2arg, Fehler) ! -> IDA3
         CALL OlfCnuP (jout, LHS, Co1arg, Un1arg, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Documentation :
         CALL IdaFuTxt (iFu, textRHS, LHS, '#')
         CALL Compose2 (doc, LHS, '='//textRHS)
         IF (iFu.ge.50) CALL Append (doc, ', #='//doc2arg)
         IF (qSR) CALL Append (doc, ' only for '//docSR)
         CALL OlfComAdd (jout, 'o', doc, Fehler)

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetSpe (j, K, nZ, Z, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (LHS.eq.'x') THEN
               iXwasSorted = irSorted (X, n)
               ENDIF

            IF (iLHS.gt.nZ) THEN
               Fehler = 'cannot operate on '//LHS//
     *                  ' since nZ only '//cl2(nZ)
               RETURN
               ENDIF

C  - Determine subrange :
            IF (qSR) CALL SelectChGet (K, n, X, Y, D, qLiSR, Fehler)
            IF (Fehler.ne.'&ff') RETURN

C  - Number of values to manipulate :
            IF (LHS.eq.'x' .or. LHS.eq.'y') THEN
               nLHS = n ! manipulate x or y
            ELSE
               nLHS = 1 ! manipulate z
               ENDIF

            IF (qKdep) THEN
C  - Determine second argument :
               IF     (i2arg.eq. 6) THEN
                  CALL NiceNum (Z(1), zahl, izahl)
                  CALL Compose3 (aus, ' spectrum '//cl3(K), ' ,z1=',
     *                           '='//zahl(1:izahl)//' :')
                  Print *, aus
                  r2arg = rAskD ('    The 2nd argument', r2arg)
                  IF (LHS.eq.'y')
     *                d2arg = rAskD ('    And its error', d2arg)
               ELSEIF (i2arg.eq. 7) THEN
                  r2arg = Z(i2arg2)
                  d2arg = 0.
               ELSEIF (i2arg.eq. 8) THEN
                  CALL GetKofKbyZX (nZ2, nK2, Z, nZ, IZXofZ,
     *                              K2, Fehler)
                  IF (Fehler.ne.'&ff') THEN
                     CALL Append (Fehler, ' '//cl3(K))
                     RETURN
                     ENDIF
                  CALL OlfGetXYD (j2, K2, n2, X2, Y2, D2, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  IF (irSorted(X2,n2).lt.2) THEN
                     Fehler = 'integral file not sorted'
                     RETURN
                     ENDIF
                  i2 = irPos (X2, n2, Z(iz2ofx), 'n')
                  r2arg = Y2(i2)
                  d2arg = D2(i2)
               ELSEIF (i2arg.eq. 9) THEN
                  r2arg = Y2(IZofKorC(K))
                  d2arg = D2(IZofKorC(K))
               ELSEIF (i2arg.eq.11) THEN
                  CALL Say2 (' enter 2nd argument for spectrum '//
     *                       cl3(K), ' :')
                  ifrom = 1
 5111             CONTINUE
                  DO i = ifrom, n
                     IF (.not.qSR .or. qLiSR(i)) THEN ! channel i not excluded
                        qPanne = .false.
                        write (aus,'(a,i4,a,g12.5,a)')
     *                     ' channel', i, ', x = ', X(i), ' --> '
                        CALL FrageC (aus, ein)
                        IF (iInpModY.eq.1) THEN
                           CALL Fi1R (ein, Y2(i))
                           IF (ein.ne.'#') qPanne = .true.
                        ELSEIF (iInpModY.eq.2 .or. iInpModY.eq.3) THEN
                           CALL FindR (ein, MRL, nRL, RL, .false.)
                           IF (nRL.ge.iMulColY) THEN
                              Y2(i) = RL(iMulColY)
                           ELSE
                              qPanne = .true.
                              ENDIF
                           IF (iInpModY.eq.3) THEN
                              IF (nRL.ge.iMulColD) THEN
                                 D2(i) = RL(iMulColD)
                              ELSE
                                 qPanne = .true.
                                 ENDIF
                              ENDIF
                           ENDIF ! iInpModY
                        ENDIF ! i is in subrange
                     ENDDO ! channels i..
                  IF (qPanne) THEN
                     Print *, ' Bad input - what now ?'
                     Print *, '   (0) break'
                     Print *, '   (1) repeat from channel 1'
                     Print *, '   (2) one channel back'
                     iOpt = iAskMu (' Option', 0, 2)
                     IF     (iOpt.eq.0) THEN
                        Fehler = ' '
                        RETURN
                     ELSEIF (iOpt.eq.1) THEN
                        ifrom = 1
                        GOTO 5111
                     ELSEIF (iOpt.eq.2) THEN
                        ifrom = max0 (1, i-2)
                        GOTO 5111
                        ENDIF
                     ENDIF
                  IF (iInpModY.ne.3) THEN
                     IF     (iInpModD.eq.1) THEN
                        Print *, ' Now enter the error :'
                        ifrom = 1
 5112                   CONTINUE
                        DO i = ifrom, n ! es fehlt was, um subrange auszulassen
                           IF (.not.qSR .or. qLiSR(i)) THEN ! i not excluded
                              qPanne = .false.
                              write (aus,'(a,i4,a,2(g12.5,a))')
     *        ' channel', i, ', x = ', X(i), ', y = ', Y2(i), ' --> '
                              CALL FrageC (aus, ein)
                              CALL Fi1R (ein, Y2(i))
                              ENDIF
                           ENDDO
                        IF (qPanne) THEN
                           Print *, ' Bad input - what now ?'
                           Print *, '   (0) break'
                           Print *, '   (1) repeat from channel 1'
                           Print *, '   (2) one channel back'
                           iOpt = iAskMu (' Option', 0, 2)
                           IF     (iOpt.eq.0) THEN
                              Fehler = ' '
                              RETURN
                           ELSEIF (iOpt.eq.1) THEN
                              ifrom = 1
                              GOTO 5112
                           ELSEIF (iOpt.eq.2) THEN
                              ifrom = max0 (1, i-2)
                              GOTO 5112
                              ENDIF
                           ENDIF
                     ELSEIF (iInpModD.eq.2) THEN
                        CALL rSet (D, 1, n, 1, 0.d0)
                        ENDIF
                     ENDIF

               ELSEIF (i2arg.eq.19) THEN
                  IF (qOneOne) THEN
                     j2r = j2
                     K2r = IZofKorC(K)
                  ELSE
                     aus = ' 2nd argument for spectrum '//cl3(K)
                     CALL Append (aus, ' (file,spectrum)')
                     CALL i2FrageD (aus, j2r, K2r, j2, K2r)
                     j2 = j2r ! default for next K
                     ENDIF
                  CALL OlfGetXYD (j2r, K2r, n2, X2, Y2, D2, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  aus = ' x-scale of spectrum '//cl3(K)
                  IF (i2col.ne.0) THEN
                     CALL GetIndex (aus, X, nLHS, X2, n2, Ix2,
     *                              tolX, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDIF
               ELSE
                  CALL Absturz (' Ida/of', 'i2arg o.o.r.')
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN
               ENDIF ! K dependent

            DO i = 1, nLHS
C  - Determine first argument :
               IF     (LHS.eq.'x') THEN
                  r1arg = X(i)
                  d1arg = 0.
               ELSEIF (LHS.eq.'y') THEN
                  r1arg = Y(i)
                  d1arg = D(i)
               ELSEIF (LHS(1:1).eq.'z') THEN
                  r1arg = Z(iLHS)
                  d1arg = 0.
                  ENDIF

               IF (qXdep) THEN
C  - Determine second argument :
                  IF     (i2arg.eq.11) THEN
                     ! #2 prepared above
                     r2arg = Y2(i)
                     d2arg = D2(i)
                  ELSEIF (i2arg.eq.19 .and. i2col.eq.0) THEN
                     r2arg = i
                     d2arg = 0.
                  ELSEIF (i2arg.eq.19 .and. i2col.eq.1) THEN
                     r2arg = X2(Ix2(i))
                     d2arg = 0.
                  ELSEIF (i2arg.eq.19 .and. i2col.eq.2) THEN
                     r2arg = Y2(Ix2(i))
                     d2arg = D2(Ix2(i))
                  ELSEIF (i2arg.eq.19 .and. i2col.eq.3) THEN
                     r2arg = D2(Ix2(i))
                     d2arg = 0.
                  ELSE
                     CALL Absturz (' ida/mf', 'i2arg o.o.r.')
                     ENDIF
                  ENDIF

               IF (.not.qSR .or. qLiSR(i)) THEN
C  - Functional transform :
                  CALL IdaFuVal (iFu, rres, dres, r1arg, d1arg,
     *                           r2arg, d2arg)
                  IF     (LHS.eq.'x') THEN
                     X(i) = rres
                  ELSEIF (LHS.eq.'y') THEN
                     Y(i) = rres
                     D(i) = dres
                  ELSEIF (LHS(1:1).eq.'z') THEN
                     Z(iLHS) = rres
                  ELSE
                     Fehler = 'SEVERE PROGRAM ERROR/ LHS oor'
                     RETURN
                     ENDIF
                  ENDIF

               ENDDO ! i

            IF (LHS.eq.'x' .and. iXwasSorted.eq.2
     *                     .and. irSorted(X,n).lt.2) THEN
               IF (.not.qSortXAsked) THEN
                  qSortX = qAskD (' Lost order of X - sort data',
     *                        intq(qSortX))
                  qSortXAsked = .true.
                  ENDIF
               IF (qSortX) THEN
                  CALL SortChannels (X, Y, D, n, .false.)
                  ENDIF
               ENDIF

            CALL OlfPutSpe (jout, K, nZ, Z, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! K

         nJListO = nJListO + 1
         JListO (nJListO) = jout
         ENDDO
C  End loop files.

c      CALL EditCnu (nJListO, JListO, LHS, Fehler)

      END ! OprPointwise

C  --------------------------------------------------------------------
      SUBROUTINE OprTensor (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 24apr92
         ! tensor product : (X,Y) x (X,Z) -> (X,Y,Z)

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER    (MRL=10)

      CHARACTER*(*) Fehler
      CHARACTER     czfrom*1
      CHARACTER*40  coZ, unZ, file2
      DIMENSION     JList(*), qLisK(MK), ZZ2(MK)

      CHARACTER*80  aus, cLisK
      CHARACTER*40  com, doc

      DATA          iFun /1/, izfrom /1/, iModZ /11/,
     *              coZ /' '/, unZ /' '/

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      Print *, ' Clone a spectrum y(x) into a file y(x,z).'
      Print *, ' Take z from'
      Print *, '    other file''s x''(11) y''(12) z''(13)'
      Print *,
     * '    regular grid, lin(31) 1/2-log(32) log(33) lin-blocks(34)'
      iModZ = iAskDMu (' Option', iModZ, 0, 50)

      IF     (qiinside(iModZ,11,13)) THEN
         ! Get z file :
         j2 = iAskDMu (' Scale z from file', j2, 0, MF)
         IF (j2.eq.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         nK2 = iOlfG(j2, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF     (iModZ.eq.11) THEN
            izfrom = 1
            czfrom ='x'
         ELSEIF (iModZ.eq.12) THEN
            izfrom = 2
            czfrom ='y'
         ELSEIF (iModZ.eq.13) THEN
            izfrom = 3
            czfrom ='z1'
         ELSE
            Fehler = ' '
            RETURN
            ENDIF
         IF (izfrom.le.2) THEN
            IF (nK2.gt.1) THEN
               Fehler = 'z-file contains more than one spectrum'
               RETURN
               ENDIF
            CALL OlfGetXYD (j2, 1, n2, X2, Y2, D2, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            DO K = 1, n2
               IF (izfrom.eq.1) THEN
                  ZZ2(K) = X2(K)
               ELSE
                  ZZ2(K) = Y2(K)
                  ENDIF
               ENDDO
         ELSE ! this option 30jun92
            CALL OlfGet1ZofK (j2, 1, n2, ZZ2, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDIF

         CALL tOlfG (j2, 'fil', file2, Fehler)
         CALL OlfCnuG (j2, czfrom, coZ, unZ, Fehler)
            Print *, ' DEBUG/ OprTens should transfer Coord Name '//coZ
            Print *, ' DEBUG/ OprTens should transfer Coord Unit '//unZ
         IF (Fehler.ne.'&ff') RETURN

         CALL Append (doc, 'from '//file2)

      ELSEIF (qiinside(iModZ,31,34)) THEN
         CALL FrageCD ('Name of new z-coordinate', coZ, coZ)
         CALL FrageCD ('And its unit', unZ, unZ)
         CALL SetGridReg (iModZ-30, sg1, sgn, nsg, ZZ2, MK, n2, doc,
     *                   Fehler)

      ELSE
         Fehler = ' '
         RETURN
         ENDIF

C  Prepare new z-coordinate (ordentlich 3jun92) :
c96/7 iz  = iParNumber(coZ)
c      IF (iz.ne.0) Print *,
c     *   'New z-coordinate replaces parameter no. '//cl3(iz)

C  Loop in j :
      DO lj = 1, nJList
      j = JList(lj)

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)

         IF (nK.gt.1) THEN
            CALL Say2 (' (X,Y)-file contains '//cl3(nK), ' spectra : ')
            CALL GetNList (' Corresponding channels of z-file',
     *         cLisK, qLisK, n2)
            IF (iqSum(qLisK,n2).ne.nK) THEN
               Fehler = 'bad number of channels given'
               RETURN
               ENDIF
         ELSE
            qLisK(1) = .true.
            CALL qSet (qLisK, 2, n2, 1, .false.)
            ENDIF

C  Documentation :
         com = 'New z = '//coZ
         IF (unZ.ne.' ') CALL Compose2 (com, '('//unZ, ')')
         CALL Append (com, ' : '//doc)
         CALL OlfComAdd (jout, 'T', com, Fehler)
         CALL OlfCnuP (jout, 'z1', coZ, unZ, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Old/New z-coordinate (19aug92, sehr unschoen) :
c96/7         izold = iPar(10)
c96/7         IF (izold.gt.0) rPar(izold) = zold
c96/7         iPar(10) = iz
c96/7         IF (iz.gt.0) rPar(iz) = 0.

C  Loop spectra :
cc         CALL OlfSpe..ZG (j, 1, nC, zold, X, Y, D, Fehler)
         Kin  = 0
         DO K = 1, n2

            IF (qLisK(K)) THEN
               ! new spectrum from (X,Y)-file :
               Kin = Kin + 1
               CALL OlfGetXYD (j, Kin, n, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDIF

            Kout = Kout + 1
            CALL OlfPutSpe (jout, Kout, 1, ZZ2(K), n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! K

         CALL OlfClos (jout, Kout, Fehler)

         ENDDO ! j

      END ! OprTensor
