C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i87   :     read raw data / light scattering
C
C  ====================================================================

C  ====================================================================
C  Raman and FPI raw data
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE RRawRam (Fehler)
C  --------------------------------------------------------------------
            ! Neues Datenformat ab 4/98

      IMPLICIT NONE

      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'
      INCLUDE      'i_wrk.f'
      INTEGER       MHLin, MWL
      PARAMETER    (MHLin=50, MWL=15)
      REAL*8        Z(MZ), c0, c0x, wvOff, wvStep, wvDel,
     *              Ti, Tf, Ts, Tsam, Treg, Tsamv, Tregv, T, TkB, beta,
     *              bg, tth, pow, dwe, fac, fact, facX,
     *              wvi, wvf, wv, wv0, rzOlfGG
      INTEGER       LRF(MK), level, iUnX, nOffs, iOffs, jLRF, nLRF,
     *              ih1, jout, Kout, K, nK, nZ, n, i, ii, n1, n2, nsc,
     *              nRF, nHLin
      CHARACTER     Fehler*(*), aux*13, TLRF*160,
     *              HLab(MHLin)*20, HDat(MHLin)*60
      CHARACTER*20  CoX, UnX, CoY, UnY, h1
      CHARACTER*80  Dir, samSearch, samFound, sli, com, pol

      DATA          iUnX / 1 /, nOffs / 1 /,
     *              c0 / 2.997925e5 /    ! in nm*THz

      CALL ExeML ('\p dir-raw-ram', Dir)

C  Preset - which level of data reduction ?
      level = iAskDMu ('Raw(0) -dark(1) ->freq(2) susc(4)',
     *                 level, 0, 5)
      IF (level.ge.1) THEN
         bg = rExeMLP ('raman-dark') ! dark counts per msec
         ENDIF

C  Preset - coordinates :
      IF (level.ge.2) THEN
         CoX = 'f'
         iUnX = iAskDMu ('Frequency unit: GHz(1) THz(2)', iUnX, 0, 2)
         IF     (iUnX.eq.0) THEN
            Fehler = ' '
            RETURN
         ELSEIF (iUnX.eq.1) THEN
            UnX = 'GHz'
            c0x = c0 * 1e3
         ELSEIF (iUnX.eq.2) THEN
            UnX = 'THz'
            c0x = c0
         ELSE
            CALL Absturz ('RRR', 'iUnX')
            ENDIF
         ! wavelength offset :
         wvOff = rAskD ('Laser line offset (in pm)', wvOff)
         nOffs = iAskDMu ('Number of different offsets', nOffs, 0, MK)
         IF     (nOffs.lt.1) THEN
            Fehler = ' '
         ELSEIF (nOffs.gt.1) THEN
            wvStep = rAskD ('Offset step (in pm)', wvStep)
            IF (wvStep.le.0) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF
      ELSE
         CoX = 'lambda'
         UnX = 'nm'
         nOffs = 1
         wvOff = 0 ! Behelf: besser gar keinen Parameter _offs abspeichern
         ENDIF

      IF (level.lt.4) THEN
         CALL Compose3 (CoY, 'I(', CoX, ')')
         UnY = 'cts/sec'
      ELSE
         CALL Compose3 (CoY, 'Chi''''(', CoX, ')')
         UnY = ' '
         ENDIF

C  Prepare output file :
c      TLRF = ' '
 3    CONTINUE
         CALL GetJList ('Numbers of raw data files',
     *                  TLRF, MK, nLRF, LRF, 1, 9999)
         IF (nLRF.eq.0) RETURN

         IF (nLRF.eq.1) THEN
            CALL OlfCreate (jout, Kout, 'r'//cl6(LRF(1)),
     *                      '&olddef', Fehler)
         ELSE
            CALL OlfCreate (jout, Kout, '&nodef', '&olddef', Fehler)
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuP (jout, 'x', CoX, UnX, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'y', CoY, UnY, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuP (jout, 'z1', 'T', 'K', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'z2', 'T_set', 'K', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'z3', '2th', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'z4', 'lambda_i', 'nm', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'z5', 'lambda_offs', 'nm', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'z6', 'numor', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         nZ = 6

         CALL iOlfP (jout, '?det-bal-sym',   0, Fehler)
         CALL iOlfP (jout, '@sam-erg-gain', -1, Fehler)
         CALL iOlfP (jout, 'plot-sy#',  0, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, ' ', 'Raman/ '//TLRF, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (level.ge.1) THEN
            CALL NiceNum (bg, h1, ih1)
            CALL OlfComAdd (jout, ' ',
     *            'y - dark counts '//h1(1:ih1)//'/msec', Fehler)
             ENDIF

C  Loop / input files :
         samSearch = ' '

         DO jLRF = 1, nLRF
            nRF = LRF(jLRF)

C  - default values for header block :
            n   = -1
            nsc = 1
            Ti  = 0
            Tf  = 0
            Ts  = 0
            tth = 0
            dwe = 0
            wvi = 0
            wvf = 0

            CALL RRawRT98 (Dir, 'ram', nRF, nHLin, MHLin, HLab, HDat,
     *                     samSearch, samFound, n, Y1, Fehler)
            IF (Fehler.eq.'&sam2') THEN
               Print *, ' WARNING: different samples ? ',
     *            samFound(1:lenU(samFound)), ' after ',
     *            samSearch(1:lenU(samSearch))
               CALL Gong (1)
               Fehler = '&ff'
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
cdeb            Print *, '  read input file'

            DO i = 1, nHLin
cdeb               Print *,i,'"',HLab(i),'" -> "',HDat(i)(1:lenU(HDat(i))),'"'
               IF     (HLab(i).eq.'set temperature') THEN
                  CALL Fi1R (HDat(i), Ts)
               ELSEIF (HLab(i).eq.'start temperature') THEN
                  CALL Fi1R (HDat(i), Ti)
               ELSEIF (HLab(i).eq.'end temperature') THEN
                  CALL Fi1R (HDat(i), Tf)
               ELSEIF (HLab(i).eq.'sample temperature') THEN
                  CALL Fi1R (HDat(i), Tsam)
               ELSEIF (HLab(i).eq.'sample temp var') THEN
                  CALL Fi1R (HDat(i), Tsamv)
               ELSEIF (HLab(i).eq.'regul. temperature') THEN
                  CALL Fi1R (HDat(i), Treg)
               ELSEIF (HLab(i).eq.'regul. temp var') THEN
                  CALL Fi1R (HDat(i), Tregv)
               ELSEIF (HLab(i).eq.'scattering angle') THEN
                  CALL Fi1R (HDat(i), tth)
               ELSEIF (HLab(i).eq.'laser wavelength') THEN
                  CALL Fi1R (HDat(i), wv0)
               ELSEIF (HLab(i).eq.'laser power') THEN
                  CALL Fi1R (HDat(i), pow)
               ELSEIF (HLab(i).eq.'no. points') THEN
                  ! n has already been set in RRawRT98 ()
               ELSEIF (HLab(i).eq.'no. sweeps' .or.
     *                 HLab(i).eq.'no. scans') THEN
                  CALL Fi1I (HDat(i), nsc)
               ELSEIF (HLab(i).eq.'dwell time') THEN
                  CALL Fi1R (HDat(i), dwe)
               ELSEIF (HLab(i).eq.'start wavelength') THEN
                  CALL Fi1R (HDat(i), wvi)
               ELSEIF (HLab(i).eq.'end wavelength') THEN
                  CALL Fi1R (HDat(i), wvf)
               ELSEIF (HLab(i).eq.'slits') THEN
                  sli = HDat(i)
               ELSEIF (HLab(i).eq.'polarization') THEN
                  IF (jLRF.eq.1) THEN
                     pol = HDat(i)
                     CALL OlfComAdd (jout, ' ',
     *                               'polarization: '//pol, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                  ELSE
                     IF (HDat(i).ne.pol) THEN
                        Print *, ' WARNING: different polarisation ? ',
     * HDat(i)(1:lenU(HDat(i))), ' vs ', pol(1:lenU(pol))
                        CALL Gong (1)
                        ENDIF
                     ENDIF
                  pol = HDat(i)
                  ENDIF
               ENDDO
cdeb            Print *, '  got all parameters'

C  - elementary checks:
            IF (n.le.0) THEN
               Fehler = 'n not found or <= 0'
               RETURN
               ENDIF
            IF (nsc.lt.1) THEN
               Fehler = 'no scans recorded'
               RETURN
               ENDIF
            IF (nsc.ne.1) THEN
               Fehler = 'nsc<>1 in raman not yet foreseen'
               RETURN
               ENDIF
            IF (wvi.eq.wvf) THEN
               Fehler =
     * 'constant wavelength incompatible with chosen output mode'
               RETURN
               ENDIF

C  - compare temperatures:
            IF (Ti.eq.0) THEN
               IF (Tf.eq.0) THEN
                  T = Tsam
               ELSE
                  T = Tf
                  ENDIF
            ELSE
               IF (Tf.eq.0) THEN
                  T = Ti
               ELSE
                  T = (Ti + Tf) / 2
                  ENDIF
               ENDIF

            IF (T.eq.0 .and. Ts.ne.0) T = Ts

C  - z coordinates:
            Z(1) = T
            Z(2) = Ts
            Z(3) = tth
            Z(4) = wv0
            Z(6) = nRF

            IF (dwe.le.0) THEN
               Fehler = 'dwell time <= 0'
               RETURN
               ENDIF
            fac  = 1.d0  / ( dwe / 1.d3 ) ! dwe is in msec

C  - input data:
            DO i = 1, n
               D1(i) = dsqrt0 (Y1(i))
               ENDDO

C  - subtract background:
            IF (level.ge.1) THEN
               DO i = 1, n
                  Y1(i) = Y1(i) - bg * dwe
                  D1(i) = D1(i) + bg * dwe
                  ENDDO
               ENDIF

C  - set frequency axis, convert count rates:
            ! loop - spectra with different wavelength offset :
            DO iOffs = 1, nOffs

               ! calculate wavelength shift in nm (preset given in pm) :
               wvDel = ( wvOff + (2*iOffs - (nOffs+1))*wvStep/2 ) *
     *                 1.d-3
               Z(5)  = wvDel

               DO i = 1, n
                  wv   = wvi + (i-1) * (wvf-wvi) / (n-1)
                  IF     (CoX.eq.'f') THEN
                     X(i) = c0x * ( 1 / wv - 1 / (wv0-wvDel) )
                  ELSEIF (CoX.eq.'lambda') THEN
                     X(i) = wv
                  ELSE
                     Fehler = 'invalid CoX: '//CoX
                     RETURN
                     ENDIF
                  Y(i) = fac * Y1(i)
                  D(i) = fac * D1(i)
                  ENDDO

               Kout = Kout + 1
               CALL OlfPutSpe (jout, Kout, nZ, Z, n, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               ENDDO ! end loop offset (iOffs)

            ENDDO ! end loop input files (jLRF)

         nK = Kout
         CALL OlfClos (jout, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfOpen (jout, 1, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL TensorCheckZ (jout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  - upgrading the data ...
         IF (level.ge.4) THEN ! convert into susceptibilities
            Kout = nK
            CALL OlfCnuP (jout, 'z+', 'iSEG', ' ', Fehler)
            DO K = 1, nK ! Increase nZ. Won't work in the next loop.
               CALL OlfGetZ (jout, K, nZ, Z, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL OlfPutZ (jout, K, nZ+1, Z, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO
            CALL UnitConv (UnX, 'meV', facX, Fehler)
            DO K = 1, nK
               CALL OlfGetSpe (jout, 1, nZ, Z, n, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               T  = rzOlfGG (jout, 1, 'T', 'K', Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (T.lt.0.1) THEN ! limit was 10
                  Fehler = 'temperature too low'
                  RETURN
                  ENDIF
               TkB   = T  * .0861733 ! K -> meV

               n1 = 0
               n2 = 0
               DO i = 1, n
                  beta = X(i)*facX / TkB
                  fact = ( dexp2(beta)-1 )
                  IF     (X(i).gt.0) THEN
                     n1 = n1 + 1
                     X1(n1) =  X(i)
                     Y1(n1) =  Y(i) * fact
                     D1(n1) =  D(i) * fact
                  ELSEIF (X(i).lt.0) THEN
                     n2 = n2 + 1
                     X2(n2) = -X(i)
                     Y2(n2) = -Y(i) * fact
                     D2(n2) = -D(i) * fact
                     ENDIF
                  ENDDO
               DO i = 1, n1
                  X(i) = X1(n1+1-i)
                  Y(i) = Y1(n1+1-i)
                  D(i) = D1(n1+1-i)
                  ENDDO

               IF (n1.gt.0) THEN
                  Z(nZ) = -1
                  Kout = Kout + 1
                  CALL OlfPutSpe (jout, Kout, nZ, Z, n1,
     *                            X, Y, D, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               IF (n2.gt.0) THEN
                  Z(nZ) = +1
                  Kout = Kout + 1
                  CALL OlfPutSpe (jout, Kout, nZ, Z, n2,
     *                            X2, Y2, D2, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               CALL OlfDelSpe (jout, 1, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               Kout = Kout - 1

               ENDDO ! K -> Kout
            ENDIF         ! level 4 -> susceptibility

         CALL OlfClos (jout, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         Print *
         TLRF = '-'
         GOTO 3

      END ! RRawRam

C  --------------------------------------------------------------------
      SUBROUTINE RRawFPI (Fehler)
C  --------------------------------------------------------------------
            ! Datenformat aehnlich wie bei Raman ab 8/98

      IMPLICIT NONE

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      INTEGER       MHLin, MWL
      PARAMETER    (MHLin=50, MWL=15)
      REAL*8        Z(MZ), Z1(MZ), Z2(MZ), YSum(MWL), ZKw(MK),
     *              c0, cn0, cwi, dfc, xCutLo, xCutHi,
     *              bg, bg2, fac, fac2, facX, fact,
     *              Tset, Tsam, tth, wv0, pow, spd, smi, fsr,
     *              Tset2, Tsam2, tth2, wv02, pow2, spd2, smi2, fsr2,
     *              smiold, T, TkB, beta,
     *              ysu, yout, varWL, varWLmax,
     *              rzOlfGG
      INTEGER       LRF(MK), JWL(MWL), nJWL, iJWL, nLRF, jLRF, jwf,
     *              nRF, nRFoff, nHLin, iwc, iZw, K, nK, Kw, nKw, nZ,
     *              level, nsc, nsc2, kWL, jout, Kout, i, ii, n, n1,
     *              n2, ndummy, ih1, ih2, ih3, ipol
      CHARACTER     Fehler*(*), TLRF*160, HLab(MHLin)*20,
     *              HDat(MHLin)*60, h1*20, h2*20, h3*20, pcw*1, Unw*20
      CHARACTER*80  Dir, date, date2, samSearch, samFound, bfi, bfi2,
     *              pol, pol2, TJWL, com, fnam
      LOGICAL       qInd, qNew, qRepeat, qwf

      DATA          c0 / 2.997925d2 /,         ! c0 in nm*GHz
     *              level / 3 /, kWL / 2 /,
     *              fnam / ' ' /, dfc / 0.d0 /, xCutLo / 0.045 /,
     *              xCutHi / 0.95 /

      qRepeat = .false.
      CALL ExeML ('\p dir-raw-fpi', Dir)

      level = iAskDMu (
     * 'Raw(0) -dark(1) ->freq(2) /white(3) susc(4) cut(5)',
     * level, 0, 5)
      IF (level.ge.1) THEN
         bg = rExeMLP ('fpi-dark')
         ENDIF
      IF (level.ge.2) THEN
         dfc = rAskD ('Frequency shift / channels', dfc)
         ENDIF
      IF (level.ge.3) THEN
         kWL = iAskDMu (
     * 'White light < enter(1) neighb(2) fit-to-neighb(3)',
     * kWL, 0, 3)
         IF (kWL.le.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         IF (kWL.eq.3) THEN ! fit to white light (9jun00)
            jwf = iAskDMu (
     * ' Fit to white light from internal file no.', jwf, 1, MF)
            qwf = qOlfGdef (jwf, '?cu', 0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (.not. qwf) THEN
               Fehler = ' file is not a curve'
               RETURN
               ENDIF
            iwc = iOlfG (jwf, 'fu#', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL pcOlfFind (jwf, 'z', 'numor', Unw, pcw, iZw, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (Unw.ne.' ') THEN
               Fehler = 'Strange unit of numor: '//Unw
               RETURN
               ENDIF
            CALL OlfGet1ZofK (jwf, iZw, nKw, ZKw, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDIF
         ENDIF
      IF (level.ge.5) THEN
         xCutLo = rAskD('Lower cutoff',xCutLo)
         xCutHi = rAskD('Upper cutoff',xCutHi)
c         xCutLo = rExeMLP ('fpi-cut-lo')
c         xCutHi = rExeMLP ('fpi-cut-hi')
         ENDIF

C  Outer loop / lists of files :
      nRFoff = 0 ! iAskD ('Raw data file number offset', nRFoff)
 3    CONTINUE

C  Kind of loop / data sets orr individual files ?
         CALL GetJList (' Numbers of raw data files',
     *                  TLRF, MK, nLRF, LRF, 1, 999999)
         IF (nLRF.eq.0) RETURN

         IF (nLRF.gt.1 .and. .not.(qRepeat .and. qInd)) THEN
            CALL FrageCD (' File name (blank: save individually)',
     *                    fnam, fnam)
            qInd = fnam.eq.' '
         ELSE
            qInd = .true.
            ENDIF

C  Loop / input files :
         samSearch = ' '
         varWLmax = 0
         DO jLRF = 1, nLRF
            nRF = nRFoff + LRF(jLRF)

            qNew = qInd .or. jLRF.eq.1
            IF (qNew) THEN

C  - create file (ask for file name, don't ask for title):
               IF (qInd) THEN
                  CALL OlfCreate (jout, Kout, 'f'//cv6(nRF),
     *                            '&noask', Fehler)
               ELSE
                  CALL OlfCreate (jout, Kout, '&noask',
     *                            '&noask', Fehler)
                  CALL tOlfP (jout, 'fil', fnam, Fehler)
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN

C  - set coordinates:
               IF (level.lt.2) THEN
                  CALL OlfCnuP (jout, 'x', 'f', 'ch.', Fehler)
               ELSE
                  CALL OlfCnuP (jout, 'x', 'f', 'GHz', Fehler)
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN

               IF     (level.le.2) THEN
                  CALL OlfCnuP (jout, 'y', 'I(f)', 'cts/msec', Fehler)
               ELSEIF (level.le.3) THEN
                  CALL OlfCnuP (jout, 'y', 'I(f)', ' ', Fehler)
               ELSE
                  CALL OlfCnuP (jout, 'y', 'Chi''''(f)', ' ', Fehler)
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN

               CALL OlfCnuP (jout, 'z1',  'lambda_i', 'nm', Fehler)
               CALL OlfCnuP (jout, 'z2',  'z0[fpi1]', 'mm', Fehler)
               CALL OlfCnuP (jout, 'z3',  'fsr', 'chs.', Fehler)
               CALL OlfCnuP (jout, 'z4',  'dz/dt', 'msec/ch', Fehler)
               CALL OlfCnuP (jout, 'z5',  '2th', ' ', Fehler)
               CALL OlfCnuP (jout, 'z6',  'pol', ' ', Fehler)
               CALL OlfCnuP (jout, 'z7',  'T', 'K', Fehler)
               CALL OlfCnuP (jout, 'z8',  'T_set', 'K', Fehler)
               CALL OlfCnuP (jout, 'z9',  'pow(in)', 'mW', Fehler)
               CALL OlfCnuP (jout, 'z10', 'numor', ' ', Fehler)
               CALL OlfCnuP (jout, 'z11', '#sweeps', ' ', Fehler)
               IF (Fehler.ne.'&ff') RETURN
               nZ = 11
               IF (level.ge.3) THEN
                  CALL OlfCnuP (jout, 'z12', 'wl-var', ' ', Fehler)
                  nZ = nZ + 1
                  ENDIF

               CALL iOlfP (jout, '?det-bal-sym',   0, Fehler)
               CALL iOlfP (jout, '@sam-erg-gain', -1, Fehler)
               CALL iOlfP (jout, 'plot-sy#',  0, Fehler)
               IF (Fehler.ne.'&ff') RETURN

C  - comment line/ original numors:
               IF (qInd) THEN
                  CALL OlfComAdd (jout, ' ', 'Tandem/ '//
     *                            cv6(nRF), Fehler)
               ELSE
                  IF (nRFoff.eq.0) THEN
                     CALL OlfComAdd (jout, ' ', 'Tandem/ '//
     *                               TLRF, Fehler)
                  ELSE
                     CALL OlfComAdd (jout, ' ',
     *                 'Tandem/ '//cv6(nRFoff)//' + '//TLRF, Fehler)
                     ENDIF
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN

               ENDIF ! qNew

C  - read input file:
            CALL RRawRT98 (Dir, 'fpi', nRF, nHLin, MHLin, HLab, HDat,
     *                     samSearch, samFound, n, Y, Fehler)
            IF (Fehler.eq.'&sam2') THEN
               Fehler = '&ff'
               IF (.not. qInd) GOTO 2
               ENDIF
            IF (Fehler.ne.'&ff') RETURN

            CALL DecParFPI (MHLin, nHLin, HLab, HDat,
     *           date, pol, bfi,
     *           Tset, Tsam, tth, wv0, pow, nsc, spd, smi, fsr,
     *           Fehler)
            IF (Fehler.ne.'&ff') RETURN

C  - store some parameters:

             IF (nLRF.eq.1 .and. date.ne.' ') THEN
                CALL OlfComAdd (jout, ' ', date, Fehler)
                IF (Fehler.ne.'&ff') RETURN
                ENDIF

             IF (qInd) THEN
                CALL tOlfP (jout, 'tit', samFound,  Fehler)
             ELSE
                IF (jLRF.eq.1) THEN
                   samSearch = samFound
                   CALL tOlfP (jout, 'tit', samFound,  Fehler)
c                   CALL OlfComAdd (jout, ' ', 'sample: '//samFound, Fehler)
c                ELSE
c                   IF (eda.ne.sam) THEN
c                      Print *, ' WARNING: different samples ? ',
c     *                   eda(1:lenU(eda)), ' vs ', sample(1:lenU(sample))
c                      CALL Gong (1)
c                      ENDIF
                   ENDIF
                ENDIF
            IF (Fehler.ne.'&ff') RETURN

            ipol = 0
            IF     (pol(1:1).eq.'V') THEN
               ipol = 10
            ELSEIF (pol(1:1).eq.'H') THEN
               ipol = 20
            ELSEIF (pol(1:1).eq.'X') THEN
               ipol = 0
            ELSE
               ipol = -1
               Print *, 'unexpected polarisation: '//pol
               GOTO 3739
               ENDIF
            IF     (pol(2:2).eq.'V') THEN
               ipol = ipol + 1
            ELSEIF (pol(2:2).eq.'H') THEN
               ipol = ipol + 2
            ELSEIF (pol(2:2).eq.'X') THEN
               ipol = ipol + 0
            ELSE
               ipol = -1
               Print *, 'unexpected polarisation: '//pol
               GOTO 3739
               ENDIF
            IF (pol(3:3).ne.' ') THEN
               ipol = -1
               Print *, 'unexpected polarisation: '//pol
               ENDIF
 3739       CONTINUE

C  - elementary checks:
            IF (nsc.lt.1) THEN
               Fehler = 'no sweeps recorded'
               RETURN
               ENDIF

C  - z coordinates:
            Z( 1) = wv0
            Z( 2) = smi
            Z( 3) = fsr
            Z( 4) = spd
            Z( 5) = tth
            Z( 6) = ipol
            Z( 7) = Tsam
            Z( 8) = Tset
            Z( 9) = pow
            Z(10) = nRF
            Z(11) = nsc

            fac  = 1.d0 / nsc / spd

C  - set X axis, convert count rates:
            DO i = 1, n
               X(i) = i
               Y(i) = fac * Y(i)
               D(i) = fac * Y(i) ! internal D = D^2
               ENDDO

C  - subtract dark counts:
            IF (level.ge.1) THEN
               DO i = 1, n
                  Y(i) = Y(i) - bg
                  D(i) = D(i) + bg*nsc
                  ENDDO
               IF (qInd .or. jLRF.eq.1) THEN
                  CALL NiceNum (bg, h1, ih1)
                  CALL OlfComAdd (jout, ' ',
     *                 'y - dark counts '//h1(1:ih1)//'/msec', Fehler)
                  ENDIF
               ENDIF

C  - convert X to frequency:
            IF (level.ge.2) THEN
               IF ((jLRF.gt.1 .or. qRepeat) .and. dfc.ne.0 .and.
     *              smi.ne.smiold) THEN
                  dfc = rAskD ('Frequency shift / channels', dfc)
                  ENDIF
               smiold = smi
               cn0 = n/2 + .5 + dfc
               cwi = c0/2/smi/fsr
               DO i = 1, n ! backwards
                  Y1(i) = Y(i)
                  D1(i) = D(i)
                  ENDDO
               DO i = 1, n
                  X (i) = cwi * (i-cn0) ! sample energy gain => change sign
                  Y (i) = Y1(n+1-i)
                  D (i) = D1(n+1-i)
                  ENDDO
               IF (qInd .or. jLRF.eq.1) THEN
                  CALL NiceNum (fsr, h1, ih1)
                  CALL NiceNum (smi, h2, ih2)
                  CALL NiceNum (dfc, h3, ih3)
                  CALL OlfComAdd (jout, ' ',
     *                 'x < fsr='//h1(1:ih1)//'chs, offs='//h3(1:ih3)//
     *                 'chs, z0='//h2(1:ih2)//'mm', Fehler)
                  ENDIF
               ENDIF ! freq axis

C  - normalize to white light:
            IF (level.ge.3) THEN
               IF (kWL.eq.1) THEN
                  TJWL = ' '
                  CALL GetJList ('White light from file nos.',
     *                 TJWL, MWL, nJWL, JWL, 1, 999999)
                  IF (nJWL.le.0) THEN
                     Fehler = ' '
                     RETURN
                     ENDIF
                  com = '/ wl < '//TJWL
               ELSEIF (kWL.eq.2 .or. kWL.eq.3) THEN ! from neighbouring scan
                  nJWL = 2
                  JWL(1) = nRF-nRFoff - 1
                  JWL(2) = nRF-nRFoff + 1
                  IF     (kWL.eq.2) THEN
                     com = '/ wl < preceeding & following scan'
                  ELSEIF (kWL.eq.3) THEN
                     com = '/ wl < fit to preceeding & following scan'
                     ENDIF
               ELSE
                  Fehler = 'this choice for wl selection not yet ready'
                  RETURN
                  ENDIF

               DO i = 1, n
                  Y1(i) = 0
                  D1(i) = 0
                  ENDDO
               DO iJWL = 1, nJWL
                  IF     (kWL.le.2) THEN ! from raw-data file
                     CALL RRawRT98 (Dir, 'fpi', nRFoff+JWL(iJWL),
     *                    nHLin, MHLin, HLab, HDat,
     *                    'white light', samFound, n2, Y3, Fehler)
                     IF (Fehler.eq.'&sam2') THEN
                        Fehler = 'norm. scan not w.l. but '//samFound
                        RETURN
                        ENDIF
                     IF (Fehler.ne.'&ff') RETURN
                     CALL DecParFPI (MHLin, nHLin, HLab, HDat,
     *                  date2, pol2, bfi2,
     *                  Tset2, Tsam2, tth2, wv02, pow2, nsc2, spd2,
     *                   smi2, fsr2, Fehler)
                     IF (wv0.ne.wv02 .or. fsr.ne.fsr .or.
     *                 smi.ne.smi2 .or. bfi.ne.bfi2 .or. n.ne.n2) THEN
                         Print *, ' n  : ', n  , n2
                         Print *, ' wv0: ', wv0, wv02
                         Print *, ' fsr: ', fsr, fsr2
                         Print *, ' smi: ', smi, smi2
                         Print *, ' bfi: ', bfi, bfi2
                         Fehler =
     * 'Essential set-up not in accord for sam vs. w.l.'
                         RETURN
                         ENDIF
                     ! bring in backward order :
                     DO i = 1, n
                        Y2(i) = Y3(n+1-i)
                        ENDDO
                     fac2  = 1.d0 / nsc2 / spd2
                     bg2   = bg
                  ELSEIF (kWL.eq.3) THEN ! from fit curve
                     DO Kw = 1, nKw
                        IF (dabs(JWL(iJWL)+nRFoff-ZKw(Kw)).lt.1d-3)
     *                   GOTO 303
                        ENDDO
                     Fehler = ' Could not find numor '//
     *                        cl6(JWL(iJWL)+nRFoff)
                     RETURN
 303                 CONTINUE
                     CALL OlfGetY (jwf, Kw, ndummy, Y4, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     CALL CuFuVal (iwc, Y4, X, Y2, n, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     fac2 = 1
                     bg2  = 0
                     nsc2 = 1 ! ????
                     ENDIF
                  ! add new Y2 to overall Y1 :
                  ysu = 0
                  DO i = 1, n
                     IF (dabs(i-cn0).gt.17 .and. dabs(i-cn0).lt.n-7)
     *                    ysu = ysu + dsqrt(dabs(i-cn0)) * Y2(i)
     *                       ! ziemlich willkuerliches Kriterium
                     Y1(i) = Y1(i) + (fac2   *Y2(i)-bg)
                     D1(i) = D1(i) + (fac2**2*Y2(i)+bg*nsc2)
                     ENDDO
                  YSum(iJWL) = ysu ! ISO 9001
                  ENDDO ! iJWL

               IF (nJWL.ge.2) THEN
                  varWL = dquot0 (YSum(nJWL)-YSum(1),
     *                            YSum(nJWL)+YSum(1))
                  Z(12) = varWL
                  varWLmax = dmax1 (dabs(varWL), varWLmax)
               ELSE
                  Z(12) = 0
                  ENDIF

               DO i = 1, n
                  yout = dquot0 (Y(i)*nJWL, Y1(i))
                  D(i) = yout**2 * ( dquot0(D (i),Y (i)**2) +
     *                               dquot0(D1(i),Y1(i)**2) )
                  Y(i) = yout
                  ENDDO

               IF (qInd .or. jLRF.eq.1) THEN
                  CALL OlfComAdd (jout, ' ', com, Fehler) ! nicht immer richtig
                  ENDIF
               ENDIF ! wl normalization


C  - restore D from D^2
            DO i = 1, n
               D(i) = dsqrt0 (D(i))
               ENDDO


C  - save:
            Kout = Kout + 1
            CALL OlfPutSpe (jout, Kout, nZ, Z, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

 2          CONTINUE

C  - collect unused z coordinates:
            IF (qInd .or. jLRF.eq.nLRF) THEN
               CALL OlfClos (jout, Kout, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               nK = Kout
               CALL OlfOpen (jout, 1, Kout, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               CALL TensorCheckZ (jout, Fehler)
               IF (Fehler.ne.'&ff') RETURN

C  - upgrading the data ...
               IF (level.ge.4) THEN ! convert into susceptibilities
                  Kout = 0
                  CALL OlfCnuP (jout, 'z+', 'iSEG', ' ', Fehler)
                  DO K = 1, nK ! Increase nZ. Won't work in the next loop.
                     CALL OlfGetZ (jout, K, nZ, Z, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     CALL OlfPutZ (jout, K, nZ+1, Z, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO
                  CALL UnitConv ('GHz', 'meV', facX, Fehler)
                  DO K = 1, nK
                     CALL OlfGetSpe (jout, 1, nZ, Z, n,
     *                               X, Y, D, Fehler)
                     IF (Fehler.ne.'&ff') RETURN

                     T  = rzOlfGG (jout, 1, 'T', 'K', Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     IF (T.lt.0.1) THEN ! limit was 10
                        Fehler = 'temperature too low'
                        RETURN
                        ENDIF
                     TkB   = T  * .0861733 ! K -> meV

                     n1 = 0
                     n2 = 0
                     DO i = 1, n
                        beta = X(i)*facX / TkB
                        fact = ( dexp2(beta)-1 )
                        IF     (X(i).gt.0) THEN
                           n1 = n1 + 1
                           X1(n1) =  X(i)
                           Y1(n1) =  Y(i) * fact
                           D1(n1) =  D(i) * fact
                        ELSEIF (X(i).lt.0) THEN
                           n2 = n2 + 1
                           X2(n2) = -X(i)
                           Y2(n2) = -Y(i) * fact
                           D2(n2) = -D(i) * fact
                           ENDIF
                        ENDDO
                     DO i = 1, n2
                        X(i) = X2(n2+1-i)
                        Y(i) = Y2(n2+1-i)
                        D(i) = D2(n2+1-i)
                        ENDDO

                     Z(nZ) = -1
                     CALL OlfPutSpe (jout,nK+K,   nZ,Z, n1,
     *                               X1,Y1,D1, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     Z(nZ) = +1
                     CALL OlfPutSpe (jout,nK+K+1, nZ,Z, n2,
     *                               X, Y ,D , Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     CALL OlfDelSpe (jout, 1, Fehler)
                     IF (Fehler.ne.'&ff') RETURN

                     ENDDO ! K -> Kout
                     Kout = 2*nK
                  ENDIF         ! level 4 -> susceptibility

C  - cut elastic peak(s):
               IF (level.ge.5) THEN
                  PRINT *, ' Cutting spectrum at ',xCutLo*c0/2./smi
     *                    ,' and ',xCutHi*c0/2./smi

                  DO K = 1, Kout
                     CALL OlfGetSpe (jout, 1, nZ, Z, n,
     *                               X, Y, D, Fehler)
                     IF (Fehler.ne.'&ff') RETURN

                     ii = 0
                     DO i = 1,n
                        IF ( X(i).gt.(xCutLo*c0/2./smi)
     *                       .and. X(i).lt.(xCutHi*c0/2./smi) ) THEN
                           ii = ii + 1
                           X1(ii) = X(i)
                           Y1(ii) = Y(i)
                           D1(ii) = D(i)
                        ENDIF
                     ENDDO

                     n = ii
                     CALL OlfPutSpe (jout,Kout+1, nZ,Z, n,
     *                               X1,Y1,D1, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     CALL OlfDelSpe (jout, 1, Fehler)
                     IF (Fehler.ne.'&ff') RETURN

                  ENDDO
               ENDIF            ! cut peak(s)

C  - at the end saving the result :
               CALL OlfClos (jout, Kout, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               ENDIF

            ENDDO ! nLRF

         Print *, ' maximum white light variation: ', varWLmax
         Print *
         qRepeat = .true.
         TLRF = '-'
         GOTO 3

      END ! RRawFPI

C  --------------------------------------------------------------------
      SUBROUTINE RRawRT98 (Dir, Ext, nRF, nHLin, MHLin, HLab, HDat,
     *                     samSearch, samFound, n, Y, Fehler)
C  --------------------------------------------------------------------
            ! as separate subroutine JWu 5jun99
         ! read a .fpi or .ram file without analysing the header data

      IMPLICIT NONE
      INCLUDE 'l_def.f'

      INTEGER       nRF, i, ii, n, nHLin, MHLin
      REAL*8        Y(*)
      CHARACTER*(*) Fehler, Dir, Ext, samSearch, samFound,
     *              HLab(MHLin), HDat(MHLin)
      CHARACTER     lab*20, dat*60, line*80, FileRaw*80

      IF (Fehler.ne.'&ff') THEN
         Print *, 'Error on entry in RRawRT98'
         Print *, '"'//Fehler(1:lenU(Fehler))//'"'
         RETURN
         ENDIF

C  Open File:
      CALL Compose2 (FileRaw, Dir, cv6(nRF))
      CALL OpenFile (11, FileRaw, Ext, 'a', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Read format line:
      Read (11, '(a)', err=91) line
      IF (line.ne.'Tandem raw data format 98' .and.
     *    line.ne.'Raman raw data format 98') THEN
         Fehler = 'bad data format: '//line(1:30)
         Close (11)
         RETURN
         ENDIF

C  Read header block:
      n = -1

      DO nHLin = 1, MHLin

         Read (11, '(a20,1x,a)', err=92) lab, dat
         IF (lab(1:3).eq.'---') GOTO 19 ! end of header block
         DO i=20,1,-1
            IF (lab(i:i).ne.'.') GOTO 12 ! deleted all trailing '.'
            lab(i:i) = ' '
            ENDDO
 12      CONTINUE ! deleted all trailing '.'

         IF     (lab.eq.'file') THEN
            CALL Fi1I (dat, ii)
            IF     (ii.le.0) THEN
               Fehler = 'no file number in input file name "'
     *                     //dat(1:lenU(dat))//'"'
               Close (11)
               RETURN
            ELSEIF (ii.ne.nRF) THEN
               Fehler = 'input file name "'
     * //dat(1:lenU(dat))//'" has other number than '//cl6(nRF)
               Close (11)
               RETURN
               ENDIF
         ELSEIF (lab.eq.'no. points') THEN
            CALL Fi1I (dat, n)
         ELSEIF (lab.eq.'sample') THEN
            samFound = dat
            IF (samSearch.ne.' ' .and. samFound.ne.samSearch) THEN
               Close (11)
               Fehler = '&sam2'
               RETURN
               ENDIF
            ENDIF

         HLab(nHLin) = lab
         HDat(nHLin) = dat
         ENDDO

      Fehler = 'too many header lines'
      RETURN

 19   CONTINUE

      IF (n.le.0) THEN
         Fehler = 'n not found or <= 0'
         RETURN
         ENDIF

C  Read data block:
      DO i = 1, (n-1)/8+1
         Read (11, *, err=93) (Y(ii), ii=i*8-7, min0(n,i*8))
         ENDDO

C  Everything read:
c      Print *, ' read ', FileRaw(1:lenU(FileRaw)), ' ...'
      Close (11)

      RETURN ! regular exit

 91   CONTINUE
      Fehler = 'Read error in first line'
      RETURN
 92   CONTINUE
      Fehler = 'Read error in header block'
      RETURN
 93   CONTINUE
      Fehler = 'Read error in data block'
      RETURN

      END ! RReadRT98

C  --------------------------------------------------------------------
      SUBROUTINE DecParFPI (MHLin, nHLin, HLab, HDat,
     *           date, pol, bfi,
     *           Tset, Tsam, tth, wv0, pow, nsc, spd, smi, fsr,
     *           Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'l_def.f'

      INTEGER       MHLin, nHLin, nsc, i
      REAL*8        Tset, Tsam, tth, wv0, pow, spd, smi, fsr, Ti, Tf
      CHARACTER*(*) date, pol, bfi, Fehler, HLab(MHLin), HDat(MHLin)

C  Default values:
      date   = ' '
      bfi    = ' '
      pol    = ' '
      wv0    = 514.5
      nsc    = -1
      Ti     = 0
      Tf     = 0
      Tset   = 0
      Tsam   = 0
      tth    = 0
      fsr    = 0

C  Decode:
      DO i = 1, nHLin
         IF (HLab(i).eq.'set temperature') THEN
            CALL Fi1R (HDat(i), Tset)
         ELSEIF (HLab(i).eq.'start temperature') THEN ! old style ?
            CALL Fi1R (HDat(i), Ti)
         ELSEIF (HLab(i).eq.'end temperature') THEN
            CALL Fi1R (HDat(i), Tf)
         ELSEIF (HLab(i).eq.'sample temperature') THEN
            CALL Fi1R (HDat(i), Tsam)
c         ELSEIF (HLab(i).eq.'sample temp var') THEN
c            CALL Fi1R (HDat(i), DTsam)
c         ELSEIF (HLab(i).eq.'regul. temperature') THEN
c            CALL Fi1R (HDat(i), Treg)
c         ELSEIF (HLab(i).eq.'regul. temp var') THEN
c            CALL Fi1R (HDat(i), DTreg)
         ELSEIF (HLab(i).eq.'scattering angle') THEN
            CALL Fi1R (HDat(i), tth)
         ELSEIF (HLab(i).eq.'laser wavelength') THEN
            CALL Fi1R (HDat(i), wv0)
         ELSEIF (HLab(i).eq.'laser power') THEN
            CALL Fi1R (HDat(i), pow)
         ELSEIF (HLab(i).eq.'no. sweeps' .or. HLab(i).eq.'no. scans')
     *      THEN
            CALL Fi1I (HDat(i), nsc)
         ELSEIF (HLab(i).eq.'speed [msec/ch]') THEN
            CALL Fi1R (HDat(i), spd)
         ELSEIF (HLab(i).eq.'mirror spacing') THEN
            CALL Fi1R (HDat(i), smi)
         ELSEIF (HLab(i).eq.'fsr [chs]') THEN
            CALL Fi1R (HDat(i), fsr)
         ELSEIF (HLab(i).eq.'band filter') THEN
            bfi = HDat(i)
         ELSEIF (HLab(i).eq.'polarization') THEN
            pol = HDat(i)
         ELSEIF (HLab(i).eq.'comment') THEN
            Print *, ' comment / '//HDat(i)(1:lenU(HDat(i)))
         ELSEIF (HLab(i).eq.'user') THEN
            ! ignore
         ELSE
            ! ignore without warning
            ENDIF

         ENDDO

C  - compare temperatures / old style :
      IF (Ti.eq.0) THEN
         IF (Tf.ne.0) THEN
            Tsam = Tf
            ENDIF
      ELSE
         IF (Tf.eq.0) THEN
            Tsam = Ti
         ELSE
            Tsam = (Ti + Tf) / 2
            ENDIF
         ENDIF

      END ! DecParFPI

C  --------------------------------------------------------------------
      SUBROUTINE NormFPI (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 17nov95, using brd.f (New York, 9oct-4nov92)
         ! Normalize the frequency axis of FPI data

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INTEGER       nJList, JList(*), lj, j, jout, K, nK, Kout, n, i,
     *              ii, ih1, ih2, ih3, npik, iCutLr, iCutCl, iCutCr,
     *              iCutRl
      LOGICAL       qOv, qOK, qCut
      REAL*8        c0m, dMirr, fsr, cn0, cwi
      CHARACTER     Fehler*(*), aux*13, cOK*1, h1*20, h2*20, h3*20,
     *              aus*80

      c0m = 3.d2 ! speed of light in mm*GHz

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      qCut = qAskD ('Cut ghosts and central peak', intq(qCut))
      IF (qCut) THEN
         iCutLr = iExeMLP ('fpi-cut-lr')
         iCutCl = iExeMLP ('fpi-cut-cl')
         iCutCr = iExeMLP ('fpi-cut-cr')
         iCutRl = iExeMLP ('fpi-cut-rl')
         ENDIF

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuP (jout, 'x', 'f', 'GHz', Fehler)
c         CALL OlfCnuP (jout, 'y', 'I(f)', 'cps/GHz', Fehler)
         CALL OlfComAdd (jout, 't', 'set x-scale for FPI', Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK

            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            dMirr = rzxOlfGG (j, K, 0, 'z0[fpi1]', 'mm', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (dMirr.le.0) THEN
               Fehler = 'nonpositive mirror spacing'
               RETURN
               ENDIF
            fsr =  rzxOlfGG (j, K, 0, 'fsr', 'chs.', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (fsr.le.0) THEN
               Fehler = 'fsr <= 0'
               RETURN
               ENDIF

C  Locate reference peaks :
            qOK = .true.
            cn0 = n/2 + .5
c            npik = idnint (fsr*.02) + 1
c           CALL LocPeak (Y, n, n/2, 5*npik/2, npik, cn0, 'center', qOK, Fehler)
c            IF (Fehler.ne.'&ff') RETURN
c            Print '(1x,a,1f8.2)', ' ...  peak at', cn0
c            IF (.not.qOK .or. cn0.lt.0.45*n .or. cn0.gt.0.55*n) THEN
c               cn0 = rAsk (' enter Rayleigh position')
c               ENDIF

            cwi = c0m/2/dMirr/fsr

C  Set energy scale, cut reference signal :
            ii = 0
            DO i = n,1,-1
               IF (.not.qCut .or.
     *             (i.gt.iCutLr .and. i.lt.iCutCl) .or.
     *             (i.gt.iCutCr .and. i.lt.iCutRl)      ) THEN
                  ii = ii + 1
                  X (ii) = cwi * (cn0-i) ! sample energy gain => change sign
                  Y1(ii) = Y(i)
                  D1(ii) = D(i)
                  ENDIF
               ENDDO
            n = ii

            CALL OlfCopZ (j, jout, K, K, Fehler)
             IF (Fehler.ne.'&ff') RETURN
            CALL OlfPutXYD (jout, K, n, X, Y1, D1, Fehler)
             IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (jout, nK, Fehler)
         ENDDO

      END ! NormFPI

C  --------------------------------------------------------------------
      SUBROUTINE LocPeak (Y, n, ncentre, nsearch, navge, c0, label,
     *                    qOK, Fehler)
C  --------------------------------------------------------------------
            ! JWu 26oct92, almost rewritten 30dec92,
            ! guess in terms of width and fsr 17nov95
            ! for automatic scans very simple new version 22sep98

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      REAL*8        Y(*)
      CHARACTER*(*) label, Fehler

C  Checks :
      IF (Fehler.ne.'&ff') RETURN
      IF (qioutside (ncentre, 1, n)) THEN
         Fehler = 'LocPeak/ nonsense ncentre'
         RETURN
         ENDIF
      IF (qioutside (nsearch, 1, n/2)) THEN
         Fehler = 'LocPeak/ nonsense nsearch'
         RETURN
         ENDIF
      IF (qioutside (navge, 1, nsearch)) THEN
         Fehler = 'LocPeak/ nonsense navge'
         RETURN
         ENDIF

C  Center of gravity :
      sxy = 0.
      sy  = 0.
      DO i = ncentre-nsearch, ncentre+nsearch
         sxy = sxy + (i-ncentre)*Y(i)
         sy  = sy  +   Y(i)
         ENDDO
      n0 = ncentre + idnint (sxy / sy)

C  Center of gravity :
      sxy = 0.
      sy  = 0.
      DO i = n0-navge, n0+navge
         sxy = sxy + (i-n0)*Y(i)
         sy  = sy  +   Y(i)
         ENDDO
      c0 = n0 + sxy / sy

      Print '(a,3i4,a,i3,a,f9.4)',
     * ' located peak with ', ncentre, nsearch, navge, ' at ',
     * n0, ' -> ', c0

      END ! LocPeak

C  --------------------------------------------------------------------
      SUBROUTINE LocDrop (Y, n, iL, iR, idir, bgmul, iDrop,
     *                    rDrop, Fehler)
C  --------------------------------------------------------------------
            ! JWu 13nov92, completely new 4jun96
         ! Locate a cut in Y(ni..nf) such as produced
         ! by suddenly dropping a filter into the beam
         ! or closing a chutter.
         ! Return channel number of the cut as idrop.

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)
      CHARACTER*(*)      Fehler

      REAL*8        Y(*)

C  Checks :
      IF (Fehler.ne.'&ff') RETURN
      IF     (iL.lt.3 .or. iR.gt.n-2 .or. iL.gt.iR-5) THEN
         Fehler = 'LocDrop/ bad limits on input'
         RETURN
      ELSEIF (idir.ne.+1 .and. idir.ne.-1) THEN
         Fehler = 'LocDrop/ no direction given'
         ENDIF

C  Estimate background :
       bg = 1.d20
       DO i = iL, iR
          IF (Y(i).lt.0) THEN
             Fehler = 'NormFPI/ unexpected count < 0'
             RETURN
             ENDIF
          IF (Y(i).gt.0 .and. Y(i).lt.bg) bg = Y(i)
          ENDDO
       bg = bg * bgmul

C  Search :
      rDrop = 0
      DO i = iL, iR
         ! average over neighbouring channels
            ! add background to handle countrate zero
         yL = dlog(Y(i-2)+bg) + dlog(Y(i-1)+bg)
         yR = dlog(Y(i+2)+bg) + dlog(Y(i+1)+bg)
         IF     (idir.eq.+1) THEN
            rD = yR / yL
         ELSE
            rD = yL / yR
            ENDIF
         IF (rD.gt.rDrop) THEN
            iDrop = i
            rDrop = rD
            ENDIF
         ENDDO
      END ! LocDrop


C  ====================================================================
C  History tables
C  ====================================================================

      SUBROUTINE RRawHistory (Fehler)
C     ----------------------------
         ! JWu 11apr91, included Any2IED 14sep94
      ! Make a completely new data file and store it in the run-time memory

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER     (MRL=20,MJ=10)
      REAL*8         RL(MRL), Z(MZ)
      INTEGER        J(MJ)
      CHARACTER*(*)  Fehler
      CHARACTER*40   CoYJ(MJ), UnYJ(MJ)
      CHARACTER*80   aus, InFile, DefFile, OutFile, Title
      CHARACTER*240  line

      DATA          Title / 'History' /, DefFile / 'shist' /

      idimneed = 1
      IF (MJ.gt.MK .or. idimneed.gt.MWrk3dim) THEN
         Fehler = 'Cannot use Wrk3dim'
         RETURN
         ENDIF

      Print *, ' read CIC history file '
      ih1 = iAsk (' starting with no.')
      ihu = iAsk ('     ... up to no.')
      IF (ih1.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF
      CALL FrageCD (' output files', DefFile, DefFile)

      nJ = 9
      CoYJ(1) = 'ps'
      UnYJ(1) = 'cts'
      CoYJ(2) = 'lhs'
      UnYJ(2) = 'cts'
      CoYJ(3) = 'rhs'
      UnYJ(3) = 'cts'
      CoYJ(4) = 'zz'
      UnYJ(4) = 'V'
      CoYJ(5) = 'dz'
      UnYJ(5) = 'V'
      CoYJ(6) = 'x1'
      UnYJ(6) = 'V'
      CoYJ(7) = 'y1'
      UnYJ(7) = 'V'
      CoYJ(8) = 'x2'
      UnYJ(8) = 'V'
      CoYJ(9) = 'y2'
      UnYJ(9) = 'V'

C  Make header :

      DO ij = 1, nJ

         CALL Compose2 (OutFile, DefFile, '-'//CoYJ(ij))
         CALL OlfCreate (J(ij), Kout, OutFile, Title, Fehler)
         IF (Fehler.ne.'&ff') RETURN
c         CALL tOlfG (J(ij), 'fil', OutFile, Fehler)  ! f"ur's n"achste Mal
         CALL tOlfG (J(ij), 'tit', Title, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuP (J(ij), 'x', 't', 'sec', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuP (J(ij), 'y', CoYJ(ij), UnYJ(ij), Fehler)
         IF (Fehler.ne.'&ff') RETURN

         nz = 1
         CALL OlfCnuP (J(ij), 'z1', 'no-of-hist', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! ij

C  Make spectra :
      x01 = 0
      DO K = 1, MK

C  Open input file(s) :
         InFile = '/pcnfs/fpi/dat-fpi/shist-'//cv2(ih1)

         CALL OpenFile (11, InFile, 'tab', 'a', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Print *, ' reading from file '//InFile(1:lenU(InFile))

C  Multicol input :
         Z(1) = K

         n = 0

         Read (11, '(a)', err=91) line
         Read (11, '(a)', err=91) line

 45      CONTINUE ! endless loop

            Read (11, '(a)', err=459) line

            CALL FindR (line, MRL, nRL, RL, .false.)

            IF (nRL.lt.1) THEN
               Fehler = ' No x value !'
               RETURN
               ENDIF
            X(n+1) = x01 + RL(1)
            IF (nRL.lt.nJ+1) THEN
               Fehler = ' trouble with y'
               RETURN
               ENDIF
            DO ij = 1, nJ
               Wrk3dim (n+1,ij,1) = RL(ij+1)
               ENDDO
            D(n+1) = 0
            n = n + 1
            IF (n.eq.MC) THEN
               Print *, ' Maximum number of channels reached'
               CALL Gong (6)
               GOTO 459
               ENDIF
         GOTO 45
 459     CONTINUE ! end of data
         Close(11)

         IF (n.le.0) THEN
            Fehler = ' empty x-range entered'
            RETURN
            ENDIF

C  Done - save spectrum :
         DO ij = 1, nJ
            CALL OlfPutSpe (J(ij), K, nZ, Z, n, X, Wrk3dim(1,ij,1),
     *                      D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         IF (ih1.eq.ihu) GOTO 5 ! regular exit from loop
         ih1 = ih1 + 1
         IF (ih1.ge.100) ih1 = 1
         x01 = x01 + n

         ENDDO ! spectra

      CALL Gong (3)
      Print *, ' maximum number of spectra created'
 5    CONTINUE

      DO ij = 1, nJ
         CALL OlfClos (J(ij), K, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      RETURN ! regular exit

 91   CONTINUE
      Fehler = ' no 1st line in file'
      RETURN

      END ! RRawHistory

C  ====================================================================
C  OKE
C  ====================================================================

C  --------------------------------------------------------------------
      SUBROUTINE RRawOke (Fehler)
C  --------------------------------------------------------------------
            ! JWu 6may95
         ! Read full output from Oke (may 02 format)

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      PARAMETER    (MLin=2500)
      CHARACTER*40  RFNam, RFSub, RFDir
      REAL*8        X(MC), R(3), Y1(MC), D1(MC), ZZ(1)
      CHARACTER     Fehler*(*)
      CHARACTER*80  Lin(MLin)

      DATA          qFirst / .true. /

      CALL OlfCreate (j1, K1, '&nodef', '&olddef', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCnuP (j1, 'x', 't', 'psec', Fehler)
      CALL OlfCnuP (j1, 'y', 'I', 'arb. units', Fehler)
      CALL OlfCnuP (j1, 'z1', 'Num', ' ', Fehler)
c      CALL OlfCnuP (j1, 'z2', 'Sign', ' ', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      temp = rAskDMu ('Temperature [K]', temp, 1.d0, 5.d2)
      CALL FrageCD ('Raw data directory', RFDir, RFDir)
c      CALL FrageC ('And subdirectory', RFSub)
      CALL rOlfP (j1, 'T', 'K', temp, Fehler)

C  Loop over input files :
      DO K = 1, MK

         iRF = iAskDu ('Raw data number [quit] ?', 0)
         IF (iRF.eq.0) GOTO 9
         CALL Compose2 (RFNam, RFDir, '/B'//cl6(iRF))
         CALL OpenFile (37, RFNam, 'Y01', 'a', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL ReadFile (37, Lin, MLin, nLin, '&eof', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         z = iRF
         n = 0

         DO il = 1, nLin
            CALL FindR (Lin(il), 3, nR, R, .false.)
            IF (nR.ne.2) THEN
               Fehler = 'RR Oke not 2 numbers in line '//cl4(il)
               RETURN
               ENDIF
               n = n+1 ! new data point
               X(n) = R(1)
               Y1(n) = R(2)
               D1(n) = 0

            ENDDO ! loop over log file lines

C  Save Results :
            K1 = K1 + 1
            ZZ(1) = z
         CALL OlfPutSpe (j1, K1, 1, ZZ, n, X, Y1, D1, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO
C  End loop over input files.

 9    CONTINUE
      CALL OlfClos (j1, K1, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! RRawOke
