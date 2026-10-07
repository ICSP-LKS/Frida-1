C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i73   :     density of states pp.
C
C  ====================================================================

C     Contents :
C            SEGconv, DOSconv, DOS
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE SEGconv (nJlist, JList, Fehler)
C  --------------------------------------------------------------------
            ! JWu 7jul93.
         ! Switch between Vercors group and rest of the world.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      CHARACTER      aus*80
      CHARACTER*40   coX, unX

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, .true., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuG (j, 'x', coX, unX, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (coX.ne.'w' .and. coX.ne.'f') THEN
            CALL Compose2 (aus, ' x-coordinate is '//
     *                     coX, ' - Continue ?')
            IF (.not.qAsk(aus)) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF

         iSEG = iOlfG (j, '@sam-erg-gain', Fehler)
         IF (Fehler.ne.'&ff' .or. iabs(iSEG).ne.1) THEN
            Fehler = '&ff'
            IF (lj.gt.1) THEN
               qSEG = qAskD(' Sign convention yet undefined: '//
     * ' does E>0 mean sample energy gain', intq(qSEG))
            ELSE
               qSEG = qAsk (' Sign convention yet undefined: '//
     * ' does E>0 mean sample energy gain ?')
               ENDIF
            iSEG = -1 + 2*intq(qSEG)
            qInv = .false.
         ELSE
            IF     (iSEG.eq.1) THEN
               Print *, ' E>0 means sample energy gain'
            ELSEIF (iSEG.eq.-1) THEN
               Print *, ' E>0 means sample energy loss'
            ELSE
               Fehler =
     * 'PROGRAM ERROR/ sign undefined/ should not arrive here'
               RETURN
               ENDIF
            IF (lj.gt.1) THEN
               qInv = qAskD(' Invert axis', intq(qInv))
            ELSE
               qInv = qAsk (' Invert axis ?')
               ENDIF
            IF (qInv) iSEGin = -iSEGin
            ENDIF

         CALL iOlfP (j, '@sam-erg-gain', iSEG, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! info only :
         qSym = qOlfG (j, '?det-bal-sym', Fehler)
         IF (Fehler.ne.'&ff') THEN
            Fehler = '&ff'
            IF (lj.le.1) THEN
               qSymPut =
     * qAsk ('Declare the spectra as det-bal symmetrized ?')
            ELSE
               qSymPut =
     * qAskD('Declare the spectra as det-bal symmetrized ?',
     *       intq(qSymPut))
               ENDIF
            CALL iOlfP (j, '?det-bal-sym', intq(qSymPut), Fehler)
         ELSEIF (qSym) THEN
            Print *, ' data are det-bal symmetrized'
         ELSE
            Print *, ' data aren''t det-bal symmetrized'
            ENDIF

         IF (qInv) THEN
            DO K = 1, nK
               CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               DO i = 1, n
                  X(i) = X(-i)
                  ENDDO
               CALL OlfPutXYD (jout, K, n, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO
            ENDIF

         CALL OlfClos (jout, nK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! lj

      END ! SEGconv

C  --------------------------------------------------------------------
      SUBROUTINE DOSconv (nJlist, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 20aug91, 28apr92. Formerly OprSpecial
            ! Rewritten with integer exponents 7/8jul94.
         ! Conversions between S,S~,g,G,Q,...
         ! Intermediate function is S(q,w)
            ! divide input by factor^ifacI
            ! then multiply by factor^ifacO

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      DIMENSION     DetE1(MK), DetE2(MK), iDetTest(512)
      CHARACTER*(*) Fehler
      INTEGER       JList(*)
      REAL*8        yQ(MC), yA(MC), yB(MC)
      CHARACTER*80  aus, doc, ord, h1, h2
      CHARACTER*40  CoXin, UnXin, CoYin, UnYin, CoYout, UnYout, Inst
      DATA          AtNu / 1.d0 /


C  Menu :
      Print *, ' Convert input y(x,z) to'
      Print *, '    (1)  Tof->E             (2)  S~(q,w)'
      Print *, '    (3)  g_1(w)             (4)  g_(1)(w)'
      Print *, '    (5)  G_1(w)             (6)  G_(1)(w)'
      Print *, '    (7)  c(w)               (8)  c(w)/T^3'
      Print *, '    (9)  u^2(w)            (10)  X''''(q,w)'
      Print *, '   (11)  q(2th,w)          (12)  exp(-2W)'
      Print *, '   (13)  1/n(w;T)          (14)  1/n~(w;T)'
      iOS = iAskDMu (' Enter option', iOS, 0, 16)
      IF (iOS.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop files :
      DO lj = 1, nJList
      j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuG (j, 'x', CoXin, UnXin, Fehler)
         CALL OlfCnuG (j, 'y', CoYin, UnYin, Fehler)
         Print *
         Print *,  ' convert file '//cr2(j)//' :  '//CoYin

         iSEG = iOlfG (j, '@sam-erg-gain', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (iabs(iSEG).ne.1) THEN
            Print *, ' ivalid iSEG = ', iSEG
            Fehler = ' sign of energy unknown/ run _pm first'
            RETURN
            ENDIF

         IF (CoXin.ne.'w' .and. CoXin.ne.'w/2pi' .and. CoXin.ne.'f'
     *       .and. CoXin.ne.'tof') THEN
            Print *, ' input x coordinate = '//CoXin
            IF (.not.qAsk(' Continue ?')) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF
         CALL UnitConv (UnXin, 'meV', facX, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Analyse output coordinate :

         ! Default settings :
         iSymO = 0
         mOrdO = 0
         idwfO = 0
         iBosO = 0
         iBetO = 0
         iQQQO = 0 ! just Q
         iCCCO = 0 ! all we need for c(w)
         iTTTO = 0 ! divide by T^3
         iUUUO = 0
         qNewY = .false.
         UnYout= UnYin
         qCtoE = .false.

         IF (iOS.eq.1) THEN
             CoYout = 'S(2th,w)'
             UnYout = 'meV-1'
             qCtoE = .true.
         ELSEIF (iOS.eq.2) THEN
            iSymO = 1
            CoYout  = 'S~(q,w)'
         ELSEIF (iOS.eq.3) THEN
            CoYout  = 'g_1(w)'
            mOrdO = 1
            idwfO = 1
            iBosO = 1
         ELSEIF (iOS.eq.4) THEN
            CoYout  = 'g_(1)(w)'
            mOrdO = 1
            iBosO = 1
         ELSEIF (iOS.eq.5) THEN
            CoYout  = 'G_1(w)'
            mOrdO = 1
            idwfO = 1
         ELSEIF (iOS.eq.6) THEN
            CoYout  = 'G_(1)(w)'
            mOrdO = 1
         ELSEIF (iOS.eq.7) THEN
            CoYout  = 'c(w)'
            mOrdO = 1
            iCCCO = 1
         ELSEIF (iOS.eq.8) THEN
            CoYout  = 'c(w)/T^3'
            mOrdO = 1
            iCCCO = 1
            iTTTO = 1
         ELSEIF (iOS.eq.9) THEN
            CoYout  = 'u^2(w)'
            iUUUO = 1
         ELSEIF (iOS.eq.10) THEN         ! chi'' as suggested by a PRL referee
            CoYout = 'Chi'''''
            UnYout = ' '
            iBosO  = 1
            iBetO  = -1
         ELSEIF (iOS.eq.11) THEN         ! supersede y by Q(x) :
            CoYout  = 'q(2th,w)'
            UnYout  = 'A-1'
            qNewY = .true.
            iQQQO = 1
         ELSEIF (iOS.eq.12) THEN         ! supersede y by DWF :
            CoYout  = 'exp(-2W)'
            UnYout  = ' '
            qNewY = .true.
            idwfO = 1
         ELSEIF (iOS.eq.13) THEN         ! supersede y
            CoYout  = '1/n(w;T)'
            UnYout  = ' '
            qNewY = .true.
            iSymO = 0
            iBosO = 1
            iBetO =-1
         ELSEIF (iOS.eq.14) THEN         ! supersede y
            CoYout  = '1/n~(w;T)'
            UnYout  = ' '
            qNewY = .true.
            iSymO = 1
            iBosO = 1
            iBetO =-1
         ELSE
            Fehler = ' Option not yet implemented'
            RETURN
            ENDIF

C  More questions about the output :
         IF (iCCCO.ne.0) THEN !  for c(T)
            AtNu = rAskD (' Number of atoms per molecule', AtNu)
            ENDIF

C  Analyse input y coordinate -> set exponents, mOrdI, ...

         ! Default setting of exponents :
         mOrdI = 0 ! no factor alpha/3
         iSymI = 0 ! not det-bal symmetrized
         idwfI = 0 ! no Debye Waller factor
         iBosI = 0 ! no Bose factor (n(beta)/beta)
         iBetI = 0 ! no beta factor

C  Auxiliary calculations :
         facI = 1.
         DO m = 1, mOrdI
            facI = facI * m
            ENDDO

C  No conversion of input y ?
         IF (qNewY) THEN
            idwfI = 0
            isymI = 0
            iBosI = 0
            aM = 1. ! needed for CalcAB (..yA..) even if yA isn't needed

C  Anaylse the input :
         ELSE
 112        CONTINUE
            CALL FindN (CoYin, 1, niIn, iIn)

            IF (CoYin(1:1).eq.'S' .or. CoYin(1:1).eq.'I' ) THEN
               iSymI = iOlfG (j, '?det-bal-sym', Fehler)
               IF (Fehler.ne.'&ff') RETURN
            ELSEIF (CoYin(1:3).eq.'X''''') THEN ! corrected 25jan00
               iBosI = 1
               iBetI = -1
            ELSEIF (CoYin(1:2).eq.'g(') THEN
               idwfI =
     * intq(qAsk(' Multiply input by dwf (not for Placzek) ?'))
               mOrdI = 1
               iBosI = 1
            ELSEIF (CoYin(1:4).eq.'g_#(') THEN
               mOrdI = 1 ! ignore iIn -> convert as if it were g_1
               idwfI = 1
               iBosI = 1
            ELSEIF (CoYin(1:6).eq.'g_(#)(') THEN
               mOrdI = 1 ! ignore iIn -> convert as if it were g_1
               iBosI = 1
            ELSEIF (CoYin(1:2).eq.'G(') THEN
               mOrdI = 1
               idwfI =
     * intq(qAsk(' Multiply input by dwf (not for Placzek) ?'))
            ELSEIF (CoYin(1:4).eq.'G_#(') THEN
               mOrdI = iIn
               idwfI = 1
            ELSEIF (CoYin(1:6).eq.'G_(#)(') THEN
               mOrdI = iIn
            ELSE
               Print *, ' Input y coordinate = '//CoYin(1:len(CoYin))
               aus = ' Is of type S(q,w), g(w), G(w), ...'
               CALL Gong (1)
               CALL FrageC (aus, CoYin)
               IF (CoYin.eq.' ') THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               GOTO 112
               ENDIF

            ! - displacement :
            IF (idwfI.ne.idwfO) THEN
               IF (lj.eq.1)
     * qDM = qAskD (' Displacement from Debye-model', intq(qDM))
               IF (qDM) THEN
                  TD = rAskD ( ' Debye temperature', TD)
                  T0 = 4. ! normalization temperature
                  CALL NiceNum (TD, h2, i2)
                  CALL Append (doc, '; TD='//h2)
               ELSE
                  u2T = rAskD (' Displacement <ux^2>(T) ', u2T)
                  u20 = rAskD (' Of which offset <ux^2>(0)', u20)
                  CALL NiceNum (u2T, h2, i2)
                  CALL Append (doc, '; u^2='//h2)
                  CALL NiceNum (u20, h2, i2)
                  CALL Append (doc, ','//h2)
                  ENDIF
               ENDIF

            ENDIF

         Print *,  ' -> '//CoYout
         CALL OlfCnuP (jout, 'y', CoYout, UnYout, Fehler)
         CALL Compose2 (doc, 'dosconv/ '//CoYin, ' -> '//CoYout)
         CALL OlfComAdd (jout, 'O', doc, Fehler)
         CALL iOlfP (jout, '?det-bal-sym', iSymO, Fehler)

C  DEBUG :
         Print '(a,2i3)', ' Conversion I->O : sym ', iSymI, iSymO
         Print '(a,2i3)', ' Conversion I->O : dwf ', idwfI, idwfO
         Print '(a,2i3)', ' Conversion I->O : ord ', mOrdI, mOrdO
         Print '(a,2i3)', ' Conversion I->O : bos ', iBosI, iBosO
         Print '(a,2i3)', ' Conversion I->O : bet ', iBetI, iBetO

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, nC, X, Y, D, Fehler)
            IF (mOrdI.ne.mOrdO .or. idwfI.ne.idwfO .or.
     *    iQQQO.eq.1 .or. iUUUO.eq.1) CALL CalcQQ (j, K, yQ,nC,Fehler)

            T  = rzOlfGG (j, K, 'T', 'K', Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (Fehler.ne.'&ff') RETURN
            IF (T.lt.0.1) THEN ! limit was 10
               Fehler = 'temperature too low'
               RETURN  ! else crash in dexp(yB/) 11mai93
               ENDIF

C  Atomic masses needed ?
            IF (mOrdI.ne.mOrdO) THEN
               aM = rOlfGG (j, 'at-mass', 'amu', Fehler)
               IF (Fehler(1:4).eq.'&pnf') THEN
                  Fehler = 'please define r-parameter "at-mass"'
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN
               IF (aM.le.0.) THEN
                  Fehler = 'at-mass <= 0'
                  RETURN
                  ENDIF
               c2M   = aM * 931.48d9  ! amu * c^2 -> meV
               chbar = 1.97329d6      ! hbar * c [meV * A]
               DO i = 1, nC
                  yA(i) = chbar**2 * yQ(i)**2 / ( 2 * c2M * TkB)
                  ENDDO
               ENDIF

            TkB   = T  * .0861733 ! K -> meV
            DO i = 1, nC
               yB(i) = X(i)*facX / TkB
               ENDDO

C  Displacements from Debye temperature ?
            IF (idwfI.ne.idwfO) THEN
               ix0 = irPos (X, nC, 0.d0, 'n')
               yQ0 = yQ(ix0)
               IF (qDM) THEN
                  u2T = u2Debye (T,  TD, aM, 1.d-3)
                  u20 = u2Debye (T0, TD, aM, 1.d-3)
                  IF (K.eq.1) THEN
                     ! notify :
                     Print '(1x,2(a,f5.1),2(a,f7.5))',
     *                  ' TD = ', TD, ' T = ', T,
     *                  ' ==> u2(0) = ', u20, ', u2(T) = ', u2T
                     ENDIF
                  ENDIF
               ENDIF

C  Read real parameters for conversion Tof/chs ->  energy
            IF (qCtoE) THEN
               FPath = rOlfGG (j,'FP', 'm', Fehler)
               IF (Fehler(1:4).eq.'&pnf') THEN
                  Fehler = 'please define r-parameter "FP"'
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN
               IF (FPath.le.0.) THEN
                  Fehler = 'FP <= 0'
                  RETURN
                  ENDIF
               Cwidth = rOlfGG (j, 'Cw', 'usec', Fehler)
               IF (Fehler(1:4).eq.'&pnf') THEN
                  Fehler = 'please define r-parameter "Cw"'
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN
               IF (Cwidth.le.0.) THEN
                  Fehler = 'Cw  <= 0'
                  RETURN
                  ENDIF
                  Eelast = rOlfGG (j, 'E0', 'meV', Fehler)
               IF (Fehler(1:4).eq.'&pnf') THEN
                  Fehler = 'please define r-parameter "E0"'
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN
               IF (Cwidth.le.0.) THEN
                  Fehler = 'E0  <= 0'
                  RETURN
                  ENDIF
               Telast = FPath / 3956.d-6 * dsqrt0(81.805/Eelast)
               CALL OlfCnuP (jout, 'x', 'w', 'meV', Fehler)
            ENDIF


C  Conversion yIn -> S(q,w) :

            DO i = 1, nC

               ! Prepare fact :
               fact = 1 ! out = in * fact

               IF (isymI.ne.isymO)
     *            fact = fact * dexp((iSymI-iSymO)*iSEG * yB(i)/2)

               IF (idwfI.ne.idwfO)
     * fact = fact *
     * dexp((idwfI-idwfO)*(-u2T*yQ(i)**2 + u20*yQ0**2))

               IF (mOrdI.ne.mOrdO)
     * fact = fact * (yA(i)/3)**(mOrdI-mOrdO) / facI

               IF (iBosI.ne.iBosO)
     * fact = fact*(-iSEG*yB(i)*(dexp2(-iSEG*yB(i))-1))**(iBosO-iBosI)

               IF (iBetI.ne.iBetO)
     * fact = fact * (-iSEG*yB(i))**(iBetO-iBetI)

               IF (iQQQO.eq.1)
     *            fact = fact * yQ(i)

               IF (iCCCO.eq.1)
     *            fact = fact * 8.3144 * AtNu * yB(i)**3 * (-iSEG) *
     *                   dquot0 (1.d0, 1-dexp2(iSEG*yB(i)))

               IF (iTTTO.eq.1)
     *            fact = fact / T**3

               IF (iUUUO.eq.1)
     *            fact = fact / yQ(i)**2 * (dexp2(-iSEG*yB(i)) + 1)

               ! Convert data (or set y) :
               IF (qNewY) THEN
                  Y(i) = fact
                  D(i) = 0

C ------------------------------------------------------------------
C      Conversion from Tof to energy analog conversion at read in
C ------------------------------------------------------------------
               ELSEIF (qCtoE) THEN
                  R1 = X(i)/(-iSEG)/Telast*FPath*1000
                  E = Eelast*(1.d0/R1/R1-1)
                  fact = .5*Telast/Eelast*R1*R1*R1/Cwidth*R1 !The last R1 for Ki/Kf
                  X(i) = -iSEG * E
                  Y(i) = fact * Y(i)
                  D(i) = fact * D(i)

               ELSE
                  Y(i) = fact * Y(i)
                  D(i) = fact * D(i)
                  ENDIF
               ENDDO ! nC

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, nC, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! K
C  -----------------------------------------------------------------
C    Detector efficiency correction analog RRT_DetEff
C    IN 5 is not possible at the moment
C  -----------------------------------------------------------------

            IF(qCtoE) THEN
                  qDetEff = qAskD('Detector efficiency correction',
     *                            intq(qDetEff))
                  IF (qDetEff) THEN
                     Print *, 'Correction for detector efficiency'
                     CALL OlfComLinG (jout, 1 , Inst, Fehler)
                  IF (Inst(1:3).eq.'IN5') THEN
                 qTakeCare = qAskD('Check det eff really working',
     *                 intq(qTakeCare))
                  ENDIF
                     DO K = 1, nK ! set parameters
                     IF (Inst(1:3).eq.'IN4') THEN
                        DetE1(K) = - .0887
                        DetE2(K) = -5.597
C  new part FK May 03 be carefull !!!
                     ELSEIF (Inst(1:3).eq.'IN5') THEN
                     IF (qTakeCare.and.(iDetTest(K).eq.0)) THEN
                                ! monitor or special block
                          !Angle(K) = -1
                        ELSEIF (qTakeCare.and.(iDetTest(K).eq.1)) THEN
                                ! IN5 type detectors (8bars,D = .9cm)
                           DetE1(K) = - .0887
                           DetE2(K) = -4.07
                        ELSEIF (qTakeCare.and.(iDetTest(K).eq.2)) THEN
                                ! IN6 type detectors
                           DetE1(K) = -.0565
                           DetE2(K) = -3.284
                        ELSEIF (qTakeCare.and.(iDetTest(K).eq.3)) THEN
                                ! IN6 type detectors
                           DetE1(K) = -.0565
                           DetE2(K) = -3.284
                        ELSEIF (qTakeCare.and.(iDetTest(K).eq.4)) THEN
                                ! IN6 type detectors
                           DetE1(K) = -.0565
                           DetE2(K) = -3.284
                        ! Multidetector not yet implemented

                        ENDIF
                     ELSEIF (Inst(1:3).eq.'IN6') THEN
                        DetE1(K) = -.0565
                        DetE2(K) = -3.284
                     ELSEIF (Inst(1:3).eq.'DCS') THEN
                        Print *, 'Check coefficients !'
                        DetE1(K) = -.0130 ! Has to be checked !!
                        DetE2(K) = -4.038
                     ELSEIF (Inst(1:3).eq.'TOF') THEN
                        DetE1(K) = -0.31
                        DetE2(K) = -9.3518
                      ENDIF
                     ENDDO ! K



                     CALL RRT_DetEff (jout, nK, Eelast, DetE1,
     *                                DetE2, Fehler)
                  ELSE
                     Print *, 'no detector efficiency correction'
                  ENDIF
            ENDIF
C  -----------------------------------------------------------------
         CALL OlfClos (jout, nK, Fehler)

         ENDDO ! lj

      END ! DOSconv

C  --------------------------------------------------------------------
      SUBROUTINE DOS (nJList, JList, Fehler)
C  --------------------------------------------------------------------
            ! JWu 14-?nov91. Revision 5-13feb92.
            ! Placzek-Expansion 30-31mar92.
            ! Arbitrary order 27-28apr92.
            ! Extrapolation from larger region 7feb95

         ! Iterative multiphonon subtraction for the
         ! density of states g(w).
         ! Before use, interpolate to constant w with X(1)=0.

            ! js        = uncorrected S(q,w)                        ! Y0
            ! jg        = 1-Phonon guess                            ! Y, Y(,1)
            ! JJ(2..mMcal) = m-Phonon contribution calculated from jg  ! Y(,m)
            ! JJ(1)     = new 1-Phonon guess = js - JJ(2) - ..      ! Y(,1)

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER        (MaxJ=10)

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      DIMENSION         JJ(MaxJ), yQ(MC), yA(MC), yB(MC), BoN(MC),
     *                  X0(MC), Y0(MC), D0(MC),
     *                  G0(-MC:MC), Gc(-MC:MC), Gm(MC),
     *                  Fak(0:MaxJ), Bin(0:MaxJ,0:MaxJ),
     *                  Gu2M(0:MaxJ), YY(MC,MaxJ), NConv(MaxJ),
     *                  qLisK(MK), qLisKold(MK), KKG(MK)

      CHARACTER*80      aus, aus2, h, cLisKg, cLisKd, ord, yGuess
      CHARACTER*40      Co, Un

      DATA              iExtra /1/, ig /1/, mMcal /5/, cofa /0.8d0/,
     *                  goalMatch /6.d0/, nloops /12/,
     *                  cLisKg /'1'/, cLisKd /'*'/,
     *                  qPlaczek /.false./, qu2Fix /.false./,
     *                  iSave /1/, qSaveExp /.false./

C  Set m! and (m over k) :
      Fak(0) = 1.
      DO i = 1, MaxJ
         Fak(i) = Fak(i-1) * i
         DO ii = 0, i
            Bin(i,ii) = Fak(i) / Fak(ii) / Fak(i-ii)
            ENDDO
         ENDDO

C  Global Setup :
      qPlaczek= qAskDi (' Multiphonon (0) or Placzek (1) expansion',
     *                     intq(qPlaczek))

      iSave   = iAskDMu (
     *   ' Save as g_(q,w) (1) or as G_(w) (2) (for m>1)', iSave, 1, 2)
      qSaveExp= qAskD (
     *   ' Save 1-phonon-term as (exp minus corr)', intq(qSaveExp))

C  Loop over input files :
      DO lj = 1, nJList
         jin = JList(lj)

         ! create new output files :
         DO m = 1, MaxJ
            CALL OlfHeadDup (jin, .false., JJ(m), nK, Kout, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

C  Input files - S(q,w) :
         Print *, ' Uncorrected S(q,w) from file '//cl2(jin)
         CALL OlfCnuG (jin, 'y', Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (Co(1:1).ne.'S') THEN
            Fehler =  'y-Coord is not S(..) but '//Co
            RETURN
            ENDIF

         CALL OlfCnuG (jin, 'x', Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (Co.ne.'w') THEN
            Print *, ' input x coordinate = '//Co
            IF (.not.qAsk(' Continue ?')) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF
         CALL UnitConv (Un, 'meV', facX, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         aM  = rOlfGG (jin, 'at-mass', 'amu', Fehler) ! for g-scale, <u^2>,...
         IF (Fehler.ne.'&ff') RETURN
         CALL NiceNum (aM, h, ih)
         Print *, ' atomic mass = '//h(1:ih)

         TS  = rOlfGG (jin, 'T', 'K', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         iSEG = iOlfG (jin, '@sam-erg-gain', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (iSEG.ne.-1) THEN
            Fehler = ' works only for neutron energy gain convention'
            RETURN
            ENDIF

         qSym = qOlfG (jin, '?det-bal-sym', Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Input files - g(w) guess ?
         aus = ' Guess for g(w) from file no. (0=from S, -1=escape)'
         jg = iAskDMu (aus, jg, -1, MemBlockInq('nF'))
         IF (jg.eq.-1) THEN
            jg = 0
            Fehler = ' '
            RETURN
            ENDIF
         IF (jg.ge.1) THEN
            ! Check parameters :
            nKold = iOlfG (jg, '#spectra', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (nKold.gt.nK) THEN
               Fehler = ' guess contains too many spectra'
               RETURN
            ELSEIF (nKold.lt.nK) THEN
 132           CONTINUE
               CALL Compose2 (aus, ' Name the '//cl3(nKold),
     *            ' spectra for which a guess is given')
               CALL GetNList (aus, cLisKg, qLisKold, nK)
               IF (iqSum(qLisKold,nK).eq.0)     GOTO  92
               IF (iqSum(qLisKold,nK).ne.nKold) GOTO 132
               ENDIF
            CALL OlfCnuG (jg, 'y', Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF     (Co(1:1).eq.'g') THEN
               yGuess = 'g'
            ELSEIF (Co(1:1).eq.'G') THEN
               yGuess = 'G'
            ELSE
               Fehler =
     * ' y-coord of 1-phonon guess is neither g(..) nor G(..)'
     * //' but '//Co
               RETURN
               ENDIF
         ELSE ! jg=0
            nKold = nK
            ENDIF

C  Loop settings :
         qWorked = .false.
 10      CONTINUE

         CALL rAskPair (
     *      ' Offset <u2(0)> (always fixed) and initial/fixed <u2(T)>',
     *                   u20, u2T, u20, u2T)
         IF (u2T.ne.0) qu2Fix = qAskD (' Keep <u2> fixed',intq(qu2Fix))

         CALL rAskRge (' Low w extrapolation from average over region',
     *        wOnset, wAvge, wOnset, wAvge)
         wDebye =
     * rAskDMu (' Fixed Debye energy (else 0)', wDebye, 0.d0, 1.d10)
         wCut = rAskDMu (' High w cut-off', wCut, wAvge, 1.d10)

         ig      = iAskDMu (' Grouping for convolution ', ig, 0, MC/4)
         IF (ig.eq.0) GOTO 92 ! emergency exit

         IF (mMsav.eq.0) mMsav = mMcal
 113     CONTINUE
         CALL i2FrageD (' Highest order internally / on exit',
     *        mMcal, mMsav, mMcal, mMsav)
         IF (.not.(1.lt.mMsav .and. mMsav.le.mMcal
     *             .and. mMcal.le.MaxJ)) THEN
            CALL Gong (3)
            GOTO 113
            ENDIF

         goalMatch = rAskDMu (' Required precision (-lg)',
     *                        goalMatch, 0.d0, 3.d1)
         aus  = ' Convergence factor (fraction of second but last g_1)'
         cofa = rAskDMu (aus, cofa, 0.d0, 1.d0)
         cu   = cofa

         aus=
     * ' Maximum number of iterations (-1=quit, 0=correct settings)'
         nLin = iAskDMu (aus, nloops, -1, 100)
         IF (nLin.eq.-1) GOTO 92 ! Notausstieg
         IF (nLin.eq. 0) GOTO 10 ! settings
         nloops = nLin

         IF (nK.gt.1) THEN
            CALL GetNList (' Iterations for spectra', cLisKd, qLisK, nK)
            nKout = iqSum (qLisK,nK)
            IF (nKout.le.0) GOTO 92 ! Notausstieg
         ELSE
            qLisK(1) = .true.
            nKout    = 1
            ENDIF

         ! assignement of spectra j <-> jg :
         IF (nKold.lt.nK) THEN
            DO K = 1, nK
               IF (qLisKold(K)) THEN
                  Kin = K
                  GOTO 134
                  ENDIF
               ENDDO
 134        CONTINUE
            DO K = 1, nK
               IF (qLisKold(K)) Kin = K
               KKG(K) = Kin
               ENDDO
         ELSE
            DO K = 1, nK
               KKG(K) = K
               ENDDO
            ENDIF

C  Now that mMsav is given we may open the output files :
         ! half a doc line about the iteration's setting :
         CALL NiceNum (wOnset, h, ih)
         aus2 = 'w = '//h(1:ih)
         CALL NiceNum (wAvge, h, ih)
         CALL Append (aus2, '/'//h(1:ih))
         CALL NiceNum (wCut, h, ih)
         IF (qPlaczek) THEN
            CALL Append (aus2, '-'//h(1:ih)//',  m <='//cl3(mMcal))
         ELSE
            CALL Append (aus2, '-'//h(1:ih)//', (m)<='//cl3(mMcal))
            ENDIF
            IF (wDebye.gt.0) THEN
               CALL NiceNum (wDebye, h, ih)
               CALL Append (aus, ' wD='//h(1:ih))
               ENDIF
         IF (ig.gt.1) CALL Append (aus2, ' ig='//cl3(ig))

         ! loop over output files :
         DO m = 1, mMsav

            ! coordinate name :
            IF     (m.eq.1 .and. qSaveExp) THEN
               ord = ' '
            ELSEIF (qPlaczek) THEN
               CALL Compose2 (ord, '_('//cl2(m), ')')
            ELSE
               ord = '_'//cl2(m)
               ENDIF
            IF     (iSave.eq.1) THEN
               CALL Compose2 (Co, 'g'//ord, '(q,w)')
            ELSE
               CALL Compose2 (Co, 'G'//ord, '(w)')
               ENDIF
            CALL OlfCnuP (JJ(m), 'y', Co, ' ', Fehler)

            ! comment line :
            IF (m.eq.1 .and. qSaveExp) THEN
               aus = 'dos/ g(w) = full - g_2 - ..'
            ELSE
               CALL Compose2 (aus, cl2(m), '-phn guess')
               ENDIF
            CALL Append (aus, '/ '//aus2)
            CALL OlfComAdd (JJ(m), 'g'//cl2(m), aus, Fehler)

            ! linetype
c96/7            i-Par4(15) = min0 (m-1,4)

            CALL OlfClos (JJ(m), 0, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

C  Loop spectra :
         Kout = 0
         DO K = 1, nK
         IF (qLisK(K)) THEN
            Kout = Kout+1

            CALL Compose2 (aus, ' iterations for spectrum '//
     *                     cl3(K), ' :')
            IF (nKout.lt.nK) CALL Append (aus, ' output spectrum '//
     *                                    cl3(Kout))
            Print *, aus

C  Get uncorrected g(w) :
            ! open complete uncorrected S(q,w) and check scale :
               Print *, ' debug 1'
            CALL OlfGetXYD (jin, K, n, X0, Y0, D0, Fehler)
            CALL CalcQQ (jin, K, yQ, n, Fehler)
            T  = rzOlfGG (jin, K, 'T', 'K', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (.not.qEqEps(T,TS)) Fehler =' problems with temperature'
            IF (Fehler.ne.'&ff') RETURN
            ! check x-scale :
               Print *, ' debug 2'
            tol = 1.d-1
            CALL CheckScale (n, X0, tol, qedX0, dX0)
            IF     (.not.qedX0) THEN
               Fehler =
     * ' x-scale not equidistant in S(q,w), spectrum '//cl3(K)
            ELSEIF (dabs(X0(1)).ge.tol*dX0) THEN
               Fehler = ' X(1)<>0 in g(w), spectrum '//cl3(K)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            ! Calculate Bose factor :
               Print *, ' debug 3'
c-tot            CALL CalcAB (n, T, X0, facX, yQ, aM, yA, yB)
            IF (T .le.0.) CALL Absturz ('CalcAB', 'T<=0')
            IF (aM.le.0.) CALL Absturz ('CalcAB', 'aM<=0')

            TkB   = T  * .0861733  ! K -> meV
            c2M   = aM * 931.48d9  ! amu * c^2 -> meV
            chbar = 1.97329d6      ! hbar * c [meV * A]

            DO i = 1, n
               yB(i) = X0(i)*facX / TkB
               yA(i) = chbar**2 * yQ(i)**2 / ( 2 * c2M * TkB)
               ENDDO

               Print *, ' debug 4'
            DO i = 1, n
               BoN(i) = yB(i) * (dexp(yB(i)) - 1.)  ! beta/n(w)
               IF (yA(i).le.0) THEN
                  Fehler = ' alpha<=0 (probably Q=0)'
                  RETURN
                  ENDIF
               ENDDO
            ! transform S(q,w) into uncorrected g(w) (except for DWF(T)):
               Print *, ' debug 5'
            DO i = 1, n
               fact = 3 / yA(i) * BoN(i) * dexp1(-u20 * yQ(1)**2)
               IF (qSym) fact = fact * dexp(-yB(i)/2)
               Y0(i) = fact * Y0(i)
               D0(i) = fact * D0(i)
               ENDDO

            ! calculate G0(0) from wDebye :
            TkB   = T  * .0861733  ! K -> meV
            IF (wDebye.gt.0) GDebye = 9 * TkB**2 / wDebye**3

C  Get 1-phonon-guess :
               Print *, ' debug 6'
            IF (jg.gt.0) THEN
               Print *, ' debug 7'
               ! open 1-phonon guess :
               CALL OlfGetXYD (jg, KKG(K), ng, X, Y1, D, Fehler) ! old: zg
               IF (Fehler.ne.'&ff') RETURN
               ! check x-scale :
               CALL CheckScale (ng, X, tol, qedX, dX)
               IF     (.not.qedX) THEN
                  Fehler =
     * ' x-scale not equidistant in g_1, spectrum '//cl3(K)
               ELSEIF (dabs(X(1)).ge.tol*dX) THEN
                  Fehler = ' X(1)<>0 in g_1, spectrum '//cl3(K)
               ELSEIF (ng.ne.n) THEN
                  Print *,
     * ' WARNING/ different length of spectra '//cl3(K)
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN
               IF     (yGuess.eq.'G') THEN
                  DO i = 1, ng
                     Y1(i) = BoN(i) * Y1(i)   ! g_1(w) <- G_1(w)
                     ENDDO

               ELSEIF (yGuess.eq.'S') THEN
                  Fehler = 'guess=S_ not yet implemented'
                  RETURN
                  ENDIF
               ENDIF

            ! wCut -> nCut :
               Print *, ' debug 8'
            nCut   = irPos (X0, n, wCut,   'r') - 1
            nOnset = irPos (X0, n, wOnset, 'r')
            nAvge  = irPos (X0, n, wAvge,  'r')
            IF     (nCut.le.3) THEN
               Fehler = ' Less than 4 channels in allowed x-range'
            ELSEIF (nCut.lt.2*nAvge) THEN
               Fehler =' Extrapolating over more than half the w-range'
            ELSEIF (nOnset.le.0) THEN
               Fehler = ' low w cut-off too small, 0. not excluded'
            ELSEIF (nAvge-nOnset.lt.1) THEN
               Fehler = ' Void region from which to extrapolate'
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            IF (nAvge-nOnset.gt.1 .and. .not.qPlaczek)
     * Print *, ' WARNING/ in calculation of <u^2> no averaging'

            ! n -> NConv :
               Print *, ' debug 9'
            DO m = 1, mMcal
                NConv(m) = min0 (nCut*m, n-1)
                ENDDO
            DO m = 3, mMcal
               NConv(m) = NConv(2)  ! Sparmassnahme
               ENDDO
            nMax = NConv(mMcal)

            ! initialize by 0. :
            DO i = 1, n
               D (i) = 0. ! no errors for m-phonon-contributions
               DO m = 1, MaxJ
                  YY (i,m) = 0. ! for m-phonon contribution outside NConv(m)
                  ENDDO
               ENDDO
            DO i = -n, n
               G0 (i) = 0.
               ENDDO

C  Loop iterations :
               Print *, ' debug 10'
            niter = nloops
            iter  = 0
 101           CONTINUE
               iter = iter + 1
               Print *, ' iter : ', iter

               IF (.not. (iter.eq.1 .and. jg.gt.0)) ! else Y1 has been set above
     *         CALL rCopy (Y1, 1, n, 1, Y,  1, 1) ! guess < result
               CALL rCopy (Y, 1, n, 1, Y0,  1, 1) ! tractandum < raw g(w)
               CALL rCopy (D, 1, n, 1, D0,  1, 1)

               IF (.not.qPlaczek) THEN ! get <u^2> (simplified 7feb95)
                  IF     (qu2Fix .or. iter.eq.1) THEN
                     u2TK = u2T
                  ELSE
                     u2TK = cu*u2TK + (1-cu)*u2 ! damping with convergence factor
                     ENDIF
cdeb                  Print *, 'u2T u2TK cu u2', u2T, u2TK, cu, u2   ! <<<
                  ! correct raw g(w) by DWF :
                  DO i = 1, n
                     fact = dexp1((u2TK)*yQ(i)**2)
                     Y(i) = fact * Y(i)
                     D(i) = fact * D(i)
                     ENDDO
                  ENDIF

               ! convergence factor : mix last and 2nd but last guess (12feb92)
               c = cofa
               IF (iter.eq.1) c = 0.

               ! guess g_1(w) -> full G0(w) :
               DO i = 2, nCut+1
                  G0( i-1) =   c  * G0( i-1) +
     * (1-c) * Y1(i) /  yB(i)  / (dexp( yB(i)) - 1.)
                  G0(-i+1) =   c  * G0(-i+1) +
     * (1-c) * Y1(i) /(-yB(i)) / (dexp(-yB(i)) - 1.)
                  ENDDO

               ! fill region around w=0 by extrapolation :
               IF (wDebye.le.0) THEN ! G0=const
                  Gleft  = 0
                  Gright = 0
                  DO i = nOnset, nAvge
                     Gleft  = Gleft  + G0(-i) / (nAvge-nOnset+1)
                     Gright = Gright + G0( i) / (nAvge-nOnset+1)
                     ENDDO
                  DO i = -nOnset+1, nOnset-1
                     G0(i) =
     * Gleft + (i+nOnset)*(Gright-Gleft)/(2*nOnset)
                     ENDDO
               ELSE ! G0 = linear
                  Gleft  = 0
                  Gright = 0
                  DO i = nOnset, nAvge
                     Gleft  =
     * Gleft  + (G0(-i)-GDebye)/i/(nAvge-nOnset+1)
                     Gright =
     * Gright + (G0( i)-GDebye)/i/(nAvge-nOnset+1)
                     ENDDO
                  DO i = 0, nOnset-1
                     G0(-i) = GDebye + i * Gleft
                     G0( i) = GDebye + i * Gright
                     ENDDO
                  ENDIF

               ! copy (for iterative convolution: Gc = G0 * G0 * ..) :
               DO i = -n, n
                  Gc(i) = G0(i)
                  ENDDO

C  Convolution -> YY (=G_m for multiphonon, =G_(m) for Placzek) :

               IF (qPlaczek) THEN
                  Gu2 =  u2TK * 3 * yQ(1)**2 / yA(1)
                  IF (Gu2.gt.1.d5) THEN
                     Print *, ' jg iter u2TK Q**2 A : ',
     *                  jg, iter, u2TK, yQ(1), yA(1)
                     Fehler = ' <u2> overflow/ G_0 > 1e5'
                     RETURN
                     ENDIF
                  Gu2M(0) = 1. ! to prevent FORTRAN error 0**0
                  DO mm = 1, mMcal
                     Gu2M(mm) = (-Gu2)**(mm)
                     ENDDO
                  DO mm = 2, mMcal
                     DO i = 1, NConv(mm)+1
                        YY(i,mm) = Bin(mm,1) * G0(i-1) * Gu2M(mm-1)
                        ENDDO
                     ENDDO
                  ENDIF

               DO m = 2, mMcal
                  ! convolution G0*Gc :
                  yc = 0.
                  DO i = 0, NConv(m)
                     IF (mod(i,ig).eq.0) THEN
                        yc = 0.
                        DO ii = max0 (-nCut, -NConv(m-1)+i),
     *                          min0 ( nCut,  NConv(m-1)+i), ig
                           yc = yc + G0(ii) * Gc(i-ii)
                           ENDDO
                        yc = yc * dX0 * ig
                        IF (dabs(yc).gt.1.d16) THEN
                           Print *, ' m, i, yc :', m, i, yc
                           Fehler =
     * ' Convolution/ Arithmetic overflow imminent'
                           RETURN
                           ENDIF
                        ENDIF
                     Gm(1+i) = yc
                     ENDDO
                  DO i = 0, NConv(m)
                     Gc( i) = Gm(1+i)
                     Gc(-i) = Gm(1+i)
                     ENDDO
                  ! -> Gc = G_m

                  IF (qPlaczek) THEN
                     DO mm = m, mMcal
                        DO i = 1, NConv(m)+1
                           YY(i,mm) = YY(i,mm) +
     *                       Bin(mm,m) * Gc(i-1) * Gu2M(mm-m)
                           ENDDO
                        ENDDO
                  ELSE
                     DO i = 1, NConv(m)+1
                        YY(i,m) = Gc(i-1)
                        ENDDO
                     ENDIF
                  ENDDO ! m

C  Conversion G -> g, subtraction from g(0) :
                         ! b / n              for restoring g from G,
                         ! 1/m! * (a/3)^(m-1) prefactor in m-phonon expansion
               DO m = 2, mMcal
                  DO i = 1, NConv(m)+1
                     IF     (dabs(YY(i,m)).gt.1.d16) THEN
                        Print *, ' NConv, m, i, YY :',
     *                            NConv(m), m, i, YY(i,m)
                        Fehler = ' G->g/ Arithmetic overflow imminent'
                        RETURN
                        ENDIF
                        Y(i) = Y(i) -
     *                      BoN(i)*(yA(i)/3)**(m-1)/Fak(m) * YY(i,m)
                     ENDDO
                  ENDDO
               DO i = 1, NConv(mMcal)+1
                  Y(i) = dmax1 (Y(i), 0.d0)
                  ENDDO

C  Save G or G_1 ?
               IF (qSaveExp) THEN
                  DO i = 1, n
                     ii = i
                     IF (i.eq.1) ii = 2
                     YY(i,1) = Y(ii) / BoN(ii) ! g -> G
                     ENDDO
               ELSE
                  DO i = 1, NConv(1) + 1
                     YY(i,1) = G0(i-1)
                     ENDDO
                  ENDIF

C  Convert the YY for saving as g_ ?
               IF     (iSave.eq.1) THEN ! contribution S_m(q,w)
                  DO m = 1, mMcal
                     DO i = 1, n  ! NConv(m)+1 but n for (m=1 and qSaveExp=true)
                        YY(i,m) =
     * BoN(i)*(yA(i)/3)**(m-1)/Fak(m)* YY(i,m)
                        ENDDO
                     ENDDO
                  ENDIF

C  Calculate new <u2> (30mar92) :
               u2 = 0.
               DO i = 1, nCut
                  ii = i
                  IF (ii.le.nOnset) ii = nOnset+1 ! extrapolation w->0
                  u2 = u2 + dX0 * Y(ii) / dtanh(yB(ii)/2) / X0(ii)
                  ENDDO
               u2 = u2 / 0.2393 / aM / 6   ! hbar^2/amu = meV A^2 / .2393

C  Estimate quality of match (6/13feb92) :
               nratio = 0
               ratio  = 0.
               box    = 0.
               squares= 0.
               igg = max0 (1, ig)
               DO i = igg, nMax, igg
                  IF (Y1(i).gt.0.) THEN
                     nratio = nratio + 1
                     ratio  = ratio + Y(i) / Y1(i)
                     ENDIF
                  box     = box     +  Y(i) * Y1(i)
                  squares = squares + (Y(i) - Y1(i))**2
                  ENDDO
               IF (nratio.le.0) THEN
                  Print *, ' old guess g_1(w) == 0.'
                  ratio = 0.
               ELSE
                  ratio = ratio / nratio
                  ENDIF
               IF (box.le.0.) THEN
                  Print *, ' new guess g_1(w) == 0. ?!'
                  squares = 0.
               ELSE
                  squares = squares / box
                  ENDIF
               squares = - dlg0 (squares)

C  Notify end of iteration :
               IF (iter.eq.1) Print *
               Print '(a,i3,a,f7.4,a,f4.1,a,f8.5)',
     * '  end of iteration', iter,
     * ', ratio =', ratio, ', -lg(match) =', squares, ', <u2>=', u2

               IF (squares.ge.goalMatch) THEN
                  ! '  found stationary solution'
                  GOTO 109
                  ENDIF

               IF (iter.ge.niter) THEN
                  aus = '  How many more iterations (-1=settings)'
                  nmore = iAskDMu (aus, niter/2, -1, 100)
                  IF     (nmore.eq.-1) THEN
                     GOTO 10
                  ELSEIF (nmore.le.0) THEN
                      GOTO 108 ! exit without solution
                  ELSE
                     niter = niter + nmore
                     ENDIF
                  ENDIF

               GOTO 101 ! loop : next iteration
C  End loop iterations.

 108        CONTINUE ! found no stationary solution
               ! tant pis
 109        CONTINUE ! save results for present K


            DO m = 1, mMsav
cdeb           IF (qAsk('Exit ?')) RETURN  ! <<<
           Print *, ' going to save in J=', JJ(m), ' K=', K, 'n=', n
               CALL OlfOpen (JJ(m), 1, Kdummy, Fehler)
                 IF (Fehler.ne.'&ff') RETURN
               DO i = 1, n ! statt YY(1,m)
                  Y4(i) = YY(i,m)
                  ENDDO
               CALL OlfCopZ (jin, JJ(m), K, K, Fehler)
                 IF (Fehler.ne.'&ff') RETURN
               CALL OlfPutXYD (JJ(m), K, n, X0, Y4, D,  Fehler)
                 IF (Fehler.ne.'&ff') RETURN
               CALL OlfClos (JJ(m), nK, Fehler)
                 IF (Fehler.ne.'&ff') RETURN
               ENDDO
            Print *, ' end saving g_()'

            IF (niter.ge.2*nloops .and. Kout.lt.nKout) THEN
               nloops =
     * iAskD ('  Increase number of iterations', nloops)
               IF (nloops.le.0) GOTO 11
               ENDIF

            ENDIF ! qLisK(K)
            ENDDO ! K
C  End loop spectra.

         qWorked = .true.
         jg = JJ(1)
         CALL qCopy (qLisKold, 1, nK, 1, qLisK, 1, 1)
         nKold = nKout
 11      CONTINUE

         IF (ig.gt.1)     Print *,
     * ' NOTE/ grouping factor still > 1'
         IF (nKout.lt.nK) Print *,
     * ' NOTE/ some spectra are still uncorrected'
         IF (qAsk (' More iterations with new settings ?')) GOTO 10
C  End loop settings.

 90      CONTINUE ! recover from error (Notausstieg)
         IF (qWorked) THEN
            ! close files :
            IF (Fehler.ne.'&ff') CALL FehlerGong (Fehler,1)
            Print *, ' m-phonon contributions ->'
            DO m = 2, mMsav
               CALL OlfClos (JJ(m), nK, Fehler)
               ENDDO
            Print *, ' 1-phonon contribution : new estimate ->'
            CALL OlfClos (JJ(1), nK, Fehler)
         ELSE
            Print *, 'did not work ??'
            ENDIF

         ENDDO ! lj

      write (aus, '(g14.7,l4)') goalMatch, qPlaczek ! exorcism of a compiler bug

C  End loop over files.
      RETURN

 92   CONTINUE ! Notausstieg
      Fehler = 'Ruecksprung leider unerwuenscht'
      RETURN
c      GOTO 90

      END ! DOS


C   ------------------------------------------------------------------
      SUBROUTINE MUPHDOS (nJList, JList, Fehler)
C   ------------------------------------------------------------------
      ! In Anlehnung an MUPHOCOR von W. Reichardt
      ! implementierte Version einer Zustandsdichteberechnung
      ! F. Kargl Feb/Mar 2003

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER        (MaxJ=14)

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      DIMENSION     JJ(MaxJ), YM0(MC), XM0(MC), DM0(MC), G(MC),
     *              Gstr(MC), GAV(MC), G0(MC), F1(MC), GI(MC),
     *              RQC(MC), EF(MC), YEQ(MC),RMAI(20), EXE(MC),
     *              EX(MC), EXA(MC), Y4n(MC), GS(MC), GD(MC),
     *              GM1(20,MC), RFC(20), GP(MC), RQD(20,MC),
     *              GAI(MC), DSIG(MC), RQUO(MC), FPI(MC), RQMI(MC),
     *              RQMA(MC), dQUOFG(MC), GGES(MC), RQUO1(MC), F(MC),
     *              AZ(MK),RSZ(MK),YMK0(MC),XMK0(MC),DMK0(MC),
     *              YMOK0(MC), YMUK0(MC)

      CHARACTER*80  h
      CHARACTER*40  Co,Un

      DATA          nPHO /5/, iVIT /10/,
     *              iDW /1/, iDifQuo /10/, iMax /20/,
     *              dW /0.05/, ngu /25/



C     Input data should be corrected S(2th,w), integrated over angle
C     over angle is done within MUPHDOS
C AmAnfang
      nPHO = iAskDMu('number of multiphonon terms: ',nPHO,1,10)
      iDW = iAskDMu('iteration of DWF (1yes/0no): ',iDW,0,1)
      iMax = iAskDMu('maximal number of iterations: ',iMax,0,20)
      iVIT = iAskDMu('change from Dif to Quot: ',iVIT,0,iMax)
      iDifQuo = iMax + 1
      IF(iVIT.gt.1) iDifQuo = iVIT
      dW = rAskDMu('Debye Waller factor (first guess): ',dW,0.d0,1.d0)
      ngu = iAskDMu('lower ch# for w^2 extrapolation: ',ngu,0,50)
      ALFI = 1.
      iPR = 0
C start not from 0 in energy !

      ngu1 = ngu + 1

C     Loop over input files:
      DO lj = 1,nJList
         jin = JList(lj)

      ! create new output files
         DO m = 1, MaxJ
            CALL OlfHeadDup (jin, .false., JJ(m), nK, Kout, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO



C     Hier fehlt noch einiges !!

C     Input files - S(2th,w) ? :
C            Print *,'Uncorrected S(2th,w) from file '//c12(jin)
            CALL OlfCnuG (jin, 'y', Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (Co(1:1).ne.'S') THEN
               Fehler = 'y-Coord is not S(..) but '//Co
               RETURN
               ENDIF

            CALL OlfCnuG (jin, 'x', Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (Co.ne.'w') THEN
               Print *, 'input x coordinate = '//Co
               IF (.not.qAsk(' Continue?')) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               ENDIF
            CALL UnitConv (UN, 'meV', facX, Fehler)
            IF (Fehler.ne.'&ff') RETURN


            AMASI = rOlfGG(jin, 'at-mass', 'amu', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL NiceNum (AMASI, h, ih)
            Print *, ' atomic mass = '//h(1:ih)
            !! have to be careful by more than one species

            DO m = 1,nPHO
               RMAI(m) = (1.0087/AMASI)**m
               ENDDO ! invers in units of neutron mass

            nK = iOlfG (jin,'#spectra',Fehler)
            nZ = iOlfG (jin,'#Z',Fehler)
            Print *,'nZ = ',nZ
            DO k = 1, nK
               CALL OlfGetZ (jin, k, nZ, Z, Fehler)
               AZ(k) = Z/2
               RSZ(k) = dsind(AZ(k))
C               Print *, 'Test 1'
C               Print *, Z, AZ(k)
C               Print *, ' '
               ENDDO
               Print *, 'nK = ', nK
! not used spectra should be deleted before
! spectra should start at x = 0
            Print *,'jin = ',jin
            Print *,'number of spectra: ',nK
            nku = iAskMu ('number min of spectrum for int', 0, nK)
            nko = iAskMu ('number max of spectrum for int', 0, nK)
            CALL OlfGetZ (jin,nku,nZ,fIMI,Fehler)
            CALL OlfGetZ (jin,nko,nZ,fIMA,Fehler)
            Print *, 'outer:', fIMI,fIMA

            !preparation of sind multiplication
C            DO k = nku,nko
C               CALL OlfGetXYD (jin,k,n,XMK0,YMK0,DMK0,Fehler)
C               DO i = 1,n
C                  YMK0(i) = YMK0(i)*RSZ(k)
C                  DMK0(i) = DMK0(i)*RSZ(k)
C                  ENDDO
C               CALL OlfOpen (jin,1,nK,Fehler)
C               CALL OlfPutXYD (jin,k,n,XMK0,YMK0,DMK0,Fehler)
C               CALL OlfClos (jin,nk,Fehler)
C            ENDDO

            ! now integration upper sum + lower sum div 2
            ! integration is witin limits accurate
            ! the higher the density the better the result
            ! you have to be sure that all spectra are
            ! on the same aequidistant grid !!
            DO i = 1,n
               YMOK0(i) = 0.
               YMUK0(i) = 0.
               ENDDO
            ! command to get n has to be modified soon
            CALL OlfGetXYD (jin,3,n,XMK0,YMK0,DMK0,Fehler)

            ! go further

            DO i = 1,n
               DO k = nku+1,nko
                  CALL OlfGetXYD (jin,k,n,XMK0,YMK0,DMK0,Fehler)
                  YMOK0(i) = YMK0(i)*RSZ(k)*(AZ(k)-AZ(k-1)) + YMOK0(i)
                  ENDDO
               DO k = nku,nko-1
                  CALL OlfGetXYD (jin,k,n,XMK0,YMK0,DMK0,Fehler)
                  YMUK0(i) = YMK0(i)*RSZ(k)*(AZ(k+1)-AZ(k)) + YMUK0(i)
                  ENDDO
                  YM0(i) = (YMOK0(i) + YMUK0(i))/2
               ENDDO

            CALL OlfOpen(JJ(1), 1, Kdummy, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL OlfCopZ(jin,JJ(1),((nko+nku)/2),1,Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL OlfPutXYD(JJ(1),1,n,XMK0,YM0,DMK0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL OlfClos(JJ(1),nK,Fehler)
            IF (Fehler.ne.'&ff') RETURN
            Print *, 'Saved sin theta summed spectrum, ', 1

            jin = JJ(1)
            Print *,'jin neu = ',jin
            IF(.not.qAsk('Continue at 1')) RETURN

C     Teilabschnitt Beginn der Berechnung nach MUPHOCOR
            K = 1
            ! OPEN S(2th,w) summed over angle weighted sin(theta)
            CALL OlfGetXYD (jin, K, n, XM0, YM0, DM0, Fehler)
C           IF (Fehler.ne.'&ff') RETURN
            ! definitions

            n1 = n+1
            n2 = 2*n
            n5 = n2
            ngoo = n
            delE = XM0(5)-XM0(4)
            E0 = rOlfGG(jin, 'E0', 'meV', Fehler)
            ! check for aequidistant grid has to be implemented

            ! check for K = 1
            IF (K.ne.1) THEN
               Print *, 'File containing more than one spectrum !!'
               RETURN
               ENDIF
            DO i = 1, n
               YEQ(i) = YM0(i) * dquot0((E0+XM0(i)),E0)**2
               YM0(i) = YM0(i) * dsqrt( dquot0( (E0 + XM0(i)),E0))
               ENDDO

C **************
            ! Call Subroutine Preparation of several Parameters
            CALL PrepSub(jin,E0,XM0,EF,n,nPHO,EXE,EXA,EX,RQD,RQC,RQMA,
     *                   RQMI,RFC,cOMI,cOMA,cOM,fIMI,fIMA,Fehler)
            IF (Fehler.ne.'&ff') RETURN



C ***************
            !prepare GNU0 from MUPHCOR
            DO i = ngu1, n5
               G0(i) = YM0(i) * EXE(i) * EX(i)/RQD(1,i)
               ENDDO



C     ***** w^2 extrapolation for small omega using mean value ***
            AFA = 0.

            DO i = 1, 3
               AFA = AFA + G0(ngu+i) / (XM0(ngu+i))**2
               ENDDO

            AFA = AFA / 3.
            DO i = 1, ngu
               G0(i) = AFA * XM0(i)**2
               ! hier ist Y0(i) = 0. in MUPHOCOR
               ENDDO

C      ****************** interpolation end ******************

            IF (iVIT .gt. 0) THEN
               DO i = 1, n5
                  IF (G0(i).lt.0.) G0(i) = 0.
                  ENDDO
               ENDIF

C       ***** Another block M^-n out in SUB Prep *************

            DO i = n1, n2
               G(i) = 0.
               ENDDO

C       ******* First guess for g(w) *************************

            DO i = 1, n
               G(i) = G0(i)*dexp(.625*RQC(i)*dW)
               ENDDO

            ! smoothed end towards 0 with (n-i)^2

            DO i = n-20 , n
               G(i) = (G(n-20)/400.)*((n-i)**2)
               ENDDO

C            DO i = 1,n5
C            Print *, 'G(',i,') = ',G(i)
C            ENDDO
C            IF(.not.qAsk('Proceed at test?')) RETURN

            ! calculation of the normalization value gnorm

               gnorm = 0.
            DO i = 1, n
               gnorm = gnorm + G(i)
               ENDDO
            gnorm = gnorm * delE
            ! normalization
            DO i = 1, n2
               G(i) = G(i)/gnorm
               Gstr(i) = 0.
               F1(i) = G(i)
               ENDDO
               Gstr(1) = G(1) * delE
            DO i = 2, n2
               Gstr(i) = Gstr(i-1) + G(i) * delE
               ENDDO
C      ** in MUPHOCOR call Sub Corrfa and Separ *****
C         CALCULATION DONE WITHIN PROGRAM only one species implemented !

            DO i = 1, n5
               dQUOFG(i) = 1.
               FPI(i) = G(i)
               ENDDO !end CORRFA and SEPAR

            iLAU = 0

            ! not yet implemented printing all settings
            Print *, 'Settings:  nPHO',nPHO,' iVIT: ',iVIT
! must be similar to inz.log in i80.f !!!!
! in order to compare all data.


C      * start iterations !!!!
 1001       CONTINUE
            iLAU = iLAU + 1

            PRINT *,'Iteration number: ', iLAU

            IF (iLAU .gt. iDifQuo) iVIT = 0

            ! normalization of calculated g(w)
            IF (iLAU .gt. 1) THEN
               gnorm = 0.
               DO i = 1, n
                  gnorm = gnorm + G(i)
                  ENDDO
               gnorm = gnorm * delE
               DO i = 1, n
                  G(i) = G(i) / gnorm
                  ENDDO
               ! start corrfa and separ
               DO i = 1, n5
                  dQUOFG(i) = 1.
                  FPI(i) = G(i)
                  ENDDO !end CORRFA and SEPAR
               ENDIF


C     ***** Calculation of T_1 = F(w) / (2hw*sinh(hw/2kT)) *****
            DO i = 1, n
               REXE = 1./EXE(i)
               GS(i) = G(i) * REXE
               FPI(i) = FPI(i) * REXE
               ENDDO



C     * Calculation of the Debye-Waller-Factor ******************
C     **** dW = Int(h^2/2M * F(w)/w * cos(hw/2kT))/sinh(hw/2kT)**

            ! Integral runs from 1 up to n as the cut off defined
            ! by the user

            DWI = 0.
            DO i = 1, n
               DWI = DWI + FPI(i) * EXA(i)
               ENDDO
            DWI = DWI * delE * 1.0087 / AMASI

            DW1 = 0.

            DO i = 1, n
               DW1 = DW1 + GS(i) * EXA(i)
               ENDDO
            DW1 = DW1 * delE * RMAI(1)

            ! Print in inz log necessary

C     * for 2 species do-Loop must be integrated at this point

            ! Storage of T_1 from -n, +n
            DO i = 1, n
               GD(i) = FPI(n+1-i)
               ENDDO
            DO i = n1, n2
               GD(i) = FPI(i-n)
               ENDDO

C      * Preparation of multiphon calculation and storage
C      ** GM1(m,i) enthaelt die T_m - Terme ***

       ! first term
       DO i = 1, n2
          GM1(1,i) = FPI(i)
          ENDDO


C      * RFC has to be calculated seperately !!!!


       DO m = 1, nPHO+1
C          Print *, 'm = ',m, 'RFC = ',RFC(m)
          CALL ViPho(Gstr,GP,GD,delE,n,m,Fehler)
            ! VIPHO has to be implemented here
            DO i = 1, n2
               GM1(m+1, i) = GP(i) * RFC(m+1)
               Gstr(i) = GP(i)
               ENDDO
          ENDDO

C ***************************************************
C  ! important testlines with saving of par in muc.log
C ***************************************************

C          IF(.not.qAsk('Teste RFC loop?'))RETURN
C          CALL OpenFile (35,'muc','log','e',Fehler)
C          IF (Fehler.ne. '&ff') RETURN

C          Write (35,'(/a)') ' XM0     G1     G2   '
C          DO i = 1,n
C             Write(35,*) -XM0(n-i),GM1(1,i),GM1(2,i)
C             ENDDO
C          DO i = n1,n2
C             Write (35,*) XM0(i-n),GM1(1,i),GM1(2,i)
C             ENDDO
C             Write (35,'(/a)') ' XM0     G3     G4   '
C          DO i = 1,n
C             Write(35,*) -XM0(n-i),GM1(3,i),GM1(4,i)
C             ENDDO
C          DO i = n1,n2
C             Write (35,*) XM0(i-n),GM1(3,i),GM1(4,i)
C             ENDDO
C             Close(35)
C             IF(.not.qAsk('Exit now after RFC ?'))RETURN



       IF (iDW.eq.1) dW = DWI


C      * Construction of the scattering law using eq 2.4 from
C      * primary report on MUPHOCOR

       DO i = 1, n5
          GAV(i) = 0.
          GGES(i) = 0.

          RQRMA = dW * RQMA(i)
          RQRMI = dW * RQMI(i)
          EXMA = dexp1(- RQRMA)
          EXMI = dexp1(- RQRMI)
          COF0 = EXMI-EXMA


          ! Calculation of I_m(x)
          DO m = 1, nPHO
             COF0 = m*COF0+EXMI*RQRMI**m-EXMA*RQRMA**m
             GAV(i) = GAV(i) + RMAI(m)*RQD(m,i)*GM1(m,i) ! check if correct
             GM1(m,i) = (1. / EX(i)) * COF0*RMAI(m)*GM1(m,i)/
     *                  (2.*dW**(m+1))
             GGES(i) = GGES(i) + GM1(m,i)
             ENDDO


          GAV(i) = GAV(i)/(RQD(1,i)*RMAI(1))
          ENDDO ! end construction part 1


          IF (IFOLD.eq.1) THEN
             ! implement folding if really used
             Print *, ' implement folding if really used'
             ENDIF

          DO i = 1, n5
             GGES(i) = GGES(i)/(RQD(1,i)*RMAI(1))
             GAI(i) = GAV(i)*EXE(i)
             GI(i) = GGES(i)*EXE(i)*EX(i)
             ENDDO



C      ******** Creation of plot file for multiphon contrib !! ***

       ! has to be implemented save phonons in new files only if
       ! ITM reached.
       IF (iLAU.ge.iMAX) THEN
          sumZ = 0.
          sumC = 0.
          Print *,'da und nicht da'
          DO i = 1, n2
             sumZ = sumZ + YM0(i)
             DO m = 1, nPHO
                sumC = sumC + GM1(m,i)
                ENDDO
             ENDDO
          sumZ = sumZ * delE
          sumC = sumC * delE

          DO i = 1, n5
             DO m = 1, nPHO
                ! if two species ALFI = ALFI(s)
                DSIG(i) = DSIG(i) + ALFI*GM1(m,i)
                ENDDO
             ENDDO

          DO l = 5,5 + nPHO
             CALL OlfOpen(JJ(l), 1, Kdummy, Fehler)
             IF (Fehler.ne.'&ff') RETURN
             IF (l.eq.5) THEN
                DO i = 1, n2
                   Y4n(i) = G(i)
C                   Y4(i) = (GM1(1,i) + GM1(2,i) + GM1(3,i) + GM1(4,i))
C     *             /sumC
                   ENDDO
                CALL OlfCopZ(jin,JJ(l),K,K,Fehler)
                IF (Fehler.ne.'&ff') RETURN
                CALL OlfPutXYD(JJ(l), K, n2, XM0, Y4n, DM0,  Fehler)
                IF (Fehler.ne.'&ff') RETURN
                CALL OlfClos(JJ(l), nK, Fehler)
                IF (Fehler.ne.'&ff') RETURN
                Print *, 'Saved multiplot 1, ', l
             ELSE
                DO i = 1, n2
                   ll = l-5
                   Y4n(i) = GM1(ll,i)/sumC
                   ENDDO
                CALL OlfCopZ(jin, JJ(l), K, K, Fehler)
                IF (Fehler.ne.'&ff') RETURN
                CALL OlfPutXYD(JJ(l),K,n2,XM0,Y4n,DM0, Fehler)
                IF (Fehler.ne.'&ff') RETURN
                CALL OlfClos(JJ(l),nK,Fehler)
                IF (Fehler.ne.'&ff') RETURN
                Print *, 'Saved multiplot 1, ', l
                ENDIF ! l save

          ENDDO ! FILE

         ENDIF ! ilau .gt. iTM


C      *******

          !end loop 2 species at this point !!



          DO i = 1, n5
             RQUO(i) = 1.
             GAV(i) = GAI(i)*ALFI
             GGES(i) = GI(i)*ALFI
             IF (GGES(i).ne.0) RQUO(i) = G0(i)/GGES(i)
             ENDDO


          PRINT *,'in it'

          gnorm = 0.
          gnorm2 = 0.
          rQUAV = 0.

          DO i = 1, ngoo
             gnorm = gnorm + GGES(i)
             gnorm2 = gnorm2 + G0(i)
             ENDDO

          DO i = 1, n
             rQUAV = rQUAV + RQUO(i)
             ENDDO

          rgnorm = gnorm / gnorm2
          rQUAV = rgnorm * rQUAV/n

          DO i = 1, n
             Gstr(i) = 0.
             GP(i) = 0.
             RQUO1(i) = 0.


C     * different new guesses for F(w) using quotient or difference
C     * method.

          IF (iVIT.eq.0) THEN !difference method
             G(i)=G(i) + GAV(i)*(rgnorm*RQUO(i)-1.)
             IF (G0(i).eq.0.) G(i) = 0.
             IF (G(i).gt.0) THEN
                GPC = G(i)/(3.*(XM0(i)**2)) ! this has to be checked
                GP(i) = 1. / GPC**.3333333
                IF(F1(i).ne.0.) THEN
                   RQUO1(i) = G(i)/F1(i)
                   ENDIF
                ENDIF
          ELSE ! quotient method
             G(i) = G(i) * RQUO(i) * rgnorm
             IF (G(i).gt.0) THEN
                GPC = G(i)/(3.*(XM0(i)**2)) ! this has to be checked
                GP(i) = 1. / GPC**.3333333
                IF(F1(i).ne.0.)  THEN
                   RQUO1(i) = G(i)/F1(i)
                   ENDIF
                ENDIF
             ENDIF ! different end
             ENDDO


        ! smooth cut off
          DO i =  n1, n2
             G(i) = 0.
             ENDDO
          DO i = n-20,n
             G(i) = (G(n-20)/400.)*(n-i)**2
             ENDDO

C     * calculation of the first moment of f(w)
          DO i = 1, n5
             Gstr(i) = 0.
             ENDDO
          Gstr(1) = G(1) * delE
          DO i = 2,n
             Gstr(i) = Gstr(i-1) + G(i)*delE
             ENDDO ! end first moment

          gnorm1 = 0.

          DO i = 1, n
             gnorm1 = gnorm1 + dabs(rgnorm*RQUO(i)- rQUAV)
             ENDDO
          gnorm1 = gnorm1/(n*rQUAV)

          IF (iPR.ne.0) GOTO 1005

          IF (iLAU.lt.iMax) GOTO 1001

          ! only for one species else CORRFA AND SEPAR to be implemented
           CONCI = 1.
 1005      DO i = 1, n5
               dQUOFG(i) = 1.
               FPI(i) = G(i)
               ENDDO !end CORRFA and SEPAR
           ! has to be changed if more than one species !!
           DO i = 1, n5
              F(i) = CONCI*FPI(i)
              ENDDO


           IF (iLAU.lt.iMax) GOTO 1001


           ! after iteration finished
           rsDSIG = 0.
           DO i = 1,n5
              rsDSIG = rsDSIG + DSIG(i)
              ENDDO
           rsDSIG = rsDSIG * delE

           DO i = 1, n2
              DSIG(i) = DSIG(i)/rsDSIG
              ENDDO


           ! implementation save routine for density of states
           ! at the moment restiriction to DOS.dat of EQUIT
           ! implementation of EQUIT at this point
           DO i=1,n5
              RQUO(i) = GGES(i) - G(i)
              ENDDO

              Print *, 'Test'
              Print *, 'rQUAV: ',rQUAV
              Print *, 'rgnorm: ',rgnorm

           CALL EQUIT(n5, RQUO, rQUAV/rgnorm, cOMI, cOMA, cOM, E0, EF,
     *                EX, EXE, YEQ, JJ, jin, delE, n, XM0, Fehler)
           ! end of call EQUIT

            ENDDO ! lj


      END ! MUPHDOS


C ----------------------------------------------------------------
      SUBROUTINE PrepSub(jin, E0, PSXM0, EF, n, nPHO, EXE, EXA,
     *                   EX, RQD, RQC, RQMA, RQMI, RFC, cOMI,
     *                   cOMA, cOM, fIMI, fIMA, Fehler)
C ----------------------------------------------------------------
      ! In Anlehnung an MUPHOCOR von W. Reichardt
      ! implementierte Version fuer Paramterberechnung
      ! F. Kargl Mai/Jun 2003

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

C      PARAMETER

      CHARACTER*(*)  Fehler

      DIMENSION     EXA(MC), EXE(MC), EX(MC), EF(MC), RQC(MC),
     *              RQMA(MC), RQMI(MC), RQD(20,MC), PSXM0(MC),RFC(20),
     *              Z(MK)


C be careful with n2 and n5 if modifying closer to original MUPHCOR
      n2 = 2*n
      n5 = n2

      T0eV = rOlfGG(jin, 'T', 'K', Fehler)/11.604 ! Temp in meV
      IF(Fehler.ne.'&ff') RETURN


      DO i = 1,n2
         EX(i) = dexp(dquot0(0.5*PSXM0(i),T0eV))
         EXE(i) = PSXM0(i)*(EX(i)-1./EX(i))
         EXA(i) = EX(i)+1./EX(i)
C         Print *, 'EXE(',i,') : ',EXE(i)
         ENDDO

      fI0  = .5*(fIMA + fIMI)
      dFI  = .5*(fIMA - fIMI)
      cOM  = dcosd(fI0)*dcosd(dFI)
      cOMA = dcosd(fIMA)
      cOMI = dcosd(fIMI)


      DO i = 1, n5
         EF(i) = E0 + PSXM0(i)
         RQC(i) = E0 + EF(i) - 2.* dsqrt0(E0*EF(i))*cOM
         RQMA(i) = E0 + EF(i) - 2.* dsqrt0(E0*EF(i))*cOMA
         RQMI(i) = E0 + EF(i) - 2.* dsqrt0(E0*EF(i))*cOMI

         DO m = 1, nPHO
            m1 = m + 1
            RQD(m,i) = (RQMA(i)**m1 - RQMI(i)**m1)/(2.*m1)
            ENDDO

         ENDDO

      RFC(1) = 1.
      DO m = 1,nPHO+1
         RFC(m+1) = RFC(m)/(m+1)
         ENDDO


      END ! subroutine PrepSub

C ----------------------------------------------------------
      SUBROUTINE ViPho(Gstr, GP, GD, delE, n, nCase, Fehler)
C ----------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'


      CHARACTER*(*)  Fehler

      DIMENSION     Gstr(MC), GP(MC), GD(MC)

      n2 = 2*n
      n1 = n + 1

C      Print *, 'Test viPho', nCase

      IF (nCase.eq.1) THEN
         DO i = 1, n2
            GP(i) = 0.
            no = n2 - i + 1
            DO k = 1, no
               GP(i) = GP(i) + GD(n2 + 1 - k)*GD(k + i - 1)
               ENDDO
            GP(i) = GP(i)*delE
            ENDDO
      ELSE
         DO i = 1, n
            GP(i) = 0.
            DO k = 1, n2
               k1 = n2 + 1 - k
               k2 = i + k - 1 - n

               IF(k2.le.0) k2 = -k2 + 1

               GP(i) = GP(i) + GD(k1)*Gstr(k2)
               ENDDO
            GP(i) = GP(i)*delE
            ENDDO

         DO i = n1,n2
            GP(i) = 0.
            no = 3*n - i + 1
            DO k = 1, no
               GP(i) = GP(i) + GD(n2+1-k)*Gstr(i+k-1-n)
               ENDDO
            GP(i) = GP(i)*delE
            ENDDO

         ENDIF

      END ! subroutine ViPho

C -----------------------------------------------------------------
      SUBROUTINE EQUIT(n5, ESQUO, rQUAV, cOMI, cOMA, cOM, E0,
     *                 ESEF, ESEX, ESEXE, ESYEQ, JJ, jin, delE, n,
     *                 XM0, Fehler)
C -----------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'


      CHARACTER*(*)  Fehler

      DIMENSION     YEQF(MC), YM1(MC), ESQUO(MC), ESEXE(MC),
     *     ESEX(MC), ESEF(MC), ESYEQ(MC), JJ(14), XM0(MC),
     *     DM0(MC)

      Print *, 'rQUAV (EQUI): ', rQUAV

      DO i = 1, n5
         YM1(i) = ESQUO(i)
         ENDDO

      rYnorm = 0.

      DO i = 1,n5
         EN = (cOMI - cOMA)*(E0 + ESEF(i) - 2*dsqrt(E0*ESEF(i))*cOM)
C         Print *, 'EN ',i,': ',EN
         YEQF(i) = ESYEQ(i)*ESEXE(i)*ESEX(i)*E0/(EN*ESEF(i)*ESEF(i))
         YM1(i) = YEQF(i)/rQUAV - YM1(i)
         IF(YM1(i).le.0) THEN
            YM1(i) = 0.
C  uncomment this     Print *, 'this is extremely bad !!',i
            ENDIF
         IF((i.gt.1).and.(i.le.n5/2)) rYnorm = rYnorm + delE*YM1(i)
         ENDDO

C *************** TEST ***********
            Print *, 'EN = ',EN
C            Print *,'i    YEQF(i)    YM1(i)'
C            DO i = 1, n5
C            Print *,i,YEQF(i),YM1(i)
C            ENDDO
            Print *,'bin hier an Partest 10'
            IF (.not.qAsk(' Continue at point 10?')) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
C ******************** END TEST

      DO i = 1, n5
         YM1(i) = YM1(i)/rYnorm
         ENDDO

      ! saving now different guesses and final result in different
      ! files
      K = 1
      Kdummy = 1
      CALL OlfOpen(JJ(2), 1, Kdummy, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCopZ(jin, JJ(2), K, K, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfPutXYD(JJ(2), K, n, XM0, ESYEQ, DM0, Fehler)
      CALL tOlfP(JJ(2), 'tit', 'Input conform MUPHCOR UN(N)', Fehler)
      CALL OlfCnuP(JJ(2), 'y', 'ESYEQ', ' ', Fehler)
      CALL OlfClos(JJ(2), nK, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      Print *, 'saved TOF distribution'
      !next file
      CALL OlfOpen(JJ(3), 1, Kdummy, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCopZ(jin, JJ(3), K, K, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfPutXYD(JJ(3), K, n, XM0, YEQF, DM0, Fehler)
      CALL tOlfP(JJ(3), 'tit', 'Input conform MUPHCOR Z0(N)', Fehler)
      CALL OlfCnuP(JJ(3), 'y', 'YEQF', ' ', Fehler)
      CALL OlfClos(JJ(3), nK, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      Print *, 'saved Input distribution'
      ! next file
      CALL OlfOpen(JJ(4), 1, Kdummy, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCopZ(jin, JJ(4), K, K, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfPutXYD(JJ(4), K, n, XM0, YM1, DM0, Fehler)
      CALL tOlfP(JJ(4), 'tit', 'Input conform MUPHCOR Z1(N)', Fehler)
      CALL OlfCnuP(JJ(4), 'y', 'YM1', ' ', Fehler)
      CALL OlfClos(JJ(4), nK, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      Print *, 'saved density of state calculated'

      END ! subroutine EQUIT
