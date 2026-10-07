C  ====================================================================
C
C      Library  IDA   :  Inelastic data treatment
C      Modul    i41   :     data manipulations / spectra
C
C  ====================================================================

C     Contents :
C        1. OrgSpectra
C              Sum, Cut, Sort, Join, Exch

C  ====================================================================
C  i41 / 1 :   OrgSpectra...
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE OrgSpectraSum (Task, qWeight, nJList, JList,
     *                          qOv, Fehler)
C  --------------------------------------------------------------------
         ! Sum spectra. A very simple task, but it becomes difficult
         ! if the spectra have different x-scales.
            ! JWu corrected 16-22feb91, routine CommonScale 3jun91,
            ! as a separate subroutine 10jul91
            ! de luxe option 19/20nov98 incorporated ..Join 5jun00
         !FK stop at - corrected

      IMPLICIT NONE

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER      Fehler*(*), aus*80, lisK*80, pc*1, un*40,
     *               cWeight*8, Task*(*)
      LOGICAL        qWeight, qOv, qList(MK), qSoJo, qReset, qJeq,
     *               qComN, qComSc, qIndivScale, qSumLarge, qRect
      INTEGER        nJList, JList(*), lj, j, jout, K, KK, nK, Kout,
     *               nCom, nC, nC1, n, nK1,
     *               i, i1, iL, iH, iLsum, iHsum, iZ, iZw, nZ, nZ1, n1,
     *               iShift(MK), mWeight, mWeightDef, mGroup, mGroupDef
      REAL*8         Z(MZ), Z1(MZ), ZK(MK), wi, Wt(MC), zTol, zLow, zHi

      DATA          lisK /' '/, mWeightDef /0/, mGroupDef /0/

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Set-up :
      ! weight
      IF (qWeight) THEN
         mWeight = iAskDMu (
     *   'Weighting: none(0) with error(1) with #scans(2) with z(3)',
     *   mWeightDef, -1, 3)
         IF (mWeight.lt.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         mWeightDef = mWeight
      ELSE
         mWeight = 0
         ENDIF
      IF     (mWeight.eq.0) THEN
         cWeight = ' '
      ELSEIF (mWeight.eq.1) THEN
         cWeight = 'err-' ! per channel / now do nothing
      ELSEIF (mWeight.eq.2) THEN
         cWeight = 'sca-'
      ELSE
         Fehler = 'weight not yet implemented'
         RETURN
         ENDIF

      ! sort ?
      IF (Task.eq.'j')
     * qSoJo = qAskD (' Sort joined spectra', intq(qSoJo))


C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! Ask for grouping :
         IF (nK.eq.1) THEN
            Fehler = 'there is only one spectrum'
            RETURN
            ENDIF
         IF (lj.eq.2) THEN
            qJeq = qAskD ('Same selection for all files', intq(qJeq))
            ENDIF
         IF (lj.eq.1 .or. .not. qJeq) THEN
            CALL Compose2 (aus, ' there are '//cl4(nK), ' spectra')
            Print *, aus
            aus = ' First spectrum of each group (z: by z-value)'
 4222       CALL FrageCD (aus, lisK, lisK)
            IF (lisK.eq.'z') THEN
               mGroup = 1
               zTol = rAskD (' Tolerance for z1',zTol)
            ELSEIF (lisK.eq.'-') THEN
               RETURN !FK inserted 27-05-07
            ELSE
               mGroup = 0
               CALL DecNList (lisK, qList, nK, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (nK.eq.0) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               IF (.not.qList(1)) THEN
                  CALL Gong(2)
                  Print *,
     * ' spectrum 1 must be retained (answer "-" to escape)'
                  GOTO 4222
                  ENDIF
               ENDIF
            ENDIF

         ! Evaluate grouping :
         IF (mGroup.eq.1) THEN
            CALL OlfGet1ZofK (j, 1, nK, ZK, Fehler)
            IF (irSorted(ZK,nK).eq.0) THEN
c               Fehler = 'grouping impossible / not sorted'
c               RETURN
               Print *, ' Spectra not sorted - check result'
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            K = 1
 4322       CONTINUE
            qList(K) = .true.
            zLow = ZK(K)
            zHi = ZK(K)
            DO KK = K+1, nK
               IF (ZK(KK).lt.zLow) zLow = ZK(KK)
               IF (ZK(KK).gt.zHi) zHi = ZK(KK)
               IF (zHi-zLow.le.zTol) THEN
                  qList(KK) = .false.
               ELSE
                  K = KK
                  GOTO 4322
                  ENDIF
               ENDDO
            ENDIF ! mGroup

         ! Documentation :
         IF     (Task.eq.'a') THEN
            IF (cWeight.eq.' ') THEN
               CALL OlfComAdd (jout, 'S', 'sptra summed '//lisK,Fehler)
            ELSE
               CALL OlfComAdd (jout, 'S',
     *              'sptra '//cWeight(1:4)//'summed '//lisK, Fehler)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (Task.eq.'j') THEN
            CALL OlfComAdd (jout, 'j', 'sptra joined '//lisK, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSE
            Fehler = 'PROGRAM ERROR/ no such task in SpectraSum'
            ENDIF

C  Execution / Add :
         IF     (Task.eq.'a') THEN
            ! Common x scale -> X1(1..nC1) :
            CALL CommonScale (j, .true., qComN, nCom, qComSc,
     *                        nC1, X1, iShift, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF     (.not.qComSc .and. .not. qComN) THEN
               Print *, ' there is neither a common scale nor'//
     *                 ' a constant number of channels :'
               Print *, ' this case is at present not allowed'
               RETURN
            ELSEIF (.not.qComSc) THEN
               aus = ' There is no common scale'
               qIndivScale = .true.
            ELSE
               qIndivScale = .false.
               ENDIF

            IF (qIndivScale) THEN
               IF (.not.qAsk(' Proceed using channel numbers ?'))RETURN
            ELSE
               qRect = qComN .and. nCom.eq.nC1
               IF (.not.qRect) THEN
                  Print *, ' the grid is not rectangular'
                  aus = ' Sum only over common(0) or over full(1) area'
                  qSumLarge = qAskDi (aus, intq(qSumLarge))
               ELSE
                  qSumLarge = .true. ! is aber egal
                  ENDIF
               ENDIF

            ! Check weight :
            IF (mWeight.eq.2) THEN
               CALL pcOlfFind (j, 'z', '#scans', un, pc, iZw, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (pc.ne.'z') THEN
                  Fehler = 'inconsistent pc from pcOlfFind'
                  RETURN
               ELSEIF (iZw.le.0) THEN
                  Fehler = 'inconsistent no from pcOlfFind'
                  RETURN
                  ENDIF
               ENDIF

            ! Loop spectra :
            DO K = 1, nK
               IF (qList(K)) THEN ! initialize sum registers :
                  n1 = 0
                  DO i = 1, nC1
                     Wt(i) = 0.
                     Y1(i) = 0.
                     D1(i) = 0.
                     ENDDO
                  IF (qSumLarge) THEN
                     iL = nC1
                     iH =  1 ! wird schon noch groesser werden
                  ELSE
                     iL = 1
                     iH = nC1
                     ENDIF

                  ! Loop over spectra to be summed :
                  DO KK = K, nK
                     IF (KK.gt.K .and. qList(KK)) GOTO 4246 ! next retainable sp.

                     CALL OlfGetSpe (j, KK, nZ, Z, nC, X, Y, D, Fehler)
                     IF (mWeight.le.1) THEN
                        wi = 1. ! default: equal weight for all spectra
                     ELSE
                        wi = Z(iZw)
                        ENDIF
                     IF (qIndivScale) THEN
                        ! Primitivfall : sum range = full spectrum
                        iL    = 1
                        iH    = nC
                        iLsum = 1
                        iHsum = nC
                     ELSE
                        IF (qSumLarge) THEN ! extend range ?
                           iL  = min0 (iL,  1+iShift(KK))
                           iH  = max0 (iH, nC+iShift(KK))
                           iLsum =  1+iShift(KK)
                           iHsum = nC+iShift(KK)
                        ELSE                ! reduce range ?
                           iL  = max0 (iL,  1+iShift(KK))
                           iH  = min0 (iH, nC+iShift(KK))
                           iLsum = iL
                           iHsum = iH
                           ENDIF
                        IF (iH.lt.iL) THEN
                           aus = ' no common x-range, spectrum '//
     * cl4(KK)
                           CALL Append (aus, ' not retained')
                           Print *, aus
                           GOTO 4249
                           ENDIF
                        ENDIF
                     DO i1 = iLsum, iHsum
                        i = i1 - iShift(KK)
                        IF (.not.qIndivScale .and.
     *              .not.qEqEps(X(i),X1(i1))) THEN
                           Print '(a,2g14.6)',
     *   'WARNING/ two points not exactly equal: ', X(i), X1(i1)
                           ENDIF
                        IF (mWeight.eq.1) THEN
                           wi = dquot0(1.d0, D(i)**2)
                           ENDIF
                        Wt (i1) = Wt(i1) + wi
                        Y1 (i1) = Y1(i1) + wi*Y(i)
                        D1 (i1) = D1(i1) +(wi*D(i))**2
                        ENDDO ! i1
                     IF (n1.eq.0) THEN
                        nZ1 = nZ
                        DO iZ = 1, nZ
                           Z1(iZ) = Z(iZ)
                           ENDDO
                     ELSE
                        IF (nZ.ne.nZ1) THEN
                           Fehler = 'nZ varies'
                           RETURN
                           ENDIF
                        DO iZ = 1, nZ
                           Z1(iZ) = Z1(iZ) + Z(iZ)
                           ENDDO
                        ENDIF
                     n1 = n1 + 1
                     ENDDO ! KK
                  ! - mean values :
 4246             CONTINUE
                  DO i = iL, iH
                     IF (Wt(i).le.0) THEN
                        Fehler = 'weight 0'
                        RETURN
                        ENDIF
                     X(i+1-iL) = X1(i)
                     Y(i+1-iL) = Y1(i) / Wt(i)
                     D(i+1-iL) = dsqrt(D1(i)) / Wt(i)
                     ENDDO
                  nC = iH+1-iL
                  DO iZ = 1, nZ
                     IF (iZ.ne.iZw) THEN
                        Z(iZ) = Z1(iZ) / n1
                     ELSE
                        Z(iZ) = rSum (Wt, iL, iH, 1) / nC ! normally Wt(i)=const
                        ENDIF
                     ENDDO
                  Kout = Kout + 1
                  CALL OlfPutSpe (jout, Kout, nZ, Z, nC,
     *                            X, Y, D, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF ! qList(K)
 4249          CONTINUE
               ENDDO ! K

C  Execution / Join :
         ELSEIF (Task.eq.'j') THEN

C  Loop spectra :
            qReset = .true.
            DO K = 1, nK
               IF (qReset) THEN
                  n1 = 0
                  nK1 = 0
                  ENDIF

               CALL OlfGetSpe (j, K, nZ, Z, n, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ! join to stored spectra :
               DO i = 1, n
                  X1(n1+i) = X(i)
                  Y1(n1+i) = Y(i)
                  D1(n1+i) = D(i)
                  ENDDO
               IF (n1.eq.0) THEN
                  nZ1 = nZ
                  DO iZ = 1, nZ
                     Z1(iZ) = Z(iZ)
                     ENDDO
               ELSE
                  IF (nZ.ne.nZ1) THEN
                     Fehler = 'nZ varies'
                     RETURN
                     ENDIF
                  DO iZ = 1, nZ
                     Z1(iZ) = Z1(iZ) + Z(iZ)
                     ENDDO
                  ENDIF
               n1  = n1  + n
               nK1 = nK1 + 1

               IF (n1.gt.MC) THEN
                  Fehler = ' Too many points when joining spectrum '//
     *                     cl4(K)
                  RETURN
                  ENDIF

               ! output, if next spectrum is retained or if nK is reached :
               IF (K.lt.nK) THEN
                  qReset = qList(K+1)
               ELSE
                  qReset = .true.
                  ENDIF

               IF (qReset) THEN
                  DO iZ = 1, nZ
                     Z1(iZ) = Z1(iZ) / nK1
                     ENDDO
                  IF (qSoJo) CALL SortChannels (X1, Y1, D1, n1, .true.)

                  Kout = Kout + 1
                  CALL OlfPutSpe (jout, Kout, nZ, Z1, n1,
     *                            X1, Y1, D1, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               ENDDO ! K

            ENDIF ! Task 'a' or 'j'
C  End Execution

         CALL TensorCheckZ (jout, Fehler)

         CALL OlfClos (jout, Kout, Fehler)

         ENDDO ! lj

      END ! OrgSpectraSum

C  --------------------------------------------------------------------
      SUBROUTINE OrgSpectraCut (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! cut spectra.
            ! JWu 1990/91, separate subroutine 10jul91

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER    (MRange=10)

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      DIMENSION      qList(MK)
      CHARACTER*80   aus, lisK, txtRange
      REAL*8         RRange(MRange), Z(MK)

      DATA          lisK /'1'/, qPerJin / .false. /, iMod / 1 /

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (nK.eq.1) THEN
            Fehler = 'there is only one spectrum'
            RETURN
            ENDIF

         IF     (lj.eq.1) THEN
            qSelAsk = .true.
            nK1 = nK
         ELSEIF (lj.eq.2) THEN
            IF (nK.eq.nK1) THEN
               qPerJin = .not. qAskD (
     *            ' Same selection for all files', 1-intq(qPerJin))
               qSelAsk = qPerJin
            ELSE
               qSelAsk = .true.
               ENDIF
         ELSEIF (lj.gt.2 .and. .not.qSelAsk .and. nK.ne.nK1) THEN
            CALL Gong (7)
            Print *,
     * ' cannot use same selection : different number of spectra'
            qSelAsk = .true.
            ENDIF
         IF (qSelAsk) THEN
            CALL Say2 (' there are '//cl4(nK), ' spectra')
            iModIn = iAskDMu (' Select by no.(1) by value(2)',
     *                        iMod, 0, 2)
            IF (iModIn.le.0) THEN
               Fehler = ' '
               RETURN
               ENDIF
            iMod = iModIn
            IF     (iMod.eq.1) THEN
               aus = ' Delete which spectra'
               CALL GetNList (aus, lisK, qList, nK)
               IF (lisK.eq.'-') THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
            ELSEIF (iMod.eq.2) THEN
               Print *,
     * ' currently only one option: select from z-range'
               iZ = iAskDMu (' Select according to which z', iZ, 0, MZ)
               IF (iZ.le.0) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               CALL rAskOnOff (' Retain range',
     *            txtRange, RRange, MRange, nRange)
               lisK = txtRange ! doc
            ELSE
               Fehler = 'PROG ERR/ iMod oor'
               RETURN
               ENDIF ! iMod
            ENDIF ! qSelAsk

C  Execute selection :
         IF (iMod.eq.2) THEN
            CALL OlfGet1ZofK (j, iZ, nK, Z, Fehler)
            DO K = 1, nK
               qList(K) = .true. ! a priori, through away
               DO ir = 1, nRange, 2
                   IF (qrinside(Z(K), RRange(ir), RRange(ir+1)))
     *                 qList(K) = .false.
                  ENDDO
               ENDDO
            ENDIF

C  Documentation :
         CALL OlfComAdd (jout, 's', 'deleted sptra '//lisK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK
            IF (.not.qList(K)) THEN
               Kout = Kout + 1
               CALL OlfCopSpe (j, jout, K, Kout, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDIF ! qList(K)
            ENDDO ! K

         CALL TensorCheckZ (jout, Fehler)

         CALL OlfClos (jout, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  End loop files :
         ENDDO

      END ! OrgSpectraCut

C  --------------------------------------------------------------------
      SUBROUTINE OrgSpectraSort (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 18feb92, 9mrz98

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f' ! using ZZofK
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*), KRang(MK)
      DATA           qLR / .true. /

C  Files :
      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      qLRset = .false.

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (nK.le.1) THEN
            Fehler = 'only one spectrum - nothing to sort'
            RETURN
            ENDIF

C  Range :
         nZ = iOlfG (j, '#Z', Fehler)
         IF (nZ.gt.1 .and. .not.qLRset) THEN
            qLR = qAskD ('Sort from left to right', intq(qLR))
            qLRset = .true.
            ENDIF
         DO iZ = 1, nZ
            IF (qLR) THEN
               CALL OlfGet1ZofK (j, iZ, nK, ZZofK(1,iZ), Fehler)
            ELSE
               CALL OlfGet1ZofK (j, iZ, nK, ZZofK(1,nZ+1-iZ), Fehler)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            ENDDO
         ifail = -1
         CALL M01DEF_local (ZZofK, MK, 1, nK, 1, nZ, 'a', KRang, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
         IF (ifail.ne.0) THEN
            Fehler = 'OrgSpectraSort/ M01DEF error'
            RETURN
            ENDIF
         ifail = 0
         CALL M01ZAF_local (KRang, 1, nK, ifail) ! invert the permutation !Artem: Replace with a self-made subroutine from lnag_local.f
         IF (ifail.ne.0) THEN
            Fehler = 'OrgSpectraSort/ M01ZAF error'
            RETURN
            ENDIF

C  Sort :
         DO KK = 1, nK
            Kout = Kout + 1
            CALL OlfCopSpe (j, jout, KRang(KK), Kout, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         CALL OlfClos (jout, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Move back :
         IF (qOv) THEN
            DO K = 1, nK
               CALL OlfCopSpe (jout, j, K, K, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO

            CALL MemFileDel (jout, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDIF

         ENDDO

      END ! OrgSpectraSort

C  --------------------------------------------------------------------
      SUBROUTINE OrgSpectraExch (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! exchange channels <-> spectra (x <-> z)
         ! JWu ID_SUM jan91, integrated feb91, new CommonScale 4jun91

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*)  Fehler
      INTEGER        JList(*)
      DIMENSION      iShift(MK), Z(MZ), Z1(MZ)
      CHARACTER      co*40, un*40, cZ*3

      DATA           iZ / 1 /

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

c      iZ = iAskDMu ('Exchange x with which z', iZ, 0, MZ)
c      IF (iZ.le.0) THEN
c         Fehler = ' '
c         RETURN
c         ENDIF
      iZ = 1
      cZ = 'z'//cl2(iZ)

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         nZ = iOlfG (j, '#Z', Fehler)
         IF (nZ.gt.1) THEN
            Fehler=' x <-> z makes no sense if there are several z''s'
            RETURN
            ENDIF

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (j,    'x', co, un, Fehler)
         CALL OlfCnuP (jout, cZ,  co, un, Fehler)
         CALL OlfCnuG (j,    cZ,  co, un, Fehler)
         CALL OlfCnuP (jout, 'x', co, un, Fehler)

         CALL OlfComAdd (jout, 'x', 'exchanged x<->'//cZ, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Common x - scale of input spectra :
         CALL CommonScale (j, .true., qComN, nCom, qCom2,
     *                     nC2, X2, iShift, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (nC2.gt.MK) THEN
            Fehler = ' '//cl4(nC2)
            CALL Append (Fehler, ' channels, but only '//cl4(MK))
            CALL Append (Fehler, ' spectra are allowed')
            RETURN
            ENDIF

         IF (.not.qCom2) THEN
            Print *, ' there is no common x-scale'
            IF (.not.qAsk(' Proceed using channel numbers ?')) RETURN
            IF (.not.qComN) THEN
               Print *,
     * ' channel number is not constant - case not allowed'
               RETURN
            ELSE
               nC2 = nCom
               ENDIF
            ENDIF

C  Proceed to create new spectra (X1,..) :
         DO i1 = 1, nC2 ! new spectrum no. i1 = Kout
            Z1(iZ) = X2(i1)
            nC1 = 0
            DO K = 1, nK
               CALL OlfGetSpe (j, K, nZ, Z, nC, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (qCom2) THEN
                  ! take channel with x=X2(i1) (= old x = new z)
                  i0 = i1 - iShift(K)
                  IF (i0.lt.1 .or. i0.gt.nC) GOTO 29 ! X2(i1) not in X
                  ! now the channel should be found, x=X(i0)
                  IF (.not.qEqEps(X(i0), z1)) THEN ! error occured 7dec95
                     Fehler =
     * 'There might be a slight difference in z values'
                     RETURN
                     ENDIF
               ELSE
                  ! take channel with number i1
                  i0 = i1
                  ENDIF
               ! new element for new spectrum
               nC1     = nC1 + 1
               X1(nC1) = Z(iZ)
               Y1(nC1) = Y(i0)
               D1(nC1) = D(i0)
 29            CONTINUE
               ENDDO
            IF (nC1.le.0) THEN
               Print *, ' i1, nC1 ', i1, nC1
               CALL Absturz ('IDA/mx', 'created empty spectrum')
               ENDIF
            IF (nC1.ge.MC)
     * CALL Absturz ('IDA/mx', 'created too long spectrum')

            CALL SortChannels (X1, Y1, D1, nC1, .true.)

            Kout = Kout + 1
            CALL OlfPutSpe (jout, Kout, nZ, Z1, nC1, X1,Y1,D1,Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO ! i1

         CALL OlfClos (jout, Kout, Fehler)

C  End loop files :
         ENDDO ! lj

      END ! OrgSpectraExch

C  --------------------------------------------------------------------
      SUBROUTINE OrgSpectraBreak (nJList, JList, Fehler)
C  --------------------------------------------------------------------
         ! Break one file into several files.
         ! JWu 28sep00 (sic) mainly for separating VV and VH, 90 and 172 deg

      IMPLICIT NONE

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER      Fehler*(*), aus*80, lisK*80, UnZ*40, CoZ*40,
     *               fil*40, tit*80, h*20, CDef(MK)*20
      LOGICAL        qList(MK), qJeq, qDef
      INTEGER        nJList, JList(*), lj, j, jout, jouti, K, KK, nK,
     *               Kout, nC, n, nK1, mGroup, mGroupDef, ih,
     *               nZ, iDef, nDef
      REAL*8         ZK(MK), Z(MZ), zout, zHi, zLow, zTol, ZDef(MK)

      DATA          lisK /' '/, mGroupDef /0/

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop files :
      nDef = 0
      DO lj = 1, nJList
         j = JList(lj)

         CALL tOlfG (j, 'fil', fil, Fehler)
         CALL tOlfG (j, 'tit', tit, Fehler)
         nK = iOlfG (j, '#spectra', Fehler)

         ! Ask for grouping :
         IF (nK.eq.1) THEN
            Fehler = 'there is only one spectrum'
            RETURN
            ENDIF
         IF (lj.eq.2) THEN
            qJeq = qAskD ('Same selection for all files', intq(qJeq))
            ENDIF
         IF (lj.eq.1 .or. .not.qJeq) THEN
            CALL Compose2 (aus, ' there are '//cl4(nK), ' spectra')
            Print *, aus
            aus = ' First spectrum of each group (z: by z-value)'
 4222       CALL FrageCD (aus, lisK, lisK)
            IF (lisK.eq.'z') THEN
               mGroup = 1
               zTol = rAskD (' Tolerance for z1',zTol)
            ELSE
               mGroup = 0
               CALL DecNList (lisK, qList, nK, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (nK.eq.0) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               IF (.not.qList(1)) THEN
                  CALL Gong(2)
                  Print *,
     * 'spectrum 1 must be retained (answer "-" to escape)'
                  GOTO 4222
                  ENDIF
               ENDIF
            ENDIF

         ! Evaluate grouping :
         IF (mGroup.eq.1) THEN
            CALL OlfGet1ZofK (j, 1, nK, ZK, Fehler)
            IF (irSorted(ZK,nK).eq.0) THEN
C               Fehler = 'grouping impossible / not sorted'
C               RETURN
               Print *, ' Spectra not sorted - check result'
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            K = 1
 4322       CONTINUE
            qList(K) = .true.
            zLow = ZK(K)
            zHi = ZK(K)
            DO KK = K+1, nK
               IF (ZK(KK).lt.zLow) zLow = ZK(KK)
               IF (ZK(KK).gt.zHi) zHi = ZK(KK)
               IF (zHi-zLow.le.zTol) THEN
                  qList(KK) = .false.
               ELSE
                  K = KK
                  GOTO 4322
                  ENDIF
               ENDDO
            ENDIF ! mGroup

         ! Documentation :
         CALL OlfCnuG (j, 'z1', CoZ, UnZ, Fehler)

         jouti = 0

         ! Loop spectra :
         DO K = 1, nK
            IF (qList(K)) THEN ! open new file
               CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (jouti.eq.0) jouti = jout

               CALL OlfGet1Z (j, K, 1, zout, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               CALL NiceNum (zout, h, ih)
               CALL Say3 ('file starting at spectrum '//cl4(K),
     *                    ', '//CoZ, ' = '//h)
               CALL tOlfP (jout, 'tit',
     * tit(1:lenU(tit))//' '//CoZ(1:lenU(CoZ))//'='//h, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               qDef = .false.
               DO iDef = 1, nDef
                  IF (zout.eq.ZDef(iDef)) THEN
                     h=CDef(iDef)
                     qDef = .true.
                     GOTO 29
                     ENDIF
                  ENDDO
 29            CONTINUE

               CALL FrageCD ('Add to file name', h, h)
               CALL tOlfP (jout, 'fil', fil(1:lenU(fil))//h, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               IF (.not.qDef .and. nDef.lt.MK) THEN
                  nDef = nDef + 1
                  ZDef(nDef) = zout
                  CDef(nDef) = h
                  ENDIF

               ENDIF

            CALL OlfGetSpe (j, K, nZ, Z, nC, X, Y, D, Fehler)
            Kout = Kout + 1
            CALL OlfPutSpe (jout, Kout, nZ, Z, nC, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (K.eq.nK .or. qList(K+1)) THEN
               CALL TensorCheckZ (jout, Fehler)
               CALL OlfClos (jout, Kout, Fehler)
               ENDIF

            ENDDO ! K

         ENDDO ! lj

      END ! OrgSpectraBreak
