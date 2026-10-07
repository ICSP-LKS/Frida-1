C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i40   :     data manipulations / channels
C
C  ====================================================================

C     Contents :
C        1. OrgCh
C              Sum, SAuto, Cut, Sort, Group, Subs
C        2. Grid
C              GridCurve, GridIntExt, GridCutSum, GridRedis, GridHist,
C              SetGrid, SetGridReg
C        3. OrgHist
C              Make

C     Aenderungsverzeichnis :
C     JWu   nov91 : reorganisation, Ida4(200) restricted to Org..
C     JWu   mai91 : Modulaufteilung
C     JWu   feb91 : Verbessert, Einbau in IDA
C     JWu   jan91 : Hauptbestandteile

C  ====================================================================
C  i40 / 1 :   OrgCh...
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE OrgChSum (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! sum channels
         ! sauber 21feb91, revu 3jun91, corr. 19jun91, separately 4jan93

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)
      DIMENSION      qList(MC), qList2(MC), iShift(MK)

      CHARACTER*80  aus, LisDoc
      CHARACTER*320 LisC, LisCofK

      DATA          LisC /' '/, LisCofK /' '/, qAllJ /.true./

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Common x scale :
         CALL CommonScale (j, .true., qComN, nCom, qCom2,
     *                     nC2, X2, iShift, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Compare subsequent files (21jun93) :
         IF (lj.eq.1) THEN
            qComN1 = qComN
            nCom1  = nCom
            qCom21 = qCom2
            nC21   = nC2
         ELSE
            IF (qComN.eqv.qComN1 .and. nCom1.eq.nCom .and.
     *          qCom21.eqv.qCom2 .and. nC21.eq.nC2) THEN
               IF (lj.eq.2) THEN
                  qAllJ = qAskD (' Same selection for all files',
     *                           intq(qAllJ))
                  qSamJ = qAllJ
                  ENDIF
            ELSE
               qSamJ = .false.
               ENDIF
            IF (qSamJ) GOTO 5 ! skip the questionary
            ENDIF

C  Select the channels once and for all ?
         IF     (nK.eq.1) THEN
            Print *, ' there is just one spectrum'
            qAskPerK = .false.
         ELSEIF (.not.qComN .and. .not.qCom2) THEN
            IF (.not.qAskD (' Select channels per spectrum', 1)) THEN
               Fehler = ' '
               RETURN
               ENDIF
            qAskPerK = .true.
         ELSE
            aus = ' The same selection for all spectra'
            qAskPerKin = .not. qAskD (aus, 1-intq(qAskPerKin))
            qAskPerK   = qAskPerKin
            ENDIF

C  When selecting them, refer to channel number or to energy ?
         IF (.not.qAskPerK) THEN
            IF     (qComN .and. qCom2 .and. nCom.ne.nC2) THEN
               aus=' Refer to common scale(0) or to channel number(1)'
               qISin = qintr (qAskDi (aus, intq(qISin)) )
               qIndivScale = qISin
            ELSEIF (qComN .and. qCom2) THEN
               ! in this case, both options have the same effect
               qIndivScale = .false.
            ELSEIF (qComN) THEN
               qIndivScale = .true.
            ELSE
               qIndivScale = .false.
               ENDIF

C  Select :
 41         CONTINUE
            IF (qIndivScale) THEN
               aus = ' there are '//cl4(nCom)
               CALL Append (aus, ' channels per spectrum')
               Print *, aus
               CALL FrageNList (' Retain which channels ', LisC, nC2) ! nCom ??
               CALL Minuskeln (LisC)
               IF     (LisC.eq.'-') THEN
                  Fehler = ' '
                  RETURN
               ELSE
                  CALL DecNList (LisC, qList2, nC2, Fehler ) ! nCom ??
                  ENDIF
            ELSE
               aus = ' the common x-scale contains '//cl4(nC2)
               CALL Append (aus, ' points')
               Print *, aus
               aus = ' Retain which channels'
               CALL GetNList (aus, LisC, qList2, nC2)
               IF (LisC.eq.'-') THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               ENDIF

            IF (.not.qList2(1)) THEN
               CALL Gong(2)
               Print *, ' channel 1 must be retained'
               GOTO 41
               ENDIF
            LisDoc = LisC

         ELSE
            LisDoc = 'as given per spectrum'
            ENDIF

C  Special case - sum over subrange :
         qSumSubrange = .not.(qCom2 .and. nCom.eq.nC2)
     *                  .and. .not.qAskPerK .and. .not.qIndivScale
               ! It will be summed using the common scale, and the
               ! individual scales are subranges of the common scale.
               ! Therefore it may occur that for a given group not all
               ! channels are contained in the individual spectrum's
               ! x-range.
               ! reformulated 19jun91
         IF (qSumSubrange) THEN
            aus=' Retain first/last group (dangerous for common grid)'
            qRetain = qAskD (aus, intq(qRetain))
            ENDIF

         CALL Insert (LisDoc, 1, 'chs summed ')

 5       CONTINUE
C  End questionary (jumped if qSamJ).

C  Save :
         CALL OlfComAdd (jout, 'C', LisDoc, Fehler)

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, nC, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (qAskPerK) THEN
C  Ask for spectra to be retained :
               CALL Compose3 (aus,  ' spectrum '//cl4(K),
     *             ' contains '//cl4(nC), ' points')
               Print *, aus
 4132          aus = ' Retain which channels'
               CALL GetNList (aus, LisCofK, qList, nC)
               IF     (LisCofK.eq.'-') THEN
                  Fehler = ' '
                  RETURN
               ELSEIF (.not.qList(1)) THEN
                  CALL Gong(2)
                  Print *, ' channel 1 must be retained'
                  GOTO 4132
                  ENDIF
            ELSE
C  Transfer common -> individual scale :
               DO i = 1, nC
                  IF (qIndivScale) THEN
                     qList(i) = qList2(i)
                  ELSE
                     IF (i+iShift(K).gt.nC2) CALL Absturz ('IDA/mc',
     *                  'major inconsistency with iShift, nCc')
                     IF (.not.qEqTol(X(i),X2(i+iShift(K)), 1.d-4)) THEN
                        ! note : CommonScale has tol = 1d-5
                        Print *, ' j K i iSh(K) = ', j, K, i, iShift(K)
                        Print *, ' X(1) X(i) X(n) =', X(1), X(i), X(nC)
                        Print *, ' X2(1+iSh), X2(i+iSh) =',
     *                       X2(1+iShift(K)), X2(i+iShift(K))
                        CALL Absturz ('IDA/mc', 'problems with iShift')
                        ENDIF
                     qList(i) = qList2(i+iShift(K))
                     ENDIF
                  ENDDO
               IF (qSumSubrange .and. qRetain) qList(1) = .true.
                  ! add a group beginning with the first channel
                  ! otherwise the group may be lost
               ENDIF

C  Select / sum :
            nC1 = 0
            DO i = 1, nC
               IF (qList(i)) THEN
                  nC1  = nC1 + 1
                  xneu = X(i)
                  yneu = Y(i)
                  dneu = D(i)**2
                  nsum = 1
                  ! sum up to the next retained channel
                  DO ii = i+1, nC
                     IF (qList(ii)) GOTO 415
                     nsum = nsum + 1
                     xneu = xneu + X(ii)
                     yneu = yneu + Y(ii)
                     dneu = dneu + D(ii)**2
                     ENDDO
 415              CONTINUE
                  X1(nC1) = xneu / nsum
                  Y1(nC1) = yneu / nsum
                  D1(nC1) = dsqrt(dneu) / nsum
                  ENDIF
               ENDDO

            IF (qSumSubrange .and. .not. qRetain) THEN
C  suppress last channel ..
               nC1 = nC1 - 1
C  .. except if it has the correct energy :
               IF (nC+1+iShift(K).le.nC2) THEN
                  IF (qList2(nC+1+iShift(K))) nC1 = nC1 + 1
                  ENDIF
               ENDIF

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, nC1, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (jout, nK, Fehler)

C  End loop files :
         ENDDO

      END ! OrgChSum

C  --------------------------------------------------------------------
      SUBROUTINE OrgChSAuto (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! sum channels such as to restrict the relative error
         ! 2mar93

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)
      CHARACTER      doc*40, h1*20, h2*20

      DATA          rErrL / 1.d-3 /, rMulX /1.d0/,
     *              rAddX /0.d0/, maxGr / 7 /

      IF (nJList.le.0) RETURN

      rErrL = rAskDMu (' Limit for dy/y', rErrL, 0.d0, 1.d0)
      rMulX = rAskDMu (' Maximum multiplicator for x-steps (1=off)',
     *                     rMulX, 0.d0, 1.d3)
      rAddX = rAskDMu (' Maximum increment for x-steps (0=off)',
     *                     rAddX, 0.d0, 1.d3)
      maxGr = iAskDMu (' Maximal # channels / group', maxGr, 1, MC)

      CALL NiceNum (rErrL, h1, ih1)
      CALL NiceNum (rMulX, h2, ih2)
      doc = 'add chs/ dy/y < '//h1(1:ih1)//'; mulX < '//h2(1:ih2)//
     *      '; maxgr='//cl4(maxGr)

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, 'C', doc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            IF (irSorted(X,n).ne.2) THEN
               Fehler = 'Spectrum not in ascending order'
               RETURN
               ENDIF

            IF (rMulX.gt.1 .and. X(1).le.0) THEN
               Fehler = 'x multiplicator incompatible with X()<0'
               RETURN
               ENDIF

            i  = 1
            n1 = 0
 1             CONTINUE
               ! new group
               n1   = n1 + 1
               xneu = X(i)
               xneu1 = xneu
               yneu = Y(i)
               dneu = D(i)**2
               nsum = 1
               ! sum up to the next retained channel
               DO ii = i+1, n
                  IF ((rErrL.gt.0 .and. dneu.le.(rErrL*yneu)**2)
     *                 .or. nsum.ge.maxGr
     *                 .or. (rMulX.gt.1 .and. X(ii).gt.xneu1*rMulX)
     *                 .or. (rAddX.gt.0 .and. X(ii).gt.xneu1+rAddX))
     *                 GOTO 4
                  nsum = nsum + 1
                  xneu = xneu + X(ii)
                  yneu = yneu + Y(ii)
                  dneu = dneu + D(ii)**2
                  ENDDO
               ! save new group
 4             CONTINUE
               i = ii
               X1(n1) = xneu / nsum
               Y1(n1) = yneu / nsum
               D1(n1) = dsqrt(dneu) / nsum
               IF (i.le.n) GOTO 1

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n1, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (jout, nK, Fehler)

C  End loop files :
         ENDDO

      END ! OrgChSAuto

C  --------------------------------------------------------------------
      SUBROUTINE OrgChCut (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! select channels
            ! revu 21feb91, 3jun91, corr. 19jun91,
            ! qVal 25jun92, from file 4jan93
            ! SelectChSet/Get -> i3: 11jul94.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*)  Fehler
      CHARACTER      LisDoc*40
      INTEGER        JList(*)
      LOGICAL        qList(MC)

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL SelectChSet ('Retain', j, lj, LisDoc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, 'c', 'chs kept '//LisDoc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         Kout = 0
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL SelectChGet (K, n, X, Y, D, qList, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            nn = 0
            DO i = 1, n
               IF (qList(i)) THEN
                  nn    = nn + 1
                  X(nn) = X(i)
                  Y(nn) = Y(i)
                  D(nn) = D(i)
                  ENDIF
               ENDDO
            IF (nn.ge.1) THEN
               Kout = Kout + 1
               CALL OlfCopZ (j, jout, K, Kout, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL OlfPutXYD (jout, Kout, nn, X, Y, D, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDIF

            ENDDO ! K

         CALL OlfClos (jout, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO ! lj

      END ! OrgChCut

C  --------------------------------------------------------------------
      SUBROUTINE OrgChSort (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! sort the channels
         ! JWu 11mar91

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      DATA          qAver /.true./

      qAver = qAskD (' Average channels with equal x', intq(qAver))

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            CALL SortChannels (X, Y, D, n, qAver)
            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         CALL OlfClos (jout, nK, Fehler)

         ENDDO
C  End loop files.

      END ! OrgChSort

C  --------------------------------------------------------------------
      SUBROUTINE OrgChGroup (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 24jun92
         ! determine optimal grouping of channels
         ! (such that the relative error remains below a limit rGmax)

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      DIMENSION     qLisX(MC)
      CHARACTER*20  h1

      DATA          rGmax /.1/, modG /1/

      IF (nJList.le.0) THEN
         Print *, '   ida> ? determine optimal grouping of channels'
         RETURN
         ENDIF

      rGmax = rAskDLu (' Maximum relative error of groups',
     *                 rGmax, 0.d0, 10.d0)
      modG  = iAskDMu (' Save no.-of-group (1), new/old (2)',
     *                 modG, 1, 2)

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL NiceNum (rGmax, h1, ih1)
         CALL OlfComAdd (jout, '#', 'groups with rel.err. < '//
     *                   h1, Fehler)
         CALL OlfCnuP   (jout, 'y', '# group', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)

C  Determine groups (for the algorithm see C2,55) :
            CALL qSet (qLisX, 1, n, 1, .false.)
            i = 1
            qLisX(1) = .true.
 101        CONTINUE
               yG = Y(i)
               dG = D(i)**2
               DO ii = i, n
                  IF (dG.lt.yG**2 * rGmax**2) GOTO 109
                  yG = yG + Y(ii)
                  dG = dG + D(ii)**2
                  ENDDO
               GOTO 119 ! no end of last group
 109           CONTINUE ! end of group
               i = ii + 1 ! start of next group
               IF (i.le.n) THEN
                  qLisX(i) = .true.
                  GOTO 101
                  ENDIF
 119        CONTINUE

C  Encode result :
            IF     (modG.eq.1) THEN
               yi = 0.
               DO ii = 1, n
                  IF (qLisX(ii)) yi = yi + 1.
                  Y(ii) = yi
                  D(ii) = 0.
                  ENDDO
            ELSEIF (modG.eq.2) THEN
               DO ii = 1, n
                  Y(ii) = dble(intq(qLisX(ii)))
                  D(ii) = 0.
                  ENDDO
               ENDIF

C  Save result :
            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! K

         CALL OlfClos (jout, nK, Fehler)

         ENDDO
C  End loop files.

      END ! OrgChGroup

C  --------------------------------------------------------------------
      SUBROUTINE OrgChSubs (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! substitute y by f(y)
         ! AMeyer Sept94; Error bars JWu Jul95

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      CHARACTER      text*40
      INTEGER        JList(*), KK(MK)

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  zweiten File erfragen und oeffnen, Anzahl der Spektren kontrollieren
         j1 = iAsk ('f(y) from File ?')
         nK1 = iOlfG (j1, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL GetK2K (j, j1, nK1, KK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Documentation :
         CALL tOlfG (j, 'tit', text, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfComAdd (jout, '_', 'y -> f(y) according to '//
     *                   text, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            CALL OlfGetXYD (j1, KK(K), n1, X1, Y1, D1, Fehler)

            DO i = 1 ,n

C  Linken und rechten Nachbarn interpolieren
C  'l' : X(irPos-1) < val <= X(irPos)

               i1 = irPos (X1, n1, Y(i), 'l')
               IF ((i1.lt.2) .or. (i1.gt.n1)) THEN
                  Fehler = ' too less values for f(y) '
                  RETURN
                  ENDIF

               CALL LinIntPol (Y(i), X1(i1-1), X1(i1), Y1(i1-1),
     *                         Y1(i1),D1(i1-1), D1(i1), yneu, dneu)
               CALL LinIntPol (Y(i)-D(i), X1(i1-1), X1(i1), Y1(i1-1),
     *                         Y1(i1),D1(i1-1), D1(i1), yneum, dneum)
               CALL LinIntPol (Y(i)+D(i), X1(i1-1), X1(i1), Y1(i1-1),
     *                         Y1(i1), D1(i1-1), D1(i1), yneup, dneup)
               Y(i) = yneu
               D(i) = ( dsqrt((yneum-yneu)**2 + dneum**2) +
     *                  dsqrt((yneup-yneu)**2 + dneup**2) ) / 2
               ENDDO
            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         CALL OlfClos (jout, nK, Fehler)

         ENDDO
C  End loop files.

      END ! OrgChSubs

C  --------------------------------------------------------------------
      SUBROUTINE OrgChExch (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! exchange x <-> y (JWu 25jul95)

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      CHARACTER      LisDoc*40
      CHARACTER*40   un, co
      DIMENSION      qList(MC)

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, 'x', 'x <-> y', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (j,    'x', un, co, Fehler)
         CALL OlfCnuP (jout, 'y', un, co, Fehler)
         CALL OlfCnuG (j,    'y', un, co, Fehler)
         CALL OlfCnuP (jout, 'x', un, co, Fehler)

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            DO i = 1, n
               a    = Y(i)
               Y(i) = X(i)
               X(i) = a
               D(i) = 0
               ENDDO

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (jout, nK, Fehler)

C  End loop files :
         ENDDO

      END ! OrgChExch

C  --------------------------------------------------------------------
      SUBROUTINE OrgChSpectra (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! break each spectrum into several spectra
            ! 23may99 for FPI histogram analysis

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*)  Fehler
      CHARACTER      LisDoc*40
      INTEGER        nJList, JList(*), lj, j, jout, K, Kout,
     *               nK, nZ, i, n, nn
      LOGICAL        qOv, qList(MC+1)
      REAL*8         Z(MZ)

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL SelectChSet ('1st channel of groups', j, lj,
     *                     LisDoc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         nZ = iOlfG (j, '#Z', Fehler)
         CALL OlfCnuP (jout, 'z'//ch1(nZ+1), 'no_spe', ' ', Fehler)
         CALL OlfComAdd (jout, 't', 'broken / channels '//
     *                   LisDoc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         Kout = 0
         DO K = 1, nK
            CALL OlfGetSpe (j, K, nZ, Z, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL SelectChGet (K, n, X, Y, D, qList, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (.not.qList(1)) THEN
               Fehler = 'channel 1 must start new group'
               RETURN
               ENDIF

            DO i = 1, n
               IF (qList(i)) THEN
                  nn = 0
                  ENDIF
               nn    = nn + 1
               X(nn) = X(i)
               Y(nn) = Y(i)
               D(nn) = D(i)
               IF (i.eq.n .or. qList(i+1)) THEN ! group completed
                  Kout = Kout + 1
                  Z(nZ+1) = Kout
                  CALL OlfPutSpe (jout, Kout, nZ+1, Z, nn,
     *                            X, Y, D, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               ENDDO ! i -> nn

            ENDDO ! K

         CALL OlfClos (jout, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO ! lj

      END ! OrgChSpectra

C  ====================================================================
C  i40 / 2 :   OrgGrid
C  ====================================================================

C  --------------------------------------------------------------------
      SUBROUTINE GridCurve (nJList, JList, Fehler)
C  --------------------------------------------------------------------
            ! JWu 1aug91, separately 19nov96

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler
      CHARACTER*80  aus, doc, h2
      CHARACTER*40  Co, Un, Co2, Un2
      INTEGER       nJList, JList(*), lj, j, jcc, jpar, jout,
     *              KccK(MK), KppK(MK), K, Kout, nK, nKdummy,
     *              n, n1, nP, i, ii, nZ, nZin, iFuNo, jCuConvAsk, iu
      REAL*8        Z(MZ)
      LOGICAL       qCurve, qConvo

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Prepare grid (outsourced) :
      CALL SetGridChoice (.false., doc, Fehler)

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)

         ! input and checks :
         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (.not.qCurve) THEN
            Fehler = ' File is no curve'
            RETURN
            ENDIF

         ! prepare evaluation and documentation :
         iFuNo  = iOlfG (j, 'fu#', Fehler)
         qConvo = qOlfGdef (j, '?conv', 0, Fehler)
         jcc = jCuConvAsk (qCurve, qConvo)
         IF (jcc.ne.0) THEN
            CALL GetK2K (j, jcc, nKdummy, KccK, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            h2 = 'convoluted/'
         ELSE
            h2 = 'evaluated/'
            ENDIF
         jpar = iOlfG (j, 'fit-par-file#', Fehler)
         IF (jpar.ne.0) THEN
            jpar = iAskDMu (' Parameter file', jpar, 1, MF)
            CALL GetK2K (j, jpar, nKdummy, KppK, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDIF

         ! output file is no curve :
         CALL iOlfP (jout, '?cu', 0, Fehler)
         CALL iOlfP (jout, 'plot-sy#', 0, Fehler)

         ! fit parameters become z :
         nP  = iOlfG (j, '#fit-par', Fehler)
         nZ  = iOlfG (j, '#Z', Fehler)
         IF (nZ+nP.gt.MZ) THEN
            Fehler = 'too many parameters (z+p)'
            RETURN
            ENDIF
         IF (Fehler.ne.'&ff') RETURN
         DO i = 1, nP
            CALL OlfCnuG (j, 'p'//cl2(i), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL OlfCnuP (jout, 'p'//cl2(i), ' ', ' ', Fehler) ! delete p*
            IF (Fehler.ne.'&ff') RETURN
            ! maybe the parameter name was also used for z :
            DO ii = 1, nZ
               CALL OlfCnuG (j, 'z'//cl2(ii), Co2, Un2, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (Co.eq.Co2) THEN
                  CALL Append (Co, '[fit]')
                  GOTO 719
                  ENDIF
               ENDDO
 719        CONTINUE
            CALL OlfCnuP (jout, 'z'//cl2(nZ+i), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         ! documentation :
         CALL Compose2 (aus, h2(1:lenU(h2)), ' '//doc)
         CALL OlfComAdd (jout, '=', aus, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! - file name :
         CALL tOlfG (j, 'fil', aus, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         iu = lenU (aus)
         IF (aus(iu:iu).eq.'-') THEN
            aus(iu:iu) = '=' ! overwrite: `-' becomes `='
         ELSE
            CALL Append (aus, '-')
            ENDIF
         CALL tOlfP (jout, 'fil', aus, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! - file title :
         CALL tOlfG (j, 'tit', aus, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL Insert (aus, 1, 'evaluated ')
         CALL tOlfP (jout, 'tit', aus, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! x-scale per file :
         CALL SetGridJ (j, Fehler)
            IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK

            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (n.ne.nP) THEN
               Fehler = 'INCONSISTENCY/ n<>nP'
               RETURN
               ENDIF

            CALL SetGridK (K, n1, Fehler) ! x grid per spectrum -> X1(1..n1)
            IF (Fehler.ne.'&ff') RETURN

            CALL CuFPFPrep (jpar, KppK, Fehler)
            CALL CuConvPrep (X1, n1, jcc, KccK(K), Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL CuConvVal (iFuNo, Y, Y1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL rSet (D1, 1, n1, 1, 0.d0)

            CALL OlfGetZ (j, K, nZin, Z, Fehler)
            IF (nZin.ne.nZ) THEN
               Fehler = 'INCONSISTENCY/ nZ<>nZ'
               RETURN
               ENDIF
            DO i = 1, nP
               Z(nZ+i) = Y(i)
               ENDDO
            CALL OlfPutZ (jout, K, nZ+nP, Z, Fehler)

            CALL OlfPutXYD (jout, K, n1, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO ! K

C  Save, end loop files :
         CALL TensorCheckZ (jout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfClos (jout, nK, Fehler)
         ENDDO ! lj

      END ! GridCurve

C  --------------------------------------------------------------------
      SUBROUTINE GridIntExt (qExt, nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 4-5jul91. Extrapolation 30jul91.
            ! Sum 19sep91. Spectrum-spectrum correpondances
            ! cleared 4feb93. IntExt separately 19nov96

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      INTEGER       MC10
      PARAMETER    (MC10=10*MC+50) !Artem: Increase the size of the working array.
      CHARACTER*(*) Fehler
      CHARACTER*80  aus, doc
      CHARACTER*40  h1, h2
      INTEGER       nJList, JList(*), lj, j, jout,
     *              K, Kout, nK, nKdummy, nC, nC1, nC2, i,
     *              iModI, iModEmax, iModEl, iModEr, iModEm, iElnn,
     *              iErnn, iOpt, iSort0, iSort1, ifail, ih1, ih2,
     *              ia, ie, i0n, i0l, i0h
      REAL*8        z, Work(MC10), tolI, valE, valEr, valEl, valEm, dx
      LOGICAL       qCurve, qExt, qOv, qAccept

C  Defaults for dialogue :
      DATA          iModI /1/, iModEl /1/, iModEr /1/,
     *              iModEm /2/, iElnn /1/, iErnn /1/

      IF     (qExt) THEN
         h2 = 'Extrapolate data from files'
      ELSE
         h2 = 'Interpolate data from files'
         ENDIF

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      CALL SetGridChoice (.true., doc, Fehler)

      Print *, ' interpolation mode :'
      Print *, '    (1) linear'
      Print *, '    (2) cubic spline'
      Print *, '    (3) don''t interpolate'
      Print *, '    (4) overwrite all'
      iModI = iAskDMu (' Option', iModI, 0, 4)
      tolI = 1.d-3
      IF (iModI.eq.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      IF      (qExt) THEN ! Extrapolation
         IF (iModI.le.2) THEN ! true extrapolation
            Print *, ' extrapolation mode :'
            Print *, '    (1) don''t extrapolate'
            Print *, '    (2) by 0'
            Print *, '    (3) by a global constant'
            Print *, '    (4) by a constant, from neighbours'
            iModEmax = 4
            CALL i2FrageD (' Option for left, right side',
     *         iModEl, iModEr, iModEl, iModEr)
            IF     (iModEl.eq.0 .or. iModEr.eq.0) THEN
               Fehler = ' '
               RETURN
            ELSEIF (qioutside(iModEl,1,iModEmax)
     *         .or. qioutside(iModEr,1,iModEmax)) THEN
               Fehler = 'invalid choice'
               RETURN
               ENDIF
         ELSE
            Print *, ' fill points :'
            Print *, '    (2) by 0'
            Print *, '    (3) by a global constant'
c            Print *, '    (5) asked inividually'
            iModEm = iAskDMu (' Option', iModEm, 0, 3)
            IF (iModEm.le.1) THEN
               Fehler = ' '
               RETURN
               ENDIF
            iModEl = iModEm
            iModEr = iModEm
            ENDIF

         IF (iModI.le.2) THEN
            IF    (iModEl.eq.3) THEN
               valEl   = rAskD
     *            (' Extrapolation: left constant', valEl)
            ELSEIF (iModEl.eq.4) THEN
               iElnn = iAskDMu
     *            (' Extrapolation: # left neighbours', iElnn, 1, MC)
               ENDIF
            IF    (iModEr.eq.3) THEN
               valEr   = rAskD
     *            (' Extrapolation : right constant', valEr)
            ELSEIF (iModEr.eq.4) THEN
               iErnn = iAskDMu
     *            (' Extrapolation : # right neighbours', iErnn, 1, MC)
               ENDIF
         ELSE
            IF     (iModEm.eq.3) THEN
               valEm = rAskD (' Extrapolation constant', valEm)
               valEl = valEm
               valEr = valEm
               ENDIF
            ENDIF

         ENDIF

C  ---------------------------------------------------------------------
C  Loop over files :
C  ---------------------------------------------------------------------

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurve) THEN
            Fehler = ' File is a curve'
            RETURN
            ENDIF

C  X-Scale per file :
         CALL SetGridJ (j, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  Documentation, open output file :
         IF (qExt) THEN
            CALL OlfComAdd (jout, 'E', 'extrapol''d/'//doc, Fehler)
         ELSE
            CALL OlfComAdd (jout, 'I', 'intrapol''d/'//doc, Fehler)
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

C  Loop spectra :
         DO K = 1, nK

            CALL OlfGetXYD (j, K, nC, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

C  X-Scale per spectrum :
            CALL SetGridK (K, nC1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

C  As a result, the new x-grid is X1.
C  Copy to X2 which may be overwritten in the following :     ! SCHLECHT
            nC2 = nC1
            DO i = 1, nC2
               X2(i) = X1(i)
               ENDDO

C  Now the second part of the K-loop, according to iModY :
            ! Check ascending order of grids :
            iSort0 = irSorted (X,  nC)
            iSort1 = irSorted (X2, nC2)
            IF (iSort0.ne.2) THEN
               Fehler = ' x-scale of file f0 is not sorted'
               RETURN
               ENDIF
            IF (iSort1.ne.2) THEN
               Fehler = ' x-scale of file f1 is not sorted'
               RETURN
               ENDIF

            ! Find the subrange of X2 which is covered by X :
            ia = irPos (X2, nC2, X(1),  'l')
            ie = irPos (X2, nC2, X(nC), 'r')
            IF (ia.gt.ie) THEN
               Fehler = ' there is no overlap of f0 and f1'
               RETURN
               ENDIF

C  And now the interpolation itself :
            IF     (iModI.eq.1) THEN
               ! simple linear interpolation :
               DO i = ia, ie
                  i0l = irPosOpt (X, nC, X2(i), 'r', i0l)
                  i0h = irPosOpt (X, nC, X2(i), 'l', i0h)
                  IF (i0l.lt.1 .or. i0h.gt.nC) THEN
                     Print *, ' i = '//cl4(i)
                     Fehler = ' Quatsch in interpolation routine'
                     RETURN
                     ENDIF
                  CALL LinIntPol (X2(i), X(i0l), X(i0h),
     *               Y(i0l), Y(i0h), D(i0l), D(i0h), Y1(i), D1(i))
                  ENDDO
            ELSEIF (iModI.eq.2) THEN ! this mode since 10sep91
               ! calculate the spline :
               ifail = 1 ! silent exit
               IF     (nC.lt.4) THEN
                  Fehler = ' Spectra (f0) too short'
                  RETURN
               ELSEIF (nC.gt.MC-4) THEN
                  Fehler = ' Spectra (f0) too long'
                  RETURN
                  ENDIF
               CALL E01BAF_local (nC, X, Y, Y2, D2, MC, Work,
     *                            MC10, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
                  ! The coefficients c(x) are returned in D2,
                  ! the knots lambda(x) in Y2.
               IF (ifail.ne.0) CALL Absturz ('OrgGrid',
     *            'E01BAF ifail='//cl4(ifail))

               ! evaluate the points :
               DO i = ia, ie
                  CALL E02BBF_local (nC+4, Y2, D2, X2(i), Y1(i), ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
                     ! point X2(i) -> value Y1(i)
                  IF (ifail.ne.0) CALL Absturz ('OrgGrid',
     *               'E01BBF ifail='//cl4(ifail))
                  D1(i) = 0. ! no error bars defined !FKwork error bars!
                  ENDDO

            ELSEIF (iModI.eq.3 .or. iModI.eq.4) THEN  ! this antimode 21jan92
               DO i = ia, ie
                  IF (iModI.eq.3) THEN ! new=old within tolerance ?
                     i0n = irPosOpt (X, nC, X2(i), 'n', i0n)
                     IF (nC.ge.2) THEN
                        dx = (X(nC)-X(1)) / (nC-1)
                     ELSE
                        dx = X(1)
                        ENDIF
 532                 CONTINUE
                     IF (qEqTol(X2(i), X(i0n), tolI)) THEN
                        qAccept = .true.
                     ELSE
                        CALL NiceNum (X2(i),  h2, ih2)
                        CALL NiceNum (X(i0n), h1, ih1)
                        Print *, ' found old x='//h1(1:ih1)//
     *                         ' for new '//h2(1:ih2)//' ?'
                        iOpt = iAskDMu (
     *                  ' Accept(1), skip(2), change tol(3)', 2, 0, 3)
                        IF     (iOpt.eq.0) THEN
                           Fehler = ' '
                           RETURN
                        ELSEIF (iOpt.eq.1) THEN
                           qAccept = .true.
                        ELSEIF (iOpt.eq.2) THEN
                           qAccept = .false.
                        ELSEIF (iOpt.eq.3) THEN
                           tolI = rAskD (' Tolerance', tolI)
                           ENDIF
                        ENDIF
                  ELSE
                     qAccept = .false.
                     ENDIF
                  IF (qAccept) THEN
                     Y1(i) = Y(i0n)
                     D1(i) = D(i0n)
                  ELSE
                     IF     (iModEm.eq.2) THEN
                        valE = 0.
                     ELSEIF (iModEm.eq.3) THEN
                        valE = valEm
                     ELSE
                        CALL Absturz ('OrgGrid', 'iModEm o.o.r.')
                        ENDIF
                     Y1(i) = valE
                     D1(i) = 0.
                     ENDIF
                  ENDDO

               ENDIF ! iModI

C  Left side :
            IF (qExt .and. iModEl.ge.2) THEN
C  - fill the rest by extrapolation :
               DO i = 1, ia-1
                  IF     (iModEl.eq.2) THEN
                     valE = 0.
                  ELSEIF (iModEl.eq.3) THEN
                     valE = valEl
                  ELSEIF (iModEl.eq.4) THEN
                     IF (ia+iElnn-1.gt.ie) THEN
                        Fehler = ' too many left neighbours required'
                        RETURN
                        ENDIF
                     valE = rSum (Y1, ia, ia+iElnn-1, 1) / iElnn
                     ENDIF
                  Y1(i) = valE
                  D1(i) = 0.
                  ENDDO
            ELSE
C  - eliminate the empty channels :
               nC2 = nC2 - ia + 1
               ie  = ie  - ia + 1
               DO i = 1, nC2
                  X2(i) = X2(ia+i-1)
                  Y1(i) = Y1(ia+i-1)
                  D1(i) = D1(ia+i-1)
                  ENDDO
               ia = 1
               ENDIF

C  Right side :
            IF (qExt .and. iModEr.ge.2) THEN
C  - fill the rest by extrapolation :
               DO i = ie+1, nC2
                  IF     (iModEr.eq.2) THEN
                     valE = 0.
                  ELSEIF (iModEr.eq.3) THEN
                     valE = valEr
                  ELSEIF (iModEr.eq.4) THEN
                     IF (ie-iErnn+1.lt.ia) THEN
                        Fehler = ' too many right neighbours required'
                        RETURN
                        ENDIF
                     valE = rSum (Y1, ie-iErnn+1, ie, 1) / iErnn
                     ENDIF
                  Y1(i) = valE
                  D1(i) = 0.
                  ENDDO
            ELSE
C  - eliminate the empty channels :
               nC2 = ie
               ENDIF

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, nC2, X2, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! loop spectra

         CALL OlfClos (jout, nK, Fehler)
         ENDDO ! lj

      END ! GridIntExt

C  --------------------------------------------------------------------
      SUBROUTINE GridSumCut (qSum, nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler
      INTEGER       nJList, JList(*), lj, j, jout, nsum,
     *              K, Kout, nK, nKdummy, nC, nC1, nC2, i, ii,
     *              i2, i2prev, iLC, iSort0
      REAL*8        z, xneu, yneu, dneu
      CHARACTER*40  doc
      LOGICAL       qSum, qCurve, qOv, qRetain, qLisC(MC)

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      CALL SetGridChoice (.false., doc, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      qRetain = qAskD (' Retain first and last group', intq(qRetain))

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurve) THEN
            Fehler = ' File is a curve'
            RETURN
            ENDIF

         CALL SetGridJ (j, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (qSum) THEN
            CALL OlfComAdd (jout, 's', 'summed chs/ '//doc, Fehler)
         ELSE
            CALL OlfComAdd (jout, 'd', 'deleted chs/ '//doc, Fehler)
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

         DO K = 1, nK

            CALL OlfGetXYD (j, K, nC, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            iSort0 = irSorted (X,  nC)
            IF (iSort0.ne.2) THEN
               Fehler = ' x-scale of file f0 is not sorted'
               RETURN
               ENDIF

            CALL SetGridK (K, nC1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

C  As a result, the new x-grid is X1.
C  Copy to X2 which may be overwritten in the following :     ! SCHLECHT
            nC2 = nC1
            DO i = 1, nC2
               X2(i) = X1(i)
               ENDDO

            ! qLisC(i) = {X(i) is first channel after a point X2} :
            i2prev = 0 ! previous entry in X2
            DO i = 1, nC
               i2 = irPosOpt (X2, nC2, X(i), 'r', i2prev)
               qLisC(i) = (i2.ne.i2prev)
               i2prev = i2
               ENDDO
               ! qLisC is determined. Now X2 is no
               ! longer used and will be overwritten.
            iLC = iqSum (qLisC, nC)
            IF (iLC.le.0) THEN
               Fehler = ' no channels retained'
               RETURN
               ENDIF
            ! Now sum/select the channels (copied from OrgChannels)
            IF (qRetain) THEN
               nC2 =  0
            ELSE
               nC2 = -1
               ENDIF
            DO i = 1, nC
               IF (qLisC(i)) THEN
                  nC2  = nC2 + 1
                  IF (nC2.eq.0) GOTO 419 ! skip the first group
                  IF (qSum) THEN
                     ! sum up to the next retained channel
                     nsum = 1
                     xneu = X(i)
                     yneu = Y(i)
                     dneu = D(i)**2
                     DO ii = i+1, nC
                        IF (qLisC(ii)) GOTO 415
                        nsum = nsum + 1
                        xneu = xneu + X(ii)
                        yneu = yneu + Y(ii)
                        dneu = dneu + D(ii)**2
                        ENDDO
 415                 CONTINUE
                     X2(nC2) = xneu / nsum
                     Y1(nC2) = yneu / nsum
                     D1(nC2) = dsqrt(dneu) / nsum
                  ELSE ! select just the channel i :
                     X2(nC2) = X(i)
                     Y1(nC2) = Y(i)
                     D1(nC2) = D(i)
                     ENDIF
                  ENDIF
 419           CONTINUE
               ENDDO
            IF (.not.qRetain) nC2 = nC2 - 1 ! skip the last group

            IF (nC2.le.0) nC2 = 0

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, nC2, X2, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! loop spectra

         CALL OlfClos (jout, nK, Fehler)

         ENDDO ! lj

      END ! GridSumCut

C  --------------------------------------------------------------------
      SUBROUTINE GridRedis (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler
      INTEGER       nJList, JList(*), lj, j, jout,
     *              K, Kout, nK, nKdummy, i, n, n1, n2, iSort0, jFill
      REAL*8        z
      CHARACTER*40  doc
      LOGICAL       qCurve, qOv

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      CALL SetGridChoice (.false., doc, Fehler)
      jFill = iAskDMu ('Fill voids: with two points(2) fully(3)',
     *     jFill, 0, 3)
      IF (jFill.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurve) THEN
            Fehler = ' File is a curve'
            RETURN
            ENDIF

         CALL SetGridJ (j, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, 'r', 'redistribute '//doc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         DO K = 1, nK

            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            iSort0 = irSorted (X, n)
            IF (iSort0.ne.2) THEN
               Fehler = ' x-scale of file f0 is not sorted'
               RETURN
               ENDIF

            CALL SetGridK (K, n1, Fehler) ! -> X1(1..n1)
            ! copy because X2,n2 may be modified in RedistrHistogr
            n2 = n1
            DO i = 1, n2
               X2(i) = X1(i)
               ENDDO
            IF (Fehler.ne.'&ff') RETURN

            CALL RedistrHistogr (MC, n,  X,  X3, Y,  D,
     *                           MC, n2, X2, X4, Y1, D1, jFill, Fehler)

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n2, X2, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! loop spectra

         CALL OlfClos (jout, nK, Fehler)

         ENDDO ! lj

      END ! GridRedis

C  --------------------------------------------------------------------
      SUBROUTINE GridHist (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler
      INTEGER       nJList, JList(*), lj, j, jout,
     *              K, Kout, nK, nKdummy, n, n1, i, i2,
     *              iHcol, iHcolIn, iSort0
      REAL*8        z, val
      CHARACTER*40  doc, co, un, coY, unY
      LOGICAL       qCurve, qOv
      DATA          iHcol / 2 /

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      CALL SetGridChoice (.false., doc, Fehler)

      iHcolIn = iAskDMu (' One-dimensional data from x(1) or y(2)',
     *                   iHcol, 0, 2)
      IF (iHcolIn.eq.0) RETURN
      iHcol = iHcolIn

      DO lj = 1, nJList
         j = JList(lj)

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurve) THEN
            Fehler = ' File is a curve'
            RETURN
            ENDIF

         CALL OlfHeadDup (j, qOv, jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF     (iHcol.eq.1) THEN
            CALL OlfCnuG (j, 'x', co, un, Fehler)
         ELSEIF (iHcol.eq.2) THEN
            CALL OlfCnuG (j, 'y', co, un, Fehler)
         ELSE
            Fehler = 'PROG ERR/ iHcol ooR'
            RETURN
            ENDIF
         CALL OlfCnuP (jout, 'x', co, un, Fehler)
         coY = '#['//co(1:lenU(co))//']'
         IF (un.ne.' ') THEN
            unY = un(1:lenU(un))//'^-1'
         ELSE
            unY = ' '
            ENDIF
         CALL OlfCnuP (jout, 'y', coY, unY, Fehler)

         CALL SetGridJ (j, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, 'h', 'group into histogra'//doc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         DO K = 1, nK

            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            iSort0 = irSorted (X, n)
            IF (iSort0.ne.2) THEN
               Fehler = ' x-scale of file f0 is not sorted'
               RETURN
               ENDIF

            CALL SetGridK (K, n1, Fehler) ! -> X1(1..n1)
            IF (Fehler.ne.'&ff') RETURN

            DO i2 = 1, n1
               Y1(i2) = 0
               ENDDO
            DO i = 1, n
               IF     (iHcol.eq.1) THEN
                  val = X(i)
               ELSEIF (iHcol.eq.2) THEN
                  val = Y(i)
               ELSE
                  Fehler = 'PROG ERR/ Hist/ iHcol oor'
                  RETURN
                  ENDIF
               i2 = irPosOpt (X1, n1, val, 'n', i2)
               Y1(i2) = Y1(i2) + 1
               ENDDO
            DO i2 = 1, n1
               D1(i2) = dsqrt (Y1(i2))
               ENDDO

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, n1, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO ! loop spectra

         CALL OlfClos (jout, nK, Fehler)

         ENDDO ! lj

      END ! GridHist

C  --------------------------------------------------------------------
      SUBROUTINE SetGrid ()
C  --------------------------------------------------------------------
            ! separated from OrgGrid 19nov96

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f' ! uses X1
      INCLUDE 'l_def.f'

      CHARACTER*(*) Fehler, doc
      CHARACTER*40  fil1, h1
      INTEGER       nsg, j, j1, K, nKdummy, nout
      INTEGER, save::iModX, iCorr, K1sel, nC1, nK1, K2K(MK), !Artem: Added "save" to persist data across entry points.
     * j1in
      REAL*8        sg1, sgn, z1
      LOGICAL       qCommens

C  Defaults for dialogue :
      DATA          iModX /11/, iCorr /2/, K1sel/1/

C  ---------------------------------------------------------------------
      ENTRY SetGridChoice (qCommens, doc, Fehler)
C  ---------------------------------------------------------------------

      IF (Fehler.ne.'&ff') THEN
         Print *, ' error on entry in SetGridChoice'
         RETURN
         ENDIF
      IF (qCommens) THEN ! old and new file have commensurable x grids
         Print *, ' new x-grid :'
         Print *,
     * '    regular grid, lin(1) 1/2-log(2) log(3) lin-blocks(4)'
         Print *, '    from a file, selected(11) 1:1(12)'
         Print *, '    set points(21)'
      ELSE
         Print *, ' x-grid :'
         Print *,
     * '    regular grid, lin(1) 1/2-log(2) log(3) lin-blocks(4)'
         Print *, '    from a file, any(11)'
         Print *, '    set points(21)'
         ENDIF
      iModX = iAskDMu (' Option', iModX, 0, 99)

      IF     (iModX.eq.0) THEN
         Fehler = ' '

      ELSEIF (qiinside(iModX,1,4)) THEN
         CALL SetGridReg (iModX, sg1, sgn, nsg, X1, MC,
     *                    nC1, doc, Fehler)

      ELSEIF (iModX.eq.11) THEN ! == new grid from file ==
         j1in = iAskD (' From file', j1in)
         IF     (j1in.le.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         nK1 = iOlfG (j1in, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         h1 = ' '
         IF (nK1.eq.1) THEN
            iCorr = 2
            K1sel = 1
         ELSE
            iCorr = iAskDMu (
     *         ' spectrum:spectrum(1) or select one(2)', iCorr, 0, 2)
            IF     (iCorr.eq.0) THEN
               Fehler = ' '
               RETURN
            ELSEIF (iCorr.eq.2) THEN
               K1sel = iAskD (' Select spectrum', K1sel)
               h1 = ' spectrum '//cl3(K1sel)
               ENDIF
            ENDIF
         CALL tOlfG (j1in, 'fil', fil1, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL Compose2 (doc, 'grid from '//fil1, h1)

      ELSEIF (iModX.eq.12) THEN ! == grid from same file ==
         j1in = 0
         iCorr = iAskDMu (' spectrum:spectrum(1) or select one(2)',
     *                    iCorr, 0, 2)
         IF     (iCorr.eq.0) THEN
            Fehler = ' '
            RETURN
         ELSEIF (iCorr.eq.2) THEN
            K1sel = iAskD (' Select spectrum', K1sel)
            h1 = ' spectrum '//cl3(K1sel)
            ENDIF
         CALL Compose2 (doc, 'grid from input x', h1)

      ELSEIF (iModX.eq.21) THEN ! == set points ==
         Print *, ' set points x :'
         CALL rAskArray (X1, MC, nC1)
         doc = 'points entered'

      ELSE
         Fehler = 'invalid option'
         ENDIF

      RETURN ! SetGridChoice

C  ---------------------------------------------------------------------
      ENTRY SetGridJ (j, Fehler)
C  ---------------------------------------------------------------------
         ! in loop over files: select X1 if not done previously

      IF (Fehler.ne.'&ff') THEN
         Print *, ' error on entry in SetGridJ'
         RETURN
         ENDIF
      IF     (qiinside(iModX,11,12)) THEN
         IF (j1in.eq.0) THEN
            j1 = j
         ELSE
            j1 = j1in
            CALL OlfCnuCheck2 (j, j1, 'x', 0, 1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDIF
         IF     (iCorr.eq.1) THEN ! set 1:1 correspondance
            CALL GetK2K (j, j1, nKdummy, K2K, Fehler)
         ELSEIF (iCorr.eq.2) THEN ! get selected spectrum
            CALL OlfGetXYD (j1, K1sel, nC1, X1, Y1, D1, Fehler)
            ENDIF
         ENDIF

      RETURN ! SetGridJ

C  ---------------------------------------------------------------------
      ENTRY SetGridK (K, nout, Fehler)
C  ---------------------------------------------------------------------
         ! in loop over files: select X1 if not done previously

      IF (Fehler.ne.'&ff') THEN
         Print *, ' error on entry in SetGridK'
         RETURN
         ENDIF
      IF     (qiinside(iModX,11,12)) THEN
         IF (iCorr.eq.1) THEN
            CALL OlfGetXYD (j1, K2K(K), nC1, X1, Y1, D1, Fehler)
            ENDIF
         ENDIF

      nout = nC1

      RETURN ! SetGridK

      END ! SetGrid

C  --------------------------------------------------------------------
      SUBROUTINE SetGridReg (iMod, s1, sn, ns, X, M, n, doc, Fehler)
C  --------------------------------------------------------------------
            ! separated from OrgGrid for use in TraFilon 27apr93
         ! set regular grid (lin log ..)

         ! Import :
         !    iMod = 1..4   type of grid
         !    M             declared dimension of X
         ! Export :
         !    X(1..n)       grid to be set here
         !    doc           some words about this grid
         ! Work space (don't modify outside) :
         !    s1,sn,ns      defaults

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      CHARACTER*(*) doc, Fehler
      CHARACTER     aus*80, h1*40, h2*40, h3*40, cl6*6
      DIMENSION     X(M)

      IF     (iMod.eq.1) THEN ! linear grid
         CALL rAskGridDef ('x', X, M, n, xL1, xLd, xLn)
         IF (n.le.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         xL1 = X(1)        ! defaults for next call
         xLd = X(2) - X(1)
         xLn = X(n)
         CALL NiceNum (X(1), h1, ih1)
         CALL NiceNum (X(n), h2, ih2)
         IF (n.gt.1) THEN
            CALL NiceNum ((X(n)-X(1))/(n-1), h3, ih3)
         ELSE
            h3 = '-'
            ih3 = 1
            ENDIF
         doc = ' '//h1(1:ih1)//'('//h3(1:ih3)//')'//h2(1:ih2)

      ELSEIF (iMod.eq.2) THEN ! log (x+1) - scale : see Protokollbuch C2,31.
         sn  = rAskDLu (' Maximum |x|', sn, 0.d0, 1.d16)
         s1 = rAskDLu (' Crossover linear - logarithmic',
     *                 s1, 0.d0, 1.d16)
         IF (sn.le.0. .or. s1.le.0.) THEN
            Fehler = ' '
            RETURN
            ENDIF
         r10 = rAskDLu (
     * ' Number of points per decade (non-integer allowed)',
     * dble(ns), 1.d-1, dble(M/2))
         ns = idnint(r10)
         ! set scale
         alpha = 2.302585 / r10
         nm = dln0 (sn/s1-1) / alpha + 1
         X(nm+1) = 0.
         DO i = 1, nm
            X(nm+1-i) = - s1 * (dexp(alpha*i) - 1)
            ENDDO
         np = nm
         DO i = 1, np
            X(nm+1+i) =   s1 * (dexp(alpha*i) - 1)
            ENDDO
         n = nm + 1 + np
         CALL NiceNum (s1, h1, ih1)
         CALL NiceNum (sn, h2, ih2)
         CALL NiceNum (r10, h3, ih3)
         doc = '1/2log '//h1(1:ih1)//'..'//h2(1:ih2)//'; '//
     *         h3(1:ih3)//'/dec'

      ELSEIF (iMod.eq.3) THEN ! logarithmic grid
 1331    CONTINUE
         s1  = rAskD (' Enter X(1)', s1)
         IF (s1.le.0) THEN
            CALL Gong (1)
            IF (qAskD('Invalid - correct', 1)) GOTO 1331
            Fehler = ' '
            RETURN
            ENDIF
         r1   = dlog10 (s1)
         sn  = rAskD (' Enter X(n)', sn)
         IF (sn.le.s1) THEN
            CALL Gong (1)
            IF (qAskD('Invalid - correct', 1)) GOTO 1331
            Fehler = ' '
            RETURN
            ENDIF
         rn   = dlog10 (sn)
         rD   = rAskDLu (' Points per decade',
     *                   dble(ns), 1.d0, dble(M-1)/(rn-r1))
         ns = idnint(rD)
         n = idint ( (rn-r1)*rD ) + 1
         IF (n.le.1) THEN
            Fehler = ' less than two points in grid f1'
            RETURN
            ENDIF
         DO i1 = 1, n
            X(i1) = 10.**(r1+(i1-1)*(rn-r1)/(n-1))
            ENDDO
         CALL NiceNum (s1, h1, ih1)
         CALL NiceNum (sn, h2, ih2)
         CALL NiceNum (rD, h3, ih3)
         doc = 'log '//h1(1:ih1)//'..'//h2(1:ih2)//'; '//
     *         h3(1:ih3)//'/dec'

      ELSEIF (iMod.eq.4) THEN ! linear blocks (for Filon integration)
         s1 = rAskDLu (' Smallest step', s1, 1.d-20, 1.d20)
         sn = rAskDLu (' Highest value', sn, s1,    1.d20)
         ns = iAskDMu (' Points per block', ns, 0, (M-1)/2)
         IF (ns.lt.1) THEN
            Fehler = ' '
            RETURN
            ENDIF

         step = s1
         DO i = 1, 2*ns+1
            X(i) = (i-1) * step
            ENDDO
         n = 2*ns+1
 1340       CONTINUE
            step = step * 2
            IF (n+ns.gt.M) THEN
               Fehler = ' too many points required'
               RETURN
               ENDIF
            DO i = 1, ns
               X(n+i) = X(n) + i*step
               ENDDO
            n = n + ns
            IF (X(n).lt.sn) GOTO 1340

            CALL NiceNum (s1, h1, ih1)
            CALL NiceNum (sn, h2, ih2)
            CALL Compose2 (aus,
     *  'lin '//h1(1:ih1)//'..'//h2(1:ih2)//'; '//cl6(ns), '/block')
            doc = aus

         ELSE

            Fehler = 'SetGridReg/ iMod not implemented'
            ENDIF

         END ! SetGridReg

C  ====================================================================
C  i40 / 3 :   OrgHist...
C  ====================================================================

C  --------------------------------------------------------------------
      SUBROUTINE OrgHistMake (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! binning of a list as a histogram
         ! JWu 25jul95

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)
      DIMENSION      qList(MC), qList2(MC), iShift(MK)

      CHARACTER*40  co, un
      CHARACTER*80  aus, LisDoc
      CHARACTER*320 LisC, LisCofK

      DATA          LisC /' '/, LisCofK /' '/, qAllJ /.true./

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop files :
      DO lj = 1, nJList
         j = JList(lj)


         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (j,    'y', co, un, Fehler)
         CALL OlfCnuP (jout, 'x', co, un, Fehler)

         CALL Compose2 (co, 'g('//co, ')')
         CALL Compose2 (un, un, '-1')
         CALL OlfCnuP (jout, 'y', co, un, Fehler)

         CALL OlfComAdd (jout, 'h', 'histogram binning, n='//
     *                   cl4(nb), Fehler)
         IF (Fehler.ne.'&ff') RETURN

         nb = iAskDMu ('Binning into how many channels', nb, 1, MC)

C  Loop spectra :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ! Min/Max ermitteln :
            ymax = Y(1)
            ymin = Y(1)
            DO i = 1, n
               ymax = dmax1 (ymax, Y(i))
               ymin = dmin1 (ymin, Y(i))
               ENDDO

            ! Kanalbreite im Histogramm :
            ychan = (ymax - ymin) / (nB - 1)

            ! Histogramm setzen :
            DO i = 1, nB
               X1(i) = ymin + ychan * (i-0.5)
               Y1(i) = 0
               D1(i) = 0
               ENDDO

            ! Histogramm hochz"ahlen :
            DO i = 1, n
               i1 = irPosOpt (X1, nB, Y(i), 'n', i1)
               Y1(i1) = Y1(i1) + 1
               ENDDO

            CALL OlfCopZ (j, jout, K, K, Fehler)
            CALL OlfPutXYD (jout, K, nB, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ENDDO

         CALL OlfClos (jout, nK, Fehler)

C  End loop files :
         ENDDO

      END ! OrgHistMake
