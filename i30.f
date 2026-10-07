C  ====================================================================
C
C      Library  IDA   :  Inelastic data treatment
C      Modul    i30   :     file manipulations, auxiliary calculations
C
C  ====================================================================

C     Contents :
C        1.  Organised data access :
C               TakeTwoSpectra, GetIndex, GetK2K, GetJ2J, CommonScale,
C               SelectCh*
C        2.  Auxiliary operations :
C               SortChannels
C        3.  Auxiliary calculations :
C               IdaMinMax, Integrate, RedistrHistogr, CoordBins, CalcQQ

C  ====================================================================
C  i30 / 1 :    Organised data access
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE TakeTwoSpectra
     * (j1, j2, K1, K2, nZ, Z, nC, X, Y1, Y2, D1, D2,tolerance,Fehler)
C     --------------------------------------------------------------------
            ! JWu 31jul91 separated from TraFourier
         ! Read two spectra that are supposed to have a
         ! common x scale. If j1 or j2 is 0, the corresponding
         ! data Y and D are set to 0.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'

      DIMENSION     X(MC), Y1(MC), D1(MC), X2(MC), Y2(MC), D2(MC),
     *              Z(MZ), Z2(MZ)

      CHARACTER     Fehler*(*)
      CHARACTER     hilf*80, cl6*6

      IF (j1.eq.0 .and. j2.eq.0) THEN
         Fehler = ' No files to be read'
         ENDIF

C  If only one spectrum is given :
      IF (j2.eq.0) THEN
         CALL OlfGetSpe (j1, K1, nZ, Z, nC, X, Y1, D1, Fehler)
         DO i = 1, MC
            Y2 (i) = 0.
            D2 (i) = 0.
            ENDDO
         RETURN ! all is done
         ENDIF
      IF (j1.eq.0) THEN
         CALL OlfGetSpe (j2, K2, nZ, Z, nC, X, Y2, D2, Fehler)
         DO i = 1, MC
            Y1 (i) = 0.
            D1 (i) = 0.
            ENDDO
         RETURN ! all is done
         ENDIF

C  Read two spectra :
      CALL OlfGetSpe (j1, K1, nZ,  Z,  nC,  X,  Y1, D1, Fehler)
      CALL OlfGetSpe (j2, K2, nZ2, Z2, nC2, X2, Y2, D2, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Prepare error message (if used or not) :
      IF (K1.eq.K2) THEN
         hilf = ' in spectrum '//cl6(K1)
      ELSE
         hilf = ' in spectra '//cl6(K1)
         CALL Append (hilf, ', '//cl6(K2))
         ENDIF

C  Checks of consistency :
      IF     (nC.ne.nC2) THEN
         Fehler = 'different number of channels'//hilf
         RETURN
      ELSEIF (nZ.ne.nZ2) THEN
         Fehler = 'different number of Z'//hilf
         RETURN
         ENDIF
      DO iZ = 1, nZ
         IF (Z(iZ).ne.Z2(iZ)) THEN
            Fehler = 'different z'//hilf
            RETURN
            ENDIF
         ENDDO
      dx = dabs(X(nC)-X(1))/(nC-1)
      DO i = 1, nC
         IF (dabs(X(i)-X2(i)).gt.tolerance*dx) THEN
            Print *, ' i, X1(i), X2(i) : ', i, X(i), X2(i)
            Fehler = 'different x-scale'//hilf
            RETURN
            ENDIF
         ENDDO

      END ! TakeTwoSpectra

      SUBROUTINE GetIndex (ein, X1, n1, X2, n2, I2, tol, Fehler)
C     ----------------------------------------------------------
            ! JWu 16sep91. Generalizing a part of OprFunctional. 4oct91.
         ! Get index of X1 in X2 : X1(i) = X2(I2(i)) with precision tol.
         ! Special case tol<0 : simply by channel, I2(i) = i.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      DIMENSION     X1(*), X2(*), I2(*) ! length can be MC or MK ...
      CHARACTER     Fehler*(*), ein*(*), aus*80, cl6*6, cr4*4, cx*16

 1    CONTINUE

      IF (tol.lt.0.) THEN ! by channel number
         IF (n2.ne.n1) THEN
            CALL Compose4 (Fehler, ein,
     *         ' different number of entries (n1,n2 = '//cl6(n1),
     *         ','//cl6(n2), ')')
            RETURN
            ENDIF
         DO i = 1, n1
            I2(i) = i
            ENDDO

c ELSEIF (n2.eq.1) THEN
c ! special case : tol must be interpreted as absolute precision.
c IF (qEqTol(X1(1), X2(1), tol)) THEN
c I2(1) = 1
c ELSE
c Fehler = ' The only channels are not equal for f1 and f2'
c ENDIF

      ELSE
         IF (n1.eq.n2) THEN
            ! are grids identical ? (26mar93)
            q1eq2 =  .true.
            DO i = 1, n1
               IF (X1(i).ne.X2(i)) THEN
                  q1eq2 = .false.
                  GOTO 121
                  ENDIF
               ENDDO
 121        CONTINUE
            IF (q1eq2) THEN
               DO i = 1, n1
                  I2(i) = i
                  ENDDO
               RETURN
               ENDIF
            ENDIF
         ! otherwise proceed :
         IF (n2.gt.1 .and. iabs(irSorted(X2,n2)).ne.2) THEN
            CALL Compose2 (Fehler, ein, '/ f2 not sorted')
            RETURN
            ENDIF
         ! find correspondences :
         qTrouble = .false.
         DO i = 1, n1
            I2(i) = irPosOpt (X2, n2, X1(i), 'n', I2(i)+1)
            ! stepwidth dx for tol-limit :
            IF (i.eq.1) THEN
               dx = X1(2)  - X1(1)
            ELSEIF (i.eq.n1) THEN
               dx = X1(n1) - X1(n1-1)
            ELSE
               dx =(X1(i+1)- X1(i-1) )/2
               ENDIF
            dx = dabs(dx)
            ! check whether |X2(i2)-X1(i)| within tolerance limit :
            IF     (dx.lt.tol*dabs(X1(i)) .or. dx.lt.1.d-14) THEN
               qTrouble = .true.
               I2(i) = -I2(i) ! use - as flag
            ELSEIF (dabs(X1(i)-X2(I2(i))).gt.tol*dx) THEN
               qTrouble = .true.
               I2(i) = -I2(i)
               ENDIF
            ENDDO
         IF (qTrouble) THEN
            ! tabulate ill correspondences (6oct93) :
            Print *, ein(1:lenU(ein)), '/ correspondence not clear'
            ilast = 0
            i = 1
 4          CONTINUE
               IF (I2(i).lt.0) THEN
                  I2(i) = -I2(i)
                  ! normally, tabulate from I2-2 to I2+2
                  iti = max0 (I2(i)-2, ilast+1)
                  itf = min0 (n2,I2(i)+2)
                  ii  = i
                  DO it = iti, itf
                     write (cx, '(g12.5)') X2(it)
                     aus = ' f2('//cr4(it)//')='//cx
                     IF (it.eq.I2(i)) THEN ! closest correspondence
                        write (cx, '(g12.5)') X1(i)
                        CALL Append (aus,
     *                     ' <-- f1('//cr4(i)//')='//cx(1:12))
                        ! are there more X1 having the same X2 correspondence ?
                        DO ii = i, n1-1
                           IF (it.ne.iabs(I2(ii+1))) GOTO 425
                           ENDDO
 425                    CONTINUE
                        IF (ii.gt.i) THEN ! yes, there are
                           write (cx, '(g12.5)') X1(ii)
                           CALL Append (aus,
     *                  ' ... f1('//cr4(ii)//')='//cx(1:12)//' ???')
                        ELSE
                           CALL Append (aus, ' ???')
                           ENDIF
                     ELSEIF (it.eq.-I2(min0(ii+1,n1))) THEN
                        ! there will be trouble again, don't show it now
                        ilast = it-1
                        GOTO 48
                        ENDIF
                     Print '(a)', aus
                     ENDDO
                  ilast = itf
 48               CONTINUE
                  I2(i) = -I2(i)
                  i = ii
                  ENDIF
               i = i + 1
               IF (i.le.n1) GOTO 4
            ! What shall we do ?
            iOpt = iAskMu (
     * 'Escape(0) Accept(1) Set tol(2) By hand(3) One-to-one(4) ?',
     * 0, 4)
            IF     (iOpt.eq.0) THEN
               Fehler = ' '
               RETURN
            ELSEIF (iOpt.eq.1) THEN
               DO i = 1, n1
                  IF (I2(i).lt.0) I2(i) = -I2(i) ! delete - flag
                  ENDDO
               RETURN ! be happy
            ELSEIF (iOpt.eq.2) THEN
               tol = rAskDLu (' tolerance', tol, 0.d0, 1.d8)
               GOTO 1 ! try everything again
            ELSEIF (iOpt.eq.3) THEN
               DO i = 1, n1
                  IF (I2(i).lt.0) THEN
                     write (cx, '(g12.5)') X1(i)
                     CALL Compose2 (aus, ' Correspondence for f1('//
     *                              cl6(i), ')='//cx)
                     I2(i) = iAskDMu (aus, -I2(i), 0, n2)
                     IF (I2(i).eq.0) THEN ! escape
                        Fehler = ' '
                        RETURN
                        ENDIF
                     ENDIF
                  ENDDO
            ELSEIF (iOpt.eq.4) THEN
               tol = -1.
               GOTO 1 ! take other case
               ENDIF
            ENDIF
         ENDIF

      END ! GetIndex

      SUBROUTINE GetK2K (j1, j2, nK1, K2K, Fehler)
C     --------------------------------------------
            ! JWu 27nov92
         ! find for each spectrum K of file j1
         ! a spectrum K2K(K) of file j2.
         ! Along with K2K, return nK1.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      CHARACTER         Fehler*(*), aus*80
      DIMENSION         K2K(*), Z1(MK), Z2(MK)

      DATA  tolZ /1.d-5/

      nK1 = iOlfG (j1, '#spectra', Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'GetK2K (file 1, nK) ')
         RETURN
         ENDIF
      nK2 = iOlfG (j2, '#spectra', Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'GetK2K (file 2, nK) ')
         RETURN
         ENDIF

      IF     (nK2.le.0) THEN
         Fehler = '2nd file is empty'
         RETURN
      ELSEIF (nK2.eq.1) THEN
         ! only possibility :
         DO K1 = 1, nK1
            K2K(K1) = 1
            ENDDO
         RETURN
         ENDIF

      nZ1 = iOlfG (j1, '#Z', Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'GetK2K (file 1, nZ) ')
         RETURN
         ENDIF
      nZ2 = iOlfG (j2, '#Z', Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'GetK2K (file 2, nZ) ')
         RETURN
         ENDIF

      IF (nZ1.eq.0) THEN
         IF (nZ2.eq.0) THEN
            IF (nK1.eq.1 .and. nK2.eq.1) THEN
               K2K(1) = 1
               RETURN
            ELSE
               Fehler = 'GetK2K/ nZ=0 but NK<>1'
               RETURN
               ENDIF
         ELSE
            Fehler = 'GetK2K/ file 1 has nZ=0, file 2 has nZ>0'
            RETURN
            ENDIF
      ELSEIF (nZ2.eq.0) THEN
         Fehler =
     * 'GetK2K/ file 2 has nZ=0, but nK>1 (should never happen)'
         RETURN
         ENDIF

      CALL OlfGet1ZofK (j1, 1, nK1, Z1, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'GetK2K (file 1) ')
         RETURN
         ENDIF
      CALL OlfGet1ZofK (j2, 1, nK2, Z2, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'GetK2K (file 2) ')
         RETURN
         ENDIF

      aus = 'setting spectrum-spectrum correspondences'
      CALL GetIndex (aus, Z1, nK1, Z2, nK2, K2K, tolZ, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! GetK2K

      SUBROUTINE AskK2K (j1, j2, nK1, K2K, Fehler)
C     --------------------------------------------
            ! JWu 4feb93 - never used ?
         ! find for each spectrum K of file j1
         ! a spectrum K2K(K) of file j2.
         ! Along with K2K, return nK1.
         ! Contrarily to GetK2K, correspondences
         ! can be freely by the user.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      CHARACTER         Fehler*(*), aus*80
      DIMENSION         K2K(*), Z1(MK), Z2(MK)

      DATA  tolZ /1.d-5/, iCorr /1/, K2 /1/

      CALL OlfGet1ZofK (j1, 1, nK1, Z1, Fehler)
      CALL OlfGet1ZofK (j2, 1, nK2, Z2, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      IF     (nK2.le.0) THEN
         Fehler = '2nd file is empty'
         RETURN
      ELSEIF (nK2.eq.1) THEN
         ! only possibility :
         DO K1 = 1, nK1
            K2K(K1) = 1
            ENDDO
      ELSE
         Print *, ' Setting spectrum-spectrum correspondences :'
         Print *, '    (1) by z-value'
         Print *, '    (2) one common scale for all spectra'
         iCorr = iAskDMu (' Option', iCorr, 0, 2)
         IF     (iCorr.eq.0) THEN
            Fehler = ' '
            RETURN
         ELSEIF (iCorr.eq.1) THEN
            aus = 'getting 1:1 spectrum-spectrum correspondences'
            CALL GetIndex (aus, Z1, nK1, Z2, nK2, K2K, tolZ, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (iCorr.eq.2) THEN
            K2 = iAskD (' Common grid from spectrum no.', K2)
            ENDIF
         ENDIF

      END ! AskK2K

      SUBROUTINE GetJ2J (Fra, JLis, nJLis, J2Lis, iKmod, Fehler)
C     ----------------------------------------------------------
            ! JWu 19mar93; consequent use of JLis 9jul93
         ! set file-file correspondences, for OprPointwise

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      CHARACTER*(*)     Fra, Fehler
      DIMENSION         JLis(*), J2Lis(*)
      CHARACTER         aus*80, cLis*80, cLis2*80, cl6*6

      IF (qErrEntry('GetJ2J', Fehler)) RETURN

      IF     (nJLis.le.0) THEN
         CALL Absturz('GetJ2J', 'nJL<=0')
      ELSE
         CALL EncJList (nJLis, JLis, cLis)
         CALL Say2 (Fra, ' for files '//cLis(1:lenU(cLis)))
         IF (nJLis.eq.1) THEN
            aus = ' .. from file'
         ELSE
            aus = ' .. from files'
            ENDIF
 11      CONTINUE
         CALL GetJList (aus, cLis2, nJLis, n2JL, J2Lis, 1, 1024)

         IF     (n2JL.le.0) THEN
            Fehler = ' '
            RETURN
         ELSEIF (n2JL.eq.1) THEN
            DO lj = 2, nJLis
               J2Lis(lj) = J2Lis(1)
               ENDDO
         ELSEIF (n2JL.eq.nJLis) THEN
            ! tout est bien
         ELSE
            CALL Gong (9)
            CALL Compose2 (aus,
     *         ' Enter 1 or '//cl6(nLis), ' filenames')
            GOTO 11
            ENDIF
         ENDIF

      END ! GetJ2J

      SUBROUTINE CommonScale (j, qInf, qComN, nCom,
     *                        qComSc, nXC, XCom, iShift, Fehler)
C     ----------------------------------------------------------
         ! JWu 16-20feb91; 3jun91
      ! this subroutine investigates whether the spectra of one
      ! file have a common x-scale, i.e. an x-scale of which
      ! the scales of the different spectra are subranges.

      ! Input :
      !  j     = number of file to be investigated
      !  qInf  = type result immediately on screen ?
      ! Output :
      !  qComN = is the number of channels the same for all spectra ?
      !  nCom  = number of channels, if qComN
      !  qComSc= is there a common x-scale ?
      !  nXC   = number of channels of the common scale, if qComSc
      !  XCom  = the common scale, if qComSc
      !  iShift= the offset in XCom, i.e. X(1,K) = XCom(1+iShift(K)).

      ! This subroutine was written for use in the IDA options
      ! >>select/sum channels<< and >>exchange x and z scales<<.

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

C  Data set :
      INCLUDE 'i_dim.f'
      DIMENSION     XCom(MC), X1(MC), iShift(MK)

C  Diverse :
      CHARACTER     Fehler*(*), cl4*4, aus*80

      IF (Fehler.ne.'&ff') THEN
         Print *, ' Error on entry in CommonScale'
         RETURN
         ENDIF

      DATA          tol / 1.d-5 /

C  Initialize the common spectrum XCom :
      nK = iOlfG (j, '#spectra', Fehler)
      CALL OlfGetX (j, 1, nC, XCom, Fehler)
      IF (nC.lt.1) Fehler =
     * ' checking x scales: first spectrum is empty'
      IF (Fehler.ne.'&ff') RETURN

      nCom      = nC
      nXC       = nC
      iShift(1) = 0
      qComN     = .true.
      qComSc    = .true.

C  Loop spectra :
      DO K = 2, nK
C  Compare XCom with the spectrum X1=X(_,K) :
         CALL OlfGetX (j, K, nC, X1, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (qComN) THEN
            IF (nC.ne.nCom) THEN
               qComN = .false.
               aus = ' spectra 1 and '//cl4(K)
               CALL Append (aus, ' have different number of channels')
               IF (qInf) Print *, aus
               ENDIF
            ENDIF

         IF (.not.qComSc) GOTO 59 ! end of loop

C  The rest of this loop only if so far there is a common scale :

C  Determine offset :
         IF ( qEqTol(X1(1), XCom(1), tol) ) THEN
C  X1 has zero offset :
            iShift(K) = 0
         ELSE
            IF (X1(1).lt.XCom(1)) THEN
C  X1 has negative offset :
               DO i = 2, nC
                  IF ( qEqTol(X1(i), XCom(1), tol) ) THEN
                     iShift(K) = 1 - i
                     GOTO 52
                     ENDIF
                  ENDDO
C  - X1 not contained in XCom :
               aus = ' spectrum '//cl4(K)
               CALL Append (aus, ' does not match common x scale')
               IF (qInf) Print *, aus
               qComSc = .false.
               GOTO 59
 52            CONTINUE
C  - Shift old channels of XCom :
               DO i = nXC, 1, -1
                  XCom(i-iShift(K)) =  XCom(i)
                  ENDDO
               nXC = nXC - iShift(K)
C  - New channels to XCom :
               DO i = 1, -iShift(K)
                  XCom(i) = X1(i)
                  ENDDO
C  - Shift the shifts :
               DO KK = 1, K
                  iShift(KK) = iShift(KK) - iShift(K)
                  ENDDO
            ELSE
C  X1 has positive offset :
               DO i = 2, nXC
                  IF ( qEqTol(X1(1), XCom(i), tol) ) THEN
                     iShift(K) = i - 1
                     GOTO 54
                     ENDIF
                  ENDDO
C  - X1 not contained in XCom :
               aus = ' spectrum '//cl4(K)
               CALL Append (aus, ' does not match common x scale')
               IF (qInf) Print *, aus
               qComSc = .false.
               GOTO 59
 54            CONTINUE
               ENDIF ! offset < or > 0
            ENDIF ! offset = or <> 0
C  Check the remaining channels :
         DO i = 2, nC
            IF (i+iShift(K).gt.nXC) THEN
C  - append new channel to XCom :
               nXC = nXC + 1
               IF (nXC.ne.i+iShift(K)) THEN
                  Fehler =  ' inconsistency while extending XCom'
                  RETURN
                  ENDIF
               XCom(nXC) = X1(i)
            ELSE
C  - check :
               IF (.not.qEqTol(X1(i), XCom(i+iShift(K)), tol) ) THEN
                  aus = ' x scales do not match at channel '//
     *                  cl4(i+iShift(K))
                  CALL Append (aus, ' of spectrum '//cl4(K))
                  IF (qInf) Print *, aus
                  qComSc = .false.
                  GOTO 59
                  ENDIF
               ENDIF
            ENDDO ! i
 59      CONTINUE
         ENDDO ! loop K

C  Loop ended

      IF (.not.qComSc) THEN
         DO K = 1, nK
            iShift(K) = 0
            ENDDO
         ENDIF

C  Info :
      IF (qInf) THEN
         IF     (qComN .and. qComSc .and. nCom.eq.nXC) THEN
            aus = ' file has a rectangular grid ('//cl4(nCom)
            CALL Append (aus, ' channels)')
         ELSEIF (qComN .and. qComSc) THEN
            aus = ' common scale, constant number of channels ('//
     *            cl4(nCom)
            CALL Append (aus, '), but no rectangular grid')
         ELSEIF (qComN) THEN
            aus = ' spectra have constant length ('//cl4(nCom)
            CALL Append (aus, ') but different scales')
         ELSEIF (qComSc) THEN
            aus = ' spectra have different length but common scale'
         ELSE
            aus = ' file has unregular scales'
            ENDIF
         CALL Insert (aus, 1, ' x-scale/')
         Print *, aus
         ENDIF

      END ! CommonScale

C  --------------------------------------------------------------------
      SUBROUTINE SelectCh ()
C  --------------------------------------------------------------------
         ! select channels
            ! from OrgChCut (revu 21feb91, 3jun91, 19jun91, 25jun92, 4jan93)
            ! for reuse in OprPointwise divided in subroutines 11jul94,
            ! restructured, simplified, new options 10jun96
            ! (first time I use ENTRY, to avoid COMMON blocks)

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      PARAMETER    (MRange=40)

      CHARACTER     Fehler*(*), LisDoc*(*), Action*(*)
      REAL*8        X(*), Y(*), D(*), ZK(MK), X1(MC), Y1(MC),
     *              X2(MC), Y2(MC), D2(MC), RRange(MRange)
      LOGICAL       qList(MC), qLisR(MC)
      INTEGER       iShift(MK), NZK(MK), K2K(MK), I2Z(MC), J2J(MF)
      CHARACTER     act*20, aus*80, ein*80, fil2*40
      CHARACTER*320 LisC, LisJ, LisCofK, txtRange

      DATA          LisC /' '/, LisCofK /' '/, iModX /1/, n00 /0/,
     *              tolXZ /1.d-5/, qFrom1 /.true./, qInPerJ /.false./

      ENTRY SelectChSet (Action, j, l, LisDoc, Fehler)
C     ------------------------------------------------
         ! Selection for file j (which is presently number l).

      act = Action

      nK = iOlfG (j, '#spectra', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      IF     (l.eq.1) THEN
         ! Questions for all files :
         Print *, ' Select channels: range from'
         Print *,
     * '    (1) channel number    (2) channel number backwards'
         Print *, '    (3) x-value           (4) y-value'
         Print *, '    (5) integral file     (6) data file''s x-range'
         Print *, '    (7) local maxima      (8) local minima'
         iModX = iAskDMu (' Option', iModX, 0, 8)
         IF (iModX.eq.0) THEN
            Fehler = ' '
            RETURN
            ENDIF

      ELSEIF (l.eq.2) THEN
         ! Make first setting global ?
         IF (qiinside(iModX,1,4)) THEN
            IF (qAskPerK) THEN
               qAskPerJ = .true.
            ELSE
               qInPerJ  = .not.qAskD (
     *          ' The same selection for all files', intq(.not.qInPerJ))
               qAskPerJ = qInPerJ
               ENDIF
         ELSE
            qAskPerJ = .true. ! J2J funktioniert noch nicht .false.
            ENDIF
         ENDIF

      IF (l.eq.1 .or. qAskPerJ) THEN
         ! Select on first call (l=1) or individually for each file :

         IF     (nK.eq.1) THEN ! there is just one spectrum
            qAskPerK = .false.
         ELSEIF (qiinside(iModX,1,4)) THEN
            aus = ' The same selection for all spectra'
            qAskPerKin = .not. qAskD (aus, 1-intq(qAskPerKin))
            qAskPerK   = qAskPerKin
            ENDIF

         ENDIF

      IF ((l.eq.1 .or. qAskPerJ) .and. .not.qAskPerK) THEN
         ! Global selection :

         IF (qiinside(iModX,1,2)) THEN

            CALL OlfGetNofK (j, nK, NZK, Fehler)
            CALL iMinMax (NZK, nK, nCmi, nCma, idumi, iduma)
            IF (nCmi.eq.nCma) THEN
               CALL Say2 (' each spectrum has '//cl4(nCmi),' channels')
            ELSE
               CALL Say3 (' the spectra have between '//cl4(nCmi),
     *                    ' and '//cl4(nCma), ' channels')
               ENDIF
            CALL Compose2 (aus, act, ' which channels')
            IF (iModX.eq.2) CALL Append (aus, ' (counting backwards)')
            CALL FrageNList (aus, LisC, nCmi)
            IF     (LisC.eq.'-') THEN
               Fehler = ' '
               RETURN
               ENDIF

            LisDoc = LisC

         ELSEIF ((iModX.eq.3.or.iModX.eq.4) .and. .not.qAskPerK) THEN
            CALL Compose2 (aus, act, ' ranges')
            CALL rAskOnOff (aus, txtRange, RRange, MRange, nRange)
            LisDoc = txtRange

         ELSEIF (iModX.eq.5) THEN
            CALL OlfGet1ZofK (j, 1, nK, ZK, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ! get 2nd file :
            aus = ' Take x-range from Y-data of file'
            j2 = iAskDMu (aus, j2, 0, MF)
            IF (j2.eq.0) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ! common scale X2 must equal scale ZK of 1st file :
            CALL CommonScale (j2, .false., qCom2, nC2,
     *                        qCom2X, nX2, X2, I2Z, Fehler)
            IF     (Fehler.ne.'&ff') THEN
               Print *, 'mcd/ 2nd file/ error in CommonScale'
               RETURN
            ELSEIF (.not.(qCom2 .and. qCom2X .and. nC2.eq.nX2)) THEN
               Fehler = '2nd file has no rectangular grid'
               ENDIF

            aus = 'x-scale of 2nd file -> z-scale of 1st file'
            CALL GetIndex (aus, ZK, nK, X2, nC2, I2Z, tolXZ, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ! inside or outside given range ?
            nK2 = iOlfG (j2, '#spectra', Fehler)
            CALL tOlfG (j2, 'fil', fil2, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            CALL Say2 (' 2nd file contains '//cl4(nK2), ' spectra')
            CALL Compose2 (aus,act,' #1..#2, #3.. (else : ..#1,#2..)')
            qFrom1 = qAskD (aus, intq(qFrom1))
            IF (qFrom1) THEN
               LisDoc = 'inside '//fil2
            ELSE
               LisDoc = 'outside '//fil2
               ENDIF

         ELSEIF (iModX.eq.6) THEN

            aus = ' Take x-range from x-range of file'
            j2 = iAskDMu (aus, j2, 0, MF)
            IF (j2.eq.0) THEN
               Fehler = ' '
               RETURN
               ENDIF

            CALL GetK2K (j, j2, nK1, K2K, Fehler)
            CALL tOlfG (j2, 'fil', fil2, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            LisDoc = 'inside x-range of '//fil2

            ENDIF ! iModX

         ENDIF! global selection

      IF (qAskPerK) THEN
         LisDoc = 'as given per spectrum'
         ENDIF

      RETURN ! SelectChSet

      ENTRY SelectChGet (K, nC, X, Y, D, qList, Fehler)
C     -------------------------------------------------
         ! Get qList(1..nC) for a given spectrum

C  Input required ?
      IF (qAskPerK) THEN

         IF     (qiinside(iModX,1,2)) THEN
            CALL Say3 (' spectrum '//cl4(K), ' contains '//
     *                 cl4(nC), ' channels')
            CALL Compose2 (aus, act, ' which channels')
            IF (iModX.eq.2) CALL Append (aus, ' (counting backwards)')
            CALL FrageNList (aus, LisC, nC)
            IF     (LisC.eq.'-') THEN
               Fehler = ' '
               RETURN
               ENDIF

         ELSEIF (qiinside(iModX,3,4)) THEN
            CALL Compose3 (aus, ' Spectrum '//cl4(K), ' : '//
     *                     act, ' ranges')
            CALL rAskOnOff (aus, txtRange, RRange, MRange, nRange)

            ENDIF

         ENDIF

C  Translate to qList :

      IF     (iModX.eq.1) THEN
         CALL DecNList (LisC, qList, nC, Fehler)

      ELSEIF (iModX.eq.2) THEN
         CALL DecNList (LisC, qLisR, nC, Fehler)
         DO i = 1, nC
            qList(i) = qLisR(nC+1-i)
            ENDDO

      ELSEIF (iModX.eq.3) THEN
         DO i = 1, nC
            qList(i) = .false.
            DO ir = 1, nRange, 2
               IF (qrinside(X(i), RRange(ir), RRange(ir+1)))
     *             qList(i) = .true.
               ENDDO
            ENDDO

      ELSEIF (iModX.eq.4) THEN
         DO i = 1, nC
            qList(i) = .false.
            DO ir = 1, nRange, 2
               IF (qrinside(Y(i), RRange(ir), RRange(ir+1)))
     *             qList(i) = .true.
               ENDDO
            ENDDO

      ELSEIF (iModX.eq.5) THEN
         ! read spectra of 2nd file :
         DO K2 = 1, nK2
            CALL OlfGetXY (j2, K2, n2, X2, Y2, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            Y1(K2) = Y2(I2Z(K))
            ENDDO
         IF (irSorted(Y1,nK2).lt.1) THEN
            CALL Compose2 (Fehler,
     *         'ranges for spectrum '//cl4(K), ' not sorted')
            RETURN
            ENDIF
         ! set qList :
         qSet = .not.qFrom1
         K2 = 1
         ii = 1
 32      CONTINUE
         DO i = ii, nC
           IF (X(i).ge.Y1(K2)) GOTO 321
           qList(i) = qSet
           ENDDO
 321     CONTINUE
         ii = i
         qSet = .not. qSet
         K2 = K2 + 1
         IF (K2.le.nK2) GOTO 32
         DO i = ii, nC
            qList(i) = qSet
            ENDDO

      ELSEIF (iModX.eq.6) THEN

         CALL OlfGetXY (j2, K2K(K), n2, X2, Y2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (n2.lt.2) THEN
            Fehler = ' less than two channels in 2nd file'
            RETURN
            ENDIF
         DO i = 1, nC
            qList(i) = qrinside (X(i), X2(1), X2(n2))
            ENDDO

      ELSEIF (iModX.eq.7) THEN ! local maxima

         IF (nC.lt.2) THEN
            Fehler = 'not enough channels'
            RETURN
            ENDIF

         qList(1) = Y(1).gt.Y(2)
         qList(nC) = Y(nC).gt.Y(nC-1)
         DO i = 2, nC-1
            qList(i) = Y(i).gt.Y(i-1) .and. Y(i).gt.Y(i+1)
            ENDDO

      ELSEIF (iModX.eq.8) THEN ! local minima

         IF (nC.lt.2) THEN
            Fehler = 'not enough channels'
            RETURN
            ENDIF

         qList(1) = Y(1).lt.Y(2)
         qList(nC) = Y(nC).lt.Y(nC-1)
         DO i = 2, nC-1
            qList(i) = Y(i).lt.Y(i-1) .and. Y(i).lt.Y(i+1)
            ENDDO

      ELSE

         CALL Absturz ('SelectChGet', 'iModX oor')

         ENDIF

      RETURN ! SelectChGet

      END ! SelectCh*

C  ====================================================================
C  i30 / 2 :   Auxiliary Operations
C  ====================================================================

      SUBROUTINE SortChannels (X, Y, D, n, qAver)
C     -------------------------------------------
         ! JWu 9feb91, 11mar91, 20mar91, 20oct94
      ! Sorts one spectrum by ascending X;
      ! qAver : if some X are equal, the corresponding Y will be summed.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE 'i_dim.f'
      DIMENSION         X(*), Y(*), D(*), iRang(MC)

C  Checks :
      IF (n.gt.MC) CALL Absturz ('SortChannels',
     *               'Auxiliary array too short')

C  Sort spectrum :
      ifail = 1
      CALL M01DAF_local (X, 1, n, 'a', iRang, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) THEN
         Print *, ' ifail: ', ifail
         Print *, ' X : ', X(1), ' ... ', X(n)
         CALL Absturz ('SortChannels/ M01DAF', 'NAG error')
         ENDIF
      ifail = 1
      CALL M01EAF_local (X, 1, n, iRang, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) CALL Absturz ('SortChannels/ M01EAF(X)',
     * 'NAG error')
      ifail = 1
      CALL M01EAF_local (Y, 1, n, iRang, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) CALL Absturz ('SortChannels/ M01EAF(Y)',
     * 'NAG error')
      ifail = 1
      CALL M01EAF_local (D, 1, n, iRang, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) CALL Absturz ('SortChannels/ M01EAF(D)',
     * 'NAG error')

C  Average over equal X :
      IF (.not.qAver) RETURN
      i = 1
 5601 ContinuE ! DO WHILE (UNIX)
      IF (i.lt.n) THEN
C  - for every i, look whether there are some j with X(j) = X(i) :
         j = i
 5602    ContinuE ! DO WHILE (UNIX)
         IF (j.lt.n .and. qEqEps(X(i), X(j+1)) ) THEN
            j = j + 1
            GOTO 5602
            ENDIF
C  - now j is the highest index with X(j) = X(i)
         IF (j.gt.i) THEN
C  - sum from i to j :
            yy = 0.
            dd = 0.
            DO ii = i, j
               yy = yy + Y(ii)
               dd = dd + D(ii)**2
               ENDDO
C  - store the average as channel i, suppress channels i+1,..,j :
            Y(i) = yy / (j+1-i)
            D(i) = dsqrt0(dd) / (j+1-i)
            n = n - (j - i)
            DO ii = i+1, n
               X(ii) = X(ii+(j-i))
               Y(ii) = Y(ii+(j-i))
               D(ii) = D(ii+(j-i))
               ENDDO
            ENDIF
C  - end of loop, next i :
         i = i + 1
         GOTO 5601
         ENDIF ! i

      END ! SortChannels

C  ====================================================================
C  i30 / 3 :   Auxiliary Calculations
C  ====================================================================

      SUBROUTINE IdaMinMax (qMin, Mod, nGru, X, Y, D, n,
     *                      i0, x0, s0, y0, d0, Fehler)
C     --------------------------------------------------
            ! JWu 6jul91, for use in OprIntegral.
         ! Calculates minimum (if qMin) or maximum of a spectrum (X, Y, D, n).
         ! Available modes are :
         ! Mod = 1 : the smallest/largest channel,
         !       2 : the average of the nGru smallest/largest channels,
         !       3 : the barycenter of the peak.
         ! The channel number i0 is calculated only for Mod=1.
         ! For Mod=3, s0 is the statistical scatter
         ! of the contributing channels.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE 'i_dim.f'
      CHARACTER         Fehler*(*)
      DIMENSION         X(MC), Y(MC), D(MC), iRang(MC)

C  Map minimum -> maximum problem :
      IF (qMin) THEN
         DO i = 1, n
            Y(i) = -Y(i)
            ENDDO
         ENDIF

      IF     (Mod.eq.1) THEN
C  Simply the maximum :
         CALL rMinMax (Y, n, dummy, y0, idummy, i0)
         x0 = X(i0)
         d0 = D(i0)

      ELSEIF (Mod.eq.2) THEN
C  Average over a group :
         IF (nGru.le.0 .or. nGru.ge.n) THEN
            Fehler =
     * 'spectrum contains less data points than your group'
            RETURN
            ENDIF
         ! ranking of Y (in DEscending order)
         ifail = 0 ! hard exit
         CALL M01DAF_local (Y, 1, n, 'd', iRang, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
         ! invert the permutation iRang :
         ifail = 0 ! hard exit
         CALL M01ZAF_local (iRang, 1, n, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
         ! now iRang(1),..,iRang(nGru) contains the indices of the largest Y
         i0 = 0
         x0 = 0.
         s0 = 0.
         y0 = 0.
         d0 = 0.
         DO i = 1, nGru
            ii = iRang(i)
            x0 = x0 + X(ii)
            s0 = s0 + X(ii)**2
            y0 = y0 + Y(ii)
            d0 = d0 + D(ii)
            ENDDO
         x0 = x0 / nGru
         s0 = dsqrt0 ( s0/nGru - x0**2 )
         y0 = y0 / nGru
         d0 = dsqrt(d0) / nGru

      ELSEIF (Mod.eq.3) THEN
C  Barycenter of Peak :
            ! M.Bee in INX_MibIn, JWu 28jul91
         ! determine peak position i0 :
         CALL rMinMax (Y, n, dummy, yPeak, idummy, i0)
         ! determine half width of the peak :
         DO ili = i0, 1, -1
            IF (Y(ili).lt..5*yPeak) GOTO 131
            ENDDO
 131     CONTINUE
         DO iri = i0, n
            IF (Y(iri).lt..5*yPeak) GOTO 132
            ENDDO
 132     CONTINUE
         ! average over peak area :
         x0 = 0.
         f0 = 0.
         DO i = ili, iri
            x0 = x0 + X(i)*Y(i)
            f0 = f0 + Y(i)
            ENDDO
         x0 = dquot0 (x0, f0)
         s0 = 0. ! no idea how to define it
         ! y0 is simply the maximum channel :
         y0 = yPeak
         d0 = 0.

      ELSE
         CALL Absturz ('IdaMinMax', 'Mod o.o.r.')
         ENDIF

C  Reverse the mapping min -> max :
      IF (qMin) THEN
         DO i = 1, n
            Y(i) = -Y(i)
            ENDDO
         y0 = -y0
         ENDIF

      END ! IdaMinMax

      SUBROUTINE Integrate (X, Y, D, n, mod, xL, xH, rInt,dInt,Fehler)
C     ------------------------------------------------------------------
         ! revu 31mai91 JWu

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER  Fehler*(*)
      INCLUDE 'i_dim.f'
      DIMENSION  X(*), Y(*), D(*), Xi(MC), Yi(MC), Di(MC)

         ! integrate the spectrum X(Y) from xL to xH.
         ! result in rInt, error calculated from D(X) in dInt.
         ! mod : 2 = polygon, 4 = four-point finite-differences
         ! special cases : xL=xH=0 means full range.

C  Check that file is sorted :
      IF (irSorted(X,n).ne.2) THEN
         Fehler = 'Integration routine/ file is not sorted'
         RETURN
         ENDIF

C  Copy into internal arrays :
      IF (n.gt.MC) CALL Absturz ('Integrate', 'n>MC')
      DO i = 1, n
         Xi(i) = X(i)
         Yi(i) = Y(i)
         Di(i) = D(i)
         ENDDO

C  Determine first and last channel to be taken into account :
      IF (xL.eq.0. .and. xH.eq.0.) THEN
         ! full range :
         iL = 1
         iH = n
      ELSEIF (xH.le.xL) THEN
         Fehler = 'integration range <= 0'
         RETURN
      ELSE
         iL = irPos (X, n, xL, 'l')  ! X(iL-1) < xL <= X(iL)
         iH = irPos (X, n, xH, 'r')  ! X(iH) <= xH < X(iH+1)
         IF (xL.lt.X(iL) .and. iL.gt.1) THEN
            ! linear interpolation --> Y(xL)
            iL = iL - 1
            Xi(iL) = xL
            CALL LinIntPol (xL, X(iL), X(iL+1), Y(iL), Y(iL+1),
     *                      D(iL), D(iL+1), Yi(iL), Di(iL))
            ENDIF
         IF (xH.gt.X(iH) .and. iH.lt.n) THEN
            ! dito --> Y(xH)
            iH = iH + 1
            Xi(iH) = xH
            CALL LinIntPol (xH, X(iH-1), X(iH), Y(iH-1), Y(iH),
     *                      D(iH-1), D(iH), Yi(iH), Di(iH))
            ENDIF
         ENDIF
      IF (iH.le.iL) THEN
         Fehler = 'range too small or negative'
         RETURN
         ENDIF

C  Now the integration :

      IF (mod.eq.2) THEN
C  - Polygon interpolation :
         stepL = Xi(iL+1)-Xi(iL)
         stepH = Xi(iH)-Xi(iH-1)
         rInt  = stepL    * Yi(iL)     + stepH    * Yi(iH)
         dInt  = stepL**2 * Di(iL)**2  + stepH**2 * Di(iH)**2
         DO i = iL+1, iH-1
            step = Xi(i+1)-Xi(i-1)
            rInt = rInt + step    * Yi(i)
            dInt = dInt + step**2 * Di(i)**2
            ENDDO
         rInt = .5 * rInt
         dInt = .5 * dsqrt(dInt)
      ELSEIF (mod.eq.4) THEN
C  - Four-point finite-differences :
c        Fehler = 'Nag routine D01GAF not available'   ! LOCAL
c        IF (0.eq.0) RETURN
         ! note the completely different definition of dInt :
         ! Di is not used here !!
         IF (iH.lt.iL+3) THEN
            Fehler = 'range too small or negative'
            RETURN
            ENDIF
         ifail = -1 ! soft exit
         CALL D01GAF_local (Xi(iL), Yi(iL), iH-iL+1, rInt, dInt, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
         dInt = dabs(dInt)
         IF (ifail.ne.0) THEN
            Fehler = 'error in NAG D01GAF'
            RETURN
            ENDIF
      ELSE
         CALL Absturz ('Integrate', 'mod o.o.r.')
         ENDIF

      END ! Integrate

      SUBROUTINE RedistrHistogr (M1, n1, X1, A1, Y1, D1,
     *                           M2, n2, X2, A2, Y2, D2, jFill, Fehler)
C     -----------------------------------------------------------------
            ! NY 17nov92, as subroutine 15mar93.

         ! redistribute histogram intensity.
         ! Input histogram :   X1,Y1,D1 (1..n1<=M1)
         ! Output histogram :  X2,Y2,D2 (1..n2<=M2)
         ! Workspace :         A1,A2

      IMPLICIT NONE

      CHARACTER*(*)     Fehler
      REAL*8            X1(*), A1(*), Y1(*), D1(*),   ! M1
     *                  X2(*), A2(*), Y2(*), D2(*),   ! M2
     *                  dx3, dx4, dsqrt0
      INTEGER           jFill, M1, n1, M2, n2, ii, if, i3, i, i4,
     *                  irPosOpt, ii0, ii2, if0, if2

      CALL CoordBins (M1, n1, X1, A1, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'RedistrHistogr/ input grid/ ')
         RETURN
         ENDIF
      CALL CoordBins (M2, n2, X2, A2, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'RedistrHistogr/ output grid/ ')
         RETURN
         ENDIF

      ! restrict A2 to ii..if+1 fully covered by A1 :
      ii = irPosOpt (A2, n2+1, A1(1), 'l', ii)
      if = irPosOpt (A2, n2+1, A1(n1),'r', if) - 1
      IF (if-ii.lt.1) THEN
         Fehler = ' No overlap of new and old grid'
         RETURN
         ENDIF
      ! redistribute intensity :
      i3 = 1
      DO i = ii, if ! outer loop new grid
         i4 = i - ii + 1
            X2(i4) = X2(i)  ! shift x : channels 1..ii-1 not used
            Y2(i4) = 0.
            D2(i4) = 0.
            dx4 = A2(i+1) - A2(i)
 461     CONTINUE   ! inner loop old grid
         IF     (i3.gt.n1+1) THEN
            Fehler = 'PROGRAM ERROR/ i3>n1+1 should never occur'
            RETURN
         ELSEIF (A1(i3+1).le.A2(i)) THEN
            i3 = i3 + 1
            GOTO 461
         ELSEIF (A1(i3).ge.A2(i+1)) THEN
            i3 = max0 (1, i3-1) ! for next i4
            GOTO 462 ! legal exit
         ELSE
            dx3 = dmin1(A2(i+1),A1(i3+1)) - dmax1(A2(i),A1(i3))
            IF (dx3.le.0.d0) THEN
               ! this time consuming check should be supressed
               ! sooner or later. - JWu 30mai94
               Fehler = 'PROGRAM ERROR/ dx3 <= 0.d0'
               RETURN
               ENDIF
            Y2(i4) = Y2(i4) +  dx3/dx4 * Y1(i3)
            D2(i4) = D2(i4) + (dx3/dx4 * D1(i3))**2 ! corr. 8mar93
            i3 = i3 + 1
            GOTO 461
            ENDIF
 462     CONTINUE
         D2(i4) = dsqrt0(D2(i4))
         ENDDO
      n2 = if - ii + 1

      IF     (jFill.eq.3) THEN
         ! accept all points
      ELSEIF (jFill.eq.2) THEN
         ! delete some X2-points if there are more than
         ! two X2-points between a pair of X1-points [11may00]
         DO i4 = 3, n2
            ii0 = irPosOpt (X1, n1, X2(i4),   'l', ii0)
            ii2 = irPosOpt (X1, n1, X2(i4-2), 'l', ii2)
            if0 = irPosOpt (X1, n1, X2(i4),   'r', if0)
            if2 = irPosOpt (X1, n1, X2(i4-2), 'r', if2)
            IF (ii0.eq.ii2 .and. if0.eq.if2) THEN
               n2 = n2 - 1
               X2(i4-1) = X2(i4)
               Y2(i4-1) = Y2(i4)
               D2(i4-1) = D2(i4)
               ENDIF
            ENDDO
      ELSE
         Fehler = 'this fill option not implemented'
         RETURN
         ENDIF

      END ! RedistrHistogr

      SUBROUTINE CoordBins (M, n, X, A, Fehler)
C     -----------------------------------------

      REAL*8        X(*), A(*) ! M
      CHARACTER*(*) Fehler

      IF     (n.gt.M) THEN
         Fehler = 'SEVERE PROGRAM ERROR: n>M'
         RETURN
      ELSEIF (n.eq.M) THEN
         Fehler = ' spectrum has maximum length - shorten by one'
         RETURN
      ELSEIF (n.le.0) THEN
         Fehler = 'SEVERE PROGRAM ERROR: n<=M'
         RETURN
      ELSEIF (n.lt.3) THEN
         Fehler = ' Input spectrum too short'
         RETURN
         ENDIF

      ! set grid A = coordinate ranges of X
      A(1)   = X(1) - (X(2)-X(1)  )/2
      A(n+1) = X(n) + (X(n)-X(n-1))/2
      DO i = 2, n
         A(i) = (X(i) + X(i-1))/2
         IF (A(i).le.A(i-1)) THEN
            Fehler = ' grid not in ascending order'
            RETURN
            ENDIF
         ENDDO

      END ! CoordBins

      SUBROUTINE CalcQQ (j, K, yQ, n, Fehler)
C     ---------------------------------------
            ! JWu 14nov91 for use in OprSpecial, DOS.
         ! get Q(2Th,w) for file j, spectrum K.

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      REAL*8        X(MC), Y(MC), D(MC), yQ(MC), z, E0, rval, tt, w,
     *              rzxOlfG, rzxOlfGG, yQ_of_W, yQ_of_0, facW
      INTEGER       j, K, i, nK, n
      CHARACTER     Fehler*(*), Test*80, un*40, UnW*40
      LOGICAL       qHasQ, qHas2th, qHasW

      Test = '&ff'
      rval = rzxOlfG (j, K, 1, 'q', un, Test) ! ACHTUNG: nicht kopieren
      qHasQ = (Test.eq.'&ff')
      IF (qHasQ .and. un.ne.'A-1') THEN
         Fehler = 'unknown unit for q: '//un
         RETURN
         ENDIF

      Test = '&ff'
      rval = rzxOlfG (j, K, 1, '2th', un, Test) ! ACHTUNG: nicht kopieren
      qHas2th = (Test.eq.'&ff')
      IF (qHas2th .and. un.ne.' ') THEN
         Fehler = 'nonsensical unit for 2th: '//un
         RETURN
         ENDIF

      Test = '&ff'
      rval = rzxOlfG (j, K, 1, 'w', UnW, Test) ! ACHTUNG: nicht kopieren
      qHasW = (Test.eq.'&ff')
      IF (qHasW) THEN
         CALL UnitConv (UnW, 'meV', facW, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF

      IF (qHasQ .and. qHas2th) THEN
         Fehler = ' Both q and 2th are defined'
         RETURN
      ELSEIF (.not.qHasQ .and. .not.qHas2th) THEN
         Fehler = ' Neither q nor 2th are defined'
         RETURN
         ENDIF

      CALL OlfGetN  (j, K, n, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Get Q :
      DO i = 1, n
         E0 = rzxOlfGG (j, K, i, 'E0', 'meV', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (qHasQ) THEN
            yQ(i) = rzxOlfGG (j, K, i, 'q', 'A-1', Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSE
            tt = rzxOlfGG (j, K, i, '2th', ' ', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (qHasW) THEN
               w = rzxOlfGG (j, K, i, 'w', UnW, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               yQ(i) = yQ_of_W (w*facW, E0, tt)
            ELSE
               yQ(i) = yQ_of_0 (E0, tt)
               ENDIF
            ENDIF
         IF (yQ(i).lt.0) THEN
            Fehler = ' Q < 0'
            RETURN
            ENDIF
         ENDDO

      END ! CalcQQ
