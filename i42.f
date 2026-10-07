C  ====================================================================
C
C      Library  IDA   :  Inelastic data treatment
C      Modul    i42   :     data manipulations / files
C
C  ====================================================================

C     Contents :
C        1. OrgFile
C              Join

C  ====================================================================
C  i42 / 1 :   OrgFile...
C  ====================================================================

C  --------------------------------------------------------------------
      SUBROUTINE tJoin (nJList, JList, jout, lab, longlab, Fehler)
C  --------------------------------------------------------------------
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler, lab, longlab
      CHARACTER*80  aus, ein
      INTEGER       nJList, JList(*), lj, jout, ii
      LOGICAL       qIdent

      IF (nJList.le.0) THEN
         Fehler = 'tJoin/ empty j-list'
         RETURN
         ENDIF

      qIdent = .true.
      CALL tOlfG (JList(1), lab, aus, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'tJoin[1]')
         RETURN
         ENDIF
      DO lj = 2, nJList
         CALL tOlfG (JList(lj), lab, ein, Fehler)
         IF (Fehler.ne.'&ff') THEN
            CALL Insert (Fehler, 1, 'tJoin[2ff]')
            RETURN
            ENDIF
         DO ii = 1, len(aus)
            IF (ein(ii:ii).ne.aus(ii:ii)) THEN
               aus(ii:ii) = ' '
               qIdent = .false.
               ENDIF
            ENDDO
         ENDDO

      IF (.not.qIdent) THEN
         CALL Kontrakt (aus)
         CALL FrageCD ('Joint '//longlab, aus, aus)
         ENDIF

      CALL tOlfP (jout, lab, aus, Fehler)

      END ! tJoin

C  --------------------------------------------------------------------
      SUBROUTINE OrgFileJoin (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
            ! JWu 17apr91, z-scale 10jul91, to Org..(Ida4) and new
            ! option Exchange 5nov91, option sum 29apr92

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Fehler
      CHARACTER*80  aus
      CHARACTER*40  Co, Un, CoZout(MZ), UnZout(MZ), CoPout(MRP), 
     *              UnPout(MRP)
      CHARACTER*20  h1, h2
      REAL*8        z, Zout(MZ), rval, rvalmin, rvalmax, rvalavg, 
     *              rzOlfGG
      INTEGER       nJList, JList(*), lj, jout, NKJ(MF), nK,
     *              nKsum, K, Kout, n, iModJ, iModJdef, iZ, nZ,
     *              iZout, nZout, iP, iPout, nPout, ih1, ih2, iModR
      LOGICAL       qSort, qOv

      DATA          iModJdef /1/, qSort /.true./

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Check that all input files exist; by the way calculate nKsum
      nKsum = 0
      DO lj = 1, nJList
         nK = iOlfG (JList(lj), '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         NKJ(lj) = nK
         nKsum   = nKsum  + nK
         ENDDO

C  Make new file :
      CALL OlfHeadDup (JList(1), .false., jout, nK, Kout, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL tJoin (nJList, JList, jout, 'fil', 'file name', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL tJoin (nJList, JList, jout, 'tit', 'title line', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Eigentlich *alle* Header miteinander vermahlen ...
      CALL OlfComAdd (jout, 'J', ' ', Fehler)

C  Make list of z :
      nZout = 0
      DO lj = 1, nJList
         nZ = iOlfG (JList(lj), '#Z', Fehler)
         DO iZ = 1, nZ
            CALL OlfCnuG (JList(lj), 'z'//cl2(iZ), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            DO iZout = 1, nZout
               IF (CoZout(iZout).eq.Co) THEN
                  IF (UnZout(iZout).ne.Un) THEN
                     Fehler = 'different units in '//Co
                     RETURN
                     ENDIF
                  GOTO 28
                  ENDIF
               ENDDO
            nZout = nZout + 1
            IF (nZout.gt.MZ) THEN
               Fehler = 'too many z'
               RETURN
               ENDIF
            CoZout(nZout) = Co
            UnZout(nZout) = Un
 28         CONTINUE
            ENDDO
         ENDDO

C  Make list of r:
      nPout = 0
      DO lj = 1, nJList
         DO iP = 1, MRP
            Co   = '&pbn '//cv3(iP)
            rval = rOlfG (JList(lj), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (Co.eq.'&empty') GOTO 38
            DO iZout = 1, nZout
               IF (Co.eq.CoZout(iZout)) GOTO 38 ! rP is Z of other file
               ENDDO
            IF (Co.eq.'&eop') GOTO 39
            DO iPout = 1, nPout
               IF (CoPout(iPout).eq.Co) THEN
                  IF (UnPout(iPout).ne.Un) THEN
                     Fehler = 'different units in '//Co
                     RETURN
                     ENDIF
                  GOTO 38
                  ENDIF
               ENDDO
            nPout = nPout + 1
            IF (nPout.gt.MRP) THEN
               Fehler = 'too many r-parameter'
               RETURN
               ENDIF
            CoPout(nPout) = Co
            UnPout(nPout) = Un
 38         CONTINUE
            ENDDO
 39      CONTINUE
         ENDDO ! lj

C  Compare values of r:
      DO iPout = 1, nPout
         Co = CoPout(iPout)
         Un = UnPout(iPout)
         rval = rOlfGG (JList(1), Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         rvalmin = rval
         rvalmax = rval
         rvalavg = rval
         DO lj = 2, nJList
            rval = rOlfGG (JList(lj), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            rvalmin = dmin1 (rvalmin, rval)
            rvalmax = dmax1 (rvalmax, rval)
c            rvalavg = rvalavg + rval
            ENDDO
c         rvalavg = rvalavg / nJList
         IF (rvalmin.ne.rvalmax) THEN
c            CALL NiceNum (rvalmin, h1, ih1)
c            CALL NiceNum (rvalmax, h2, ih2)
c            Print *, '  r-parameter values for '//Co(1:lenU(Co))//
c     *         ' vary from '//h1(1:ih1)//' to '//h2(1:ih2)//' '//Un(1:lenU(Un))
c            iModR = iAskMu (' Average(1), make to z(2), throw away(3)', 0, 4)
c            IF     (iModR.eq.0) THEN
c               Fehler = ' '
c               RETURN
c            ELSEIF (iModR.eq.1) THEN
c               CALL rOlfP (jout, Co, Un, rvalavg, Fehler)
c            ELSEIF (iModR.eq.2) THEN
               CALL rOlfDel (jout, Co, Fehler)
               nZout = nZout + 1
               IF (nZout.gt.MZ) THEN
                  Fehler = 'too many z - eliminate some z or r'
                  RETURN
                  ENDIF
               CoZout(nZout) = Co
               UnZout(nZout) = Un
c            ELSEIF (iModR.eq.3) THEN
c               CALL rOlfDel (jout, Co, Fehler)
c               ENDIF
            IF (Fehler.ne.'&ff') RETURN
         ELSE
            CALL rOlfP (jout, Co, Un, rvalavg, Fehler)
            ENDIF
         ENDDO ! r-par

C  Z's to header:
         DO iZout = 1, nZout
            CALL OlfCnuP (jout, 'z'//cl2(iZout),
     *                    CoZout(iZout), UnZout(iZout), Fehler)
            ENDDO
         nZout = nZout + 1
         IF (nZout.gt.MZ) THEN
            Fehler = 'too many z - eliminate some z or r'
            RETURN
            ENDIF
         CALL OlfCnuP (jout, 'z'//cl2(nZout), 'no.', ' ', Fehler)

C  Main loop: join spectra
         DO lj = 1, nJList
            DO K = 1, NKJ(lj)
               Kout = Kout + 1
               DO iZout = 1, nZout-1
                  Zout(iZout) = rzOlfGG (JList(lj), K, CoZout(iZout),
     *                                   UnZout(iZout), Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDDO
               Zout(nZout) = lj
               CALL OlfPutZ (jout, Kout, nZout, Zout, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               CALL OlfCopXYD (JList(lj), jout, K, Kout, Fehler)
               IF (Fehler.ne.'&ff') RETURN

               ENDDO ! K
            ENDDO ! lj

      CALL OlfClos (jout, Kout, Fehler)

      ! Delete input files ?
      IF (qOv) CALL FileKill (nJList, JList, Fehler)

      END ! OrgFileJoin
