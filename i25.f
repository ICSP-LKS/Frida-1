C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i25   :     file delete, copy, make
C
C  ====================================================================

C     Contents :
C        1.  FileKill, FileCopy
C        2.  FileMake
C        3.  SaveArray

C  ====================================================================
C  i25 / 1 :   file delete, copy
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE FileKill (nJList, JList, Fehler)
C     -------------------------------------------
            ! JWu 30jul91
         ! Delete some files

      INCLUDE 'i_dim.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF
      DO lj = nJList, 1, -1
         CALL MemFileDel (JList(lj), Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      END ! FileKill

      SUBROUTINE FileCopy (nJList, JList, Fehler)
C     -------------------------------------------
            ! JWu 23oct91. Reduction > MemFileDup 7nov96
         ! Duplicate files in the run-time memory

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF
      DO lj = 1, nJList
         CALL MemFileDup (JList(lj), jout, 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      END ! FileCopy

C  ====================================================================
C  i25 / 2 :   create file from nothing
C  ====================================================================

      SUBROUTINE FileMake (Fehler)
C     ----------------------------
         ! JWu 11apr91, included Any2IED 14sep94
      ! Make a completely new data file and store it in the run-time memory

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER     (MRL=40, MRein=80)
      DIMENSION      qEnt(3), iOM(2), RL(MRL), Rein(MRein)
      CHARACTER*(*)  Fehler
      CHARACTER*80   aus, File, Title
      CHARACTER*80   ein
      CHARACTER*40   CoX, UnX, CoY, UnY, CoZ, UnZ

      DATA          iOx / 1 /, iOz / 1 /, iOM / 1, 3 /,
     *              Title / ' ' /, File / ' ' /

      Print *, ' make your own IED file :'

C  Make header :
      CALL OlfCreate (jout, Kout, File, Title, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL tOlfG (jout, 'fil', File, Fehler)  ! f"ur's n"achste Mal
      CALL tOlfG (jout, 'tit', Title, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL FrageTD (' x Coordinate name', CoX, CoX)
      CALL FrageCD (' Unit', UnX, UnX)
      CALL OlfCnuP (jout, 'x', CoX, UnX, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL FrageTD (' y Coordinate name', CoY, CoY)
      CALL FrageCD (' Unit', UnY, UnY)
      CALL OlfCnuP (jout, 'y', CoY, UnY, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL FrageCD (' z Coordinate name', CoZ, CoZ)
      qZ = (CoZ.ne.' ')
      IF (qZ) THEN
         nZ = 1
         CALL FrageCD (' Unit', UnZ, UnZ)
         CALL OlfCnuP (jout, 'z1', CoZ, UnZ, Fehler)
         IF (Fehler.ne.'&ff') RETURN
      ELSE
         nZ = 0
         Print *, ' no z given - only one spectrum can be read'
         ENDIF

 32   CONTINUE
      Print *, ' Input mode :'
      Print *, '   (1) enter x coordinates first'
      Print *, '   (2) set x on equidistant grid'
      Print *, '   (3) enter x-y pairs'
      Print *, '   (4) enter x-y-d triples'
      Print *, '   (5) select from multicolumn input'
      IF (qZ) Print *, '   (6) dito, including z' ! 8mar97
      iOx = iAskDMu (' Option', iOx, 0, 6)
      qEnt(1) = .false. ! y will be entered along with x
      qEnt(2) = .false. ! d will be entered along with x
      qEnt(3) = .false. ! z will be entered along with x
      IF     (iOx.eq.0) THEN
         Fehler = ' '
         RETURN
      ELSEIF (iOx.eq.3) THEN
         qEnt(1) = .true.
      ELSEIF (iOx.eq.4) THEN
         qEnt(1) = .true.
         qEnt(2) = .true.
         imcx = 1
         imcy = 2
         imcd = 3
      ELSEIF (iOx.eq.5) THEN
         CALL i3FrageD (' Columns of x,y,d (0=no input)',
     *                  imcx, imcy, imcd, imcx, imcy, imcd)
         IF (imcx.le.0 .or. imcy.le.0) THEN
            Print *, ' x and y must be read from multicolumn input'
            CALL Gong (3)
            GOTO 32
            ENDIF
         qEnt(1) = .true.
         qEnt(2) = (imcd.gt.0)
      ELSEIF (iOx.eq.6) THEN
         IF (.not.qZ) THEN
            Fehler = 'inconsistent choices'
            RETURN
            ENDIF
         imcx = iAskD (' x from column no.', imcx)
         imcy = iAskD (' y from column no.', imcy)
         imcd = iAskD (' d from column no. (0=no input)', imcd)
         imcz = iAskD (' z from column no.', imcz)
         IF (imcx.le.0 .or. imcy.le.0 .or. imcz.le.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         qEnt(1) = .true.
         qEnt(2) = (imcd.gt.0)
         qEnt(3) = .true.
         ENDIF

      IF (.not.qEnt(1)) THEN
         Print *, ' Input mode for y :'
         Print *, '   (1) enter individually'
         Print *, '   (2) common value'
         Print *, '   (3) common value 0.0'
         Print *,
     * '   (4) enter, not necessarily with CR LF between data'
         iOpt = iAskDMu (' Option', iOM(1), 0, 4)
         IF (iOpt.eq.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         iOM(1) = iOpt
         ENDIF

      IF (.not.qEnt(2)) THEN
         Print *, ' Input mode for error d :'
         Print *, '   (1) enter individually'
         Print *, '   (2) common value'
         Print *, '   (3) common value 0.0'
         Print *, '   (5) sqrt(y)'
 377     iOpt = iAskDMu (' Option', iOM(2), 0, 5)
         IF (iOpt.eq.4) THEN
            CALL Gong (1)
            GOTO 377
            ENDIF
         IF (iOpt.eq.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         iOM(2) = iOpt
         ENDIF

 38   CONTINUE
      IF (qZ .and. .not.qEnt(3)) THEN
         Print *, ' Input mode for z :'
         Print *, '   (1) for each spectrum, just enter z'
         Print *, '   (2) take z from a header block'
         iOz = iAskDMu (' Option', iOz, 0, 2)
         IF (iOz.eq.0) THEN
            Fehler = ' '
            RETURN
            ENDIF
         IF     (iOz.eq.1) THEN
            iHead = 1
            iLz = 1
            iCz = 1
         ELSEIF (iOz.eq.2) THEN
            Print *,
     *  '   please note: data blocks must be followed by '//
     *  'an empty line - '
            Print *,
     *  '   this line does *not* count as part of the '//
     *  'following header block.'
            Print *,
     *  '   the first line of the header block itself '//
     *  'may *not* be empty -'
            Print *, '   otherwise, no more spectra are read in'

            iHead = iAskDMu (' Number of header lines', iHead, 1, 4096)
            CALL i2FrageD (' Position of z (line, column)',
     *                     iLz, iCz, iLz, iCz)
            IF (iLz.le.0 .or. iLz.gt.iHead .or. iCz.le.0) THEN
               CALL Gong (3)
               GOTO 38
               ENDIF
            ENDIF
         ENDIF

      Print *, '   now everything is ready for reading the data;'
      Print *, '   for reading from a file, use "\i <file-name>".'

      IF (iOx.eq.6) THEN ! separat erledigen
         CALL FileMakeXYZ (jout, imcx, imcy, imcz, imcd,iOM(2),Fehler)
         RETURN
         ENDIF

C  Make spectra :
      DO K = 1, MK

C  Get z :
 41      CONTINUE
         IF (qZ) THEN
            Print *, ' make spectrum no. '//cl6(K)

c            IF     (iOz.eq.1) THEN
c               z = rAsk (' Enter z')
c            ELSEIF (iOz.eq.2) THEN
               DO jh = 1, iHead
                  aus = ' Enter header line'
                  IF (jh.eq.iLz) CALL Append (aus, ' containing z')
                  IF (jh.eq.1)   CALL Append (aus, ' (or " " to stop)')
                  CALL FrageC (aus, ein)
                  IF (jh.eq.1 .and. ein.eq.' ') GOTO 5 ! eoi
                  IF (jh.eq.iLz) THEN
                     CALL FindR (ein, MRL, nRL, RL, .false.)
                     IF (nRL.lt.iCz) THEN
                        CALL Gong(7)
                        Print *,
     * ' cannot find z in that line, setting z=0'
                        z = 0
                     ELSE
                        z = RL(iCz)
                        ENDIF
                     ENDIF
                  ENDDO ! jh
c               ENDIF
            Print *, ' and now the data themselves :'
            ENDIF ! z ?

C  Get x (and more ?) :
         IF (iOx.eq.1) THEN
            CALL rAskArray (X, MC, n)
         ELSEIF (iOx.eq.2) THEN
            CALL rAskGrid  ('x', X, MC, n)
         ELSEIF (iOx.eq.3) THEN
            CALL rAskPairs (X, Y, MC, n)
         ELSEIF (iOx.eq.4 .or. iOx.eq.5) THEN ! multicol input (ex Any2IED)
            Print *, ' Enter data lines; enter " " to stop:'
            n = 0
 45         CONTINUE ! endless loop
               CALL FrageC (' > ', ein)
               IF (ein.eq.' ') GOTO 459 ! end of data block - regular exit
               CALL FindR (ein, MRL, nRL, RL, .false.)
               IF (nRL.ge.imcx) THEN
                  X(n+1) = RL(imcx)
               ELSE
                  Print *, ' No x value !'
                  CALL Gong (1)
                  GOTO 45
                  ENDIF
               IF (nRL.ge.imcy) THEN
                  Y(n+1) = RL(imcy)
               ELSE
                  Print *, ' No y value !'
                  CALL Gong (1)
                  GOTO 45
                  ENDIF
               IF (imcd.gt.0) THEN
                  IF (nRL.ge.imcd) THEN
                     D(n+1) = RL(imcd)
                  ELSE
                     Print *, ' No d value !'
                     CALL Gong (1)
                     GOTO 45
                     ENDIF
                  ENDIF
               n = n + 1
               IF (n.eq.MC) THEN
                  Print *, ' Maximum number of channels reached'
                  CALL Gong (6)
                  GOTO 459
                  ENDIF
               GOTO 45
 459        CONTINUE ! end of data

            ENDIF ! iOx

         IF (n.le.0) THEN
            IF (qZ) THEN
               Print *,
     * ' empty x-range entered -- try it again or stop :'
               CALL Gong (2)
               GOTO 41
            ELSE
               Fehler = ' empty x-range entered'
               RETURN
               ENDIF
            ENDIF

C  Get y and d, if not done before :
         DO m = 1, 2
            IF (qEnt(m)) GOTO 49
            IF (m.eq.1) Print *, ' enter the y-values :'
            IF (m.eq.2) Print *, ' enter the error bars :'

            IF     (iOM(m).eq.1) THEN
               DO i = 1, n
                  write (aus,'(a,i4,a,g12.5,a)')
     *               ' channel', i, ', x = ', X(i), ' --> '
                  IF (m.eq.1) Y(i) = rAsk (aus)
                  IF (m.eq.2) D(i) = rAsk (aus)
                  ENDDO
            ELSEIF (iOM(m).le.3) THEN
               IF (iOM(m).eq.2) val = rAsk (' common value ?')
               IF (iOM(m).eq.3) val = 0.
               IF (m.eq.1) CALL rSet (Y, 1, n, 1, val)
               IF (m.eq.2) CALL rSet (D, 1, n, 1, val)
            ELSEIF (iOM(m).eq.4) THEN
               Print *, ' now feed in the data :'
               i = 0
               DO WHILE (i.lt.n)
                  CALL FrageC ('.. ', ein)
                  CALL FindR (ein, MRein, nRein, Rein, .false.)
                  DO ii = 1, nRein
                     Y(i+1) = Rein(ii)
                     i = i + 1
                     ENDDO
                  ENDDO
               IF (i.gt.n) Print *,
     * ' too many data given > will be ignored'
            ELSEIF (iOM(m).eq.5) THEN
               IF (m.ne.2) CALL Absturz ('FileMake',
     * 'Option 5 only for D')
               qWarn = .false.
               DO i = 1, n
                  IF (Y(i).lt.0) THEN
                     qWarn = .true.
                     D(i) = 0
                  ELSE
                     D(i) = sqrt(Y(i))
                     ENDIF
                  ENDDO
               IF (qWarn) THEN
                  CALL Gong (7)
                  Print *, 'WARNING/ some Y(i)<0  ==> cannot set error'
                  ENDIF
               ENDIF

 49         CONTINUE
            ENDDO ! m : Y..D

C  Done - save spectrum :
         Kout = Kout + 1
         CALL OlfPutSpe (jout, Kout, nZ, z, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (.not.qZ) GOTO 5 ! only one spectrum
         ENDDO ! spectra

      CALL Gong (3)
      Print *, ' maximum number of spectra created'
 5    CONTINUE

      Print *, '   new file completed'
      Print *, '   to modify parameters, use "ed", "er", ...'
      Print *, '   to store the file on disk, use "fs"'
      CALL OlfClos (jout, Kout, Fehler)

      END ! FileMake

      SUBROUTINE FileMakeXYZ (jout, imcx, imcy, imcz, imcd,iOd,Fehler)
C     ------------------------------------------------------------------
            ! 8mar97 (f"ur TS-Scans)

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER     (MRL=40, MRein=80)
      DIMENSION      RL(MRL), Rein(MRein)
      CHARACTER*(*)  Fehler
      CHARACTER*240  ein

      Print *, ' Enter data lines; enter " " to stop:'
      nK = 0

 1    CONTINUE                  ! endless loop

         CALL FrageC (' > ', ein)
         IF (ein.eq.' ') GOTO 9 ! end of data block - regular exit
         CALL FindR (ein, MRL, nRL, RL, .false.)

         IF (nRL.lt.imcx .or. nRL.lt.imcy .or.
     * (imcd.gt.0 .and. nRL.lt.imcd) .or. nRL.lt.imcz) THEN
            Print *, ' Not enough columns given'
            CALL Gong (1)
            GOTO 1
            ENDIF

         xn = RL(imcx)
         yn = RL(imcy)
         zn = RL(imcz)
         IF (imcd.gt.0) THEN
            dn = RL(imcd)
         ELSEIF (iOd.eq.3) THEN
            dn = 0
         ELSE
            Fehler =
     * 'implementing this d-mode is trivial but not yet done'
            RETURN
            ENDIF

         ! compare zn to the z given so far
         qNewZ = .false.
         IF (nK.ge.1) THEN
            ENDIF

         IF (qNewZ) THEN
            ! new z -> new spectrum

         ELSE
            ! existing z -> new entry


            n = n + 1
            IF (n.eq.MC) THEN
               Print *, ' Maximum number of channels reached'
               CALL Gong (6)
               GOTO 9
               ENDIF
            ENDIF

         GOTO 1

 9    CONTINUE
      Print *, '   new file completed'
      Print *, '   to modify parameters, use "dt", "dr", ...'
      Print *, '   to store the file on disk, use "fw"'
      CALL OlfClos (jout, Kout, Fehler)

      END ! FileMakeXYZ

      SUBROUTINE FileMakeHist (Fehler)
C     --------------------------------
            ! JWu 28jan00
         ! Make a completely new histogram file from single-events log
      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER     (MRL=10, MRein=20)
      REAL*8         RL(MRL), Rein(MRein), Z1(MK)
      CHARACTER*(*)  Fehler
      CHARACTER*80   File, Title, ein, doc
      CHARACTER*40   CoX, UnX, CoY, UnY, CoZ, UnZ

      DATA           Title / ' ' /, File / ' ' /


C  Make header :
      CALL OlfCreate (jout, Kout, File, Title, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL tOlfG (jout, 'fil', File, Fehler)  ! f"ur's n"achste Mal
      CALL tOlfG (jout, 'tit', Title, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL FrageTD (' x Coordinate name', CoX, CoX)
      CALL FrageCD (' Unit', UnX, UnX)
      CALL OlfCnuP (jout, 'x', CoX, UnX, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CoY = '#['//CoX(1:lenU(CoX))//']'
      IF (UnX.ne.' ') THEN
         UnY = UnX(1:lenU(UnX))//'^-1'
      ELSE
         UnY = ' '
         ENDIF
      CALL OlfCnuP (jout, 'y', CoY, UnY, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL FrageCD (' z Coordinate name', CoZ, CoZ)
      qZ = (CoZ.ne.' ')
      IF (qZ) THEN
         nZ = 1
         CALL FrageCD (' Unit', UnZ, UnZ)
         CALL OlfCnuP (jout, 'z1', CoZ, UnZ, Fehler)
         IF (Fehler.ne.'&ff') RETURN
      ELSE
         nZ = 0
         Print *, ' no z given - only one spectrum can be read'
         ENDIF

      IF (qZ) THEN

         CALL i3FrageD (' x,weight,z', ix, iw, iz, ix, iw, iz)
         CALL i3FrageD (' Columns of x,weight,z', imcx, imcw, imcz,
     *                  imcx, imcw, imcz)
         nK = 0
      ELSE
         CALL i2FrageD (' Columns of x,weight', imcx, imcw, imcx, imcw)
         nK = 1
         ENDIF

      CALL SetGridChoice (.false., doc, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL SetGridJ (j, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL SetGridK (1, n, Fehler) ! -> X1(1..n)
      IF (Fehler.ne.'&ff') RETURN

C  Reset arrays :
      DO i = 1, n
         D(i) = 0
         DO K = 1, MK
            Wrk3dim(i,K,1) = 0
            ENDDO
         ENDDO
      totweight = 0

C  Get input line (endless loop) :
 41      CONTINUE

         CALL FrageC (' > ', ein)
         IF (ein.eq.' ') GOTO 49 ! end of data block - regular exit
         CALL FindR (ein, MRL, nRL, RL, .false.)

         IF (nRL.ge.imcx) THEN
            val = RL(imcx)
         ELSE
            Print *, ' No x value !'
            CALL Gong (1)
            GOTO 49
            ENDIF
         ix = irPosOpt (X1, n, val, 'n', ix)

         IF (.not.qZ) THEN
            K = 1
         ELSEIF (nRL.ge.imcz) THEN
            zin = RL(imcz)
            DO K = 1, nK
               IF (zin.eq.Z1(K)) GOTO 4219
               ENDDO
            ! new value of z => new spectrum
            nK = nK + 1
            IF (nK.gt.MK) THEN
               Fehler = 'too many different z'
               RETURN
               ENDIF
            K = nK
            Z1(K) = zin
 4219       CONTINUE
         ELSE
            Print *, ' No z value !'
            CALL Gong (1)
            GOTO 49
            ENDIF

         IF (imcw.eq.0) THEN
            weight = 1
         ELSEIF (nRL.ge.imcw) THEN
            weight = RL(imcw)
         ELSE
            Print *, ' No weight !'
            CALL Gong (1)
            GOTO 49
            ENDIF
         IF (weight.lt.0) THEN
            Fehler = 'negative weight'
            RETURN
            ENDIF

         Wrk3dim(ix,K,1) = Wrk3dim(ix,K,1) + weight
         totweight = totweight + weight

         GOTO 41 ! endless loop
 49   CONTINUE

      IF (totweight.le.0) THEN
         Fehler = 'total weight not > 0'
         RETURN
         ENDIF

C  Done - save spectra :
      DO K = 1, nK
         ! normalise :
         DO i = 1, n
            Wrk3dim(i,K,1) = Wrk3dim(i,K,1) / totweight
            ENDDO

         CALL OlfPutSpe (jout, K, 1, Z1(K), n, X1,
     *                   Wrk3dim(1,K,1), D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO ! spectra

      CALL OlfClos (jout, nK, Fehler)

      END ! FileMake

C  ====================================================================
C  i25 / 3 :   save array
C  ====================================================================

      SUBROUTINE SaveMatrix (jout, Fil, Tit,
     *                       CoX, UnX, CoZ, UnZ, CoY, UnY,
     *                       MCin, MKin, nC, nK, Xin, Zin, Yin, Fehler)
C     -----------------------------------------------------------------
            ! JWu 11jan00 for use in i94/mscat
         ! Save a rectangular array Y(X,Z1) without error bars

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'
      INTEGER       jout, MCin, MKin, MY, nC, nK, Kout, K
      REAL*8        Xin(MCin), Zin(MKin), Yin(MCin,MKin)
      CHARACTER*(*) Fil, Tit, Fehler, CoX, UnX, CoZ, UnZ, CoY, UnY

      IF     (nK.gt.MK) THEN
         CALL Compose2 (Fehler, 'SaveMatrix/ cannot store '//cl6(nK),
     *        'spectra; maximum is '//cl6(MK))
         RETURN
      ELSEIF (nC.gt.MC) THEN
         CALL Compose2 (Fehler, 'SaveMatrix/ cannot store '//cl6(nC),
     *        'channels; maximum is '//cl6(MC))
         RETURN
         ENDIF

      CALL OlfCreate (jout, Kout, '&noask','&noask', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL tOlfP (jout, 'fil', Fil, Fehler)
      CALL tOlfP (jout, 'tit', Tit, Fehler)
      CALL OlfCnuP (jout, 'x',  CoX, UnX, Fehler)
      CALL OlfCnuP (jout, 'y',  CoY, UnY, Fehler)
      CALL OlfCnuP (jout, 'z1', CoZ, UnZ, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO K = 1, nK
         CALL OlfPutZ (jout, K, 1, Zin(K), Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfPutXY0 (jout, K, nC, Xin, Yin(1,K), Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO ! K

      CALL OlfClos (jout, K, Fehler)

      END ! SaveArray
