C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i23   :     directories, editing, r/z-handling
C
C  ====================================================================

C     Contents :
C        1.  Directories and Editing :
C               MemInfoF/Z/K/Y
C        2.  Editing :
C               EditCnu, EditDoc, EditZ, EditRPar, EditIPar, EditGPar
C        3.  Handling of z/r :
C               ..

C  ====================================================================
C  i23 / 1 :   directories
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE MemInfoF (Fehler)
C     ----------------------------
            ! JWu 16may95
         ! directory of on-line files

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER   Fehler*(*), aus*80
      CHARACTER   h1*10, cMod*1, Fil*80, Doc*40, Tit*80,
     *            xCo*40, yCo*40, xUn*40, yUn*40
      DIMENSION   NofK(MK)

      iRest = MemBlockInq ('fE')
      iTot  = MemBlockInq ('ME')

C  Head line :
      write (h1, '(f5.1)') 100 * dble(iRest)/iTot
      CALL DelLeft (h1)
      CALL Say4 (' '//cl8(iRest),' lines free in run-time memory (',
     *           h1,'%)')

C  Table of files :
      DO j = 1, MemBlockInq ('nF')

         CALL MemFileStatG (j, iStat, Fehler)
         IF (Fehler.ne.'&ff') GOTO 80
         IF (qBitGet(iStat,1)) THEN
            cMod = ' '
         ELSE
            cMod = '='
            ENDIF

         CALL tOlfG (j, 'fil', Fil, Fehler)
         IF (Fehler.ne.'&ff') GOTO 80
         CALL tOlfG (j, 'tit', Tit, Fehler)
         IF (Fehler.ne.'&ff') GOTO 80
         CALL tOlfG (j, 'doc', Doc, Fehler)
         IF (Fehler.ne.'&ff') GOTO 80

         CALL OlfCnuG (j, 'x', xCo, xUn, Fehler)
         IF (Fehler.ne.'&ff') GOTO 80
         CALL OlfCnuG (j, 'y', yCo, yUn, Fehler)
         IF (Fehler.ne.'&ff') GOTO 80

         CALL OlfGetNofK (j, nK, NofK, Fehler)
         IF (Fehler.ne.'&ff') GOTO 80
         n = NofK(1)
         DO K = 2, nK
            IF (NofK(K).ne.n) THEN
               n = 0
               GOTO 29
               ENDIF
            ENDDO
 29      CONTINUE

         IF (qBitGet(iStat,2)) Doc = 'read-only'

         idoc = lenU(Doc)
         IF (idoc.lt.11) THEN
            idoc = 11
            ENDIF
         IF (idoc.gt.11) THEN
            Doc(idoc-10:idoc-9) = '..'
            ENDIF
         aus = ' '//cr2(j)//cMod//Fil(1:15)//' '//Doc(idoc-10:idoc)//
     *            ' '//Tit(1:26)//' '//xCo(1:3)//
     *            ' '//yCo(1:6)//' '//cr3(nK)//'*'
         IF (n.eq.0) THEN
            CALL Append (aus, '  ?')
         ELSE
            CALL Append (aus, cr4(n))
            ENDIF
         GOTO 90

 80      CONTINUE ! error occured
         aus = ' '//cr2(j)//' ILLISIBLE: '//Fehler
         Fehler = '&ff'

 90      CONTINUE
         Print '(a80)', aus
         ENDDO

      END ! MemInfoF

      SUBROUTINE MemInfoZ (nJlist, JList, Fehler)
C     -------------------------------------------
            ! JWu 17may95

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER      Fehler*(*), Co*40, Un*40, h1*4, h2*20, aus*80
      INTEGER        nJList, JList(*), lj, j, nK, K, iOlfG, nZ, iZ
      REAL*8         X(MC), Z(MZ)

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      DO lj = 1, nJList
      j = JList(lj)

         nK = iOlfG (j, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         nZ = iOlfG (j, '#Z', Fehler)
         IF (nZ.ge.1) THEN
            ! List header :
            Print '(a)', ' z coordinates of file '//cl3(j)
            IF (Fehler.ne.'&ff') RETURN
            aus = ' '
            DO iZ = 1, nZ
               CALL OlfCnuG (j, 'z'//cl2(iZ), Co, Un, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL Compose3 (aus(4+(iZ-1)*12:4+iZ*12-1), Co, ' ['//
     *                        Un, ']')
               ENDDO
            Print '(a80)', aus

            ! List of spectra :
            DO K = 1, nK
               CALL OlfGetZ (j, K, nZ, Z, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               aus = cr3(K)
               DO iZ = 1, nZ
                  write (aus(4+(iZ-1)*12:4+iZ*12-1),'(2x,g10.4)') Z(iZ)
                  ENDDO
               Print '(a80)', aus
               ENDDO ! K

         ELSE
            IF (nK.eq.1) THEN
               Print '(a)',
     *  ' only 1 spectrum and no z coordinates in file '//cl3(j)
            ELSE
               Print '(a)',
     *  ' WARNING/ no z coordinates though several spectra in file '//
     *  cl3(j)
               ENDIF
            ENDIF

         ENDDO ! lj
      Print *

      END ! MemInfoZ

      SUBROUTINE MemInfoK (nJList, Jlist, Fehler)
C     -------------------------------------------
            ! JWu 17may95

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER      Fehler*(*), Co*40, Un*40, h1*4, h2*20, aus*80
      INTEGER        nJList, JList(*), lj, j, nK, K, iOlfG, nZ, is, n
      REAL*8         X(MC), Z(MZ), dx
      LOGICAL        qed

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      DO lj = 1, nJList
      j = JList(lj)

         nK = iOlfG (j, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         Print '(a)', ' spectra of file '//cl3(j)

         ! List of spectra :
         DO K = 1, nK

            CALL OlfGetZ (j, K, nZ, Z, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            CALL OlfGetX (j, K, n, X, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            is = irSorted (X, n)
            IF     (is.eq.-2) THEN
               h1 = ' >> '
            ELSEIF (is.eq.-1) THEN
               h1 = ' >= '
            ELSEIF (is.eq. 0) THEN
               h1 = ' <> '
            ELSEIF (is.eq. 1) THEN
               h1 = ' <= '
            ELSEIF (is.eq. 2) THEN
               h1 = ' << '
               ENDIF
            IF (iabs(is).eq.2) THEN
               CALL CheckScale (n, X, 1.d-2, qed, dx)
            ELSE
               qed = .false.
               ENDIF
            IF (qed) THEN
               write (h2, '(a,g10.4)') ' dX = ', dx
            ELSE
               h2 = ' '
               ENDIF

            Print '(i3,a,g10.4,a,i4,a,g10.4,a4,g10.4,a17)',
     *             K, '#  z = ', Z(1), ' n = ', n, ' X = ',
     *             X(1), h1, X(n), h2

            ENDDO ! K
         ENDDO ! lj
      Print *

      END ! MemInfoK

      SUBROUTINE MemInfoY (nJList, JList, Fehler)
C     -------------------------------------------
            ! JWu 17may95
         ! directory of on-line data points

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      INTEGER        JList(*)

      CHARACTER      h1*4, h2*12, lisCh*80, lisChold*80, aus*80
      REAL*8         X(MC), Y(MC), D(MC), Z(MZ)
      LOGICAL        qList(MC)

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      DO lj = 1, nJList
      j = JList(lj)

         Print '(a)', ' listing x-y for file '//cl3(j)
         nK = (MemBlockNum(j) - MBH) / 4
         Kdef = Kold

 150     CONTINUE
         IF (nK.gt.1) THEN
            K = iAskDMu (' Spectrum no. (0=quit)', Kdef, -1, nK)
            IF (K.eq.0) RETURN
            Kold = K
         ELSE
            Print *, ' there is just one spectrum'
            K = 1
            ENDIF

         CALL OlfGetSpe (j, K, nZ, Z, nC, X, Y, D, Fehler)
         IF (nC.eq.0) THEN
            Print *, '  Spectrum is empty !'
            GOTO 150
            ENDIF

         lisCh = lisChold
 151     CONTINUE
         IF (nC.gt.12) THEN
            aus = '  Show which channels'
            CALL GetNList (aus, lisCh, qList, nC)
            IF (lisCh.eq.'-') GOTO 158
            lisChold = lisCh
         ELSE
            ! don't bore the user with that question,
            ! show all those few channels :
            DO i = 1, nC
               qList(i) = .true.
               ENDDO
            ENDIF

         Print '(i3,a,g10.4)', K,'#  z1 = ', Z(1)

         DO i = 1, nC
            IF (qList(i)) THEN
               Print '(a,i4,a,g12.6,a,g12.6,a,g10.4,a,g10.4)',
     *      ' ch', i, '  x = ', X(i), ', y = ', Y(i), ' +-', D(i),
     *      '; dx = ', rStepAt (X, nC, i)
               ENDIF
            ENDDO

         IF (nC.gt.12) THEN
            lisCh = '-'
            GOTO 151
            ENDIF
 158     CONTINUE
         Kdef = 0
         IF (nK.gt.1) GOTO 150
 159     CONTINUE
         ENDDO

      END ! MemInfoY

C  ====================================================================
C  i23 / 2 :   editing
C  ====================================================================

      INTEGER FUNCTION iLabOfCnu (Lab)
C     --------------------------------
            ! 14jan99

      IMPLICIT NONE
      CHARACTER     Lab*(*), cl6*6
      INTEGER       lenU, ichar1, il

      IF     (Lab.eq.'x') THEN
         iLabOfCnu = -1
         RETURN
      ELSEIF (Lab.eq.'y') THEN
         iLabOfCnu = 0
         RETURN
      ELSEIF (Lab(1:1).eq.'z') THEN
         IF     (lenU(Lab).gt.2) THEN
            CALL Gong (3)
            Print *, 'iLabOfCnu / unexpected label "',
     *              Lab(1:lenU(Lab)), '"'
            iLabOfCnu = -2
            RETURN
         ELSEIF (Lab.eq.'z') THEN
            CALL Gong (3)
            Print *,
     * 'iLabOfCnu / the use of "z" instead of "z1" is not recommended'
            iLabOfCnu = 1
            RETURN
            ENDIF
         il = ichar1(Lab(2:2))
         IF (il.lt.1 .or. il.gt.9) THEN
            Print *, 'iLabOfCnu / unexpected label "',Lab(1:lenU(Lab)),
     *           '", apparently no. '//cl6(il)
            iLabOfCnu = -3
            RETURN
            ENDIF
         iLabOfCnu = il
         RETURN
         ENDIF
         Print *, 'iLabOfCnu / unexpected label ', Lab
         iLabOfCnu = -4

      END ! iLabOfCnu

      SUBROUTINE EditCnu (nJList, JList, whext, Fehler)
C     -------------------------------------------------
            ! JWu 3jun91, completely renewed (FileTPar) 9mar93,
            ! new again (EditCnu) 18nov96/4mar98

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      INCLUDE 'i_dim.f'
      CHARACTER*(*) Fehler, whext
      CHARACTER*80  aus
      CHARACTER*40  UnJ, CoJ, UnE, CoE
      CHARACTER*4   which
      INTEGER       nJList, JList(*), lj, lji, iCo, nZ, iOlfG

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      ! whext z.Zt. ignoriert

C  Uebersichtstabelle
      nZ = iOlfG (JList(1), '#Z', Fehler)
 10   CONTINUE
      DO iCo = -1, nZ
         IF     (iCo.eq.-1) THEN
            which = 'x'
         ELSEIF (iCo.eq.0) THEN
            which = 'y'
         ELSE
            which = 'z'//cl2(iCo)
            ENDIF

         lj = 0
 11      CONTINUE
         lj = lj + 1
         CALL OlfCnuG (JList(lj), which, CoJ, UnJ, Fehler)
         IF (Fehler(1:4).eq.'&pnf') GOTO 19
         IF (Fehler.ne.'&ff') RETURN
         lji = lj
         DO lj = lji, nJList-1
            CALL OlfCnuG (JList(lj+1), which, CoE, UnE, Fehler)
            IF (Fehler(1:4).eq.'&pnf') GOTO 19
            IF (Fehler.ne.'&ff') RETURN
            IF (CoE.ne.CoJ .or. UnE.ne.UnJ) GOTO 15
            ENDDO
         lj = nJList
 15      CONTINUE

         Print '(a2,a,a20,a,a15,a,i3,a,i3)',
     * which, ' is ', CoJ, ' in ', UnJ, ' for files ', lji, ' .. ', lj

         IF (lj.lt.nJList) GOTO 11

 19      CONTINUE
         Fehler = '&ff'
         ENDDO ! iCo / Uebersichtstabelle

 20   CONTINUE
      Print *
      CALL FrageC (' Edit which coordinate (x,y,..)', which)
      IF (which.eq.' ') RETURN

      DO lj = 1, nJList
         CALL OlfCnuG (JList(lj), which, CoJ, UnJ, Fehler)
         IF (Fehler(1:4).eq.'&pnf') THEN
            Print *,
     * ' do not use edit_coordinate to add new coordinates'
            Fehler = '&ff'
            GOTO 20
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

         IF (lj.eq.1) THEN
            CoE = CoJ
            UnE = UnJ
         ELSE
            IF (CoJ.ne.CoE) CoE = '&?'
            IF (UnJ.ne.UnE) UnE = '&?'
            ENDIF
         ENDDO

      IF (nJList.eq.1) THEN
         CALL FrageTD (' Coordinate name of '//which, CoJ, CoJ)
         CALL FrageCD (' And its unit', UnJ, UnJ)
         CALL OlfCnuP (JList(1), which, CoJ, UnJ, Fehler)
      ELSE
         IF (CoE.eq.'&?') THEN
            CALL Compose2 (aus, ' Coordinate name of '//which,
     *                     ' [enter indivdually]')
            CALL FrageC (aus, CoE)
         ELSE
            CALL Compose2 (aus, ' Coordinate name of '//which,
     *                     ' (answer ''\ '' to enter indivdually) ')
            CALL FrageCD (aus, CoE, CoE)
            ENDIF
         IF (CoE.eq.' ') THEN ! enter individually
            DO lj = 1, nJList
               CALL OlfCnuG (JList(lj), which, CoJ, UnJ, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL FrageTD (' Coordinate name', CoJ, CoJ)
               CALL FrageCD (' And its unit', UnJ, UnJ)
               CALL OlfCnuP (JList(lj), which, CoJ, UnJ, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO
         ELSE
            IF(UnE.eq.'&?') THEN
               CALL FrageT (' And its unit', UnE)
            ELSE
               CALL FrageCD (' And its unit', UnE, UnE)
               ENDIF
            DO lj = 1, nJList
               CALL OlfCnuP (JList(lj), which, CoE, UnE, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO
            ENDIF
         ENDIF

      GOTO 10
c      IF (whext.eq.'?') GOTO 1

      END ! EditCnu

      SUBROUTINE EditDoc (nJList, JList, Fehler)
C     ------------------------------------------
            ! JWu 3jun91. Revised (Type as parameter, List) 2sep91.
            ! Completely renewed (FileTPar) 9mar93.
            ! Completely rewritten 9oct97

      IMPLICIT NONE
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      CHARACTER*(*) Fehler
      CHARACTER*80  line, tit
      CHARACTER*40  fil, doc, ein1
      INTEGER       nJList, JList(*), lj, j,
     *              iP, iPList(MTP), niPList, iPL, iiPL

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Loop over files/ no attempt is made to collate doc's form different files :
      DO lj = 1, nJList
         j = JList(lj)

         Print *, ' file '//cl3(j)

C  Show doc :
 10      CONTINUE
         Print *
         CALL tOlfG(j, 'fil', fil, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL tOlfG(j, 'tit', tit, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL tOlfG(j, 'doc', doc, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Print '(a3,2x,a40)', 'f', fil
         Print '(a3,2x,a72)', 't', tit
         Print '(a3,2x,a40)', 's', doc
         iP = 0
 11      CONTINUE
            CALL OlfComLinG (j, iP+1, line, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (line.eq.'&eoc') GOTO 19
            iP = iP + 1
            Print '(i3,2x,a72)', iP, line
            GOTO 11
 19      CONTINUE
         Print *

C  Modify :
 30      CONTINUE
         CALL FrageCD (' Modify lines (f,t,s,<list>,a,-)', ein1, '-')
         IF (ein1.eq.'a') THEN
            CALL FrageC ('  Add line:', line)
            IF (line.eq.' ') THEN
               CALL Gong (1)
               GOTO 30
               ENDIF
            IF (iP.ge.MTP) THEN
               CALL Gong (11)
               Print *, ' too many lines'
               GOTO 30
               ENDIF
            iP = iP + 1
            CALL OlfComLinP (j, iP, line, Fehler)
         ELSEIF (ein1.eq.'?') THEN
            Print *, '    f       modify file name'
            Print *, '    t       modify title line'
            Print *, '    s       modify short documentation'
            Print *, '    <list>  modify existing lines'
            Print *, '    a       add new line'
            Print *, '    -       quit'
         ELSEIF (ein1.eq.'f') THEN
            CALL FrageCD ('  File name', fil, fil)
            CALL tOlfP (j, 'fil', fil, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (ein1.eq.'t') THEN
            CALL FrageCD ('  Title line', tit, tit)
            CALL tOlfP (j, 'tit', tit, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (ein1.eq.'s') THEN
            CALL FrageCD ('  Short documentation', doc, doc)
            CALL tOlfP (j, 'doc', doc, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (ein1.eq.'-') THEN
            GOTO 39
         ELSEIF (ein1.eq.' ') THEN ! re-display
         ELSE ! list given or bad input
            CALL DecJList (ein1, MTP, niPList, iPList, 1, iP, Fehler)
            IF (Fehler.ne.'&ff') THEN
               CALL FehlerGong (Fehler, 1)
               GOTO 30
               ENDIF
            DO iPL = 1, niPList
               CALL OlfComLinG (j, iPList(iPL), line, Fehler)
               CALL FrageCD ('  Line '//cl3(iPList(iPL)), line, line)
               IF (line.eq.' ') THEN
                  CALL OlfComLinDel (j, iPList(iPL), Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  iP = iP - 1
                  DO iiPL = iPL+1, niPList
                     IF (iPList(iiPL).gt.iPList(iPL))
     *                   iPList(iiPL) = iPList(iiPL) -1
                     ENDDO
               ELSE
                  CALL OlfComLinP (j, iPList(iPL), line, Fehler)
                  ENDIF
               ENDDO
            ENDIF
         GOTO 10 ! redisplay

 39      CONTINUE

         Print *
         ENDDO ! lj

      END ! EditDoc

      SUBROUTINE EditZ (nJList, JList, Fehler)
C     ----------------------------------------
            ! 4mar98

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER      Fehler*(*), Co*40, Un*40, Co2*40, Un2*40,
     *               action*8, aus*80
      INTEGER        nJList, JList(*), lj, j, nK, K, iOlfG,
     *               nZ, iZ, iZ2, Kdum
      REAL*8         Z(MZ), zav, zval

 10   CONTINUE
      CALL MemInfoZ (nJList, JList, Fehler)

 20   CONTINUE
      nZ = iOlfG (JList(1), '#spectra', Fehler)
      CALL FrageC (' Modify (a,d,r,s,x)', action)
      IF     (action.eq.'h' .or. action.eq.'?') THEN
         Print *, ' the following modifications are accessible:'
         Print *, '    a  add one more z-coordinate'
         Print *, '    d  delete one z-coordinate'
         Print *, '    r  average one z-coordinate and save as r-param'
         Print *, '    s  sort: leftmost columns varying most rapidly'
         Print *, '    x  exchange order of z-coordinates'
         GOTO 20
      ELSEIF (action.eq.' ') THEN
         GOTO 99
      ELSEIF (action.eq.'a') THEN
         CALL FrageC ('Add a coordinate named', Co)
         IF (Co.eq.' ') GOTO 20
         CALL FrageC ('And its unit', Un)
         DO lj = 1, nJList
            j = JList(lj)
            CALL OlfOpen (j, 1, Kdum, Fehler)
            nK = iOlfG (j, '#spectra', Fehler)
            nZ = iOlfG (j, '#Z', Fehler)
            CALL OlfCnuP (j, 'z'//cl2(nZ+1), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            Print *, ' set new z for file '//cl4(j)
            zval = 0
            DO K = 1, nK
               CALL OlfGetZ (j, K, nZ, Z, Fehler)
               IF (Fehler.ne.'&ff') GOTO 99
               zval = rAskD ('Value for spectrum '//cl4(K), zval)
               Z(nZ+1) = zval
               CALL OlfPutZ (j, K, nZ+1, Z, Fehler)
               IF (Fehler.ne.'&ff') GOTO 99
               ENDDO
            CALL OlfClos (j, 0, Fehler)
            ENDDO
      ELSEIF (action.eq.'d') THEN
         iZ = iAskMu (' delete which z-coordinate', 0, MZ)
         IF (iZ.le.0) GOTO 20
         DO lj = 1, nJList
            j = JList(lj)
            CALL OlfOpen (j, 1, Kdum, Fehler)
            CALL OlfDel1Z (j, j, iZ, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL OlfClos (j, 0, Fehler)
            ENDDO
      ELSEIF (action.eq.'r') THEN
         iZ = iAskMu (' average which z-coordinate', 0, MZ)
         IF (iZ.le.0) GOTO 20
         DO lj = 1, nJList
            j = JList(lj)
            CALL OlfOpen (j, 1, Kdum, Fehler)
            nK = iOlfG (j, '#spectra', Fehler)
            zav = 0
            DO K = 1, nK
               CALL OlfGet1Z (j, K, iZ, zval, Fehler)
               zav = zav + zval
               ENDDO
            zav = zav / nK
            CALL OlfCnuG (j, 'z'//cl2(iZ), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL OlfDel1Z (j, j, iZ, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL rOlfP (j, Co, Un, zav, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL OlfClos (j, 0, Fehler)
            ENDDO
      ELSEIF (action.eq.'s') THEN
         CALL OrgSpectraSort (nJList, JList, .true., Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         GOTO 10
      ELSEIF (action.eq.'x') THEN
         IF (nZ.gt.2) THEN
            iZ = iAskMu (' move which z-coordinate', 0, MZ)
            IF (iZ.le.0) GOTO 20
         ELSEIF (nZ.lt.2) THEN
            CALL Gong (4)
            Print *, ' nothing to move'
            GOTO 20
         ELSE
            iZ = 2
            ENDIF
         iZ2 = iAskDMu (' new position', 1, 0, MZ)
         IF (iZ2.eq.0 .or. iZ2.eq.iZ) THEN
            CALL Gong (1)
            GOTO 20
            ENDIF
         DO lj = 1, nJList
            j = JList(lj)
            CALL OlfOpen (j, 1, Kdum, Fehler)
            nK = iOlfG (j, '#spectra', Fehler)
            DO K = 1, nK
               CALL OlfGetZ (j, K, nZ, Z, Fehler)
               IF (Fehler.ne.'&ff') GOTO 99
               zav = Z(iZ)
               Z(iZ) = Z(iZ2)
               Z(iZ2) = zav
               CALL OlfPutZ (j, K, nZ, Z, Fehler)
               IF (Fehler.ne.'&ff') GOTO 99
               ENDDO
            CALL OlfCnuG (j, 'z'//cl2(iZ), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL OlfCnuG (j, 'z'//cl2(iZ2), Co2, Un2, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL OlfCnuP (j, 'z'//cl2(iZ), Co2, Un2, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL OlfCnuP (j, 'z'//cl2(iZ2), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            CALL OlfClos (j, 0, Fehler)
            ENDDO
      ELSE
         CALL Gong (4)
         Print *, ' unknown option - type h for help'
         GOTO 20
         ENDIF

      GOTO 10

 99   CONTINUE
      IF (Fehler.ne.'&ff') CALL FehlerGong (Fehler, 3)

      END ! EditZ

      SUBROUTINE EditRPar (nJList, JList, Fehler)
C     -------------------------------------------
            ! JWu 10jul91, renewed 18mar92, completely new 24nov96, 8oct97
         ! List and change rPar

      IMPLICIT NONE
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      CHARACTER*(*) Fehler
      CHARACTER*40  file, Co, Un, CoPout(MRP), UnPout(MRP),
     *              h3, ein, ein1
      LOGICAL       qPdef(MRP), qPuni(MRP)
      CHARACTER*16  h1, h2
      INTEGER       nJList, JList(*), lj, j, np, iP, iPout, nPout,
     *              niPList, iPList(MRP), Kdum, nK, K, nZ, iOlfG
      REAL*8        rval, rvalmin, rvalmax, rvalavg, rOlfG,
     *              rOlfGG, zval, Z(MZ)

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Make list of r:
 10   CONTINUE
      nPout = 0
      DO lj = 1, nJList
         DO iP = 1, MRP
            Co   = '&pbn '//cv3(iP)
            rval = rOlfG (JList(lj), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (Co.eq.'&empty') GOTO 18
            IF (Co.eq.'&eop') GOTO 19
            DO iPout = 1, nPout
               IF (CoPout(iPout).eq.Co) THEN
                  IF (UnPout(iPout).ne.Un) THEN
                     UnPout(iPout) = '&diff' ! 'different units in '//Co
                     ENDIF
                  GOTO 18 ! parameter already in list
                  ENDIF
               ENDDO
            nPout = nPout + 1 ! new parameter
            IF (nPout.gt.MRP) THEN
               Fehler = 'too many r-parameter for being listed here'
               RETURN
               ENDIF
            CoPout(nPout) = Co
            UnPout(nPout) = Un
 18         CONTINUE
            ENDDO
 19      CONTINUE
         ENDDO ! lj

C  Compare values of r:
      DO iPout = 1, nPout
         Co = CoPout(iPout)
         Un = UnPout(iPout)
         qPdef(iPout) = .true.
         qPuni(iPout) = .false.
         rval = rOlfGG (JList(1), Co, Un, Fehler)
         IF (Fehler.ne.'&ff') THEN
            Fehler = '&ff'
            qPdef(iPout) = .false.
            GOTO 228
            ENDIF
         rvalmin = rval
         rvalmax = rval
         rvalavg = rval
         DO lj = 2, nJList
            rval = rOlfGG (JList(lj), Co, Un, Fehler)
            IF (Fehler.ne.'&ff') THEN
               Fehler = '&ff'
               qPdef(iPout) = .false.
               GOTO 228
               ENDIF
            rvalmin = dmin1 (rvalmin, rval)
            rvalmax = dmax1 (rvalmax, rval)
            rvalavg = rvalavg + rval
            ENDDO
         rvalavg = rvalavg / nJList
         IF (rvalmin.ne.rvalmax) THEN
            write (h1, '(g16.8)') rvalmin
            write (h2, '(g16.8)') rvalmax
            h3 = 'varies from '//h1//' to '//h2
         ELSE
            write (h3, '(g16.8)') rvalavg
            qPuni(iPout) = .true.
            ENDIF
         GOTO 229
 228     CONTINUE ! arrived here if qPdef = false
         h3 = 'for some files undefined'
 229     CONTINUE
         Print '(i3,3x,a16,2x,a12,2x,a)', iPout, Co, Un, h3
         ENDDO ! r-par
      Print *

C  Menu :
 30      CONTINUE
         CALL FrageCD (' Modify parameters (<list>,a,-)', ein1, '-')
         Print *
 301     CONTINUE
         IF (ein1.eq.'a') THEN
            Print *, '  Add parameter:'
            nPout = nPout + 1
            iPout = nPout
            CALL FrageC ('  Coordinate', CoPout(iPout))
            CALL FrageC ('  And its unit', UnPout(iPout))
            DO lj = 1, nJList
               CALL rOlfP (JList(lj), CoPout(iPout), UnPout(iPout),
     *                     0.d0, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO
            ein1 = cl4(iPout)
            GOTO 301 ! goto `ELSE'
         ELSEIF (ein1.eq.'?') THEN
            Print *, '    <list>  modify existing parameters'
            Print *, '    a       add new parameter'
            Print *, '    -       quit'
            ein1 = ' '
         ELSEIF (ein1.eq.'-') THEN
            RETURN
         ELSEIF (ein1.eq.' ') THEN ! re-display
         ELSE ! list given or bad input
            CALL DecJList (ein1, MRP, niPList, iPList, 1,nPout,Fehler)
            ein1 = ' '
            IF (Fehler.ne.'&ff') THEN
               CALL FehlerGong (Fehler, 1)
               Print *
               GOTO 10
               ENDIF
            DO iP = 1, niPList
               iPout = iPList(iP)
 31            CONTINUE
               Print *, '  parameter: '//CoPout(iPout)
               ! display current values :
               IF (.not.qPuni(iPout)) THEN
                  DO lj = 1, nJList
                     rval = rOlfG (JList(lj), CoPout(iPout),Un,Fehler)
                     IF (Fehler.ne.'&ff') THEN
                        Fehler = '&ff'
                        h3 = ' ** undefined'
                     ELSE
                        write (h3, '(g16.8,1x,a20)') rval, Un
                        ENDIF
                     Print '(a,i2,a,a)', '    present value for file ',
     *                  JList(lj), ': ', h3
                     ENDDO
                  Print *
                  ENDIF
               CALL FrageC ('  Modification (val,f,n,u,z,d) [none]',
     *                      ein)
               CALL Fi1R (ein, rval)
               IF (ein.eq.'#') THEN
                  DO lj = 1, nJList
                     CALL rOlfP (JList(lj), CoPout(iPout),
     *                 UnPout(iPout), rval, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO
               ELSEIF (ein.eq.'f') THEN
                  DO lj = 1, nJList
                     rval = rOlfG (JList(lj), CoPout(iPout),Un,Fehler)
                     rval = rAskD ('   Value for file '//
     *                             cl4(JList(lj)), rval)
                     CALL rOlfP (JList(lj), CoPout(iPout), Un, rval,
     *                           Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO
               ELSEIF (ein.eq.'n') THEN
                  CALL FrageCD ('   Coordinate', Co, CoPout(iPout))
                  CALL FrageCD ('   And its unit', Un, UnPout(iPout))
                  DO lj = 1, nJList
                     rval = rOlfG (JList(lj),
     *                           CoPout(iPout), UnPout(iPout), Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     IF (CoPout(iPout).ne.Co)
     *                 CALL rOlfDel (JList(lj), CoPout(iPout), Fehler)
                     CALL rOlfP (JList(lj), Co, Un, rval, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     CoPout(iPout) = Co
                     UnPout(iPout) = Un
                     ENDDO
               ELSEIF (ein.eq.'u') THEN
                  Print *, '   unfertig '
               ELSEIF (ein.eq.'z') THEN
                  Print *, ' bricolage pour mfj'
                  DO lj = 1, nJList
                     Co = CoPout(iPout)
                     Un = UnPout(iPout)
                     j = JList(lj)
                     CALL OlfOpen (j, 1, Kdum, Fehler)
                     nK = iOlfG (j, '#spectra', Fehler)
                     nZ = iOlfG (j, '#Z', Fehler)
                     CALL OlfCnuP (j, 'z'//cl2(nZ+1), Co, Un, Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     zval = 0
                     DO K = 1, nK
                        CALL OlfGetZ (j, K, nZ, Z, Fehler)
                        IF (Fehler.ne.'&ff') GOTO 999
                        zval = rOlfGG(j, Co, Un, Fehler)
                        Z(nZ+1) = zval
                        CALL OlfPutZ (j, K, nZ+1, Z, Fehler)
                        IF (Fehler.ne.'&ff') GOTO 999
                     ENDDO
                     CALL OlfClos (j, 0, Fehler)
                     CALL rOlfDel (JList(lj), CoPout(iPout), Fehler)
                  ENDDO
               ELSEIF (ein.eq.'d') THEN
                  DO lj = 1, nJList
                     CALL rOlfDel (JList(lj), CoPout(iPout), Fehler)
                     IF (Fehler.ne.'&ff') Fehler = '&ff'
                     ENDDO
               ELSEIF (ein.eq.'?') THEN
                  Print *, '    <real value>  global parameter value'
                  Print *, '    f             different value for'//
     *                     ' each file'
                  Print *, '    n             coordinate name and unit'
                  Print *, '    u             unit conversion'
                  Print *, '    d             delete'
                  Print *, '    z             promote r to z'
                  Print *, '    <RETURN>      no modification'
                  GOTO 31
               ELSEIF (ein.eq.' ') THEN
                  ! regular exit
               ELSE
                  CALL Gong (3)
                  GOTO 31
                  ENDIF
               ENDDO ! iP
            ENDIF
            Print *
            GOTO 10

 999        CONTINUE
            IF (Fehler.ne.'&ff') CALL FehlerGong (Fehler, 3)


      END ! EditRPar

      SUBROUTINE EditIPar (nJList, JList, Fehler)
C     -------------------------------------------
         ! List and change iPar

      IMPLICIT NONE
      INCLUDE      'l_def.f'
      CHARACTER*(*) Fehler
      CHARACTER     co*24, un*24, file*40
      INTEGER       nJList, JList(*), lj, j, np
      INTEGER       ival, iOlfG

      DO lj = 1, nJList
         j = JList(lj)

         CALL tOlfG(j, 'fil', file, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL Say2 (' integer parameters of file '//cl3(j),' : '//file)

         np = 0
 11      CONTINUE
            co   = '&pbn '//cv3(np+1)
            ival = iOlfG (j, co, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (co.eq.'&eop') GOTO 19
            np = np + 1
            IF (co.eq.'&empty') GOTO 11
            Print '(i3,2x,a16,2x,i8)', np, co, ival
            GOTO 11

 19      CONTINUE

         ENDDO ! lj

      END ! EditIPar

      SUBROUTINE EditGPar (nJList, JList, Fehler)
C     -------------------------------------------
            ! JWu 4-6nov91
         ! Manipulate the iPar and rPar used for plotting

      IMPLICIT NONE
      INCLUDE      'l_def.f'
      CHARACTER*(*) Fehler
      INTEGER       nJList, JList(*), lj, j, icu, ips, ncp, iOlfGdef
      REAL*8        rpi, rpf, rOlfGGdef

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF
      DO lj = 1, nJList
         j = JList(lj)

         icu = iOlfGdef (j, '?cu',      0, Fehler)
         ips = iOlfGdef (j, 'plot-sy#', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ips = iAskD (
     * 'Linestye(>0)  or plotsymbol(<0) or automatic(0)', ips)
         CALL iOlfP (j, 'plot-sy#', ips, Fehler)

         IF (qintr(icu)) THEN
            ncp = iOlfGdef  (j, 'plot-#pts',   0,    Fehler)
            rpi = rOlfGGdef (j, 'plot-i', ' ', 0.d0, Fehler)
            rpf = rOlfGGdef (j, 'plot-f', ' ', 0.d0, Fehler)

            ncp = iAskDMu ('Number of points', ncp, 2, 10000)
            CALL rAskRgeFull ('Plot range', rpi, rpf, rpi, rpf)

            CALL iOlfP (j, 'plot-#pts',   ncp, Fehler)
            CALL rOlfP (j, 'plot-i', ' ', rpi, Fehler)
            CALL rOlfP (j, 'plot-f', ' ', rpf, Fehler)
            ENDIF

         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      END ! EditGPar

C  ====================================================================
C  i23 / 3 :   handling of z/r
C  ====================================================================

      SUBROUTINE AskNoZ (quest, j, iZ, cZ, qMustHaveZ, Fehler)
C     --------------------------------------------------------
         ! ask for number of z-parameter

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      CHARACTER      quest*(*), cZ*(*), Fehler*(*)
      INTEGER        j, iZ, nZ, iOlfG
      LOGICAL        qMustHaveZ

      IF (Fehler.ne.'&ff') THEN
         Print *, ' Fehler on entry in AskNoZ'
         RETURN
         ENDIF

      nZ = iOlfG (j, '#Z', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF (nZ.gt.1) THEN
         iZ = iAskDMu (quest, iZ, 1, nZ)
      ELSEIF (nZ.gt.0) THEN
         iZ = 1
      ELSE
         IF (qMustHaveZ) THEN
            Fehler = ' no z defined'
         ELSE
            iZ = 0
            ENDIF
         cZ = 'z-'
         RETURN
         ENDIF

      cZ = 'z'//cl2(iZ)

      END ! AskNoZ

      SUBROUTINE TensorCheckZ (j, Fehler)
C     -----------------------------------
            ! 4mar98

      IMPLICIT NONE
      CHARACTER     Fehler*(*), Co*40, Un*40, cl2*2
      INTEGER       j, nZ, iOlfG, iZ, K, nK
      REAL*8        z, z2

      nZ = iOlfG (j, '#Z', Fehler)
      nK = iOlfG (j, '#spectra', Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

      DO iZ = nZ, 1, -1

         CALL OlfGet1Z (j, 1, iZ, z, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         DO K = 2, nK
            CALL OlfGet1Z (j, K, iZ, z2, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            IF (z2.ne.z) GOTO 19 ! z remains useful
            ENDDO

         ! z is constant => becomes rPar
         CALL OlfCnuG (j, 'z'//cl2(iZ), Co, Un, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         CALL OlfDel1Z (j, j, iZ, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         CALL rOlfP (j, Co, Un, z, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99

 19      CONTINUE
         ENDDO

      RETURN
 99   CONTINUE
      Print *, 'error passed through TensorCheckZ'

      END ! TensorCheckZ

      SUBROUTINE TensorSaveInt (j, jout, nK, nKout, Yin, Din, Fehler)
C     ---------------------------------------------------------------
            ! JWu 6mrz98
         ! save result of integral operation (used by oi and ci)

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER      Fehler*(*), aus*80, Co*40, Un*40,
     *               CoG*40, UnG*40       ! letzte 2 DEBUG
      INTEGER        j, jout, nK, K, KK, nKout, KNew(MK), nZ, iZ,
     *               iZsel, n, iOlfG
      REAL*8         Z(MZ), ZZ(MZ), Yin(*), Din(*), Xout(MC), Yout(MC),
     *               Dout(MC)
      LOGICAL        qMultiDim

C  special case : one single data point
      nZ = iOlfG (j, '#Z', Fehler)
      IF (nZ.eq.0) THEN
         IF (nK.gt.1) THEN
            Fehler = 'TensorSaveInt/ nZ=0 but nK>1'
            RETURN
            ENDIF
         Xout(1) = 0
         CALL OlfPutSpe (jout, 1, 0, Z, 1, Xout, Yin, Din, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'x', '&no x', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         RETURN
         ENDIF

C  determine groups :
      DO K = 1, nK
         KNew(K) = K
         ENDDO
      DO K = 1, nK
         IF (KNew(K).eq.K) THEN ! begin of new group
            CALL OlfGetZ (j, K, nZ, Z, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            DO KK = K+1, nK
               IF (KNew(KK).eq.KK) THEN ! not yet regrouped
                  CALL OlfGetZ (j, KK, nZ, ZZ, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  DO iZ = 2, nZ
                     IF (ZZ(iZ).ne.Z(iZ)) GOTO 180
                     ENDDO
                  ! Z2, Z3, ... are all equal
                  KNew(KK) = K
 180              CONTINUE
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

C  save :
      IF (j.eq.jout) THEN
         Fehler = 'TensorSaveInt/ cannot overwrite'
         RETURN
         ENDIF

C  user choice :
      IF (KNew(nK).gt.1) THEN ! there are indeed several groups
         DO iZ = 1, nZ
            CALL OlfCnuG (j, 'z'//cl2(iZ), Co, Un, Fehler)
            IF (Fehler(1:4).eq.'&pnf') THEN
               Print *, 'PROGRAM ERROR/ CnuG in TensorCheckInt'
               RETURN
               ENDIF
            Print '(a,a,a20,a,a15)', 'z'//cl2(iZ), ' is ',
     *            Co, ' in ', Un
            ENDDO
         CALL Compose2 (aus,
     *     'Save in several spectra (0) or retain only one z (1-'//
     *     cl2(nZ), ')')
         iZsel = iAskDMu (aus, iZsel, 0, nZ)
         qMultiDim = (iZsel.eq.0)
      ELSE
         qMultiDim = .true. ! eigentlich Quatsch (z identisch -> ist r, nicht z)
         ENDIF

      IF (qMultiDim) THEN
         nKout = 0
         DO K = 1, nK
            IF (KNew(K).eq.K) THEN
               nKout = nKout + 1
               CALL OlfCopZ (j, jout, K, nKout, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               n = 0
               DO KK = K, nK
                  IF (KNew(KK).eq.K) THEN
                     n = n + 1
                     CALL OlfGet1Z (j, KK, 1, Xout(n), Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     Yout(n) = Yin(KK)
                     Dout(n) = Din(KK)
                     ENDIF
                  ENDDO
               CALL OlfPutXYD (jout, nKout, n, Xout, Yout,Dout,Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDIF
            ENDDO
         ! now z1 has become useless :
         CALL OlfCnuG (j, 'z1', Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'x', Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfDel1Z (jout, jout, 1, Fehler)
         IF (Fehler.ne.'&ff') RETURN

      ELSE ! only one z survives: everything becomes very simple (11jan99)
         DO K = 1, nK
            CALL OlfGet1Z (j, K, iZsel, Xout(K), Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO
         CALL OlfPutSpe (jout, 1, 0, Z, nK, Xout, Yin, Din, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (j, 'z'//cl2(iZsel), Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout, 'x', Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDIF

      END ! TensorSaveInt

      SUBROUTINE SetZXofZ (j1, j2, IZXofZ, iz2ofx, nZ2, Fehler)
C     ---------------------------------------------------------
            ! JWu 9mrz98
         ! for all iZ1 of file 1:
         !    search for corresponding iZ2 of file 2
         ! return results as IZXofZ(iZ1) := iZ2
         ! `corresponding' means that Co(iZ) and Un(iZ) agree
         ! special case:
         !    iZ2=0 stands for CoX instead of CoZ
         !    i2ofx denotes the iZ1 for which this special case happens

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      CHARACTER     Fehler*(*), Co1*40, Un1*40, Co2*40, Un2*40
      INTEGER       j1, j2, IZXofZ(*), iz2ofx, iOlfG, nZ1, nZ2,iZ1,iZ2

      nZ1 = iOlfG (j1, '#Z', Fehler)
      nZ2 = iOlfG (j2, '#Z', Fehler)
      iz2ofx = 0
      IF (Fehler.ne.'&ff') RETURN
      DO iZ1 = 1, nZ1
         CALL OlfCnuG (j1, 'z'//cl2(iZ1), Co1, Un1, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         DO iZ2 = 1, nZ2
            CALL OlfCnuG (j2, 'z'//cl2(iZ2), Co2, Un2, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (Co1.eq.Co2) GOTO 15
            ENDDO
         iZ2 = 0
         CALL OlfCnuG (j2, 'x', Co2, Un2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (Co1.eq.Co2) GOTO 15

         ! coordinate not found
         CALL Say2 ('  GARDEZ/ coordinate '//
     *              Co1, ' not found in 2nd file')
         iZ2 = -1
         GOTO 18
         RETURN

 15      IF (Un1.ne.Un2) THEN
            Fehler = ' different units for '//Co1
            RETURN
            ENDIF

 18      IZXofZ(iZ1) = iZ2
         IF (iZ2.eq.0) iz2ofx = iZ1

         ENDDO

      END ! SetZXofZ

      SUBROUTINE GetKofKbyZX (nZ2, nK2, Z1, nZ1, IZXofZ, K2, Fehler)
C     --------------------------------------------------------------
         !

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f' ! ZZofK
      CHARACTER    Fehler*(*)
      REAL*8       Z1(MZ)
      INTEGER      nZ2, nK2, nZ1, K, K2, IZXofZ(MZ), iZ1, KK,
     *             iZ2, nKfound
      LOGICAL      qK(MK)

      DO KK = 1, nK2
         qK(KK) = .true.
         ENDDO
      DO iZ1 = 1, nZ1
         iZ2 = IZXofZ(iZ1)
         IF (iZ2.ge.1) THEN
            DO KK = 1, nK2
               IF (ZZofK(KK,iZ2).ne.Z1(iZ1)) qK(KK) = .false.
               ENDDO
            ENDIF
         ENDDO

      nKfound = iqSum (qK, nK2)
      IF (nKfound.eq.0) THEN
         Fehler = 'no correspondence found for spectrum '
         RETURN
      ELSEIF (nKfound.gt.1) THEN
         Fehler = 'integral point not unique for spectrum '
         RETURN
      ELSE
         DO K2 = 1, nK2
            IF (qK(K2)) RETURN
            ENDDO
         ENDIF

      END ! GetKofKbyZX
