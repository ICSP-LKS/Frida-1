C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i20   :     on-line memory
C
C  ====================================================================

C     General information :
C        Module written by J.Wuttke, jan/feb91
C        Here I have reinvented the use of pointers and the dynamical
C        allocation of memory. A little bit complicated, but it works.

C        Completely new data format; almost everything rewritten mar/may95.

C     Contents :
C        1.  On-line memory / block level :
C               MemBlockPut/Get/Inq/Num/Siz,
C               MemBlSubAdd/Del/Get/Put, MemRestInfo,
C               MemDims, qMemFill, MemFileDel
C                  ! these routines just handle a number of files consisting
C                  ! of different blocks; they ignore everything about the
C                  ! contents of the blocks.
C        2.  On-line memory / status shell :
C               MemFileStatP/G, MemFileBitP, MemBlockOvr
C                  ! block(1) contains the status bits; the routines of this
C                  ! level provide protected access to on-line files
C        3.  On-line memory / data translation :
C        4.  On-line memory / user shell :
C               OlfOpen/Clos, OlfParP/G, OlfSpeP/G, OlfGet.., OlfPtrG, FileClean

C     Aenderungsverzeichnis :
C        JWu  6feb97 : Restructured in view of poly-z
C        JWu 24nov96 : *Par replaced by direct labelling
C        JWu  9nov96 : OlfParG/P in outer routines replaced by i/r/tOlfG/P
C        JWu 17may95 : New structure has old functionality
C        JWu 15mar95 : Block structure (preparations for major revision)
C        JWu 25jan95 : Status bits (since feb91 only minor changes)
C        JWu 26nov91 : Parameter qOpen [now iPar(20)(1)
C        JWu 16sep91 : MemGet.. to replace SpectrumTake
C        JWu 31jul91 : Interactive manipulations -> IDA3
C        JWu 19jun91 : iOnext(0:MF)
C        JWu   apr91 : Modulaufteilung
C        JWu   feb91 : FileModif, SpectrumModif
C        JWu   jan91 : on-line memory

C  ====================================================================
C     1.  On-line memory / low-level : access blocks
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      BLOCK DATA MemBlockPreset
C     -------------------------

      IMPLICIT NONE
      INCLUDE       'i_dim.f'
      INTEGER        nF, nFB, nFA
      REAL*8         FMem

      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

      DATA           nF / 0 /, nFB(0) / 0 /, nFA(0) / 0 /

      END ! MemBlockPreset

      SUBROUTINE MemBlockPut (j, k, n, R, Fehler)
C     -------------------------------------------
            ! JWu 15-17mar95 (reusing OlfOpen 23.1.91, 12.2.91)
         ! Save R(1..n) as block k of file j in on-line memory

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

C  Parameters :
      INTEGER       j, k, n
      REAL*8        R(*)
      CHARACTER     Fehler*(*)

C  Internals :
      INTEGER       jj, kk, kmx, kdel, i, ndel, nold

C  Data / on-line memory / :
      INTEGER       nF, nFB, nFA
      REAL*8        FMem
      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

C  Checks :
      IF     (Fehler.ne.'&ff') THEN
         CALL Gong (3)
         Print *, 'error on entry in MemBlockPut'
         RETURN
         ENDIF

c      Print '(15x,15i4)', (nFA(i), i=1,15) ! <<< DEBUG
c      Print '(a,3i4,2f12.2)', 'MBP: ', j, k, n, R(1), R(n) ! <<< DEBUG

C  New File ?
      IF (j.lt.0 .or. j.gt.MF) THEN
         Fehler = 'PROG ERR/ MemBlockPut/ j oor'
         RETURN
      ELSEIF (j.gt.nF) THEN
         Fehler = 'file not load'
         RETURN
      ELSEIF (j.le.0 .and. nF.eq.MF) THEN
         Fehler = 'too many files in memory - new file cannot be saved'
         RETURN
      ELSEIF (j.eq.0) THEN
         nF = nF + 1
         j = nF
         ! New file is empty :
         nFB(j) = nFB(j-1)
         ENDIF

C  New Block ?
c      Print '(40x,a,2i5)', '-> ', nFB(j-1)+k, nFA(nFB(j-1)+k-1) ! <<< DEBUG
      kmx = nFB(j) - nFB(j-1)
      IF (k.le.0 .or. k.gt.kmx+1) THEN ! delete file (k=0 or error)
         IF (k.lt.0 .or. k.gt.MB) THEN
            Fehler = 'PROG ERR/ MemBlockPut/ k oor'
         ELSEIF (k.gt.kmx+1) THEN
            CALL Compose2 (Fehler, 'MemBlockPut/ P-ERR/ k='//
     * cl6(k),' exceeds no. of stored blocks kmx+1='//cl6(kmx+1))
            ENDIF
         ndel = nFA(nFB(j)) - nFA(nFB(j-1))
         DO i = nFA(nFB(j))+1, nFA(nFB(nF))
            FMem(i-ndel) = FMem(i)
            ENDDO
         kdel = nFB(j) - nFB(j-1)
         DO kk = nFB(j), nFB(nF)
            nFA(kk-kdel) = nFA(kk) - ndel
            ENDDO
         DO jj = j, nF
            nFB(jj-1) = nFB(jj) - kdel
            ENDDO
         nF = nF - 1
         RETURN
      ELSEIF (k.eq.kmx+1) THEN ! insert new block
         IF (nFA(nFB(nF))+n.gt.Mmem) THEN
            CALL Compose3 (Fehler,
     *         'MBP/ not enough space to write file '//cl3(j),
     *         ' block '//cl3(k), ' in on-line-memory')
            RETURN
            ENDIF
         DO jj = nF, j, -1
            nFB(jj) = nFB(jj) + 1
            ENDDO
         DO kk = nFB(nF), nFB(j), -1 ! nachgetragen 27jul95
            nFA(kk) = nFA(kk-1)
            ENDDO
         !         nFA(nFB(j-1)+k) = nFA(nFB(j-1)+k-1) ! block is empty
         ENDIF

C  Insert data :
      nold = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)
      IF (n.lt.0) THEN ! delete block
         ndel = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)
         DO i = nFA(nFB(j-1)+k)+1, nFA(nFB(nF)) ! error until 6mrz98
            FMem(i-ndel) = FMem(i)
            ENDDO
         DO kk = nFB(j-1)+k, nFB(nF)
            nFA(kk-1) = nFA(kk) - ndel
            ENDDO
         DO jj = j, nF
            nFB(jj) = nFB(jj) - 1
            ENDDO
         RETURN
      ELSEIF (n.eq.nold) THEN ! simply overwrite
         DO i = 1, n
            FMem(nFA(nFB(j-1)+k-1)+i) = R(i)
            ENDDO
         RETURN
      ELSEIF (n.lt.nold) THEN ! overwrite and compress
         DO i = 1, n
            FMem(nFA(nFB(j-1)+k-1)+i) = R(i)
            ENDDO
         ndel = nold - n
         DO i = nFA(nFB(j-1)+k)+1, nFA(nFB(nF)) ! error until 6mrz98
            FMem(i-ndel) = FMem(i)
            ENDDO
         DO kk = nFB(j-1)+k, nFB(nF)
            nFA(kk) = nFA(kk) - ndel
            ENDDO
         RETURN
      ELSEIF (n.gt.nold) THEN ! expand and fill in
         ndel = n - nold
         IF (nFA(nFB(nF))+ndel.gt.Mmem) THEN
            CALL Compose3 (Fehler,
     *         'MBP/ not enough space to write file '//cl3(j),
     *         ' block '//cl3(k), ' in on-line-memory')
            RETURN
            ENDIF
         DO i = nFA(nFB(nF)), nFA(nFB(j-1)+k)+1, -1
            FMem(i+ndel) = FMem(i)
            ENDDO
         DO i = 1, n
            FMem(nFA(nFB(j-1)+k-1)+i) = R(i)
            ENDDO
         DO kk = nFB(j-1)+k, nFB(nF)
            nFA(kk) = nFA(kk) + ndel
cdeb            Print *, '>>  kk ndel nFA(kk)_neu = ', kk, ndel, nFA(kk)
            ENDDO
         RETURN
      ELSE
         Fehler = 'PROG ERR/ MemBlockPut/ case n = ??'
         ENDIF

      END ! MemBlockPut

      SUBROUTINE MemBlSubAdd (j, k, noldcheck, nadd, Fehler)
C     ------------------------------------------------------
            ! JWu 24nov96
         ! extend block by nadd new entries
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      INTEGER       j, k, noldcheck, nadd
      CHARACTER     Fehler*(*)
      INTEGER       jj, kk, kmx, kdel, i, ndel, nold

C  Data / on-line memory / :
      INTEGER       nF, nFB, nFA
      REAL*8        FMem
      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

C  Checks :
      IF     (Fehler.ne.'&ff') THEN
         CALL Gong (3)
         Print *, 'error on entry in MemBlSubAdd'
         RETURN
         ENDIF

      IF (j.lt.1 .or. j.gt.nF) THEN
         Fehler = 'PROG ERR/ MemBlSubAdd/ j oor'
         RETURN
         ENDIF
      kmx = nFB(j) - nFB(j-1)
      IF (k.le.0 .or. k.gt.kmx) THEN ! delete file (k=0 or error)
         Fehler = 'PROG ERR/ MemBlSubAdd/ k oor'
         ENDIF
      ! hier k"onnte man evtl auch neuen Block erzeugen lassen

C  Insert data :
      nold = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)
      IF (nold.ne.noldcheck) THEN
         Print *, 'j k nold[intern] nold[extern] ', j, k, nold,
     *            noldcheck
         Fehler = 'PROG ERR/ MemBlSubAdd/ nold deviates'
         RETURN
         ENDIF
      IF (nadd.le.0) THEN
         Fehler = 'PROG ERR/ MemBlSubAdd/ nadd<=0'
         RETURN
         ENDIF
      ! expand and fill in
      IF (nFA(nFB(nF))+nadd.gt.Mmem) THEN
         CALL Compose3 (Fehler,
     *      'MBSA/ not enough space to write file '//cl3(j),
     *      ' block '//cl3(k), ' in on-line-memory')
         RETURN
         ENDIF
      DO i = nFA(nFB(nF)), nFA(nFB(j-1)+k)+1, -1
         FMem(i+nadd) = FMem(i)
         ENDDO
      DO i = 1, nadd
         FMem(nFA(nFB(j-1)+k-1)+nold+i) = 0
         ENDDO
      DO kk = nFB(j-1)+k, nFB(nF)
         nFA(kk) = nFA(kk) + nadd
         ENDDO

      END ! MemBlSubAdd

      SUBROUTINE MemBlSubDel (j, k, npos, ndel, Fehler)
C     -------------------------------------------------
            ! JWu 3/4feb97
         ! delete some lines from a block
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      INTEGER       j, k, npos, ndel
      CHARACTER     Fehler*(*)
      INTEGER       kk, kmx, kdel, i, nold

C  Data / on-line memory / :
      INTEGER       nF, nFB, nFA
      REAL*8        FMem
      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

C  Checks :
      IF     (Fehler.ne.'&ff') THEN
         CALL Gong (3)
         Print *, 'error on entry in MemBlSubDel'
         RETURN
         ENDIF

      IF (j.lt.1 .or. j.gt.nF) THEN
         Fehler = 'PROG ERR/ MemBlSubDel/ j oor'
         RETURN
         ENDIF
      kmx = nFB(j) - nFB(j-1)
      IF (k.le.0 .or. k.gt.kmx) THEN
         Fehler = 'PROG ERR/ MemBlSubDel/ k oor'
         RETURN
         ENDIF

      nold = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)
      IF (npos+ndel.gt.nold) THEN
         Print *, 'j k nold[intern] npos ndel ', j, k, nold, npos, ndel
         Fehler = 'PROG ERR/ MemBlSubDel/ cannot delete so many lines'
         RETURN
         ENDIF

C  Delete :
      DO i = nFA(nFB(j-1)+k-1)+npos+ndel+1, nFA(nFB(nF))
         FMem(i-ndel) = FMem(i)
         ENDDO
      DO kk = nFB(j-1)+k, nFB(nF)
         nFA(kk) = nFA(kk) - ndel
         ENDDO

      END ! MemBlSubDel

      SUBROUTINE MemBlockGet (j, k, n, nmax, R, Fehler)
C     -------------------------------------------------
            ! JWu 17mar95
         ! Read R(1..n) from block k of file j in on-line memory

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       j, k, n, nmax, i
      REAL*8        R(*)
      CHARACTER     Fehler*(*)

      INTEGER       nF, nFB, nFA
      REAL*8        FMem
      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

      IF     (j.le.0 .or. j.gt.nF) THEN
         Fehler = 'PROG ERR/ MemBlockGet/ j oor'
         RETURN
      ELSEIF (k.le.0 .or. k.gt.nFB(j)-nFB(j-1)) THEN
         Fehler = 'PROG ERR/ MemBlockGet/ k oor'
         RETURN
         ENDIF

      n = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)
      IF (n.lt.0) THEN
         CALL Compose2 (Fehler, 'PROG ERR/ MemBlockGet/ k='//
     *                  cl3(k), 'n < 0')
         RETURN
c      ELSEIF (n.eq.0) THEN
c         CALL Compose2 (Fehler, 'PROG ERR/ MemBlockGet/ k='//cl3(k), 'n = 0')
c         RETURN
      ELSEIF (n.gt.nmax) THEN
         CALL Compose4 (Fehler, 'PROG ERR/ MemBlockGet/ file '//cl3(j),
     *      '/ block '//cl4(k), '/ n='//cl4(n), ' > nmax='//cl4(nmax))
         RETURN
         ENDIF
      DO i = 1, n
         R(i) = FMem(nFA(nFB(j-1)+k-1)+i)
         ENDDO

      END ! MemBlockGet

      SUBROUTINE MemBlSubGet (j, k, i0, nsub, R, Fehler)
C     --------------------------------------------------
            ! JWu 17mar95
         ! read entries i0+1...i0+nsub from block k of on-line-file j
         ! into R(1...nsub)

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       j, k, i0, nsub, i, n
      REAL*8        R(*)
      CHARACTER     Fehler*(*)

      INTEGER       nF, nFB, nFA
      REAL*8        FMem

      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

      IF     (j.le.0 .or. j.gt.nF) THEN
         Fehler = 'PROG ERR/ MemBlSubGet/ j oor'
         RETURN
      ELSEIF (k.le.0 .or. k.gt.nFB(j)-nFB(j-1)) THEN
         Print *, ' j, k nK ', j, k, nFB(j)-nFB(j-1)
         Fehler = 'PROG ERR/ MemBlSubGet/ k oor'
         RETURN
         ENDIF

      n = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)
      IF     (i0.lt.0 .or. nsub.lt.1) THEN
         CALL Compose3 (Fehler, 'PROG ERR/ MemBlSubGet/ k='//cl3(k),
     *                  '; i0='//cl6(i0), '; nsub='//cl6(nsub))
         RETURN
      ELSEIF (i0+nsub.gt.n) THEN
         CALL Compose4 (Fehler, 'PROG ERR/ MemBlSubGet/ k='//cl3(k),
     *         '; n='//cl6(n),  '; i0='//cl6(i0), '; nsub='//cl6(nsub))
         RETURN
         ENDIF

      DO i = 1, nsub
         R(i) = FMem(nFA(nFB(j-1)+k-1)+i0+i)
         ENDDO

      END ! MemBlSubGet

      SUBROUTINE MemBlSubPut (j, k, i0, nsub, R, Fehler)
C     --------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       j, k, i0, nsub, i, n
      REAL*8        R(*)
      CHARACTER     Fehler*(*)
      INTEGER       nF, nFB, nFA
      REAL*8        FMem

      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

      IF     (j.le.0 .or. j.gt.nF) THEN
         Fehler = 'PROG ERR/ MemBlSubPut/ j oor'
         RETURN
      ELSEIF (k.le.0 .or. k.gt.nFB(j)-nFB(j-1)) THEN
         Fehler = 'PROG ERR/ MemBlSubPut/ k oor'
         RETURN
         ENDIF

      n = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)
      IF     (i0.lt.0 .or. nsub.lt.1) THEN
         CALL Compose3 (Fehler, 'PROG ERR/ MemBlSubPut/ k='//cl3(k),
     *                  '; i0='//cl6(i0), '; nsub='//cl6(nsub))
         RETURN
      ELSEIF (i0+nsub.gt.n) THEN
         CALL Compose4 (Fehler, 'PROG ERR/ MemBlSubPut/ k='//cl3(k),
     *         '; n='//cl6(n),  '; i0='//cl6(i0), '; nsub='//cl6(nsub))
         RETURN
         ENDIF

      DO i = 1, nsub
         FMem(nFA(nFB(j-1)+k-1)+i0+i) = R(i)
         ENDDO

      END ! MemBlSubPut

      INTEGER FUNCTION MemBlockInq (what)
C     -----------------------------------
            ! JWu 1feb93, 11mar93, 16may95
         ! inquire memory status

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      CHARACTER     what*(*)

      INTEGER       nF, nFB, nFA
      REAL*8        FMem

      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

      IF     (what.eq.'nF') THEN
         MemBlockInq = nF
      ELSEIF (what.eq.'MF') THEN
         MemBlockInq = MF
      ELSEIF (what.eq.'fF') THEN ! free files
         MemBlockInq = MF - nF
c      ELSEIF (what.eq.'nE') THEN
c         MemBlockInq = iOnext(nF) - 1
      ELSEIF (what.eq.'ME') THEN
         MemBlockInq = Mmem
      ELSEIF (what.eq.'fE') THEN ! free entries
         MemBlockInq = Mmem - nFA(nFB(nF)) ! ????
      ELSE
         CALL Absturz('MemBlockInq', 'Option not implemented :'//what)
         ENDIF

      END ! MemBlockInq

      INTEGER FUNCTION MemBlockNum (j)
C     --------------------------------
            ! 17may95
         ! return number of blocks of file j

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       nF, nFB, nFA, j
      REAL*8        FMem

      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

      IF (j.lt.1 .or. j.gt.nF) THEN
         MemBlockNum = 0 ! must be tested in the calling routing
         RETURN
         ENDIF
      MemBlockNum = nFB(j) - nFB(j-1)

      END ! MemBlockNum

      INTEGER FUNCTION MemBlockSiz (j, k)
C     -----------------------------------
            ! 17may95
         ! return number of lines of block k of file j

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       nF, nFB, nFA, j, k
      REAL*8        FMem

      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF

      IF     (j.lt.0 .or. j.gt.nF) THEN
         MemBlockSiz = -1
c         type *, ' WARNUNG/ MemBlockSiz/ j oor'
         RETURN
      ELSEIF (k.le.0 .or. k.gt.nFB(j)-nFB(j-1)) THEN
         MemBlockSiz = 0 ! war -1, aber warum nicht ganz normal 0 abfragen ?
c         type *, 'UNGLAUBLICH/ MemBlockSiz/ k oor : ', k
         RETURN
         ENDIF

      MemBlockSiz = nFA(nFB(j-1)+k) - nFA(nFB(j-1)+k-1)

      END ! MemBlockSiz

      SUBROUTINE MemRestInfo ()
C     -------------------------

      IMPLICIT NONE
      INTEGER       iRest, iAlt, MCmem, MemBlockInq
      CHARACTER     cl8*8

      iRest = MemBlockInq ('fE')
      MCmem = MemBlockInq ('ME')
      IF (iRest.ne.iAlt .and. iRest.le.MCmem/7) THEN
         CALL Say2 (' '//cl8(iRest), ' lines free in run-time memory')
         iAlt = iRest
         ENDIF

      END ! MemRestInfo

      SUBROUTINE MemDims ()
C     ---------------------
            ! JWu FileDims 23nov92
         ! Print current array dimensions

      INCLUDE 'i_dim.f'

      Print '(a)',
     *   ' Current array dimensions (data files and on-line memory) :'
      Print '(a,i8)', ' # files             ', MF
      Print '(a,i8)', ' # spectra/file      ', MK
      Print '(a,i8)', ' # channels/spectrum ', MC
      Print '(a,i8)', ' # tPar/file         ', MP
      Print '(a,i8)', ' # total channels    ', Mmem

      END ! MemDims

      LOGICAL FUNCTION qMemFill (newF, newE, Fehler)
C     ----------------------------------------------
            ! 12mar93
         ! Is there not enough place to accomodate
         ! newF new files with a total of newE entries ?

      IMPLICIT LOGICAL (q)
      CHARACTER         Fehler*(*), cl8*8

      nF = MemBlockInq ('fF')
      nE = MemBlockInq ('fE')
      IF     (newF.gt.nF) THEN
         CALL Compose3 (Fehler,
     *      ' not enough place in memory : '//cl8(nF),
     *      ' more files allowed, '//cl8(newF), ' required')
         qMemFill = .true.
      ELSEIF (newE.gt.nE) THEN
         CALL Compose3 (Fehler,
     *      ' not enough place in memory : '//cl8(nE),
     *      ' lines free, '//cl8(newE), ' required')
         qMemFill = .true.
      ELSE
         qMemFill = .false.
         ENDIF

      END ! qMemFill

      SUBROUTINE MemFileDel (j, Fehler)
C     ---------------------------------
         ! 16may95, 17jan91

      ! Delete one file from OLM.

      IMPLICIT NONE
      INCLUDE 'l_def.f'

      CHARACTER     Fehler*(*)
      INTEGER       j, n, iAlt, iNew, MemBlockInq
      REAL*8        R(1)

      iAlt = MemBlockInq ('fE')
      CALL MemBlockPut (j, 0, n, R, Fehler)
      iNew = MemBlockInq ('fE')

      CALL Say3 ('  file '//cl3(j), ' deleted ('//
     *           cl6(iNew-iAlt), ' lines)')
      CALL MemRestInfo ()

      END ! MemFileDel

      SUBROUTINE MemFileDup (jin, jout, nBdup, Fehler)
C     ------------------------------------------------
            ! JWu 7nov96
         ! duplicate the first nBdup blocks of an on-line-file
         ! without any notion of its semantics, preserving the status bits.
         ! nBdup=0 means: copy *all* blocks

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      INCLUDE 'i_dim.f'

      CHARACTER     Fehler*(*)
      INTEGER       jin, jout, nB, nBdup, iB, n, iAlt, iNew,
     *              MemBlockNum, MemBlockInq
      REAL*8        R(MmemBlo)

      iAlt = MemBlockInq ('fE')

      nB = MemBlockNum (jin)
      IF (nB.lt.1) THEN
         Fehler = 'MemFileDup/ empty or nonexisting file'
         RETURN
         ENDIF
      IF (nBdup.ge.1) THEN
         IF (nB.lt.nBdup) THEN
            Fehler = 'MemFileDup/ too many blocks required'
            RETURN
            ENDIF
         nB = nBdup
         ENDIF
      jout = 0

      DO iB = 1, nB
         CALL MemBlockGet (jin, iB, n, MmemBlo, R, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL MemBlockPut (jout, iB, n, R, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      iNew = MemBlockInq ('fE')
      CALL Say4 ('  file '//cl3(jin), ' -> '//cl3(jout),
     *           ' duplicated ('//cl6(iAlt-iNew), ' lines)')
      CALL MemRestInfo ()

      END ! MemFileDup

C  ====================================================================
C     2.  Status shell : access to files restricted by their status
C  ====================================================================

         ! The status bits (introduced 25jan95) encode in particular
         ! information about the allowed access mode.

         ! bit (0) = 0 / 1 : online / being modified
         ! bit (1) = 0 / 2 : as on disc / modified
         ! bit (2) = 0 / 4 : changeable / protected

      SUBROUTINE MemFileStatP (j, iStat, Fehler)
C     ------------------------------------------
            ! JWu 16may95 reusing OlfOpen (23jan91)

      ! Open or close a file for modifications

      IMPLICIT NONE
      INCLUDE       'l_def.f'

      INTEGER        j, iStat, iR
      REAL*8         R(1)
      CHARACTER*(*)  Fehler
      EQUIVALENCE   (R(1), iR)

      IF (iStat.lt.0 .or. iStat.ge.8) THEN
         Fehler = 'MFSP/ Illegal status bits'
         RETURN
         ENDIF

      IF (iBitGet(iStat,0).eq.1 .and. iBitGet(iStat,2).eq.1
     *    .and. j.ne.0) THEN
         Fehler = 'MFSP/ Attempt to modify protected file'
         RETURN
         ENDIF

      iR = iStat
      CALL MemBlockPut (j, 1, 1, R, Fehler)
      IF (Fehler.ne.'&ff') Print *, 'error passed through MFSP'

      END ! MemFileStatP

      SUBROUTINE MemFileStatG (j, iStat, Fehler)
C     ------------------------------------------

      ! Get status information about file j

      IMPLICIT NONE
      INCLUDE       'l_def.f'

      INTEGER        j, iStat, iR, n
      REAL*8         R(1)
      CHARACTER*(*)  Fehler
      EQUIVALENCE   (R(1), iR)

      CALL MemBlockGet (j, 1, n, 1, R, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through MFSG'
         RETURN
         ENDIF

      IF (n.ne.1) THEN
         Fehler =
     * 'MFSG/ Block 1 did not contain 1 R*8 but #elements='//cl6(n)
         RETURN
         ENDIF
      iStat = iR

      IF (iStat.lt.0 .or. iStat.ge.8) THEN
         Fehler = 'MFSG/ Illegal status bits'
         RETURN
         ENDIF

      END ! MemFileStatG

      SUBROUTINE MemFileBitP (j, nPos, iVal, Fehler)
C     ----------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*)  Fehler
      INTEGER        j, nPos, iVal, iStat

      CALL MemFileStatG (j, iStat, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL BitSet (iStat,nPos,iVal)

      CALL MemFileStatP (j, iStat, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! MemFileBitP

      SUBROUTINE MemBlockOvr (j, k, n, R, Fehler)
C     -------------------------------------------
         ! Save R(1..n) as block k of file j in on-line memory
         ! but only if the present status of file j allows it

         ! From outer levels, calls to MemBlockPut should always
         ! pass through MemBlockOvr

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       j, k, n, iStat
      REAL*8        R(*)
      CHARACTER     Fehler*(*)

      CALL MemFileStatG (j, iStat, Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

      IF (iBitGet(iStat,0).eq.0) THEN
         CALL Compose2 (Fehler, 'MBO/ File '//cl3(j),
     *                  ' is not open for modifications')
         RETURN
         ENDIF

      IF (iBitGet(iStat,1).eq.0) THEN
         CALL BitSet(iStat,1,1) ! file will be different from version on disc
         CALL MemFileStatP (j, iStat, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         ENDIF

      CALL MemBlockPut (j, k, n, R, Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

      RETURN
 99   CONTINUE
      Print *, 'error passed through MBO'

      END ! MemBlockOvr

C  ====================================================================
C     3. Open, close, and copy files
C  ====================================================================

      ! Block(1)   : status      MFStat (see above)
      ! Block(2)   : integer parameters
      ! Block(3)   : real parameters
      ! Block(4)   : text parameters
      ! Block(5)   : documentation
      ! Block(6)   : coordinate names
      ! Block(7..) : four blocks for each spectrum (z X Y D)

      ! For modifying an on-line-file, proceed as follows
      ! (new organisation from nov96) :

      ! for creating a new file starting from nothing, start with
      !    CALL OlfCreate (jout, Kout, FileDef, TitlDef, Fehler)
      ! for modifying an existing file, use
      !    CALL OlfHeadDup (jin, qOv, jout, Kout, Fehler)
      ! except it shall always be overwritten, in which case you use
      !    CALL OlfOpen (j, 1, Kout)          % if j=0, new j on exit


      ! * old procedure : *
      !    CALL OlfOpen (j, 1, Kout)          % if j=0, new j on exit
      !    CALL OlfParP (j, iPar, rPar, tPar) % required if file is new
      !    Loop
      !       CALL OlfSpeP (j, Kout,...)      % Kout+=1; overwrite spectrum
      !       IF Fehler GOTO 1
      !       LoopEnd
      !  1 CALL OlfClos (j, Kout)             % Kout is new nK

      SUBROUTINE OlfOpen (j, iStat, Kout, Fehler)
C     -------------------------------------------
               ! JWu FileModif 23.1.91, 12.2.91, OlfOpen 17may95

         ! open an existing on-line-file for modification

      IMPLICIT NONE
      CHARACTER     Fehler*(*)
      INTEGER       j, iStat, Kout

      IF (j.le.0) THEN
         Fehler = 'PROG ERR/ OlfOpen now restricted to existing files'
         RETURN
         ENDIF
      Kout  = 0 ! Counter for modified spectra
      CALL MemFileStatP (j, iStat, Fehler)
      IF (Fehler.ne.'&ff') CALL Insert (Fehler, 1, 'OO/ ')

      END ! OlfOpen

      SUBROUTINE OlfClos (j, Kout, Fehler)
C     ------------------------------------

      IMPLICIT NONE
      INCLUDE 'l_def.f'
      CHARACTER     Fehler*(*)
      INTEGER       j, Kout

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in OlfClos'
         RETURN
         ENDIF

      IF     (Kout.ge.2) THEN
         CALL Say2 ('  stored '//cl3(Kout), ' spectra as file '//
     *              cl3(j))
      ELSEIF (Kout.ge.1) THEN
         Print *, '  stored spectrum as file ', j
         ENDIF
      CALL MemRestInfo()

      CALL MemFileBitP (j, 0, 0, Fehler) ! close for modif
      IF (Fehler.ne.'&ff') Print *, 'error passed through OlfClos'

      END ! OlfClos

      SUBROUTINE OlfHeadDup (jin, qOv, jout, nK, Kout, Fehler)
C     --------------------------------------------------------
            ! JWu 8nov96
         ! qOv=true:  open existing file for modification
         ! qOv=false: copy header blocks to new file and open
         !            the new file for modification

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      CHARACTER     Fehler*(*)
      INTEGER       jin, jout, Kout, nK, nB, MemBlockNum
      LOGICAL       qOv

      IF (Fehler.ne.'&ff') THEN
         Print *, ' OlfHeadDup/ error on entry'
         RETURN
         ENDIF

      nB = MemBlockNum(jin)
      IF (nB.lt.1) THEN
         Print *, 'jin = ', jin
         Fehler = 'OlfHeadDup/ input file does not exist'
         RETURN
      ELSEIF (nB.le.MBH) THEN
         Fehler = 'OlfHeadDup/ input File has uncomplete header'
         RETURN
         ENDIF
      nK = (nB - MBH) / 4
      Kout = 0
      IF (qOv) THEN
         jout = jin
      ELSE
         jout = 0
         CALL MemFileDup (jin, jout, MBH, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF
      CALL MemFileStatP (jout, 1, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! OlfHeadDup

      SUBROUTINE OlfCreate (jout, Kout, FilDef, TitDef, Fehler)
C     ---------------------------------------------------------
            ! JWu 8nov96
         ! qOv=true:  open existing file for modification
         ! qOv=false: copy header blocks to new file and open
         !            the new file for modification

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      CHARACTER*(*) FilDef, TitDef, Fehler
      CHARACTER*40  Fil, Tit
      INTEGER       jout, Kout, iBH, MemBlockNum
      REAL*8        R(1)
      SAVE          Fil, Tit
      DATA          Fil /' '/, Tit /' '/

      IF (Fehler.ne.'&ff') THEN
         Print *, ' OlfCreate/ error on entry'
         RETURN
         ENDIF

      ! create the file and open it for modifications :
      jout = 0 ! ignore input value
      Kout = 0 ! Counter for modified spectra
      CALL MemFileStatP (jout, 1, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      ! initialize the header blocks :
      DO iBH = 2, 6
         CALL MemBlockPut (jout, iBH, 0, R, Fehler)
         ENDDO

      ! set filename and title :
      IF     (FilDef.eq.'&noask') THEN
         Fil = ' '
      ELSEIF (FilDef.eq.'&nodef') THEN
         CALL FrageT  ('File name', Fil)
      ELSEIF (FilDef.eq.'&olddef') THEN ! last input -> new default
         CALL FrageTD ('File name', Fil, Fil)
      ELSE
         CALL FrageTD ('File name', Fil, FilDef)
         ENDIF
      IF     (TitDef.eq.'&noask') THEN
         Tit = ' '
      ELSEIF (TitDef.eq.'&nodef') THEN
         CALL FrageT  ('and Title', Tit)
      ELSEIF (TitDef.eq.'&olddef') THEN ! last input -> new default
         CALL FrageTD ('and Title', Tit, Tit)
      ELSE
         CALL FrageTD ('and Title', Tit, TitDef)
         ENDIF
      CALL tOlfP (jout, 'fil', Fil, Fehler)
      CALL tOlfP (jout, 'tit', Tit, Fehler)
      CALL tOlfP (jout, 'doc', ' ', Fehler)
      CALL tOlfP (jout, 'dir', ' ', Fehler)

      END ! OlfCreate

C  ====================================================================
C     4. Access to header blocks
C  ====================================================================

      SUBROUTINE OlfLabG (j, k, lab, RL, nRL, Fehler)
C     -----------------------------------------------
            ! JWu 23/24nov96
         ! get any labelled parameter from block k
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) lab, Fehler
      CHARACTER     name*24
      INTEGER       j, k, nRL, MRlab, MRmax, nR, i, iP, nP, MemBlockSiz
      PARAMETER     (MRlab=24/8, MRmax=MRlab+80/8)
      REAL*8         R(MRmax), RL(*)
      EQUIVALENCE   (R(1), name)

      nR = nRL + MRlab
      IF (nR.gt.MRmax) THEN
         Fehler = 'OlfLabG/ nR>MRmax'
         RETURN
         ENDIF
      nP = MemBlockSiz(j, k) / nR
      IF (nP.lt.0) THEN
         Fehler = 'OlfLabG/ block not accessible'
         ENDIF

      IF (lab(1:4).eq.'&pbn') THEN ! parameter by number
         CALL Fi1N (lab, iP)
         IF (lab.ne.'&pbn #') THEN
            Fehler = 'OlfLabG/ invalid Macro '//lab
            RETURN
         ELSEIF (iP.lt.1) THEN
            Fehler = 'OlfLabG/ invalid iP<1'
            RETURN
         ELSEIF (iP.eq.nP+1) THEN
            lab = '&eop' ! end of parameters
            RETURN
         ELSEIF (iP.gt.nP+1) THEN
            Fehler = 'OlfLabG/ invalid iP>>nP'
            RETURN
            ENDIF
         CALL MemBlSubGet (j, k, (iP-1)*nR, nR, R, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         lab = name
         DO i = 1, nRL
            RL(i) = R(MRlab+i)
            ENDDO
         RETURN
         ENDIF

      DO iP = 0, nP-1
         CALL MemBlSubGet (j, k, iP*nR, nR, R, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (name.eq.lab) THEN
            DO i = 1, nRL
               RL(i) = R(MRlab+i)
               ENDDO
            RETURN ! success
            ENDIF
         ENDDO

      ! parameter-not-found :
      CALL Compose3 (Fehler, '&pnf ('//cl2(k), ') "'//lab,
     *   '" by OlfLabG in j='//cl3(j))

      END ! OlfLabG

      SUBROUTINE OlfLabDel (j, k, lab, nRL, Fehler)
C     ---------------------------------------------
            ! JWu 3/4feb97
         ! remove any labelled parameter from block k
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) lab, Fehler
      CHARACTER     name*24
      INTEGER       j, k, nRL, MRlab, MRmax, nR, i, iP, nP, MemBlockSiz
      PARAMETER     (MRlab=24/8, MRmax=MRlab+80/8)
      REAL*8         R(MRmax)
      EQUIVALENCE   (R(1), name)

      nR = nRL + MRlab
      nP = MemBlockSiz(j, k) / nR

      IF (nP.lt.0) THEN
         Fehler = 'OlfLabDel/ block not accessible'
         ENDIF

      DO iP = 0, nP-1
         CALL MemBlSubGet (j, k, iP*nR, nR, R, Fehler)
         IF (name.eq.lab) THEN
            CALL MemBlSubDel (j, k, iP*nR, nR, Fehler)
            RETURN ! success
            ENDIF
         ENDDO

      ! parameter-not-found :
      CALL Compose3 (Fehler, '&pnf ('//cl2(k), ') "'//
     *               lab, '" in OlfLabDel')

      END ! OlfLabDel

      SUBROUTINE OlfLabP (j, k, lab, RL, nRL, Fehler)
C     -----------------------------------------------
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) lab, Fehler
      CHARACTER     name*24
      INTEGER       j, k, nRL, MRlab, MRmax, nR, i, iP, nP, nPadd,
     *              MemBlockSiz
      PARAMETER     (MRlab=24/8, MRmax=MRlab+80/8)
      REAL*8         R(MRmax), RL(*)
      EQUIVALENCE   (R(1), name)

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in OlfLabP'
         RETURN
         ENDIF
      nR = nRL + MRlab
      IF (nR.gt.MRmax) THEN
         Fehler = 'OlfLabP/ nR>MRmax'
         RETURN
         ENDIF
      nP = MemBlockSiz(j, k) / nR
      IF (nP.lt.0) THEN
         Fehler = 'tOlfP/ block not accessible'
         RETURN
         ENDIF

      DO iP = 0, nP-1
         CALL MemBlSubGet (j, k, iP*nR, nR, R, Fehler)
         IF (Fehler.ne.'&ff') THEN
            CALL Insert (Fehler, 1, 'OlfLabP (get-old)/ ')
            RETURN
            ENDIF
         IF (name.eq.'&empty') name = lab ! take free slot
         IF (name.eq.lab) THEN
            DO i = 1, nRL
               R(MRlab+i) = RL(i)
               ENDDO
            CALL MemBlSubPut (j, k, iP*nR, nR, R, Fehler)
            IF (Fehler.ne.'&ff') THEN
               CALL Insert (Fehler, 1, 'OlfLabP (put-ovr) / ')
               RETURN
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            RETURN ! overwrite
            ENDIF
         ENDDO

 2    CONTINUE ! no more slots free

      nPadd = nP + 8
      CALL MemBlSubAdd (j, k, nP*nR, nPadd*nR, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'OlfLabP (add many)/ ')
         RETURN
         ENDIF
      name = lab
      DO i = 1, nRL
         R(MRlab+i) = RL(i)
         ENDDO
      CALL MemBlSubPut (j, k, nP*nR, nR, R, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'OlfLabP (put-new)/ ')
         RETURN
         ENDIF

      name = '&empty'
      DO i = 1, nRL
         R(MRlab+i) = 0
         ENDDO
      DO iP = nP+1, nP+nPadd-1
         CALL MemBlSubPut (j, k, iP*nR, nR, R, Fehler)
         IF (Fehler.ne.'&ff') THEN
            CALL Insert (Fehler, 1, 'OlfLabP (put-void)/ ')
            RETURN
            ENDIF
         ENDDO

      END ! OlfLabP

      INTEGER FUNCTION iOlfG (j, lab, Fehler)
C     ---------------------------------------
         ! return one integer parameter from the header of file j.
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      CHARACTER*(*) lab, Fehler
      INTEGER       j, nB, ival, MemBlockNum, MemBlockSiz, nZ
      REAL*8        R(1)
      EQUIVALENCE  (R(1), ival)

      iOlfG = 0

      IF (lab.eq.'#spectra') THEN
         nB = MemBlockNum(j)
         IF (nB.lt.1) THEN
            Fehler = 'iOlfG/ File does not exist'
            RETURN
         ELSEIF (nB.le.MBH) THEN
            Fehler = 'iOlfG/ File has uncomplete header'
            RETURN
            ENDIF
         iOlfG = (nB - MBH) / 4 ! bei "Anderung auch lokal: nK = .. "andern
      ELSEIF (lab.eq.'#Z') THEN
         nZ = MemBlockSiz(j, MBH+1)
         IF (nZ.lt.0) THEN
            Fehler = 'iOlfG/ PROGRAM ERROR/ MemBlockSize=-1'
         ELSE
            iOlfG = nZ
            ENDIF
      ELSE
         CALL OlfLabG (j, 3, lab, R, 1, Fehler)
         iOlfG = ival
         ENDIF

      END ! iOlfG

      SUBROUTINE iOlfP (j, lab, iin, Fehler)
C     --------------------------------------
         ! overwrite one integer parameter in the header of file j.
      IMPLICIT NONE
      CHARACTER*(*) lab, Fehler
      INTEGER       j, iin, ival, MemBlockNum
      REAL*8        R(1)
      EQUIVALENCE  (R(1), ival)

      IF (lab.eq.'#spectra') THEN
         Fehler = 'PROGR ERR/ cannot modify #spectra via iOlfP'
         RETURN
         ENDIF

      ival = iin
      CALL OlfLabP (j, 3, lab, R, 1, Fehler)

      END ! iOlfP

      REAL*8 FUNCTION rOlfG (j, co, un, Fehler)
C     -----------------------------------------
         ! read one real parameter from the header of file j.
      IMPLICIT NONE
      CHARACTER*(*) co, un, Fehler
      CHARACTER*24  uval
      INTEGER       j
      REAL*8        rval, R(4)
      EQUIVALENCE  (R(1), rval)
      EQUIVALENCE  (R(2), uval)

      CALL OlfLabG (j, 4, co, R, 4, Fehler)
      IF (Fehler.ne.'&ff') THEN
         rOlfG = 0
         RETURN
         ENDIF
      rOlfG = rval
      un    = uval

      END ! rOlfG

      SUBROUTINE rOlfDel (j, co, Fehler)
C     ----------------------------------
            ! JWu 3feb97
         ! delete one real parameter from the header of file j.
      IMPLICIT NONE
      CHARACTER*(*) co, Fehler
      INTEGER       j

      CALL OlfLabDel (j, 4, co, 4, Fehler)

      END ! rOlfDel

      SUBROUTINE rOlfP (j, co, un, rin, Fehler)
C     -----------------------------------------
            ! JWu 9nov96
         ! overwrite one real parameter in the header of file j.

      IMPLICIT NONE
      CHARACTER*(*) co, un, Fehler
      CHARACTER*24  uval
      INTEGER       j
      REAL*8        rin, rval, R(4)
      EQUIVALENCE  (R(1), rval)
      EQUIVALENCE  (R(2), uval)

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in rOlfP'
         RETURN
         ENDIF

      rval = rin
      uval = un
      CALL OlfLabP (j, 4, co, R, 4, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'rOlfP'
         RETURN
         ENDIF

      END ! rOlfP

      SUBROUTINE tOlfG (j, lab, text, Fehler)
C     ---------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) lab, text, Fehler
      CHARACTER     line*80
      INTEGER       j
      REAL*8        R(8)
      EQUIVALENCE  (R(1), line)

      CALL OlfLabG (j, 2, lab, R, 8, Fehler)
      text = line

      END ! tOlfG

      SUBROUTINE tOlfP (j, lab, text, Fehler)
C     ---------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) lab, text, Fehler
      CHARACTER     line*80
      INTEGER       j
      REAL*8        R(8)
      EQUIVALENCE  (R(1), line)

      line = text
      CALL OlfLabP (j, 2, lab, R, 8, Fehler)

      END ! tOlfP

C  --------------------------------------------------------------------
C     4.3 user interface / access to header block 5 (coordinate names)
C  --------------------------------------------------------------------

      SUBROUTINE OlfCnuG (j, lab, co, un, Fehler)
C     -------------------------------------------
            ! JWu  8nov96 using old tPar
            ! JWu 24nov96 direct access to block
         ! search coordinate-name-and-unit

      IMPLICIT NONE
      CHARACTER*(*) lab, co, un, Fehler
      CHARACTER*24  cval, uval
      INTEGER       j
      REAL*8        R(6)
      EQUIVALENCE  (R(1),cval)
      EQUIVALENCE  (R(4),uval)

      CALL OlfLabG (j, 5, lab, R, 6, Fehler)
      co = cval
      un = uval
      IF (Fehler.ne.'&ff') Print *, 'error passed through OlfCnuG ('//
     *                             lab//')'

      END ! OlfCnuG

      SUBROUTINE OlfCnuP (j, lab, co, un, Fehler)
C     -------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) lab, co, un, Fehler
      CHARACTER*24  cval, uval
      CHARACTER     cl2*2
      INTEGER       j, nZ, iOlfG
      REAL*8        R(6)
      EQUIVALENCE  (R(1),cval)
      EQUIVALENCE  (R(4),uval)

      IF (Fehler.ne.'&ff') THEN
         Print *, 'BAD PROGRAMMING STYLE/ error on entry in OlfCnuP'
         RETURN
         ENDIF
      IF (co.eq.' ') THEN
         CALL OlfLabDel (j, 5, lab, 6, Fehler)
         IF (Fehler.ne.'&ff') Print *,
     *       'error passed through OlfCnuP ('' '')'
      ELSE
         cval = co
         uval = un
         IF (lab.eq.'z+') THEN
            nZ = iOlfG (j, '#Z', Fehler)
            IF (Fehler.ne.'&ff') THEN
               Print *, 'OlfCnuP/ error in iOlfG'
               RETURN
               ENDIF
            CALL OlfLabP (j, 5, 'z'//cl2(nZ+1), R, 6, Fehler)
            IF (Fehler.ne.'&ff') Print *,
     *          'error passed through OlfCnuP (z+)'
         ELSE
            CALL OlfLabP (j, 5, lab, R, 6, Fehler)
            IF (Fehler.ne.'&ff')
     *           Print *, 'error passed through OlfCnuP ('//lab//')'
            ENDIF
         ENDIF

      END ! OlfCnuP

C  --------------------------------------------------------------------
C     4.4 user interface / access to header block 6 (comment)
C  --------------------------------------------------------------------

      SUBROUTINE OlfComAddFull (j, adoc, acom, idat, Fehler)
C     ------------------------------------------------------
            ! renewed JWu 8/24nov96
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) adoc, acom, Fehler
      CHARACTER*80  doc, line
      CHARACTER     sz*8, Zeit*8, sd*10, Datum*10
      INTEGER       j, nP, MR, MemBlockSiz, idat
      PARAMETER     (MR=80/8)
      REAL*8         R(MR)
      EQUIVALENCE   (R(1), line)

      IF (adoc.ne.' ') THEN
         CALL tOlfG (j, 'doc', doc, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL Append (doc, adoc)
         CALL tOlfP (j, 'doc', doc, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF

      nP = MemBlockSiz(j, 6) / MR
      IF (nP.lt.0) THEN
         Fehler = 'OlfAddCom/ block not accessible'
         ENDIF

      CALL MemBlSubAdd (j, 6, nP*MR, MR, Fehler)
      IF (idat.ge.1) THEN
         sd = Datum(7)
         sz = Zeit(5)
         line = sd(1:7)//' '//sz(1:5)//' '//acom
      ELSE
         line = acom
         ENDIF
      CALL MemBlSubPut (j, 6, nP*MR, MR, R, Fehler)

      END ! OlfComAddFull

      SUBROUTINE OlfComAdd (j, adoc, acom, Fehler)
C     --------------------------------------------
      CHARACTER*(*) adoc, acom, Fehler
      CALL OlfComAddFull (j, adoc, acom, 0, Fehler)
      END ! OlfComAdd

      SUBROUTINE OlfComLinP (j, iP, lin, Fehler)
C     -------------------------------------------
         ! only for use in EditDoc
      IMPLICIT NONE
      CHARACTER*(*) Fehler, lin
      CHARACTER     line*80
      INTEGER       j,  nP, iP, MR, MemBlockSiz
      PARAMETER    (MR=80/8)
      REAL*8        R(MR)
      EQUIVALENCE  (R(1), line)

      nP = MemBlockSiz(j, 6) / MR
      IF (nP.lt.0) THEN
         Fehler = 'OlfComLinP/ block not accessible'
         RETURN
         ENDIF

      IF (iP.lt.1 .or. iP.gt.nP+1) THEN
         Fehler = 'OlfComLinP/ iP oor'
      ELSE
         IF (iP.eq.nP+1) CALL MemBlSubAdd (j, 6, nP*MR, MR, Fehler)
         line = lin
         CALL MemBlSubPut (j, 6, (iP-1)*MR, MR, R, Fehler)
         ENDIF

      END ! OlfComLinP

      SUBROUTINE OlfComLinG (j, iP, lout, Fehler)
C     -------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) Fehler, lout
      CHARACTER     line*80
      INTEGER       j,  nP, iP, MR, MemBlockSiz
      PARAMETER    (MR=80/8)
      REAL*8        R(MR)
      EQUIVALENCE  (R(1), line)

      nP = MemBlockSiz(j, 6) / MR
      IF (nP.lt.0) THEN
         Fehler = 'OlfComLinG/ block not accessible'
         RETURN
         ENDIF

      IF (iP.eq.nP+1) THEN
         lout = '&eoc' ! end-of-comment
      ELSEIF (iP.lt.1 .or. iP.gt.nP+1) THEN
         Fehler = 'OlfComLinG/ iP oor'
      ELSE
         CALL MemBlSubGet (j, 6, (iP-1)*MR, MR, R, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         lout = line
         ENDIF

      END ! OlfComLinG

      SUBROUTINE OlfComLinDel (j, iP, Fehler)
C     ---------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) Fehler
      INTEGER       j,  nP, iP, MR, MemBlockSiz
      PARAMETER    (MR=80/8)

      nP = MemBlockSiz(j, 6) / MR
      IF (nP.lt.0) THEN
         Fehler = 'OlfComLinDel/ block not accessible'
         RETURN
         ENDIF

      IF (iP.lt.1 .or. iP.gt.nP) THEN
         Fehler = 'OlfComLinG/ iP oor'
      ELSE
         CALL MemBlSubDel (j, 6, (iP-1)*MR, MR, Fehler)
         ENDIF

      END ! OlfComLinDel

      SUBROUTINE OlfComLinDelAll (j, Fehler)
C     --------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) Fehler
      INTEGER       j,  nP, iP, MR, MemBlockSiz
      PARAMETER    (MR=80/8)

      nP = MemBlockSiz(j, 6) / MR
      IF (nP.lt.0) THEN
         Fehler = 'OlfComLinDel/ block not accessible'
         RETURN
         ENDIF

      DO iP = nP, 1, -1
         CALL MemBlSubDel (j, 6, (iP-1)*MR, MR, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      END ! OlfComLinDelAll

C  --------------------------------------------------------------------
C     4.Appendix / indirect access
C  --------------------------------------------------------------------

      SUBROUTINE pcOlfFind (j, sel, co, un, pc, no, Fehler)
C     -----------------------------------------------------
            ! JWu 19nov98
         ! search for co in header (selected by sel) blocks of file j
         ! returns un, pc, sel

      IMPLICIT NONE
      INCLUDE    'l_def.f'
      INTEGER     j, no, is, iZ, nZ, iOlfG
      CHARACTER   sel*(*), co*(*), un*(*), pc*(*), Fehler*(*),
     *            coI*40, cs*1

      DO is = 1, lenU(sel)
         cs = sel(is:is)
         IF (cs.eq.'z') THEN
            nZ = iOlfG (j, '#Z', Fehler)
            DO iZ = 1, nZ
               CALL OlfCnuG (j, 'z'//cl2(iZ), coI, un, Fehler)
               IF (Fehler.ne.'&ff') THEN
                  CALL Insert (Fehler, 1, 'pcOlfFind/')
                  RETURN
                  ENDIF
               IF (coI.eq.co) THEN ! yes, it is 'z'
                  pc = 'z'
                  no = iZ
                  RETURN
                  ENDIF
               ENDDO

         ELSE
            Fehler = ' this search option not yet implemented'
            RETURN
            ENDIF
         ENDDO

      pc = '-'
      no =  0
      Fehler = '&pnf'

      END ! pcOlfFind

      LOGICAL FUNCTION qOlfG (j, lab, Fehler)
C     ---------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) lab, Fehler
      LOGICAL       qintr
      INTEGER       iOlfG, j, ival

      ival  = iOlfG(j, lab, Fehler)
      qOlfG = qintr(ival)

      END ! qOlfG

      INTEGER FUNCTION iOlfGdef (j, lab, idef, Fehler)
C     ------------------------------------------------
      CHARACTER*(*) lab, Fehler
      iOlfGdef = iOlfG (j, lab, Fehler)
      IF (Fehler(1:4).eq.'&pnf') THEN
         iOlfGdef = idef
         Fehler = '&ff'
         ENDIF
      END ! iOlfGdef

      LOGICAL FUNCTION qOlfGdef (j, lab, idef, Fehler)
C     ------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) lab, Fehler
      LOGICAL       qintr
      INTEGER       iOlfGdef, j, ival, idef

      ival = iOlfGdef (j, lab, idef, Fehler)
      qOlfGdef = qintr(ival)

      END ! qOlfGdef

      REAL*8 FUNCTION rOlfGdef (j, co, un, rdef, Fehler)
C     --------------------------------------------------
         ! read one real parameter from the header of file j.
      IMPLICIT NONE
      CHARACTER*(*) co, un, Fehler
      INTEGER       j
      REAL*8        rdef, rOlfG

      rOlfGdef = rOlfG (j, co, un, Fehler)
      IF (Fehler(1:4).eq.'&pnf') THEN
         rOlfGdef = rdef
         Fehler = '&ff'
         ENDIF

      END ! rOlfGdef

      REAL*8 FUNCTION rzOlfG (j, K, co, un, Fehler)
C     ---------------------------------------------
            ! JWu 9nov96
         ! co is either 'z', or a parameter

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'
      CHARACTER*(*) co, un, Fehler
      CHARACTER*40  coI
      INTEGER       j, K, nZ, iZ, iOlfG
      REAL*8        rval, rOlfG

      rzOlfG = 0
      nZ = iOlfG (j, '#Z', Fehler)
      DO iZ = 1, nZ
         CALL OlfCnuG (j, 'z'//cl2(iZ), coI, un, Fehler)
         IF (Fehler.ne.'&ff') THEN
            Print *, 'rzOlfG[CnuG]'
            RETURN
            ENDIF
         IF (coI.eq.co) THEN ! yes, it is 'z'
            CALL OlfGet1Z (j, K, iZ, rval, Fehler)
            IF (Fehler.ne.'&ff') THEN
               Print *, 'rzOlfG[1Z]'
               RETURN
               ENDIF
            rzOlfG = rval
            RETURN
            ENDIF
         ENDDO
      rzOlfG = rOlfG (j, co, un, Fehler) ! either it's a param, or Fehler..
      IF (Fehler.ne.'&ff' .and. Fehler(1:4).ne.'&pnf')
     *     Print *, 'rzOlfG -> rOlfG'

      END ! rzOlfG

      REAL*8 FUNCTION rzxOlfG (j, K, i, co, un, Fehler)
C     -------------------------------------------------
            ! JWu 9nov96
         ! co is either 'z', or 'x', or a parameter

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'
      CHARACTER*(*) co, un, Fehler
      CHARACTER*40  coI
      INTEGER       j, K, i, n, nZ, iZ, iOlfG
      REAL*8        rval, rOlfG, X(MC), Y(MC)

      rzxOlfG = 0
      CALL OlfCnuG (j, 'x', coI, un, Fehler)
cdeb      Print *, '1/ ', coI, ":", Fehler
      IF (Fehler.ne.'&ff') RETURN
      IF (coI.eq.co) THEN ! yes, it is 'x'
         CALL OlfGetXY (j, K, n, X, Y, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (i.lt.1 .or. i.gt.n) THEN
            Fehler = 'search for x/z/par: invalid index of x'
            RETURN
            ENDIF
         rzxOlfG = X(i)
         RETURN
         ENDIF
      nZ = iOlfG (j, '#Z', Fehler)
      DO iZ = 1, nZ
         CALL OlfCnuG (j, 'z'//cl2(iZ), coI, un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (coI.eq.co) THEN ! yes, it is 'z#'
            CALL OlfGet1Z (j, K, iZ, rval, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            rzxOlfG = rval
            RETURN
            ENDIF
         ENDDO
      rzxOlfG = rOlfG (j, co, un, Fehler) ! either it's a param, or Fehler..
      IF (Fehler.ne.'&ff' .and. Fehler(1:4).ne.'&pnf')
     *     Print *, 'rzxOlfG -> rOlfG'

      END ! rzxOlfG

      REAL*8 FUNCTION rOlfGG (j, co, un, Fehler)
C     ------------------------------------------
            ! JWu 9nov96
         ! befor returning the value of parameter (co),
         ! check that (un) agrees with the stored unit.

      IMPLICIT NONE
      CHARACTER*(*) co, un, Fehler
      CHARACTER*40  unI
      INTEGER       j
      REAL*8        rval, rOlfG

      rOlfGG = 0
      rval = rOlfG (j, co, unI, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF (un.ne.unI) THEN
         CALL Compose3 (Fehler, 'coordinate '//co, ' has unit '//unI,
     *                  ' instead of requested '//un)
         RETURN
         ENDIF
      rOlfGG = rval

      END ! rOlfGG

      REAL*8 FUNCTION rOlfGGdef (j, co, un, rdef, Fehler)
C     ---------------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) co, un, Fehler
      CHARACTER*40  unI
      INTEGER       j
      REAL*8        rval, rdef, rOlfG

      rOlfGGdef = 0
      rval = rOlfG (j, co, unI, Fehler)
      IF (Fehler(1:4).eq.'&pnf') THEN
         rOlfGGdef = rdef
         Fehler = '&ff'
         RETURN
         ENDIF
      IF (Fehler.ne.'&ff') RETURN
      IF (un.ne.unI) THEN
         CALL Compose3 (Fehler, 'coordinate '//co, ' has unit '//unI,
     *                  ' instead of requested '//un)
         RETURN
         ENDIF
      rOlfGGdef = rval

      END ! rOlfGGdef

      REAL*8 FUNCTION rzOlfGG (j, K, co, un, Fehler)
C     ----------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) co, un, Fehler
      CHARACTER*40  unI
      INTEGER       j, K
      REAL*8        rval, rzOlfG

      rzOlfGG = 0
      rval = rzOlfG (j, K, co, unI, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF (un.ne.unI) THEN
         CALL Compose3 (Fehler, 'coordinate '//co, ' has unit '//unI,
     *                  ' instead of requested '//un)
         RETURN
         ENDIF
      rzOlfGG = rval

      END ! rzOlfGG

      REAL*8 FUNCTION rzxOlfGG (j, K, i, co, un, Fehler)
C     --------------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) co, un, Fehler
      CHARACTER*40  unI
      INTEGER       j, K, i
      REAL*8        rval, rzxOlfG

      rzxOlfGG = 0
      rval = rzxOlfG (j, K, i, co, unI, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF (un.ne.unI) THEN
         CALL Compose3 (Fehler, 'coordinate '//co, ' has unit '//unI,
     *                  ' instead of requested '//un)
         RETURN
         ENDIF
      rzxOlfGG = rval

      END ! rzxOlfGG

      SUBROUTINE iOlfCopy (jin, jout, what, Fehler)
C     ---------------------------------------------
            ! JWu 9nov96

      CHARACTER*(*) what, Fehler

      ival = iOlfG (jin, what, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL iOlfP (jout, what, ival, Fehler)

      END ! iOlfCopy

      SUBROUTINE rOlfCopy (jin, jout, Co, Fehler)
C     -------------------------------------------
            ! JWu 9nov96

      IMPLICIT NONE
      CHARACTER*(*) Co, Fehler
      CHARACTER*40  Un
      REAL*8        rval, rOlfG
      INTEGER       jin, jout

      rval = rOlfG (jin, Co, Un, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL rOlfP (jout, Co, Un, rval, Fehler)

      END ! rOlfCopy

      SUBROUTINE OlfCnuCheck2 (j1, j2, which, lCo, lUn, Fehler)
C     ---------------------------------------------------------
            ! JWu 10nov96
         ! Check agreement between coord's of two files
         ! Levels: 0=check nothing, 1=ask user, 2=always Fehler

      IMPLICIT NONE
      CHARACTER*(*) which, Fehler
      INTEGER       j1, j2, lCo, lUn
      CHARACTER*40  Co1, Un1, Co2, Un2
      CHARACTER*80  aus
      LOGICAL       qAsk

      IF (Fehler.ne.'&ff') RETURN

      CALL OlfCnuG (j1, which, Co1, Un1, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCnuG (j2, which, Co2, Un2, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      IF (Co1.ne.Co2) THEN
         IF     (lCo.ge.2) THEN
            CALL Compose2 (Fehler,
     *         'files have different '//which//' coordinates '//Co1,
     *         ' and '//Co2)
            RETURN
         ELSEIF (lCo.ge.1) THEN
            CALL Compose3 (aus,
     *         'files have different '//which//' coordinates '//Co1,
     *         ' and '//Co2, ' - continue ?')
            IF (.not.qAsk(aus)) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF
         ENDIF

      IF (Un1.ne.Un2) THEN
         IF     (lUn.ge.2) THEN
            CALL Compose2 (Fehler,
     *         'files have different '//which//' units '//Un1,
     *         ' and '//Un2)
            RETURN
         ELSEIF (lUn.ge.1) THEN
            CALL Compose3 (aus,
     *         'files have different '//which//' units '//Un1,
     *         ' and '//Un2, ' - continue ?')
            IF (.not.qAsk(aus)) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF
         ENDIF

      END ! OlfCnuCheck2

C  ====================================================================
C     5. Access to data blocks
C  ====================================================================

C  --------------------------------------------------------------------
C     5.1 data put
C  --------------------------------------------------------------------

      SUBROUTINE OlfPutZ (j, K, nZ, Z, Fehler)
C     ----------------------------------------
         ! JWu 23.1.91 (MemGetXY), 17may95
      ! Overwrite a spectrum in / MEM / by another one.

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, nZ, nZ1, MemBlockNum, MemBlockSiz
      REAL*8        Z(*)
      CHARACTER     Fehler*(*), cr3*3

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in OlfPutZ'
         RETURN
      ELSEIF (MemBlockNum(j).lt.MBH) THEN
         Print *, ' j K MBH MemBlNum(j) ', j, K, MBH, MemBlockNum(j)
         Fehler = 'OPZ/ cannot save spectrum before'//
     *            'header blocks are saved'
         RETURN
      ELSEIF (K.lt.0) THEN ! ??
         Fehler = 'OPZ/ call with K<0'
         RETURN
      ELSEIF (K.gt.MK) THEN
         Fehler = 'OPZ/ file shall not contain more than '//
     *            cr3(MK)//' spectra'
         RETURN
         ENDIF

      IF (K.eq.1) THEN
         IF (nZ.lt.0) THEN
            Fehler = 'OPZ/ call with nZ<0'
            RETURN
         ELSEIF (nZ.gt.MZ) THEN
            Fehler = 'OPZ/ call with nZ>MZ'
            RETURN
            ENDIF
      ELSE
         nZ1 = MemBlockSiz(j, MBH+1)
         IF (nZ1.lt.0) THEN
            Fehler = 'OPZ/ nZ(#1) invalid'
            RETURN
         ELSEIF (nZ.ne.nZ1) THEN
            Print *, ' nZ(in) nZ1(<MemBlock) j K ', nZ, nZ1, j, K
            Fehler = 'OPZ/ #Z different from spectrum 1'
            RETURN
            ENDIF
         ENDIF

      CALL MemBlockOvr (j, MBH+4*(K-1)+1, nZ, Z, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'OlfPutZ'
         RETURN
         ENDIF

      END ! OlfPutZ

      SUBROUTINE OlfPutXYD (j, K, n, X, Y, D, Fehler)
C     -----------------------------------------------
         ! JWu 23.1.91 (MemGetXY), 17may95
      ! Overwrite a spectrum in / MEM / by another one.

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       i, j, K, n, MemBlockNum
      REAL*8        X(*), Y(*), D(*), z
                    ! for semi-external input routines, X,Y,D may be declared
                    ! with a size MCin < MC : therefore we use (*).
      CHARACTER     Fehler*(*), cr3*3

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in OlfPutXYD'
         RETURN
      ELSEIF (MemBlockNum(j).lt.MBH) THEN
         Fehler = 'OSP/ cannot save spectrum before header'//
     *            ' blocks are saved'
         RETURN
      ELSEIF (n.lt.0) THEN
         Fehler = 'OSP/ n<0'
         RETURN
      ELSEIF (n.gt.MC) THEN
         Fehler = 'OSP/ call with n>MC'
         RETURN
      ELSEIF (n*3+1.gt.MmemSpe) THEN
         Fehler = 'PROGR ERR/ MmemSpe inconsistent with MC'
         RETURN
      ELSEIF (K.lt.0) THEN ! ??
         Fehler = 'OSP/ call with K<0'
         RETURN
      ELSEIF (K.gt.MK) THEN
         Fehler = 'OSP/ file shall not contain more than '//
     *            cr3(MK)//' spectra'
         RETURN
      ELSEIF (n.le.0) THEN
         Fehler = 'OSP/ call with n=0 replaced by OlfDelSpe'
         RETURN
      ELSEIF (MemBlockNum(j).lt.MBH+4*(K-1)+1) THEN
         Print *, ' j K MBH MemBlNum(j) ', j, K, MBH, MemBlockNum(j)
         Fehler = 'OlfPutXYD/ cannot save X,Y,D before Z'
         RETURN
         ENDIF

      CALL MemBlockOvr (j, MBH+4*(K-1)+2, n, X, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfPutXYD / X'
         RETURN
         ENDIF
      CALL MemBlockOvr (j, MBH+4*(K-1)+3, n, Y, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfPutXYD / Y'
         RETURN
         ENDIF
      CALL MemBlockOvr (j, MBH+4*(K-1)+4, n, D, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfPutXYD / D'
         RETURN
         ENDIF

      END ! OlfPutXYD

C  --------------------------------------------------------------------
C     .. indirect calls
C  --------------------------------------------------------------------

      SUBROUTINE OlfPutXY0 (j, K, n, X, Y, Fehler)
C     ---------------------------------------------
         ! Overwrite z, X and Y within a spectrum, set D=0
         ! (needed for consistency with GetXY)
         ! Vorsicht: soll nicht dazu dienen, die Fehlerrechnung auszuhebeln ..

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, n
      REAL*8        X(*), Y(*), D(MC)
      CHARACTER     Fehler*(*)
      DATA          D / MC * 0.d0 /

      CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)

      END ! OlfPutXY0

      SUBROUTINE OlfPutSpe (j, K, nZ, Z, n, X, Y, D, Fehler)
C     -------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, nZ, n
      REAL*8        X(*), Y(*), D(*), Z(*)
      CHARACTER     Fehler*(*)

      IF (Fehler.ne.'&ff') THEN
         Print *, 'Error on entry in OlfPutSpe'
         RETURN
         ENDIF
      CALL OlfPutZ (j, K, nZ, Z, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'j K nZ n ', j, K, nZ, n
         Print *, 'error passed through OPS(Z)'
         RETURN
         ENDIF
      CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OPS(XYD)'
         RETURN
         ENDIF

      END ! OlfPutSpe

C  --------------------------------------------------------------------
C     5.2  data get
C  --------------------------------------------------------------------

      SUBROUTINE OlfGetZ (j, K, nZ, Z, Fehler)
C     ----------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, nZ
      REAL*8        Z(MZ)
      CHARACTER     Fehler*(*)

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in OlfGetZ'
         RETURN
         ENDIF

      CALL MemBlockGet (j, MBH+4*(K-1)+1, nZ, MZ, Z, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfGetZ'
         RETURN
         ENDIF

      END ! OlfGetZ

      SUBROUTINE OlfGetN (j, K, nC, Fehler)
C     -------------------------------------
            ! JWu 6feb97

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INTEGER       j, K, nC, nCY, nCD, MemBlockSiz
      CHARACTER     Fehler*(*)

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in OlfGetN'
         RETURN
         ENDIF

      nC  =  MemBlockSiz (j, MBH+4*(K-1)+2)
      nCY =  MemBlockSiz (j, MBH+4*(K-1)+3)
      nCD =  MemBlockSiz (j, MBH+4*(K-1)+4)

      IF (nC.ne.nCY .or. nC.ne.nCD) THEN
         Print *, ' j K nX nY nD ', nC, nCY, nCD
         Fehler = 'OlfGetN/ nX inconsistent with nY or nD'
         RETURN
         ENDIF

      END ! OlfGetN

      SUBROUTINE OlfGetXYD (j, K, nC, X, Y, D, Fehler)
C     ------------------------------------------------
            ! JWu 16sep91, 17/19may95
         ! Get spectrum K of file j from memory.

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, nC, nCY, nCD
      REAL*8        X(MC), Y(MC), D(MC)
      CHARACTER     Fehler*(*)

      CALL MemBlockGet (j, MBH+4*(K-1)+2, nC,  MC, X, Fehler)
        IF (Fehler.ne.'&ff') GOTO 99
      CALL MemBlockGet (j, MBH+4*(K-1)+3, nCY, MC, Y, Fehler)
        IF (Fehler.ne.'&ff') GOTO 99
      CALL MemBlockGet (j, MBH+4*(K-1)+4, nCD, MC, D, Fehler)
        IF (Fehler.ne.'&ff') GOTO 99

      IF (nC.ne.nCY .or. nC.ne.nCD) THEN
         Print *, 'j K #X #Y #D', j, K, nC, nCY, nCD
         Fehler = 'OlfGetXYD/ nC inconsistent'
         RETURN
         ENDIF

      RETURN
 99   CONTINUE
      Print *, 'error passed through OlfGetXYD'

      END ! OlfGetXYD

      SUBROUTINE OlfGetXY (j, K, nC, X, Y, Fehler)
C     --------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, nC, nCY
      REAL*8        X(MC), Y(MC)
      CHARACTER     Fehler*(*)

      CALL MemBlockGet (j, MBH+4*(K-1)+2, nC,  MC, X, Fehler)
        IF (Fehler.ne.'&ff') GOTO 99
      CALL MemBlockGet (j, MBH+4*(K-1)+3, nCY, MC, Y, Fehler)
        IF (Fehler.ne.'&ff') GOTO 99

      IF (nC.ne.nCY) THEN
         Print *, 'j K #X #Y ', j, K, nC, nCY
         Fehler = 'OlfGetXY/ nC inconsistent'
         RETURN
         ENDIF

      RETURN
 99   CONTINUE
      Print *, 'error passed through OlfGetXY'

      END ! OlfGetXY

      SUBROUTINE OlfGetX (j, K, n, X, Fehler)
C     ---------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, n
      REAL*8        X(MC)
      CHARACTER     Fehler*(*)

      CALL MemBlockGet (j, MBH+4*(K-1)+2, n, MC, X, Fehler)
      IF (Fehler.ne.'&ff') Print *, 'error passed through OlfGetX'

      END ! OlfGetX

      SUBROUTINE OlfGetY (j, K, n, Y, Fehler)
C     ---------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       j, K, n
      REAL*8        Y(MC)
      CHARACTER     Fehler*(*)

      CALL MemBlockGet (j, MBH+4*(K-1)+3, n, MC, Y, Fehler)
      IF (Fehler.ne.'&ff') Print *, 'error passed through OlfGetY'

      END ! OlfGetY

C  --------------------------------------------------------------------
C     .. indirect access
C  --------------------------------------------------------------------

      SUBROUTINE OlfGetSpe (j, K, nZ, Z, nC, X, Y, D, Fehler)
C     -------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       i, j, K, nC, nZ, MemBlockSiz
      REAL*8        X(MC), Y(MC), D(MC), Z(MZ)
      CHARACTER     Fehler*(*)

      CALL OlfGetZ (j, K, nZ, Z, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfGetXYD (j, K, nC, X, Y, D, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! OlfGetSpe

      SUBROUTINE OlfGet1Z (j, K, iZ, zout, Fehler)
C     --------------------------------------------
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INTEGER       j, K, nZ, iZ
      REAL*8        Z(MZ), zout
      CHARACTER     Fehler*(*), cl4*4

      IF (Fehler.ne.'&ff') THEN
         Print *, 'error on entry in OlfGet1Z'
         RETURN
         ENDIF

      CALL OlfGetZ (j, K, nZ, Z, Fehler)
      IF (Fehler.ne.'&ff') THEN
         IF (Fehler.ne.'&ff') Print *, 'error passed through OlfGet1Z'
         RETURN
         ENDIF

      IF (iZ.lt.0) THEN
         Fehler = 'OlfGet1Z/ iz<1 requested'
         RETURN
      ELSEIF (iZ.gt.nZ) THEN
         CALL Compose2 (Fehler, 'OlfGet1Z/ iz='//cl4(iZ),
     *        ' requested while nz='//cl4(nZ))
         RETURN
         ENDIF
      zout = Z(iZ)

      END ! OlfGet1Z

      SUBROUTINE OlfGet1ZofK (j, iZ, nK, Z, Fehler)
C     ---------------------------------------------
            ! JWu 16sep91 (MemGetZ), 17may95
         ! Get nK, and for K=1,..,nK : Z(iZ;K) from memory.

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      CHARACTER     Fehler*(*)
      INTEGER       j, nB, iZ, K, nK, MemBlockNum
      REAL*8        Z(MK)

      nB = MemBlockNum(j)
      IF (nB.lt.1) THEN
         Fehler = 'OlfGet1ZofK/ File does not exist'
         RETURN
      ELSEIF (nB.le.MBH) THEN
         Fehler = 'OlfGet1ZofK/ File has uncomplete header'
         RETURN
         ENDIF
      nK = (nB - MBH) / 4
      DO K = 1, nK
         CALL OlfGet1Z (j, K, iZ, Z(K), Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      END ! OlfGet1ZofK

      SUBROUTINE OlfGetNofK (j, nK, Nz, Fehler)
C     -----------------------------------------
         ! Get nK, and for K=1,..,nK n(K) from memory.

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      CHARACTER     Fehler*(*)
      INTEGER       j, nB, K, nK, nC, Nz(MK), MemBlockNum

      nB = MemBlockNum(j)
      IF (nB.lt.1) THEN
         Fehler = 'OlfGetNofK/ File does not exist'
         RETURN
      ELSEIF (nB.le.MBH) THEN
         Fehler = 'OlfGetNof/ File has uncomplete header'
         RETURN
         ENDIF
      nK = (nB - MBH) / 4
      DO K = 1, nK
         CALL OlfGetN  (j, K, Nz(K), Fehler)
         IF (Fehler.ne.'&ff') THEN
            Print *, 'error passed through OlfGetNofK'
            RETURN
            ENDIF
         ENDDO

      END ! OlfGetNofK

      SUBROUTINE OlfDel1Z (j, jout, iZ, Fehler)
C     -----------------------------------------
            ! JWu 4/6mrz98

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      CHARACTER         Fehler*(*), Co*40, Un*40, cl2*2
      INTEGER           j, jout, iZ, nB, MemBlockNum, nK, K, nZ, iiZ
      REAL*8            Z(MZ)

      nB = MemBlockNum(j)
      nK = (nB - MBH) / 4
      DO K = 1, nK
         CALL OlfGetZ (j, K, nZ, Z, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (iZ.gt.nZ) THEN
            Fehler = ' cannot delete coordinate z'//cl2(iZ)
            RETURN
            ENDIF
         DO iiZ = iZ, nZ-1
            Z(iiZ) = Z(iiZ+1)
            ENDDO
         CALL OlfPutZ (jout, K, nZ-1, Z, Fehler)
         ENDDO
         CALL OlfCnuG (jout, 'y', Co, Un, Fehler)

      DO iiZ = iZ, nZ-1
         CALL OlfCnuG (j, 'z'//cl2(iiZ+1), Co, Un, Fehler)
         CALL OlfCnuP (jout, 'z'//cl2(iiZ), Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      CALL OlfCnuP (jout, 'z'//cl2(nZ), ' ', ' ', Fehler)

      END ! OlfDel1Z

C  --------------------------------------------------------------------
C     5.3 spectra copy
C  --------------------------------------------------------------------

      SUBROUTINE OlfCopZ (jin, jout, Kin, Kout, Fehler)
C     -------------------------------------------------
            ! JWu 6feb97 corr. 4mrz98
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INTEGER        jin, jout, Kin, Kout, nZ
      REAL*8         Z(MZ)
      CHARACTER*(*)  Fehler

      CALL OlfGetZ (jin, Kin, nZ, Z, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfPutZ (jout, Kout, nZ, Z, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! OlfCopZ

      SUBROUTINE OlfCopXYD (jin, jout, Kin, Kout, Fehler)
C     ---------------------------------------------------
            ! JWu 6feb97
      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INTEGER        jin, jout, Kin, Kout, nC
      REAL*8         X(MC), Y(MC), D(MC)
      CHARACTER*(*)  Fehler

      CALL OlfGetXYD (jin, Kin, nC, X, Y, D, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfPutXYD (jout, Kout, nC, X, Y, D, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! OlfCopXYD

      SUBROUTINE OlfCopSpe (jin, jout, Kin, Kout, Fehler)
C     ---------------------------------------------------
            ! JWu 6feb97
      IMPLICIT NONE
      INTEGER        jin, jout, Kin, Kout
      CHARACTER*(*)  Fehler

      CALL OlfCopZ   (jin, jout, Kin, Kout, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCopXYD (jin, jout, Kin, Kout, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! OlfCopSpe

C  --------------------------------------------------------------------
C     5.4 spectra delete
C  --------------------------------------------------------------------

      SUBROUTINE OlfDelSpe (j, K, Fehler)
C     -----------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'

      INTEGER       i, j, K, nC, MemBlockNum
      REAL*8        R(1) ! dummy
      CHARACTER     Fehler*(*)

      CALL MemBlockOvr (j, MBH+4*(K-1)+4, -1, R, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfDelSPe/ ~D'
         RETURN
         ENDIF
      CALL MemBlockOvr (j, MBH+4*(K-1)+3, -1, R, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfDelSPe/ ~Y'
         RETURN
         ENDIF
      CALL MemBlockOvr (j, MBH+4*(K-1)+2, -1, R, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfDelSPe/ ~X'
         RETURN
         ENDIF
      CALL MemBlockOvr (j, MBH+4*(K-1)+1, -1, R, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'error passed through OlfDelSPe/ ~z'
         RETURN
         ENDIF

      END ! OlfDelSpe

      SUBROUTINE FileClean ()
C     -----------------------
            ! JWu 26nov91
         ! Delete files which are not closed or which are empty.

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      REAL*8 dummy(1)

      CHARACTER*80 Fehler

      Fehler = '&ff'
      DO j = MemBlockInq('nF'), 1, -1
         nK = (MemBlockNum(j) - MBH) / 4
         DO K = MemBlockNum(j), MBH+1, -4
            IF (MemBlockSiz(j,K-2).le.0) THEN
               CALL Say2 (' .. cleaning up/ delete empty spectrum '//
     *              cl4((K-MBH)/4),
     *              ' in file '//cl3(j))
               DO KK = K, K-3, -1
                  CALL MemBlockPut (j, KK, -1, dummy, Fehler)
                  IF (Fehler.ne.'&ff') CALL FehlerGong (Fehler,1)
                  ENDDO
               ENDIF
            ENDDO
         nK = (MemBlockNum(j) - MBH) / 4
         IF (nK.le.0) THEN
            Print '(a)', ' .. cleaning up/ delete empty file '//cl3(j)
            CALL MemBlockPut (j, 0, -1, dummy, Fehler)
            IF (Fehler.ne.'&ff') CALL FehlerGong (Fehler,1)
            ENDIF
         ENDDO

      END ! FileClean
