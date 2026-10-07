C  ====================================================================
C
C      Program  IDA   :  Inelastic Data Analysis
C      Modul    i10   :     input / output
C
C  ====================================================================

C     Contents :
C
C        1.  Load / save files:
C               FileLoad, FileReadOld, FileSave
C
C        2.  File access:
C               AskPath, OpenIdaFile
C
C        3.  Read / write in '96 format:
C               FileWrite96, FileRead96
C
C        4.  Read old data files:
C               ['92 formats] LoadSpectrum92, i/r/tOlfPold
C               [IED format]  ReadIED(sub)

C  ====================================================================
C  i10 / 1 :   Read / write file (interface to on-line-memory)
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE FileLoad (Object, Fehler)
C     ------------------------------------
         ! Read data files and store the contents

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      CHARACTER*(*) Object, Fehler

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       iPar(MP)
      REAL*8        rPar(MP)
      CHARACTER*40  tPar(MP), FilExt, DirExt, FilInt, DirInt
      REAL*8        X(MC), Y(MC), D(MC)
      CHARACTER*79  aus, Path
      DATA          DirExt /' '/

      IF (Fehler.ne.'&ff') THEN
         Print *, 'FileLoad/ error on entry'
         RETURN
         ENDIF

C  Loop files :
      qLoop = .false.
      iProtect = 0
 1    CONTINUE
         IF     (Object.eq.' ') THEN ! interactive loop
            aus = ' Load file'
            FilExt= ' '
            qLoop = .true.
         ELSEIF (Object(1:5).eq.'&int ') THEN ! internal (e.g. for numtab)
            aus = '&noq'
            DirExt = ' '
            FilExt = Object(6:lenU(Object))
            iProtect = 4
         ELSE
            aus = '&noq'
            FilExt = Object
            ENDIF
 2       CALL AskPath (aus, FilExt, DirExt, Path, Fehler)
         IF (Path.eq.' ') RETURN

         CALL OpenIdaFile (31, Path, iFormat, Fehler)
         IF (Fehler.eq.'&fnf') THEN
            CALL Gong(1)
            Print *, 'Cannot find data file '//Path
            Fehler = '&ff'
            aus = ' Load file'
            GOTO 2
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCreate (jout, Kout, '&noask', '&noask', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (iFormat.eq.960) THEN
            CALL FileRead96 (31, jout, nK, Fehler)
         ELSE
            CALL FileReadOld (31, iFormat, Path, jout, nK, Fehler)
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

         ! update internal file name and directory :
         CALL tOlfG (jout, 'fil', FilInt, Fehler)
         CALL tOlfG (jout, 'dir', DirInt, Fehler)
         IF (Fehler.ne.'&ff') THEN
            Fehler = '&ff' ! no file or dir name given -> na und ?
         ELSE
            ! c96/7 add a comment
            ENDIF
         CALL tOlfP (jout, 'fil', FilExt, Fehler)
         CALL tOlfP (jout, 'dir', DirExt, Fehler)

         CALL OlfClos (jout, nK, Fehler)
         CALL MemFileStatP (jout, iProtect, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  End of loop - next file ?
      IF (qLoop) GOTO 1

      END ! FileLoad

      SUBROUTINE FileReadOld (NU, iFormat, Path, jout, Kout, Fehler)
C     --------------------------------------------------------------
         ! Read data files and store the contents

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Path, Fehler
      INTEGER       iPar(MP)
      REAL*8        rPar(MP)
      CHARACTER*40  tPar(MP), null
      REAL*8        X(MC), Y(MC), D(MC), ZZ(1)
      CHARACTER*79  aus

C  Loop spectra :
      nK = 0
      DO K = 1, MK+1
         CALL LoadSpectrum92 (NU, Path, iFormat, K, nK, iPar, rPar,
     *                        tPar, z, MC, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') THEN
            RETURN
         ELSEIF (nK.gt.0) THEN
            RETURN ! eof - regular exit
         ELSEIF (K.gt.MK) THEN
            Fehler = 'File contains more than '//cr3(MK)//' spectra'
            RETURN
            ENDIF

         IF (K.eq.1) THEN
            CALL iOlfPold (jout, iPar, Fehler)
            CALL rOlfPold (jout, rPar, Fehler)
            CALL tOlfPold (jout, tPar, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            DO ip = 11, iPar(11)
               IF (tPar(ip).ne.' ' .and. tPar(ip).ne.null)
     *             CALL OlfComAddFull (jout, ' ', tPar(ip), 0, Fehler)
               ENDDO
            ENDIF

         ZZ(1) = z
         CALL OlfPutSpe (jout, K, 1, ZZ, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      Fehler = 'FileReadOld/ unexpectedly K>MK'

      END ! FileReadOld

      SUBROUTINE FileSave (nJList, JList, Fehler)
C     -------------------------------------------
         ! Completely new 11-nov96

      IMPLICIT NONE
      INCLUDE      'l_def.f'
      CHARACTER*(*) Fehler*(*)
      CHARACTER*40  FilInt, DirInt, FilExt, DirExt, FilJ, DirJ
      CHARACTER*80  Path
      CHARACTER     stat*2
      INTEGER       nJList, JList(*), lj, j, jj, MemBlockInq
      LOGICAL       qOvFD, qOvKWTD
      DATA          stat /'n!'/

      IF (nJList.le.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      qOvKWTD = .false. ! overwrite: know-what-to-do
      DO lj = 1, nJList
         j = JList(lj)

         CALL tOlfG (j, 'fil', FilInt, Fehler)
         CALL tOlfG (j, 'dir', DirInt, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! default :
         FilExt = FilInt
         DirExt = DirInt

         CALL AskOpenDatFile (' Write to file ', 32, FilExt, DirExt,
     *                        'i96',stat, 'seq', 'for', 0)
         IF (FilExt.eq.' ') THEN ! legal emergency exit
            Fehler = ' '
            RETURN
            ENDIF

         IF (FilExt.ne.FilInt .or. DirExt.ne.DirInt) THEN
            IF (.not.qOvKWTD) qOvFD = .true.
c ausser Betrieb : = qAskD (' Overwrite internal file/dir', intq(qOvFD))
            qOvKWTD = .true.
            IF (qOvFD) THEN
               CALL tOlfP (j, 'fil', FilExt, Fehler)
               CALL tOlfP (j, 'dir', DirExt, Fehler)
               ENDIF
            ENDIF

         ! here a switch between data formats ?
         CALL FileWrite96 (32, j, Fehler)
         Close (32) ! close the file also in case of write-out error (feb99)
         IF (Fehler.ne.'&ff') RETURN

         ! Update status bit (4mar98) :
         CALL MemFileBitP (j, 1, 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         DO jj = 1, MemBlockInq ('nF')
            IF (jj.ne.j) THEN
               CALL tOlfG (jj, 'fil', FilJ, Fehler)
               CALL tOlfG (jj, 'dir', DirJ, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (FilJ.eq.FilExt .and. DirJ.eq.DirExt)
     *            CALL MemFileBitP (jj, 1, 1, Fehler)
               ENDIF
            ENDDO

         ENDDO ! lj

      END ! FileWrite

C  ====================================================================
C  i10 / 2 :   File access
C  ====================================================================
            ! old LoadSpectrum splitted in subroutines 10nov96

      SUBROUTINE AskPath (Quest, File, Dir, Path, Fehler)
C     ---------------------------------------------------
         ! ask for File, optionally change Dir, construct Path=Dir/File

      CHARACTER*(*) Quest, File, Dir, Path, Fehler
      CHARACTER*80  aus

      IF (Quest.ne.'&noq') THEN
 21      CONTINUE
         IF (Dir.ne.' ' .and. Dir.ne.'&nod') THEN
            CALL Compose2 (aus, Quest,
     *         ' (dir=' // Dir(1:lenU(Dir)) // ') ?')
         ELSE
            CALL Compose2 (aus, Quest, ' ?')
            ENDIF
         CALL FrageC (aus, File)
         ! directory command given ?
         IF (File(1:3).eq.'cd ' .and. Dir.ne.'&nod') THEN
            Dir = File(4:len(File))
            CALL DelLeft (Dir)
            GOTO 21
            ENDIF
         ENDIF
      IF (File.eq.' ') THEN
         Path = ' ' ! no file wanted
         RETURN
         ENDIF
      ! construct full path name (2sep93, 3dec93 here) :
      IF (Dir.ne.' ' .and. Dir.ne.'&nod') THEN
         CALL Compose2 (Path, Dir, '/'//File)
      ELSE
         Path = File
         ENDIF

      END ! AskPath

      SUBROUTINE OpenIdaFile (NU, Path, iFormat, Fehler)
C     --------------------------------------------------
         ! try to open ida file and determine its format
            ! from 1996 on, a new format means always a new extension,
            ! which simplifies this routine as well as file maintenance.

      CHARACTER*(*) Path, Fehler
      CHARACTER     Code*8

      IF (Fehler.ne.'&ff') RETURN
      Code = '&eof'

      ! first attempt: try extension .i96 (must be ASCII-96)
      CALL OpenFile (NU, Path, 'i96', 'a', Fehler)
      IF (Fehler.eq.'&ff') THEN
         Read (NU,  '(a8)', err=20) Code
         IF (Code.eq.'ASCII-96') THEN ! successfully opened
            iFormat = 960
            RETURN ! success
            ENDIF
 20      Close (NU)
         Fehler = ' .i96-file has illegal format code '//Code
         RETURN
         ENDIF
      Fehler = '&ff'

      ! second attempt: try extension .dat, code binary92
      CALL OpenDatFile (NU, Path, 'dat', 'l', 'seq', '-cc', -1, Fehler)
      IF (Fehler.eq.'&ff') THEN
         Read (NU, err=30) Code
         IF (Code.eq.'binary92') THEN ! successfully opened
            iFormat = 921
            RETURN ! success
            ENDIF
 30      Close (NU) ! .dat, but not binary92
         ! in that case, .dat is a very old ascii file
         CALL OpenFile (NU, Path, 'dat', 'a', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         GOTO 80
         ENDIF
      Fehler = '&ff'

      ! last attempt : extension .asc
      CALL OpenFile (NU, Path, 'asc', 'a', Fehler)
      IF (Fehler.ne.'&ff') THEN
         Fehler = '&fnf' ! file-not-found
         RETURN
         ENDIF

 80   CONTINUE ! now an ascii-file is open:
      Read (NU, '(a8)', err=1012) Code
      IF (Code.eq.'ASCII-92') THEN
         iFormat = 922
         RETURN
      ELSE
         Print *, ' Code = ', Code
         ENDIF
 1012 CONTINUE  ! 21jul92 : error 1011 occured in some old files
      iFormat = 890 ! or 893 or ...
      Close(NU) ! there is no code : rewind, and reopen later (must be .dat)
      RETURN

      END ! OpenIdaFile

C  ====================================================================
C  i10 / 3 :   Write / read files in '96 format
C  ====================================================================

      SUBROUTINE FileWrite96 (NU, j, Fehler)
C     --------------------------------------
            ! Completely new 11-nov96
         ! Write internal file j to external unit 32

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'
      CHARACTER*(*) Fehler*(*)
      CHARACTER     co*24, un*24, lab*12, line*80,
     *              format*24, format2*24, format3*24
      INTEGER       NU, j, np, n, i, ival, iOlfG, nK, K, nZ
      REAL*8        rval, rOlfG, Z(MZ), z1, X(MC), Y(MC), D(MC)

      Write (NU, '(a8)', err=90) 'ASCII-96'

      ! block 2: t-par's
      format = '(a24,a56)'
      Write (NU, '(a)', err=91) format
      np = 0
 21   CONTINUE
         co   = '&pbn '//cv3(np+1)
         CALL tOlfG (j, co, line, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (co.eq.'&eop') GOTO 29
         np = np + 1
         IF (co.eq.'&empty') GOTO 21
         Write (NU, format, err=91) co, line
         GOTO 21
 29   CONTINUE
      Write (NU, '(a)',err=91)
     * '&eob 2 (end of block)  ------------------------'

      ! block 3: i-par's
      format = '(a24,i16)'
      Write (NU, '(a)', err=91) format
      np = 0
 31   CONTINUE
         co   = '&pbn '//cv3(np+1)
         ival = iOlfG (j, co, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (co.eq.'&eop') GOTO 39
         np = np + 1
         IF (co.eq.'&empty') GOTO 31
         Write (NU, format,err=91) co, ival
         GOTO 31
 39   CONTINUE
      Write (NU, '(a)')
     * '&eob 3 (end of block)  ------------------------'

      ! block 4: r-par's
      format = '(a24,a24,g20.10)'
      Write (NU, '(a)',err=91) format
      np = 0
 41   CONTINUE
         co   = '&pbn '//cv3(np+1)
         rval = rOlfG (j, co, un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (co.eq.'&eop') GOTO 49
         np = np + 1
         IF (co.eq.'&empty') GOTO 41
         Write (NU, format,err=91) co, un, rval
         GOTO 41
 49   CONTINUE
      Write (NU, '(a)')
     * '&eob 4 (end of block)  ------------------------'

      ! block 5: coord's
      format = '(a4,a24,a24)'
      Write (NU, '(a)',err=91) format
      np = 0
 51   CONTINUE
         lab   = '&pbn '//cv3(np+1)
         CALL OlfCnuG (j, lab, co, un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (lab.eq.'&eop') GOTO 59
         np = np + 1
         IF (lab.eq.'&empty') GOTO 51
         Write (NU, format,err=91) lab, co, un
         GOTO 51
 59   CONTINUE
      Write (NU, '(a)',err=91)
     * '&eob 5 (end of block)  ------------------------'

      ! block 6: long-doc
      format = '(a80)'
      Write (NU, '(a)',err=91) format
      np = 1
 61   CONTINUE
         CALL OlfComLinG (j, np, line, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (line.eq.'&eoc') GOTO 69
         Write (NU, format,err=91) line
         np = np + 1
         GOTO 61
 69   CONTINUE
      Write (NU, '(a)',err=91)
     * '&eob 6 (end of block)  ------------------------'

      ! data blocks
      nK = iOlfG (j, '#spectra', Fehler)
      nZ = iOlfG (j, '#Z', Fehler)
      format = '(2i8)'
      format2= '(i16,4g16.8/5g16.8)'
      !format3= '(3(2x,g16.8))'
      format3= '(3(2x,g24.16))' !Artem: increase precision
      Write (NU, '(a)',err=91) format
      Write (NU, '(a)',err=91) format2
      Write (NU, '(a)',err=91) format3
      Write (NU, format,err=91) nK, nZ
      DO K = 1, nK
         Write (NU, '(a,i3)',err=91) '&spectrum ', K
         CALL OlfGetSpe (j, K, nZ, Z, n, X, Y, D, Fehler)
         Write (NU, format2, err=91) n, (Z(i), i=1,nZ)
         Write (NU, format3, err=91) (X(i),Y(i),D(i), i=1,n)
         ENDDO

      RETURN

 90   Fehler = 'write error in 1st line '//
     *         '(disk full or something worse ??)'
      RETURN
 91   Fehler = 'write error (disk full ?)'
      RETURN

      END ! FileWrite96

      SUBROUTINE FileRead96 (NU, j, nK, Fehler)
C     -----------------------------------------
            ! Completely new 11-nov96
         ! Get internal file j from external unit 32

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      INCLUDE      'l_def.f'
      CHARACTER*(*) Fehler*(*)
      CHARACTER     co*24, un*24, lab*12, line*80,
     *              format*24, format2*24, format3*24
      INTEGER       NU, j, n, i, np, nK, K, nZ, nZlab, ival, iOlfG
      REAL*8        rval, rOlfG, Z(MZ), z1, X(MC), Y(MC), D(MC)

      ! block 2: t-par's
      Read (NU, '(a)', err=921) format
      np = 0
 21   CONTINUE
         Read (NU, format, err=922) co, line
         IF (co(1:4).eq.'&eob') GOTO 29
         CALL tOlfP (j, co, line, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         GOTO 21
 29   CONTINUE

      ! block 3: i-par's
      Read (NU, '(a)', err=931) format
 31   CONTINUE
         Read (NU, '(a)', err=932) line
         IF (line(1:4).eq.'&eob') GOTO 39
         read (line, format, err=933) co, ival
         CALL iOlfP (j, co, ival, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         GOTO 31
 39   CONTINUE

      ! block 4: r-par's
      Read (NU, '(a)', err=941) format
 41   CONTINUE
         Read (NU, '(a)', err=942) line
         IF (line(1:4).eq.'&eob') GOTO 49
         read (line, format, err=943) co, un, rval
         CALL rOlfP (j, co, un, rval, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         GOTO 41
 49   CONTINUE

      ! block 5: coord's
      Read (NU, '(a)', err=951) format
      nZlab = -2
 51   CONTINUE
         Read (NU, format, err=952) lab, co, un
         IF (lab(1:4).eq.'&eob') GOTO 59
         IF (lab.eq.'z') lab = 'z1'
         IF (co.eq.' ') THEN
            co = '? [empty coordinate name in ASCII-96 file]'
            ENDIF
         CALL OlfCnuP (j, lab, co, un, Fehler)
         nZlab = nZlab + 1
         IF (Fehler.ne.'&ff') RETURN
         GOTO 51
 59   CONTINUE

      ! block 6: long-doc
      Read (NU, '(a)', err=961) format
 61   CONTINUE
         Read (NU, format, err=962) line
         IF (line(1:4).eq.'&eob') GOTO 69
         CALL OlfComAddFull (j, ' ', line, 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         GOTO 61
 69   CONTINUE

      ! data blocks
      Read (NU, '(a)',err=971) format
      Read (NU, '(a)',err=971) format2
      Read (NU, '(a)',err=971) format3
      Read (NU, format,err=972) nK, nZ
      IF (nZ.gt.nZlab) THEN
         Print *,
     * 'incorrect data format: adding name(s) for z coordinate(s)'
         DO i = nZlab+1, nZ
            CALL OlfCnuP (j, 'z'//cl2(i),
     * '? [unlabelled coordinate in ASCII-96 file]', ' ', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO
      ELSEIF (nZ.lt.nZlab) THEN
         Print *, 'WARNING/ unused z-labels in input file/ nZ nZlab ',
     *            nZ, nZlab
         ENDIF
      DO K = 1, nK
         Read (NU, '(a)',err=973) line
         IF (line(1:9).ne.'&spectrum') THEN
            Fehler = 'missed beginning of spectrum'
            RETURN
            ENDIF
         Read (NU, format2, err=975) n, (Z(i), i=1,nZ)
         IF (n.gt.MC) THEN
            Print *, ' nC, MC: ', n, MC
            Fehler = ' too many channels/ recompile with bigger MC'
            RETURN
            ENDIF
         Read (NU, format3, err=976) (X(i),Y(i),D(i), i=1,n)
         CALL OlfPutSpe (j, K, nZ, Z, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') THEN
            Print *, Fehler
            Fehler = 'FileRead: cannot save spectrum '//cl6(K)
            RETURN
            ENDIF
         ENDDO

      Close(NU)

      RETURN

 921  Fehler = 'read err 921'
      RETURN
 922  Fehler = 'read err 922'
      RETURN
 931  Fehler = 'read err 931'
      RETURN
 932  Fehler = 'read err 932'
      RETURN
 933  Fehler = 'read err 933'
      RETURN
 941  Fehler = 'read err 941'
      RETURN
 942  Fehler = 'read err 942'
      RETURN
 943  Fehler = 'read err 943'
      RETURN
 951  Fehler = 'read err 951'
      RETURN
 952  Fehler = 'read err 952'
      RETURN
 961  Fehler = 'read err 961'
      RETURN
 962  Fehler = 'read err 962'
      RETURN
 971  Fehler = 'read err 971'
      RETURN
 972  Fehler = 'read err 972'
      RETURN
 973  Fehler = 'read err 973'
      RETURN
 975  Fehler = 'read err 975'
      RETURN
 976  Fehler = 'read err 976'
      RETURN

      END ! FileRead96

C  ====================================================================
C  i10 / 4 :   Read spectra from old data files
C  ====================================================================

C  --------------------------------------------------------------------
C     read 92' formats
C  --------------------------------------------------------------------

      SUBROUTINE LoadSpectrum92 (NU, Path, iForm, K, nK,
     *                           iPar, rPar, tPar,
     *                           z, M, n, X, Y, D, Fehler)
C     ----------------------------------------------------
            ! JWu 4may92
         ! Read one spectrum from a file in the old formats
         !    921 = binary92
         !    922 = ASCII-92
         !    890 = IED / SQW
         ! If K = 1 then open the file
         !   (if Quest<>'&noq' then ask for its name),
         ! if eof then close the file, and return nK>0.

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      CHARACTER*40  tPar(MP)
      DIMENSION     iPar(MP), rPar(MP), X(*), Y(*), D(*)
      CHARACTER*(*) Path, Fehler
      CHARACTER     aus*80, cl3*3, cl6*6

C  Open File :
      IF (K.eq.1) THEN

C  Read header ?
         IF     (iForm.eq.921) THEN
            Read (NU, err=102) niP, (iPar(i), i=1, niP)
            Read (NU, err=102) nrP, (rPar(i), i=1, nrP)
            Read (NU, err=102) ntP, (tPar(i), i=1, ntP)
         ELSEIF (iForm.eq.922) THEN
            Read (NU, '(3(i4,2x))', err=1120) niP, nrP, ntP
            DO i = 1, niP
               Read (NU, '(i10)',  err=1121) iPar(i)
               ENDDO
            DO i = 1, nrP
               Read (NU, '(g16.8)',err=1122) rPar(i)
               ENDDO
            DO i = 1, ntP
               Read (NU, '(a40)',  err=1123) tPar(i)
               ENDDO
            ENDIF

         ENDIF ! K=1

      IF (nK.ne.0)
     * CALL Absturz ('LoadSpectrum', 'someone manipulated nK')

C  Read one spectrum :
      IF (iForm.eq.921 .or. iForm.eq.922) THEN
         IF (iForm.eq.921) THEN
            Read (NU, err=103) n
         ELSE
            Read (NU, '(i8)', err=1131) n
            ENDIF
         IF (n.eq.-1) THEN
            Close (NU)
            GOTO 100 ! eof
         ELSEIF (n.le.0) THEN
            Fehler = 'Wrong entry : n<-2 or n=0 in spectrum '//cl3(K)
            Close (NU)
            RETURN
         ELSEIF (n.gt.M) THEN ! variable M since 30oct92
            Fehler = 'Spectrum too long : n = '//cl6(n)
            Close (NU)
            RETURN
            ENDIF

         IF (iForm.eq.921) THEN
            Read (NU, err=103) z
            Read (NU, err=104) (X(i), i=1, n)
            Read (NU, err=105) (Y(i), i=1, n)
            Read (NU, err=106) (D(i), i=1, n)
         ELSE
            Read (NU, '(g16.8)',      err=1132)  z
            Read (NU, '(3(2x,g16.8))',err=1133) (X(i),Y(i),D(i), i=1,n)
            ENDIF

      ELSEIF (iForm.eq.890) THEN
         CALL ReadIED (Path, K, iPar, rPar, tPar, z, M, n, X, Y, D,
     *                 Fehler)
         IF (Fehler.eq.'&eof') THEN
            Fehler = '&ff'
            GOTO 100
         ELSEIF (Fehler.ne.'&ff') THEN
            RETURN
            ENDIF

      ELSE
         Fehler = ' LoadSpectrum92 with invalid format '//cl3(iForm)
         RETURN
         ENDIF

      IF (tPar(7).eq.' ' .and. tPar(10).eq.' ') z = 0 ! to prevent z=..E-312
      RETURN
C  Regular exit.

C  EOF condition :
 100  CONTINUE
      nK = K-1
      IF (nK.le.0) CALL Absturz (
     *   'LoadSpectrum', 'eof in first spectrum')
      RETURN

C  Errors :
c 101  CONTINUE    ! ausser Betrieb -- zur Zeit keine Fehlermeldung
c      Fehler = ' File code could not be read'
c      Close (NU)
c      RETURN
 102  CONTINUE
      Fehler = ' Error in header of binary file'
      Close (NU)
      RETURN
 103  CONTINUE
      CALL Compose2 (Fehler,
     *    ' Error while reading z in spectrum '//cl3(K),
     *    ' of binary file')
      Close (NU)
      RETURN
 104  CONTINUE
      CALL Compose2 (Fehler,
     *    ' Error while reading X in spectrum '//cl3(K),
     *    ' of binary file')
      Close (NU)
      RETURN
 105  CONTINUE
      CALL Compose2 (Fehler,
     *    ' Error while reading Y in spectrum '//cl3(K),
     *    ' of binary file')
      Close (NU)
      RETURN
 106  CONTINUE
      CALL Compose2 (Fehler,
     *    ' Error while reading D in spectrum '//cl3(K),
     *    ' of binary file')
      RETURN
 1120 CONTINUE
      Fehler = ' Error in parameter block header of ASCII file'
      Close (NU)
      RETURN
 1121 CONTINUE
      Fehler = ' Error in integer parameter block of ASCII file'
      Close (NU)
      RETURN
 1122 CONTINUE
      Fehler = ' Error in real parameter block of ASCII file'
      Close (NU)
      RETURN
 1123 CONTINUE
      Fehler = ' Error in text parameter block of ASCII file'
      Close (NU)
      RETURN
 1131 CONTINUE
      CALL Compose2 (Fehler,
     * ' Error while reading n in spectrum '//cl3(K), ' of ASCII file')
      Close (NU)
      RETURN
 1132 CONTINUE
      CALL Compose2 (Fehler,
     * ' Error while reading z in spectrum '//cl3(K), ' of ASCII file')
      Close (NU)
      RETURN
 1133 CONTINUE
      CALL Compose2 (Fehler,
     * ' Error while reading X Y D in spectrum '//cl3(K),
     * ' of ASCII file')
      Close (NU)
      RETURN

      END ! LoadSpectrumOld

C  --------------------------------------------------------------------
C     translation <-> old data format
C  --------------------------------------------------------------------

      SUBROUTINE iOlfPold (j, iPar, Fehler)
C     -------------------------------------

      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      CHARACTER*(*) Fehler
      INTEGER       j, iPar(MP)

      IF (Fehler.ne.'&ff') RETURN
                         CALL iOlfP (j, '?cu', iPar(3), Fehler)
      IF (iPar( 4).ne.0) CALL iOlfP (j, 'fu#', iPar(4), Fehler)
      IF (iPar( 5).ne.0) CALL iOlfP (j, '#fit-par', iPar(5), Fehler)
      IF (iPar( 6).ne.0) CALL iOlfP (j, '?weight-stp-x', iPar(6),
     *                               Fehler)
      IF (iPar( 7).ne.0) CALL iOlfP (j, '?weight-err-y', iPar(7),
     *                               Fehler)
      IF (iPar( 8).ne.0) CALL iOlfP (j, 'fit-dat-file#', iPar(8),
     *                               Fehler)
      IF (iPar( 9).ne.0) CALL iOlfP (j, '@fixed', iPar(9), Fehler)
      IF (iPar(12).ne.0) CALL iOlfP (j, '?det-bal-sym', iPar(12),
     *                               Fehler)
      IF (iPar(13).ne.0) CALL iOlfP (j, '@sam-erg-gain', iPar(13),
     *                               Fehler)
      IF (iPar(14).ne.0) CALL iOlfP (j, 'plot-#pts', iPar(14), Fehler)
      IF (iPar(15).ne.0) CALL iOlfP (j, 'plot-sy#', iPar(15), Fehler)
      IF (iPar(16).ne.0) CALL iOlfP (j, '?weight-log-y', iPar(16),
     *                               Fehler)
      IF (iPar(17).ne.0) CALL iOlfP (j, '?conv', iPar(17), Fehler)
      IF (iPar(18).ne.0) CALL iOlfP (j, 'fit-par-file#', iPar(18),
     *                               Fehler)

      END ! iOlfPold

      SUBROUTINE rOlfPold (j, rPar, Fehler)
C     -------------------------------------
      IMPLICIT NONE
      INCLUDE      'i_dim.f'
      CHARACTER*(*) Fehler
      INTEGER       j
      REAL*8        rPar(MP)

      IF (Fehler.ne.'&ff') RETURN
      IF (rPar( 1).ne.0) CALL rOlfP (j, '2th', ' ', rPar(1), Fehler)
      IF (rPar( 2).ne.0) CALL rOlfP (j, 'E0', 'meV', rPar(2), Fehler)
      IF (rPar( 3).ne.0) CALL rOlfP (j, 'q', 'A-1', rPar(3), Fehler)
      IF (rPar( 4).ne.0) CALL rOlfP (j, 'L[fpi1]', 'mm', rPar(4),
     *                               Fehler)
      IF (rPar( 6).ne.0) CALL rOlfP (j, 'P', 'kbar', rPar(6), Fehler)
      IF (rPar( 7).ne.0) CALL rOlfP (j, 'T', 'K', rPar(7), Fehler)
      IF (rPar( 8).ne.0) CALL rOlfP (j, 'at-mass', 'amu', rPar(8),
     *                               Fehler)
      IF (rPar( 9).ne.0) CALL rOlfP (j, 'cts[mon]', ' ', rPar(9),
     *                               Fehler)
      IF (rPar(10).ne.0) CALL rOlfP (j, 't[scan]', 'sec/10 ?',
     *                               rPar(10), Fehler)
      IF (rPar(14).ne.0) CALL rOlfP (j, 'plot-i', ' ', rPar(14),Fehler)
      IF (rPar(15).ne.0) CALL rOlfP (j, 'plot-f', ' ', rPar(15),Fehler)
      IF (rPar(16).ne.0) CALL rOlfP (j, 'fit-i', ' ', rPar(16), Fehler)
      IF (rPar(17).ne.0) CALL rOlfP (j, 'fit-f', ' ', rPar(17), Fehler)

      END ! rOlfPold

      SUBROUTINE tOlfPold (j, tPar, Fehler)
C     -------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      CHARACTER*(*) Fehler
      CHARACTER     tPar(MP)*40
      INTEGER       j

      CALL tOlfP (j, 'fil', tPar(1), Fehler)
      CALL tOlfP (j, 'tit', tPar(2), Fehler)
      CALL tOlfP (j, 'doc', tPar(3), Fehler)
      CALL tOlfP (j, 'dir', tPar(4), Fehler)
      CALL OlfCnuP (j, 'x', tPar(5), tPar( 8), Fehler)
      CALL OlfCnuP (j, 'y', tPar(6), tPar( 9), Fehler)
      IF (tPar(7).ne.' ') CALL OlfCnuP (j, 'z1', tPar(7), tPar(10),
     *                                  Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! tOlfPold

C  --------------------------------------------------------------------
C  read old IED format
C  --------------------------------------------------------------------
            ! The IED format was used until May 92.
            ! It extended the ILL's SQW format (R.Gosh) to
            ! which it was kept two-way-compatible.

      SUBROUTINE ReadIED (File, K, iPar, rPar, tPar,
     *                    z, M, n, X, Y, D, Fehler)
C     ----------------------------------------------
            ! (JWu)15.1.91, 13.2.91 as ReadSpectrum (-> ReadCompact)
            ! abbreviated nov96
         ! If K=1, the parameter arrays will be set;
         ! if K>1, they will be checked for consistency.
         ! The reading of the data blocks is done in ReadIEDsub

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)
      INCLUDE      'i_dim.f'
      INTEGER       iPar(MP), iParIn(MP)
      REAL*8        rPar(MP), rParIn(MP)
      CHARACTER*40  tPar(MP), tParIn(MP)
      REAL*8        X(*), Y(*), D(*)
      CHARACTER*(*) File, Fehler
      CHARACTER     cr2*2, cl2*2

      IF (Fehler.ne.'&ff') CALL Absturz ('ReadSpectrum','err on entry')

      IF (K.lt.1) CALL Absturz ('ReadSpectrum', 'K < 1 on entry')
      qSet = (K.eq.1)

      IF (qSet) THEN
C  Open file :
         CALL OpenFile (31, File, 'dat', 'l', Fehler)
         IF (Fehler.ne.'&ff') RETURN
C  Set Par to 0 :
         DO i = 1, MP
            iParIn(i) = 0
            rParIn(i) = 0.
            tParIn(i) = '&-'
            ENDDO
         ENDIF

      tPar(1) = File
      CALL ReadIEDsub  (31, iParIn(9) .ne. 0, iParIn(1), iParIn(11),
     *                  iParIn(12),
     *                  rParIn(1), rParIn(2), rParIn(3), rParIn(9),
     *                  rParIn(4), rParIn(5), rParIn(6), rParIn(7),
     *                  rParIn(8),
     *                  tParIn(2), tParIn(3), tParIn(5), tParIn(6),
     *                  tParIn(7),
     *                  tParIn(8), tParIn(9), tParIn(10),
     *                  tPar(11), M, X, Y, D, Fehler)

      IF (Fehler.eq.'&eof') THEN
         IF (qSet) THEN
            CALL Insert (Fehler, 1, 'file is completely empty')
         ELSE
            Close (31) ! Reached end of file
            ENDIF
         RETURN
      ELSEIF (Fehler.ne.'&ff') THEN
         ! Error occured :
         CALL Insert (Fehler, 1, 'spectrum '//cr2(K)//' ')
         RETURN
         ENDIF

      n = iParIn(1)

      IF (qSet) THEN
c96/7         iz = iParNumber(tParIn(7))
c      ELSE
         iz = iPar(10)
         ENDIF
      iParIn(10) = iz
      z = rParIn(iz)
      rParIn(iz) = 0.0

      IF (qSet) THEN
         DO i = 1, MP
            iPar (i) = iParIn (i)
            rPar (i) = rParIn (i)
            tPar (i) = tParIn (i)
            ENDDO
      ELSE
         IF (iParIn(1).ne.iPar(1)) iPar(1) = 0  ! nC
         DO i = 3,MP
            IF (iParIn(i).ne.iPar(i)) THEN
               Print *, 'irregularity in iPar at i,K = ',
     *                   i, K, ' : value = ', iParIn(i)
               ENDIF
            ENDDO
         DO i = 1,MP
            IF (rParIn(i).ne.rPar(i)) THEN
               Print *, 'irregularity in rPar '//cl2(i)//
     *               ' at K ='//cr2(K)//', value = ', rParIn(i)
               rPar(i) = 0.
               ENDIF
            ENDDO
         ENDIF

      iPar ( 2) = K
      iPar (10) = iz

      END ! ReadIED

      SUBROUTINE ReadIEDsub (NR, qToF, nC, ntpar, iSym,
     *                angle, E0, Q0, z, deltaE, deltaTau, deltaK,
     *                Temp, AMasse, HeadLine, Doc,
     *                xCoord, yCoord, zCoord, xUnit, yUnit, zUnit,
     *                Comment, M, X, Y, D, Fehler )
C     -------------------------------------------------------------------
            ! renewed 6.1.91, error messages 21.1.91

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      CHARACTER     muell*40, Fehler*(*), cr3*3
      CHARACTER*40  xUnit, yUnit, zUnit, xCoord, yCoord, zCoord,
     *              Headline, Doc, Comment(15)
      REAL *8       X(*), Y(*), D(*), Q0
      INTEGER       nKopf(8)

      nTrial = 0
 11   CONTINUE
      nTrial = nTrial + 1
      IF (nTrial.gt.2) THEN
         Fehler = 'ReadIED/ trapped in ToF-BS bubble'
         RETURN
         ENDIF

C     Zone 0 :
      Read (NR,'(8i5)', end=80, err=90)  nKopf

C     Zone 1 :
      Read (NR,'(1x,a39)', err=91)  Headline
      IF (nKopf(2).gt.1) THEN
      Read (NR,'(1x,a39)', err=91)  Doc
      ELSE
         Doc = ' '
         ENDIF

C     Zone 2 :
      z = 0.0
      Read (NR,'(1x,f6.2,f8.3,f8.4,f9.3,f6.1,i2,/'//
     *          '3x,e13.5,3f8.4)', err=92)
     *                   angle, E0, Q0, Temp, AMasse, iSym,
     *                   z, DeltaE, DeltaTau, DeltaK

C     Zone 3 :
      IF     (nKopf(4).eq.0) THEN
C        Special procedure for SQW/CrossX data :
         xCoord = 'w'
         xUnit  = 'meV'
         yCoord = 'S(q,w)'
         yUnit  = 'meV-1'
         zCoord = '2Th'
         zUnit  = ' '
      ELSEIF (nKopf(4).eq.3) THEN
C        For all other data : JWu-format :
         Read (NR,'(1x,a27,a12)', err=93)  xCoord, xUnit
         Read (NR,'(1x,a27,a12)', err=93)  yCoord, yUnit
         Read (NR,'(1x,a27,a12)', err=93)  zCoord, zUnit
      ELSE
         GOTO 934
         ENDIF

C     Zone 4 :
      DO i = 1, nKopf(5)
         IF (i.le.10) THEN
            Read (NR,'(1x,a39)', err=94)  Comment(i)
         ELSE
            Read (NR,'(1x,a39)', err=94)  muell
            ENDIF
         ENDDO
      ntpar = 10 + min0 (nKopf(5), 15)

C     Zone 5, 6 :
      DO i = 1, nKopf(6)
         Read (NR,'(1x,a39)', err=95)  muell
         ENDDO
      DO i = 1, nKopf(7)
         Read (NR,'(1x,a39)', err=96)  muell
         ENDDO

C     Zone 7 :
      nC = nKopf(8)
      IF (nC.gt.M) GOTO 970

         ! for standard SQW/CrossX as well as for own data :
         ! by order of mufti :

         ! qToF      <=> X /   1 <=> 'f9.5' <=> descending X
         ! .not.qToF <=> X /1000 <=> 'f9.6' <=> ascending X

C   DIVISION / 1000 NICHT PROGRAMMIERT !

      IF (.not.qTOF) THEN
 103        FORMAT(6X,F9.6,E13.5,E12.4)
         Read (NR,103, err=971) X(1),Y(1),D(1)
         IF (nC.gt.1) THEN
            Read (NR,103, err=972) X(2),Y(2),D(2)
            qUpwards=(X(1).LE.X(2))
            IF (.NOT.qUpwards) THEN
               qTOF = .true.
c               Print *, ' apparently reading TOF-data'
               Rewind (NR)
               GOTO 11
               ENDIF
            ENDIF
      ELSE
 104        FORMAT(6X,F9.5,E13.5,E12.4)
         Read (NR,104, err=971) X(1),Y(1),D(1)
         IF (nC.gt.1) THEN
            Read (NR,104, err=972) X(2),Y(2),D(2)
            qUpwards=(X(1).LE.X(2))
            IF (qUpwards) THEN
               qTOF = .false.
c               Print *, ' apparently reading backscattering-data'
               Rewind (NR)
               GOTO 11
               ENDIF
            ENDIF
         ENDIF

      IF (.not.qUpwards .and. nC.gt.1) THEN
         X (nC+1-1)=X (1)
         Y (nC+1-1)=Y (1)
         D(nC+1-1)=D(1)
         X (nC+1-2)=X (2)
         Y (nC+1-2)=Y (2)
         D(nC+1-2)=D(2)
         ENDIF

      DO iCh=3,nC
         IF (qUpwards) THEN
            i=iCh
         ELSE
            i=nC+1-iCh
            ENDIF

         IF (qTOF) THEN
            Read (NR,104, err=97) X(i), Y(i), D(i)
         ELSE
            Read (NR,103, err=97) X(i), Y(i), D(i)
            ENDIF

         ENDDO

C  For some units restoration of factor 100 or 1000 :
      IF (xUnit.eq.'hK') THEN
         xUnit = 'K'
         DO i = 1, nC
            X(i) = X(i) * 100
            ENDDO
         ENDIF
      IF (xCoord.eq.'2Th/100') THEN
         xCoord = '2Th'
         DO i = 1, nC
            X(i) = X(i) * 100
            ENDDO
         ENDIF

      DO i = 1, nC
         IF (D(i).lt.0.) THEN
            D(i) = 0.
            IF (.not.qWarn) THEN
               Print *, ' data have negative error bars'
               Print *, ' minderwertige Software benutzt ?'
               CALL Gong (7)
               ENDIF
            qWarn = .true.
         ENDIF
      ENDDO
      RETURN

 80   CONTINUE
      Fehler = '&eof'
      RETURN
 90   CONTINUE ! check whether format<>8i5 or file empty
      Read (NR, '(a)', err=901) muell
      Fehler = ' is no IED file (bad format in zone 0)'
      RETURN
 901  Fehler = ' is empty'
      RETURN
 91   Fehler = ' bad format in zone 1'
      RETURN
 92   Fehler = ' bad format in zone 2'
      RETURN
 93   Fehler = ' bad format in zone 3'
      RETURN
 934  Fehler = ' has bad format : #lines (zone 3) <> 0,3'
      RETURN
 94   Fehler = ' bad format in zone 4'
      RETURN
 95   Fehler = ' bad format in zone 5'
      RETURN
 96   Fehler = ' bad format in zone 6'
      RETURN
 970  Fehler = ' spectrum too long'
      RETURN
 971  Fehler = ' bad format in first data line'
      RETURN
 972  Fehler = ' bad format in 2nd data line'
      RETURN
 97   Fehler = ' bad format in zone 7 (data line'//cr3(iCh)//')'
      RETURN
      END ! ReadIEDsub
