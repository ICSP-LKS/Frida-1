C  ====================================================================
C
C     Library  WuL :  General FORTRAN Library
C     Module   L5  :       File Access, Help
C
C  ====================================================================

C     Contents :
C        5.1.  Open File Commands :
C                       these commands ultimately call OpenExe from
C                       l0xxx which depends on the operating system.
C                 AskOpenFile, AskDatOpenFile, OpenFile, OpenDatFile,
C                 FileNameMake
C        5.2.  File Read-in :
C                 ReadFile
C        5.4.  Help system :
C                 ErsteHilfe
C        5.5.  Account for CPU time (noch hier) :
C                 ChaUnit, ListCPS

C     History :
C        aug93 Help
C        jan93 System calls now in l0xxx
C        oct92 File access for UNIX(bsd)
C        90/92 File access for VMS
C        89/90 File access, copy, delete, timer for NOS/VE
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  ====================================================================
C  L5.1.  Open File Commands :
C  ====================================================================

       ! Parameters for file access :
       ! DefExt       : Filename will be <Name>.<DefExt>
       !                if Name doesn't contain a '.'
       ! stat = 'a'   : alt :  muss schon existieren
       !        'l'   : ?   :  the same as 'a' ??
       !        'n'   : neu :  darf noch nicht existieren
       !        'n?'  :        ueberschreiben nur nach Rueckfrage
       !        'n!'  :        ueberschreiben, aber Warnung
       !        'e'   : egal : in jedem Fall oeffnen
       ! acc  = 'seq' : sequential
       !        'app' : sequential, append
       !        'dir' : direct (or random) access
       ! cc   = '-cc' : no carriage control
       !        'for' : fortran interpretation of first character
       !        'lis' : list carriage control
       !       (ignored for old files)
       ! irec = 1-2499: recordlength=irec, recordtype=fixed, formatted
       !        0     : recordlength and -type undetermined,
       !                  form(Def)=formatted(if seq), unformatted(if dir)
       !        -1    : ... undetermined, unformatted

      SUBROUTINE AskOpenFile (Frage, Nr, Name, Dir, DefExt, stat)
C     -----------------------------------------------------------
         CHARACTER *(*) Frage, Name, Dir, DefExt, stat
         CALL AskOpenDatFile (Frage, Nr, Name, Dir,
     *                        DefExt, stat, 'seq', 'lis', 0)
         END ! AskOpenFile

      SUBROUTINE AskOpenDatFile (Frage, Nr, Name, Dir, DefExt,
     *                           stat, acc, cc, irec)
C     --------------------------------------------------------
            ! JWu13dec90
         ! ask for a file and try to open it.
      IMPLICIT LOGICAL (q)
      CHARACTER *(*) Frage, Name, Dir, DefExt, stat, acc, cc
      CHARACTER *80  File, aus, quest, ein, Fehler

      Fehler = '&ff'
c     Print *, '  ' ! :: NY - Absturzstelle

      aus = Frage
 11   CONTINUE
      IF (aus.ne.'&noq') THEN
         quest = aus
         IF (Dir.ne.' ' .and. Dir.ne.'&nod') THEN
            CALL Append (quest, ' (dir='//Dir(1:lenU(Dir))//')')
            ENDIF
         CALL FrageHD (quest, ein, Name)
         IF     (ein.eq.'?') THEN
            Print *, ' INPUT HELP/'
            Print *,
     * '    the program is going to open an external data file'
            Print *,'    by default, the external file name will be "',
     *                       Name(1:lenU(Name)), '"'
            Print *, '    enter any name to overwrite this default'
            Print *,
     * '    enter "\" if you don''t want to open a file at all'
            IF (Dir.ne.'&nod') THEN
            IF (Dir.eq.' ') THEN
            Print *,
     * '    by default, files will be opened in the directory ',
     * 'from which this program was started'
            ELSE
            Print *,
     * '    by default, files will be opened in directory ',
     * Dir(1:lenU(Dir))
            ENDIF
            Print *, '    enter "cd <dir_name>" to change this default'
            ENDIF
            GOTO 11
      ! directory command given ?
         ELSEIF (ein(1:3).eq.'cd ' .and. Dir.ne.'&nod') THEN
            Dir = ein(4:len(ein))
            CALL DelLeft (Dir)
            GOTO 11
         ELSEIF (ein.eq.'y' .or. ein.eq.'1') THEN ! haeufiger Eingabefehler
            CALL Gong (2)
            IF (.not.qAsk(' Indeed '//ein(1:lenU(ein))//' ?')) GOTO 11
            ENDIF
         Name = ein
         ENDIF
      IF (Name.eq.' ') RETURN ! user decided to open no file

      IF (Dir.ne.' ' .and. Dir.ne.'&nod') THEN
         File = Dir(1:lenU(Dir)) // '/' // Name
      ELSE
         File = Name
         ENDIF
      CALL OpenDatFile (Nr, File, DefExt, stat, acc, cc, irec, Fehler)

      IF (Fehler.ne.'&ff') THEN
         IF (Fehler.ne.'&rifusato') THEN
            CALL Gong(3)
            Print *, Fehler
            ENDIF
         Fehler = '&ff'
         IF (aus.eq.'&noq') aus = ' Enter other filename'
         Name = ' ' ! no default for second attempt, 3jun92
         GOTO 11
         ENDIF

      END ! AskOpenDatFile

      SUBROUTINE OpenFile (Nr, Name, DefExt, status, Fehler)
C     ------------------------------------------------------
         ! abbreviated call to OpenDatFile
      CHARACTER *(*) Name, DefExt, Fehler, status

      CALL OpenDatFile (Nr, Name, DefExt, status, 'seq',
     *                  'lis', 0, Fehler)

      END ! OpenFile

      SUBROUTINE OpenDatFile (Nr, Name, DefExt, stat, acIn,
     *                        ccIn, irec, Fehler)
C     -------------------------------------------------------------------------
            ! J.Wuttke 23.11.89, split 10dec92.

      IMPLICIT LOGICAL (q)
      CHARACTER *(*) Name, DefExt, Fehler, stat, acIn, ccIn
      CHARACTER *20 acc, cc
      CHARACTER *80 fnam

      CALL FileNameMake (fnam, Name, DefExt)

      IF     (acIn.eq.'seq') THEN
         acc = 'sequential'
      ELSEIF (acIn.eq.'app') THEN
         acc = 'append'
      ELSEIF (acIn.eq.'dir') THEN
         acc = 'direct'
      ELSE
         CALL Absturz ('OpenDatFile', 'acIn o.o.r.')
         ENDIF

      IF     (ccIn.eq.'-cc') THEN
         cc = 'none'
      ELSEIF (ccIn.eq.'for') THEN
         cc = 'fortran'
      ELSEIF (ccIn.eq.'lis') THEN
         cc = 'list'
      ELSE
         CALL Absturz ('OpenDatFile', 'cc o.o.r.')
         ENDIF

      IF (irec.lt.-1 .or. irec.gt.2499)
     *   CALL Absturz ('OpenDatFile', 'irec o.o.r.')

      CALL OpenExe (Nr, fnam, stat, acc, cc, irec, Fehler)

      END ! OpenDatFile

      SUBROUTINE FileNameMake (fnam, Name, DefExt)
C     --------------------------------------------
            ! JWu 5nov92, mainly for ".dat" under UNIX.
         ! Compose filename fnam = Name // '.' // DefExt
         ! except if Name contains already a '.'
         ! or is ending with '\' which means no-extension
         ! or if DefExt = '&noext'
      ! to be called exclusively by AskOpenDatFile

      IMPLICIT LOGICAL (q)
      CHARACTER*(*) fnam, Name, DefExt

      lf = len(fnam)
      ln = lenU(Name)

      IF (Name(ln:ln).eq.'\') THEN
         fnam = Name(1:ln-1) ! no extension
         RETURN
         ENDIF

      i  = 0
      qExt = (DefExt.eq.'&noext')
      qDir = .false.
 1    CONTINUE ! search for '.' and '/'
         i = i + 1
         IF (Name(i:i).eq.'[') THEN
            DO ii=i+1, ln
               IF (Name(ii:ii).eq.']') THEN
                  i = ii
                  GOTO 1
                  ENDIF
               ENDDO
            ENDIF
         IF     (Name(i:i).eq.'.') THEN
            IF     (qSubStrEq(Name,i-1,'\')) THEN ! UNIX-escaped
            ELSEIF (qSubStrEq(Name,i+1,'.')) THEN ! .. gilt nicht (23jun93)
               i = i + 1
            ELSE ! no exceptions : found valid '.'
               qExt = .true.
               ENDIF
            ENDIF
         IF (i.lt.ln) GOTO 1

C  Extension :
      IF (qExt) THEN
         fnam(1:lf) = Name
      ELSE
         fnam(1:lf) = Name(1:ln) // '.' // DefExt
         ENDIF

      END ! FileNameMake

C  ====================================================================
C  L5.2.  File Read-in, Write-Out :
C  ====================================================================

      SUBROUTINE ReadFile (nu, line, ML, nL, ceoi, Fehler)
C     ----------------------------------------------------
            ! JWu  7.12.89, revu 3mar93
         ! Read from unit nu into array line(1..nL).
         ! If ceoi='&eof' read until end-of-file,
         ! else read until a line starts with ceoi;
         ! this line will not be written to line.

      IMPLICIT LOGICAL (q)
      CHARACTER*(*) Fehler, ceoi, line(ML)
      CHARACTER     cl7*7

      IF (qErrEntry('ReadFile', Fehler)) RETURN

      DO i = 1, ML
         Read (nu, '(a)', end=8, err=9) line(i)
         IF (line(i)(1:len(ceoi)).eq.ceoi) THEN
            nL = i - 1
            IF (ceoi.eq.'&eof') CALL Compose2 (Fehler,
     *         'line '//cl7(i), ' contains reserved string &eof')
            RETURN
            ENDIF
         ENDDO
      CALL Compose2 (Fehler, 'file too long: '//cl7(i), ' lines read')
      RETURN

 8    CONTINUE ! eof
      IF (ceoi.eq.'&eof') THEN
         nL = i - 1
      ELSE
         Fehler = 'reached end-of-file without finding '//ceoi
         ENDIF
      RETURN
 9    CONTINUE ! error condition
      Fehler = 'error while reading line '//cl7(i)
      RETURN

      END ! ReadFile

C  ====================================================================
C  L5.4.  The Help System :
C  ====================================================================

      SUBROUTINE ErsteHilfe (CalFile, CalKey)
C     ---------------------------------------
            ! JWu 13aug93
      IMPLICIT LOGICAL (q)
      EXTERNAL      RunPreset

      PARAMETER    (MHF=20,MHL=2000)
      CHARACTER*(*) CalFile, CalKey
      CHARACTER*40  HFile(MHF), DirAux
      CHARACTER*80  HLine(MHL), fname, blabla, Fehler
      DIMENSION     NHL(0:MHF)

      DATA  nHF / 0 /

C  Initialisations :
      NHL(0) = 0 ! NHL(j)=last line of file j
      Fehler = '&ff'

C  If necessary, open help-file :
      DO iHF = 1, nHF
         IF (CalFile.eq.HFile(iHF)) GOTO 11
         ENDDO
      IF (nHF.ge.MHF) THEN
         Print *, ' Help/', CalFile
         Print *, ' sorry, too many help-files open'
         RETURN
         ENDIF
      ! open new help-file :
      CALL ExeML ('\p dir-hlp', DirAux)
      CALL Compose2 (fname, DirAux, CalFile)
      CALL OpenFile (37, fname, 'hlp', 'a', Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, ' Help/', CalFile
         Print *, ' sorry, cannot open help-file ', fname
         RETURN
         ENDIF
      HFile(nHF+1) = CalFile
      ! read in :
      CALL ReadFile (37, HLine(NHL(nHF)+1), MHL, nin, '&eof', Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL FehlerGong (Fehler,1)
         RETURN
         ENDIF
      NHL(nHF+1) = NHL(nHF)+nin
      nHF = nHF + 1
      iHF = nHF
      Close (37)

C  Help file has been read, now use it :
 11   CONTINUE

      ! find node CalKey :
      DO i = NHL(iHF-1)+1, NHL(iHF)
         IF (HLine(i)(1:1).eq.'*') THEN
            IF (HLine(i).eq.'*'//CalKey) THEN
               iNode = i
               GOTO 31
               ENDIF
            ENDIF
         ENDDO
      ! nothing found :
      Print *, ' Help/', CalFile
      Print *, ' sorry, there is no documentation on ', CalKey
      RETURN

      ! output one node :
 31   CONTINUE
      DO i = iNode+1, NHL(iHF)
         IF (HLine(i)(1:1).eq.'*') THEN ! begin of next node
            RETURN
            ENDIF
         Print '(a79)', HLine(i)
         ENDDO

      END ! ErsteHilfe

C  ====================================================================
C  L5.5.  Account for CPU time :
C  ====================================================================

      BLOCK DATA UnitPreset
C     ---------------------

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      CHARACTER      UnitNam*40

      COMMON /CPUZZ/ CPS(1:100), UnitNam(1:100), ialt, iuralt

      DATA           iuralt /0/, ialt /0/,
     *               UnitNam / 100* '?' /, CPS / 100* 0. /

      END ! UnitPreset

      SUBROUTINE ChaUnit (ineuIn, Name)
C     ---------------------------------
           ! JW 26.10.89, renovation 13apr92
       ! Wechselt das Konto, auf das die CPU-Sekunden gezaehlt werden.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      EXTERNAL       RunPreset

      CHARACTER      Name*(*), UnitNam*40
      COMMON /CPUZZ/ CPS(1:100), UnitNam(1:100), ialt, iuralt

C  Initialization ?
      IF (ialt.eq.0) THEN
         SecIst     = CPUseconds (0) ! reset

C  Increment last account :
      ELSE
         SecIst     = CPUseconds (1)
         delta      = SecIst - SecAlt
         CPS (ialt) = CPS (ialt) + delta  !  add to account
         ENDIF

C  Back in history ?
      IF (ineuIn.eq.-1) THEN
         ineu = iuralt
      ELSE
         ineu = ineuIn
         ENDIF
      iuralt     = ialt
      ialt       = ineu

      IF (ineu.eq.0) RETURN ! stop or interrupt

C  New account ?
      IF (ineu.lt.1 .or. ineu.gt.100) THEN
         Print *, ' iuralt ialt ineu ineuIn :', iuralt,ialt,ineu,ineuIn
         CALL Absturz ( 'ChaUnit', 'new unit number out of range' )
         ENDIF
      IF (UnitNam(ineu).eq.'?') THEN
         ! open new account :
         UnitNam(ineu) = Name
      ELSEIF (ineuIn.gt.0) THEN
         ! check identity of old account :
         IF (Name.ne.UnitNam(ineu)) THEN
            Print *, ' old name : ', UnitNam(ineu)
            Print *, ' new name : ', Name
            CALL Absturz ('ChaUnit','Name geaendert: ')
            ENDIF
         ENDIF

C  Restart counting :
      SecAlt     = CPUseconds (1)

      END ! ChaUnit

      SUBROUTINE ListCPS (NU)
C     ----------------------- ! JWu 26.10.89, 14apr92
         ! Listet verbrauchte CPU-Zeit auf

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      EXTERNAL       RunPreset

      CHARACTER      UnitNam *40
      COMMON /CPUZZ/ CPS(100), UnitNam(100), ialt, iuralt
      COMMON / System   / qCPUT

      IF (.not.qCPUT) RETURN

      CALL ChaUnit (0,' ') ! ausschalten
      SecTot = CPUseconds (1)

      IF (NU.le.0) RETURN

      write (NU,*)
      write (NU,*) ' CPU-Zeit-Verbrauch : '
      write (NU,*)
      CPsum = 0.
      DO iu = 1,100
         IF (UnitNam(iu).ne.'?') THEN
            write (*,'(a,i2,a,a17,a,f8.2)')
     *      ' Im Programmteil ', iu,' = ',UnitNam(iu),' : ',CPS(iu)
            CPsum = CPsum + CPS(iu)
            ENDIF
         ENDDO
      write (NU,*)
      write (NU,'(a,17x,a,f8.2)')
     *      ' Summe                ',' : ',CPsum
      write (NU,'(a,17x,a,f8.2)')
     *      ' Total                ',' : ',SecTot
      write (NU,*)

      END ! ListCPS
