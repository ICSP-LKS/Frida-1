C  ====================================================================
C
C     Library  WuLib :  General FORTRAN Library
C     Module   L0sos :  Settings and Calls Specific to SUN-OS
C
C  ====================================================================

C     Outline :
C        This module contains block data and executable subroutines
C        that must be rewritten for each operating system

C     Contents :
C        0.1.  Settings for Library Usage (interactive mode) :
C                 RunPreset
C        0.2.  Settings for graphics :
C                 GraPreset, FormPreset
C        0.4.  File access :
C                 OpenExe
C        0.5.  Account for CPU time :
C                 CPUseconds
C        0.6.  Get real time :
C                 Datum, Zeit
C        0.8.  Arithmetic functions :
C                 dpow0(REAL*8)

C  ====================================================================
C  L0sun.1.  Settings for Library Usage :
C  ====================================================================

      BLOCK DATA RunPreset
C     -------------------- ! 20may94 (ex TermPreset)
         ! - DirAux :   Directory for ErsteHilfe (l5),
         ! - Terminal : linefeed character "+" ?
         !              Gong ?
         ! - Impltn :   OS -> Absturz (l3)

      IMPLICIT LOGICAL (q)
      CHARACTER  DirAux*80, OS*4

      COMMON / RunDir   / DirAux
      COMMON / Terminal / qPlus, qDroehn
      COMMON / Impltn   / OS

      DATA  DirAux / '/home/llb2/wuttke/exe' / ! Saclay (~ doesn't work)
C      DATA  DirAux / '~jwuttke/for' / ! E13

      DATA  qPlus, nDroehn / .false., 1 / ! Falco

C      DATA  OS / 'X11' /
      DATA  OS / 'SOS' / ! Sun OS, vt100 mode

      END ! RunPreset

C  ====================================================================
C  L0sun.2.  Settings for graphics :
C  ====================================================================

      BLOCK DATA GraPreset
C     -------------------- ! 23feb93
         ! Presets for g1 :
	 !    FileAux   =  location of PostScript macros
	 !    qWindow   =  graphic has its own window
	 !    qGTOverlay=  graphic and dialog in different sections
         !                    of one window (Pericom)
	 !    qGToldAdr =  reduction Bytes->Bits works
	 !    iGraP     =  0..1 : less or more '+' in graphic commands
	 !    iEsc      =  0..1 : switch between windows with ESC ?

      IMPLICIT LOGICAL (q)
      CHARACTER  FileAux*80

      COMMON / GraDir / FileAux   ! <-> g1
      COMMON / GraTerm/ qWindow, qGTOverlay, qGToldAdr
      COMMON / GraPlus / iGraP, iEsc

      DATA  FileAux / '/home/llb2/wuttke/for/g3.ps' / ! Saclay
C      DATA  FileAux / '~jwuttke/for/g3.ps' / ! E13

      DATA  qWindow /.false./, qGTOverlay /.false./,
     *      qGToldAdr /.false./, iGraP / 0 /, iEsc / 0 /    ! Falco
C      DATA  qWindow /.true./, qGTOverlay /.false./,
C     *      qGToldAdr /.false./, iGraP / 0 /, iEsc / 1 /   ! X-windows


      END ! GraPreset

      BLOCK DATA FormPreset
C     --------------------- ! 23feb93/rep. 6jun93
         ! Default choice for graphic windows in g2

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h, o-p, r-z)

      COMMON / Scroll / iScroIst, iScroGra
      COMMON / Format / kScrF, kLasF, kPS_F

      DATA    iScroIst, iScroGra / 0, 20 /
      DATA    kScrF, kLasF, kPS_F / 9, 1, 3 /     ! 1 window
c      DATA    kScrF, kLasF, kPS_F / 5, 1, 3 /     ! multiwindow

      END ! FormPreset

C  ====================================================================
C  L5sun.4.  File access :
C  ====================================================================

      SUBROUTINE OpenExe (Nr, fnam, stat, acc, cc, irec, Fehler)
C     -----------------------------------------------------------
            ! oct-dec92 UNIX(bsd) version.
         ! voir l5vax.for pour une version commentee

         IMPLICIT LOGICAL (q)
         CHARACTER *(*) fnam, Fehler, stat, acc, cc
         CHARACTER *80 aus

         IF (irec.lt.-1 .or. irec.gt.2499)
     *      CALL Absturz ('OpenDatFile', 'irec o.o.r.')

         IF     (stat.eq.'a' .or. stat.eq.'l') THEN
            ! Open existing file. READONLY doesn't exist.
            IF     (irec.ge.0)  THEN
               Open (unit=Nr, file=fnam, status='old', err=11,
     *            access=acc)
            ELSEIF (irec.eq.-1) THEN
               Open (unit=Nr, file=fnam, status='old', err=11,
     *            access=acc, form='unformatted')
               ENDIF
            RETURN
 11         Open (unit=Nr, file=fnam, status='old', err=12)
            Fehler = 'wrong access keyword for old file '//fnam
            RETURN
 12         Fehler = 'cannot open old file '//fnam
            RETURN
         ELSEIF (stat.eq.'n') THEN
            ! Open new file - first make sure it doesnt exist.
            Open (unit=Nr, file=fnam, status='old', err=21)
               Fehler = 'File ' // fnam
               CALL Append (Fehler, ' already exists')
               Close (Nr)
               RETURN
 21         CONTINUE
            IF     (irec.eq.0)  THEN
               Open (unit=Nr, file=fnam, status='new', err=22,
     *            access=acc)
            ELSEIF (irec.eq.-1) THEN
               Open (unit=Nr, file=fnam, status='new', err=22,
     *            access=acc, form='unformatted')
            ELSE
               Open (unit=Nr, file=fnam, status='new', err=22,
     *            access=acc, form='formatted',
     *            recl=irec)
               ENDIF
            RETURN
 22         Fehler = 'cannot create file ' // fnam
            RETURN
         ELSEIF (stat.eq.'n?' .or. stat.eq.'n!' .or. stat.eq.'e') THEN
            ! Open as new, or overwrite, but take care :
            Open (unit=Nr, file=fnam, status='old', err=41, access=acc)
               ! does already exist :
               ! Because there are no version numbers,
               ! 'e' under UNIX acts as 'n!' under VMS,
               ! and 'n!' as 'n?' :
               IF     (stat.eq.'n?' .or. stat.eq.'n!') THEN
                  aus = ' File '//fnam
                  CALL Append (aus, ' exists already - overwrite ?')
                  qOver = qAsk(aus)
               ELSE
                  CALL Say2 (' WARNING/ file '//fnam,
     *               ' will be overwritten')
                  qOver = .true.
                  ENDIF
               IF (.not.qOver) THEN
                  Fehler = '&rifusato'
                  Close(Nr)
                  ENDIF
            Close(Nr)
 41         CONTINUE
            ! Open as with specific access mode and format :
            IF     (irec.eq.0)  THEN
               Open (unit=Nr, file=fnam, err=43, status='unknown',
     *            access=acc)
            ELSEIF (irec.eq.-1) THEN
               Open (unit=Nr, file=fnam, err=43, status='unknown',
     *            access=acc, form='unformatted')
            ELSE
               Open (unit=Nr, file=fnam, err=43, status='unknown',
     *            access=acc, form='formatted',
     *            recl=irec)
               ENDIF
               RETURN
 43         CONTINUE
            Fehler = 'cannot open file ' // fnam
            RETURN
         ELSE
            Print *, ' stat = ', stat
            CALL Absturz ('OpenExe','stat o.o.r.')
            RETURN
            ENDIF

         END ! OpenExe

C  ====================================================================
C  L5sun.5.  Account for CPU time :
C  ====================================================================

      REAL*8 FUNCTION CPUseconds (iReset)
C     -----------------------------------

      RETURN
      CALL Absturz ('CPUseconds', 'Under UNIX, use other tools')

C  Reset the timer ?
      IF (iReset.eq.0) THEN

C  else return the elapsed time :
      ELSE
         IF (Code.ne.2) CALL Absturz ('CPUseconds',
     *      'The timer hasn''t been initialized')
         ENDIF

      END ! CPUseconds

C  ====================================================================
C  L5sun.6.  Get real time
C  ====================================================================


      SUBROUTINE SysDat (id, im, iy)
C     ------------------------------
      INTEGER IDAT(3)
      CALL idate (IDAT) ! SPARC
c      CALL idate (im,id,iy) ! AUSPROBIEREN !!!
      id = IDAT(2)
      im = IDAT(1)
      iy = IDAT(3)
      END

      SUBROUTINE SysTim (ch, cm, cs)
C     ------------------------------
      CHARACTER s*8, ch*2, cm*2, cs*2
      CALL Time (s)
      ch = s(1:2)
      cm = s(4:5)
      cs = s(7:8)
      END

C  ====================================================================
C  L0sun.8.  Arithmetic functions
C  ====================================================================

      REAL*8 FUNCTION dpow0 (x,c)          ! x**c
	 IMPLICIT REAL*8 (a-h,o-p,r-z)
         IF ((c.le.0).and.(x.le.0)) THEN
            dpow0 = 0.
         ELSEIF (x.le.0) THEN
            IF (c.eq.dfloat(int(c))) THEN
               dpow0 = x ** int(c)
            ELSE
               dpow0 = 0.
               ENDIF
         ELSE
            h = c * dlog(x)
            dpow0 = dexp1(h)
            ENDIF
         END ! dpow0

      REAL*8 FUNCTION dpowi (x,m)
C     ---------------------------
         ! with integer argument m
      IMPLICIT REAL*8 (a-h,o-p,r-z)

      IF (x.eq.0) THEN
         dpowi = 0
      ELSE
         isig = 1
         IF (x.lt.0 .and. mod(m,2).eq.1) isig = -1
         dpowi = isig * dexp1(m*dlog(dabs(x)))
         ENDIF
      END ! dpowi