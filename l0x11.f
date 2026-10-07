C  ====================================================================
C
C     Library  WuLib :  General FORTRAN Library
C     Module   L0x11 :  Settings and Calls Specific to X11 (Ultrix/OSF)
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
C        0.7.  Site specific calls :
C                 Idol_Fil
C        0.8.  Arithmetic functions :
C                 dpow0(REAL*8), dpowi

C  ====================================================================
C  L0*.1.  Settings for Library Usage :
C  ====================================================================

      BLOCK DATA RunPreset
C     -------------------- ! 20may94 (ex TermPreset)

      IMPLICIT LOGICAL (q)
      CHARACTER         OS*4

      COMMON / Terminal / qPlus, nDroehn
      COMMON / Impltn   / OS

      DATA  qPlus, nDroehn / .false., 1 / ! Falco, X11
      DATA  OS / 'X11' /

      END ! RunPreset

C  ====================================================================
C  L5*.4.  File access :
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
C  L5*.5.  Account for CPU time :
C  ====================================================================

      REAL*8 FUNCTION CPUseconds (iReset)
C     -----------------------------------
         ! Under UNIX, use other tools
      CPUseconds = 0
      END ! CPUseconds

C  ====================================================================
C  L5*.6.  Get real time
C  ====================================================================

      SUBROUTINE SysDat (id, im, iy)
C     ------------------------------
C      CALL idate (im,id,iy) ! SPARC : american logic (month-day-year)
      id=12
      im=12
      iy=2000
      END

      SUBROUTINE SysTim (ch, cm, cs)
C     ------------------------------
      CHARACTER s*8, ch*2, cm*2, cs*2
      CHARACTER DateCH*8, TimeCH*10, ZoneCH*5
      INTEGER VALUES(8)
C      CALL Time (s)
C      ch = s(1:2)
C      cm = s(3:4)
C      cs = s(5:6)
      CALL DATE_AND_TIME(DateCH, TimeCH, ZoneCH, VALUES) !Artem replace CALL Time (s) with DATE_AND_TIME
      ch = TimeCH(1:2)
      cm = TimeCH(3:4)
      cs = TimeCH(5:6)
      END


C  ====================================================================
C  l0 (7) Site specific calls
C  ====================================================================

      SUBROUTINE Idol_Fil (Inst, irun, nC, nK, LongTit, HelgaPar, iC,
     *                     HelgaSca, MData, IData, istat)
C     ----------------------------------------------------------------
         ! exists only at the ILL
         !    [ link with /usr/local/idol/idol_sub.o -lrpcsvc -lsun ]

      REAL*4    HelgaPar(512), HelgaSca(512)
      INTEGER   IData(MData)
c      BYTE      HelgaByt(512) ! non-standard Fortran ??
      CHARACTER*(*) Inst, LongTit

      istat = -1000 ! subroutine does not exist at other sites

      END

C  ====================================================================
C  L0*.8.  Arithmetic functions
C  ====================================================================

      REAL*8 FUNCTION dpow0 (x,c)          ! x**c
C     ---------------------------
      IMPLICIT NONE
      REAL*8         x,c,h,dexp1

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
