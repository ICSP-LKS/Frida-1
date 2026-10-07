C  ====================================================================
C
C     Library  WuL :  General FORTRAN Library
C     Module   L3  :     Input/Output and MetaLanguage
C
C  ====================================================================

C     J.Wuttke - 1991ff.

C     Outline :
C        This module provides a unified user interface
C        underlying the various routines of L4.
C        All input passes through the routine Lies.
C        It is possible to escape from normal dialogue by
C        entering MetaLanguage commands which are executed by ExeML.
C        ExeML can also be directly called from a program,
C        e.g. for initializations.
C        All programs using this dialogue system should handle
C        fatal errors by calling Absturz; this routine dumps
C        all previous input before terminating program execution.

C     Contents :
C        3.1.  ANSI screen commands :
C                 ScrollOn/Off, PurgeSc, Cursor, PosLin, Gong
C        3.2.  Input/Output interface :
C                 ML_Ini, ExeML, i/rExeMLP, Lies, LinEdi, Sage
C        3.8.  Output Schnickschnack :
C                 Warte, Counter, Say2-7
C        3.9   Control :
C                 qErrEntry, FehlerGong, Absturz

C     File Access :
C              uses unit 61 : inout.ml

C     History :
C        oct97  scroll, prompt, Schreib eliminated
C        mar96  $-variables
C        may95  parameters
C        mar93  division L3/L4
C        jul92  formal help on all questions
C        91/92  first MetaLanguage elements

C  ====================================================================
C  L3.1  ANSI screen commands :
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE ScrollOn
C     -------------------
            ! Scroll, renamed 27aug91 because of interference
         CHARACTER seq *16
         seq = ' ' // char(27) // 'PMK1' // char(27) // char(92) ! reset scroll
         Write (6,'(a)') seq
         END ! Scroll

      SUBROUTINE ScrollOff
C     --------------------
         CHARACTER seq *16
         seq = ' ' // char(27) // 'PMK0' // char(27) // char(92)
      ! set NoScroll mode
         Write (6,'(a)') seq
         END ! NoScroll

      SUBROUTINE PurgeSc
C     ----------------
      !  purge screen
         Write (6,'(3a)') ' ', char(27), '[2J'      ! erase complete screen
         END ! PurgeSc

      SUBROUTINE Cursor (irel, ix, iy)
C     --------------------------------
         ! reposition cursor (JWu 25oct93)
         CHARACTER cv2 *2

         RETURN ! doesn't work on X-terminals and might even cause core dump

c         IF (irel.eq.0) THEN ! absolute
c            Print '(a)', char(27)//'['//cv2(iy)//';'//cv2(ix)//'H'
c         ELSE ! relative
c            IF (ix.gt.0) Print '(a)', char(27)//'['//cv2( ix)//'C' ! forward
c            IF (ix.lt.0) Print '(a)', char(27)//'['//cv2(-ix)//'D' ! backward
c            IF (iy.gt.0) Print '(a)', char(27)//'['//cv2( iy)//'B' ! down
c            IF (iy.lt.0) Print '(a)', char(27)//'['//cv2(-iy)//'A' ! up
c            ENDIF

         END ! Cursor

      SUBROUTINE PosLin (lz)
C     ----------------------
      !  Positionierung des Cursors nach Linie (lz) :
         CHARACTER cv2 *2
         Write (6,'(a)') '+' // char(27) // '[' // cv2(lz) // ';0H'
         END ! PosLin

      SUBROUTINE Gong (n)
C     -------------------
      IMPLICIT LOGICAL (q)
      EXTERNAL         RunPreset

      INCLUDE         'l_dim.f'
      CHARACTER        MLParNam*16, MLParVal*48, MLVarVal*16

      CHARACTER*80     Task(MLT), Done                          !Artem Task -> Task(MLT)

      COMMON / Terminal / qPlus, nDroehn
      COMMON /  ML  /  InUnit, qSave, Task, nTask, Done(MLD),   !Artem Task(MLT) -> Task
     *                 nDone, MLParNam(MLP), MLParVal(MLP), nMLPar,
     *                 MLVarVal(MLV), qEcho

      IF (nDroehn.le.0) RETURN

      IF (InUnit.eq.0) THEN ! 12jul94, 25oct95
         Print '(2a)', ' >>> GONG <<< ', char(7)

      ELSE
         IF (qPlus) THEN
            DO i = 1, n / nDroehn
               Print '(a)',  '+' // char(7)
               ENDDO
         ELSE
            DO i = 1, n / nDroehn
               Print '(2a)', char(0), char(7) ! NUL BEL ! NUL => no line feed
               ENDDO
            ENDIF

         ENDIF

      END ! Gong

C  ====================================================================
C  L3.2  Input/Output interface :
C  ====================================================================

      BLOCK DATA ML_Ini
C     -----------------

      IMPLICIT LOGICAL (q)

      INCLUDE         'l_dim.f'
      CHARACTER        MLParNam*16, MLParVal*48, MLVarVal*16
      CHARACTER*80     Task(MLT), Done  !Artem Task -> Task(MLT)

      COMMON /  ML  /  InUnit, qSave, Task, nTask, Done(MLD),  !Artem Task(MLT) -> Task
     *                 nDone, MLParNam(MLP), MLParVal(MLP),
     *                 nMLPar, MLVarVal(MLV), qEcho

      DATA   InUnit /5/, qSave /.true./, nTask /0/, nDone /0/,
     *       nMLPar /0/, qEcho /.true./, MLVarVal /MLV*' '/

      END ! ML_Ini

      SUBROUTINE ExeML (Text, Report)
C     -------------------------------
            ! JWu 14jan92
         ! Execute a metalanguage command. To be called from
         ! either the application program or Lies(Text).

      IMPLICIT LOGICAL (q)

      INCLUDE        'l_dim.f'
      CHARACTER       Text*(*), Report*(*),
     *                MLParNam*16, nam*16, MLParVal*48, MLVarVal*16
      CHARACTER       Datum*10, Zeit*8, cr3*3, cr4*4, cl4*4, cr5*5,
     *                cl6*6, cl8*8
      CHARACTER*80    Fehler, Word, Object, Task(MLT), Done, aus, ein !Artem Task -> Task(MLT)
      DIMENSION       qExe(MLT), qRecall(MLD)

      COMMON /  ML  / InUnit, qSave, Task, nTask, Done(MLD), !Artem Task(MLT) -> Task
     *                nDone, MLParNam(MLP), MLParVal(MLP), nMLPar,
     *                MLVarVal(MLV), qEcho

      Fehler = '&ff'
      Report = ' '

      IF (Text(1:1).ne.'\' .and. Text(1:1).ne.'$
     *    ') THEN  !  program error
         DO i = 1, min0 (10, len(Text))
            Print *, ' Text (i) = ichar : i, ichar ',
     *               i, ichar(Text(i:i))
            ENDDO
         Print *, ' again, as string :', Text
         CALL Absturz1 ('ExeML', 'No "\" command')
      ELSEIF (Text(1:3).eq.'\..' .or. Text(1:2).eq.'\%' .or.
     *        Text(1:2).eq.'\$' .or.  Text(1:2).eq.'\\') THEN
         CALL DelVonBis (Text, 1, 1) ! escape to original meaning
         RETURN
         ENDIF

      Object = Text
      CALL TakeVorDel (Object, Word, ' ')

      IF     (Word.eq.'\?') THEN ! help
         Print *,
     * '\ML MetaLanguage summary (type RETURN after each section) :'
         Print *
         Print *,
     * '1/ the following characters have special meaning whereever'//
     * ' they are used :'
         Print *, '   %            rest of line is comment'
         Print *, '   $0 $1        variables $0 $1 and so on'
         Print *, '   \            linebreak'//
     *            ' (backslash followed by blank or eol)'
         Print *, '   ..           continue input in next line'
         Print *, '   \% \$ \\ \..  escape to % $ \ ..'
         Read (*,*)
         Print *,
     * '2/ all other commands can be given at the beginning'//
     * ' of any input line'
         Print *, '   %             comment line'
         Print *, '   ?             request help on input syntax'
         Print *, '   ??            obtain more specific help'//
     *            '  or escape to ?'
         Print *, '   \?                this help'
         Print *, '   \?a               array dimensions'
         Print *, '   \* <n> <..>       n times ..'
         Print *,
     * '   \*=<list> <..>    n times, replacing * by text from list'
         Print *, '   \b <n>            bell'
         Print *, '   \dt               display time'
         Print *, '   \qui              quit'//
     *            ' (the same as ^C but more civilised)'
         Print *, '   \p <pnam> <val>   store parameter'
         Print *, '   \p <pnam>         recall parameter'
         Print *, '   \v <n>            video on/off'
         Read (*,*)
         Print *,
     * '3/ the following commands allow to save and to recall'//
     * ' previous input'
         Print *, '   \ls [<list>]      list session'
         Print *, '   \x <list>         recall from session'
         Print *, '   \ws <file>        write session to file (.ML)'
         Print *, '   \wsd              - file = ~/Dump.ML'
         Print *,
     * '   \i [-<n>] <file>  include task file (.ML) start at line n'
         Print *, '   \ie <file>        include and expand'          ! 7apr92
         Print *, '   \ied              - file = ~/Dump.ML'         ! 2jul92
         Print *, '   \lt [<list>]      list task'
         Print *, '   \st <list>        select task'
         Print *, '   \e                execute task'
         Print *, '   \et               empty task'
         Read (*,*)
         Print *,
     * '4/ the following commands handle the variables $0 to $9'
         Print *, '   $3=<text>         assign value to $3'
         Print *, '   \for $3=<list>    start loop'
         Print *, '   \rof              end loop'
         Read (*,*)
         Print *,
     * '5/ the following commands are of interest only in .ml files'
         Print *, '   \reset            next input from console'
         Print *

      ELSEIF (Word.eq.'\?a') THEN ! display array dimensions (24aug94)
         Print *, ' Maximum number of include lines = ', cr5(MLT)
         Print *, ' Size of session stack           = ', cr5(MLD)
         nDone = nDone - 1

      ELSEIF (Word.eq.'\dt') THEN ! display time (31oct91)
         Print *, Datum(9), ', ', Zeit(8)
         nDone = nDone - 1

      ELSEIF (Word.eq.'\b') THEN ! bell (7feb92)
        IF (Object.eq.' ') THEN
            n = 1
         ELSEIF (Object.eq.'#') THEN
            Fehler = 'abuse of # !!'
            n = 0
         ELSE
            CALL Fi1N (Object, n)
            IF     (Object.ne.'#') THEN
               Fehler = '\b : no positive integer given'
               n = 0
               ENDIF
            ENDIF
         CALL Gong (n)

      ELSEIF (Word.eq.'\qui') THEN ! quit (Reformationstag 1991)
         ! Overlay out, full scroll, goto last line :
         CALL ClearScroll ()
         STOP '\ML quit'

      ELSEIF (Word.eq.'\e') THEN ! execute
         IF (nTask.ge.1) THEN
            InUnit = 0
         ELSE
            Fehler = '\e : no tasks on stack'
            ENDIF

      ELSEIF (Word.eq.'\lt') THEN ! list task
         IF (Object.ne.' ') THEN
            CALL DecNList (Object, qExe, MLT, Fehler)
            IF (Fehler.ne.'&ff') GOTO 900
         ELSE
            CALL qSet (qExe, 1, nTask, 1, .true.)
            ENDIF
         IF     (nTask.lt.0) THEN
            Print *, ' PROGRAM ERROR/ nTask<0'
            nTask = 0
         ELSEIF (nTask.gt.MLT) THEN
            Print *, ' PROGRAM ERROR/ nTask>MLT'
            nTask = 0
         ELSEIF (nTask.eq.0) THEN
            Print *, '\lt : no tasks on stack'
         ELSE
            DO i = 1, nTask
               IF (qExe(i)) THEN
                  aus = '\lt '//cr3(i)//' : '//Task(i)
                  Print *, aus
                  ENDIF
               ENDDO
           ENDIF

      ELSEIF (Word.eq.'\et') THEN ! empty task
         nTask = 0

      ELSEIF (Word.eq.'\st') THEN ! select task
         CALL DecNList (Object, qExe, MLT, Fehler)
         IF (Fehler.ne.'&ff') GOTO 900
         DO i = 1, MLT
            IF (.not.qExe(i) .and. i.le.nTask) THEN
               nTask = nTask - 1
               DO ii = i, nTask
                  Task(ii) = Task(ii+1)
                  ENDDO
               ENDIF
            ENDDO

      ELSEIF (Word.eq.'\i' .or. Word(1:3).eq.'\ie') THEN ! include
         qExpand = (Word(1:3).eq.'\ie')
         IF (Word.eq.'\ied') Object = '~/Dump'
         IF (Object(1:1).eq.'-') THEN
            CALL DelVonBis (Object, 1, 1)
            CALL TakeVorDel (Object, Word, ' ')
            CALL DelLeft (Object)
            CALL Fi1N (Word, linestart)
            IF (Word.ne.'#') THEN
               Fehler = '\i -: no line number given'
               GOTO 900
               ENDIF
         ELSE
            linestart = 1
            ENDIF
         IF (Object.eq.' ') THEN
            Fehler = '\i : no include file given'
            GOTO 900
            ENDIF
         ! read input file, write in unchanged order on top of task-stack :
         IF (Object.ne.'>') THEN
            CALL OpenFile (61, Object, 'ml', 'a', Fehler)
            IF (Fehler.ne.'&ff') GOTO 900
            inline = 0
            DO ii = 1, linestart-1
               inline = inline + 1
               Read (61,'(a)',err=218,end=219) ein
               ENDDO
            ENDIF
         nTaskOld = nTask
 200     CONTINUE
            Read (61,'(a)',err=218,end=22) ein
            inline = inline + 1
            nTask = nTask + 1
            ! shift old task :
            DO ii = nTask, nTask-nTaskOld+1, -1
               Task(ii) = Task(ii-1)
               ENDDO
            IF (nTask.ge.MLT) THEN
c               Task(MLT) = '\i -'//cl6(linestart+nTask-1)//' '//Object
               Task(MLT) = '\i >'
c              Print *, '\i : include file too long, rest buffered: '//
c     *                  Task(MLT)
               Print *, '\i : include file too long, buffered from line'
     *                  //cl8(inline)
               GOTO 24
               ENDIF
            ! add new line :
            Task(nTask-nTaskOld) = ein
            GOTO 200
 218     CONTINUE            ! error
         Fehler = '\i : error while reading ML file'
         Close(61)
         GOTO 900
 219     CONTINUE               ! error
         Fehler = '\i -: eof before reaching start line'
         Close(61)
         GOTO 900
 22      CONTINUE ! regular end of file
         Close(61)
 24      CONTINUE
         Report = '\i + ML file read in'
         ! now read from Task
         InUnit = 0
         IF (.not.qExpand) THEN
            qSave = .false.
         ELSE
            nDone = nDone - 1  !  don't save the '\ie' command
            ENDIF

      ELSEIF (Word(1:3).eq.'\ws') THEN ! write session   (15jan92)
         nDone = nDone-1 ! don't save session commands
         IF     (Word.eq.'\wsd') THEN
            Object = '~/Dump'
         ELSEIF (Object.eq.' ') THEN
            Fehler = '\ws : no output file given'
            GOTO 900
            ENDIF
         CALL OpenDatFile (61, Object, 'ml', 'e', 'seq', 'lis',
     *                     0, Fehler)
         IF (Fehler.ne.'&ff') GOTO 900
         DO i = 1, nDone
            Write (61,'(a)',err=251) Done(i)(1:lenU(Done(i)))
            ENDDO
            GOTO 252
 251     CONTINUE ! error
         Fehler = '\ws : error while writing to session file'
         Close(61)
         GOTO 900
 252     CONTINUE ! regular end of information
         Report = '\ws + session file written'
         Close(61)

      ELSEIF (Word.eq.'\ls') THEN ! list session
         nDone = nDone-1 ! don't save session commands
         IF (Object.eq.' ') THEN
            Print '(a)', '\ls : the last entry was no. '//cl4(nDone)
            Print '(a)', '      use \ls <list> to see selected entries'
         ELSE
            CALL DecNList (Object, qRecall, nDone, Fehler)
            IF (Fehler.ne.'&ff') GOTO 900
            Print '(a)', '\ls : listing entries '//
     *            Object(1:lenU(Object))
            DO i = 1, nDone
               IF (qRecall(i)) THEN
                  Print '(a80)', cr4(i)//'/ '//Done(i)
                  ENDIF
               ENDDO
            ENDIF

      ELSEIF (Word.eq.'\v') THEN ! video
         IF (Object.eq.' ') THEN
            IF (qEcho) THEN
               Print *, '\ echo is ON --- type \v + to switch it off'
            ELSE
               Print *, '\ echo is OFF --- type \v - to switch it off'
               ENDIF
         ELSEIF (Object.eq.'-') THEN
            qEcho = .false.
         ELSEIF (Object.eq.'+') THEN
            qEcho = .true.
         ELSE
            Fehler = '\v : use + or - to switch echo on or off'
            ENDIF

      ELSEIF (Word.eq.'\x') THEN ! recall
         nDone = nDone-1 ! don't save session commands
         CALL DecNList (Object, qRecall, MLD, Fehler)
         IF (Fehler.ne.'&ff') GOTO 900
         DO i = 1, MLD
            IF (qRecall(i)) THEN
               nTask = nTask + 1
               Task(nTask) = Done(i)
               ENDIF
            ENDDO
         InUnit = 0

      ELSEIF (Word.eq.'\reset') THEN
         IF (InUnit.ne.5) THEN
            InUnit   = 5
            Print *, '\ ML reset'
            Print *, '\ control returned to console'
            Print *, '\ type "\et"   to void the stack'
            Print *, '\ type "\e"    to continue execution'
            ENDIF

      ELSEIF (Word.eq.'\*') THEN
         nDone = nDone-1 ! expand => don't save
         CALL Fi1N (Object, n)
         IF     (Object(1:1).ne.'#') THEN
            Fehler = '\* : no positive integer in "'//
     *           Object(1:lenU(Object))//'"'
            GOTO 900
            ENDIF
         ! shift the stack by n entries :
         nTask = nTask + n
         DO ii = nTask, n+1, -1
            Task(ii) = Task(ii-n)
            ENDDO
         ! now write the new command n times on the stack :
         CALL DelVonBis (Object, 1, 1) ! remove '#'
         CALL DelLeft (Object)
         DO i = 1, n
            Task(i) = Object
            ENDDO
         InUnit   = 0

      ELSEIF (Word(1:3).eq.'\*=') THEN ! 9jun95
         nDone = nDone-1 ! expand => don't save
         CALL DelVonBis (Word, 1, 3)
         IF (Word.eq.' ') THEN
            Fehler = '\*= no list given'
            GOTO 900
            ENDIF
         nTaskNew = 0
 821     CONTINUE
         nTaskNew = nTaskNew + 1
         ! shift the stack by 1 entry :
         nTask = nTask + 1
         DO ii = nTask, nTaskNew+1, -1
            Task(ii) = Task(ii-1)
            ENDDO
         ! new stack line :
         Task(nTaskNew) = Object
         CALL TakeVorDel (Word, aus, ',')
         CALL ReplaceC (Task(nTaskNew), '*', aus(1:lenU(aus)))
         IF (Word.ne.' ') GOTO 821
         InUnit   = 0

      ELSEIF (Word.eq.'\p') THEN
            ! ML parameter register (JWu 31may95)

         IF     (Object.eq.' ') THEN
            CALL Say2 ('There are '//cl4(nMLPar),
     *         ' ML parameters - type "\p ?" for a listing')
         ELSEIF (Object.eq.'?') THEN
            Print *, ' The following ML parameters are stored: '
            DO i = 1, nMLPar
               Print '(i3,2x,a16,2x,a50)', i, MLParNam(i), MLParVal(i)
               ENDDO
         ELSE
            CALL TakeVorDel (Object, nam, ' ')
            CALL DelLeft (Object)
            IF (Object.eq.' ') THEN
               ! display value
               DO i = 1, nMLPar
                  IF (MLParNam(i).eq.nam) THEN
                     Report = MLParVal(i)
                     GOTO 589
                     ENDIF
                  ENDDO
 589           CONTINUE
            ELSE
               DO i = 1, nMLPar ! par already stored ?
                  IF (MLParNam(i).eq.nam) THEN
                     MLParVal(i) = Object ! overwrite
                     GOTO 599
                     ENDIF
                  ENDDO
               IF (nMLPar.ge.MLP) THEN
                  Fehler =  'ML parameter register is full'
                  GOTO 900
                  ENDIF
               ! add new ML par :
               nMLPar = nMLPar + 1
               MLParNam(nMLPar) = nam
               MLParVal(nMLPar) = Object
 599           CONTINUE
               ENDIF
            ENDIF

      ELSEIF (Word(1:1).eq.'$') THEN ! assignment $i=anything
         IF (Word(3:3).ne.'=') THEN
            Fehler = 'PROGRAM ERROR/ $_!='
            ENDIF
         iV = ichar1(Word(2:2)) + 1
         IF (iV.gt.10) THEN
            Report = '\ML-noc'
            RETURN
            ENDIF
         MLVarVal(iV) = Word(4:lenU(Word))

      ELSE
         Fehler = '\ ML command not implemented'
         ENDIF

      IF (Fehler.ne.'&ff') GOTO 900
      RETURN ! command successfully executed

C  Errors :
 900  CONTINUE
      Report = '\ML-err: '//Fehler
      RETURN

      END ! ExeML

      INTEGER FUNCTION iExeMLP (par)
C     ------------------------------
            ! JWu 1jun95 extract integer ML parameter

      CHARACTER  par*(*), ein*60, aus*60

      ein = '\p '//par
      CALL ExeML (ein, aus)
      iExeMLP = 0
      CALL Fi1I (aus, iExeMLP)
cdeb      Print *, ' iExeMLP/ found ', iExeMLP, ' in "', aus(1:lenU(aus)),'"'

      END ! iExeMLP

      REAL*8 FUNCTION rExeMLP (par)
C     -----------------------------
      IMPLICIT NONE
      CHARACTER  par*(*), ein*60, aus*60

      ein = '\p '//par
      CALL ExeML (ein, aus)
      rExeMLP = 0.
      CALL Fi1R (aus, rExeMLP)

      END ! rExeMLP

      SUBROUTINE Lies (questIn, answIO)
C     ---------------------------------
            ! JWu 1989 for NumGra. Oct91 first MetaLanguage elements.
            ! 21jan92 to include question.
         ! ALLE Eingaben sollen ueber diese Routine erfolgen.
         ! Von dieser Routine aus ist die MetaLanguage zugaenglich.

      IMPLICIT LOGICAL (q)

      INCLUDE         'l_dim.f'
      CHARACTER     questIn*(*), answIO*(*),
     *              MLParNam*16, MLParVal*48, MLVarVal*16
      CHARACTER*80  quest, aus, Task(MLT), Done, exereport !Artem: Task -> Task(MLT);
      CHARACTER*83  answ !Artem: CHARACTER*80 -> CHARACTER*83 Increase ANSW size to 83 because ANSW(i:i+3) now falls inside the array.
      CHARACTER     cl6*6

      EXTERNAL      ML_Ini ! initialization of /screen/ and /ML/.

      COMMON /  ML  /  InUnit, qSave, Task, nTask, Done(MLD), !Artem Task(MLT) -> Task
     *                 nDone, MLParNam(MLP), MLParVal(MLP), nMLPar,
     *                 MLVarVal(MLV), qEcho

      la = len(answ)-3 ! now fixed ! Artem: Increase ANSW size to 83 because ANSW(i:i+3) now falls inside the array.

C  Question and answer :
 1    CONTINUE
      quest = questIn
 2    CONTINUE
      answ = ' '

      IF (InUnit.eq.5) THEN ! read from terminal
         ! write question, followed by *no* carriage return
         IF (quest.eq.'\-') THEN
            Write (6,'(a,$)') ' '
         ELSE
            Write (6,'(a,$)') quest(1:lenU(quest)+2)
            ENDIF
         ! get answer
         IF (answIO.ne.'&noq') THEN
            Read (*,'(a)') answ
         ELSE
            Write (6,'(x)') ! just text display, no question
            ENDIF

      ELSEIF (InUnit.eq.0) THEN ! read from stack
         IF (nTask.le.0) THEN
            CALL Gong (6)
            Print *, quest
            Print *, ' PROGRAM ERROR/ Lies/ nTask < 0'
            Print *, ' control returned'
            nTask = 0
            InUnit = 5
            qSave = .true.
            RETURN
            ENDIF
         IF (answIO.ne.'&noq') THEN
            answ = Task(1)
            nTask = nTask - 1
            DO ii = 1, nTask
               Task(ii) = Task(ii+1)
               ENDDO
            ! suppress the comment that contains the old question :
            iq = jPos1 (answ, '%?')
            IF (iq.le.la .and. iq.ge.45) answ(iq:la) = ' '
            ENDIF
         IF (qEcho .and. answ.ne.'\v -') THEN
            ! show answ :
            IF     (answIO.eq.'&noq') THEN ! call just for text display :
               Print '(a)', quest
            ELSEIF (answ.eq.'%') THEN      ! empty line as comment line
               Print *
            ELSEIF (answ(1:1).eq.'%') THEN ! comment line
               Print '(a)', answ
            ELSEIF (answ(1:1).eq.'\') THEN ! Meta-Command: do NOT show question
               IF (lenU(answ).gt.75) THEN
                  CALL Compose2 (aus, ' /'//cl6(nTask+1), '/ ..')
                  Print '(a)', aus
                  Print '(a)', ' .. ', answ
               ELSE
                  CALL Compose2 (aus, ' /'//cl6(nTask+1), '/ '//answ)
                  Print '(a)', aus
                  ENDIF
            ELSE                           ! standard case
               IF (lenU(quest)+lenU(answ).gt.72) THEN
                  CALL Compose3 (aus, quest,' /'//cl6(nTask+1),'/ ..')
                  Print '(a)', aus
                  Print '(a)', ' .. ', answ
               ELSE
                  CALL Compose3 (aus, quest, ' /'//
     *                           cl6(nTask+1), '/ '//answ)
                  Print '(a)', aus
                  ENDIF
               ENDIF
            ENDIF
         IF (nTask.eq.0) THEN ! stack voided
            InUnit = 5
            qSave = .true.
            ENDIF
      ELSE
         CALL Absturz ('Lies', 'InUnit o.o.r.')

         ENDIF

C  Continuation line requested ?
      IF (answ.eq.'..') THEN
         quest = ' .. '
         GOTO 2
         ENDIF

C  Save at the end of the session stack :
      IF (qSave) THEN
         IF (answIO.eq.'&noq') THEN
            ! save as a comment line :
            aus = '%'
            aus(45:80) = '  ! '//questIn
         ELSE
            ! append the question as a comment :
            aus = answ
            IF (lenU(aus).ge.45) THEN ! don't shorten the answer
cUNIX : 2*     CALL Append (aus, ' %? '//questIn)
               CALL Append (aus, ' %? ')
               CALL Append (aus, questIn)
            ELSE
               aus(45:80) = ' %? '//questIn
               ENDIF
            ENDIF

         CALL SaveSessLine (aus)
         ENDIF

      IF (answIO.eq.'&noq') RETURN ! 23apr92  use Lies in Sage for output only

C  Treat MetaLanguage commands :
      IF ((answ(1:1).eq.'\' .and. answ.ne.'\-' .and. answ.ne.'\') .or.
     *    (answ(1:1).eq.'$' .and. answ(3:3).eq.'=')) THEN
         CALL ExeML (answ, exereport)
         IF (exereport.eq.'\ML-noc') GOTO 79 ! was no ML command
         IF (exereport.ne.' ') THEN
            IF (exereport(1:7).eq.'\ML-err') CALL Gong(3)
            IF (exereport(1:7).eq.'\ML-err' .or. qEcho)
     *         Print *, exereport(1:lenU(exereport))
            ENDIF
         GOTO 1       ! ML command executed, now get next input
      ELSEIF (answ(1:1).eq.'%') THEN
         GOTO 1       ! comment line, get next input
         ENDIF
 79   CONTINUE

C  Loop : in-line macros ? (20jan92)
      i = 1
 80   CONTINUE
         IF     (answ(i:i).eq.'%') THEN
            answ(i:la) = ' ' ! rest of line is comment
            GOTO 90
         ELSEIF (answ(i:i+1).eq.'\\' .or.
     *           answ(i:i+1).eq.'\$' .or.
     *           answ(i:i+1).eq.'\%') THEN ! means \,$,%
            CALL DelVonBis (answ, i, i)
            i = i+1
         ELSEIF (answ(i:i).eq.'$') THEN
            iV = ichar1(answ(i+1:i+1)) + 1
            IF (iV.gt.10) GOTO 89 ! no variable => no macro
            CALL DelVonBis (answ, i, i+1)
            la = lenU(MLVarVal(iV))
            CALL Insert (answ, i, MLVarVal(iV)(1:la))
            i = i+la-1
         ELSEIF (answ(i:i+2).eq.'\pi') THEN ! inline macros
            CALL DelVonBis (answ, i, i+2)
            CALL Insert (answ, i, '3.1415926536')
            i = i+12-3
         ELSEIF (answ(i:i+3).eq.'\2pi') THEN
            CALL DelVonBis (answ, i, i+3)
            CALL Insert (answ, i, '6.2831853072')
            i = i+12-3
         ELSEIF (i.ne.1 .and. answ(i:i+1).eq.'\ ') THEN ! break line
            IF (i+2.le.la) THEN
               nTask = nTask + 1
               DO ii = nTask, 2, -1
                  Task(ii) = Task(ii-1)
                  ENDDO
               Task(1) = answ(i+2:la)
               answ(i:la) = ' '
               InUnit = 0
               GOTO 90
            ELSE
               CALL Gong (5)
               Print *,'STRANGE INPUT/ break "\ " at end of input line'
               ENDIF
            ENDIF
 89      CONTINUE
         i = i+1
         IF (i.lt.la) GOTO 80   ! last letter of line is not examined

 90   CONTINUE
      answIO = ' '
      answIO = answ

      END ! Lies

      SUBROUTINE SaveSessLine (aus)
C     -----------------------------
         ! 11mar96 as subroutine. Used by Lies, FehlerGong.

      IMPLICIT LOGICAL (q)

      INCLUDE         'l_dim.f'
      CHARACTER     aus*(*)
      CHARACTER*80  Task(MLT), Done, MLParNam*16, MLParVal*48, !Artem Task -> Task(MLT)
     *              MLVarVal*16

      COMMON /  ML  /  InUnit, qSave, Task, nTask, Done(MLD), !Artem Task(MLT) -> Task
     *                 nDone, MLParNam(MLP), MLParVal(MLP), nMLPar,
     *                 MLVarVal(MLV), qEcho

      IF (nDone.ge.MLD) THEN ! kill them en bloc (16oct93)
         Print '(a)',
     * '\.. session has become very long: deleting first entries'
         nDP = MLD / 4
         nDone = nDone-nDP
         DO ii = 1, nDone
            Done(ii) = Done(ii+nDP)
            ENDDO
         ENDIF
      nDone = nDone + 1
      Done(nDone) = aus

      END ! SaveSessLine

      SUBROUTINE LinEdi (line)
C     ------------------------
         ! 28apr93 : line editor.

      IMPLICIT LOGICAL (q)
      CHARACTER*(*) line
      CHARACTER     com*80, txt*80

         ! the recursive use of FrageC must be eliminated before this
         ! can be called by a ML command \e.

 1    CONTINUE
      Print *, line
      CALL FrageC ('>', com)
      IF     (com.eq.' ') THEN
         RETURN
      ELSEIF (com.eq.'?') THEN
         Print *, 'line editor commands :'
         Print *,
     * 'position the following under characters to be changed :'
         Print *, '   i  insert (at location of i)'
         Print *,
     * '   r  replace (one character by string of any length)'
         Print *, '   d  delete (use rddd to replace four characters)'
      ELSE
         DO i = 77, 1, -1
            IF     (com(i:i).eq.'d') THEN
               CALL DelVonBis (line, i+3, i+3)
            ELSEIF (com(i:i).eq.'i') THEN
               CALL FrageC ('Insert text >', txt)
               CALL Insert (line, i+3, txt(1:lenU(txt)))
            ELSEIF (com(i:i).eq.'r') THEN
               CALL FrageC ('Replace text >', txt)
               CALL DelVonBis (line, i+3, i+3)
               CALL Insert (line, i+3, txt(1:lenU(txt)))
               ENDIF
            ENDDO
         ENDIF
      GOTO 1

      END ! LinEdi

      SUBROUTINE Sage (text)
C     ----------------------
            ! JWu 23apr92
         ! Write text to terminal if echo; save as comment.
      CHARACTER  text*(*)

      CALL Lies (text, '&noq')

      END ! Sage

C  ====================================================================
C  L3.8  Output Schnickschnack
C  ====================================================================

      SUBROUTINE Warte ()
C     -------------------
         ! hold text in scroll area until user allows to continue
         CHARACTER muell *80
         CALL FrageC (' .. to continue, press RETURN', muell)
         IF (muell.ne.' ') CALL Gong (1)
         END ! Warte

      SUBROUTINE Counter (i, j, n, str)
C     ---------------------------------
            ! JWu 12aug91. Text 24feb93
         ! show every j-th step the new counter value i,
         ! the final value is n.

         CHARACTER form*20, cn*1, ch1*1, str*(*)

         IF (i.eq.1) Print *

         ls = min0 (len(str), 60)

         IF (n.gt.0) THEN
            IF ((mod(i,j).eq.0.or.i.eq.n) .and. n.ge.2*j ) THEN
               cn   = ch1(nDigits(n))
               form = '(i'//cn//',a,i'//cn//',1x,a)'
               Print form, i,'/', n, str(1:ls)
               ENDIF
         ELSE ! n unknown
            IF (mod(i,j).eq.0) THEN
               cn   = ch1(nDigits(i)+1)
               form = '(i'//cn//',a,a)'
               Print form, i,'/ ', str(1:ls)
               ENDIF
            ENDIF

         CALL Cursor (1, 0, -1) ! up one line

         END ! Counter

      SUBROUTINE Say2 (t1, t2)
C     ------------------------
         CHARACTER*(*) t1, t2
         CHARACTER*80  aus

         CALL Compose2 (aus, t1, t2)
         Print '(a)', aus

         END ! Say2

      SUBROUTINE Say3 (t1, t2, t3)
C     --------------------------------
         CHARACTER*(*) t1, t2, t3
         CHARACTER*80  aus

         CALL Compose3 (aus, t1, t2, t3)
         Print '(a)', aus

         END ! Say3

      SUBROUTINE Say4 (t1, t2, t3, t4)
C     --------------------------------
         CHARACTER*(*) t1, t2, t3, t4
         CHARACTER*80  aus

         CALL Compose4 (aus, t1, t2, t3, t4)
         Print '(a)', aus

         END ! Say4

      SUBROUTINE Say5 (t1, t2, t3, t4, t5)
C     ------------------------------------
         CHARACTER*(*) t1, t2, t3, t4, t5
         CHARACTER*80  aus

         CALL Compose5 (aus, t1, t2, t3, t4, t5)
         Print '(a)', aus

         END ! Say5

C  ====================================================================
C  L3.9  Control : Absturz
C  ====================================================================

      LOGICAL FUNCTION qErrEntry (wo, Fehler)
C     ---------------------------------------
            ! JWu 3mar93
         ! Check Fehler on entry to subroutine wo
         ! usage : IF (qErrEntry) RETURN
      CHARACTER*(*) Fehler, wo

      IF (Fehler.ne.'&ff') THEN
         Print *, ' old error = ', Fehler(1:lenU(Fehler))
         CALL Gong (1)
         Fehler = 'Error on entry in module '//wo
         qErrEntry = .true.
      ELSE
         qErrEntry = .false.
         ENDIF

      END ! qErrEntry

      SUBROUTINE FehlerGong (Fehler, n)
C     ---------------------------------
            ! ainsi a partir du 22jun94
      CHARACTER Fehler*(*), blabla*80

      IF (Fehler.eq.' ') THEN
         Print *
      ELSE
         CALL Gong (n)
         Print *, Fehler
         CALL ExeML ('\reset % Fehler: '//
     *               Fehler(1:lenU(Fehler)), blabla)
         ENDIF
      Fehler = '&ff'

      END ! FehlerGong

      SUBROUTINE ClearScroll
C     ----------------------
         ! Overlay out, full scroll, goto last line
         CHARACTER     OS*4
         COMMON /Impltn / OS

         IF (OS.eq.'VMS') THEN
            Print '(1x,5a1,$)',Char(29),Char(27),
     *            Char(92),Char(51),Char(24)
            Print '(1x,a1,a)', Char(27), '[0;25r'
            Print '(1x,a1,a)', Char(27), '[25;1H'
         ELSE
            ! nothing to do
            ENDIF

         END ! ClearScroll

      SUBROUTINE Absturz (wo, warum)
C     ------------------------------
         EXTERNAL      RunPreset
         CHARACTER*(*) wo, warum
         CHARACTER     OS*4, blabla*80
         COMMON /Impltn / OS

         CALL ClearScroll ()

C  Error message :
         Print *, 'CRASH in : ', wo
         Print *, 'MESSAGE  : ', warum

C  Save session :
         Print *, 'DUMP  -->  ~/Dump.ml'
         CALL ExeML ('\ws ~/Dump', blabla)

C  Crash :
         IF (OS.eq.'VMS') THEN
            Print *, 'Now provocing another error in order'//
     *              ' to obtain a system traceback :'
            read (wo, '(i3)') wo ! a guaranteed error
         ELSEIF (OS.eq.'X11') THEN
            Print *, ' Hit RETURN to stop'
            Read (*,*)
            STOP
         ELSE
            STOP
            ENDIF

         END ! Absturz
