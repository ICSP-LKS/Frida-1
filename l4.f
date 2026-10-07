C  ====================================================================
C
C     Library  WuL :  General FORTRAN Library
C     Module   L4  :     Get and Decode Input
C
C  ====================================================================

C     J.Wuttke - 1989ff.

C     Outline :
C        This module provides a variety of routines to read
C        input in interactive programs. Every input passes by
C        the routine Lies in module L3 which allows to escape
C        to the MetaLanguage level.
C        All programs using this dialogue system should handle
C        fatal errors by calling Absturz; this routine dumps
C        all previous input before terminating program execution.

C     Contents :
C        4.0.  Interface to L3 :
C                 FragLies
C        4.1.  String input :
C                 FrageC(D),H(D)
C        4.2.  Integer input :
C                 iAsk_
C        4.3.  Real input :
C                 rAsk_
C        4.5.  Boolean input :
C                 qAsk_
C        4.6.  Integer List Input :
C                 FrageNList, GetNList, FrageJList, GetJList
C        4.7.  Integer Multiinput :
C                 i2Frage_
C        4.8.  Real Multiinput :
C                 rAskPair/Trip, rAskRge(Full,Txt), rAskArray/Pairs/Grid/OnOff

C     History :
C        mar93  division L3/L4
C        jul92  formal help
C        1989   first systematic input interface

C  ====================================================================
C  L4.0.  Interface to L3 :
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE FragLies (quest, answ, defin, qhelp)
C     -----------------------------------------------
            ! FrageC 1989, FrageCD 31okt89, 6/21jun91, FragLies 1jul92

         ! Only for internal use by FrageC/H(D).
         ! Modifications of parameters reserved.

         ! Ask for general text.
         ! A default is proposed if it is not '&nodef'
         ! If qhelp then a general help is given.

         IMPLICIT LOGICAL (q)
         CHARACTER*(*)    quest, answ, defin
         CHARACTER        ein*80, aus*80, def*400
         CHARACTER        cl3*3

         lam= len (answ)
         ldm= min0 (lam,len(def))
         lf = lenU(quest)

         qDef = (defin.ne.'&nodef')
         IF (qDef) THEN
            ld = max0 (1, lenU(defin))
            IF (ld.gt.ldm) THEN ! kein Absturz mehr / 13jan00
               def = defin(1:max0(1,ldm-4))//' pp.'
            ELSE
               def = defin ! FORTRAN : defin might occupy the same
                           !           storage as answ, and therefore
                           !           it could unvoluntarily be changed.
               ENDIF
            ENDIF

C  Display the full question :
         IF (.not.qDef) THEN
            aus = quest
         ELSEIF (ichar(def(1:1)).eq.0 .or. def.eq.' ') THEN
            ! there is no default (18may92) :
            aus = quest(1:lf)//' ?'
            qDef = .false.
         ELSEIF (lf+ld+5 .gt. 80) THEN
            ! type question+default in more than one line :
            Print *, quest(1:lf)
            IF (ld.gt.74) THEN
               Print *, '  [', def(1:74)
               DO i = 75, ld-74, 74
                  Print *, '   ', def(i:i+73)
                  ENDDO
               aus = '   '//def(i:ld)//'] ?'
            ELSE
               aus = '  ['//def(1:ld)//'] ?'
               ENDIF
         ELSE
            aus = quest(1:lf)//' ['//def(1:ld)//'] ?'
            ENDIF
         IF (lenU(aus).ge.68) THEN
            ! long question
            Print *, aus
            aus = ' .. '
            ENDIF

C  Loop (1) : try new answer :
  1      CONTINUE
         answ = ' '
         ia = 1   ! first line of answer
C  Loop (2) : continuation of answer :
  2      CONTINUE
            CALL Lies (aus, ein)
            ls = lenU (ein)
            la = len  (answ)

C  Treat continuation demand :
            qCont = .false.
            IF (ls.ge.3) THEN
               IF (ein(ls-1:ls).eq.'..') THEN
                  IF (ein(ls-2:ls-2).eq.char(92)) THEN
                     ! restore .. from \..
                     ein(ls-2:ls) = '..'
                  ELSE
                     ! continue after the present line
                     ls = ls-2
                     qCont = .true.
                     ENDIF
                  ENDIF
               ENDIF
            IF (ia+ls-1.gt.la) THEN
               CALL Gong (3)
               aus = 'maximum string length '//cl3(la)
               GOTO 1
               ENDIF
            IF (ls.ge.1) answ(ia:ia+ls-1) = ein
            IF (qCont) THEN
               ia = ia+ls
               aus = ' .. '
               GOTO 2
               ENDIF
C  End loop (2) : a full answer is read in.

         qConfirm = .false.
         IF     (qhelp .and. answ.eq.'?') THEN
            Print *, 'INPUT HELP/'
            Print *, '   required input : any text'
            IF (qdef) THEN
               CALL Say2 ('   default is "'//def, '"')
               Print *, '   say RETURN to use the default value'
               Print *, '   say "\"    to enter a blank (" ")'
               Print *, '   say "\\"   to enter "\"'
               Print *, '   say "//"   followed by any text '//
     *            'to concatenate with the default'
               Print *, '   say "--//" followed by any text '//
     *            'to delete 2 last letters'
               Print *, '              and concatenate the rest '//
     *            'with the default'
               ENDIF
            Print *, '   use ".." to demand a continuation line'
            Print *, '   try "??" to obtain more specific help, '//
     *                    ' or to answer "?"'
            aus = ' ..'
            GOTO 1
         ELSEIF (qhelp .and. answ(1:2).eq.'??') THEN
            CALL DelVonBis (answ,1,1)
         ELSEIF (qDef .and. answ.eq.' ') THEN
            answ = def
         ELSEIF (answ.eq.char(92)) THEN ! to replace RETURN for answering ' '
            answ = ' '
            Print *, ' understood : '' '''
         ELSEIF (answ.eq.'\\') THEN
            answ =  char(92)
            Print *, ' understood : ''', char(92), ''''
         ELSEIF (qDef .and. answ(1:2).eq.'//') THEN ! concatenation with default
            CALL DelVonBis (answ, 1, 2)
            CALL Insert (answ, 1, def(1:ld))
            Print *, ' understood : ', answ(1:lenU(answ))
         ELSEIF (qDef .and. answ(1:1).eq.'-') THEN ! delete part of default ?
            DO i = 1, len(answ)
               IF (answ(i:i).ne.'-') GOTO 389
               ENDDO
 389        CONTINUE
            IF (answ(i:i+1).eq.'//') THEN ! delete def(1:i-1), append answ(i+2:)
               CALL DelVonBis (answ, 1, i+1)
               IF (i-1.lt.ld) CALL Insert (answ, 1, def(1:ld-(i-1))) !corr aug94
               Print *, ' understood : ', answ(1:lenU(answ))
               ENDIF
            ENDIF ! special cases treated

         IF (lenU(answ).gt.lam) THEN
            CALL Gong (3)
            aus = 'maximum string length '//cl3(lam)
            GOTO 1
            ENDIF
C  End loop (1) : the answer is valid.

         END ! FragLies

C  ====================================================================
C  L4.1.  String input :
C  ====================================================================

      SUBROUTINE FrageC (quest, answ)
C     -------------------------------
      CHARACTER*(*) quest, answ
      CALL FragLies (quest, answ, '&nodef', .true.)
      END ! FrageC

      SUBROUTINE FrageCD (quest, answ, defin)
C     ---------------------------------------
      CHARACTER*(*) quest, answ, defin
      CALL FragLies (quest, answ, defin, .true.)
      END ! FrageCD

      SUBROUTINE FrageT (quest, answ)
C     -------------------------------
      CHARACTER*(*) quest, answ
 1    CALL FragLies (quest, answ, '&nodef', .true.)
      IF (answ.eq.' ') THEN
         CALL Gong(1)
         Print *, ' answer must be nonblank - please repeat input :'
         GOTO 1
         ENDIF
      END ! FrageT

      SUBROUTINE FrageTD (quest, answ, defin)
C     --------------------------------------- ! 9nov96
      CHARACTER*(*) quest, answ, defin
 1    CALL FragLies (quest, answ, defin, .true.)
      IF (answ.eq.' ') THEN
         CALL Gong(1)
         Print *, ' answer must be nonblank - please repeat input :'
         GOTO 1
         ENDIF
      END ! FrageTD

      SUBROUTINE FrageH (quest, answ)
C     -------------------------------
      CHARACTER*(*) quest, answ
      CALL FragLies (quest, answ, '&nodef', .false.)
      END ! FrageH

      SUBROUTINE FrageHD (quest, answ, defin)
C     ---------------------------------------
      CHARACTER*(*) quest, answ, defin
      CALL FragLies (quest, answ, defin, .false.)
      END ! FrageHD

C  ====================================================================
C  L4.2.  Integer Input :
C  ====================================================================

C     Einzelne Funktionen iAsk... :
C     -----------------------------

      INTEGER FUNCTION iAsk (cFra )
         CHARACTER  *(*) cFra
         iAsk = iAskGen (cFra,0,.false.,.false.,0,0,.false.,.false.)
         END

      INTEGER FUNCTION iAskD (cFra, idef)
         CHARACTER  *(*) cFra
         iAskD = iAskGen (cFra,idef,.true.,.true.,0,0,.false.,.false.)
         END

      INTEGER FUNCTION iAskDu (cFra, idef)
         CHARACTER  *(*) cFra
         iAskDu = iAskGen(cFra,idef,.true.,.false.,0,0,.false.,.false.)
         END

      INTEGER FUNCTION iAskM (cFra, imin, imax)
         CHARACTER  *(*) cFra
         iAskM = iAskGen (cFra,0,.false.,.false.,
     *                    imin,imax,.true.,.true.)
         END

      INTEGER FUNCTION iAskMu (cFra, imin, imax)
         CHARACTER  *(*) cFra
         iAskMu = iAskGen (cFra,0,.false.,.false.,
     *                     imin,imax,.true.,.false.)
         END
      INTEGER FUNCTION iAskDM (cFra, idef, imin, imax)
         CHARACTER  *(*) cFra
         iAskDM = iAskGen (cFra,idef,.true.,.true.,
     *                     imin,imax,.true.,.true.)
         END

      INTEGER FUNCTION iAskDMu (cFra, idef, imin, imax)
         CHARACTER  *(*) cFra
         iAskDMu = iAskGen (cFra,idef,.true.,.true.,
     *                      imin,imax,.true.,.false.)
         END

      INTEGER FUNCTION iAskDuM (cFra, idef, imin, imax)
         CHARACTER  *(*) cFra
         iAskDuM = iAskGen (cFra,idef,.true.,.false.,
     *                      imin,imax,.true.,.true.)
         END

      INTEGER FUNCTION iAskDuMu (cFra, idef, imin, imax)
         CHARACTER  *(*) cFra
         iAskDuMu = iAskGen (cFra,idef,.true.,.false.,
     *                       imin,imax,.true.,.false.)
         END

      INTEGER FUNCTION iAskGen (cFra, idefin, qdefin, qdefdisin,
     *                            imin, imax, qmm, qmmdis )
C     -----------------------------------------------------
            ! JWu ca.1989. LiesIC incorporated jan92.
         !  Mutterfunktion fuer iAskM/D/DM/...
         !  Falls qdef/qmm = on, steht Default zur Verfuegung und
         !  wird auf Einhaltung der Grenzen [imin,imax] getestet;
         !  falls qdefdis/qmmdis = on, werden idef/imin,imax in
         !  die Anzeige cFra aufgenommen.

         !  ' '    means default
         !  '=+#'  means default plus number #    (22jan92)

      IMPLICIT LOGICAL (q)
      PARAMETER        (MI=6)
      CHARACTER         cFra*(*), aus*80, ein*80, cl8*8
      INTEGER           iIn(MI)

C  Check parameters :
      qdef = qdefin
      idef = idefin
      IF (qdef) THEN
         IF (idef.gt.(1.E19)) THEN
            aus = 'idef > 1E19, also wahrscheinlich INDEFINITE'
            CALL Absturz ('iAskGen',aus)
            ENDIF
         IF (qmm .and. (idef.lt.imin .or. idef.gt.imax)) THEN ! 25aug94
            qdef = .false.
            ENDIF
         ENDIF
      qdefdis = qdef .and. qdefdisin

      IF (qmm) THEN
         IF (imax.lt.imin) THEN
            aus = 'imax < imin'
            CALL Absturz ('iAskGen',aus)
            ENDIF
         IF (imax.gt.(1.E19)) THEN
            aus = 'imax > 1E19, also wahrscheinlich INDEFINITE'
            CALL Absturz ('iAskGen',aus)
            ENDIF
         ENDIF

C  Construct question line :
 3    CONTINUE
      aus = cFra
      IF     (     qdefdis.and.     qmmdis)  THEN
         CALL Append (aus, ' [ ' // cl8(imin) )
         CALL Append (aus, '..')
         CALL Append (aus, cl8(imax) )
         CALL Append (aus, ' ; def=')
         CALL Append (aus, cl8(idef) )
         CALL Append (aus, ' ] ?')
      ELSEIF (     qdefdis.and..not.qmmdis)  THEN
         CALL Append (aus, ' ['//cl8(idef) )
         CALL Append (aus, '] ?')
      ELSEIF (.not.qdefdis.and.     qmmdis)  THEN
         CALL Append (aus, ' [' // cl8(imin) )
         CALL Append (aus, '..' // cl8(imax) )
         CALL Append (aus, '] ?')
      ELSE
         ! do *not* append a '?' - may or may not be contained in cFra
         ENDIF

C  Question and answer :
 10   CONTINUE
      CALL Lies (aus, ein)
      CALL FindI (ein, MI, ni, iIn)
      CALL Kontrakt (ein)

      IF     (ein.eq.'?') THEN
         Print *, 'INPUT HELP/'
         Print *, '   required input : an integer value'
         IF (qmm) THEN
            Print *, '   minimal value = '//cl8(imin)
            Print *, '   maximal value = '//cl8(imax)
            ENDIF
         IF (qdef) THEN
            Print *, '   default value = '//cl8(idef)
            Print *, '      answer RETURN to obtain the default value'
            Print *, '      answer "=+5", "=-1", ... to modify it'
            ENDIF
         GOTO 10
      ELSEIF (ein.eq.' ' .and. qdef) THEN
         i = idef
      ELSEIF (ein.eq.'#') THEN
         IF (ni.ne.1)
     * CALL Absturz ('iAskGen', 'ni=1 expected, found '//cl8(ni))
         i = iIn(1)
      ELSEIF (ein.eq.'=#' .and. qdef) THEN
         IF (ni.ne.1)
     * CALL Absturz ('iAskGen', 'ni=1 expected, found '//cl8(ni))
         i = idef + iIn(1)
         Print *, ' understood : '//cl8(i)
      ELSEIF (ein.eq.char(27)//'[A' .and. qdef) THEN ! arrow upwards
         idef = idef - 1
         GOTO 3
      ELSEIF (ein.eq.char(27)//'[B' .and. qdef) THEN ! arrow dnwards
         idef = idef + 1
         GOTO 3
      ELSE
         CALL Gong (5)
         GOTO 10
         ENDIF

      IF (qmm  .and. (i.lt.imin .or. i.gt.imax)) THEN
         CALL Gong (3)
         aus = 'BAD INPUT/ allowed range is '//cl8(imin)
         CALL Append (aus, '...'//cl8(imax))
         CALL Append (aus, '.')
         Print *, aus
         GOTO 3
         ENDIF

      iAskGen = i

      END ! iAskGen

C  ====================================================================
C  L4.3.  Real Input :
C  ====================================================================

      REAL*8 FUNCTION rAskGen (cFra, def, rmi, rma, jd, jm, jl)
C     ----------------------------------------------------------
         ! only for internal use,
         ! to be called by the specific functions rAsk.. below.
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER         cFra*(*), cd*16, cmi*16, cma*16,
     *                  ein*80, aus*80, c*1
      ! jd = 0 : no default
      !    = 1 : hidden default
      !    = 2 : show default

      ! jm =...: rmi =< rAskGen =< rma
      ! jl =...: rmi <  rAskGen <  rma

         jjd = jd
         jjm = jm
         jjl = jl

         IF (jjm.ge.1) THEN ! check of consistency
            IF (rmi.gt.rma) CALL Absturz ('rAskGen', 'min >= max')
            IF (jjd.ge.1 .and. (def.lt.rmi .or. def.gt.rma)) jjd = 0
               ! default out of range (allowed since 5jul91) : proceed without
            ENDIF
         IF (jjl.ge.1) THEN ! check of consistency
            IF (rmi.ge.rma) CALL Absturz ('rAskGen', 'min >= max')
            IF (jjd.ge.1 .and. (def.le.rmi .or. def.ge.rma)) jjd = 0
            ENDIF

         CALL NiceNum (def, cd,  id)
         CALL NiceNum (rmi, cmi, imi)
         CALL NiceNum (rma, cma, ima)

 3       aus = cFra
         IF     (jjd.ge.2 .and. jjm.ge.2) THEN
            CALL Append (aus, ' ['//cmi(1:imi)//'..'//cma(1:ima)//
     *                        ' ; def='//cd(1:id)//'] ?')
         ELSEIF (jjd.ge.2 .and. jjl.ge.2) THEN
            CALL Append (aus, ' ['//cmi(1:imi)//'<..<'//cma(1:ima)//
     *                        ' ; def='//cd(1:id)//'] ?')
         ELSEIF (jjd.ge.2) THEN
            CALL Append (aus, ' ['//cd(1:id)//'] ?')
         ELSEIF (jjd.lt.2 .and. jjm.ge.2) THEN
            CALL Append (aus,' ['//cmi(1:imi)//'..'//cma(1:ima)//'] ?')
         ELSEIF (jjd.lt.2 .and. jjl.ge.2) THEN
            CALL Append (aus, ' ['//cmi(1:imi)//'<..<'//cma(1:ima)//
     *                        '] ?')
         ELSE
            ! do NOT append (aus, ' ?')
            ENDIF

 2       CONTINUE
         CALL Lies (aus, ein)
         CALL DelLeft (ein)
         CALL Fi1R (ein, xneu)

         IF     (ein.eq.'?') THEN
            Print *, 'INPUT HELP/'
            Print *, '   required input : a real value'
            IF     (jjm.ge.1) THEN
               Print *, '   minimal value = '//cmi
               Print *, '   maximal value = '//cma
            ELSEIF (jjl.ge.1) THEN
               Print *, '   value must be > '//cmi
               Print *, '             and < '//cma
               ENDIF
            IF (jjd.ge.1) THEN
               Print *, '   default value = '//cd
               Print *,
     *          '      answer RETURN to obtain the default value'
               ENDIF
            Print *, '   format samples :'
            Print *, '      1, -1., .3, 5e-6, -90.6E2, ...'
            GOTO 2
         ELSEIF (ein.eq.' ' .and. jjd.ge.1) THEN
            rAskGen = def
            RETURN
         ELSEIF (ein.eq.'#') THEN
            ! value given
         ELSE
            ! Texteingabe nicht vorgesehen
            CALL Gong (3)
            GOTO 2
            ENDIF

         IF (jjm.ge.1 .and. (xneu.lt.rmi .or. xneu.gt.rma)) THEN
            CALL Gong (1)
            jjm = 2
            GOTO 3
            ENDIF
         IF (jjl.ge.1 .and. (xneu.le.rmi .or. xneu.ge.rma)) THEN
            CALL Gong (1)
            jjl = 2
            GOTO 3
            ENDIF
         rAskGen = xneu
         RETURN

         END ! rAskGen

      REAL*8 FUNCTION rAsk (cFra)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAsk = rAskGen (cFra, 0.d0, 0.d0, 0.d0, 0, 0, 0)
         END

      REAL*8 FUNCTION rAskD (cFra, def)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskD = rAskGen (cFra, def, 0.d0, 0.d0, 2, 0, 0)
         END

      REAL*8 FUNCTION rAskDu (cFra, def)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDu = rAskGen (cFra, def, 0.d0, 0.d0, 1 , 0, 0)
         END

      REAL*8 FUNCTION rAskM (cFra, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskM = rAskGen (cFra, 0.d0, rmin, rmax, 0, 2, 0)
         END

      REAL*8 FUNCTION rAskMu (cFra, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskMu = rAskGen (cFra, 0.d0, rmin, rmax, 0, 1, 0)
         END

      REAL*8 FUNCTION rAskL (cFra, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskL = rAskGen (cFra, 0.d0, rmin, rmax, 0, 0, 2)
         END

      REAL*8 FUNCTION rAskLu (cFra, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskLu = rAskGen (cFra, 0.d0, rmin, rmax, 0, 0, 1)
         END

      REAL*8 FUNCTION rAskDM (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDM = rAskGen (cFra, def, rmin, rmax, 2, 2, 0)
         END

      REAL*8 FUNCTION rAskDMu (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDMu = rAskGen (cFra, def, rmin, rmax, 2, 1, 0)
         END

      REAL*8 FUNCTION rAskDuM (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDuM = rAskGen (cFra, def, rmin, rmax, 1, 2, 0)
         END

      REAL*8 FUNCTION rAskDuMu (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDuMu = rAskGen (cFra, def, rmin, rmax, 1, 1, 0)
         END

      REAL*8 FUNCTION rAskDL (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDL = rAskGen (cFra, def, rmin, rmax, 2, 0, 2)
         END

      REAL*8 FUNCTION rAskDLu (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDLu = rAskGen (cFra, def, rmin, rmax, 2, 0, 1)
         END

      REAL*8 FUNCTION rAskDuL (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDuL = rAskGen (cFra, def, rmin, rmax, 1, 0, 2)
         END

      REAL*8 FUNCTION rAskDuLu (cFra, def, rmin, rmax)
         IMPLICIT REAL*8 (a-h,o-p,r-z)
         CHARACTER  *(*) cFra
         rAskDuLu = rAskGen (cFra, def, rmin, rmax, 1, 0, 1)
         END

C  ====================================================================
C  L4.5. Logical Input :
C  ====================================================================

      LOGICAL FUNCTION qAsk (cFra)
         IMPLICIT LOGICAL (q)
         CHARACTER  *(*) cFra
         qAsk = qAskGen (cFra, 0, 0)
         END

      LOGICAL FUNCTION qAskD (cFra, idef)
         IMPLICIT LOGICAL (q)
         CHARACTER  *(*) cFra
         qAskD = qAskGen (cFra, idef, 2)
         END

      LOGICAL FUNCTION qAskDi (cFra, idef)
         IMPLICIT LOGICAL (q)
         CHARACTER  *(*) cFra
         qAskDi = qAskGen (cFra, idef, 3)
         END

      LOGICAL FUNCTION qAskDu (cFra, idef)
         IMPLICIT LOGICAL (q)
         CHARACTER  *(*) cFra
         qAskDu = qAskGen (cFra, idef, 1)
         END

      LOGICAL FUNCTION qAskGen (cFra, idef, jd )
C     -------------------------------------------
         IMPLICIT LOGICAL (q)
         CHARACTER cFra*(*), ein*80, aus*80, ch1*1
      ! jd   = 0 : no default
      !      = 1 : hidden default
      !      = 2 : show default
      !      = 3 : show default as integer value

      ! idef = 0 : off=false=0=n=N
      !      = 1 : on =true =1=y=Y=j=J=o=O=s=S

 3       aus = cFra
         IF     (jd.eq.2 .and. qintr(idef)) THEN
            CALL Append (aus, ' [y] ?')
         ELSEIF (jd.eq.2) THEN
            CALL Append (aus, ' [n] ?')
         ELSEIF (jd.eq.3) THEN
            CALL Append (aus, ' ['//ch1(idef)//'] ?')
         ELSE
            ! CALL Append (aus, ' ?')
            ENDIF
 2       CALL Lies (aus, ein)

         IF (ein.eq.'?') THEN
            Print *, 'INPUT HELP/'
            Print *, '   required input : a logical value'
            Print *, '   format :  1,+,y,Y,j,J  or  0,-,n,N '
            IF (jd.ge.1) THEN
               Print *, '   default value = '//ch1(idef)
               Print *,
     *   '      answer RETURN to obtain the default value'
               ENDIF
            GOTO 2
            ENDIF
         IF (ein.eq.' ') THEN
            IF (jd.ge.1) THEN
               qAskGen = qintr (idef)
               RETURN
            ELSE
               CALL Gong (5)
               GOTO 2
               ENDIF
            ENDIF
         CALL Kontrakt (ein)
         IF (lenU(ein).ne.1) THEN
            CALL Gong (3)
            GOTO 2
            ENDIF
         IF     (jPos1('1yYjJ+',ein(1:1)).le.6) THEN
            qAskGen = .true.
         ELSEIF (jPos1('0nN-',ein(1:1)).le.4) THEN
            qAskGen = .false.
         ELSE
            CALL Gong (3)
            GOTO 2
            ENDIF

         END ! qAskGen

C  ====================================================================
C  L4.6.  Integer Lists :
C  ====================================================================

      SUBROUTINE FrageNList (quest, answer, nL)
C     -----------------------------------------
            ! JWu 7feb91. Ein Sorgenkind. 21jun91.
         ! ask for input in list format, help facility, no checks
         ! if answer is not empty (' '), it is also proposed as default

         CHARACTER  *(*) quest, answer
         CHARACTER       seq*100, cl6*6

         seq = quest

 1       CONTINUE
         CALL FrageHD (seq, answer, answer)

         IF     (answer.eq.' ') THEN
            CALL Gong (1)
            Print *, ' for an empty list, answer ''-'''
            GOTO 1
         ELSEIF (answer.eq.'?') THEN
            Print *, 'INPUT HELP/'
            Print *, '   required input : '//
     *               ' an list of integers in ascending order'
            Print *, '   allowed range  : from 1 to '//cl6(nL)
            Print *, '   format samples :'
            Print *, '      17'
            Print *, '      1,3,5'
            Print *, '      1-5i2          means  1,3,5'
            Print *, '      /17-19         means  all but 17,18,19'
            Print *, '      *              means  whole range'
            Print *, '      6-*i3          means  6,9,12,..'
            Print *, '      -              means  empty list'
            GOTO 1
            ENDIF

         END ! FrageNList

      SUBROUTINE GetNList (quest, answer, qList, nL)
C     ----------------------------------------------
         ! JWu 18apr91, default answer 4jul91.
      ! ask for input in list format and decode it.
      ! brings together the subroutines FrageNList and DecNList.

      CHARACTER quest*(*), answer*(*), Fehler*80
      LOGICAL   qList(nL)

C  Check parameter nL (fast ein Fall fuer Absturz) :
      IF (nL.le.0) THEN
         CALL Gong (12)
         Print *, ' WARNING/ WuL4/ GetNList/ n<1, returning empty list'
         answer = '-'
         RETURN
         ENDIF

C  Check default answer :
      Fehler = '&ff'
      CALL DecNList (answer, qList, nL, Fehler)
      IF (Fehler.ne.'&ff') answer = ' ' ! dans ce cas, pas de default

C  Get new answer and decode it :
 1    CONTINUE
      Fehler = '&ff'
      CALL FrageNList (quest, answer, nL)
      CALL DecNList (answer, qList, nL, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL FehlerGong (Fehler,1)
         answer = ' ' ! plus de default
         GOTO 1
         ENDIF

      END ! GetNList

      SUBROUTINE FrageJList (quest, answer, ji, jf, MJ)
C     -------------------------------------------------
            ! JWu 9jul93 copied from FrageNList.
         ! ask for input in j-list format, help facility, no checks
         ! if answer is not empty (' '), it is also proposed as default

      CHARACTER*(*) quest, answer
      CHARACTER     cl6*6

 1    CONTINUE
      CALL FrageHD (quest, answer, answer)

      IF     (answer.eq.' ') THEN
         CALL Gong (1)
         Print *, ' for an empty list, answer ''-'''
         GOTO 1
      ELSEIF (answer.eq.'?') THEN
         Print *, 'INPUT HELP/'
         Print *, '   required input : a list of integers'
         CALL Say2 ('   allowed range  : from '//
     *              cl6(ji), ' to '//cl6(jf))
         CALL Say2 ('   not more than '//cl6(MJ), ' elements')
         Print *, '   format samples :'
         Print *, '      17'
         Print *, '      1,5,2,6'
         Print *, '      1-5i2          means  1,3,5'
         Print *, '      /17-19         means  all but 17,18,19'
         Print *, '      *              means  whole range'
         Print *, '      6-*i3          means  6,9,12,..'
         Print *, '      -              means  empty list'
         GOTO 1
         ENDIF

      END ! FrageJList

      SUBROUTINE GetJList (ques, answ, MJL, nJL, JL, ji, jf)
C     ------------------------------------------------------
         ! JWu 9jul93
      ! ask for input in j-list format and decode it.
      ! brings together the subroutines FrageJList and DecJList.
         ! import : question text ques,
         !          limits : <= MJL elements ranging from ji to jf.
         ! export : literal answer answ,
         !          # elements nJL, elements JL(1..nJL)

      IMPLICIT NONE
      INTEGER   ji, jf, MJL, nJL, JL(*), lenU
      CHARACTER ques*(*), answ*(*), Fehler*80

C  Check parameter :
      IF     (MJL.le.0) THEN
         CALL Absturz ('GetJList', 'MJL<=0')
      ELSEIF (ji.gt.jf) THEN
         CALL Absturz ('GetJList', 'ji>jf')
         ENDIF

C  Check default answer :
      Fehler = '&ff'
      CALL DecJList (answ, MJL, nJL, JL, ji, jf, Fehler)
      IF (Fehler.ne.'&ff') answ = ' ' ! dans ce cas, pas de default

C  Get new answer and decode it :
 1    CONTINUE
      Fehler = '&ff'
      CALL FrageJList (ques, answ, ji, jf, MJL)
      CALL DecJList (answ, MJL, nJL, JL, ji, jf, Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL FehlerGong (Fehler,1)
         answ = ' ' ! plus de default
         GOTO 1
         ENDIF

      END ! GetJList

C  ====================================================================
C  L4.7.  Integer Multiinput :
C  ====================================================================

      SUBROUTINE immFrageD (cFra, imin, imax, imiD, imaD)
C     ---------------------------------------------------
         CHARACTER  cFra*(*), cl8*8, ein*80, aus*80, aus2*80, c*1
         INTEGER    imin, imax, imiD, imaD, iSta(3)

         aus = cFra

         CALL Append (aus, ' ['//cl8(imiD) )
         CALL Append (aus, ','//cl8(imaD) )
         CALL Append (aus, '] ?')

 2       CALL FrageH (aus, ein)
         CALL FindI (ein, 3, ni, iSta)
         CALL Kontrakt (ein)

         IF     (ni.eq.0 .and. ein.eq.' ') THEN
            imin = imiD
            imax = imaD
         ELSEIF (ni.eq.1 .and. ein.eq.'#,') THEN
            imin = iSta(1)
            imax = imaD
         ELSEIF (ni.eq.1 .and. ein.eq.',#') THEN
            imin = imiD
            imax = iSta(1)
         ELSEIF (ni.eq.2 .and.(ein.eq.'#,#' .or. ein.eq.'##'))THEN
            imin = iSta(1)
            imax = iSta(2)
         ELSE
            Print *, ' BAD FORMAT/'
            Print *, ' answer "x, y" or "x y" or "x, " or " ,y" or " "'
            CALL Gong(1)
            GOTO 2
            ENDIF

         IF (imin.ge.imax) THEN
            CALL Gong (4)
            Print *, ' BAD INPUT/ please let min < max'
            GOTO 2
            ENDIF

         END ! immFrageD

      SUBROUTINE i2FrageD (cFra, imin, imax, imiD, imaD)
C     --------------------------------------------------
         CHARACTER  cFra*(*), cl8*8, ein*80, aus*80
         INTEGER    imin, imax, imiD, imaD, iSta(3)

         aus = cFra

         CALL Append (aus, ' [ '//cl8(imiD) )
         CALL Append (aus, ' , '//cl8(imaD) )
         CALL Append (aus, ' ] ?')

 2       CALL FrageH (aus, ein)
         CALL FindI (ein, 3, ni, iSta)
         CALL Kontrakt (ein)

         IF     (ni.eq.0 .and. ein.eq.' ') THEN
            imin = imiD
            imax = imaD
         ELSEIF (ni.eq.1 .and. ein.eq.'#,') THEN
            imin = iSta(1)
            imax = imaD
         ELSEIF (ni.eq.1 .and. ein.eq.',#') THEN
            imin = imiD
            imax = iSta(1)
         ELSEIF (ni.eq.2 .and.(ein.eq.'#,#' .or. ein.eq.'##'))THEN
            imin = iSta(1)
            imax = iSta(2)
         ELSE
            Print *, ' BAD FORMAT/'
            Print *, ' answer "x, y" or "x y" or "x, " or " ,y" or " "'
            CALL Gong(1)
            GOTO 2
            ENDIF

      END ! i2FrageD

      SUBROUTINE i2Frage (cFra, imin, imax)
C     -------------------------------------
         CHARACTER  cFra*(*), cl8*8, ein*80, aus*80
         INTEGER    imin, imax, iSta(3)

         aus = cFra

 2       CALL FrageH (aus, ein)
         CALL FindI (ein, 3, ni, iSta)
         CALL Kontrakt (ein)

         IF (ni.eq.2 .and.(ein.eq.'#,#' .or. ein.eq.'##'))THEN
            imin = iSta(1)
            imax = iSta(2)
         ELSE
            Print *, ' BAD FORMAT/'
            Print *, ' answer "x, y" or "x y" or "x, " or " ,y" or " "'
            CALL Gong(1)
            GOTO 2
            ENDIF

         END ! i2Frage

      SUBROUTINE i3FrageD (cFra, i1, i2, i3, i1d, i2d, i3d)
C     -----------------------------------------------------
            ! JWu 14sep94
         IMPLICIT NONE
         CHARACTER  cFra*(*), cl8*8, ein*80, aus*80
         INTEGER    i1, i2, i3, i1d, i2d, i3d, iSta(4), ni

         CALL Compose5 (aus, cFra,
     *        ' ['//cl8(i1d), ','//cl8(i2d), ','//cl8(i3d), '] ?')

 2       CALL FrageH (aus, ein)
         CALL FindI (ein, 4, ni, iSta)
         CALL Kontrakt (ein)

         i1 = i1d
         i2 = i2d
         i3 = i3d
         IF     (ni.eq.0 .and. ein.eq.' ') THEN
            ! defaults accepted
         ELSEIF (ni.eq.1 .and. ein.eq.'#,,') THEN
            i1 = iSta(1)
         ELSEIF (ni.eq.1 .and. ein.eq.',#,') THEN
            i2 = iSta(1)
         ELSEIF (ni.eq.1 .and. ein.eq.',,#') THEN
            i3 = iSta(1)
         ELSEIF (ni.eq.2 .and.(ein.eq.'##,' .or. ein.eq.'#,#,'))THEN
            i1 = iSta(1)
            i2 = iSta(2)
         ELSEIF (ni.eq.2 .and. ein.eq.'#,,#') THEN
            i1 = iSta(1)
            i3 = iSta(2)
         ELSEIF (ni.eq.2 .and.(ein.eq.',##' .or. ein.eq.',#,#'))THEN
            i2 = iSta(1)
            i3 = iSta(2)
         ELSEIF (ni.eq.3 .and.(ein.eq.'#,#,#' .or. ein.eq.'###'))THEN
            i1 = iSta(1)
            i2 = iSta(2)
            i3 = iSta(3)
         ELSEIF (ein.eq.'?') THEN
            Print *, ' INPUT HELP:'
            Print *,
     * ' answer "x,y,z" or leave out a number to accept the default'
         ELSE
            CALL Gong(1)
            GOTO 2
            ENDIF

      END ! i3FrageD

C  ====================================================================
C  L4.8.  Real Multiinput :
C  ====================================================================

      SUBROUTINE rAskPair (cFra, r1, r2, r1D, r2D)
C     --------------------------------------------!WARNING: rAskPairs<>rAskPair
            ! Restored from rAskRge 17nov93
         ! Ask for a real pair r1, r2.
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER        cFra*(*), c1*16, c2*16, ein*80, aus*80
      DIMENSION        rSta(3)

C  Build question :
         CALL NiceNum (r1D, c1, i1)
         CALL NiceNum (r2D, c2, i2)
         CALL Compose2 (aus,cFra,' ['//c1(1:i1)//' '//c2(1:i2)//'] ?')

C  Input and decoding :
 2       CALL FrageH (aus, ein)
         CALL FindR (ein, 3, nr, rSta, .false.)
         CALL Kontrakt (ein)

C  Default or new values ?
         r1 = r1D
         r2 = r2D
         IF     (nr.eq.0 .and. ein.eq.' ') THEN
            ! both defaults accepted
         ELSEIF (nr.eq.1 .and. ein.eq.'#,') THEN
            r1 = rSta(1)
         ELSEIF (nr.eq.1 .and. ein.eq.',#') THEN
            r2 = rSta(1)
         ELSEIF (nr.eq.2 .and.(ein.eq.'#,#' .or. ein.eq.'##'))THEN
            r1 = rSta(1)
            r2 = rSta(2)
         ELSEIF (ein.eq.'?') THEN
            Print *, ' INPUT HELP:'
            Print *, ' answer "x, y" or "x y" or "x, " or " ,y" or " "'
         ELSE
            CALL Gong(1)
            GOTO 2
            ENDIF
         END ! rAskPair

      SUBROUTINE rAskTrip (cFra, r1, r2, r3, r1D, r2D, r3D)
C     -----------------------------------------------------
            ! Copied from rAskPair 25jan94.
         ! Ask for a triplet of reals.
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER        cFra*(*), c1*16, c2*16, c3*16, ein*80, aus*80
      DIMENSION        rSta(4)

C  Build question :
         CALL NiceNum (r1D, c1, i1)
         CALL NiceNum (r2D, c2, i2)
         CALL NiceNum (r3D, c3, i3)
         CALL Compose2 (aus, cFra,
     *      ' ['//c1(1:i1)//' '//c2(1:i2)//' '//c3(1:i3)//'] ?')

C  Input and decoding :
 2       CALL FrageH (aus, ein)
         CALL FindR (ein, 4, nr, rSta, .false.)
         CALL Kontrakt (ein)

C  Default or new values ?
         r1 = r1D
         r2 = r2D
         r3 = r3D
         IF     (nr.eq.0 .and. ein.eq.' ') THEN
            ! all defaults accepted
         ELSEIF (nr.eq.1 .and. ein.eq.'#,,') THEN
            r1 = rSta(1)
         ELSEIF (nr.eq.1 .and. ein.eq.',#,') THEN
            r2 = rSta(1)
         ELSEIF (nr.eq.1 .and. ein.eq.',,#') THEN
            r3 = rSta(1)
         ELSEIF (nr.eq.2 .and.(ein.eq.'##,' .or. ein.eq.'#,#,'))THEN
            r1 = rSta(1)
            r2 = rSta(2)
         ELSEIF (nr.eq.2 .and. ein.eq.'#,,#') THEN
            r1 = rSta(1)
            r3 = rSta(2)
         ELSEIF (nr.eq.2 .and.(ein.eq.',##' .or. ein.eq.',#,#'))THEN
            r2 = rSta(1)
            r3 = rSta(2)
         ELSEIF (nr.eq.3 .and.(ein.eq.'#,#,#' .or. ein.eq.'###'))THEN
            r1 = rSta(1)
            r2 = rSta(2)
            r3 = rSta(3)
         ELSEIF (ein.eq.'?') THEN
            Print *, ' INPUT HELP:'
            Print *,
     * ' answer "x,y,z" or leave out a number to accept the default'
         ELSE
            CALL Gong(1)
            GOTO 2
            ENDIF
         END ! rAskTrip

      SUBROUTINE rAskRge (cFra, r1, r2, r1D, r2D)
C     -------------------------------------------
            ! old r2/mmFrageD. Renamed 3mar93.
         ! Ask for a real range r1<r2.
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER        cFra*(*), c1*16, c2*16, ein*80, aus*80
      DIMENSION        rSta(3)

C  Check default :
      IF (rmiD.gt.rmaD) THEN  !CALL Absturz('rmmFrageD','Def:mi>ma')
         CALL Gong (8)
         Print *, ' BAD USE of WuL/ rAskRge/ bad defaults given :'
         Print *, ' min, max = ', rmiD, rmaD
         rmiD = 0.d0
         rmaD = 0.d0
         ENDIF

C  Build question :
         CALL NiceNum (r1D, c1, i1)
         CALL NiceNum (r2D, c2, i2)
         CALL Compose2 (aus, cFra,
     *      ' ['//c1(1:i1)//' '//c2(1:i2)//'] ?')

C  Input and decoding :
 2       CALL FrageH (aus, ein)
         CALL FindR (ein, 3, nr, rSta, .false.)
         CALL Kontrakt (ein)

C  Default or new values ?
         r1 = r1D
         r2 = r2D
         IF     (nr.eq.0 .and. ein.eq.' ') THEN
            ! both defaults accepted
         ELSEIF (nr.eq.1 .and. ein.eq.'#,') THEN
            r1 = rSta(1)
         ELSEIF (nr.eq.1 .and. ein.eq.',#') THEN
            r2 = rSta(1)
         ELSEIF (nr.eq.2 .and.(ein.eq.'#,#' .or.
     *                         ein.eq.'##'))THEN
            r1 = rSta(1)
            r2 = rSta(2)
         ELSEIF (ein.eq.'?') THEN
            Print *, ' INPUT HELP:'
            Print *, ' answer "x, y" or "x y" or "x, " or " ,y" or " "'
         ELSE
            CALL Gong(1)
            GOTO 2
            ENDIF
         IF (r1.ge.r2) THEN
            CALL Gong (4)
            Print *, ' BAD INPUT/ min < max required'
            GOTO 2
            ENDIF

         END ! rAskRge

      SUBROUTINE rAskRgeTxt (cFra, ein, r1, r2, r1D, r2D)
C     ---------------------------------------------------
            !  JWu 18jan93, generalized 3mar93
         ! Ask for real range r1<r2 but accept also text
         ! Default = ein, or [r1D,r2D] if ein=' '.
         ! User must provide help of is own.
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER  cFra*(*), ein*(*), c1*16, c2*16, aus*80, ans*80
      DIMENSION  rSta(3)

      IF (r1D.gt.r2D) THEN
         CALL Gong (8)
         Print *, ' BAD USE of WuL/ rAskRgeTxt/ bad defaults given :'
         Print *, ' min, max = ', r1D, r2D
         r1D = 0.d0
         r2D = 0.d0
         ENDIF

 2    CONTINUE
      IF (ein.eq.' ') THEN
         ! numerical default :
         CALL NiceNum (r1D, c1, i1)
         CALL NiceNum (r2D, c2, i2)
         CALL Compose2 (aus,cFra,' ['//c1(1:i1)//' '//c2(1:i2)//'] ?')
         CALL FrageH (aus, ans)
      ELSE
         CALL FrageHD (cFra, ans, ein)
         ENDIF

      CALL FindR (ans, 3, nr, rSta, .false.)
      CALL Kontrakt (ans)

      r1 = r1D
      r2 = r2D
      IF     (nr.eq.0 .and. ans.eq.' ') THEN
         ! both defaults accepted
      ELSEIF (nr.eq.0) THEN
         ein = ans ! text given
         RETURN    ! don't test whether r1<r2
      ELSEIF (nr.eq.1 .and. ans.eq.'#,') THEN
         r1 = rSta(1)
         ein = ' '
      ELSEIF (nr.eq.1 .and. ans.eq.',#') THEN
         r2 = rSta(1)
         ein = ' '
      ELSEIF (nr.eq.2 .and.(ans.eq.'#,#' .or.
     *                      ans.eq.'##'))THEN
         r1 = rSta(1)
         r2 = rSta(2)
         ein = ' '
      ELSE
         CALL Gong(1)
         GOTO 2
         ENDIF
      IF (r1.ge.r2) THEN
         CALL Gong (4)
         Print *, ' BAD INPUT/ min < max required'
         GOTO 2
         ENDIF

      END ! rAskRgeTxt

      SUBROUTINE rAskRgeFull (cFra, rL, rH, rLD, rHD)
C     -----------------------------------------------
            ! JWu 31mai91 RangeFrageH. Reduced to call to rAskRgeTxt 3mar93.

         ! ask for a range rL...rH.
         ! rL=rH=0.0 means full range.
         ! rH<rL is rejected
         ! default is rLD...rHD

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      CHARACTER  cFra*(*), aus*80, ein*80

      qFull = (rLD.eq.0. .and. rHD.eq.0.)

      aus = cFra(1:lenU(cFra))//' (* = full range)'
 2    CONTINUE
      IF (qFull) THEN
         rLDi = 0.
         rHDi = 0.
         ein  = '*'
      ELSE
         IF (rLD.gt.rHD) THEN
            ! correct bad default :
            rLDi = rHD
            rHDi = rLD
         ELSE
            rLDi = rLD
            rHDi = rHD
            ENDIF
         ein = ' '
         ENDIF
      CALL rAskRgeTxt (aus, ein, rL, rH, rLDi, rHDi)

      IF     (ein.eq.'?') THEN
         Print *, 'INPUT HELP :'
      ELSEIF (ein.eq.'*') THEN
         rL = 0.
         rH = 0.
      ELSEIF (ein.ne.' ') THEN
         CALL Gong (3)
         GOTO 2
         ENDIF

      END ! rAskRgeFull

      SUBROUTINE rAskOnOff (fra, ein, RR, MR, nR)
C     -------------------------------------------
            !  JWu 3mar93 without format; 5oct93 with Parser

      IMPLICIT REAL*8 (a-h,o-p,r-z)

      PARAMETER (MA=32)
      CHARACTER  fra*(*), ein*(*), aux*80, def*80, sub*8, cl3*3
      DIMENSION  RR(*), RA(32)

C  Encode default:
 2    CONTINUE
      def = ' '
      IF     (mod(nR,2).ne.0) THEN
         CALL Gong (3)
         Print *, ' rAskOnOff/ invalid default'
         nR = 0
      ELSE
         DO i = 2, nR, 2
            CALL NiceNum (RR(i-1), aux, ia)
            IF (i.gt.2) THEN
               CALL Append (def, ' & '//aux)
            ELSE
               CALL Append (def, aux)
               ENDIF
            CALL NiceNum (RR(i), aux, ia)
            CALL Append (def, ' '//aux)
            ENDDO
         ENDIF

      CALL FrageHD (fra, ein, def)

      IF (ein.eq.'?') THEN
         Print '(a)',' INPUT HELP/'
         CALL Say2 ( '    Specify up to '//cl3(MR/2),
     *        ' subranges, separated by a ''&'' character')
         Print '(a)','    Specify a subrange by two real numbers, '//
     *         'separated by a comma or a blank'
         GOTO 2
         ENDIF

      aux = ein
      CALL FindR (aux, MA, nH, RA, .false.)
      CALL Kontrakt (aux)
      la = lenU(aux)
      IF (aux(la:la).eq.'&') aux(la+1:la+1) = ',' ! accept def. for last subr.

C  Loop/ decode subranges :
      iA = 0 ! taken from RA
      nA = 0 ! transcribed into RR
 4    CONTINUE
      CALL TakeVorDel (aux, sub, '&')
cdb      Print *, ' subrange "', sub
cdb      Print *, ' restrang "', aux
      nA = nA + 2
      IF     (nA.gt.MR) THEN
         Print *, ' Too many subranges given'
         GOTO 2
      ELSEIF (sub.eq.'##' .or. sub.eq.'#,#') THEN
         iA = iA + 2
         RR(nA-1) = RA(iA-1)
         RR(nA)   = RA(iA)
      ELSEIF (sub.eq.',#' .and. nA.le.nR) THEN
         iA = iA + 1
         RR(nA)   = RA(iA)
      ELSEIF (sub.eq.'#,' .and. nA.le.nR) THEN
         iA = iA + 1
         RR(nA-1) = RA(iA)
      ELSEIF ((sub.eq.' ' .or. sub.eq.',') .and. nA.le.nR) THEN
      ELSE
         Print *, ' Illegal subrange : '//sub
         GOTO 2
         ENDIF
cdb      Print *, ' => nA, RR ', nA, RR(nA-1), RR(nA)
      IF (aux.ne.' ') GOTO 4
      nR = nA
C  End loop/ decode subranges.

      IF (irSorted(RR, nR).lt.2) THEN
         CALL Gong (3)
         Print *, 'Values must be in ascending order'
         GOTO 2
         ENDIF

cdb debug :
cdb      Print *, ' nR ===== ', nR
cdb      DO i = 1, nR
cdb         Print *, i, ' ----> ', RR(i)
cdb         ENDDO
      END ! rAskOnOff

      SUBROUTINE rAskArray (X, M, n)
C     ------------------------------
            ! JWu 10sep91, help 10mar93
         ! ask for a whole array X(1..n), n<=M

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      CHARACTER         cl6*6, aus*80, ein*80
      DIMENSION         X(*)

      IF (M.lt.1) CALL Absturz('rAskArray', 'M<1')
      n = 0
 1       CONTINUE
         aus = ' Entry '//cl6(n+1)
         CALL Append (aus, ' [end]')
 12      CALL FrageH (aus, ein)
         CALL Fi1R (ein, r)
         IF     (ein.eq.'#') THEN
            n = n + 1
            X(n) = r
         ELSEIF (ein.eq.' ') THEN
            RETURN
         ELSEIF (ein.eq.'?') THEN
            Print *, ' INPUT HELP/'
            Print *, '     enter array elements (real values)'
            Print *, '     after each element, type RETURN'
            Print *, '     type RETURN to end the input'
         ELSE
            CALL Gong (3)
            GOTO 12
            ENDIF
         IF (n.eq.M) THEN
            Print *, ' array completely filled'
            RETURN
            ENDIF
         GOTO 1

      END ! rAskArray

      SUBROUTINE rAskPairs (X, Y, M, n)
C     ---------------------------------  !  WARNING: rAskPairs<>rAskPair
            ! JWu 15jan93
         ! ask for a pair of arrays X,Y(1..n), n<=M

      CHARACTER         cl6*6, aus*60, ein*60
      REAL*8            X(*), Y(*), Rein(2)

      DO i = 1, M
         aus = ' Entry '//cl6(i)
         CALL Append (aus, ' : enter X, Y [end]')
 12      CALL FrageH (aus, ein)

         CALL FindR (ein, 2, nr, Rein, .false.)
         CALL Kontrakt (ein)
         IF     (nr.eq.0 .and. ein.eq.' ') THEN
            n = i - 1
            RETURN
         ELSEIF (nr.eq.2 .and. (ein.eq.'##' .or. ein.eq.'#,#')) THEN
            X(i) = Rein(1)
            Y(i) = Rein(2)
         ELSE
            CALL Gong(3)
            GOTO 12
            ENDIF
         ENDDO

      END ! rAskPairs

      SUBROUTINE rAskGrid (cona, X, M, n)
C     -----------------------------------
            ! obsolete 15dec95
         ! ask for an equidistant grid and set
         ! the whole array X(1..n), n<=M.
         ! cona is the coordinate's name.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      CHARACTER         cona*(*)
      DIMENSION         X(M)

      IF (M.lt.2) CALL Absturz ('rAskGrid', 'M<2')
      IF (n.gt.1 .and. n.le.M) THEN
         xn = X(n)
      ELSE
         xn = X(1) ! no default will be given
         ENDIF
      CALL rAskGridDef (cona, X, M, n, X(1), X(2)-X(1), xn)

      END ! rAskGrid

      SUBROUTINE rAskGridDef (cona, X, M, n, x1d, dxd, xnd)
C     -----------------------------------------------------
            ! JWu 27feb92; defaults improved 10feb95, again (rAGD) 15dec95
         ! ask for an equidistant grid and set
         ! the whole array X(1..n), n<=M.
         ! cona is the coordinate's name.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      CHARACTER         cona*(*)
      DIMENSION         X(M)

      IF (M.lt.2) CALL Absturz ('rAskGridDef', 'M<2')
      lcn = lenU(cona)
 1    CONTINUE
         X(1) = rAskD (' Enter '//cona(1:lcn)//'(1)', x1d)
         IF (dxd.gt.0) THEN
            dx = rAskD (' Step in '//cona(1:lcn), dxd)
         ELSE
            dx = rAsk (' Step in '//cona(1:lcn))
            ENDIF
         IF (dx.le.1.d-30 .or. dx.gt.1.d30) GOTO 2
         IF (qroutside(xnd, x1d+dxd, x1d+(M-1)*dxd)) THEN
            Xmax  = rAsk (' Enter '//cona(1:lcn)//'(n)')
         ELSE
            Xmax  = rAskD (' Enter '//cona(1:lcn)//'(n)', xnd)
            ENDIF
         IF (Xmax.lt.dx) GOTO 2

         n   = (Xmax - X(1) + dX/1.d5) / dX + 1
         Print *, ' => number of points = ', n
         IF (n.gt.M) THEN
            CALL Gong (1)
            Print *, ' allowed maximum is ', M
            GOTO 1
            ENDIF

         DO i1 = 1, n ! set equidistant scale :
            X(i1) = X(1) + (i1-1)*dx
            ENDDO
         RETURN

C  Error treatment / escape :
 2    CONTINUE
         CALL Gong (3)
         Print *, ' BAD INPUT'
         IF (qAsk(' Correct ?')) GOTO 1
         n = 0
         RETURN

      END ! rAskGridDef
