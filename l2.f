C  ====================================================================
C
C     Library  WuL :  General FORTRAN Library
C     Module   L2  :     String Manipulation
C
C  ====================================================================

C     Contents :
C        2.1.  Integer and Logical Functions :
C                 jPos1, lenU, qSubStrEq
C        2.2.  Concatenation and Insertion :
C                 Append, Compose_, Insert
C        2.3.  Cutting and Deleting :
C                 StaTake, DelVonBis,
C                 TakeVonBis, TakeVor, TakeVorDel
C        2.4.  Replacing :
C                 ReplaceC/T, Minuskeln, Majuskeln, DelLeft, Kontrakt
C        2.5.  General Decoding :
C                 Parser
C        2.6.  Search for numbers :
C                 FindN, FindI, FindR, Fi1N, Fi1I, Fi1R
C        2.7.  Integer Lists :
C                 DecNRange, DecNList, ExpNList,
C                 DecJList, EncJList
C        2.9.  Parentheses :
C                 jPos1Kla, TakeVorDelKla, TakeVorKla
C
C     Notes :
C        Maximum string length *2048 may be
C        implementation dependent.

C     History :
C        JWu  mar93      clean-up, L2 and L3 joined together
C        JWu  oct90      Revision, VMS, WuLib, Module WuL2
C        JWu  oct/nov89  Revision, Lib A, Module A_
C        JWu  1989       fuer NumGraph

C  ====================================================================
C  L2.1.  Integer and Logical Functions :
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      INTEGER FUNCTION jPos1 (stri, substri)
C     --------------------------------------                              JW 89
      !  Position des ersten Zeichens von substri
      !  beim ersten Auftreten von substri in stri.
      !  Wenn substri nicht in stri enthalten ist,
      !  wird jPos1 = len(stri) + 1 gesetzt.

         CHARACTER *(*) stri, substri
         ls  = len(stri)
         IF (ls.gt.512) THEN ! NY : detect abnormous len's
            Print *, ' ls = ', ls
            Print *, ' stri(1..320) = "', stri(1:320), '".'
            CALL Absturz1 ('jPos1', 'stri>512 z.Z. ungewollt')
            ENDIF
         lss = len(substri)
         DO 1 i = 1,ls-(lss-1)
            IF (stri(i:i+(lss-1)).eq.substri(1:lss)) THEN
               jPos1 = i
               RETURN
               ENDIF
 1          CONTINUE
         jPos1 = len(stri) + 1
         END

      INTEGER FUNCTION lenU (stri)
C     ----------------------------
            ! JWu 1989, rev. 2nov89.
         ! length of stri without blanks on the right.
      CHARACTER *(*) stri
      i  = len(stri)
 1    CONTINUE ! DO-Schleife hat nicht funtioniert
         ic = ichar(stri(i:i))
         IF (ic.ne.32 .and. ic.ne.0) GOTO 2
         i = i - 1
         IF (i.gt.0) GOTO 1
 2    lenU = i

      END ! lenU

      LOGICAL FUNCTION qSubStrEq (stri, i1, sub)
C     ------------------------------------------
         ! JWu 23jun93
         ! test if stri(i1..)==sub; false if i1 out of range.

      CHARACTER*(*) stri, sub

      l1 = len(stri)
      i2 = i1-1+len(sub)
      IF     (i1.lt.1) THEN
         qSubStrEq = .false.
      ELSEIF (i2.gt.l1) THEN
         qSubStrEq = .false.
      ELSE
         qSubStrEq = ( stri(i1:i2).eq.sub )
         ENDIF

      END ! qSubStrEq

C  ====================================================================
C        2.2.  Concatenation and Insertion :
C  ====================================================================

      SUBROUTINE Append (Stack, In)
C     -----------------------------
            ! JWu 1989, corr oct91.
         ! Appends In at the position of the last blank in Stack.
      CHARACTER Stack *(*), In *(*)
      le = lenU (Stack)
      IF (le.gt.0) THEN
         Stack = Stack(1:le)//In
      ELSE ! Stack is empty
         Stack = In
         ENDIF
      END ! Append

      SUBROUTINE Compose2 (out, in1, in2)
C     -----------------------------------
            ! ComposeN : JWu 30sep91
         ! append in2 on in1 and write it to out
         CHARACTER *(*) out, in1, in2
         out = in1
         CALL Append (out, in2)
         END ! Compose2

      SUBROUTINE Compose3 (out, in1, in2 , in3)
C     -----------------------------------------
         CHARACTER *(*) out, in1, in2, in3
         out = in1
         CALL Append (out, in2)
         CALL Append (out, in3)
         END ! Compose3

      SUBROUTINE Compose4 (out, in1, in2 , in3, in4)
C     ----------------------------------------------
         CHARACTER *(*) out, in1, in2, in3, in4
         out = in1
         CALL Append (out, in2)
         CALL Append (out, in3)
         CALL Append (out, in4)
         END ! Compose4

      SUBROUTINE Compose5 (out, in1, in2 , in3, in4, in5)
C     ---------------------------------------------------
         CHARACTER *(*) out, in1, in2, in3, in4, in5
         out = in1
         CALL Append (out, in2)
         CALL Append (out, in3)
         CALL Append (out, in4)
         CALL Append (out, in5)
         END ! Compose5

      SUBROUTINE Compose6 (out, in1, in2 , in3, in4, in5, in6)
C     --------------------------------------------------------
         CHARACTER *(*) out, in1, in2, in3, in4, in5, in6
         out = in1
         CALL Append (out, in2)
         CALL Append (out, in3)
         CALL Append (out, in4)
         CALL Append (out, in5)
         CALL Append (out, in6)
         END ! Compose6

      SUBROUTINE Insert (stri, ji, substri)
C     -------------------------------------
            ! JWu 23oct89, refait 16sep91.
       ! Setzt substri als stri(ji:..) ein.

         CHARACTER *(*) stri, substri

         ls    = len (stri)
         lss   = len (substri)
         IF (lss.eq.0) RETURN
         IF (ji.gt.ls) RETURN

         DO i = ls, ji+lss, -1
            stri(i:i) = stri(i-lss:i-lss)
            ENDDO
         stri(ji:ji+lss-1) = substri
         ! stri(1:ji-1) unchanged

         END ! Insert

C  ====================================================================
C  2.3.  Cutting and Deleting :
C  ====================================================================

      SUBROUTINE StaTake (n, stri, out)
C     ---------------------------------
         !  JWu 89, rev 20.10.89
       ! nimmt die n ersten Stellen aus stri
       ! und schreibt sie nach out
         CHARACTER *(*) stri, out
         IF (n.lt.1) THEN
            out = ' '
            RETURN
            ENDIF
         out   = stri (1:n)
         ls    = len (stri)
         IF (n.lt.ls) THEN
            stri = stri ((n+1):ls)
         ELSE
            stri = ' '
            ENDIF
         END ! StaTake

      SUBROUTINE DelVonBis (stri, j1, j2)
C     -----------------------------------
            ! JWu 20.10.89, erg. 23.10.89, corr. 1. 2.90
         ! Kuerzt stri um stri(j1:j2)

         CHARACTER *(*) stri

         l  = len (stri)
         IF (j1.gt.j2) CALL Absturz1 ('DelVonBis', 'j1>j2')
         IF (j2.gt. l) CALL Absturz1 ('DelVonBis', 'j2>l')

         ld = j2 - j1 + 1
         IF (j2.lt.l) stri (j1:l-ld) = stri (j2+1:l)
         DO i = l-ld+1, l
            stri (i:i) = ' '
            ENDDO

         END ! DelVonBis

      SUBROUTINE TakeVonBis (stri, j1, j2, substri)
C     ---------------------------------------------                JW 23.10.89
       ! Kuerzt stri um stri(j1:j2), welchselbiges nach
       ! substri geschrieben wird

         CHARACTER *(*) stri, substri

         l  = len (stri)
         ls = len (substri)
         jdel = j2 - j1 + 1
         IF (j1.gt.j2) CALL Absturz1 ('TakeVonBis', 'j1>j2')
         IF (j2.gt. l) CALL Absturz1 ('TakeVonBis', 'j2>l')
         IF (jdel.gt.ls) CALL Absturz1 ('TakeVonBis','substri zu kurz')

         substri = stri (j1:j2)
         IF (j2.lt.l) THEN
            stri (j1:l-jdel) = stri (j2+1:l)
            ENDIF
         DO i = l-jdel+1, l
            stri (i:i) = ' '
            ENDDO

         END ! TakeVonBis

      SUBROUTINE TakeVor (Stack, Out, Substri)
C     ----------------------------------------                           JW 89
       ! Nimmt alles vor dem ersten Auftreten von Substri
       ! aus Stack und schreibt es nach Out

         CHARACTER *(*) Stack, Out, Substri
         n     = jPos1 (Stack, Substri) - 1
         CALL StaTake (n, Stack, Out)
         END ! TakeVor

      SUBROUTINE TakeVorDel (Stack, Out, Substri)
C     -------------------------------------------
         !  JWu 89, rev 8.11.89
      !  Nimmt alles vor dem ersten Auftreten von Substri
      !  aus Stack und schreibt es nach Out;
      !  entfernt ausserdem Substri aus Stack.

         CHARACTER *(*) Stack, Out, Substri
         CHARACTER *2048 muell
         n     = jPos1 (Stack, Substri) - 1
         CALL StaTake (n, Stack, Out)
         CALL StaTake (len(Substri), Stack, muell)

         END ! TakeVorDel

C  ====================================================================
C     2.4. Replacing :
C  ====================================================================

      SUBROUTINE ReplaceC (stri, alt, neu)
C     ------------------------------------
            ! JWu 89 "Ersetz", jun91 "ReplaceC"
         ! Ersetzt in stri alt durch neu
         CHARACTER *(*) stri, alt, neu
         CHARACTER *2048 Hilf

         ls = len (stri)
         la = len (alt)
         ln = len (neu)
         l  = 0
         lh = 0
 1       CONTINUE ! copy (modified) stri to Hilf
            l = l + 1 ! search for alt in stri(l:_)
            lh= lh+ 1
            IF (l+la-1.gt.ls) GOTO 2 ! no more place for alt - stop searching
            IF (stri(l:l+la-1).eq.alt) THEN ! found alt
               Hilf(lh:lh+ln-1) = neu
               lh = lh+ln-1
               l  = l+la-1
            ELSE
               Hilf(lh:lh) = stri(l:l)
               ENDIF
            GOTO 1
 2       CONTINUE
         IF (l.le.ls) THEN ! copy the rest
            lrest = ls-l
            Hilf(lh:lh+lrest) = stri(l:ls)
            Hilf(lh+lrest+1:len(Hilf)) = ' ' ! 17sep95
            ENDIF
         stri = Hilf ! copy back to stri
         END ! ReplaceC

      SUBROUTINE ReplaceThisC (stri, alt, neu)
C     ----------------------------------------
            ! JWu 6sep91 for use in FPP.
         ! if stri(1:_)=alt, then replace this occurence of alt by new.
         CHARACTER *(*) stri, alt, neu
         CHARACTER *2048 Hilf

         ls = len (stri)
         la = len (alt)
         ln = len (neu)

         IF (stri(1:la).eq.alt) THEN
            Hilf = neu
            Hilf(ln+1:ls+ln-la) = stri(la+1:ls)
            stri = Hilf
            ENDIF

         END ! ReplaceThisC

      SUBROUTINE ReplaceT (stri, alt, neu)
C     ------------------------------------
              ! JWu 11feb91, corr. 30sep91, alt(1:lenU) 25oct91.
         ! Replace alt by neu in stri. Contrarily to ReplaceC,
         ! ending blanks in alt and neu are ignored.

      CHARACTER *(*) stri, alt, neu

      ln = lenU(neu)
      la = lenU(alt)
      IF (la.lt.1) THEN
         CALL Absturz1 ('ReplaceT', 'Empty string cannot be replaced')
         ENDIF
      IF (ln.ge.1) THEN
         CALL ReplaceC (stri, alt(1:la), neu(1:ln) )
      ELSE
         DO i = 1, len(stri)-la+1
            IF (stri(i:i+la-1).eq.alt)
     *         CALL DelVonBis (stri, i, i+la-1)
            ENDDO
         ENDIF
      END ! ReplaceT

      SUBROUTINE Minuskeln (stri)
C     ---------------------------
            ! JWu 3nov89
         ! Ersetzt in stri Majuskeln durch Minuskeln
         CHARACTER *(*) stri
         ica = ichar ('A')
         icz = ichar ('Z')
         id  = ichar ('a') - ichar ('A')

         DO 1 j = 1,len(stri)
            ic = ichar(stri(j:j))
            IF (ic.ge.ica .and. ic.le.icz) stri(j:j) = char(ic+id)
 1          CONTINUE
         END ! Minuskeln

      SUBROUTINE Majuskeln (stri)
C     ---------------------------
            ! JWu 3.11.89
       ! Ersetzt in stri Minuskeln durch Majuskeln
         CHARACTER stri *(*)
         ica = ichar ('a')
         icz = ichar ('z')
         id  = ichar ('A') - ichar ('a')

         DO 1 j = 1,len(stri)
            ic = ichar(stri(j:j))
            IF (ic.ge.ica .and. ic.le.icz) stri(j:j) = char(ic+id)
 1          CONTINUE

         END ! Majuskeln

      SUBROUTINE DelLeft (stri)
C     -------------------------
               ! JWu 30jul91
            ! Delete blanks on the left side of the text contained in string.
         CHARACTER *(*) stri
         l = len(stri)
         DO i = 1, l
            IF (stri(i:i).ne.' ') GOTO 1
            ENDDO
         RETURN ! all is blank
 1       CONTINUE
         IF (i.eq.1) RETURN ! no blanks
         CALL DelVonBis (stri, 1, i-1)

         END ! DelLeft

      SUBROUTINE Kontrakt (stri)
C     --------------------------
            ! JWu 20.10.89, rev  2.11.89
       ! Kontrahiert stri durch Auslassen aller Leerzeichen

         CHARACTER *(*) stri

         l    = lenU (stri)
         IF (l.eq.0) RETURN

         i = 1
 1       CONTINUE
            IF (stri(i:i).eq.' ') THEN
               IF (i.lt.l) THEN
                  stri (i:l) = stri (i+1:l) // ' '
                  l = l - 1
                  GOTO 1
               ELSE
                  stri (i:l) = ' '
                  RETURN
                  ENDIF
            ELSE
               i = i + 1
               ENDIF
            IF (i.le.l) GOTO 1

         END ! Kontrakt

C  ====================================================================
C  L2.5.  General Decoding
C  ====================================================================

      SUBROUTINE Parser (ein, aus, M, n, tAus, iAus, rAus, dir, Fehler)
C     -----------------------------------------------------------------
            ! JWu 25nov92.
         ! Input :   ein    : any text
         !           M      : declared length of arrays _Aus
         !           dir    : directives
         ! Output :  aus    : compact version of ein
         !           n      : used length of aus
         !           _Aus   : text and numbers found in ein
         !           Fehler : <>'&ff' if errors occured
         ! Directives : n   : convert non-negativ numbers
         !              i   : convert integer numbers
         !              r   : convert real numbers
         !              m   : make letters small
         !              w   : bundle words
         !              _   : eliminate blanks
         !              p   : parentheses uncoded into tAus
         !              q   : quotations uncoded into tAus
         ! Codings :    n   : number
         !              c   : character or word
         !              p   : contents of parenthesis
         !              q   : contents of quote

         ! Example :
         ! when directives are 'npw', the text ein='pi=3(oder 4)'
         ! is parsed as aus='c=np', thus n=4,
         ! and the arrays contain tAus='pi','=','3','oder 4'
         ! and iAus=-,-,3,-.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      CHARACTER*(*) ein, aus, tAus, dir, Fehler
      DIMENSION     tAus(*), rAus(*), iAus(*)

      CHARACTER     ci*1

      DATA          iNa,iNz, iKa,iKz, iGa,iGz /48,57, 97,122, 65,90/
                    ! ichar of 0,9, a,z, A,Z

C  Checks :
      IF     (Fehler.ne.'&ff') THEN
         Fehler = 'error on entry in Parser'
         RETURN
      ELSEIF (M.lt.1) THEN
         CALL Absturz ('Parser', 'M<1')
      ELSEIF (len(aus).lt.M) THEN
         Fehler = 'Parser/ aus too short'
         ENDIF

C  Presets : don't decode anything :
      idecnum = 0
      qSmalls = .false.
      qParent = .false.
      qQuotes = .false.
      qWords  = .false.
      qEliBla = .false.

C  Decode directives :
      DO i = 1, len(dir)
         IF     (dir(i:i).eq.'n') THEN
             idecnum = 1
         ELSEIF (dir(i:i).eq.'i') THEN
             idecnum = 2
         ELSEIF (dir(i:i).eq.'r') THEN
             idecnum = 3
         ELSEIF (dir(i:i).eq.'m') THEN
             qSmalls = .true.
         ELSEIF (dir(i:i).eq.'p') THEN
             qParent = .true.
         ELSEIF (dir(i:i).eq.'q') THEN
             qQuotes = .true.
         ELSEIF (dir(i:i).eq.'w') THEN
             qWords  = .true.
         ELSEIF (dir(i:i).eq.'_') THEN
             qEliBla = .true.
         ELSEIF (dir(i:i).eq.' ') THEN
         ELSE
             Fehler = 'invalid parser directive : '//dir
             RETURN
             ENDIF
         ENDDO


      IF (idecnum.ne.1) THEN
         Fehler = 'Parser option ''?'' not yet implemented'
         RETURN
         ENDIF
      IF (qWords) THEN
         Fehler = 'Parser option ''w'' not yet implemented'
         RETURN
         ENDIF
      IF (qParent) THEN
         Fehler = 'Parser option ''p'' not yet implemented'
         RETURN
         ENDIF
      IF (qQuotes) THEN
         Fehler = 'Parser option ''q'' not yet implemented'
         RETURN
         ENDIF

C  Now go on :
      aus = ' ' ! preset
      n = 0     ! index(aus)
      i = 0     ! index(ein)
 20   CONTINUE ! next character
         IF (qParserGet(ein, i, ci, ii)) THEN
c             Print *, ' ein : "', ein(1:lenU(ein)), '", i=',i
c             Print *, ' aus : "', aus(1:lenU(aus)), '", n=',n
c             DO i = 1, n
c                Print *, i, iAus(i), rAus(i), '"', tAus(i), '"'
c                ENDDO
            RETURN ! regular exit
            ENDIF
 21      CONTINUE

         IF     (ci.eq.' ') THEN
            IF (qEliBla) GOTO 20 ! eliminate all blanks
            IF (qParserPut(aus, M, n, ' ', Fehler)) RETURN
            tAus(n) = ' '
            iAus(n) = 0
            rAus(n) = 0.
            ! eliminate following blanks :
 301        CONTINUE
               IF (qParserGet(ein, i, ci, ii)) RETURN
               IF (ci.eq.' ') GOTO 301
               GOTO 21
         ELSEIF (qiinside(ii,iKa,iKz) .or.
     *           qiinside(ii,iGa,iGz)) THEN   ! it's a letter
            IF (qParserPut(aus, M, n, 'c', Fehler)) RETURN
            IF (qSmalls .and. ii.le.iGz) THEN ! make small
               tAus(n) = char(ii+iKa-iGa)
            ELSE
               tAus(n) = ci
               ENDIF
            iAus(n) = 1
            rAus(n) = 0.
         ELSEIF (qiinside(ii,iNa,iNz)) THEN ! it's a digit
            IF (idecnum.le.0) GOTO 391 ! no treatment
            istnum = ii - iNa
            ia = i ! first digit
 321        CONTINUE
               IF (.not.qParserGet(ein, i, ci, ii)) THEN
                  IF (qiinside(ii,iNa,iNz)) GOTO 327 ! kommen noch mehr Ziffern
                  ENDIF
               ! integer ended in ein(i-1)
               ie = i-1 ! last digit
               IF (qParserPut(aus, M, n, 'n', Fehler)) RETURN
               tAus(n) = ein(ia:ie)
               iAus(n) = istnum
               rAus(n) = dble(iAus(n))
c               Print *, 'SAVE/ ', istnum
               GOTO 21 ! das war's, fuer diese Zahl
 327           CONTINUE
               IF (istnum.gt.1.d8) THEN ! maxint=2.1E9
                  Fehler = 'Integer value too large'
                  RETURN
                  ENDIF
               istnum = 10*istnum + (ii-iNa)
c                  Print *, 'ADD/ ', ii, ' -> ', istnum
               GOTO 321
         ELSE ! anything else
            GOTO 391
            ENDIF
         GOTO 20

 391  CONTINUE
      IF (qParserPut(aus, M, n, ci, Fehler)) RETURN
      tAus(n) = ci
      iAus(n) = 0
      rAus(n) = 0.
      GOTO 20

      END ! Parser

      LOGICAL FUNCTION qParserGet (ein, i, ci, ii)
C     --------------------------------------------
            ! used EXCLUSIVELY by Parser.
      CHARACTER ein*(*), ci*1
      i = i + 1
      IF (i.gt.len(ein)) THEN
         qParserGet = .true.
      ELSE
         ci = ein(i:i)
         ii = ichar(ci)
c               Print *, 'PG/ ', i, ': "', ci, '" ->', ii
         qParserGet = .false.
         ENDIF
      END ! qParserGet

      LOGICAL FUNCTION qParserPut (aus, M, n, cn, Fehler)
C     ---------------------------------------------------
            ! used EXCLUSIVELY by Parser.
      CHARACTER aus*(*), cn*1, Fehler*(*)
      IF (n.ge.M) THEN
         Fehler = 'Parser overflow'
         qParserPut = .true.
      ELSE
         n = n + 1
         aus(n:n) = cn
         qParserPut = .false.
         ENDIF
      END ! qParserPut

C  ====================================================================
C  L2.6.  Search for numbers :
C  ====================================================================

      SUBROUTINE FindN (stri, Mi, ni, Inte)
C     -------------------------------------
            ! JWu 89, slight modif 20jun91
       ! Findet Integer >= 0 in stri, schreibt sie nach Inte, und
       ! ersetzt sie in stri durch '#'.
       ! Es koennen maximal Mi Integer gefunden werden, die
       ! tatsaechlich gefundene Anzahl steht in ni.

         IMPLICIT LOGICAL (q)
         CHARACTER stri *(*), c *1
         INTEGER   Inte (*)

         IF (Mi.lt.1) CALL Absturz1 ('FindN','Mi<1')

         l    = len (stri)
         ni   = 0
         i    = 0

C  Check absence of '#' : (22jan92)
         IF (jPos1(stri,'#').le.l) THEN
            CALL Gong (7)
            Print *, 'FindN/ stri = ', stri(1:lenU(stri)), '"'
            Print *, 'FORBIDDEN INPUT/ "#" where numbers expected'
            stri = ' '
            RETURN
            ENDIF

C  Search for digits :
         qZahl= .false.
 1       CONTINUE
            i  = i + 1
            IF  (i.gt.l) RETURN ! legal EXIT
            c  = stri(i:i)
            ia = ichar(c)-ichar('0')
            IF (qiinside(ia,0,9)) THEN ! Ziffer gefunden
               IF (qZahl) THEN ! an aufgebaut werdende Zahl anhaengen
                  IF (Inte(ni).gt.2.d8) THEN ! MaxInt = 2.1E9 (27jan92)
                     Print *, 'FORBIDDEN INPUT/ Number too long'
                     stri(i:l) = ' '
                     ni = ni - 1
                     RETURN
                     ENDIF
                  Inte(ni) = Inte(ni) * 10 + ia
                  IF (i.lt.l) THEN
                     stri(i:l) = stri(i+1:l) // ' '  ! stri aufruecken
                  ELSE
                     stri(i:i) = ' ' ! end of string 4may95
                     ENDIF
                  l = l - 1
                  i = i - 1
               ELSE ! neue Zahl aufbauen
                  ni = ni + 1
                  Inte (ni) = ia
                  stri (i:i) = '#' ! Marke in stri
                  qZahl = .true.
                  ENDIF
            ELSE
               IF (qZahl) THEN
                  qZahl = .false.
                  IF (ni.eq.Mi) RETURN ! there is no place for more numbers
                  ENDIF
               ENDIF
            GOTO 1

         END ! FindN

      SUBROUTINE FindI (stri, Mi, ni, Ista)
C     -------------------------------------
            ! JWu 1989/90
      !   finds signed integers

         IMPLICIT LOGICAL (q)
         CHARACTER        stri *(*), c *1, c2 *2
         INTEGER          Ista (*)

         CALL FindN ( stri, Mi, ni, Ista )
         lS = len(stri)

         is = 0
         DO 1 ii = 1,ni
 11         CONTINUE ! searching the single '#'
               is = is + 1
               IF (is.gt.lS) THEN
                  Print *, ' stri = "', stri, '"'
                  Print *, ' ii ni is lS ',ii,ni,is,lS
                  CALL Absturz1 ('FindI', 'Nonsense result from FindN')
                  ENDIF
               IF (stri(is:is).eq.'#') GOTO 12
               GOTO 11
 12         CONTINUE ! found '#'
            IF (is.eq.1) GOTO 1
            c2 = stri(is-1:is)
            IF (c2.eq.'-#') Ista(ii) = - Ista(ii)
            IF (c2.eq.'-#' .or. c2.eq.'+#') THEN
               CALL TakeVonBis (stri, is-1, is-1, c)
               ENDIF
 1          CONTINUE

         END ! FindI

      SUBROUTINE FindR (stri, Mr, nr, Rsta, qAuchKomma)
C     -------------------------------------------------
            ! JWu 1989(?), renewed oct90 for VMS
      !   find Real in stri.
      !   Format : # #. #.# #e# #E# #.#e# #m# ...
      !   Wenn qAuchKomma auch "," statt "."

      !   Ist leider etwas kompliziert geworden. Liegt aber nur daran,
      !   dass die NOS/VE - Version unter VMS nicht lief.


         IMPLICIT LOGICAL (q)
         CHARACTER        stri *(*), c *1, aux*80
         INTEGER          Nexp (1)
         REAL *8          Rsta (*), r

         lS = len(stri)

C  Check absence of '#' : (27jan92)
         IF (jPos1(stri,'#').le.lS) THEN
            CALL Gong (7)
            Print *, 'FORBIDDEN INPUT/ "#" where numbers expected'
            stri = ' '
            RETURN
            ENDIF

         is = 0
         nr = 0
 1       CONTINUE
C  search for begin of number :
            qNeu  = .false.
            r     = 0.
 12         CONTINUE
            qNega = .false.
            is = is + 1        ! position in stri to be read
            i1 = is            ! 1st position of number, if found
            IF (is.gt.lS) RETURN
            c = stri(is:is)
            IF ((c.eq.'+' .or. c.eq.'-') .and. is.lt.lS) THEN
C  a sign was found, use only if number follows :
               qNega = (c.eq.'-')
               is = is + 1
               IF (is.gt.lS) RETURN
               c = stri(is:is)
               ENDIF
            IF (jPos1('0123456789',c).le.10) THEN
C  a digit found - decode integer part :
               CALL Fi1N (stri(is:lS), nip)
               IF (stri(is:is).ne.'#') THEN
                  Print *, 'FindR/ bad integer from Fi1N'
                  stri = ' '
                  nr = 0
                  RETURN
                  ENDIF
               r   = dfloat (nip)     ! value of integer part
               qNeu= .true.           ! number is being read
               is  = is + 1
               IF (is.gt.lS) GOTO 19
               c   = stri(is:is)
               ENDIF
            IF (c.eq.'.' .or. (qAuchKomma .and. c.eq.',')) THEN
C  decimal fraction expected :
               id = 0
 16            CONTINUE
               is = is + 1
               id = id + 1
               IF (is.gt.lS) GOTO 19
               c = stri(is:is)
               IF (jPos1('0123456789',c).le.10) THEN
                     ! 13sep94 completely new to allow for Bartsch bandworms
                  ic = ichar(c) - ichar('0')
                  IF (ic.gt.0) qNeu = .true.
                  r  = r + dble(ic) / 10.d0**id
                  GOTO 16
                  ENDIF
               ENDIF
            IF     (qNeu .and. (c.eq.'e' .or. c.eq.'E')) THEN
C  decadic exponent expected :
               is = is + 1
               IF (is.gt.lS) GOTO 19
               CALL FindI (stri(is:lS), 1, ni, Nexp)
               IF (ni.lt.0 .or. stri(is:is).ne.'#') THEN
C  no, it was just a letter E, nothing more :
                  is = is - 1
                  GOTO 19
                  ENDIF
C  decadic exponent found :
               is = is + 1
               r = r * 10.d0**Nexp(1)
            ELSEIF (qNeu .and. c.eq.'m') THEN
C  decadic exponent expected :
               is = is + 1
               IF (is.gt.lS) GOTO 19
               CALL FindI (stri(is:lS), 1, ni, Nexp)
               IF (ni.lt.0 .or. stri(is:is).ne.'#') THEN
C  no, it was just a letter m, nothing more :
                  is = is - 1
                  GOTO 19
                  ENDIF
C  decadic exponent found :
               is = is + 1
               r = r * 10.d0**(-Nexp(1))
               ENDIF
  19        CONTINUE
C  if there was a number, all of it is read :
            IF (qNeu) THEN
               IF (qNega) r = -r
               nr = nr + 1
               IF (nr.gt.Mr) RETURN
               Rsta (nr) = r
C  now what is left of the number is replaced by '#' :
               CALL TakeVonBis (stri, i1, is-1, aux)
               CALL Insert (stri, i1, '#')
C  continue search for the next number :
               is = i1 + 1  ! this line added 09-11-90
               ENDIF
            GOTO 1
         END ! FindR

      SUBROUTINE Fi1N (stri, i)
C     -------------------------
            ! JWu 1989, reduced 1mar93
       ! Findet das erste Integer in stri, schreibt es
       ! nach i und ersetzt es in stri durch '#'.

         CHARACTER stri *(*)
         INTEGER   Inte (1)

         CALL FindN (stri,1,n,Inte)
         IF (n.ge.1) THEN
            i = Inte(1)
            ENDIF

         END ! Fi1N

      SUBROUTINE Fi1I (stri, i)
C     -------------------------
            ! JWu 14jan92, reduced 1mar93
         ! Find the first signed integer in stri and return
         ! its value i; in stri the number is replaced by '#'.

      CHARACTER stri*(*)
      INTEGER   Inte(1)

      CALL FindI (stri,1,n,Inte)
      IF (n.ge.1) THEN
         i = Inte(1)
         ENDIF

      END ! Fi1I

      SUBROUTINE Fi1R (stri, r)
C     -------------------------
            ! JWu 3mar93
         ! Find the first real number in stri and return
         ! its value i; in stri the number is replaced by '#'.

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      CHARACTER        stri*(*)
      DIMENSION         RR(1)

      CALL FindR (stri, 1, n, RR, .false.)
      IF (n.ge.1) THEN
         r = RR(1)
         ENDIF

      END ! Fi1R

C  ====================================================================
C  L2.7.  Integer Lists
C  ====================================================================

      SUBROUTINE DecNRange (stri, jmx, ji, jf, Fehler )
C     -------------------------------------------------
            ! JWu 20.10.89, 12.2.91, underlay of dppIFS 27nov92
         ! decodes Wuttke n-range format as described in
         ! FrageNRange, see WuL4, section 4.3

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      CHARACTER*(*)     stri, Fehler
      PARAMETER        (Mcomp=22)
      CHARACTER         comp*22, tComp*20
      DIMENSION         tComp(Mcomp), iComp(Mcomp), rComp(Mcomp)

C  Parse :
      CALL Parser (stri, comp, Mcomp, nc,
     *             tComp, iComp, rComp, 'nm_', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Dec post parser :
      CALL dppIFS (1, nc, nc, comp, tComp, iComp,
     *             1, jmx, ji, jf, js, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'DecNRange/ "', stri(1:lenU(Stri)), '"'
         RETURN
         ENDIF

      IF     (js.ne.1) THEN
         Fehler = 'range with step=1 expected'
      ELSEIF (ji.gt.jf) THEN
         Fehler = 'range ji-jf with ji<=jf expected'
         ENDIF

      END ! DecNRange

      SUBROUTINE dppIFS (ia, if, nc, comp, tComp, iComp,
     *                   MI, MF, ji, jf, js, Fehler)
C     --------------------------------------------------
            ! JWu 27nov92
         ! decode post parser - initial final step
         ! format : n-nin or n-n or n or * or n-* or *in or ...

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      CHARACTER*(*) comp, tComp, Fehler
      DIMENSION     tComp(*), iComp(*)
      CHARACTER     cl6*6

C  Checks :
      IF     (nc.lt.1) THEN
         Fehler = 'dppIFS/ nothing to decode'
         RETURN
      ELSEIF (ia.lt.1 .or. if.gt.nc .or. ia.gt.if) THEN
         Print *, ' ia if nc ', ia, if, nc
         Fehler = 'dppIFS/ invalid subrange of comp'
         RETURN
      ELSEIF (MI.gt.MF) THEN
         Fehler = 'dppIFS/ no integer subrange allowed'
         RETURN
         ENDIF

C  Decode :
      js = 1 ! default
      IF     (comp(ia:ia).eq.'n') THEN
         ji = iComp(ia)
         IF (qioutside(ji,MI,MF)) THEN
            CALL Compose3 (Fehler, 'range starting with '//cl6(ji),
     *           ' outside allowed '//cl6(MI), '...'//cl6(MF))
            RETURN
            ENDIF
         jf = ji ! default
      ELSEIF (comp(ia:ia).eq.'*') THEN
         ji = MI
         jf = MF ! default
      ELSE
         Fehler = 'expected number or "*" in range'
         RETURN
         ENDIF
      IF (if.eq.ia) RETURN

      IF     (comp(ia+1:ia+1).eq.'c') THEN
         ic = ia+1
         GOTO 30
      ELSEIF (comp(ia+1:ia+1).eq.'-') THEN
         IF     (comp(ia:ia).eq.'*') THEN
            Fehler = 'format "*-.." is undefined'
            RETURN
         ELSEIF (if.lt.ia+2) THEN
            Fehler = 'range ended after "-"'
            RETURN
         ELSEIF (comp(ia+2:ia+2).eq.'n') THEN
            jf = iComp(ia+2)
            IF (qioutside(jf,MI,MF)) THEN
               Fehler = 'range ends outside allowed integers'
               RETURN
               ENDIF
         ELSEIF (comp(ia+2:ia+2).eq.'*') THEN
            jf = MF
         ELSE
            Fehler = 'expected number after "-"'
            RETURN
            ENDIF
         ic = ia+3
         ! now two regular possibilities :
         IF (ic.gt.if) RETURN            ! eoi
         IF (comp(ic:ic).eq.'c') GOTO 30 ! js follows
         Fehler = 'range followed by invalid char'
         RETURN
      ELSE ! added 22aug94 on remark by O.Randl
         Fehler = 'invalid character in range'
         RETURN
         ENDIF
 30   CONTINUE ! found 'c'
      IF (tComp(ic).ne.'i') THEN
         Fehler = 'valid range followed by invalid letter'
         RETURN
         ENDIF
      IF (ic+1.gt.if) THEN
         Fehler = 'range ended after "i"'
         RETURN
         ENDIF
      js = iComp(ic+1)

      END ! dppIFS

      SUBROUTINE DecNList (stri, qLis, mL, Fehler)
C     ---------------------------------------------
            ! JW 20.10.89; increment 7.2.91; all but .. 24oct91.
            ! fast version using Parser completely new 27nov92
         ! Returns qLis (il) = { il ist in der Liste stri enthalten }.
         ! For an explanation of the Wuttke N-list format, see the
         ! subroutine FrageNList in WuL4, section 4.3.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      CHARACTER *(*)    stri, Fehler
      DIMENSION         qLis(*)

      PARAMETER        (Mcomp=320)
      CHARACTER         comp*320, tComp*20
      DIMENSION         tComp(Mcomp), iComp(Mcomp), rComp(Mcomp)

C  Parse :
      CALL Parser (stri, comp, Mcomp, nc, tComp, iComp,
     *             rComp, 'nm_', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Special cases and preset :
      IF     (comp.eq.' ') THEN
         Fehler = ' DecNList : empty list is now "-"'
         RETURN
      ELSEIF (comp.eq.'-') THEN
         CALL qSet (qLis, 1, mL, 1, .false.)
         RETURN  !  added 13-11-90, changed 12-2-91
      ELSEIF (comp(1:1).eq.'/') THEN
         CALL qSet (qLis, 1, mL, 1, .true.)
         qNew = .false.
         ia = 2
      ELSE
         CALL qSet (qLis, 1, mL, 1, .false.)
         qNew = .true.
         ia = 1
         ENDIF

C  Loop : decode ranges (separated by , )
      i = ia+1
 2    CONTINUE
         qEoR = i.gt.nc
         IF (.not.qEoR) qEoR = comp(i:i).eq.','
         IF (qEoR) THEN
            CALL dppIFS (ia, i-1, nc, comp, tComp, iComp,
     *                   1, mL, ki, kf, ks, Fehler)
            IF (Fehler.ne.'&ff') THEN
               Print *, 'DecNList/ "', stri(1:lenU(stri)), '"'
               RETURN
               ENDIF
            CALL qSet (qLis, ki, kf, ks, qNew)
            IF (i.gt.nc) RETURN ! regular exit
            ia = i+1
            ENDIF
         i  = i+1
         GOTO 2

      END ! DecNList

      SUBROUTINE ExpNList (qLis, nLis, J, nJ)
C     ---------------------------------------
            ! JWu 4oct91
         ! expand the list array qLis(nLis) into an index array J(nj),
         ! such that qLis(J(1..nj))=.true.

         IMPLICIT LOGICAL (q)
         DIMENSION         qLis(*), J(*)

         nJ = 0
         DO i = 1, nLis
            IF (qLis(i)) THEN
               nJ    = nJ + 1
               J(nJ) = i
               ENDIF
            ENDDO

         END ! ExpNList

      SUBROUTINE DecJList (stri, MJL, nJLis, JLis, ji, jf, Fehler)
C     ------------------------------------------------------------
            ! JWu 1feb93 indirect, 6jul93 direct
         ! decode stri, set JLis(1..nJLis), values must be in ji..jf

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      CHARACTER*(*)  stri, Fehler
      DIMENSION      JLis(*)

      PARAMETER     (Mcomp=320)
      CHARACTER      comp*320, tComp*20, cl6*6
      DIMENSION      tComp(Mcomp), iComp(Mcomp), rComp(Mcomp)

C  Parse :
      CALL Parser (stri, comp, Mcomp, nc,
     *             tComp, iComp, rComp, 'nm_', Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Insert (Fehler, 1, 'Dec/JList')
         RETURN
         ENDIF

C  Special cases and preset :
      IF     (comp.eq.' ') THEN
         Fehler = ' DecJList : empty list is now "-"'
         RETURN
      ELSEIF (comp.eq.'-') THEN
         nJLis = 0
         RETURN
      ELSEIF (comp(1:1).eq.'/') THEN
         Fehler = 'DecJList/ "/" not allowed'
         RETURN
         ENDIF
      qNew = .true.
      ia = 1
      nJLis = 0

C  Loop : decode ranges (separated by , )
      i = ia+1
 2    CONTINUE
         qEoR = i.gt.nc
         IF (.not.qEoR) qEoR = comp(i:i).eq.','
         IF (qEoR) THEN
            CALL dppIFS (ia, i-1, nc, comp, tComp, iComp,
     *                   ji, jf, ki, kf, ks, Fehler)
            IF (Fehler.ne.'&ff') THEN
               Print *, 'DecJList/ "', stri(1:lenU(Stri)), '"'
               RETURN
               ENDIF
            DO k = ki, kf, ks
               nJLis = nJLis + 1
               IF (nJLis.gt.MJL) THEN
                  CALL Compose2 (Fehler,
     *  'list input/ calling routine expects less than '//
     *  cl6(MJL), ' elements')
                  RETURN
                  ENDIF
               JLis(nJLis) = k
               ENDDO
            IF (i.gt.nc) RETURN ! regular exit
            ia = i+1
            ENDIF
         i  = i+1
         GOTO 2

      END ! DecJList

      SUBROUTINE EncJList (nJLis, JLis, stri)
C     ---------------------------------------
            ! JWu 1feb93 for IDA-main
         ! encode JLis(1..nJLis) into stri

      CHARACTER stri*(*), cl6*6
      DIMENSION JLis(*)

      IF     (nJLis.lt.0) THEN
         stri = '?'
      ELSEIF (nJLis.eq.0) THEN
         stri = '-'
      ELSE
         stri = cl6(JLis(1))
         DO lj = 2, nJLis
            IF (JLis(lj).eq.JLis(lj-1)+1) THEN
               jM = 1
            ELSE
               jM = 0
               ENDIF
            IF (lj.eq.nJLis) THEN
               jP = 0
            ELSE
               IF (JLis(lj+1).eq.JLis(lj)+1) THEN
                  jP = 1
               ELSE
                  jP = 0
                  ENDIF
               ENDIF
            IF     (jM.eq.1 .and. jP.eq.0) THEN
               CALL Append (stri, '-'//cl6(JLis(lj)))
            ELSEIF (jM.eq.1 .and. jP.eq.1) THEN
               ! dont write ( included in '-' )
            ELSE
               CALL Append (stri, ','//cl6(JLis(lj)))
               ENDIF
            ENDDO
         ENDIF

      END ! EncJList

C  ====================================================================
C  L2.9.  Parentheses :
C  ====================================================================

      INTEGER FUNCTION jPos1Kla ( stri, substri, Fehler)
C 331 --------------------------------------------------             JW 13.11.89
      !  Position des ersten Zeichens von substri
      !  beim ersten Auftreten von substri in stri.
      !  Wenn substri nicht in stri enthalten ist,
      !  wird jPos1 = len(stri) + 1 gesetzt.
      !  Klammern werden uebersprungen


         CHARACTER *(*) stri, substri, Fehler
         CHARACTER *1   c
         ls  = len(stri)
         lss = len(substri)
         nKla = 0

         DO 1 i = 1,ls-(lss-1)
            c = stri(i:i)
            IF     (c.eq.'(' .or. c.eq.'[' .or. c.eq.'{') THEN
               nKla = nKla + 1
            ELSEIF (c.eq.')' .or. c.eq.']' .or. c.eq.'}') THEN
               nKla = nKla - 1
               IF (nKla.lt.0) THEN
                  Fehler = '3311 too many ")"'
                  GOTO 2
                  ENDIF
            ELSEIF (nKla.eq.0 .and.
     *              stri(i:i+(lss-1)).eq.substri(1:lss)) THEN
               jPos1Kla = i
               RETURN
               ENDIF
 1          CONTINUE
         IF (nKla.gt.0) THEN
            Fehler = '3312 ")" expected'
            ENDIF
 2       jPos1Kla = len(stri) + 1

         END ! jPos1Kla

      SUBROUTINE TakeVorDelKla ( Stack, Out, Substri, Fehler)
C 332 -------------------------------------------------------        JW 13.11.89
      !  Nimmt alles vor dem ersten Auftreten von Substri
      !  aus Stack und schreibt es nach Out;
      !  entfernt ausserdem Substri aus Stack.
      !  Klammern werden uebersprungen.

         CHARACTER *(*) Stack, Out, Substri, Fehler
         CHARACTER *99 muell
         n     = jPos1Kla (Stack, Substri,Fehler) - 1
            IF (Fehler.ne.'&ff') RETURN
         CALL StaTake (n, Stack, Out)
         CALL StaTake (len(Substri), Stack, muell)

         END ! TakeVorDelKla

      SUBROUTINE TakeVorKla ( Stack, Out, Substri, Fehler)
C 333 ----------------------------------------------------           JW 19. 2.90
      !  Nimmt alles vor dem ersten Auftreten von Substri
      !  aus Stack und schreibt es nach Out.
      !  Klammern werden uebersprungen.

         CHARACTER *(*) Stack, Out, Substri, Fehler
         n     = jPos1Kla (Stack, Substri,Fehler) - 1
            IF (Fehler.ne.'&ff') RETURN
         CALL StaTake (n, Stack, Out)

         END ! TakeVorKla
