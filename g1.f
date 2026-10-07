C  ====================================================================
C
C     Library  WuG :  Graphics
C     Module   G1  :  Commands to specific stations
C
C  ====================================================================

C     J.Wuttke  -  Version Jun91

C     Contents :
C
C        1.1.   Pericom 4014 :   Scrollmode etc.
C                  OverlayOn/Off, GMode, Scrollarea, ClearGraphic/Screen
C
C        1.2.   Tektronix :      Direct graphic driver
C                  TekBytes, ResetOldAdress, TekLinType/GoTo/DrawTo,
C                  TekCharSize/Area/-
C
C        1.3.   Windows :        Routines for both Tek and PS
C                  SetWindow, SetCoord, Coord, PS_Coord, Ticks
C
C        1.4.   TekGraphics :    Graphics via 1.2.
C                  TekFrameClear, TekAxis, Ticks, TekPlotCS,
C                  TekSetSymbol, TekPlotSymbol, TekPoint, TekCurve, TekText,
C                  TekSetDevice
C
C        1.5.   PS_Graphics :    Graphics with postscript
C                  PS_Numbers, PS_PlotCS,
C                  PS_Point, PS_Curve, PS_Text,
C                  OpenPS, ClosePS

C     File Access :
C        uses units 51 : out.tek
C                   52 : out.ps
C                   53 : in.ps
C               ioUnit : terminal or 51 as set by TekSetDevice

C     Major modifications :
C        92/93  made work with Falco and X11 display
C        Jun91  Postscript driver from HP.Schildberg
C        Okt90  Library WuGra; complete revision
C        Okt87  Version Deorie, module DeGra
C        Aug87  Tektronix direct driver from HP.Schildberg
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  ====================================================================
C  g1.1.  Pericom 4014
C  ====================================================================

      SUBROUTINE OverlayOn
C     --------------------
         Print '(1x,5a1,$)',
     *      Char(29),Char(27),Char(92),Char(53),Char(24) ! GS ESC B/SL 5 CAN
         END

      SUBROUTINE OverlayOff
C     ---------------------
         Print '(1x,5a1,$)',
     *      Char(29),Char(27),Char(92),Char(51),Char(24) ! GS ESC B/SL 3 CAN
         END

      SUBROUTINE GMode (iWn)
C     ----------------------
      COMMON / GraPlus / iGraP, iEsc
         DATA  iWo /0/ ! instead of 99=undef'd: 23jun94
               ! renewed 22oct92
            ! iWn = 1 : go to Graph window
            ! iWn = 0 : go to Text window
         IF (iWn.ne.iWo) THEN
            IF     (iWn.eq.0) THEN ! goto Text
               IF     (iEsc.ge.1) THEN
                  Print '(2a1,$)', Char(27), Char(3) ! ESC ETX
               ELSEIF (iGraP.ge.1) THEN
                  Print '(2a1,$)', '+', Char(24)  ! CAN
               ELSE
                  Print '(a1,$)', Char(24)
                  ENDIF
            ELSEIF (iWn.eq.1) THEN
               IF     (iEsc.ge.1) THEN
                  Print '(2a,$)', Char(27),'[?38h'
               ELSEIF (iGraP.ge.1) THEN
                  Print '(2a1,$)', '+', Char(29)  ! GS
               ELSE
                  Print '(a1,$)', Char(2)
                  ENDIF
            ELSE
               CALL Absturz ('GMode', 'iWn o.o.r.')
               ENDIF
            iWo = iWn
            ENDIF

         END ! GMode

      SUBROUTINE ScrollArea (iAnf, iEnd)
C     ----------------------------------

      Print '(1x,a1,a)',         Char(27), '[2J'                    !erase Scree
      Print '(1x,2a1,2(i2,a1))', Char(27), '[', iAnf, ';', iEnd, 'r'!scroll area
      Print '(1x,2a1,i2,a)',     Char(27), '[', iAnf, ';1H'         !goto scroll

      END ! ScrollArea

      SUBROUTINE ClearGraphic
C     -----------------------
         !  clears the graphic and performs some usefull
         ! reset operations on the terminal setup mode.
      COMMON / GraPlus / iGraP, iEsc

      IF (iEsc.ge.1) THEN ! xterm*tek
         print *, 'CLEAR CLEAR CLEEEEEEEAAAAAAAAr??'
         Print '(2a)', Char(27), Char(12) ! ESC LF(10)
      ELSE ! Pericom (HPS)
         Print '(1x,33a1,3(a1,a1,a2),a1)',
     *      Char(29),Char(27),Char(12), ! GS ESC NP
     *     (Char(22),I=1,30),        ! SYN do nothing, let terminal recover
     *      Char(27),Char(92),'E1',  ! select PERICOM 4014 GRAPHICS  ?
     *      Char(27),Char(92),'s0',  ! select LARGE GIN CURSOR       ?
     *      Char(27),Char(92),'t0',  ! select SOLID GIN CURSOR       ?
     *      Char(24)                                      ! CAN
         ENDIF

      END ! ClearGraphic

      SUBROUTINE ClearScreen
C     ----------------------
      Print '(1x,a1,a3)', Char(27),'[2J'    ! erase text
      Print '(1x,a1,a5)', Char(27),'[1;1H'  ! goto line 1
      END

C  ====================================================================
C  g1.2.  Tektronix
C  ====================================================================
C               The routines of this section are due to H.P.Schildberg

C  --------------------------------------------------------------------
C  g1.2.1.  Screen coordinates as optimized strings
C  --------------------------------------------------------------------

      SUBROUTINE TekBytes (IX,IY,OutChar,LEN)
C     ---------------------------------------
         ! Transforms the screen coordinates IX,IY into a string of 5 bytes,
         ! (starting from element OutChar(9)). An outputstring optimized with
         ! respect to length can be found starting from OutChar(1). LEN is the
         ! length of this string. Lateron either String can be sent to the
         ! Tektronics graphic area.
         ! Input:  IX,IY
         ! Output: OutChar, LEN

      IMPLICIT LOGICAL (q)
      CHARACTER*1 OutChar(*), OldChar(16)

      COMMON / GraTerm / qWindow, qGTOverlay, qGToldAdr
      COMMON / oldie   / OldChar

C calculate the characters for the point adress.
C order: highY, lsbYX, lowY, highX, lowX
        OutChar( 9) = Char( mod(IY/128,32)+32 )
        OutChar(10) = Char( mod(IY,4)*4+mod(IX,4)+96 )
        OutChar(11) = Char( mod(IY/4,32)+96 )
        OutChar(12) = Char( mod(IX/128,32)+32 )
        OutChar(13) = Char( mod(IX/4,32)+64 )

C  "abspecken" of output, doesn't work on all terminals :
      IF (qGToldAdr) THEN
        LEN=0
        IF (OldChar(9) .ne. OutChar(9)) THEN
           ! high Y needed
           LEN=1
           OldChar(9)=OutChar(9)
           OutChar(1)=OutChar(9)
           ENDIF
        q10 = OldChar(10) .ne. OutChar(10)
        q11 = OldChar(11) .ne. OutChar(11)
        q12 = OldChar(12) .ne. OutChar(12)
        IF (q10) THEN
           ! LSBYX needed
           LEN=LEN+1
           OldChar(10)=OutChar(10)
           OutChar(LEN)=OutChar(10)
           ENDIF
        IF (q10 .or. q11 .or. q12) THEN
           ! LOW Y needed
            LEN=LEN+1
           OldChar(11)=OutChar(11)
           OutChar(LEN)=OutChar(11)
           ENDIF
        IF (q12) THEN
           ! HIGH X needed
           LEN=LEN+1
           OldChar(12)=OutChar(12)
           OutChar(LEN)=OutChar(12)
           ENDIF
        ! take lowX in any case.
        LEN=LEN+1
        OldChar(13)=OutChar(13)
        OutChar(LEN)=OutChar(13)
C now the first LEN bytes in OutChar contain the optimized outputstring,
C from OutChar(9)  to OutChar(13) one finds the full string.

      ELSE
         DO i = 1, 5
            OutChar(i) = OutChar(8+i)
            ENDDO
         LEN = 5

         ENDIF ! Abspecken or not

      END ! TekBytes

      SUBROUTINE ResetOldAdress
C     ------------------------------
         ! will load the register for the old point adress with hex. FF = 255,
         ! namely a bit pattern, which never appears in the graphic output.
      CHARACTER*1 OldChar(16)
      COMMON / Oldie   / OldChar
      DO i = 1,16
         OldChar(i) = Char(255)
         ENDDO
      END ! ResetOldAdress

C  --------------------------------------------------------------------
C  g1.2.2.  Lines
C  --------------------------------------------------------------------

      SUBROUTINE TekLinPrint (LinTyp)
C     ------------------------------
         ! sets the linetype for the Tektronics graphic
         !     1 : solid
         !     2 : dotted
         !     3 : dot-dash
         !     4 : short dash
         !     5 : long dash
         !     6 : points

      COMMON / GraPlus / iGraP, iEsc
      COMMON / output  / ioUnit

      iTyp = mod(LinTyp-1,6) + 1 ! result is in 1..6 (23jun94)

      IF (LinTyp.ne.6) THEN ! select linetype
         IF (iGraP.ge.1) THEN
            Write(ioUnit,'(4a1,$)') '+',Char(29),Char(27),Char(95+iTyp)
         ELSE
            Write(ioUnit,'(3a1,$)') Char(29),Char(27),Char(95+iTyp)
            ENDIF
      ELSE ! linetype is points
         IF (iGraP.ge.1) THEN
            Write(ioUnit,'(3a1,$)') '+',Char(29),Char(28)
         ELSE
            Write(ioUnit,'(2a1,$)') Char(29),Char(28)
            ENDIF
         ENDIF

      END ! TekLinPrint

      SUBROUTINE TekLinGoTo ( ix, iy)
C     -------------------------------
         ! Go to (ix, iy), and invoke the vector plot module

      CHARACTER    OutChar(16)*1, OutFormat*8, cv2*2
      COMMON / GraPlus / iGraP, iEsc
      COMMON / output  / ioUnit

      ixi = iinside (ix, 0, 4080)
      iyi = iinside (iy, 0, 3060)
      CALL TekBytes(ixi,iyi,OutChar,LEN)
      IF (iGraP.ge.1) THEN
         OutFormat = '('//cv2(2+LEN)//'a1,$)'
         Write(ioUnit,OutFormat) '+', Char(29), (OutChar(I),I=1,LEN)
      ELSE
         OutFormat = '('//cv2(1+LEN)//'a1,$)'
         Write(ioUnit,OutFormat) Char(29), (OutChar(I),I=1,LEN)
         ENDIF

      END ! TekLinGoTo

      SUBROUTINE TekLinDrawTo ( ix, iy)
C     ---------------------------------
         ! Draw a line to (ix, iy).
         ! The vector plot module has to be invoked by a preceeding
         ! call of TekLinGoTo

      CHARACTER    OutChar(16)*1, OutFormat*8, cv2*2
      COMMON / GraPlus / iGraP, iEsc
      COMMON / output  / ioUnit

      ixi = iinside (ix, 0, 4080)
      iyi = iinside (iy, 0, 3060)
      CALL TekBytes(ixi,iyi,OutChar,LEN)
      IF (iGraP.ge.1) THEN
         OutFormat = '('//cv2(1+LEN)//'a1,$)'
         Write(ioUnit,OutFormat) '+', (OutChar(I),I=1,LEN)
      ELSE
         OutFormat = '('//cv2(0+LEN)//'a1,$)'
         Write(ioUnit,OutFormat) (OutChar(I),I=1,LEN)
         ENDIF

      END ! TekLinDrawTo

      SUBROUTINE TekLin ( ixlo, iylo, ixhi, iyhi)
C     -------------------------------------------
         ! Draw a line  from (ixl, iyl) to (ixh, iyh).

      CALL TekLinGoTo   (ixlo, iylo)
      CALL TekLinDrawTo (ixhi, iyhi)

      END ! TekLin

      SUBROUTINE TekErase (ixl,iyl,ixh,iyh)
C     -------------------------------------
         ! performs a selective erase of the graphic area within the
         ! rectangle defined by ixl,iyl,ixh,iyh.
         ! This routine does only work properly, if the CEFTI-Pericom was first
         ! put into Graphics 2 state and has been reset to PERICOM 4014 graphics
         ! after the first graphic operation was performed by the program

      Print '(3a1,$)', '+', Char(27), 'x'
      Print '(3a1,$)', '+', Char(27), Char(2)
      CALL TekLin (ixl, iyl, ixh, iyh)
      Print '(3a1,$)', '+', Char(27), Char(1)
      Print '(3a1,$)', '+', Char(27), '`'

      END ! TekErase

C  --------------------------------------------------------------------
C  g1.2.3.  Character strings
C  --------------------------------------------------------------------

      SUBROUTINE TekCharSize (i)
C     --------------------------
         ! Select a character size in Tektronics graphic area.

      CHARACTER cl6*6
      COMMON / output  / ioUnit

      IF (i.lt.1 .or. i.gt.4)
     *    CALL Absturz ('TekCharSize', 'i[size] = '//cl6(i) )

      Write(ioUnit,'(1x,4a1,$)') Char(29),Char(27),Char(55+i),Char(24)

      END ! TekCharSize

      BLOCK DATA TekCharArea
C     ----------------------
         ! JWu 5jun91, as block data 1apr93
         ! TCAx(i)*TCAy(i) is the area of a character of size no. i
         ! calculated from total area = 4080*3060 :
            ! i=1:  74 Characters, 35 Lines
            ! i=2:  81 Characters, 38 Lines
            ! i=3: 121 Characters, 58 Lines
            ! i=4: 133 Characters, 64 Lines

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      COMMON / CharArea / TCAx(4), TCAy(4)
      DATA     TCAx       / 55.14d0, 50.37d0, 33.72d0, 30.68d0 /,
     *         TCAy       / 87.43d0, 80.53d0, 52.76d0, 47.81d0 /

      END ! TekCharArea

      SUBROUTINE TekChar (ix, iy, String)
C     -----------------------------------
         ! Write String starting from position (ix, iy)
         ! The charactersize is assumed to be set by TekCharSize

      CHARACTER    PosByt(16)*1, Form*20, String*(*), cv2*2
      COMMON / output  / ioUnit

      IF (ix.lt.0 .or. iy.lt.0) RETURN ! bad limits -> simply ignore the text

      CALL TekBytes (ix, iy, PosByt, LPB) !create byte sequence for coordinates
      Form = '(1x,'//cv2(LPB+2)//'a1,a)'
      Write (ioUnit, Form) Char(29), (PosByt(i),i=1,LPB),
     *                                Char(31), String

      END ! TekChar

C  ====================================================================
C  g1.3.   Windows
C  ====================================================================

C  --------------------------------------------------------------------
C  g1.3.1.  Window set up
C  --------------------------------------------------------------------

      SUBROUTINE SetWindow (iW)
C     -------------------------
         !  Set the viewport (ixl..ixh, iyl..iyh),
         !  the size of ticks, numbers, and symbols,
         !  and the size and position of text lines.

         !  Units are pixels (3060*4080) for Tek. No longer called for PS.

      IMPLICIT LOGICAL (q)

      PARAMETER    (MW=9) ! # implemented formats
      INTEGER      iFrame(4,MW), iSize(4,MW), iLabel(5,MW), iText(4,MW)

      COMMON / Viewport/ iTOTx, iTOTy
      COMMON / graph   / ixl, ixh, iyl, iyh,
     *                   iTack, itick, jNumber, jSymbol
      COMMON / InfPos /  infx(200), infy(200), infSiz(200)

C  Frame co-ordinates : ixl,ixh,iyl,iyh
      DATA  iFrame /   590, 2000, 1750, 3030, ! 1: TEK, half window on the left
     *                2590, 4000, 1750, 3030, ! 2: TEK, half window on the right
     *                 800, 2350,  720, 2920, ! 3: TEK, large, Hochformat
     *                 800, 3000, 1000, 2550, ! 4: TEK, large, Querformat
     *                 695, 4000,  940, 3060, ! 5: TEK, scroll, 20 lines
     *                1300, 2700, 2200, 3060, ! 6: TEK, scroll, 8 lines
     *                 500, 1000,  500, 1240, ! 7: PS, Hochformat OBSOLET
     *                 500, 1240,  500, 1000, ! 8: PS, Querformat "
     *                 695, 4080,  135, 3060/ ! 9: TEK, 25 lines (full screen)
C  Lines : length of iTack, itick; size of numbers, symbols
      DATA  iSize /  30, 12, 3, 120,   ! 1
     *               30, 12, 3, 120,   ! 2
     *               40,  0, 1, 300,   ! 3
     *               40,  0, 1, 300,   ! 4
     *               36, 20, 2, 160,   ! 5
     *               36, 20, 2, 160,   ! 6
     *               36, 20, 2,   0,   ! 7
     *               36, 20, 2,   0,   ! 8
     *               54, 18, 2, 175/   ! 9
C  Label : position of x-,y- label; size
      DATA  iLabel /  2000, 1610,  500, 3080, 3, ! 1
     *                4000, 1610, 2500, 3080, 3, ! 2
     *                2350,  600,  800, 2880, 1, ! 3
     *                2000,  750,  800, 2620, 1, ! 4
     *                2700, 2900, 2700, 2780, 2, ! 5  ! z.Zt. im Innern der Gr.
     *                3000, 2880, 3000, 3000, 2, ! 6
     *                 700,  900,  700,  780, 2, ! 7
     *                 700,  900,  700,  780, 2, ! 8
     *                4050, 2960, 4050, 2840, 2/ ! 9  ! z.Zt. im Innern der Gr.
C  Text : x, y-offset, y-incr, size
      DATA  iText /    540, 1540, 48,  4, ! 1
     *                2540, 1540, 48,  4, ! 2
     *                   0,    0,  0,  0, ! 3
     *                   0,    0,  0,  0, ! 4
     *                   0, 3080, 44,  3, ! 5
     *                   0, 3080, 44,  3, ! 6
     *                   0,    0,  0, 10, ! 7
     *                   0,    0,  0, 10, ! 8
     *                   0, 3070, 44,  3/ ! 9

C  Check :
      IF (qiOutside(iW,1,MW)) THEN
          Print *, 'iW = ', iW
          CALL Absturz ('SetWindow', 'iW o.o.r.')
         ENDIF

C  Set the common block / Viewport / :
      IF (iW.le.6) THEN
         ! TEK
         iTOTx = 4080
         tTOTy = 3060
      ELSE
         ! PS-DIN A4
         iTOTx = 2950
         iTOTy = 2200
         ENDIF

C  Set the common block /graph/ :
      ixl    = iFrame (1, iW)
      ixh    = iFrame (2, iW)
      iyl    = iFrame (3, iW)
      iyh    = iFrame (4, iW)

      iTack  = iSize (1, iW)
      itick  = iSize (2, iW)
      jNumber= iSize (3, iW)
      jSymbol= iSize (4, iW)

C  Set the common block /infpos/ :

      infx   (2) = iLabel (1, iW)
      infy   (2) = iLabel (2, iW)
      infSiz (2) = iLabel (5, iW)

      infx   (3) = iLabel (3, iW)
      infy   (3) = iLabel (4, iW)
      infSiz (3) = iLabel (5, iW)

      DO i = 11, 100
         infx   (i) = iText (1, iW)
         infy   (i) = iText (2, iW) - (i-11) * iText (3, iW)
         infSiz (i) = iText (4, iW)
         ENDDO

      END ! SetWindow

C  --------------------------------------------------------------------
C  g1.3.2.  Coordinate system : points
C  --------------------------------------------------------------------

      SUBROUTINE SetCoord (mode)
C     --------------------------
            ! TekCoord JWu 20mar91, Coord/SetCoord 26jan94
         ! set / gcoord / for use in Coord (..)

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      CHARACTER          mode*(*)

      REAL*8             RD(2)

      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)
      COMMON / graph   / ixl, ixh, iyl, iyh, iTack, itick,
     *                   jNumber, jSymbol
      COMMON / gcoord  / R0(2), RC(3,2)

C  Seize of plot RD and position of origin R0 in relative units :
      RD(1) = 0
      RD(2) = 0
      R0(1) = 0
      R0(2) = 0
      DO j = 1, nDim
         RD(1) = RD(1) + GAREL(j)*dabs(dcosd(GANG(j)))
         RD(2) = RD(2) + GAREL(j)*dabs(dsind(GANG(j)))
         R0(1) = R0(1) - dmin1 (0.d0, GAREL(j)*dcosd(GANG(j)))
         R0(2) = R0(2) - dmin1 (0.d0, GAREL(j)*dsind(GANG(j)))
         ENDDO

C  Dito in absolute coordonates :
      IF     (mode.eq.'TEK') THEN
         RD(1) = (ixh - ixl) / RD(1)
         RD(2) = (iyh - iyl) / RD(2)
         R0(1) = ixl + RD(1)*R0(1)
         R0(2) = iyl + RD(2)*R0(2)
      ELSEIF (mode.eq.'PS' ) THEN
         RD(1) = 10000 / RD(1)
         RD(2) = 10000 / RD(2)
         R0(1) = RD(1)*R0(1)
         R0(2) = RD(2)*R0(2)
         ENDIF

C  Transformation of x,y,z :
      DO j = 1, nDim
         IF (LLG(j).eq.0) THEN
            RC(j,2) = GAREL(j) / (GMM(j,2)-GMM(j,1))
         ELSE
            RC(j,2) = GAREL(j) / dlog10(GMM(j,2)/GMM(j,1))
            ENDIF
         RC(j,1) = RD(1) * RC(j,2) * dcosd(GANG(j))
         RC(j,2) = RD(2) * RC(j,2) * dsind(GANG(j))
         ENDDO

      END ! SetCoord

      SUBROUTINE Coord (x, y, z, ix, iy)
C     ----------------------------------
            ! TekCoord JWu 20mar91, Coord/SetCoord 26jan94
         ! transform real co-ordinates (x,z,y) -> graphic point (ix,iy)

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      REAL*8             Pt(3)
      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)
      COMMON / gcoord  / R0(2), RC(3,2)

      Pt(1) = x
      Pt(2) = y
      Pt(3) = z

C  Origin of Coordonate system :
      RX = R0(1)
      RY = R0(2)

C  Transformation of x,y,z :
      DO j = 1, nDim
         IF (LLG(j).eq.0) THEN
            RX = RX + RC(j,1) * (Pt(j)-GMM(j,1))
            RY = RY + RC(j,2) * (Pt(j)-GMM(j,1))
         ELSE
            IF (Pt(j).gt.0.) THEN ! zweiteinfachste Absturzsicherung
               RX = RX + RC(j,1) * dlog10(Pt(j)/GMM(j,1))
               RY = RY + RC(j,2) * dlog10(Pt(j)/GMM(j,1))
               ENDIF
            ENDIF
         ENDDO

      ix = idnint(RX)
      iy = idnint(RY)

      END ! Coord

      SUBROUTINE PS_Coord (x, y, z, rx, ry)
C     -------------------------------------
            ! JWu 11/13jun91; reduction -> Coord 26jan94

      IMPLICIT REAL*8    (a-h,o-p,r-z)

      CALL Coord (x, y, z, ix, iy)

      rx = dble(ix) / 1000          ! -> units of 0.1axis
      ry = dble(iy) / 1000

      END ! PS_Coord

C  --------------------------------------------------------------------
C  g1.3.3.  Frame : ticks
C  --------------------------------------------------------------------

      SUBROUTINE Ticks (lilo, rmin, rmax, rTick, MTick, nTick, rTack,
     *                  MTack, nTack, rTLim, nTpT, Fehler)
C     ----------------------------------------------------------------
         !  Calculates the positions of small and large ticks
         !  (ticks and tacks) to be drawn on the co-ordinate axes.
         !  Now both linear and logarithmic axes are possible (lilo=0/1).
         !  All co-ordinates are real co-ordinates.
         !  the routine returnes rTick(1..nTick), rTack(1..nTack),
         !  rTLim(2), and nTpT (ticks per tack).

         !  JWu 15mar91, adapted from H.P.Schildberg's new version
         !  "ticks_array" in HPSGRA3.FOR.

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)
      DIMENSION          rTick(MTick), rTack(MTack), rTLim(2)
      CHARACTER          Fehler*(*)

         !  Int*4 <-> Real*8 conversion : dFlotJ, idNint, idInt, iPint.

C  Check limits :
      IF (rmin.ge.rmax) THEN
         Fehler = 'PROGRAM ERROR/ Ticks/ Bad Limits'
         RETURN
         ENDIF

C  Initialize tick counters :
      nTick = 0
      nTack = 0

      IF (lilo.eq.0) THEN
C  Linear scale :

C  - R = logarithm to base 10 of plot range :
         R  = dlog10 (rmax-rmin)
         IR = idint (R)
         IF (R.lt.0.d0) IR = idint (R-1.d0)

C  - RD = fractional part of R :
         RD = R - dble(IR)

C  - Calculate TL = spacing between large ticks : (revised JWu 2jul92,12oct93)
         RDR = 10.**RD
         IF     (RDR.gt.10.-1.d-5) THEN ! allowing for eps=1d-5 (12oct93)
            TL = 2.5* 10.0**IR
            nTpT = 5
         ELSEIF (RDR.gt.5+1.d-5) THEN
            TL = 2. * 10.0**IR
            nTpT = 4 ! Hinweis von Hanne
         ELSEIF (RDR.gt.2.5+1.d-5) THEN
            TL = 1. * 10.0**IR
            nTpT = 5
         ELSEIF (RDR.gt.1.6+1.d-5) THEN
            TL = .5 * 10.0**IR
            nTpT = 5
         ELSEIF (RDR.gt.1.25+1.d-5) THEN
            TL = .4 * 10.0**IR
            nTpT = 4
         ELSE
            TL = .25* 10.0**IR
            nTpT = 5
            ENDIF
C - TL0 = startposition for first TL :
         r0 = TL * iPint (rmin/TL)  ! preceeding int, -> WuL1
         IF ( dAbs(rmin-r0) .lt. 1.d-6 *TL) THEN
            TL0 = r0 !  Should equal rmin
         ELSE
            TL0 = r0 + TL
            ENDIF
C - Number of large ticks, set array :
         nTack = 1 + idInt( (rmax-TL0)/TL + .1 )
         IF (nTack.gt.MTack) THEN
            Print *, ' nTack, Mtack : ', nTack, MTack
            Fehler = 'PROGRAM ERROR/ Ticks/ not enough tacks foreseen'
            RETURN
            ENDIF
         DO i = 1, nTack
            rTack(i) = TL0 + (i-1)*TL
            ENDDO
         IF (rTack(nTack).gt.rmax+1.d-6*TL) nTack = nTack-1 ! reverse the + .
         rTLim(1) = TL0 - TL
         rTLim(2) = rTack(nTack) + TL

C - The same for the small ticks :
         TS = TL / nTpT ! usually 5, sometimes 4 (16nov93)
         TS0 = TL0 - TS * idInt( (TL0-rmin)/TS )
         nTick = 1 + idInt( (rmax-TS0)/TS )
         IF (nTick.gt.MTick) THEN
            Print *, ' nTick, Mtick : ', nTick, MTick
            Fehler = 'PROGRAM ERROR/ Ticks/ not enough ticks foreseen'
            RETURN
            ENDIF
C - - For simplicity, at each large tick we will also have a small tick :
         DO i = 1, nTick
            rTick(i) = TS0 + (i-1)*TS
            ENDDO

      ELSE
C  Logarithmic scale :

C  - Check limits :
         IF (rmin.le.0.) THEN
            Fehler = 'PROGRAM ERROR/ Ticks/ negative Log'
            RETURN
            ENDIF

C  - Determine smallest and largest exponent for large ticks :
         rlgmin = dlog10 (rmin)
         rlgmax = dlog10 (rmax)

         rlgrel  = rlgmax - rlgmin
         IF (rlgrel.gt.40.) THEN
            Fehler = 'Ticks/ Range exceeded in LOG_TICKS'
            RETURN
            ENDIF

         IF (dabs(dmod(rlgmin,1.d0)) .lt. 1.d-6) THEN
            minexp = idNint (rlgmin)
         ELSE
            minexp = iPint (rlgmin) + 1
            ENDIF
         IF (dabs(dmod(rlgmax,1.d0)) .lt. 1.d-6) THEN
            maxexp = idNint (rlgmax)
         ELSE
            maxexp = iPint (rlgmax)
            ENDIF

C  - Increment = number of decades per large tick :
         Increment = 1 + rlgrel/7
         IF (Increment.gt.1) THEN
            minexp = minexp - mod(minexp,Increment)
            maxexp = maxexp - mod(maxexp,Increment)
            ENDIF

C  - Set large ticks array :
         eins = 1. + 1.d-6 ! Unity plus arithmetic tolerance
         nTack = 0
         DO i = minexp, maxexp, Increment
            r0 = 1.d1 ** i
            IF ( qrinside(r0,rmin/eins,rmax*eins) ) THEN
               nTack = nTack + 1
               rTack(nTack) = r0
               ENDIF
            ENDDO
         IF (nTack.ge.1) THEN
            rTLim(1) = rTack(1)     / 1.d1 ** Increment
            rTLim(2) = rTack(nTack) * 1.d1 ** Increment
         ELSE
            ! some lines got lost, Grenoble, may95
            ENDIF

C  - Set small ticks array (rewritten 13jan93) :
         IF     (Increment.eq.1 .and. rlgrel.le.6.) THEN
            nTpT = 9
            DO i = minexp-1,maxexp+1
               DO j = 2, 9
                  r0 = 10.d0**i * j
                  IF ( qrinside(r0,rmin,rmax) ) THEN
                     nTick = nTick + 1
                     rTick(nTick) = r0
                     ENDIF
                  ENDDO
               ENDDO
         ELSEIF (Increment.eq.1 .and. rlgrel.le.12.) THEN
            nTpT = 3
            ! 1-2-5-10 - Schritte : haesslich ?
            DO i = minexp-Increment,maxexp+Increment
               DO j = 1,3
                  IF (j.eq.1) r0 = 10.d0**i * 1
                  IF (j.eq.2) r0 = 10.d0**i * 2
                  IF (j.eq.3) r0 = 10.d0**i * 5
                  IF ( qrinside(r0,rmin,rmax) ) THEN
                     nTick = nTick + 1
                     rTick(nTick) = r0
                     ENDIF
                  ENDDO
               ENDDO
         ELSE
            IF (Increment.eq.1) THEN
               nTpT = 1
            ELSE
               nTpT = -Increment ! means 1/Increment means 1 tick per decade
               ENDIF
            DO i = minexp-Increment,maxexp+Increment
               r0 = 10.d0**i * j
               IF ( qrinside(r0,rmin,rmax) ) THEN
                  nTick = nTick + 1
                  rTick(nTick) = r0
                  ENDIF
               ENDDO
            ENDIF

         ENDIF ! lilo

      END ! Ticks

C  ====================================================================
C  g1.4.   TekGraphics
C  ====================================================================

C  --------------------------------------------------------------------
C  g1.4.1.  Frame
C  --------------------------------------------------------------------

      SUBROUTINE TekFrameClear
C     ------------------------
         !  clear the box

      COMMON / graph   / ixl, ixh, iyl, iyh, iTack, itick,
     *                   jNumber, jSymbol

      i = 1
      CALL TekErase (ixl+i, iyl+i, ixh-i, iyh-i)

      END ! TekFrameClear

      SUBROUTINE TekAxis (j, x0, y0, z0, qTic, qLab, Fehler)
C     ------------------------------------------------------
            !  TekFrameTicks, renewed JWu 15mar91
            !  rewritten 27jan94
         !  Draw an axis from (x0,y0,z0) in direction j.
         !  If qTic, draw ticks and tacks on axis.
         !  If qLab, write labels on tacks.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8    (a-h,o-p,r-z)
      PARAMETER         (MTick=100, MTack=10)
      DIMENSION          rTick(MTick), rTack(MTack),
     *                   DummyL(2), Pi(3), Pf(3)
      CHARACTER          Fehler*(*), String*20

      EXTERNAL TekCharArea

      COMMON / graph   / ixl, ixh, iyl, iyh, iTack, itick, jNumber,
     *       jSymbol
      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)
      COMMON / CharArea / TCAx(4), TCAy(4)

      IF (Fehler.ne.'&ff') RETURN

C  Starting point :
      Pi(1) = x0
      Pi(2) = y0
      Pi(3) = z0
      DO jj = 1, nDim
         Pf(jj) = Pi(jj)
         ENDDO

C  This axis is no. j :
      lilo = LLG(j)
      gmi  = GMM(j,1)
      gma  = GMM(j,2)

C  Draw axis :
      Pf(j) = GMM(j,2)
      CALL SetCoord ('TEK')
      CALL TekLinPrint (1)
      CALL Coord (Pi(1), Pi(2), Pi(3), jxl, jyl)
      CALL Coord (Pf(1), Pf(2), Pf(3), jxh, jyh)
      CALL TekLin (jxl, jyl, jxh, jyh)

      IF (.not.qTic) RETURN ! done.

C  Angle of axis in drawing plane :
      jxd = jxh - jxl
      jyd = jyh - jyl
      IF (jxd.eq.0 .and. jyd.eq.0) RETURN ! should not occur
      angax = datan2d (dble(jyd), dble(jxd))
c      write (9,*) ' j x y ang ', j, angax, jxd, jyd

C  Direction of ticks :
      IF (angax.lt.-135 .or. angax.gt.45) THEN
         utix = -dsind(angax)
         utiy =  dcosd(angax)
         irl  = -1 ! ticks are on left of axis
      ELSE
         utix =  dsind(angax)
         utiy = -dcosd(angax)
         irl  = +1 ! ticks are on right of axis
         ENDIF

C  Get positions of ticks and tacks :
      IF (gmi.ge.gma) THEN
         Fehler = 'PROGRAM ERROR/ TekAxis/ GMM(1,1)>=GMM(1,2)'
         RETURN
         ENDIF
      CALL Ticks (lilo, gmi, gma, rTick, MTick, nTick, rTack,
     *            MTack, nTack, DummyL, nTpT, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Draw small ticks ("ticks") :
      jtix = idnint(utix*itick)
      jtiy = idnint(utiy*itick)
      DO i = 1, nTick
         Pf(j) = rTick(i)
         CALL Coord (Pf(1), Pf(2), Pf(3), jx, jy)
         CALL TekLin (jx, jy, jx+jtix, jy+jtiy)
         ENDDO

C  Draw large ticks ("tacks") :
      jtix = idnint(utix*iTack)
      jtiy = idnint(utiy*iTack)
      DO i = 1, nTack
         Pf(j) = rTack(i)
         CALL Coord (Pf(1), Pf(2), Pf(3), jx, jy)
         CALL TekLin (jx, jy, jx+jtix, jy+jtiy)
         ENDDO

      IF (.not.qLab) RETURN ! done.

C  Set character size :
      CALL TekCharSize (jNumber)
      rx = TCAx(jNumber)
      ry = TCAy(jNumber)

      jtix = idnint(utix*rx*0.7)
      jtiy = idnint(utiy*ry*1.4)

C  Draw Labels :
      DO i = 1, nTack
         Pf(j) = rTack(i)
         CALL Coord (Pf(1), Pf(2), Pf(3), jx, jy)
         CALL NiceNum(rTack(i), String, nNum)
         jlabx = idnint (nNum * rx)   ! size of label
         jlaby = idnint (ry)
         xl = jx+1.35*jtix-.5*jlabx        ! x-position for centered label
         yl = jy+1.*jtiy-.33*jlaby        ! y-position
         IF (dabs(dsind(angax)).gt.0.03) THEN ! shift x-position
            xfree = 0. !dabs((dabs(1.2d0*jtiy)+.5*jlaby)/dsind(angax)) + rx*.7
            xl = xl + dsign(1.d0,utix) * dmax1(0.d0,jlabx/2.-xfree)
            ENDIF
         CALL TekChar(iinside(idnint(xl),0,4080-jlabx),
     *                iinside(idnint(yl),0,3060-jlaby),
     *                String(1:nNum))
         ENDDO

      END ! TekAxis

      SUBROUTINE TekPlotCS (LL, GG, AxCro, AxAng, AxLen, ndi,
     *                      qBox, Fehler)
C     ---------------------------------------------------------------------
         ! Plot coordinate system (axes j=1,2,3)
         ! LL(j)    : scale (=1,2 : lin,log)
         ! GG(j,k)  : co-ordinates of edges (k=1,2 : min,max)

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)
      DIMENSION          LL(3), GG(3,2), AxCro(3), AxAng(2), AxLen(3)
      CHARACTER          Fehler*(*)

      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)
      COMMON / graph   / ixl, ixh, iyl, iyh, iTack, itick,
     *                   jNumber, jSymbol

C  Set COMMON variables :
      nDim = ndi
      DO j = 1,nDim
         IF (GG(j,1).ge.GG(j,2)) CALL Absturz ('TekPlotC', 'bad range')
         IF (LL(j).eq.1 .and. GG(j,1).le.0.) CALL Absturz ('TekPlotCS',
     *      'logarithmic scale incompatibel with negative lower limit')
         IF (.not.qBox .and. qroutside(AxCro(j), GG(j,1), GG(j,2)))
     * CALL Absturz ('TekPlotCS',
     *               'axes cross outside coordinate range')
         IF (AxLen(j).le.0)
     * CALL Absturz ('TekPlotCS', ' axis length not > 0')
         GMM(j,1) = GG(j,1)
         GMM(j,2) = GG(j,2)
         LLG(j)   = LL(j)
         GAREL(j) = AxLen(j)
         ENDDO
      IF (nDim.ge.3) THEN
         GANG(1) = AxAng(1)
         GANG(2) = 90
         GANG(3) = AxAng(2)
      ELSE
         GANG(1) = 0
         GANG(2) = 90
         ENDIF

C  Begin to plot :
      CALL ResetOldAdress
      CALL SetCoord ('TEK')

      IF (nDim.eq.2) THEN
         IF (qBox) THEN
            CALL TekAxis
     * (1, GG(1,1), GG(2,1), 0.d0, .true., .true., Fehler)
            CALL TekAxis
     * (1, GG(1,1), GG(2,2), 0.d0, .true., .false.,Fehler)
            CALL TekAxis
     * (2, GG(1,1), GG(2,1), 0.d0, .true., .true., Fehler)
            CALL TekAxis
     * (2, GG(1,2), GG(2,1), 0.d0, .true., .false.,Fehler)
         ELSE
            CALL TekAxis
     * (1, GG(1,2), AxCro(2), 0.d0, .true., .true., Fehler)
            CALL TekAxis
     * (2, AxCro(1), GG(2,1), 0.d0, .true., .true., Fehler)
            ENDIF
      ELSE
         CALL TekAxis
     * (1, GG(1,1), AxCro(2), AxCro(3), .true., .true., Fehler)
         CALL TekAxis
     * (2, AxCro(1), GG(2,1), AxCro(3), .true., .true., Fehler)
         CALL TekAxis
     * (3, AxCro(1), AxCro(2), GG(3,1), .true., .true., Fehler)
         ENDIF

      END ! TekPlotCS

C  --------------------------------------------------------------------
C  g1.4.2.  Data points
C  --------------------------------------------------------------------

      SUBROUTINE TekSetSymbol (iSymbol, rSize)
C     ----------------------------------------                       JWu 22nov90
         ! set relative co-ordinates to draw
         ! symbol iSymbol with magnification rSize/100.
         ! The symbol consists of nCurv curves, each curve
         ! consisting of nLin straight lines.

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      DIMENSION        nLines(11), nCurves(11)
      REAL*4           Offsets(12,11), Radius(11)   ! strictly local

      COMMON / TekSymbol / nCurv, nLin, iOff(12)    !  = Output

      DATA  nCurves / 1,1,1,1,1,1,2,2,3,1,1 /
     *      nLines  / 4,4,4,4,3,3,1,1,1,3,3 /
     *      Offsets /
     *  -1., 1.,  1., 1.,  1.,-1., -1.,-1., -1., 1.,  0., 0.,   ! square
     *   0., 1.,  1., 0.,  0.,-1., -1., 0.,  0., 1.,  0., 0.,   ! karo
     *  -1., 1.,  1.,-1., -1.,-1.,  1., 1., -1., 1.,  0., 0.,   ! eieruhr
     *  -1., 1.,  1.,-1.,  1., 1., -1.,-1., -1., 1.,  0., 0.,   ! valve
     *  -7.,-4.,  0., 8.,  7.,-4., -7.,-4.,  0., 0.,  0., 0.,   ! 3angle
     *   7., 4.,  0.,-8., -7., 4.,  7., 4.,  0., 0.,  0., 0.,   ! cedez passage
     *  -1., 0.,  1., 0.,  0.,-1.,  0., 1.,  0., 0.,  0., 0.,   ! +
     *  -1., 1.,  1.,-1., -1.,-1.,  1., 1.,  0., 0.,  0., 0.,   ! x
     *  -5., 0.,  5., 0., -4.,-3.,  4., 3., -4., 3.,  4.,-3.,   ! *
     *  -4.,-7.,  8., 0., -4., 7., -4.,-7.,  0., 0.,  0., 0.,   ! |>
     *   4.,-7., -8., 0.,  4., 7.,  4.,-7.,  0., 0.,  0., 0./   ! <|
     *      Radius /
     *   .08, .12, .09, .09, .016, .016, .1, .08, .02, .016, .016 /

      CALL TekLinPrint(1)   !symbols always as solid lines

      isy = mod (iSymbol-1, 11) + 1  ! Provisorischer Standort
      nLin = nLines (isy)
      nCurv= nCurves(isy)
      DO j = 1, 2*nCurv*(nLin+1)
         iOff(j) = nint(Offsets(j,isy)*Radius(isy)*rSize)
         ENDDO

      END ! TekSetSymbol

      SUBROUTINE TekPlotSymbol (ixc, iyc)
C     -----------------------------------
         ! Plots a symbol at position ixc,iyc.
         ! The symbol has to be selected by a preceeding
         ! call of TekSetSymbol
                  ! GoTo/DrawTo: JWu,28nov90

      COMMON / TekSymbol / nCurv, nLin, iOff(12)

      DO ic = 0, (nCurv-1) * 2*(nLin+1), 2*(nLin+1)
            CALL TekLinGoTo   ( ixc+iOff(ic+1), iyc+iOff(ic+2) )
         DO il = ic + 2, ic + nLin*2, 2
            CALL TekLinDrawTo ( ixc+iOff(il+1), iyc+iOff(il+2) )
            ENDDO
         ENDDO

      END ! TekPlotSymbol

      SUBROUTINE TekPoint (X, Y, D, n, z, iSymb, rSyMag, qErr)
C     --------------------------------------------------------
         !  Plot Y(i) vs. X(i) as data points
         !  Eventually error bars D(i) are added.

               ! J.Wu  7.11.90 : X has no longer to be in ascending order
               !       6. 3.91 : line drawing renewed
               !      20. 3.91 : logarithmic scale
               !      11.10.91 : TekPaint -> GraPaint, TekPoint, TekCurve

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)

      DIMENSION          X(*),Y(*), D(*)

      COMMON / output  / ioUnit
      COMMON / graph   / ixl, ixh, iyl, iyh, iTack,
     *                   itick, jNumber, jSymbol
      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)

      IF (iSymb.le.0) CALL Absturz ('TekPoint', 'iSymb <= 0')

      CALL TekSetSymbol (iSymb, jSymbol*rSyMag)
      CALL SetCoord('TEK')

      DO i = 1, n
         IF (qrinside(X(i), GMM(1,1), GMM(1,2)) .and.
     *       qrinside(Y(i), GMM(2,1), GMM(2,2))      ) THEN

            CALL Coord (X(i), Y(i), z, ixc, iyc)
            CALL TekPlotSymbol (ixc, iyc)
            IF (qErr) THEN ! error bar, restricted to graph range
               RYU   = dinside ( Y(i)+D(i), GMM(2,1), GMM(2,2) )
               RYD   = dinside ( Y(i)-D(i), GMM(2,1), GMM(2,2) )
               CALL Coord (X(i), RYU, z, ixc, iycU)
               CALL Coord (X(i), RYD, z, ixc, iycD)
               CALL TekLin ( ixc,iycD, ixc,iycU )
               ENDIF ! qErr
            ENDIF ! point inside
         ENDDO ! i

      END ! TekPoint

      SUBROUTINE TekCurve (X, Y, n, z, iLine, rSyMag)
C     -----------------------------------------------
            ! JWu 11oct91 separate subroutine.
         ! Plot Y(i) vs. X(i) as a polygon.
         ! Lines intersecting the graphic frame must be
         ! treated on entry.

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)

      DIMENSION          X(*), Y(*)

      IF (iLine.le.0) CALL Absturz ('TekCurve', 'iLine <= 0')

      CALL TekLinPrint(iLine) ! rSyMag nicht waehlbar
      CALL SetCoord ('TEK')

      CALL Coord (X(1), Y(1), z, ix, iy)
      CALL TekLinGoTo (ix, iy)
      DO i = 2, n
         CALL Coord (X(i), Y(i), z, ix, iy)
         CALL TekLinDrawTo (ix, iy)
         ENDDO

      END ! TekCurve

C  --------------------------------------------------------------------
C  g1.4.3.  Text
C  --------------------------------------------------------------------

      SUBROUTINE TekText (iNr, Text)
C     ------------------------------
         ! write Text at position iNr

      CHARACTER  Text*(*), aus*256

      COMMON / Viewport / iTOTx, iTOTy
      COMMON / InfPos /  infx(200), infy(200), infSiz(200)

      IF (iNr.lt.1 .or. iNr.gt.200) RETURN  ! (9mar93) kein Grund fuer Absturz

      isize = infSiz(iNr)
      IF (isize.eq.0) RETURN

      la = max0 (1, min0 (lenU(Text), 256)) ! continue in case Text=' ' (?)
      aus = Text

C  Auxiliary, for macros :
      IF     (isize.eq.1) THEN
         width =  55.1
         iySh  =  60
         jSy   = 500
      ELSEIF (isize.eq.2) THEN
         width =  50.4
         iySh  =  40
         jSy   = 360
      ELSEIF (isize.eq.3) THEN
         width =  33.7
         iySh  =  20
         jSy   = 150
      ELSEIF (isize.eq.4) THEN
         width =  30.7
         iySh  =  12
         jSy   = 100
         ENDIF

C  Determine position of first character :
         ! JWu 3jul91
      ix = infx(iNr)
      iy = infy(iNr)

      IF     (aus(1:6).eq.'&cent ') THEN
         ! middle centered
         CALL DelVonBis (aus, 1, 6)
         zeile = width * lenU(aus)
         ix = max0 (ix - nint(zeile/2), 0)
      ELSEIF (aus(1:7).eq.'&right ') THEN
         ! right centered
         CALL DelVonBis (aus, 1, 7)
         zeile = width * lenU(aus)
         ix = max0 (ix - nint(zeile), 0)
      ELSE
         ! left centered
         ! force into graphic range
         zeile = width * lenU(aus)
         ix = min0 (ix, iTOTx-nint(zeile))
         ENDIF

C  Search macros :
      DO i = 1, la
         IF (aus(i:i).eq.'&') THEN
            IF     (aus(i+1:i+3).eq.'sy=') THEN
               CALL Fi1I (aus(i+4:i+5), isy)
               IF (aus(i+4:i+4).ne.'#') THEN
                  CALL Gong (13)
                  Print *, 'ERROR/ macro &sy in TekText'
                  RETURN
                  ENDIF
               CALL DelVonBis (aus, i, i+4)
               CALL Insert (aus, i, '  ')
               CALL TekSetSymbol (isy, jSy*1.d0)
               ixs = ix + nint (width * (i + 0.5))
               CALL TekPlotSymbol (ixs, iy+iySh)
            ELSEIF (aus(i+1:i+3).eq.'li=') THEN
               CALL Fi1I (aus(i+4:i+5), ili)
               IF (aus(i+4:i+4).ne.'#') THEN
                  CALL Gong (13)
                  Print *, 'ERROR/ macro &li in TekText'
                  RETURN
                  ENDIF
               CALL DelVonBis (aus, i, i+4)
               CALL Insert (aus, i, '    ')
               CALL TekLinPrint (ili)
               ix1 = ix + nint (width * (i + 0.))
               ix2 = ix + nint (width * (i + 3.))
               CALL TekLin (ix1, iy+iySh, ix2, iy+iySh)
               ENDIF
            ENDIF ! '&'
         ENDDO ! i

      CALL TekCharSize (isize)
      CALL TekChar (ix, iy, aus(1:la))

      END ! TekText

C  --------------------------------------------------------------------
C  g1.4.4.  Master : set unit, laserfile
C  --------------------------------------------------------------------

      SUBROUTINE TekSetDevice (mm)
C     ----------------------------
          !  Set the unit mm the next output will be sent to

      COMMON / Output  / ioUnit

      IF (mm.eq.0) CALL ClearGraphic
      ioUnit = mm
      IF (ioUnit.ne.0) CALL ResetOldAdress

      END ! TekSetDevice

C  ====================================================================
C  g1.5.   PS: PostScript driver
C  ====================================================================

C  --------------------------------------------------------------------
C  g1.5.1. Frame
C  --------------------------------------------------------------------

      SUBROUTINE PS_Numbers (ixy, XY, nXY)
C     ------------------------------------
            ! HPS. JWu 11jun91, 4dec91.
         ! Large ticks and numbers attached to the axes:
         ! PostScript line : 0.200 (0.5) yN

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      DIMENSION         XY(*)
      CHARACTER         num*20, pref*9

      CALL SetCoord('PS')
      z = 0
      Write (52, '(a)') '['
      DO i = 1,nXY
         IF     (ixy.eq.1) THEN
            CALL PS_Coord (XY(i), 1.d0, z, r, dummy)
         ELSEIF (ixy.eq.2) THEN
            CALL PS_Coord (1.d0, XY(i), z, dummy, r)
         ELSE
            CALL Absturz ('PS_Numbers', 'ixy o.o.r.')
            ENDIF
         CALL NiceNum (XY(i), num, ih)
         Write (pref,'(f9.5)') r
         CALL PS_Text (num(1:ih), pref, ' ')
         ENDDO
      Write (52, '(a)') '] SetTacVec'

      END ! PS_Numbers

      SUBROUTINE PS_PlotCS (LL, GG, labelX, labelY, Fehler)
C     -----------------------------------------------------
            ! JWu 11jun91
         ! Plot co-ordinate system (axes j=1,2,3)
         ! LL(j)  = 0 : linear scale
         !          1 : logarithmic scale
         ! GG(j,k)    : co-ordinates of edges (k=1,2 : min,max)

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)

      PARAMETER    (MTick=100, MTack=10)
      DIMENSION     rTick(MTick), rTack(MTack), rTLim(2),
     *              GG(3,2), LL(3)
      CHARACTER*(*) Fehler, labelX, labelY
      CHARACTER     cllx*3, clly*3

      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)

      CALL SetCoord('PS')
      z = 0

C  Write to file :
      Write (52, '(a/)') '%: Coordinate system : '
      IF (LLG(1).eq.0) THEN
         cllx = 'Lin'
      ELSE
         cllx = 'Log'
         ENDIF
      write (52, '(3a,g14.7,a,g14.7)')
     *   '%  ', cllx, ' x-axis from ', GMM(1,1), ' ', GMM(1,2)
      IF (LLG(2).eq.0) THEN
         clly = 'Lin'
      ELSE
         clly = 'Log'
         ENDIF
      write (52, '(3a,g14.7,a,g14.7)')
     *   '%  ', clly, ' y-axis from ', GMM(2,1), ' ', GMM(2,2)

      Write (52, '(/a/)') 'Resets'

      CALL Ticks (LLG(1), GMM(1,1), GMM(1,2),
     * rTick, MTick, nTick, rTack, MTack, nTack, rTLim, nTpT, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL PS_Coord (rTLim(1), 1.d0, z, rTx1, dummy)
      CALL PS_Coord (rTLim(2), 1.d0, z, rTxn, dummy)
      Write (52, '(2(f9.5,1x),i2,1x,i2,2a)')
     *   rTx1, rTxn, nTack+2, nTpT, ' SetTicVec', cllx
      CALL PS_Numbers(1, rTack, nTack)
      Write (52, '(a)')
     * '0 10   0  0     0  90 OneAxx Axx Tic xTacL xNumL % low x axis'
      Write (52, '(a)')
     * '0 10   0 10     0 270 OneAxx Axx Tic xTacH       % top x axis'

      Write (52, '(1x)')
      CALL Ticks (LLG(2), GMM(2,1), GMM(2,2),
     * rTick, MTick, nTick, rTack, MTack, nTack, rTLim, nTpT, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL PS_Coord (1.d0, rTLim(1), z, dummy, rTy1)
      CALL PS_Coord (1.d0, rTLim(2), z, dummy, rTyn)
      Write (52, '(2(f9.5,1x),i2,1x,i2,2a)')
     *   rTy1, rTyn, nTack+2, nTpT, ' SetTicVec', clly
      CALL PS_Numbers(2, rTack, nTack)
      Write (52, '(a)')
     * '0 10   0  0    90   0 OneAxx Axx Tic yTacL yNumL % left y axis'
      Write (52, '(a)')
     * '0 10  10  0    90 180 OneAxx Axx Tic '//
     * 'yTacH % yNumH % right y axis'

      Write (52, '(1x)')
      CALL PS_Text (labelX, ' ', 'xCL')
      CALL PS_Text (labelY, ' ', 'yCL')

      Write (52, '(1x)')

      END ! PS_PlotCS

C  --------------------------------------------------------------------
C  g1.5.2. Data
C  --------------------------------------------------------------------

      SUBROUTINE PS_Point (X, Y, D, n, z, iSymb, rSyMag, qErr)
C     --------------------------------------------------------
            ! HPS; JWu 11jun91, 11oct91.
         !  Plot Y(i) vs. X(i), eventually error bars D(i) are added.

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)
      DIMENSION           X(*),Y(*), D(*)
      CHARACTER           h1*20, h2*20, text*80, cl2*2

      COMMON / output  / ioUnit
      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)

      IF (iSymb.le.0) CALL Absturz ('PS_Point', 'iSymb .le. 0')
      CALL SetCoord('PS')

      Write (52, '(a,g12.4)') '%: one spectrum / z = ', z
      Write (52, '(i2,a)')     iSymb, ' pstyle'
      DO i = 1, n
         IF (qrinside(X(i), GMM(1,1), GMM(1,2)) .and.
     *       qrinside(Y(i), GMM(2,1), GMM(2,2))      ) THEN
c            RYU   = dinside ( Y(i)+D(i), GMM(2,1), GMM(2,2) )
c            RYD   = dinside ( Y(i)-D(i), GMM(2,1), GMM(2,2) )
            CALL PS_Coord (X(i), Y(i),      z, xc, yc)
            CALL PS_Coord (X(i), Y(i)+D(i), z, xc, ycU)
            CALL PS_Coord (X(i), Y(i)-D(i), z, xc, ycD)
            dc = dmin1 ((ycU-ycD)/2, 99.999d0) ! prevent overflow
            IF     (i.eq.1) THEN
               Write (52, '(3(1x,f7.3),a)') xc, yc, dc, ' ti'
            ELSEIF (i.eq.n) THEN
               Write (52, '(3(1x,f7.3),a/)') xc, yc, dc, ' tf'
            ELSE
               Write (52, '(3(1x,f7.3),a)') xc, yc, dc, ' t'
               ENDIF
            ENDIF ! point
         ENDDO ! i

      END ! PS_Point


      SUBROUTINE PS_PointLog (X, Y, D, n, z, iSymb, rSyMag, qErr)
C     --------------------------------------------------------
         ! FK added due to problems with correct error bars
         ! for logarithmic plot similar to PS_Point aug05

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)
      DIMENSION           X(*),Y(*), D(*)
      CHARACTER           h1*20, h2*20, text*80, cl2*2

      COMMON / output  / ioUnit
      COMMON / glimit  / GMM(3,2), LLG(3), nDim, GANG(3), GAREL(3)

      IF (iSymb.le.0) CALL Absturz ('PS_Point', 'iSymb .le. 0')
      CALL SetCoord('PS')

      Write (52, '(a,g12.4)') '%: one spectrum / z = ', z
      Write (52, '(i2,a)')     iSymb, ' pstyle'
      DO i = 1, n
         IF (qrinside(X(i), GMM(1,1), GMM(1,2)) .and.
     *       qrinside(Y(i), GMM(2,1), GMM(2,2))      ) THEN
c            RYU   = dinside ( Y(i)+D(i), GMM(2,1), GMM(2,2) )
c            RYD   = dinside ( Y(i)-D(i), GMM(2,1), GMM(2,2) )
            CALL PS_Coord (X(i), Y(i),      z, xc, yc)
            CALL PS_Coord (X(i), Y(i)+D(i), z, xc, ycU)
            CALL PS_Coord (X(i), Y(i)-D(i), z, xc, ycD)
            IF (yCD.lt.0.0) THEN
               yCD = 0.0
            ENDIF
            yCU = dmin1 (ycU, 99.999d0) ! prevent overflow
            yCD = dmin1 (yCD, 99.999d0) ! prevent overflow
            IF     (i.eq.1) THEN
               Write (52, '(3(1x,f7.3),a)') xc, yc, yCU, ' ti'
               Write (52, '(3(1x,f7.3),a)') xc, yc, yCD, ' t'
            ELSEIF (i.eq.n) THEN
               Write (52, '(3(1x,f7.3),a/)') xc, yc, yCU, ' t'
               Write (52, '(3(1x,f7.3),a/)') xc, yc, yCD, ' tf'
            ELSE
               Write (52, '(3(1x,f7.3),a)') xc, yc, yCU, ' t'
               Write (52, '(3(1x,f7.3),a)') xc, yc, yCD, ' t'
               ENDIF
            ENDIF ! point
         ENDDO ! i

      END ! PS_PointLog

      SUBROUTINE PS_Curve (X, Y, n, z, iLine, rSyMag)
C     -----------------------------------------------
            ! JWu 11oct91.
         ! Plot Y(i) vs. X(i) as a polygon.
         ! The points have to be already forced
         ! into the graphic frame.

      IMPLICIT REAL*8    (a-h,o-p,r-z)
      IMPLICIT LOGICAL   (q)
      DIMENSION           X(*),Y(*)
      CHARACTER           h1*10, h2*10, text*80, cl2*2

      COMMON / output  / ioUnit

      IF (iLine.le.0) CALL Absturz ('PS_Curve', 'iLine .le. 0')

      Write (52, '(a,g12.4)')  '%: one curve / z = ', z
      Write (52, '(i2,a)')     iLine, ' cstyle'

      CALL SetCoord ('PS')
      CALL PS_Coord (X(1), Y(1), z, xc, yc)
      Write (52, '(3(1x,f7.3),a)') xc, yc, 0., ' ti'
      DO i = 2, n-1
         CALL PS_Coord (X(i), Y(i), z, xc, yc)
         Write (52, '(3(1x,f7.3),a)') xc, yc, 0., ' t'
         ENDDO ! i
      CALL PS_Coord (X(n), Y(n), z, xc, yc)
      Write (52, '(3(1x,f7.3),a/)') xc, yc, 0., ' tf'

      END ! PS_Curve

C  --------------------------------------------------------------------
C  g1.5.3.  Text
C  --------------------------------------------------------------------

      SUBROUTINE PS_Text (Text, Prefix, Suffix)
C     -----------------------------------------
            ! JWu 11jun91 dummy, 28feb92, macros activated 2jul92,
            ! Prefix, Suffix to include PS_Label 21jul93
         ! write a text to position no. iNr

      CHARACTER  Text*(*), Prefix*(*), Suffix*(*),
     *           aus*256, aux*256, h1*80, cl3*3, cr2*2

      COMMON / Viewport / iTOTx, iTOTy

      IF (Text.eq.'&start_text') THEN
         Write (52, '(a)')  'black 0 -4 13 newlist'
         RETURN
         ENDIF

      laus = len(aus)
      IF (lenU(Text).gt.laus+2) THEN
         Print *, 'PS_Text/ the following line is too long:'
         Print *, Text
         CALL Compose2 (aus, '('//Text(1:laus-2), ')')
      ELSE
         CALL Compose2 (aus, '('//Text, ')')
         ENDIF

C  Search macros :
      i = 0
      iKla = 0
 21   CONTINUE
      IF (i.lt.lenU(aus) .and. i.lt.laus) THEN
         i = i + 1
         IF     (aus(i:i).eq.'&') THEN
            IF     (aus(i+1:i+3).eq.'sy=') THEN
               CALL Fi1I (aus(i+4:i+5), isy)
               IF (aus(i+4:i+4).ne.'#') THEN
                  CALL Gong (13)
                  Print *, 'ERROR/ macro &sy in PS_Text'
                  Print *, aus
                  RETURN
                  ENDIF
               ! replace macro by PS command :
               CALL DelVonBis (aus, i, i+4)
               CALL Insert (aus, i,') 8 spce '//cr2(isy)
     *                      //' pstyle pins (')
               i = i-1
            ELSEIF (aus(i+1:i+3).eq.'li=') THEN
               CALL Fi1I (aus(i+4:i+5), ili)
               IF (aus(i+4:i+4).ne.'#') THEN
                  CALL Gong (13)
                  Print *, 'ERROR/ macro &li in PS_Text'
                  RETURN
                  ENDIF
               ! replace macro by PS command :
               CALL DelVonBis (aus, i, i+4)
               CALL Insert (aus, i, ') '//cr2(ili)//' cstyle cins (')
               i = i-1
            ELSE ! '&' means '&', nothing else 24oct96
               ENDIF
         ELSEIF (aus(i:i).eq.'(') THEN
            iKla = iKla + 1
         ELSEIF (aus(i:i).eq.')') THEN
            iKla = iKla - 1
            IF (iKla.lt.0) THEN
               CALL Gong (9)
               Print *, ' too many ''))'' in the following line:'
               Print *, aus
               CALL DelVonBis (aus, i, i)
               ENDIF
         ELSEIF (aus(i:i).eq.' ') THEN
            ! if there are many ' ', replace them by /spce (10jun94)
            DO ii = i+1, lenU(aus) ! find first character <> ' '
               IF (aus(ii:ii).ne.' ') GOTO 243
               ENDDO
 243           CONTINUE
            IF (ii-i.gt.4) THEN ! found many ' ' indeed
               CALL DelVonBis (aus, i, ii-1)
               CALL Compose2 (aux,  ') '//cl3(ii-i), ' spce (')
               laux = lenU(aux)
               CALL Insert (aus, i, aux(1:laux))
               i = i+laux
            ELSE
               i = ii-1
               ENDIF
            ENDIF ! special character

         GOTO 21 ! end loop i
         ENDIF

      IF (iKla.gt.0) THEN
         CALL Gong (9)
         Print *, ' too many ''(('' in the following line :'
         Print *, aus
         laus = len(aus)
         DO i = 1, iKla
            IF (aus(laus:laus).ne.' ') THEN
               Print *, 'line not written to PostScript'
               RETURN
               ENDIF
            CALL Append (aus, ')')
            ENDDO
         ENDIF

      Write (52,'(5a)')  Prefix, ' {', aus(1:lenU(aus)), '} ', Suffix

      END ! PS_Text

C  --------------------------------------------------------------------
C  g1.5.4.  Master : laserfile
C  --------------------------------------------------------------------

      SUBROUTINE OpenPS (FileExt, FileInt, Fehler)
C     --------------------------------------------
            ! JWu 11jun91
         ! open new postscript file as station 52, with default name l#.ps

      IMPLICIT LOGICAL   (q)

      CHARACTER*(*) FileExt, FileInt, Fehler
      CHARACTER*80  aus
      CHARACTER     cl3*3, cv2*2

      DATA          iPS /0/

      IF (Fehler.ne.'&ff') THEN
         CALL Gong (17)
         Print *, ' Error on entry in OpenPS'
         RETURN
         ENDIF

      IF (FileExt.eq.' ') THEN
         ! create default file name :
 21      CONTINUE
         Fehler = '&ff'
         iPS = iPS+ 1
         FileInt = 'l'//cl3(iPS)

         IF (mod(iPS,50).eq.0) THEN
            ! ask whether all right :
            CALL Compose2 (aus, 'Already opening file '//FileInt,
     *           ' - continue ?')
            IF (.not.qAsk (aus)) THEN
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF

         !  Try to open the file :
         CALL OpenDatFile (52, FileInt, 'ps', 'n', 'seq',
     *                     'lis', 0, Fehler)
         IF (Fehler.ne.'&ff') GOTO 21 ! try next
         Print *, ' postscript file '//FileInt

      ELSE
         FileInt = FileExt
         CALL OpenDatFile (52, FileInt, 'ps', 'n!', 'seq',
     *                     'lis', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF

      END ! OpenPS

      SUBROUTINE CopyPS (From, Fehler)
C     --------------------------------
            ! JWu 11jun91 in OpenPS, 13sep98 separately
         ! include initialization file (g3.ps or other)

      CHARACTER*(*) From, Fehler
      CHARACTER*80  FileAux, aus, line, form
      CHARACTER     cl3*3, cv2*2

C  Copy the auxiliary file :
      CALL ExeML ('\p '//From, FileAux)
 33   CONTINUE
      CALL OpenFile (53, FileAux, '&noext', 'l', Fehler)
      IF (Fehler.ne.'&ff') THEN
         CALL Gong (3)
         CALL Say2 (' Cannot open PostScript definition file "'//
     *              FileAux(1:lenU(FileAux)), '"')
         Fehler = '&ff'
         CALL FrageC ('Try new file name ?', FileAux)
         IF (FileAux.ne.' ') GOTO 33
         Fehler = ' '
         RETURN
         ENDIF

C  loop - copy from 53 to 52 :
      DO i = 1, 1200
         read (53, '(a80)', end=98, err=99) line
         ll = lenU(line)
         IF (ll.le.0) THEN
            write (52, '( )')
         ELSE
            form = '(a'//cv2(ll)//')'
            write (52, form) line
            ENDIF
         ENDDO ! loop copy

C  reached MaxLines :
      Fehler = 'CopyPs/ too many lines/ something is wrong'
      RETURN

C  end-of file :
 98   CONTINUE
      IF (i.le.3) THEN
         Fehler = 'gra-su-file too short: eof in line '//cl3(i)
         RETURN
         ENDIF

C  all right :
      Close (53)
      RETURN

C  read error :
 99   CONTINUE
      Fehler = 'SEVERE/ CopyPS/ reading the auxiliary file'
      RETURN

      END ! CopyPS

      SUBROUTINE ClosePS
C     ------------------
            ! JWu 12jun91
         ! close last postscript file

C  Write end commands to postscript :
      Write (52, '(/a)') 'EndFrame'
      ! no longer needed ? Write (52, '(a1)')  char(4)

C  Close the file :
      Close (Unit=52)

      END ! ClosePS
