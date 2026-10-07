C  ====================================================================
C
C     Library  WuGra :  Graphics
C     Modul    g2.f  :  Interface
C
C  ====================================================================

C     J. Wuttke

C     Contents :
C
C        outer shell / setup, info :
C                 WdwPreset, GraSetup, GraDims, GraChoice,
C                 GraSetWdw, GraWdwList, GraSetAxs
C        outer shell / set CS :
C                 GraSetSca, GraInquireCS
C        outer shell / plot :
C                 GraQuit, GraPlotCS, GraPaint, GraText, GraSoftCopy
C        Tektronix screen driver :
C                 SetTek, SetScroll
C        auxiliary / math for line drawing :
C                 SectBoxLin, qSectLinLin
C        auxiliary / compose labels :
C                 GraLabel

C     Major modifications :
C        Sep98  TEK file and LIN output suppressed, restructured
C        Jan94  3d begun
C        Jul93  Text register revised
C        Apr93  Window switch
C        Nov92  Some optimization for CCNY
C        Oct91  Module WuG2 restructured
C        Jun91  Option PostScript
C        Nov90  Options TEK/ILL
C        Okt90  Library WuGra
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  ====================================================================
C  outer shell / setup ...
C  ====================================================================

      BLOCK DATA WdwPreset
C     --------------------
            ! 23feb, 12apr93
         ! Default choice for window-setups

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h, o-p, r-z)

      INCLUDE 'g_dim.f'

      COMMON / GraWdw / iGW, nGW
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSpe / qErrBar(MGW), qStandSymb(MGW), rSyMag(MGW),
     *                  qForce(3,MGW)
      COMMON / GraSca / LinLog(3,MGW)
      COMMON / GraAx  / AxCro(3,MGW), AxAng(2,MGW), AxLen(3,MGW),
     *                  nWdim(MGW), qBox2(MGW)

      DATA    iGW, nGW / 1, 0 /
      DATA    qErrBar, qStandSymb, rSyMag / MGW*.false.,
     *        MGW*.true., MGW*1.d0 /
      DATA    qForce / MGW*.false., MGW*.true., MGW*.false. /
      DATA    LinLog / MGW*0, MGW*0, MGW*0 /
      DATA    Gmm    / MGW6*0.d0 /
      DATA    AxCro, AxAng, AxLen, nWdim, qBox2 /
     *                 MGW*0.d0, MGW*0.d0, MGW*0.d0, MGW*0.d0,
     *                 MGW*0.d0, MGW*1.d0, MGW*1.d0, MGW*1.d0,
     *                 MGW*2, MGW*.true. /

      END ! WdwPreset

      SUBROUTINE GraSetup ()
C     ----------------------
            ! as Block Data 23feb93, read from ML-par 1jun95
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

      COMMON / GraTerm/ qWindow, qGTOverlay, qGToldAdr
      COMMON / GraPlus/ iGraP, iEsc
      COMMON / Scroll / iScroIst, iScroGra
      COMMON / Format / kScrF, kLasF, kPS_F

      qWindow    = qintr (iExeMLP('su-gra-window'))
      qGTOverlay = qintr (iExeMLP('su-gra-ovrlay'))
      qGToldAdr  = qintr (iExeMLP('su-gra-oldAdr'))

      iGraP    = iExeMLP ('su-gra-iGraP')
      iEsc     = iExeMLP ('su-gra-Escape')
      iScroIst = iExeMLP ('su-scroll-ist')
      iScroGra = iExeMLP ('su-scroll-gra')

      kScrF    = iExeMLP ('su-tekwdw-Scr')
      kLasF    = iExeMLP ('su-tekwdw-Las')

      END ! GraSetup

      SUBROUTINE GraDims ()
C     ---------------------
            ! JWu 23nov92
         ! Print current array dimensions

      INCLUDE 'g_dim.f'

      Print '(a)', ' Current array dimensions (graphic register) :'
      Print '(a,i8)', ' # spectra           ', MKreg
      Print '(a,i8)', ' # channels/spectrum ', MCgra
      Print '(a,i8)', ' # total channels    ', MGreg
      Print '(a,i8)', ' # text lines        ', MTreg
      Print '(a,i8)', ' # labels            ', MSreg

      END ! GraDims

      SUBROUTINE GraChoice (What, Object, Fehler)
C     -------------------------------------------
            ! renove 9jul93

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      CHARACTER         What*(*), Object*(*), Fehler*(*), aus*80, cl3*3

      INCLUDE 'g_dim.f'

      EXTERNAL          WdwPreset

      COMMON / Scroll / iScroIst, iScroGra
      COMMON / GraTerm/ qWindow, qGTOverlay, qGToldAdr
      COMMON / Format / kScrF, kLasF, kPS_F

      COMMON / GraWdw / iGW, nGW
      COMMON / GraSpe / qErrBar(MGW), qStandSymb(MGW), rSyMag(MGW),
     *                  qForce(3,MGW)
      COMMON / GraSca / LinLog(3,MGW)
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraAx  / AxCro(3,MGW), AxAng(2,MGW), AxLen(3,MGW),
     *                  nWdim(MGW), qBox2(MGW)

      IF     (What.eq.' ') THEN  ! help
         Print *, 'graphics commands :'
           ! --- interpretation of these commands must be user-written :
         Print *, '   [files] p [spectra]   : plot'
         Print *, '   [files] a [spectra]   : add to plot'
         Print *,
     * '   gs [filename]         : copy plot to PostScript file'
         Print *,
     * '   gp [filename]         : copy plot to PS file inlc head'
         Print *, '   ga [filename]         : - without definitions'
           ! --- the following commands are implemented in this routine :
         Print *, '   g-                    : close graphic display'
         Print *, '   gg                    : display size'
         Print *, '   gw [window-no]        : set window'
         Print *, '   g:                    : table of the following :'
         Print *,
     * '   glx/y/z               : logarithmic x/y/z scale ?'
         Print *,
     * '   gox/y/z [fu-no [arg]] : functional operation on x/y/z'
         Print *, '   gfx/y/z               : force x/y/z into frame ?'
         Print *, '   gd                    : default symbols ?'
         Print *, '   ge                    : error bars ?'
         Print *, '   gr [size]             : symbol radius'
         Print *, '   gb                    : box ?'

      ELSEIF (What.eq.'-') THEN
         CALL SetScroll (0)

      ELSEIF (What.eq.'g') THEN
         Print '(a)',  ' modify graphic format :'
         Print '(a,i3)', '    ( 1) # scroll lines       =', iScroGra
         Print '(a,i3)', '    ( 2) laser format         =', kLasF
         Print '(a,i3)', '    ( 3) PostScript format    =', kPS_F
 12      CONTINUE
         iMod = iAskDuMu (' Modify option [exit] ?', 0, 0, 3)
         IF     (iMod.eq. 0) THEN ! exit
            RETURN
         ELSEIF (iMod.eq. 1) THEN
            IF (qWindow)THEN
               CALL Gong (2)
               Print *, ' Not applicable for multi-window terminal'
            ELSE
               iScroGra = iAskDM (' Lines for graphics',
     *                            iScroGra, 8, 25)
               IF (iScroGra.ge.20) THEN
                  kScrF = 5 ! Grossformat
               ELSE
                  kScrF = 6 ! Kleinformat
                  ENDIF
               IF (iScroGra.lt.iScroIst) CALL SetScroll (iScroGra)
               ENDIF
         ELSEIF (iMod.eq. 2) THEN
            aus = ' Laser format : Einklebe(1), Hoch(2), Quer(3)'
            kLasF = iAskMu (aus, 1, 3)
         ELSEIF (iMod.eq. 3) THEN
            aus = ' Postscript format : Hoch(2), Quer(3)'
            kPS_F = iAskMu (aus, 2, 3)
            ENDIF
         GOTO 12

      ELSEIF (What.eq.'w') THEN

         CALL Fi1I (Object, iO)
         IF     (Object.eq.' ') THEN
            CALL GraWdwList ()
            iO = iAskDMu (' Switch to graphic window no. ',
     *                    iGW, 1, nGW+1)
         ELSEIF (Object.ne.'#') THEN
            Print *, ' usage : gw (window-number)'
            RETURN
         ELSEIF (qioutside(iO,1,MGW)) THEN
            Fehler = ' Window number outside allowed range 1..'//
     *               cl3(MGW)
            RETURN
            ENDIF
         iGW = iO

      ELSEIF (What.eq.':') THEN
         CALL Say2 (' setup of window '//cl3(iGW), ' :')
         Print '(a,i3)',   '    (e)  error bars         =',
     *         intq(qErrBar(iGW))
         Print '(a,i3)',   '    (d)  default symbols    =',
     *         intq(qStandSymb(iGW))
         Print '(a,f6.3)', '    (r)  symbol radius      = ',
     *         rSyMag(iGW)
         Print '(a,i3)',   '    (b)  box                =',
     *         intq(qBox2(iGW))

         Print '(a,i3)',   '    (lx) log x scale        =',
     *         LinLog(1,iGW)
         Print '(a,i3)',   '    (ly) log y scale        =',
     *         LinLog(2,iGW)
         Print '(a,i3)',   '    (lz) log z scale        =',
     *         LinLog(3,iGW)

         Print '(a,i3)',   '    (fx) force x into frame =',
     *         intq(qForce(1,iGW))
         Print '(a,i3)',   '    (fy) force y into frame =',
     *         intq(qForce(2,iGW))
         Print '(a,i3)',   '    (fz) force z into frame =',
     *         intq(qForce(3,iGW))

         Print '(a,i3,2x,g9.3)','    (ox) operation on x     =',
     *                                  iGFu(1,iGW), rGFu(1,iGW)
         Print '(a,i3,2x,g9.3)','    (oy) operation on y     =',
     *                                  iGFu(2,iGW), rGFu(2,iGW)
         Print '(a,i3,2x,g9.3)','    (oz) operation on z     =',
     *                                  iGFu(3,iGW), rGFu(3,iGW)

      ELSEIF (What.eq.'e') THEN
         qErrBar(iGW) = .not. qErrBar(iGW)
      ELSEIF (What.eq.'d') THEN
         qStandSymb(iGW) = .not.qStandSymb(iGW)
      ELSEIF (What.eq.'r') THEN
         aus = ' Symbol radius ?'
         rSyMag(iGW) = rAskMu (aus, 0.d0, 1.d1)
      ELSEIF (What.eq.'b') THEN
         qBox2(iGW) = .not. qBox2(iGW)

      ELSEIF (What.eq.'lx') THEN
         LinLog(1,iGW) = 1 - LinLog(1,iGW)
      ELSEIF (What.eq.'ly') THEN
         LinLog(2,iGW) = 1 - LinLog(2,iGW)
      ELSEIF (What.eq.'lz') THEN
         LinLog(3,iGW) = 1 - LinLog(3,iGW)

      ELSEIF (What.eq.'fx') THEN
         qForce(1,iGW) = .not. qForce(1,iGW)
      ELSEIF (What.eq.'fy') THEN
         qForce(2,iGW) = .not. qForce(2,iGW)
      ELSEIF (What.eq.'fz') THEN
         qForce(3,iGW) = .not. qForce(3,iGW)

      ELSEIF (What.eq.'ox') THEN
         CALL FuAsk (' Function for x-scale',
     *               Object, iGFu(1,iGW), rGFu(1,iGW))
      ELSEIF (What.eq.'oy') THEN
         CALL FuAsk (' Function for y-scale',
     *               Object, iGFu(2,iGW), rGFu(2,iGW))
      ELSEIF (What.eq.'oz') THEN
         CALL FuAsk (' Function for z-scale',
     *               Object, iGFu(3,iGW), rGFu(3,iGW))

      ELSE
         Fehler = ' this option not accessible'
         ENDIF

      END ! GraChoice

      SUBROUTINE GraSetWdw (ndim, Coor, Unit, Fehler)
C     -----------------------------------------------
            ! JWu 12mar96
         ! Change window if dim has changed or if Co/Un don't fit

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      INCLUDE 'g_dim.f'
      CHARACTER*(*)      Coor(3), Unit(3), Fehler
      CHARACTER*40       GCoor, GUnit

      COMMON / GraWdw / iGW, nGW
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSca / LinLog(3,MGW)
      COMMON / GraLab / GCoor(3,MGW), GUnit(3,MGW)
      COMMON / GraAx  / AxCro(3,MGW), AxAng(2,MGW), AxLen(3,MGW),
     *                  nWdim(MGW), qBox2(MGW)

      IF (Fehler.ne.'&ff') RETURN
      iGWold = iGW
 1    CONTINUE
      IF (iGW.gt.nGW) THEN ! open new window
         nWdim(iGW) = ndim
         DO idim = 1, ndim
            GUnit(idim,iGW) = Unit(idim)
            ENDDO
         nGW = min0(iGW, MGW-1)
         IF (iGWold.ge.1) THEN ! inherit setup from previous window
            DO id = 1, 3
               IF (Unit(id).eq.GUnit(id,iGWold)) THEN
                  LinLog(id,iGW) = LinLog(id,iGWold)
                  Gmm(id,1,iGW)  = Gmm(id,1,iGWold)
                  Gmm(id,2,iGW)  = Gmm(id,2,iGWold)
                  ENDIF
               ENDDO
            ENDIF
      ELSE
         IF (ndim.ne.nWdim(iGW) .or.
     *       Unit(1).ne.GUnit(1,iGW) .or. Unit(2).ne.GUnit(2,iGW) .or.
     *       (ndim.eq.3 .and. Unit(3).ne.GUnit(3,iGW))) THEN
            ! circular search for other window :
            iGW = iGW-1
            IF (iGW.le.0) iGW = nGW
            IF (iGW.eq.iGWold) iGW = nGW + 1
            GOTO 1
            ENDIF
         ENDIF

      DO idim = 1, ndim
         GCoor(idim,iGW) = Coor(idim)
         ENDDO

      END ! GraSetWdw

      SUBROUTINE GraWdwList ()
C     ------------------------
         ! List contents of activated graph windows (9aug93)

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      INCLUDE 'g_dim.f'
      CHARACTER         Activ*1, LX*1, LY*1, GCoor*40, GUnit*40

      COMMON / GraWdw / iGW, nGW
      COMMON / GraLab / GCoor(3,MGW), GUnit(3,MGW)
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSca / LinLog(3,MGW)

      DO i = 1, MGW
         IF (Gmm(1,1,i).ne.Gmm(1,2,i)) THEN ! else the window hasn't been used
            IF (i.eq.iGW) THEN
               Activ = '*'
            ELSE
               Activ = ' '
               ENDIF
            IF (LinLog(1,i).eq.1) THEN
               LX    = 'L'
            ELSE
               LX    = ' '
               ENDIF
            IF (LinLog(2,i).eq.1) THEN
               LY    = 'L'
            ELSE
               LY    = ' '
               ENDIF

            Print '(i2,a1,1x,2(a1,1x,a7,a6,1x,g10.2,g10.2,3x))',
     *       i, Activ,
     *       LX, GCoor(1,i), '('//GUnit(1,i)(1:lenU(GUnit(1,i)))//')',
     *       Gmm(1,1,i), Gmm(1,2,i),
     *       LY, GCoor(2,i), '('//GUnit(2,i)(1:lenU(GUnit(2,i)))//')',
     *       Gmm(2,1,i), Gmm(2,2,i)

            ENDIF
         ENDDO

      END ! GraWdwList

      SUBROUTINE GraSetAxs (Fehler)
C     -----------------------------
            ! JWu 12mar96 separated from GraSetCS
         ! Determine situation of coordinate axis in the drawing

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      INCLUDE 'g_dim.f'
      CHARACTER          Fehler*(*)
      CHARACTER*40       CoUn(6)

      COMMON / GraWdw / iGW, nGW
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraAx  / AxCro(3,MGW), AxAng(2,MGW), AxLen(3,MGW),
     *                  nWdim(MGW), qBox2(MGW)

 41   CONTINUE
      IF (nWdim(iGW).eq.3) THEN
         CALL rAskTrip (' x/y/z Axes cross in point',
     *      AxCro(1,iGW), AxCro(2,iGW), AxCro(3,iGW),
     *      AxCro(1,iGW), AxCro(2,iGW), AxCro(3,iGW))
         CALL rAskTrip (' Relative length of axes',
     *      AxLen(1,iGW), AxLen(2,iGW), AxLen(3,iGW),
     *      AxLen(1,iGW), AxLen(2,iGW), AxLen(3,iGW))
         DO j = 1, 3
            IF (qroutside(AxCro(j,iGW), Gmm(j,1,iGW), Gmm(j,2,iGW)))
     *         THEN
               CALL Gong(6)
               Print *, ' Axes cross must be within coordinate limits'
               GOTO 41
               ENDIF
            ENDDO
         CALL rAskPair (' And x/z axes have angles',
     *      AxAng(1,iGW), AxAng(2,iGW), AxAng(1,iGW), AxAng(2,iGW))
      ELSEIF (.not.qBox2(iGW)) THEN
         CALL rAskPair (' x/y Axes cross in point',
     *      AxCro(1,iGW), AxCro(2,iGW), AxCro(1,iGW), AxCro(2,iGW))
         DO j = 1, 2
            IF (qroutside(AxCro(j,iGW), Gmm(j,1,iGW), Gmm(j,2,iGW)))
     *         THEN
               CALL Gong(6)
               Print *, ' Axes cross must be within coordinate limits'
               GOTO 41
               ENDIF
            AxLen(j,iGW) = 1
            ENDDO
         ENDIF

      END ! GraSetAxs

C  ====================================================================
C  outer shell / set scale
C  ====================================================================

      SUBROUTINE GraSetSca (iD, cD, si, sf, action, qLin, igf, rgf)
C     -------------------------------------------------------------
            ! JWu oct90. revisions 19sep/11oct91, 18jan93.
            ! Simplified version, recursive interaction with i0: 12mar96.
         ! Determine graphic range

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      INCLUDE 'g_dim.f'
      CHARACTER          aus*80, ein*80, Fehler*80, what*20, cD*1,
     *                   action*1
      CHARACTER*40       CoUn(6)

      COMMON / GraWdw / iGW, nGW
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSca / LinLog(3,MGW)

      Fehler = '&ff'

      IF (si.ge.sf) THEN ! no default limits given => take last ones
         si = Gmm(iD,1,iGW)
         sf = Gmm(iD,2,iGW)
         ENDIF

 10   CONTINUE

      ! Menu (except in some special cases) :
      IF (si.ge.sf) THEN ! should happen only on 1st call
         ein = 'n'
      ELSEIF (ein.eq.'a') THEN
         ein = ' ' ! accept default, don't ask for confirmation
      ELSE
         aus = ' Rescale '//cD//' (h=help)'
         ein = ' '
         CALL rAskRgeTxt (aus, ein, si, sf, si, sf)
         ENDIF

 11   CONTINUE
      IF     (ein.eq.'h' .or. ein.eq.'?') THEN
         Print *, ' rescale '//cD//'-range of plot :'
         Print *, '    RETURN = accept default range'
         Print *, '    r1 r2  = new range r1..r2'
         Print *, '    r1,    = overwrite r1, accept other default'
         Print *, '    ,r2    = overwrite r2, accept other default'
         Print *, '    a      = automatic calculation'
         Print *, '    n      = calculate new default'
         Print *, '   <g-opt> = setup commands d,e,r,lx/y,..,'
         Print *, '    :      = setup list'
         IF (iD.gt.1)
     *   Print *, '    x      = correct x-range'
         Print *, '    -      = exit, no plot'
         GOTO 10

      ELSEIF (ein.eq.' ') THEN ! accept default
         ! test whether new limits are consistent with log plot :
         CALL FuVal (iGFu(iD,iGW),
     *               siF, dummy, si, 0.d0, rGFu(iD,iGW), 0.d0)
         CALL FuVal (iGFu(iD,iGW),
     *               sfF, dummy, sf, 0.d0, rGFu(iD,iGW), 0.d0)
         IF (LinLog(iD,iGW).eq.1 .and.
     *       (siF.le.0 .or. sfF.le.0)) THEN
            Print *, ' limits inconsistent with log scale'
            CALL Gong (3)
            GOTO 10
         ELSE
            Gmm(iD,1,iGW) = si
            Gmm(iD,2,iGW) = sf
            action = ' '
            ENDIF

      ELSEIF (jPos1('n-x',ein(1:1)).le.3) THEN
         action = ein(1:1)

      ELSEIF (ein.eq.'a') THEN
         action = 'n'

      ELSE
         CALL TakeVorDel (ein, what, ' ')
         CALL GraChoice (what, ein, Fehler)
         IF (Fehler.ne.'&ff') THEN
            CALL FehlerGong (Fehler, 3)
            GOTO 10
         ELSEIF (what(1:1).eq.'l') THEN
            IF (what(2:2).eq.cD) THEN
               action = 'n'
            ELSEIF (ichar(what(2:2)).lt.ichar(cD)) THEN
               action = what(2:2)
               ENDIF
         ELSEIF (what(1:1).eq.'w') THEN
            action = 'x'
         ELSE
            GOTO 10
            ENDIF

         ENDIF

      ! the following parameters are needed for calculating new limits :
      qLin = .not. qintr(LinLog(iD,iGW))
      igf  = iGFu(iD,iGW)
      rgf  = rGFu(iD,iGW)

      END ! GraSetSca

      SUBROUTINE GraInquireCS (cxy, rmi, rma, iFu, rFuPar, lilo)
C     ----------------------------------------------------------
         ! Access to COMMON / GraCS / for application program

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      INCLUDE 'g_dim.f'
      CHARACTER         cxy*1

      COMMON / GraWdw / iGW, nGW
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSca / LinLog(3,MGW)

      j = 0
      IF (cxy.eq.'x') j=1
      IF (cxy.eq.'y') j=2
      IF (cxy.eq.'z') j=3
      IF (j.eq.0) CALL Absturz ('GraInquire', 'cxy o.o.r.')

      rmi    = Gmm   (j,1,iGW)
      rma    = Gmm   (j,2,iGW)
      iFu    = iGFu  (j,iGW)
      rFuPar = rGFu  (j,iGW)
      lilo   = LinLog(j,iGW)

      END ! GraInquireCS

C  ====================================================================
C  outer shell / plot
C  ====================================================================

      SUBROUTINE GraQuit ()
C     ---------------------
         ! Quit the graph mode to continue with text.

      IMPLICIT LOGICAL (q)
      CHARACTER  ein *80

      COMMON / GraTerm/ qWindow, qGTOverlay, qGToldAdr

      IF ((.not.qGTOverlay) .and. (.not.qWindow)) THEN
         CALL FrageC (' Say something to continue', ein)
         ENDIF

      CALL GMode (0)

      END ! GraQuit

      SUBROUTINE GraPlotCS (Fehler)
C     -----------------------------
         ! Open or clear the graphics, plot scales, labels and title.
         ! The coordinate window should be set by a previous call
         ! to GraSetCS.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)
      INCLUDE 'g_dim.f'

      CHARACTER         RegTX*80, RegSX*80, GCoor*40, GUnit*40,
     *                  Labl(3)*40, Fehler*(*)

      COMMON / GraWdw / iGW, nGW
      COMMON / GraLab / GCoor(3,MGW), GUnit(3,MGW)
      COMMON / GraReg / NReg(MKreg), PReg(MGreg,3), ZReg(MKreg), iReg,
     *                  iRLS(MKreg), RegTX(MTreg), RegSX(MSreg,6),
     *                  iRTS(7)
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSca / LinLog(3,MGW)
      COMMON / GraAx  / AxCro(3,MGW), AxAng(2,MGW), AxLen(3,MGW),
     *                  nWdim(MGW), qBox2(MGW)

      CALL GraSetup ()

      DO idim = 1, 3
         CALL GraLabel (idim, Labl(idim), GCoor(idim,iGW),
     *                  GUnit(idim,iGW))
         ENDDO

C  Clear GraReg :
      iReg = 0
      DO j=1,7
         iRTS(j) = 0
         ENDDO

C  Enregister :
      RegSX(1,1) = ' ' ! title
      RegSX(1,2) = Labl(1)
      RegSX(1,3) = Labl(2)
      iRTS(1)    = 1
      iRTS(2)    = 1
      iRTS(3)    = 1
      CALL GraText ('&start_text', 3)

C  Open graphics :
      CALL SetTek ()
      CALL GMode (1) ! added for xterm
      CALL ClearGraphic

C  Plot co-ordinate system and text :
      CALL GMode (1)
      qBox = (qBox2(iGW) .and. nWdim(iGW).eq.2)
      CALL TekPlotCS (LinLog(1,iGW), Gmm(1,1,iGW),
     *                AxCro(1,iGW), AxAng(1,iGW),
     *                AxLen(1,iGW), nWdim(iGW), qBox, Fehler)
      CALL TekText (2, '&right '//Labl(1))
      CALL TekText (3, '&right '//Labl(2))

      END ! GraPlotCS

      SUBROUTINE GraPaint (X, Y, D, n, z, iLSin, iLS)
C     -----------------------------------------------
         ! Choose linestyle/symbols, transform data, call ExePaint.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      INCLUDE 'g_dim.f'
      DIMENSION  X(*), Y(*), D(*), Xp(MCgra), Yp(MCgra), Dp(MCgra),
     *           XS(4), YS(4) ! should be (2)
      CHARACTER  RegTX*80, RegSX*80

      COMMON / GraWdw / iGW, nGW
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSpe / qErrBar(MGW), qStandSymb(MGW), rSyMag(MGW),
     *                  qForce(3,MGW)
      COMMON / GraReg / NReg(MKreg), PReg(MGreg,3), ZReg(MKreg), iReg,
     *                  iRLS(MKreg), RegTX(MTreg), RegSX(MSreg,6),
     *                  iRTS(7)

C  Checks :
      IF (n.gt.MCgra) THEN
         CALL GMode(0)
         Print *, ' Graphic overflow/ too many points'
         RETURN
         ENDIF

C  Choice of linestyle/symbol (help suppressed 28jan92) :
      IF (qStandSymb(iGW)) THEN
         iLS = iLSin
      ELSE
         ! construct default :
         IF     (iReg.le.0) THEN
            iLSdef = iLSin
         ELSEIF (iLS.gt.0) THEN
            iLSdef = iLS + 1
         ELSE
            iLSdef = iLS - 1
            ENDIF
         ! ask :
         CALL GMode (0)
         iLS = iAskD (' Linestyle(>0) or plotsymbol(<0)', iLSdef)
         ENDIF
      IF (iLS.eq.0) RETURN ! break
      qPoints = (iLS.lt.0)

C  Loop over data points (all operations here since 6nov92) :
      np = 0
      qDrawing = .false. ! old point (xa,ya) is still undefined
      DO i = 1, n
C  - Functional transforms :
         IF (iGFu(1,iGW).eq.0) THEN
            xi = X(i)
         ELSE
            IF (iGFu(1,iGW).eq.1 .or. iGFu(1,iGW).eq.2) THEN
               IF (X(i).le.0.) GOTO 19
               ENDIF
            CALL FuVal (iGFu(1,iGW), xi, dummy, X(i), 0.d0,
     *                  rGFu(1,iGW), 0.d0)
            ENDIF
         IF (iGFu(2,iGW).eq.0) THEN
            yi = Y(i)
            di = D(i)
         ELSE
            IF (iGFu(2,iGW).eq.1 .or. iGFu(2,iGW).eq.2) THEN
               IF (Y(i).le.0.) GOTO 19
               ENDIF
            CALL FuVal (iGFu(2,iGW), yi, di, Y(i), D(i),
     *                  rGFu(2,iGW), 0.d0)
            ENDIF

         IF (qPoints) THEN
C  - Force into graphic range ?
            IF (qForce(1,iGW)) THEN
               xi = dinside (xi, Gmm(1,1,iGW), Gmm(1,2,iGW))
            ELSEIF (qroutside(xi, Gmm(1,1,iGW), Gmm(1,2,iGW))) THEN
               GOTO 19
               ENDIF
            IF (qForce(2,iGW)) THEN
               yi = dinside (yi, Gmm(2,1,iGW), Gmm(2,2,iGW))
            ELSEIF (qroutside(yi, Gmm(2,1,iGW), Gmm(2,2,iGW))) THEN
               GOTO 19
               ENDIF
         ELSE ! draw line
C  - Handle intersections with frame :
            IF (np.gt.MCgra-2) THEN
               CALL GMode (0)
               Print *, ' Graphic overflow/ too many points'
               GOTO 90
               ENDIF
            IF (qDrawing) THEN
               CALL SectBoxLin (Gmm(1,1,iGW), Gmm(1,2,iGW),
     *                          Gmm(2,1,iGW), Gmm(2,2,iGW),
     *                          xa, ya, xi, yi,
     *                          qIn1, qIn2, nSect, XS, YS)
               DO jSect = 1, nSect ! ainsi simplifie le 22mai94
                  Xp(np+jSect) = XS(jSect)
                  Yp(np+jSect) = YS(jSect)
                  Dp(np+jSect) = -1000-100*intq(qIn1)-
     *                           10*intq(qIn2)-jSect
                  ENDDO
               np = np + nSect
               IF (.not.qIn2) GOTO 18
            ELSE
               IF (.not.(qrinside(xi,Gmm(1,1,iGW),Gmm(1,2,iGW)) .and.
     *             qrinside(yi,Gmm(2,1,iGW),Gmm(2,2,iGW)))) GOTO 18
               ENDIF
            ENDIF ! points/line
C  - Retain entry :
         np = np + 1
         Xp(np) = xi
         Yp(np) = yi
         Dp(np) = di
C  - End of loop :
 18      CONTINUE
         qDrawing = .true.
         xa = xi ! old points
         ya = yi
 19      CONTINUE
         ENDDO

C  Enregister (one-dimensional storage 29oct92) :
      IF (iReg.eq.0) THEN
         iOff = 0
      ELSE
         iOff = NReg(iReg)
         ENDIF
      IF     (iOff+np.gt.MGreg) THEN
         Print *, 'GraReg overflow/ too many points'
      ELSEIF (iReg.ge.MKreg) THEN
         Print *, 'GraReg overflow/ too many spectra'
      ELSE
         iReg = iReg + 1
         NReg(iReg) = iOff + np ! last filled line
         iRLS(iReg) = iLS
         DO i = 1, np
            PReg(iOff+i,1) = Xp(i)
            PReg(iOff+i,2) = Yp(i)
            PReg(iOff+i,3) = Dp(i)
            ENDDO
         ZReg(iReg) = z
         ENDIF

C  Plot points :
 90   CONTINUE
      CALL GMode (1)
      IF (iLS.lt.0) THEN
         CALL TekPoint (Xp, Yp, Dp, np, z, -iLS, rSyMag(iGW),
     *                  qErrBar(iGW))
      ELSE
         CALL TekCurve (Xp, Yp, np, z, iLS, rSyMag(iGW))
         ENDIF

      END ! GraPaint

      SUBROUTINE GraText (text, idev)
C     -------------------------------
             ! rewritten JWu 22jul93
         ! decode macros, enregister, call TekText.

         ! idev = 1   : TEK on screen
         !        2   : register (and later on file)
         !        3   : both (1 and 2)

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE 'g_dim.f'
      CHARACTER         RegTX*80, RegSX*80, text*(*)
      SAVE              iTEKpos

      COMMON / GraReg / NReg(MKreg), PReg(MGreg,3), ZReg(MKreg), iReg,
     *                  iRLS(MKreg), RegTX(MTreg), RegSX(MSreg,6),
     *                  iRTS(7)

C  Enregister :
      IF (idev.eq.2 .or. idev.eq.3) THEN
         iRTS(7) = iRTS(7) + 1
         IF (iRTS(7).gt.MTreg) GOTO 19
         RegTX(iRTS(7)) = text
         DO ii = 81, lenU(text), 75
            iRTS(7) = iRTS(7) + 1
            IF (iRTS(7).gt.MTreg) GOTO 19
            RegTX(iRTS(7)) = '... '//text(ii:min0(lenU(text),ii+74))
            ENDDO
         ENDIF
 19   CONTINUE

C  Plot :
      IF (idev.eq.1 .or. idev.eq.3) THEN
         IF     (text.eq.'&start_text') THEN
            iTEKpos = 11
         ELSE
            CALL TekText (iTEKpos, text)
            iTEKpos = iTEKpos + 1
            ENDIF
         ENDIF

      END ! GraText

      SUBROUTINE GraSoftCopy (FileExt, FileIniMac, Fehler)
C     ----------------------------------------------------
         ! A copy of the current graphics (as saved in / GraReg / )
         ! is written to an output file optionally given by Object.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)
      INCLUDE 'g_dim.f'
      CHARACTER*(*) FileExt, FileIniMac, Fehler
      CHARACTER     FileInt*80, RegTX*80, RegSX*80, aus*240, Datum*10,
     *              Zeit*10,Datuma*24
      REAL*8        X(MCgra), Y(MCgra), D(MCgra)

      COMMON / GraReg / NReg(MKreg), PReg(MGreg,3), ZReg(MKreg), iReg,
     *                  iRLS(MKreg), RegTX(MTreg), RegSX(MSreg,6),
     *                  iRTS(7)
      COMMON / GraWdw / iGW, nGW
      COMMON / GraSpe / qErrBar(MGW), qStandSymb(MGW),
     *                  rSyMag(MGW), qForce(3,MGW)
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)
      COMMON / GraSca / LinLog(3,MGW)
      COMMON / GraAx  / AxCro(3,MGW), AxAng(2,MGW), AxLen(3,MGW),
     *                  nWdim(MGW), qBox2(MGW)

C  Checks :
      IF (iReg.le.0) THEN
         Fehler = 'graphic register is empty'
         RETURN
         ENDIF

C  Open output file :
      CALL OpenPS (FileExt, FileInt, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL CopyPS (FileIniMac, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL PS_Text (FileInt, ' ',
     *     ' /filename exch def 10 -3 18 showfilename')

C  Copy plot from register :

c     nDi = nWdim(iGW)                   ! unused
c     qBox = (qBox2(iGW) .and. nDi.eq.2)

      CALL PS_PlotCS (LinLog(1,iGW), Gmm(1,1,iGW),
     *                RegSX(1,2), RegSX(1,3), Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Data points and curves :

      DO ir = 1, iReg
         IF (ir.eq.1) THEN
            iOff = 0
         ELSE
            iOff = NReg(ir-1)
            ENDIF
         n = NReg(ir) - iOff
         DO i = 1, n
            X(i) = PReg(iOff+i,1)
            Y(i) = PReg(iOff+i,2)
            D(i) = PReg(iOff+i,3)
            ENDDO

         IF (iRLS(ir).lt.0.and.LinLog(2,iGW).eq.0) THEN
            CALL PS_Point (X, Y, D, n, ZReg(ir), -iRLS(ir),
     *                     rSyMag(iGW), qErrBar(iGW))
         ELSEIF (iRLS(ir).lt.0.and.LinLog(2,iGW).eq.1) THEN
            CALL PS_PointLog (X, Y, D, n, ZReg(ir), -iRLS(ir),
     *                     rSyMag(iGW), qErrBar(iGW))
         ELSE
            CALL PS_Curve (X, Y, n, ZReg(ir), iRLS(ir), rSyMag(iGW))
            ENDIF

         ENDDO ! ir

C  Text lines :
      i = 0
      iPos = 0
 11   CONTINUE
         i = i+1
         iPos = iPos + 1
         IF (i.gt.iRTS(7)) GOTO 19
         aus = RegTX(i)
 12      IF (i.lt.iRTS(7) .and. RegTX(i+1)(1:3).eq.'&cd') THEN
            i = i + 1
            CALL Append (aus, RegTX(i)(4:len(RegTX(i))))
            GOTO 12
            ENDIF
         CALL PS_Text (aus, ' ', 'infline') ! position 10+i only for TEK
         GOTO 11
 19   CONTINUE


C Fixed date/time output (by Christian Geisler, Sep04)
C modified by FK Oct06 now fully operational (System dependent!!)
      CALL  fdate(Datuma)
      CALL  Compose5 (Datum(1:10),Datuma(9:10),'-',Datuma(5:7),'-',
     *                Datuma(23:24))
      CALL  Compose2 (Zeit(1:9),Datuma(12:19),'-' )
      CALL Compose3 (aus(1:21),Datum(1:10), ', ', Zeit(1:8))
      CALL PS_Text (aus,  ' ', 'infline')

C  Close output file :
      CALL ClosePS ()

      END ! GraSoftCopy

C  ====================================================================
C  Tektronix screen driver
C  ====================================================================

      SUBROUTINE SetTek ()
C     --------------------
         ! set device = TEKTRONIX

      IMPLICIT LOGICAL (q)

      COMMON / GraTerm/ qWindow, qGTOverlay, qGToldAdr
      COMMON / Scroll / iScroIst, iScroGra
      COMMON / Format / kScrF, kLasF, kPS_F

      IF (qGTOverlay .and. iScroIst.ne.iScroGra)
     *    CALL SetScroll (iScroGra)
      CALL TekSetDevice (6)
      CALL SetWindow (kScrF)

      END ! SetTek

      SUBROUTINE SetScroll (iScroNew)
C     -------------------------------
         ! Set the scroll area.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      COMMON / Scroll / iScroIst, iScroGra

      IF (qioutside(iScroNew, 0, 25))
     *    CALL Absturz ('SetScroll', 'iScroNew o.o.r.')

      CALL ScrollArea (iScroNew,25)
      IF (iScroNew.gt.0) THEN
         CALL OverlayOn
      ELSE
         CALL OverlayOff
         ENDIF

      iScroIst = iScroNew

      END ! SetScroll

C  ====================================================================
C  auxiliary / math for line drawing
C  ====================================================================

      SUBROUTINE SectBoxLin (xf1, xf2, yf1, yf2, xa, ya, xb, yb,
     *                       qIn1, qIn2, nSect, XS, YS)
C     ----------------------------------------------------------
            ! JWu 6mar91 in real coordinates, 20mar91 as TekFrameEdge
            ! in graphic coordinates, 20feb92 again in real coordinates,
            ! 30aug93 new as SectBoxLin
         ! for use in GraPaint :
         ! calculates the intersections (xs_,ys_) of the
         ! line (xa,ya)-(xb,yb) with the box (xf_,yf_)

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)
      DIMENSION         XS(4), YS(4) ! at least (3) : qSectLinLin writes
                                     ! to XS even if there is no intersection

      IF (xf1.ge.xf2 .or. yf1.ge.yf2) CALL Absturz (
     *     'SectBoxLin', 'box ill defined')

      qIn1 = qrinside(xa, xf1, xf2) .and. qrinside(ya, yf1, yf2)
      qIn2 = qrinside(xb, xf1, xf2) .and. qrinside(yb, yf1, yf2)
      nSect= 0 ! # intersections

      IF (qIn1 .and. qIn2) RETURN ! points inside > ignore border case

C  Check for intersections with all four lines limiting the box:
      IF (qSectLinLin(xf1, yf1, xf1, yf2, xa, ya, xb, yb,
     *                XS(nSect+1), YS(nSect+1))) nSect=nSect+1
      ! the right order should be conserved :
      ! if ya<yb, then look first at the lower boundary yf1, later at yf2
      IF (ya.lt.yb) THEN
         IF (qSectLinLin(xf1, yf1, xf2, yf1, xa, ya, xb, yb,
     *                   XS(nSect+1), YS(nSect+1))) nSect=nSect+1
         IF (qSectLinLin(xf1, yf2, xf2, yf2, xa, ya, xb, yb,
     *                   XS(nSect+1), YS(nSect+1))) nSect=nSect+1
      ELSE
         IF (qSectLinLin(xf1, yf2, xf2, yf2, xa, ya, xb, yb,
     *                   XS(nSect+1), YS(nSect+1))) nSect=nSect+1
         IF (qSectLinLin(xf1, yf1, xf2, yf1, xa, ya, xb, yb,
     *                   XS(nSect+1), YS(nSect+1))) nSect=nSect+1
         ENDIF
      IF (qSectLinLin(xf2, yf1, xf2, yf2, xa, ya, xb, yb,
     *                XS(nSect+1), YS(nSect+1))) nSect=nSect+1

C  Security checks (to be replaced later by correct arithmetic error handling) :
      IF (nSect.gt.2) CALL Absturz ('SectBoxLin', 'nSect>2')
c      IF ((qIn1.ne.qIn2).and.(nSect.ne.1)) THEN
c    *     CALL Absturz ('SectBoxLin', 'nSect<>1')
      IF (nSect.ge.7) THEN ! f"ur Fehlersuche - au"ser Betrieb
         CALL GMode (0)
         Print *, 'SectBoxLin/ qIn1, qIn2, nSect = ', qIn1, qIn2, nSect
         Print *, ' xa ya = ', xa, ya
         Print *, ' xb yb = ', xb, yb
         DO i = 1, nSect
            Print *, ' XS YS = ', XS(i), YS(i)
            ENDDO
         ENDIF

      END ! SectBoxLin

      LOGICAL FUNCTION qSectLinLin (xA1, yA1, xA2, yA2,
     *                              xB1, yB1, xB2, yB2, xs, ys)
C     ---------------------------------------------------------
            ! JWu 31aug93
         ! Is there an intersection (xs,ys) of two lines A, B ?
         ! First application: (A) are the lines delimiting a Box.
         ! A junction is not counted as an Intersection (21may94)

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      IF     (xA1.eq.xA2) THEN  ! vertical case
         xs = xA1
         IF (qrInOpen (xs, xB1, xB2)) THEN
            CALL LinIntPol (xs, xB1, xB2, yB1, yB2, 0.d0, 0.d0,ys,dum)
            qSectLinLin = qrInOpen (ys, yA1, yA2)
         ELSE
            qSectLinLin = .false.
            ENDIF
      ELSEIF (yA1.eq.yA2) THEN  ! horizontal case
         ys = yA1
         IF (qrInOpen (ys, yB1, yB2)) THEN
            CALL LinIntPol (ys, yB1, yB2, xB1, xB2, 0.d0, 0.d0,xs,dum)
            qSectLinLin = qrInOpen (xs, xA1, xA2)
         ELSE
            qSectLinLin = .false.
            ENDIF
      ELSE                      ! general case
         CALL Absturz ('qSectLinLin', 'gen. case not yet implemented')
         ENDIF

      END ! qSectLinLin

C  ====================================================================
C  auxiliary / compose labels
C  ====================================================================

      SUBROUTINE GraLabel (j, label, var, unit)
C     -----------------------------------------
            ! refait pour la n-ieme fois le 11 mars 91 JWu
         ! returns the label, composed of var and unit
         ! according to the function chosen for axis j.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8 (a-h,o-p,r-z)

      INCLUDE 'g_dim.f'
      CHARACTER *(*) label, var, unit
      CHARACTER *40  text2, hilf

      COMMON / GraWdw / iGW, nGW
      COMMON / GraCS  / iGFu(3,MGW), rGFu(3,MGW), Gmm(3,2,MGW)

      IF (qioutside(j,1,3)) CALL Absturz ('GraLabel', 'axis o.o.r.')

      hilf = var
      IF (unit.ne.' ') CALL Compose4 (hilf, hilf, ' (', unit, ')')

      CALL NiceNum (rGFu(j,iGW), text2, i2)
      CALL FuTxt (iGFu(j,iGW), label, hilf, text2(1:i2))

      END ! GraLabel

C  ====================================================================
C  g2.f / eof
C  ====================================================================
