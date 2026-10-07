C  ====================================================================
C
C      Library  IDA   :  Ingenious Data Analysis
C      Modul    i01   :     plot (interface to g2.f)
C
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

      SUBROUTINE IdaPlot (Word, nJList, JList, Object, Fehler)
C     --------------------------------------------------------
         ! call plot routines to create a frame and to show data
            ! as a subroutine JWu 30jul91.
            ! 12mar96: clear division: first setup, then plot.

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler, Object, Word
      DIMENSION      JList(*)

      INTEGER        K1K(MK), K2K(MK), K3K(MK)
      REAL*8         Zp(MK), Z(MZ)

      CHARACTER*1    action
      CHARACTER*40   Co3(3), Un3(3), CoZ(MZ), UnZ(MZ),
     *               CoP, UnP, Un, file, title, doc,hilf
      CHARACTER*80   aus, lisK, lisKold, line
      CHARACTER*240  tline0, tline1

      DATA   lisK /'*'/

      IF (nJList.eq.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

C  Input commands :
      IF     (Word.eq.'p' ) THEN
         qPlotNew  = .true.
         ndim      = 2
      ELSEIF (Word.eq.'a' ) THEN
         qPlotNew  = .false.
      ELSEIF (Word.eq.'pp' ) THEN
         qPlotNew  = .true.
         ndim      = 3
      ELSE
         Fehler = ' Unknown command : '//Word
         RETURN
         ENDIF

C  New coordinate system ?
      IF (qPlotNew) THEN
         j = JList(1)
         nK = iOlfG (j, '#spectra', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (qCurve) THEN
            ifc = iOlfG (j, 'fu#', Fehler)
            ENDIF
         IF (Fehler.ne.'&ff') RETURN

         ! which spectra ?
         IF (Object.ne.' ') THEN
            CALL DecJList (Object, MK, nK1, K1K, 1, nK, Fehler)
            lisK = Object
            IF (nK1.le.0 .or. Fehler.ne.'&ff') RETURN
         ELSEIF (nK.eq.1) THEN
            K1K(1) = 1
            nK1 = 1
         ELSE
            lisKold = lisK
            CALL Compose2 (aus, ' Plot which of the '//
     *                          cl3(nK), ' spectra')
            CALL GetJList (aus, lisK, MK, nK1, K1K, 1, nK)
            IF (nK1.le.0) THEN
               lisK = lisKold ! old default
               Fehler = ' '
               RETURN
               ENDIF
            ENDIF

         ! get coordinate names & units > determine appropriate window
         CALL OlfCnuG (j, 'x', Co3(1), Un3(1), Fehler)
         CALL OlfCnuG (j, 'y', Co3(2), Un3(2), Fehler)
         IF (ndim.ge.3) CALL OlfCnuG (j, 'z1', Co3(3), Un3(3), Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL GraSetWdw (ndim, Co3, Un3, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! (re)set x y (and z) scales of the chosen window
               ! der Ablauf ist ganz einfach : man ruft die Routine
               ! {GraSetSca}, die ihrerseits "uber {action} eine
               ! bestimmte Reaktion verlangt, die im folgenden
               ! ausgef"uhrt wird und gegebenenfalls nach {GraSetSca}
               ! zur"uckspringt.
         xi = 0
         xf = 0
 22      CALL GraSetSca (1, 'x', xi, xf, action, qLin, igfx, rgfx)
         IF     (action.eq.'-') THEN ! escape
            Fehler = ' '
            RETURN
         ELSEIF (action.eq.'n') THEN ! calculate new defaults
            IF (qCurve) THEN
               IF (qLin) THEN
                  xi = -1
                  xf =  1
               ELSE
                  xi = .1
                  xf = 10
                  ENDIF
            ELSE
               ! search for minimum and maximum
               xi =  1d20
               xf = -1d20
               ! loop spectra
               DO jK1 = 1, nK1
                  K = K1K(jK1)
                  CALL OlfGetX (j, K, n, X, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  DO i = 1, n
                     CALL FuVal (igfx, xx, dummy, X(i), 0d0, rgfx, 0d0)
                     IF (xx.lt.xi .and. (qLin .or. xx.gt.0)) xi = xx
                     xf = dmax1 (xf, xx)
                     ENDDO
                  ENDDO
               IF (qLin) THEN
                  CALL RoundLin (xi, xf, .052d0, 2)
               ELSE
                  CALL RoundLog (xi, xf, .1d0, 1)
                  ENDIF
               ENDIF ! qCurve or not
            IF (xi.ge.xf) THEN ! avoid endless loop 24mrz98
               xi = xi - 1
               xf = xf + 1
               ENDIF
            GOTO 22

         ELSEIF (action.eq.'x') THEN ! revise 'x'
            GOTO 22
            ENDIF

         yi = 0
         yf = 0
 23      CALL GraSetSca (2, 'y', yi, yf, action, qLin, igfy, rgfy)
         IF     (action.eq.'-') THEN ! escape
            Fehler = ' '
            RETURN
         ELSEIF (action.eq.'n') THEN ! calculate new defaults
            yi =  1d49
            yf = -1d49
            ! loop spectra
            DO jK1 = 1, nK1
               K = K1K(jK1)
               IF (qCurve) THEN
                  ! calculate value without convolution (and without par(18))
                  n = 2
                  X(1) = xi
                  X(2) = xf
                  CALL OlfGetY (j, K, ndummy, Y1, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  CALL CuFuVal (ifc, Y1, X, Y, n, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
               ELSE
                  CALL OlfGetXY (j, K, n, X, Y, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               DO i = 1, n
                  CALL FuVal (igfx, xx, dummy, X(i), 0d0, rgfx, 0d0)
                  IF (xi.le.xx .and. xx.le.xf) THEN
                     CALL FuVal (igfy, yy, dummy, Y(i), 0d0, rgfy, 0d0)
                     IF (yy.lt.yi .and. (qLin .or. yy.gt.0)) yi = yy
                     yf = dmax1 (yf, yy)
                     ENDIF
                  ENDDO
               ENDDO
            IF ((.not.qLin .and. yi.le.0 .or. yf.le.yi)) THEN ! kommt vor
               CALL Gong (2)
               Print *, 'curve starts or ends below 0: yi yf = ',yi,yf
               yi = 1d-10
               yf = 1d10
               ENDIF
            IF (qLin) THEN
               CALL RoundLin (yi, yf, .052d0, 2)
            ELSE
               CALL RoundLog (yi, yf, .1d0, 1)
               ENDIF
            IF (yi.ge.yf) THEN ! avoid endless loop 12nov96
               yi = yi - 1
               yf = yf + 1
               ENDIF
            GOTO 23

         ELSEIF (action.eq.'x') THEN ! revise 'x'
            GOTO 22
            ENDIF

         IF (ndim.lt.3) GOTO 25
         zi = 0
         zf = 0
 24      CALL GraSetSca (3, 'z', zi, zf, action, qLin, igfz, rgfz)
         IF     (action.eq.'-') THEN ! escape
            Fehler = ' '
            RETURN
         ELSEIF (action.eq.'n') THEN ! calculate new defaults
            CALL OlfGet1ZofK (j, 1, nKz, Zp, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            zi =  1d20
            zf = -1d20
            DO K = 1, nKz
               CALL FuVal (igfz, z1, dummy, Zp(K), 0d0, rgfz, 0d0)
               IF (z1.lt.zi .and. (qLin .or. z1.gt.0)) zi = z1
               zf = dmax1 (zf, z1)
               ENDDO
            IF (qLin) THEN
               CALL RoundLin (zi, zf, .052d0, 2)
            ELSE
               CALL RoundLog (zi, zf, .1d0, 1)
               ENDIF
            GOTO 24

         ELSEIF (action.eq.'x') THEN ! revise 'x'
            GOTO 22

            ENDIF

 25      CONTINUE

         ! Now install the chosen axes and proceed to plot them:
         CALL GraSetAxs (Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL GraPlotCS (Fehler)
         IF (Fehler.ne.'&ff') RETURN
         iK = 0

         ENDIF ! qPlotNew

C  Proceed to plot (loop files / spectra) :
      DO lj = 1, nJList
         j = JList(lj)
         qFileNew = .true.

         nK = iOlfG (j, '#spectra', Fehler)
         CALL tOlfG (j, 'fil', file,  Fehler)
         CALL tOlfG (j, 'tit', title, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         qCurve = qOlfGdef (j, '?cu', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qCurve) THEN
            ifc     = iOlfG (j, 'fu#', Fehler)
            qConv   = qOlfG (j, '?conv', Fehler)
            jparin  = iOlfGdef (j, 'fit-par-file#', 0, Fehler)
            nplopts = iOlfGdef (j, 'plot-#pts', 100, Fehler)
            ri      = rOlfGdef (j, 'plot-i', Un, 0.d0, Fehler)
            rf      = rOlfGdef (j, 'plot-f', Un, 0.d0, Fehler)
            doc = '96/7'
            IF (Fehler.ne.'&ff') RETURN
            IF (qConv) THEN ! possibly with convolution
               CALL GraQuit ()
               j2 = jCuConvAsk (qCurve, qConv)
               IF (j2.ne.0) THEN
                  CALL GetK2K (j, j2, nK, K2K, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
            ELSE
               j2 = 0
               ENDIF
            IF (jparin.ne.0) THEN
               jpar = iAskDMu (' Parameter file', jparin, 0, MF)
               IF (jpar.le.0) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               IF (jpar.ne.jparin) THEN
                  CALL OlfOpen (j, 1, Fehler)
                  CALL iOlfP (j, 'fit-par-file#', jpar, Fehler)
                  CALL OlfClos (j, nK, Fehler)
                  ENDIF
               CALL GetK2K (j, jpar, nK, K3K, Fehler)
               IF (Fehler.ne.'&ff') RETURN
            ELSE
               jpar = 0 !Artem: Otherwise, the values may be random
               K3K = 0 !Artem: Otherwise, the values may be random
               ENDIF
            ENDIF ! qCurve

         IF (.not.qPlotNew .or. lj.gt.1) THEN
            ! which spectra :
            IF (Object.ne.' ') THEN
               CALL DecJList (Object, MK, nK1, K1K, 1, nK, Fehler)
               lisK = Object
               IF (nK1.le.0 .or. Fehler.ne.'&ff') RETURN
            ELSEIF (nK.eq.1) THEN
               K1K(1) = 1
               nK1 = 1
            ELSE
               lisKold = lisK
               CALL GraQuit ()
               CALL Compose2 (aus, ' Plot which of the '//
     *                             cl3(nK), ' spectra')
               CALL GetJList (aus, lisK, MK, nK1, K1K, 1, nK)
               IF (nK1.le.0) THEN
                  lisK = lisKold ! old default
                  Fehler = ' '
                  RETURN
                  ENDIF
               ENDIF
            ENDIF

         ! loop spectra :
         isy = iOlfGdef (j, 'plot-sy#', 0, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         DO jK1 = 1, nK1
            ! to avoid excessive waiting :
            IF (isy.le.0 .and. jK1.eq.2 .and. nK1.gt.18 .and.
     *          ndim.ne.3) THEN
               CALL GraQuit ()
               CALL Compose2 (aus, ' Plot '//
     *                        cl3(nK1-1), ' more spectra')
               IF (.not.qAskD(aus, 1)) THEN
                  Fehler = ' '
                  RETURN
                  ENDIF
               ENDIF
            K = K1K(jK1)
            CALL OlfGetXYD (j, K, nC, X, Y, D, Fehler)
            CALL OlfGetZ (j, K, nZ, Z, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (nC.eq.0) THEN
               Fehler = ' spectrum is empty'
               RETURN
               ENDIF

            IF (qCurve) THEN ! plot curve
               ! x-range :
               CALL GraInquireCS ('x', xmi, xma, igfx, gfpx, liloX)
               IF (ri.ne.rf) THEN
                  xmi = dmax1 (xmi, ri)
                  xma = dmin1 (xma, rf)
                  ENDIF
               ! y-stepwidth :
               CALL GraInquireCS ('y', ymi, yma, igfy, gfpy, liloY)
               relY = 6.d0 / nplopts ! maximum y-step
               ! do everything in i6.f :

               CALL CuNice (ifc, Y, j2, K2K(K), jpar, K3K(K),
     *                      xmi, xma, igfx, gfpx, liloX,
     *                      ymi, yma, igfy, gfpy, liloY,
     *                      nplopts, relY,
     *                      MC, nC1, X1, Y1, D1, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL GraPaint (X1, Y1, D1, nC1, Z(1), isy, iLS)

            ELSE ! plot data
               iK = iK + 1
               IF (iK.ge.100) iK = 1 ! to avoid overflow of &sy
               IF (isy.eq.0) THEN ! default
                  IF (ndim.eq.3) THEN
                     iLSdef =  1   ! always straight lines
                  ELSE
                     iLSdef = -iK  ! incremental plot symbols
                     ENDIF
               ELSE ! plotsymbol is file-parameter
                  iLSdef = isy
                  ENDIF
               CALL GraPaint (X, Y, D, nC, Z(1), iLSdef, iLS)
               ENDIF

            IF (iLS.eq.0) RETURN ! user's choice : draw nothing

            IF (Fehler.ne.'&ff') RETURN
            DO iZ = 1, nZ
               CALL OlfCnuG (j, 'z'//cl2(iZ), CoZ(iZ), UnZ(iZ), Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDDO

            ! == display some text lines on the screen == !
            IF (qFileNew) THEN ! display name of file
               CALL GraText (file, 1)
               CALL GraText (' ', 1)
               ENDIF
            IF (ndim.eq.2) THEN ! some info for each spectrum
               aus = cr3(j)//'#'//cl3(K)
               IF (iLS.lt.0) THEN
                  CALL GraText ('&sy='//cv2(-iLS)//aus, 1)
               ELSEIF (iLS.gt.0) THEN
                  CALL GraText ('&li='//cv2( iLS)//aus, 1)
                  ENDIF
               IF (nK.gt.1) THEN ! display z
                  DO iZ = 1, nZ
                     CALL NiceNum (Z(iZ), hilf, ia)
                     aus = CoZ(iZ)
                     CALL Append (aus, ' = '//hilf(1:ia)//' '//UnZ(iZ))
                     CALL GraText (aus, 1)
                     ENDDO
                  ENDIF
               ENDIF
            CALL GraText (' ', 1) ! empty line

            ! == text for register == !
            ! head line : fname, doc, title
            IF (qFileNew) THEN
               CALL Compose3 (aus, file, ',  '//title, '    '//doc)
               CALL GraText (aus, 2)
               ENDIF
            ! encode graphic symbol :
            IF (iLS.lt.0) THEN
               aus = '&sy='//cv2(-iLs )
            ELSE
               aus = '&li='//cv2( iLS )
               ENDIF
            ! list of spectra, z, parameters if curve :
            IF (qCurve) THEN ! print parameters (11sep91, wurde auch Zeit)
               tline0 = ' ' ! head line (1/file)
               tline1 = ' ' ! data line (1/spectrum)
               ! First 3 columns : linestyle, K, z :
               tline1(2:9) = aus ! linestyle
               IF (nK.gt.1) THEN ! K, z
                  tline1(10:12) = cl3(K)
                  Write(tline1(13:24), '(g12.4e2)') Z(1)
                  CALL Compose3 (aus, '# '//CoZ(1), '('//UnZ(1), ')')
                  tline0(15:80) = aus ! (1:14) remain empty
                  jipar = 24
                  nCoutmax = 18
               ELSE
                  jipar = 12
                  nCoutmax = 19
                  ENDIF
               ! Following columns : curve-parameters
               nCout = min0 (nCoutmax, nC) ! # parameters that can be printed
               DO i = 1, nCout
                  ! could be fine-tuned (more place if nCout<5)
                  ji = jipar + i*12-11
                  jf = jipar + i*12
                  IF (qFileNew) THEN ! names of parameters
c                     IF (qintr(i-Par(20+i))) THEN
c                        tline0(ji:ji) = '='  ! fix
c                     ELSE
c                        tline0(ji:ji) = '!'  ! free
c                        ENDIF
                     CALL OlfCnuG (j, 'p'//cl2(i), CoP, UnP, Fehler)
                     tline0(ji+1:jf) = CoP
                     ENDIF
                  Write (hilf,'(g12.4e2)') Y(i)
                  tline1(ji:jf) = hilf
               ENDDO
               IF (qFileNew) CALL GraText (tline0, 2)
               CALL GraText (tline1, 2)
            ELSE ! no curve
               IF (nK.gt.1) THEN ! name and value of z
                  CALL Append (aus, '# '//cl3(K))
                  CALL Append (aus, ' :')
                  DO iZ = 1, nZ
                     CALL NiceNum (Z(iZ), hilf, ia)
                     IF (iZ.gt.1) CALL Append (aus, ',')
                     CALL Append (aus, ' '//CoZ(iZ))
                     CALL Append (aus, ' = '//hilf(1:ia)//' '//UnZ(iZ))
                     ENDDO
                  ENDIF
               CALL GraText (aus, 2)
               ENDIF ! qCurve

            qFileNew = .false.
 629        CONTINUE
            ENDDO ! loop spectra

            ! == file-info for register == !
            i = 1
 71         CONTINUE
               CALL OlfComLinG (j, i, line, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               IF (line.eq.'&eoc') GOTO 79
               CALL GraText ('  '//line(1:lenU(line)), 2) ! comment
               i = i + 1
               GOTO 71
 79         CONTINUE

         ENDDO ! loop files

      CALL GraQuit ()

      END ! IdaPlot
