C  ====================================================================
C
C      Library  IDA   :  Inelastic data treatment
C      Modul    i74   :     calculate self absorption coefficients
C
C  ====================================================================

C      Contents :
C         AbsCoeffs, AskGeometry, CoCaSlab, CoCaCyl, AbsCyl, TraCyl

C      History :
C         9nov96   2th-w-grid just as given
C         2sep96   recovered from sqw1
C         22jun94   geometry parameters in output comment lines
C         7apr93   more comments, SQW1 corrected
C         12. 8.91, comments, units, renamed COCA.
C         16. 5.91, as test program ACT
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE AbsCoeffs (nJList, JList, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'l_def.f'
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INTEGER       nJList, JList(*)
      CHARACTER*40  File, FileOut
      CHARACTER*(*) Fehler
      CHARACTER*80  ein, aus
      CHARACTER*40  doc(2)
      LOGICAL       qCont

      DATA          nT / 40 /, nR / 20 / ! had been ill-adjusted in A1,38ff
                                         ! does no longer depend on size
      IF (nJlist.eq.0) THEN
         Fehler = ' '
         RETURN
         ENDIF

      DO lj = 1, nJList
         j = JList(lj)

         ! Description of sample cell :
         CALL AskGeometry (wB, cB, jGeo, back, SamD, alpha, SamR,
     *                     ConT1, ConD1, ConT2, ConD2, qCont,
     *                     abs_s, abs_a, doc, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (jGeo.eq.2 .or. jGeo.eq.3) THEN
            nT = iAskD ('# angular meshes', nT)
            nR = iAskD ('# radial meshes', nR)
            ENDIF

         IF (abs_a.ne.0) THEN
            Ei = rOlfGG (j, 'E0', 'meV', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            absi = abs_s + dsqrt (25.305/Ei) * abs_a ! 25.3meV = 2200m/sec
         ELSE
            absi = abs_s
            ENDIF

         ! Output files :
         CALL OlfHeadDup (j, .false., j1, nK, K1, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (j1, 'y', 'Assc', ' ', Fehler)
         CALL tOlfG (j, 'fil', File, Fehler)
         CALL Compose2 (FileOut, File, '_A1')
         CALL tOlfP (j1, 'fil', FileOut, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (qCont) THEN
            CALL OlfHeadDup (j, .false., j2, nK, K2, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            CALL OlfCnuP (j2, 'y', 'Arel', ' ', Fehler)
            CALL Compose2 (FileOut, File, '_A2')
            CALL tOlfP (j2, 'fil', FileOut, Fehler)
            ENDIF

         ! the calculation :
         DO K = 1, nK
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)

            DO i = 1, n
               tt = rzxOlfGG (j, K, i, '2th', ' ', Fehler)

               IF (abs_a.ne.0) THEN
                  Ef = Ei + rzxOlfGG (j, K, i, 'w', 'meV', Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  absf = abs_s + dsqrt (25.305/Ef) * abs_a
               ELSE
                  absf = abs_s
                  ENDIF

               IF     (jGeo.eq.1) THEN
                  CALL CoCaSlab (tt, back, SamD, alpha, ConT1, ConT2,
     *                           absi, absf, Assc, Arel, Fehler)
               ELSEIF (jGeo.eq.2 .or. jGeo.eq.3) THEN
                  CALL CoCaCyl (tt, wB, cB, jGeo, back, SamD, SamR,
     *                          ConT1, ConD1, ConT2, ConD2,
     *                          absi, absf, Assc, Arel, trans,
     *                          nT, nR, area, Fehler)
               ELSE
                  Fehler = 'PROG ERR/ bad j'
                  RETURN
                  ENDIF
               IF (Fehler.ne.'&ff') RETURN

               Y1(i) = Assc
               Y2(i) = Arel
               D(i)  = 0

               ENDDO
c            CALL Counter (K, 1, nK, '... working hard')

            CALL OlfCopZ (j, j1, K, K, Fehler)
            CALL OlfPutXY0 (j1, K, n, X, Y1, Fehler)
            IF (qCont) THEN
               CALL OlfCopZ (j, j2, K, K, Fehler)
               CALL OlfPutXY0 (j2, K, n, X, Y2, Fehler)
               ENDIF
            IF (Fehler.ne.'&ff') RETURN
            ENDDO ! K

         CALL OlfClos (j1, nK, Fehler)
         IF (qCont) CALL OlfClos (j2, nK, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         Print *, ' target section area ', area
         Print *, ' transmission        ', trans, ' IS WRONG if wB>0'

         ENDDO ! lj

      END ! AbsCoeffs

C  --------------------------------------------------------------------
      SUBROUTINE AskGeometry (wB, cB, jGeo, back, d, alpha, R,
     *                        t1, d1, t2, d2, qCont,
     *                        abs_s, abs_a, doc, Fehler)
C  --------------------------------------------------------------------
                                                     ! JWuttke, May 1991

         ! Questionary for geometry and scattering intensity.
         ! For absorption corrections, one needs the product abso
         ! and complete information on the experimental geometry.

         ! Output :
         !  wB    = beam width
         !  cB    = beam centre
         !  jGeo  = code for geometry
         !  back  = ratio of back-scattered beam passing through the sample
         !  d     = thickness of slab or hollow cylinder
         !  alpha = angle (2 theta) of slab
         !  R     = outer radius of cylinder
         !  t1    = transmission of frontside/inner container (C1)
         !  d1    = thickness of C1
         !  t2    = transmission of backside/outer container (C2)
         !  d2    = thickness of C2
         !  qCont = is there a container ?
         !  abs_s = (number density) * (scattering(i+c) cross section)
         !  abs_a = (number density) * (absorption cross section at 2200m/sec)

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      CHARACTER        aus*80, Fehler*(*)
      CHARACTER*(*)    doc(2)

 1    CONTINUE

      wB = rAskD ('Width of beam (0:fully illuminated)', wB)
      IF (wB.gt.0) cB = rAskD ('And its center', cB)

C  Which kind of geometry ?
      Print *, 'Geometry : ' ! (0) unspecified'
      Print *, '           (1) flat slab,'
      Print *, '           (2) solid cylinder,'
      Print *, '           (3) hollow cylinder,'
      jGeo = iAskDMu ('Case', jGeo, 0, 3)

C  Absorption correction ?
      IF     (jGeo.eq.0) THEN
         Fehler = ' '
         RETURN ! 'no absorption correction for unspecified geometry.'
         ENDIF

C  Information on geometry :
      qCont = .true.
      IF (jGeo.eq.1) THEN
C  Enter slab parameters :
         d    = rAskDMu ('Thickness of the slab (cm)', d, 1.d-4, 1.d1)
 111     CONTINUE
         alpha= rAskDMu
     *      ('Angle between plane of the slab and transmitted beam',
     *        alpha, 0.d0, 180.d0)
         IF (alpha.lt.3. .or. alpha.gt.177.) THEN
            CALL Gong (3)
            Print *, '  Angles should be given in degrees.'
            Print *,
     *   '  Angles outside 3...177 are rejected for physical reasons.'
            GOTO 111
            ENDIF

         t = rAskDMu (
     * 'Container transmission (for perpendicular beam)',
     * t, 1.d-3, 1.d0)
            ! frontside and backside not different
         t1   = dsqrt (t)
         t2   = t1
      ELSEIF (jGeo.eq.2) THEN
C  Enter solid cylinder parameters :
         R    = rAskMu ('Radius of the sample (cm) ?', 1.d-4, 1.d0)
         t2   = rAskMu ('Container transmission ?', 1.d-3, 1.d0)
         IF (t2.lt.1.) THEN
            d2 = rAskMu ('Thickness of container ?', 1.d-6, 1.d1)
         ELSE
            CALL Gong (1)
            Print *, '  WARNING :'
            Print *, '  The relative correction factor'
     *           //' (container/sample) will be set to A2/A1 = 1'
            d2 = 0.
            ENDIF
C  In the following, the solid cylinder will be treated as
C  a special hollow cylinder :
         d  = R
         t1 = 1.
         d1 = 0.
      ELSEIF (jGeo.eq.3) THEN
         R  =
     * rAskDMu ('Outer radius of the sample (cm) ?', R, 1.d-3, 1.d3)
         d  =
     * rAskDMu ('Thickness of the sample layer (cm) ?', d, 1.d-4, R)
         t1 =
     * rAskDMu ('Transmission of the inner container?',
     *          t1, 1.d-3, 1.d0)
         IF (t1.lt.1.) THEN
            d1 = rAskMu ('Thickness of the inner container (cm) ?',
     *             1.d-6, 1.d1)
         ELSE
            d1 = 0
            ENDIF
         t2 =
     * rAskDMu ('Transmission of the outer container?',
     *          t2, 1.d-3, 1.d0)
         IF (t2.lt.1.) THEN
            d2 = rAskMu ('Thickness of the outer container (cm) ?',
     *             1.d-6, 1.d1)
         ELSE
            d2 = 0
            ENDIF
         IF (t1.eq.1 .and. t2.eq.1) qCont = .false.
         ENDIF

      back = rAskDMu (
     * 'Fraction of back-scattered beam transversing the sample',
     * back, 0.d0, 1.d0)

C  Information on the sample :
      dens21 = rAskDMu ('Number density of scatterers (10^21 cm-3)',
     *                  dens21, 1.d-4, 1.d8)
      dens = 1.d-3 * dens21  ! density in 10^24 cm-3
      sigscat = rAskDMu ('Scattering cross section (coh+incoh) (barn)',
     *                   sigscat, 0.d0, 1.d5)
      sigabs = rAskDMu (
     * 'Absorption cross section at 2200 m/sec (barn)',
     * sigabs, 0.d0, 1.d5)

      IF (sigabs+sigscat.le.0) THEN
         Fehler = 'no scattering, no absorption -> no coefficients'
         RETURN
         ENDIF

      abs_s = dens * sigscat
      abs_a = dens * sigabs

C  Documentation :
      IF     (jGeo.eq.0) THEN
         aus = '   geometry unspecified'
      ELSEIF (jGeo.eq.1) THEN
            write (aus, '(a,f5.2,a,f5.1,a,f5.2,a,f5.2)')
     * '  slab : d=', d, 'cm alpha=', alpha, ' t_C=', t, ' fB=', back
      ELSEIF (jGeo.eq.2) THEN
         aus = '* cyl'
      ELSEIF (jGeo.eq.3) THEN
            write (aus, '(7(a,f5.2))')
     *       '  cyl : d=', d*10, 'mm Ra=', R*10,
     *       'mm ti=', t1, ' ta=', t2,
     *       ' di=', d1*10, 'mm da=', d2*10, 'mm fB=', back
         ENDIF
      doc(1) = aus
      Print *, aus

      write (aus, '(a,e9.2,a,f7.2,a,f7.2,a)')
     *    '          n=', dens*1.d24, 'cm-3 sig_s=', sigscat,
     *    'barn sig_a=', sigabs, 'barn'
      doc(2) = aus
      Print *, aus

      IF (.not.qAskD('All right', 1)) GOTO 1

      END ! AskGeometry

C  --------------------------------------------------------------------
      SUBROUTINE CoCaSlab (angle, back, d, alpha, t1, t2, absi, absf,
     *                    Assc, Arel, Fehler)
C  --------------------------------------------------------------------

         ! calculates the two absorption correction coefficients A1, A2

         ! angle = scattering angle (2 theta)
         ! ...  see above (AskGeometry)
         ! absi  = mu_i = (density of scatterers) * (total cross section at Ei)
         ! absf  = mu_f = (density of scatterers) * (total cross section at Ef)
              ! for back-scattering, we have always absi = absf.
              ! however, both quantities are kept separated for
              ! the sake of conceptual clarity.

      IMPLICIT NONE
      CHARACTER       Fehler*(*)
      REAL*8          angle, back, d, alpha, t1, t2, absi, absf, Assc,
     *                Arel, beta, angmin, Ass, xsi, xsf, tsi, tsf,
     *                SigC, f, Bs, Bc, Bsc, tci, tcf

      angmin = 1.d-2

      beta = alpha - angle

      IF (dabs(beta).lt.angmin .or.  ! scattering in slab direction
     *    alpha.lt.angmin .or. alpha.gt.180-angmin) THEN
         Assc = 0.
         Arel = 0.
         RETURN
         ENDIF

C  - absorption by sample :
      xsi = absi * d / dsind(alpha)
      xsf = absf * d / dsind(beta)
      tsi = dexp ( - xsi )      ! transmission of incident beam through sample
      tsf = dexp ( -dabs(xsf) ) ! transmission of scattered beam through sample

C  - integral : scattering by sample, absorption by sample :
      IF (dabs(xsi-xsf).lt.1.d-8) THEN ! to avoid division by 0.
         IF (beta.gt.0.) THEN
            Ass = dexp(-xsi) * (1 - .5*(xsf-xsi) )
         ELSE
            Ass = 1 + .5*(xsf-xsi) ! of course the second summand is negative
            ENDIF
      ELSE
         IF (beta.gt.0.) THEN   ! transmission case :
            Ass = (dexp(-xsi) - dexp(-xsf)) / (xsf - xsi)
         ELSE                   ! reflection case :
            Ass = (dexp(xsf-xsi) - 1) / (xsf - xsi)
            ENDIF
         ENDIF

C  - absorption by one of the two container walls :
      IF (t1.ne.t2) THEN
         Fehler = ' PROGRAM ERROR / Transmission different for C1, C2'
         RETURN
         ENDIF
      SigC = - dlog (t1)
      tci  = dexp ( - SigC / dsind(alpha) )
      tcf  = dexp ( - SigC / dabs(dsind(beta)) )

C  - absorption of the back-scattered beam :

      f = back * dabs(dsind(beta))
         ! fraction that passes a second time through the sample

      Bsc = (1 - f) + f * tcf * tsf * tcf  ! transmission through C and S
      Bc  = (1 - f) + f * tcf    *    tcf  ! .. through C alone
      Bs  = Bsc / Bc

C  - now calculate the absorption coefficients :
      ! A(S,SC) :
      Assc = Ass * tci * tcf * Bsc

      ! Arel = A(C,SC) / A(C,C) :
      IF (beta.gt.0) THEN   ! transmission case
         Arel = (tsf*tcf + tci*tsi) / (tcf + tci) * Bs
      ELSE                  ! reflection case
         Arel = (1 + tci*tsi*tsf*tcf) / (1 + tci*tcf) * Bs
         ENDIF

      END ! CoCaSlab

C  --------------------------------------------------------------------
      SUBROUTINE CoCaCyl (angle, wB, cB, jGeo, back,
     *                    d, R, t1, d1, t2, d2, absi, absf,
     *                    Assc, Arel, trans, nT, nR, area, Fehler)
C  --------------------------------------------------------------------

         ! calculates the two absorption correction coefficients A1, A2

         ! angle = scattering angle (2 theta)
         ! ...  see above (AskGeometry)
         ! absi  = mu_i = (density of scatterers) * (total cross section at Ei)
         ! absf  = mu_f = (density of scatterers) * (total cross section at Ef)
              ! for back-scattering, we have always absi = absf.
              ! however, both quantities are kept separated for
              ! the sake of conceptual clarity.
         ! A1    = multiplicator for N(S+C)
         ! A2    = multiplicator for N(C)
         ! trans = transmission
         ! area  = mean sample thickness

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      REAL*8            Rad(4), Rabsi(4), Rabsf(4)
      CHARACTER         Fehler*(*)

C  Cylinder correction :

C  - Set arrays :
      Rad (1) = R - d - d1 ! air
      Rad (2) = R - d      ! inner container (C1)
      IF (d1.gt.0.) THEN
         Rabsi(2) = -dlog(t1) / (2 * d1)
         Rabsf(2) = -dlog(t1) / (2 * d1)
         n1 = 2
      ELSE
            ! nota bene : d1=0 is interpreted as t1=1
         n1 = 3 ! no inner container cylinder
         ENDIF
      Rad  (3) = R          ! sample (S)
      Rabsi(3) = absi
      Rabsf(3) = absf
      Rad  (4) = R + d2     ! outer container (C2)
      IF (d2.gt.0.) THEN
         Rabsi(4) = -dlog(t2) / (2 * d2)
         Rabsf(4) = -dlog(t2) / (2 * d2)
         nU = 4
      ELSE
         nU = 3
         ENDIF

C  - Calculate A(S,SC) :
      CALL AbsCyl (angle, wB, cB, n1, nU, 3,
     *             Rad, Rabsi, Rabsf, nT, nR, ts, area)
      Assc = ts / area

      CALL TraCyl (n1, nU, Rad, Rabsf, R, trans)
      IF (back.gt.0.) THEN
         Assc = Assc * ( (1-back) + back * trans )
         ENDIF

C  - Calculate Arel :
      IF (n1.eq.nU) THEN
         ! both containers are transparent
         Arel = 1.
      ELSE
C  - - To calculate A(C,SC), integrate over C1 and C2 :
         Acsc = 0.
         IF (n1.eq.2) THEN
            CALL AbsCyl (angle, wB, cB, n1, nU, 2,
     *                   Rad, Rabsi, Rabsf, nT, nR, ts, area)
            Acsc = Acsc + ts
            ENDIF
         IF (nU.eq.4) THEN
            CALL AbsCyl (angle, wB, cB, n1, nU, 4,
     *                   Rad, Rabsi, Rabsf, nT, nR, ts, area)
            Acsc = Acsc + ts
            ENDIF
               ! no normalization by area because only the
               ! quotient Acsc/Acc is required
         IF (back.gt.0.) THEN
            Acsc = Acsc * ( (1-back) + back * tb )
               ! the absorption coefficient for the back-scattered
               ! beam is the same for Assc and Acsc.
            ENDIF
C  - - To calculate A(C,C), do the same with air in the sample cylinder :
         Rabsi(3) = 0.
         Rabsf(3) = 0.
         Acc = 0.
         IF (n1.eq.2) THEN
            CALL AbsCyl (angle, wB, cB, n1, nU, 2,
     *                   Rad, Rabsi, Rabsf, nT, nR, ts, area)
            Acc = Acc + ts
            ENDIF
         IF (nU.eq.4) THEN
            CALL AbsCyl (angle, wB, cB, n1, nU, 4,
     *                   Rad, Rabsi, Rabsf, nT, nR, ts, area)
            Acc = Acc + ts
            ENDIF
         IF (back.gt.0.) THEN
            CALL TraCyl (n1, nU, Rad, Rabsf, R, tb)
            Acc = Acc * ( (1-back) + back * tb )
            ENDIF

         Arel = Acsc / Acc
         ENDIF

      END ! CoCaCyl

C  --------------------------------------------------------------------
      SUBROUTINE AbsCyl (ang, wB, cB, n1, nU, nInt, Rad, Absi, Absf,
     *                   nT, nR, ts, area)
C  --------------------------------------------------------------------
                                     ! F.Rieutord jan91, J.Wuttke may91

         ! Determination of the self shielding factor
         ! for an assembly of concentrical cylinders.

         ! Input :
         !  ang  = scattering angle (2 theta)
         !  wB, cB = beam width and center
         !  n1   = first cylinder to be taken into account
         !  nU   = last cylinder ...
                      ! the cylinder are counted from the center
         !  nInt = cylinder in which the beam is scattered
                      ! nInt = 1 is forbidden - for a solid cylinder
                      ! let Rad(1) = 0 and use nInt = 2
         !  Rad  = outer radii of the cylinders
         !  Absi = density * cross section (at Ei)
         !  Absf = dito (at Ef)
         !  nT   = # angular points
         !  nR   = # radial points
         ! Output :
         !  ts   = transmission for the scattered beam (scattered by 2theta)
         !  area = illuminated arealine

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      ! input :
      DIMENSION        Rad(4), Absi(4), Absf(4)
      ! internal :
      DIMENSION        Rad2(4)

C  Auxiliary stores for acceleration of the integration :
      DO i = 1, nU
         Rad2 (i) = Rad(i)**2
         ENDDO

C  Set integration step width :
      d = Rad(nInt) - Rad(nInt-1)
      Rstep = d / nR
      Tstep = 360.d0 / nT

C  Start integration :
      ts    = 0.
      area  = 0.

C  Loop over angle T :
         ! T is defined the same way as the scattering angle (2 theta)
      DO iT = 1, nT
         T = (iT-.5) * Tstep
         ! calculate as much as possible outside the inner loop :
         sTi2 = dsind(T)**2
         cTi  = dcosd(T)
         sTf2 = dsind(T+180-ang)**2
         cTf  = dcosd(T+180-ang)
         IF (wB.gt.0 .and.
     * ((Rad(nInt-1)*dsind(T).lt.cB-wB/2 .and.
     *   Rad(nInt)*dsind(T).lt.cB-wB/2) .or.
     *  (Rad(nInt-1)*dsind(T).gt.cB+wB/2 .and.
     *   Rad(nInt)*dsind(T).gt.cB+wB/2)))
     *      GOTO 799

C  Loop over radius R :
         R = Rad(nInt-1) - .5 * Rstep
         DO iR = 1, nR
            R  = R + Rstep
            IF (wB.gt.0 .and.
     * (R*dsind(T).lt.cB-wB/2 .or. R*dsind(T).gt.cB+wB/2)) GOTO 79
                      ! outside beam
            R2 = R**2

C  Incident beam :
            ySi2 =  R2 * sTi2
            xSi  =  R  * cTi
               ! scattering at point (xSi,ySi);
               ! coordinate system : x in direction of incident beam.
            Sig = 0.
C  - first, the way to |xSi| :
            DO i = nInt+1, nU
               xa  = dsqrt0 (Rad2(i)   - ySi2)
               xi  = dsqrt0 (Rad2(i-1) - ySi2)
               Sig = Sig + (xa-xi) * Absi(i)
               ENDDO
            xa  = dsqrt0 (Rad2(nInt) - ySi2)
            Sig = Sig + (xa-dabs(xSi)) * Absi(nInt)
C  - then, if necessary, the way from -xSi to xSi :
            IF (xSi.gt.0) THEN
               DO i = n1, nInt-1
                  xa  = dsqrt0 (Rad2(i)   - ySi2)
                  xi  = dsqrt0 (Rad2(i-1) - ySi2)
                  Sig = Sig + 2 * (xa-xi) * Absi(i)
                  ENDDO
               xi  = dsqrt0 (Rad2(nInt-1) - ySi2)
               Sig = Sig + 2 * (xSi-xi) * Absi(nInt)
               ENDIF

C  Scattered beam :
            ySf2 =  R2 * sTf2
            xSf  =  R  * cTf
               ! Just the same procedure in the coordinate system
               ! of the scattered beam.
            DO i = nInt+1, nU
               xa  = dsqrt0 (Rad2(i)   - ySf2)
               xi  = dsqrt0 (Rad2(i-1) - ySf2)
               Sig = Sig + (xa-xi) * Absf(i)
               ENDDO
            xa  = dsqrt0 (Rad2(nInt) - ySf2)
            Sig = Sig + (xa-dabs(xSf)) * Absf(nInt)
            IF (xSf.gt.0) THEN
               DO i = n1, nInt-1
                  xa  = dsqrt0 (Rad2(i)   - ySf2)
                  xi  = dsqrt0 (Rad2(i-1) - ySf2)
                  Sig = Sig + 2 * (xa-xi) * Absf(i)
                  ENDDO
               xi  = dsqrt0 (Rad2(nInt-1) - ySf2)
               Sig = Sig + 2 * (xSf-xi) * Absf(nInt)
               ENDIF

C  Now the new element is added to the integral :
            ts   = ts   + R * dexp(-Sig)
            area = area + R
               ! factor dy = R |cosT| dT = xSi
               ! is the non-constant part of the integral measure
               ! (wrongly corrected 7apr93, eliminated 14jan00)
 79         CONTINUE
            ENDDO ! R
 799        CONTINUE
         ENDDO ! T

C  Yet to come : the constant part of the integral measure :
      ts   = ts   * Rstep * (Tstep*3.14159/180)
      area = area * Rstep * (Tstep*3.14159/180)
      IF (area.le.0) CALL Absturz ('AbsCyl', 'area .le. 0')

      END ! AbsCyl

C  --------------------------------------------------------------------
      SUBROUTINE TraCyl (n1, nU, Rad, Abs, B, tr)
C  --------------------------------------------------------------------
                                                     ! J.Wuttke 18may91

         ! Determination of the transmission
         ! for an assembly of concentrical cylinders.

         ! Input :
         !  n1   = Innerster Zylinder to be taken into account
         !  nU   = Aeusserster Zylinder ...
         !  Rad  = outer radii of the cylinders
         !  Abs  = density * cross section (at Ei)
         !  B    : integration from y=-B to y=+B
                     ! by symmetry, one has to integrate only from y=0 to y=B
         ! Output :
         !  tr   = transmission

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      DIMENSION        Rad(1), Abs(1)

C  Arbitrary choice of integration step width :
      ny    = 100
      dy    = B / ny

C  Start integration :
      tr    = 0.

C  Loop over y = distance of beam from the center :
      DO iy = 1, ny
         y = (iy-.5) * dy

         Sig = 0.
         DO i = n1, nU
            xa   = dsqrt0 (Rad(i)  **2 - y**2)
            xi   = dsqrt0 (Rad(i-1)**2 - y**2)
            Sig  = Sig + (xa-xi) * Abs(i)
            ENDDO

         tr = tr + dexp (-2*Sig)

         ENDDO

      tr = tr / ny

      END ! TraCyl
