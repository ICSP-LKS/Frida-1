C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i95   :     simulation for optics
C
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks  

C  --------------------------------------------------------------------
      SUBROUTINE OptikFPI (nJList, JList, qOv, Fehler)
C  --------------------------------------------------------------------
         ! Dokumentation -> I1,86

      IMPLICIT NONE

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*) Fehler

      INTEGER       JList(*), nJList, lj, j, jout, K, nK, i, n, np,
     *              j1, m1, j2, m2, j3, m3
      LOGICAL       qOv
      REAL*8        pi, th, fi, si, ri, l, so, ro, rf, s, rzOlfGG, k0
      COMPLEX*16    im, e, ejj

      pi = 4 * dAtan (1.d0)
      im = (0.d0,1.d0)

      DO lj = 1, nJList
         j = JList(lj)

         CALL OlfHeadDup (j, qOv, jout, nK, K, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, '_', 'FPI Optik', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         DO K = 1, nK

            CALL OlfGetX (j, K, n, X, Fehler) ! X = k in m^-1
            CALL OlfCopZ (j, jout, K, K, Fehler)

            k0 = rzOlfGG (j, K, 'k0', 'cm-1', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            fi = rzOlfGG (j, K, 'f-in', 'cm', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (fi.le.0) THEN
               Fehler = 'focus fi <= 0'
               RETURN
               ENDIF
            ri = rzOlfGG (j, K, 'R-in', 'cm', Fehler) / fi
            IF (Fehler.ne.'&ff') RETURN
            ro = rzOlfGG (j, K, 'R-out', 'cm', Fehler) / fi
            IF (Fehler.ne.'&ff') RETURN
            IF (ri.le.0 .or. ro.le.0) THEN
               Fehler = 'radii must be > 0'
               RETURN
               ENDIF
            l  = rzOlfGG (j, K, 'd-mirr', 'cm', Fehler)
            IF (Fehler.ne.'&ff') RETURN
            rf = rzOlfGG (j, K, 'refl', ' ', Fehler)
            IF (Fehler.ne.'&ff') RETURN

            np = idnint( rzOlfGG (j, K, '#-pass', ' ', Fehler) )
            IF (Fehler.ne.'&ff') RETURN
            IF (np.lt.1) THEN
               Fehler = 'np < 1'
               RETURN
               ENDIF

            m1 = idnint( rzOlfGG (j, K, '#-th', ' ', Fehler) )
            IF (Fehler.ne.'&ff') RETURN
            m2 = idnint( rzOlfGG (j, K, '#-ri', ' ', Fehler) )
            IF (Fehler.ne.'&ff') RETURN
            m3 = idnint( rzOlfGG (j, K, '#-ro', ' ', Fehler) )
            IF (Fehler.ne.'&ff') RETURN

            IF (m1.lt.1 .or. m2.lt.1 .or. m3.lt.1) THEN
               Fehler = '# steps < 1'
               RETURN
               ENDIF

            DO i = 1, n
               Y(i) = 0

               ! numeric integration
               DO j3 = 1, m3
                  so = (j3-0.5)*ro/m3
                  e = (0.d0,0.d0)

                  DO j2 = 1, m2
                     si = (j2-0.5)*ri/m2
                     DO j1 = 1, m1
                        th = (j1-0.5)*pi/m1

                        s = sqrt ( 1 + si**2 + so**2 - 
     *                             2*si*so*dcos(th) )
                        ejj = rf * exp ( im * l * ( X(i) * s + 
     *                                   k0 * (s-1) ) )
                        e = e + si * dquot0 ( ejj, ( 1 - ejj**2 ))

                        ENDDO
                     ENDDO

                  Y(i) = Y(i) + so * ( 2 / ro * m3 ) *
     * ( abs ( e * conjg(e) ) / ( ri * m1 * m2 )**2 ) ** np
                  ENDDO ! j3

               ENDDO ! i

            CALL OlfPutXY0 (jout, K, n, X, Y, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO ! K

         CALL OlfClos (jout, nK, Fehler)
         ENDDO ! lj

      END ! OptikFPI
