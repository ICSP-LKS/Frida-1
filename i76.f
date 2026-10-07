C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i76   :     mode coupling integration
C
C  ====================================================================

C     Contents :
C         MCT_FuT, MCT_Fit

C     History :
C         implemented with help by A.Singh 11/97
C         dPhi_s(t)_dt and transformation GHz to ps by S.Wiebel 05/02
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE MCT_FuT (ifu, Name, Formula, ParDef, ParUni, nP)
C  --------------------------------------------------------------------

      CHARACTER*(*)   Name, Formula, ParDef, ParUni

      Name = 'MCT/'
      ParDef = '#dezi;blocksize/2;tol;ht;W1;W2;F1;F2;'//
     *         'v1;v2;vs;debug;amp;gaussW;time'
C     old line before change
C     *         'lambda;epsilon;V12;debug;amp;gaussW'
      ParUni = '0,0;0,0;0,0;1,0;-1,0;-1,0;-1,0;-1,0;0,0;0,0;0,0;0,0;//
     *          0,1;1,0;0,0'
      nP = 15

      IF     (ifu.eq.1) THEN
         CALL Append (Name,    'F12/ PhiM')
         Formula = 'Phi_m(t)'
c         CALL Append (ParDef,  'relAmp;w_sep;Gwid;rGwid')
c         CALL Append (ParUni,  '0,0;?;?;0,0')
c         nP = nP + 4
      ELSEIF (ifu.eq.2) THEN
         CALL Append (Name,    'F12/ PhiS')
         Formula = 'Phi_s(t)'
      ELSEIF (ifu.eq.3) THEN
         CALL Append (Name,    'F12/ ChiM''')
         Formula = 'X_m''(w)'
      ELSEIF (ifu.eq.4) THEN
         CALL Append (Name,    'F12/ ChiS''')
         Formula = 'X_s''(w)'
      ELSEIF (ifu.eq.5) THEN
         CALL Append (Name,    'F12/ ChiM''''')
         Formula = 'X_m''''(w)'
      ELSEIF (ifu.eq.6) THEN
         CALL Append (Name,    'F12/ ChiS''''')
         Formula = 'X_s''''(w)'
      ELSEIF (ifu.eq.7) THEN
         CALL Append (Name,    'F12/ S_M''''')
         Formula = 'S_m''''(w)'
      ELSEIF (ifu.eq.8) THEN
         CALL Append (Name,    'F12/ S_S''''')
         Formula = 'S_s''''(w)'
      ELSEIF (ifu.eq.9) THEN
         CALL Append (Name,    'F12/ dPhiS_dt')
         Formula = 'dPhi_s(t)_dt'
      ELSE
         Name    = '&undefined'
         ENDIF

      END ! MCT_FuT

C  --------------------------------------------------------------------
      SUBROUTINE MCT_Fit (ifuext, P, nf, XF, YF, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       MCO, MBL, MSV, ifu, ideb, ifuext, nf, if, is,
     *              it, i2, iL,nBL, nCh, nC4, id, nd, ii, jcout,
     *              jfout, n1, nn, ip
      PARAMETER    (MCO=2,MBL=MC,MSV=8000)
      LOGICAL       Time
      REAL*8        P(*), XF(*), YF(*), AX(MBL), Psi(MBL,MCO),
     *              Phi(MBL,MCO), Ern(MBL,MCO), dPhi(MBL,MCO),
     *              dErn(MBL,MCO), SavX(MSV), SavY(MSV), Fric(MCO),
     *              Omeg(2), falt1, falt2, RedF(MCO), RedO(MCO),
     *              AuxA(MCO), AuxB(MCO), AuxC(MCO), Xs(MSV), Ys(MSV),
     *              lambda, epsilon, V1sq, V1li, V12, ph1A, ph1B,
     *              ph2A, ph2B, ht, sub1, sub2, sub3, sub4, genau,
     *              Amp, GaussW, Phinfty, Phigma, sumY, alpha, theta,
     *              delta, twopi
      CHARACTER     Fehler*(*), aux*20
      COMMON / MathConst / twopi

      IF (Fehler.ne.'&ff') RETURN

C  Parameter :
      nd      = idnint (P(1)) ! # Dezimierungen
      nC4     = idnint (P(2)) ! 1/2 blocksize
      genau   = P(3)
      ht      = P(4) ! Startwert: wird beim Dezimieren jeweils verdoppelt
      Omeg(1) = P(5)
      Omeg(2) = P(6)
      Fric(1) = P(7)
      Fric(2) = P(8)
      ! F12-Modell-spezifisch
      V1sq    = dmax1 (0.d0, P(10))
      V1li    = dmax1 (0.d0, P(9))
! two old lines !!
C      lambda  = dinside (P(9), .5d0, 1.d0)
C      epsilon = P(10)
      V12     = dmax1 (0.d0, P(11))
      ideb    = idnint (P(12)) ! Debug ?
      Amp     = P(13)
      GaussW  = P(14)
      Time    = P(15).GT.0.0 ! wandelt GHz in ps, falls true.
                             ! FK 29-Jul-2009 added GT

C  Funktionsauswahl :
      IF (ideb.ge.1) THEN
         ifu = ideb
      ELSE
         ifu = ifuext
         ENDIF
      IF     (ifu.eq.1) THEN
         jcout = 1
         jfout = 1
      ELSEIF (ifu.eq.2) THEN
         jcout = 2
         jfout = 1
      ELSEIF (ifu.eq.3) THEN
         jcout = 1
         jfout = 2
      ELSEIF (ifu.eq.4) THEN
         jcout = 2
         jfout = 2
      ELSEIF (ifu.eq.5) THEN
         jcout = 1
         jfout = 3
      ELSEIF (ifu.eq.6) THEN
         jcout = 2
         jfout = 3
      ELSEIF (ifu.eq.7) THEN
         jcout = 1
         jfout = 4
      ELSEIF (ifu.eq.8) THEN
         jcout = 2
         jfout = 4
      ELSEIF (ifu.eq.9) THEN
         jcout = 2
         jfout = 5
      ELSE
         Fehler = 'MCT/ few-corr-models/ funtion not implemented'
         RETURN
         ENDIF

C  "Ubernahme der MCT-Parameter :
      ! Numerik-Setup :
      nBL     = 2 * nC4 + 1
      nCh     = 4 * nC4 + 1
      IF (nCh.gt.MBL) THEN
         Fehler = 'MCT/ too much work space required'
         RETURN
         ENDIF
      IF ((nd+2)*2*nC4+1.gt.MSV) THEN
         Fehler = 'MCT/ too much internal storage required'
         RETURN
         ENDIF
      ! Modell :
      IF (Omeg(1).le.0 .or. Omeg(2).le.0) THEN
         Fehler = 'Invalid Omega .le. 0'
         RETURN
         ENDIF

      ! Sjoegren-Modell-spezifisch
! two old lines have been uncommented before change !!
C      V1sq = dmax1 (0.d0, 1 / lambda**2)
C      V1li = dmax1 (0.d0, 2*dsqrt(V1sq) - V1sq + epsilon)
c       Print *, 'Coupling coefficients V1li V1sq : ', V1li, V1sq

C  Integration im Zeitraum
      ! Anfangswerte fuer t=0 :
      Phi(1,1) = 1
      Phi(1,2) = 1
      Ern(1,1) = V1li * Phi(1,1) + V1sq * Phi(1,1)**2        ! Sjoegren-Modell
      Ern(1,2) = V12  * Phi(1,1) * Phi(1,2)

      ! "startwerte"
      sub1= 1 - Fric(1)*ht
      sub2= Omeg(1)**2 * ht
      sub3= 1 - Fric(2)*ht
      sub4= Omeg(2)**2 * ht
      falt1=0
      falt2=0
      DO is = 1, nCh
         AX(is) = (is-1)*ht
         ENDDO
      DO is = 2, nBL
         Psi(is,1) = sub1*Psi(is-1,1) - sub2*Phi(is-1,1) - sub2*falt1
         Psi(is,2) = sub3*Psi(is-1,2) - sub4*Phi(is-1,2) - sub4*falt2
         Phi(is,1) = ht*Psi(is,1)+Phi(is-1,1)
         Phi(is,2) = ht*Psi(is,2)+Phi(is-1,2)
         Ern(is,1) = V1li * Phi(is,1) + V1sq * Phi(is,1)**2
         Ern(is,2) = V12  * Phi(is,1) * Phi(is,2)
         falt1 = 0
         falt2 = 0
         DO it = 1, is-1
            falt1 = falt1 + ht/2 *
     * (Ern(1+is-it,1)*Psi(it,1) + Ern(is-it,1)*Psi(it+1,1))
            falt2 = falt2 + ht/2*
     * (Ern(1+is-it,2)*Psi(it,2) + Ern(is-it,2)*Psi(it+1,2))
            ENDDO
         ENDDO
      ! "initialisiere"
      DO is = 2, nBL
         dPhi(is,1) = (Phi(is-1,1) + Phi(is,1)) / 2 ! schon durch ht geteilt
         dErn(is,1) = (Ern(is-1,1) + Ern(is,1)) / 2
         dPhi(is,2) = (Phi(is-1,2) + Phi(is,2)) / 2
         dErn(is,2) = (Ern(is-1,2) + Ern(is,2)) / 2
         ENDDO
cdeb      Print *, ' DEB ', Psi(3,1), Phi(3,1), Ern(3,1)

      ! Speichern :
      DO is = 1, nBL ! sozusagen id=0
         SavX(is) = AX(is)
         SavY(is) = Phi(is,jcout)
         ENDDO

      ! Hauptschleife:
      id = 1
 100  CONTINUE
         ! "hauptalgorithmus"
         RedF(1) = 0.5*Fric(1) / ht / Omeg(1)**2
         RedF(2) = 0.5*Fric(2) / ht / Omeg(2)**2
         RedO(1) = 1.d0 / ht**2 / Omeg(1)**2
         RedO(2) = 1.d0 / ht**2 / Omeg(2)**2
         AuxA(1) = 1 + dErn(2,1) + 3*RedF(1) + 2*RedO(1)
         AuxA(2) = 1 + dErn(2,2) + 3*RedF(2) + 2*RedO(2)
         IF (AuxA(1).eq.0 .or. AuxA(2).eq.0) THEN
            Fehler = 'unexpectedly AuxA=0 --> division by 0 imminent'
            RETURN
            ENDIF
         AuxB(1) = (1 - dPhi(2,1)) / AuxA(1)
         AuxB(2) = (1 - dPhi(2,2)) / AuxA(2) ! war zeitweise AuxA(1) !??
         DO is = nBL+1, nCh
            Phi(is,1) = 0
            Ern(is,1) = 0
            Phi(is,2) = 0
            Ern(is,2) = 0
            ! Faltungsintegral :
            i2 = (is+1)/2
            AuxC(1) = Ern(i2,1) * Phi(1+is-i2,1)
            AuxC(2) = Ern(i2,2) * Phi(1+is-i2,2)
            DO it = 2, i2
               AuxC(1) = AuxC(1) +
     *              dErn(it,1) * (Phi(2+is-it,1) - Phi(1+is-it,1)) +
     *              dPhi(it,1) * (Ern(2+is-it,1) - Ern(1+is-it,1))
               AuxC(2) = AuxC(2) +
     *              dErn(it,2) * (Phi(2+is-it,2) - Phi(1+is-it,2)) +
     *              dPhi(it,2) * (Ern(2+is-it,2) - Ern(1+is-it,2))
               ENDDO ! it
            DO it = i2+1, is-i2+1
               AuxC(1) = AuxC(1) +
     *              dPhi(it,1) * (Ern(2+is-it,1) - Ern(1+is-it,1))
               AuxC(2) = AuxC(2) +
     *              dPhi(it,2) * (Ern(2+is-it,2) - Ern(1+is-it,2))
               ENDDO ! it
            AuxC(1) = AuxC(1) +
     * RedO(1) * (-5*Phi(is-1,1)+4*Phi(is-2,1)-Phi(is-3,1)) +
     * RedF(1) * (-4*Phi(is-1,1)+  Phi(is-2,1))
            AuxC(2) = AuxC(2) +
     * RedO(2) * (-5*Phi(is-1,2)+4*Phi(is-2,2)-Phi(is-3,2)) +
     * RedF(2) * (-4*Phi(is-1,2)+  Phi(is-2,2))
            AuxC(1) = AuxC(1) / AuxA(1)
            AuxC(2) = AuxC(2) / AuxA(2)

            ! Jetzt die Iteration (im Sjoegren-Modell hintereinander moeglich)
            ! ---> innerste Schleife (bottleneck) <--- !
            ! Zuerst Korrelator 1:
            ph1B = Phi(is-1,1)
 102        ph1A = ph1B
            ph1B = -AuxC(1) + AuxB(1) *
     *           ( V1li * ph1A + V1sq * ph1A**2 )
            IF (dabs(ph1B-ph1A).gt.genau) GOTO 102

            ! Dann Korrelator 2:
            ph2B = Phi(is-1,2)
 104        ph2A = ph2B
            ph2B = -AuxC(2) + AuxB(2) * ( V12 * ph1B * ph2A )
            IF (dabs(ph2B-ph2A).gt.genau) GOTO 104

            ! Ergebnis :
            Phi(is,1) = ph1B
            Phi(is,2) = ph2B
            Ern(is,1) = V1li * ph1B + V1sq * ph1B**2        ! Sjoegren-Modell
            Ern(is,2) = V12 * ph1B * ph2B
            ENDDO ! is

         ! "speichere"
         DO is = nBL+1, nCh ! fuer id=1..nd
            SavX((id-1)*(nBL-1)+is) = AX(is)
            SavY((id-1)*(nBL-1)+is) = Phi(is,jcout)
            ENDDO

         ! "dezimierungen"
         IF (id.gt.nd) GOTO 199
         ht = 2*ht
         DO is = 2, nCh
            AX(is) = (is-1) * ht
            ENDDO
         DO is = 2, nC4+1
            dPhi(is,1) = (dPhi(2*is-2,1) + dPhi(2*is-1,1)) / 2
            dErn(is,1) = (dErn(2*is-2,1) + dErn(2*is-1,1)) / 2
            dPhi(is,2) = (dPhi(2*is-2,2) + dPhi(2*is-1,2)) / 2
            dErn(is,2) = (dErn(2*is-2,2) + dErn(2*is-1,2)) / 2
            ENDDO
         DO is = nC4+2, nBL
            dPhi(is,1) =
     * (Phi(2*is-3,1) + 4*Phi(2*is-2,1) + Phi(2*is-1,1)) / 6
            dErn(is,1) =
     * (Ern(2*is-3,1) + 4*Ern(2*is-2,1) + Ern(2*is-1,1)) / 6
            dPhi(is,2) =
     * (Phi(2*is-3,2) + 4*Phi(2*is-2,2) + Phi(2*is-1,2)) / 6
            dErn(is,2) =
     * (Ern(2*is-3,2) + 4*Ern(2*is-2,2) + Ern(2*is-1,2)) / 6
            ENDDO
         DO is = 2, nBL
            Phi(is,1) = Phi(2*is-1,1)
            Ern(is,1) = Ern(2*is-1,1)
            Phi(is,2) = Phi(2*is-1,2)
            Ern(is,2) = Ern(2*is-1,2)
            ENDDO
         id = id + 1
         GOTO 100
 199  CONTINUE

C  Output selection :

      IF (Time) THEN
         ! cambia la scala da 1/omega(GHz^(-1)) a t(picosecondi) ;-)
         DO ip = 1, nd*(nBL-1)+nCh
            SavX(ip) = SavX(ip)*1000/twopi
            ENDDO
            ENDIF


      IF     (jfout.eq.1) THEN
         ! interpolation to experimental grid:
         DO if = 1, nf
            CALL LinIntPolArray (SavX, SavY, nd*(nBL-1)+nCh,
     *                           XF(if), YF(if), 0.d0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

      ELSEIF (jfout.eq.5) THEN
         ! Amp * Ableitung nach der Zeit:
        DO il = 1, nd*(nBL-1)+nCh-1
           Xs(il) = (SavX(il+1) + SavX(il)) / 2
           Ys(il) = -Amp*(SavY(il+1) - SavY(il))/(SavX(il+1)-SavX(il))
C          Ds(il) = dsqrt (SavY(il+1)**2 +
C     *              SavY(il)**2) / (SavX(il+1) - SavX(il))
             ENDDO

         ! interpolation to experimental grid:
         DO if = 1, nf-1
            CALL LinIntPolArray (Xs, Ys, nd*(nBL-1)+nCh-1,
     *                           XF(if), YF(if), 0.d0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

C         DO if = 1, nf-1
C            XF(if) = Xs(if)
C            YF(if) = Ys(if)
C            ENDDO



      ELSEIF (jfout.ge.2 .and. jfout.le.4) THEN
         ! get Phi(infty)
         Phinfty = 0
         DO is = nBL+1, nCh
            Phinfty = Phinfty + SavY(nd*(nBL-1)+is)
            ENDDO
         Phinfty = Phinfty / (nCh-nBL)
         Phigma  = 0
         DO is = nBL+1, nCh
            Phigma = Phigma + (SavY(nd*(nBL-1)+is)-Phinfty)**2
            ENDDO
         Phigma = Phigma / (nCh-nBL)
         IF (Phigma.gt.1d-12) THEN
            Print *, ' Phi(infty) Phi(sigma) ', Phinfty, Phigma
            CALL NiceNum (Phigma, aux, i2)
            Fehler = 'Phi(t) converges only within +- '//aux(1:i2)
            RETURN
            ENDIF
         DO is = 1, (nd+1)*(nBL-1)+1
            SavY(is) = SavY(is) - Phinfty
            ENDDO

         ! blockwise Fourier transform
         DO if = 1, nf
            YF(if) = 0
            ENDDO
         DO id = 0, nd
            n1 = id*(nBL-1)+1
            nn = id*(nBL-1)+nBL
            ht = (SavX(nn) - SavX(n1)) / (nBL-1)
            DO if = 1, nf
               ! lowest-order Filon integration
               theta = XF(if) * ht
               IF (dabs(theta).lt.1d-40) THEN
                  alpha = 1
                  delta = 0
               ELSE
                  alpha = ( dsin(theta/2) / (theta/2) )**2
                  delta = (theta-dsin(theta)) / theta**2
                  ENDIF
               sumY = 0
               IF (jfout.eq.3 .or. jfout.eq.4) THEN ! real part of S(w) = X''/w
                  DO is = n1+1, nn-1
                     sumY = sumY + SavY(is) * dcos(SavX(is)*XF(if))
                     ENDDO ! is
                  sumY = alpha * ( sumY +
     *               ( SavY(n1) * dcos(SavX(n1)*XF(if))
     *               + SavY(nn) * dcos(SavX(nn)*XF(if)) ) / 2 )
     *                    + delta *
     *               ( - SavY(n1) * dsin(SavX(n1)*XF(if))
     *                 + SavY(nn) * dsin(SavX(nn)*XF(if)) )
               ELSEIF (jfout.eq.2) THEN  ! imaginary part of S(w) = X'/w
                  DO is = n1+1, nn-1
                     sumY = sumY + SavY(is) * dsin(SavX(is)*XF(if))
                     ENDDO ! is
                  sumY = alpha * ( sumY +
     *               ( SavY(n1) * dsin(SavX(n1)*XF(if))
     *               + SavY(nn) * dsin(SavX(nn)*XF(if)) ) / 2 )
     *                    + delta *
     *               ( + SavY(n1) * dcos(SavX(n1)*XF(if))
     *                 - SavY(nn) * dcos(SavX(nn)*XF(if)) )
               ELSE
                  Fehler = 'PROGRAM ERROR / Re vs Im ???'
                  RETURN
                  ENDIF ! Re vs. Im
               YF(if) = YF(if) + ht * sumY
               ENDDO ! if
            ENDDO ! id

         IF (jfout.eq.2 .or. jfout.eq.3) THEN ! susceptibility
            DO if = 1, nf
               YF(if) = Amp * YF(if) * XF(if)
                  ! without factor 1/pi: is still in Goetze convention
               ENDDO
         ELSEIF (jfout.eq.4) THEN ! scattering law
            DO if = 1, nf
               YF(if) = Amp * YF(if) * 2 / twopi ! Faktor 1/pi 22okt98
               ENDDO
            IF (GaussW.gt.0) THEN ! add pseudoelastic curve {18jan00}
               DO if = 1, nf
                  YF(if) = YF(if) + ( Amp / dsqrt(twopi/2) / GaussW ) *
     * Phinfty * dexp (- (XF(if)/GaussW)**2 )
                  ENDDO
               ENDIF
         ELSE
            Fehler = 'Scattering law or susceptibility ?'
            RETURN
            ENDIF

      ELSE
         Fehler = 'PROGRAM ERROR/ unexpected output choice'
         ENDIF

      END ! MCT_Fit

C  --------------------------------------------------------------------
      SUBROUTINE MCT_SFT (ifu, Name, Formula, ParDef, ParUni, nP)
C  --------------------------------------------------------------------

      CHARACTER*(*)   Name, Formula, ParDef, ParUni

      Name = 'MCT/'
      ParDef =
     * '#dezi;blocksize/2;tol;ht;t0;lambda;sigma;delta;h;f;alg;deb'
      ParUni = '0,0;0,0;0,0;1,0;1,0;0,0;0,0;1,0;0,1;0,1;0,0;0,0'
      nP = 12

      IF     (ifu.eq.11) THEN
         CALL Append (Name,    'Phi(t,sig,del)')
         Formula = 'Phi(t,sig,del)'
      ELSEIF (ifu.eq.12) THEN
         CALL Append (Name,    'Chi''(w,sig,del)')
         Formula = 'Chi''(w,sig,del)'
      ELSEIF (ifu.eq.13) THEN
         CALL Append (Name,    'Chi''''(w,sig,del)')
         Formula = 'Chi''''(w,sig,del)'
      ELSEIF (ifu.eq.14) THEN
         CALL Append (Name,    'S(w,sig,del)')
         Formula = 'S(w,sig,del)'
      ELSE
         Name    = '&undefined'
         ENDIF

      END ! MCT_SFT

C  --------------------------------------------------------------------
      SUBROUTINE MCT_SkaFu (ifu, P, nf, XF, YF, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      INTEGER       MBL, MSV, ifu, nf, if, is, it, i2, ntot, nBL, nCh,
     *              nC4, id, nd, ii, jcout, jfout, n1, nn, ialg, ideb
      PARAMETER    (MBL=MC,MSV=18000)
      REAL*8        P(*), XF(*), YF(*), AX(MBL),
     *              Phi(MBL), dPhi(MBL), SavX(MSV), SavY(MSV),
     *              lambda, aex, bex, sigma, hopp, ampl, offs,
     *              ht0, t0, genau, bexhopp, phinfty, powft,
     *              Sg, Cg, Rg, Qg, phA, phB, phC, fuB, fuC, ht, xll,
     *              sumY, alpha, theta, delta, twopi
      LOGICAL       qSavSubA, qSavDiff, qSavTime, qAutRang, qTraSine
      CHARACTER     Fehler*(*), aux*20
      COMMON / MathConst / twopi

      IF (Fehler.ne.'&ff') RETURN

C  "Ubernahme der Parameter :
      nd      = idnint (P(1)) ! # Dezimierungen
      nC4     = idnint (P(2)) ! 1/2 blocksize
      genau   = P(3)
      ht0     = dinside (P( 4), 1.d-24, 1.d12) ! beide Grenzen willkuerlich
      t0      = dinside (P( 5), 1.d-12, 1.d24) ! beide Grenzen willkuerlich
      lambda  = dinside (P( 6), .5d0, 1.d0)
      sigma   = dinside (P( 7),-1.d2, 1.d2) ! beide Grenzen willkuerlich
      hopp    = dinside (P( 8), 0.d0, 1.d3) ! Obergrenze willkuerlich
      ampl    = dinside (P( 9), 0.d0, 1.d3) ! Obergrenze willkuerlich
      offs    = dinside (P(10), 0.d0, 1.d3) ! Obergrenze willkuerlich
      ialg    = idnint (P(11)) ! welchen Algorithmus ?
      ideb    = idnint (P(12)) ! Debug ?

C  First checks :
      nBL     = 2 * nC4 + 1
      nCh     = 4 * nC4 + 1 ! = 2*nBL - 1
      IF (nCh.gt.MBL) THEN
         Fehler = 'MCT/ too much work space required'
         RETURN
         ENDIF
      IF ((nd+2)*2*nC4+1.gt.MSV) THEN
         Print *, (nd+2)*2*nC4+1, MSV
         Fehler = 'MCT/ too much internal storage required'
         RETURN
         ENDIF

      IF (ialg.lt.1 .or. ialg.gt.3) THEN
         Fehler =
     * 'Algorithm must be reg.falsi(1) simple iter(2) sqrt(3)'
         RETURN
         ENDIF

C  Funktionsauswahl :
      IF (ideb.gt.0) THEN
         jfout = ideb
      ELSE
         jfout = ifu-10
         ENDIF

      IF (jfout.lt.1 .or. jfout.gt.6) THEN
         Print *, 'function no. ', jfout, ' < ', ifu, ideb
         Fehler = 'MCT/ SkaFu/ funtion not implemented'
         RETURN
         ENDIF

C  Aktuelle Funktionenliste :
      !  1 = Phi(t)
      !  2 = Chi'(w)
      !  3 = Chi''(w)
      !  4 = S(w)
      !  5 = Phi(t)-t^-a
      !  6 = d(Phi(t)-t^-a)/dt

      qSavSubA = (jfout.gt.1)
      qSavDiff = (jfout.gt.1 .and. jfout.ne.5)
      qSavTime = (jfout.eq.1 .or. jfout.ge.5)
      qAutRang = (jfout.ge.2 .and. jfout.le.4)
      qTraSine = (jfout.eq.3 .or. jfout.eq.4)

C  Langzeitgrenze setzen :
      xll     = ht0*((nBL-1)*2.d0**(nd+1)-2)
      IF (qSavTime .and. XF(nf).gt.xll) THEN
         Print *, ' longest time ', xll
         Fehler = ' iteration does not extend far enough'
         RETURN
         ENDIF

      CALL GoetzeExp (lambda, aex, bex)
c      Print *, ' l -> a b ', lambda, aex, bex

      ! sehr-Langzeit-Limes :
      IF (hopp.gt.0) THEN ! universell fuer Glas / Uebergang / Fluessigkeit
         phinfty = 0
         bexhopp = dmax1 (bex, 0.5d0) ! fuer lambda>pi/4 1/2 statt b
      ELSEIF (sigma.lt.0) THEN ! ideale (hoppingfreie) Fluessigkeit
         phinfty = 0
         bexhopp = bex
      ELSE ! ideales (hoppingfreies) Glas
         phinfty = dsqrt(sigma/(1-lambda))
         bexhopp = 0
         ENDIF

C  Beginn der Integration / Setze Anfangswerte :
      ht = ht0
      DO is = 1, nCh
         AX(is) = (is-1)*ht
         ENDDO

      Phi(1)  = 0 ! undefiniert, wird aber auch gebraucht
      dPhi(1) = 0 ! erst recht undefiniert, wird aber auch nicht gebraucht
      DO is = 2, nBL
         Phi(is)  = dpow0 (t0/AX(is), aex)
         dPhi(is) = (t0/ht/(1-aex)) *
     *                (dpow0(AX(is)/t0,1-aex)-dpow0(AX(is-1)/t0,1-aex))
         ENDDO

      ! Speichern :
      DO is = 1, nBL-1 ! sozusagen id=0
         SavX(is) = AX(is)
         IF (qSavSubA) THEN
            SavY(is) = 0 ! t^-a subtracted
         ELSE
            SavY(is) = Phi(is)
            ENDIF
         ENDDO

C  Hauptschleife: Integration auf logarithmischer Zeitachse
      id = 1
 100  CONTINUE
         ! "hauptalgorithmus"
         DO is = nBL+1, nCh
            ! Faltungsintegral :
            i2 = is/2 ! geaendert < war (is+1)/2
            Sg = 0
            DO it = 2, i2
               Sg = Sg + (Phi(is-it+1) - Phi(is-it)) * dPhi(it+1)
               ENDDO ! it
            DO it = 2, is-1-i2
               Sg = Sg + (Phi(is-it+1) - Phi(is-it)) * dPhi(it+1)
               ENDDO ! it
            Cg = Sg + Phi(is-i2)*Phi(i2+1) - 2*Phi(is-1)*dPhi(2)
     *              - sigma + hopp * AX(is)
            Rg = Cg / (2*dPhi(2))
            Qg = lambda / (2*dPhi(2))

            ! Jetzt die Iteration ---> innerste Schleife (bottleneck) <---
            IF     (ialg.eq.1) THEN ! regula falsi
               phB = Phi(is-1)
               phC = Phi(is-1) - ht/AX(is)
 104           CONTINUE
               fuB = -Rg + Qg * phB**2 - phB
               fuC = -Rg + Qg * phC**2 - phC
               phA = phB - (phB-phC)*fuB/(fuB-fuC)
               IF (dabs(phA-phB).gt.genau) THEN
                  phC = phB
                  phB = phA
                  GOTO 104
                  ENDIF

               Phi(is) = phA ! Ergebnis

            ELSEIF (ialg.eq.2) THEN ! simple iteration

               phB = Phi(is-1)
 102           phA = phB
               phB = -Rg + Qg * phA**2
               IF (dabs(phB-phA).gt.genau) GOTO 102

               Phi(is) = phB

            ELSE ! exact solution
               Phi(is) = dPhi(2)/lambda
     *                   - dsqrt0 ((dPhi(2)/lambda)**2 + Cg/lambda)

               ENDIF ! ialg

c            Print '(a,f13.5,1x,2i3,g13.5)', '>',
c     *              Phi(is), is, id, (Phi(is)-Phi(is-1))/Phi(is)

            ENDDO ! is

         ! "speichere" und berechne Durchschnittsamplitude
         DO is = nBL, nCh-1 ! fuer id=1..nd
            SavX((id-1)*(nBL-1)+is) = AX(is)
            IF (qSavSubA) THEN
               IF (qSavDiff) THEN
                  SavY((id-1)*(nBL-1)+is) = (Phi(is+1) - Phi(is-1)) /
     * 2 / ht + aex * dpow0 (t0/AX(is), aex) / AX(is)
               ELSE
                  SavY((id-1)*(nBL-1)+is) = Phi(is) -
     *                                      dpow0 (t0/AX(is), aex)
                  ENDIF
            ELSE
               SavY((id-1)*(nBL-1)+is) = Phi(is)
               ENDIF
            ENDDO

         IF (id.gt.nd) GOTO 199

         ! "dezimierungen"
         ht = 2*ht
         DO is = 2, nCh
            AX(is) = (is-1) * ht
            ENDDO
         DO is = 2, nC4+1
            dPhi(is) = (dPhi(2*is-2) + dPhi(2*is-1)) / 2
            ENDDO
         DO is = nC4+2, nBL
            dPhi(is) = (Phi(2*is-3) + 4*Phi(2*is-2) + Phi(2*is-1)) / 6
            ENDDO
         DO is = 2, nBL
            Phi(is) = Phi(2*is-1)
            ENDDO
         id = id + 1
         GOTO 100
 199  CONTINUE

C  Ausgabe je nach gewuenschter Darstellung :

      nd = id - 1
      ntot = nd*(nBL-1)+nCh-2

      IF     (qSavTime) THEN
         ! interpolation to experimental grid:
         DO if = 1, nf
            CALL LinIntPolArray (SavX, SavY, ntot,
     *                           XF(if), YF(if), 0.d0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            YF(if) = offs + ampl * YF(if)
            ENDDO

      ELSE

         IF (.not. (qSavSubA .and. qSavDiff)) THEN
            Fehler = 'not prepared for FT'
            RETURN
            ENDIF

         ! blockwise Fourier transform
         DO if = 1, nf
            YF(if) = 0
            ENDDO
         DO id = 0, nd
            n1 = id*(nBL-1)+1
            nn = id*(nBL-1)+nBL
            ht = (SavX(nn) - SavX(n1)) / (nBL-1)
            DO if = 1, nf
               ! lowest-order Filon integration
               theta = XF(if) * ht
               IF (dabs(theta).lt.1d-40) THEN
                  alpha = 1
                  delta = 0
               ELSE
                  alpha = ( dsin(theta/2) / (theta/2) )**2
                  delta = (theta-dsin(theta)) / theta**2
                  ENDIF
               sumY = 0
               IF (.not. qTraSine) THEN ! real part of S(w) = X''/w
                  ! cosine transform
                  DO is = n1+1, nn-1
                     sumY = sumY + SavY(is) * dcos(SavX(is)*XF(if))
                     ENDDO ! is
                  sumY = alpha * ( sumY +
     *               ( SavY(n1) * dcos(SavX(n1)*XF(if))
     *               + SavY(nn) * dcos(SavX(nn)*XF(if)) ) / 2 )
     *                    + delta *
     *               ( - SavY(n1) * dsin(SavX(n1)*XF(if))
     *                 + SavY(nn) * dsin(SavX(nn)*XF(if)) )
               ELSE                     ! imaginary part of S(w) = X'/w
                  ! sine transform
                  DO is = n1+1, nn-1
                     sumY = sumY + SavY(is) * dsin(SavX(is)*XF(if))
                     ENDDO ! is
                  sumY = alpha * ( sumY +
     *               ( SavY(n1) * dsin(SavX(n1)*XF(if))
     *               + SavY(nn) * dsin(SavX(nn)*XF(if)) ) / 2 )
     *                    + delta *
     *               ( + SavY(n1) * dcos(SavX(n1)*XF(if))
     *                 - SavY(nn) * dcos(SavX(nn)*XF(if)) )
                  ENDIF ! Re vs. Im
               YF(if) = YF(if) + ht * sumY
c               IF (XF(if).gt.15 .and. XF(if).lt.17) THEN
c               Print '(2i4,4g14.5)', id, if, XF(if), alpha, delta,
c     *                 sumY*ht
c               ENDIF
               ENDDO ! if
            ENDDO ! id

         ! Potenzgesetz t^-a wieder draufaddieren :
         IF (qTraSine) THEN
            powft = - twopi/4/dgamma1(aex)/dcos(twopi*aex/4)
         ELSE
            powft = twopi/4/dgamma1(aex)/dsin(twopi*aex/4)
            ENDIF
         DO if = 1, nf
            YF(if) = YF(if) + powft * dpow0 (XF(if)*t0, aex)
            ENDDO

         IF (jfout.eq.2 .or. jfout.eq.3) THEN
            DO if = 1, nf
               YF(if) = - ampl * YF(if) ! susceptibility
               ENDDO
         ELSEIF (jfout.eq.4) THEN
            Fehler = ' riskiere div / 0'
            RETURN
c            DO if = 1, nf
c               YF(if) = ampl * YF(if) * 2 / twopi ! scattering law
c                  ! Faktor 1/pi 22okt98
c               ENDDO
         ELSE
            Fehler = 'Scattering law or susceptibility ?'
            RETURN
            ENDIF

         ENDIF

      END ! MCT_SkaFu
