C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i75   :     nuclear resonant scattering
C
C  ====================================================================

C     Contents :
C            NFS_FuT, NFS_Fit

         ! Smirnow-Kohn-Fits JWu apr96
         ! Nichtexponentielle Relaxation JWu okt96
C  16.02.2026 Artem Panchenko: Corrected several line breaks
C  --------------------------------------------------------------------
      SUBROUTINE NFS_FuT (ifu, Name, Formula, ParDef, ParUni, nP)
C  --------------------------------------------------------------------

      CHARACTER*(*)   Name, Formula, ParDef, ParUni

      Name = 'NFS/'
      Formula = 'I(t) < D(w)='
      ParDef =
     * 'A;n/2;range;L;dL;Lwid;Lcut;cutoff;bgr;no-bunch;dt-bunch;'
      ParUni = '0,1;0,0;1,0;0,0;0,0;0,0;0,0;0,0;0,1;0,0;1,0;'
      nP = 11

      IF (ifu.eq.2) THEN
         CALL Append (Name,    ' 2 Gauss/ elastic')
         CALL Append (Formula, '2Gauss ')
         CALL Append (ParDef,  'relAmp;w_sep;Gwid;rGwid')
         CALL Append (ParUni,  '0,0;?;?;0,0')
         nP = nP + 4
      ELSEIF (ifu.eq.7) THEN
         CALL Append (Name,  ' 2 delta lines/ exponential relaxation')
         CALL Append (Formula, '2delta')
         CALL Append (ParDef,  'relAmp;w_sep;1/tau;_rel')
         CALL Append (ParUni,  '0,0;?;-1,0;0,0')
         nP = nP + 4
      ELSEIF (ifu.eq.17) THEN
         CALL Append (Name,    ' 2 delta lines/ Kohlrausch relaxation')
         CALL Append (Formula, '2delta')
         CALL Append (ParDef,  'relAmp;w_sep;1/<tau>;beta')
         CALL Append (ParUni,  '0,0;?;-1,0;0,0')
         nP = nP + 4
      ELSE
         Name    = '&undefined'
         ENDIF

      END ! NFS_FuT

C  --------------------------------------------------------------------
      SUBROUTINE NFS_Fit (ifu, P, nf, XF, YF, Fehler)
C  --------------------------------------------------------------------
         ! Front end for call of NFS-Conv-Traf by i6 (JWu 16apr96)
         ! parameters are :
         !    ifu  = number of energy distribution function
         !    P( 1) = A = linear prefactor in I(t)
         !    P( 2) = n/2 = number of positive (or negative) internal channels
         !    P( 3) = dx = energy step (unit = G_0/2 = 2.34 neV)
         !    P( 4) = L = mean effective thickness
         !    P( 5) = step in L
         !    P( 6) = Lsig = mean deviation of L
         !    P( 7) = cutoff in L (if relative weight < cutoff, stop loop)
         !    P( 8) = cutoff in phi : omit channels with Y(i) < cutoff * Y_max
         !    P( 9) = background
         !    P(10) = no. of bunches
         !    P(11) = time between bunches
         !    P(12-)= parameters of distribution function

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      PARAMETER    (MET=21, M2=MC)
      DIMENSION     P(*), XF(*), YF(*), X(MC), X3(MC), Y(MC),
     *              Y1(MC), Y2(MC), Y3(MC), Y4(MC), Y5(MC), Y6(MC),
     *              Work(M2+15), D3(MC), ET(MET), WT(MET) !Artem Work(MC2) -> Work(MC2+15) for rfftf()
      COMPLEX       YC(MC) !Artem: add for cfftf()
      CHARACTER     Fehler*(*), aus*20

      DATA          G02 /2.3322869d0/,    ! 57Fe line width Gamma_0/2 in neV
     *              hbar /0.658218d3/,    ! in neV * nsec
     *              twopi /6.2831853d0/, tol /1.d-2/

      IF (Fehler.ne.'&ff') RETURN

C  Set equidistant scale (revised 18oct96, cf. A5,21) :
      nerg = 2 * idnint (P(2)) ! must be even, should be 2^..
      ntim = nerg/2 + 1
      IF     (nerg.lt.5) THEN
         Fehler = 'bad parameter: not enough channels'
         RETURN
      ELSEIF (nerg.gt.MC) THEN
         Fehler = 'bad parameter: too many channels'
         RETURN
         ENDIF

      resowidth = G02 / hbar
      IF (P(3).le.0) THEN
         Fehler = 'bad parameter: nonpositive stepwidth'
         RETURN
         ENDIF
      dx = hbar / nerg / P(3) ! input in nsec -> energy step in neV
      Ampl = P(1) / nerg / P(3)**2 * hbar**2 / twopi**2
      DO i = 1, nerg/2
         ! positive halfwave : erg = 0 .. (nerg/2-1)*dx
         X(i) = (i-1) * dx
         ! negative halfwave : erg = -(nerg/2)*dx .. -dx
         X(nerg/2+i) = (i-1-nerg/2) * dx
         ENDDO

C  Set sample effective thickness distribution :
      ETmean  = dmax1 (P(4), 0.d0)
      IF (P(7).ge.1 .or. P(5).le.0 .or. P(6).le.0) THEN
         ! *no* averaging
         ET(1) = ETmean
         WT(1) = 1
         nT = 1
      ELSE
         Fehler = 'Dickenverteilung -> koennte Fehler enthalten'
         RETURN ! 11nov96 : Mittelung uebr Ergebnisse funktioniert nicht ??
         ETcut   = P(7)
         ETstep  = P(5)
         ETsig   = P(6)
         nT = 0 ! schlecht programmiert, NFS_SetGauss braucht Vorgabe 0
         ETsum = 0
         CALL NFS_SetGauss(1.d0, ETmean, ETsig, ETstep, ETcut,
     *                     MET, nT, ET, WT, ETsum)
         ! normieren :
         IF (ETsum.le.0) THEN
            Fehler = 'program error: thickness not normalizable'
            RETURN
            ENDIF
         DO iT = 1, nT
            WT(iT) = WT(iT) / ETsum
            ENDDO
         ENDIF

C  Calculate the nuclear resonance function phi(w) for different models :
      cut = P(8)
      IF (ifu.le.9) THEN ! analytic calculation in w-space
         sum = 0
         nR = 0
         ! Set nuclear resonance energy distribution
         IF     (ifu.eq.2) THEN ! 2 Gaussians
            CALL NFS_SetGauss (P(12),   P(13)/2/hbar, P(14)/hbar,
     *                         dx,cut,MC,nR,X3,Y3,sum)
            CALL NFS_SetGauss (1-P(12),-P(13)/2/hbar, P(14)/hbar*P(15),
     *                         dx,cut,MC,nR,X3,Y3,sum)
            CALL rSet (D3, 1, nR, 1, 0.d0)
         ELSEIF (ifu.eq.7) THEN    ! 2 delta lines
            X3(1) = -P(13)/2/hbar
            X3(2) =  P(13)/2/hbar
            Y3(1) =  P(12)
            Y3(2) =  1-P(12)
            D3(1) =  dabs(P(14))
            D3(2) =  dabs(P(14)) * dabs(P(15))
            nR = 2
            sum = 1
         ELSE
            Fehler = 'NFS/ function undefined'
            RETURN
            ENDIF

         ! Normalize distribution :
         IF (sum.le.0) THEN
            Fehler =
     * 'bad parameters: resonance distribution nonnormalizable'
            RETURN
            ENDIF
         DO i = 1, nR
            Y3(i) = Y3(i) / sum
            ENDDO

         !Set Re and Im of phi(w) :
         DO i = 1, nerg
            Y1(i) = 0 ! Re
            Y2(i) = 0 ! Im
            DO ii = 1, nR ! over resonances
               ynen = (X(i)-X3(ii))**2 + (D3(ii)+resowidth)**2
               Y1(i) = Y1(i) + Y3(ii) * (-resowidth) *
     *                 (X(i)-X3(ii)) / ynen
               Y2(i) = Y2(i) + Y3(ii) * (D3(ii)+resowidth) *
     *                 resowidth  / ynen
               ENDDO
c           YF(i) = Y1(i) ! DEB
            ENDDO
c        RETURN ! DEB

      ELSE ! numeric FT from t- to w-space

         IF     (ifu.eq.17) THEN   ! 2 delta lines / Kohlrausch
            dt   = twopi / (nerg * dx)
            fac  = resowidth * dt * dsqrt(dble(nerg))
            sepp = P(13)/2/hbar  ! separation of positions
            dfac = (1-2*P(12))*fac   ! amplitude1 - amplitude2
            beta =  P(15)
            ! 1/tau from 1/<tau> :
            tau_inv  = dquot0 (dgamma1(dquot0(1.d0,beta)) *
     *                                 dabs(P(14)), beta)

            ! set S(qt) * e^{-Gammma/2*t} * sum[e^{+-i*w_sep/2*t}]/2 :
            DO it = 1, nerg/2
               tim = (it-1) * dt
               sg  = dexp(-resowidth*tim -(dpow0(tau_inv*tim,beta)))
               wt2 = sepp * tim
               Y2(it) =-sg *  fac * dcos(wt2) ! re
               Y1(it) = sg * dfac * dsin(wt2) ! im
               ENDDO
            IF (sg.gt.cut) THEN
               Fehler = 'S(qt)*e^{-G/2*t} does not decay below cut'
               RETURN
               ENDIF
            ! halfwave corresponding to times t<0 vanishes
            DO it = 1, nerg/2
               Y1(nerg/2+it) = 0
               Y2(nerg/2+it) = 0
               ENDDO

            ! Fourier transform :
            ifail = 0
            !Artem: Replace with cfftf function from slatec/fishfft/cfftf.f: CALL C06FCF (Y1, Y2, nerg, Work, ifail)
            DO i = 1, nerg
               YC(i) = CMPLX(Y1(i), Y2(i))
            ENDDO
            CALL CFFTI (nerg, Work)
            CALL CFFTF (nerg, YC, Work)
            DO i = 1, nerg
               Y1(i) = REAL(YC(i))
               Y2(i) = AIMAG(YC(i))
            ENDDO
            IF (ifail.eq.1 .or. ifail.eq.2) THEN
               Fehler = ' FFT t->w/ #chs has bad primefactors'
               RETURN
            ELSEIF (ifail.ne.0) THEN
               Fehler = ' FFT t->w/ ifail = '//cl4(ifail)
               RETURN
               ENDIF
            ifail = 0
            CALL RFFTB (nerg, Y2, Work)
            !Artem: Replace with rfftb function from slatec/fishfft/rfftb.f: CALL C06GCF (Y2, nerg, ifail)
c           DO i= 1,nerg
c              YF(i) = Y1(i)
c              ENDDO
c           RETURN  ! DEB

         ELSE
            Fehler = 'NFS/ function undefined'
            RETURN
            ENDIF

         ENDIF

C  Calculate time profile of transmitted beam as function of thickness (ET) :
      IF (nT.gt.1) THEN
         DO i = 1, nerg
            Y5(i) = 0
            Y6(i) = 0
            ENDDO
      ELSEIF (nT.lt.1) THEN
         CALL Absturz ('NFS_ConTra', 'nT invalid')
         ENDIF

      ! loop over ET :
      DO iT = 1, nT
         et2 = ET(iT)/2

         ! Set Re and Im of e^{-L*phi(w)} :
         DO i = 1, nerg
            Y3(i) = dcos ( et2*Y1(i)) * dexp (-et2*Y2(i)) - 1  ! Re e^{i..}
            Y4(i) = dsin ( et2*Y1(i)) * dexp (-et2*Y2(i))      ! Im e^{i..}
            ENDDO
               ! der Zusatzterm -2 entspricht einer delta-Funktion bei t=0
               ! (2 und nicht 1 ad-hoc korrigiert: Punkt z"ahlt nur halb ?)

         ! Now the y-Trafo :
         ifail = 1 ! silent exit
         !Artem: Replace with cfftf function from slatec/fishfft/cfftf.f: CALL C06FCF (Y3, Y4, nerg, Work, ifail)
         DO i = 1, nerg
            YC(i) = CMPLX(Y3(i), Y4(i))
         ENDDO
         CALL CFFTI (nerg, Work)
         CALL CFFTF (nerg, YC, Work)
         DO i = 1, nerg
            Y3(i) = REAL(YC(i))
            Y4(i) = AIMAG(YC(i))
         ENDDO

         IF (ifail.eq.1 .or. ifail.eq.2) THEN
            Fehler = ' FFT/ #chs has bad primefactors'
            RETURN
         ELSEIF (ifail.ne.0) THEN
            Fehler = ' FFT/ ifail = '//cl4(ifail)
            RETURN
            ENDIF

         IF (nT.eq.1) THEN
            ! abbreviated procedure, if there is no distribution of ET's
            ! Correct for prefactor and calculate the squared modulus;
            ! Retain only results for t>0
            DO i = 1, ntim
               Y(i) = Y3(i)**2 + Y4(i)**2
                  ! factor hbar**2 eliminated > unit of I(t) is now nsec^-2
               ENDDO
            GOTO 80
            ENDIF

         ! sum Re and Im weigted with WT :
         DO i = 1, ntim
            Y5(i) = Y5(i) + WT(iT) * Y3(i)
            Y6(i) = Y6(i) + WT(iT) * Y4(i)
            ENDDO

         ENDDO ! end loop over ET.

            print *, 'DEB nach FT : ', ntim, nerg/2,
     *               Y5(nerg/2), Y6(nerg/2)

      ! Correct for prefactor and calculate the squared modulus;
      DO i = 1, ntim
         Y(i) = Y5(i)**2 + Y6(i)**2
         ENDDO
 80   CONTINUE ! arrive here if there is no distribution of ET's (nT=1)

C  New x-scale :
      DO i = 1, ntim
         X(i) = (i-1) *  twopi / (dx * nerg)
         ENDDO

C  Loop over bunches and interpolation to experimental time grid :
      nbunches = idnint (P(10))
      IF     (nbunches.lt.1) THEN
         Fehler = 'bad parameter: # bunches < 1'
         RETURN
      ELSEIF (nbunches.gt.20) THEN
         Fehler = 'bad parameter: too many bunches required'
         RETURN
      ELSEIF (nbunches.gt.1 .and. P(11).le.0) THEN
         Fehler = 'bad parameter: time between bunches nonpositive'
         RETURN
      ELSEIF (XF(1).lt.X(1)) THEN
         Fehler = 'unphysical data set: experimental times < 0'
         RETURN
      ELSEIF (XF(nf)+(nbunches-1)*P(11).gt.X(ntim)) THEN
         CALL NiceNum((XF(nf)+(nbunches-1)*P(11))/X(ntim), aus, ia)
         Fehler = 't-step too short by a factor of '//aus
         RETURN
         ENDIF

      DO if = 1, nf
         YF(if) = P(9) ! background
         DO ibunch = 1, nbunches ! loop over bunches
            t = XF(if) + (ibunch-1)*P(11)
            CALL LinIntPolArray (X, Y, ntim, t, yi, 0.d0, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            YF(if) = YF(if) + Ampl*yi
            ENDDO
         ENDDO

      END ! NFS_Fit

C  --------------------------------------------------------------------
      SUBROUTINE NFS_SetGauss (A, Pos, Wid, step, cut, nmax,
     *                         n1, X1, Y1, sum)
C  --------------------------------------------------------------------

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      DIMENSION      X1(*), Y1(*)

      X1(n1+1) = Pos
      Y1(n1+1) = A
      sum = sum + A
      n1 = n1 + 1

      imax = (nmax-n1)/2
      DO i = 1, imax
         xg = i * step
         gauss = A * dexp (-(xg/Wid)**2)
         IF (gauss.lt.A*cut) RETURN
         X1(n1+1) = Pos - xg
         X1(n1+2) = Pos + xg
         Y1(n1+1) = gauss
         Y1(n1+2) = gauss
         n1 = n1 + 2
         sum = sum + 2 * gauss
         ENDDO

      END ! NFS_SetGauss
