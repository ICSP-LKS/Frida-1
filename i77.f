C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i77   :     mode coupling integration / full model
C
C  ====================================================================

C     Link only when needed: huge memory consumption

C     Contents :
C        FullMCT
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE FullMCT (nJList, JList, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      INTEGER       MCO, MBL, MKK
      PARAMETER    (MCO=2,MBL=MC,MKK=MCO*MK3V)
      INTEGER       ifu, nf, ns, istep, if, is, it, i2, niter, maxiter,
     *              nBL, nCh, nC4, id, nd, ii, jcout, jfout, n1, nn,
     *              nout, nJList, JList(*), nK, nKin, K, KK, KKK,
     *              Kin, Kout, nK1, nK2, lj, j, jout, nZ
      REAL*8        AX(MBL), Psi(MKK,MBL),
     *              Vtx1(MK3V,MK3V,MK3V), Vtx2(MK3V,MK3V,MK3V),
     *              Phi(MKK,MBL), Ern(MKK,MBL), dPhi(MKK,MBL),
     *              dErn(MKK,MBL), SavX(MS3V,MK3V), SavY(MS3V,MK3V),
     *              Fric(MKK), Omeg2(MKK), RedF(MKK), RedO(MKK),
     *              OmRel(MCO), FrRel(MCO),
     *              AuxA(MKK), AuxB(MKK), AuxC(MKK),
     *              AuxZ(MKK), AuxY(MKK),
     *              ht, ht0, genau, ungenau, ualt, xmi, xma, xll,
     *              delQ, eta, diam, q0, dens, Tm, qs, ks, ps,
     *              Phinfty, Phigma, sumY, alpha, theta, delta, twopi,
     *              Z(MZ), rzOlfGG, PercYevDiam
      LOGICAL       qCurve, qEqui
      CHARACTER     Fehler*(*), aux*20, docGrid*80, Co*40, Un*40
      COMMON / MathConst / twopi

      IF (Fehler.ne.'&ff') RETURN

C  Output selection :
      jcout = iAskDMu ('Save correlator (1=density, 2=self)',
     *                 jcout, 1, 2)
      jfout = iAskDMu ('Save Phi(1), Chi''(2), Chi''''(3), S(4), f(5)',
     *                 jfout, 1, 5)
      nf = 1024

C  Numeric Set-up :
      ! Numerik-Setup :
      genau = rAskD ('precision', genau)
      IF (jfout.le.4) THEN
         nd      = iAskDMu ('# decimations', nd, 1, 1000)
         nC4     = iAskDMu ('1/2 blocksize', nC4, 1, 1000)
         nBL     = 2 * nC4 + 1
         nCh     = 4 * nC4 + 1
         IF (nCh.gt.MBL) THEN
            Fehler = 'MCT/ too much work space required'
            RETURN
            ENDIF
         IF ((nd+2)*2*nC4+1.gt.MS3V) THEN
            Fehler = 'MCT/ too much internal storage required'
            RETURN
            ENDIF
         ht0 = rAskD ('initial time step (psec)', ht0)
         xll = ht0*(nBL-1)*2.d0**(nd+1)
         Print '(a,g12.6)', ' longest time will be ', xll

         IF (jfout.gt.1) THEN
            xmi =  rAskDMu ('Lowest frequency', xmi,1/xll,1/twopi/ht0)
            nout = iAskDMu ('Number of data points', nout, 1, MC)
         ELSE
            nout = MC
            ENDIF

      ELSE
         maxiter = iAskD ('max # iterations', maxiter)
         IF (maxiter.le.0) RETURN
         ENDIF ! full dynamics

C  Input files :

      DO lj = 1, nJList
      j = JList(lj)

      IF (jfout.le.4) THEN
         ENDIF ! full dynamics

      CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      qCurve = qOlfGdef (j, '?cu', 0, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF (qCurve) THEN
         Fehler = ' File is a curve'
         RETURN
         ENDIF

      ! Input x becomes output z :
      Print *, ' resetting xyz labels'
      CALL OlfCnuG (j, 'x', Co, Un, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF     (Co.ne.'q') THEN
         Fehler = ' x-coordinate not q but '//Co
         RETURN
      ELSEIF (Un.ne.'A-1') THEN
         Fehler = ' x-unit not A-1 but '//Un
         RETURN
         ENDIF

      ! Check input y :
      CALL OlfCnuG (j, 'y', Co, Un, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF     (Co.ne.'S(q)') THEN
         Fehler = ' y-coordinate not S(q) but '//Co
         RETURN
      ELSEIF (Un.ne.' ') THEN
         Fehler = ' y not dimensionless but in '//Un
         RETURN
         ENDIF

      ! Set output coordinates :
      IF (jfout.eq.5) THEN
         ! CALL OlfCnuP (jout, 'x', 'q', 'A-1', Fehler)
      ELSE
         CALL OlfCnuP (jout, 'z+', 'q', 'A-1', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (jfout.eq.1) THEN
            CALL OlfCnuP (jout, 'x', 't', 'psec', Fehler)
         ELSE
            CALL OlfCnuP (jout, 'x', 'w/2pi', 'THz', Fehler)
            ENDIF
         ENDIF

      IF     (jfout.eq.1) THEN
         CALL OlfCnuP (jout, 'y', 'Phi(q,t)', ' ', Fehler)
      ELSEIF (jfout.eq.2) THEN
         CALL OlfCnuP (jout, 'y', 'Chi''(q,w)', ' ', Fehler)
      ELSEIF (jfout.eq.3) THEN
         CALL OlfCnuP (jout, 'y', 'Chi''''(q,w)', ' ', Fehler)
      ELSEIF (jfout.eq.4) THEN
         CALL OlfCnuP (jout, 'y', 'S~(q,w)', 'THz-1', Fehler)
      ELSEIF (jfout.eq.5) THEN
         CALL OlfCnuP (jout, 'y', 'f(q)', ' ', Fehler)
      ELSE
         Fehler = 'PROGRAM ERROR/ jfput oor'
         ENDIF
      IF (Fehler.ne.'&ff') RETURN

      CALL iOlfP (jout, '?det-bal-sym',   0, Fehler)
      CALL iOlfP (jout, '?sam-erg-gain', -1, Fehler)

      ! Documentation, open output file :
      CALL OlfComAdd (jout, '_', 'mode-coupling model', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C Loop: input spectra :
      Print *, ' input spectra'
      nKin = iOlfG (j, '#spectra', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO Kin = 1, nKin

      CALL OlfGetZ (j, Kin, nZ, Z, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      eta = rzOlfGG (j, Kin, 'eta', ' ', Fehler)
      Print *, ' found eta = ', eta
      IF (Fehler.ne.'&ff') RETURN

      ! diameter can be given either explicitely or through q0 :
      diam = rzOlfGG (j, Kin, 'diam', 'A', Fehler)
      IF     (Fehler(1:4).eq.'&pnf') THEN
         Fehler = '&ff'
         q0 = rzOlfGG (j, Kin, 'q0', 'A-1', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         diam = PercYevDiam (q0, eta)
      ELSEIF (Fehler.ne.'&ff') THEN
         RETURN
      ELSE
         q0 = rzOlfGG (j, Kin, 'q0', 'A-1', Fehler)
         IF (Fehler.eq.'&ff') THEN
            Fehler = 'found both diam and q0'
            RETURN
         ELSEIF (Fehler(1:4).ne.'&pnf') THEN
            RETURN
            ENDIF
         Fehler = '&ff'
         ENDIF
      Print *, ' found diam[A] = ', diam

      Tm = rzOlfGG (j, Kin, '<v^2>', '(A/psec)^2', Fehler)
      Print *, ' found T/m[(A/psec)^2] = ', Tm
      IF (Fehler.ne.'&ff') RETURN

      IF (diam.le.0 .or. eta.le.0) THEN
         Fehler = ' invalid hard spheres parameters'
         RETURN
         ENDIF

      IF (jcout.ne.5) THEN
         OmRel(1) = rzOlfGG (j, Kin, 'Omega2_rel', ' ', Fehler)
         Print *, ' found Omega_rel = ', OmRel(1)
         IF (Fehler.ne.'&ff') RETURN
         FrRel(1) = rzOlfGG (j, Kin, 'Friction', 'psec', Fehler)
         Print *, ' found Friction = ', FrRel(1)
         IF (Fehler.ne.'&ff') RETURN

         IF (jcout.ge.2) THEN
         OmRel(2) = rzOlfGG (j, Kin, 'Omega2_rel_self', ' ', Fehler)
         Print *, ' found Omega_rel_self = ', OmRel(2)
         IF (Fehler.ne.'&ff') RETURN
         FrRel(2) = rzOlfGG (j, Kin, 'Friction_self', 'psec', Fehler)
         Print *, ' found Friction_self = ', FrRel(2)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF

         ENDIF

C Input: static structure factor S(k):
      Print *, ' read S(q)'
      CALL OlfGetXY (j, Kin, nK, X3, Y2, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      IF     (nK.gt.MK) THEN
         Fehler = 'too many spectra required'
         RETURN
      ELSEIF (nK.gt.MK3V) THEN
         Fehler = 'too many q-points / recompile with bigger MK3V'
         RETURN
      ELSEIF (jfout.ne.5 .and. nKin*nK.gt.MK3V) THEN
         Fehler =
     * 'too many spectra required / recompile with bigger MK3V'
         RETURN
         ENDIF
      CALL CheckScale (nK, X3, 1d-6, qEqui, delQ)
      IF (.not. qEqui) THEN
         Fehler = 'static structure must be given on equidistant grid'
         RETURN
         ENDIF
      IF (delQ.le.0 .or. X3(1).lt.0) THEN
         Fehler = 'S(q) must be given on positive grid'
         RETURN
         ENDIF
      IF (.not. qEqTol (X3(1), delQ/2, 1.d-6)) THEN
         Fehler = 'S(q): grid must start at 0.5 stepwidth'
         RETURN
         ENDIF

      nK1 = nK + 1
      IF (jcout.gt.1) THEN
         nK2 = 2*nK
      ELSE
         nK2 = nK
         ENDIF

C  Vertizes vorbereiten :
      Print *, ' set frequencies'
      DO K = 1, nK
         IF (Y2(K).le.0) THEN
            Fehler = 'S(q) <= 0'
            RETURN
            ENDIF
         ! nC(q): Ornstein-Zernike direct correlation
         Y3(K) = (Y2(K)-1) / Y2(K)
         ! Omega^2 = T * q^2 / m / S(q) :
         Omeg2(K) = OmRel(1) * Tm * X3(K)**2 / Y2(K)
         Fric(K)  = FrRel(1) * Tm * X3(K)**2 / Y2(K)
         ENDDO
      DO K = nK1, nK2
         ! Omega_s ohne S(q)
         Omeg2(K) = OmRel(2) * Tm * X3(K-nK)**2
         Fric(K)  = FrRel(2) * Tm * X3(K-nK)**2
         ENDDO

      dens = 12 * eta / twopi / diam**3
      Print *, ' found density[A-3] = ', dens

      Print *, ' set vertices'
      ! Franosch et al 1997b (7) & Fuchs et al 1998 (11)
      DO K = 1, nK
         qs = K - 0.5
         DO KK = 1, nK
            ks = KK - 0.5
            DO KKK = 1, nK
               ps = KKK - 0.5
               IF (dabs(qs-ks)+0.5.le.ps .and. ps.le.qs+ks-0.5) THEN

                  Vtx1(K,KK,KKK) = (delQ**3 / 8 / twopi**2 / dens) *
     * Y2(K)*Y2(KK)*Y2(KKK) * (ks*ps/qs**5) *
     * ((qs**2+ks**2-ps**2)*Y3(KK) + (qs**2+ps**2-ks**2)*Y3(KKK))**2

                  Vtx2(K,KK,KKK) = (delQ**3 / 4 / twopi**2 / dens) *
     * Y2(KK) * (ks*ps/qs**5) * ((qs**2+ks**2-ps**2)*Y3(KK))**2

               ELSE
                  Vtx1(K,KK,KKK) = 0
                  Vtx2(K,KK,KKK) = 0
                  ENDIF
               ENDDO ! KKK
            ENDDO ! KK
         ENDDO ! K

C  Nur DWF berechnen ?
      IF (jfout.eq.5) THEN ! Franosch et al 1997b III.A.
         Print *, ' now calculate f(q)'
         Kout = Kin
         DO K = 1, nK
            AuxY(K) = 1 ! Startwert
            ENDDO
         niter = 0
         ualt = 1
 31      CONTINUE
            niter = niter + 1
            IF (niter.gt.maxiter) THEN
               Fehler = ' did not converge'
               RETURN
               ENDIF
            ungenau = 0
            CALL Vertex (nK, Vtx1, AuxY, AuxY, AuxZ)
            DO K = 1, nK
               IF (AuxZ(K).le.-1) THEN
                  Fehler = 'F(f)<-1'
                  RETURN
                  ENDIF
               AuxZ(K) = AuxZ(K) / (1 + AuxZ(K))
               ungenau = dmax1 (ungenau, dabs(AuxZ(K)-AuxY(K)))
               ENDDO
            IF ((niter.le.100 .or.
     *          (niter.le.1000 .and. mod(niter,10).eq.0) .or.
     *          mod(niter,100).eq.0) .and. nK.ge.7)
     * Print *, niter, 'f(q) -> ', AuxZ(2), dquot0(ungenau,ualt)
            IF (ungenau.lt.genau) GOTO 39
            ualt = ungenau
            DO K = 1, nK
               AuxY(K) = AuxZ(K)
               ENDDO
            GOTO 31
 39      CONTINUE
         ! speichern :
         CALL OlfPutZ (jout, Kout, nZ, Z, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfPutXY0 (jout, Kout, nK, X3, AuxZ, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ! schon fertig !!
         GOTO 900
         ENDIF ! case f(q)

C  Integration im Zeitraum
      ! Anfangswerte fuer t=0 :
      Print *, ' set initial values'
      ht = ht0
      DO K = 1, nK2
         IF (Omeg2(K).le.0) THEN
            Fehler = 'Omega<=0 incompatible with Newtonian dynamics'
            RETURN
            ENDIF
         RedF(K) = 1 - Fric(K)*ht
         RedO(K) = Omeg2(K) * ht
         AuxC(K) = 0
         Phi(K,1) = 1
         ENDDO
      CALL Vertex (nK, Vtx1, Phi(1,1), Phi(1,1), Ern(1,1))
      IF (jcout.gt.1)
     *CALL Vertex (nK, Vtx2, Phi(1,1), Phi(nK1,1), Ern(nK1,1))

      ! "startwerte"
      Print *, ' set first block'
      DO is = 1, nCh
         AX(is) = (is-1)*ht
         ENDDO
      DO is = 2, nBL
         DO K = 1, nK2
            Psi(K,is) = RedF(K)*Psi(K,is-1) -
     *                  RedO(K)*(Phi(K,is-1)+AuxC(K))
            Phi(K,is) = ht*Psi(K,is)+Phi(K,is-1)
            ENDDO ! K

         CALL Vertex (nK, Vtx1, Phi(1,is), Phi(1,is), Ern(1,is))
         IF (jcout.gt.1)
     *   CALL Vertex (nK, Vtx2, Phi(1,is), Phi(nK1,is), Ern(nK1,is))

         DO K = 1, nK2
            AuxC(K) = 0
            DO it = 1, is-1
               AuxC(K) = AuxC(K) + ht/2 *
     *            (Ern(K,1+is-it)*Psi(K,it) + Ern(K,is-it)*Psi(K,it+1))
               ENDDO
            ENDDO ! K

         ENDDO ! is

      DO K = 1, nK2
         ! "initialisiere"
         DO is = 2, nBL
            dPhi(K,is) = (Phi(K,is-1) + Phi(K,is)) / 2
            dErn(K,is) = (Ern(K,is-1) + Ern(K,is)) / 2
            ENDDO
         ENDDO ! K

      DO K = 1, nK
         ! Speichern :
         DO is = 1, nBL ! sozusagen id=0
            SavX(is,K) = AX(is)
            SavY(is,K) = Phi(K+(jcout-1)*nK,is)
            ENDDO
         ENDDO ! K

      ! Hauptschleife:
      Print *, ' main loop'
      id = 1
 100  CONTINUE
         Print *, '  id = ', id
         ! "hauptalgorithmus"
         DO K = 1, nK2
            RedF(K) = 0.5*Fric(K) / ht / Omeg2(K)
            RedO(K) = 1.d0 / ht**2 / Omeg2(K)
            AuxA(K) = 1 + dErn(K,2) + 3*RedF(K) + 2*RedO(K)
            AuxB(K) = (1 - dPhi(K,2)) / AuxA(K)
            ENDDO
         DO is = nBL+1, nCh
c            Print *, '   is = ', is
            DO K = 1, nK2
               Phi(K,is) = 0
               Ern(K,is) = 0
               ! Faltungsintegral :
               i2 = (is+1)/2
               AuxC(K) = Ern(K,i2) * Phi(K,1+is-i2)
               DO it = 2, i2
                  AuxC(K) = AuxC(K) +
     *                 dErn(K,it) * (Phi(K,2+is-it) - Phi(K,1+is-it)) +
     *                 dPhi(K,it) * (Ern(K,2+is-it) - Ern(K,1+is-it))
                  ENDDO ! it
               DO it = i2+1, is-i2+1
                  AuxC(K) = AuxC(K) +
     *                 dPhi(K,it) * (Ern(K,2+is-it) - Ern(K,1+is-it))
                  ENDDO ! it
               AuxC(K) = AuxC(K) +
     *          RedO(K) * (-5*Phi(K,is-1)+4*Phi(K,is-2)-Phi(K,is-3)) +
     *          RedF(K) * (-4*Phi(K,is-1)+  Phi(K,is-2))
               AuxC(K) = AuxC(K) / AuxA(K)
               AuxY(K) = Phi(K,is-1)
               ENDDO ! K

            ! innermost loop (the bottleneck) :
 101        CONTINUE ! alt -> Z
            DO K = 1, nK
               AuxZ(K) = AuxY(K)
               ENDDO ! K
            CALL Vertex (nK, Vtx1, AuxZ(1), AuxZ(1), AuxY(1))
            ungenau = 0
            DO K = 1, nK ! neu -> Y
               AuxY(K) = -AuxC(K) + AuxB(K) * AuxY(K)
               ungenau = dmax1 (ungenau, dabs(AuxY(K)-AuxZ(K)))
               ENDDO ! K
            IF (ungenau.gt.genau) GOTO 101

            ! dito for the second set of correlators (coupled to the first) :
            IF (jcout.gt.1) THEN
 103        CONTINUE            ! alt -> Z
            DO K = nK1, nK2
               AuxZ(K) = AuxY(K)
               ENDDO ! K
            CALL Vertex (nK, Vtx2, AuxY(1), AuxZ(nK1), AuxY(nK1))
            ungenau = 0
            DO K = nK1, nK2 ! neu -> Y
               AuxY(K) = -AuxC(K) + AuxB(K) * AuxY(K)
               ungenau = dmax1 (ungenau, dabs(AuxY(K)-AuxZ(K)))
               ENDDO ! K
            IF (ungenau.gt.genau) GOTO 103
            ENDIF ! Sjoegren

            ! Ergebnis :
            DO K = 1, nK2
               Phi(K,is) = AuxY(K)
               ENDDO ! K
            CALL Vertex (nK, Vtx1, AuxY(1), AuxY(1), Ern(1,is))
            IF (jcout.gt.1)
     *      CALL Vertex (nK, Vtx2, AuxY(1), AuxY(nK1), Ern(nK1,is))

            ENDDO ! is

         DO K = 1, nK
            ! "speichere"
            DO is = nBL+1, nCh ! fuer id=1..nd
               SavX((id-1)*(nBL-1)+is,K) = AX(is)
               SavY((id-1)*(nBL-1)+is,K) = Phi(K+(jcout-1)*nK,is)
               ENDDO
            ENDDO ! K

         ! "dezimierungen"
         IF (id.gt.nd) GOTO 199
         ht = 2*ht
         DO is = 2, nCh
            AX(is) = (is-1) * ht
            ENDDO
         DO K = 1, nK2
            DO is = 2, nC4+1
               dPhi(K,is) = (dPhi(K,2*is-2) + dPhi(K,2*is-1)) / 2
               dErn(K,is) = (dErn(K,2*is-2) + dErn(K,2*is-1)) / 2
               ENDDO
            DO is = nC4+2, nBL
               dPhi(K,is) = (Phi(K,2*is-3)+4*Phi(K,2*is-2)+
     *                       Phi(K,2*is-1))/6
               dErn(K,is) = (Ern(K,2*is-3)+4*Ern(K,2*is-2)+
     *                       Ern(K,2*is-1))/6
               ENDDO
            DO is = 2, nBL
               Phi(K,is) = Phi(K,2*is-1)
               Ern(K,is) = Ern(K,2*is-1)
               ENDDO
            ENDDO ! K
         id = id + 1
         GOTO 100
 199  CONTINUE
      Print *, ' end of iterations'

C  Loop output spectra :
      Print *, ' now convert to output function'
      ns = nd*(nBL-1)+nCh
      IF (jfout.eq.1) THEN
         istep = max0 (1, (ns-1)/nout+1)
         nf = ns / istep
      ELSEIF (jfout.lt.4) THEN
         nf = min0 (nout, ns)
      ELSE
         nf = min0 ((nout-1)/2, ns)
         ENDIF

      DO K = 1, nK
         Kout = (Kin-1)*nK + K

C  Conversion to output function :
         IF     (jfout.eq.1) THEN
            DO if = 1, nf
               X1(if) = SavX(if*istep,K)
               Y1(if) = SavY(if*istep,K)
               ENDDO

         ELSEIF (jfout.ge.2 .and. jfout.le.4) THEN
            ! set log grid
            xma = 1/ht0 ! war 1/SavX(2, K)
            xll = dlog (xma/xmi)
            DO if = 1, nf
               X2(if) = xmi * dexp((if-1)*xll/(nf-1))
               ENDDO

            ! get Phi(infty)
            Phinfty = 0
            DO is = nBL+1, nCh
               Phinfty = Phinfty + SavY(nd*(nBL-1)+is,K)
               ENDDO
            Phinfty = Phinfty / (nCh-nBL)
            Phigma  = 0
            DO is = nBL+1, nCh
               Phigma = Phigma + (SavY(nd*(nBL-1)+is,K)-Phinfty)**2
               ENDDO
            Phigma = Phigma / (nCh-nBL)
            IF (Phigma.gt.1d-12) THEN
               Print *, ' Phi(infty) Phi(sigma) ', Phinfty, Phigma
               CALL NiceNum (Phigma, aux, i2)
               Fehler = 'Phi(t) converges only within +- '//aux(1:i2)
               RETURN
               ENDIF
            DO is = 1, (nd+1)*(nBL-1)+1
               SavY(is,K) = SavY(is,K) - Phinfty
               ENDDO

            ! blockwise Fourier transform
            DO if = 1, nf
               Y2(if) = 0
               ENDDO
            DO id = 0, nd
               n1 = id*(nBL-1)+1
               nn = id*(nBL-1)+nBL
               ht = (SavX(nn,K) - SavX(n1,K)) / (nBL-1)
               DO if = 1, nf
                  ! lowest-order Filon integration
                  theta = X2(if) * ht
                  IF (dabs(theta).lt.1d-40) THEN
                     alpha = 1
                     delta = 0
                  ELSE
                     alpha = ( dsin(theta/2) / (theta/2) )**2
                     delta = (theta-dsin(theta)) / theta**2
                     ENDIF
                  sumY = 0
                  IF (jfout.eq.3 .or. jfout.eq.4) THEN ! real part S(w) = X''/w
                     DO is = n1+1, nn-1
                        sumY = sumY + SavY(is,K) *
     *                         dcos(SavX(is,K)*X2(if))
                        ENDDO ! is
                     sumY = alpha * ( sumY +
     *                  ( SavY(n1,K) * dcos(SavX(n1,K)*X2(if))
     *                  + SavY(nn,K) * dcos(SavX(nn,K)*X2(if)) ) / 2 )
     *                       + delta *
     *                  ( - SavY(n1,K) * dsin(SavX(n1,K)*X2(if))
     *                    + SavY(nn,K) * dsin(SavX(nn,K)*X2(if)) )
                  ELSEIF (jfout.eq.2) THEN  ! imaginary part of S(w) = X'/w
                     DO is = n1+1, nn-1
                        sumY = sumY + SavY(is,K) *
     *                         dsin(SavX(is,K)*X2(if))
                        ENDDO ! is
                     sumY = alpha * ( sumY +
     *                  ( SavY(n1,K) * dsin(SavX(n1,K)*X2(if))
     *                  + SavY(nn,K) * dsin(SavX(nn,K)*X2(if)) ) / 2 )
     *                       + delta *
     *                  ( + SavY(n1,K) * dcos(SavX(n1,K)*X2(if))
     *                    - SavY(nn,K) * dcos(SavX(nn,K)*X2(if)) )
                  ELSE
                     Fehler = 'PROGRAM ERROR / Re vs Im ???'
                     RETURN
                     ENDIF ! Re vs. Im
                  Y2(if) = Y2(if) + ht * sumY
                  ENDDO ! if
               ENDDO ! id

            IF (jfout.eq.2 .or. jfout.eq.3) THEN
               DO if = 1, nf
                  X1(if) = X2(if)
                  Y1(if) = Y2(if) * X2(if) ! susceptibility
                  ENDDO
               nout = nf
            ELSEIF (jfout.eq.4) THEN
               DO if = 1, nf
                  X1(nf+1-if) = -X2(if)
                  X1(nf+1+if) =  X2(if)
                  Y1(nf+1-if) = Y2(if) * 2 / twopi ! scattering law
                  Y1(nf+1+if) = Y2(if) * 2 / twopi ! scattering law
                     ! Faktor 1/pi 22okt98
                  ENDDO
               ! elastic line (Debye-Waller factor or interpolation) :
               X1(nf+1) = 0
               Y1(nf+1) = dmax1 (Phinfty * 2 / (X1(nf+2) - X1(nf)),
     *                           (Y1(nf+2) + Y1(nf))/2 )
               nout = 2*nf+1
            ELSE
               Fehler = 'Scattering law or susceptibility ?'
               RETURN
               ENDIF

         ELSE
            Fehler = 'PROGRAM ERROR/ unexpected output choice'
            ENDIF

         ! Save result :
         Z(nZ+1) = X3(K) ! Q[A-1]
         CALL OlfPutZ (jout, Kout, nZ+1, Z, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfPutXY0 (jout, Kout, nout, X1, Y1, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! K

 900  ENDDO ! Kin
      CALL OlfClos (jout, Kout, Fehler)
      ENDDO ! J

      END ! MCT_Full

      SUBROUTINE Vertex (nK, Vtx, Phi1, Phi2, Ern)
C     --------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INTEGER        K, KK, KKK, nK
      REAL*8         Vtx(MK3V,MK3V,MK3V), Phi1(MK3V), Phi2(MK3V),
     *               Ern(MK3V)

      DO K = 1, nK
         Ern(K) = 0
         DO KK = 1, nK
            DO KKK = 1, nK
               Ern(K) = Ern(K) + Vtx(K,KK,KKK) * Phi1(KK) * Phi2(KKK)
               ENDDO ! KKK
            ENDDO ! KK
         ENDDO ! K

      END ! Vertex
