C  ====================================================================
C
C      Program  IDA   :  Inelastic Data Analysis
C      Modul    i72   :     tof kinematics
C
C  ====================================================================

C     Contents :
C
C        1.  Neutron kinematics:
C               yQ_of_0, yQ_of_W, Th2_of_W, E_of_l,
C               Tau_of_W, W_of_Tau, T_of_W, W_of_T, yKi_of_W
C
C        2.  Constant q interpolation:
C               CoQ
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  ====================================================================
C  i72 / 1 :   Neutron kinematics
C  ====================================================================

      REAL*8 FUNCTION yQ_of_0 ( E0, Th2)
C     ----------------------------------
         ! Calculate elastic wavenumber q for given scattering angle 2th
      REAL*8  K0, E0, Th2
      IF (E0.le.0) THEN
         Print *, 'yQ_of_0/ SEVERE ERROR/ E0<=0'
         yQ_of_0 = 0
         RETURN
         ENDIF
      K0 = dsqrt( E0 / 2.0723 )
      IF (Th2.lt.0) THEN
         Print *, 'yQ_of_0/ WARNING/ Th2<0'
         ENDIF
      yQ_of_0 = 2 * K0 * dsind ( Th2/2 )
      END ! yQ_of_0

      REAL*8 FUNCTION Th2_of_Q0 ( E0, Q)
C     ----------------------------------
         ! Calculate scattering angle 2th for given elastic wavenumber q
      REAL*8  K0, E0, Q, args
      K0 = dsqrt( E0 / 2.0723 )
      args = Q / 2 / K0
      IF (args.lt.-1 .or. args.gt.1) THEN
         Print *, 'Th2_ofQ0/ WARNING/ q out of elastic range'
         Th2_of_Q0 = 0
      ELSE
         Th2_of_Q0 = 2 * dasind ( args )
         ENDIF
      END ! Th2_of_Q0

      REAL*8 FUNCTION yQ_of_W ( W, E0, Th2 )
C     --------------------------------------
         ! Calculate inelastic q for scattering with 2th,w
      REAL*8  K, K0, W, E0, Th2
      IF (E0+W.lt.0.) THEN
         yQ_of_W = 0.
         Print *, 'invalid parameters in yQ_of_W : E0, W ', E0, W
         RETURN
         ENDIF
      K    = dsqrt( (E0 + W) / 2.0723 )
      K0   = dsqrt(  E0      / 2.0723 )
      yQ_of_W = dsqrt ( K**2 + K0**2 - 2*K*K0*dcosd(Th2) )
      END ! yQ_of_W

      REAL*8 FUNCTION Th2_of_W ( W, E0, Q )
C     -------------------------------------
         ! Calculate 2th for inelastic scattering with q,w
      REAL*8  Q, K, K0, W, E0, arg
      IF (E0+W.lt.0.) THEN
         Print *, ' Absturz in Th2_of_W : E0, W ', E0, W
         STOP
         ENDIF
      K    = dsqrt( (E0 + W) / 2.0723 )
      K0   = dsqrt(  E0      / 2.0723 )
      arg  = ( K**2 + K0**2 - Q**2 ) / (2*K*K0)
      IF (dabs(arg).gt.1.) THEN
         Print *, ' K K0 Q : ', K, K0, Q
         CALL Absturz ('Th2_of_W', '|arg|>+1')
      ELSE
         Th2_of_W = dacosd (arg)
         ENDIF
      END ! Th2_of_W

      REAL*8 FUNCTION E_of_l (lambda)
C     -------------------------------
         ! Convert neutron wavelength to energy
      REAL*8 lambda
      E_of_l = 81.805d0 / lambda**2
      END ! E_of_l

      REAL*8 FUNCTION Tau_of_W (W, E0)
C     --------------------------------
         ! time of flight for a neutron with energy E0+hquer*w
         ! since 21jun91 in milliseconds/meter, not microseconds !
      REAL*8 W, E0
      IF (W.le.-E0) THEN
         CALL Gong (17)
         Print *, ' BAD INPUT in Tau_of_W : w <= -E0'
         Tau_of_W = 1.
         RETURN
         ENDIF
      Tau_of_W = dsqrt ( 5.2271d0 / ( E0 + W ) )
      END ! Tau_of_W

      REAL*8 FUNCTION W_of_Tau (Tau, E0)
C     ----------------------------------
         ! energy gain E-E0 for a neutron with time of flight tau.
         ! msec/m -> meV
      REAL*8 Tau, E0
      IF (E0.eq.0.) THEN
         CALL Gong (17)
         CALL Absturz ('W_of_Tau', 'E0 = 0.0')
         ENDIF
      W_of_Tau =  5.2271d0 / Tau**2 - E0
      END ! W_of_Tau

      REAL*8 FUNCTION T_of_W (W)
C     --------------------------
         ! just t = hvoll/w [nsec / meV]
      REAL*8 W
      IF (W.le.0.) CALL Absturz ('T_of_W', 'W =< 0')
      T_of_W = 4.13570d-3 / W
      END ! T_of_W

      REAL*8 FUNCTION W_of_T (T)
C     --------------------------
         ! just w = hvoll/t [meV / nsec]
      IF (T.le.0.) CALL Absturz ('W_of_T', 'T =< 0')
      W_of_T = 4.13570d-3 / T
      END ! W_of_T

      REAL*8 FUNCTION DKi_of_W ( W, E0 )
C     ----------------------------------
         ! DKi = Ki(W) - Ki(W=0)
      REAL*8  K0, W, E0
      IF (E0+W.lt.0.) THEN
         Print *, ' Absturz in DKi_of_W : E0, W ', E0, W
         STOP
         ENDIF
      K0   = dsqrt(  E0      / 2.0723 )
      DKi_of_W = dsqrt (W/2.0723+ K0**2 )-K0
      END ! DKi_of_W

C  ====================================================================
C  i72 / 2 :   Constant q interpolation
C  ====================================================================

            ! FK  15feb07 : option 4 activated
            ! JWu  7nov96 : Integrated in Ida (429 -> 328 lines)
            ! JWu   feb92 : Public version for Saclay
            ! JWu 24aug92 : Finite sample width (see C2,57)
            ! JWu 16may92 : Format binary92
            ! JWu  5sep91 : Double linear interpolation (see C2,28)
            ! JWu  2sep91 : Real*8, extended IED format, dialogue.
            ! FRi   jun91 : Version InGridM
            ! FRi   jan91 : Program InGrid by Francois Rieutord

      SUBROUTINE CoQ (nJList, JList, Fehler)
C     --------------------------------------

      IMPLICIT NONE

      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      REAL*8        yQ_of_0, yQ_of_W

      INCLUDE      'i_wrk.f'

      INTEGER       maxQsample
      PARAMETER    (maxQsample=100)

      CHARACTER*(*) Fehler
      INTEGER       nJList, JList(*), lj, jin, jout

      REAL*8        Qsoll(MK), Z(MZ), ZZ(MK)
      INTEGER       iOptQ, iOptE, nCmin, n2, nZ,
     *              nK, nK2, nKsoll, K, K2, Ka, Kb, Kout, KPS,
     *              i, n1, iaa, iab, iba, ibb, iout,
     *              nQsample, iP, ih1, ih2
      REAL*8        E0, Emin, Emax, Esa, Esb, Erel, Qstep, Qextn,
     *              Qcut, Qmax, Qwidth, wtot, facX, Xs, Ys, Ds, Trel,
     *              Ts, Yas, Ybs, Das, Dbs, Qi, Qf, Qs, arg
      REAL*8        DQsample(maxQsample), WQsample(maxQsample)
      CHARACTER*40  Co, Un, UnX, UnZ, com, comQ, comE, docE, h1, h2

      DATA          iOptQ /1/, iOptE /2/, nCmin /7/, Qcut /.01d0/

C  --------------------------------------------------------------------
C     Initializations, Begin dialogue, Input
C  --------------------------------------------------------------------

      IF (nJList.lt.1) THEN
         Fehler = ' '
         RETURN
         ENDIF

      Print *, ' choose energy grid :'
      CALL SetGridChoice (.true., docE, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO lj = 1, nJList
         jin = JList(lj)

         CALL OlfHeadDup (jin, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL SetGridJ (jin, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         E0 = rOlfGG (jin, 'E0', 'meV', Fehler)

         CALL OlfCnuG (jin, 'z1', Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (Co.ne.'2th') THEN
            Fehler = ' WARNING/ z-coordinate not 2th but '//Co
            RETURN
            ENDIF
         CALL OlfCnuP (jout, 'z1', 'q', 'A-1', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfCnuG (jin, 'y', Co, Un, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF     (Co.eq.'S(2th,w)') THEN
            CALL OlfCnuP (jout, 'y', 'S(q,w)', Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSEIF (Co.eq.'S~(2th,w)') THEN
            CALL OlfCnuP (jout, 'y', 'S~(q,w)', Un, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ELSE
            Fehler = 'y-coordinate not S(2th,w) but '//Co
            RETURN
            ENDIF

         CALL OlfCnuG (jin, 'x', Co, UnX, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL UnitConv (UnX, 'meV', facX, Fehler)
         IF (Fehler.ne.'&ff') RETURN

C  --------------------------------------------------------------------
C     Prepare q,w - grid
C  --------------------------------------------------------------------

C  Set 2-dimensional arrays :

c  Compiler-Absturz unter DEC Unix Nov96 :
c         IF (MWrk3dim.lt.3) THEN
c            Fehler = 'Wrk3dim insufficient for CoQ'
c            RETURN
c            ENDIF

         Qmax = 0
         DO K = 1, nK

            CALL OlfGetSpe (jin, K, nZ, Z, NofK(K), X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (nZ.ne.1) THEN
               Fehler = 'PROVISORISCH/ nZ<>1 verboten'
               RETURN
               ENDIF
            ZZ(K) = Z(1)

            CALL SetGridK (K, n1, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            DO i = 1, NofK(K)
               Wrk3dim(i,K,1) = X(i)
               Wrk3dim(i,K,2) = Y(i)
               Wrk3dim(i,K,3) = D(i)
               Qmax = dmax1 (Qmax, yQ_of_W (X(i)*facX, E0, ZZ(K)))
               ENDDO
            ENDDO

C  Set Q - grid :
         Print '(a,g10.4)', ' elast. q min = ', yQ_of_0 (E0, ZZ(1))
         Print '(a,g10.4)', ' elast. q max = ', yQ_of_0 (E0, ZZ(nK))
         Print '(a,g10.4)', ' inel.  q max = ', Qmax

         IF (lj.eq.1) THEN
            Print *
            Print *, ' interpolation in q :'
            Print *, '    (1) equidistant q'
            Print *, '    (2) q(w=0)'
            Print *, '    (3) q(w=0) + inelastic extension'
            Print *, '    (4) enter individually'
            iOptQ = iAskDMu (' Option', iOptQ, 0, 4)

            IF     (iOptQ.eq.1) THEN
               CALL rAskGrid ('q', Qsoll, MK, nKsoll)
               comQ = ' q = eqd.'
            ELSEIF (iOptQ.eq.2 .or. iOptQ.eq.3) THEN ! 10jan00
               CALL OlfGet1ZofK (jin, 1, nKsoll, ZZ, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               DO K = 1, nKsoll
                  Qsoll(K) = yQ_of_0 (E0, ZZ(K))
                  ENDDO
               IF (iOptQ.eq.2) THEN
                  comQ = ' q(w=0)'
               ELSE ! extend the q-range for inelastic scattering
                  IF (nKsoll.lt.2) THEN
                     Fehler = ' nKsoll < 3'
                     RETURN
                     ENDIF
                  IF (Qextn.le.0) Qextn = Qmax ! default for first call
                  Qextn = rAskD (' Extend q-range up to', Qextn)
                  Qstep = Qsoll(nKsoll) - Qsoll(nKsoll-1)
                  IF (Qstep.le.0) THEN
                     Fehler = 'q-step <= 0'
                     RETURN
                     ENDIF
                  DO WHILE (Qsoll(nKsoll)+Qstep.lt.Qextn
     *                      .and. nKsoll.lt.MK)
                     nKsoll = nKsoll + 1
                     Qsoll(nKsoll) = Qsoll(nKsoll-1)+Qstep
                     ENDDO
                  comQ = ' q(w=0) + inelastic extension'
                  ENDIF
            ELSEIF (iOptQ.eq.4) THEN
               nKsoll = iAskDmu ('number of points', nKsoll, 0, 25)
               CALL rAskArray (Qsoll, MK, nKsoll)
               comQ = ' q = entered'
            ELSE
               Fehler = 'invalid choice'
               RETURN
               ENDIF

            nQsample = iAskDMu (' No. of slices', nQsample,
     *                          1, maxQsample)
            IF (nQsample.gt.1) THEN
               Qwidth = rAskDMu (' Radius of arc for Q sampling',
     *              Qwidth, 0.d0, Qmax)
               IF (Qwidth.le.0) THEN
                  Fehler =
     * 'if no_of_slices > 0 then sampling_width > 0 required'
                  RETURN
                  ENDIF
               wtot = 0
               DO ip = 1, nQsample
                  DQsample(ip) = (2*ip-nQsample-1) * Qwidth /
     *                           (nQsample-1)
                  ! weigth: arc with radius Qwidth
                  WQsample(ip) = dsqrt0 (1 - (DQsample(ip)/Qwidth)**2)
                  wtot = wtot + WQsample(ip)
                  ENDDO
               DO ip = 1, nQsample
                  WQsample(ip) = WQsample(ip) / wtot
                  ENDDO

                     ! S(q,w) will be calculated as a weighted average
                     ! <WQsample(i) * S(Q_i,w)> with Q_i = Q + DQsample(i)

               CALL NiceNum (Qwidth,  h1, ih1)
               CALL Append (comQ, ' '//h1(1:ih1)//', #='//
     *                      cl4(nQsample) )
            ELSE
               DQsample(1) = 0.
               WQsample(1) = 1.
               ENDIF

            ENDIF ! first file

C  Further set-up :
         IF (lj.eq.1) THEN
            nCmin = iAskDMu (' Minimal number of channels per Q',
     *                       nCmin, 1, n1)
            ENDIF

C  Documentation :
         CALL Compose2 (com, 'coq/ '//comQ, comE)
         CALL OlfComAdd (jout, 'Q', com, Fehler)

C  --------------------------------------------------------------------
C     Interpolation
C  --------------------------------------------------------------------

            ! The interpolation routine is envelopped by a double
            ! running over the new (q,w) grid. Nevertheless, the
            ! interpolation algorithm itself (C2, 29) is formulated
            ! in (2th,w), not in (q,w).

         Kout = 0
         DO K = 1, nKsoll
            iout = 0
            DO i = 1, n1

               Ys = 0
               Ds = 0

               DO KPS = 1, nQsample       ! loop over Q - samples
                  Qs = Qsoll(K) + DQsample(KPS)
                  Xs = X1(i)

C  Find the mesh that contains the point (K,i) :

c                       2Th  Kb :  iba    ibb
c                        ^
c                        |   Ks        is
c                        |
c                        |   Ka :  iaa    iab
c                           - - - - - - - - - - >  E


C  Determine Ka so that 2Th(Ka) <= 2Th(i,K) < 2Th(Ka+1) :

                  IF (Xs*facX.le.-E0) THEN
                     Print *, ' K i E0 Xs ', K, i, E0, Xs
                     Fehler = ' Impossible energy scale'
                     RETURN
                     ENDIF

                  ! determine 2Theta(Qs) [incorporate Th2_of_W] :
                  Qi   = dsqrt(  E0          / 2.0723 )
                  Qf   = dsqrt( (E0+Xs*facX) / 2.0723 )
                  arg  = ( Qi**2 + Qf**2 - Qs**2 ) / (2*Qi*Qf)
                  IF (dabs(arg).gt.1.) GOTO 28 ! outside the dynamical range
                  Ts = dacosd (arg)

                  Ka = irPosOpt (ZZ, nK, Ts, 'r', Ka)
                  Kb = Ka + 1
                  IF (Ka.lt.1 .or. Kb.gt.nK) GOTO 28 ! outside the old area

C  Now for Ka,Kb determine ia,ib so that X(i_,K_) <= X1(i) < X(i_+1,K) :

                  iaa= irPosOpt (Wrk3dim(1,Ka,1), NofK(Ka), Xs,
     *                           'r', iaa+1)
                  iab= iaa + 1
                  IF (iaa.lt.1 .or. iab.gt.NofK(Ka)) GOTO 28
                  iba= irPosOpt (Wrk3dim(1,Kb,1), NofK(Kb), Xs,
     *                           'r', iba+1)
                  ibb= iba + 1
                  IF (iba.lt.1 .or. ibb.gt.NofK(Kb)) GOTO 28

C  The four neighbours are found, the point is inside the old grid :

                  ! determine Erel :
                  Trel = (Ts-ZZ(Ka)) / (ZZ(Kb)-ZZ(Ka))
                  Esa  = (1-Trel)*Wrk3dim(iaa,Ka,1) +
     *                   Trel*Wrk3dim(iba,Kb,1)
                  Esb  = (1-Trel)*Wrk3dim(iab,Ka,1) +
     *                   Trel*Wrk3dim(ibb,Kb,1)
                  Erel = (Xs-Esa) / (Esb-Esa)

                  ! now use Erel to interpolate in E :
                  Yas  = (1-Erel)*Wrk3dim(iaa,Ka,2) +
     *                   Erel*Wrk3dim(iab,Ka,2)
                  Das  = (1-Erel)*Wrk3dim(iaa,Ka,3) +
     *                   Erel*Wrk3dim(iab,Ka,3)

                  Ybs  = (1-Erel)*Wrk3dim(iba,Kb,2) +
     *                   Erel*Wrk3dim(ibb,Kb,2)
                  Dbs  = (1-Erel)*Wrk3dim(iba,Kb,3) +
     *                   Erel*Wrk3dim(ibb,Kb,3)

                  ! finally interpolate in 2Th :
                  Ys = Ys + WQsample(KPS) * ( (1-Trel)*Yas + Trel*Ybs )
                  Ds = Ds + WQsample(KPS) * ( (1-Trel)*Das + Trel*Dbs )

                  ENDDO ! KPS

C  The point was inside the old grid for all KPS :
               iout       = iout + 1
               X2(iout) = Xs
               Y2(iout) = Ys
               D2(iout) = Ds

 28            CONTINUE
               ENDDO ! i

C  Save the new spectrum :

            IF (iout.ge.nCmin) THEN
               Kout = Kout + 1
               CALL OlfPutSpe (jout, Kout, 1,Qsoll(K), iout,
     *                         X2, Y2, D2, Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDIF

            ENDDO ! outer loop (K)

         CALL OlfClos (jout, Kout, Fehler)

         ENDDO ! lj

      END ! CoQ
