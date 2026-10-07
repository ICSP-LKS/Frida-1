C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i94   :     simulation / multiple scattering
C
C  ====================================================================

C  History :
C      JWu 7-8/99 :  integrated in IDA as i94.f
C      JWu   3/96 :  real*4 suppressed, resolution from file
C      JWu   2/96 :  backscattering: arrays eps<>beta
C      JWu  11/94 :  geometry-dependent parts rewritten
C      JWu   7/94 :  made it running again, simplified, splitted
C      JWu 4-6/92 :  dialogue, file in/output, listing,
C                     array dimensions, variable names,
C                     program structure, arithmetic exceptions, ...
C      Program MSCAT by J.D.R.Copley et al.
C         (Comput.Phys.Commun. 7,289(1974),...)
C      Version MSCAT84 copied from J.B.Suck [C6.VISUC.SIERK]MSCAT.TXT;94
C         (as of 2-apr-92)
C      16.02.2026 Artem Panchenko: Corrected several line breaks

C  ====================================================================
      SUBROUTINE MScat (nJList, JList, Fehler)
C  ====================================================================

C  --------------------------------------------------------------------
C     Some explanations, definitions, ..
C  --------------------------------------------------------------------

      ! Types of Collisions :
      !     IE = inelastic
      !     SC = scattering = EL + IE
      !     AB = absorption
      !     TT = total = SC + AB

      ! Coordinate systems :
      !     T-CS   = Target CS,
      !     DE-CS  = CS for Bragg scattering
      !     Det-CS = detector CS, used only in CrossDE.

      ! Units :
      !     cm    for lengths; cross sections are converted to cm-1
      !     usec  for time of flight
      !     A-1   for momenta, and A^2 for mean square displacement <u^2>
      !     meV   for enErgIEs, and internally for Temperature,
      !               widely also for momenta

C  --------------------------------------------------------------------
C     Array dimensions :
C  --------------------------------------------------------------------

      IMPLICIT NONE ! 11jan00

      INCLUDE           'i_dim.f'
      INCLUDE           'i_wrk.f'
      INCLUDE           'l_def.f'

      INTEGER    MAC, MBC, MTC, MEC, MMat, Mscor, Mcol, Msqw
      PARAMETER (
     *           MAC   = 120,       ! max # internal alpha meshes
     *           MBC   = 480,       ! max # internal beta  meshes
     *           MTC   = 128,       ! max # angles for output
     *           MEC   = MBC,       ! max # enErgIEs for output
     *           Mmat  =   3,       ! max # Target sections DON'T CHANGE !!
     *              ! ATTENTION : a fixed Mmat=3 is assumed in some SUBROUTINEs
     *           Mscor = 2*MMat+1,  ! # Scoring classes
     *                              !     (singl/mat, mult/mat, mult/mix)
     *           Mcol  =   9,       ! # classes for collision statistics
     *           Msqw  = Mscor+MMat)! # Mscor + Mmat (for ideal scattering law)

C  --------------------------------------------------------------------
C    Local variables :
C  --------------------------------------------------------------------

      INTEGER        JList(*),           ! IDA -> input file no.
     *               JOut(Mscor),        ! output -> IDA file nos.
     *               JOutF(9),           ! auxiliary output -> IDA file nos.
     *               NUnscored(MEC)      ! accu: # neutrons lost from score
      REAL*8         Xout(MC),           ! output -> IDA
     *               SQW(MC,Msqw,2),     ! original scoring groups
     *               SIO(MC,Mscor,2),    ! output scoring groups
     *               AC(MAC), BC(MBC),   ! internal alpha,beta-mesh
     *               AI(MAC), BI(MBC),   ! input alpha,beta-mesh
     *               Sint(MAC),          ! dummy -> GetIE
     *               Eext(MBC),          !
     *               SAB(MAC,MBC,MMat),  ! scattering law
     *               FS(MAC,MBC,MMat),   ! integrated ..
     *               PMAR(MBC,MEC,MMat), ! marginal ..
     *               SiMatSC(MMat),      ! bound scattering section
     *               SiMatAB(MMat),      ! absorption section for 2200km/s
     *               SiErgIE(MBC,MMat),  ! scattering cross section vs. erg
     *               SiErgTT(MBC,MMat),  ! total cross section ...
     *               SigSC(MMat),        ! scattering section at given erg
     *               SigTT(MMat),        ! total section ...
     *               Travel(MMat),       ! distance remains to travel
     *               DistSect(MMat)      ! distance through section
      REAL*8         RSG(MEC,MTC,Mscor), ! accu: response per group
     *               RSA(MEC,MTC,Mscor), ! accu: total response
     *               RSQ(MEC,MTC,Mscor), ! accu: squared ..
     *               TC(MTC),            ! detector angles
     *               DirDet(3,MTC),      ! detector position cosines
     *               TimOut(MEC),        ! output time mesh
     *               DTimOut(MEC),       ! widths of ...
     *               EC(MEC),            ! output? energy mesh
     *               DEC(MEC),           ! widths of ... (?)
     *               THIGH(MEC),         ! dummy -> IndexTim
     *               WgtCol(Mcol)        ! accu: weight vs collision-order

      CHARACTER*20   hnum, hnu2, NamMat(MMat)
      CHARACTER*40   UnitErg, UnitSqw, fil, tit
      CHARACTER*80   aus
      CHARACTER      Fehler*(*)

      ! variables for COMMON blocks :
      REAL*8         TempE, ConvVof, ErgN, TimN, PosN, DirN
      INTEGER        nPRT, jTrace
      LOGICAL        qHasSectn

      ! subroutines :
      INTEGER        Index1, Index2
      REAL*8         Alpha_BEC, Alpha_Base !Artem: Replace with rand function from slatec/fnlib/rand.f:  G05CAF, G05DAF
      INTEGER        count, rate, maxc !Artem add for rand()
      REAL           ran, rseed !Artem add for rand()
      real, external :: rand !Artem add for rand()

      ! local variables :
      INTEGER        j, nJList, Kout, jField, iRanStart, nHi, nGr, iHi,
     *               iGr, nScEvMax, nScEvRc,
     *               I, iBC, nAI, nBI, nAC, nBC, iB, iB1, iHig, iLow,
     *               IEfirst, nEC, iE, iE0, nTC, iA, iT, iC, IS, ir1,
     *               ir2, iStyle, KTim0, KEL, IndexTim, M1, M2min, M2,
     *               jMat, nMat, iRefMat, InMat, InOneMat,
     *               nLossA, nLossBL, nLossBH, nImp, nTrc, NnoHit,
     *               IScor, II, Iscout, isco, ij, nScor,
     *               ih, ih1, ih2, lf
      REAL*8         twopi, BOLTZ, WgtCut, Wgt, Erg0, Erg0width, TempK,
     *               DistSD, TARV, VL, XL, CO, SI,
     *               wMin, wMax, ErgMin, ErgMax, E, Erg1, PSR, Eps,
     *               Epsnew,
     *               Pint, t0, DTim0, E0, DistCol, DisImp, DisThru,
     *               aLow1, aHig1, aLow2, aHig2, aHig, aLow, alpha,
     *               A, B, w, bl, br, a1, a2, beta, prob, wr, wr1, wr2,
     *               fB, b1, b2, fr1, fr2, fE, s12, s1, s2, Fmax, Fmin,
     *               Bran, Fran, aziS, aziC, exneg, RelRed,
     *               frac, survive, Transm, TarThick, wLossA, wLossBL,
     *               wLossBH, PE, FacY, SCO, FacX, dummy, SIE,
     *               FactA, FactR, FactQ,
     *               aux
      LOGICAL        qBetaEqProb, qIgnToF, qSymOut, qHit

C  --------------------------------------------------------------------
C    Global variables (COMMON blocks) :
C  --------------------------------------------------------------------

      COMMON /Convert/   TempE, ConvVof
      COMMON /Current/   ErgN, TimN, qHasSectn(MMat)
      COMMON /NeutPos/   PosN(3), DirN(3)
      COMMON /Trace/     nPRT, jTrace

C  --------------------------------------------------------------------
C     MS /  Initializations ( DATA and first instructions ) :
C  --------------------------------------------------------------------

      DATA NamMat / 'sample', 'container', 'absorber' /
      DATA          twopi   /6.2831853/,
     *              BOLTZ   /86.1733d-3/          ! meV / K
      DATA WgtCut / 3d-2 /

      ConvVof = 522.71                        ! [meV] = C / [usec/cm]^2

C  --------------------------------------------------------------------
C     MS /  Parameters from model_file :
C  --------------------------------------------------------------------

      IF     (nJList.lt.1) THEN
         Fehler =' no file given/ input must contain scattering law(s)'
         RETURN
      ELSEIF (nJList.gt.MMat) THEN
         Fehler = ' too much files given/ input are scattering laws of'
     *        //' sample sections'
         RETURN
         ENDIF

      j = JList(1)

      Erg0 = rOlfGG (j, 'E0', 'meV', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      TempK = rOlfGG (j, 'T', 'K', Fehler)
      TempE  = TempK * BOLTZ
      IF (Fehler.ne.'&ff') RETURN

      CALL tOlfG (j, 'fil', fil, Fehler)
      CALL tOlfG (j, 'tit', tit, Fehler)
      lf = lenU(fil)
      IF (Fehler.ne.'&ff') RETURN

C  --------------------------------------------------------------------
C     MS /  Interactive input :
C  --------------------------------------------------------------------

      IF (Fehler.ne.'&ff') RETURN

      CALL Sage (' ')
      CALL Sage (' 1/   Program control and Monte-Carlo setup')

C  Output files :
      CALL Sage (' 1/1  Output')

      nPRT = 16
      CALL OpenDatFile (nPRT, 'MScat.prt', 'lis', 'e',
     *                  'seq', 'for', 0, Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'cannot open lineprinter'
         RETURN
         ENDIF

      ! Output header :
      Write (nPRT,'(a/)') 'IDA / MScat'
      Write (nPRT,'(76(1h=)/a/76(1h=)/)') 'SETUP'
      Write (nPRT,'(a/17(1h-)/)') 'Program Control :'

      jField = iAskDMu (' Level for Array printing', jField, 0, 5)
      jTrace = iAskDMu (' Level for Neutron tracing', jTrace, 0, 9)

C  Monte-Carlo :
      CALL Sage (' 1/2  Monte-Carlo setup')

      ! Initialize NAG random number generator G05CAF/DAF/.. :
      iRanStart = iAskD (
     * 'Start value for random generator (0-> from clock)?', iRanStart)
      IF (iRanStart.eq.0) THEN
         call system_clock(count, rate, maxc) !Artem: add for rand()
         rseed = ( real(mod(count, 4194304)) + 0.5 ) !Artem: add for rand()
     *           / 4194304.0
         if (rseed <= 0.0) rseed = 0.5/4194304.0 !Artem: add for rand()
         ran = rand (rseed) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CCF()
      ELSE
         rseed = ( real(mod(abs(iRanStart), 4194304)) + 0.5 ) / 4194304.0 !Artem: add for rand()
         ran = rand (rseed) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CBF (iRanStart)

         ENDIF
      Write (nPRT, '(a,i8)') ' Random generator start value ',iRanStart

      nHi = iAskDMu  ('# neutrons per group ?', nHi, 1, 100000)
      nGr = iAskDMu  ('# groups ?', nGr, 2, 10000)

      nScEvMax = iAskDMu (
     * 'Russian roulette(0) or maximum # Collisions(>0)',
     * nScEvMax, 0, Mcol)
      IF (nScEvMax.eq.0) THEN
         WgtCut  = rAskDMu  ('Cut-off weight for Monte-Carlo',
     *                       WgtCut, 1d-6, 1d0)
         nScEvRc = Mcol ! max # scattering events that can be recorded
      ELSE
         nScEvRc = nScEvMax
         ENDIF

      qBetaEqProb = qAskD ('Equal probability for every beta',
     *                     intq(qBetaEqProb))

      CALL Sage (' ')
      CALL Sage (' 2/   Experimental setup')
      CALL GetInst (qIgnToF, DistSD, Erg0, Erg0width, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL GetGeom(Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL Sage (' ')
      CALL Sage (' 3/   Scattering properties')

C  Scattering properties of target components :

      nMat = 0
      DO jMat = 1, MMat
         IF (qHasSectn(jMat)) THEN
            nMat = nMat + 1
            IF (nJList.lt.nMat) THEN
               Fehler = 'not enough ideal scattering laws given'
               RETURN
               ENDIF
             ! scattering lengths
            CALL GetMProp (jMat, NamMat(jMat),
     *                     SiMatSC(jMat), SiMatAB(jMat), Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ! scattering law :
            CALL GetIE (JList(nMat), MAC, nAI, AI, MBC, nBI, BI,
     *                  SAB(1,1,jMat), Sint, UnitErg, Fehler)
            IF (Fehler.ne.'&ff') RETURN

            ! the following checks should better be done with general IDA fct's
            IF (nMat.eq.1) THEN
               ! set internal (a,b) grid :
               nAC = nAI
               DO I = 1, nAC
                  AC(I) = AI(I)
                  ENDDO
               nBC = nBI
               DO I = 1, nBC
                  BC(I) = BI(I)
                  ENDDO
               IEfirst = jMat

            ELSE
               ! The grid must be identical for all Target sections :
               DO I=1,nAC
                  IF (AI(I).NE.AC(I)) THEN
                     Fehler='different materials have different A grid'
                     RETURN
                     ENDIF
                  ENDDO
               DO I=1,nBC
                  IF (BI(I).NE.BC(I)) THEN
                     Fehler='different materials have different B grid'
                     RETURN
                     ENDIF
                  ENDDO
               ENDIF

            IF (4*Erg0/TempE.gt.AC(nAC)) THEN
               Print *, ' ratio ', 4*Erg0/TempE/AC(nAC)
               Fehler = 'Erg0 too high or given Q values too small'
               RETURN
               ENDIF
            ENDIF
         ENDDO

C  Set energy grid :

      ! energy transfer limits :
      wMax = BC(nBC)*TempE * (1-1d-8)
      wMin = BC(1)  *TempE * (1+1d-8)

      Write (nPRT,'(a,i3,3g13.6)')
     * 'nBC, BC*T(1), BC*T(nBC), Erg0width',
     *  nBC, BC(1)*TempE, BC(nBC)*TempE, Erg0width
      IF (Erg0width.gt.0.5*dmin1(wMax,-wMin)) THEN
         Print *, ' nBC, BC(1), BC(nBC), T, Erg0width',
     *  nBC, BC(1), BC(nBC), TempE, Erg0width
         Fehler =
     * '(spread in incident energy) > 0.5 * (maximum energy transfer)'
         RETURN
         ENDIF

      ErgMax = Erg0 + wMax
      ErgMin = dmax1 (Erg0 + wMin, Erg0*1.d-2)

c      IF (Erg0+wMin.lt.ErgMin) THEN ! add one point to user-defined grid
c         EC(1) = ErgMin
c         nEC = 1
c      ELSE
         nEC = 0
c         ENDIF
      DO iBC = 1, nBC
         E = Erg0+BC(iBC)*TempE
         IF (ErgMin.lt.E .and. E.lt.ErgMax) THEN
            nEC = nEC + 1
            IF (nEC.gt.MEC) THEN
               Fehler = 'MEC exceeded'
               RETURN
               ENDIF
            EC(nEC) = E
            ENDIF
         ENDDO
c      IF (EC(nEC).lt.ErgMax) THEN ! add one point to user-defined grid
c         nEC = nEC + 1
c         EC(nEC) = ErgMax
c         ENDIF
      IF (nEC.lt.3) THEN
         Fehler = '(almost) no points in energy grid'
         RETURN
         ENDIF
      iE0 = Index1 (EC, nEC, Erg0)
      Write (nPRT, '(a,2i4,5g13.6)')
     * ' internal energy grid: iE0, nEC, Erg0, '//
     * 'wMin, wMax, EC(1), EC(n) ',
     * iE0, nEC, Erg0, wMin, wMax, EC(1), EC(nEC)

C  Normalize to which cross section ?
      IF (qHasSectn(1) .and. qHasSectn(2)) THEN
         CALL Sage (' 3/-  There are several Target sections ...')
         iRefMat = iAskMu ('Normalize output to material ?', 1, 2) ! MMat)
      ELSE
         IF (qHasSectn(1)) THEN
            iRefMat = 1
         ELSE
            iRefMat = 2
            ENDIF
         ENDIF

C  Set angles :
      CALL Sage (' ')
      CALL Sage (' 4/   Output')
      CALL Sage (' 4/1  Enter detector angles (in degrees)')

      iStyle = iAskDMu ('Enter(0), equidist(1)', iStyle, 0, 1)
      IF     (iStyle.eq.0) THEN
         CALL rAskArray (TC, MTC, nTC)
      ELSEIF (iStyle.eq.1) THEN
         nTC = iAskDMu ('Number of detectors', nTC, 1, MTC)
         DO iA = 1, nTC
            TC(ia) = (iA-.5)*180 / nTC
            ENDDO
         ENDIF

      ! set detector direction cosines :
      DO j = 1, nTC
         DirDet(1,j) =  dsind(TC(j))
         DirDet(2,j) =  dcosd(TC(j))
         DirDet(3,j) =  0
         ENDDO

C  Read in and set time channels and channel widths (was SUBROUTINE GetErgOut) :
      ! Set TimOut :
      Print *, '   set time grid'
      DO iE = 1,nEC
         TimOut(iE) = DistSD * sqrt (ConvVof/EC(iE))
         ENDDO

      ! Set channel width :
      Print *, '   set channel width'
      DO iE = 2, nEC-1
         DEC(iE) = ( EC(iE+1) - EC(iE-1) ) / 2
         IF (DEC(iE).le.0) THEN
            Fehler = 'GetEC / Energy must be increasing'
            RETURN
            ENDIF
         ENDDO
      DEC(1)   = 0 ! dmin1 (DEC(2), dabs(EC(1)-EoutMin))
      DEC(nEC) = 0 ! DEC(nEC-1)
      DO iE = 1,nEC
         DTimOut(iE) = DistSD * sqrt (ConvVof/(EC(iE)-DEC(iE)/2)) -
     *                DistSD * sqrt (ConvVof/(EC(iE)+DEC(iE)/2))
         ENDDO

      ! Elastic channel width :
      Print *, '   determine elastic channel'
      t0 = sqrt (ConvVof/Erg0) * DistSD
      Print '(/a38,g10.3)', ' Elastic time-of-flight (usec) = ', t0

      ! find the channel number for elastic scattering :
      KTim0 = IndexTim (MEC, THIGH, TimOut, DTimOut, EC, DEC, nEC, t0)
      IF (KTim0.eq.0) THEN
         Print *, ' WARNUNG/ No elastic channel found'
         KEL = 1
         DTim0 = 1.
      ELSE
         DTim0 = DEC(KTim0)
         Print '(a,i7,a,g10.3)', ' Elastic channel number = ', KTim0,
     *        ' with width (usec) = ', DTim0
         ENDIF

      ! Print out :
      IF (jField.ge.1) THEN
         Write (nPRT,'(/1x,a)') 'Time [usec] from sample to detector :'
         Write (nPRT,'(1x,12e10.3)')  (TimOut(iE), iE=1, nEC) !Artem '(1x,12e10.3))' -> '(1x,12e10.3)'
         ENDIF

      IF (jField.ge.2) THEN
         Write (nPRT,'(/1x,a)') 'Channel widths [usec] :'
         Write (nPRT,'(1x,12e10.3)')  (DTimOut(iE), iE=1, nEC) !Artem '(1x,12e10.3))' -> '(1x,12e10.3)'

         Write (nPRT,'(/1x,a)') 'Corresponding enErgIEs [meV] :'
         Write (nPRT,'(1x,12e10.3)')  (EC(iE), iE=1, nEC) !Artem '(1x,12e10.3))' -> '(1x,12e10.3)'

         Write (nPRT,'(/1x,a)') 'Channel widths [meV] :'
         Write (nPRT,'(1x,12e10.3)')  (DEC(iE), iE=1, nEC) !Artem '(1x,12e10.3))' -> '(1x,12e10.3)'

         ENDIF

      CALL Sage (' 4/3  Scoring and files')
      M2min  = 1 !iAskDMu ('Lowest scattering order for Scoring', 1, 1, nScEvRc)
      qSymOut = qAskD('Symmetrized S~(2th,w) on output', intq(qSymOut))

C  --------------------------------------------------------------------
C     MS /  Prepare Arrays :
C  --------------------------------------------------------------------

C  Calculate auxiliary scattering functions (was SUBROUTINE PrepIE) :
      DO jMat = 1, MMat
         IF (qHasSectn(jMat)) THEN

            ! the alpha-integrated scattering function FS :
            DO iB = 1, nBC
               FS(1,iB,jMat) = 0
               DO iA = 2, nAC
                  FS(iA,iB,jMat) = FS(iA-1,iB,jMat) +
     * (AC(iA)-AC(iA-1)) * (SAB(iA,iB,jMat) + SAB(iA-1,iB,jMat)) / 2
                  IF (FS(iA,iB,jMat).lt.FS(iA-1,iB,jMat)) THEN
                     Fehler = 'PROGRAM ERROR/ FS non monotonic in A'
                     RETURN
                     ENDIF
                  ENDDO
               ENDDO

            ! the beta-integrated functions PMAR, SiErgIE :
            DO iE = 1, nEC ! loop over neutron energy epsilon_0

               E0 = EC(iE) / TempE

               iB1 = irPos (BC, nBC, -E0, 'l')
               aLow2 = Alpha_BEC (BC(iB1), E0, +1.d0)
               aHig2 = Alpha_BEC (BC(iB1), E0, -1.d0)

               PMAR (1,iE,jMat) = 0 ! integral in beta and alpha

               DO iB = iB1+1, nBC ! loop over scattering energy beta

                  Pint = 0 ! integral in alpha

                  ! integration limits :
                  aLow1 = aLow2 ! from last value of B
                  aHig1 = aHig2
                  aLow2 = Alpha_BEC (BC(iB), E0, +1.d0)
                  aHig2 = Alpha_BEC (BC(iB), E0, -1.d0)

                  DO iA = 2, nAC
                     Pint = Pint +
     *                    Alpha_Base (AC(iA-1), AC(iA),
     *                    aLow1, aLow2, aHig1, aHig2, Fehler) *
     *                    (SAB(iA,iB,jMat) + SAB(iA-1,iB,jMat)) / 2
                     IF (Fehler.ne.'&ff') THEN
                        Print *, 'args :', iE, E0, iB, BC(iB)*TempE,
     * iA, AC(iA-1), AC(iA), aLow1, aLow2, aHig1, aHig2
                        RETURN
                        ENDIF
                     ENDDO ! iA

                  PMAR (iB,iE,jMat) = PMAR(iB-1,iE,jMat)+
     *                                Pint*(BC(iB)-BC(iB-1))

                  ENDDO ! iB

               ! Inelastic scattering probability :
               SiErgIE(iE,jMat) = PMAR(nBC,iE,jMat) *
     *                            SiMatSC(jMat) / E0 / 4

               ENDDO ! iE

            IF (jField.ge.3) THEN
               CALL SaveMatrix (JOutF(1), fil(1:lf)//
     *                      'sc', 'SC(a,b) < '//tit,
     *                      'b', ' ', 'a', ' ', 'SC(a,b)', ' ',
     *                      MBC, MAC, nBC, nAC, BC, AC, SAB(1,1,jMat),
     *                      Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL SaveMatrix (JOutF(2), fil(1:lf)//
     *                      'fs', 'FS(a,b) < '//tit,
     *                      'b', ' ', 'a', ' ', 'FS(a,b)', ' ',
     *                      MBC, MAC, nBC, nAC, BC, AC, FS(1,1,jMat),
     *                      Fehler)
               IF (Fehler.ne.'&ff') RETURN
               CALL SaveMatrix (JOutF(3), fil(1:lf)//
     *             'pm', 'PM(b,E) < '//tit,
     *             'b', ' ', 'E', '?', 'PM(b,E) not normalised', ' ',
     *             MBC, MEC, nBC, nEC, BC, EC, PMAR(1,1,jMat),
     *             Fehler)
               IF (Fehler.ne.'&ff') RETURN
               ENDIF ! save auxiliary fields

            DO iE = 1, nEC ! loop over neutron energy epsilon_0
               ! Normalize scattering probability :
               DO iB = 1, nBC
                  PMAR(iB,iE,jMat) = PMAR(iB,iE,jMat) /
     *                               PMAR(nBC,iE,jMat)
                  ENDDO
               ENDDO ! iE

            ENDIF ! has material
         ENDDO ! jMat

C  Compute total cross sections (for use in ScorIE) :
      DO jMat = 1, MMat
         IF (jField.ge.2) THEN
            Write (nPRT, '(/a,i2/2(a,f10.4))')
     *        ' Total cross section for material', jMat,
     *        '    calculated with sigma_sc=', SiMatSC(jMat),
     *        ' and sigma_abs(2200km/sec)=', SiMatAB(jMat)
            Write (nPRT, '(a3,2x,6(a8,2x))')
     *        'K', 'E', 'S_IE', 'S_AB', 'S_TT'
            ENDIF

         IF (qHasSectn(jMat)) THEN
            DO iE = 1, nEC
               ! Cross sections are function of neutron's energy :
               E = EC(iE)
               ! total = scattering + absorption cross section :
               SiErgTT(iE,jMat) = SiErgIE(iE,jMat)
     *                          + SiMatAB(jMat)*dsqrt(25.30d0/E)

               IF (jField.ge.2) Write (nPRT,'(i3,2x,4(g11.6,2x))')
     *  iE, E, SiErgIE(iE,jMat), SiMatAB(jMat) *
     *  dsqrt(25.30d0/E), SiErgTT(iE,jMat)
               ENDDO ! iE
            ENDIF
         ENDDO ! jMat

C  Initialize counters :
      NnoHit  = 0   ! bad impact -> target not hit
      nLossBL = 0   ! # neutrons lost out of energy range
      wLossBL = 0   ! their weight
      nLossBH = 0   ! # neutrons lost out of energy range
      wLossBH = 0   ! their weight
      nLossA  = 0   ! # neutrons lost out of q range
      wLossA  = 0   ! ..
      Transm  = 0   ! sample transmission, averaged over all incident neutrons
      TarThick= 0   ! average thickness of target (cm)

      DO IS = 1,Mscor
         DO iT = 1, nTC
            DO iE = 1, nEC
               RSA (iE,iT,IS) = 0
               RSQ (iE,iT,IS) = 0
               ENDDO
            ENDDO
         ENDDO
      DO iE=1,nEC
         NUnscored(iE) = 0
         ENDDO
      DO iC = 1, Mcol
         WgtCol(iC) = 0
         ENDDO

C  --------------------------------------------------------------------
C     MS /  The Monte-Carlo Loop :
C  --------------------------------------------------------------------
         ! Outer main loop : groups of neutrons treated together
         !                   for the purpose of calculating standard deviations

      Print *, ' MS enters the main loop'
      Write (nPRT,'(//76(1h=)/a/76(1h=)//)') 'MONTE CARLO LOG'

      DO iGr = 1,nGr

C  Reset group response :
         DO IS = 1,Mscor
            DO iT = 1, nTC
               DO iE = 1, nEC
                  RSG (iE,iT,IS) = 0
                  ENDDO
               ENDDO
            ENDDO

C  Inner main loop : single neutron histories :
         DO iHi = 1,nHi

C  Initialize the neutron :
            M1 = 1  ! binary event indicator
            M2 = 0  ! # Collisions
            Wgt= 1. ! weight

 2100       CALL RanImp (DisImp, DisThru, qHit)
            IF (.not.qHit) THEN
               NnoHit = NnoHit + 1
               IF (NnoHit.gt.nGr*nHi*10) THEN
                  Fehler ='ini neutron / seldom if ever hitting target'
                  RETURN
                  ENDIF
               GOTO 2100
               ENDIF
            CALL ImpactTOF (qIgnToF, Erg0,
     *                      TimN, ErgN, Erg0width, EC(1), EC(nEC))

            IF (.not.qIgnToF) TimN = TimN + DisImp * sqrt(ErgN/ConvVof)

C  One-scattering-event loop :
 2200       CONTINUE

            IF (jTrace.ge.3) THEN
               Write (nPRT, '(a,f7.4,x,f8.5,x,3(f5.3,1x))')
     *            '>>  fly with w, E, D                  ',
     *                       Wgt, ErgN, DirN(1), DirN(2), DirN(3)
               ENDIF

            ! distances SSAM, SCAN to travel through @s and @c :
            CALL CalcDist (DirN, DistSect)

C  Scattering SigSC and total SigTT cross section :
            DO jMat = 1, MMat
               IF (qHasSectn(jMat)) THEN
                  IF (ErgN.gt.EC(nEC)) THEN
                     Fehler = 'PROGRAM ERROR/  cannot calculate SigSC'
                     RETURN
                     ENDIF
                  iE = Index1 (EC, nEC, ErgN)
                  SigSC(jMat) = SiErgIE(iE,jMat)
                  SigTT(jMat) = SigSC(jMat) +
     *                           SiMatAB(jMat)*dsqrt(25.30d0/ErgN) ! total
                  IF (SigTT(jMat).le.0) THEN
                     Fehler = 'PROGRAM ERROR/ sigma = 0'
                     RETURN
                     ENDIF
               ELSE
                  SigSC(jMat) = 0
                  SigTT(jMat) = 0
                  ENDIF
               ENDDO

C  Modify Wgt according to the chance the neutron has to leave the sample :
            survive = dexp(-DistSect(1)*SigTT(1)
     *                     -DistSect(2)*SigTT(2)
     *                     -DistSect(3)*SigTT(3))
            Wgt = Wgt * (1-survive)

            IF (jTrace.ge.3) Write (nPRT, '(a,3f9.4,a,f7.4)')
     * '>> distances thru S,C,A : ', DistSect, ' => reduced w = ', Wgt
            IF (jTrace.ge.4) Write (nPRT, '(a,3f9.4)')
     *         '& SigTT(SCA) : ', SigTT

            IF (M2.eq.0) THEN
               Transm   = Transm   + survive
               TarThick = TarThick + DistSect(iRefMat)
               ENDIF

C  Choose a collision point :
            CALL RanCol (SigTT, DistCol, InMat)
            IF (.not.qHasSectn(InMat)) THEN
               Print *, InMat
               Fehler =
     * 'RanCol finds collision in invalid sample section'
               RETURN
               ENDIF

            ! time needed to reach collsion point :
            IF (.not.qIgnToF) TimN = TimN + DistCol*sqrt(ErgN/ConvVof)

            ! force scattering (reduce weight instead of allowing absorption) :
            Wgt=Wgt*SigSC(InMat)/SigTT(InMat)
            IF (Wgt.le.0) THEN
               Fehler =
     * 'PROGRAM ERROR/ SigSC<=0 (sample section absorbs all)'
               RETURN
               ENDIF

            ! Statistics :
            M2 = M2 + 1
            IF (M2.le.Mcol) WgtCol(M2) = WgtCol(M2) + Wgt

            IF (jTrace.ge.3) Write (nPRT, '(a,i2,a,f9.4,a,i1,a,f7.4)')
     * '>> collision no. ', M2, ' after dist ', DistCol, ' in mat ',
     * InMat, ' => reduced w = ', Wgt

C  Perform the appropriate Scoring :
            ! Determine class IScor of event to be scored :
            IF (M2.lt.2) THEN
               IScor = InMat
               InOneMat = InMat
            ELSEIF (InMat.eq.InOneMat) THEN
               IScor = MMat + InMat
            ELSE
               IScor=2*MMat+1
               InOneMat = 0
               ENDIF

            IF (M2.lt.M2min) GOTO 2600 ! skip low-order scoring to save time

C  Perform inelastic scoring :
            PE   = dsqrt(ErgN)
            ! relative cross section = bound / scattering cross section :
            FacY = Wgt * SiMatSC(InMat)/SigSC(InMat) / TempE / PE

            iA = 1
            iB = 1
            ! Outer loop over detector angles :
            DO iT=1,nTC

               CALL CalcDist (DirDet(1,iT), Travel)

               CO = DirN(1)*DirDet(1,iT) + DirN(2)*DirDet(2,iT) +
     *              DirN(3)*DirDet(3,iT)
               IF (jTrace.ge.5) Write (nPRT,
     *          '(a,i3,a,f6.2,a,f6.3,a,3f8.4)')
     *          '  score to det.', iT, ' at ang ', TC(iT), ' CO ', CO,
     *          ' dist(mat) ', Travel

               ! Distance to detector :
               IF (.not.qIgnToF) THEN
                  XL = DistSD - ( PosN(1)*DirDet(1,iT) +
     *                 PosN(2)*DirDet(2,iT) + PosN(3)*DirDet(3,iT) )
                  VL = dsqrt(ConvVof) * XL
                  ENDIF

               ! Inner loop over time channels :

                  ! NOTE : this is the innermost loop of the program,
                  !        here some optimization has been done.

               DO iE = 1, nEC
                  IF (qIgnToF) THEN
                     PSR = dsqrt(EC(iE))
                  ELSE
                     ! Time of arrival at detector :
C                     TARV = G05DAF (TimOut(iE)-DTimOut(iE)/2,
C     *                              TimOut(iE)+DTimOut(iE)/2) !Artem: Replace with rand function from slatec/fnlib/rand.f
                      ran = rand(0.0)
                      TARV = TimOut(iE)-DTimOut(iE)/2 +
     * (TimOut(iE)+DTimOut(iE)/2 - TimOut(iE)-DTimOut(iE)/2) *
     * dble(ran)
                     PSR =  VL / (TARV-TimN)
                     ENDIF
                  ! Scattered energy :
                  Erg1 = PSR**2
                  w   = Erg1 - ErgN
                  IF (w.le.wMin .or. w.ge.wMax .or.
     *                Erg1.lt.EC(1) .or. Erg1.gt.EC(nEC)) THEN
                     NUnscored(iE)=NUnscored(iE)+1
                     GOTO 259
                     ENDIF
                  ! calculate A,B :
                  A = (Erg1+ErgN-2*PE*PSR*CO)/TempE
                  B = w / TempE
                  IF (A.gt.AC(nAC)) THEN
                     NUnscored(iE)=NUnscored(iE)+1
                     GOTO 259
                     ENDIF
                  ! interpolate S(a,b) {war jahrelang falsch 4aug99}
                     iA = Index2 (AC, nAC, A, iA)
                     iB = Index2 (BC, nBC, B, iB)
                     br = (B-BC(iB-1)) /  (BC(iB)-BC(iB-1))
                     bl = (BC(iB  )-B) /  (BC(iB)-BC(iB-1))
                     a2 = (A-AC(iA-1)) /  (AC(iA)-AC(iA-1))
                     a1 = (AC(iA  )-A) /  (AC(iA)-AC(iA-1))
                     s1 = bl * SAB(iA-1,iB-1,InMat) +
     *                    br * SAB(iA-1,iB,InMat)
                     s2 = bl * SAB(iA  ,iB-1,InMat) +
     *                    br * SAB(iA  ,iB,InMat)
                     s12= a1 * s1 + a2 * s2
                  ! find transmission
                  exneg = - Travel(1)*SiErgTT(iE,1)
     *                    - Travel(2)*SiErgTT(iE,2)
     *                    - Travel(3)*SiErgTT(iE,3)
                  ! Contribution to inelastic Score :
                  RSG(iE,iT,IScor) = RSG(iE,iT,IScor) +
     *                 FacY * s12 * PSR * dexp(exneg)
c               IF (M2.eq.2 .and. iE.eq.iE0)
c     *              write (nPRT, '(a,4g15.8)'), 'HEUTE: ',
c     *              TC(iT), dacosd(aux), dacosd(CO),
c     *              FacY * s12 * PSR * dexp(exneg)
 259              CONTINUE

                  ENDDO ! EnErgIEs
               ENDDO ! Angles

 2600       CONTINUE

C  Determine fate of neutron :
            IF (nScEvMax.gt.0) THEN
               ! # Collisions fixed :
               IF (M2.ge.nScEvMax) GOTO 2790 ! exit

            ELSE
               ! Compare Wgt with cutoff weight :
               DO WHILE (Wgt.le.WgtCut) ! Important correction 13jul94
                  ! .. play Russian roulette :
                  IF (dble(rand(0.0)).gt.0.5) GOTO 2790 ! exitus  !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
                  Wgt = 2*Wgt
                  IF (jTrace.ge.2) Write (nPRT, '(a,f7.4)')
     * '>>>>   survived Russian roulette, doubling weight to ', Wgt
                  ENDDO
               ENDIF

C  Select inelastic scattering event :
            IF (ErgN.gt.EC(nEC)) THEN
               Fehler =
     * 'PROGRAM ERROR/ ErgN>EC(nEC) on input to new collision'
               RETURN
               ENDIF
            iE = Index1 (EC, nEC, ErgN)
            IF (iE.lt.2 .or. iE.gt.nEC) THEN
               Fehler = 'PROGRAM ERROR/ iE oor'
               RETURN
               ENDIF
            fE = (ErgN-EC(iE-1)) / (EC(iE)-EC(iE-1))

            ! Select beta :
            IF (qBetaEqProb) THEN ! this option since 27jan00

               Fehler = 'beta-equiprobable does not yet work'
               RETURN

               ! beta_min(E) nur vorlaeufig hier:
               b1 = - ErgN/TempE
               DO iB1 = 1, nBC
                  IF (BC(iB1).gt.b1) GOTO 2619
                  ENDDO
               Fehler = 'no channels with Enew>0'
               RETURN
 2619          CONTINUE
               ! beta equally probable within user-defined mesh
               prob = iB1 + (nBC-iB1) * dble(rand(0.0)) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
               iB = idint(prob) + 1 ! (will also be used for alpha-selection)
               fB = prob - (iB-1)
               beta = BC(iB-1) + fB*(BC(iB)-BC(iB-1))

               ! now calculate probability for this event
               wr1 = PMAR(iB,iE-1,InMat) - PMAR(iB-1,iE-1,InMat)
               wr2 = PMAR(iB,iE,  InMat) - PMAR(iB-1,iE  ,InMat)
               wr  = wr1 + fE*(wr2-wr1)
               IF (wr.le.0) THEN
                  Fehler = 'PROGRAM ERROR/ SEl -> weight=0'
                  RETURN
                  ENDIF

               IF (jTrace.ge.6) Write (nPRT, *)  ' SEL ',
     *                ErgN, prob, beta, wr, nBC-iB1

               Wgt = Wgt * wr * (nBC-iB1)

            ELSE ! conventional choice of beta according to probability
               Bran = dble(rand(0.0)) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
               ir1 = Index1 (PMAR(1,iE-1,InMat), nBC, Bran)
               ir2 = Index1 (PMAR(1,iE,  InMat), nBC, Bran) ! Index2 ??

               fr1 = (Bran-PMAR(ir1-1,iE-1,InMat))
     *              / (PMAR(ir1,iE-1,InMat)-PMAR(ir1-1,iE-1,InMat))
               fr2 = (Bran-PMAR(ir2-1,iE,  InMat))
     *              / (PMAR(ir2,iE,  InMat)-PMAR(ir2-1,iE,  InMat))

               b1  = BC(ir1-1) + fr1*(BC(ir1)-BC(ir1-1))
               b2  = BC(ir2-1) + fr2*(BC(ir2)-BC(ir2-1))
               beta = b1 + fE*(b2-b1)

               IF (jTrace.ge.6) THEN
                  Write (nPRT, '(a,g12.7,i3,3g11.5)') ' current erg: ',
     * ErgN, iE, EC(iE-1), EC(iE), fE
                  Write (nPRT, '(a,g11.6,i4,3g11.6)') ' ran(iE-1): ',
     * Bran, ir1, PMAR(ir1-1,iE-1,InMat), PMAR(ir1,iE-1,InMat), fr1
                  Write (nPRT, '(a,g11.6,i4,3g11.6)') ' ran(iE)  : ',
     * Bran, ir2, PMAR(ir2-1,iE,  InMat), PMAR(ir2,iE,  InMat), fr2
                  Write (nPRT, *) EC(iE),   ':',
     * PMAR(ir2-1,iE-1,InMat), PMAR(ir2,iE,InMat)
                  Write (nPRT, '(a,4g12.6)') ' B(iE-1): ',
     * BC(ir1-1), BC(ir1), fr1, b1
                  Write (nPRT, '(a,4g12.6)') ' B(iE)  : ',
     * BC(ir2-1), BC(ir2), fr2, b2
                  Write (nPRT, '(a,4g12.6)') ' --> ',
     * fE, beta, beta*TempE, ErgN
                  ENDIF

               iB = ir1 ! for later use (the approximations become cruder)
               ENDIF ! qBetaEqProb

            ! beta implies new neutron energy :
            Eps    = ErgN / TempE
            Epsnew = Eps + beta
            ErgN   = Epsnew * TempE
c            IF     (ErgN.le.0) THEN
c               Fehler = 'PROGRAM ERROR/ scattered energy is nonpositive'
c               RETURN
c               ENDIF
            IF     (ErgN.le.EC(1)) THEN
               wLossBL = wLossBL + Wgt
               nLossBL = nLossBL + 1    ! neutron is lost out of energy range
               GOTO 2790
            ELSEIF (ErgN.ge.EC(nEC)) THEN
               Write (nPRT, '(a,i2,4g12.6)') 'loss E>EC(n): ',
     *              M2, Eps*TempE, ErgN, beta*TempE, EC(nEC)
               wLossBH = wLossBH + Wgt
               nLossBH = nLossBH + 1    ! neutron is lost out of energy range
               GOTO 2790
               ENDIF

            ! Prepare selection of alpha (accessible range depends on beta) :
            aHig = Alpha_BEC (beta, Eps, -1.d0)
            aLow = Alpha_BEC (beta, Eps, +1.d0)

            IF (aHig.gt.AC(nAC)) THEN ! not covered by given scattering law
               wLossA = wLossA + Wgt
               nLossA = nLossA + 1
               GOTO 2790
               ENDIF
            iHig = Index1 (AC, nAC, aHig)
            IF (aLow.lt.0) THEN
               Fehler = 'PROGRAM ERROR/ alpha selection/ aLow < 0'
               RETURN
            ELSEIF (aLow.eq.0) THEN
               iLow = 1
            ELSE
               iLow = Index1 (AC, nAC, aLow) - 1
               ENDIF
            IF (iLow.lt.1) THEN
               Fehler = 'PROGRAM ERROR/ alpha selection/ iLow < 1'
               RETURN
               ENDIF

            CALL LinIntPol1 (aLow, AC(iLow), AC(iLow+1),
     * FS(iLow,iB,InMat), FS(iLow+1,iB,InMat), Fmin)
            CALL LinIntPol1 (aHig, AC(iHig-1), AC(iHig),
     * FS(iHig-1,iB,InMat), FS(iHig,iB,InMat), Fmax)

            ! Now select alpha :
C            Fran = G05DAF (Fmin, Fmax) !Artem: Replace with rand function from slatec/fnlib/rand.f
            ran = rand(0.0)
            Fran = Fmin + (Fmax - Fmin) * dble(ran)
            IF (Fran.gt.FS(nAC,iB,InMat)) THEN
               Fehler = 'PROGRAM ERROR/ vor Index 1 abgefangen/ FS'
               RETURN
               ENDIF

            iA   = Index1 (FS(1,iB,InMat), nAC, Fran)
            CALL LinIntPol1 (Fran, FS(iA-1,iB,InMat), FS(iA,iB,InMat),
     *                       AC(iA-1), AC(iA), alpha)
            IF (alpha.gt.aHig .or. alpha.lt.aLow) THEN
               Print *, 'a.. ', alpha, aHig, aLow, iA
               Fehler = 'PROGRAM ERROR/ alpha o.o.r.'
               RETURN
               ENDIF

            ! alpha implies scattering cosine CO :
            CO=(Epsnew+Eps-alpha) / 2 / dsqrt(Epsnew*Eps)
            IF (dabs(CO).gt.1) THEN
               Fehler = 'PROGRAM ERROR/ Illegal cattering cosinus'
               RETURN
               ENDIF

            ! choose azimutal angle at random :
            CALL RandAng (aziS,aziC) ! sine and cosine
            SI = dsqrt(1-CO*CO)      ! scattering sine
            CALL NewDir (DirN, CO, SI, aziC, aziS)
            IF (jTrace.ge.4) Write (nPRT, '(a,2f8.5,a,3f8.5)')
     *            '> scattered by ', CO, aziC, ' in new dir ', DirN

            aux = CO

            GOTO 2200
C  End of one-scattering-event loop.

 2790       CONTINUE
            ENDDO ! iHi
C  End of inner main loop (neutron history).

C  Increment the response arrays :
         DO IS = 1,Mscor
            DO iT=1,nTC
               DO iE = 1, nEC
                  RSA (iE,iT,IS) = RSA (iE,iT,IS) + RSG(iE,iT,IS)
                  RSQ (iE,iT,IS) = RSQ (iE,iT,IS) + RSG(iE,iT,IS)**2
                  ENDDO
               ENDDO
            ENDDO ! iT

         CALL Counter (iGr, 1, nGr, 'end of group')

         ENDDO ! iGr
C  End of outer main loop (neutron group).

C  --------------------------------------------------------------------
C     MS /  Normalization and Output :
C  --------------------------------------------------------------------

      Print *, ' ... end of the main loop, now printing results ...'
      Print *

C  History statistics :

      Print '(a)', 'Neutron history statistics :'

      nImp = nGr * nHi

      Print '(a,i9)',
     *     ' neutrons leaving the choppers     ', nImp + NnoHit !Artem '(a,i9))' -> '(a,i9)'
      Print '(a,i9,a,f5.2,a)',
     *     ' of which missing the sample       ', NnoHit,
     *     ', i.e. ', 100.d0*NnoHit/(nImp+NnoHit), ' %' !Artem '(a,i9,a,f5.2,a))' -> '(a,i9,a,f5.2,a)'
      Print '(a,i9)',
     *     ' => hitting the sample                  ', nImp + NnoHit !Artem '(a,i9))'' -> '(a,i9)''
      Print '(a,i9,a,f5.2,a,f6.4)',
     *   ' of which lost below energy range ', nLossBL,
     *   ', i.e. ', 100.d0*nLossBL/nImp, ' % with <w> ',
     *    dnquot (wLossBL, nLossBL) !Artem '(a,i9,a,f5.2,a,f6.4))'' -> '(a,i9,a,f5.2,a,f6.4)''
      Print '(a,i9,a,f5.2,a,f6.4)',
     *   ' of which lost above energy range ', nLossBH,
     *   ', i.e. ', 100.d0*nLossBH/nImp, ' % with <w> ',
     *    dnquot (wLossBH, nLossBH) !Artem '(a,i9,a,f5.2,a,f6.4))' -> '(a,i9,a,f5.2,a,f6.4)'
      Print '(a,i9,a,f5.2,a,f6.4)',
     *   '                and out of q range ', nLossA,
     *   ', i.e. ', 100.d0*nLossA/nImp, ' % with <w> ',
     *    dnquot (wLossA, nLossA) !Artem '(a,i9,a,f5.2,a,f6.4))'' -> '(a,i9,a,f5.2,a,f6.4)''
      nTrc = nImp - nLossBL -nLossBH - nLossA
      Print '(a,i9)',
     *   ' => traced until death             ', nTrc !Artem '(a,i9))' -> '(a,i9)'

      DO iE= 1, nEC
         IF (NUnscored(iE).ne.0) THEN
            Print '(a,i4,a,i9)', '  in channel ', iE,
     *           ' we lost ', NUnscored(iE), ' neutrons'
            ENDIF
         ENDDO

      Print '(a)', 'Neutron weight statistics:'
      DO iC = 1, Mcol
         Print '(a,i2,a,f10.6)', '   coll. no. ', iC, ': <wgt> = ',
     *        WgtCol(iC) / nImp
         ENDDO

C  Transmission :

      Transm    = Transm  / nImp              ! average transmission
      TarThick = TarThick / nImp              ! average thickness of target

      Print '(a,f10.5)', 'Transmission          ', Transm
      Print '(a,f10.5)', 'Target thickness (cm) ', TarThick

C  Scattering law -> output files :

      Print *
      Print *, ' ... now save scattering laws ...'
      Print *

      IF (nMat.le.1) THEN
         nScor = 4
      ELSE
         nScor = 7
         ENDIF
      IF (nScor.gt.Mscor) THEN
         Fehler = 'PROGRAM ERROR/ Mscor oor'
         RETURN
         ENDIF

      ! open internal files (IDA-online-memory) :
      CALL OlfCreate (JOut(1), Kout, fil(1:lf)//'i', 'sim (ideal) '//
     *                tit, Fehler)
      CALL OlfCreate (JOut(2), Kout, fil(1:lf)//'s', 'sim (singl) '//
     *                tit, Fehler)
      CALL OlfCreate (JOut(3), Kout, fil(1:lf)//'m', 'sim (multi) '//
     *                tit, Fehler)
      CALL OlfCreate (JOut(4), Kout, fil(1:lf)//'t', 'sim (total) '//
     *                tit, Fehler)
      IF (nScor.gt.4) THEN
      CALL OlfCreate (JOut(5), Kout, fil(1:lf)//'e', 'sim (sampl) '//
     *                tit, Fehler)
      CALL OlfCreate (JOut(6), Kout, fil(1:lf)//'c', 'sim (contr) '//
     *                tit, Fehler)
      CALL OlfCreate (JOut(7), Kout, fil(1:lf)//'a', 'sim (absrp) '//
     *                tit, Fehler)
         ENDIF

      ! Save coordinate names and some documentation :

      UnitSqw = UnitErg(1:lenU(UnitErg))//'-1'
      CALL UnitConv ('meV', UnitErg, FacX, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO ij = 1, nScor

         CALL OlfCnuP (JOut(ij), 'x', 'w/2pi', UnitErg, Fehler)
         IF (qSymOut) THEN
            CALL OlfCnuP (JOut(ij), 'y', 'S~(2th,w)', UnitSqw, Fehler)
            CALL iOlfP (JOut(ij), '?det-bal-sym',   1, Fehler)
         ELSE
            CALL OlfCnuP (JOut(ij), 'y', 'S(2th,w)', UnitSqw, Fehler)
            CALL iOlfP (JOut(ij), '?det-bal-sym',   0, Fehler)
            ENDIF
         CALL OlfCnuP (JOut(ij), 'z1', '2th', ' ', Fehler)

         CALL rOlfP (JOut(ij), 'E0', 'meV', Erg0, Fehler)
         CALL rOlfP (JOut(ij), 'T', 'K', TempK, Fehler)
         CALL iOlfP (JOut(ij), '@sam-erg-gain', -1, Fehler)

         CALL Compose3 (aus, 'Simulated '//cl6(nGr), '*'//
     *                  cl6(nHi), ' neutrons')
         IF (nScEvMax.gt.0) THEN
            CALL Compose2 (aus, aus, ' up to order '//cl3(nScEvMax))
         ELSE
            CALL NiceNum (WgtCut, hnum, ih)
            CALL Compose2 (aus, aus, ' down to weight '//hnum)
            ENDIF
         CALL OlfComAdd (JOut(ij), ' ', aus, Fehler)

         CALL NiceNum (Transm, hnum, ih1)
         CALL NiceNum (TarThick, hnu2, ih2)
         aus = 'transmission = '//hnum(1:ih1)//
     *         ', eff. thickness = '//hnu2(1:ih2)//'cm'
         CALL OlfComAdd (JOut(ij), ' ', aus, Fehler)

         IF (Fehler.ne.'&ff') RETURN
         ENDDO

C  Prefactors for normalization of response arrays :
      IF (TarThick.le.0) THEN
         Fehler = 'PROGRAM ERROR ? TarThick <= 0'
         RETURN
         ENDIF
      RelRed = 1/TarThick ! see C2,137
      FactA = dnquot (RelRed, nGr * nHi)
      FactQ = dnquot (RelRed**2, (nGr-1)*nHi**2)
      FactR = dnquot (1.d0, nGr)

         ! mean(R) = FactA * sum (R)
         ! var (R) = FactQ *(sum(R**2) - FactR*sum(R)**2)

C  Main loop : scattering angles :

      DO iT = 1, nTC
         SCO = dcosd(TC(iT)) ! the cosine

C  Loop over time channels :
         ! first guess for Index2 :
         iA = 1
         iB = 1

         DO iE = 1, nEC

            Erg1 = EC(iE)
            B    = (Erg1-Erg0)/TempE
            Xout(iE) = (Erg1-Erg0) * FacX ! w>0 = neutron energy gain

            ! Conversion factor (cross section -> S(q,w)) :
            FacY = dsqrt(Erg0/Erg1) / SiMatSC(iRefMat) / FacX

            ! Detailed balance correction (ScorIE produces asymmetric S(q,w)) :
            IF (qSymOut) FacY = FacY * dexp(B/2)

C  Calculate ideal inelastic cross section :
            IF (qHasSectn(1) .or. qHasSectn(2)) THEN
               A=(Erg1+Erg0-2*dsqrt(Erg0*Erg1)*SCO)/TempE
               IF (A.gt.AC(nAC)) THEN ! happens with big Erg0

                  DO jMat = 1, MMat
                     SQW(iE,Mscor+jMat,1) = 0
                     SQW(iE,Mscor+jMat,2) = 0
                     ENDDO
                  GOTO 3267
                  ENDIF
               iA = Index2 (AC, nAC, A, iA)
               iB = Index2 (BC, nBC, B, iB)
               ENDIF

            DO jMat = 1, MMat
               ! inelastic contribution :
               IF (qHasSectn(jMat)) THEN
                  ! interpolate S(a,b)
                  br = (B-BC(iB-1)) /  (BC(iB)-BC(iB-1))
                  bl = (BC(iB  )-B) /  (BC(iB)-BC(iB-1))
                  a2 = (A-AC(iA-1)) /  (AC(iA)-AC(iA-1))
                  a1 = (AC(iA  )-A) /  (AC(iA)-AC(iA-1))
                  s1 = bl * SAB(iA-1,iB-1,jMat) + br*SAB(iA-1,iB,jMat)
                  s2 = bl * SAB(iA  ,iB-1,jMat) + br*SAB(iA  ,iB,jMat)
                  s12= a1 * s1 + a2 * s2
                     ! convert b-1 -> hbarw-1 and reverse trivial prefactors
                 SIE = s12 / TempE / FacX
                  IF (qSymOut) SIE = SIE * dexp(B/2)
               ELSE
                  SIE = 0.
                  ENDIF
               SQW(iE,Mscor+jMat,1) = SIE   ! ideal S(q,w)
               SQW(iE,Mscor+jMat,2) = 0.    ! its variance
               ENDDO

 3267       CONTINUE
                 !  Normalized means RSA, squared standard errors RSQ,
            !  and joined response SQW from all scattering types :
            DO IS = 1,Mscor
               SQW(iE,IS,1) = FacY    * FactA * RSA(iE,iT,IS)
               SQW(iE,IS,2) = FacY**2 * FactQ *
     *              (RSQ(iE,iT,IS) - FactR * RSA(iE,iT,IS)**2)
               SQW(iE,IS,2) = dmax1 (SQW(iE,IS,2), 0.d0) ! precaution
               ENDDO

            !  New grouping of scattering types for output (assuming MMat==3) :
            DO II = 1,2
               SIO(iE,1,II) =  SQW(iE,8,II) + SQW(iE,9,II)      ! ideal
               SIO(iE,2,II) =  SQW(iE,1,II) + SQW(iE,2,II)
     *                                      + SQW(iE,3,II)      ! single
               SIO(iE,3,II) =  SQW(iE,4,II) + SQW(iE,5,II)
     *                       + SQW(iE,6,II) + SQW(iE,7,II)      ! multi
               SIO(iE,4,II) =  SIO(iE,2,II) + SIO(iE,3,II)      ! total
               SIO(iE,5,II) =  SQW(iE,4,II) + SQW(iE,1,II)      ! @s
               SIO(iE,6,II) =  SQW(iE,5,II) + SQW(iE,2,II)      ! @c
               SIO(iE,7,II) =  SQW(iE,6,II) + SQW(iE,3,II)      ! @a
               ENDDO

            DO Iscout = 1,nScor
               SIO(iE,Iscout,2) = dsqrt(SIO(iE,Iscout,2)) ! variance
               ENDDO

            ENDDO ! iE

         ! write spectra to internal storage (IDA-online-memory) :
         DO isco = 1, nScor
            CALL OlfPutSpe (JOut(isco), iT, 1, TC(iT),
     *              nEC, Xout, SIO(1,isco,1), SIO(1,isco,2), Fehler )
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

C  End loop in angle.
         ENDDO ! iT

      DO isco = 1, nScor ! close internal files
         CALL OlfClos (JOut(isco), nTC, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      Print *, ' simulation terminated successfully'

      END ! MS$main

C  ====================================================================
C
C      Module IE (inelastic scattering)
C
C  ====================================================================

C  Contents of this module : GetIE, ScorIE

C  ====================================================================
      SUBROUTINE GetIE (j, MAC, nAC,AC,MBC,nBC,BC,SAB,Sint,UnX,Fehler)
C  ====================================================================
         ! get the external scattering function S~(a,b).

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)
      INCLUDE           'i_dim.f'

      REAL*8        AC(MAC), BC(MBC), SAB(MAC,MBC), Sint(MAC), SinMin,
     *              SintMax
      CHARACTER*(*) UnX, Fehler
      CHARACTER*40  Co, Un, UnSoll
      REAL*8        Xin(MC), Yin(MC), Din(MC), X0(MC)

      COMMON /Convert/   TempE, ConvVof
      COMMON /Trace/     nPRT, jTrace

      IF (Fehler.ne.'&ff') THEN
         Print *, 'ERROR on entry in GetIE'
         RETURN
         ENDIF

C  Open input file :

      CALL OlfCnuG (j, 'x', Co, UnX, Fehler)
      IF (Fehler.ne.'&ff') GOTO 99
      IF     (Co.ne.'w/2pi') THEN
         IF (Co(1:1).ne.'w') THEN
            Fehler = 'x-coordinate not w/2pi but '//Co
            GOTO 99
         ELSE
            CALL Say2 (' assuming x-coordinate '//Co, ' means w/2pi')
            ENDIF
         ENDIF
      CALL UnitConv (UnX, 'meV', FacX, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL OlfCnuG (j, 'y', Co, Un, Fehler)
      IF (Fehler.ne.'&ff') GOTO 99
      IF     (Co.ne.'S~(q,w)') THEN
         Fehler = 'y-coordinate not S~(q,w) but '//Co
         GOTO 99
      ENDIF
      UnSoll = UnX(1:lenU(UnX))//'-1'
      IF (Un.ne.UnSoll) THEN
         CALL Compose2 (Fehler, 'y-unit is not '//UnSoll, ' but '//Un)
         GOTO 99
         ENDIF

      CALL OlfCnuG (j, 'z1', Co, Un, Fehler)
      IF (Fehler.ne.'&ff') GOTO 99
      IF     (Co.ne.'q') THEN
         Fehler = 'z1-coordinate not q but '//Co
         GOTO 99
      ELSEIF (Un.ne.'A-1') THEN
         Fehler = 'z1-unit not A-1 but '//Un
         GOTO 99
         ENDIF

C  Loop over q/ read spectra:

      nAC = iOlfG (j, '#spectra', Fehler) + 1
      IF (Fehler.ne.'&ff') THEN
         GOTO 99
      ELSEIF (nAC.le.2) THEN
         Fehler = 'nK<=2'
         GOTO 99
      ELSEIF (nAC.gt.MAC) THEN
         Print *, ' nAC = ', nAC, ', MAC = ', MAC
         Fehler = 'nK>MAC'
         GOTO 99
         ENDIF

      AC(1) = 0 ! scattering law for q=0 will be set later
      Sint(1) = 1
      SintMin = 1
      SintMax = 1

      DO iA = 2, nAC

         CALL OlfGet1Z (j, iA-1, 1, z, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         ! conversion A-1 -> meV -> dimensionless
         AC(iA) = (z/.694838)**2 / TempE ! checks below (iA=1 / iA>1)

         CALL OlfGetXY (j, iA-1, nC, Xin, Yin, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99

         IF (iA.eq.2) THEN ! obtain grid from spectrum 1 (->2)

            ! check w-scale :
            IF (irSorted(Xin,nC).ne.2) THEN
               Fehler ='Input energy transfer w not in ascending order'
               GOTO 99
               ENDIF
            IF (.not.(Xin(1).lt.0 .and. Xin(nC).gt.0)) THEN
               Fehler = 'input S~(qw) should cover w<0 and w>0'
               GOTO 99
               ENDIF

            ! obtain b-scale :
            IF (nC.gt.MBC) THEN
               Print *, ' nC = ', nC, ', MBC = ', MBC
               Fehler = 'External S(q,w) has too many channels'
               GOTO 99
               ENDIF
            nBC = nC
            DO iB = 1, nBC
               X0(iB) = Xin(iB)       ! for later comparison
               BC(iB) = X0(iB) * FacX / TempE
               ENDDO
            IF (BC(nBC).gt.80) THEN
               Print *, ' BC(nBC) = ', BC(nBC)
               Fehler =
     * 'too large energy transfers / cannot calculate e^-b/2'
               GOTO 99
               ENDIF
            SMALL = 1.d-8 * (X0(nBC)-X0(1))

            IF (AC(2).le.0) THEN
               Fehler =
     * 'Q<=0 in spectrum 1 (which becomes internal spectrum 2)'
               GOTO 99
               ENDIF

         ELSE ! spectra 2,..,nAC
            ! check b-scale :
            IF (nC.ne.nBC) THEN
               Fehler = '# channels of external S(q,w) is not constant'
               GOTO 99
               ENDIF
            DO iB = 1, nBC
               IF (dabs(Xin(iB)-X0(iB)).gt.SMALL) THEN
                  Fehler = 'different energy scales'
                  GOTO 99
                  ENDIF
               ENDDO
            IF (AC(iA).le.AC(iA-1)) THEN
               Fehler = 'Q not sorted'
               GOTO 99
               ENDIF

            ENDIF

         ! integrate :
         Sint(iA) = 0
         ymin = Yin(iB)
         ymax = Yin(iB)
         DO iB = 2, nBC
            Sint(iA) = Sint(iA) + (Xin(iB)-Xin(iB-1))*
     *                            (Yin(iB)+Yin(iB-1))/2
            ymin = dmin1 (ymin, Yin(iB))
            ymax = dmax1 (ymax, Yin(iB))
            ENDDO
         SintMin = dmin1 (SintMin, Sint(iA))
         SintMax = dmax1 (SintMax, Sint(iA))
         IF (ymin.lt.0 .and. -ymin.gt.1.d-4*ymax) THEN
            Print *, iA, ymin, ymax
            Fehler = 'negative S(q,w) given'
            RETURN
            ENDIF

         ! desymmetrize : S~ -> S
         DO iB = 1, nBC
            SAB(iA,iB) = dmax1 (0.d0,
     *           Yin(iB) / FacX * TempE * dexp (-BC(iB)/2) )
            ENDDO

         ENDDO ! iA

C  Set spectrum 1:
      iB0 = irPos (BC, nBC, 0.d0, 'n')
      DO iB = 1, nBC
         SAB(1,iB) = 0
         ENDDO
      SAB(1,iB0) = 2 / (BC(iB0+1)-BC(iB0-1))

C  Check a-scale :
      IF (irSorted(AC,nAC).ne.2) THEN
         Fehler = 'input q''s not in ascending order'
         GOTO 99
         ENDIF

C  Normalisation ?
      tol = 1.03
      IF (SintMax.gt.tol) THEN
         Print '(a,g10.3,a,g10.3)',
     *     ' integral ranges from ', SintMin, ' to ', SintMax
         Fehler = 'unphysical: Int[S(q,w)] > 1'
         RETURN
      ELSEIF (SintMin.lt.1/tol) THEN
         Print '(a,g10.3,a,g10.3)',
     *     ' integral ranges from ', SintMin, ' to ', SintMax
         Print *, 'WARNING / many neutrons will be lost'
c         DO iA = 1, nAC
c            DO iB = 1, nBC
c               SAB(iA,iB) = SAB(iA,iB) / Sint(iA)
c               ENDDO
c            ENDDO
      ELSE
         Print *, ' Fine. Ideal scattering law was normalised.'
         ENDIF

C  Check Alpha vs. Beta :
      Bmax = dmax1 (-BC(1), BC(nBC))
      IF (AC(nAC).lt.(3+2*dsqrt(2.0d0))*Bmax) THEN
         Print *, ' AC(nAC), |BC_max| : ', AC(nAC), Bmax
         Print *, ' please provide Q_max^2 > 2.81 |w| (in A-1, meV)'
         Fehler ='criterion  AC(nAC) >= BC(nBC)*(3+2*sqrt(2)) violated'
         GOTO 99
         ENDIF

      RETURN
 99   CONTINUE
      Print *, 'error in GetIE'
      RETURN

      END ! GetIE

C  ====================================================================
C
C      Module TG (target geometry)
C
C  ====================================================================

C  Contents of this module :
C      Geometry (GetGeom, RanImp, CalcDist, RanCol)
C      CutCyl, NewDir

C  This module has been completely rewritten (JWuttke November 1994)
C  in order to more flexibly accomodate new sample geometries

C  ====================================================================
      SUBROUTINE Geometry()
C  ====================================================================

C     Entries are :
C        GetGeom       Input target geometry, initializations
C        RanImp        Initialize neutron (random impact)
C        CalcDist      Calculate distance from collision point out of target
C        RanCol        Choose subsequent collision points

      IMPLICIT NONE                                           ! 13jan00
      INTEGER       MMat, MItem, MFragm
      PARAMETER    (MMat=3, MItem=100, MFragm=2*MItem)

      INCLUDE      'l_def.f'

      ! geometry
      INTEGER       ItemMat(MItem),FragmMat(MFragm),FragmRank(MFragm),
     *              iGeoTyp, iCylTyp, nTarBlock, iArrTyp, nTarItem
      REAL*8        ItemH2(MItem), ItemRo(MItem), ItemRi(MItem),
     *              ItemPos(MItem,3),
     *              Fragm0(MFragm), FragmD(MFragm),
     *              DirBeam(3),
     *              BeamW, BeamH, BeamC, BeamZ,
     *              TarH, TarR4, TarD4, TarW4, TarSep, TarAng

      ! local
      INTEGER       iBlock, nFragm, iFragm, iItem, ifail, i
      REAL*8        d, x, y, z, PP(3), XS1(3), XS2(3), P1,P1i,P2,P2i,
     *              Pmin, DEPS
      LOGICAL       qCut, qCuti

      ! parameters for calls :
      INTEGER       InMat
      REAL*8        AbsCo(3), DistSect(MMat), Dir(3), Dist, Thru,
     *              aux, survival, event, evlog
      LOGICAL       qHit
      CHARACTER*(*) Fehler

      ! external :
      !REAL*8        G05DAF !Artem: Replace with rand function from slatec/fnlib/rand.f
      REAL           ran !Artem add for rand()
      real, external :: rand !Artem add for rand()

      ! common :
      INTEGER       nPRT, jTrace
      REAL*8        ErgN, TimN, PosN, DirN
      LOGICAL       qHasSectn

      COMMON /Current/   ErgN, TimN, qHasSectn(MMat)
      COMMON /NeutPos/   PosN(3), DirN(3)
      COMMON /Trace/     nPRT, jTrace

      DATA DEPS /1d-10/, qHasSectn / .true., .false., .false. /

C  --------------------------------------------------------------------
      ENTRY GetGeom(Fehler)
C  --------------------------------------------------------------------

C        In this section, the user must specify the sample geometry.
C        The parameters are then transferred to internal arrays.

C  Dialogue :

      CALL Sage (' 2/2  Beam geometry')

      BeamW = rAskDMu ('Width  of incoming beam (cm)',
     *                 BeamW, 1.d-2, 1.d2)
      BeamC = rAskDMu ('Centered at', BeamC, -1.d2, +1.d2)
      BeamH = rAskDMu ('Height of incoming beam (cm)',
     *                 BeamH, 1.d-2, 1.d2)
      BeamZ = 0

      CALL Sage (' 2/3  Target geometry')

      Print *,' The target may contain sample, '//
     * 'container and absorber sections:'
      qHasSectn(1) = qAskD ('Is there sample material (S)',
     *     intq(qHasSectn(1)))
      qHasSectn(2) = qAskD ('Is there container material (C)',
     *     intq(qHasSectn(2)))
      qHasSectn(3) = qAskD ('Is there absorber material (A)',
     *     intq(qHasSectn(3)))

      Print *, ' Kind of geometry :'
      Print *,
     * '   (1) cylindric (axes perpendicular to scattering plane)'
      Print '(a)', '=> at present, the only choice is (1)'
      iGeoTyp = 1 ! iAskMu ('Option ?', 1, 1)
      Print *

      IF (iGeoTyp.eq.1) THEN
         Print *, ' Type of cylinders :'
         Print *, '   (1) solid rods (SC)'
         Print *, '   (2) hollow tubes (-CSC)'
         iCylTyp = iAskDMu ('Option ?', iCylTyp, 1, 2)
         Print *

         nTarBlock = iAskDMu ('Number of rods/tubes', 1, 1, 1000)
         Print *

         Print *, ' Spatial arrangement :'
         Print *, '   (1) linear, one row'
         Print *, '   (2) circular, one band'
         iArrTyp = 1 ! iAskMu ('Option ?', 1, 2)
         Print *, ' => at present, the only choice is (1)'
         Print *

         TarH  = rAskD    ('Heigth of rods/tubes (cm) ?', TarH)

         IF     (iCylTyp.eq.1) THEN
            TarR4 = rAsk    ('Outer radius of sample volume (cm) ?')
            IF (qHasSectn(2))
     *         TarD4 = rAsk    ('Wall thickness of container (cm) ?')
         ELSEIF (iCylTyp.eq.2) THEN
            TarR4 = rAskD    ('Outer radius of sample volume (cm) ?',
     *                        TarR4)
            TarW4 = rAskDMu  ('Wall thickness of sample (cm) ?',
     *                        TarW4, 0.d0, TarR4)
            IF (qHasSectn(2)) TarD4 = rAskMu  (
     * 'Wall thickness of container (cm) ?', TarD4, TarR4-TarW4)
            ENDIF

         IF (nTarBlock.gt.1) THEN
            TarSep = rAsk    ('Tube separation (cm) ?')
            TarAng = rAskM (
     * 'Angle between transmitted beam and target plane',
     * TarAng, -180.d0, 180.d0)
            ENDIF

         ENDIF ! iGeoTyp / cylindric

C  Prepare arrays for later use :

C  - Beam direction:
      DirBeam(1) = 0
      DirBeam(2) = 1
      DirBeam(3) = 0

C  - Sample items :
      IF (iGeoTyp.eq.1) THEN
         nTarItem = 0
         DO iBlock = 1, nTarBlock
            ! coordinates of block center :
            d = (iBlock-1-float(nTarBlock)/2) * TarSep
            x = d * dsind(TarAng) ! a verifier
            y = d * dcosd(TarAng) ! a verifier
            ! items in block :
            IF     (iCylTyp.eq.1) THEN  ! solid cylinder
               IF (qHasSectn(1)) ! sample section
     *            CALL AddItemCyl (MItem, nTarItem,
     * ItemMat, ItemH2, ItemRo, ItemRi, ItemPos,
     * 1, TarH/2, TarR4, 0d0, x, y, 0d0, Fehler)  !Artem add Fehler
               IF (qHasSectn(2)) ! container around solid cylinder
     *            CALL AddItemCyl (MItem, nTarItem,
     * ItemMat, ItemH2, ItemRo, ItemRi, ItemPos,
     * 2, TarH/2, TarR4+TarD4, TarR4, x, y, 0d0, Fehler) !Artem add Fehler

            ELSEIF (iCylTyp.eq.2) THEN  ! hollow cylinder
               IF (qHasSectn(1)) THEN ! sample section
                  CALL AddItemCyl (MItem, nTarItem,
     * ItemMat, ItemH2, ItemRo, ItemRi, ItemPos,
     * 1, TarH/2, TarR4, TarR4-TarW4,
     * x, y, 0d0, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               IF (qHasSectn(2)) THEN ! outer container
                  CALL AddItemCyl (MItem, nTarItem,
     * ItemMat, ItemH2, ItemRo, ItemRi, ItemPos,
     * 2, TarH/2, TarR4+TarD4, TarR4,
     * x, y, 0d0, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               IF (qHasSectn(2)) THEN ! inner container
                  CALL AddItemCyl (MItem, nTarItem,
     * ItemMat, ItemH2, ItemRo, ItemRi, ItemPos,
     * 2, TarH/2, TarR4-TarW4, TarR4-TarW4-TarD4,
     * x, y, 0d0, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDIF
               ENDIF
            ENDDO

         ENDIF ! iGeoTyp / cylindric

C  Print-out :
      Write (nPRT, '(/a/17(1h-)/)') 'Target geometry :'

      Write (nPRT, '(a,f10.3)') ' tube length (cm)          ', TarH
      Write (nPRT, '(a,f10.3)') ' tube outer radius (cm)    ', TarR4
      IF (iCylTyp.eq.2)
     *Write (nPRT, '(a,f10.3)') ' sample thickness (cm)     ', TarW4
      IF (qHasSectn(2))
     *Write (nPRT, '(a,f10.3)') ' container thickness (cm)  ', TarD4
      Write (nPRT, '(a,f10.3)') ' tube separation (cm)      ', TarSep
      Write (nPRT, '(a,i10  )') ' number of tubes           ',nTarBlock
      Write (nPRT, '(a,i10  )') ' number of items           ', nTarItem
      Write (nPRT, '(a,f10.3)') ' angle beam><tube-plane    ', TarAng
      Write (nPRT, '(a,3f10.3)')' incident beam cosines     ', DirBeam
      Write (nPRT, '(2(a,f10.3))')' beam dimensions (cm)      ', BeamH,
     *                            '    by',                      BeamW
      Write (nPRT, '(2(a,f10.3))')' centered in               ', BeamZ,
     *                            '    by',                      BeamC

      RETURN ! End *GetGeom

C  --------------------------------------------------------------------
      ENTRY RanImp (Dist, Thru, qHit)
C  --------------------------------------------------------------------

C        In this entry, a new neutron is given a direction,
C        and an impact position on the surface of the target is chosen.
C        Output :     Dist : distance to impact position, relative to the
C                            distance from the neutron source to point 0,0,0.
C                     Thru : pathlength thru sample section (for normalisation)
C                     qHit  : whether the target has been qHit at all.
C        OverWrites : PosN : Present neutron position

C  Set incident direction = beam direction :
         ! this could be modified including beam divergence
      DirN(1) = DirBeam(1)
      DirN(2) = DirBeam(2)
      DirN(3) = DirBeam(3)

C  Choose impact in a plane perpendicular to DirN :
C      z = G05DAF (BeamZ-BeamH/2, BeamZ+BeamH/2)!Artem: Replace with rand function from slatec/fnlib/rand.f
      ran = rand(0.0)
      z = BeamZ-BeamH/2 + (BeamZ+BeamH/2 - BeamZ-BeamH/2) * dble(ran)
C      d = G05DAF (BeamC-BeamW/2, BeamC+BeamW/2)!Artem: Replace with rand function from slatec/fnlib/rand.f
      ran = rand(0.0)
      d = BeamC-BeamW/2 + (BeamC+BeamW/2 - BeamC-BeamW/2) * dble(ran)

C  Transform to T-CS :
      IF (DirN(3).ne.0) CALL Absturz ('Impact',
     *     'General case Dir(3)<>0 not allowed')
      PP(1) =   d * DirN(2)
      PP(2) = - d * DirN(1)
      PP(3) =   z

C  Find the target item the neutron qHits first :
      qHit = .false.
      Thru = 0
      DO iItem = 1, nTarItem
         ! search for item with minimal P1 :
         CALL CutCyl (PP, DirN, ItemPos(iItem,1), ItemH2(iItem),
     *                ItemRo(iItem), qCut, P1, XS1, P2, XS2)
         IF (qCut) THEN
            IF (.not.qHit) THEN
               Pmin = P1
               ENDIF
            IF (.not.qHit .or. P1.lt.Pmin) THEN
               Pmin = P1
               ! overWrite neutron position
               DO i = 1, 3
                  PosN(i) = XS1(i)
                  ENDDO
               ENDIF
            qHit = .true.
            ! increment pathlength thru sample (for normalization) :
            !  IF mat=S ??????????????????
            CALL CutCyl (PP, DirN,
     *                 ItemPos(iItem,1), ItemH2(iItem), ItemRi(iItem),
     *                 qCuti, P1i, XS1, P2i, XS2) ! XS* is dummy, not needed
            IF (qCuti) THEN
               Thru = Thru + P2-P2i + P1i-P1
            ELSE
               Thru = Thru + P2 - P1
               ENDIF
            ENDIF
         ENDDO
      Dist = Pmin

C  Info / Debug :
      IF (jTrace.ge.1) THEN
         Write (nPRT, '(//a/a,3f11.6)') '>>>  NEW NEUTRON  <<<',
     *                                  '     impact: ', PosN
         IF (.not.qHit) THEN
            Write (nPRT, '(a)') '    missed sample'
            ENDIF
         ENDIF

      RETURN ! End *RanImp

C  --------------------------------------------------------------------
      ENTRY CalcDist (Dir, DistSect)
C  --------------------------------------------------------------------

            ! completely rewritten (JWu 17/18nov94)
C        In this entry, the distances to travel thru the different
C          sample sections are calculated. The intersections are saved
C          (and later used in RanCol).
C        Input :      Dir(3)     : a direction (REAL*8)
C        Output :     DistSect   : distances travelled thru S, C, A (R*4)
C        OverWrites :

      IF (jTrace.ge.5) Write (nPRT, '(2(a,3f8.4))')
     *      '    calc dist from ', PosN, ' in dir ', Dir

      DO i = 1, 3
         DistSect(i) = 0
         ENDDO
      nFragm = 0
      DO iItem = 1, nTarItem
         ! Is there an intersection with the item's outer cylinder ?
         CALL CutCyl (PosN, Dir, ItemPos(iItem,1), ItemH2(iItem),
     *                ItemRo(iItem), qCut, P1, XS1, P2, XS2)
         IF (.not.qCut) GOTO 9 ! next item
         IF (P2.lt.0) GOTO 9 ! item is lying behind the neutron
         IF (P1.lt.0) P1 = 0 ! neutron is *inside* the item
         ! Yes, the item has been hit - maybe, the inner cylinder is hit, too ?
         CALL CutCyl (PosN, Dir, ItemPos(iItem,1), ItemH2(iItem),
     *                ItemRi(iItem), qCuti, P1i, XS1, P2i, XS2) ! XS* is dummy, not needed
         IF (qCuti .and. P2i.lt.0) qCuti = .false. ! hole is lying behind
         IF (qCuti .and. P1i.lt.0) P1i = 0 ! we are *in* the hole
         ! distance travelled thru one item :
         IF (qCuti) THEN
            d = P2-P2i + P1i-P1
         ELSE
            d = P2 - P1
            ENDIF
         ! Info / Debug :
         IF (jTrace.ge.4) THEN
            IF (qCuti) THEN
               Write (nPRT, '(a,i2,a,4f8.4)') '    -> line cuts item ',
     *      iItem, ' twice/ paths: ', P1, P1i, P2i, P2
               ! line cuts twice - but one section may already lie behind
            ELSE
               Write (nPRT, '(a,i2,a,2f8.4)')
     *      '    -> line cuts item ', iItem, ' once/ path is ', P1, P2
               ENDIF
            ENDIF
         ! Increment the corresponding output array :
         DistSect(ItemMat(iItem)) = DistSect(ItemMat(iItem)) + d
         ! Save intersection points and path lengths for use in RanCol :
         nFragm = nFragm + 1
         Fragm0(nFragm) = P1
         FragmMat(nFragm) = ItemMat(iItem)
         IF (qCuti) THEN
            FragmD(nFragm) = P1i - P1
            nFragm = nFragm + 1
            Fragm0(nFragm) = P2i
            FragmMat(nFragm) = ItemMat(iItem)
            FragmD(nFragm) = P2 - P2i
         ELSE
            FragmD(nFragm) = P2 - P1
            ENDIF
 9       CONTINUE
         ENDDO

      IF (nFragm.le.0) THEN
         CALL Absturz ('CalcDist',
     *                 'No target section intersected (nFragm=0)')
         ENDIF

      RETURN ! End *CalcDist

C  --------------------------------------------------------------------
      ENTRY RanCol (AbsCo, Dist, InMat)
C  --------------------------------------------------------------------

            ! completely rewritten (JWu 18nov94)
C        In this entry, a new collision point is choosen.
C        Input :      AbsCo(Mat)   : total cross sections
C        Output :     Dist         : travelled distance
C                     InMat        : it happened in S, C, A
C        OverWrites : PosN(3)      : neutron position

C  Rank fragments according to the intersections Fragm0() :
      ifail = 0
      CALL M01DAF_local (Fragm0, 1, nFragm, 'a', FragmRank, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) CALL Absturz ('RanCol', 'NAG M01DAF error')

C  Rearrange them according to their rank :
      ifail = 0
      CALL M01EAF_local (Fragm0, 1, nFragm, FragmRank, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) CALL Absturz ('RanCol', 'NAG M01EAF error')

      ifail = 0
      CALL M01EAF_local (FragmD, 1, nFragm, FragmRank, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) CALL Absturz ('RanCol', 'NAG M01EAF error')

      ifail = 0
      CALL M01EBF_local (FragmMat, 1, nFragm, FragmRank, ifail) !Artem: Replace with a self-made subroutine from lnag_local.f
      IF (ifail.ne.0) CALL Absturz ('RanCol', 'NAG M01EBF error')

      IF (jTrace.ge.3) THEN
         DO iFragm = 1, nFragm
            Write (nPRT, '(a,2i3,3f9.4)') !Artem '(a,2i3,3f9.4))' -> '(a,2i3,3f9.4))'
     *       '& Fragm no. mat from thru mu*d',
     *       iFragm, FragmMat(iFragm), Fragm0(iFragm), FragmD(iFragm),
     *       FragmD(iFragm)*AbsCo(FragmMat(iFragm))
            ENDDO
         ENDIF

C  Calculate survival probability :
      aux = 0
      DO iFragm = 1, nFragm
         aux = aux + AbsCo(FragmMat(iFragm))*FragmD(iFragm)
         ENDDO
      survival = dexp (-aux)

C  Random choice :
C      event = G05DAF (survival, 1.d0) !Artem: Replace with rand function from slatec/fnlib/rand.f
      ran = rand(0.0)
      event = survival + (1.d0 - survival) * dble(ran)
      IF (jTrace.ge.3) Write (nPRT, '(2(a,f7.4))')
     *     '& RanCol/ random x = ', event, ' where min = ', survival
      aux = 0
      Dist = 0 ! needed for tof increment
      evlog = -dlog (event)
      DO iFragm = 1, nFragm
         aux = aux + AbsCo(FragmMat(iFragm))*FragmD(iFragm)
         IF (aux.ge.evlog) GOTO 32
         Dist = Dist + FragmD(iFragm)
         ENDDO
 32      CONTINUE ! it happened in fragment no. iFragm
      IF (jTrace.ge.3) Write (nPRT, '(a,i3,4(a,f8.4))')
     *     '& RanCol/ in Fragm.', iFragm, ' aux=', aux, ', evlog=',
     *     evlog, ', dist(before)=', Dist, ' dist(here)=', d
      d = FragmD(iFragm) - (aux-evlog)/AbsCo(FragmMat(iFragm))
      Dist = Dist + d

C  New position :
      DO i = 1, 3
         PosN(i) = PosN(i) + DirN(i) * (Fragm0(iFragm) + d)
         ENDDO
      InMat = FragmMat(iFragm)

C  Info / Debug :
      IF (jTrace.ge.2) THEN
         Write (nPRT, '(a,3f9.4)') !Artem '(a,3f9.4))' -> '(a,3f9.4)'
     *      '>>>     scattered in ', PosN(1), PosN(2), PosN(3)
         ENDIF

      RETURN ! End *RanCol

      END ! Geometry

C  ====================================================================
      SUBROUTINE AddItemCyl (MItem, nTarItem,
     *                       ItemMat, ItemH2, ItemRo, ItemRi, ItemPos,
     *                       jMat, H2, Ro, Ri, Px, Py, Pz, Fehler)
C  ====================================================================

      IMPLICIT NONE

      INTEGER    MItem, nTarItem, ItemMat(MItem), jMat
      REAL*8     ItemH2(MItem), ItemRo(MItem), ItemRi(MItem),
     *           ItemPos(MItem,3), H2, Ro, Ri, Px, Py, Pz
      CHARACTER  Fehler*(*)

      nTarItem = nTarItem + 1
      IF (nTarItem.gt.MItem) THEN
         Fehler = 'Too many target items'
         RETURN
         ENDIF
      ItemMat(nTarItem) = jMat
      ItemH2(nTarItem) = H2
      ItemRo(nTarItem) = Ro
      ItemRi(nTarItem) = Ri
      ItemPos(nTarItem,1) = Px
      ItemPos(nTarItem,2) = Py
      ItemPos(nTarItem,3) = Pz

      END ! AddItemCyl

C  ====================================================================
      SUBROUTINE CutCyl (X0, D0, XM, Z2, R, qCut, P1, XS1, P2, XS2)
C  ====================================================================

            ! Completely new (JWuttke 17nov94)
         ! Calculate the intersection of a line with a cylinder
         !   (of finite height).
         ! Input :  X0, D0  define the line
         !          XM      center of the cylinder
         !          Z2      half height
         !          R       radius
         ! Output : Cut     Are there intersections ?
         !          P1      path length from X0 to XS1
         !          XS1     first intersection point
         !          P2, XS2 dito, second intersection
         ! Note that P* may be negative. Intersections are ordered
         !   with respect to the direction D0. We have always P1<=P2.

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      REAL*8   X0(3), D0(3), XM(3), XS1(3), XS2(3)
      REAL*8   XX(3)

      DATA     ZEROCOS / 1d-6 /

      qCut = .false.
      IF (R .le. 1.d-14) CALL Absturz (
     * 'CutCyl', 'Cylindre infinement mince')
      ZERODIST = R * ZEROCOS

C  Intersection with curved surface ?
      aabb = D0(1)**2 + D0(2)**2        ! in-plane-cos**2
      DO i = 1, 3
         XX(i) = XM(i) - X0(i)
         ENDDO
      IF (dabs(aabb).lt.ZEROCOS) THEN ! perpendicular direction
C  Exception : perpendicular qCut, intersection with front planes :
         IF (XX(1)**2 + XX(2)**2 .lt. R**2) THEN ! X0 inside infinite cylinder
            D0(3) = D0(3)/dabs(D0(3)) ! D0(3) is +1 or -1
            P1 = XX(3)/D0(3) - Z2
            P2 = XX(3)/D0(3) + Z2
            XS1(1) = 0
            XS1(2) = 0
            XS1(3) = XM(3) - Z2/D0(3)
            XS2(1) = 0
            XS2(2) = 0
            XS2(3) = XM(3) + Z2/D0(3)
            qCut = .true.
         ! ELSE : there is no intersection
            ENDIF
         RETURN ! everything done
         ENDIF
C  Continue with normal case :
      xxab = XX(1)*D0(1) + XX(2)*D0(2)  ! in-plane-distance XM-X0
      P0 = xxab / aabb ! distance X0-(closest point to cyl axis)
      aux = xxab**2 - aabb * (XX(1)**2 + XX(2)**2 - R**2)
      IF (aux .le. ZERODIST**2) RETURN ! no useful pair of intersections

      qCut = .true.
      delta = dsqrt(aux) / aabb ! half path thru sample
      P1 = P0 - delta
      P2 = P0 + delta

      DO i = 1, 3 ! set intersection points :
         XS1(i) = X0(i) + P1 * D0(i)
         XS2(i) = X0(i) + P2 * D0(i)
         ENDDO

C  Now check for the z coordinate :
      q1low = XS1(3) .lt. XM(3)-Z2
      q1hig = XS1(3) .gt. XM(3)+Z2
      q2low = XS2(3) .lt. XM(3)-Z2
      q2hig = XS2(3) .gt. XM(3)+Z2

      IF (.not. (q1low .or. q1hig .or. q2low .or. q2hig)) RETURN ! everything oK

C  At least one intersection point not in the cylindrical wall:
      IF ((q1low .and. q2low) .or. (q1hig .and. q2hig)) THEN
         ! neutron passes below or above the sample
         qCut = .false.
         RETURN
         ENDIF
      ! Now we can be sure that D0(3)<>0.
      IF (D0(3).eq.0.d0) CALL Absturz ('CutCyl', 'z-cos<>0')

C  Now solve the problems one by one :
      IF (q1low) P1 = (XX(3)-Z2) / D0(3)
      IF (q1hig) P1 = (XX(3)+Z2) / D0(3)
      IF (q2low) P2 = (XX(3)-Z2) / D0(3)
      IF (q2hig) P2 = (XX(3)+Z2) / D0(3)

      DO i = 1, 3 ! reset intersection points :
         XS1(i) = X0(i) + P1 * D0(i)
         XS2(i) = X0(i) + P2 * D0(i)
         ENDDO

      END ! CutCyl

C  ====================================================================
      SUBROUTINE NewDir (Dir, scaC, scaS, aziC, aziS)
C  ====================================================================
         ! change direction Dir by angles (sca,azi)
         ! for each angle, cosine and sine are given

      IMPLICIT NONE
      REAL*8      Dir(3), Dir2(3), scaC, scaS, aziC, aziS, sumsq, aux2

      sumsq=Dir(1)*Dir(1)+Dir(2)*Dir(2)
      IF (sumsq.gt.1e-12) THEN
         aux2 = dsqrt(sumsq)
         Dir2(1)=( Dir(2)*aziC+Dir(3)*Dir(1)*aziS)*scaS/aux2 +
     *             Dir(1)*scaC
         Dir2(2)=(-Dir(1)*aziC+Dir(3)*Dir(2)*aziS)*scaS/aux2 +
     *             Dir(2)*scaC
         Dir2(3)=    -aux2*scaS*aziS             + Dir(3)*scaC
      ELSE ! old direction was (0,0,+-1).
         Dir2(1)=scaS*aziC
         Dir2(2)=scaS*aziS
         Dir2(3)=scaC*Dir(3)
         ENDIF
      ! ensure normalisation of Dir :
c      sumsq = Dir2(1)*Dir2(1) + Dir2(2)*Dir2(2) + Dir2(3)*Dir2(3)
c      aux2 = 1 / dsqrt(sumsq)
c      Print '(a,f14.10)', ' normalisation of new direction: ', aux2
      Dir(1)=Dir2(1) !*aux2
      Dir(2)=Dir2(2) !*aux2
      Dir(3)=Dir2(3) !*aux2

      END ! NewDir

C  ====================================================================
C
C      Module TF (time-of-flight, random-distributions, misc.)
C
C  ====================================================================

C  Contents of this module :
C      IndexTim
C      GetInst, ImpactTOF
C      GetMProp
C      Index1, Index2
C      Alpha_BEC, Alpha_Base, RandAng, TruncGauss

C  ====================================================================
      INTEGER FUNCTION IndexTim (MEC,THIGH,TimOut,DTimOut,EC,DEC,nEC,T)
C  ====================================================================

            ! Find index IndexTim=K of T in TimOut(1..nEC).
            ! Return K=0 if T is not contained in TimOut+-DTimOut/2.
            ! Used to determine the elastic channel.

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      REAL*8         THIGH (MEC), TimOut(*), DTimOut(*), EC(*), DEC(*)

      COMMON /Trace/ nPRT, jTrace

C  NEEDED ONLY ON FIRST CALL :::
      ! It is assumed that TimOut is in DEcreasing order.
      IF (jTrace.ge.7) Write (nPRT, '(a,i3)')
     *     '&&& IndexTim : Inverting I = 1, ...,', nEC
      DO I = 1, nEC
         THIGH(nEC+1-I) = TimOut(I) + 0.5*DTimOut(I)
         ENDDO

C  Find index K of T in THIGH :
      IF (jTrace.ge.7) Write (nPRT, '(a,f9.4,a,f9.4,a,i3,a,f9.4)')
     *      '&&& IndexTim/ going to search T=', T, ' between TH(1)=',
     *      THIGH(1), ' and TH(', nEC, ')=', THIGH(nEC)

      IF     (T.lt.THIGH(1)) THEN
         K = 1 ! maybe T.gt.TLOW ...
      ELSEIF (T.gt.THIGH(nEC)) THEN
         IndexTim = 0
         IF (jTrace.ge.5) Write (nPRT, '(a)') '& neutron lost'
         RETURN     ! neutron lost
      ELSE
         K = Index1 (THIGH, nEC, T) ! => T<=THIGH(K)
         ENDIF

      IndexTim = nEC+1-K
      IF (jTrace.ge.7) Write (nPRT, '(a,i3)') '&&& => IndexTim=',
     *                 IndexTim

      IF (T.lt.TimOut(IndexTim)-0.5*DTimOut(IndexTim)) IndexTim = 0 ! lost

      IF (jTrace.ge.7) Write (nPRT, '(a,i3)') '&&& => returning ',
     *                 IndexTim

      END ! IndexTim

C  ====================================================================
      SUBROUTINE GetInst (qIgnToF, DistSD, Erg0, Erg0width, Fehler)
C  ====================================================================
            ! JWu 5mar96 : backscattering or time-of-flight
         ! ask for instrumental setup, read resolution function

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      PARAMETER    (MRes=256)
      REAL*8        XRes(MRes), YRes(MRes), PRes(MRes)
      CHARACTER     Fehler*(*), FileRes*80

      COMMON /Convert/   TempE, ConvVof
      COMMON /Trace/     nPRT, jTrace
      COMMON /ErgIn/     Vof0, Vof0del, Vof0min, Vof0max,
     *                   DistCS, Tim0del, Tim0min, Tim0max,
     *                   XRes, PRes, nRes, qRdRes

      CALL Sage (' 2/1  Instrumental setup')
      Vof0 = sqrt (ConvVof/Erg0)

      iInstTyp = iAskDMu (
     * 'Perfect Resolution (0), Backscattering (1), '//
     * 'Time-of-flight (2)',
     * iInstTyp, 0, 2)

      IF     (iInstTyp.eq.0) THEN
         qIgnToF = .true.
         qRdRes  = .false.
         Erg0width = 0
      ELSEIF (iInstTyp.eq.1) THEN
         qIgnToF = .true.
         qRdRes  = .true.
      ELSEIF (iInstTyp.eq.2) THEN
         qIgnToF = .false.
         qRdRes  = qAskDi (
     * 'Resolution from normal distribution(0) or from file (1)', 1)
         IF (.not.qRdRes) THEN
            CALL Sage (
     * ' Standard deviations for Gaussian distributions :')

            Vof0rel = rAsk ('Spread in v-1 of incident beam (%) ?')
            Tim0del = rAsk ('Time spread at the chopper (usec) ?')
            DistCS  = rAsk ('Distance chopper to sample (cm) ?')

            ! Mean inverse velocity in usec/cm :
            Vof0del = Vof0 * Vof0rel/100
            ! Truncation for the inverse velocity Vof :
            Vof0max = Vof0 + 5*Vof0del
            Vof0min = Vof0 - 5*Vof0del
            ! Truncation for time TINC at sample :
            Tim0max =  5*Tim0del
            Tim0min = -5*Tim0del
            IF (Vof0min.le.0.) THEN
               Fehler = 'Vof_variance too large => Vof<0 not excuded'
               RETURN
               ENDIF
            ! Largest permitted incident energy (for check against grids) :
            Erg0width = dmax1 (Erg0-ConvVof/Vof0max**2,
     *                         ConvVof/Vof0min**2-Erg0)

            ENDIF
         ENDIF

      DistSD  = rAskDMu ('Distance sample to detector (cm) ?',
     *                   DistSD, 1.d0, 1.d5)

      IF (qRdRes) THEN
         CALL FrageC ('Read resolution from file', FileRes)
c         CALL ReadArray (FileRes, 'w', 'meV', 'S', 'meV-1',
c     *                   MRes, nRes, XRes, YRes, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ! integrate and normalise resolution function :
         PRes(1) = 0
         DO i = 2, nRes
            ym = (YRes(i)+YRes(i-1))/2
            IF (ym.le.0) THEN
               Fehler = 'Resolution contains S<0'
               RETURN
               ENDIF
            PRes(i) = PRes(i-1) + (XRes(i)-XRes(i-1)) * ym
            ENDDO
         DO i = 1, nRes
            PRes(i) = PRes(i) / PRes(nRes)
            ENDDO

         ENDIF

      END ! GetInst

C  ====================================================================
      SUBROUTINE ImpactTOF (qIgnToF, Erg0, Tim, Erg,
     *                      Erg0width, ErgMin, ErgMax)
C  ====================================================================
         ! Initialize a neutron : set Tim and Erg (fix or random).

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INTEGER       MRes, nRes
      PARAMETER    (MRes=256)
      REAL*8        XRes(MRes), PRes(MRes)
      real, external :: rand !Artem add for rand()

      COMMON /Convert/   TempE, ConvVof
      COMMON /Trace/     nPRT, jTrace
c      COMMON /Resoltn/   qIgnToF, DistSD
      COMMON /ErgIn/     Vof0, Vof0del, Vof0min, Vof0max,
     *                   DistCS, Tim0del, Tim0min, Tim0max,
     *                   XRes, PRes, nRes, qRdRes

 1    CONTINUE
      IF (qRdRes) THEN              ! ---- resolution was read from file
         Tim = 0
         ran = dble(rand(0.0)) !Artem: Replace with rand function from slatec/fnlib/rand.f: G05CAF(dummy)
         ir  = irPos (PRes, nRes-1, ran, 'r')
         Erg = Erg0 + XRes(ir) +
     * (ran-PRes(ir))/(PRes(ir+1)-PRes(ir))*(XRes(ir+1)-XRes(ir))
         IF (Erg.le.ErgMin .or. Erg.ge.ErgMax) GOTO 1
      ELSEIF (Erg0width.eq.0) THEN ! ---- no energy distibution
         Tim = 0
         Erg = Erg0
      ELSE                         ! ---- Gaussian energy distribution
         CALL TruncGauss (Vof, Vof0, Vof0del, Vof0min, Vof0max)
         Erg = ConvVof / Vof**2
         IF (Erg.le.ErgMin .or. Erg.ge.ErgMax) GOTO 1

         CALL TruncGauss (Tim, 0.d0, Tim0del, Tim0min, Tim0max)
         Tim = Tim + (Vof-Vof0)*DistCS
         ENDIF

      END ! ImpactTOF

C  ====================================================================
      SUBROUTINE GetMProp (jMat, Name, SiMatSC, SiMatAB, Fehler)
C  ====================================================================
            ! As a subroutine JWu 25jul94 (loop over materials)
         ! Get material properties.

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)
      INTEGER       nPRT, jTrace
      CHARACTER*(*) Name, Fehler
      CHARACTER     ch1*1, Under*80

      COMMON /Trace/     nPRT, jTrace

      DATA   Avogadro / 6.022045d23 /

      IF (Fehler.ne.'&ff') THEN
         Print *, 'ERROR on entry in GetMProp'
         RETURN
         ENDIF

      lName = lenU(Name)
      CALL Sage (' 3/'//ch1(jMat)//' Scattering properties of the '//
     *           Name(1:lName)//' :')
      Write (nPRT, '(/2a/)') 'Properties of the ', Name(1:lName)

      ! Re-convert old cross sections for use as default :
      IF (RhoNum.gt.0) THEN
         SiMatAB = SiMatAB / ( 1d-24 * RhoNum )
         SiMatSC = SiMatSC / ( 1d-24 * RhoNum )
         ENDIF

      SiMatSC = rAskD ('Bound cross section per molecule (barns) ?',
     * SiMatSC)
      SiMatAB = rAskD (
     * 'Absorption cross section for 2200m/sec neutrons ?',
     * SiMatAB)
      WtMol = rAskDMu('Molecular weight (0=use number density) ?',
     * WtMol, 0.d0, 1.d6)
      IF (WtMol.le.0.) THEN
         RhoNum = rAskD ('Number density (1/cm3) ?', RhoNum)
      ELSE
         RhoNum = rAsk ('Density (g/cm3) ?')
         RhoNum = RhoNum / WtMol*Avogadro
         ENDIF

      Write (nPRT, '(a,f10.5)') ' bound cross section (barn)    ',
     *       SiMatSC
      Write (nPRT, '(a,f10.5)') ' absorption cross section      ',
     *       SiMatAB
      Write (nPRT, '(a,g10.5)') ' number density (cm-3)         ',
     *       RhoNum
      Write (nPRT, '(a,f10.5)') ' molecular weight              ',
     *       WtMol

      ! Convert cross sections from barn/cm3 to cm-1
      SiMatAB = SiMatAB * 1d-24 * RhoNum
      SiMatSC = SiMatSC * 1d-24 * RhoNum

      IF (SiMatSC.gt.0)
     * Write (nPRT, '(a,f10.5)') ' scattering free path (cm)     ',
     *        1./SiMatSC
      IF (SiMatAB.gt.0)
     * Write (nPRT, '(a,f10.5)') ' absorption free path          ',
     *        1./SiMatAB

      END ! GetMProp

C  ====================================================================
      INTEGER FUNCTION Index1 (X, N, xi)
C  ====================================================================
         ! Find the index Index1 of xi in the array X(1..N),
         ! such that  X(Index1-1) < xi <= X(Index1).
         ! Absturz if xi < X(1) or xi > X(N).

         ! Original version : interval halving technique

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      REAL*8  X(*)

      IF     (xi.lt.X(1))  THEN
         Print '(a,i3,2x,3(1x,g10.4),3x,g10.4)',
     *       ' N X(1,2,N) x : ', N, X(1), X(2), X(N), xi
         CALL Absturz ('Index1', 'xi<X(1)')
      ELSEIF (xi.gt.X(N)) THEN
         Print '(a,i3,5g12.4)',
     *       ' N X(1,2,N-1,N) x : ', N, X(1), X(2), X(N-1), X(N), xi
         CALL Absturz ('Index1', 'xi>X(N)')
         ENDIF

      KLO=0
      KHI=N
  200 CONTINUE
      IF (KHI-KLO.le.1) THEN
         Index1 = KHI
         RETURN
         ENDIF
      KCHK = (KHI+KLO)/2
      IF (xi.le.X(KCHK)) THEN
         KHI=KCHK
      ELSE
         KLO=KCHK
         ENDIF
      GOTO 200

      END ! Index1

C  ====================================================================
      INTEGER FUNCTION Index2 (X, N, xind, ILAST)
C  ====================================================================

         ! Find the index Index2 of xind in the array X(1..N),
         ! such that  X(Index2-1) < xind <= X(Index2).

         ! Accelerated version : start search at ILAST.
         ! Absturz if xind < X(0) or xind > X(N).
         ! No checks of N or ILAST.

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      REAL*8  X(*)

      IF (xind.gt.X(ILAST)) THEN
         DO I = ILAST+1, N
            IF (xind.le.X(I)) THEN
               Index2 = I
               RETURN
               ENDIF
            ENDDO
         IF (xind.lt.X(N)*(1+1.d-8)) THEN
            Index2 = N
            Print *, 'WARNING/ Index2 above limits'
            RETURN
            ENDIF
         Print '(a,2i3,4(1x,g10.4))', ' N ilast X(1,N) X(ilast) x0 = ',
     *           N, ILAST, X(1), X(N), X(ILAST), xind
         CALL Absturz ('Index2', 'xind>XN')

      ELSE
         DO I = ILAST-1, 1, -1
            IF (xind.gt.X(I)) THEN
               Index2 = I+1
               RETURN
               ENDIF
            ENDDO
         IF (xind.ge.X(1)) THEN
            Index2 = 2
         ELSE
         Print '(a,i3,3(1x,g10.4),3x,g10.4)',
     *       ' N ilast X(1,2,N) x = ', N, ILAST, X(1), X(2), X(N), xind
            CALL Absturz ('Index2', 'xind<X(1)')
            ENDIF

         ENDIF

      END ! Index2

C  ====================================================================
      REAL*8 FUNCTION Alpha_BEC (beta, eps, scos)
C  ====================================================================
         ! calculate alpha for given
         !     energy transfer / temperature = beta,
         !     incoming energy / temperature = eps,
         !     scattering angle cosine       = scos.

      IMPLICIT REAL*8   (a-h,o-p,r-z)

      IF (eps+beta.lt.0) CALL Absturz ('Alpha_BEC', 'arg(sqrt)<0')
      Alpha_BEC = dmax1 (0.d0, 2*eps + beta -
     *                   2 * dsqrt(eps*(eps+beta)) * scos)

      END ! Alpha_BEC

C  ====================================================================
      REAL*8 FUNCTION Alpha_Base (rL, rH, gL1, gL2, gH1, gH2, Fehler)
C  ====================================================================
            ! reestablished 13jan00

         !  <--- gL2 --->        <--- gH2 --->
         !       rL                    rH
         !  <--- gL1 --->        <--- gH1 --->

         ! calculate fraction of r-rectangle that is filled by g-parallelogram
         ! multiplied by its baseline (rH-rL).

      IMPLICIT NONE
      CHARACTER     Fehler*(*)
      REAL*8        rL, rH, gL1, gL2, gH1, gH2,
     *              gLi, gHi, gLo, gHo, wL, wH

      Alpha_Base = 0
      IF     (rH .le.rL ) THEN
         Fehler = 'FracMS: rH < rL'
         RETURN
      ELSEIF (gH1.lt.gL1) THEN
         Fehler = 'FracMS: gH1 < gL1'
         RETURN
      ELSEIF (gH2.lt.gL2) THEN
         Fehler = 'FracMS: gH2 < gL2'
         RETURN
         ENDIF

      gLi = dmax1 (gL1, gL2) ! inner bound
      gHi = dmin1 (gH1, gH2)
      gLo = dmin1 (gL1, gL2) ! outer bound
      gHo = dmax1 (gH1, gH2)

      wL  = dmax1 (gLi, rL)
      wH  = dmin1 (gHi, rH)

      IF     (rH.lt.gLo) THEN
         Alpha_Base = 0
      ELSEIF (rL.gt.gHo) THEN
         Alpha_Base = 0
      ELSEIF (gLi.gt.gHi) THEN
         Print *, 'gLi gHi : ', gLi, gHi
         Fehler = 'Alpha_Base: gLi>gHi not forseen'
         RETURN
      ELSEIF (rH.lt.gLi) THEN
         Alpha_Base = (rH-gLo)**2 / (gLi-gLo) / 2
         IF (rL.gt.gLo) THEN
            Alpha_Base = Alpha_Base - (rL-gLo)**2 / (gLi-gLo) / 2
            ENDIF
      ELSEIF (rL.gt.gHi) THEN
         Alpha_Base = (gHo-rL)**2 / (gHo-gHi) / 2
         IF (rH.lt.gHo) THEN
            Alpha_Base = Alpha_Base - (gHo-rH)**2 / (gHo-gHi) / 2
            ENDIF
      ELSE
         Alpha_Base = wH - wL ! rectangular part
         IF (gLi.gt.rL) THEN
            IF (gLo.lt.rL) THEN ! g-border cuts r-border
               Alpha_Base = Alpha_Base + (gLi-rL)**2 / (gLi-gLo) / 2
            ELSE
               Alpha_Base = Alpha_Base + (gLi-gLo) / 2
               ENDIF
            ENDIF
         IF (gHi.lt.rH) THEN
            IF (gHo.gt.rH) THEN ! g-border cuts r-border
               Alpha_Base = Alpha_Base + (rH-gHi)**2 / (gHo-gHi) / 2
            ELSE
               Alpha_Base = Alpha_Base + (gHo-gHi) / 2
               ENDIF
            ENDIF
         ENDIF

c      Print '(a,3(2f10.3,2x),f8.4)', 'X ', rL, rH, gL1, gL2, gH1, gH2,
c     *                          Alpha_Base/(rH-rL)

      END ! Alpha_Base

C  ====================================================================
      SUBROUTINE RandAng (S, C)
C  ====================================================================
               ! J.von Neumann, NBS Appl.Math.Ser.12,36(1951)

         ! Calculate sine S and cosine C of a random angle.
         ! (without calling the SIN, COS, or just sqrt function)

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      real, external :: rand !Artem add for rand()

  100 CONTINUE
C      X  = G05DAF (-1.d0, 1.d0)  ! arbitrary x from (-1,1) !Artem: Replace with rand function from slatec/fnlib/rand.f
      ran = rand(0.0)
      X = -1.d0 + (1.d0 + 1.d0) * dble(ran)
C      Y  = G05CAF (dummy)        ! arbitrary y from (0, 1) !Artem: Replace with rand function from slatec/fnlib/rand.f
      ran = rand(0.0)
      Y = dble(ran)
      R2 = X*X + Y*Y             ! corresponding r^2
      IF (R2.ge.1.  .or.         ! if (x,y) outside the unit sphere
     *    R2.lt.1.E-6) GOTO 100  !   or too close to the origin

                                 ! now same probability for all angles phi/2

      C = (X*X - Y*Y) / R2       ! cos phi = cos^2 phi/2 - sin^2 phi/2
      S =   2*X*Y     / R2       ! sin phi = 2 sin phi/2 cos phi/2

      END ! RandAng

C  ====================================================================
      SUBROUTINE TruncGauss (E, EPO, EWD, Emin, Emax)
C  ====================================================================
         ! Choose a random value E from a Gaussian (EPO,EWD)
         ! with truncation (Emin,Emax).

      IMPLICIT REAL*8   (a-h,o-p,r-z)

 100  CONTINUE
      E = G05DDF_local (EPO, EWD) !Artem: Replace with a self-made subroutine from lnag_local.f: G05DDF (EPO, EWD)

         ! I did not check that EWD is defined the same way as before
      IF (E.le.Emin .or. E.ge.Emax) GOTO 100 ! try it again

      END ! TruncGauss
