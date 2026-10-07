C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i80   :     read raw data / neutron scattering
C
C  ====================================================================

C     Contents :
C        (1) NRSE-Saclay : [->store]
C               RRawNRSE, DCorrNRSE
C        (2) NSE/ IN11 :
C               RRawIN11
C        (3) Backscattering/ elastic scans :
C               RRawIN10e
C        (4) ILL inelastic data format :
C               Medir
C        (5) Backscattering/ full spectra (old SQW) :
C               RRawBS
C        (6) Time-of-flight (old INZ) :
C               RRawTOF
C        (7) HFBS data files
C               RRT_In_Hfbs
C               RRT_In_Hfbso
C        (8) FANS data files
C               RRT_In_Fans
C        (9) BT2 NIST data files + PSI DMC data
C               RRT_In_BT2
C        (10) IN3 Read In still preliminary thus not called
C               RRT_In_IN3

C  ====================================================================
C  i8 (2) NSE / IN11
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE RRawIN11 (Fehler)
C  --------------------------------------------------------------------
            ! JWu 6may95
         ! Read full output from IN11 (march95 format)

      IMPLICIT REAL *8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      PARAMETER    (MLin=1500, MBB=30)
      CHARACTER*40  RFNam, aux
      DIMENSION     X(MC), R(7), Y1(MC), D1(MC), FIP(MBB,2), CBB(MBB)
      CHARACTER     Fehler*(*)
      CHARACTER*80  Lin(MLin)

      DATA          qFirst / .true. /

      CALL OlfCreate (j1, K1, '&nodef', '&olddef', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCnuP (j1, 'x', 't', 'psec', Fehler)
      CALL OlfCnuP (j1, 'y', 'Echo/Aver', 'Cts', Fehler)
      CALL OlfCnuP (j1, 'z1', 'T', 'K', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      wavel = rAskDMu ('Wavelength [A]', wavel, 1.d0, 3.d1)
      E0 = E_of_l(wavel)
      CALL rOlfP (j1, 'E0', 'meV', E0, Fehler)

C  Geometry :
      IF (qFirst) THEN
         DO iBB = 1, MBB
            CALL rAskPair (
     *           'Field integral prefactors (Oe cm) for coil # '//
     *           cl2(iBB), FIP(iBB,1), FIP(iBB,2), 0.d0, 0.d0)
            CBB(iBB) = 0
            ENDDO
      ELSE
         Print *, 'Using field integral prefactors as before'
         ENDIF

      CBB(29) = 16.7 ! Polarizer
      CBB(30) = 16.7 ! Analyzer

C  Loop over input files :
      DO K = 1, MK

         CALL FrageC ('Raw data file [quit] ?', RFNam)
         IF (RFNam.eq.' ') GOTO 9
         CALL OpenFile (37, RFNam, 'i11', 'a', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL ReadFile (37, Lin, MLin, nLin, '&eof', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         z = rAskMu ('Temperature', 0.d0, 1.d4)

C  Decode log :
         il = 0
         dwn = 0
         upp = 0
 21      CONTINUE
            cdet = 0
            cmon = 0
            npts = 0
 22         CONTINUE
            iL = iL + 1
            IF (iL.gt.nLin) THEN
               Fehler = 'RR IN11/ did not find C-lines'
               RETURN
               ENDIF
            IF (Lin(iL)(1:2).eq.'C ') THEN
               CALL FindR (Lin(il), 7, nR, R, .false.)
               IF (nR.lt.3) THEN
                  Fehler = 'RR IN11/ not 3 numbers in B-line '//cl4(il)
                  RETURN
                  ENDIF
               cdet = cdet + R(3)
               cmon = cmon + R(2)
               npts = npts + 1
               GOTO 22
               ENDIF
            IF (Lin(iL)(1:2).ne.'B ') GOTO 22
            IF (npts.ge.1) THEN
               IF     (dwn.eq.0) THEN
                  dwn = cdet / cmon
                  ddwn= dsqrt (cdet/cmon**2 + (cdet**2/cmon**3))
               ELSEIF (upp.eq.0) THEN
                  upp = cdet / cmon
                  dupp= dsqrt (cdet/cmon**2 + (cdet**2/cmon**3))
               ELSE
                  GOTO 29
                  ENDIF
               ENDIF
            GOTO 21
 29      CONTINUE

         X (1) = 0
         Y1(1) = (upp - dwn) / (upp + dwn) ! Monitor k"urzt sich raus
         D1(1) = 4*(upp**2 * ddwn**2 + dwn**2 * dupp**2) / (upp+dwn)**4
         n = 1

         DO il = iL, nLin
            IF     (Lin(il)(1:2).eq.'B ') THEN
               CALL FindR (Lin(il), 7, nR, R, .false.)
               IF (nR.lt.4) THEN
                  Fehler = 'RR IN11/ not 4 numbers in B-line '//cl4(il)
                  RETURN
                  ENDIF
               iBB = idnint (R(1)) ! # of coil
               IF (qioutside(iBB,1,MBB-2)) THEN
                  Fehler = 'RR IN11/ bad coil # in B-line '//cl4(il)
                  RETURN
                  ENDIF
               CBB(iBB) = R(2)     ! set value

            ELSEIF (Lin(il)(1:7).eq.'EC aver') THEN
               CALL FindR (Lin(il), 7, nR, R, .false.)
               IF (nR.lt.4) THEN
                  Fehler = 'RR IN11/ not 4 numbers in EC_av-line '//
     *                     cl4(il)
                  RETURN
                  ENDIF
               IF (R(1).gt.0) THEN
                  n = n+1 ! new data point
                  Y1(n) = R(3) / R(1)
                  D1(n) = dsqrt ((R(4)/R(1))**2 +
     *                          (R(3)*R(2)/R(1)/R(1))**2)

c                  ! Log :
c                  Print '(i2,1x,6f8.2,2x,2f6.0)',
c    *            n, CBB(2), CBB(17), CBB(9), CBB(30), CBB(18),
c    *             CBB(5), R(1), R(3)

                  IF     (CBB(4).gt.0 .and. CBB(7).eq.0) THEN
                     iKG = 1 ! small spin echo
                  ELSEIF (CBB(4).eq.0 .and. CBB(7).gt.0) THEN
                     iKG = 2 ! large spin echo
                  ELSE
                     Fehler = 'RR IN11/ small or large ?'
                     RETURN
                     ENDIF

                  ! Calculate spin echo time :
                  X(n) = 0
                  DO iBB = 1, MBB
                     X(n) = X(n) + FIP(iBB,iKG) * CBB(iBB)
                     ENDDO
                  X(n) = 1.863d-4 * X(n) * wavel**3

               ELSE
                  Fehler = 'RR IN11/ no counts in EC_av-line '//
     *             cl4(il)
                  RETURN
                  ENDIF

               ENDIF
            ENDDO ! loop over log file lines

C  Save Results :
         K1 = K1 + 1
         CALL OlfPutSpe (j1, K1, 1, z, n, X, Y1, D1, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         qFirst = .false.
         ENDDO
C  End loop over input files.

 9    CONTINUE
      CALL OlfClos (j1, K1, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! RRawIN11

C  ====================================================================
C  i8 (3) Backscattering/ elastic scans
C  ====================================================================

C  --------------------------------------------------------------
      SUBROUTINE RRawIN10e (Fehler)
C  --------------------------------------------------------------
         ! read elastic scans from IN10 dump file

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      PARAMETER    (MLin=MC)
      CHARACTER*40  RFNam, aux
      DIMENSION     X (MC), Y (MC), D (MC), RR(11), ZZ(8)
      CHARACTER     Fehler*(*)
      CHARACTER*80  ein, Lin(MLin)
      DATA          qFirst / .true. /

      IF (qFirst) THEN
         Print *, ' Setup for IN10 elastic scan:'
         wavel = rAsk ('Wavelength (A) ?')
         nK = iAskMu ('Number of detectors used ?', 1, 8)
         DO K = 1, nK
            ZZ(K) = rAsk ('Enter detector angle')
            ENDDO
         qFirst = .false.
      ELSE
         Print *,
     * 'Re-using previous setup (for other setup, restart IDA)'
         ENDIF

 1    CONTINUE

      nRF = iAskDMu ('Number of raw data file (0=qui)', nRF+1, 0, 9999)
      IF (nRF.eq.0) RETURN

      RFNam   = cv4(nRF)
      CALL OpenFile (37, RFNam, 'col', 'a', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL ReadFile (37, Lin, MLin, nLin, '&eof', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Output file :
      CALL OlfCreate (jout, Kout, 'r'//cv4(nRF), '&olddef', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCnuP (jout, 'x', 'T', 'K', Fehler)
      CALL OlfCnuP (jout, 'y', 'Cts', ' ', Fehler)
      CALL OlfCnuP (jout, 'z1', '2th', ' ', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO K = 1, nK
         n = 0
         DO i = 5, nLin ! FOUR lines HEADER assumed <<<<<<<< GEAENDERT !!!!
            n = n + 1
            ein = Lin(i)
            CALL FindR (ein, 11, nR, RR, .false.)
            X(n) = RR(2)
            Y(n) = RR(K+3) / RR(3)
            D(n) = dsqrt0(RR(K+3)) / RR(3)
            ENDDO

         CALL OlfPutSpe (jout, K, 1, ZZ(K), n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      CALL OlfClos (jout, nK, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      GOTO 1

      END ! RRawIN10e

C  ====================================================================
C  i8 (4) Inelastic data access/ 1995 ILL format
C  ====================================================================

C  --------------------------------------------------------------------
      LOGICAL FUNCTION Medir (Instru, iCycle, irun, iBlock, exp, date,
     *                        medpar,pp2,p1,p2,iAdd)
C  --------------------------------------------------------------------
               ! From INZ; copy for IN10sgi may95
         !  Read input as produced from the ILL program SPECTRA
         !  VMS compatibility still maintained, but obsolete 15dec95

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      PARAMETER (MCin=1024, MM=156, ML=450, MN=128)

C  Variables for use in Medir :
      REAL*4    p1(128), p2(128)
      INTEGER   medpar(156), iD(MCin), iM(MCin), iAdd(MCin)
      BYTE      pp2(512)
      CHARACTER exp*10, date*18

      CHARACTER*80  Fehler, LongTit, DirRaw, FileRaw, line
      Character*(*) Instru
      CHARACTER     cl3*3, cv3*3, cv4*4, cv6*6

      Fehler = '&ff'
      Medir = .false.

C  Open the input files :
      IF (irun.lt.1 .or. irun.gt.99999) THEN
         Print *, ' ILL_In/ Bad number of run : ',irun
         RETURN
         ENDIF

      IF     (iCycle.eq.-1) THEN
         CALL ExeML ('\p dir-raw-n-lst', DirRaw)
      ELSEIF (iCycle.eq. 0) THEN
         CALL ExeML ('\p dir-raw-n-new', DirRaw)
      ELSEIF (iCycle.eq. 1) THEN
         CALL ExeML ('\p dir-raw-n-dac', DirRaw)
      ELSE
         CALL ExeML ('\p dir-raw-n-old', DirRaw)
         ENDIF
      IF (DirRaw(1:4).eq.'&err') THEN
         Print *, 'This cycle is not accessible with the present setup'
         RETURN
         ENDIF
      FileRaw = DirRaw

      IF     (Instru.eq.'IN10') THEN
         CALL ReplaceT (FileRaw, '&inst', 'in10')
      ELSEIF (Instru.eq.'IN16') THEN
         CALL ReplaceT (FileRaw, '&inst', 'in16')
      ELSEIF (Instru.eq.'IN13') THEN
         CALL ReplaceT (FileRaw, '&inst', 'in13')
      ELSE
         Print *, 'Invalid instrument : ', Instru
         RETURN
         ENDIF

      IF (iCycle.lt.1000) THEN
         CALL ReplaceT (FileRaw, '&cycle', cv3(iCycle))
      ELSE
         CALL ReplaceT (FileRaw, '&cycle', cv4(iCycle))
         ENDIF

      CALL Append (FileRaw, cv6(irun))

      CALL OpenFile (11, FileRaw, '&noext', 'a', Fehler)
      IF (Fehler.ne.'&ff') THEN
         Print *, 'ILL_In/ Could not open file ', FileRaw
         RETURN
         ENDIF

C  Read header blocks :
      Read(11, '(/i8)') jrun

      Read(11, '(/////16(10(i8)/))', end=11) (medpar(i), i=1,156)
 11   CONTINUE

      nK = medpar(1)

      Read (11,'(a)') line
      IF (line(1:3).ne.'AAA') THEN
         Print *, line
         CALL Absturz ('ILL_In', 'Block 4/ expected AAAAAA')
         ENDIF
      Read(11, '(/a80,6(/))') LongTit

      Read (11,'(a)') line
      IF (line(1:3).ne.'FFF') THEN
         Print *, line
         CALL Absturz ('ILL_In', 'Par1/ expected FFFFFF')
         ENDIF
      Read(11, '(/40(5(2x,d14.8)/))', end=12) (p1(i), i=1,128)
 12   CONTINUE

      Read (11,'(a)') line
      IF (line(1:3).ne.'FFF') THEN
         Print *, line
         CALL Absturz ('ILL_In', 'Par2/ expected FFFFFF')
         ENDIF
      Read(11, '(/40(5(2x,d14.8)/))', end=13) (p2(i), i=1,128)
 13   CONTINUE

      IF (Instru.eq.'IN6')
     *     Read (11,'(53(/),x)') ! a new, nonsensical integer field on IN6

C  Read data blocks :
      DO K=1, iBlock
         Read (11,'(a)') line
         IF (line(1:3).ne.'SSS') THEN
            Print *, 'spectrum '//cl3(K)
            Print *, line
            CALL Absturz ('ILL_In', 'expected SSSSSSS')
            ENDIF
         Read (11,'(//i8)') iBlo
         IF (iBlo.gt.MCin) CALL Absturz (
     *       'ILL_In', '#channels > MCin in spectrum '//cl3(K))
         Read (11, '(10i8)')  (iAdd(i), i=1,iBlo)
         ENDDO

      Close(11)

C  Checks :
      IF    (jrun.ne.irun) THEN
         Print *, '#run : given/read ', irun, jrun
         CALL Absturz ('ILL_In',' Bad number of run')
c      ELSEIF (nC.ne.511) THEN
c         CALL Absturz ('ILL_In', '#channels<>511)')
         ENDIF

      Medir = .true. ! successfully executed

      END ! Medir

C  ====================================================================
C  i8 (5) Backscattering (old SQW)
C  ====================================================================

C   General information :
C        Completely rewritten March-May 1991 by M.Wendel, J.Wuttke
C
C   Table of major modifications :
C
C           27. 7.07  FK: included choice S(q,w) S(2th,w) IN16/IN10
C           27. 2.96  IN16 version (from O.Randl's version of SQW)
C            3. 5.95  Implementaion for SG-Workstation; read in ASCII-File
C           26. 9.91  Debye integral
C           23. 9.91  Log file, outer loop, defaults
C           14. 8.91  Dialogue revised, some repairs
C           21. 5.91  Absorption correction added

C   Structure (location of modules given if <> sqw0) :
C        Main
C           GetData        ! sum numors, set common x-scale
C              ReadILL(2)  ! read 1 numor, determine energies
C                 Medir(SpeLib)
C           GetSpectrum    ! take spectra from array, group them

C   Implementation :
C        Link sqw0-2, i1, l0-l6, (ComData:SpeLib/Li), NAGLIB2/L

C  --------------------------------------------------------------------
      SUBROUTINE RRawBS (Instru, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      PARAMETER (MCin=1024, MKin=35)
      DIMENSION  SDY(MCin), SDX(MCin), SDYE(MCin), SDAngleA(MKin),
     *           EshiftA(MKin), SDYA(MCin,MKin), SDYEA(MCin,MKin)

      REAL*8     SDAngle

      CHARACTER*(*) Fehler, Instru
      CHARACTER*40  Title, TitOut
      DATA          TitOut / ' ' /, iz /3/

C  --------------------------------------------------------------------
C        Preset
C  --------------------------------------------------------------------

      IF     (Instru.eq.'IN10' .or. Instru.eq.'IN16') THEN
         iErgPow = -6 ! ueV
         IF (qAskD('output S(2th,w)',merge(1, 0, qiz))) THEN !Artem qiz -> merge(1, 0, qiz)
            iz   =  1 !2th
         ELSE
            iz      =  3   ! q0
            ENDIF
      ELSEIF (Instru.eq.'IN13') THEN
         iErgPow = -3 ! meV
         iz      =  3 ! q0
      ELSE
         Fehler = 'RawBS/ no BS instrument'
         RETURN
         ENDIF

C  --------------------------------------------------------------------
C        Read raw data
C  --------------------------------------------------------------------

 1    CONTINUE ! new run
      CALL GetData(SDX1, SDXstep, SDYA, SDYEA, SDAngleA, EshiftA,
     *             SDTemp, nSD, nKSD, Ef, Title, Instru)
      IF (nKSD.le.0) RETURN
      nK = nKSD
      IF (nSD.le.0) THEN
         Fehler = 'RRawBS/ read empty spectra ???'
         RETURN
         ENDIF

      IF (SDTemp.le.0.) THEN
            ! for IN10, the sample temperature is not recorded,
            ! GetData returns Temp=0.
         SDTemp = rAskMu ('Sample temperature [K] ?',0.d0,8000.d0)
         ENDIF

      CALL OlfCreate (jout, Kout, '&nodef', '&olddef', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      IF     (iErgPow.eq.-3) THEN
         CALL OlfCnuP (jout, 'x', 'w', 'meV', Fehler)
         CALL OlfCnuP (jout, 'y', 'S(q,w)', 'meV-1', Fehler)
      ELSEIF (iErgPow.eq.-6) THEN
         CALL OlfCnuP (jout, 'x', 'w', 'ueV', Fehler)
         IF (iz.eq.1) THEN
            CALL OlfCnuP (jout, 'y', 'S(2th,w)', 'ueV-1', Fehler)
         ELSEIF (iz.eq.3) THEN
            CALL OlfCnuP (jout, 'y', 'S(q,w)', 'ueV-1', Fehler)
         ENDIF
      ELSE
         Fehler = 'iErgPow o.o.r.'
         ENDIF
      IF (Fehler.ne.'&ff') RETURN

      IF     (iz.eq.1) THEN
         CALL OlfCnuP (jout, 'z1', '2th', ' ', Fehler)
      ELSEIF (iz.eq.3) THEN
         CALL OlfCnuP (jout, 'z1', 'q', 'A-1', Fehler) ! in guter Naeherung
      ELSE
         Fehler = 'iz o.o.r.'
         ENDIF
      IF (Fehler.ne.'&ff') RETURN

      CALL iOlfP (jout, '?det-bal-sym',   0, Fehler)
      CALL iOlfP (jout, '@sam-erg-gain', -1, Fehler)

      CALL rOlfP (jout, 'E0', 'meV', Ef,     Fehler)
      CALL rOlfP (jout, 'T',  'K',   SDTemp, Fehler)

      CALL OlfComAdd (jout, ' ',
     *     'raw data '//Instru(1:lenU(Instru))//': '//Title, Fehler)

C  --------------------------------------------------------------------
C        Loop over grouped spectra
C  --------------------------------------------------------------------

      DO K = 1, nK

C  - Get sample data :
         DO i = 1,nSD
            SDY(i)  = SDYA(i,K)
            SDYE(i) = SDYEA(i,K)
            ENDDO
         SDAngle = SDAngleA(K)

C  - Write out :
         DO i = 1,nSD
            SDX(i) = SDX1+(i-1)*SDXstep
            ENDDO
         IF (iErgPow.eq.-6) THEN ! meV -> ueV
            DO i = 1, nSD
               SDX (i) = SDX (i) * 1000
               SDY (i) = SDY (i) * 1000
               SDYE(i) = SDYE(i) * 1000
               ENDDO
            ENDIF

         IF     (iz.eq.1) THEN
            z = SDAngle
         ELSEIF (iz.eq.3) THEN
            z = 2 * dsqrt (Ef / 2.0723) * dsind (SDAngle/2) ! Q = kf-ki
            ENDIF

         CALL OlfPutSpe (jout, K, 1, z, nSD, SDX, SDY, SDYE, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! K

      CALL OlfClos (jout, nK, Fehler)
      Print *

      GOTO 1

      END ! SQW_main

C  --------------------------------------------------------------------
      SUBROUTINE GetData (X1,Xstep,YA,YEA,AngleA,EshiftA,Temp,
     *                    n,nK,Ef,Title,Instru)
C  --------------------------------------------------------------------
         ! reads all data of all angles in a 2-dim array
         ! sets energy scale

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      PARAMETER (MCin=1024, MKin=35)
      DIMENSION YM(MCin),YMa(MCin),Xo(MCin),Xu(MCin),
     *          Xao(MCin),Xau(MCin)
      DIMENSION YA(MCin,MKin),YEA(MCin,MKin),YAa(MCin,MKin)
      DIMENSION AngleA(MKin),AngleAa(MKin),EshiftA(MKin),EshiftAa(MKin)
      DIMENSION Coder(30),Codera(30)
      CHARACTER*(*) Instru,Title
      CHARACTER exp*10, date*18, Titlea*80 !Artem add Titlea*80
c     EQUIVALENCE (YAE,YAa) !! to save storage capacity if necessary

C     external character functions :
      CHARACTER     ch1*1, cr2*2, cl2*2, cr3*3, cr4*4, cl6*6

C  Read first numor :
      CALL ReadILL (X1,Xstep,YA,YM,AngleA,Coder,EshiftA,Temp,
     *                   n,nK,Ef,Title,Instru)
      IF (nK.eq.0) RETURN ! no sample data are given
      IF (n.le.0) THEN
         CALL Gong (9)
         Print *, 'GetData/ got n<=0 from ReadILL'
         RETURN
         ENDIF

C  Add further numors as long as wanted :
      DO WHILE (1.ne.0)

C  - Read in :
200      CONTINUE
         Print *, 'Add further numors ?'
         CALL ReadILL (X1a,Xstepa,YAa,YMa,AngleAa,Codera,EshiftAa,
     *                  Tempa,na,nKa,Efa,Titlea,Instru)
         IF (nKa.eq.0) GOTO 300

C  - Check for consistency :
         IF (nka.ne.nk) THEN
            CALL Gong (3)
            Print *, 'WARNING/ This numor ignored.'
            Print *, 'Number of detectors inconsistent.'
            GOTO 200
            ENDIF
         CALL Append (Title, Titlea)
C  - - Check consistency of all angles and energy shifts :
         qAngleIncons = .false.
         qEshiftIncons = .false.
         DO k = 1,nK
            IF ((dabs(AngleA(k)-AngleAa(k))).gt.(0.1)) THEN
               qAngleIncons = .true.
               ENDIF
            IF ((dabs(EshiftA(k)-EshiftAa(k))).gt.(0.1*Xstep)) THEN
               qEshiftIncons = .true.
               ENDIF
            ENDDO
         IF (qAngleIncons) THEN
            Print *, 'WARNING/ Angle inconsistent.'
            ENDIF
         IF (qEshiftIncons) THEN
            Print *,
     * 'WARNING/ Energy shift (analyser offset) inconsistent.'
            ENDIF
C  - - Check consistence of temperature :
         IF ((dabs(Temp-Tempa)).gt.(1+0.01*Temp)) THEN
            Print *, 'WARNING/ Temperature inconsistent.'
            ENDIF
C  - - Check consistence of neutron energy :
         IF ((dabs(Ef-Efa)).gt.(0.01*Ef)) THEN
            Print *, 'WARNING/ Neutron energy inconsistent.'
            ENDIF
C  - - Check consistence of coder :
         DO l = 1,30
            IF ((dabs(Coder(l)-Codera(l))).gt.(0.01*Coder(l))) THEN
               Print *, 'WARNING/ Coder '//cr2(l)//' inconsistent.'
               ENDIF
            ENDDO

C  - Add the new numor (the non-trivial part of this subroutine) :

         ! calculate number of channels added at the beginning :
         IF (X1a-Xstepa/2.lt.x1-Xstep/2) THEN
            dXu = (X1-Xstep/2)-(X1a-Xstepa/2)
            nu = idint(dXu/Xstep)
         ELSE
            nu = 0
            ENDIF

         ! calculate number of channels added at the end :
         IF (X1a+Xstepa*(na+0.5).gt.X1+Xstep*(n+0.5)) THEN
            dXo = (X1a+Xstepa*(na+0.5))-(X1+Xstep*(n+0.5))
            no = idint(dXo/Xstep)
         ELSE
            no = 0
            ENDIF

        ! warns if channel number of extended energy scale to large

         IF (n+nu+no.gt.MCin) THEN
            Print *, 'WARNING/ number of channels exceeds '//cl6(MCin)
            Print '(2(a,g10.4))', ' old : X1 = ', X1,  ', dX = ',Xstep
            Print '(2(a,g10.4))', ' new : X1 = ', X1a, ', dX = ',Xstepa
            Print '(3(a,i4))', ' #old =', n,
     *         ', #new(<1) = ', nu, ', #new(>N) = ', no
            Print *, 'last numor not added'
            GOTO 210
            ENDIF

        ! shifts the data nu channels

         IF (nu.gt.0) THEN
            DO i = n,1,-1
               YM(i+nu) = YM(i)
               DO k = 1,nK
                  YA(i+nu,k) = YA(i,k)
                  YEA(i+nu,k) = YEA(i,k)
                  ENDDO
               ENDDO
            DO i = 1,nu
               YM(i) = 0.
               DO k = 1,nK
                  YA(i,k) = 0.
                  YEA(i,k) = 0.
                  ENDDO
               ENDDO
            ENDIF ! shift of data

         ! calculates new channel number and new X1

         n = n+nu+no
         X1 = X1-nu*Xstep

         ! add other channels

         DO i = 1,n                ! calculates boundaries of each
            Xo(i) = X1+(i-0.5)*Xstep ! channel of old energy scale
            Xu(i) = Xo(i)-Xstep
            ENDDO
         DO j = 1,na+1                ! calculates boundaries of each
            Xao(j) = X1a+(j-0.5)*Xstepa ! channel of new energy scale
            Xau(j) = Xao(j)-Xstepa
            ENDDO

         j = 1                           ! j: channel counter of new scale
         DO i = 1,n                      ! i:   "         "    " old   "
            DO WHILE (Xao(j).lt.Xu(i))   ! while no overlap of channels i
               j = j+1                   ! of old scale and j of new scale
               IF (j.gt.n) GOTO 210
               ENDDO
            DO WHILE (Xau(j).lt.Xo(i))   !while overlap of channels i and j
               part = (dmin1(Xao(j),Xo(i))-dmax1(Xau(j),Xu(i)))/Xstepa
                                         ! part: relative extend of overlap
cdeb         print '(a,2i4,5g10.2)', 'i j p YM.. YA(8)..', i, j, part,
cdeb     *      YM(i), YMa(j), YA(i,8), YAa(j,8)
               YM(i) = YM(i)+part*YMa(j) ! add monitor data
               DO k = 1,nK               ! k: angle counter
                  YA(i,k) = YA(i,k)+part*YAa(j,k) ! add count rates
                  ENDDO
               j = j+1
               IF (j.gt.na) GOTO 220   ! if all new channels added
               ENDDO
            IF (j.gt.1) j = j-1
220         CONTINUE
            ENDDO

C  - The new numor is added, end of loop :
210      CONTINUE
         ENDDO
300   CONTINUE ! exit from loop : no further numors are given


C  Normalize detector to monitor counts :

      DO k = 1,nK
         DO i = n,1,-1
            ymr = YM(i)
            ydr = YA(i,k)
            YA(i,k) = dquot0(ydr,ymr)
            YEA(i,k) = YA(i,k)*dsqrt(dquot0(1.d0,ydr)+dquot0(1.d0,ymr))
            ENDDO ! i
         ENDDO   ! K

      END ! GetData

C  --------------------------------------------------------------------
      SUBROUTINE ReadILL (X1,Xstep,YA,YM,AngleA,Coder,EshiftA,Temp,
     *                    n,nK,Ef,Title,Instru)
C  --------------------------------------------------------------------

      ! reads all data of all angles in a 2-dim array
      ! calculates error of all values
      ! calculates X1 and Xstep

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      PARAMETER (MCin=1024, MKin=35)

C  Global variables :
      DIMENSION YA(MCin,MKin), YM(MCin)
      DIMENSION AngleA(MKin), Coder(30), EshiftA(MKin)
      CHARACTER Instru*4, Title*40

C  Local variables :
      REAL*8    X(MCin)
      CHARACTER udb*1, udb0*1

C  Variables for use in Medir :
      REAL*4    p1(128), p2(128)
      INTEGER   medpar(156), iD(MCin), iM(MCin), iAdd(MCin), iE(MCin)
      BYTE      pp2(512)
      CHARACTER exp*10, date*18, TitleIn*40

      EQUIVALENCE (pp2(1),TitleIn)

C  External character functions :
      CHARACTER cr3*3, cr6*6, cl8*8

C  External function :
      LOGICAL   Medir

C  Defaults for dialogue :
      DATA      udb /'b'/, udb0 /' '/

C  Ask for numor identification :

10    CONTINUE
      numor = iAskDu('Numor (or RETURN) ? ',0)
      IF (numor.le.0) THEN
         nK = 0
         RETURN
         ENDIF
      iCycle = iAskDMu ('Cycle (1=bus, 0=present, -1=last)',
     *                  iCycle,0,9999)
      Title = cl8(Numor)

C  read header blocks :
      iBlock = 1
      IF (.not.Medir(Instru,iCycle,numor,iBlock,exp,date,medpar,pp2,
     *               p1,p2,iD)) THEN
         CALL Gong (2)
         GOTO 10
         ENDIF

C  read number of blocks and number of channels :
      ! default : no double blocks (corr. JWu 5may93)
      Koffset = 0
      qboth   = .false.
      udb0 = ' '
      IF (Instru.eq.'IN10') THEN      ! IN10
         qboth = .false.
         nK = int(p1(20))
         nKMon = nK+1
         IF (p1(22).eq.0) THEN         ! if IN10 Doppler
            qMono = .false.
            n = medpar(2)-1
         ELSEIF (p1(22).eq.12 .or. p1(22).eq.13) THEN   ! if IN10 Monochromator
                    ! Que'est-ce que c'est que p1(22) ? Ajoute le cas =13 ?
            qMono = .true.
            n = medpar(2)
         ELSE
            Print *, ' ReadILL/ IN10 A or B ?? p1(22) = ', p1(22)
            nK = 0
            RETURN
            ENDIF
         ENDIF                          ! end IN10

      IF (Instru.eq.'IN16') THEN      ! IN16
         qboth = .false.
         nK = int(p1(8))
         nKMon = nK+1
         qMono = .false.
         n = int(p1(7)-1.)
         ENDIF                          ! end IN16

      IF (Instru.eq.'IN13') THEN      ! IN13
         qMono = .true.
         n = int(p1(7))
         nK = int(p1(8))-2
         nKMon = int(p1(8)) - 1 ! (apr00) monitor now in block 36
         nBlock = medpar(1)
         Print *, 'ReadILL/ IN13/ ', n, nK, nBlock
         IF (nBlock.eq.2*nK+4) THEN
            qDouble = .true.              ! if up- and down-scan
         ELSE IF (nBlock.eq.nK+2) THEN
            qDouble = .false.
         ELSE
            Print *, 'FATAL ERROR :'
            Print *, 'Problems with data format.'
            Print *, 'Number of blocks fits neither to up - scan'//
     *              'nor to up and down - scan.'
            nK = 0
            RETURN
            ENDIF
         IF (qDouble) THEN   ! sets variables for up and down scan
132         CONTINUE
            CALL FrageCD
     *         ('Take up (u), down (d) or both (b) scans', udb, udb)
            CALL Minuskeln(udb)
            IF (udb.eq.'u') THEN
            ELSE IF (udb.eq.'d') THEN
               Koffset = nK+2
            ELSE IF (udb.eq.'b') THEN
               qboth = .true.
            ELSE
               CALL Gong(3)
               GOTO 132
               ENDIF
            udb0 = udb ! for doc
            ENDIF ! qDouble
         ENDIF ! IN13

C  read sample temperature, ...

      IF (Instru.eq.'IN10') THEN
         Temp = 0.
      ELSE IF (Instru.eq.'IN13') THEN
         Temp = dble(p1(9))
      ELSE IF (Instru.eq.'IN16') THEN
         Temp = dble(p1(10))
         ENDIF

C  read coders, angles

      IF ((Instru.eq.'IN10').or. (Instru.eq.'IN13')) THEN

         DO i = 1,30                   !  IN10 and IN13
            Coder(i) = dble(p1(50+i))
            ENDDO
         DO k = 1,nK                   !  IN10 and IN13
            AngleA(k) = dble(p2(k))
            ENDDO

      ELSE IF (Instru.eq.'IN16') THEN

         DO i = 1,30
            Coder(i) = dble(p1(30+i))
            ENDDO

         IF (nK.gt.20) THEN

            DO k = 1, 20
               AngleA(k) = dble(p2(k))+(p1(65)-1.)*5.0d0 ! multidetector
               ENDDO

            DO k = 21,nK
               AngleA(k) = dble(p2(k)) !  single detectors
               ENDDO

          ELSEIF (nK.eq.20) then

             DO k = 1,20
                AngleA(k) = dble(p2(k))+(p1(65)-1.)*5.0d0 ! multidetector
                ENDDO

          ELSE
             DO k = 1,nK
                AngleA(k) = dble(p2(k+20)) ! single detectors
                ENDDO

             ENDIF ! K > or = or < 20

         ENDIF ! Instru

C  determine lattice constant

      IF (Instru.eq.'IN10') THEN
         dlattice = dble(p1(84))
         ENDIF
      IF (Instru.eq.'IN16') THEN
         dlattice = dble(p1(80))
         ENDIF
      IF (Instru.eq.'IN13') THEN
         deltaTemp = dble(p1(91))-28.  ! p1(91):analyser temperature
          ! temperature difference to reference temperature 28 C
         dlattice = dble(p1(82)) *
     *             ( 1. + dble(p1(89))*deltaTemp +
     *                    dble(p1(90))*(deltaTemp**2)/2 )
          ! p1(82) : lattice constant of analyzer-crystal at 28 C
          ! p1(89),p1(90) material constants for thermal expansion
         ENDIF

C  calculate energy Ef, energy shift EshiftA(k)

      Ef = 81.8055/((2.*dlattice)**2.)
        ! dlattice : lattice constant of analyzer crystal
        ! lattice-const-->wavelength-->energy of detected neutrons
      DO k = 1,nK
         EshiftA(k) = Ef*(dtand(dble(p2(50+k))/2.)**2.)
          ! p2(50+k) : analyser offsets (angle)
         ENDDO

C  calculate X1, Xstep

      IF (Instru.eq.'IN10') THEN     ! if IN10
         IF (.not.qMono) THEN         ! if IN10 Doppler
            EM = 81.8055 / (2*dble(p1(82)))**2
             ! p1(82) : lattice constant of monochromator crystal
             ! lattice-const = .5 wavelength-->energy of incident neutr.
            ED = 19.656769d8 * dble(p1(2)) * 4.135701327d-12 /
     *                                          (2*dble(p1(82)))
             ! 19.65..d8 : length of Doppler-drive [A] (d8:cm-->A)
             ! p1(2) : Doppler-frequency [s^-1]
             ! 4.1357..(-12) : h [meV*s]
             ! calculates maximal energy difference by doppler effect
             !INFO: Xstep and X1 calculated after analysing monitor rates
         ELSE IF (qMono) THEN   ! if IN10 Monochromator
            iBlock = nK+2
            IF (.not.Medir(Instru,iCycle,numor,iBlock,exp,date,medpar,
     *         pp2,p1,p2,iD)) THEN
               Print *, 'FATAL ERROR :'
               Print *, 'Medir cannot read the energy block.'
               nK = 0
               RETURN
               ENDIF
            fact = dble(p1(23))      ! converts energy scale to microeV
             ! Exclude channels with X=0. :

c  das ist noch gar nicht schoen;
c  hier muessen schon nL und nR eingefuehrt werden.
            j = 1
            DO WHILE (iD(j).eq.0)   ! first energy value not equal 0
               j = j+1
               ENDDO
            m = n
            DO WHILE (iD(m).eq.0)   ! last energy value not equal 0
               m = m-1
               ENDDO
c  auch schlecht :
            Xone = dble(iD(j))   / (fact*1000)   ! first written energy value
            Xtwo = dble(iD(j+1)) / (fact*1000)   ! second   "       "     "
            Xstep = Xtwo - Xone
            X1 = Xone - (j-1) * Xstep
            DO i = j+1,m  ! control if energy steps are equidistant
               X(i)   = dble (iD(i)  ) / (fact*1000)
               X(i-1) = dble (iD(i-1)) / (fact*1000)
                  ! Energies are stored in nano(!)eV as integer values.
                  ! 1000 converts neV to ueV
                  ! fact is usually 1000 and converts ueV to meV
               IF (dabs(X(i)-X(i-1)-Xstep).gt.(0.001*Xstep)) THEN
c                  Print *, 'FATAL ERROR :'
                  CALL Gong (5)
                  Print *, 'Energy scale is not equidistant: i = ',i
c                  STOP
                  Print *, 'Notreparatur :'
                  X1    = rAsk ('E(#1) [ueV]') / 1000
                  Xstep = rAsk ('dE    [ueV]') / 1000
                  GOTO 209 ! exit loop
                  ENDIF
               ENDDO
 209        CONTINUE
         ELSE     ! if no inelastic scan
            Print *, 'ERROR :'
            Print *, 'this numor is no inelastic scan : p1(22) = ',
     *               p1(22)
            GOTO 10
            ENDIF                      ! IN10 Monochromator
         ENDIF                        ! IN10

      IF (Instru.eq.'IN16') THEN

            EM = 81.8055 / (2*dble(p1(70)))**2
             ! p1(70) : lattice constant of monochromator crystal
             ! lattice-const = .5 wavelength-->energy of incident neutr.

            ED = 15.939917d8 * dble(p1(3)) * 4.135701327d-12 /
     *                                          (2*dble(p1(70)))
             ! 15.93..d8 max.speed (A/s) for Doppler-freq. = 1Hz
             ! p1(3) : Doppler-frequency [s^-1]
             ! 4.1357..(-12) : h [meV*s]
             ! p1(70) lattice parameter monochromator (A)
             ! calculates maximal energy difference by doppler effect
             !INFO: Xstep and X1 calculated after analysing monitor rates
         ENDIF ! IN16

      IF (Instru.eq.'IN13') THEN
c         X0 = dble(p1(2)-p1(3))/1000
c          ! p1(2) = centre of energy range (microeV)
c          ! p1(3) = half width of energy range
c         Xstep = dble(p1(12))/1000
c         X1 = X0+Xstep/2
         ! apr00: energies are in block 38
         iBlock = nKMon + 2
         IF (.not.Medir(Instru,iCycle,numor,iBlock,exp,date,medpar,pp2,
     *                  p1,p2,iE)) THEN
            Print *, 'problem with Medir reading energies '
            GOTO 10
            ENDIF
         X1    = iE(1) / 1d5
         Xstep = (iE(2)-iE(1)) / 1d5
         Xn    = iE(n) / 1d5
         IF (Xn.ne.X1 + (n-1)*Xstep) THEN
            Print *, 'energy scale not equidistant ????'
            RETURN
            ENDIF
         ENDIF

C  read monitor counts

      iBlock = nKMon+Koffset
      IF (.not.Medir(Instru,iCycle,numor,iBlock,exp,date,medpar,pp2,
     *               p1,p2,iM)) THEN
         Print *, 'problems with Medir reading monitor-data '
         GOTO 10
         ENDIF

C  add up and down monitor rates (IN13) if wanted

      IF (qboth) THEN
         iBlock = nKMon+(nK+2)
         IF (.not.Medir(Instru,iCycle,numor,iBlock,exp,date,medpar,pp2,
     *                 p1,p2,iAdd)) THEN
            Print *, 'problems with Medir reading monitor-data of'
            Print *, 'second scan of up and down scan'
            GOTO 10
            ENDIF
         DO i = 1,n
            iM(i) = iM(i)+iAdd(i)
            ENDDO
         ENDIF

C  determine significant channels nL,...,nR (revised JWu 14aug91) :
      nL  = 1
      nR  = n
      rel = .1 ! arbitrary cutoff condition

 320  CONTINUE
      ! mean monitor counts :
      YmMean = 0
      DO i = nL,nR
         YmMean = YmMean+dble(iM(i))
         ENDDO
      nFull  = nR - nL + 1
      YmMean = YmMean / nFull

      ! throw away empty or poor channels on the left :
      DO WHILE (iM(nL).le.(rel*YmMean))
         nL = nL+1
         ENDDO
      ! the same on the right :
      DO WHILE (iM(nR).le.(rel*YmMean))
         nR = nR-1
         ENDDO

      Print *, 'ReadILL/ nR = ', n, ' -> ',nR,'; nL = ', 1, ' -> ',nL
      n = nR - nL + 1                   ! number of significant channels
      IF (n.lt.nFull) GOTO 320          ! repeat the operation with new YmMean

C  adjust energy scale :
      IF (((Instru.eq.'IN10').and.(.not.qMono)) .or.
     *     (Instru.eq.'IN16')) THEN                   ! doppler scan
         Xstep = 2*ED/n
          ! energy difference between each channel
         X1 = -Xstep*(n-1)/2.+Ef-EM ! corrected 29aug94 (error found by ORa)
          ! minimal energy (compared to incoming neutron energy)
      ELSE
         X1 = X1+(nL-1)*Xstep                         ! monochromator scan
         ENDIF

C  save monitor data as YM(i) :
      DO i = 1,n
         YM(i) = dble(iM(i+(nL-1)))
         ENDDO

C  read all spectra and save as YA(i,k)

      DO k = 1,nK
         iBlock = k+Koffset            ! Koffset for IN13 up and down scan
         IF (.not.Medir(Instru,iCycle,numor,iBlock,exp,date,medpar,pp2,
     *                 p1,p2,iD)) THEN
            Print *, 'problems with Medir reading spectrum',k
            GOTO 10
            ENDIF
         IF (qboth) THEN              ! if adding up and down scan IN13
            iBlock = k+nK+2
            IF (.not.Medir(Instru,iCycle,numor,iBlock,exp,date,medpar,
     *                   pp2,p1,p2,iAdd)) THEN
               Print *, 'problems with Medir reading second spectrum',k
               GOTO 10
               ENDIF
            DO i = 1,n
               iD(i) = iD(i)+iAdd(i)
               ENDDO
            ENDIF
         DO i = n,1,-1                 ! +(nL-1) : begins with first
            YA(i,k) = dble(iD(i+(nL-1))) ! significant channel
            ENDDO ! i
         ENDDO   ! K

      END ! ReadILL

C  ====================================================================
C  i8 (6) Time-of-flight (old INX/INZ)
C  ====================================================================

C   Structure :
C      RRawTOF
C        RRT_InSub
C          RRT_InSum
C            RRT_Input
C              RRT_In_ILL
C              RRT_In_Frm
C              RRT_In_Mib
C              RRT_In_Sil
C              RRT_In_Fcs
C              RRT_In_Foc
C              RRT_In_Foco
C              RRT_EPeak_Fcs
C        RRT_Check
C        RRT_EPeak
C        RRT_Energ
C        RRT_DetEff
C        RRT_Norm
C        RRT_Sort
C        RRT_Add

C  History :
C     S.Busch apr08 consider wavelength-dependence of the monitor
C     S.Busch dec07 corrected error of vanadium-normalization,
C                   changed defaults to suit TOFTOF
C     F.Kargl mar06 implement T. Unruh's and F. Juranyis FOCUS ReadIn
C     F.Kargl jul05 changes for TOFTOF at FRM2 (frame overlap etc.)
C     F.Kargl jun05 add TOFTOF read in routine
C     F.Kargl jul02 changed vandium,sample, other file elastic channel
C                   treatment.
C     F.Kargl oct01 #frame overlap#
C     A.Meyer aug01 # of det. IN6
C     A.Meyer may01 RRT_In_Foc
C     A.Meyer may01 RRT_In_Mib new data format
C     A.Meyer feb98 RRT_In_Fcs
C     JWu 5dec95 2dim arrays eliminated
C     JWu 15aug95 as subroutine RRawTOF
C     A.Meyer + JWu feb95 : RRT_In_Sil, spec sum+del under output
C     JWu jan93 : installation Saclay, entruempelt, REAL*8
C     JWu jun92 : frame overlap, output, input, entruempelt
C     JWu sep91 : DetEff for IN6
C     JWu jul91 : RRT_In_Mib version 2, DetEff for Mibemol, Text in EppCon
C     J.Wuttke aug91 : version INZ without menus
C     M.Bee '91 : Mibemol version 1
C     Program INX by Fr. Rieutord 1990

C  --------------------------------------------------------------------
      SUBROUTINE RRawTOF (Inst, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      PARAMETER   (MNuSu=100)

      DIMENSION    AngleS(MK), AngleV(MK),AngleF(MK), FPath(MK),
     *             DetE1(MK), DetE2(MK),
     *             qLisK(MK), qKnoAnal(MK),
     *             XepS(MK), Xepfcs(MK), YepS(MK), XepV(MK), YepV(MK),
     *             XepF(MK), YepF(MK)
      CHARACTER*1  BadS(MK), BadV(MK), BadF(MK)

      CHARACTER*(*) Inst, Fehler
      CHARACTER*40  FilOut, Title
      CHARACTER*80  ListSD, ListSC, ListVD, ListVC, ListFD, ListFC,
     *              LKnoAnal, LisKdel

      DIMENSION     JNumorSD(MNuSu), JNumorSC(MNuSu),
     *              JNumorVD(MNuSu), JNumorVC(MNuSu),
     *              JNumorFD(MNuSu), JNumorFC(MNuSu)

C  --------------------------------------------------------------------
C  Default values :
C  --------------------------------------------------------------------

      DATA     nLoop /0/, ListSC, ListVC, ListFC /3*'-'/,
     *         qVan /.true./, qOutE /.true./, q0fromS /.false./,
     *         q0fromF /.false./, iFmod /2/, EFmax /1.d3/,
     *         EFminRel /-.99d0/, qMon /.true./, qDetEff /.true./,
     *         iDel /1/, qAngEqAdd / .true. /

C  --------------------------------------------------------------------
C  Checks and Initializations :
C  --------------------------------------------------------------------

      IF (.not. (Inst.eq.'IN4' .or. Inst.eq.'IN5' .or.
     *           Inst.eq.'IN6' .or. Inst.eq.'MIB' .or.
     *           Inst.eq.'SIL' .or. Inst.eq.'FCS' .or.
     *           Inst.eq.'FOCUS' .or. Inst.eq.'NEAT' .or.
     *           Inst.eq.'FOCUSO'.or. Inst.eq.'DCS'.or.Inst.eq.'TOF'))
     * THEN
         Fehler =
     *     'Instruments are IN4,IN5,IN6,TOF,MIB,SIL,FCS,FOC,DCS,NEAT -
     *      do not know '//Inst
         GOTO 99
         ENDIF

C  Open Log-File :
      CALL OpenFile (35, 'inz', 'log', 'e', Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

C  --------------------------------------------------------------------
C  Begin dialogue :
C  --------------------------------------------------------------------

C      IF (Inst(1:2).eq.'IN')
      iCycle = iAskD (
     * 'Data from instrument(-1), current cycle(0), archive(>0)',
     * iCycle)
      iNOff = iAskD ('Run number offset', iNOff)

 110  CONTINUE ! outer loop
      nLoop = nLoop + 1

      ListSD = ' '
      CALL GetJList ('Sample run numbers', ListSD,
     *               MNuSu, nNumorSD, JNumorSD, 0, 9999999)
      IF (nNumorSD.le.0) GOTO 99
      DO i = 1, nNumorSD
         JNumorSD(i) = JNumorSD(i) + iNOff
         ENDDO

      IF (nLoop.gt.1) THEN
         qTheSame = qAskD('All the rest as before',intq(qTheSame))
         IF (qTheSame) GOTO 200
      ELSE
         qTheSame = .false.
         ENDIF

      IF (Inst(1:3).eq.'DCS') THEN
         iCycle = iAskDMu('Year of the experiment?',iCycle,1997,2010)
         ENDIF

      CALL GetJList ('Background run numbers', ListSC,
     *               MNuSu, nNumorSC, JNumorSC, 0, 9999999)
      DO i = 1, nNumorSC
         JNumorSC(i) = JNumorSC(i) + iNOff
         ENDDO

      CALL GetJList ('Vanadium run numbers', ListVD,
     *               MNuSu, nNumorVD, JNumorVD, 0, 9999999)
      DO i = 1, nNumorVD
         JNumorVD(i) = JNumorVD(i) + iNOff
         ENDDO

      IF (nNumorVD.ne.0) THEN
         CALL GetJList ('Vanadium background run numbers', ListVC,
     *                   MNuSu, nNumorVC, JNumorVC, 0, 9999999)
         DO i = 1, nNumorVC
            JNumorVC(i) = JNumorVC(i) + iNOff
            ENDDO
         qMon = .true.
      ELSE
         qMon = qAskD ('Normalization to monitor', merge(1, 0, qMon)) !Artem qMon -> merge(1, 0, qMon)
         ENDIF

      Print *
      Print *, 'Sample treatment :'

C      qOutE   = .true. ! qAskD ('Conversion to energy', intq(qOutE))
      qOutE   = qAskD ('Conversion to energy', intq(qOutE))

      IF (nNumorVD.ne.0) THEN
         q0fromS= qAskDi (
     *      'Elastic channel from vanadium(0) or sample/other file(1)',
     *       intq(q0fromS))
         q0fromF = .false.  !FK jul02
      ELSE
         q0fromS = .true.
         ENDIF
      IF (q0fromS) THEN
         q0fromF= qAskDi (
     * 'Elastic channel from sample(0) or other file(1)',
     * intq(qofromF))
      ENDIF

      IF (q0fromF) THEN
         CALL GetJList ('Other file number', ListFD,
     *               MNuSu, nNumorFD, JNumorFD, 0, 9999999)
         DO i = 1, nNumorFD
         JNumorFD(i) = JNumorFD(i) + iNOff
         ENDDO
      ENDIF

      IF (nNumorFD.gt.0 .and. q0fromF) THEN !changed FK jul02
         CALL GetJList ('Other file background number', ListFC,
     *                 MNuSu, nNUmorFC, JNumorFC, 0, 9999999)
         DO i = 1, nNumorFC
            JNumorFC(i) = JNumorFC(i) + iNOff
            ENDDO
      ENDIF

      IF (qOutE) THEN
         iSEG    = -1 ! neutron energy gain is positive
         qDetEff = qAskD ('Detector efficiency correction',
     *                    intq(qDetEff))
      ELSE
         iSEG    = -1 ! not used
         qDetEff = .false.
         ENDIF

      IF (nNumorVD.ne.0) THEN
         Print *
         Print *, 'Vanadium treatment :'
         qVan = qAskD ('Take DWF for vanadium', intq(qVan))
         IF (.not.qVan) THEN
            u2xV = rAskDMu
     * ('Mean square displacement <(u_x)^2> [A^2]', u2xV, 0.d0, 10.d0)
            ENDIF

         ENDIF

      ! for use in EppCon (handle frame-overlap) :
      Print *
      Print *, 'Handling of frame overlap :'
      Print *, '   (0) do nothing'
C      Print *, '   (1) set energies to limit the frame,'
      Print *, '   (2) proper handling of frame overlap'
      Print *, '   (4) proper handling special TOFTOF'
C      Print *, '   (2) explicitly enter channel numbers,'
C      Print *, '   (3) determine from global minimum of scattering'


      iFmod = iAskDMu ('Choose option', iFmod, 0, 4)
C      IF    (iFmod.eq.1) THEN
C         EFmax    = rAskDMu (
C     * 'Maximal neutron energy gain in frame (meV)', EFmax, 0.d0, 1.d6)
C         EFminRel = - rAskDLu (
C     * 'Maximal energy loss in frame (unit=E0)', -EFminRel, 0.d0, 1.d0)
      IF (iFmod.eq.0) THEN
         NooF = 0
         NooD = 0
      ELSEIF (iFmod.eq.2.or.iFmod.eq.4) THEN
         EFmax    = rAskDMu (
     * 'Maximal neutron energy gain in frame (meV)', EFmax, 0.d0, 1.d6)
         EFminRel = - rAskDLu (
     * 'Maximal energy loss in frame (unit=E0)', -EFminRel, 0.d0, 1.d0)
C      ELSEIF (iFmod.eq.2) THEN
C         NooF = iAskDMu (
C     * 'Number of fast neutron channels to move around', NooF, 0,
C      * MCin)
C         NooD = iAskDMu (
C     * 'Number of intermediate channels to delete', NooD, 0, MCin)
      ELSE ! case iFmod=3
         ! nothing to ask
         ENDIF

      qAngEqAdd = qAskD (
     *   'Sum spectra at equal angle (before background subtraction)',
     *   intq(qAngEqAdd))

      IF (iFmod.eq.2.or.iFmod.eq.4) THEN  ! (.or. iFmod.eq.1) THEN
         IF (EWmax.le.0. .or. EWmax.gt.EFmax) EWmax = EFmax ! default
      ELSE
         EFmax = 1.d19 ! no limitation
         IF (EWmax.le.0. .or. EWmax.gt.EFmax) EWmax = 1.d3 ! default
         ENDIF

      EWmax = EFmax ! rAskDMu ('Maximal energy on output', EWmax, 0.d0, EFmax)

      CALL GetNList ('Exclude spectra from analysis', LKnoAnal,
     *               qKnoAnal, MK)

 200  CONTINUE

C  --------------------------------------------------------------------
C  Get vanadium data / prepare normalization :
C  --------------------------------------------------------------------

      Write (35,'(/3a)') 'going to read data from instrument ', Inst,
     *     ', cycle '//cl6(iCycle)

      IF (qTheSame) THEN ! this option new 28aug91
         Print *
         Print *, 'vanadium as in last run'
         Write (35, '(/a)') 'Vanadium as above'

      ELSEIF (nNumorVD.le.0) THEN
         WRITE (35,'(/a)') 'No vanadium given - no normalization.'

      ELSE
         Print *
         Print *, 'getting vanadium runs :'

         Write (35,'(/a)') 'going to read vanadium run(s):'
         CALL RRT_InSub (qMon, Inst, iCycle, qKnoAnal, nNumorVD,
     *               JNumorVD, nNumorVC, JNumorVC, j, VTemp, EelastV,
     *               WaveLV, CwidthV, PeriodV, qAngEqAdd,
     *               nKV, AngleV, FPath, DetE1, DetE2, BadV, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99

         Print *, 'processing the vanadium :'

         IF (Inst.eq.'FCS') THEN
            CALL RRT_EPeak (j, nKV, AngleV, Xepfcs, YepV, BadV, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99

            CALL RRT_EPeak_FCS (j, nKV, AngleV, FPath, Xepfcs,
     *        WaveLV, EelastV, CwidthV, PeriodV, Fehler)

            CALL RRT_EPeak (j, nKV, AngleV, XepV, YepV, BadV, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99

         ELSE
            CALL RRT_EPeak (j, nKV, AngleV, XepV, YepV, BadV, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
            ENDIF
         Write (35,'(/a)') 'analyse vanadium data :'

C  Correct for Debye-Waller factor of vanadium :
         IF (qVan) THEN
         Write (*, '(a,f6.2)') 'Vanadium temperature: ',VTemp
           IF (285.0.gt.VTemp .OR. VTemp.gt.310.0) THEN
             VTemp = rAsk ('Enter Vanadium temperature in K: ')
             ENDIF
            u2xV = u2Debye (VTemp, 359.d0, 50.94d0, 2.5d-5) ! WuL6
            Write (*, '(a,f6.3,a)') ' <u_x^2> = ', u2xV*1000,
     *                '*10-3 A^2'
            ENDIF
         Write (35, '(a,f6.3,a)') 'Van : <u_x^2> = ', u2xV*1000,
     *              '*10-3 A^2'
         pi = 3.14159d0
         DWFmin = 1.d0
         DO K = 1, nKV
            elaQ     = 4*pi/WaveLV * dsind(AngleV(K)/2)
            DWF      = dexp (-u2xV * elaQ**2)
            YepV(K)  = YepV(K)  / DWF
            DWFmin   = dmin1(DWF,DWFmin)
            ENDDO

         Write (*, '(a,f6.4,a,f6.2,a)') 'DWF at highest angle: ',
     *             DWFmin, ' for Vanadium at ' ,VTemp, 'K'

         CALL MemFileDel (j, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99

         ENDIF

      IF (.not.qMon) WRITE (35,'(/a)') 'No normalization to monitor.'

C  ===================================================================
C     Get other file data for elastic peak
C  ===================================================================
      ! Florian Kargl Oct01


      IF (q0fromF) THEN

       Print *
       Print *, 'going to read runs from other file'

       Write (35,*) 'going to read runs from other file'

       CALL RRT_InSub(qMon,Inst,iCycle,qKnoAnal, nNumorFD, JNumorFD,
     *               nNumorFC, JNumorFC, j, FTEMP, EelastF, WaveLF,
     *               CwidthF, PeriodF, qAngEqAdd, nKF, AngleF, FPath,
     *               DetE1, DetE2, BadF, Fehler)
       IF (Fehler.ne.'&ff') GOTO 99


       Print *, 'Processing other file'

       IF (Inst.eq.'FCS') THEN
          CALL RRT_EPeak (j, nKF, AngleF, Xepfcs, YepF, BadF, Fehler)
          IF(Fehler.ne.'&ff') GOTO 99

          CALL RRT_EPeak_FCS (j, nKF, AngleF, FPath, Xepfcs, WaveLF,
     *                        EelastF, CwidthF, PeriodF, Fehler)

          CALL RRT_EPeak (j, nKF, AngleF, XepF, YepF, BadF, Fehler)
          IF(Fehler.ne.'&ff') GOTO 99

       ELSE
          CALL RRT_EPeak (j, nKF, AngleF, XepF, YepF, BadF, Fehler)
          IF(Fehler.ne.'&ff') GOTO 99
       ENDIF

       CALL MemFileDel (j, Fehler)
       IF (Fehler.ne.'&ff') GOTO 99

      ENDIF
C  =================================================================




C  --------------------------------------------------------------------
C  Get sample data / normalize and correct them :
C  --------------------------------------------------------------------

      Print *
      Print *, 'getting sample runs :'

      Write (35,*) 'going to read sample run(s):'
      CALL RRT_InSub (qMon, Inst, iCycle, qKnoAnal,
     *               nNumorSD, JNumorSD, nNumorSC, JNumorSC, j,
     *               STemp, Eelast, WaveL, Cwidth, Period, qAngEqAdd,
     *               nK, AngleS, FPath, DetE1, DetE2, BadS, Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

C     write real parameters if necessary for conversion to energy
      IF (qOutE .eqv. .false.) THEN  !Artem .eq. -> .eqv.
      CALL rOlfP (j, 'FP', 'm', FPath(1), Fehler)    !Artem FPath -> FPath(1)
      CALL rOlfP (j, 'Cw', 'usec', Cwidth, Fehler)
      ENDIF


      IF (nNumorVD.gt.0) THEN
         ! check whether parameters from sample and vanadium are consistent :
         CALL RRT_Check ('sample vs vanadium',
     *                STemp, Eelast, WaveL, Cwidth, Period,
     *                nK, AngleS, FPath, DetE1, DetE2,
     *                VTemp, EelastV, WaveLV, CwidthV, PeriodV,
     *                nKV, AngleV, FPath, DetE1, DetE2,
     *                Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         ENDIF

      Write (35, '(/a)')        'TOF parameters (sample) :'
      Write (35, '(a28,g12.5)') '  wavelength (A) ', WaveL
      Write (35, '(a28,g12.5)') '  energy (meV) ', Eelast
      Write (35, '(a28,g12.5)') '  channel width (usec) ', Cwidth
      Write (35, '(a28,g12.5)') '  period (chs) ', Period

      Print *
      Print *, 'processing the sample :'

      Write (35,'(/a)') 'analyse sample data :'

      IF (Inst.eq.'FCS') THEN
         IF (nNumorVD.ne.0) THEN
            DO K = 1, nK
               XepS(K) = Xepfcs(K)
               ENDDO
            ELSE
            CALL RRT_EPeak (j, nK, AngleS, XepS, YepS, BadS, Fehler)
            ENDIF

         CALL RRT_EPeak_FCS (j, nK, AngleS, FPath, XepS,
     *       WaveL, Eelast, Cwidth, Period, Fehler)

         CALL RRT_EPeak (j, nK, AngleS, XepS, YepS, BadS, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99


      ELSE
         CALL RRT_EPeak (j, nK, AngleS, XepS, YepS, BadS, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         ENDIF

C   ---------------------------------------------------------------------
C     check if data from sample and file are consistent
C   ---------------------------------------------------------------------
      IF (q0fromF) THEN
         Print *
         Print *, 'check if data from file and sample are consistent'

       CALL RRT_Check ('file vs sample',
     *                STEMP, Eelast, WaveL, Cwidth, Period, nK, AngleS,
     *                FPath, DetE1, DetE2, FTEMP, EelastF, WaveLF,
     *                CwidthF, PeriodF, nKF, AngleF, FPath, DetE1,
     *                DetE2, Fehler)
       IF (Fehler.ne.'&ff') GOTO 99
       ENDIF




      IF (q0fromS .and. (.not.q0fromF)) THEN
         Write (35, '(/a)') '  elastic peak from sample'
         Print *, '  elastic peak from sample'
      ELSEIF ((.not.q0fromS) .and. (.not.qofromF)) THEN
         Write (35, '(/a)') '  elastic peak from vanadium'
         Print *, '  elastic peak from vanadium'
      ELSE
         Write (35, '(/a)') '  elastic peak from other file'
         Print *, '  elastic peak from other file'
         ENDIF

      Write (35,*)
      WRITE (35,*) 'Parameters used for normalization :'
      qBad0 = .false.
      qBada = .false.
      qBadx = .false.
      qBadm = .false.
      qBade = .false.
      qBadv = .false.
      qBadg = .false.
      qBadn = .false.
      WRITE (35,*)
     *  '          |     EP Intensity    |    EP Position | Trouble |'
      WRITE (35,*)
     *  ' Spec|Ang.| Sample   |   Vana   | Sample | Vana  | S  |  V |'
      DO K=1, nK
         WRITE (35,435) K, AngleS(K),
     * YepS(K), YepV(K), XepS(K), XepV(K), BadS(K), BadV(K)
         IF     (BadV(K).eq.'0' .or. BadS(K).eq.'0') THEN
            qBad0 = .true.
         ELSEIF (BadV(K).eq.'x' .or. BadS(K).eq.'x') THEN
            qBadx = .true.
         ELSEIF (BadV(K).eq.'a' .or. BadS(K).eq.'a') THEN
            qBada = .true.
         ELSEIF (BadV(K).eq.'-' .or. BadS(K).eq.'-') THEN
            qBadm = .true.
         ELSEIF (BadV(K).eq.'e' .or. BadS(K).eq.'e') THEN
            qBade = .true.
         ELSEIF (BadV(K).eq.'n' .or. BadS(K).eq.'n') THEN
            qBadn = .true.
         ELSEIF (BadV(K).eq.'v' .or. BadS(K).eq.'v') THEN
            qBadv = .true.
         ELSEIF (BadV(K).eq.'g' .or. BadS(K).eq.'g') THEN
            qBadg = .true.
            ENDIF
         ENDDO
435   FORMAT (2X,i3,X,F5.1,X,G10.4,X,G10.4,X,F6.2,2x,G10.4,x,a3,2x,a3)
      Write (35, *)
      ! explain error marks :
      IF (qBadx) Write(35, *)
     * '   trouble x : spectrum excluded a priori'
      IF (qBad0) Write(35, *) '   trouble 0 : empty data block'
      IF (qBadm) Write(35, *) '   trouble - : negative count rates'
      IF (qBada) Write(35, *)
     *        '   trouble a : invalid angle (monitor or special block)'
      IF (qBadv) Write(35, *) '   trouble v : nearly vanishing peak'
      IF (qBadg) Write(35, *) '   trouble g : giant peak'
      IF (qBade) Write(35, *) '   trouble e : excentric peak'
      IF (qBadn) Write(35, *) '   trouble n : normalization impossible'
      qBad = qBad0 .or. qBadm .or. qBada .or. qBadv .or. qBadg .or.
     *       qBade .or. qBadn .or. qBadx

      IF (qOutE) THEN
      Print *
      Print *, 'conversion to energy :'
      Write (35,'(/a)') ' conversion to energy :'
      IF     (iSEG.eq.+1) THEN
         Write (35, '(a)') ' E>0 means sample energy gain'
      ELSEIF (iSEG.eq.-1) THEN
         Write (35, '(a)') ' E>0 means neutron energy gain'
         ENDIF
      ELSE
      Print *
      Print *, 'no conversion to energy'
      ENDIF

      IF (q0fromS .and. (.not. q0fromF)) THEN
         CALL RRT_Energ (qOutE, iSEG, XepS, BadS, AngleS, FPath,
     *                iFmod, NooW, NooD,
     *                EFmax, EWmax, EFminRel,
     *                WaveL, Eelast, Cwidth, Period,
     *                j, nK, Fehler)
      ELSEIF ((.not. q0fromS) .and. (.not.q0fromF)) THEN
         CALL RRT_Energ (qOutE, iSEG, XepV, BadV, AngleV, FPath,
     *                iFmod, NooW, NooD,
     *                EFmax, EWmax, EFminRel,
     *                WaveL, Eelast, Cwidth, Period,
     *                j, nK, Fehler)
       ELSE
         CALL RRT_Energ (qOutE, iSEG, XepF, BadF, AngleF, FPath,
     *                iFmod, NooW, NooD,
     *                EFmax, EWmax, EFminRel,
     *                WaveL, Eelast, Cwidth, Period,
     *                j, nK, Fehler)
         ENDIF

      IF (qDetEff) THEN
         WRITE (35,*) 'detector efficiency correction for ', Inst
         Print *, 'detector efficiency correction'
         CALL RRT_DetEff (j, nK, Eelast, DetE1, DetE2, Fehler)
      ELSE
         WRITE (35,*) 'no detector efficiency correction'
         Print *, 'no detector efficiency correction'
         ENDIF
      IF (Fehler.ne.'&ff') GOTO 99

C  Normalisation to the vanadium
       ! The normalisation to vanadium should be done after the grouping
       ! otherwise one should take an average over Angle instead of a sum

      IF (nNumorVD.eq.0) THEN
         Print *, 'no vanadium - no normalization'
         DO K=1, nK
            BadV(K) = ' '
            ENDDO
      ELSE
C         Print *, 'normalization to vanadium'
         CALL RRT_Norm (j, nK, YepV, BadV, BadS, Fehler)
         IF (Fehler.ne.'&ff') GOTO 99
         CALL OlfComAdd (j, ' ', 'normalized to '//ListVD, Fehler)
         ENDIF

C  --------------------------------------------------------------------
C     Output :
C  --------------------------------------------------------------------

C  default output filename :
      CALL FrageCD ('File name', FilOut, 'r'//cl6(JNumorSD(1)))
      CALL FrageCD ('And title', Title, Title)
      CALL tOlfP (j, 'fil', FilOut, Fehler)
      CALL tOlfP (j, 'tit', Title,  Fehler)
      CALL tOlfP (j, 'sub', Inst , Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

      IF (qTheSame) THEN
         Print *, ' delete some spectra according to option ', iDel
      ELSE
         Print *
         Print *, 'Output :'
         Print *, 'Delete spectra :'
         Print *,
     * '   (0) none                    (1) if some data are bad'
         Print *,
     * '   (2) if sample data are bad  (3) if vanadium data are bad'
         Print *, '   (4) enter list yourself  '
         iDel = iAskDMu ('Option ', iDel, 0, 4)
         ENDIF
      Print *

C  For bad data blocks, set AngleS=0 :
      IF     (iDel.eq.1) THEN
         DO K = 1, nK
            IF (BadS(K).ne.' ' .or. BadV(K).ne.' ') AngleS(K) = 0.
            ENDDO
      ELSEIF (iDel.eq.2) THEN
         DO K = 1, nK
            IF (BadS(K).ne.' ') AngleS(K) = 0.
            ENDDO
      ELSEIF (iDel.eq.3) THEN
         DO K = 1, nK
            IF (BadV(K).ne.' ') AngleS(K) = 0.
            ENDDO
      ELSEIF (iDel.eq.4) THEN
         CALL GetNList ('Delete spectra nos.', LisKdel, qLisK, nK)
         DO K = 1, nK
            IF (qLisK(K)) AngleS(K) = 0.
            ENDDO
         ENDIF

      Write (35,'(/2a)') 'Output: ', FilOut
      Write (35,'(2a/)') 'Title:  ', Title

      CALL RRT_Sort (j, nK, AngleS, STemp, qOutE, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Outermost loop: next run ?
      Print *
      Print *, 'Next data set:'
      GOTO 110

C  Errors :
 99   CONTINUE
      Close(35)

      END ! RRawTOF

C  --------------------------------------------------------------------
      SUBROUTINE RRT_InSub (qMon, Inst, iCycle, qKnoAnal,
     *                nNumor, JNumor, nNumor1, JNumor1, j,
     *                Temp, Eelast, WaveL, Cwidth, Period, qAngEqAdd,
     *                nK, Angle, FPath, DetE1, DetE2, Status, Fehler)
C  --------------------------------------------------------------------
         ! read in, subtract empty can

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     JNumor(*), JNumor1(*)
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*),
     *              Angle1(MK), FPath1(MK), DetE11(MK), DetE21(MK)
      CHARACTER*1   Status(*)
      CHARACTER*(*) Inst, Fehler
      CHARACTER*80  LisK       !Artem add LisK
      LOGICAL       qKnoAnal(*)

      DATA          LisK /'-'/

C  Read first file:
      CALL RRT_InSum (qMon, Inst, iCycle, nNumor, JNumor, j,
     *                Temp, Eelast, WaveL, Cwidth, Period, qAngEqAdd,
     *                nK, Angle, FPath, DetE1, DetE2, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO K = 1, nK
         IF (Angle(K).le.0) THEN
            Status(K) = 'a'
         ELSEIF (qKnoAnal(K)) THEN
            Status(K) = 'x'
         ELSE
            Status(K) = ' '
            ENDIF
         ENDDO

C  Read second file ?
      IF (nNumor1.le.0) RETURN

      Write (35, '(a)') 'subtract background:'
      CALL RRT_InSum (qMon, Inst, iCycle, nNumor1, JNumor1, j1,
     *             Temp1, Eelast1, WaveL1, Cwidth1, Period1, qAngEqAdd,
     *             nK1, Angle1, FPath1, DetE11, DetE21, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL RRT_Check ('Subtract background/ ',
     *                Temp, Eelast, WaveL, Cwidth, Period,
     *                nK, Angle, FPath, DetE1, DetE2,
     *                Temp1, Eelast1, WaveL1, Cwidth1, Period1,
     *                nK1, Angle1, FPath1, DetE11, DetE21, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO K = 1, nK
         CALL OlfGetXYD (j,  K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfGetXYD (j1, K, n1, X1, Y1, D1, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (n.ne.n1) THEN
            Fehler = 'RRT_InSub/ Inconsistent length of spectra '//
     *               cl4(K)
            RETURN
            ENDIF
         DO i = 1, n
            Y(i) = Y(i) - Y1(i)
            D(i) = dsqrt( D(i)**2 + D1(i)**2 )
            ENDDO
         CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      CALL MemFileDel (j1, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      END ! RRT_InSub

C  --------------------------------------------------------------------
      SUBROUTINE RRT_InSum (qMon, Inst, iCycle, nNumor, JNumor, j,
     *                  Temp, Eelast, WaveL, Cwidth, Period, qAngEqAdd,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
         ! sum numors
         ! normalize to monitor, set error

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     JNumor(*)
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*),
     *              Angle1(MK), FPath1(MK), DetE11(MK), DetE21(MK)
      CHARACTER*(*) Inst, Fehler

         Write (35,*) nNumor, ' run(s) to be read:'
         Write (35,'(20i5)') (JNumor(iNu),iNu=1, nNumor)

C  First File:
      CALL RRT_Input (Inst, iCycle, JNumor(1), j,
     *                rMon, Temp, Eelast, WaveL, Cwidth, Period,
     *                nK, Angle, FPath, DetE1, DetE2, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Add more files ?
      DO iNu = 2, nNumor

         CALL RRT_Input (Inst, iCycle, JNumor(iNu), j1,
     *               rMon1, Temp1, Eelast1, WaveL1, Cwidth1, Period1,
     *               nK1, Angle1, FPath1, DetE11, DetE21, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL RRT_Check ('Sum numors/ ',
     *                   Temp, Eelast, WaveL, Cwidth, Period,
     *                   nK, Angle, FPath, DetE1, DetE2,
     *                   Temp1, Eelast1, WaveL1, Cwidth1, Period1,
     *                   nK1, Angle1, FPath1, DetE11, DetE21, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         rMon = rMon + rMon1

         DO K = 1, nK
            CALL OlfGetXYD (j,  K, n, X, Y, D, Fehler)
            CALL OlfGetXYD (j1, K, n1, X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (n.ne.n1) THEN
               Fehler = 'RRT_InSum/ Inconsistent length of spectra '//
     *                  cl4(K)
               RETURN
               ENDIF
            DO i = 1, n
               Y(i) = Y(i) + Y1(i)
               ENDDO
            CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         CALL MemFileDel (j1, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! sum numors

C  Correct the wavelength-dependent sensitivity of the monitor
C  sbusch 2008-APR
C      Ei = rOlfGG (j, 'E0', 'meV', Fehler)
C      IF (Fehler.ne.'&ff') RETURN
C      rMon = rMon / dsqrt (25.305/Ei) ! 25.3meV = 2200m/sec



C  Save summed monitor:
      CALL rOlfP (j, 'cts[mon]', ' ', rMon, Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Add at-mass to real par
      Call rOlfP (j, 'at-mass', 'amu', 1.d0, Fehler)

C  Normalize to monitor and set error :
      IF (qMon .and. rMon.le.0.d0) THEN
         Print *, 'WARNING/ monitor = 0'
         rMon = rAsk ('Enter monitor')
         ENDIF
      DO K = 1, nK
         CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (qMon) THEN
            DO i = 1, n
               D(i) = dsqrt(Y(i)) / rMon
               Y(i) = Y(i) / rMon
               ENDDO
         ELSE
            DO i = 1, n
               D(i) = dsqrt(Y(i))
               ENDDO
            ENDIF
         CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      IF (qAngEqAdd) CALL RRT_Add (j, nK, Angle, Fehler)

      END ! RRT_InSum

C  --------------------------------------------------------------------
      SUBROUTINE RRT_Input (Inst, iCycle, iNumor, j,
     *                      rMon, Temp, Eelast, WaveL, Cwidth, Period,
     *                      nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
         ! call institute specific input routines

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*)
      CHARACTER*(*) Inst, Fehler
      CHARACTER*80  DirRaw, LongTit, LongDat

C  Open on-line file :
      CALL OlfCreate (j, Kout, '&noask', '&noask', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Get pathname from setup file :
      IF     (iCycle.eq.-1) THEN
         CALL ExeML ('\p dir-raw-n-dac', DirRaw)
      ELSEIF (iCycle.eq. 0) THEN
         CALL ExeML ('\p dir-raw-n-new', DirRaw)
      ELSE
         CALL ExeML ('\p dir-raw-n-old', DirRaw)
         ENDIF

C  Call specific input routine:
      IF     (Inst(1:2).eq.'IN') THEN
         CALL RRT_In_ILL (Inst, DirRaw, iCycle, iNumor, j, rMon,
     *               Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSEIF (Inst.eq.'TOF') THEN
         CALL RRT_In_Frm (Inst, DirRaw, iCycle, iNumor, j, rMon,
     *               Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSEIF (Inst.eq.'MIB') THEN
         CALL RRT_In_Mib (DirRaw, iCycle, iNumor, j,
     *               rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSEIF (Inst.eq.'FCS') THEN
        CALL RRT_In_Fcs (DirRaw, iNumor, j,
     *               rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSEIF (Inst.eq.'FOCUS') THEN
        CALL RRT_In_Foc (DirRaw, iNumor, j,
     *               rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSEIF (Inst.eq.'FOCUSO') THEN
         CALL RRT_In_Foco (DirRaw, iNumor, j,
     *               rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSEIF (Inst.eq.'DCS') THEN
         CALL RRT_In_DCS (DirRaw, iNumor, j, iCycle,
     *               rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSEIF (Inst.eq.'NEAT') THEN
        CALL RRT_In_NEAT (DirRaw, iNumor, j,
     *               rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *               nK, Angle, FPath, DetE1, DetE2, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         Eelast  = 81.805 / WaveL**2             ! in meV
      ELSE
         Fehler = 'RRT_Input/ Unknown instrument'
         RETURN
         ENDIF

      Write (35, *) 'read numor ', iNumor, ' - Monitor = ', rMon

      CALL OlfComAdd (j, ' ', Inst(1:lenU(Inst))//': '//
     *                        LongDat, Fehler)
      CALL OlfComAdd (j, ' ', '"'//LongTit(1:lenU(LongTit))//
     *                        '"', Fehler)

      CALL iOlfP (j, '?det-bal-sym',   0, Fehler)
      CALL iOlfP (j, '@sam-erg-gain', -1, Fehler)
      CALL iOlfP (j, 'plot-sy#',  0, Fehler)

      CALL rOlfP (j, 'E0', 'meV', Eelast, Fehler)
      CALL rOlfP (j, 'T',  'K',   Temp,   Fehler)

      END ! RRT_Input

C  --------------------------------------------------------------------
      SUBROUTINE RRT_In_ILL (Inst, DirRaw, iCycle, iNumor, j, rMon,
     *                  Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
            !  Von K.Bauszus eingetippt, 1993.
            !  Modifications after ILL-upstart mrz95
            !  Idol by Helga Schwab, removed 6mai98

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      PARAMETER    (MM=156, ML=450, MN=128, MD=1024)
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*),
     *              Par0(MM), Par1(ML), Par2(MN), iDet(512)

      CHARACTER*80  FileRaw, line
      CHARACTER*(*) Inst, Fehler, DirRaw, LongTit

C  Direct access or not ?
      IF (DirRaw.eq.'&error') THEN
         Fehler = 'Directory not accessible (&error)'
         RETURN
         ENDIF

      FileRaw = DirRaw

      IF (iCycle.lt.1000) THEN
         CALL ReplaceT (FileRaw, '&cycle', cv3(iCycle))
      ELSE
         CALL ReplaceT (FileRaw, '&cycle', cv4(iCycle))
         ENDIF

      IF     (Inst.eq.'IN6') THEN
         CALL ReplaceT (FileRaw, '&inst', 'in6')
      ELSEIF (Inst.eq.'IN5') THEN
         CALL ReplaceT (FileRaw, '&inst', 'in5')
      ELSEIF (Inst.eq.'IN4') THEN
         CALL ReplaceT (FileRaw, '&inst', 'in4')
      ELSE
         Fehler = 'Invalid instrument : '//Inst
         RETURN
         ENDIF

C  Construct file name :
      CALL Append (FileRaw, cv6(iNumor))

C  Open the input files :
      CALL OpenFile (11, FileRaw, '&noext', 'a', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Read header blocks :
C  - Block 1: Numor
      Read (11, '(/a80)', err=91) line
      CALL Fi1I (line, jrun)
      IF    (jrun.ne.iNumor) THEN
         Print *, 'header line "', line(1:lenU(line)), '"'
         Print *, '#run : given/read ', iNumor, jrun
         Fehler = 'inconsistent number of run or invalid data format'
         Close(11)
         RETURN
         ENDIF

C  - Block 2: Short title
C  - Block 3: Some integer parameters (Par0)
      Read (11, '(/////16(10(f8.0)/))', err=92, end=11)
     *      (Par0(i), i=1,156)
 11   CONTINUE

      nK = Par0(1)
      nCh = Par0(2)

C  - Block 4: Long title
      Read(11, '(//60x,a80//////)') LongTit
      Print *, LongTit

C  - Block 5: Real parameters (Par1), mainly scattering angles
      Read (11,'(a)') line
      IF (line(1:3).ne.'FFF') THEN
         Print *, line
         Fehler = 'RRT_In_ILL/ Par2/ expected FFFFFF'
         RETURN
         ENDIF
      Read (11, '(i8)', err=93) nPar1
      Read (11, '(92(5(2x,d14.8)/))', err=93, end=12)
     *      (Par1(i),i=1,nPar1)
 12   CONTINUE

C  - Block 6: Other real parameters (Par2)
      Read (11,'(a)') line
      IF (line(1:3).ne.'FFF') THEN
         Print *, line
         Fehler = 'RRT_In_ILL/ Par2/ expected FFFFFF'
         RETURN
         ENDIF
      Read (11, '(/80(5(2x,d14.8)/))', err=94, end=13)
     *      (Par2(i), i=1,128)
 13   CONTINUE

      IF (Inst.eq.'IN5') THEN
         Read (11,'(27(/),x)') ! a new, nonsensical field (128*real)
         ENDIF ! (vgl. Mail von Andreas Meyer 26aug97)

C  A new field, intended to contain a code for the detector type :
      Read (11,'(a)') line ! dieser Test 1okt96
      IF (line(1:3).ne.'III') THEN
         Print *, line
         Fehler = 'RRT_In_ILL/ Det-Codes/ expected III'
         RETURN
         ENDIF
      Read (11, '(i8)', err=95) niDet
      Read (11, '(53(10(i8)/))', err=96, end=14) (iDet(i), i=1,niDet)
 14   CONTINUE

C  Read data blocks :
      IF (nK.gt.MK) THEN
         CALL Compose2 (Fehler, 'RRT_In_ILL/ found '//cl6(nK),
     *      ' data blocks while MK='//cl6(MK))
         RETURN
         ENDIF

      DO K = 1, nK
         Read (11,'(a)') line
         IF (line(1:3).ne.'SSS') THEN
            Print *, 'spectrum '//cl3(K)
            Print *, line
            Fehler = 'RRT_In_ILL/ expected SSSSSSS'
            RETURN
            ENDIF
         Read (11,'(//i8)', err=97) iBlo
         IF (iBlo.ne.nCh) THEN
            Fehler = 'RRT_In_ILL/ #channels <> nCh in spectrum '//
     *               cl3(K)
            RETURN
            ENDIF
         Read (11, '(10f8.0)', err=98)  (Y(i), i=1,nCh)
         IF (K.eq.1) THEN
            rMon = 0
            DO i = 1, nCh-1
               rMon = rMon + Y(i)
            ENDDO
         ENDIF

         CALL OlfPutSpe (j, K, 1, 0.d0, nCh-1, X, Y, D, Fehler)
            ! ch. nCh is used for no. of block
         IF (Fehler.ne.'&ff') RETURN
         ENDDO
      Close (11)

C  Checks:
      nCnom = NINT(Par2(2))  ! nominelle Zahl der verwendeten Kanaele !Artem: Replace with a standard NINT function: jidnnt(Par2(2))
      IF (nCnom.ne.nCh) THEN
         Fehler = 'RRT_In_ILL/ Decode Par2/ unexpected #chs. = '//
     *            cl6(nCnom)
         RETURN
         ENDIF

C  Par2 enth"alt SetUp-Daten :
         Temp   = Par2(11)          ! Probentemp
         Period = Par2(13)          ! Periode der Messung
         Cwidth = Par2(18)          ! Kanalbreite
         WaveL  = Par2(21)          ! Wellenlaenge
         Eelast = 81.805 / WaveL**2             ! in meV

C  K dependent parameters :

         DO K = 1, nK
            Angle (K) = Par1(K+31) ! die ersten 30 Zahlen sind keine Winkel
            FPath (K) = Par2(27)   ! Abstand Probe Detektor

            DetE1(K) = 1001 ! preset = nonsense
            DetE2(K) = 1002 ! preset = nonsense

            IF     (Inst.eq.'IN4') THEN
               DetE1(K) = - .0887
               DetE2(K) = -5.597
            ELSEIF (Inst.eq.'IN5') THEN
               IF     (iDet(K).eq.0) THEN ! monitor or special block
                  Angle(K) = -1
               ELSEIF (iDet(K).eq.1) THEN ! IN5 type detectors (8bars,D = .9cm)
                  DetE1(K) = - .0887
                  DetE2(K) = -4.07
               ELSEIF (iDet(K).eq.2) THEN ! IN6 type detectors
                  DetE1(K) = -.0565
                  DetE2(K) = -3.284
               ELSEIF (iDet(K).eq.3) THEN ! IN6 type detectors
                  DetE1(K) = -.0565
                  DetE2(K) = -3.284
               ELSEIF (iDet(K).eq.4) THEN ! IN6 type detectors
                  DetE1(K) = -.0565
                  DetE2(K) = -3.284
               ELSEIF (iDet(K).eq.5) THEN ! Multidetector
                  DetE1(K) = -.0565
                  DetE2(K) = -3.284
                  FPath(K) =  Par2(27) - .30
                  ENDIF
            ELSEIF (Inst.eq.'IN6') THEN
               DetE1(K) = -.0565
               DetE2(K) = -3.284
C               IF (K.le.6) Angle(K) = 0 ! monitor blocks
               ENDIF
         ENDDO

C  Errors :
      RETURN
 91   CONTINUE
      CALL Compose2 (Fehler,
     *               'Error while reading 1st lines of "'//FileRaw, '"')
      Close(11)
      RETURN
 92   CONTINUE
      Fehler = 'Read from data file / error 92'
      RETURN
 93   CONTINUE
      Fehler = 'Read from data file / error 93'
      RETURN
 94   CONTINUE
      Fehler = 'Read from data file / error 94'
      RETURN
 95   CONTINUE
      Fehler = 'Could not read new integer block (with detector code)'
      RETURN
 96   CONTINUE
      Fehler = 'Read from data file / error 96'
      RETURN
 97   CONTINUE
      Fehler = 'Read from data file / error 97'
      RETURN
 98   CONTINUE
      Fehler = 'Read from data file / error 98'
      RETURN

      END ! RRT_In_ILL

C  -------------------------------------------------------------------
      SUBROUTINE RRT_In_Frm (Inst, DirRaw, iCycle, iNumor, j, rMon,
     *                  Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
            !  FK12 modified to enable read-in upgraded detector#
            !  FK06 modified for faster detector eff correction
            !  FK05 based on the ILL version of Bauszus with
            !  with modified detector numbers
            !  Von K.Bauszus eingetippt, 1993.
            !  Modifications after ILL-upstart mrz95
            !  Idol by Helga Schwab, removed 6mai98

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      PARAMETER    (MM=156, ML=1025, MN=128, MD=1025)
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*),
     *              Par0(MM), Par1(ML), Par2(MN), iDet(1025)

      CHARACTER*80  FileRaw, line
      CHARACTER*(*) Inst, Fehler, DirRaw, LongTit

C  Direct access or not ?
      IF (DirRaw.eq.'&error') THEN
         Fehler = 'Directory not accessible (&error)'
         RETURN
         ENDIF

      FileRaw = DirRaw

      IF (iCycle.lt.1000) THEN
         CALL ReplaceT (FileRaw, '&cycle', cv3(iCycle))
      ELSE
         CALL ReplaceT (FileRaw, '&cycle', cv4(iCycle))
         ENDIF

      IF     (Inst.eq.'TOF') THEN
         CALL ReplaceT (FileRaw, '&inst', 'tof')
      ELSE
         Fehler = 'Invalid instrument : '//Inst
         RETURN
         ENDIF

C  Construct file name :
      CALL Append (FileRaw, cv6(iNumor))

C  Open the input files :
      CALL OpenFile (11, FileRaw, '&noext', 'a', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Read header blocks :
C  - Block 1: Numor
      Read (11, '(/a80)', err=91) line
      CALL Fi1I (line, jrun)
      IF    (jrun.ne.iNumor) THEN
         Print *, 'header line "', line(1:lenU(line)), '"'
         Print *, '#run : given/read ', iNumor, jrun
         Fehler = 'inconsistent number of run or invalid data format'
         Close(11)
         RETURN
         ENDIF

C  - Block 2: Short title
C  - Block 3: Some integer parameters (Par0)
      Read (11, '(/////16(10(f8.0)/))', err=92, end=11)
     *     (Par0(i), i=1,156)
 11   CONTINUE

      nK = Par0(1)
      nCh = Par0(2)

C  - Block 4: Long title
      Read(11, '(//60x,a80//////)') LongTit
      Print *, LongTit

C  - Block 5: Real parameters (Par1), mainly scattering angles
      Read (11,'(a)') line
      IF (line(1:3).ne.'FFF') THEN
         Print *, line
         Fehler = 'RRT_In_Frm/ Par2/ expected FFFFFF'
         RETURN
         ENDIF
      Read (11, '(i8)', err=93) nPar1
      IF (iNumor.lt.34300) THEN
         Read (11, '(132(5(2x,d14.8)/))', err=93, end=12)(Par1(i),i=1,
     *      nPar1)
      ELSE
         Read (11, '(205(5(2x,d14.8)/))', err=93, end=12)(Par1(i),i=1,
     *      nPar1)
         ENDIF
 12   CONTINUE

C  - Block 6: Other real parameters (Par2)
      Read (11,'(a)') line
      IF (line(1:3).ne.'FFF') THEN
         Print *, line
         Fehler = 'RRT_In_Frm/ Par2/ expected FFFFFF'
         RETURN
         ENDIF
      Read (11, '(/80(5(2x,d14.8)/))', err=94, end=13)
     * (Par2(i), i=1,128)
 13   CONTINUE

C  the true number of channels that is saved in the files
      nChS = int(Par2(13)/Par2(18))

C  A new field, intended to contain a code for the detector type :
      Read (11,'(a)') line ! dieser Test 1okt96
      IF (line(1:3).ne.'III') THEN
         Print *, line
         Fehler = 'RRT_In_Frm/ Det-Codes/ expected III'
         RETURN
         ENDIF
      Read (11, '(i8)', err=95) niDet
      IF (niDet.lt.609) THEN
      Read (11, '(67(10(i8)/))', err=96, end=14) (iDet(i), i=1,niDet)
      ELSE
            IF (niDet.lt.981) THEN
            Read (11, '(98(10(i8)/))', err=96, end=14) (iDet(i),
     *            i=1,niDet)
            ELSE
            Read (11, '(99(10(i8)/))', err=96, end=14) (iDet(i),
     *            i=1,niDet)
            ENDIF
      ENDIF
 14   CONTINUE

C  A new field, intended to contain a code for the detector type :
C      Read (11,'(a)') line ! dieser Test 1okt96
C      IF (line(1:3).ne.'III') THEN
C         Print *, line
C         Fehler = 'RRT_In_Frm/ Det-Codes/ expected III'
C         RETURN
C         ENDIF
C      Read (11, '(i8)', err=95) niDet
C      IF (iNumor.lt.34000) THEN
C      Read (11, '(67(10(i8)/))', err=96, end=14) (iDet(i), i=1,niDet)
C      ELSE
C            Read (11, '(98(10(i8)/))', err=96, end=14) (iDet(i),
C     *            i=1,niDet)
C            ENDIF
C 14   CONTINUE

C  Read data blocks :
      IF (nK.gt.MK) THEN
         CALL Compose2 (Fehler, 'RRT_In_Frm/ found '//cl6(nK),
     *      ' data blocks while MK='//cl6(MK))
         RETURN
         ENDIF

      DO K = 1, nK
         Read (11,'(a)') line
         IF (line(1:3).ne.'SSS') THEN
            Print *, 'spectrum '//cl3(K)
            Print *, line
            Fehler = 'RRT_In_Frm/ expected SSSSSSS'
            RETURN
            ENDIF
         Read (11,'(//i8)', err=97) iBlo
         IF (iBlo.ne.nCh) THEN
            Fehler = 'RRT_In_Frm/ #channels <> nCh in spectrum '//
     *               cl3(K)
            RETURN
            ENDIF
         Read (11, '(10f8.0)', err=98)  (Y(i), i=1,nCh)
         IF (K.eq.1) THEN
            rMon = 0
            DO i = 1, nCh-1
               rMon = rMon + Y(i)
            ENDDO
         ENDIF

         CALL OlfPutSpe (j, K, 1, 0.d0, nChS, X, Y, D, Fehler)
            ! ch. nCh is used for no. of block
         IF (Fehler.ne.'&ff') RETURN
         ENDDO
      Close (11)

C  Checks:
      nCnom = NINT(Par2(2))  ! nominelle Zahl der verwendeten Kanaele !Artem: Replace with a standard NINT function: jidnnt(Par2(2))
      IF (nCnom.ne.nCh) THEN
         Fehler = 'RRT_In_Frm/ Decode Par2/ unexpected #chs. = '//
     *            cl6(nCnom)
         RETURN
         ENDIF

C  Par2 enth"alt SetUp-Daten :
         Temp   = Par2(11)          ! Probentemp
         Period = Par2(13)          ! Periode der Messung
         Cwidth = Par2(18)          ! Kanalbreite
         WaveL  = Par2(21)          ! Wellenlaenge
         Eelast = 81.805 / WaveL**2             ! in meV

C  K dependent parameters :

         DO K = 1, nK
            Angle (K) = Par1(K+31) ! die ersten 30 Zahlen sind keine Winkel
            FPath (K) = Par2(27)   ! Abstand Probe Detektor

            DetE1(K) = 1001 ! preset = nonsense
            DetE2(K) = 1002 ! preset = nonsense

           IF   (Inst.eq.'TOF') THEN
               DetE1(K) = -.31
               DetE2(K) = -9.3518
C      Coefficients modified according to T. Unruh Jun05
C               IF (K.le.3) Angle(K) = 0 ! monitor blocks
               ENDIF
         ENDDO

C  Errors :
      RETURN
 91   CONTINUE
      CALL Compose2 (Fehler,
     *          'Error while reading 1st lines of "'//FileRaw, '"')
      Close(11)
      RETURN
 92   CONTINUE
      Fehler = 'Read from data file / error 92'
      RETURN
 93   CONTINUE
      Fehler = 'Read from data file / error 93'
      RETURN
 94   CONTINUE
      Fehler = 'Read from data file / error 94'
      RETURN
 95   CONTINUE
      Fehler = 'Could not read new integer block (with detector code)'
      RETURN
 96   CONTINUE
      Fehler = 'Read from data file / error 96'
      RETURN
 97   CONTINUE
      Fehler = 'Read from data file / error 97'
      RETURN
 98   CONTINUE
      Fehler = 'Read from data file / error 98'
      RETURN

      END ! RRT_In_Frm

C  --------------------------------------------------------------------
      SUBROUTINE RRT_In_NEAT (DirRaw, iNumor, j, rMon,
     *                  Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      PARAMETER    (MM=156, ML=450, MN=128, MD=512)
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*),
     *              Par0(MM), Par2(MN), iDet(MD)

      CHARACTER*80  FileRaw, line
      CHARACTER*(*) Fehler, DirRaw, LongTit

      INTEGER       j,iang,nbheaderblock
      INTEGER       nbchan_SD,nbspec_SD,int_ty_SD,nbblock_SD
      INTEGER*2     buffer1(256)
      INTEGER*4     buffer3(4096)
      REAL*8        buffer2(128)
      CHARACTER*20  LongDat

C  Direct access or not ?
      IF (DirRaw.eq.'&error') THEN
         Fehler = 'Directory not accessible (&error)'
         RETURN
         ENDIF

C  Construct file name and open input file :
      FileRaw = DirRaw
      CALL Append (FileRaw, cv4(iNumor))
      CALL Append (FileRaw, 'neat.out')
      CALL OpenFile (11, FileRaw, '&noext', 'a', Fehler)
      IF (Fehler.ne.'&ff') THEN
         Fehler = '&ff'
         FileRaw = DirRaw
         CALL Append (FileRaw, cv4(iNumor))
         CALL Append (FileRaw, 'NEAT.OUT')
         CALL OpenFile (11, FileRaw, '&noext', 'a', Fehler)
         ENDIF
      IF (Fehler.ne.'&ff') RETURN

C  Read header blocks :

C.......Read Header-Block No 1

      Read (11, '(a80)', err=91) line
      IF (line(1:3).ne.'RRR') THEN
         Print *, line
         Fehler = 'missed beginning of block 1'
         RETURN
         ENDIF
      DO i = 1, 16
         READ (11, '(a)', Err=91) line ! (buffer1(i),i=1,256)
         ENDDO
      nbheaderblock = 10        ! buffer1(256) ! nach Konversion i2 -> i4 !
cdeb      Print *, 'read block 1'

C.......Read Header-Block No 2

      Read (11, '(a80)') line
      IF (line(1:3).ne.'RRR') THEN
         Print *, line
         Fehler = 'missed beginning of block 2'
         RETURN
         ENDIF
      Read (11, '(a80)') LongTit
      Print '(2a)', ' Run : ', LongTit(1:20) ! , 'started at :', LongDat(1:20)
cdeb      Print *, 'read block 2'

C.......Read Header-Block No 3

      Read (11, '(a80)') line
      IF (line(1:3).ne.'RRR') THEN
         Fehler = 'missed beginning of block 3'
         RETURN
         ENDIF
      Read (11, '(15(8(d10.2)/),8d10.2)') (buffer2(i), i=1,128)

        nbchan_SD=buffer2(56)
        nbspec_SD=buffer2(70)
        int_ty_SD=buffer2(14)
        nbblock_SD=nbchan_SD/256*(1+int_ty_SD)

      Period = buffer2(53)
      Cwidth = buffer2(61) * .125   ! Kanalbreite (usec) < time step 125nsec
      WaveL  = buffer2(27)          ! Wellenlaenge
      FPathSD= buffer2(21) / 1000   ! Abstand Probe - Detektor (mm -> m)

      Period = Period / Cwidth ! (#chs)

      nC = nbchan_SD
      nK = nbspec_SD

      Print *, '#spec, #block, #chan :', nbspec_SD, nbblock_SD,
     *                                   nbchan_SD
      IF (nC.gt.MC) THEN
         Fehler = 'channels/spectrum > MC'
         RETURN
         ENDIF
      IF (nK.gt.MK) THEN
         CALL Compose2 (Fehler, 'RRT_In_NEAT/ found '//cl6(nK),
     *      ' data blocks while MK='//cl6(MK))
         RETURN
         ENDIF
cdeb      Print *, 'read block 3'

C.......Read single detector angles
        iang=0
        DO i=1,16
           iang=iang+1
c          IF (iang .GT. nbspec_SD) GOTO 50
           Angle(i)=buffer2(112+i)
           ENDDO
        Read (11, '(a80)') line
        IF (line(1:3).ne.'RRR') THEN
           Print *, '::', line(1:76)
           Fehler = 'missed beginning of block 4'
           RETURN
           ENDIF
        READ (11, '(15(8(d10.2)/),8d10.2)', err=94)
     *            (buffer2(i),i=1,128)
        DO i=17,144
           iang=iang+1
c          IF (iang .GT. nbspec_SD) GOTO 50
           Angle(i)=buffer2(i-16)
           ENDDO
        Read (11, '(a80)') line
        IF (line(1:3).ne.'RRR') THEN
           Print *, '::', line(1:76)
           Fehler = 'missed beginning of block 5'
           RETURN
           ENDIF
        READ (11, '(15(8(d10.2)/),8d10.2)', err=94)
     *            (buffer2(i),i=1,128)
        DO i=145,272
           iang=iang+1
c           IF (iang .GT. nbspec_SD) GOTO 50
           Angle(i)=buffer2(i-144)
           ENDDO
        Read (11, '(a80)') line
        IF (line(1:3).ne.'RRR') THEN
           Fehler = 'missed beginning of block 6'
           RETURN
           ENDIF
        READ (11, '(15(8(d10.2)/),8d10.2)', err=94)
     *            (buffer2(i),i=1,128)
        DO i=273,400
           iang=iang+1
           IF (iang .GT. nbspec_SD) GOTO 50
           Angle(i)=buffer2(i-272)
           ENDDO
50    CONTINUE
cdeb      Print *, 'read angles, iang = ', iang

C.......Read single detector spectra

      Read (11, '(a80)') line
      IF (line(1:3).ne.'RRR') THEN
           Print *, '::', line(1:76)
         Fehler = 'missed beginning of data blocks'
         RETURN
         ENDIF
      DO KK=1, nbspec_SD
         IF (nbchan_SD.ne.512) THEN
            Fehler = 'data format requires fixed nC'
            RETURN
            ENDIF
         READ (11, '(51(10i8/),2i8)', end=51) (buffer3(i),i=1,512)
 51      CONTINUE
         DO i=1,nbchan_SD
            Y(i)=buffer3(i)
            ENDDO
         IF (KK.eq.1) THEN
            rMon = 0
            DO i = 1, nC
               rMon = rMon + Y(i)
               ENDDO
            ENDIF
         CALL OlfPutSpe (j, KK, 1, 0.d0, nbchan_SD, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO
cdeb      Print *, 'read single detector spectra'

C  skip multidetector data ...

C  end of input file
      Close (11)

      TempIn = rAskD (' Temperature', TempIn)         ! Probentemp
      Temp   = TempIn
      Eelast = 81.805 / WaveL**2             ! in meV

C  K dependent parameters :
      Print *, ' flight path SD (m) = ', FPathSD
      DO K = 1, nK
         FPath (K) = FPathSD  ! Abstand Probe Detektor (m)
         DetE1(K) = -.0565
         DetE2(K) = -3.284
         ENDDO

C  Errors :
      RETURN
 91   CONTINUE
      CALL Compose2 (Fehler,
     *          'Error while reading 1st lines of "'//FileRaw, '"')
      Close(11)
      RETURN
 94   CONTINUE
      Fehler = 'Error while reading integer block'
      RETURN

      END ! RRT_In_NEAT

C  --------------------------------------------------------------------
      SUBROUTINE RRT_In_Mib (DirRaw, iCycle, irun, j,
     *                  rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
         !  Lecture des donnees du spectrometre temps de vol MIBEMOL de Saclay
         !  Version 1 by M.Bee, 25jan91
         !  Adapted to new raw data format by J.Wuttke, 25jul91
         !  New data format 10apr01 by A.Meyer
         !  Reads two input files d++++.ust and d++++.asc,
         !  where ++++ is a four-digit integer.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*)

      CHARACTER*(*)  Fehler, LongTit, DirRaw
      CHARACTER*80   inline, aus, FileUST, FileASC, nam
      CHARACTER      cl5*5, cl4*4, cl3*3, cl2*2, cv4*4

      IF (Fehler.ne.'&ff') THEN
         Print *, 'RRT_In_Mib/ Fehler on entry'
         RETURN
         ENDIF

C  Open the input files :
      IF (DirRaw.eq.' ') THEN
         nam = 'd'
      ELSE
         CALL Compose2 (nam, DirRaw, 'd')
         ENDIF
      CALL ReplaceT (nam, '&cycle', cv4(iCycle))

      IF (irun.ge.1 .and. irun.lt.9999) THEN
         isub = 0 ! accumulated data
      ELSEIF (irun.ge.100000 .and. irun.le.999999) THEN
         isub = mod(irun, 100)
         irun = (irun - isub) / 100
      ELSE
         Fehler = 'RRT_In_Mib/ bad number of run'
         RETURN
         ENDIF

      CALL Compose3 (FileUST, nam, cl4(irun), '.ust')
      IF (isub.le.0) THEN
         CALL Compose3 (FileASC, nam, cl4(irun), '.asc')
      ELSE
         CALL Compose3 (FileASC, nam, cl4(irun), '.pro'//cl2(isub))
         ENDIF
      CALL OpenFile (11, FileUST, '&noext', 'a', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OpenFile (12, FileASC, '&noext', 'a', Fehler)
      IF (Fehler.ne.'&ff') RETURN

C  Read *.ust :
      IF (irun.ge.4000) THEN ! new format encountered 2/98
         Read (11, '(x)') ! -/-
         IF (irun.ge.6000) THEN ! warum muessen die immer was aendern ?
            Read (11, '(19x,i4)') jrun
         ELSE
            Read (11, '(20x,i4)') jrun
            ENDIF
         Print *, ' reading ust for d'//cl4(jrun)
         Read (11, '(x)') ! date
         Read (11, '(x)') ! -/-
         Read (11, '(x)') ! subr max
         Read (11, '(x)') ! subr is
         Read (11, '(x)') ! n-count control
         Read (11, '(x)') ! T control
         Read (11, '(a)') inline
            CALL Fi1R (inline, Temp)
            print *, ' read T = ', Temp
         Read (11, '(x)') ! delta T
         Read (11, '(x)') ! waiting time
         Read (11, '(x)') ! grad
         Read (11, '(a)') inline
            CALL Fi1R (inline, delay)
            print *, ' read delay = ', delay
         Read (11, '(a)') inline
            CALL Fi1R (inline, ch_width)
            print *, ' read width = ', ch_width
         Read (11, '(x)') ! time preset
         Read (11, '(x)') ! Exp pars
         Read (11, '(x)') ! --------
         Read (11, '(19x,a17)') LongTit
            print *, ' read long title : ', LongTit(1:40)
         Read (11, '(x)') ! experimentateurs
         Read (11, '(23x,d8.0)') WaveL
            print *, ' read wavelength = ', WaveL
         Read (11, '(x)') ! geom
         Read (11, '(x)') ! angle
         Read (11, '(23x,d3.0)') fre_chops
         Read (11, '(23x,d1.0)') ratio_ch4
         Read (11, '(23x,d3.0)') ratio2
         channels = 512
         Read (11, '(30x,i2)')   nKentry
            print *, ' read nK = ', nKentry
         nK = nKentry
         Read (11, '(x)')
      ELSE
         ! d10.0 reads any real number which contains a decimal point
      Read (11, '(x)')                              ! three empty lines
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(37x,i5)')       jrun ! normally, it's i4.
      Read (11, '(x)')
      Read (11, '(x)')                              ! start time
      Read (11, '(6x,24x,x)')                       ! # subruns
      Read (11, '(6x,24x,a)')     LongTit           ! sample name
c      Print *, ' reading .ust titled ', LongTit(1:lenU(LongTit))
      Read (11, '(6x,24x,x)')                       ! transmission
      Read (11, '(6x,24x,d10.0)') WaveL             ! lambda
      Read (11, '(6x,24x,x)')                       ! container
      Read (11, '(6x,24x,d10.0)') Temp              ! sample temperature
      Read (11, '(6x,24x,x)')                       ! sample angle
      Read (11, '(6x,24x,d10.0)') ch_width          ! time channel width
      Read (11, '(6x,24x,d10.0)') delay   !(10-7sec)! time offset T_0
      Read (11, '(6x,24x,x)')                       ! timer preset
      Read (11, '(6x,24x,x)')                       ! monitor preset
      Read (11, '(6x,24x,d10.0)') channels          ! # tof channels
c      Print *, ' #channels is ', channels
      Read (11, '(6x,24x,i10)')   nKentry           ! # det groups
      nK = nKentry ! Ablesen aus Winkeltabelle waere auch nicht besser
c      Print *, ' #spectra is ', nK
      Read (11, '(6x,24x,d10.0)') fre_chops         ! frequency choppers
         ! 1/2 RPM choppers 1,2,5,6 (2 trous)
      Read (11, '(6x,24x,d10.0)') ratio_ch4         ! ratio ch4 / others
      Read (11, '(6x,24x,x)')                       ! user name
      Read (11, '(x)')
      ENDIF

C  Read the last block which contains the scattering angles :
      DO il = 1, MK/10+1
         Read (11, '(10f6.2)', end=211)
     *        (Angle(ii), ii=il*10-9, min0(nK,il*10))
         ENDDO
c      GOTO 212
 211  CONTINUE
c      Fehler = 'unexpected end-of-list in angle table'
c      RETURN
c 212  CONTINUE
      Close (11)

C  There is always a meaningless first entry in Angle :
      nK = nK - 1
      DO K = 1, nK
         Angle (K) = Angle (K+1)
         ENDDO

C  There may be empty angles at the end :
      DO WHILE (Angle(nK).le.0.d0)
         Print *, 'WARNING: unexpected zero angle at nK = ', nK
         IF (nK.lt.0) STOP
         nK = nK - 1
         ENDDO

C  Tell :
      Print *, 'PROVISORISCH nK = ', nK
      CALL Say2 (' Found '//cl3(nK), ' spectra as set in table')

C  Checks :
      IF    (jrun.ne.irun) THEN
         Print *, '#run : given/read ', irun, jrun
         Fehler = 'RRT_In_Mib/.ust/ Bad number of run'
         RETURN
      ELSEIF (nint(channels).ne.512) THEN
         Print *, 'channels = ', channels
         Fehler = 'RRT_In_Mib/.ust/ #channels<>512)'
         RETURN
      ELSEIF (qroutside(fre_chops,3.d2,3.d4)) THEN
         IF (fre_chops.ge.3.d4) THEN
            Fehler = 'RRT_In_Mib/.ust/ Chopper frequency must be in Hz'
            RETURN
         ELSE
            Print *, 'Chopper frequency times 60 ...'
            fre_chops = fre_chops * 60
            ENDIF
      ELSEIF (qroutside(ratio_ch4,1.d0,1.d1)) THEN
         CALL Gong (4)
         Print *, ' definition of "frequency chopper 4" has changed'
         Print *, ' it is now the ratio of frequencies (ch4/others)'
         ratio_ch4 = fre_chops / ratio_ch4
         IF (qroutside(ratio_ch4,1.d0,1.d1)) THEN
            Fehler = 'RRT_In_Mib/.ust/ value ratio_ch4 outside 1..10'
            RETURN
            ENDIF
         Print *, ' assuming it is meant ', ratio_ch4
         ENDIF

C  Get ToF parameters :
      Cwidth  = ch_width / 10                 ! in usec
      Period  = 60.d6 * ratio_ch4 / (2*fre_chops) / Cwidth
                                              ! (usec/rot)/(usec/ch)

C  Read *.asc :
      DO K = 1,nK
         Read(12,'(8f10.0)')  (Y(i),i=1,512)
            ! En realite le format d'ecriture est 8i10
         IF (K.eq.1) rMon = Y (509)
         CALL OlfPutSpe (j, K, 1, 0.d0, 508, X, Y, D, Fehler)
            ! ch. 509-512 are used for monitors, time, ?
         IF (Fehler.ne.'&ff') RETURN
         ENDDO
      Close (12)

C  Instrument parameters : detector efficiency, flight path
      DO K = 1, nK
         DetE1(K) =  0.
         DetE2(K) = -2.822
         ! dist_4_s   = 2.81d0 ! chopper4 - sample   ! not used at all
         FPath(K) = 3.58d0     ! sample - detectors
         ENDDO

      END ! RRT_In_Mib

C  --------------------------------------------------------------------
      SUBROUTINE RRT_In_Foc (DirRaw, irun, j,
     *                  rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
      ! currently not working uncomment for working !!! FKmar06
         !  by Fanni Juranyi, Nov. 2001
         !  after Mark Konnecke and Andreas Meyer
         !  Focus at SINQ
         !  adapted to new FOCUS file format by T.Unruh 21.12.2005
         !  adapted to new FOCUS file format by T.Unruh 30.10.2023
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_FLOAT, C_INT,
     * c_intptr_t
      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE "nexus/napif_iface.inc"
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*)

      CHARACTER*(*) LongTit, Fehler
      CHARACTER*80  inline, aus
      CHARACTER*80  Fileida, DirRaw !Artem add DirRaw, CHARACTER*40 -> CHARACTER*80
      CHARACTER     cl5*5
      CHARACTER*1   dettyp
      INTEGER       detno,detgrp
      REAL          detdist

C  Definitions for reading NeXus file:
      INCLUDE      'nexus/napif.inc'

      CHARACTER     cl4*4,cv5*5,cv6*6
      CHARACTER     titel*80, StartTime*20
      INTEGER(C_INT)     FILEID(NXHANDLESIZE), STATUS, iBank, iyear !Artem add iyear
      INTEGER       monitor, counts(4000*383), aaa, nZ
      DOUBLE PRECISION temperature
      DOUBLE PRECISION wavelength, Drpm, Period1, Z
      REAL          detang(384), t_o_f(4001)

      DATA  t_o_f/4001*0.0/, detang /384*0.0/, counts/1532000*0/,
     *      iBank /0/, iyear /2005/

      Fehler = '&ff' !

C  Open the input files :

      IF (irun.lt.1 .or. irun.gt.99999) THEN
         CALL Absturz ('RRT_In_Foc', 'Bad number of run : '//cl5(irun))
         ENDIF
      !CALL FrageCD ('measurement path ', path, '/home/FOCUS/data/')

      iyear=iAskDMu('Year of the measurement',iyear,1998,2093)

      IF(iyear.lt.2004) THEN
         CALL Compose3 (Fileida,trim(DirRaw)// !Artem: replace '/home/FOCUS/data/' -> trim(DirRaw)
     *                  cl4(iyear),'/focus'//
     *                  cv5(irun)//cl4(iyear),'.hdf')
      ELSE
         CALL Compose3 (Fileida,trim(DirRaw)//'focus'//cl4(iyear),'n'// !Artem add trim(DirRaw)
     *                  cv6(irun),'.hdf')
      ENDIF
      IF(NXOPEN(Fileida,NXACC_READ,FILEID) .NE. NX_OK) THEN
         RETURN
      ENDIF
      IF(NXOPENGROUP(FILEID,'entry1','NXentry') .NE. NX_OK) THEN
         CALL Absturz ('RRT_In_Foc', 'Could not open entry'//Fileida)
      ENDIF
      print *, 'Which detector bank do you want to use?'
      print *, '   (0) Merged'
      print *, '   (1) Upper'
      print *, '   (2) Middle'
      print *, '   (3) Lower'
      iBank = iAskDMu ('Choose option', iBank, 0, 3)
      IF(iBank.eq.0) THEN
       IF(NXOPENGROUP(FILEID,'merged','NXdata').NE.NX_OK) THEN
        IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
         CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//
     *                  Fileida)
        ELSE
         print *, 'old file, open middle bank'
        ENDIF
       ENDIF
      ELSEIF(iBank.eq.1) THEN
       IF(NXOPENGROUP(FILEID,'upperbank','NXdata').NE.NX_OK) THEN
        IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
         CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//
     *                 Fileida)
        ELSE
         print *, 'old file, open middle bank'
        ENDIF
       ENDIF
      ELSEIF(iBank.eq.2) THEN
       IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
        CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//
     *                Fileida)
       ENDIF
      ELSEIF(iBank.eq.3) THEN
       IF(NXOPENGROUP(FILEID,'lowerbank','NXdata').NE.NX_OK) THEN
        IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
         CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//
     *                 Fileida)
        ELSE
         print *, 'old file, open middle bank'
        ENDIF
       ENDIF
      ENDIF
         STATUS=NXOPENDATA(FILEID,'time_binning')
            STATUS=NXGETDATA(FILEID,t_o_f)

            nC=1
            DO WHILE(t_o_f(nC).NE.0)
               nC=nC+1
            ENDDO
            nC=nC-1
C ********** added y T. Unruh, 31.10.2023
C            IF(iyear.gt.2022) THEN
C               t_o_f(nC) = 0.0
C               nC=nC-1
C               DO aaa = 1, nC
C                  t_o_f(aaa) = t_o_f(aaa) + (t_o_f(2)-t_o_f(1))/2.0
C               ENDDO
C            ENDIF
C **********
            Cwidth=REAL(t_o_f(2)-t_o_f(1))
            IF(iyear.gt.2022) THEN
                t_o_f(nC) = 0.0
                nC=nC-1
                DO aaa = 1, nC
                    t_o_f(aaa) = t_o_f(aaa) + Cwidth*0.5 !Artem: replaced 03.03.2026
                ENDDO
            ENDIF
            print *, 'time bins:',t_o_f(1),t_o_f(2),t_o_f(nC),nC
         STATUS=NXCLOSEDATA(FILEID)
         STATUS=NXOPENDATA(FILEID,'theta')
            STATUS=NXGETDATA(FILEID,detang)
            nK=1
            DO WHILE(detang(nK).NE.0.0)
               nK=nK+1
            ENDDO
            nK=nK-1
            DO k=1,nK
               Angle(k)=detang(k)
            ENDDO
            print *, 'angles:',Angle(1),Angle(nK),nK
            STATUS=NXCLOSEDATA(FILEID)
         IF (iyear.lt.2004) THEN
            STATUS=NXOPENDATA(FILEID,'counts')
            STATUS=NXGETDATA(FILEID,counts)
            STATUS=NXCLOSEDATA(FILEID)
            STATUS=NXOPENDATA(FILEID,'monitor')
            STATUS=NXGETDATA(FILEID,monitor)
            STATUS=NXCLOSEGROUP(FILEID)
         ELSE
            STATUS=NXOPENDATA(FILEID,'counts')
            STATUS=NXGETDATA(FILEID,counts)
            STATUS=NXCLOSEDATA(FILEID)
            STATUS=NXCLOSEGROUP(FILEID)
            IF(NXOPENGROUP(FILEID,'FOCUS','NXinstrument').NE.NX_OK)
     *      THEN
               CALL Absturz (
     * 'RRT_In_Foc', 'Failed to open instrument'//Fileida)
            ELSE
               IF(NXOPENGROUP(FILEID,'counter','NXmonitor').NE.NX_OK)
     *         THEN
                  CALL Absturz (
     * 'RRT_In_Foc', 'Could not open counter'//Fileida)
               ELSE
                  STATUS=NXOPENDATA(FILEID,'monitor')
                  STATUS=NXGETDATA(FILEID,monitor)
                  STATUS=NXCLOSEDATA(FILEID)
                  STATUS=NXCLOSEGROUP(FILEID)
                  STATUS=NXCLOSEGROUP(FILEID)
               ENDIF
            ENDIF
         ENDIF
            print *,'Monitor:',monitor
            rMon=REAL(monitor)

C Fileida/entry1/: sample/name,temperature
      IF(NXOPENGROUP(FILEID,'sample','NXsample') .EQ. NX_OK) THEN
         IF(NXOPENDATA(FILEID,'name') .EQ. NX_OK) THEN
            IF(NXGETCHARDATA(FILEID,titel) .NE. NX_OK) THEN !Artem NXGETDATA -> NXGETCHARDATA
               print *, 'failed to read sample name'
            ELSE
               LongTit=titel
               print *, 'title:',LongTit
            ENDIF
         STATUS=NXCLOSEDATA(FILEID)
         ENDIF
         IF(NXOPENDATA(FILEID,'temperature') .EQ. NX_OK) THEN
            IF(NXGETDATA(FILEID,temperature) .NE. NX_OK) THEN
               print *, 'failed to read temperature'
            ELSE
               IF(temperature.LE.0.0) THEN
                  Temp=295
                  print *,
     * 'temperature is not measured, it is set to 295K.'
               ELSE
                  Temp=temperature
                  print *, 'temperature: ',Temp,'K'
               ENDIF
            ENDIF
           STATUS=NXCLOSEDATA(FILEID)
C         ELSE
C            Temp=295
C            print *,
C     * 'temperature is not measured, it is set to 295K.' !Artem: Add handling for cases where the temperature field is missing from the file.
         ENDIF
      STATUS=NXCLOSEGROUP(FILEID)
      ENDIF
C Fileida/entry1/: start_time
      IF(NXOPENDATA(FILEID,'start_time').EQ. NX_OK) THEN
         IF(NXGETCHARDATA(FILEID,StartTime).NE. NX_OK) THEN !Artem NXGETDATA -> NXGETCHARDATA
            print *, 'failed to read start time'
         ELSE
            print *, 'start time: ',StartTime
         ENDIF
      STATUS=NXCLOSEDATA(FILEID)
      ENDIF

      IF(NXOPENGROUP(FILEID,'FOCUS','NXinstrument').NE.NX_OK) THEN
         CALL Absturz ('RRT_In_Foc', 'Failed to open instrument'//
     *                 Fileida)
      ENDIF

C Fileida/entry1/FOCUS: monochromator/lambda
      IF(
     * NXOPENGROUP(FILEID,'monochromator','NXmonochromator').EQ.NX_OK)
     * THEN
         IF(NXOPENDATA(FILEID,'lambda') .EQ. NX_OK) THEN
            IF(NXGETDATA(FILEID,wavelength).NE. NX_OK) THEN
               CALL Absturz ('RRT_In_Foc', 'Failed to read wavelength'
     1//Fileida)
            ELSE
               WaveL=wavelength
               print *, 'wavelength: ',WaveL,'Angstrom'
            ENDIF
         STATUS=NXCLOSEDATA(FILEID)
         ENDIF
      STATUS=NXCLOSEGROUP(FILEID)
      ELSE
         CALL Absturz ('RRT_In_Foc', 'Failed to read wavelength'//
     *                 Fileida)
      ENDIF
      WaveL=rAskDMu('change wavelength?',WaveL,0.d6,50.d6)
      print *, 'wavelength: ',WaveL,'Angstrom'

C Fileida/entry1/FOCUS: disk_chopper/rotation_speed -> Period
      IF(NXOPENGROUP(FILEID,'disk_chopper','NXchopper').EQ.NX_OK) THEN
         IF(NXOPENDATA(FILEID,'rotation_speed').EQ.NX_OK) THEN
            IF(NXGETDATA(FILEID,Drpm).NE.NX_OK) THEN
               CALL Absturz ('RRT_In_Foc', 'Failed to read rotation
     1speed'//Fileida)
            ELSE
C two holes on the chopper, Period in microsec:
               Period = 30000000./Drpm ! Disk chopper speed value may be buggy
               Period1 = t_o_f(nC) - t_o_f(1)
               IF (dabs(Period-Period1).gt.Cwidth) THEN
                  Period = Period1
               ENDIF
               print *, 'Disk chopper speed: ',Drpm,'rpm'
               print *, 'Duty cycle: ',Period,'microsec'
            ENDIF
         STATUS=NXCLOSEDATA(FILEID)
         ELSE
            CALL Absturz (
     * 'RRT_In_Foc', 'Failed to read rotation speed'//Fileida)
         ENDIF
      STATUS=NXCLOSEGROUP(FILEID)
      ELSE
         CALL Absturz ('RRT_In_Foc', 'Failed to read rotation speed'//
     *                 Fileida)
      ENDIF

      STATUS=NXCLOSEGROUP(FILEID)
      STATUS=NXCLOSEGROUP(FILEID)
      STATUS=NXCLOSE(FILEID)


      Nmax = nC
      DO K = 1,nK
C        Detector effeciency
         DetE1(K) = -0.017            !  rectangular 30*10mm^2, 6bar
         DetE2(K) = -5.818            !  deff=10mm

         Fpath(K) = 2.50              !  in m
         DO i = 1,nC
C ********** added y T. Unruh, 31.10.2023 !Artem: I didn’t see a reason
C           IF(iyear.gt.2022) THEN
C               IF(K.eq.330) THEN
C                  Y(i)=0.0
C                  D(i)=0.0
C                  CONTINUE
C               ENDIF
C           ENDIF
C **********
           Y(i)=REAL(counts(i+((nC+1)*(K-1))))
           D(i)=sqrt(Y(i))
         ENDDO
         CALL OlfPutSpe (j, K, 1, 0.d0, Nmax, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
      ENDDO

      END ! RRT_In_Foc

C  --------------------------------------------------------------------
      SUBROUTINE RRT_In_Foco (DirRaw, irun, j,
     *                  rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
         !  by Andreas Meyer, May 2001
         !  Focus at SINQ

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*)

      CHARACTER*(*) DirRaw, LongTit, Fehler
      CHARACTER*80  FileRaw, FileFoc
      CHARACTER     cl5*5
      INTEGER       detno, detgrp
      REAL          detangle, detdist


      Fehler = '&ff' !

C  The following variables are no longer read in :
C  New file header necessary
      nK = 383          ! # of spectra
      nC = 631          ! # of channels
      WaveL = 4.09992   !AA
      Cwidth = 5.0      !mus
      Period = 3155.0   !mus

C  Open the input files :
      IF (irun.lt.1 .or. irun.gt.99999) THEN
         CALL Absturz ('RRT_In_Foc', 'Bad number of run : '//cl5(irun))
         ENDIF

      CALL Compose2 (FileFoc, 'f'//cl5(irun), '.inx')

C  Construct file name :
      FileRaw = DirRaw
      CALL Append (FileRaw, FileFoc)

      CALL OpenFile (11, FileRaw, '&noext', 'a', Fehler)

      IF (Fehler.ne.'&ff') RETURN

C  Read Usertable--File Header :
         ! d10.0 reads any real number which contains a decimal point
         ! f10.0 reads any integer number and transfer into a real one

      Read (11, '(x)')
      Read (11, '(a)') LongTit              ! sample name and comment

         Print *, ' sample: ',LongTit
         Print *, ' wavelength     = ',WaveL

C  Read *.ida :
      IF (MK.lt.nK) THEN
         CALL Absturz ('RRT_In_Foc', 'Recompile with MK >= nK')
         ENDIF
C  Data
      IF (MC.lt.nC) THEN
         CALL Absturz ('RRT_In_Foc', 'Recompile with MC >= nC')
      ENDIF

      Nmax = nC
      DO K = 1,nK
C        Detector effeciency
         DetE1(K) = -0.017            !  rectangular 30*10mm^2, 6bar
         DetE2(K) = -5.818            !  deff=10mm

         Fpath(K) = 2.50              !  in m
         Read (11, '(d10.0)') Angle(K)
         IF (K.eq.1) THEN
            Read (11, '(d10.0,4x,d10.0)') Temp,rMon
         ELSE
            READ (11, '(x)')
         ENDIF
         DO i = 1,nC
           Read(11,'(18x,d10.0)') Y(i)
           D(i)=sqrt(Y(i))
         ENDDO
         CALL OlfPutSpe (j, K, 1, 0.d0, Nmax, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (K.eq.nK) GOTO 222
         Read (11, '(x)')
         Read (11, '(x)')
         ENDDO
 222     CONTINUE
         Print *, ' Temperature: ',Temp
         Print *, ' Monitor:     ',rMon
         Print *

      Close (11)

      END ! RRT_In_Foco

C  --------------------------------------------------------------------
      SUBROUTINE RRT_In_DCS (DirRaw, irun, j, iCycle,
     *                  rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
      ! by F. Kargl 2003
      ! DCS @ NIST

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*),Par1(6),
     *              iPar2(1024)

      CHARACTER*(*) DirRaw, LongTit, Fehler
      CHARACTER*80  FileRaw, FileDCS,irTest
      CHARACTER     cl4*4,cl7*7

      REAL          rFPATH
      DATA          iAmod / 0 /


      Print * ,'------------------------------'
      Print * ,' options for spectra handling'
      Print * ,'------------------------------'
      Print * ,' (0) all angles'
      Print * ,' (1) exclude negative angles'
      Print * ,' (2) only negative angles'

      iAmod  = iAskDMu('Choose option',iAmod,0,2)

C     Create Filename
      IF(iabs(ndigits(irun)).eq.6) THEN
         FileRaw = '0'//cl7(irun)
      ELSEIF(iabs(ndigits(irun)).eq.7) THEN
       ! IF(cl7(irun)(1:1).lt.'1'.or.cl7(irun)(1:2).gt.'12')
         FileRaw = cl7(irun)
      ELSE
         CALL Absturz('irun is not a valid file number',Fehler)
      ENDIF

      CALL Compose4(FileDCS,cl4(iCycle),FileRaw(1:4),'_',FileRaw(5:7))
      CALL Append(FileDCS,'.asc')
      FileRaw = DirRaw
      CALL Append(FileRaw,FileDCS)


C     Open file
      CALL OpenFile(11,FileRaw,'&noext','a',Fehler)

      IF (Fehler.ne.'&ff') RETURN


C      Read in necessary data !!!
      Read(11,'(x)')
      Read(11,'(x)')
      Read(11,'(x)')
      Read(11,'(x)')
      Read(11,'(a3)') irTest
      IF (irTest(1:3).ne.'DCS') THEN
         CALL Absturz ('not the right instrument:',irTest)
         ENDIF
      Read (11,'(x)')
      Read (11,'(x)')
      Read (11,'(x)')
      Read (11,'(a20)')  irTest
      IF (irTest(3:6).eq.'rows') nK = ichar1(irTest(9:9))*100 +
     *                                 ichar1(irTest(10:10))*10 +
     *                                 ichar1(irTest(11:11))
      Read (11,'(a20)') irTest
      IF (irTest(3:9).eq.'columns') nC = ichar1(irTest(12:12))*1000 +
     *                                 ichar1(irTest(13:13))*100 +
     *                                 ichar1(irTest(14:14))*10 +
     *                                 ichar1(irTest(15:15))
      DO K = 1,nK
            Read (11,'(6(6x,d7.0))') (Y(i), i = 1,nC)
            CALL OlfPutSpe(j, K, 1, 0.d0, nC, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ENDDO
C     input of data is done
C     reading additional information

C -------------------------------------------------------
       Print *, 'read in accomplished !! '
C -------------------------------------------------------

      Read(11,'(x)')
      Read(11,'(a20)') irTest
C -------------------------------------------------------

C -------------------------------------------------------
      IF(irTest(9:14).eq.'ch_inp') THEN
         Read(11,'(x//)')
         Read(11,'(4(7x,d9.2))') (Par1(i), i = 1,6)
         ENDIF
         DO i = 1,6
            Print *, 'Par(',i,')',Par1(i)
            ENDDO
      Read(11,'(x)')
      Read(11,'(a20)') irTest
C -------------------------------------------------------
       Print *,'det dist ?',irTest ,'iPar field:',Par1(1),
     *                           Par1(6)
C -------------------------------------------------------
      IF(irTest(9:15).eq.'det_dis') THEN
         Read(11,'(x)')
         Read(11,'(7x,d9.0)') rFPath
         DO K = 1,nK
            FPath(K) = rFPath/1000.
            ENDDO
         ENDIF

C     Read monitor (field 4 in HistoHigh)

      Read(11,'(x)')
      Read(11,'(a20)') irTest
C -------------------------------------------------------

C -------------------------------------------------------
      IF (irTest(9:17).eq.'histohigh') THEN
         Read(11,'(x//)')
         DO K = 1,18
         Read(11,'(6(3x,i9))') (iPar2(i), i = 1,nC)
         IF (K.eq.4) THEN
            rMon = 0.
            DO i = 1,1024
               rMon = rMon + iPar2(i)
               ENDDO
            ENDIF
         ENDDO
         ENDIF
C -------------------------------------------------------
       Print *, 'Monitor value:', rMon,iPar2(800)
C -------------------------------------------------------
C     Read Temperature
      Read(11,'(x)')
      Read(11,'(a20)') irTest
      IF (irTest(9:16).eq.'temp_set') THEN
         Read(11,'(x)')
         Read(11,'(7x,d9.2)') Temp
         ENDIF

C     Read Angles
      Read(11,'(x)')
      Read(11,'(a20)') irTest
      IF (irTest(9:15).eq.'angle_p') THEN
         Read(11,'(x//)')
         Read(11,'(6(5x,d8.2))') (Angle(K), K = 1,nK)
         ENDIF


C sum over same negative and positive scattering angles
      IF (iAmod.eq.0) THEN
         Print *, 'negative not exluded'
         Write (35,'(/a27)') 'negative angles not excluded'
         DO K = 1,nK
            Angle(K) = dabs(Angle(K))
            ENDDO

C only positive scattering angles will remain
      ELSEIF (iAmod.eq.1) THEN
         Print *, 'negative angles excluded'
         Write (35,'(/a23)') 'negative angles excluded'
         DO K = 1,nK
            IF(Angle(K).le.0.) THEN
               Angle(K) = 0.
               Write(35,*) 'negative at K =', K
               ENDIF
            ENDDO
C now calcualtion of shifted spectra field
            KK = 1
         IF (Angle(1).eq.0) THEN
            KK = 2
            ENDIF

         DO K = 2,nK
            IF (Angle(K).eq.0) THEN
               DO K2 = 1,K-KK
                  CALL OlfGetXYD(j,K-K2,nC,X,Y,D,Fehler)
                  CALL OlfPutXYD(j,K-K2+1,nC,X,Y,D,Fehler)
                  Angle(K-K2+1)= Angle(K-K2)
                  ENDDO
                  Angle(KK)=0.
                  KK = KK+1
               ENDIF

            ENDDO

C only negative scattering angles will remain

      ELSEIF (iAmod.eq.2) THEN
         Print *, 'only negative angles'
         Write (35,'(/a20)') 'only negative angles'
         DO K = 1,nK
            IF (Angle(K).ge.0) THEN
               Angle(K) = 0.
            ELSE
               Angle(K)=dabs(Angle(K))
               Write(35,*) 'angle negative K =', K
               ENDIF
            ENDDO

! now calcualtion of shifted spectra field
            KK = 1
         IF (Angle(1).eq.0) THEN
            KK = 2
            ENDIF

         DO K = 2,nK
            IF (Angle(K).eq.0) THEN
               DO K2 = 1,K-KK
                  CALL OlfGetXYD(j,K-K2,nC,X,Y,D,Fehler)
                  CALL OlfPutXYD(j,K-K2+1,nC,X,Y,D,Fehler)
                  Angle(K-K2+1)= Angle(K-K2)
                  ENDDO
                  Angle(KK)=0.
                  KK = KK+1
               ENDIF

            ENDDO

         ENDIF


      DO K = 1,nK
         DetE1(K) = -.013  ! dummy value (has still to be checked)
         DetE2(K) = -4.039 ! calculated from values given for DCS
         ENDDO

      WaveL = Par1(1)
      Cwidth = Par1(6)/1000
      Period = Par1(6)
      Eelast = 81.805 / WaveL**2

      Close(11)

      END ! RRT_In_DCS

C  --------------------------------------------------------------------
      SUBROUTINE RRT_In_Fcs (DirRaw, irun, j,
     *                  rMon, Temp, WaveL, Cwidth, Period, LongTit,
     *                  nK, Angle, FPath, DetE1, DetE2, Fehler)
C  --------------------------------------------------------------------
         !  by Andreas Meyer, Feb. 1998
         !  NIST fermi choper tof

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      DIMENSION     Angle(*), FPath(*), DetE1(*), DetE2(*)

      CHARACTER*(*) LongTit, Fehler, DirRaw !Artem add DirRaw
      CHARACTER*80  inline, aus
      CHARACTER*40  Fileida
      CHARACTER     cl5*5
      CHARACTER*1   dettyp
      INTEGER       detno,detgrp
      REAL          detangle, detdist


      Fehler = '&ff' !

C  The following variables are no longer read in :
      nC         = 1024       ! # of tof channels
      nK         = 63         ! # of spectra

C  Open the input files :

      IF (irun.lt.1 .or. irun.gt.99999) THEN
         CALL Absturz ('RRT_In_Fcs', 'Bad number of run : '//cl5(irun))
         ENDIF

      CALL Compose2 (Fileida, 'f'//cl5(irun), '.ida')

      CALL OpenFile (11, Fileida, '&noext', 'a', Fehler)

      IF (Fehler.ne.'&ff')
     *   CALL Absturz ('RRT_In_Fcs', 'Could not open file '//Fileida)


C  Read Usertable--File Header :
         ! d10.0 reads any real number which contains a decimal point
         ! f10.0 reads any integer number and transfer into a real one

      Read (11, '(x)')                          ! three lines not read in
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(14x,a)')     LongTit          ! sample name and comment
      Read (11, '(x)')
      Read (11, '(20x,d10.0)') rMon             ! Monitor
      Read (11, '(x)')
      Read (11, '(20x,d10.0)') WaveL            ! wavelength in AA
      Read (11, '(20x,d10.0)') Cwidth           ! channel width in musec
      Read (11, '(20x,d10.0)') Period           ! Period in Hz
      Read (11, '(20x,d10.0)') Temp             ! sample temperature (K)
      Read (11, '(x)')
      Read (11, '(x)')

      Period = 1.0 / Period * 1000000.0  ! Period in musec

         Print *, ' sample: ',LongTit
         Print *, ' temperature    = ',Temp
         Print *, ' monitor        = ',rMon
         Print *, ' wavelength     = ',WaveL
         Print *

C  Read *.ida :
      IF (MK.lt.nK) THEN
         CALL Absturz ('RRT_In_Fcs', 'Recompile with MK >= 63')
         ENDIF

C  Detector angles and efficiency
      DO K = 1,nK
         Angle(K) = 0.0
      ENDDO

      DO K = 1, nK
         Read(11, '(1x,i3,1x,f8.4,1x,f8.1,2x,a1,3x,i2)')
     *       detno,detangle,detdist,dettyp,detgrp
         IF (detgrp.lt.37) THEN
           Angle(K) = detangle - 0.6363
           Read (11, '(x)')
         ELSEIF (detgrp.eq.37) THEN
           Angle(K) = detangle
         ELSEIF (detgrp.gt.37 .and. detgrp.lt.47) THEN
           Angle(K) = detangle - 0.6363
           Read (11, '(x)')
         ELSE
           Angle(K) = detangle
         ENDIF
         FPath(K) = detdist * 0.001   !  in meter
C  Detector effeciency
         DetE1(K) = -0.0806           !  cylindric d=25.4mm, 4bar
         DetE2(K) = -5.306            !  deff=20mm

      ENDDO

C  Data
      IF (MC.lt.nC) THEN
         CALL Absturz ('RRT_In_Fcs', 'Recompile with MC >= 1024')
      ENDIF

      Read (11,'(129(/),x)') ! Monitor in first block

      Nmax = min0(int(Period /Cwidth), nC)
      DO K = 1,nK
         Read(11,'(x)')
         DO l = 0,127
           Read(11,'(7(d11.0,1x),d11.0)')  (Y(i+l*8), i=1,8)
         ENDDO
         CALL OlfPutSpe (j, K, 1, 0.d0, Nmax, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      Close (11)

      END ! RRT_In_Fcs

C  --------------------------------------------------------------------
      SUBROUTINE RRT_EPeak_FCS (j, nK, Angle, FPath, Xep,
     *           WaveL, Eelast, Cwidth, Period, Fehler)
C  --------------------------------------------------------------------
         !  by A. Meyer, Feb. 1998
         !  special treatment for NIST fermi chopper tof
         !  shift elastic peak to correct position

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE          'l_def.f'
      INCLUDE          'i_dim.f'
      INCLUDE          'i_wrk.f'
      CHARACTER*(*)    Fehler
      DIMENSION        Angle(*), FPath(*), Xep(*)

      nC = 1024

      DO K = 1, nK
         Nelist = int(Xep(K))

         CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         Telast = FPath(K) / 3956.d-6 * WaveL   ! in usec

         Nmax = min0(int(Period /Cwidth), nC)
         Nfend= min0(Nmax - int(Telast / Cwidth) + Nelist, Nmax-1)

C  Wrap around :
         DO i = 1, nC
            Y2(i) = Y(i)
            D2(i) = D(i)
            ENDDO
         DO i = 1,Nfend
            Y(i + Nmax - Nfend) = Y2(i)
            D(i + Nmax - Nfend) = D2(i)
            ENDDO
         DO i = Nfend+1, Nmax
            Y(i - Nfend) = Y2(i)
            D(i - Nfend) = D2(i)
            ENDDO

         CALL OlfPutXYD (j, K, Nmax, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! K

      END ! RRT_EPeak_FCS

C  --------------------------------------------------------------------
      SUBROUTINE RRT_Check (what,
     *                Temp, Eelast, WaveL, Cwidth, Period,
     *                nK, Angle, FPath, DetE1, DetE2,
     *                Temp1, Eelast1, WaveL1, Cwidth1, Period1,
     *                nK1, Angle1, FPath1, DetE11, DetE21,
     *                Fehler)
C  --------------------------------------------------------------------

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)

      INCLUDE      'l_def.f'
      CHARACTER*(*)     Fehler, what
      DIMENSION         Angle(*), FPath(*), DetE1(*), DetE2(*),
     *                  Angle1(*), FPath1(*), DetE11(*), DetE21(*)

      IF (nK1.ne.nK) THEN
         CALL Compose2 (Fehler,
     *         what//'/ no. of data blocks '//cl4(nK), '<>'//cl4(nK1))
         RETURN
         ENDIF
      DO K = 1, nK
         IF (dabs(Angle1(K)-Angle(K)).ge.1.d0) THEN
            Fehler = what//'/ angles too different'
            RETURN
            ENDIF
         ENDDO

      IF (.not.qEqTol(Eelast1, Eelast, 5.d-2)) THEN
         Fehler = what//'/ wavelengths too different'
         RETURN
         ENDIF

      END ! RRT_Check

C  --------------------------------------------------------------------
      SUBROUTINE RRT_EPeak (j, nK, Angle, Xep, Yep, Bad, Fehler)
C  --------------------------------------------------------------------
             ! JWu 30jan93 : rewritten, using M.Bee's RRT_In_Mib

         ! calculates elastic peak area (mainly for use with vanadium),
         ! and peak position (mainly for use with sample).

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'
      CHARACTER*(*) Fehler
      DIMENSION     Xep(*), Yep(*), Angle(*), SX0(MK)
      CHARACTER*1   Bad(*)
      CHARACTER     aus*80

C  First pass : rough estimate of elastic channel :
      Write (35, '(a)') '  first pass'
      nK1    = 0
      xmep1  = 0.
      ymep1  = 0.
      DO K = 1, nK
         CALL OlfGetXYD (j,  K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         IF (Bad(K).eq.' ') THEN
            IF (Angle(K).le.0.d0) Fehler =
     *   'EPeak/ Invalid angle (should have been detected earlier), K ='
     *   //cl4(K)
            ! find maximum :
            ima = 1
            yma = 0.
            DO i = 1, n
               IF (Y(i).gt.yma) THEN
                  ima = i
                  yma = Y(i)
                  ENDIF
               ENDDO

            IF (yma.eq.0.) THEN
               Write (*, *) '  > no counts in data block ', K
               Write (35,*) '  > no counts in data block ', K
               Bad(K) = '0' ! empty
            ELSEIF (yma.le.0.) THEN
               Write (*, *) '  > negative counts in data block ', K
               Write (35,*) '  > negative counts in data block ', K
               Bad(K) = '-' ! negativ
            ELSE
               Xep(K) = dble(ima)
               Yep(K) = Y(ima)
               xmep1 = xmep1 + ima * Y(ima)
               ymep1 = ymep1 +       Y(ima)
               nK1  = nK1  + 1
               ENDIF
            ENDIF
         ENDDO
      xmep1 = xmep1 / ymep1 ! weighted average
      ymep1 = ymep1 / nK1

C  Second pass : check maximum, sum over region around it :
      ! initial tolerance (10*final) :
      xtol  = iExeMLP('tof-epeak-xtol-i')
      ytol  = iExeMLP('tof-epeak-ytol-i')
      xtolf = iExeMLP('tof-epeak-xtol-f')
      ytolf = iExeMLP('tof-epeak-ytol-f')
      IF (xtol.lt.1 .or. ytol.lt.1 .or. xtolf.lt.1 .or.ytolf.lt.1) THEN
         Fehler = 'Please set tolerances in ida.su'
         RETURN
         ENDIF
      Write (35, '(a)') '  second pass'
 21   CONTINUE
      nK2    = nK1
      nK1    = 0
      xmep2  = 0.
      ymep2  = 0.
      DO K = 1, nK
         IF (Bad(K).eq.' ') THEN
            ! check height of max :
            IF (Yep(K).lt.ymep1/ytol) THEN
               Write (aus, '(a,i3,a,g8.3,a,g8.3)')
     *            '  > peak too small in spectrum ', K, ': ',
     *            Yep(K), 'vs mean ', ymep1
               Write (* , '(a)') aus
               Write (35, '(a)') aus
               Bad(K) = 'v'  !  vanishing
               GOTO 29
               ENDIF

            IF (Yep(K).gt.ymep1*ytol) THEN
               Write (aus, '(a,i3,a,g8.3,a,g8.3)')
     *            '  > elefantic peak in spectrum ', K, ': ',
     *            Yep(K), 'vs mean ', ymep1
               Write (* , '(a)') aus
               Write (35, '(a)') aus
               Bad(K) = 'g'  !  giant
               GOTO 29
               ENDIF

            ! compare individual to average peak position :
            IF (dabs(Xep(K)-xmep1).gt.xtol) THEN
               Write (aus, '(a,i3,a,g8.3,a,g8.3)')
     *            '  > excentric peak in spectrum ', K, ': ',
     *            Xep(K), 'vs mean ', xmep1
               Write (* , '(a)') aus
               Write (35, '(a)') aus
               Bad(K) = 'e' !  excentric
               GOTO 29
               ENDIF

            ! determine resolution width :
            CALL OlfGetXYD (j,  K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            yhalf = Yep(K) / 2
            imax = idnint (Xep(K))
            DO ich=imax-1,1,-1
            IF (Y(ich).lt.yhalf .or. ich.eq.1) then
               ileft=ich
               goto 1                     ! Sortie mi-hauteur a gauche
               ENDIF
            ENDDO
 1          CONTINUE
            DO ich=imax+1, n
               IF (Y(ich).lt.yhalf .or. ich.eq.n) then
               iright=ich
               goto 2                      ! Sortie mi-hauteur a droite
               ENDIF
            ENDDO
 2          CONTINUE
            nfwhm = iright-ileft+1      ! Nombre de points a mi-hauteur
            nfwB  = min0 (  nfwhm, imax-1, n-imax)   ! pour barycentre
            nfwC  = min0 (2*nfwhm, imax-1, n-imax)   ! pour integrale
            ilB = iinside (imax-nfwB    , 1, imax)   ! pour barycentre
            irB = iinside (imax+nfwB    , imax, n)
            ilC = iinside (imax-nfwC    , 1, imax)   ! pour integrale
            irC = iinside (imax+nfwC    , imax, n)
            IF (irB.lt.ilB .or. irC.lt.ilC) THEN
               Print *, ' problems with peak: K imax ileft iright : ',
     *                 K, imax, ileft, iright
               STOP
               ENDIF

            ! Barycenter :
            x_peak   = 0.
            s_peak   = 0.
            DO ich=ilB, irB                ! Boucle sur canaux du pic
            x_peak   = x_peak   + Y(ich) * ich
            s_peak   = s_peak   + Y(ich)
            ENDDO
            IF (s_peak.le.0.d0) THEN ! encountered by AToelle jun95
               Write (*,  '(a)') '  > negative peak in spectrum ', K
               Write (35, '(a)') '  > negative peak in spectrum ', K
               Bad(K) = '-' ! negativ
               ENDIF
            Xep(K) = x_peak/s_peak
            xmep2  = xmep2 + Xep(K)
            ! sum over peak :
            ! sbusch 2007-dec: before, SX0 was a float containing S(0) of
            ! the current spectrum. It was written into Yep(K) to have access
            ! to it. However, doing so already here fills Yep with content that
            ! is unexpected at the beginning of this routine and will therefore
            ! crash the vanadium normalization if (after spectrum deletion)
            ! this routine is re-started.
            ! Solution: Write S(0) into an array SX0(K) and only after the
            ! routine is finished for sure, copy this array into Yep(K) to
            ! ensure full operationality.
            !SX0 = 0.
            SX0(K) = 0.


            DO I = ilC, irC
!              SX0 = SX0+Y(I)
               SX0(K) = SX0(K) + Y(I)
               ENDDO
!           Yep(K) = SX0
            ymep2  = ymep2 + Yep(K)
            nK1    = nK1  + 1
 29         CONTINUE
            ENDIF
         ENDDO

      IF (nK1.le.0) THEN
         Write (* , '(a)') '  > no spectrum has good shape'
         Write (35, '(a)') '  > no spectrum in good shape/'//
     *                     ' mean peak position undefined'
         RETURN
         ENDIF

      xmep1 = xmep2 / nK1
      ymep1 = ymep2 / nK1

      IF (nK1.lt.nK2) GOTO 21    ! averages changed -> repeat
      ! sbusch 2007-dec: only now, after the last GOTO, it is safe to
      ! write S(0) into Yep
      DO K = 1, nK
        Yep(K) = SX0(K)
        ENDDO

      ! decrease tolerance step by step :

      xtol = xtol/1.2
      ytol = ytol/1.2
      IF (xtol.lt.xtolf .and. ytol.lt.ytolf) RETURN

      xtol = dmax1 (xtol, xtolf)
      ytol = dmax1 (ytol, ytolf)

      END ! Analysis

C  --------------------------------------------------------------------
      SUBROUTINE RRT_Energ (qOutE, iSEG, Xep, Bad, Angle, FPath,
     *                   iFmod, NooW, NooD,
     *                   EFmax, EWmax, EFminRel,
     *                   WaveL, Eelast, Cwidth, Period,
     *                   j, nK, Fehler)
C  --------------------------------------------------------------------
         ! FK oct06 modified TOFTOF option
         ! FK jul05 add special for TOFTOF
         ! proper handling of frame overlap FK oct01
                 ! rewritten JWu 22/23jun92, X=tof option 30apr93
         ! shift elastic peak   (from Xep = XepS or XepV)
         ! shift frame overlap  (par. EFmax, EWmax, EFminRel)
         ! convert to energy    (if qOutE)

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*(*) Fehler
      DIMENSION     Xep(*), Angle(*), FPath(*)
      CHARACTER*1   Bad(*)

C   ----------------------------------------------------------------------
C      set counters zero
C   ----------------------------------------------------------------------
      NooF   = 0
      NooD   = 0
      rmNooF = 0.0
      imNFC  = 0
      rmNooD = 0.0
      imNDC  = 0
      rChS   = 0.0

C   ----------------------------------------------------------------------
C     FILE HANDLING
C   ----------------------------------------------------------------------

      DO K = 1, nK
         CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         Telast = FPath(K) / 3956.d-6 * WaveL ! in usec ! TBM
         rDelay = (Period - (n+1)*Cwidth)/Cwidth
         ! rDelay could be an error source for Focus
         ! could be an error source for DCS too
         ! appended FK Nov01
         IF (rDelay.lt.0) THEN
            rDelay = 0.0
            ENDIF ! because of negative rDelay

C   ----------------------------------------------------------------------
C   special TOFTOF
C   ----------------------------------------------------------------------
      IF (iFmod.eq.4.and.K.eq.1) THEN
         rChS = Period/Cwidth - int(Period/Cwidth)
         IF (rChS.gt.1E-3*Cwidth) THEN
         WRITE (35,*) 'rChS for TOFTOF eq. ', rChS
         ELSE
            rChS = 0.0
            ENDIF
       ENDIF

C   ----------------------------------------------------------------------
C   calculate dead-time between two frames
      IF (iFmod.eq.2 .and. K.eq.1) THEN
         IF (int(rDelay).gt.5) THEN !if dead time > 5*Cwidth
            Print *
            Print *, 'Check Period in Raw Data File !!!'
            Print *
         ENDIF
      ENDIF
C ------------------------------------------------------------------------

C  Prepare determination of the frame :

C   Prepare necessary parameters
      IF (iFmod.eq.2.or.iFmod.eq.4) THEN
         EFmin  = EFminRel * Eelast
         TFmax  = Telast/(dsqrt(EFmin/Eelast+1))
         TFmin  = Telast/(dsqrt(EFmax/Eelast+1))


      ELSEIF (iFmod.eq.3) THEN
         Fehler = 'RRT_Energ/ Option 3 ausser Betrieb'
         RETURN

         ENDIF ! iFmod

C  channel number of elastic peak :
         REPP = Xep(K)

C  Set the frame :


         IF (iFmod.eq.2.or.iFmod.eq.4) THEN
            NooFL = int((Telast-TFmin)/Cwidth)+1
            NooFR = int((TFmax-Telast)/Cwidth)+1
            NooF  = int(REPP-(Telast-TFmin)/Cwidth)+1
         ELSE
            NooW = NooF + NooD
            ENDIF

C  Move around :
         IF (iFmod.eq.0) THEN
            ! save first NooF time channels :
            DO i = 1,NooF
               Y2(i) = Y(i)
               D2(i) = D(i)
               ENDDO
            ! delete first NooW channels :
            DO i = 1, n
               X(i) =  i-REPP ! distance from elastic peak, in channels
C               D(i) =  D(i)
C               Y(i) =  Y(i)
               ENDDO
            ! append the overlapping NooF channels at the end :
            DO i = 1, NooF
               X (n-NooW+i) = n + i-REPP ! + Period
               Y (n-NooW+i) = Y2(i)
               D (n-NooW+i) = D2(i)
               ENDDO
               nout = n - NooW + NooF

         ELSEIF (iFmod.eq.2.or.iFmod.eq.4) THEN
            IF(NooF.ge.0 .and. (Bad(K).eq.' ')) THEN !FK10/02 Bad(K) added
! move only left channels and cut channels on the right
              ! save first NooF time channels :
               DO i = 1,NooF
                  Y2(i) = Y(i)
                  D2(i) = D(i)
                  ENDDO
              ! delete first NooW channels :
               DO i = 1+NooF, n
                  X(i-NooF) =  i-REPP ! distance from elastic peak, in channels
                  D(i-NooF) =  D(i)
                  Y(i-NooF) =  Y(i)
                  ENDDO
              ! append the overlapping NooF channels at the end :
               DO i = 1, NooF
                  X (n-NooF+i) = n + i - REPP + rDelay
                  Y (n-NooF+i) = Y2(i)
                  D (n-NooF+i) = D2(i)
                  ENDDO
C  ------------------------------------------------------------------------
C  cut channels if REPP+NooFR < n + NooF + rDelay
C  ------------------------------------------------------------------------
                  IF ((n-REPP).gt.NooFR) THEN
                     NooD = NooF+n-REPP-NooFR
                     nout = n - NooD
                     rmNooD = rmNooD + NooD
                     imNDC = imNDC + 1
                  ELSEIF(((n-REPP).lt.NooFR)
     *                .and.(n-REPP+NooF+int(rDelay)).gt.NooFR) THEN
                     NooD = min0(int(n-REPP+NooF+
     *                           rDelay-NooFR),NooF) !number of channels to be cut
                     nout = n - NooD
!FK10/02                     rmNooD = rmNooD + NooD
!FK10/02                     imNDC = imNDC + 1
                   ELSE
                      nout = n
                     ENDIF !cut channels

            ELSEIF ((NooF+int(rDelay).lt.0) .and. (Bad(K).eq.' ')) THEN
                  NooF = NooF + int(rDelay)
              ! move necessary right channels to the left and cut both sides
               ! save ncs channels
               DO i = n+NooF+1,n
                  Y2(i-n-NooF) = Y(i)
                  D2(i-n-NooF) = D(i)
                  ENDDO
              ! move channels around
               DO i = 1, n+NooF
                  X(n-i+1) = n+NooF+1-i-REPP
                  Y(n-i+1) = Y(n+NooF+1-i)
                  D(n-i+1) = D(n+NooF+1-i)
                  ENDDO
              ! append right NooF channels to left
                  rChS = rChS + 2 * int(rDelay)
               DO i = 1, -NooF
                  X(i) = i - REPP + NooF - rChS
                  Y(i) = Y2(i)
                  D(i) = D2(i)
                  ENDDO
C  ----------------------------------------------------------------------
C  cut channels if REPP+NooFR < n+int(rDelay)+NooF
C  ----------------------------------------------------------------------
                  IF((n-REPP+NooF).gt.NooFR) THEN
                     NooD = n-REPP+NooF-NooFR !number of channels to be cut
                     nout = n - NooD
!FK10/02                     rmNooD = rmNooD + NooD
!FK10/02                     imNDC = imNDC + 1
                  ELSE
                     nout = n
                     ENDIF !cut channels
C  ------------------------------------------------------------------------
            ELSE
               ! if Bad(K).ne.' '
               NooF = 0
                ! save first NooF time channels :
C               DO i = 1,NooF
C                  Y2(i) = Y(i)
C                  D2(i) = D(i)
C                  ENDDO
              ! delete first NooW channels :
               DO i = 1+NooF, n
                  X(i-NooF) =  i-REPP ! distance from elastic peak, in channels
C                  D(i-NooF) =  D(i)
C                  Y(i-NooF) =  Y(i)
                  ENDDO
              ! append the overlapping NooF channels at the end :
C               DO i = 1, NooF
C                  X (n-NooF+i) = n + i-REPP
C                  Y (n-NooF+i) = Y2(i)
C                  D (n-NooF+i) = D2(i)
C                  ENDDO
               nout = n
            ENDIF ! NooF-Loop

         ENDIF  !iFmod move around

C  Calculate average number of channels moved around
         IF (Bad(K).eq.' ') THEN
            rmNooF = rmNoof + NooF
            imNFC = imNFC + 1
            rmNooD = rmNooD + NooD !FK10/02 added
            imNDC = imNDC + 1      !FK10/02 added
        ENDIF

C Set energies, and convert to S(2Th,w) :
         DO i = 1, nout
            IF (qOutE) THEN
               R1 = 1+Cwidth/Telast*X(i)
               E = Eelast*(1.d0/R1/R1-1)
               fac = .5*Telast/Eelast*R1*R1*R1/Cwidth*R1 !The last R1 for Ki/Kf
            ELSE ! channel number -> time-of-flight
               E = ( Telast + Cwidth*X(i) ) / FPath(K) / 1000.
                  ! units: (usec/m) -> (msec/m)
               fac = 1. ! no conversion, counts remain counts
               ENDIF
            X(i) = -iSEG * E      ! 7jul93. Vercors group has -(-1)
            D(i) = fac * D(i)
            Y(i) = fac * Y(i)
            ENDDO

         CALL OlfPutXYD (j, K, nout, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! K


      imNooF = int(rmNooF/imNFC)
      IF(imNDC.ne.0) THEN
      imNooD = int(rmNooD/imNDC)
      ELSE
      imNooD = 0
      ENDIF
C  If channels were moved around :
      IF (iFmod.gt.0) THEN
C         CALL Say2 ('  throw away '//cl6(NooW-NooF), ' high energy channels')
         CALL Say2 ('  move around '//cl6(imNooF), ' channels')
         CALL Say2 ('  throw away '//cl6(imNooD),
     *                   ' of them')
         CALL Say3 ('  thus retaining '//cl6(n-imNooD),
     *                   ' of '//cl6(n), ' original channels')
         ENDIF

      END ! EppCon

C  --------------------------------------------------------------------
      SUBROUTINE RRT_DetEff (j, nK, E0, C1, C2, Fehler)
C  --------------------------------------------------------------------
            ! Correct for the variation of detector efficiency with
            !   neutrons energy, according to the coefficients C1, C2
            !   which have to be preset in the input routine.
            ! The influence of the detector stainless steel walls
            !   which produce a cutoff around 4 Angstroms is neglected.

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     C1(*), C2(*)
      CHARACTER*(*) Fehler

      write (35,'(a3,6a10)') 'K', 'C1(K)', 'C2(K)', 'E0', 'X(1)',
     *                       'X(n)', 'eff0'

      DO K = 1, nK
         IF (C1(K).ge.0 .or. C2(K).ge.0) THEN
            write (35, '(i3,a,2(2x,g8.3))') K, C1(K), C2(K),
     *       ' no detector efficiency correction (parameters not set)'
         ELSE
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            eff0 = dexp1(C1(K)/dsqrt(E0)) * (1-dexp1(C2(K)/dsqrt(E0)))
            write (35, '(i3,6(2x,g8.3))') K, C1(K), C2(K), E0, X(1),
     *                                    X(n), eff0
            DO i = 1, n
               E = E0 + X(i)
               IF (E.gt.1.d-8) THEN
                  eff = dexp1(C1(K)/dsqrt(E)) *
     *                       (1-dexp1(C2(K)/dsqrt(E)))
               ELSE
                  eff = 1
                  ENDIF
               Y(i) = Y(i)*eff0/eff
               D(i) = D(i)*eff0/eff
               ENDDO
            CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ENDIF

         ENDDO ! K

      END! DetEff

C  --------------------------------------------------------------------
      SUBROUTINE RRT_Norm (j, nK, YepV, BadV, BadS, Fehler)
C  --------------------------------------------------------------------
              !Normalization to Vanadium

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*1   BadV(*), BadS(*)
      DIMENSION     YepV(*)
      CHARACTER*(*) Fehler

      DO K = 1, nK
         CALL OlfGetXYD (j,  K, n, X, Y, D, Fehler)

         IF     (BadV(K).ne.' ') THEN
            BadS(K) = 'n'
         ELSEIF (YepV(K).le.0) THEN
            Fehler = 'RRT_Norm/ van-int=0 not recognized in _EPeak'
            RETURN
         ELSE
            DO i = 1, n
               Y(i) = Y(i) / YepV(K)
               D(i) = D(i) / YepV(K)
               Y(i) = Y(i)
               D(i) = D(i)
               ENDDO
            ENDIF

         CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         ENDDO ! K
      END ! RRT_Norm

C  --------------------------------------------------------------------
      SUBROUTINE RRT_Add (j, nK, Angle, Fehler)
C  --------------------------------------------------------------------
         ! Add spectra if angles are equal
            ! mrz95 by Andreas Meyer for data treatment of SIL tof spectra
            ! JWu 3mrz95 generalized for IN6 (3 angles may be equal)

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     Angle(*)
      CHARACTER*(*) Fehler

      K = 1
 300  CONTINUE
         KK = K+1 ! next K (if nothing happens..)
         CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
         IF (Angle(K).le.0) GOTO 309
         DO KK = K+1, nK
            IF (dabs(Angle(K)-Angle(KK)) .gt. Angle(K)/1000) GOTO 301
            ENDDO
 301     CONTINUE ! found out that AngleS(K)=..=AngleS(KK-1)
         IF (KK.le.K+1) GOTO 309 ! angles are *not* equal
         ! sum
         DO i = 1, n
            D(i) = D(i)**2
            ENDDO
         DO KKK = K+1, KK-1
            CALL OlfGetXYD (j, KKK, n, X1, Y1, D1, Fehler)
            DO i = 1, n
               Y(i) = Y(i) + Y1(i)
               D(i) = D(i) + D1(i)**2
               ENDDO
            Angle(KKK) = 0 ! spectrum will be deleted later
            ENDDO ! KKK
         DO i = 1, n
            Y(i) =       Y(i)  / (KK-K)
            D(i) = dsqrt(D(i)) / (KK-K)
            ENDDO
         Write (35,'(i1,a,2i4,g9.3)') KK-K, ' spectra summed: ',
     *      K, KK, Angle(K)
         CALL OlfPutXYD (j, K, n, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
 309     CONTINUE
         IF (KK.le.nK) THEN
            K = KK
            GOTO 300
            ENDIF

      END ! RRT_Add

C  --------------------------------------------------------------------
      SUBROUTINE RRT_Sort (j, nK, Angle, Temp, qOutE, Fehler)
C  --------------------------------------------------------------------
         ! sort energies, eliminate spectra with angle=0
         ! set coordinate names

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h,o-p,r-z)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*(*) Fehler
      DIMENSION     Angle(*)

      Write (35, '(a)') ' RRT_Sort/ begin'
      jout = 0
      CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      IF (qOutE) THEN
         CALL OlfCnuP (jout, 'x', 'w', 'meV', Fehler)
         CALL OlfCnuP (jout, 'y', 'S(2th,w)', 'meV-1', Fehler)
      ELSE
         CALL OlfCnuP (jout, 'x', 'tof', 'msec/m', Fehler)
         CALL OlfCnuP (jout, 'y', 'I(2th,tof)', 'cts', Fehler)
         ENDIF
      CALL OlfCnuP (jout, 'z1', '2th', ' ', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      Kout = 0
      DO K = 1, nK
         IF (Angle(K).ne.0.) THEN
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ! sort spectrum :
            DO i = 1, n
               X1(n+1-i) = X(i)
               Y1(n+1-i) = Y(i)
               D1(n+1-i) = D(i)
               ENDDO
            Kout = Kout + 1
            CALL OlfPutSpe (jout, Kout, 1, Angle(K), n,
     *                      X1, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDIF
         ENDDO

      CALL OlfClos (jout, Kout, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL MemFileDel (j, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      Write (35, '(a)') ' RRT_Sort/ end'

      END ! RRT_Sort

C  ====================================================================
C  i8  (7) :   read HFBS data
C  ====================================================================

      SUBROUTINE RRT_In_Hfbs (Fehler)
C     ----------------------------
      ! rebuild Florian Kargl 08/03 (new data format)
        ! Andreas Meyer 12/99
      ! Read converted HFBS (NIST) data and normalize to monitor (FC)

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      CHARACTER*5    crun
      CHARACTER*40   Filehfbs,cTest
      CHARACTER*30   aus, File, Title, LongTit

      DIMENSION rPM(4096),rPD(4096),XE(4096),rMon(4096),rWB(500),
     *          rTemp(500),rDoppl(500),rTB(500),zq(16),dBG(16)

      DATA Title / ' ' /, File / ' ' /,iym/200308/,qTheSame/.false./,
     *     nLoop/0/,qMon/.true./,qEnerg/.true./,qBG/.false./

      Fehler = '&ff' !

         Print *, ' read raw HFBS data files '

 105  CONTINUE               ! outer loop
      nLoop = nLoop+1


      CALL FrageH('sample run number (mm_nn) ',crun)
      IF (crun.eq.'-1') GOTO 151

      IF ((crun(1:2).lt.'01'.or.crun(1:2).gt.'31') .or. (crun(4:5).lt.
     *              '01' .or.  crun(4:5).gt.'99')) THEN
C         CALL Absturz ('RRT_In_Hfbs' , 'Bad number of run : ' //crun)
         Print *, ('Not a valid Filenumber !')
         GOTO 151
         ENDIF



C  The following variables are fixed :

      nK = 16   ! # of spectra: detectors + monitor


C  Select between normalization to monitor, normalization
C  and background subtraction according to algorithm
C  or none of them.

      IF (nLoop.gt.1) THEN
      qTheSame = qAskD ('The same as before ?',intq(qTheSame))
      IF (qTheSame) GOTO 110
      ENDIF

      iym  = iAskD('offset:',iym)

      qEnerg = qAskD('Transformation channel to energy ?',intq(qEnerg))
      qMon = qAskD ('Normalization to monitor ?',intq(qMon))
      IF (qMon) THEN
      qBG  = qAskD ('Background subtraction ?',intq(qBG))
      ELSE
      qBG = .false.
      ENDIF
      qThetaq = qAskD ('Transformation 2th to q?',intq(qThetaq))

C  Open the input files :

 110  CONTINUE

      CALL Compose2 (Filehfbs, cl6(iym),crun)
      CALL Compose2 (Filehfbs, Filehfbs,'.hfbs')
      Print *, 'Filehfbs: ',Filehfbs


      CALL OpenFile (11, Filehfbs, '&noext', 'a', Fehler)

      IF (Fehler.ne.'&ff')
     *   CALL Absturz ('RRT_In_Hfbs', 'Could not open file '//Filehfbs)


      Read (11, '(21x,a)')        LongTit   ! sample name and comment
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(21x,i3)') Ncycl            ! number of cycles
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      DO i = 1,Ncycl
         Read (11, '(2x,d14.11)') rDoppl(i)
         ENDDO
      Read (11, '(x)')
      DO i = 1,Ncycl
         Read (11, '(2x,d9.1)') rWB(i)
         ENDDO
      Read (11, '(x)')
      DO i = 1,Ncycl
         Read (11, '(2x,d9.1)') rTB(i)
         ENDDO
      Read (11, '(x)')
      DO i = 1,Ncycl
         Read (11, '(2x,d6.1)') rTemp(i)
         ENDDO
! temperature check
         sTemp = 0.
      DO i = 1, Ncycl
         sTemp = sTemp + rTemp(i)/Ncycl
         ENDDO
         aTemp = 0.
      DO i = 1, Ncycl
         aTemp = aTemp + (rTemp(i)-sTemp)**2/Ncycl
         ENDDO
         aTemp = dsqrt0(aTemp)
         ! arbitrary value for mean deviation
         IF (aTemp.gt.2.) THEN
            Print *,('Check File for Temperatures !!!, aTemp = '),aTemp
            ENDIF
! This field is the fission chamber monitor
      Read (11, '(x)')
      DO i = 1, 4096
         Read (11,'(5x,d10.2)') rMon(i)
         ENDDO
! The first 16 Detectors will now be read in



         Print *, ' sample:    ',LongTit
         Print *, ' Cycles: ',Ncycl
         Print *

      CALL OlfCreate (jout, Kout, File, Title, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL tOlfG (jout, 'fil', File, Fehler)
      CALL tOlfG (jout, 'tit', Title, Fehler)
      IF (Fehler.ne.'&ff') GOTO 120

! setting axes labels
      IF (qEnerg) THEN
         CALL OlfCnuP (jout, 'x', 'w', 'ueV', Fehler)
         IF (Fehler.ne.'&ff') RETURN
      ELSE
         CALL OlfCnuP (jout, 'x', 'ch#', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF

      IF (qMon) THEN
         CALL OlfCnuP (jout, 'y', 'norm. Int.', 'meV-1', Fehler)
         IF (Fehler.ne.'&ff') RETURN
      ELSEIF (qMon .and. qBG) THEN
         CALL OlfCnuP (jout, 'y', 'I(q)', 'meV-1', Fehler)
         IF (Fehler.ne.'&ff') RETURN
      ELSE
         CALL OlfCnuP (jout, 'y', 'Int.', 'cts', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF

      IF(.not.qThetaq) THEN
         CALL OlfCnuP (jout, 'z1', '2th', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN
      ELSE
         CALL OlfCnuP (jout, 'z1', 'q', 'A-1',Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDIF
      Einc = 2.080000
      CALL rOlfP (jout, 'E0', 'meV', Einc , Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL rOlfP (jout, 'T', 'K', sTemp, Fehler)
      IF (Fehler.ne.'&ff') RETURN

! 2th-values for all detectors

         zq(1)=14.46
         zq(2)=20.98
         zq(3)=27.08
         zq(4)=32.31
         zq(5)=36
         zq(6)=43.77
         zq(7)=51.5
         zq(8)=59.25
         zq(9)=67
         zq(10)=74.75
         zq(11)=82.5
         zq(12)=90.25
         zq(13)=98
         zq(14)=105.75
         zq(15)=113.15
         zq(16)=121.25

! If qThetaq true then recalculate q values
         IF(qThetaq) THEN
            DO K = 1,nK
               zq(K) = 2.0 * dsind(zq(K) / 2.0)
               IF((roundN(zq(K),2)).lt.1.) THEN
                  zq(K) = roundN(zq(K),2)
                  ELSE
                  zq(K) = roundN(zq(K),3)
                  ENDIF
               Print *, ('z('),K,(') = '),zq(K)
               ENDDO
            ENDIF

C -------  getting raw data ------------------------------------
      DO K = 1,nK
         Read (11,'(x)')
         DO i = 1,4096
           Read(11,'(d4.0,d10.2)')  X(i),Y(i)
           D(i) = sqrt(Y(i))
         ENDDO

         CALL OlfPutSpe (jout, K, 1, zq(K), 4096, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

         Read (11,'(x)')
         DO i = 1,4096
            Read (11,'(4x,d10.2)') XE(i)
            ENDDO

C create own XE(i)
         DO i = 1,4096
            XE(i) = -200. + (400.0/4095.0*(i-1.))
C            Print *, ('XE('),i,(') = '),XE(i)
            ENDDO

C read in correction factors
         Read (11,'(x)')
         DO i = 1,4096
            Read (11,'(4x,d10.2)') rPD(i)
            ENDDO
         Read (11,'(x)')
         DO i = 1,4096
            Read (11,'(4x,d10.2)') rPM(i)
            ENDDO


C calculate cutoff ranges !!
            iSE = 2
            iEE = 0
         DO i = 2,4096
            IF ((rPD(i).eq.0.0) .and. (XE(i).lt.-1.0)) THEN
               iSE = iSE + 1
               ENDIF
            IF ((rPD(i).gt.0.0) .and. (XE(i).gt.1.0)) THEN
               iEE = i
               ENDIF
            ENDDO
            Print *,('iSE = '),iSE
            Print *,('iEE = '),iEE


         ! writing energy in channel nummbers
         IF (qEnerg) THEN

            IF (qMon) THEN !normalization to monitor

               IF (qBG) THEN ! constant BG subtr a la Dave


                  DO K = 1, nK
                     CALL OlfGetSpe (jout,K,1,zq(K),4096,X,Y,D,Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     rBGm = 0.
                     rBGp = 0.
                     rBGMm = 0.
                     rBGMp = 0.
                     rBGg = 0.
                     drBGm = 0.
                     drBGp = 0.
                     drBGg = 0.
                     DO i = iSE,iEE
                        IF (XE(i).le.-5.0) THEN
                           rBGm = rBGm + (Y(i)/rPD(i))*rPM(i)/rMon(i)*
     *                          (XE(i+1)-XE(i))
                           rBGMm = rBGMm + (rPM(i)/rMon(i))*(XE(i+1)-
     *                             XE(i))
                           drBGm = drBGm + ((Y(i)/rPD(i))/
     *                              (rMon(i)/rPM(i))**2)*
     *                              (XE(i+1)-XE(i))**2
                        ELSEIF (XE(i).ge.5.0) THEN
                           rBGp = rBGp + (Y(i)/rPD(i))*rPM(i)/rMon(i)*
     *                          (XE(i+1)-XE(i))
                           rBGMp = rBGMp + (rPM(i)/rMon(i))*(XE(i+1)-
     *                             XE(i))
                           drBGp = drBGp + ((Y(i)/rPD(i))/
     *                              (rMon(i)/rPM(i))**2)*
     *                              (XE(i+1)-XE(i))**2
                           ENDIF
                        ENDDO
                     rBGg = (rBGp-rBGm)/(rBGMp-rBGMm)
                     drBGg = dsqrt0(drBGp + drBGm)/(rBGMp-rBGMm)
C                     Print *,('rBGp('),K,('): '),rBGp,('rBGm('),K,('): '),
C     *                    rBGm,('rBGMp('),K,('): '),rBGMp,('rBGMm('),
C     *                    K,('): '),rBGMm,('rBGg('),K,('): '),rBGg,
C     *                    ('drBGg('),K,('): '),drBGg,('drBGp('),
C     *                    K,('): '),drBGp,('drBGm('),K,('): '),drBGm

                     DO i = iSE,iEE
                        IF (rMon(i).ne.0.0 .and. rPD(i).ne.0.0) THEN
! error takes also monitor error into account
                          D(i) = dsqrt0(Y(i)+dquot0(Y(i)**2,rMon(i)))*
     *                          dquot0(rPM(i),rPD(i)*rMon(i))
C old                         D(i) = dsqrt0(Y(i))/rPD(i)*rPM(i)/rMon(i)
                           Y(i) = (Y(i)/rPD(i)-rBGg)*rPM(i)/rMon(i)
                        ELSE
                           Y(i) = 0.
                           D(i) = 0.
                           ENDIF
                        ENDDO

C --- cut data at end positions ---------------------------------
                     nE = iEE - iSE + 1
                     DO i = 1,nE
                        XE(i) = XE(i-1+iSE)
                        Y(i) = Y(i-1+iSE)
                        D(i) = D(i-1+iSE)
                        ENDDO

                     CALL OlfPutSpe (jout, K, 1, zq(K), nE, XE,Y,
     *                       D,Fehler)
                     IF (Fehler.ne.'&ff') RETURN
                     ENDDO !K

               ELSE ! without BG


                  DO K = 1, nK
                     CALL OlfGetSpe (jout,K,1,zq(K),4096,X,Y,D,Fehler)
                     IF (Fehler.ne.'&ff') RETURN

                     DO i = iSE,iEE
                        IF (rMon(i).ne.0.0 .and. rPD(i).ne.0.0) THEN
! error takes also monitor error into account
                          D(i) = dsqrt0(Y(i)+dquot0(Y(i)**2,rMon(i)))*
     *                          dquot0(rPM(i),rPD(i)*rMon(i))
C old                         D(i) = dsqrt0(Y(i))/rPD(i)*rPM(i)/rMon(i)
                           Y(i) = Y(i)*rPM(i)/rMon(i)/rPD(i)
                        ELSE
                           Y(i) = 0.
                           D(i) = 0.
                           ENDIF
                        ENDDO

C --- cut data at end positions ---------------------------------
                     nE = iEE - iSE + 1
                     DO i = 1,nE
                        XE(i) = XE(i-1+iSE)
                        Y(i) = Y(i-1+iSE)
                        D(i) = D(i-1+iSE)
                        ENDDO


                        CALL OlfPutSpe (jout, K, 1, zq(K), nE, XE,Y,
     *                       D,Fehler)
                        IF (Fehler.ne.'&ff') RETURN
                     ENDDO !K
                  ENDIF !BG


            ELSE !only transformation to energy, but no normalization
                 ! to monitor !!!


               DO K = 1, nK
                  CALL OlfGetSpe (jout,K,1,zq(K),4096,X,Y,D,Fehler)
                  IF (Fehler.ne.'&ff') RETURN

                  DO i = iSE,iEE
                     IF(rPD(i).gt.0.0) THEN
                        Y(i) = Y(i)/rPD(i)
                        D(i) = dsqrt0(Y(i)/rPD(i))
                     ELSE
                        Y(i) = 0.
                        D(i) = 0.
                        ENDIF
                     ENDDO

C --- cut data at end positions ---------------------------------
                     nE = iEE - iSE + 1
                     DO i = 1,nE
                        XE(i) = XE(i-1+iSE)
                        Y(i) = Y(i-1+iSE)
                        D(i) = D(i-1+iSE)
                        ENDDO


                  CALL OlfPutSpe (jout, K, 1, zq(K), nE, XE,Y,D,Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDDO ! K

                  ENDIF !qMon

         ELSE  ! implemented special option for test purposes

                  dBG(1)=0.000319
                  dBG(2)=0.000608
                  dBG(3)=0.000335
                  dBG(4)=0.000124
                  dBG(5)=0.0000429
                  dBG(6)=0.0000313
                  dBG(7)=0.0000296
                  dBG(8)=0.0000298
                  dBG(9)=0.0000253
                  dBG(10)=0.0000259
                  dBG(11)=0.0000243
                  dBG(12)=0.0000229
                  dBG(13)=0.0000208
                  dBG(14)=0.0000232
                  dBG(15)=0.0000282
                  dBG(16)=0.0000264



           DO K = 1, nK
                  CALL OlfGetSpe (jout,K,1,zq(K),4096,X,Y,D,Fehler)
                  IF (Fehler.ne.'&ff') RETURN



                  DO i = iSE,iEE
                     IF(rPD(i).gt.0.0) THEN
                        D(i) = dsqrt0(Y(i))/rPD(i)*rPM(i)/rMon(i)
                        Y(i) = (Y(i)/rPD(i)-dBG(K))*rPM(i)/rMon(i)
                     ELSE
                        Y(i) = 0.
                        D(i) = 0.
                        ENDIF
                     ENDDO

C --- cut data at end positions ---------------------------------
                     nE = iEE - iSE + 1
                     DO i = 1,nE
                        XE(i) = XE(i-1+iSE)
                        Y(i) = Y(i-1+iSE)
                        D(i) = D(i-1+iSE)
                        ENDDO


                  CALL OlfPutSpe (jout, K, 1, zq(K), nE, XE,Y,D,Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  ENDDO ! K




            ENDIF !qEnerg

            !Save Monitor

      CALL OlfCreate (jout+1, Kout, File, Title, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL tOlfG (jout, 'fil', File, Fehler)
      CALL tOlfG (jout, 'tit', Title, Fehler)
      IF (Fehler.ne.'&ff') GOTO 120

         CALL OlfCnuP (jout+1, 'x', 'w', 'ueV', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout+1, 'y', 'I', 'cts ', Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfCnuP (jout+1, 'z1', '2th', ' ', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         DO i=iSE,iEE
            D(i-iSE+1) = dsqrt(rMon(i))/rPM(i)
            rMon(i-iSE+1) = rMon(i)/rPM(i)
            ENDDO

      CALL OlfPutSpe (jout+1, 1, 1, 0.0d0, nE, XE,rMon,D,Fehler)  !Artem 0 -> 0.0d0
                  IF (Fehler.ne.'&ff') RETURN


 120  CONTINUE

      Close (11)

      GOTO 105

 151  CONTINUE

      END ! RRT_In_Hfbs

C     -----------------------------------------------------------
      SUBROUTINE RRT_In_Hfbso (Fehler)
C     -----------------------------------------------------------
      ! old implementation of HFBS input routine (before Aug 2003)
        ! Andreas Meyer 12/99
      ! Read converted HFBS data and normalize to monitor (FC)

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      CHARACTER*(*)  Fehler
      CHARACTER*12   Filehfbs
      CHARACTER*80   Title, File, LongTit  !Artem add LongTit

      DATA Title / ' ' /, File / ' ' /

      Fehler = '&ff' !

         Print *
         Print *, ' read energy converted HFBS data files '
         Print *


 105  Continue               ! outer loop

      CALL OlfCreate (jout, Kout, File, Title, Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL FrageTD (' run number', irun, irun)

      IF (irun.lt.0) GOTO 151

      IF (irun.lt.1 .or. irun.gt.999999) THEN
         CALL Absturz ('RRT_In_Hfbs', 'Bad number of run : '//
     *                 cl6(irun))
         ENDIF

      CALL OlfCnuP (jout, 'x', 'w', 'ueV', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL OlfCnuP (jout, 'y', 'norm. Int.', 'Cnts', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      CALL OlfCnuP (jout, 'z', 'q', 'A-1', Fehler)
      IF (Fehler.ne.'&ff') RETURN


C  The following variables are fixed :

      nK = 17   ! # of spectra: detectors + monitor
      nC = 1161 ! for Dec99 Dibut

C  Open the input files :


      CALL Compose2 (Filehfbs, 'd'//cl6(irun), '.con')

      CALL OpenFile (11, Filehfbs, '&noext', 'a', Fehler)

      IF (Fehler.ne.'&ff') THEN
        CALL Absturz ('RRT_In_Hfbs', 'Could not open file '//Filehfbs)
        ENDIF

      Read (11, '(21x,a)')     LongTit          ! sample name and comment
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(21x,f10.0)') Ncycl            ! number of cycles
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(x)')
      Read (11, '(23x,d10.0)') rMon             ! Monitor
      Read (11, '(x)')
      Read (11, '(x)')

      DO K = 1, Ncycl
       Read (11, '(x)')
       ENDDO

         Print *, ' sample: ',LongTit
         Print *, ' monitor        = ',rMon
         Print *
C  Data
      IF (MC.lt.nC) THEN
         CALL Absturz ('RRT_In_Hfbs', 'Recompile with MC >= nC')
      ENDIF

      DO K = 1,nK
      Read (11, '(21x,d10.0)') zAng             ! Angle
         zq = 2.0 * dsin(zAng / 2.0)            ! Approximation
         DO i = 1,nC
           Read(11,'(2x,d10.0,4x,d10.0)')  X(i),Y(i)
           Y(i) = Y(i) / rMon
           D(i) = sqrt(Y(i))
         ENDDO

         CALL OlfPutSpe (j, K, 1, zq, nC, X, Y, D, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO

      Close (11)

      GOTO 105

 151  CONTINUE

      END ! RRT_In_Hfbso


C ==================================================================
C   i8  (8)  Input Routine FANS (NIST BT4 sub)
C ==================================================================

      SUBROUTINE RRT_In_Fans (Fehler)
C     ----------------------------
      ! Florian Kargl 07/03
      ! Read data from instrument file.
      ! instrument files are named, you have first to
      ! transform filename to a number with 4 digits
      ! and no fileextension.
      ! read in is simple done as counts over angles
      ! z-component is incoming neutron energy
      ! in order to get S(w) one has to sum over all
      ! angles. (except channel one and 2)
      ! summing would be integrated later on with
      ! various options.
      ! Concerning the density of states calculations
      ! please be advised that energy grid is not constant !!
      ! FK 29-Jul_2009 change D1 dimension from 512 to MC

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)


      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER (MNuSu = 100)

      CHARACTER*(*)  Fehler
      CHARACTER*40   Filefans, ListFans, Fansext
      CHARACTER*80   aus, File, Title, LongTit
      CHARACTER*80   cTest

      DIMENSION rQX(512),rQY(512),rQZ(512),rCT(512),rMin(512),Z(512),
     *         irun(MNuSu)
C       ,D1(MC)
      DATA Title / ' ' /, File / ' ' /,ListFans /'-'/,qInter /.false./

      Fehler = '&ff' !

         Print *, ' read Fans data summing up over angles '
         Print *, ' Be careful: two types of raw data, choose first'

 205  CONTINUE               ! outer loop

      ListFans = ' '
      CALL GetJList ('Sample run number', ListFans,
     *               MNuSu, nirun, irun, 0, 999)
      IF (nirun.le.0) GOTO 251

      DO j = 1, nirun
       IF (irun(j).lt.1 .or. irun(j).gt.999) THEN
          CALL Absturz ('RRT_In_FANS', 'Bad number of run :
     *     '//cl6(irun(j)))
          ENDIF

       Print *,('\Please insert raw data name without extension
     *        and number:')
       Read (*,'(A)') Filefans !Artem '(a\)' -> '(A)'

       CALL OlfCreate (jout, Kout, File, Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL tOlfG (jout, 'fil', File, Fehler)
       CALL tOlfG (jout, 'tit', Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN

       CALL OlfCnuP (jout, 'x', 'ch#', ' ', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'y', 'I(2th,w)', 'Cnts', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'z1', 'E', 'meV ', Fehler)
       IF (Fehler.ne.'&ff') RETURN

       ! in second file all times will be stored versus energy
       CALL OlfCreate (jout+1, Kout, File, Title, Fehler)
       CALL OlfCnuP (jout+1, 'x', 'E', 'meV', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout+1, 'y', 'time', 'min', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout+1, 'z1', 't', 'min ', Fehler)
       IF (Fehler.ne.'&ff') RETURN

       Fansext = ' '
      ! obey all information about Xstal orientation !!
       CALL Append(Fansext, cv3(irun(j)))
       Print *, ('Filefans = '), Filefans
       Print *, ('Fansext = '), Fansext
       CALL Compose2(Filefans, Filefans, Fansext)
       CALL OpenFile (11, Filefans,'bt4','a',Fehler)
       IF (Fehler.ne.'&ff') THEN
          Print *,('Could not open file')//Filefans
          RETURN
          ENDIF
C     *    CALL Absturz('RRT_In_FANS', 'Could not open file '//Filefans)


       qInter = qAskD('Had the experiment been interrupted ?',
     *               intq(qInter))

       ! This takes care about the iPt problem: iPt is fixed at start
       ! independent of actual number of energy steps counted
       ! implementation of end of file reading is planed
       ! !! at present you have to insert a blank in the last line
       IF(qInter) THEN
C      estimating the number of lines in file
        ilin = 0
        idlin = 0
        ilin2 = 0
        rlin3 = 0.
        i = 0
        cTest = '1 try'
        DO WHILE (cTest.ne.' ')
           ilin = ilin + 1
           Read (11,'(a80)') cTest
            IF (cTest(1:4).eq.' 3,0') THEN
               IF (i.eq.0) THEN
                  ilin2 = ilin
C                  qTest = qAsk ('Should we proceed in 1 ?')
               ELSE
                  idlin = ilin - ilin2
                  ilin2 = ilin
                  rlin3 = rlin3 + idlin
                  ENDIF
               i = i + 1
            ENDIF
           ENDDO
        Close(11)
        rlin3 = rlin3/(i-1)
        Print *,('idlin = '),idlin
        Print *,('rlin3 = '),rlin3
        Print *,('If rlin3 not equal idlin: Check file !!')
        nE = (ilin - 13)/idlin
        Print *,('Fields in file: '), nE
       ENDIF


       CALL OpenFile (11, Filefans,'bt4','a',Fehler)
       IF (Fehler.ne.'&ff') THEN
          Print *,('Could not open file')//Filefans
          RETURN
          ENDIF
C     *    CALL Absturz('RRT_In_FANS', 'Could not open file '//Filefans)
       Read (11,'(43x,d9.2,x,d2.0,10x,i3/)')  rMon, rPrf, iPt
       Print *, ('rMon = '),  rMon !number of monitor counts until run stops
       Print *, ('rPrf = '),  rPrf !repetition rate at on certain energy
       Print *, ('iPt  = '),  iPt  ! number of total points over energy
       Read (11,'(a80////)')  LongTit ! Long title as given in data file
       Print *, ('Longtitle: '),LongTit
       Read (11,'(3x,d7.3,3x,d7.3,6x,d7.3////)') rES,delE,rEF
       ! rES: starting energy , delE: stepwidth in fraction of resolution
       ! rEF: final neutron energy considered
       Print *,('Start energy: '),rES
       Print *,('stepwidth E: '),delE
       Print *,('final energy: '),rEF

       !store real parameters
       CALL rOlfP(jout,'rMon','cts',rMon,Fehler)
       CALL rOlfP(jout,'dPrf','t*cts',dPrf,Fehler)
       ! scaling by summing up different files is done by
       ! rMon and iPrf taken into account.
       ! normalization to monitor is later done to a ratio
       ! of 30000 counts.

       IF (.not.qInter) nE = iPt


      ! now start to read in all fields
       DO i = 1, nE
        Read (11,'(3(d8.4),x,d8.4,3x,d6.3,5x,d6.0)') rQX(i),
     *   rQY(i),rQZ(i),Z(i),rMin(i),rCT(i)
         !rQX ... is not used, Z(i) is energy of the incoming beam
         !rMin: minutes counted for this spectrum, rCT: total
         ! number of neutrons counted in spectrum
         ! several option concerning background subtraction
         ! would be calculated next week. including the time weighting !!
         ! one has to be careful at that point about shifting of
         ! random statistics.
         DO k = 1,128
              Y(k) = 0
              ENDDO
        Read (11,*) (Y(k), k=1,128)
        DO k = 1,128
          X(k) = k-2
          ! angle doesen't matter for this type of instrument
          ! random distribution of counts over all detector numbers !
          ! have a look on spectra: random distribution -> no real
          ! q-information.
          D(k) = dsqrt(Y(k))
          ENDDO

        CALL OlfPutSpe (jout, i, 1, Z(i), 128, X , Y, D, Fehler)
          IF (Fehler.ne.'&ff') RETURN
        ENDDO

        DO i = 1, nE
            D1(i) = 0.
            ENDDO
        CALL OlfPutSpe (jout+1,1,1,10.0d0,nE, Z, rMin,D1,Fehler)  !Artem 10 -> 10.0d0
          IF (Fehler.ne.'&ff') RETURN
      Close (11)

      ENDDO !j

      GOTO 205

 251   CONTINUE

      END ! End FANS

C ==================================================================
C   i8  (9) Input Routine BT2 NIST + DMC data
C ==================================================================

      SUBROUTINE RRT_In_BT2 (Fehler)
C     ----------------------------
      ! Florian Kargl 07/03
      ! Read data from instrument file.
      ! instrument files are named, you have first to
      ! transform filename to a number with 4 digits
      ! and no fileextension.
      ! read in is simple done as counts over angles
      ! z-component is incoming neutron energy
      ! in order to get S(w) one has to sum over all
      ! angles. (except channel one and 2)
      ! summing would be integrated later on with
      ! various options.
      ! Concerning the density of states calculations
      ! please be advised that energy grid is not constant !!

      ! reads also PSI DMC data

      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)


      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER (MNuSu = 100)

      CHARACTER*(*)  Fehler
      CHARACTER*40   Filebt2, Listbt2, Bt2ext, FilePSI, ListFans !Artem add ListFans
      CHARACTER*80   aus, File, Title, LongTit
      CHARACTER*80   cTest
      CHARACTER*1    cType

      DIMENSION rQX(512),rQY(512),rQZ(512),rE(512),rT(512),rMin(512),
     *         irun(MNuSu),rQ(512),iY(512)
      DATA Title / ' ' /, File / ' ' /,ListFans /'-'/,qInter /.false./,
     *      qMon/.true./,nLoop/0/


      Fehler = '&ff' !

         Print *, ' read BT2 data'
         Print *, ' Be careful: two types of raw data, choose first'

      Print *, 'special version for PSI data !!'
      qPSI = qAsk('Read in PSI data ?')
      IF(qPSI) GOTO 405

 205  CONTINUE               ! outer loop
      nLoop = nLoop+1

      Listbt2 = ' '
      CALL GetJList ('Sample run number', Listbt2,
     *               MNuSu, nirun, irun, 0, 999)
      IF (nirun.le.0) GOTO 251

      DO j = 1, nirun
       IF (irun(j).lt.1 .or. irun(j).gt.999) THEN
          CALL Absturz ('RRT_In_BT2', 'Bad number of run :
     *     '//cl6(irun(j)))
          ENDIF

       IF (nLoop.gt.1) THEN
          qTheSame = qAskD('The rest as before ?',intq(qTheSame))
          IF(qTheSame) GOTO 230
          ENDIF

       Print *,('\Please insert raw data name without extension
     *        and number:')
       Read (*,'(A)') Filebt2  !Artem '(a\)' -> '(A)'

       qMon = qAskD('Normalization to Monitor ?',intq(qMon))

 230   CONTINUE

       ! in second file all times will be stored versus energy
C       CALL OlfCreate (jout+1, Kout, File, Title, Fehler)
C       CALL OlfCnuP (jout+1, 'x', 'E', 'meV', Fehler)
C       IF (Fehler.ne.'&ff') RETURN
C       CALL OlfCnuP (jout+1, 'y', 'time', 'min', Fehler)
C       IF (Fehler.ne.'&ff') RETURN
C       CALL OlfCnuP (jout+1, 'z1', 't', 'min ', Fehler)
C       IF (Fehler.ne.'&ff') RETURN

       Bt2ext = ' '
      ! obey all information about Xstal orientation !!
       CALL Append(Bt2ext, cv3(irun(j)))
       Print *, ('Filebt2 = '), Filebt2
       Print *, ('Bt2ext = '), Bt2ext
       CALL Compose2(Filebt2, Filebt2, Bt2ext)
       CALL OpenFile (11, Filebt2,'bt2','a',Fehler)
       IF (Fehler.ne.'&ff') THEN
          Print *,('Could not open file')//Filebt2
          RETURN
          ENDIF
C     *    CALL Absturz('RRT_In_FANS', 'Could not open file '//Filefans)


       qInter = qAskD('Had the experiment been interrupted ?',
     *               intq(qInter))



       CALL OpenFile (11, Filebt2,'bt2','a',Fehler)
       IF (Fehler.ne.'&ff') THEN
          Print *,('Could not open file')//Filebt2
          RETURN
          ENDIF

C     *    CALL Absturz('RRT_In_FANS', 'Could not open file '//Filefans)
       Read (11,'(36x,a1,6x,d9.2,x,d2.0,10x,i3/)')  cType,rMon,rPrf,iPt

       IF (cType.eq.'Q'.or.cType.eq.'E') THEN
         Print *, ('rMon = '),  rMon !number of monitor counts until run stops
         Print *, ('rPrf = '),  rPrf !repetition rate at on certain energy
         Print *, ('iPt  = '),  iPt  ! number of total lines
         Read (11,'(a80////)')  LongTit ! Long title as given in data file
         Print *, ('Longtitle: '),LongTit
         Read (11,'(3x,d7.3,3x,d7.3,6x,d7.3/)') rES,delE,rEF
         ! rES: energy transfer, delE: stepwidth in meV
         ! rEF: energy of monochromator in meV
         Print *,('energy transfer: '),rES
         Print *,('stepwidth E: '),delE
         Print *,('energy of mono: '),rEF
         Read (11,'(3(d8.4),x,3(d8.4)//)') raqx,
     *              raqy,raqz,rdqx,rdqy,rdqz

       ELSEIF (cType.eq.'I') THEN
          Print *, ('rMon = '),  rMon !number of monitor counts until run stops
          Print *, ('rPrf = '),  rPrf !repetition rate at on certain energy
          Print *, ('iPt  = '),  iPt ! number of total lines

          IF (qMon) THEN
             rMonD=rAskDMu('Base Number of Monitor to normalize on',
     *                     rMonD,1.0d0,rMon*rPrf) !Artem 1 -> 1.0d0
             rMon = rMon/rMonD
             ENDIF

          Read (11,'(a80//////////)')  LongTit
          Print *, LongTit

          ENDIF

       ! scaling by summing up different files is done by
       ! rMon and iPrf taken into account.
       ! normalization to monitor is later done to a ratio
       ! of 30000 counts.

       IF (.not.qInter) nE = iPt


      ! now start to read in all fields
! reading in constant energy scans

       IF(cType.eq.'Q') THEN
        DO i = 1, nE
         Read (11,'(3(d8.4),x,d8.4,2x,d7.3,3x,d5.2,3x,i9)') rQX(i),
     *         rQY(i),rQZ(i),rE(i),rT(i),rMin(i),iY(i)


         ! normalization to monitor
          IF(qMon) THEN
          Y(i) = iY(i)/1.
          D(i) = dsqrt0(Y(i))
          Y(i) = Y(i) / rMon / rPrf
          D(i) = D(i)/ rMon/ rPrf
          ELSE
          Y(i) = iY(i)/1.
          D(i) = dsqrt0(Y(i))
          ENDIF

         IF(rqdx.ne.0.0) THEN
            rQ(i) = rQX(i)
            rdq = 1.0
         ELSEIF (rqdy.ne.0.0) THEN
            rQ(i) = rQY(i)
            rdq = 2.0
          ELSE
            rQ(i) = rQZ(i)
            rdq = 3.0
            ENDIF
            ENDDO

       CALL OlfCreate (jout, Kout, File, Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL tOlfG (jout, 'fil', File, Fehler)
       CALL tOlfG (jout, 'tit', Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN

       CALL OlfCnuP (jout, 'x', 'q', 'A-1', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'y', 'I(q)', 'Cnts/s', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'z1', 'E', 'meV ', Fehler)
       IF (Fehler.ne.'&ff') RETURN

C !store real parameters
       CALL rOlfP(jout,'rMon','cts',rMon,Fehler)
       CALL rOlfP(jout,'dPrf','t*cts',rPrf,Fehler)
       CALL rOlfP(jout,'delE','meV',rES,Fehler)

          CALL OlfPutSpe (jout, 1, 1, rdq, nE, rQ , Y, D, Fehler)
          IF (Fehler.ne.'&ff') RETURN

! reading in constant q scans !!
        ELSEIF (cType.eq.'E') THEN
           DO i = 1, nE
             Read (11,'(3(d8.4),x,d8.4,2x,d7.3,3x,d5.2,3x,i9)') rQX(i),
     *             rQY(i),rQZ(i),rE(i),rT(i),rMin(i),iY(i)


           ! normalization to monitor
             IF(qMon) THEN
             Y(i) = iY(i)/1.
             D(i) = dsqrt0(Y(i))
             Y(i) = Y(i) / rMon / rPrf
             D(i) = D(i)/ rMon/ rPrf
             ELSE
             Y(i) = iY(i)/1.
             D(i) = dsqrt0(Y(i))
             ENDIF

             ENDDO
             rdq = 0.


       CALL OlfCreate (jout, Kout, File, Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL tOlfG (jout, 'fil', File, Fehler)
       CALL tOlfG (jout, 'tit', Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN

       CALL OlfCnuP (jout, 'x', 'E', 'meV', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'y', 'I(E)', 'Cnts/s', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'z1', 'q', 'A-1', Fehler)
       IF (Fehler.ne.'&ff') RETURN

C !store real parameters
       CALL rOlfP(jout,'rMon','cts',rMon,Fehler)
       CALL rOlfP(jout,'dPrf','t*cts',rPrf,Fehler)
       CALL rOlfP(jout,'QX',' ',raqx,Fehler)
       CALL rOlfP(jout,'QY',' ',raqy,Fehler)
       CALL rOlfP(jout,'QZ',' ',raqz,Fehler)

       CALL OlfPutSpe (jout, 1, nK, rdq, nE, rE , Y, D, Fehler)
       IF (Fehler.ne.'&ff') RETURN


       ELSEIF (cType.eq.'I') THEN !only diffraction 5 angles fixed

          DO i = 1, nE
             Read (11,'(4x,d7.3,x,d8.4,2x,d5.2,3x,i9)') rQX(i),
     *          rT(i),rMin(i),iY(i)

           ! normalization to monitor
             IF(qMon) THEN
             Y(i) = iY(i)/1.
             D(i) = dsqrt0(Y(i))
             Y(i) = Y(i) / rMon / rPrf
             D(i) = D(i)/ rMon/ rPrf
             ELSE
             Y(i) = iY(i)/1.
             D(i) = dsqrt0(Y(i))
             ENDIF

             ENDDO
             rdq = 0.

             File = Filebt2
             Title = Longtit
             CALL OlfCreate (jout, Kout, File, Title, Fehler)
             IF (Fehler.ne.'&ff') RETURN
             CALL tOlfG (jout, 'fil', File, Fehler)
             CALL tOlfG (jout, 'tit', Title, Fehler)
             IF (Fehler.ne.'&ff') RETURN

             CALL OlfCnuP (jout, 'x', '2th', ' ', Fehler)
             IF (Fehler.ne.'&ff') RETURN
             CALL OlfCnuP (jout, 'y', 'I(2th)', 'Cts', Fehler)
             IF (Fehler.ne.'&ff') RETURN
             CALL OlfCnuP (jout, 'z1', 'E', 'meV', Fehler)
             IF (Fehler.ne.'&ff') RETURN

C !store real parameters
             CALL rOlfP(jout,'rMon','cts',rMon,Fehler)
             CALL rOlfP(jout,'dPrf','t*cts',rPrf,Fehler)
C             CALL rOlfP(jout,'QX',' ',raqx,Fehler)
C             CALL rOlfP(jout,'QY',' ',raqy,Fehler)
C             CALL rOlfP(jout,'QZ',' ',raqz,Fehler)

             CALL OlfPutSpe (jout, 1, 1, 1.0d0, nE, rQX , Y, D, Fehler) !Artem 1 -> 1.0d0
             IF (Fehler.ne.'&ff') RETURN




       ENDIF

          Close (11)

      ENDDO !j

      GOTO 205

 405  Print *, 'PSI data input programm DMC'

         Print *,('\Please insert raw data name')
       Read (*,'(A)') FilePSI  !Artem '(a\)' -> '(A)'
       Print *, 'FilePSI = ', FilePSI
       CALL OpenFile (11, FilePSI,'dat','a',Fehler)
       IF (Fehler.ne.'&ff') RETURN

       Read (11,'(a4)') FilePSI
       Print *, FilePSI
       IF (FilePSI(1:4).ne.'HRPT') THEN
          Print *, 'Bin 1'
       Read (11,'(x)')
       Read (11,'(x,d7.3,x,d7.3,x,d7.3)') rtanf,rtstep,rtend
       Print *, 'rtanf rstep rtend ',rtanf,rstep,rtend

       n = ((rtend-rtanf)/rtstep + 1)/10
       Print *, 'n = ',n

       DO i = 1,n
          Read (11,'(10(d8.0))') Y(i*10-9),Y(i*10-8),Y(i*10-7),
     *                           Y(i*10-6),Y(i*10-5),Y(i*10-4),
     *                           Y(i*10-3),Y(i*10-2),Y(i*10-1),
     *                           Y(i*10)
          ENDDO
       DO i = 1,n
          Read (11,'(10(d8.0))') D(i*10-9),D(i*10-8),D(i*10-7),
     *                           D(i*10-6),D(i*10-5),D(i*10-4),
     *                           D(i*10-3),D(i*10-2),D(i*10-1),
     *                           D(i*10)
          ENDDO
       nE = 10*n
       DO i = 1,nE
          X(i) = rtanf + (i-1)*rtstep
C          D(i) = dsqrt0(Y(i))
          ENDDO

       ELSE

          Read (11,'(x)')
          Read (11,'(x,d7.3,x,d7.3,x,d7.3)') rtanf,rtstep,rtend
          Print *, 'rtanf rstep rtend ',rtanf,rstep,rtend

          n = ((rtend-rtanf)/rtstep + 1)/10
          Print *, 'n = ',n

          DO i = 1,n
             Read (11,'(10(d8.0))') Y(i*10-9),Y(i*10-8),Y(i*10-7),
     *                              Y(i*10-6),Y(i*10-5),Y(i*10-4),
     *                              Y(i*10-3),Y(i*10-2),Y(i*10-1),
     *                              Y(i*10)
             ENDDO
          Read (11,'(8(d8.0))')  Y(n*10+1),Y(n*10+2),Y(n*10+3),
     *                           Y(n*10+4),Y(n*10+5),Y(n*10+6),
     *                           Y(n*10+7),Y(n*10+8)

          DO i = 1,n
             Read (11,'(10(d8.0))') D(i*10-9),D(i*10-8),D(i*10-7),
     *                              D(i*10-6),D(i*10-5),D(i*10-4),
     *                              D(i*10-3),D(i*10-2),D(i*10-1),
     *                              D(i*10)
             ENDDO
          Read (11,'(8(d8.0))')  D(n*10+1),D(n*10+2),D(n*10+3),
     *                           D(n*10+4),D(n*10+5),D(n*10+6),
     *                           D(n*10+7),D(n*10+8)
          nE = 10*n+8
          DO i = 1,nE
             X(i) = rtanf + (i-1)*rtstep
C            D(i) = dsqrt0(Y(i))
             ENDDO

          ENDIF

       CALL OlfCreate (jout, Kout, File, Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL tOlfG (jout, 'fil', File, Fehler)
       CALL tOlfG (jout, 'tit', Title, Fehler)
       IF (Fehler.ne.'&ff') RETURN

       CALL OlfCnuP (jout, 'x', '2th', ' ', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'y', 'I(2th)', 'Cnts', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfCnuP (jout, 'z1', 'E', 'meV ', Fehler)
       IF (Fehler.ne.'&ff') RETURN
       CALL OlfPutSpe (jout, 1, 1, 1.0d0, nE , X , Y, D, Fehler) !Artem 1 -> 1.0d0

       Close (11)

 251   CONTINUE

      END ! End BT2

C================================================================
C  i8  (10)              Input Routine IN3 (ILL)
C================================================================

C ---------------------------------------------------------------
      SUBROUTINE RRT_In_IN3 (Fehler)
C ---------------------------------------------------------------
C     FK Jan 2005 at the moment only isotropic scatterer accepted
C     FK Mar 2006 still preliminary
C
      IMPLICIT REAL*8   (a-h,o-p,r-z)
      IMPLICIT LOGICAL  (q)


      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'

      PARAMETER (MNuSu = 100)

      CHARACTER*(*)  Fehler
      CHARACTER*40   Filein3, Listin3
      CHARACTER*80   aus, File, Title, LongTit
      CHARACTER*80   cTest
      CHARACTER*1    cType

      DIMENSION rQX(512),rQY(512),rQZ(512),rE(512),rT(512),rMin(512),
     *         irun(MNuSu),rQ(512),iY(512)
      DATA Title / ' ' /, Filein3 / ' ' /,Listin3 /'-'/,
     *          qInter /.false./, qMon/.true./,nLoop/0/

      PRINT *, 'read data of IN3 and store S(2th,w) for isotropic'

 205  CONTINUE               ! outer loop
      nLoop = nLoop+1

      Listin3 = ' '
      CALL GetJList ('Sample run number', Listin3,
     *               MNuSu, nirun, irun, 0, 999999)
      IF (nirun.le.0) GOTO 251

      DO j = 1, nirun
       IF (irun(j).lt.1 .or. irun(j).gt.999999) THEN
          CALL Absturz ('RRT_In_IN3', 'Bad number of run :
     *     '//cl6(irun(j)))
          ENDIF
          ENDDO
       IF (nLoop.gt.1) THEN
          qTheSame = qAskD('The rest as before ?',intq(qTheSame))
          IF(qTheSame) GOTO 230
          ENDIF

       qMon = qAskD('Normalization to Monitor ?',intq(qMon))


 230   CONTINUE


 251  CONTINUE

      END ! End IN3
