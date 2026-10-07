C This file contains IRIS specific subroutines
C The different programs are written by F.K.
C IRIS Data read in is based on WS Howells Fortran
C Routines for the IRIS instrument implemented in the
C ISIS Modes Package.
C
C ----------------------------------------------------
C              History
C ----------------------------------------------------
C   May06 F. Kargl 1st built
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C
C ----------------------------------------------------
C  1. Mean square displacement calculation
C  2. Data read in IRIS at ISIS
C
C
C
C
C *****************************************************
C 1.) MSD calculus
C -----------------------------------------------------
C
C -----------------------------------------------------
      SUBROUTINE MSDCalc(Fehler)
C -----------------------------------------------------
      ! FK 23may06
      ! MSD are taken from OpenGenie data


      IMPLICIT REAL*8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      INCLUDE 'i_wrk.f'


      PARAMETER (MNuSu=100)
      CHARACTER*80 ListMD, ListBGMD,DirRaw,FileRaw
      CHARACTER*40 Title, FilOut, Inst
      INTEGER    iTemp

      CHARACTER Fehler*(*)


      DIMENSION JNumorBGMD(MNuSu), JNumorMD(MNuSu),iTemp(100),
     *          Temp(100)

C  ------------------------------------------------------
C     Default varables
C  ------------------------------------------------------
      DATA qBGsub /.false./
C  ------------------------------------------------------

C  Open log file
      CALL OpenFile (30, 'msd', 'log', 'e', Fehler)
       IF (Fehler.ne.'&ff') GOTO 99


C  Begin dialogue for read in of the data:
C  --------------------------------------------------

      !ask for path variable

 110  CONTINUE
       nLoop = nLoop + 1

      ListMD = ' '
      CALL GetJList ('MSD run numbers', ListMD, MNuSu,
     *               nNumorMD, JNumorMD, 0, 99999)
      IF (nNumorMD.le.0) GOTO 99

      qBGsub = qAskD('Background subtraction?', intq(qBGsub))
      IF(qBGsub) THEN
         CALL GetJList ('MSD background run numbers', ListBGMD, MNuSu,
     *                  nNumorBGMD, JNumorBGMD, 0, 99999)
         ENDIF


      iTemps = iAskD ('Start temperature', iTemps)
      iTempi = iAskD ('Temperature increment',iTempi)

      DO i=1,nNumorMD
         iTemp(i) = iTemps + (i-1)*iTempi
         ENDDO
      DO i=1,nNumorMD
         Temp(i) = iTemp(i)/1.0
         ENDDO

      WRITE (35,'(/a)') 'going to read sample runs'
      Print *, 'reading in the msd files'


      CALL ExeML ('\p dir-raw-n-dac', DirRaw)

      CALL OlfCreate (j,Kout, '&noask','&noask', Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

      CALL OlfCnuP (j, 'x', 'q^2', 'A-2', Fehler)
      CALL OlfCnuP (j, 'y', 'S(q,w=0)', ' ', Fehler)
      CALL OlfCnuP (j, 'z1', 'T', 'K', Fehler)

      DO K = 1, nNumorMD

         FileRaw = DirRaw

      CALL Append (FileRaw, 'irs')
      Print *,'step1 dir raw: ', FileRaw
      CALL Append (FileRaw, cv5(JNumorMD(K)))
      Print *,'step2 dir raw: ', FileRaw

      CALL OpenFile (11,FileRaw, 'msd', 'a', Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO i = 1,51
C      READ (11, '(3(10f8.0))') X(i),Y(i),D(i)
      READ (11, *) X(i),Y(i),D(i)
       ENDDO

       Print *, 'iTemp(',K,'): ', iTemp(K), Temp(K)
       CALL OlfPutSpe (j,K,1,Temp(K),51,X,Y,D,Fehler)

       Close(11)
      ENDDO !K

      CALL FrageCD ('File name: ', FilOut, 'irs'//cl5(JNumorMD(1)))
      CALL FrageCD ('And Title', Title, Title)
      CALL tOlfP (j, 'fil', FilOut, Fehler)
      CALL tOlfP (j, 'tit', Title, Fehler)
      CALL tOlfP (j, 'sub', 'IRIS', Fehler)


      GOTO 110

      Close(30)

 99    Write (30,'(a)') 'end of read in'
      Close (30)


      END !MSD Calc


C *****************************************************
C 2.) Data read in IRIS at ISIS
C -----------------------------------------------------
C
C -----------------------------------------------------
      SUBROUTINE RRaw_IRS(Inst, Fehler)
C -----------------------------------------------------
      ! FK 05jun06

      IMPLICIT REAL*8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      INCLUDE 'i_wrk.f'

      PARAMETER   (MNuSu=100)

      DIMENSION   AngleS(MK),XepV(MK),YepV(MK)

      CHARACTER*1 BadS(MK), BadV(MK), BadF(MK)

      CHARACTER*(*) Inst, Fehler

      CHARACTER*40 FilOut, Title
      CHARACTER*80 ListSD,ListSC,ListVD,ListVC,LKnoAnal,LisKdel

      DIMENSION    JNumorSD(MNuSu), JNumorSC(MNuSu),
     *             JNumorVD(MNuSu), JNumorVC(MNuSu)


C  Open Log-File :
      CALL OpenFile (35, 'inzrev', 'log', 'e', Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

      iCycle = iAskD (
     * 'Data from instrument(-1), current cycle(0), archive(>0)',
     * iCycle)
      iNOff = iAskD ('Run number offset', iNOff)

 110  CONTINUE !outer loop
      nLoop = nLoop + 1

      ListSD = ' '
      CALL GetJList('Sample run numbers', ListSD, MnuSu, nNumorSD,
     *                JNumorSD, 0, 99999)

      IF (nNumorSD.le.0) GOTO 99
      DO i = 1, nNumorSD
         JNumorSD(i) = JNumorSD(i) + iNOff
         ENDDO

      IF (nLoop.gt.1) THEN
         qTheSame = qAskD('All the rest as before', intq(qTheSame))
         IF (qTheSame) GOTO 200
      ELSE
         qTheSame = .false.
         ENDIF

      qMon = qAskD ('Normalization to monitor', intq(qMon))

C Ask for output of data options are time in usec,
C wavelength in A or energy transfer in meV.
      Print *
      Print *, 'Sample treatment :'
      Print *, ' '
      Print *, 'Output of data: Time(0), WaveL(1), Energy(2)'
      Print *, ' '

      iOmod = iAskDMu ('Choose option', iOmod, 0, 2)

      IF (iOmod.eq.1.or.iOmod.eq.2) THEN
        iSEG = -1 ! neutron energy gain is positive
        qDetEff = qAskD('Detector efficiency correction',
     *                  intq(qDetEff))
       ELSE
         iSEG    = -1 ! not used
         qDetEff = .false.
         ENDIF

      IF (qDetEff) THEN
         ListVD = ' '
         CALL GetJList('Vanadium run numbers', ListVD, MnuSu, nNumorVD,
     *                 JNumorVD, 0, 99999)

         IF (nNumorVD.le.0) THEN
            GOTO 99
            ELSE
             VTemp = 285.0
             VTemp = rAskDMu('Vanadium temperature:',VTemp,
     *                       285.d0,310.d0)
            ENDIF

         DO i = 1, nNumorVD
            JNumorVD(i) = JNumorVD(i) + iNOff
            ENDDO
         ENDIF

C choose now option how to treat the frame
C (0) means regular no monitor unwrap
C (1) means monitor is unwraped
      IF (qMon) THEN
         Print *
         Print *, 'Handling of sample frame:'
         Print *, '   (0) do nothing'
         Print *, '   (1) unwrap sample to monitor'

         iFmod = iAskDMu ('Choose option', iFmod, 0, 1)
         ENDIF

C choose type of analyser
      Print *, ' '
      Print *, '---------------------'
      Print *, 'Choose analyser bank:'
      Print *, '(0) PG bank'
      Print *, '(1) Mica'
      Print *, '(2) both'
      Print *, '---------------------'

      iAnaTyp = iAskDMu ('Which type of analysers', iAnaTyp, 0, 2)

      IF (iAnaTyp.eq.0) THEN
         Print *, ' '
         Print *, '-------------------'
         Print *, '(0) 002 reflection'
         Print *, '(1) 004 reflection'
         Print *, '-------------------'
         iRefTyp = iAskDMu ('Which reflection',iRefTyp,0,1)
      ELSEIF (iAnaTyp.eq.1) THEN
                  Print *, ' '
         Print *, '-------------------'
         Print *, '(0) 002 reflection'
         Print *, '(1) 004 reflection'
         Print *, '(2) 006 reflection'
         Print *, '-------------------'
         iRefTyp = iAskDMu ('Which reflection',iRefTyp,0,2)
      ELSE
         iRefTyp = 3
         ENDIF

C choose output whether conform to OpenGenie or exact

      IF (iOmod.ne.0) THEN
         Print *, ' '
         qMod = qAskD('Output a la OpenGenie (non exact treatment)?',
     *               intq(qMod))
         ENDIF

C      CALL GetNList ('Exclude spectra from analysis (monitor bad
C     *                detectors, banks)', LKnoAnal, qKnoAnal, MK)

 200  CONTINUE

      Write (35,'(/3a)') 'going to read data from instrument ', Inst,
     *     ', cycle '//cl6(iCycle)

C  --------------------------------------------------------------------
C  Get sample data / normalize and correct them :
C  --------------------------------------------------------------------



      IF (qDetEff) THEN
         Print *
         Print *, 'getting Vanadium runs'

         CALL RRaw_IRSSub(Inst,iCycle,nNumorVD,JNumorVD,j,nK,XepV,
     *   BadV,iAnaTyp,iRefTyp,.false.,0,.false.,0,Fehler)
         CALL RRaw_IRSVan(j,nK,VTemp,XepV,YepV,iAnaTyp,iRefTyp,Fehler)
         CALL MemFileDel(j,Fehler)
         Print *, 'Vanadium succesffuly read: ', YepV(1)
         IF (Fehler.ne.'&ff') GOTO 99
      ELSE
         Print *, 'No Vanadium no detector efficiency'
         ENDIF


      Print *
      Print *, 'getting sample runs :'

      Write (35,*) 'going to read sample run(s):'

      CALL RRaw_IRSSub(Inst,iCycle,nNumorSD,JNumorSD,j,nK,AngleS,BadS,
     *                 iAnaTyp,iRefTyp,qMod,iOmod,qMon,iFmod,Fehler)
C      Print *, 'j = ',j,nNumorSD
      IF (Fehler.ne.'&ff') GOTO 99
      IF (qDetEff) THEN
         CALL RRaw_IRSNor(j,nNumorSD,nK,XepV,YepV,Fehler)
         ENDIF


C  Outermost loop: next run ?
      Print *, ' '
      Print *, 'Next data set:'
      qTheSame = .true.
      GOTO 110


 99   Close (35)


       END !RRaw_IRS

C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSSub(Inst,iCycle,nNumor,JNumor,j,nNOS,Angle,
     *                      Status,iAna,iRef,qMod,iOmod,qMon,iFmod,
     *                      Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     JNumor(*)
      DIMENSION     Angle(*)
      REAL*8        AngleM(MK)
      REAL*8        rAng(MK),Xout(2001)
      INTEGER*4     Yout(2001*115)

      CHARACTER*(*) Inst, Fehler
      CHARACTER*1   Status(*)
      CHARACTER*40  FilOut, Title
      CHARACTER*80  FileRaw,DirRaw, LongTit, LongDat

      INTEGER*4     nTc,nK,nUSE,iFehler


C  decide whether to add all files (case non MSD)
C  or not to add files (MSD calculation)
      IF (nNumor.gt.1) THEN
         qAdd = qAskD('Add all files to a common one?', intq(qAdd))
      ELSE
         qAdd = .false.
         ENDIF

C  Open on-line file :
      DO iol = 1,nNumor

      IF (.not.qAdd) THEN
         IF (qMon) THEN
            CALL RRaw_IRSIn(Inst,iCycle,JNumor(iol),j,nNOS,Angle,
     *                   Status,iAna,iRef,qMod,iOmod,Fehler)
            IF (iFmod.eq.1) THEN
               Write (35,*) 'Parameter after Read in: (nNOS)',nNOS
               CALL RRaw_IRSuns (j,nNOS,iAna,iRef,iOmod,Fehler)
               ENDIF
            nNOS1 = nNOS
C  default output filename :
            CALL FrageCD ('File name', FilOut, 'r'//cl6(JNumor(iol)))
            CALL FrageCD ('And title', Title, Title)
            CALL tOlfP (j, 'fil', FilOut, Fehler)
            CALL tOlfP (j, 'tit', Title,  Fehler)
            CALL tOlfP (j, 'sub', Inst , Fehler)

            CALL OlfClos (j, nNOS, Fehler)
            CALL RRaw_IRSIn(Inst,iCycle,JNumor(iol),j1,nNOS,AngleM,
     *                   Status,3,iRef,qMod,iOmod,Fehler)
            Title = 'monitor of r'
            CALL Append (Title,cl6(JNumor(iol)))
            CALL Append (FilOut,'m')
            CALL tOlfP (j1, 'fil', FilOut, Fehler)
            CALL tOlfP (j1, 'tit', Title,  Fehler)
            CALL tOlfP (j1, 'sub', Inst , Fehler)
            CALL OlfClos (j1, nNOS, Fehler)

            IF (iOmod.eq.1) THEN
               CALL RRaw_IRSmon(j,j1,nNOS1,Fehler)
               CALL MemFileDel (j,Fehler)
               CALL MemFileDel (j,Fehler)
               CALL MemFileDel (j,Fehler)
               IF (qAskD('Output in energy transfer', intq(qET))) THEN
                  CALL RRaw_IRSmcw(j,nNOS1,iAna,iRef,Fehler)
                  ENDIF
               nNOS = nNOS1
               ENDIF

         ELSE
            CALL RRaw_IRSIn(Inst,iCycle,JNumor(iol),j,nNOS,Angle,
     *                   Status,iAna,iRef,qMod,iOmod,Fehler)
C  default output filename :
            CALL FrageCD ('File name', FilOut, 'r'//cl6(JNumor(iol)))
            CALL FrageCD ('And title', Title, Title)
            CALL tOlfP (j, 'fil', FilOut, Fehler)
            CALL tOlfP (j, 'tit', Title,  Fehler)
            CALL tOlfP (j, 'sub', Inst , Fehler)

            CALL OlfClos (j, nNOS, Fehler)

            ENDIF !qMon

      ELSE ! add spectra

         IF (qMon) THEN
            CALL RRaw_IRSIn(Inst,iCycle,JNumor(iol),j,nNOS,Angle,
     *                   Status,iAna,iRef,qMod,iOmod,Fehler)
            IF (iFmod.eq.1) THEN
               Write (35,*) 'Parameter after Read in: (nNOS)',nNOS
               CALL RRaw_IRSuns (j,nNOS,iAna,iRef,iOmod,Fehler)
               ENDIF
            nNOS1 = nNOS

            CALL OlfClos (j, nNOS, Fehler)
            CALL RRaw_IRSIn(Inst,iCycle,JNumor(iol),j1,nNOS,AngleM,
     *                   Status,3,iRef,qMod,iOmod,Fehler)

            CALL OlfClos (j1, nNOS, Fehler)
            IF (iOmod.eq.1) THEN
               CALL RRaw_IRSmon(j,j1,nNOS1,Fehler)
               CALL MemFileDel (j,Fehler)
               CALL MemFileDel (j,Fehler)
               CALL MemFileDel (j,Fehler)
               IF (qAskD('Output in energy transfer', intq(qET))) THEN
                  CALL RRaw_IRSmcw(j,nNOS1,iAna,iRef,Fehler)
                  ENDIF
               nNOS = nNOS1
               IF(iol.gt.1) THEN
                  CALL RRaw_IRSAdd(j-1,j,nNOS1,Fehler)
                  CALL MemFileDel (j,Fehler)
                  ENDIF
            ELSE
               IF(iol.gt.1) THEN
                  CALL RRaw_IRSAdd(j-2,j,nNOS1,Fehler)
                  CALL RRaw_IRSAdd(j1-2,j1,nNOS1,Fehler)
                  CALL MemFileDel(j,Fehler)
                  CALL MemFileDel(j,Fehler)
                  ENDIF
               ENDIF

         ELSE
            CALL RRaw_IRSIn(Inst,iCycle,JNumor(iol),j,nNOS,Angle,
     *                   Status,iAna,iRef,qMod,iOmod,Fehler)
            CALL OlfClos (j, nNOS, Fehler)

            IF(iol.gt.1) THEN
               CALL RRaw_IRSAdd(j-1,j,nNOS,Fehler)
               CALL MemFileDel(j,Fehler)
               ENDIF
            ENDIF !qMon

C case definition for sample description
            IF (qMon.and.iOmod.ne.1.and.iol.eq.nNumor) THEN
C  default output filename :
               CALL OlfOpen(j-1,1,iK,Fehler)
               CALL OlfOpen(j,1,iK1,Fehler)
               CALL RRaw_IRSAdN(j-1,nNOS1,nNumor,Fehler)
               CALL RRaw_IRSAdN(j,nNOS1,nNumor,Fehler)
               CALL FrageCD ('File name', FilOut, 'r'//cl6(JNumor(1)))
               CALL FrageCD ('And title', Title, Title)
               CALL tOlfP (j-1, 'fil', FilOut, Fehler)
               CALL tOlfP (j-1, 'tit', Title,  Fehler)
               CALL tOlfP (j-1, 'sub', Inst , Fehler)
               Title = 'monitor of'
               CALL Append (FilOut,'m')
               CALL Append (Title,FilOut)
               CALL tOlfP (j, 'fil', FilOut, Fehler)
               CALL tOlfP (j, 'tit', Title,  Fehler)
               CALL tOlfP (j, 'sub', Inst , Fehler)
               CALL OlfClos(j-1, nNOS1,Fehler)
               CALL OlfClos(j,nNOS1,Fehler)
               j=j-1
               nNumor=1
            ELSEIF (iol.eq.nNumor) THEN
               CALL OlfOpen(j-1,1,iK,Fehler)
               CALL RRaw_IRSAdN(j-1,nNOS1,nNumor,Fehler)
               CALL FrageCD ('File name', FilOut, 'r'//cl6(JNumor(1)))
               CALL FrageCD ('And title', Title, Title)
               CALL tOlfP (j-1, 'fil', FilOut, Fehler)
               CALL tOlfP (j-1, 'tit', Title,  Fehler)
               CALL OlfClos(j-1, nNOS1,Fehler)
               j=j-1
               nNumor=1
               ENDIF ! description of file !!


         ENDIF !qAdd

      ENDDO !iol


      END !RRaw_IRSSub

C  -------------------------------------------------------------
      SUBROUTINE RRaw_IRSIn(Inst,iCycle,iJNumor,j,nNOS,Angle,
     *                      Status,iAna,iRef,qMod,iOmod,Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

C      DIMENSION     JNumor(*)
      DIMENSION     Angle(*)
      REAL*4        rAng(MK),Xout(2001)
      INTEGER*4     Yout(2001*115)

      CHARACTER*(*) Inst, Fehler
      CHARACTER*1   Status(*)
      CHARACTER*80  FileRaw,DirRaw, LongTit, LongDat

      INTEGER*4 nTc,nK,nUSE,iFehler,iJNumor

C  Determine number of spectra to read depending
C  on chosen analyser
      IF (iAna.eq.0) THEN
         nNOS = 51
         iFIR = 3
         Eelast = 1.845
      ELSEIF (iAna.eq.1) THEN
         nNOS = 51
         iFIR = 54
      ELSEIF (iAna.eq.2) THEN
         nNOS = 113
      ELSE
         nNOS = 1
         iFIR = 1
         ENDIF
      nMN = 2001


      CALL OlfCreate (j, Kout, '&noask', '&noask', Fehler)
      IF (Fehler.ne.'&ff') GOTO 99

C  Get pathname from setup file :
      IF     (iCycle.eq.-1) THEN
         CALL ExeML ('\p dir-raw-n-dac', DirRaw)
      ELSEIF (iCycle.eq. 0) THEN
         CALL ExeML ('\p dir-raw-n-new', DirRaw)
      ELSE
         CALL ExeML ('\p dir-raw-n-old', DirRaw)
         ENDIF

C      Print *, DirRaw, iJ
      FileRaw = DirRaw

      CALL Append (FileRaw, 'irs')
      CALL Append (FileRaw, cv5(iJNumor))
      IF (iCycle.ne.0) THEN
      CALL Append (FileRaw, '.raw')
      ELSE
      CALL Append (FileRaw, '.sav')
      ENDIF
      Print *, 'Opening File :', FileRaw
C IFIR start detector NOS number of detectors used

      CALL OPEN_DATA_FILE(FileRaw,nTC,nK,nUSE,iFehler)
C      PRINT *, iFehler, nTC, nDET, nUSE
       IF (iFehler.ne.0) THEN
          Print *, 'Error (',iFehler,') in OPEN DATA FILE'
          GOTO 99
          ELSE
             Fehler = '&ff'
          ENDIF
      Print *, 'successfull opening file'

        CALL GETPARR(FileRaw,'TCB1',Xout,nMN,nMOUT,iFehler)
        Print *, 'Par 1: ', nMN, nMOUT, iFehler
        CALL GETPARR(FileRaw,'TTHE',rAng,nMN,nMOUTA,iFehler)
        Print *, 'Par 2: ', rAng(3), nMN, nMOUTA, iFehler
        CALL GETDAT(FileRaw, iFIR, nNOS, Yout, ILENGTH,iFehler)
        Print *, 'Par 3: ', iFIR, nNOS, ILENGTH, iFehler
        IF (iFehler.ne.0) THEN
          Print *, 'Error (',Fehler,') in READING DATA'
          GOTO 99
          ELSE
             Fehler = '&ff'
          ENDIF
      Print *, 'sucessfull extracting parameters'

C write data to IDA workspace
      DO i = 1,nNOS
         Angle(i) = rAng(i+2)
         ENDDO
      DO i = 1,nMN
         X(i) = Xout(i)
         ENDDO
      IF (nNOS.eq.1) THEN
         Angle(1) = rAng(1)
         ENDIF

       DO i = 1,nNOS
          DO ii = 1,2001
             Y(ii) = 0.000
C             WRITE (35,*) 'Yout(',(i-1)*2001+ii,')=',Yout((i-1)*2001+ii)
             Y(ii) = Yout((i-1)*2001+ii)
             D(ii) = dsqrt(Y(ii))
             ENDDO
             CALL OlfPutSpe (j,i,1,Angle(i),nMOUT,X,Y,D,Fehler)
C             CALL OlfPutXYD (j, nK, nMOUT, X, Y, D, Fehler)
             IF (Fehler.ne.'&ff') GOTO 99
       ENDDO

       CALL RRaw_IRSOut (Inst,j,nNOS,iOmod,iAna,iRef,qMod,
     *                   Eelast,Fehler)
C      CALL OlfComAdd (j, ' ', Inst(1:lenU(Inst))//': '//LongDat, Fehler)
C      CALL OlfComAdd (j, ' ', '"'//LongTit(1:lenU(LongTit))//'"', Fehler)

      CALL iOlfP (j, '?det-bal-sym',   0, Fehler)
      CALL iOlfP (j, '@sam-erg-gain', -1, Fehler)
      CALL iOlfP (j, 'plot-sy#',  0, Fehler)

      CALL rOlfP (j, 'E0', 'meV', Eelast, Fehler)
C      CALL rOlfP (j, 'T',  'K',   Temp,   Fehler)

      IF (iOmod.eq.0) THEN
         CALL OlfCnuP (j, 'x', 'tof', 'usec', Fehler)
         CALL OlfCnuP (j, 'y', 'I(2th,tof)', 'cts', Fehler)
      ELSEIF (iOmod.eq.1) THEN
         CALL OlfCnuP (j, 'x', 'l', 'A', Fehler)
         CALL OlfCnuP (j, 'y', 'I(2th,l)', 'cts/A', Fehler)
      ELSE
         CALL OlfCnuP (j, 'x', 'w', 'meV', Fehler)
         CALL OlfCnuP (j, 'y', 'S(2th,w)', 'meV-1', Fehler)
         ENDIF

      CALL OlfCnuP (j, 'z1', '2th', ' ', Fehler)


 99   CONTINUE


      END !RRaw_IRSIn

C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSOut(Inst,j,nK,iOmod,iAn,iRef,qMod,Eel,Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'


      CHARACTER*(*) Inst, Fehler
      PARAMETER    (HBPAR1=3.9560345693E-3,HBPAR2=81.804251686,
     *              FPPAR1=36.41,FPPAR2=1.45,FPPAR3=1.47,
     *              FPPAR4=-0.37,FPPAR5=2.12,FPPAR6=0.8576,
     *              DEPAR1=1.276E-3,DEPAR2=0.025)

C      CHARACTER*1   Status(*)

C      Print *, 'Parameters subout: '
C      Print *, 'iAn: ',iAn,'iRef: ',iRef,'iOmod: ', iOmod
C      Print *, 'qMod = ',qMod

      Write (35,*) 'Parameter in IRIS Sub'
      Write (35,*) 'HBPAR1: ', HBPAR1, 'HBPAR2: ', HBPAR2
      Write (35,*) 'FPPAR1: ', FPPAR1, 'FPPAR2: ', FPPAR2,
     *             'FPPAR3: ', FPPAR3
      Write (35,*) 'FPPAR4: ', FPPAR4, 'FPPAR5: ', FPPAR5,
     *             'FPPAR6: ', FPPAR6

C determine the final energies (meV)(s. dialogue.gcl in genie modes)
      IF (iAn.eq.0) THEN
         IF (iRef.eq.0) THEN
            rFPA = 1.845
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR2,HBPAR1)
C 2440.606014207
C            Print *, '1, rFPA, tFPA: ', rFPA, tFPA
         ELSEIF (iRef.eq.1) THEN
            rFPA = 7.3812
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR2,HBPAR1)
C            Print *, '2, rFPA, tFPA: ', rFPA, tFPA
            ENDIF
      ELSEIF (iAn.eq.1) THEN
         IF (iRef.eq.0) THEN
            rFPA = 0.2067
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
C            Print *, '3, rFPA, tFPA: ', rFPA, tFPA
         ELSEIF (iRef.eq.1) THEN
            rFPA = 0.8255
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
C            Print *, '4, rFPA, tFPA: ', rFPA, tFPA
         ELSEIF (iRef.eq.2) THEN
            rFPA = 1.8567
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
C            Print *, '5, rFPA, tFPA: ', rFPA, tFPA
            ENDIF
      ELSEIF (iAn.eq.2.and.iRef.eq.3) THEN
            rFPA = 1.845
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR2,HBPAR1)
            rFPA2 = 1.8567
            tFPA2 = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
      ELSEIF (iAn.eq.3) THEN
            GOTO 103
            Write (35,*) 'iAnalyzer = ', iAn
         ENDIF

         Eel = rFPA
C         Print *, 'Eel: ', Eel

      IF (.not.qMod) THEN
         GOTO 101
      ELSE
         GOTO 102
         ENDIF

C correct energy axis output
 101  IF (iOmod.eq.0) THEN
         GOTO 199
      ELSEIF (iOmod.eq.1) THEN
         CALL OlfGetXYD (j,1,n,X,Y,D,Fehler)
C         Print *, 'xtest X(1),X(599),X(n)',X(1),X(599),X(n)
C         Print *, 'test dquot0', dquot0((X(599)-tFPA),FPPAR1)
C         Print *, 'test multiplication: ', dquot0((X(599)-tFPA),
C     *            FPPAR1)*HBPAR1
         DO i = 1,n
           X1(i) = dquot0((X(i)-tFPA),FPPAR1)*HBPAR1
C 0.000108652
C 0.000109768
           ENDDO
C         Print *, 'test 1 out, X1(599)', X1(599)
         DO K = 1, nK
           CALL OlfGetXYD (j,K,n,X,Y,D,Fehler)
           dX = X(2)-X(1)
           dL = dsqrt(dquot0(HBPAR2,rFPA))
           DO i = 2, n
              dX1  = X1(i)-X1(i-1)
              Y(i-1) = (Y(i)*dX)/(dX1*DEPAR1*(1.-dexp(-8.3*DEPAR2*dL)))
              D(i-1) = (D(i)*dX)/(dX1*DEPAR1*(1.-dexp(-8.3*DEPAR2*dL)))
              ENDDO
           CALL OlfPutXYD (j,K,n-1,X1,Y,D,Fehler)
           ENDDO
         GOTO 299
      ELSE
         CALL OlfGetXYD (j,1,n,X,Y,D,Fehler)
         DO i = 1,n
           X1(i) = rFPA-dquot0(HBPAR2,(dquot0((X(i)-tFPA),FPPAR1)*
     *                     HBPAR1)**2)
           ENDDO
         DO K = 1, nK
           CALL OlfGetXYD (j,K,n,X,Y,D,Fehler)
           dX = X(2)-X(1)
           DO i = 2, n
              dX1 = X1(i)-X1(i-1)
              Y(i-1) = Y(i)*dX/dX1
              D(i-1) = D(i)*dX/dX1
              ENDDO
           CALL OlfPutXYD (j,K,n-1,X1,Y,D,Fehler)
           ENDDO
         GOTO 299
         ENDIF

C Open Genie like transformation of energy axis
 102  IF (iOmod.eq.1) THEN
         CALL OlfGetXYD (j,1,n,X,Y,D,Fehler)
          DO i = 1,n
             X1(i) = dquot0(X(i),(FPPAR1 + FPPAR2))*HBPAR1
             ENDDO
          DO K = 1, nK
          CALL OlfGetXYD (j,K,n,X,Y,D,Fehler)
          dX = X(2)-X(1)
          DO i = 2, n
              dX1 = X1(i)-X1(i-1)
              Y(i-1) = Y(i)*dX/dX1
              D(i-1) = D(i)*dX/dX1
              ENDDO
           CALL OlfPutXYD (j,K,n-1,X1,Y,D,Fehler)
           ENDDO
         GOTO 399
      ELSEIF (iOmod.eq.2) THEN
          CALL OlfGetXYD (j,1,n,X,Y,D,Fehler)
          DO i = 1,n
             X1(i) = rFPA-dquot0(HBPAR2,(dquot0(X(i),(FPPAR1 +
     *               FPPAR2))*HBPAR1)**2)
             ENDDO
          DO K = 1, nK
          CALL OlfGetXYD (j,K,n,X,Y,D,Fehler)
          dX = X(2)-X(1)
          DO i = 2, n
              dX1 = X1(i)-X1(i-1)
              Y(i-1) = Y(i)*dX/dX1
              D(i-1) = D(i)*dX/dX1
              ENDDO
           CALL OlfPutXYD (j,K,n-1,X1,Y,D,Fehler)
           ENDDO
         GOTO 399
         ENDIF

 103  Print *, 'monitor handling'
      IF (iOmod.eq.0) THEN
         CALL RRaw_IRSunw(j,Fehler)
         GOTO 499
      ELSEIF (iOmod.eq.1.or.iOmod.eq.2) THEN
         CALL RRaw_IRSunw(j,Fehler)
         CALL OlfGetXYD(j,1,n,X,Y,D,Fehler)
         Write (35,*) 'first number: ',X(1)
         DO i = 1,n
            X1(i) = dquot0(X(i),(FPPAR1 + FPPAR4))*HBPAR1
            ENDDO
            dX = X(2) - X(1)
         DO i = 2,n
            dX1 = X1(i) - X1(i-1)
            Y(i-1) = (Y(i)*dX)/(dX1*DEPAR1*
     *               (1.-dexp(-8.3*DEPAR2*X1(i))))
            D(i-1) = (D(i)*dX)/(dX1*DEPAR1*
     *               (1.-dexp(-8.3*DEPAR2*X1(i))))
            ENDDO
         CALL OlfPutXYD (j,1,n-1,X1,Y,D,Fehler)
         GOTO 499
         ENDIF


 199  CONTINUE
      Write (35,*) 'Output is in time of flight'
      Print *, 'Save spectra in ToF'
      RETURN

 299  CONTINUE
      Write (35,*) 'X-Axis conversion successfully accomplished'
      RETURN

 399  CONTINUE
      Write (35,*) 'X-Axis conversion similar to Modes in OpenGenie'
      RETURN

 499  CONTINUE
      Write (35,*) 'Monitor handling'
      RETURN

      END

C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSmon(j,j2,nK,Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*80  doc
      CHARACTER*(*) Fehler

! work here need file i40.f !!!!
      K1sel = 1
      iCorr = 2
      iModI = 1
      nK1=nK
      iModX = 11
      j1in = j
      j1 = j
      CALL OlfOpen (j,1,nK,Fehler)
      IF (Fehler.ne.'&ff') RETURN
      !CALL OlfOpen (j1,1,iKi,Fehler)

      CALL OlfHeadDup (j2, .false., j1out, nK2, Kout1, Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfCnuCheck2 (j2, j1, 'x', 0, 1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
      CALL OlfGetXYD (j1, K1sel, nC1, X1, Y1, D1, Fehler)



      DO K = 1, nK2

         ! do interpolation
         CALL OlfGetXYD (j2,K,nC,X,Y,D,Fehler)
         IF (Fehler.ne.'&ff') RETURN
         nC2 = nC1
         DO i = 1,nC2
            X2(i) = X1(i)
            ENDDO

         iSort0 = irSorted (X, nC)
         iSort1 = irSorted (X2, nC2)
          IF (iSort0.ne.2) THEN
               Fehler = ' x-scale of file f0 is not sorted'
               RETURN
               ENDIF
            IF (iSort1.ne.2) THEN
               Fehler = ' x-scale of file f1 is not sorted'
               RETURN
               ENDIF

         ia = irPos (X2, nC2, X(1),  'l')
         ie = irPos (X2, nC2, X(nC), 'r')
         IF (ia.gt.ie) THEN
            Fehler = ' there is no overlap of f0 and f1'
            RETURN
            ENDIF

         ! simple linear interpolation :
         DO i = ia, ie
            i0l = irPosOpt (X, nC, X2(i), 'r', i0l)
            i0h = irPosOpt (X, nC, X2(i), 'l', i0h)
            IF (i0l.lt.1 .or. i0h.gt.nC) THEN
            Print *, ' i = '//cl4(i)
            Fehler = ' Quatsch in interpolation routine'
            RETURN
            ENDIF
            CALL LinIntPol (X2(i), X(i0l), X(i0h),
     *      Y(i0l), Y(i0h), D(i0l), D(i0h), Y1(i), D1(i))
            ENDDO
C  - eliminate the empty channels :
            nC2 = nC2 - ia + 1
            ie  = ie  - ia + 1
            DO i = 1, nC2
               X2(i) = X2(ia+i-1)
               Y1(i) = Y1(ia+i-1)
               D1(i) = D1(ia+i-1)
               ENDDO
            ia = 1
C  - eliminate the empty channels :
            nC2 = ie
            CALL OlfCopZ (j2, j1out, K, K, Fehler)
            CALL OlfPutXYD (j1out, K, nC2, X2, Y1, D1, Fehler)
            IF (Fehler.ne.'&ff') RETURN
         ENDDO

         CALL OlfClos (j,nK,Fehler)

         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         DO K = 1,nK
            CALL OlfGetXYD (j,K,nC,X,Y,D,Fehler)
            DO i = 1,nC2
               Y(i) = Y(i) / Y1(i)
               D(i) = D(i) / Y1(i)
               ENDDO
            CALL OlfCopZ (j,jout,K,K,Fehler)
            CALL OlfPutXYD (jout,K,nC2,X,Y,D,Fehler)
            ENDDO

         CALL OlfClos (jout,nK,Fehler)
         CALL OlfClos (j1out,1,Fehler)


      END ! RRaw_IRSmon

C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSunw(j,Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*(*) Fehler

C  --- unwrap first determining minimum and shifting channels ---
      WRITE (35,*) 'Bin hier an 1 unwrap'
      CALL OlfGetXYD(j,1,n,X,Y,D,Fehler)
      WRITE (35,*) 'Bin hier an 2 unwrap'
      WRITE (35,*) 'Unwrap parameter ',j,n
        rMiny = Y(2)
        rMinx = X(2)
        iMin = 1
        DO i = 3,n
           IF (Y(i).lt.rMiny) THEN
              rMiny = Y(i)
              rMinx = X(i)
              iMin  = i
              ENDIF
           ENDDO
        WRITE (35,*) 'Minimum unwrap iMin = ',iMin
        WRITE (35,*) 'Minimum X(',iMin-iMin+1,',) = ',X(1)
        WRITE (35,*) 'Minima X(',iMin,',) = ',rMinx
        WRITE (35,*) 'Maximum X(',n,') =',X(n)
        WRITE (35,*) 'Minimum Y(',iMin,',) = ',rMiny
        DO i = 1,n-iMin-1
           X1(i+1) = X(1)-(X(n)-X(i+iMin))
           Y1(i+1) = Y(i+iMin)
           D1(i+1) = D(i+iMin)
           ENDDO
        WRITE (35,*) 'Min X1(',iMin-iMin+1,',) = ',X1(1)
        DO i = 2,iMin
           X1(n-iMin+i) = X(i)
           Y1(n-iMin+i) = Y(i)
           D1(n-iMin+i) = D(i)
           ENDDO
        rMinlow = 0.
        dMinlow = 0.
        rMinhigh = 0.
        dMinhigh = 0.
C      -- remove last point by linear interpolation --
        DO i = 1,5
           rMinlow = Y1(n-iMin-5+i) + rMinlow
           dMinlow = D1(n-iMin-5+i) + dMinlow
           rMinhigh = Y1(n-iMin+7-i) + rMinhigh
           dMinhigh = D1(n-iMin+7-i) + dMinhigh
           ENDDO
           X1(n-iMin+1) = X(1)
           Y1(n-iMin+1) = (rMinlow + rMinhigh)/10.
           D1(n-iMin+1) = (dMinlow + dMinhigh)/10.
C      -- first non sential point is later removed --
           X1(1) = X(1)-(X(n)-X(iMin))
           Y1(1) = Y(1)
           D1(1) = D(1)
           WRITE (35,*) 'Max X1(',n-1,',) = ',X1(n-1)
        CALL OlfPutXYD(j,1,n,X1,Y1,D1,Fehler)
C  ---------------- end unwrapping monitor -----------------------

      END ! end unwrap

C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSuns(j,nNOS,iAna,iRef,iOmod,Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*(*) Fehler

      PARAMETER    (HBPAR1=3.9560345693E-3,HBPAR2=81.804251686,
     *              FPPAR1=36.41,FPPAR2=1.45,FPPAR3=1.47,
     *              FPPAR4=-0.37,FPPAR5=2.12,FPPAR6=0.8576,
     *              DEPAR1=1.276E-3,DEPAR2=0.025)

C     this function unwraps the sample data
C     presently working for PG002 only


C determine the final energies (meV)(s. dialogue.gcl in genie modes)
      IF (iAn.eq.0) THEN
         IF (iRef.eq.0) THEN
            rFPA = 1.845
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR2,HBPAR1)
            iCut = 1880
            Print *, '1, rFPA, tFPA, iCut: ', rFPA, tFPA, iCut
         ELSEIF (iRef.eq.1) THEN
            rFPA = 7.3812
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR2,HBPAR1)
            Print *, '2, rFPA, tFPA: ', rFPA, tFPA
            ENDIF
      ELSEIF (iAn.eq.1) THEN
         IF (iRef.eq.0) THEN
            rFPA = 0.2067
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
            Print *, '3, rFPA, tFPA: ', rFPA, tFPA
         ELSEIF (iRef.eq.1) THEN
            rFPA = 0.8255
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
            Print *, '4, rFPA, tFPA: ', rFPA, tFPA
         ELSEIF (iRef.eq.2) THEN
            rFPA = 1.8567
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
            Print *, '5, rFPA, tFPA: ', rFPA, tFPA
            ENDIF
      ELSEIF (iAn.eq.2.and.iRef.eq.3) THEN
            rFPA = 1.845
            tFPA = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR2,HBPAR1)
            rFPA2 = 1.8567
            tFPA2 = dsqrt(dquot0(HBPAR2,rFPA))*dquot0(FPPAR3,HBPAR1)
C      ELSEIF (iAn.eq.3) THEN
C            GOTO 103
C            Write (35,*) 'iAnalyzer = ', iAn
         ENDIF
C     --------------depending on input calculate back or directly grid

      IF (iOmod.eq.0) THEN
         CALL OlfGetXYD (j,1,n,X,Y,D,Fehler)
            DO i = iCut+1,n
               X1(i-iCut) = X(1)-(X(n)-X(i-1))
               ENDDO
            DO i = 1,iCut
               X1(n-iCut+i) = X(i)
               ENDDO
         DO i = 1,nNOS
            CALL OlfGetXYD (j,i,n,X,Y,D,Fehler)
            DO k = iCut+1,n
               IF (i.eq.1) THEN
               Print *,'k = ',k
               ENDIF
               Y1(k-iCut) = Y(k)
               D1(k-iCut) = D(k)
               ENDDO
            DO k = 1,iCut
               Y1(n-iCut+k) = Y(k)
               D1(n-iCut+k) = D(k)
               ENDDO
         CALL OlfPutXYD (j,i,n,X1,Y1,D1,Fehler)
            ENDDO
         Print *,'n = ',n
      ELSEIF (iOmod.eq.1) THEN
         Print *, 'not yet available'
      ELSEIF (iOmod.eq.2) THEN
         Print *, 'not yet available 2'
         ENDIF


C     ---------------   ----------------------------

      END ! End RRaw_IRSuns

C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSVan(j,nK,VTemp,Xep,Yep,iAna,iRef,Fehler)
C  --------------------------------------------------------------------

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     Xep(*),Yep(*)
      CHARACTER*(*) Fehler


      IF (iAna.eq.0) THEN
         IF (iRef.eq.0) THEN
            rbl = 5.9d4
            rbh = 6.15d4
            rl  = 6.25d4
            rh  = 6.5d4
            rlamel = 6.6587
            GOTO 344
         ELSEIF (iRef.eq.1) THEN
            rbl = 2.5d4
            rbh = 2.7d4
            rl  = 3.15d4
            rh  = 3.25d4
            rlamel = 3.329
            GOTO 344
         ELSE
            Fehler = 'This is not allowed'
            GOTO 355
            ENDIF

      ELSEIF (iAna.eq.1) THEN
         IF (iRef.eq.0) THEN
            rbl = 1.86d5
            rbh = 1.88d5
            rl  = 1.89d5
            rh  = 1.92d5
            rlamel = 19.894
            GOTO 344
         ELSEIF (iRef.eq.1) THEN
            rbl = 1.0d5
            rbh = 1.015d5
            rl  = 9.45d4
            rh  = 9.65d4
            rlamel = 9.9548
            GOTO 344
         ELSEIF (iRef.eq.2) THEN
            rbl = 5.9d4
            rbh = 6.15d4
            rl  = 6.25d4
            rh  = 6.5d4
            rlamel = 6.6377
            GOTO 344
         ELSE
            Fehler = 'This is not allowed'
            GOTO 355
            ENDIF
      ELSE
         Fehler = 'This is not allowd'
         Print *, 'Probably diffraction setup check parameters'
         GOTO 355
         ENDIF

 344     CALL OlfGetXYD (j,1,n,X,Y,D,Fehler)
         DO i = 1,n!watch out:always same nr of det for bg and YepV
            IF (X(i).lt.rbl) THEN
               ibl = i + 1
            ELSEIF (X(i).lt.rbh) THEN
               ibh = i + 1
            ELSEIF (X(i).lt.rl) THEN
               il = i + 1
            ELSEIF (X(i).lt.rh) THEN
               ih = i + 1
               ENDIF
            ENDDO
         IF ((ibh-ibl).ne.(ih-il)) THEN ! to be sure !!
            ibl = ibh - (ih-il)
            ENDIF
         DO i = 1, nK
            CALL OlfGetXYD (j,i,n,X,Y,D,Fehler)
            rbg = 0.
            Yep(i) = 0.
            DO ii = ibl,ibh
               rbg = rbg+Y(ii)
               Yep(i) = Yep(i) + Y(il+ii-ibl)
               ENDDO
            Yep(i) = Yep(i) - rbg
C     Write (35,*) 'X/Yep(',i,') = ', Xep(i), Yep(i)
         ENDDO
         u2xV =  u2Debye (VTemp, 359.d0, 50.94d0, 2.5d-5)
         Write (*, '(a,f6.3,a)') '<u_x^2> =', u2xV*1000, '*10-3 A^2'
         pi = 3.14159d0
         DWFmin = 1.d0
         rsum = 0.
         DO i = 1, nK
            elaQ     = 4*pi/rlamel * dsind(Xep(i)/2)
            DWF      = dexp (-u2xV * elaQ**2)
            Yep(i)  = Yep(i)  / DWF
            rsum = rsum + Yep(i)
            DWFmin   = dmin1(DWF,DWFmin)
            ENDDO
         rsum = rsum / nK
         Write (35,*) 'After summation detector efficiency'
         DO i = 1,nK
            Yep(i) = Yep(i) / rsum
            Write (35,*) 'X/Yep(',i,') = ', Xep(i), Yep(i)
            ENDDO
         Write (*,'(a,f6.4,a,f6.2,a)')'DWF at highest angle:',DWFmin,
     *           ' for Vanadium at ' ,VTemp, 'K'

            RETURN

 355        Print *, 'Break in Vannorm' , Fehler
            Write (35,*) 'Break in Vannorm, ', Fehler
            RETURN

      END !End RRaw_IRSVan

C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSNor(j,nNumorSD,nK,Xep,Yep,Fehler)
C  --------------------------------------------------------------------
      ! This subroutine is used for normalization to Vanadium
      ! determined detector efficiencies


      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      DIMENSION     Xep(*),Yep(*)
      REAL*8        Angle(MK)
      REAL*8        Z(MK)
      CHARACTER*(*) Fehler


      Print *, 'In Nor, j,Xep(1),Yep(1)',j,Xep(1),Yep(1)
      DO i = j-nNumorSD+1,j !outer loop
         CALL OlfOpen (i,1,iK,Fehler)
         IF (Fehler.ne.'&ff') RETURN
         DO ii = 1,nK
            CALL OlfGetZ(i,ii,nZ,Z(ii),Fehler)
            IF (Fehler.ne.'&ff') RETURN
            !Check sample and normalisation file angles
            IF (dabs(Z(ii)-Xep(ii)).ge.1.d0) THEN
            Fehler = 'angles too different'
              RETURN
              ENDIF
            ENDDO
         WRITE (35,*) 'nK in loop ', nK
         DO K = 1,nK
            CALL OlfGetXYD (i,K,n,X,Y,D,Fehler)
C            WRITE (35,*) 'Yep(K) in loop',Yep(K)
            IF (Fehler.ne.'&ff') RETURN
            DO ii = 1,n
               Y(ii) = Y(ii) / Yep(K)
               D(ii) = D(ii) / Yep(K)
               ENDDO
            CALL OlfPutXYD (i,K,n,X,Y,D,Fehler)
            IF (Fehler.ne.'&ff') RETURN
            ENDDO

         WRITE (35,*) 'Normalization of file ',i,'of ',
     *                nNumorSD,'accomplished'
         CALL OlfClos (i,nK,Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO ! end outer loop

      END !End RRaw_IRSNor


C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSmcw(j,nK,iAna,iRef,Fehler)
C  --------------------------------------------------------------------
      ! This subroutine is used for transform of lambda to del E

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*(*) Fehler
      PARAMETER    (HBPAR2=81.804251686)


      IF (iAna.eq.0) THEN
         IF (iRef.eq.0) THEN
            rFPA = 1.845
         ELSE
            rFPA = 7.3812
            ENDIF
      ELSEIF (iAna.eq.1) THEN
         IF (iRef.eq.0) THEN
            rFPA = 0.2067
         ELSEIF (iRef.eq.1) THEN
            rFPA = 0.8255
         ELSE
            rFPA = 1.8567
            ENDIF
      ELSE
         Print *, 'This option is not allowed'
         RETURN
         ENDIF


      CALL OlfOpen(j,1,iK,Fehler)
      IF (Fehler.ne.'&ff') RETURN
      DO K = 1,nK
         CALL OlfGetXYD (j,K,n,X,Y,D,Fehler)
         IF (Fehler.ne.'&ff') RETURN
C         dX1=X(1)-(X(2)-X(1))
C         dX1 =  -1.*(dquot0(HBPAR2,dX1**2)-rFPA)
         X(1) = -1.*(dquot0(HBPAR2,X(1)**2)-rFPA)
C         dX1 = dabs(X(1)-dX1)
         Y(1) = Y(1)*dsqrt(1-dquot0(X(1),rFPA))
         D(1) = D(1)*dsqrt(1-dquot0(X(1),rFPA))
         DO i = 2,n
            X(i) = -1.*(dquot0(HBPAR2,X(i)**2)-rFPA)
            Y(i) = Y(i)*dsqrt(1-dquot0(X(i),rFPA))
            D(i) = D(i)*dsqrt(1-dquot0(X(i),rFPA))
            ENDDO !i
         CALL OlfPutXYD (j,K,n,X,Y,D,Fehler)
         ENDDO !K

      CALL OlfCnuP (j, 'x', 'w', 'meV', Fehler)
      CALL OlfCnuP (j, 'y', 'S(2th,w)', 'meV-1', Fehler)
      IF (Fehler.ne.'&ff') RETURN
      CALL OlfClos (j,nK,Fehler)

      END !RRaw_IRSmcw


C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSAdd(j,j1,nK,Fehler)
C  --------------------------------------------------------------------
      ! This subroutine is used for adding files

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*(*) Fehler

      ! opening first file
      CALL OlfOpen(j,1,iK,Fehler)
      IF (Fehler.ne.'&ff') RETURN
      ! opening 2nd file
      CALL OlfOpen(j1,1,iK1,Fehler)
      IF (Fehler.ne.'&ff') RETURN

      DO K=1,nK
         CALL OlfGetXYD (j,K,n,X,Y,D,Fehler)
         IF (Fehler.ne.'&ff') RETURN
         CALL OlfGetXYD (j1,K,n1,X1,Y1,D1,Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ! Test for same number of channels
         IF (n.ne.n1) THEN
            Print *, 'not the same #ch for file1 and file2'
            Write (35,*) 'adding of spectra wrong channel number'
            RETURN
            ENDIF

         DO i = 1,n
            Y(i) = Y(i)+Y1(i)
            D(i) = dsqrt0(D(i)**2+D1(i)**2)
            ENDDO
            CALL OlfPutXYD (j,K,n,X,Y,D,Fehler)
         ENDDO !K
         CALL OlfClos(j,nK,Fehler)
         CALL OlfClos(j1,nK,Fehler)

      END !RRaw_IRSAdd


C  --------------------------------------------------------------------
      SUBROUTINE RRaw_IRSAdN(j,nK,nNumor,Fehler)
C  --------------------------------------------------------------------
      ! This subroutine normalizes the added files

      IMPLICIT REAL*8  (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)
      INCLUDE      'l_def.f'
      INCLUDE      'i_dim.f'
      INCLUDE      'i_wrk.f'

      CHARACTER*(*) Fehler

      DO K = 1,nK
         CALL OlfGetXYD(j,K,n,X,Y,D,Fehler)
         DO i = 1,n
            Y(i) = Y(i) / dble(nNumor)
            D(i) = D(i) / dble(nNumor)
            ENDDO
         CALL OlfPutXYD(j,K,n,X,Y,D,Fehler)

         ENDDO

      END ! RRaw_IRSAdN
