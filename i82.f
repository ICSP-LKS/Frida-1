C Subroutines for BASIS spectrometer data read in.
C
C ----------------------------------------------------
C              History
C ----------------------------------------------------
C
C   Feb. 2010 F. Yang 1st version
C
C ----------------------------------------------------

C -----------------------------------------------------
      SUBROUTINE RRaw_BSS(Inst, Fehler)
C -----------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      IMPLICIT LOGICAL (q)

      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      INCLUDE 'i_wrk.f'

      PARAMETER   (MNuSu=100)

      DIMENSION   AngleS(*), FPath(*), DetE1(*), DetE2(*)

      CHARACTER*1 BadS(MK), BadV(MK), BadF(MK)


      CHARACTER*(*) Inst, Fehler

      CHARACTER*40 FilOut, Title
      CHARACTER*40  Fileida

      CHARACTER*80 ListSD,ListSC,ListVD,ListVC,LKnoAnal,LisKdel

      DIMENSION    JNumorSD(MNuSu), JNumorSC(MNuSu),
     *             JNumorVD(MNuSu), JNumorVC(MNuSu)






C
C      CHARACTER     cl5*5
C      CHARACTER*1   dettyp
C      INTEGER       detno,detgrp
C      REAL          detdist
C
CC  Definitions for reading NeXus file:
C      INCLUDE      'napif.inc'
C
C      CHARACTER     cl4*4,cv5*5,cv6*6
C      CHARACTER     titel*80, StartTime*20
C      INTEGER       FILEID(NXHANDLESIZE), STATUS,iBank
C      INTEGER       monitor, counts(4000*383)
C      REAL          temperature, wavelength, Drpm, Period1
C      REAL          detang(384), t_o_f(4001)
C
C      DATA  t_o_f/4001*0.0/, detang /384*0.0/, counts/1532000*0/, iBank /0/
C     1iyear /2005/
C
C      Fehler = '&ff' !
C
CC  Open the input files :
C
C      IF (irun.lt.1 .or. irun.gt.99999) THEN
C         CALL Absturz ('RRT_In_Foc', 'Bad number of run : '//cl5(irun))
C         ENDIF
C
C      iyear=iAskDMu('Year of the measurement',iyear,1998,2093)
C
C      IF(iyear.lt.2004) THEN
C         CALL Compose3 (Fileida,'/home/FOCUS/data/'//cl4(iyear),'/focus'
C     1//cv5(irun)//cl4(iyear),'.hdf')
C      ELSE
C         CALL Compose3 (Fileida,'focus'//cl4(iyear),'n'//cv6(irun),'.hdf')
C      ENDIF
C      print *, fileida
C
C      IF(NXOPEN(Fileida,NXACC_READ,FILEID) .NE. NX_OK) THEN
C         RETURN
C      ENDIF
C      IF(NXOPENGROUP(FILEID,'entry1','NXentry') .NE. NX_OK) THEN
C         CALL Absturz ('RRT_In_Foc', 'Could not open entry'//Fileida)
C      ENDIF
C
C      print *, 'Which detector bank do you want to use?'
C      print *, '   (0) Merged'
C      print *, '   (1) Upper'
C      print *, '   (2) Middle'
C      print *, '   (3) Lower'
C      iBank = iAskDMu ('Choose option', iBank, 0, 3)
C      IF(iBank.eq.0) THEN
C       IF(NXOPENGROUP(FILEID,'merged','NXdata').NE.NX_OK) THEN
C        IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
C         CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//Fileida)
C        ELSE
C         print *, 'old file, open middle bank'
C        ENDIF
C       ENDIF
C      ELSEIF(iBank.eq.1) THEN
C       IF(NXOPENGROUP(FILEID,'upperbank','NXdata').NE.NX_OK) THEN
C        IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
C         CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//Fileida)
C        ELSE
C         print *, 'old file, open middle bank'
C        ENDIF
C       ENDIF
C      ELSEIF(iBank.eq.2) THEN
C       IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
C        CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//Fileida)
C       ENDIF
C      ELSEIF(iBank.eq.3) THEN
C       IF(NXOPENGROUP(FILEID,'lowerbank','NXdata').NE.NX_OK) THEN
C        IF(NXOPENGROUP(FILEID,'bank1','NXdata').NE.NX_OK) THEN
C         CALL Absturz ('RRT_In_Foc', 'Could not open detector bank'//Fileida)
C        ELSE
C         print *, 'old file, open middle bank'
C        ENDIF
C       ENDIF
C      ENDIF
C         STATUS=NXOPENDATA(FILEID,'time_binning')
C            STATUS=NXGETDATA(FILEID,t_o_f)
C            nC=1
C            DO WHILE(t_o_f(nC).NE.0)
C               nC=nC+1
C            ENDDO
C            nC=nC-1
C            Cwidth=REAL(t_o_f(2)-t_o_f(1))
C            print *, 'time bins:',t_o_f(1),t_o_f(2),t_o_f(nC),nC
C         STATUS=NXCLOSEDATA(FILEID)
C         STATUS=NXOPENDATA(FILEID,'theta')
C            STATUS=NXGETDATA(FILEID,detang)
C            nK=1
C            DO WHILE(detang(nK).NE.0.0)
C               nK=nK+1
C            ENDDO
C            nK=nK-1
C            DO k=1,nK
C               Angle(k)=detang(k)
C            ENDDO
C            print *, 'angles:',Angle(1),Angle(nK),nK
C            STATUS=NXCLOSEDATA(FILEID)
C         IF (iyear.lt.2004) THEN
C            STATUS=NXOPENDATA(FILEID,'counts')
C            STATUS=NXGETDATA(FILEID,counts)
C            STATUS=NXCLOSEDATA(FILEID)
C            STATUS=NXOPENDATA(FILEID,'monitor')
C            STATUS=NXGETDATA(FILEID,monitor)
C            STATUS=NXCLOSEGROUP(FILEID)
C         ELSE
C            STATUS=NXOPENDATA(FILEID,'counts')
C            STATUS=NXGETDATA(FILEID,counts)
C            STATUS=NXCLOSEDATA(FILEID)
C            STATUS=NXCLOSEGROUP(FILEID)
C            IF(NXOPENGROUP(FILEID,'FOCUS','NXinstrument').NE.NX_OK) THEN
C               CALL Absturz ('RRT_In_Foc', 'Failed to open instrument'//Fileida)
C            ELSE
C               IF(NXOPENGROUP(FILEID,'counter','NXmonitor').NE.NX_OK) THEN
C                  CALL Absturz ('RRT_In_Foc', 'Could not open counter'//Fileida)
C               ELSE
C                  STATUS=NXOPENDATA(FILEID,'monitor')
C                  STATUS=NXGETDATA(FILEID,monitor)
C                  STATUS=NXCLOSEDATA(FILEID)
C                  STATUS=NXCLOSEGROUP(FILEID)
C                  STATUS=NXCLOSEGROUP(FILEID)
C               ENDIF
C            ENDIF
C         ENDIF
C            print *,'Monitor:',monitor
C            rMon=REAL(monitor)
C
CC Fileida/entry1/: sample/name,temperature
C      IF(NXOPENGROUP(FILEID,'sample','NXsample') .EQ. NX_OK) THEN
C         IF(NXOPENDATA(FILEID,'name') .EQ. NX_OK) THEN
C            IF(NXGETDATA(FILEID,titel) .NE. NX_OK) THEN
C               print *, 'failed to read sample name'
C            ELSE
C               LongTit=titel
C               print *, 'title:',LongTit
C            ENDIF
C         STATUS=NXCLOSEDATA(FILEID)
C         ENDIF
C         IF(NXOPENDATA(FILEID,'temperature') .EQ. NX_OK) THEN
C            IF(NXGETDATA(FILEID,temperature) .NE. NX_OK) THEN
C               print *, 'failed to read temperature'
C            ELSE
C               IF(temperature.LE.0.0) THEN
C                  Temp=295
C                  print *, 'temperature is not measured, it is set to 295K.'
C               ELSE
C                  Temp=temperature
C                  print *, 'temperature: ',Temp,'K'
C               ENDIF
C            ENDIF
C         STATUS=NXCLOSEDATA(FILEID)
C         ENDIF
C      STATUS=NXCLOSEGROUP(FILEID)
C      ENDIF
CC Fileida/entry1/: start_time
C      IF(NXOPENDATA(FILEID,'start_time').EQ. NX_OK) THEN
C         IF(NXGETDATA(FILEID,StartTime).NE. NX_OK) THEN
C            print *, 'failed to read start time'
C         ELSE
C            print *, 'start time: ',StartTime
C         ENDIF
C      STATUS=NXCLOSEDATA(FILEID)
C      ENDIF
C
C      IF(NXOPENGROUP(FILEID,'FOCUS','NXinstrument').NE.NX_OK) THEN
C         CALL Absturz ('RRT_In_Foc', 'Failed to open instrument'//Fileida)
C      ENDIF
C
CC Fileida/entry1/FOCUS: monochromator/lambda
C      IF(NXOPENGROUP(FILEID,'monochromator','NXmonochromator').EQ.NX_OK) THEN
C         IF(NXOPENDATA(FILEID,'lambda') .EQ. NX_OK) THEN
C            IF(NXGETDATA(FILEID,wavelength).NE. NX_OK) THEN
C               CALL Absturz ('RRT_In_Foc', 'Failed to read wavelength'
C     1//Fileida)
C            ELSE
C               WaveL=wavelength
C               print *, 'wavelength: ',WaveL,'Angstrom'
C            ENDIF
C         STATUS=NXCLOSEDATA(FILEID)
C         ENDIF
C      STATUS=NXCLOSEGROUP(FILEID)
C      ELSE
C         CALL Absturz ('RRT_In_Foc', 'Failed to read wavelength'//Fileida)
C      ENDIF
C      WaveL=rAskDMu('change wavelength?',WaveL,0.d6,50.d6)
C      print *, 'wavelength: ',WaveL,'Angstrom'
C
CC Fileida/entry1/FOCUS: disk_chopper/rotation_speed -> Period
C      IF(NXOPENGROUP(FILEID,'disk_chopper','NXchopper').EQ.NX_OK) THEN
C         IF(NXOPENDATA(FILEID,'rotation_speed').EQ.NX_OK) THEN
C            IF(NXGETDATA(FILEID,Drpm).NE.NX_OK) THEN
C               CALL Absturz ('RRT_In_Foc', 'Failed to read rotation
C     1speed'//Fileida)
C            ELSE
CC two holes on the chopper, Period in microsec:
C               Period = 30000000./Drpm ! Disk chopper speed value may be buggy
C               Period1 = t_o_f(nC) - t_o_f(1)
C               IF (dabs(Period-Period1).gt.Cwidth) THEN
C                  Period = Period1
C               ENDIF
C               print *, 'Disk chopper speed: ',Drpm,'rpm'
C               print *, 'Duty cycle: ',Period,'microsec'
C            ENDIF
C         STATUS=NXCLOSEDATA(FILEID)
C         ELSE
C            CALL Absturz ('RRT_In_Foc', 'Failed to read rotation speed'//Fileida)
C         ENDIF
C      STATUS=NXCLOSEGROUP(FILEID)
C      ELSE
C         CALL Absturz ('RRT_In_Foc', 'Failed to read rotation speed'//Fileida)
C      ENDIF
C
C      STATUS=NXCLOSEGROUP(FILEID)
C      STATUS=NXCLOSEGROUP(FILEID)
C      STATUS=NXCLOSE(FILEID)
C
C
C      Nmax = nC
C      DO K = 1,nK
CC        Detector effeciency
C         DetE1(K) = -0.017            !  rectangular 30*10mm^2, 6bar
C         DetE2(K) = -5.818            !  deff=10mm
C
C         Fpath(K) = 2.50              !  in m
C         DO i = 1,nC
C           Y(i)=REAL(counts(i+(nC*(K-1))))
C           D(i)=sqrt(Y(i))
C         ENDDO
C         CALL OlfPutSpe (j, K, 1, 0.d0, Nmax, X, Y, D, Fehler)
C         IF (Fehler.ne.'&ff') RETURN
C      ENDDO
C
C      END ! RRT_In_Foc