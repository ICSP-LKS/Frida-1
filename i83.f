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


