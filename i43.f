C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i43   :     special data manipulations
C
C  ====================================================================

C     Contents :
C        RetainMaster

C  ====================================================================
C  i43 / 1 :   RetainMaster and auxiliary routines
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  --------------------------------------------------------------------
      SUBROUTINE RetainMaster (nJList, JList, Fehler)
C  --------------------------------------------------------------------
            ! JWu 4dec98. Zweiter Anlauf 25feb99.
         ! retain channels that belong to a master curve

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'i_wrk.f'
      INCLUDE 'l_def.f'
      CHARACTER     Fehler*(*)
      LOGICAL       qComN, qCom2, qReverse, qKSav(MK), qKIns(MK)
      INTEGER       nJList, JList(*), lj, j, jout, nK, Kout, K, KK, nZ,
     *              nCom, n2, n, nn, iShift(MK), ii, i2, IKdel(MK)
      REAL*8        ymean, dmean, dma, dl0, dl1, YK(MK), DK(MK), Z(MZ)

      DATA          dl0 /0.d0/, dl1 /1.d0/

      qReverse = qAskDi (' Operate from left(0) or right(1) end',
     *     intq(qReverse))

      DO lj = 1, nJList
         j = JList(lj)
         CALL OlfHeadDup (j, .false., jout, nK, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL OlfComAdd (jout, '_', 'retained master curve ', Fehler)
         IF (Fehler.ne.'&ff') RETURN

         CALL CommonScale (j, .false., qComN, nCom, qCom2,
     *                     n2, X2, iShift, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         IF (.not.qCom2) THEN
            Fehler = ' there is no common x-scale'
            RETURN
            ENDIF

         DO K = 1, nK
            qKSav(K) = .false.
            IKdel(K) = 0
            ENDDO

         DO ii = 1, n2
            IF (qReverse) THEN
               i2 = n2+1-ii
            ELSE
               i2 = ii
               ENDIF

            ! set y section:
            DO K = 1, nK
               qKIns(K) = .false.
               IF (i2.ge.1+iShift(K)) THEN
                  CALL OlfGetY (j, K, n, Y, Fehler)
                  IF (Fehler.ne.'&ff') RETURN
                  IF (i2+iShift(K).le.n) THEN
                     YK(K) = Y(i2-iShift(K))
                     ! wenn log weight, dann hier !!
                     qKIns(K) = .true.
                     ENDIF
                  ENDIF
               ENDDO

 20         CONTINUE
            KK = 0 ! Anzahl
            ymean = 0
            dmean = 0
            DO K = 1, nK
               IF (qKIns(K)) THEN
                  KK = KK + 1
                  ymean = ymean + YK(K)
                  dmean = YK(K)**2
                  ENDIF
               ENDDO
            IF (KK.le.0) GOTO 90
            ymean = ymean / KK
            dmean = dsqrt0 ( dmean / KK - ymean**2 )
            dma = 0
            DO K = 1, nK
               IF (qKIns(K) .and. .not. qKSav(K)) THEN
                  DK(K) = dabs (YK(K)-ymean)
                  IF (DK(K).ge.dma) THEN
                     KK = K ! pos. max.
                     dma = DK(K)
                     ENDIF
                  ENDIF
               ENDDO
            IF (dma.gt.dl0+dl1*dmean) THEN
               qKIns(KK) = .false.
               IKdel(KK) = ii ! Pkt. mit max. Abw. wegbeissen
c               Print *, ' deleting ii KK ', ii, KK
               GOTO 20
               ENDIF
 90         CONTINUE
            DO K = 1, nK
               IF (qKIns(K)) qKSav(K) = .true.
               ENDDO

            ENDDO ! ii

         Kout = 0
         DO K = 1, nK
            CALL OlfGetSpe (j, K, nZ, Z, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            nn = n - IKdel(K)
            IF (nn.le.0) GOTO 190
            Kout = Kout + 1
            IF (qReverse) THEN
               ii = 1
            ELSE
               ii = 1 + IKdel(K)
               ENDIF
            CALL OlfPutSpe (jout, Kout, nZ, Z, nn,
     *                      X(ii), Y(ii), D(ii), Fehler)
            IF (Fehler.ne.'&ff') RETURN
 190        CONTINUE
            ENDDO

         CALL OlfClos (jout, Kout, Fehler)
         IF (Fehler.ne.'&ff') RETURN
         ENDDO ! lj

      END ! RetainMaster

C  --------------------------------------------------------------------
      SUBROUTINE xFindPos (j, nK, val, NKdn, NKup, qKval, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER       Fehler*(*)
      LOGICAL         qKval(MK)
      INTEGER         j, nK, NKdn(MK), NKup(MK), K, n
      REAL*8          val, X(MC)

      DO K = 1, nK
         CALL OlfGetX (j, K, n, X, Fehler)
         NKdn(K) = irPos (X, n, val, 'r')
         NKup(K) = irPos (X, n, val, 'l')
         qKval(K) = qiinside(NKdn(K),1,n) .and. qiinside(NKup(K),1,n)
         ENDDO

      END ! xFindPos

C  --------------------------------------------------------------------
      SUBROUTINE OlfGetDatOfK (j, nK, IKin, QKin, NKC,
     *                         XK, YK, DK, Fehler)
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER       Fehler*(*)
      LOGICAL         QKin(MK)
      INTEGER         j, nK, IKin(MK), NKC(MK), n, K
      REAL*8          XK(MK), YK(MK), DK(MK), X(MC), Y(MC), D(MC)

      DO K = 1, nK
         IF (QKin(K)) THEN
            CALL OlfGetXYD (j, K, n, X, Y, D, Fehler)
            IF (Fehler.ne.'&ff') RETURN
            IF (IKin(K).lt.1 .or. IKin(K).gt.n) THEN
               Fehler = 'GetDatOfK/ IKin o.o.r.'
               RETURN
               ENDIF
            NKC(K) = n
            XK(K) = X(IKin(K))
            YK(K) = Y(IKin(K))
            DK(K) = D(IKin(K))
            ENDIF
         ENDDO

      END ! OlfGetDatOfK
