C  ====================================================================
C
C     Program  IDA   :  Inelastic Data Analysis
C     Modul    I_WRK :     Working Space for all IDA Subroutines
C
C  ====================================================================
         ! shared working space since 16may95

      ! i_dim must be included before including i_wrk

      ! note for multidim. arrays: leftmost variable varies most rapidly

      REAL*8            X, X1, X2, X3, X4, Y, Y1, Y2, Y3, Y4,
     *                  D, D1, D2, D3, D4, Wrk3dim, ZZofK
      INTEGER           NofK, NofK1, NofK2, NofJ

      REAL*8            rOlfG, rOlfGdef, rOlfGG, rzxOlfGG
      INTEGER           iOlfG, iOlfGdef
      LOGICAL           qOlfG, qOlfGdef

      COMMON  / iWrk /  X(MC), X1(MC), X2(MC), X3(MC), X4(MC),
     *                  Y(MC), Y1(MC), Y2(MC), Y3(MC), Y4(MC),
     *                  D(MC), D1(MC), D2(MC), D3(MC), D4(MC),
     *                  Wrk3dim(MC,MK,MWrk3dim),
     *                  NofK(MK), NofK1(MK), NofK2(MK), ZZofK(MK,MZ),
     *                  NofJ(MF)