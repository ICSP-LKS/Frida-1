C  ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i67   :     very special fit functions
C
C  ====================================================================

C     Contents :
C        GoetzeCoeff/G/X, HavNeg, Voigt, Kohlrausch, Percus-Yevick
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C  ====================================================================
C  Goetze Coefficients
C  ====================================================================

      SUBROUTINE GoetzeCoeff (yl, a, b, A1, A2, A3, B0, B1, tStar,
     *                        A0ss, A1ss, A2ss, A3ss, B0ss, B1ss,
     *                        wCro, wMin, xMin)
C     --------------------------------------------------------------
            ! JWu 23jan92
         ! Exponents and expansion coefficients from Goetze, J.Phys., 90a
         ! wCro selbst ermittelt, siehe C2,59

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h, o-p, r-z)

      REAL*4            Table1( 9,12), Table2( 9,12),
     *                  Table3(10,12), Table4(10,12)
      DIMENSION         Table(18,24)

      DATA  q1stCall / .true. /
C  Table1/2 : lambda, a, b, eA1, A2, eA3, B, B1, t*
      DATA  Table1 /
     *   .50, .395, 1.,    .616, -.053, .0040, .228, 0.,   1.571,
     *   .52, .390, .961,  .640, -.057, .0042, .255, .021, 1.508,
     *   .54, .384, .922,  .666, -.061, .0045, .284, .044, 1.443,
     *   .56, .377, .885,  .695, -.065, .0048, .317, .069, 1.376,
     *   .58, .371, .848,  .726, -.070, .0051, .353, .096, 1.308,
     *   .60, .364, .812,  .761, -.076, .0054, .393, .126, 1.237,
     *   .62, .358, .777,  .799, -.083, .0058, .438, .158, 1.164,
     *   .64, .350, .742,  .841, -.090, .0061, .488, .194, 1.088,
     *   .66, .343, .708,  .888, -.098, .0064, .545, .234, 1.011,
     *   .68, .335, .674,  .940, -.108, .0067, .609, .278,  .932,
     *   .70, .327, .641, 1.000, -.120, .0069, .681, .327,  .851,
     *   .72, .318, .608, 1.067, -.133, .0069, .765, .384,  .767 /
      DATA  Table2 /
     *   .74, .309, .575, 1.145, -.150, .0068, .861, .448,  .683,
     *   .76, .300, .542, 1.235, -.170, .0061, .974, .523,  .597,
     *   .78, .290, .509, 1.342, -.194, .0048,1.107, .610,  .512,
     *   .80, .279, .476, 1.469, -.225, .0022,1.266, .715,  .427,
     *   .82, .267, .443, 1.624, -.264,-.0027,1.459, .843,  .344,
     *   .84, .255, .409, 1.816, -.316,-.0119,1.700,1.002,  .264,
     *   .86, .241, .375, 2.062, -.387,-.0292,2.006,1.207,  .1898,
     *   .88, .226, .339, 2.388, -.489,-.0635,2.410,1.480,  .1243,
     *   .90, .209, .303, 2.842, -.643,-.1367,2.967,1.864,  .0704,
     *   .91, .200, .284, 3.143, -.753,-.2040,3.333,2.120,  .0488,
     *   .92, .190, .264, 3.518, -.898,-.3112,3.784,2.442,  .0313,
     *   .93, .180, .244, 3.997,-1.096,-.4907,4.355,2.857,  .0180 /
C  Table3/4 : lambda, A0'', eA1'', A2'', eA3'', B0'', B1'', wCro, wmi, Xmi
      DATA  Table3 /
     * .50,.860, -.318,.056,-.0003, .228,2.193,.0042,   .598,1.194,
     * .52,.842, -.326,.059,-.0007, .250,2.032,.0095,   .637,1.198,
     * .54,.823, -.335,.063,-.0011, .273,1.888,.0155,   .679,1.201,
     * .56,.804, -.345,.068,-.0015, .298,1.747,.0220,   .725,1.204,
     * .58,.785, -.356,.073,-.0020, .324,1.618,.0311,   .775,1.207,
     * .60,.765, -.367,.078,-.0026, .351,1.499,.0422,   .831,1.209,
     * .62,.745, -.379,.085,-.0031, .381,1.386,.0582,   .893,1.211,
     * .64,.725, -.392,.092,-.0037, .411,1.281,.0767,   .966,1.213,
     * .66,.704, -.406,.100,-.0043, .445,1.182,.1016,  1.050,1.215,
     * .68,.682, -.422,.109,-.0049, .480,1.091,.1318,  1.151,1.217,
     * .70,.660, -.439,.119,-.0055, .517,1.007,.1660,  1.271,1.218,
     * .72,.637, -.458,.131,-.0059, .558, .926,.1972,  1.420,1.219 /
      DATA  Table4 /
     * .74,.613, -.479,.145,-.0061, .602, .851,.2331,  1.607,1.220,
     * .76,.589, -.503,.161,-.0058, .651, .780,.2692,  1.850,1.221,
     * .78,.563, -.530,.181,-.0047, .704, .713,.3144,  2.119,1.222,
     * .80,.537, -.561,.205,-.0022, .762, .650,.3673,  2.624,1.223,
     * .82,.509, -.598,.234, .0028, .828, .590,.4193,  3.273,1.224,
     * .84,.480, -.641,.272, .0124, .903, .533,.4926,  4.282,1.224,
     * .86,.449, -.693,.321, .0306, .990, .479,.5432,  5.982,1.225,
     * .88,.416, -.757,.387, .0661,1.093, .426,.7,     9.178,1.225,
     * .90,.379, -.841,.481, .1392,1.218, .374,1.0593,16.315,1.225,
     * .91,.360, -.892,.545, .2040,1.293, .349,1.8,   23.517,1.225,
     * .92,.340, -.954,.625, .3042,1.378, .324,3.1081,36.761,1.225,
     * .93,.318,-1.028,.729, .4655,1.478, .299,4.,    64.049,1.225 /

C  Because only 19 continuation lines are allowed :
C  And conversion *4 -> *8 !
      IF (q1stCall) THEN
         DO j = 1, 12
            DO i = 1, 9
               Table(i,j)    = Table1(i,j)
               Table(i,j+12) = Table2(i,j)
               ENDDO
            DO i = 2, 10
               Table(8+i,j)    = Table3(i,j)
               Table(8+i,j+12) = Table4(i,j)
               ENDDO
            ENDDO
         q1stCall = .false.
         ENDIF

C  Find interpolation position (channel ja + fraction dl) :
      IF    (yl.le..50) THEN
         ja = 1
         dl = 0.
      ELSEIF (yl.le..90) THEN
         ja = (yl-.50)/.02 +  1
         dl = (yl-Table(1,ja))/.02
      ELSEIF (yl.lt..93) THEN
         ja = (yl-.90)/.01 + 21
         dl = (yl-Table(1,ja))/.01
      ELSE
         ja = 23
         dl = 1.
         ENDIF

      IF (dl.lt.-.001 .or. dl.gt.1.001 .or.
     *    ja.lt.1     .or. ja.gt.24) THEN
         Print *, ' yl dl ja T(1,ja) = ', yl, dl, ja, Table(1,ja)
         CALL Absturz ('GoetzeCoeff', 'Index o.o.r.')
         ENDIF

C  Interpolate :
      a    = (1-dl)*Table( 2,ja) + dl*Table( 2,ja+1)
      b    = (1-dl)*Table( 3,ja) + dl*Table( 3,ja+1)
      A1   = (1-dl)*Table( 4,ja) + dl*Table( 4,ja+1)
      A2   = (1-dl)*Table( 5,ja) + dl*Table( 5,ja+1)
      A3   = (1-dl)*Table( 6,ja) + dl*Table( 6,ja+1)
      B0   = (1-dl)*Table( 7,ja) + dl*Table( 7,ja+1)
      B1   = (1-dl)*Table( 8,ja) + dl*Table( 8,ja+1)
      tStar= (1-dl)*Table( 9,ja) + dl*Table( 9,ja+1)
      A0ss = (1-dl)*Table(10,ja) + dl*Table(10,ja+1)
      A1ss = (1-dl)*Table(11,ja) + dl*Table(11,ja+1)
      A2ss = (1-dl)*Table(12,ja) + dl*Table(12,ja+1)
      A3ss = (1-dl)*Table(13,ja) + dl*Table(13,ja+1)
      B0ss = (1-dl)*Table(14,ja) + dl*Table(14,ja+1)
      B1ss = (1-dl)*Table(15,ja) + dl*Table(15,ja+1)
      wCro = (1-dl)*Table(16,ja) + dl*Table(16,ja+1)
      wMin = dexp((1-dl)*dlog(Table(17,ja)) + dl*dlog(Table(17,ja+1)))
      xMin = (1-dl)*Table(18,ja) + dl*Table(18,ja+1)

      END ! GoetzeCoeff

      SUBROUTINE GoetzeExp (yl, a, b)
C     -------------------------------
         ! abbreviated call: lambda -> a, b

      IMPLICIT NONE
      REAL*8       yl, a, b, A1, A2, A3, B0, B1, tStar,
     *             A0ss, A1ss, A2ss, A3ss, B0ss, B1ss, wCro, wMin, xMin

      CALL GoetzeCoeff (yl, a, b, A1, A2, A3, B0, B1, tStar,
     *           A0ss, A1ss, A2ss, A3ss, B0ss, B1ss, wCro, wMin, xMin)

      END ! GoetzeCoeff

      REAL*8 FUNCTION GoetzeG (tt, yl, qFull)
C     ---------------------------------------
            ! JWu 23jan92
         ! The scaling function G(tt); yl=lambda;
         ! if not qFull, only the first order terms.

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h, o-p, r-z)

      CALL GoetzeCoeff (yl, a, b, A1, A2, A3, B0, B1, tStar,
     *           A0ss, A1ss, A2ss, A3ss, B0ss, B1ss, wCro, wMin, xMin)

      IF (qFull) THEN
         IF (tt.le.tStar) THEN ! formula 12b
            GoetzeG =        dpow0(tt, -a) - A1 * dpow0(tt,  a)
     *                + A2 * dpow0(tt,3*a) - A3 * dpow0(tt,5*a)
         ELSE ! formula (6b)
            GoetzeG = - B0 * dpow0(tt, b) + dquot0(B1,B0)*dpow0(tt,-b)
            ENDIF
      ELSE
         IF (tt.le.tStar) THEN
            GoetzeG =        dpow0(tt, -a)
         ELSE
            GoetzeG = - B0 * dpow0(tt, b)
            ENDIF
         ENDIF

      END ! GoetzeG

      REAL*8 FUNCTION GoetzeX (wIn, yl, qFull, qResc)
C     -----------------------------------------------
            ! JWu 4mar93
         ! The scaling function wg(w); yl=lambda;
         ! if not qFull, only the first order terms.
         ! if qResc, rescale minimum to (1,1)

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h, o-p, r-z)

      CALL GoetzeCoeff (yl, a, b, A1, A2, A3, B0, B1, tStar,
     *           A0ss, A1ss, A2ss, A3ss, B0ss, B1ss, wCro, wMin, xMin)

      IF (qResc) THEN
         ww = wIn * wMin
      ELSE
         ww = wIn
         ENDIF

      IF (qFull) THEN
         IF (ww.gt.wCro) THEN ! formula 14b
            xOut = A0ss * dpow0(ww,   a) - A1ss * dpow0(ww,  -a)
     *              + A2ss * dpow0(ww,-3*a) - A3ss * dpow0(ww,-5*a)
         ELSE ! formula (18)
            xOut = B0ss * dpow0(ww,-b) + B1ss * dpow0(ww,b)
            ENDIF
      ELSE ! Goetze Sjoegren '89
         xOut = A0ss * dpow0(ww,a) + B0ss * dpow0(ww,-b)
         ENDIF

      IF (qResc) THEN
         GoetzeX = xOut / xMin
      ELSE
         GoetzeX = xOut
         ENDIF

      END ! GoetzeX

C  ====================================================================
C  Empirical susceptibility fits
C  ====================================================================

      SUBROUTINE HavNeg (i, w, tau, A0, Ai, alpha, gamma, HNR, HNI)
C     -------------------------------------------------------------
            ! as subroutine JWu 11aug93
         ! calculate Re and Im of (1-(iwt)^a)^-gamma.
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      COMMON / MathConst / twopi

      IF (i.eq.1) THEN ! first call
         phea = twopi/4 * (1-alpha)
         cp   = dcos(phea)
         sp   = dsin(phea)
         ENDIF
      wta = dpow0(w*tau,alpha)
      pg  = gamma * datan ( dquot0 (wta*cp, wta*sp+1) ) !sp corrected SWi 3/02
      ar  = dpow0 ( 1 + 2*wta*sp + wta**2, gamma/2 )
      HNR = Ai - (Ai-A0) * dquot0(dcos(pg),ar) ! VZW 3/02
      HNI =      (Ai-A0) * dquot0(dsin(pg),ar)

      END ! HavNeg

      SUBROUTINE ColDav (w, tau, beta, deltaqudr, CDR, CDI)
C     -------------------------------------------------------------
         ! calculate Re and Im of Cole-Davidson
      IMPLICIT REAL*8 (a-h,o-p,r-z)
      COMMON / MathConst / twopi
      delta=deltaqudr*(twopi**2)
      wta = w * tau
      pg  = dpow0 ((1+(wta)**2),-beta / 2)
      ar  = beta*datan(-wta)
      co  = dquot0(delta,(w))
      CDI = -co * pg * dsin(ar)
      CDR = co * (pg * dcos(ar) - 1)
      END ! ColDav

      REAL*8 FUNCTION AlterC (T, alpha, aJ, qdi)
C     --------------------------------------------------------------
         ! provide second term of alternatin chain model
      IMPLICIT  LOGICAL (q)
      IMPLICIT  REAL*8 (a-h,o-p,r-z)
      REAL*8    RPARHEIS(6)

      IF (qdi) THEN
C     real parameter Heisenberg chain
C     1: A 2: B 3: C 4: D 5: E 6: F
         RPARHEIS(1) = 0.25
         RPARHEIS(2) =  -0.12587 + 0.22752*alpha
         RPARHEIS(3) = 0.019111 - .13307*alpha + .50967*alpha**2 -
     *           1.3167 * alpha**3 + 1.0081 * alpha**4
         RPARHEIS(4) = .10772 + 1.4192*alpha
         RPARHEIS(5) = -2.8521E-3 - 0.42346 * alpha + 2.1953*alpha**2-
     *           .82412 * alpha**3
         RPARHEIS(6) = 0.37754 - 6.7022E-2 * alpha + 6.9805*alpha**2 -
     *           21.678 * alpha**3 + 15.838*alpha**4



         xZahl = RPARHEIS(1) + RPARHEIS(2)*aJ/2/T +
     *           RPARHEIS(3)*aJ**2 / 4 / T**2
         xNenn = 1. + (RPARHEIS(4) * aJ/ 2/ T) + (RPARHEIS(5) * aJ**2 /
     *           4 / T**2) + RPARHEIS(6) * aJ**3 / 8 / T**3

         AlterC = xZahl/xNenn
C         PRINT *, 'xZahl, xNEnn, AlterC', xZahl, xNenn, AlterC

      ELSE
C         Print *, 'Bin hier',alpha, Aalt

         RPARHEIS(1) = 0.25
         RPARHEIS(2) = -.13695 + 0.26387 * alpha
         RPARHEIS(3) = 0.017025 - 0.12668 * alpha + 0.49113*alpha**2 -
     *          1.1977 * alpha**3 + 0.87257 * alpha**4
         RPARHEIS(4) = 0.070509 + 1.3042 * alpha
         RPARHEIS(5) = -3.5767E-3 - 0.40837 * alpha + 3.4862*alpha**2-
     *          0.73888 * alpha**3
         RPARHEIS(6) = 0.36184 - 0.065528 * alpha + 6.65875*alpha**2 -
     *          20.945 * alpha**3 + 15.425 * alpha**4



         xZahl = RPARHEIS(1) + RPARHEIS(2) * aJ / 2 / T +
     *           RPARHEIS(3) * aJ**2 / 4 / T**2
         xNenn = 1 + RPARHEIS(4)*aJ / 2 / T + RPARHEIS(5) * aJ**2 /
     *           4 / T**2 + RPARHEIS(6) * aJ**3 / 8 / T**3
         AlterC = xZahl/xNenn
         ENDIF

      END ! AlterC

      REAL*8 FUNCTION AlterS (T, alpha, aJ, TN, qdi)
C     --------------------------------------------------------------
         ! provide second term of alternating chain model
         ! for special case J-J_s
      IMPLICIT  LOGICAL (q)
      IMPLICIT  REAL*8 (a-h,o-p,r-z)
      REAL*8    RPARCHAIN(6)


         fT = 5.8
         fT = roundN(fT,7)
         AT = 0.32
         AT = roundN(AT,7)

         aJs   = dquot0(TN,4 * AT * sqrt(dlog(fT * aJ/
     *          (TN))))
         aJ    = aJ - aJs

      IF (qdi) THEN
C kopiert aus ach no 54 funktioniert
         RPARCHAIN(1) = 0.25
         RPARCHAIN(2) =  -0.12587 + 0.22752*alpha
         RPARCHAIN(3) = 0.019111 - .13307*alpha + .50967*alpha**2 -
     *           1.3167 * alpha**3 + 1.0081 * alpha**4
         RPARCHAIN(4) = .10772 + 1.4192*alpha
         RPARCHAIN(5) = -2.8521E-3 - 0.42346 * alpha + 2.1953*alpha**2-
     *           .82412 * alpha**3
         RPARCHAIN(6) = 0.37754 - 6.7022E-2 * alpha + 6.9805*alpha**2 -
     *           21.678 * alpha**3 + 15.838*alpha**4

C alt
C         RPARCHAIN(1) = 0.25
C         RPARCHAIN(2) = -.12587 + 0.22752*alpha
C         RPARCHAIN(3) = 0.019111 - 0.13307*alpha + 0.50967*alpha**2 -
C     *        1.3167*alpha**3 + 1.0081*alpha**4
C         RPARCHAIN(4) = 0.10772 + 1.4192*alpha
C         RPARCHAIN(5) = -2.8521E-3 - 0.42346*alpha + 2.1953*alpha**2-
C     *        0.82412*alpha**3
C         RPARCHAIN(6) = 0.37753 - 0.067022*alpha + 6.9805*alpha**2 -
C     *          21.678*alpha**3 + 15.838*alpha**4



         xZahlS = RPARCHAIN(1) + RPARCHAIN(2) * aJ/ 2 / T +
     *           RPARCHAIN(3) * aJ**2 / 4 / T**2
         xNennS = 1 + RPARCHAIN(4) * aJ / 2 / T + RPARCHAIN(5)*aJ**2 /
     *           4 / T**2 + RPARCHAIN(6) * aJ**3 / 8 / T**3

         AlterS = xZahlS/xNennS

      ELSE
         RPARCHAIN(1) = 0.25
         RPARCHAIN(2) = -.13695 + 0.26387 * alpha
         RPARCHAIN(3) = 0.017025 - 0.12668 * alpha + 0.49113*alpha**2 -
     *          1.1977 * alpha**3 + 0.87257 * alpha**4
         RPARCHAIN(4) = 0.070509 + 1.3042 * alpha
         RPARCHAIN(5) = -3.5767E-3 - 0.40837 * alpha + 3.4862*alpha**2-
     *          0.73888 * alpha**3
         RPARCHAIN(6) = 0.36184 - 0.065528 * alpha + 6.65875*alpha**2 -
     *          20.945 * alpha**3 + 15.425 * alpha**4

C         RPARCHAIN(1) = 0.25
C         RPARCHAIN(2) = -.13695 + 0.26387*alpha
C         RPARCHAIN(3) = 0.017025 - 0.12668*alpha + 0.49113*alpha**2 -
C     *        1.1977 * alpha**3 + 0.87257*alpha**4
C         RPARCHAIN(4) = 0.070509 + 1.3042*alpha
C         RPARCHAIN(5) = -3.5767E-3 - 0.40837*alpha + 3.4862*alpha**2 -
C     *        0.73888*alpha**3
C         RPARCHAIN(6) = 0.36184 - 0.065528*alpha + 6.65875*alpha**2 -
C     *        20.945*alpha**3 + 15.425*alpha**4


         xZahlS = RPARCHAIN(1) + RPARCHAIN(2) * aJ / 2 / T +
     *           RPARCHAIN(3) * aJ**2 / 4 / T**2
         xNennS = 1 + RPARCHAIN(4) * aJ / 2 / T + RPARCHAIN(5)*aJ**2 /
     *           4 / T**2 + RPARCHAIN(6) * aJ**3 / 8 / T**3

         AlterS = xZahlS/xNennS

         ENDIF

      END ! AlterS

      REAL*8 FUNCTION AlterS1 (T, alpha, aJ, aJs, qdi)
C     --------------------------------------------------------------
         ! provide second term of alternating chain model
         ! for special case J-J_s
      IMPLICIT  LOGICAL (q)
      IMPLICIT  REAL*8 (a-h,o-p,r-z)
      REAL*8    RPARCHAIN(6)



      IF (qdi) THEN
         RPARCHAIN(1) = 0.25
         RPARCHAIN(2) = -.12587 + 0.22752*alpha
         RPARCHAIN(3) = 0.019111 - 0.13307*alpha + 0.50967*alpha**2 -
     *        1.3167*alpha**3 + 1.0081*alpha**4
         RPARCHAIN(4) = 0.10772 + 1.4192*alpha
         RPARCHAIN(5) = -2.8521E-3 - 0.42346*alpha + 2.1953*alpha**2-
     *        0.82412*alpha**3
         RPARCHAIN(6) = 0.37753 - 0.067022*alpha + 6.9805*alpha**2 -
     *          21.678*alpha**3 + 15.838*alpha**4


         xZahlS1 = RPARCHAIN(1) + RPARCHAIN(2) * (aJ - aJs)/ 2 / T +
     *           RPARCHAIN(3) * (aJ - aJs)**2 / 4 / T**2
         xNennS1 = 1 + RPARCHAIN(4) * (aJ - aJs) / 2 / T + RPARCHAIN(5)
     *           * (aJ - aJs)**2 / 4 / T**2 + RPARCHAIN(6) *
     *           (aJ - aJs)**3 / 8 / T**3

         AlterS1 = xZahlS1/xNennS1

      ELSE
         RPARCHAIN(1) = 0.25
         RPARCHAIN(2) = -.13695 + 0.26387*alpha
         RPARCHAIN(3) = 0.017025 - 0.12668*alpha + 0.49113*alpha**2 -
     *        1.1977 * alpha**3 + 0.87257*alpha**4
         RPARCHAIN(4) = 0.070509 + 1.3042*alpha
         RPARCHAIN(5) = -3.5767E-3 - 0.40837*alpha + 3.4862*alpha**2 -
     *        0.73888*alpha**3
         RPARCHAIN(6) = 0.36184 - 0.065528*alpha + 6.65875*alpha**2 -
     *        20.945*alpha**3 + 15.425*alpha**4



         xZahlS1 = RPARCHAIN(1) + RPARCHAIN(2) * (aJ - aJs) / 2 / T +
     *           RPARCHAIN(3) * (aJ - aJs)**2 / 4 / T**2
         xNennS1 = 1 + RPARCHAIN(4) * (aJ - aJs) / 2 / T + RPARCHAIN(5)
     *           * (aJ - aJs)**2 / 4 / T**2 + RPARCHAIN(6) *
     *           (aJ - aJs)**3 / 8 / T**3

         AlterS1 = xZahlS1/xNennS1

         ENDIF

      END ! AlterS1


C  ====================================================================
C  Voigt function (Gauss-Lorentz convolution ?)
C  ====================================================================

      REAL*8 FUNCTION Voigt (bp, bg, bl, ww)
C     ---------------------------------------
            ! AMeyer apr96
         ! analytical approximation from
         ! D.G.Rancourt, Nucl. Instr. Meth. B44, 199-210 (1989)

      IMPLICIT LOGICAL (q)
      IMPLICIT REAL*8  (a-h, o-p, r-z)

      DIMENSION        T(4,4)

      DATA T /
     *   -1.2150d0, -1.3509d0, -1.2150d0, -1.3509d0,
     *    1.2359d0,  0.3786d0, -1.2359d0, -0.3786d0,
     *   -0.3085d0,  0.5906d0, -0.3085d0,  0.5906d0,
     *    0.0210d0, -1.1858d0, -0.0210d0,  1.1858d0 /

      IF (bg.eq.0) THEN
         Voigt = 0
         RETURN
         ENDIF

      x = (ww - bp) / ( bg * 2 * sqrt(2.d0) )
      y =     bl    / ( bg * 2 * sqrt(2.d0) )

      sum = 0
      DO i = 1, 4
         sum = sum + dquot0((T(i,3) *
     *      (y - T(i,1)) + T(i,4) * (x - T(i,2))),
     *     ((y - T(i,1)) * (y - T(i,1)) + (x - T(i,2)) * (x - T(i,2))))
         ENDDO
      Voigt = sum * bl / bg

      END ! Voigt

C  ====================================================================
C  Kohlrausch stretched exponential in the frequency domain
C  ====================================================================

         ! Stand 25aug98: Bisher nur fuer Im Chi getestet.
         ! Im groben gute Uebereinstimmung mit selbstberechneter Tabelle.
         ! Einzelne Ausreisser im Mittelbereich: Integration noch ganz falsch ?
         ! Zu langsam: Schleife in z nach innen verlegen, Koeffizienten vorab
         ! berechnen, besseres Abbruch/Genauigkeitskriterium
         ! Absolute Hoehe weicht um Faktor beta von selbstberechnetem ab:
         ! hier richtiger ??

      REAL*8 FUNCTION RauschInteg (z, beta, flag)
C     -------------------------------------------
         ! numerical integration

      REAL*8   z, beta, twopi, si, ul, h, a, ab, xx, tol, RauschFunc
      INTEGER  flag,n , k
      COMMON / MathConst / twopi

c ! old style (Chung) : Simpson's rule (to be restored later ?)
c si=50. ! step interval
c ll=0.0     ! lower limit
c ul=twopi/4 ! upper limit
c h=(ul-ll)/si
c a=-RauschFunc(flag,ll,beta,z)
c DO k=1, si/2
c DO m=1, 2
c yy = RauschFunc(flag,ll+(2*k+m-3)*h,beta,z)
c a = a + 2*m*yy
c c            IF (k.gt.0.996*si/2 .and. dabs(yy).lt.1e-3*a/si) THEN
c c               RauschInteg = a*h/3
c c               RETURN
c c               ENDIF
c ENDDO
c ENDDO
c RauschInteg = a*h/3

      ul= twopi/4 ! upper limit
      tol = 1.d-4
      a = 0
      n = 20
 1    CONTINUE
         h = ul / (n+1)
         ab= a
         a = 0
         DO k = 1, n
            xx = (k-0.5)*h
            a  = a + RauschFunc (flag,xx,beta,z)
            ENDDO
         a = a*h
         IF ((ab.gt.0 .and. dabs(a-ab).lt.tol*a) .or. n.gt.50000) THEN
            RauschInteg = a
            RETURN
            ENDIF
         n = 2*n
         GOTO 1

      END ! RauschInteg

      REAL*8 FUNCTION RauschFunc(flag,x,beta,z)
C     -----------------------------------------
         ! Integrand for numeric Fourier transform of the Kohlrausch function
      IMPLICIT NONE
      INTEGER flag
      REAL*8  x, beta, z, etx

      IF (x.eq.0) THEN
         etx = 1
      ELSE
         etx = dexp(-dexp(beta*dlog(dtan(x))))
         ENDIF
      IF     (flag.eq.0) THEN
         RauschFunc = etx * dcos(z*dtan(x)) / dcos(x) / dcos(x)
      ELSEIF (flag.eq.1) THEN
         RauschFunc = etx * dsin(z*dtan(x)) / dcos(x) / dcos(x)
      ELSE
         CALL Absturz ('RauschFunc', 'flag oor')
         ENDIF

      END ! RauschFunc

      REAL*8 FUNCTION RauschReihe (z, beta, flagg)
C     --------------------------------------------
         ! Series evaluation of Fourier transform of Kohlrausch function
      IMPLICIT NONE
      INTEGER  flagg, n, nn,  isig
      REAL*8   beta, z, twopi, halfpi, dgamma1, dpowii,
     *         x, v, vv, result, result1, rr, tt1, tt2, tolerance

      COMMON / MathConst / twopi
      halfpi = twopi / 4
      tolerance = 1e-5

      IF (flagg.eq.7) THEN
         nn=0
      ELSE
         nn=1
         ENDIF
      result = 0
      isig = 1

      DO n = nn, 150
         result1 = result
         IF     (flagg.eq.5) THEN ! Q, convergent
            x = n*beta+1
            rr = dgamma1(x)/dgamma1(n+1.d0)/dexp(x*dlog(z))
            result = result + isig*rr*dsin(halfpi*n*beta)
         ELSEIF (flagg.eq.6) THEN ! Q, divergent
            x = 2*n-1
            rr = dgamma1(x/beta)*dexp((x-1)*dlog(z))/dgamma1(x)/beta
            result = result + isig*rr
        ELSEIF (flagg.eq.7) THEN ! V, convergent
            x = n*beta+1
            rr = dgamma1(x)/dgamma1(n+1.d0)/dexp(x*dlog(z))
            result = result + isig*rr*dcos(halfpi*n*beta)
         ELSEIF (flagg.eq.8) THEN ! V, divergent
            x = 2*n
            rr = dgamma1(x/beta)*dexp((x-1)*dlog(z))/dgamma1(x)/beta
            result = result + isig*rr
         ENDIF
         IF (result.le.0) THEN
            Print *, z, beta, flagg, result
            Print *, 'SEVERE : result <= 0'
            RauschReihe = 0
            RETURN
            ENDIF
         tt1=dabs(rr/result)
         IF     (tt1.lt.tolerance .and. n.gt.40) THEN
            RauschReihe = result
            RETURN
         ELSEIF ((flagg.eq.6 .or. flagg.eq.8) .and.
     *           n.ge.3 .and. tt1.gt.tt2) THEN
            RauschReihe = result1
            RETURN
            ENDIF
         tt2 = tt1
         isig = -isig
         ENDDO
      ! not converged
      RauschReihe = result
      END ! RauschReihe

      REAL*8 FUNCTION RauschKohl (qImag, z, beta)
C     -------------------------------------------
         ! the Fourier transform of the Kohlrausch stretched exponential/

         ! References: Dishon et al. J. Res. Natl. Bur. Stand. 90, 27 (1985).
         ! Using qvi.c by SH Chung, revised August, 1995
         ! JWu 24aug98

      IMPLICIT NONE
      REAL*8        z, zl, beta, twopi, q, v,
     *              Qinf(8), Qsup(8), Vinf(8), Vsup(8), Blim(8),
     *              RauschReihe, RauschInteg, dgamma1
      INTEGER       i
      LOGICAL       qImag
      COMMON / MathConst / twopi

      !# Aenderungen bei Qinf, Qsup. Originalversion der Daten:
      !  DATA Qinf / .001d0, .03d0, .02d0, .1d0, .15d0, .35d0, .65d0, 2.d0 /
      !  DATA Qsup / .001d0, .0002d0, .01d0, .05d0, .1d0, .25d0, .5d0, 1.d0 /

      DATA Qinf / .0010d0,.030d0,.03d0,.15d0,.25d0,.45d0,.85d0,2.d0 /
      DATA Qsup / .0001d0,.002d0,.01d0,.05d0,.1d0,.25d0,.4d0,.5d0 /
      DATA Vinf / .0021d0,.01d0,.025d0,.06d0,.15d0,.35d0,.65d0,2.d0 /
      DATA Vsup / .001d0,.001d0,.01d0,.045d0,.1d0,.25d0,.5d0,1.d0 /
      DATA Blim / .3d0,.4d0, .5d0, .6d0, .7d0, .8d0, .9d0, 1.d0 /

      IF (z.le.0 .or. beta.le.0.2) THEN
         RauschKohl = 0
         RETURN
         ENDIF

      ! determine i from which critical z values will be determined
      DO i = 1, 8
         IF (beta.le.Blim(i)) GOTO 2
         ENDDO
 2    CONTINUE

      zl = z * dgamma1(1/beta) / beta ! to use original limits {24jan00}

      IF (qImag) THEN ! Imag of X'' -> cos transform of dPhi/dt -> evaluate Q
         IF     (zl.ge.Qinf(i)) THEN ! use Eqn (5) in Dishon et al
            RauschKohl = RauschReihe (z,beta,5)
         ELSEIF (zl.lt.Qsup(i)) THEN ! use Eqn (6)
            RauschKohl = RauschReihe (z,beta,6)
         ELSE                       ! use numerical integration
            RauschKohl = RauschInteg (z,beta,0)
            ENDIF

      ELSE           ! Re -> sin transform   -> evaluate V
         IF     (zl.ge.Vinf(i)) THEN ! use Eqn (7)
            Print *, 'sin trafo (7) of ', zl
            RauschKohl = RauschReihe (z,beta,7)
         ELSEIF (zl.lt.Vsup(i)) THEN ! use Eqn (8)
            Print *, 'sin trafo (8) of ', zl
            RauschKohl = RauschReihe (z,beta,8)
         ELSE                       ! use numerical integration to
            Print *, 'sin trafo (1) of ', zl
            RauschKohl = RauschInteg (z,beta,1)
            ENDIF

         ENDIF

      END ! Kohlrausch

C  ====================================================================
C  Statische Strukturfaktoren
C  ====================================================================

      REAL*8 FUNCTION PercYevC (qs, eta)
C     ----------------------------------
            ! JWu 29jul99
         ! nC(q*sigma;eta) from Boon Yip (2.3.15)

      IMPLICIT NONE
      REAL*8          qs, eta, a0, a1, a2

      IF (eta.le.0 .or. eta.ge.1) THEN
         PercYevC = 1
      ELSEIF (qs.lt.0) THEN
         PercYevC = 1
      ELSE
         a0 = (1+2*eta)**2 / (1-eta)**4
         a1 = -6*eta*(1+eta/2)**2/(1-eta)**4
         a2 = eta*a0/2
         IF (qs.lt.0.02) THEN ! C2,129 (3aug99)
            PercYevC = - 24*eta *
     *               ( a0 * ( 1.d0/3 - qs**2/6/5 + qs**4/120/ 7 )
     *               + a1 * ( 1.d0/4 - qs**2/6/6 + qs**4/120/ 8 )
     *               + a2 * ( 1.d0/6 - qs**2/6/8 + qs**4/120/10 ) )
         ELSE
            PercYevC = - (24*eta/qs**3) *
     *     ( a0 * (dsin(qs)-qs*dcos(qs)) +
     *     (a1/qs) * ( 2*qs*dsin(qs) - (qs**2-2)*dcos(qs) - 2 ) +
     *     (a2/qs**3) *( (4*qs**3-24*qs)*dsin(qs) - (qs**4-12*qs**2+24)
     *     *dcos(qs) + 24 ) )
            ENDIF
         ENDIF

      END ! PercYevC

      REAL*8 FUNCTION PercYevS (qs, eta)
C     ----------------------------------
         ! S(q*sigma;eta) from Boon Yip (2.3.15)

      IMPLICIT NONE
      REAL*8          qs, eta, pyc, PercYevC

      IF (eta.le.0 .or. eta.ge.1) THEN
         PercYevS = 0
      ELSE
         pyc = PercYevC(qs,eta)
         IF (pyc.ge.1) THEN
            PercYevS = 0
         ELSE
            PercYevS = 1 / ( 1 - pyc )
            ENDIF
         ENDIF

      END ! PercYevS

      REAL*8 FUNCTION PercYevDiam (q0, eta)
C     -------------------------------------

      IMPLICIT NONE
      REAL*8          q0, eta

      IF (eta.le.0 .or. eta.ge.1 .or. q0.le.0) THEN
         PercYevDiam = 0
      ELSE
         PercYevDiam = (7.00637 + (eta-.5)*(3.99569 +
     *       (eta-.5)*(3.04330 +
     *       (eta-.5)*(-3.37543 + (eta-.5)*(-11.6586))))) / q0
         ENDIF

      END ! PercYevDiam

C  ====================================================================
C  fin de fichier
C  ====================================================================
