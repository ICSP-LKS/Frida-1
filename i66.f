C ====================================================================
C
C      Library  IDA   :  Inelastic Data Analysis
C      Modul    i66   :     ever-growing collection of fit functions
C
C  ====================================================================
C  16.02.2026 Artem Panchenko: Corrected several line breaks

C     Contents :
C        CuFuText, CuFuVal

      SUBROUTINE CuFuText (ifc, Name, Formula, ParDef, ParUni, nP)
C     ------------------------------------------------------------

      CHARACTER *(*)  Name, Formula, ParDef, ParUni

      IF     (ifc.eq. 1) THEN
         Name    = 'linear'
         Formula = 'p1 + p2 * x'
         ParDef  = 'a;b;'
         ParUni  = '0,1;-1,1'
         nP      = 2
      ELSEIF (ifc.eq. 2) THEN
         Name    = 'cubic polyn.'
         Formula = 'p1 + .. + p4 * x^3'
         ParDef  = 'a0;a1;a2;a3;'
         ParUni  = '0,1;-1,1;-2,1;-3,1'
         nP      = 4
      ELSEIF (ifc.eq. 3) THEN
         Name    = 'polynomial'
         Formula = 'p1 + .. + p7 * x^7'
         ParDef  = 'a0;a1;a2;a3;a4;a5;a6;a7'
         ParUni  = '0,1;-1,1;-2,1;-3,1;-4,1;-5,1;-6,1;-7,1'
         nP      = 8
      ELSEIF (ifc.eq. 4) THEN
         Name    = 'linear'
         Formula = 'p1 * (x-p2)'
         ParDef  = 'a;xc;'
         ParUni  = '-1,1;1,0'
         nP      = 2
      ELSEIF (ifc.eq. 5) THEN
         Name    = 'Taylor'
         Formula = 'p2+p3*(x-p1)+p4*()^2+...'
         ParDef  = 'x0;a0;a1;a2;a3;a4;a5;a6'
         ParUni  = '1,0;0,1;-1,1;-2,1;-3,1;-4,1;-5,1;-6,1'
         nP      = 8
      ELSEIF (ifc.eq. 6) THEN
         Name    = 'rational 3/3'
         Formula = 'p1+p2*x+p3*x^2+p4*x^3 / p5+p6*x+p7*x^2+p8*x^3'
         ParDef  = 'a0;a1;a2;a3;b0;b1;b2;b3'
         ParUni  = '0,1;-1,1;-2,1;-3,1;0,0;-1,0;-2,0;-3,0'
         nP      = 8
      ELSEIF (ifc.eq. 7) THEN
         Name    = 'rational 3+3/3'
         Formula =
     *'p1+p2*x+p3*x^2+p4*x^3 + p5+p6*x+p7*x^2+p8*x^3 /'//
     *'p9+p10*x+p11*x^2+p12*x^3'
         ParDef  = 'c0;c1;c2;c3;a0;a1;a2;a3;b0;b1;b2;b3'
         ParUni  =
     * '0,1;-1,1;-2,1;-3,1;0,1;-1,1;-2,1;-3,1;0,0;-1,0;-2,0;-3,0'
         nP      = 12
      ELSEIF (ifc.eq. 8) THEN
         Name    = 'FANS fast bg Cu (p5 static dspace = 1.27)'
         Formula = 'p1+p2/x+p3/x^2+p4*sin-1(sqrt(20.45/x)/p5)'
         ParDef  = 'a0;a1;a2;a3;a4'
         ParUni  = ''
         nP      = 5
      ELSEIF (ifc.eq. 9) THEN
         Name    = 'sqrt(1+x)'
         Formula = 'p1*sqrt(1+x/p2)'
         ParDef  = 'a;b'
         ParUni  = ' '
         nP      = 2
      ELSEIF (ifc.eq.10) THEN !fit a la jump relaxation model
         Name    = 'jumprelmodColmenero'
         Formula = 'p1 + p1 / p2*x^2'
         ParDef  = 'taub;l02'
         ParUni  = ' '
         nP      = 2
      ELSEIF (ifc.eq.11) THEN
         Name    = 'power law'
         Formula = 'p1*x^p2'
         ParDef  = 'h;s;'
         ParUni  = '?;0,0'
         nP      = 2
      ELSEIF (ifc.eq.12) THEN
         Name    = 'fix + power law'
         Formula = 'p1+p2*(x/p3)^p4'
         ParDef  = 'f;h;tau;s;'
         ParUni  = '0,1;1,0;0,1;0,0'
         nP      = 4
      ELSEIF (ifc.eq.13) THEN
         Name    = 'critical'
         Formula = 'p1 * |(x-p2)/p2|^p3'
         ParDef  = 'a;xc;g'
         ParUni  = '0,1;1,0;0,0'
         nP      = 3
      ELSEIF (ifc.eq.14) THEN
         Name    = 'Type A polynom ts'
         Formula = 'p1 * (x^2-2xp2 + p2^2)'
         ParDef  = 'A;Tc'
         ParUni  = '0,0;0,0'
         nP      = 2
      ELSEIF (ifc.eq.15) THEN
         Name    = 'shifted power law'
         Formula = 'p1 * |x-p2|^p3'
         ParDef  = 'a;xc;g'
         ParUni  = '0,1;1,0;0,0'
         nP      = 3
      ELSEIF (ifc.eq.16) THEN
         Name    = 'beam spec ISIS (IRIS)'
         Formula = 'p1*(p2/x)^2*e(-(p2/(x*p3))^2)*(1-e(-p4*x))'
         ParDef  = 'a;l;b;af'
         ParUni  = '0,0;0,0;0,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.17) THEN
         Name    = 'sigmoidal'
         Formula = 'p1-p1/(1+p2*e(-p3*(x-p4)))'
         ParDef  = 'a;b;c;x0'
         ParUni  = '0,0;0,0;0,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.18) THEN
         Name    = 'Singwi Sjoelander jump diff'
         Formula = '1/p1*(1-DWF/(1+p2*p1*q^2))'
         ParDef  = 'tau_0;D;<u^2>;T'
         ParUni  = '0,0;0,0;0,0;0,0;'
         nP      = 4
      ELSEIF (ifc.eq.19) THEN
         Name    = 'Diffusion law'
         Formula = 'c1+c2/2 + c1-c2/2*erf(z-z0/sqrt(4Dt)) '
         ParDef  = 'C1;C2;Z0;D;t'
         ParUni  = '0,0;0,0;0,0;0,0;0,0'
         nP      = 5
      ELSEIF (ifc.eq.21) THEN
         Name    = 'exponential'
         Formula = 'p1*exp(-p2*x)'
         ParDef  = 'a;k'
         ParUni  = '0,1;-1,0'
         nP      = 2
      ELSEIF (ifc.eq.22) THEN
         Name    = 'Exponential'
         Formula = 'p1*exp(-x/p2)'
         ParDef  = 'a;b'
         ParUni  = '0,1;1,0'
         nP      = 2
      ELSEIF (ifc.eq.23) THEN
         Name    = 'Fix + Exponential'
         Formula = 'p1 + p2*exp(-p3*x)'
         ParDef  = 'c;a;k'
         ParUni  = '0,1;0,1;-1,0'
         nP      = 3
      ELSEIF (ifc.eq.24) THEN
         Name    = 'Exponentials'
         Formula = 'sum a_*e(-x/b_)'
         ParDef  = 'aTOT;b1;a2;b2;a3;b3;a4;b4;a5;b5;a6;b6;'
         ParUni  = ' '
         nP      = 12
      ELSEIF (ifc.eq.25) THEN
         Name    = 'Arrhenius'
         Formula = 'a*e(-b/x)'
         ParDef  = 'D0;b'
         ParUni  = ' '
         nP      = 2
      ELSEIF (ifc.eq.26) THEN
         Name    = 'special function exponential'
         Formula = '1-e(-x^2 D c1) / 1-e(-x^2 D c2)'
         ParDef  = 'D;c1;c2'
         ParUni  = ' '
         nP      = 3
      ELSEIF (ifc.eq.27) THEN
         Name    = 'Bose special T'
         Formula = '1/1-e(-E/kb T)'
         ParDef  = 'E'
         ParUni  = '0,1'
         nP      = 1
      ELSEIF (ifc.eq.28) THEN
         Name    = 'Bose'
         Formula = 'p1/(e^-x/p2 - 1)'
         ParDef  = 'A;kT'
         ParUni  = '0,1;1,0'
         nP      = 2
      ELSEIF (ifc.eq.29) THEN
         Name    = 'Fermi'
         Formula = 'p1/(e^-x/p2 + 1)'
         ParDef  = 'A;kT'
         ParUni  = '0,1;1,0'
         nP      = 2
      ELSEIF (ifc.eq.30) THEN
         Name    = 'linear + Gauss'
         Formula = 'p1+ p2*x + p3*exp(-((x-p5)/p4)^2)'
         ParDef  = 'bg;m;a;b;d'
         ParUni  = '0,1;-1,1;0,1;1,0;1,0'
         nP      = 5
      ELSEIF (ifc.eq.31) THEN
         Name    = 'Gauss'
         Formula = 'p1*exp(-((x-p3)/p2)^2)'
         ParDef  = 'a;b;d;'
         ParUni  = '0,1;1,0;1,0'
         nP      = 3
      ELSEIF (ifc.eq.32) THEN
         Name    = 'Gauss 2'
         Formula = 'p1*exp(-p2*x^2)'
         ParDef  = 'a;u^2;'
         ParUni  = '0,1;2,0'
         nP      = 2
      ELSEIF (ifc.eq.33) THEN
         Name    = 'Gaussians'
         Formula = 'sum a_*e-((x-d_)/b_)^2'
         ParDef  = 'aTOT;b1;d1;a2;b2;d2;a3;b3;d3;'//
     *             'a4;b4;d4;a5;b5;d5;a6;b6;d6;'
         ParUni  = ' '
         nP      = 18
      ELSEIF (ifc.eq.34) THEN
         Name    = 'Gaussians 2'
         Formula = 'sum a_*e(-u_*x^2)'
         ParDef  = 'aTOT;u1;a2;u2;a3;u3;a4;u4;a5;u5;a6;u6;'
         ParUni  = ' '
         nP      = 12
      ELSEIF (ifc.eq.35) THEN
         Name    = 'Gauss + bg'
         Formula = 'p1+p2*exp(-((x-p4)/p3)^2)'
         ParDef  = 'bg;a;b;d;'
         ParUni  = '0,1;0,1;1,0;1,0'
         nP      = 4
      ELSEIF (ifc.eq.36) THEN
         Name    = 'Polynom*Gauss'
         Formula = '(p1+p2*x+p3*x^2)*exp(-((x-p4)/p5)^2)'
         ParDef  = 'a0;a1;a2;d;b'
         ParUni  = '0,1;-1,1;-2,1;1,0;1,0'
         nP      = 5
      ELSEIF (ifc.eq.37) THEN
         Name    = '2*Gauss+bg)'
         Formula = 'p1+p2*exp(-((x-p4)/p3)^2)+p5*exp(-((x-p7)/p6)^2)'
         ParDef  = 'bg;a1;b1;d1;a2;b2;d2'
         ParUni  = '0,1;0,1;1,0;1,0;0,1;1,0;1,0'
         nP      = 7
      ELSEIF (ifc.eq.38) THEN
         Name    = 'Gauss + Voigt + bg'
         Formula = 'b0+ag*exp(-((x-d)/bg)^2)+av*Voigt(d,bg,bl,x)'
         ParDef  = 'b0;ag;bg;d;bl;av'
         ParUni  = ''
         nP      = 6
      ELSEIF (ifc.eq.39) THEN
         Name    = 'Viscosity Cohen and Grest'
         Formula = 'no*exp(2*c1/(T-T0+sqrt[(T-T0)^2+c2*T])'
         ParDef  = 'no;T0;c1;c2'
         ParUni  = ''
         nP      = 4
      ELSEIF (ifc.eq.40) THEN
         Name    =  'dimer + Phononenfit'
         Formula = '1.5*N_d*R*exp(-2J/kBT)/(1+3*exp(-2J/kBT))+ph'
         Pardef  = 'N_d;J;b1;b2'
         ParUni  = '0,1;0,1;1,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.41) THEN
         Name    = 'ln Scaling'
         Formula = 'p1 - p2 ln (|x-p3| / p3)'
         ParDef  = 'a;c;T_0;'
         ParUni  = '0,1;0,1;1,0'
         nP      = 3
      ELSEIF (ifc.eq.42) THEN
         Name    = 'reso IN6'
         Formula = 'bg + SumGauss(a,b,d1,d2)+SumLor(aL,w,d1,d2)'
         ParDef  = 'bg;aG;bG;aL;bL;d1;d2;'
         ParUni  = '0,1;0,1;1,0;0,1;1,0;1,0;1,0'
         nP      = 7
      ELSEIF (ifc.eq.43) THEN
         Name    = 'reso 2 IN6'
         Formula = 'bg + SumGauss(a,b,d1,d2)+SumLor(aL,w,d3,d4)'
         ParDef  = 'bg;aG;bG;aL;bL;d1;d2;d3;d4:'
         ParUni  = '0,1;0,1;1,0;0,1;1,0;1,0;1,0;1,0;1,0'
         nP      = 9
      ELSEIF (ifc.eq.44) THEN
         Name    = 'dispersion'
         Formula = 'sqrt(p1^2-(p2*x)^2)-p1'
         ParDef  = 'mc2;hbarc;eps'
         ParUni  = '0,1;-1,1;0,0'
         nP      = 3
      ELSEIF (ifc.eq.45) THEN
         Name    = 'oscillation'
         Formula = 'p4+p2*cos(p2*x+p3)'
         ParDef  = 'a;omega;phi;cst'
         ParUni  = '1,0;0,-1;0,0;0,1'
         nP      = 4
      ELSEIF (ifc.eq.46) THEN
         Name    = 'LA'
         Formula = 'p1*|sin(p2*x+p3)|'
         ParDef  = 'a;omega;phi'
         ParUni  = '1,0;0,-1;0,0'
         nP      = 3
      ELSEIF (ifc.eq.47) THEN
         Name    = 'S(Q)_dimer'
         Formula = 'FFcu2p^2**A*(1-sin(Qd)/(Qd))+bg+BQ^2'
         ParDef  = 'A;d;bg;B'
         ParUni  = '0,1;1,0;0,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.48) THEN
         Name    = 'S(Q)_ach'
         Formula = 'FFcu2p^2*A*(1-sin(Qd)/(Qd)+alpha*f(q,d,b))+bg+BQ^2'
         ParDef  = 'A;d;alpha;b;bg;B'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0'
         nP      = 6
       ELSEIF (ifc.eq.49) THEN
         Name    = 'FFcu2p'
         Formula = 'A*FF(s)'
         ParDef  = 'A'
         ParUni  =  '0,0'
         nP      =  1
      ELSEIF (ifc.eq.50) THEN
         Name    = 'chain+Phononenfit'
         Formula = '2/3*N_c*R/J*T+b1*T^3+b2*T^5'
         ParDef  = 'N_c;J;b1;b2'
         ParUni  = '0,0;0,0;1,0;-1,1'
         nP      = 4
      ELSEIF (ifc.eq.51) THEN
         Name    = 'Curie Weiss mit Potenzgesetz (cgs)'
         Formula = 'p1*p2^2*m_b^2*N_A/3*k_B * p3(p3+1)/(x-p4)^p6 + p5'
         ParDef  = 'N;g;S;theta;Chi0;gamma'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0'
         nP      = 6
      ELSEIF (ifc.eq.52) THEN
         Name    = 'Bonner Fischer alpha = 1 x = p3/k_B T + CW'
         Formula = 'p1*p2^2*const/T*A+B*x+Cx^2/1+Dx+Ex^2+Fx^3+p4+CW'
         ParDef  = 'N_c;g;J;Chi0;N_cw;g_cw;S;theta;gamma'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0'
         nP      = 9
      ELSEIF (ifc.eq.53) THEN
         Name    = 'Dimer (alpha = 0)'
         Formula = 'p1*p2^2*const/x * 1/3+exp(p3/k_B*x) + p4 + CW'
         ParDef  = 'N_c;g;J;Chi0;N_cw;g_cw;S;theta;gamma'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0'
         nP      = 9
      ELSEIF  (ifc.eq.54) THEN
         Name    = 'ACH (0<=alpha<=1)'
         Formula = 'p2*p3^2*m_B*N_A/k_B*x * AlterC(x,p1,p4,qvari) + p5'
         ParDef  = 'alpha;N_c;g;J;Chi0;N_cw;g_cw;S;theta;gamma'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0'
         nP      = 10
      ELSEIF  (ifc.eq.55) THEN
         Name    = 'ACH N_tot'
         Formula =
     * 'p2*p3^2*m_B*N_A/k_B*x * AlterS1(x,p1,p4,qvari) + p5'
         ParDef  = 'alpha;N_ges;g;J;Chi0;N_cw;g_cw;S;theta;gamma'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0'
         nP      = 10
      ELSEIF (ifc.eq.56) THEN
         Name    = 'Curie Weiss inv mit Potenzgesetz (cgs)'
         Formula =
     * '1/(p1*p2^2*m_b^2*N_A/3*k_B* p3(p3+1))*(x-p4)^p6 +p5'
         ParDef  = 'N;g;S;theta;Chi0;gamma'
         nP      =  6
      ELSEIF  (ifc.eq.57) THEN
         Name    = 'ACH J-J_s (Js)'
         Formula =
     * 'p2*p3^2*m_B*N_A/k_B*x*AlterS1(x,p1,p4,p5,qvari) +p5'
         ParDef  = 'alpha;N_c;g;J;J_s;Chi0;N_cw;g_cw;S;theta;gamma'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0'
         nP      = 11
      ELSEIF (ifc.eq.60) THEN
         Name    = 'bg + Lorentz + elastic'
         Formula = 'p1 + p2[p3*|p5|/pi/((x-p4)^2+p5^2)+(1-p3)d(x-p6)]'
         ParDef  = 'bg;a;aL;dL;wL;dD'
         ParUni  = '0,0;1,1;0,0;1,0;1,0;1,0'
         nP      = 6
      ELSEIF (ifc.eq.61) THEN
         Name    = 'Lorentz'
         Formula = 'p1*|p3|/pi/((x-p2)^2+p3^2)'
         ParDef  = 'a;d;w;'
         ParUni  = '1,1;1,0;1,0'
         nP      = 3
      ELSEIF (ifc.eq.62) THEN
         Name    = '2 Lorentz'
         Formula = 'p1*|p2|/pi((x-p5)^2+p2^2)+dito(p3,p4)'
         ParDef  = 'a1;w1;a2;w2;d'
         ParUni  = '1,1;1,0;1,0;1,0;1,0'
         nP      = 5
      ELSEIF (ifc.eq.63) THEN
         Name    = 'Lorentz + elastic'
         Formula = 'p1[p2*|p4|/pi/((x-p3)^2+p4^2)+(1-p2)d(x-p5)]'
         ParDef  = 'a;aL;dL;wL;dD'
         ParUni  = '1,1;0,0;1,0;1,0;1,0'
         nP      = 5
      ELSEIF (ifc.eq.64) THEN
         Name    = 'bg + Lorentz'
         Formula = 'p4+p5*x+..p10*x^6+p1/((x-p2)^2+p3^2)'
         ParDef  = 'a;xc;xw;b0;b1;b2;b3;b4;b5;b6;'
         ParUni  = ' '
         nP      = 10
      ELSEIF (ifc.eq.65) THEN
         Name    = 'bg + SUM Lorentz'
         Formula = 'p1+p2*x+p3*x^2  + SUM p4/((x-p5)^2+p6^2)'
         ParDef  =
     * 'b0;b1;b2;a1;c1;w1;a2;c2;w2;a3;c3;w3;a4;c4;w4;a5;c5;w5;'
         ParUni  = ' '
         nP      = 18
      ELSEIF (ifc.eq.66) THEN
         Name    = 'Power laws and Lorentzian'
         Formula = 'p1+p2*((x/p3)^p4+(p3/x)^p5)+p6/((x-p7)^2+p8^2)'
         ParDef  = 'bg;H;ws;a;b;L;xc;xw;'
         ParUni  = ' '
         nP      =  8
      ELSEIF (ifc.eq.67) THEN
         Name    = 'Lorentz^n'
         Formula = 'p1*(p3^2/(x-p2)^2+p3^2)^p4'
         ParDef  = 'a;c;w;n;'
         ParUni  = '0,1;1,0;1,0;0,0'
         nP      =  4
      ELSEIF (ifc.eq.68) THEN
         Name    = 'Polynome bg + Lorentz'
         Formula = 'p6+p7*x+..p10*x^4+(p1+..+p3*x^2)/((x-p4)^2+p5^2)'
         ParDef  = 'a0;a1;a2;xc;xw;b0;b1;b2;b3;b4;'
         ParUni  = ' '
         nP      = 10
      ELSEIF (ifc.eq.69) THEN
         Name    = 'bg + Gauss + Lorentz'
         Formula =
     * 'p1+p2*|p4|/pi/((x-p3)^2+p4^2)+p5*exp(-((x-p6)/p7)^2)'
         ParDef  = 'bg;aL;dL;bL;aG;dG;bG;'
         ParUni  = ' '
         nP      =  7
      ELSEIF (ifc.eq.70) THEN
         Name    = 'bg0+ bg1*x + Gauss + Lorentz'
         Formula =
     * 'p1+p2*x+p3*|p5|/pi/((x-p4)^2+p5^2)+p6*exp(-((x-p7)/p8)^2)'
         ParDef  = 'bg;bg1;aL;dL;bL;aG;dG;bG;'
         ParUni  = ' '
         nP      =  8
      ELSEIF  (ifc.eq.71) THEN
         Name    = '2 * Lorentz doublet'
         Formula = 'p1 + Lorentz(p2,p3,p4,p5)+ dito(p6,p7,p8,p9)'
         ParDef  = 'bg;a1;c1;c1*;w1;a2;c2;c2*;w2'
         ParUni  = ' '
         nP      = 9
      ELSEIF (ifc.eq.72) THEN
         Name    = 'Lorentz (a=max)'
         Formula = 'p1*p3^2/((x-p2)^2+p3^2)'
         ParDef  = 'a;d;w;'
         ParUni  = '0,1;1,0;1,0'
         nP      = 3
      ELSEIF (ifc.eq.73) THEN
         Name    = 'S(Q)_dimer'
         Formula = 'A*(1-sin(Qd)/(Qd))+bg'
         ParDef  = 'A;d;bg'
         ParUni  = '0,1;1,0;0,0'
         nP      = 3
      ELSEIF (ifc.eq.74) THEN
         Name    = 'S(Q)_ach'
         Formula = 'A*(1-sin(Qd)/(Qd)+alpha*f(q,d,b))+bg'
         ParDef  = 'A;d;alpha;b;bg'
         ParUni  = '0,0;0,0;0,0;0,0;0,0'
         nP      = 5
      ELSEIF (ifc.eq.78) THEN
         Name    = 'Pd catalyst THF'
         Formula = 'bg+Gauss+2*Lorentz'
         ParDef  = 'bg;aG;bG;aL1;bL1;aL2;bL2;d;'
         ParUni  = '0,1;0,1;1,0;1,0;1,0;1,0;1,0;1,0'
         nP      = 8
      ELSEIF (ifc.eq.79) THEN
         Name    = 'special Pd catalyst'
         Formula = 'DWF*(L_T+SUM L_i)+const'
         ParDef  = 'zi;<u2>;R;bg;A;c1;g1;g2;g3;cq;wT;w1;w2;w3;d'
         ParUni  = '0,1;0,1;0,1;0,1;1,0;0,1;0,1;0,1;0,1;0,1;1,0;
     *              1,0;1,0;1,0;0,1'
         nP      = 15
      ELSEIF (ifc.eq.81) THEN
         Name    = 'Kohlrausch + bg'
         Formula = 'p1*exp(-(x/p2)^p3)'
         ParDef  = 'a;tau;beta;bg'
         ParUni  = '0,1;1,0;0,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.82) THEN
         Name    = 'improved Kohlrausch'
         Formula = 'p1*exp(-(x/t[p2])^p3)'
         ParDef  = 'a;<tau>;beta;'
         ParUni  = '0,1;1,0;0,0'
         nP      = 3
      ELSEIF (ifc.eq.83) THEN
         Name    = '-d/dx Kohlrausch'
         Formula = '-d/dx p1*exp(-(x/p2)^p3)'
         ParDef  = 'a;tau;beta;'
         ParUni  = '1,1;1,0;0,0'
         nP      = 3
      ELSEIF (ifc.eq.84) THEN
         Name    = 'Kohlrausch * exponential'
         Formula = 'p1*exp(-(x/t[p2])^p3)*exp(-x/p4)'
         ParDef  = 'a;<tau>;beta;tau_exp'
         ParUni  = '0,1;1,0;0,0;1,0'
         nP      = 4
      ELSEIF (ifc.eq.85) THEN
         Name    = 'impr. Kohlrausch + background'
         Formula = 'p1+p2*exp(-(x/t[p3])^p4)'
         ParDef  = 'b;a;<tau>;beta'
         ParUni  = '0,1;0,1;1,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.86) THEN
         Name    = 'double Kohlrausch'
         Formula = 'p1*exp(-(x/t[p2])^p3)+p4*exp(-(x/t[p5])^p6)'
         ParDef  = 'b;<tau1>;beta1;a;<tau2>;beta2;'
         ParUni  = '0,1;1,0;0,0;0,1;1,0;0,0'
         nP      = 6
      ELSEIF (ifc.eq.87) THEN
         Name    = 'special Ft 123, 2*Kohlrausch, 1Exp'
         Formula = 'p1*exp(-x/p2)+ p3 + p4*G(t) +
     *              p7*exp(-(x/t[p8])^p9 + p10 * exp(-(x/t[p11]))^p12'
         ParDef  = 'pho;tp;fq;h*s;te;lambda;a;tau1;beta1;b;tau2;beta2'
         ParUni  = '0,1;1,0;0,1;0,1;1,0;0,0;0,1;1,0;0,0;0,1;1,0;0,0'
         nP      = 12
      ELSEIF (ifc.eq.88) THEN
         Name    = 'impr. Kohl + bg with a+b fixed'
         Formula = '(p1-p2)+p2*exp(-(x/t[p3])^p4)'
         ParDef  = 'fq;a;<tau>;beta'
         ParUni  = '0,1;0,1;1,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.89) THEN
         Name    = 'impr. Kohlrausch + pos.  background'
         Formula = '|p1|+p2*exp(-(x/t[p3])^p4)'
         ParDef  = 'b;a;<tau>;beta'
         ParUni  = '0,1;0,1;1,0;0,0'
         nP      = 4
      ELSEIF (ifc.eq.91) THEN
         Name    = 'ln Vogel-Fulcher'
         Formula = 'p1 - ( p2 / (x-p3) )'
         ParDef  = 'A;E_0;T_0;'
         ParUni  = '0,1;1,1;1,0'
         nP      = 3
      ELSEIF (ifc.eq.92) THEN
         Name    = 'ln VF+..'
         Formula = 'p1 - p2/(x-p3) - (p4/(x-p5))^2'
         ParDef  = 'A;E_0;T_0;E_1;T_1'
         ParUni  = '0,1;1,1;1,0;1,1;1,0'
         nP      = 5
      ELSEIF (ifc.eq.95) THEN
         Name    = 'order parameter T_c eval'
         Formula = 'M_s*(1-T/T_c)^b'
         ParDef  = 'M_s;T_c;b'
         ParUni  = '0,0;1,0;0,0'
         nP      = 3
      ELSEIF (ifc.eq.96) THEN
         Name    = 'Bloch T^3/2 law'
         Formula = 'M_s*(1-a*T^3/2)'
         ParDef  = 'M_s;a'
         ParUni  = '0,0;0,0'
         nP      = 2
      ELSEIF (ifc.eq.101) THEN
         Name    = 'tof-background'
         Formula = 'p1'
         ParDef  = 'bg;E0;'
         ParUni  = '0,1;?'
         nP      =  2
      ELSEIF (ifc.eq.102) THEN
         Name    = 'Delta tof-background'
         Formula = 'bg(w)-bg(-w)'
         ParDef  = 'bg;E0;'
         ParUni  = '0,1;?'
         nP      =  2
      ELSEIF (ifc.eq.103) THEN
         Name    = 'INS detector efficiency'
         Formula = 'DetEff[p1,p2,p5;p3;p4]'
         ParDef  = 'C1;C2;E0;E_b;f_b'
         ParUni  = '?;?;1,0;1,0;?'
         nP      =  5
      ELSEIF (ifc.eq.106) THEN
         Name    = 'fpi-reflexion'
         Formula = 'p1-L(p2,p4,p6,p8)-L(p3,p5,p6*p7,p8*p9)'
         ParDef  = 'bg;a1;a2;w1;w2;d1;dr2;p1;dr2;nord'
         ParUni  = '0,1;0,1;0,1;1,0;1,0;1,0;0,0;1,0;0,0;0,0'
         nP      = 10
      ELSEIF (ifc.eq.107) THEN
         Name    = 'fpi-transmission'
         Formula = 'p1*A(p3*(x-p4),p2)^3*A(p3*p5*(x-p4),p2)^3'
         ParDef  = 'A;F;k;x0;costh'
         ParUni  = '0,1;0,0;-1,0;1,0;0,0'
         nP      = 5
      ELSEIF (ifc.eq.108) THEN
         Name    = 'Lorentz^n * Saturation'
         Formula = '[p1*(p3/(x-p2)^2+p3^2)^p4]*e-p5*[..]'
         ParDef  = 'a;c;w;n;tau'
         ParUni  = '?;1,0;1,0;0,0;0,-1'
         nP      =  5
      ELSEIF (ifc.eq.109) THEN
         Name    = 'resolution FPI'
         Formula = 'p1*R(x-p2-n*p4;p3)^3*R(..-..*p5-p7,p5*p6*p3)'
         ParDef  = 'A;w0;G;P;rG;rF;dw'
         ParUni  = '?;1,0;1,0;0,0;0,0;0,0;1,0'
         nP      =  7
      ELSEIF (ifc.eq.110) THEN
         Name    = 'resolution FPI * saturation'
         Formula =
     * 'p1*R(x-p2-n*p4;p3)^3*R(..-..*p5-p7,p5*p6*p3)*e-p8[..]'
         ParDef  = 'A;w0;G;P;rG;rF;dw;tau'
         ParUni  = '?;1,0;1,0;0,0;0,0;0,0;1,0;0,-1'
         nP      =  8
      ELSEIF (ifc.eq.111) THEN
         Name    = 'echo (rectang)'
         Formula = 'p1/2 (1+p2*cos 2pi*x/p3*(sin/arg) 2pi*x/p3/p4)'
         ParDef  = 'A;P;xp;env;x0'
         ParUni  = ' '
         nP      = 5
      ELSEIF (ifc.eq.112) THEN
         Name    = 'echo (triang)'
         Formula = 'p1/2 (1+p2*cos 2pi*x/p3*(sin/arg)^2 2pi*x/p3/p4)'
         ParDef  = 'A;P;xp;env;x0'
         ParUni  = ' '
         nP      = 5
      ELSEIF (ifc.eq.113) THEN
         Name    = 'diffraction by grid'
         Formula = '[(1/p1)(sin p2 x p1 / sin p2 x) sinc p3 x]^2'
         ParDef  = 'N;line-halfsep;line-halfwidth'
         ParUni  = '0,0;-1,0;-1,0'
         nP      = 3
      ELSEIF (ifc.eq.121) THEN  !  in der Form Li et al 92a, eqn. 4.14
         Name    = 'Sjoegren in w'
         Formula = 'p1 (p4*x/p2^p3+p3*p2/x^p4)/p3+p4'
         ParDef  = 'c;ws;a;b;'
         ParUni  = ' '
         nP      =  4
      ELSEIF (ifc.eq.122) THEN
         Name    = 'Two power laws'
         Formula = 'p1+p2*p3/x^p5-p2*(x/p4)^p6)'
         ParDef  = 'f;h;A;tau;a;b'
         ParUni  = ' '
         nP      =  6
      ELSEIF (ifc.eq.123) THEN
         Name    = 'Goetze f+h*G(t) liquid'
         Formula = 'p1+p2*G(t)'
         ParDef  = 'f;h*s;te;lambda'
         ParUni  = '0,1;0,1;1,0;0,0'
         nP      =  4
      ELSEIF (ifc.eq.124) THEN
         Name    = 'Goetze f+h*G(t) .. short'
         Formula = 'p1+p2*G(t)..short'
         ParDef  = 'f;h*s;te;lambda'
         ParUni  = '0,1;0,1;1,0;0,0'
         nP      =  4
      ELSEIF (ifc.eq.125) THEN
         Name    = 'Goetze X(w)'
         Formula = 'p1*w/p2*g_p3(w/p2)'
         ParDef  = 'X0;w0;lambda'
         ParUni  = '0,1;1,0;0,0'
         nP      =  3
      ELSEIF (ifc.eq.126) THEN
         Name    = 'Goetze X(w) .. short'
         Formula = 'p1*w/p2*g''_p3(w/p2)'
         ParDef  = 'X0;w0;lambda'
         ParUni  = '0,1;1,0;0,0'
         nP      =  3
      ELSEIF (ifc.eq.127) THEN
         Name    = 'Goetze S(w)'
         Formula = 'p1/p2*g_p3(w/p2)'
         ParDef  = 'X0;w0;lambda'
         ParUni  = ' '
         nP      =  3
      ELSEIF (ifc.eq.128) THEN
         Name    = 'Goetze f+h*G(t) liquid t-fix'
         Formula = 'p1+p2*G(t)'
         ParDef  = 'f;h*s;t_ampl;lambda'
         ParUni  = ' '
         nP      =  4
      ELSEIF (ifc.eq.129) THEN
         Name    = 'Goetze X(w/wmin)/Xmin'
         Formula = 'p1*w/p2*g_p3(w/p2)'
         ParDef  = 'Xmin;wmin;lambda'
         ParUni  = '0,1;1,0;0,0'
         nP      =  3
      ELSEIF (ifc.eq.131) THEN
         Name    = 'square root law'
         Formula = 'p1+Th(si)p2*sqrt(si)+p4*si, si=p3-x/p3'
         ParDef  = 'fc;h;Tc;a'
         ParUni  = ' '
         nP      =  4
      ELSEIF (ifc.eq.132) THEN
         Name    = 'square root law (two lin slopes)'
         Formula =
     * 'p1+Th(si)p2*sqrt(si)+Th(si)p4*si+Th(-si)p5*si, si=p3-x/p3'
         ParDef  = 'fc;h;Tc;a_cold;a_hot'
         ParUni  = ' '
         nP      =  5
      ELSEIF (ifc.eq.133) THEN
         Name    = 'square root law'
         Formula = 'p1*e^(p4(p3-x)) + (if x>p3)p2*sqrt((p3-x)/p3)'
         ParDef  = 'fc;h;Tc;a'
         ParUni  = ' '
         nP      =  4
      ELSEIF (ifc.eq.141) THEN
         Name    = 'Debye <u^2>'
         Formula = '<u^2 (x,p1,p2)>'
         ParDef  = 'TD;m;eps;'
         ParUni  = '1,0;0,0;1,0'
         nP      =  3
      ELSEIF (ifc.eq.142) THEN
         Name    = 'Debye <u^2>'
         Formula = '<u^2 (x,p1,p2)> - <u^2 (p3,..)>'
         ParDef  = 'TD;m;T0;eps;'
         ParUni  = '1,0;0,0;1,0;0,0'
         nP      =  4
      ELSEIF (ifc.eq.143) THEN
         Name    = 'Debye-Waller factor from Debye <u^2>'
         Formula = 'p6*exp-p5^2*(<u^2 (x,p1,p2)> - <u^2 (p3,..)>)'
         ParDef  = 'TD;m;T0;eps;q;A'
         ParUni  = ' '
         nP      =  6
      ELSEIF (ifc.eq.144) THEN
         Name    = 'm-phonon Q-integral'
         Formula = 'sum^m (u2 x^2)^m/m! * exp(-u2 x^2)'
         ParDef  = 'u2;m;'
         ParUni  = ' '
         nP      = 2
      ELSEIF (ifc.eq.145) THEN
         Name    = 'm-phonon contribution'
         Formula = 'a (u2 x^2)^m/m! * exp(-u2 x^2)'
         ParDef  = 'a;u2;m;'
         ParUni  = ' '
         nP      = 3
      ELSEIF (ifc.eq.146) THEN
         Name    = 'Debye C(T)'
         Formula = '3R*n*c(T/TD)'
         ParDef  = 'TD;n;nint;'
         ParUni  = ' '
         nP      =  3
      ELSEIF (ifc.eq.147) THEN
         Name    = 'Recoil'
         Formula = 'p1*S_rec(p2,x;p3,p4)'
         ParDef  = 'a;q;m;T'
         ParUni  = ' '
         nP      = 4
      ELSEIF (ifc.eq.151) THEN
         Name    = '1-phonon Q-dependance'
         Formula = 'p1 x^2 * exp(-p2^2 x^2)'
         ParDef  = 'a;u;'
         ParUni  = ' '
         nP      = 2
      ELSEIF (ifc.eq.152) THEN
         Name    = 'for g(w) - i'
         Formula = '(p2*u^2+p3*u^3+p4*u^4+p5*u^5)e(-u^2/2),u=x/p1'
         ParDef  = 'w0;a2;a3;a4;a5;'
         ParUni  = ' '
         nP      = 5
      ELSEIF (ifc.eq.153) THEN
         Name    = 'for g(w) - ii'
         Formula = '(p2*u^2+p3*u^3+p4*u^4+p5*u^5)e(-u),u=x/p1'
         ParDef  = 'w0;a2;a3;a4;a5;'
         ParUni  = ' '
         nP      = 5
      ELSEIF (ifc.eq.154) THEN
         Name    = 'Debye-Gauss components'
         Formula = 'SUM a(w^2/u^3)exp(-w^2/2u^2)'
         ParDef  = 'a1;w1;a2;w2;a3;w3;a4;w4;a5;w5;'//
     *             'a6;w6;a7;w7;a8;w8;a9;w9;a10;w10;'
         ParUni  = ' '
         nP      = 20
      ELSEIF (ifc.eq.155) THEN
         Name    = 'Debye-Lorentz components'
         Formula = 'SUM a(w^2/u^3)exp(-w/u)'
         ParDef  = 'a1;w1;a2;w2;a3;w3;a4;w4;a5;w5;'//
     *             'a6;w6;a7;w7;a8;w8;a9;w9;a10;w10;'
         ParUni  = ' '
         nP      = 20
      ELSEIF (ifc.eq.156) THEN
         Name    = 'sharp Debye components'
         Formula = 'SUM a(3w^2/u^3) for w<=u'
         ParDef  = 'a1;w1;a2;w2;a3;w3;a4;w4;a5;w5;'//
     *             'a6;w6;a7;w7;a8;w8;a9;w9;a10;w10;'
         ParUni  = ' '
         nP      = 20
      ELSEIF (ifc.eq.157) THEN
         Name    = 'SUM(Lorentz) + elastic + bg'
         Formula = 'a[SUM(an*L(wn,dn))+(1-SUM(an))*d(x-dD)+bg]'
         ParDef  = 'a;a1;d1;w1;a2;d2;w2;a3;d3;w3;a4;d4;w4;b1;b2;'//
     *        'b3;b4;dD;dp;'
         ParUni  = '1,1;0,0;1,0;1,0;0,0;1,0;1,0;0,0;1,0;1,0;0,0;'//
     *        '1,0;1,0;1,0;1,0;1,0;1,0;1,0;1,0'
         nP      = 19
      ELSEIF (ifc.eq.158) THEN
         Name    = '2 * Lorentz (normalized)'
         Formula = 'a*[|a1|*L(w1,d1)+(1-|a1|)*L(w2,d2)]+bg'
         ParDef  = 'a;a1;d1;w1;d2;w2;bg;'
         ParUni  = '1,1;0,0;1,0;1,0;0,0;1,0;1,0;'
         nP      = 7
      ELSEIF (ifc.eq.159) THEN
         Name    = 'long range diff. & iso. rot. diff.'
         Formula = 'a*[A0(|aq|)*L(d1,w1)+SUM[AN(|aq|)*L(d2,Dr)]+bg'
         ParDef  = 'a;aq;d1;w1;d2;Dr;bg;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 7
      ELSEIF (ifc.eq.160) THEN
         Name    = 'long range diff. & uniaxial rot. diff.'
         Formula = 'a*[A0(|aq|)*L(d1,w1)+SUM[AN(|aq|)*F(|tr|)]+bg'
         ParDef  = 'a;aq;d1;w1;tr;bg;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 6
      ELSEIF (ifc.eq.161) THEN
         Name    = 'long range diff. & 2 uniaxial rot. diff.'
         Formula = 'a*[Tr(D)+Rot1(a1,t1)+Rot2(a2,t2)]+bg'
         ParDef  = 'a;D;a1;t1;a2;t2;bg;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 7
      ELSEIF (ifc.eq.170) THEN
         Name    = 'sph. bess. first kind J_0(x)'
         Formula = 'P0*J_0(P1*X(i))'
         ParDef  = 'A;R;'
         ParUni  = '0,0;0,0;'
         nP      = 2
      ELSEIF (ifc.eq.171) THEN
         Name    = 'Zorn et al. JCP 2002'
         Formula = 'P2*((1-P0)*3fold jump + P0)'
         ParDef  = 'cf;R;A;'
         ParUni  = '0,0;0,0;0,0;'
         nP      = 3
      ELSEIF (ifc.eq.172) THEN ! yf 12.08.2010
         Name    = 'sph. bess. first kind J_0^2(x)+bg'
         Formula = 'A*[(1-bg)*J_0^2(R*X(i))+bg]'
         ParDef  = 'A;R;bg'
         ParUni  = '0,0;0,0;0,0'
         nP      = 3
      ELSEIF (ifc.eq.180) THEN
         Name    = 'special small angle DLM'
         Formula = 'I(q)*S(q)'
         ParDef  = 'G;B;Rg;P;k;xsi;A;n;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 8
      ELSEIF (ifc.eq.181) THEN
         Name    = 'special small angle II'
         Formula = 'mod. I(Q)*S(Q)'
         ParDef  = 'G1;B1;Rg1;P1;k1;ksi1;G2;B2;Rg2;P2;k2;xsi2;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 12
      ELSEIF (ifc.eq.190) THEN
         Name    = 'brioullin rayleigh'
         Formula = 'delta + 2*Lorentz + bg'
         ParDef  = 'A;B;w;d;bg;lim;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 6
      ELSEIF (ifc.eq.191) THEN
         Name    = 'brioullin lorentz + rayleigh'
         Formula = 'Lorentz + 2*Lorentz + bg'
         ParDef  = 'A;B;w1;d1;w2;d2;bg;bg1;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 8
      ELSEIF (ifc.eq.192) THEN
         Name    = '190 amplitude lorentz decoupled'
         Formula = 'delta + 2*Lorentz + bg'
         ParDef  = 'A;B;w;d;bg;lim;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 6
      ELSEIF (ifc.eq.193) THEN
         Name    = '191 amplitude lorentz decouple'
         Formula = 'Lorentz + 2*Lorentz + bg'
         ParDef  = 'A;B;w1;d1;w2;d2;bg;bg1;'
         ParUni  = '0,0;0,0;0,0;0,0;0,0;0,0;0,0;0,0;'
         nP      = 8
      ELSEIF (ifc.eq.195) THEN ! special function NMR ratios
         Name    = 'special NMR ratioed Gaussians'
         Formula = 'sum a_*e-((x-d_)/b_)^2'
         ParDef  = 'aTOT;p1;b1;d1;p2;b2;d2;p3;b3;d3;'//
     *             'p4;b4;d4;p5;b5;d5;p6;b6;d6;'
         ParUni  = ' '
         nP      = 19
      ELSEIF (ifc.eq.201) THEN ! Lorentz-Amplitude = -Omega / (x + i*Omega)
         Name    = 'Re Lorentz'
         Formula = '-p1*p2*x/[(x-p3)^2+p2^2]'
         ParDef  = 'A;G;dw'
         ParUni  = '0,1;1,0;1,0'
         nP      = 3
      ELSEIF (ifc.eq.202) THEN
         Name    = 'Im Lorentz'
         Formula = 'p1*p2^2/[(x-p3)^2+p2^2]'
         ParDef  = 'A;G;dw'
         ParUni  = '0,1;1,0;1,0'
         nP      = 3
      ELSEIF (ifc.eq.211) THEN
         Name    = 'S(qw) Brill. Liq.'
         Formula = 'p1*[p2/((x+p3)^2+p2^2) + ..]'
         ParDef  = 'A;d;w4'
         ParUni  = '1,-1;1,0;1,0'
         nP      = 3
      ELSEIF (ifc.eq.221) THEN
         Name    = 'Re Havriliak Negami'
         Formula = 'p2+(p1-p2)*Re(1-ixp2^p3)^-p4'
         ParDef  = 'A0;Ai;tau;a;g'
         ParUni  = ' '
         nP      = 5
      ELSEIF (ifc.eq.222) THEN
         Name    = 'Im Havriliak Negami'
         Formula = 'p1*Im(1-ixp2^p3)^-p4'
         ParDef  = 'A;tau;a;g;'
         ParUni  = ' '
         nP      = 4
      ELSEIF (ifc.eq.223) THEN
         Name    = 'S(w) Havriliak Negami'
         Formula = 'p1/x*Re(1-ixp2^p3)^-p4 + p5'
         ParDef  = 'A;tau;a;g;bg'
         ParUni  = ' '
         nP      = 5
      ELSEIF (ifc.eq.224) THEN
         Name    = 'Brill-S(qw) of HN'
         Formula = 'p1/x*ImX-1(p2,p3,p4):HN(p5,p6,p7)+fri(p8,p9,p10)'
         ParDef  = 'A;q;c0;ci;tau;a;g;G0;G1;G2'
         ParUni  = ' '
         nP      = 10
      ELSEIF (ifc.eq.225) THEN
         Name    = 'Brill-X''(qw) of HN'
         Formula = 'p1*ImX-1(p2,p3,p4):HN(p5,p6,p7)+fri(p8,p9,p10)'
         ParDef  = 'A;q;c0;ci;tau;a;g;G0;G1;G2'
         ParUni  = ' '
         nP      = 10
      ELSEIF (ifc.eq.226) THEN
         Name    = 'X'' of HN and boson-peak'
         Formula = 'Int_W W^2 I0(W) (2p5/p6) X"(x)/x of W^2-w^2-d2*HN'
         ParDef  = '&j;tau;a;g;Oo;T;T0;G;Oc;G2;debug'
         ParUni  = ' '
         nP      = 10
      ELSEIF (ifc.eq.227) THEN
         Name    = 'X'' of dos and HN'
         Formula = 'Int_W W^2 I0(W) (2p5/p6) X"(x)/x of W^2-w^2-d2*HN'
         ParDef  = '&j;tau;a;g;Del;T;T0;pow;fri'
         ParUni  = ' '
         nP      = 8
      ELSEIF (ifc.eq.228) THEN
         Name    = 'KWW Gauss bg hybrid'
         Formula = 'bg + S_KWW#234 + Gauss'
         ParDef  = 'bg0;A;<tau>;b;ag;d;bg'
         ParUni  = '0,1;0,1;-1,0;0,0;0,1;0,1;0,1'
         nP      = 7
      ELSEIF (ifc.eq.229) THEN
         Name    = 'KWW S(q,w) like 234 plus sloping bg'
         Formula = '234 + p6*x'
         ParDef  = 'A;<tau>;b;bg0;ad;bg1'
         ParUni  = '0,1;-1,0;0,0;0,1;0,1;0,0'
         nP      = 6
      ELSEIF (ifc.eq.230) THEN
         Name    = 'KWW S(q,w) like 233 plus sloping bg'
         Formula = '233 + p5*x'
         ParDef  = 'A;<tau>;b;bg0;bg1'
         ParUni  = '0,1;-1,0;0,0;0,1;0,0'
         nP      = 5
      ELSEIF (ifc.eq.231) THEN
         Name    = 'Re KWW (w) [ft-tab]'
         Formula = 'p1 + p2 * ReX(x*p3;p4)'
         ParDef  = 'A_infty;Delta;<tau>;b'
         ParUni  = ' '
         nP      =  4
      ELSEIF (ifc.eq.232) THEN
         Name    = 'Im KWW (w)'
         Formula = 'p1 * ImX(x*p2;p3)'
         ParDef  = 'A;<tau>;b'
         ParUni  = ' '
         nP      =  3
      ELSEIF (ifc.eq.233) THEN
         Name    = 'S_KWW (w)'
         Formula = 'p1 * S(x*p2;p3) + p4'
         ParDef  = 'A;<tau>;b;bg'
         ParUni  = '0,1;-1,0;0,0;0,1'
         nP      =  4
      ELSEIF (ifc.eq.234) THEN
         Name    = 'S_KWW (w) + elast'
         Formula = 'p1 * S(x*p2;p3) + p4 + p5*delta'
         ParDef  = 'A;<tau>;b;bg;ad'
         ParUni  = '0,1;-1,0;0,0;0,1;0,1'
         nP      =  5
      ELSEIF (ifc.eq.235) THEN
         Name    = 'KWW: e''''(e'')'
         Formula = 'p1+p2*Im(Re-p3/p4-p3)'
         ParDef  = 'offs;A;eps0;epsI;beta'
         ParUni  = '0,1;0,1;1,0;1,0;0,0'
         nP      =  5
      ELSEIF (ifc.eq.236) THEN
         Name    = 'Re KWW (w) [series]'
         Formula = 'p1 + p2 * ReX(x*p3;p4)'
         ParDef  = 'A_infty;Delta;<tau>;b'
         ParUni  = ' '
         nP      =  4
      ELSEIF (ifc.eq.237) THEN
         Name    = 'Im KWW (w)'
         Formula = 'p1 * ImX(x*p2;p3)'
         ParDef  = 'A;<tau>;b'
         ParUni  = ' '
         nP      =  3
      ELSEIF (ifc.eq.238) THEN
         Name    = 'S_KWW (w)'
         Formula = 'p1 * S(x*p2;p3) + p4'
         ParDef  = 'A;<tau>;b;bg'
         ParUni  = '0,1;-1,0;0,0;0,1'
         nP      =  4
      ELSEIF (ifc.eq.239) THEN
         Name    = '--> debug'
         Formula = 'integrand'
         ParDef  = 'b;z'
         ParUni  = '0,0;1,0'
         nP      =  2
      ELSEIF (ifc.eq.240) THEN
         Name    = 'Delta + bg'
         Formula = 'p1*delta + p2'
         ParDef  = 'a;bg'
         ParUni  = '0,1;0,1'
         nP      =  2
      ELSEIF (ifc.eq.241) THEN
         Name    = 'X'' of dos and KWW [tab]'
         Formula = 'Int_W W^2 I0(W) (2p5/p6) X"(x)/x of W^2-w^2-d2*KWW'
         ParDef  = '&j;tau;beta;Del;T;T0;pow;fri'
         ParUni  = ' '
         nP      = 7
      ELSEIF (ifc.eq.243) THEN
         Name    = 'S_KWW (w) + positive bg'! yf 09.08.2010
         Formula = 'p1 * S(x*p2;p3) + |p4|'
         ParDef  = 'A;<tau>;b;bg'
         ParUni  = '0,1;-1,0;0,0;0,1'
         nP      =  4
      ELSEIF (ifc.eq.251) THEN
         Name    = 'Martin:Brill.m.CD'
         Formula = 'p1(p2+CDI)/(x^2-p3^2+xCDR)^2...'
         ParDef  = 'A;go;Temp;tau;beta;deltaqudr;w1;w2'
         ParUni  = ' '
         nP      = 8
      ELSEIF (ifc.eq.252) THEN
         Name    = 'Martin:Brill.m.hybrid'
         Formula = 'CD1+ iwP(1/(9)-iw)^(P(10)-1)'
         ParDef  = 'A;go;Temp;tau;beta;deltaqudr;w1;w2;B;a'
         ParUni  = ' '
         nP      = 10
      ELSEIF (ifc.eq.261) THEN
         Name    = 'Percus-Yevick direct correlation'
         Formula = 'nC(k)'
         ParDef  = 'diam;eta'
         ParUni  = '-1,0;0,0'
         nP      =  2
      ELSEIF (ifc.eq.262) THEN
         Name    = 'Percus-Yevick structure factor'
         Formula = 'S(k)'
         ParDef  = 'diam;eta'
         ParUni  = '-1,0;0,0'
         nP      =  2
      ELSEIF (ifc.eq.263) THEN
         Name    =
     * 'Percus-Yevick structure factor with variable diameter'
         Formula = 'S(k)'
         ParDef  = 'q0;eta'
         ParUni  = '1,0;0,0'
         nP      =  2
      ELSEIF (ifc.eq.264) THEN
         Name    = 'Percus-Yevick vertex'
         Formula = 'V(q;k,p;d,e)'
         ParDef  = 'diam;eta;k;p'
         ParUni  = '-1,0;0,0;1,0;1,0'
         nP      =  4
      ELSEIF (ifc.eq.265) THEN
         Name    =
     * 'PY structure factor with distributed and variable diameter'
         Formula = 'S(k)'
         ParDef  = 'q0;eta;drelrad'
         ParUni  = '1,0;0,0;0,0'
         nP      =  3
      ELSEIF (ifc.eq.272) THEN
         Name    = 'tube <y^2>'
         Formula = 'see C2, 145ff.'
         ParDef  = 'eps'
         ParUni  = '0,0'
         nP      = 1
      ELSEIF (ifc.eq.273) THEN
         Name    = '4 Voigt + bg'
         Formula = 'Rancourt/Ping'
         ParDef  = 'd0;d1;q1;b1;q2;b2;a1;a2;a3;a4;gl;bg'
         ParUni  = ' '
         nP      = 12
      ELSEIF (ifc.eq.274) THEN
         Name    = 'Voigt + bg'
         Formula = 'a*Voigt(d,bg,bl,x) + b0'
         ParDef  = 'd;bg;bl;a;b0'
         ParUni  = ''
         nP      = 5
      ELSEIF (ifc.eq.275) THEN
         Name    = 'damped oscillation'
         Formula = 'c+a*cos[2pi*f(x-x0)]*exp[-x*b]'
         ParDef  = 'c;a;f;x0;b'
         ParUni  = ''
         nP      = 5
      ELSEIF (ifc.eq.276) THEN
         Name    = 't dep. damped osci.'
         Formula = 'bg+A*cos[(f+dx)*(x-x0)]*exp[-x*(b+cx)]'
         ParDef  = 'bg;A;f;d;x0;b;c'
         ParUni  = ''
         nP      = 7
      ELSEIF (ifc.eq.277) THEN
         Name    = 'damped osci.+additional osci.'
         Formula =
     * 'c+a*cos[2pi*f(x-x0)]*exp[-x*b]+d*cos[2pi*f1*(x-x1)]'
         ParDef  = 'c;a;f;x0;b;d;f1;x1'
         ParUni  = ''
         nP      = 8
      ELSEIF (ifc.eq.278) THEN
         Name    = 'damped osci.+ additional osci. coupled'
         Formula = 'c+a*cos[f(x-x0)]*exp[-x*b]+d*cos[2pi*f*(x-x0)]'
         ParDef  = 'c;a;f;x0;b;d;'
         ParUni  = ''
         nP      = 6
      ELSEIF (ifc.eq.279) THEN
         Name    = 'Free radiative cooling'
         Formula = '-a(t-b)^(-1/3)'
         ParDef  = 'a;b;'
         ParUni  = ''
         nP      = 2
      ELSEIF (ifc.eq.280) THEN
         Name    = 'NFS: FT{e^(i*Lorentz)}'
         Formula = 'p1*decay(x/p5)'
         ParDef  = 'A;L0f;bgr;tau_0;1/t_diff;no.bunches;t_bunchsep'
         ParUni  = '0,1;0,0;0,1;1,0;-1,0;0,0;1,0'
         nP      =  7
      ELSEIF (ifc.eq.281) THEN
         Name    = 'NFS: zwei FT{e^(i*Lorentz)} (Ferrocene)'
         Formula = 'p1*decay(x/p5)*interf(p8)'
         ParDef  =
     * 'A;L0f;bgr;tau_0;1/t_diff;no.bunches;t_bunchsep;w_sep'
         ParUni  = '0,1;0,0;0,1;1,0;-1,0;0,0;1,0;?'
         nP      =  8
      ELSEIF (ifc.eq.282) THEN
         Name    = 'Fe3Si N"aherung'
         Formula = '2 exponentials'
         ParDef  = 'A;q;L;V;B'
         ParUni  = ' '
         nP      =  5
      ELSEIF (ifc.eq.283) THEN
         Name    = 'density, time oscillation'
         Formula = 'p0+p*(-a(t-t0)^(-1/3))+d*sin(w(t-t1)-f0)'
         ParDef  = 'p0;p;a;t0;d;w;t1;f0'
         ParUni  = ''
         nP      =  8
      ELSEIF (ifc.eq.284) THEN
         Name    = 'density, Temp oscillation'
         Formula = 'p0+p*T+d*sin(w(-(a/T)^3+t0-t1)-f0)'
         ParDef  = 'p0;p;d;w;a;t0;t1;f0'
         ParUni  = ''
         nP      =  8
      ELSEIF (ifc.eq.291) THEN
         Name    = 'Debye-S(t) [inelastic part]'
         Formula = 'e^-2W(0) (e^2W(p1,x;a=p2,wD=p3)-1)'
         ParDef  = 'q;r0^2;wD'
         ParUni  = '0,0;0,0;-1,0'
         nP      =  3
      ELSEIF (ifc.eq.292) THEN
         Name    = 'Gauss-Debye-S(t) [inelastic part]'
         Formula = 'e^-2W(0) (e^2W(p1,x;a=p2,wD=p3)-1)'
         ParDef  = 'q;r0^2;wD'
         ParUni  = '0,0;0,0;-1,0'
         nP      =  3
      ELSEIF (ifc.eq.293) THEN
         Name    = 'provis. C2,156'
         Formula = 'f(tau1)'
         ParDef  = '2th;n'
         ParUni  = '0,0;0,0'
         nP      =  2
      ELSEIF (ifc.ge.300 .and. ifc.le.399) THEN ! NFS special routine
         CALL NFS_FuT (ifc-300, Name, Formula, ParDef, ParUni, nP)
      ELSEIF (ifc.ge.400 .and. ifc.le.409) THEN ! MCT special routine
         CALL MCT_FuT (ifc-400, Name, Formula, ParDef, ParUni, nP)
      ELSEIF (ifc.ge.410 .and. ifc.le.419) THEN ! MCT special routine
         CALL MCT_SFT (ifc-400, Name, Formula, ParDef, ParUni, nP)

      ELSEIF (ifc.gt.400) THEN
         Name    = '&lastentry'
      ELSE
         Name    = '&undefined'
         ENDIF

      END ! CuFuText

      SUBROUTINE CuFuVal (ifc, P, X, Y, n, Fehler)
C     --------------------------------------------
         ! Calculate Y(1..n)

      IMPLICIT REAL*8  (a-h, o-p, r-z)
      INCLUDE 'i_dim.f'

      DIMENSION         P(*), X(*), Y(*), Y1(MC), Y2(MC)
      CHARACTER         Fehler*(*)

      ! Linux: gegen vermeintlichen type mismatch in dsqrt,... :
      REAL*8        E0, E

      COMMON / MathConst / twopi
      COMMON / FitParFil / XFPF(MC), YFPF(MC), nFPF
      DATA     hbar /0.658218/    ! in 10^-15 eV sec

C     T. Unruh added field for more complex model fit functions
C     157-161 incorp. standard FRIDA by FK05
      DOUBLE PRECISION  FNU, ITAU1, ITAU2, Q10RATIO
      DOUBLE PRECISION  JJJ0(100), AAA(100)
      DOUBLE PRECISION  AAA1(100), AAA2(100), JJJ1(100), JJJ2(100)
      INTEGER           IFAIL, NUM, NZ, NUM1, NUM2, i, j, k
      CHARACTER*1       SCALE
      COMPLEX*16        Z
      COMPLEX*16        CY(100)
      REAL*8            bZ, baY(100) !Artem add for dbesi

      IF      (ifc.eq. 1) THEN
         DO i = 1, n
            Y(i) = P(1) + P(2) * X(i)
            ENDDO
      ELSEIF (ifc.eq. 2) THEN
         DO i = 1, n
            Y(i) = P(1) + P(2)*X(i) + P(3)*X(i)**2 + P(4)*X(i)**3
            ENDDO
      ELSEIF (ifc.eq. 3) THEN
         DO i = 1, n
            Y(i) = P(1) + X(i)*( P(2) + X(i)*( P(3) + X(i)*( P(4) +
     *                    X(i)*( P(5) + X(i)*( P(6) + X(i)*( P(7) +
     *                    X(i)*  P(8)  ))))))
            ENDDO
      ELSEIF (ifc.eq. 4) THEN ! linear
         DO i = 1, n
            Y(i) = P(1) * (X(i)-P(2))
            ENDDO
      ELSEIF (ifc.eq. 5) THEN ! Taylor 'p2+p3*(x-p1)+p4*()^2+...'
         DO i = 1, n
            xa = X(i)-P(1)
            Y(i) = P(2) + xa*( P(3) + xa*( P(4) + xa*( P(5) +
     *                    xa*( P(6) + xa*( P(7) + xa*P(8))))))
            ENDDO
      ELSEIF (ifc.eq. 6) THEN ! p1+p2*x+p3*x^2+p4*x^3 / p5+p6*x+p7*x^2+p8*x^3
         DO i = 1, n
            Y(i) = dquot0 (
     * P(1) + P(2)*X(i) + P(3)*X(i)**2 + P(4)*X(i)**3,
     * P(5) + P(6)*X(i) + P(7)*X(i)**2 + P(8)*X(i)**3)
            ENDDO
      ELSEIF (ifc.eq. 7) THEN ! polynom + rational function
         DO i = 1, n
            Y(i) =
     *         P( 1) + P( 2)*X(i) + P( 3)*X(i)**2 + P( 4)*X(i)**3 +
     * dquot0 (P( 5) + P( 6)*X(i) + P( 7)*X(i)**2 + P( 8)*X(i)**3,
     *         P( 9) + P(10)*X(i) + P(11)*X(i)**2 + P(12)*X(i)**3)
            ENDDO
      ELSEIF (ifc.eq. 8) THEN
         DO i = 1,n
           Y(i)  =  P(1) + dquot0(P(2),X(i)) + dquot0(P(3),X(i)**2) +
     *              P(4) * dasin0(dsqrt0(20.45/X(i))/P(5))
            ENDDO
      ELSEIF (ifc.eq. 9) THEN ! constant * sqrt
         DO i = 1, n
            sigma = dquot0((P(2)+X(i)),P(2))
            IF (sigma.gt.0) THEN
               Y(i) = P(1) * dsqrt(sigma)
            ELSE
               Y(i) = 0
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.10) THEN ! Colmenero jump relax model
         DO i = 1, n
            Y(i) = P(1) + dquot0(P(1),P(2)*X(i)**2)
            ENDDO
      ELSEIF (ifc.eq.11) THEN ! better power law (numerically more stable)
         DO i = 1, n
            Y(i) = P(1) * dpow0 (X(i),P(2))
            ENDDO
      ELSEIF (ifc.eq.12) THEN ! fix + power law
         DO i = 1, n
            Y(i) = P(1) + P(2) *
     *             dpow0 (dquot0(X(i),P(3)),P(4))
            ENDDO
      ELSEIF (ifc.eq.13) THEN ! 'critical'
         DO i = 1, n
            Y(i) = P(1) * dpow0 (dabs(dquot0 (X(i)-P(2), P(2)) ),P(3))
            ENDDO
      ELSEIF (ifc.eq.14) THEN ! Type a for ts
         DO i = 1, n
            Y(i) = P(1) * (X(i)**2 - 2*P(2)*X(i) + P(2)**2)
            ENDDO
      ELSEIF (ifc.eq.15) THEN ! shifted power law
         DO i = 1, n
            Y(i) = P(1) * dpow0 ( dabs (X(i)-P(2)), P(3) )
            ENDDO
      ELSEIF (ifc.eq.16) THEN ! ISIS IRIS beam spec
         DO i = 1, n
            Y(i) = (P(1) * dquot0(P(2),X(i))**2) *
     *             dexp1(- dquot0(P(2),X(i)*P(3))**2)*
     *             (1. - dexp1(-P(4)*X(i)))
            ENDDO
      ELSEIF (ifc.eq.17) THEN ! sigmoidal for el. int.
         DO i= 1, n
             Y(i) = P(1) - dquot0(P(1),(1.+P(2)*dexp1(P(3)*(X(i)-
     *                P(4)))))
             ENDDO
      ELSEIF (ifc.eq.18) THEN ! Singwi Sjoelander
         DO i=1, n
            Y(i) = 1./P(1)*(1. - dquot0(dexp1(-P(3)*P(4)*X(i)),
     *             (1. + P(1) * P(2) * X(i))))
            ENDDO
      ELSEIF (ifc.eq.19) THEN ! Diffusion law
         DO i = 1,n
            Y(i) = (P(1)+P(2))/2. + (P(1)-P(2))/2. *
     *         ERF(dquot0(X(i)-P(3),dsqrt(4*P(4)*P(5)))) !Artem: Replace with a standard ERF function: S15AEF(dquot0(X(i)-P(3),dsqrt(4*P(4)*P(5))),IFAIL)
            ENDDO
      ELSEIF (ifc.eq.21) THEN ! exponential
         DO i = 1, n
            Y(i) = P(1) * dexp1 ( - P(2) * X(i) )
            ENDDO
      ELSEIF (ifc.eq.22) THEN ! exponential
         DO i = 1, n
            Y(i) = P(1) * dexp1 ( - dquot0(X(i),P(2)) )
            ENDDO
      ELSEIF (ifc.eq.23) THEN ! fix + exponential
         DO i = 1, n
            Y(i) = P(1) + P(2) * dexp1 ( - P(3) * X(i) )
            ENDDO
      ELSEIF (ifc.eq.24) THEN ! 6 Exponentials
         DO i = 1, n
            p1 = P(1) - P(3) - P(5) - P(7) - P(9) - P(11)         ! P(1) = aTOT
            Y(i) =    p1    * dexp1 ( -dquot0(X(i),P( 2)) )
     *             +  P( 3) * dexp1 ( -dquot0(X(i),P( 4)) )
     *             +  P( 5) * dexp1 ( -dquot0(X(i),P( 6)) )
     *             +  P( 7) * dexp1 ( -dquot0(X(i),P( 8)) )
     *             +  P( 9) * dexp1 ( -dquot0(X(i),P(10)) )
     *             +  P(11) * dexp1 ( -dquot0(X(i),P(12)) )
            ENDDO
      ELSEIF (ifc.eq.25) THEN ! Arrhenius
         DO i = 1, n
            Y(i) = P(1) * dexp1 (-P(2) / X(i))
            ENDDO
      ELSEIF (ifc.eq.26) THEN ! special function exponential
         DO i = 1, n
            Y(i) = dquot0((1 - dexp(-X(i)**2 * P(1) * P(2))),
     *               (1 - dexp(-X(i)**2 * P(1) * P(3))))
            ENDDO
      ELSEIF (ifc.eq.27) THEN ! Bose special T
         DO i=1,n
            Y(i) = 1./ (1.-dexp(-P(1)/(8.61739E-2 * X(i))))
            ENDDO
      ELSEIF (ifc.eq.28) THEN ! Bose
         DO i = 1, n
            Y(i) = P(1) / ( 1. - dexp( dquot0(X(i),P(2)) ) )
            ENDDO
      ELSEIF (ifc.eq.29) THEN ! Fermi
         DO i = 1, n
            Y(i) = P(1) / ( 1. + dexp( dquot0(X(i),P(2)) ) )
            ENDDO
      ELSEIF (ifc.eq.30) THEN ! linear + Gaussian
         DO i = 1, n
            Y(i) = P(1) + P(2)*X(i) + P(3) *
     *                   dexp1 ( - dquot0(X(i)-P(5),P(4))**2 )
            ENDDO
      ELSEIF (ifc.eq.31) THEN ! Gaussian
         DO i = 1, n
            Y(i) = P(1) * dexp1 ( - dquot0(X(i)-P(3),P(2))**2 )
            ENDDO
      ELSEIF (ifc.eq.32) THEN ! Gaussian for use as DWF
         DO i = 1, n
            Y(i) = P(1) * dexp1 ( - P(2) * X(i)**2 )
            ENDDO
      ELSEIF (ifc.eq.33) THEN ! 6 Gaussians
         DO i = 1, n
            p1 = P(1) - P(4) - P(7) - P(10) - P(13) - P(16)       ! P(1) = aTOT
            Y(i) =    p1    * dexp1 ( -dquot0(X(i)-P( 3), P( 2))**2 )
     *             +  P( 4) * dexp1 ( -dquot0(X(i)-P( 6), P( 5))**2 )
     *             +  P( 7) * dexp1 ( -dquot0(X(i)-P( 9), P( 8))**2 )
     *             +  P(10) * dexp1 ( -dquot0(X(i)-P(12), P(11))**2 )
     *             +  P(13) * dexp1 ( -dquot0(X(i)-P(15), P(14))**2 )
     *             +  P(16) * dexp1 ( -dquot0(X(i)-P(18), P(17))**2 )
            ENDDO
      ELSEIF (ifc.eq.34) THEN ! 6 Gaussians for DWF
         DO i = 1, n
            p1 = P(1) - P(3) - P(5) - P(7) - P(9) - P(11)       ! P(1) = aTOT
            Y(i) =    p1    * dexp1 ( - P( 2) * X(i)**2)
     *             +  P( 3) * dexp1 ( - P( 4) * X(i)**2)
     *             +  P( 5) * dexp1 ( - P( 6) * X(i)**2)
     *             +  P( 7) * dexp1 ( - P( 8) * X(i)**2)
     *             +  P( 9) * dexp1 ( - P(10) * X(i)**2)
     *             +  P(11) * dexp1 ( - P(12) * X(i)**2)
            ENDDO
      ELSEIF (ifc.eq.35) THEN ! Gaussian
         DO i = 1, n
            Y(i) = P(1) + P(2) * dexp1 ( - dquot0(X(i)-P(4),P(3))**2 )
            ENDDO
      ELSEIF (ifc.eq.36) THEN ! Polynom*Gaussian
         DO i = 1, n
            Y(i) = (P(1) + P(2)*X(i) + P(3)*X(i)**2) *
     *             dexp1(- dquot0(X(i),P(5))**2)
            ENDDO
      ELSEIF (ifc.eq.37) THEN ! 2*Gaussian + bg
         DO i = 1,n
            Y(i) = P(1) + P(2) * dexp1 ( - dquot0(X(i)-P(4),P(3))**2 )+
     *             P(5) * dexp1 ( - dquot0(X(i)-P(7),P(6))**2 )
            ENDDO
      ELSEIF (ifc.eq.38) THEN ! Gaussian + Voigt
         DO i = 1, n
            ww=X(i)
            Y(i) = P(1) + P(2) * dexp1 ( - dquot0(X(i)-P(4),P(3))**2 )
     *           + P(6) * Voigt(P(4),dabs(P(3)),P(5),ww)
            ENDDO
      ELSEIF (ifc.eq.39) THEN ! Viscosity Cohen and Grest
         DO i = 1, n
            Y(i) =
     * P(1)*dexp1(2.0*P(3)/(X(i)-P(2)+dsqrt((X(i)-P(2))**2+P(4)*X(i))))
            ENDDO
      ELSEIF (ifc.eq.40) THEN ! dimer+Phononenfit
         f = 8314.51
         f = roundN(f,7)
         h1 = 1.
         h2 = 3.
         h3 = 6.
         DO i = 1, n
            Y(i) = h3 * f * P(1) *  P(2)**2 / X(i)**2 *
     *             dquot0 ( dexp1( -2*P(2)/X(i) ), h1 + h2
     *             * dexp1(-2*P(2)/X(i)))
     *             + P(3)*X(i)**3 + P(4)*X(i)**5
            ENDDO
      ELSEIF (ifc.eq.25) THEN ! Polynom*Gaussian
         DO i = 1, n
            Y(i) = (P(1) + P(2)*X(i) + P(3)*X(i)**2) *
     *             dexp1 ( - dquot0(X(i)-P(4),P(5))**2 )
            ENDDO
      ELSEIF (ifc.eq.41) THEN
         DO i = 1, n
            a = dquot0 ( dabs(X(i)-P(3)), P(3) )
            Y(i) = P(1) - P(2) * dln0(a)
            ENDDO
      ELSEIF (ifc.eq.42) THEN ! IN6 reso1
         pi = twopi/2
         DO i = 1, n
            Y(i) = P(1) + P(2)*(dexp1(- dquot0(X(i)-P(6),P(3))**2 ) +
     *             dexp1(- dquot0(X(i)-P(7),P(3))**2)) +
     *             P(4)*(dquot0(dabs(P(5))/pi,(X(i)-P(6))**2+P(5)**2) +
     *             dquot0(dabs(P(5))/pi,(X(i)-P(7))**2+P(5)**2))
            ENDDO
      ELSEIF (ifc.eq.43) THEN ! IN6 reso2
         pi = twopi/2
         DO i = 1, n
            Y(i) = P(1) + P(2)*(dexp1(- dquot0(X(i)-P(6),P(3))**2 ) +
     *             dexp1(- dquot0(X(i)-P(7),P(3))**2)) +
     *             P(4)*(dquot0(dabs(P(5))/pi,(X(i)-P(8))**2+P(5)**2) +
     *             dquot0(dabs(P(5))/pi,(X(i)-P(9))**2+P(5)**2))
            ENDDO
      ELSEIF (ifc.eq.44) THEN ! 'dispersion'
         DO i = 1, n
            IF (dabs(P(2)*X(i)).lt.P(1)*P(3)) THEN ! zwengs der Numerik
               Y(i) = (P(2)*X(i))**2/P(1)/2 - (P(2)*X(i))**4/P(1)**3/8
            ELSE
               Y(i) = dsqrt(P(1)**2 + (P(2)*X(i))**2) - P(1)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.45) THEN ! 'oscillation'
         DO i = 1, n
            Y(i) = P(4) + P(1) * dcos (P(2)*X(i) + P(3))
            ENDDO
      ELSEIF (ifc.eq.46) THEN ! 'LA'
         DO i = 1, n
            Y(i) = P(1) * dabs( dsin (P(2)*X(i) + P(3)) )
            ENDDO
      ELSEIF (ifc.eq.47) THEN ! S(Q)_dimer * FFCu2p^2
         PP = 3.141592654
         PP = roundN(PP,7)
C         Print *, 'PP=',PP
         a = 0.0232
         a = roundN(a,7)
C         Print *, 'a=', a
         b = 34.969
         b = roundN(b,7)
C         Print *, 'b=', b
         c = 0.4023
         c = roundN(c,7)
C         Print *, 'c=', c
         d = 11.564
         d = roundN(d,7)
         e = 0.5882
         e = roundN(e,7)
         f = 3.843
         f = roundN(f,7)
         g = -0.0137
         g = roundN(g,7)
C         Print *, 'D=', D
         r = 1.5189
         r = roundN(r,7)
         s = 10.478
         s = roundN(s,7)
         t = 1.1512
         t = roundN(t,7)
         u = 3.813
         u = roundN(u,7)
         v = 0.2918
         v = roundN(v,7)
         w = 1.398
         w = roundN(w,7)
         z = 0.0017
         z = roundZN(z,7)
         DO i = 1, n
         Y(i) = (( a * dexp1(-b*X(i)**2/PP**2/16)
     *         + c * dexp1(-d*X(i)**2/PP**2/16)
     *         + e * dexp1(-f*X(i)**2/PP**2/16) + g)
     *         + (1-2/1.2)
     *     * ( r * X(i)**2/PP**2/16 * dexp1(-s*X(i)**2/PP**2/16)
     *      +  t * X(i)**2/PP**2/16 * dexp1(-u*X(i)**2/PP**2/16)
     *      +  v * X(i)**2/PP**2/16 * dexp1(-w*X(i)**2/PP**2/16)
     *      +  z * X(i)/PP/4))**2 *
     *       P(1) * (1 - dquot0( dsin(X(i)*P(2)), X(i)*P(2)))
     *            +  P(3) + P(4) * X(i)**2
            ENDDO
      ELSEIF (ifc.eq.48) THEN ! S(Q)_dimer_ach *FFCu2p^2
         PP = 3.141592654
         PP = roundN(PP,7)
C         Print *, 'PP=',PP
         a = 0.0232
         a = roundN(a,7)
C         Print *, 'a=', a
         b = 34.969
         b = roundN(b,7)
C         Print *, 'b=', b
         c = 0.4023
         c = roundN(c,7)
C         Print *, 'c=', c
         d = 11.564
         d = roundN(d,7)
         e = 0.5882
         e = roundN(e,7)
         f = 3.843
         f = roundN(f,7)
         g = -0.0137
         g = roundN(g,7)
C         Print *, 'D=', D
         r = 1.5189
         r = roundN(r,7)
         s = 10.478
         s = roundN(s,7)
         t = 1.1512
         t = roundN(t,7)
         u = 3.813
         u = roundN(u,7)
         v = 0.2918
         v = roundN(v,7)
         w = 1.398
         w = roundN(w,7)
         z = 0.0017
         z = roundZN(z,7)
         DO i = 1, n
         Y(i) = (( a * dexp1(-b*X(i)**2/PP**2/16)
     *         + c * dexp1(-d*X(i)**2/PP**2/16)
     *         + e * dexp1(-f*X(i)**2/PP**2/16) + g)
     *         + (1-2/1.2)
     *     * ( r * X(i)**2/PP**2/16 * dexp1(-s*X(i)**2/PP**2/16)
     *      +  t * X(i)**2/PP**2/16 * dexp1(-u*X(i)**2/PP**2/16)
     *      +  v * X(i)**2/PP**2/16 * dexp1(-w*X(i)**2/PP**2/16)
     *      +  z * X(i)/PP/4))**2 *
     *       P(1)*(1 - dquot0( dsin(X(i)*P(2)), X(i)*P(2)) +
     *           0.5 * P(3) *
     *           dquot0( dsin(X(i)*P(4)) , X(i)*P(4))  -
     *           0.25 * P(3) *
     *           dquot0( dsin(X(i)*(P(2) + P(4))) , X(i)*(P(2)+P(4)))
     *           - 0.25 * P(3) *
     *           dquot0( dsin(X(i)*(P(4)- P(2))), X(i)*(P(4) - P(2))))
     *           + P(5)+ P(6) * X(i)**2
            ENDDO
      ELSEIF (ifc.eq.49) THEN ! FFCu2p
         PP = 3.141592654
         PP = roundN(PP,7)
C         Print *, 'PP=',PP
         a = 0.0232
         a = roundN(a,7)
C         Print *, 'a=', a
         b = 34.969
         b = roundN(b,7)
C         Print *, 'b=', b
         c = 0.4023
         c = roundN(c,7)
C         Print *, 'c=', c
         d = 11.564
         d = roundN(d,7)
         e = 0.5882
         e = roundN(e,7)
         f = 3.843
         f = roundN(f,7)
         g = -0.0137
         g = roundN(g,7)
C         Print *, 'D=', D
         r = 1.5189
         r = roundN(r,7)
         s = 10.478
         s = roundN(s,7)
         t = 1.1512
         t = roundN(t,7)
         u = 3.813
         u = roundN(u,7)
         v = 0.2918
         v = roundN(v,7)
         w = 1.398
         w = roundN(w,7)
         z = 0.0017
         z = roundZN(z,7)
         DO i = 1, n
         Y(i) = P(1) * (( a * dexp1(-b*X(i)**2/PP**2/16)
     *         + c * dexp1(-d*X(i)**2/PP**2/16)
     *         + e * dexp1(-f*X(i)**2/PP**2/16) + g)
     *         + (1-2/1.2)
     *     * ( r * X(i)**2/PP**2/16 * dexp1(-s*X(i)**2/PP**2/16)
     *      +  t * X(i)**2/PP**2/16 * dexp1(-u*X(i)**2/PP**2/16)
     *      +  v * X(i)**2/PP**2/16 * dexp1(-w*X(i)**2/PP**2/16)
     *      +  z * X(i)/PP/4))
            ENDDO
      ELSEIF (ifc.eq.50) THEN ! chain+Phonoenfit
         f = 8314.51
         f = roundN(f,7)
         h1 = 1.
         DO i = 1, n
            Y(i) = 2 * P(1) * f * X(i) / P(2) / 3
     *             + P(3)*X(i)**3 + P(4)*X(i)**5
            ENDDO
      ELSEIF (ifc.eq.51) THEN ! Curie Weiss
         g1 = 0.1250486
         g1 = roundN(g1,7)
         h1 = 1.
         DO i = 1, n
            Y(i) = P(1)* P(2)**2 * g1*P(3)*dquot0(P(3)+h1,(X(i)-P(4))**
     *      P(6)) +  P(5)
            ENDDO
      ELSEIF (ifc.eq.52) THEN ! Bonner Fischer a = 1
         a = .25
         b = .07498
         b = roundN(b,7)
         c = 0.07524
         c = roundN(c,7)
         d = 0.99310
         d = roundN(d,7)
         e = 0.17214
         e = roundN(e,7)
         f = 0.75783
         f = roundN(f,7)
         g = 0.3751458
         g = roundN(g,7)
         g1 = 0.1250486
         g1 = roundN(g1,7)
         h1 = 1.
         h2 = 2.
         h3 = 3.
         DO i = 1, n
            Y(i) = P(1) * P(2)**2 * dquot0(g,X(i)) *
     *      dquot0(a + b * dquot0(P(3),X(i)) +
     *      c * dpow0 (dquot0(P(3), X(i)),h2),(
     *      h1 + d * dquot0(P(3),X(i)) + e *
     *      dpow0 (dquot0(P(3), X(i)),h2) + f *
     *      dpow0 (dquot0(P(3), X(i)),h3))) + P(4) +
     *      P(5)* P(6)**2 * g1*P(7)*dquot0(P(7)+h1,(X(i)-P(8))**P(9))
            ENDDO
      ELSEIF (ifc.eq.53) THEN ! Dimer alpha = 0
 1       a = 0.3751458
         a = roundN(a,7)
         b = 1.
         c = 3.
         g1 = 0.1250486
         g1 = roundN(g1,7)
         h1 = 1.
         DO i = 1,n
            Y(i) = P(1) * P(2)**2 * dquot0(a,X(i))*
     *      dquot0 (b, c + dexp1(P(3)/X(i))) + P(4) +
     *      P(5)* P(6)**2 * g1 * P(7)*dquot0(P(7)+h1,(X(i)-P(8))**P(9))
            ENDDO
      ELSEIF (ifc.eq.54) THEN ! alternating chain Hall et al. ..
C         aa = 0.3751458
C         aa = roundN(a,7)
C         Print *, 'aa =',aa
         g1 = 0.1250486
         g1 = roundN(g1,7)
         h1 = 1.
         IF (P(1).ge.0.and.P(1).le.0.4) THEN
            DO i = 1,n
               Y(i) = P(2) * P(3)**2 * 0.3751458/X(i) *
     *           AlterC(X(i),P(1),P(4),.true.) + P(5) +
     *          P(6)* P(7)**2 * g1*P(8)*dquot0(P(8)+h1,(X(i)-P(9))
     *         **P(10))
               ENDDO
         ELSEIF (P(1).gt.0.4.and.P(1).le.1.0) THEN
            DO i = 1,n
               Y(i) = P(2) * P(3)**2 * 0.3751458/X(i) *
     *           AlterC(X(i),P(1),P(4),.false.) + P(5) +
     *          P(6)* P(7)**2 * g1*P(8)*dquot0(P(8)+h1,(X(i)-P(9))
     *          **P(10))
               ENDDO
            ENDIF
       ELSEIF (ifc.eq.55) THEN ! alternating chain Hall et al. ..
C         aa = 0.3751458
C         aa = roundN(a,7)
C         Print *, 'aa =',aa
         g1 = 0.1250486
         g1 = roundN(g1,7)
         h1 = 1.
         IF (P(1).ge.0.and.P(1).le.0.4) THEN
            DO i = 1,n
               Y(i) = (P(2)-P(6))* P(3)**2 * 0.3751458/X(i) *
     *           AlterC(X(i),P(1),P(4),.true.) + P(5) +
     *          P(6)* P(7)**2 * g1*P(8)*dquot0(P(8)+h1,(X(i)-P(9))
     *         **P(10))
               ENDDO
         ELSEIF (P(1).gt.0.4.and.P(1).le.1.0) THEN
            DO i = 1,n
               Y(i) = (P(2)-P(6)) * P(3)**2 * 0.3751458/X(i) *
     *           AlterC(X(i),P(1),P(4),.false.) + P(5) +
     *          P(6)* P(7)**2 * g1*P(8)*dquot0(P(8)+h1,(X(i)-P(9))
     *          **P(10))
               ENDDO
            ENDIF
      ELSEIF (ifc.eq.56) THEN ! Curie Weiss invers
         g1 = 0.1250486
         g1 = roundN(g1,7)
         h1 = 1.
         DO i = 1, n
            Y(i) = 1/P(1)/ P(2)**2 / g1/P(3)*dquot0((X(i)-P(4))** P(6),
     *      P(3)+h1) +  P(5)
            ENDDO
      ELSEIF (ifc.eq.57) THEN ! alternating chain Hall et al. special
C                               J-J_s
C         a = 0.3751458
C         a = roundN(a,7)
         g1 = 0.1250486
         g1 = roundN(g1,7)
         h1 = 1.
         IF (P(1).ge.0.and.P(1).le.0.4) THEN
            DO i = 1,n
               Y(i) = P(2) * P(3)**2 * 0.3751458/X(i) *
     *           AlterS1(X(i),P(1),P(4),P(5),.true.) + P(6) +
     *          P(7)* P(8)**2 * g1*P(9)*dquot0(P(9)+h1,(X(i)-P(10))
     *          **P(11))
               ENDDO
         ELSEIF (P(1).gt.0.4.and.P(1).le.1.0) THEN
            DO i = 1,n
               Y(i) = P(2) * P(3)**2 * 0.3751458/X(i) *
     *           AlterS1(X(i),P(1),P(4),P(5),.false.) + P(6) +
     *          P(7)* P(8)**2 * g1*P(9)*dquot0(P(9)+h1,(X(i)-P(10))
     *         **P(11))
               ENDDO
            ENDIF
      ELSEIF (ifc.eq.60) THEN ! bg+Lorentz+elastic
         pi = twopi/2
         DO i = 1,n
            Y(i) = P(1) + P(2) * dquot0(P(3)*
     *        dabs(P(5))/pi,(X(i)-P(4))**2+P(5)**2)
            ENDDO
         iXL = irPosOpt (X, n, P(6), 'r', iXL)
         iXR = iXL + 1
         IF (iXL.lt.0 .or. iXR.gt.n) RETURN
         dX = X(iXR) - X(iXL)
         IF (dX.lt.1d-20) RETURN
         reR = (P(6)-X(iXL))/dX
         Y(iXL) = Y(iXL) + P(2) * (1-P(3)) * (1-reR) / dX
         Y(iXR) = Y(iXR) + P(2) * (1-P(3)) *   reR   / dX
      ELSEIF (ifc.eq.61) THEN ! Lorentzian
         pi = twopi/2
         DO i = 1, n
            Y(i) = dquot0(P(1)*dabs(P(3))/pi,(X(i)-P(2))**2+P(3)**2)
            ENDDO
      ELSEIF (ifc.eq.62) THEN ! 2 Lorentzians
         pi = twopi/2
         DO i = 1, n
            Y(i) = dquot0(P(1)*dabs(P(2))/pi,(X(i)-P(5))**2+P(2)**2)
     *            +dquot0(P(3)*dabs(P(4))/pi,(X(i)-P(5))**2+P(4)**2)
            ENDDO
      ELSEIF (ifc.eq.63) THEN ! Lorentzian + elastic
         pi = twopi/2
         DO i = 1, n
            Y(i) =
     * P(1) * dquot0(P(2)*dabs(P(4))/pi,(X(i)-P(3))**2+P(4)**2)
            ENDDO
         iXL = irPosOpt (X, n, P(5), 'r', iXL)
         iXR = iXL + 1
         IF (iXL.lt.0 .or. iXR.gt.n) RETURN
         dX = X(iXR) - X(iXL)
         IF (dX.lt.1d-20) RETURN
         reR = (P(5)-X(iXL))/dX
         Y(iXL) = Y(iXL) + P(1) * (1-P(2)) * (1-reR) / dX
         Y(iXR) = Y(iXR) + P(1) * (1-P(2)) *   reR   / dX
      ELSEIF (ifc.eq.64) THEN ! bg + Lorentzian
         DO i = 1, n
            Y(i) =
     *   P(4) + P(5)*X(i) + P(6)*X(i)**2 + P(7)*X(i)**3 + P(8)*X(i)**4
     *                    + P(9)*X(i)**5 + P(10)*X(i)**6
     *   + dquot0 ( P(1), (X(i)-P(2))**2 + P(3)**2 )
            ENDDO
      ELSEIF (ifc.eq.65) THEN ! bg + 5 Lorentzians
         DO i = 1, n
            Y(i) = P(1) + P(2)*X(i) + P(3)*X(i)**2
     *             + dquot0 ( P( 4), (X(i)-P( 5))**2 + P( 6)**2 )
     *             + dquot0 ( P( 7), (X(i)-P( 8))**2 + P( 9)**2 )
     *             + dquot0 ( P(10), (X(i)-P(11))**2 + P(12)**2 )
     *             + dquot0 ( P(13), (X(i)-P(14))**2 + P(15)**2 )
     *             + dquot0 ( P(16), (X(i)-P(17))**2 + P(18)**2 )
            ENDDO
      ELSEIF (ifc.eq.66) THEN  ! Two power laws + Lorentzian
         DO i = 1, n
            xred  = X(i) / P(3)
            Y(i)  = P(1)
     *     + P(2) * ( dpow0 (xred,  P(4)) + dpow0 (xred, -P(5)) )
     *     + dquot0 ( P(6), (X(i)-P(7))**2 + P(8)**2 )
            ENDDO
      ELSEIF (ifc.eq.67) THEN  ! Lorentzian^n
         DO i = 1, n
            Y(i)  = P(1) *  dpow0 (dquot0 (P(3)**2,
     *                      (X(i)-P(2))**2 + P(3)**2 ), P(4) )
            ENDDO
      ELSEIF (ifc.eq.68) THEN ! Polynome bg + Lorentzian
         DO i = 1, n
            Y(i) =
     *   P(6) + P(7)*X(i) + P(8)*X(i)**2 + P(9)*X(i)**3 + P(10)*X(i)**4
     *   + dquot0 ( P(1)+P(2)*X(i)+P(3)*X(i)**2,
     *              (X(i)-P(4))**2 + P(5)**2 )
            ENDDO
      ELSEIF (ifc.eq.69) THEN ! bg + Gauss + Lorentz
         pi = twopi / 2
         DO i = 1,n
             Y(i) =
     *       P(1) + dquot0( P(2)*dabs(P(4))/pi,(X(i)-P(3))**2 +
     *       P(4)**2) + P(5) * dexp1(- dquot0(X(i)-P(6),P(7))**2)
            ENDDO
      ELSEIF (ifc.eq.70) THEN ! sloping bg + Gauss + Lorentz
         pi = twopi / 2
         DO i = 1,n
             Y(i) =
     *       P(1) + P(2)*X(i) + dquot0( P(3)*
     *       dabs(P(5))/pi,(X(i)-P(4))**2 + P(5)**2)
     *        + P(6) * dexp1(- dquot0(X(i)-P(7),P(8))**2)
            ENDDO
      ELSEIF (ifc.eq.71) THEN ! 2 * Lorentz doublet
C         pi = twopi / 2
         DO i = 1,n
            Y(i) = P(1) + P(2) / ((X(i)-P(3))**2 + P(5)**2) + P(2) /
     *            ((X(i)-P(4))**2 + P(5)**2) + P(6) / ((X(i)-P(7))**2 +
     *            P(9)**2) + P(6) / ((X(i)-P(8))**2 + P(9)**2)
            ENDDO
      ELSEIF (ifc.eq.72) THEN ! Lorentzian (amplitude = value at max)
         pi = twopi/2.
         DO i = 1, n
            Y(i) = dquot0(P(1)*P(3)**2,(X(i)-P(2))**2+P(3)**2)
            ENDDO
      ELSEIF (ifc.eq.73) THEN ! S(Q)_dimer
         DO i = 1, n
         Y(i) = P(1) * (1 - dquot0( dsin(X(i)*P(2)), X(i)*P(2)))
     *            +  P(3)
            ENDDO
      ELSEIF (ifc.eq.74) THEN ! S(Q)_dimer_ach
         DO i = 1, n
         Y(i) = P(1)*(1 - dquot0( dsin(X(i)*P(2)), X(i)*P(2)) +
     *           0.5 * P(3) *
     *           dquot0( dsin(X(i)*P(4)) , X(i)*P(4))  -
     *           0.25 * P(3) *
     *           dquot0( dsin(X(i)*(P(2) + P(4))) , X(i)*(P(2)+P(4)))
     *           - 0.25 * P(3) *
     *           dquot0( dsin(X(i)*(P(4)- P(2))), X(i)*(P(4) - P(2))))
     *           + P(5)
            ENDDO
      ELSEIF (ifc.eq.78) THEN ! Pd catalyst and THF
         pi = twopi/2.
         DO i = 1, n
            Y(i) = P(1) + P(2)/dsqrt0(pi) / P(3) *
     *             dexp1(- dquot0(X(i)-P(8),P(3))**2) +
     *             dquot0(P(4)*dabs(P(5))/pi,(X(i)-P(8))**2+P(5)**2) +
     *             dquot0(P(6)*dabs(P(7))/pi,(X(i)-P(8))**2+P(7)**2)
            ENDDO
      ELSEIF (ifc.eq.79) THEN ! special function Pd catalyst
         pi = twopi/2.
         FNU = P(1)*P(3)
         EISF = 0.333*(1. + 2.*dbesi0(FNU)) !Artem: Replace with dbesi0 function from slatec/fnlib/dbesi0.f:  S17AEF (FNU,IFAIL)
         DWF  = dexp1(- P(1)**2 * 0.3333 * P(2))
         DO i = 1,n
            Y(i) = DWF * ( dquot0( (P(5)*EISF+P(6)+P(10)) *
     *      dabs(P(11))/pi,(X(i)-P(15))**2+P(11)**2) + P(5)*(1.-EISF)*
     *      (dquot0(P(7)*dabs(P(12))/pi,(X(i)-P(15))**2+P(12)**2)+
     *      dquot0(P(8)*dabs(P(13))/pi,(X(i)-P(15))**2+P(13)**2)+
     *      dquot0(P(9)*dabs(P(14))/pi,(X(i)-P(15))**2+P(14)**2)
     *      )) + P(4)
            ENDDO
      ELSEIF (ifc.eq.81) THEN ! Kohlrausch
         DO i = 1, n
            Y(i) =
     * P(4) + P(1) * dexp1 ( - dpow0 (dquot0(X(i),P(2)),P(3)) )
            ENDDO
      ELSEIF (ifc.eq.82) THEN ! improved Kohlrausch 15dec95
         tt = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         DO i = 1, n
            Y(i) = P(1) * dexp1 ( - dpow0 (dquot0(X(i),tt),P(3)) )
            ENDDO
      ELSEIF (ifc.eq.83) THEN ! -d/dt Kohlrausch
         DO i = 1, n
            Y(i) =
     * P(1) * dquot0(P(3),X(i)) * dpow0(dquot0(X(i),P(2)),P(3)) *
     *                  dexp1 ( - dpow0 (dquot0(X(i),P(2)),P(3)) )
            ENDDO
      ELSEIF (ifc.eq.84) THEN ! improved Kohlrausch * exp.
         tt = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         DO i = 1, n
            Y(i) = P(1) * dexp1 ( - dpow0 (dquot0(X(i),tt),P(3)) )
     *             * dexp1 ( -dquot0(X(i),P(4)) )
            ENDDO
      ELSEIF (ifc.eq.85) THEN ! improved Kohlrausch + background
         tt = dquot0 (P(4) * P(3), dgamma1(dquot0(1.d0,P(4))))
         DO i = 1, n
            Y(i) =
     * P(1) + P(2) * dexp1 ( - dpow0 (dquot0(X(i),tt),P(4)) )
            ENDDO
      ELSEIF (ifC.eq.86) THEN !improved Kohlrausch * 2
         tt1 = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         tt2 = dquot0 (P(6) * P(5), dgamma1(dquot0(1.d0,P(6))))
         DO i = 1,n
            Y(i) = P(1) * dexp1 ( - dpow0 (dquot0(X(i),tt1),P(3)) )
     *             + P(4) * dexp1 ( - dpow0 (dquot0(X(i),tt2),P(6)) )
            ENDDO
      ELSEIF (ifc.eq.87) THEN !special function Exp + Gotze + 2xKohl
         tt1 = dquot0 (P(8) * P(9), dgamma1(dquot0(1.d0,P(9))))
         tt2 = dquot0 (P(11) * P(12), dgamma1(dquot0(1.d0,P(12))))
          DO i = 1, n
            IF (P(5).le.0.) THEN
               Y(i) = 1.d4 ! punish Newton
            ELSE
               tt = dquot0(X(i), P(5))
               Y(i)  = P(1) * dexp1( - dquot0(X(i),P(2)) ) + P(3) +
     *            P(4) * GoetzeG (tt, P(6), .true.) + P(7) *
     *            dexp1 ( - dpow0 (dquot0(X(i),tt1),P(9))) + P(10)*
     *            dexp1 ( - dpow0 (dquot0(X(i),tt2),P(12)))
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.88) THEN ! impr Kohl + bg constant a+b
         tt = dquot0 (P(4) * P(3), dgamma1(dquot0(1.d0,P(4))))
         DO i = 1, n
            Y(i) = P(1)-P(2)+P(2)*dexp1( - dpow0(dquot0(X(i),tt),P(4)))
            ENDDO
      ELSEIF (ifc.eq.89) THEN ! improved Kohlrausch + pos. background
         tt = dquot0 (P(4) * P(3), dgamma1(dquot0(1.d0,P(4))))
         DO i = 1, n
            Y(i) =
     * dabs (P(1)) + P(2) * dexp1 ( - dpow0 (dquot0(X(i),tt),P(4)) )
            ENDDO
      ELSEIF (ifc.eq.91) THEN ! ln VFT
         DO i = 1, n
            Y(i) = P(1) - dquot0 ( P(2), X(i)-P(3) )
            ENDDO
      ELSEIF (ifc.eq.92) THEN ! ln VFT + ...
         DO i = 1, n
            Y(i) = P(1) - dquot0 (P(2), X(i)-P(3))
     *                  -(dquot0 (P(4), X(i)-P(5)))**2
            ENDDO
      ELSEIF (ifc.eq.95) THEN ! order parameter critical
         DO i = 1, n
            Y(i) = P(1)*dpow0((1.0 - dquot0(X(i),P(2))),P(3))
            ENDDO
      ELSEIF (ifc.eq.96) THEN ! Bloch T^3/2 law
         DO i = 1,n
            Y(i) = P(1)*(1.0 - P(2) * (X(i)**1.5))
            ENDDO
      ELSEIF (ifc.eq.101) THEN  ! bg
         halfmass = 5.2271  ! 0.5*mass_of_neutron [meV,millisec,m]
         fact     = halfmass * dsqrt (halfmass/P(2))
         DO i = 1, n
            tof  = Tau_of_W (X(i), P(2))
            IF (tof.gt.0) THEN
               dwdt = 2*fact/tof**4
               Y(i) = P(1) / dwdt
            ELSE
               Y(i) = 0
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.102) THEN  ! bg(w) - bg(-w)
         halfmass = 5.2271  ! 0.5*mass_of_neutron [meV,millisec,m]
         fact     = halfmass * dsqrt (halfmass/P(2))
         DO i = 1, n
            tofP  = Tau_of_W (X(i), P(2))
            dwdtP = 2*fact/tofP**4
            tofM  = Tau_of_W (-X(i), P(2))
            dwdtM = 2*fact/tofM**4
            Y(i)  = P(1) / dwdtP - P(1) / dwdtM
            ENDDO
      ELSEIF (ifc.eq.103) THEN  ! INS Detector Efficiency
         E0 = P(3)
         IF (E0.gt.1.d-8) THEN
            eff0 = dexp1(P(1)/dsqrt(E0)) * (1-dexp1(P(2)/dsqrt(E0)))
         ELSE
            eff0 = 1
            ENDIF
         DO i = 1, n
            E = E0 + X(i)
            IF (E.gt.P(4)) THEN
               eff = P(5) * (1-dexp1(P(2)/dsqrt(E)))
            ELSEIF (E.gt.1.d-8) THEN
               eff = dexp1(P(1)/dsqrt(E)) * (1-dexp1(P(2)/dsqrt(E)))
            ELSE
               eff = 1
               ENDIF
            Y(i) = eff0/eff
            ENDDO
      ELSEIF (ifc.eq.106) THEN ! FPI-Reflexion
         DO i = 1, n
            Y(i) = P(1)
            nord = idnint(P(10))
            wid1 = dabs(P(6))
            wid2 = dabs(P(6)*P(7))
            DO ii = -nord, nord
               x1 = X(i) + ii*P(8)      - P(4)
               x2 = X(i) + ii*P(8)*P(9) - P(5)
               Y(i) = Y(i) - dquot0(P(2)*wid1, x1**2 + wid1**2)
     *                     - dquot0(P(3)*wid2, x2**2 + wid2**2)
               ENDDO
            ENDDO
      ELSEIF (ifc.eq.107) THEN ! FPI-Transmission
         fin = 16 * P(2)**2 / twopi**2
         DO i = 1, n
            Y(i) = dquot0 (P(1),
     * (1+fin*dsin(P(3)*(X(i)-P(4))     )**2)**3 *
     * (1+fin*dsin(P(3)*(X(i)-P(4))*P(5))**2)**3 )
            ENDDO
      ELSEIF (ifc.eq.108) THEN  ! Lorentzian^n * Saturation
         DO i = 1, n
            Y(i)  = P(1) *  dpow0 (dquot0 (dabs(P(3)),
     *                      (X(i)-P(2))**2 + P(3)**2 ), P(4) )
            Y(i)  = Y(i) * dexp1 (-P(5)*Y(i))
            ENDDO
      ELSEIF (ifc.eq.109) THEN  ! resolution FPI
                          ! 'p1*R(x-p2-n*p4;p3)^3*R(..-..*p5-p7,p5*p6*p3)'
                          ! 'A;w0;G;P;rG;rF;dw'
         DO i = 1, n
            Y(i)  = P(1)
     *      * (  dLorentz(X(i)-P(2),               P(3)          )**3
     *         + dLorentz(X(i)-P(2)     -P(4),     P(3)          )**3
     *         + dLorentz(X(i)-P(2)     +P(4),     P(3)          )**3 )
     *      * (  dLorentz(X(i)-P(2)-P(7),          P(3)*P(5)*P(6))**3
     *         + dLorentz(X(i)-P(2)-P(7)-P(4)*P(5),P(3)*P(5)*P(6))**3
     *         + dLorentz(X(i)-P(2)-P(7)+P(4)*P(5),P(3)*P(5)*P(6))**3 )
            ENDDO
      ELSEIF (ifc.eq.110) THEN  ! resolution FPI * saturation
         DO i = 1, n
            Y(i)  = P(1)
     *      * (  dLorentz(X(i)-P(2),               P(3)          )**3
     *         + dLorentz(X(i)-P(2)     -P(4),     P(3)          )**3
     *         + dLorentz(X(i)-P(2)     +P(4),     P(3)          )**3 )
     *      * (  dLorentz(X(i)-P(2)-P(7),          P(3)*P(5)*P(6))**3
     *         + dLorentz(X(i)-P(2)-P(7)-P(4)*P(5),P(3)*P(5)*P(6))**3
     *         + dLorentz(X(i)-P(2)-P(7)+P(4)*P(5),P(3)*P(5)*P(6))**3 )
            Y(i)  = Y(i) * dexp1 (-P(8)*Y(i))
            ENDDO
      ELSEIF (ifc.eq.111) THEN ! echo (rectang lambda -> sinx / x)
         DO i = 1, n
            dx = X(i) - P(5)
            argosc = dquot0 (twopi*dx, P(3))
            argenv = dquot0 (twopi*dx, P(3)*P(4))
            Y(i) = P(1)/2 * ( 1 + P(2)*dcos (argosc)*
     *                dquot1(dsin(argenv),argenv))
            ENDDO
      ELSEIF (ifc.eq.112) THEN ! echo (triang lambda -> (sinx / x)^2 )
         DO i = 1, n
            dx = X(i) - P(5)
            argosc = dquot0 (twopi*dx, P(3))
            argenv = dquot0 (twopi*dx, P(3)*P(4))
            Y(i) = P(1)/2 * ( 1 + P(2)*dcos (argosc)*
     *                dquot1(dsin(argenv),argenv)**2)
            ENDDO
      ELSEIF (ifc.eq.113) THEN
         ! Name    = 'diffraction by grid'
         ! Formula = '[(1/p1)(sin p2 x p1 / sin p2 x) sinc p3 x]^2'
         nn = idnint(P(1))
         DO i = 1, n
            IF (nn.lt.1 .or. (nn.gt.1 .and. P(2).le.0) .or.
     *         P(3).le.0) THEN
               Y(i) = 0
            ELSE
               ak = P(2)*X(i)
               IF (nn.eq.1) THEN
                  fac2 = 1
               ELSEIF (dabs(dsin(ak)).lt.1d-30) THEN
                  fac2 = nn * dcos(nn*ak) / dcos(ak)
               ELSE
                  fac2 = dsin(nn*ak) / dsin(ak)
                  ENDIF
               ENDIF
            bk = P(3)*X(i)
            IF (dabs(bk).lt.1d-30) THEN
               fac3 = 1
            ELSE
               fac3 = dsin(bk) / bk
               ENDIF
            Y(i) = ((1.d0/nn) * fac2 * fac3)**2
            ENDDO
      ELSEIF (ifc.eq.121) THEN  ! Sjoegren interpolation
         pref = dquot0 (P(1), P(3)+P(4))
         DO i = 1, n
            w     = dquot0 (X(i), P(2))
            Y(i)  = pref * ( P(4)*dpow0(w,P(3)) + P(3)*dpow0(w,-P(4)) )
            ENDDO
      ELSEIF (ifc.eq.122) THEN  ! The two power laws in t
         DO i = 1, n
            Y(i)  = P(1) + P(2) * ( P(3) * dpow0 (X(i), -P(5))
     *                             - dpow0 (dquot0(X(i),P(4)), P(6)) )
            ENDDO
      ELSEIF (ifc.eq.123) THEN  ! Goetze f+h*G(t)
         DO i = 1, n
            IF (P(3).le.0.) THEN
               Y(i) = 1.d4 ! punish Newton
            ELSE
               tt = dquot0(X(i), P(3))
               Y(i)  = P(1) + P(2) * GoetzeG (tt, P(4), .true.)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.124) THEN  ! Goetze f+h*G(t) .. short
         DO i = 1, n
            IF (P(3).le.0.) THEN
               Y(i) = 1.d4 ! punish Newton
            ELSE
               tt = dquot0(X(i), P(3))
               Y(i)  = P(1) + P(2) * GoetzeG (tt, P(4), .false.)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.125) THEN  ! Goetze X(w)
         DO i = 1, n
            IF (P(2).le.0.) THEN
               Y(i) = 1.d4
            ELSE
               ww = dquot0(X(i), P(2))
               Y(i)  = P(1) * GoetzeX (ww, P(3), .true., .false.)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.126) THEN  ! Goetze X(w) short
         DO i = 1, n
            IF (P(2).le.0.) THEN
               Y(i) = 1.d4
            ELSE
               ww = dquot0(X(i), P(2))
               Y(i)  = P(1) * GoetzeX (ww, P(3), .false., .false.)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.127) THEN  ! Goetze S(w)
         DO i = 1, n
            IF (P(2).le.0.) THEN
               Y(i) = 1.d4
            ELSE
               ww = dquot0(X(i), P(2))
               Y(i)  =
     * dquot0(P(1), X(i)) * GoetzeX (ww, P(3), .true., .false.)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.128) THEN  ! Goetze f+h*G(t) mit fixem t_sig zu h_q
         CALL GoetzeCoeff (P(4), a, b, A1, A2, A3, B0, B1, tStar,
     *         A0ss, A1ss, A2ss, A3ss, B0ss, B1ss, wCro, wMin, xMin)
         tsig = P(3) * dpow0 (P(2), -dquot0(1.d0, a))
         DO i = 1, n
            IF (tsig.le.0.) THEN
               Y(i) = 1.d4 ! punish Newton
            ELSE
               tt = dquot0(X(i), tsig)
               Y(i)  = P(1) + P(2) * GoetzeG (tt, P(4), .true.)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.129) THEN  ! Goetze X(w)
         DO i = 1, n
            IF (P(2).le.0.) THEN
               Y(i) = 1.d4
            ELSE
               ww = dquot0(X(i), P(2))
               Y(i)  = P(1) * GoetzeX (ww, P(3), .true., .true.)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.131) THEN  ! sqrt law
         DO i = 1, n
            sigma = dquot0((P(3)-X(i)),P(3))
            IF (sigma.gt.0.) THEN
               Y(i) = P(1) + P(2) * dsqrt(sigma) + P(4) * sigma
            ELSE
               Y(i) = P(1) + P(4) * sigma
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.132) THEN  ! sqrt law / two slopes
         DO i = 1, n
            sigma = dquot0((P(3)-X(i)),P(3))
            IF (sigma.gt.0.) THEN
               Y(i) = P(1) + P(2) * dsqrt(sigma) + P(4) * sigma
            ELSE
               Y(i) = P(1) + P(5) * sigma
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.133) THEN  ! sqrt law
                               ! 'p1*e^(p4(p3-x)) + (if x>p3)p2*sqrt((p3-x)/p3)'
         DO i = 1, n
            sigma = dquot0((P(3)-X(i)),P(3))
            IF (sigma.gt.0.) THEN
               Y(i) = P(1)*dexp1(P(4)*(P(3)-X(i))) + P(2)*dsqrt(sigma)
            ELSE
               Y(i) = P(1)*dexp1(P(4)*(P(3)-X(i)))
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.141) THEN  ! Debye <u^2> (from WuL6)
         DO i = 1, n
            Y(i) =   u2Debye (X(i), P(1), P(2), P(3))
            ENDDO
      ELSEIF (ifc.eq.142) THEN  ! Debye <u^2> (from WuL6)
         DO i = 1, n
            Y(i) =   u2Debye (X(i), P(1), P(2), P(4))
     *             - u2Debye (P(3), P(1), P(2), P(4))
            ENDDO
      ELSEIF (ifc.eq.143) THEN  ! DWF from Debye <u^2> (from WuL6)
         DO i = 1, n
            Y(i) =   u2Debye (X(i), P(1), P(2), P(4))
     *             - u2Debye (P(3), P(1), P(2), P(4))
            Y(i) = P(6) * dexp1(-P(5)**2 * Y(i))
            ENDDO
      ELSEIF (ifc.eq.144) THEN  ! 'm-phonon Q-integral'
         DO i = 1, n
            uQ2 = X(i)**2 * P(1)
            dwf = dexp (-uQ2)
            sum = dwf
            fac = 1.
            DO m = 1, idnint(P(2))
               fac = fac * m
               sum = sum + uQ2**m / fac * dwf
               ENDDO
            Y(i) = sum
            ENDDO
      ELSEIF (ifc.eq.145) THEN  ! 'm-phonon contribution'
         fac = 1.
         m = idnint(P(3))
         DO mm = 1, m
            fac = fac * mm
            ENDDO
         DO i = 1, n
            uQ2 = X(i)**2 * P(2)
            Y(i)= P(1) * dexp1(-uQ2) * uQ2**m / fac
            ENDDO
      ELSEIF (ifc.eq.146) THEN  ! 'Debye C(T)'
         ! programmiert wie Sau, den 23jul92
         nIntegr = idnint(P(3))
         DO i = 1, n
            IF (X(i).le.0.) THEN
               Y(i) = 0.
            ELSE
               xx = P(1) / X(i)     ! TD / T
               ! integrate => D(x)
               DD = 0.
               dx = xx / nIntegr
               DO ii = 1, nIntegr
                  xi = ii * xx / nIntegr
                  DD = DD + dx * 3/xx**3 * dquot0 (xi**3, dexp1(xi)-1)
                  ENDDO
               Y(i) = 3 * 8.3144 * P(2) *
     *            ( 4 * DD - 3 * dquot0 (xx, dexp1(xx)-1) )
               ! see LL V (66,8)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.147) THEN
            ! q,m,T in A-1,amu,K; X in meV
         ER = dquot0 (389.387d10*P(2), 931.48d9*P(3)) ! hbar^2*q^2/2M
         sig2 = dquot0 (389.387d10*P(2) * .08617*P(4), 931.48d9*P(3))
                   ! hbar^2*q^2*kT/M
         IF (sig2.le.0) sig2 = 1
         B = 2 * sig2
         A = P(1) / dsqrt(twopi) / dsqrt(sig2)
         DO i = 1, n
            Y(i) = A * dexp1(-(X(i)-ER)**2/B)
            ENDDO
      ELSEIF (ifc.eq.151) THEN
         DO i = 1, n
            Y(i) = P(1) * X(i)**2 * dexp1 ( - P(2)**2 * X(i)**2 )
            ENDDO
      ELSEIF (ifc.eq.152) THEN
         DO i = 1, n
            u    = X(i)/P(1)
            Y(i) = (P(2)*u**2 + P(3)*u**3 + P(4)*u**4 + P(5)*u**5)
     *              * dexp (-u**2/2)
            ENDDO
      ELSEIF (ifc.eq.153) THEN
         DO i = 1, n
            u    = X(i)/P(1)
            Y(i) = (P(2)*u**2 + P(3)*u**3 + P(4)*u**4 + P(5)*u**5)
     *              * dexp (-u)
            ENDDO
      ELSEIF (ifc.eq.154) THEN  ! Debye-Gauss
         DO i = 1, n
            xi   = X(i)
            yi   = 0.
            DO m = 1, 10
               IF (P(2*m).ne.0.)
     *         yi = yi + P(2*m-1)*dexp(-xi**2/2/P(2*m)**2)/P(2*m)**3
               ENDDO
            Y(i) = yi * xi**2
            ENDDO
      ELSEIF (ifc.eq.155) THEN  ! Debye-Lorentz
         DO i = 1, n
            xi   = X(i)
            yi   = 0.
            DO m = 1, 10
               IF (P(2*m).ne.0.)
     *         yi = yi + P(2*m-1) * dexp(-xi/P(2*m)) / P(2*m)**3
               ENDDO
            Y(i) = yi * xi**2
            ENDDO
      ELSEIF (ifc.eq.156) THEN  ! sharp Debye
         DO i = 1, n
            xi   = X(i)
            yi   = 0.
            DO m = 1, 10
               IF (xi.le.P(2*m) .and. P(2*m).gt.0.)
     *            yi = yi + P(2*m-1)*xi**2/P(2*m)**3
               ENDDO
            Y(i) = yi
            ENDDO
      ELSEIF (ifc.eq.157) THEN ! Lorentzian + elastic
         pi = twopi/2
         DO i = 1, n
            Y(i) = P(1)
     *   * ( dquot0(dabs(P(2))*dabs(P(4))/pi,(X(i)-P(3))**2+P(4)**2)
     *   + dquot0(dabs(P(5))*dabs(P(7))/pi,(X(i)-P(6))**2+P(7)**2)
     *   + dquot0(dabs(P(8))*dabs(P(10))/pi,(X(i)-P(9))**2+P(10)**2)
     *   + dquot0(dabs(P(11))*dabs(P(13))/pi,(X(i)-P(12))**2+P(13)**2)
     *   + P(14) + P(15) * X(i) + P(16) * X(i)**2 + P(17) *X(i)**3)
         ENDDO
         iXL = irPosOpt (X, n, P(18), 'r', iXL)
         iXR = iXL + 1
         IF (iXL.lt.0 .or. iXR.gt.n) RETURN
         dX = X(iXR) - X(iXL)
         IF (dX.lt.1d-20) RETURN
         reR = (P(18)-X(iXL))/dX
         Y(iXL) = Y(iXL)
     *          + P(1)*P(19)*(1-dabs(P(2))-dabs(P(5))-dabs(P(8))
     *          - dabs(P(11))) * (1-reR) / dX
         Y(iXR) = Y(iXR)
     *           + P(1)*P(19)*(1-dabs(P(2))-dabs(P(5))-dabs(P(8))
     *          - dabs(P(11))) *   reR / dX
      ELSEIF (ifc.eq.158) THEN ! 2 * Lorentzian (normalized)
         pi = twopi/2
         DO i = 1, n
            Y(i) = P(1)
     *  * ( dquot0(dabs(P(2))*dabs(P(4))/pi,(X(i)-P(3))**2+P(4)**2)
     *  + dquot0((1-dabs(P(2)))*dabs(P(6))/pi,(X(i)-P(5))**2+P(6)**2))
     *  + P(7)
         ENDDO
      ELSEIF (ifc.eq.159) THEN !
         pi = twopi/2
         FNU = 1
         NUM = 10
         SCALE = 'U'
         DO i = 1, n
C            Z = CMPLX (dabs(P(2)),0.0) !Artem: Replace with dbesi function from slatec/src/dbesi.f
            bZ =  dabs(P(2))
            CALL dbesi(bZ, FNU, 1, NUM, baY, NZ)
C            Y(i) = S17DEF (FNU,Z,NUM,SCALE,CY,NZ,IFAIL) !Artem: Replace with dbesi function from slatec/src/dbesi.f:  S17DEF (FNU,Z,NUM,SCALE,CY,NZ,IFAIL)
            Y(i) = dquot0(dabs(P(4))/pi,(X(i)-P(3))**2 + P(4)**2)*
     *            baY(1)**2  ! Y(i) = dquot0(dabs(P(4))/pi,(X(i)-P(3))**2 + P(4)**2)*dble(CY(1))**2
            Do j = 2, NUM
               Y(i) = Y(i) +
     *    dquot0(j*(j+1)*dabs(P(6))/pi, ((X(i)-P(5))**2+(j*(j+1)*
     *    dabs(P(6)))**2)) * (2*j+1) * baY(j)**2 ! Y(i) = Y(i) + dquot0(j*(j+1)*dabs(P(6))/pi, ((X(i)-P(5))**2+(j*(j+1)* dabs(P(6)))**2)) * (2*j+1) * dble(CY(j))**2
            ENDDO
            Y(i) = Y(i) * P(1) + P(7)
         ENDDO
      ELSEIF (ifc.eq.160) THEN !
         pi = twopi/2
         FNU = 1
         NUM = 12
         DO i = 1, NUM
            FNU = 2*dabs(P(2))*dsin(pi*i/NUM)
            JJJ0(i) = DBESI0(FNU) !Artem: Replace with dbesi0 function from slatec/fnlib/dbesi0.f:  S17AEF (FNU,IFAIL)
         ENDDO
         DO i = 1, NUM
            AAA(i) = 0
            DO j = 1, NUM
               AAA(i) = AAA(i) + JJJ0(j) * dcos(twopi*i*j/NUM)
            ENDDO
            AAA(i) = AAA(i) / NUM
         ENDDO
         DO i = 1, n
            Y(i) = dquot0(dabs(P(4))/pi,(X(i)-P(3))**2+P(4)**2)*
     *            AAA(NUM)
            Do j = 1, NUM-1
               FNU = 1.0 / (4*dabs(P(5))*(dsin(pi*j/NUM))**2)
!               print *, 'FNU: ', FNU
               Y(i) = Y(i) + AAA(j)/pi * dquot0(FNU,1+(X(i)*FNU)**2)
            ENDDO
            Y(i) = Y(i) * P(1) + P(6)
         ENDDO
      ELSEIF (ifc.eq.161) THEN !
         pi = twopi/2
         FNU = 1
         NUM1 = 12
         NUM2 = 12
         Q10RATIO = 13/9
         DO i = 1, NUM1
            FNU = 2*dabs(P(3))*dsin(pi*i/NUM1)
            JJJ1(i) = dbesi0 (FNU) !Artem: Replace with dbesi0 function from slatec/fnlib/dbesi0.f:  S17AEF (FNU,IFAIL)
         ENDDO
         DO i = 1, NUM2
            FNU = 2*dabs(P(5))*dsin(pi*i/NUM2)
            JJJ2(i) = dbesi0 (FNU) !Artem: Replace with dbesi0 function from slatec/fnlib/dbesi0.f:  S17AEF (FNU,IFAIL)
         ENDDO
         DO i = 1, NUM1
            AAA1(i) = 0
            DO j = 1, NUM1
               AAA1(i) = AAA1(i) + JJJ1(j) * dcos(twopi*i*j/NUM1)
            ENDDO
            AAA1(i) = AAA1(i) / NUM1
!            print *, 'B1n: ', i, AAA1(i)
         ENDDO
         DO i = 1, NUM2
            AAA2(i) = 0
            DO j = 1, NUM2
               AAA2(i) = AAA2(i) + JJJ2(j) * dcos(twopi*i*j/NUM2)
            ENDDO
            AAA2(i) = AAA2(i) / NUM2
!            print *, 'B2n: ', i, AAA2(i)
         ENDDO
         DO i = 1, n
            Y(i) = dquot0(dabs(P(2)),pi*((X(i))**2+P(2)**2))*
     *            AAA1(NUM1) * AAA2(NUM2)
            Do j = 1, NUM1-1
            ITAU1 = dquot0((dsin(pi*j/NUM1))**2,(dabs(P(4))*
     *      (dsin(pi/NUM1))**2))
!               print *, 'ITAU1: ', ITAU1
               Do k=1, NUM2-1
            ITAU2 = dquot0((dsin(pi*k/NUM2))**2, (dabs(P(6))*
     *              (dsin(pi/NUM2))**2))
                  FNU = ITAU1+ITAU2+dabs(P(2))
                  Y(i) = Y(i) + (AAA1(j) * AAA2(k) *
     *                 dquot0(FNU,FNU**2+(X(i))**2)) / 2 / pi**2
               ENDDO
            ENDDO
            FNU = 0
            Do j = 1, NUM1-1
               ITAU1 = (dsin(pi*j/NUM1))**2 / (dabs(P(4))*
     *                 (dsin(pi/NUM1))**2)
               FNU = FNU + AAA1(j) * dquot0(ITAU1+dabs(P(2)),(ITAU1+
     *               dabs(P(2)))**2+(X(i))**2)
            ENDDO
            Y(i) = Y(i) + AAA2(NUM2)/dsqrt(2*pi**3) * FNU
            FNU = 0
            Do j = 1, NUM2-1
               ITAU2 = (dsin(pi*j/NUM2))**2 / (dabs(P(6))*
     *                 (dsin(pi/NUM2))**2)
               FNU = FNU + AAA2(j) * dquot0(ITAU2+dabs(P(2)),(ITAU2+
     *         dabs(P(2)))**2+(X(i))**2)
            ENDDO
            Y(i) = Y(i) + AAA1(NUM1)/dsqrt(2*pi**3) * FNU
            Y(i) = Y(i) / Q10RATIO
            FNU = 0
            Do j = 1, NUM1-1
               ITAU1 = (dsin(pi*j/NUM1))**2 / (dabs(P(4))*
     *                 (dsin(pi/NUM1))**2)
               FNU = FNU + AAA1(j) * dquot0(ITAU1+dabs(P(2)),(ITAU1+
     *         dabs(P(2)))**2+(X(i))**2)
            ENDDO
            Y(i) = Y(i) + Q10RATIO * (AAA1(NUM1)/dsqrt(2*pi**3)+
     *      dquot0(dabs(P(2)),(P(2))**2+(X(i))**2)+1/2/pi**2*FNU)
            Y(i) = Y(i) * twopi * P(1) + P(7)
         ENDDO
      ELSEIF (ifc.eq.170) THEN ! sph. bessel ftn. of first kind J_0(x)
         IFAIL = 0
         DO i = 1,n
            FNU =  X(i)*P(2)
            Y(i) = P(1)*dbesi0(FNU) !Artem: Replace with dbesi0 function from slatec/fnlib/dbesi0.f:  S17AEF (FNU,IFAIL)
            ENDDO
      ELSEIF (ifc.eq.171) THEN ! 3-fold jump a la Zorn
         IFAIL = 0
         DO i = 1,n
            FNU =  X(i)*P(2)
            Y(i) = P(3) * (((1. - P(1)) * 0.333 * (1. +
     *                  2. * dbesi0(FNU))) + P(1)) !Artem: Replace with dbesi0 function from slatec/fnlib/dbesi0.f:  S17AEF (FNU,IFAIL)
            ENDDO
      ELSEIF (ifc.eq.172) THEN ! sph. bessel ftn. of first kind J_0(x)*J_0(x) ! yf 12.08.2010
         IFAIL = 0
         DO i = 1,n
            FNU =  X(i)*P(2)
            Y(i) =
     * P(1)*((1-P(3))*dbesi0(FNU)*dbesi0(FNU)+P(3)) !Artem: Replace with dbesi0 function from slatec/fnlib/dbesi0.f:  S17AEF (FNU,IFAIL)
            ENDDO
      ELSEIF (ifc.eq.180) THEN ! special small angle scattering
         IFAIL = 0
         DO i = 1,n
            FNU=6.0
            FNU=dquot0(X(i)*P(3),dsqrt(FNU))
            RFNU = ERF(FNU) !Artem: Replace with a standard ERF function: S15AEF(FNU,IFAIL)
            Y(i) = (P(1) * dexp(- (P(3)**2 / 3.0)* X(i)**2 )  + P(2) *
     *           (dquot0(RFNU**3, X(i)))**P(4))*
     *           dquot0((X(i)*P(6))**3,((X(i)*P(6))**3 + 3.0*P(5)*
     *            (dsin(X(i)*P(6)) - X(i)*P(6)*dcos(X(i)*P(6)))))+
     *           P(7)*X(i)**P(8)
            ENDDO
      ELSEIF (ifc.eq.181) THEN ! modified special small angle
         IFAIL = 0
         IFAIL1 = 0
         DO i = 1,n
            FNU=6.0
            FNU1=dquot0(X(i)*P(9),dsqrt(FNU))
            FNU=dquot0(X(i)*P(3),dsqrt(FNU))
            RFNU = ERF(FNU) !Artem: Replace with a standard ERF function: S15AEF(FNU,IFAIL)
            RFNU1 = ERF (FNU) !Artem: Replace with a standard ERF function: S15AEF(FNU,IFAIL1)
            Y(i) = (P(1) * dexp(- (P(3)**2 / 3.0)* X(i)**2 )  + P(2) *
     *           (dquot0(RFNU**3, X(i)))**P(4))*
     *           dquot0((X(i)*P(6))**3,((X(i)*P(6))**3 + 3.0*P(5)*
     *            (dsin(X(i)*P(6)) - X(i)*P(6)*dcos(X(i)*P(6)))))+
     *           (P(7) * dexp(- (P(9)**2 / 3.0)* X(i)**2 )  + P(8) *
     *           (dquot0(RFNU1**3, X(i)))**P(10))*
     *           dquot0((X(i)*P(12))**3,((X(i)*P(12))**3 + 3.0*P(11)*
     *            (dsin(X(i)*P(12)) - X(i)*P(12)*dcos(X(i)*P(12)))))
            ENDDO
      ELSEIF (ifc.eq.190) THEN ! Brioullin plus Rayleigh
               pi = twopi / 2
                  DO i = 1,n
               Y(i) = P(2) * (dquot0(dabs(P(3))/pi,(X(i)-P(4))**2 +
     *               P(3)**2) + dquot0(dabs(P(3))/pi,(X(i)+P(4))**2 +
     *               P(3)**2)) + P(5)
            IF (X(i).lt. P(6) .and. X(i).gt. -P(6)) THEN
               Y(i) = Y(i) + P(1)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.191) THEN ! Brioullin plus Rayleigh
               pi = twopi / 2
                  DO i = 1,n
               Y(i) = P(1) * dquot0(dabs(P(3))/pi,(X(i)-P(4))**2 +
     *                P(3)**2) + P(2) * (dquot0(dabs(P(5))/pi,(X(i)-
     *                P(6))**2 + P(5)**2) + dquot0(dabs(P(5))/pi,(X(i)+
     *                P(6))**2 + P(5)**2)) + P(7) + P(8)*X(i)
            ENDDO
      ELSEIF (ifc.eq.192) THEN ! 190 amplitude lorentz decoupled
               pi = twopi / 2
                  DO i = 1,n
               Y(i) = P(2) * (dquot0(P(3)**2,(X(i)-P(4))**2 +
     *               P(3)**2) + dquot0(P(3)**2,(X(i)+P(4))**2 +
     *               P(3)**2)) + P(5)
            IF (X(i).lt. P(6) .and. X(i).gt. -P(6)) THEN
               Y(i) = Y(i) + P(1)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.193) THEN ! 191 amplitude lorentz decoupled
               pi = twopi / 2
                  DO i = 1,n
               Y(i) = P(1) * dquot0(P(3)**2,(X(i)-P(4))**2 +
     *                P(3)**2) + P(2) * (dquot0(P(5)**2,(X(i)-
     *                P(6))**2 + P(5)**2) + dquot0(P(5)**2,(X(i)+
     *                P(6))**2 + P(5)**2)) + P(7) + P(8)*X(i)
            ENDDO

      ELSEIF (ifc.eq.195) THEN ! 6 special ratioed Gauss NMR
         DO i = 1, n
            p1 = P(1)*P(2)
            p2 = P(1)*P(5)
            p3 = P(1)*P(8)
            p4 = P(1)*P(11)
            p5 = P(1)*P(14)
            p6 = P(1)*P(17)   ! P(1) = aTOT (rest percentages 0:1)
            Y(i) =    p1 * dexp1 ( -dquot0(X(i)-P( 4), P( 3))**2 )
     *             +  P2 * dexp1 ( -dquot0(X(i)-P( 7), P( 6))**2 )
     *             +  P3 * dexp1 ( -dquot0(X(i)-P(10), P( 9))**2 )
     *             +  P4 * dexp1 ( -dquot0(X(i)-P(13), P(12))**2 )
     *             +  P5 * dexp1 ( -dquot0(X(i)-P(16), P(15))**2 )
     *             +  P6 * dexp1 ( -dquot0(X(i)-P(19), P(18))**2 )
            ENDDO
       ELSEIF (ifc.eq.201) THEN ! Re Lorentz
         DO i = 1, n
            Y(i) = P(1) * dquot0(-P(2)*X(i), P(2)**2 + (X(i)-P(3))**2)
            ENDDO
      ELSEIF (ifc.eq.202) THEN ! Im Lorentz
         DO i = 1, n
            Y(i) = P(1) * dquot0(P(2)**2, P(2)**2 + (X(i)-P(3))**2)
            ENDDO
      ELSEIF (ifc.eq.211) THEN ! Brillouin-S(qw) of Lorentz (30nov95)
         ! liquid-like limiting case w*tau<<1 (Montrose Solovyev Litovitz 68)
         DO i = 1, n
            xp   = X(i) + P(3)
            xm   = X(i) - P(3)
            Y(i) = P(1) * ( dquot0 ( P(2), xp**2 + P(2)**2 )
     *             + dquot0 ( P(2), xm**2 + P(2)**2 )
     *             + dquot0 ( P(2)*xp, P(3)*(xp**2 + P(2)**2) )
     *             - dquot0 ( P(2)*xm, P(3)*(xm**2 + P(2)**2) ) )
            ENDDO
      ELSEIF (ifc.eq.221) THEN ! Re Havriliak Negami. See Alvaraez..91, eqn.15
            ! 'A0;Ai;tau;a;g'
         DO i = 1, n
            CALL HavNeg (i, X(i), P(3), P(1), P(2),
     *                   P(4), P(5), Y(i), dummy)
            ENDDO
      ELSEIF (ifc.eq.222) THEN ! Im Havriliak Negami. See Alvaraez..91, eqn.16
         DO i = 1, n
            CALL HavNeg (i, X(i), P(2), 0.d0, P(1),
     *                   P(3), P(4), dummy, Y(i))
            ENDDO
      ELSEIF (ifc.eq.223) THEN ! S(w) Havriliak Negami = Im / w + bg
         DO i = 1, n
            CALL HavNeg (i, X(i), P(2), 0.d0, P(1),
     *                   P(3), P(4), dummy, HNI)
            Y(i) = dquot0(HNI,X(i)) + P(5)
            ENDDO
      ELSEIF (ifc.eq.224 .or. ifc.eq.225) THEN
         ! Brillouin-S(qw) or X''(qw) of Havriliak-Negami (6aug93)
            ! 'A;q;c0;ci;tau;a;g;G0;G1;G2'
            ! units : q in nm-1, x in GHz, c0 in m/sec, tau in GHz-1.
         A0  = (P(3)*P(2)/twopi)**2      ! (q*(speed of sound at w=0))**2
         Ai  = (P(4)*P(2)/twopi)**2      ! dito at w=infty
         DO i = 1, n
            wi  = dabs(X(i))
            ahf = P(8) + P(9)*wi ! friction
            CALL HavNeg (i, wi, P(5), A0, Ai, P(6), P(7), HNR, HNI)
            HNI = HNI - ahf
            HND = HNR - wi**2 ! - w^2
            IF (ifc.eq.224) THEN
               Y(i)= dquot0 (-P(1)*HNI, wi*(HND**2 + HNI**2))
            ELSE
               Y(i)= dquot0 (-P(1)*HNI, (HND**2 + HNI**2))
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.226) THEN ! Novosibirsk + Havriliak Negami
            ! '&j;tau;a;g;d;T;T0;G0;G1;G2;debug'
         A  = 4 / twopi * dquot0 (P(5), P(6))
         ! check par-file :
         DO ii = 1, nFPF
            IF (XFPF(ii).le.0) THEN
               Print *, 'CuFuVal/ par-file contains x<= 0'
               CALL rSet (Y, 1, n, 1, -31.d0)
               RETURN
               ENDIF
            ENDDO
         DO i = 1, n
            CALL HavNeg (i, X(i), P(1), 0.d0, 1.d0,
     *                   P(2), P(3), HNR, HNI)
            yInt = 0.
            DO ii = 1, nFPF-1 ! Integration over Eigenmodes Omega
               Oma  = ( XFPF(ii+1) + XFPF(ii) ) / 2      ! Omega
c              IF (Oma.gt.P(8)) THEN
               OmaD =   XFPF(ii+1) - XFPF(ii)            ! Integration step
               fri  =  P(7) * (OmaD*X(I)/OMA)**2         ! aux. friction
c              d2 = dmax1 (Oma**2 - P(4)**2, 0.d0)
               d2   = P(4) * Oma**2
               HNIY = -d2 * HNI
               HNRY = Oma**2 + d2 * (HNR-1)
               xim  = dquot0 (fri+HNIY,
     *                        (X(i)**2-HNRY)**2 + (HNIY+fri)**2)
               yInt = yInt + OmaD * Oma * YFPF(ii) * xim
c              ENDIF
               ENDDO ! Omega
            Y(i) = A * yInt !dquot0 (yInt, X(i))
            ENDDO ! w
      ELSEIF (ifc.eq.227) THEN ! Convolute Havriliak-Negami (Neufassung 23jan95)
            ! '&j:tau;a;g;d;T;T0;pow;G0'
               ! pow : curve is g(w) / w^pow
         A  = 4 / twopi * dquot0 (P(5), P(6))
         ! check par-file :
         DO ii = 1, nFPF
            IF (XFPF(ii).le.0) THEN
               Print *, 'CuFuVal/ par-file contains x<= 0'
               CALL rSet (Y, 1, n, 1, -31.d0)
               RETURN
               ENDIF
            ENDDO
         DO i = 1, n
            CALL HavNeg (i, X(i), P(1), 0.d0, 1.d0,
     *                   P(2), P(3), HNR, HNI)
            yInt = 0.
            DO ii = 1, nFPF-1 ! Integration over Eigenmodes Omega
               Oma  = ( XFPF(ii+1) + XFPF(ii) ) / 2      ! Omega
               OmaD =   XFPF(ii+1) - XFPF(ii)            ! Integration step
               fri  = P(8) * 2*OmaD*X(i)         ! aux. friction
               d2   = P(4) * Oma**2
               HNIY = -d2 * HNI
               HNRY = Oma**2 + d2 * (HNR-1)
               Bruch= dquot0 (fri+HNIY,
     *                        (X(i)**2-HNRY)**2 + (HNIY+fri)**2)
               yInt = yInt + OmaD * YFPF(ii) * dpow0(Oma,P(7))* Bruch
               ENDDO ! Omega
            Y(i) = A * yInt * dpow0 (X(i), 1-P(7))
            ENDDO ! w
      ELSEIF (ifc.eq.228) THEN  ! hybrid Gauss KWW bg
         tau = dquot0 (P(4) * P(3), dgamma1(dquot0(1.d0,P(4))))
         CALL ExtTabXZ ('kww_sqw', X, 1, n,
     *                  P(4), 0.d0, tau/hbar, P(1), P(2), Y, Fehler)
         DO i = 1,n
              Y(i) = Y(i) + P(5) * dexp1(- dquot0(X(i)-P(6),P(7))**2)
            ENDDO
      ELSEIF (ifc.eq.229) THEN  ! function 234 + sloping bg
         tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         CALL ExtTabXZ ('kww_sqw', X, 1, n,
     *                  P(3), 0.d0, tau/hbar, P(4), P(1), Y, Fehler)
         DO i = 1,n
            IF (X(i).lt. 0.001 .and. X(i).gt. -0.001) THEN
               Y(i) = Y(i) + P(5)
               ENDIF
               Y(i) = Y(i) + P(6) * X(i)
            ENDDO
      ELSEIF (ifc.eq.230) THEN  ! function 233 + sloping bg
            tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         CALL ExtTabXZ ('kww_sqw', X, 1, n,
     *                  P(3), 0.d0, tau/hbar, P(4), P(1), Y, Fehler)
         DO i = 1, n
            Y(i) = Y(i) + P(5) * X(i)
            ENDDO
      ELSEIF (ifc.eq.231) THEN  ! Re KWW in w
         tau = dquot0 (P(4) * P(3), dgamma1(dquot0(1.d0,P(4))))
         CALL ExtTabXZ ('kww_real', X, 1, n,
     *                  P(4), 0.d0, tau, 0.d0, P(2), Y, Fehler)
         DO i = 1, n
            Y(i) = P(1) + Y(i)
            ENDDO
      ELSEIF (ifc.eq.232) THEN  ! Im KWW in w
         tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         CALL ExtTabXZ ('kww_imag', X, 1, n,
     *                  P(3), 0.d0, tau, 0.d0, P(1), Y, Fehler)
      ELSEIF (ifc.eq.233) THEN  ! S_KWW in w
         tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         CALL ExtTabXZ ('kww_sqw', X, 1, n,
     *                  P(3), 0.d0, tau/hbar, P(4), P(1), Y, Fehler)
      ELSEIF (ifc.eq.234) THEN  ! S_KWW in w + elast
         tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         CALL ExtTabXZ ('kww_sqw', X, 1, n,
     *                  P(3), 0.d0, tau/hbar, P(4), P(1), Y, Fehler)
         DO i = 1,n
            IF (X(i).lt. 0.001 .and. X(i).gt. -0.001) THEN
               Y(i) = Y(i) + P(5)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.235) THEN  ! KWW: eps'' vs. eps'
         ymul = P(2) * (P(4)-P(3))
         xmul = dquot0 (1.d0, (P(4)-P(3)))
         CALL ExtTabXZ ('kww_im_re', X, 1, n,
     *                  P(5), P(3), xmul, P(1), ymul, Y, Fehler)
      ELSEIF (ifc.eq.236) THEN  ! Re KWW in w
         tau = dquot0 (P(4) * P(3), dgamma1(dquot0(1.d0,P(4))))
         DO i = 1, n
            Y(i) = P(1) + P(2) * ! 2 / twopi *
     *       X(i) * tau * RauschKohl (.false., tau*X(i), P(4))
            ENDDO
      ELSEIF (ifc.eq.237) THEN  ! Im KWW in w
         tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         DO i = 1, n
            Y(i) = P(1) * ! 2 / twopi *
     *       X(i) * tau * RauschKohl (.true., tau*X(i), P(3))
            ENDDO
      ELSEIF (ifc.eq.238) THEN  ! S_KWW in w
         tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         Print *, 'CHECK WHETHER TO ADD 2 / twopi'
         DO i = 1, n
            Y(i) = P(1) * tau *
     *       RauschKohl (.true., tau*X(i), P(3)) + P(4)
            ENDDO
      ELSEIF (ifc.eq.239) THEN  ! S_KWW in w
         DO i = 1, n
            Y(i) = RauschFunc (0, X(i), P(1), P(2))
            ENDDO
      ELSEIF (ifc.eq.240) THEN !delta + bg
         DO i = 1,n
            IF (X(i).lt. 0.001 .and. X(i).gt. -0.001) THEN
               Y(i) = Y(i) + P(1)
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.241) THEN ! Convolute Kohlrausch
            ! '&j:tau;beta;Delta;T;T0;pow;G0'
               ! pow : curve is g(w) / w^pow
         A  = 4 / twopi * dquot0 (P(4), P(5))
         ! check par-file :
         DO ii = 1, nFPF
            IF (XFPF(ii).le.0) THEN
               Print *, 'CuFuVal/ par-file contains x<= 0'
               CALL rSet (Y, 1, n, 1, -31.d0)
               RETURN
               ENDIF
            ENDDO
         IF (n.gt.MC) CALL Absturz ('CuFuVal', 'Y1, Y2 <= MC')
         CALL ExtTabXZ ('kww_imag', X, 1, n, P(2),
     *                  P(1), P(3), 0.d0, 1.d0, Y2, Fehler)
         CALL ExtTabXZ ('kww_real', X, 1, n, P(2),
     *                  P(1), P(3), 0.d0, 1.d0, Y1, Fehler)
         DO i = 1, n
            yInt = 0.
            DO ii = 1, nFPF-1 ! Integration over Eigenmodes Omega
               Oma  = ( XFPF(ii+1) + XFPF(ii) ) / 2      ! Omega
               OmaD =   XFPF(ii+1) - XFPF(ii)            ! Integration step
               fri  = P(7) * 2*OmaD*X(i)         ! aux. friction
               HNIY = Oma**2 * Y2(i)
               HNRY = Oma**2 * ( 1 - P(3) +  Y1(i) ) ! Y1(0..infty)=P(3)..0
               Bruch= dquot0 (fri+HNIY,
     *                        (X(i)**2-HNRY)**2 + (HNIY+fri)**2)
               yInt = yInt + OmaD * YFPF(ii) * dpow0(Oma,P(6))* Bruch
               ENDDO ! Omega
            Y(i) = A * yInt * dpow0 (X(i), 1-P(6))
            ENDDO ! w
      ELSEIF (ifc.eq.243) THEN  ! S_KWW in w
         tau = dquot0 (P(3) * P(2), dgamma1(dquot0(1.d0,P(3))))
         CALL ExtTabXZ ('kww_sqw', X, 1, n,
     * P(3), 0.d0, tau/hbar, dabs( P(4) ), P(1), Y, Fehler)
      ELSEIF (ifc.eq.251) THEN
            ! Brillouin of Cole Davidson [M.Frank 2000]
            ! 'A;go;Temp;tau;beta;deltaqudr;w1;w2'
         DO i = 1, n
            w=X(i)*twopi
            g=P(2)*twopi
            wo=twopi*(P(7)-P(8)*P(3))
            CALL ColDav ( w, P(4), P(5), P(6), CDR, CDI)
            ob = g+CDI
            una = (w**2-wo**2+w*CDR)**2
            unb = (w*g+w*CDI)**2
            Y(i)= P(1)*dquot0(ob,(una+unb))
            ENDDO
      ELSEIF (ifc.eq.252) THEN
            ! Brillouin of Hybrid model [M.Frank 2000]
            ! 'A;go;Temp;tau;beta;deltaqudr;w1;w2;B;a'
         DO i = 1, n
            w=X(i)*twopi
            g=P(2)*twopi
            wo=twopi*(P(7)-P(8)*P(3))
            bet=P(9)*P(6)*dpow0(P(4)**(-2)+w**2,(P(10)-1)/2)
            agm=(P(10)-1)*datan(-w*P(4))
            CALL ColDav ( w, P(4), P(5), P(6), CDR, CDI)
            agi = CDI+bet*dcos(agm)
            REA = CDR-bet*dsin(agm)
            ob = g+agi
            una = (w**2-wo**2+w*REA)**2
            unb = (w*g+w*agi)**2
            Y(i)= P(1)*dquot0(ob,(una+unb))
            ENDDO
      ELSEIF (ifc.eq.261) THEN  ! Percus-Yevick nC(k) [Boon Yip 2.3.15]
         DO i = 1, n
            Y(i) = PercYevC (P(1)*X(i), P(2))
            ENDDO
      ELSEIF (ifc.eq.262) THEN  ! Percus-Yevick S(k) [Boon Yip 2.3.14,15]
         DO i = 1, n
            Y(i) = PercYevS (P(1)*X(i), P(2))
            ENDDO
      ELSEIF (ifc.eq.263) THEN  ! Percus-Yevick S(k) with diam=diam(q0,eta)
         diam = PercYevDiam (P(1), P(2))
         DO i = 1, n
            Y(i) = PercYevS (diam*X(i), P(2))
            ENDDO
      ELSEIF (ifc.eq.264) THEN  ! Percus-Yevick MC vertex
         diam = P(1) ! PercYevDiam (P(1), P(2))
         IF (diam.le.0) THEN
            CALL rSet (Y, 1, n, 1, -164.1d0)
            RETURN
            ENDIF
         dens = 12 * P(2) / twopi / diam**3
         DO i = 1, n
            wq = X(i)*diam
            wk = P(3)*diam
            wp = P(4)*diam
            Sq = PercYevS (wq, P(2))
            Sk = PercYevS (wk, P(2))
            Sp = PercYevS (wp, P(2))
            IF (Sq.le.0 .or. Sk.le.0 .or. Sp.le.0 .or. wq.le.0) THEN
               Y(i) = 0
            ELSE
               Y(i) = ( 1.d0 / 8 / twopi**2 / dens) *
     *  Sq*Sk*Sp * (wk*wp/wq**5) *
     *  ((wq**2+wk**2-wp**2)*(Sk-1)/Sk
     * + (wq**2+wp**2-wk**2)*(Sp-1)/Sp)**2
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.265) THEN  ! PYS(k) with diam(q0,eta) and distribution
         diam = PercYevDiam (P(1), P(2))
         DO i = 1, n
            Y(i) = 0
            ENDDO
         DO ii = 1, 99
            ang = twopi * ii / 200
            rel = 1 + P(3) * dcos(ang)
            wgt = dsin(ang) * twopi / 4 / 99
            DO i = 1, n
               Y(i) = Y(i) + wgt * PercYevS (diam*rel*X(i), P(2))
               ENDDO
            ENDDO
      ELSEIF (ifc.eq.272) THEN ! gleich wieder loeschen 18jan00
         eps = dmin1 (1.d0, dmax1 (P(1),0.d0))
         rho = 1 - eps
         amp = 2 * dquot0 (1.d0, twopi/8 * (1 - rho**2) * eps)
         DO i = 1, n
            xarg = dabs(X(i))
            IF (xarg.gt.1) THEN
               Y(i) = 0
            ELSEIF (xarg.gt.rho) THEN
               Y(i) = amp * (1-xarg**2)
            ELSE
               Y(i) = amp * (dsqrt(1-xarg**2)-dsqrt(rho**2-xarg**2))**2
               ENDIF
            ENDDO
      ELSEIF (ifc.eq.273) THEN ! 4 Voigt + bg
         DO i = 1, n
            ww = X(i)
            Y(i) = P(12)
     *           + P(7) * Voigt((P(1)+P(3)*(P(2)+0.5)),
     *           (dabs(P(2)+0.5)*P(4)), P(11), ww)
     *           + P(8) * Voigt((P(1)+P(3)*(P(2)-0.5)),
     *           (dabs(P(2)-0.5)*P(4)), P(11), ww)
     *           + P(9) * Voigt((P(1)+P(5)*(P(2)+0.5)),
     *           (dabs(P(2)+0.5)*P(6)), P(11), ww)
     *           + P(10) * Voigt(P(1)+P(5)*(P(2)-0.5),
     *           (dabs(P(2)-0.5)*P(6)), P(11), ww)
             ENDDO
      ELSEIF (ifc.eq.274) THEN ! Voigt + bg
         DO i = 1, n
            ww = X(i)
            Y(i) = P(5) + P(4)*Voigt(P(1),P(2),P(3),ww)
             ENDDO
      ELSEIF (ifc.eq.275) THEN ! damped oscillation
         DO i = 1, n
            Y(i) = P(1) + dabs(P(2)) * dcos(twopi * P(3) * (
     *             X(i) - P(4))) * dexp (-X(i) * P(5))
             ENDDO
      ELSEIF (ifc.eq.276) THEN ! t dep. damped osci.
         DO i = 1, n
            Y(i) = P(1) + dabs(P(2)) *
     *       dcos(twopi * (P(3) + P(4) * X(i))
     *       * (X(i) - P(5))) * dexp (-X(i) * (P(6) + X(i) * P(7)))
             ENDDO
      ELSEIF (ifc.eq.277) THEN ! damped oscillation+additional oscillation
         DO i = 1, n
            Y(i) = P(1) + dabs(P(2)) * dcos(twopi * P(3) * (
     * X(i) - P(4))) * dexp (-X(i) * P(5)) + dabs (P(6)) *
     * dcos( twopi * P(7) * (X(i) - P(8)))
             ENDDO
      ELSEIF (ifc.eq.278) THEN ! damped oscillation+ additional oscillations coupled
         DO i = 1, n
            Y(i) = P(1) + dabs(P(2)) * dcos(twopi * P(3) * (
     * X(i) - P(4))) * dexp (-X(i) * P(5)) + dabs (P(6)) *
     * dcos( twopi * P(3) * (X(i) - P(4)))
             ENDDO
      ELSEIF (ifc.eq.279) THEN ! damped oscillation+ additional oscillations coupled
         cbrt = twopi/twopi/3
         DO i = 1, n
            Y(i) = (-P(1)) * ((X(i)-P(2)) ** (-cbrt))
             ENDDO
      ELSEIF (ifc.eq.280 .or. ifc.eq.281) THEN
         nbunches = idnint (P(6))
         IF     (nbunches.lt.1) THEN
            Fehler = 'bad parameter: # bunches < 1'
            RETURN
         ELSEIF (nbunches.gt.20) THEN
            Fehler = 'bad parameter: too many bunches required'
            RETURN
         ELSEIF (nbunches.gt.1 .and. P(7).le.0) THEN
            Fehler = 'bad parameter: time between bunches nonpositive'
            RETURN
            ENDIF
         tau0 =  P(4)
         IF (tau0.le.0) THEN ! should be 141.11
            Fehler = 'natural lifetime tau_0 <= 0'
            RETURN
            ENDIF
         IF     (ifc.eq.280) THEN ! one Line (Kagan-Afa..-Kohn)
            DO i = 1, n
               ifail = 0
               Y(i) = P(3) ! background
               DO ibu = 1, nbunches ! serie of bunches
                  time = X(i) + (ibu-1) * P(7)
                  Y(i) = Y(i) + P(1) * P(2)**2 / 16 / tau0**2 * (
     *                   dexp1( -time/tau0 - 2*P(5)*time )
     *                 * 4 * dquot0(tau0,P(2)*time)
     *                 * (dbesi1(dsqrt(P(2)*time/tau0)))**2 !Artem: Replace with dbesi1 function from slatec/fnlib/dbesi1.f:  S17AFF(dsqrt(P(2)*time/tau0), ifail)
     *                 + intq(time.eq.0) )
                  ENDDO
               ENDDO
         ELSEIF (ifc.eq.281) THEN ! two widely separated lines (A5,13)
            DO i = 1, n
               ifail = 0
               Y(i) = P(3) ! background
               sepp = P(8) / (2*hbar*1d3)
               DO ibu = 1, nbunches ! serie of bunches
                  time = X(i) + (ibu-1) * P(7)
                  Y(i) = Y(i) + P(1) * P(2)**2 / 16 / tau0**2 * (
     *                   dexp1( -time/tau0 - 2*P(5)*time )
     *                 * (dcos(sepp*time) )**2
     *                 * 4 * dquot0(2*tau0,P(2)*time)
     *                 * (dbesi1(dsqrt(P(2)*time/(2*tau0))))**2 !Artem: Replace with dbesi1 function from slatec/fnlib/dbesi1.f:  S17AFF(dsqrt(P(2)*time/(2*tau0)), ifail)
     *                 + intq(time.eq.0) )
                  ENDDO
               ENDDO
         ELSE
            Fehler = 'invalid subcase'
            RETURN
            ENDIF
      ELSEIF (ifc.eq.282) THEN  ! Fe3Si N"aherung
         DO i = 1, n
               ifail = 0
               tlife = 141.11
               tbunch = 176.0
               X(i) = dquot0 (X(i),tlife)
               Y(i) =  P(1)
     *              * (dexp1( - X(i)) + P(4) * dexp1( - P(2) * X(i)))
     *              * dquot0 (P(3),X(i))
     *              * (dbesi1 (dsqrt(P(3) * X(i))))**2 !Artem: Replace with dbesi1 function from slatec/fnlib/dbesi1.f:  S17AFF (dsqrt(P(3) * X(i)), ifail)
     *              + P(5)
               X(i) = X(i) + dquot0 (tbunch,tlife)
               Y(i) =  Y(i) + P(1)
     *              * (dexp1( - X(i)) + P(4) * dexp1( - P(2) * X(i)))
     *              * dquot0 (P(3),X(i))
     *              * (dbesi1 (dsqrt(P(3) * X(i))))**2 !Artem: Replace with dbesi1 function from slatec/fnlib/dbesi1.f:  S17AFF (dsqrt(P(3) * X(i)), ifail)
               X(i) = X(i) + dquot0 (tbunch,tlife)
               Y(i) =  Y(i) + P(1)
     *              * (dexp1( - X(i)) + P(4) * dexp1( - P(2) * X(i)))
     *              * dquot0 (P(3),X(i))
     *              * (dbesi1 (dsqrt(P(3) * X(i))))**2 !Artem: Replace with dbesi1 function from slatec/fnlib/dbesi1.f:  S17AFF (dsqrt(P(3) * X(i)), ifail)
               X(i) =  X(i) * tlife - 2 * tbunch
            ENDDO
      ELSEIF (ifc.eq.283) THEN
         cbrt = twopi/twopi/3
         pi = twopi/2
         DO i = 1, n
            Y(i) = P(1) + P(2) * (-P(3)*( (X(i)-P(4))**(-cbrt) )) +
     *             P(5) * dsin ( P(6)*(X(i)-P(7)) - P(8) )
            ENDDO
      ELSEIF (ifc.eq.284) THEN
         DO i = 1, n
            Y(i) = P(1) + P(2) * X(i) + P(3) *
     * dsin ( P(4)*(-(P(5)/X(i))**3 + P(6) - P(7)) - P(8))
            ENDDO
      ELSEIF (ifc.eq.291) THEN  ! e^-2W(q,t) from Debye model
         alphusw = P(2) * P(1)**2 ! amp*q^2
         DO i = 1, n
            wt = P(3)*X(i)
            IF (dabs(wt).lt.1.d-20) THEN
               siwt = 1
            ELSE
               siwt = dsin(wt) / wt
               ENDIF
            Y(i) = dexp(-alphusw) * (dexp (alphusw*siwt) - 1)
            ENDDO
      ELSEIF (ifc.eq.292) THEN  ! e^-2W(q,t) from Gauss-Debye model C2,152
         alphusw = P(2) * P(1)**2 ! amp*q^2
         DO i = 1, n
            wt = P(3)*X(i) ! omega in 1/s, nicht in Hz !
            Y(i) = dexp(-alphusw) * (dexp (alphusw*dexp(-wt**2/4)) - 1)
            ENDDO
      ELSEIF (ifc.eq.293) THEN  ! WEG DAMIT
         ni = NINT(P(2)) !Artem: Replace with a standard NINT function: jidnnt(P(2))
         tau = P(1)
         DO i = 1, n
            tau1 = X(i)
            f = 0
            DO ii = 1, ni
               phi = (ii-0.5)*360.d0/ni
               f = f + dsind(tau1)*dsqrt0(1-
     * (dcosd(tau1)*dcosd(tau)+dsind(tau1)*dsind(tau)*dcosd(phi))**2)
               ENDDO
            Y(i) = f / ni
            ENDDO
      ELSEIF (ifc.ge.300 .and. ifc.le.399) THEN ! NFS special routine
         CALL NFS_Fit (ifc-300, P, n, X, Y, Fehler)
      ELSEIF (ifc.ge.400 .and. ifc.le.409) THEN ! MCT special routine
         CALL MCT_Fit (ifc-400, P, n, X, Y, Fehler)
      ELSEIF (ifc.ge.410 .and. ifc.le.419) THEN ! MCT special routine
         CALL MCT_SkaFu (ifc-400, P, n, X, Y, Fehler)
      ELSE
         Print *, 'CuFuVal / ifc :', ifc
         CALL Absturz ('CuFuVal', 'ifc o.o.r.')
         ENDIF

      END ! CuFuVal
