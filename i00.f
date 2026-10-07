C  ====================================================================
C
C      Library  IDA   :  Ingenious Data Analysis
C      Modul    i00   :     main program
C
C  ====================================================================

C  --------------------------------------------------------------------
C  --------------------------------------------------------------------
C      LICENSE TERMS (not to be modified except by the authors)
C  --------------------------------------------------------------------
C  FRIDA (fast reliable inelastic data analysis)
C  <http://frida.sourceforge.net> is a program for generic spectral
C  analysis, with many specialized routines for inelastic neutron
C  scattering. The FORTRAN version Frida-1 is an updated version of
C  Joachim Wuttke's IDA, with contributions from the community. The
C  maintainer is Florian Kargl <f_kargl@users.sourceforge.net>.
C  FRIDA is released under the GNU public license.
C  (C) Joachim Wuttke 1990-2001
C  (C) Florian Kargl 2006
C  --------------------------------------------------------------------
C
C  The package includes the following maintained modules:
C  i00.f - i99.f,l1.f-l6.f,g1.f,g2.f
C  i_*.f,l_*.f,g_*.f contain definitions
C  l0*.f contain system specific information
C  --------------------------------------------------------------------
C  --------------------------------------------------------------------

C      General information :
C         Link i00-i99,l0*,l1-l6,g1-g2,NAGLIB/DP
C         Home page for detailed documentation (to be announced
C         via sourceforge mailinglist)

C      Contents of the IDA modules :
C         i00 : main program
C         i01 : plot
C         i10 : input/output
C         i20 : on-line memory
C         i23 : directories, editing, r/z-handling
C         i25 : file copy, delete, make
C         i30 : file manipulations, auxiliary calculations
C         i32 : general functions and symbolic calculation
C         i40 : manipulations on data / per channel
C         i41 : manipulations on data / per spectrum
C         i42 : manipulations on data / per file
C         i43 : special manipulations on data
C         i50 : operations on data
C         i60 : curves and fit
C         i66 : collection of fit functions
C         i67 : very special fit functions
C         i70 : rescaling operations
C         i71 : Fourier transforms
C         i72 : neutron kinematics, constant-q
C         i73 : density of states
C         i74 : self absorption
C         i75 : nuclear forward scattering
C         i76 : mode coupling integration
C         i77 : full mode coupling model             ! LINK ONLY WHEN NEEDED
C         i80 : raw data read / neutron scattering
C         i87 : raw data read / light scattering
C         i94 : simulation / multiple scattering     ! LINK ONLY WHEN NEEDED
C         i95 : simulation / optics
C  16.02.2026 Artem Panchenko: Corrected several line breaks
C  ====================================================================
C  The Main Program :
C  ====================================================================

      PROGRAM  IDA

C  --------------------------------------------------------------------
C     Declarations :
C  --------------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'

      CHARACTER*80  Fehler, Object, Word, aus
      CHARACTER     Stream*400
      INTEGER       MemBlockInq, JList(MF), nJList, lj, iWord !Artem: not used: nF, nFold
      LOGICAL       qOv
      INTEGER       JLMem(MLM,MF), NJLMem(MLM), i, j !Artem: Add global variables for the stack.
      common /JLMEMORY/ JLMem, NJLMem                !Artem: Add global variables for the stack.

C  --------------------------------------------------------------------
C     Initializations :
C  --------------------------------------------------------------------
      NJLMem=0 !Artem: Add global variables for the stack.
      JLMem=0  !Artem: Add global variables for the stack.
      Print *
      Print *,'*******************************************************'
      Print *
      Print *,'This is FRIDA-1 (fast reliable inelastic data analysis)'
      Print *,'FORTRAN release by Florian Kargl and Joachim Wuttke'
      Print *,'Linux Version 1.4 - May 2010 '
      Print *
      Print *,'*******************************************************'
      Print *

      Fehler = '&ff'
C      nF = 0 !Artem: nF not used
      nJList = 0
      Stream   = ';!me'

      ! search for the setup file :
      CALL ExeML ('\i ida.su', aus)
      IF (aus(1:4).ne.'\i +') CALL ExeML ('\i ../ida.su', aus)
      IF (aus(1:4).ne.'\i +') CALL ExeML ('\i ../../ida.su', aus)
      IF (aus(1:4).ne.'\i +') CALL ExeML ('\i /home/yang_fa/FRIDA1/'//
     *                          'ida.su', aus)
      IF (aus(1:4).ne.'\i +') CALL ExeML ('\i /usr/local/lib/ida.su',
     *                                    aus)
      IF (aus(1:4).ne.'\i +') THEN
         Print *
         Print *, ' Could not load setup file ida.su'
         Print *
         ENDIF

      Object = ' '
      CALL FileLoad  (Object, Fehler)
      IF (Fehler.ne.'&ff') CALL FehlerGong (Fehler, 3)
      CALL LiGetDef  (nJList, JList)
      Print *
      Print *, 'You reach the IDA command level. For help, type h'
      Print *

C  --------------------------------------------------------------------
C     Main loop : execute the commands contained in Stream
C  --------------------------------------------------------------------

 1    CONTINUE

      CALL NextCommand (Stream, Word, Object, iWord, nJList,JList,qOv)

C  --------------------------------------------------------------------
      IF (Word.eq.'h') THEN                                     !  Help
C  --------------------------------------------------------------------

      Print *, 'help :'
      Print *, '   hc  = commands'
      Print *, '   hd  = array dimensions'

      ELSEIF (Word.eq.'hc') THEN
         Print *, 'command groups / type first letter to obtain '//
     *                                'full lists:'
         Print *, '   f*  = files in on-line-memory (load,save,'//
     *                                'make,delete)'
         Print *, '   d*  = directories of files in on-line-memory'
         Print *, '   e*  = edit files (in particular file headers)'
         Print *, '   g*  = graphics'
         Print *, '   c*  = fit curves'
         Print *, '   m*  = manipulate data'
         Print *, '   o*  = operate on data'
         Print *, '   t*  = transform on data'
         Print *, '   r*  = raw data input'
         Print *, '   _*  = very special data treatment'
         Print *, 'further commands:'
         Print *, '   p   = plot'
         Print *, '   a   = add to plot'
         Print *, '   qui = quit IDA'

      ELSEIF (Word.eq.'hd') THEN
         CALL MemDims ()
         Print *
         CALL GraDims ()

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'f') THEN                                !  Files
C  --------------------------------------------------------------------

      Print *, 'files :'
      Print *, '   fl  = load file from disk into on-line-memory'
      Print *, '   fs  = save file from on-line-memory on disk'
      Print *, '   fc  = copy (duplicate) on-line-file'
      Print *, '   fx  = move file to last position in on-line-memory'
      Print *, '   fdel= delete file from on-line-memory'
      Print *, '   fm  = make a new file starting from nothing'
      Print *, '   fmh = make a new histogram file from an event log'

      ELSEIF (Word.eq.'fl') THEN
         CALL FileLoad  (Object, Fehler)
      ELSEIF (Word.eq.'fc') THEN
         CALL FileCopy (nJList, JList, Fehler)
      ELSEIF (Word.eq.'fdel') THEN
         CALL FileKill (nJList, JList, Fehler)
         CALL LiStackIO (0, nJList, JList)
      ELSEIF (Word.eq.'fx') THEN
         CALL FileCopy (nJList, JList, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99
         CALL FileKill (nJList, JList, Fehler)
         CALL LiStackIO (0, nJList, JList)
      ELSEIF (Word.eq.'fs') THEN
         CALL FileSave (nJList, JList, Fehler)
      ELSEIF (Word.eq.'fm') THEN
         CALL FileMake (Fehler)
      ELSEIF (Word.eq.'fmh') THEN
         CALL FileMakeHist (Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'d') THEN                            !  Directory
C  --------------------------------------------------------------------

      Print *, 'directories :'
      Print *, '   df  = directory of files in on-line-memory'
      Print *, '   ds  = directory of spectra in selected files'
      Print *, '   dz  = directory of z-values in selected files'
      Print *, '   dp  = listing of data points in selected spectra'

      ELSEIF (Word.eq.'df') THEN
         CALL MemInfoF (Fehler)
      ELSEIF (Word.eq.'ds') THEN
         CALL MemInfoK (nJList, JList, Fehler)
      ELSEIF (Word.eq.'dz') THEN
         CALL MemInfoZ (nJList, JList, Fehler)
      ELSEIF (Word.eq.'dp') THEN
         CALL MemInfoY (nJList, JList, Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'e') THEN                            !  Edit
C  --------------------------------------------------------------------

      Print *, 'edit :'
      Print *, '   ec  = edit coordinate names and units'
      Print *, '   ed  = edit documentation'
      Print *, '   ez  = edit z coordinates'
      Print *, '   er  = edit real parameters'
      Print *, '   ei  = edit integer parameters'
      Print *, '   eg  = edit graphics parameters'

      ELSEIF (Word.eq.'ec') THEN
         CALL EditCnu (nJList, JList, '?', Fehler)
      ELSEIF (Word.eq.'ed') THEN
         CALL EditDoc (nJList, JList, Fehler)
      ELSEIF (Word.eq.'ez') THEN
         CALL EditZ (nJList, JList, Fehler)
      ELSEIF (Word.eq.'er') THEN
         CALL EditRPar (nJList, JList, Fehler)
      ELSEIF (Word.eq.'ei') THEN
         CALL EditIPar (nJList, JList, Fehler)
      ELSEIF (Word.eq.'eg') THEN
         CALL EditGPar (nJList, JList, Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word(1:1).eq.'g') THEN                        !  Graphics
C  --------------------------------------------------------------------


         IF     (Word.eq.'gs') THEN
            CALL GraSoftCopy (Object, 'fil-gra-def.ps', Fehler)
         ELSEIF (Word.eq.'gp') THEN
            CALL GraSoftCopy (Object, 'fil-gra-ful.ps', Fehler)
         ELSEIF (Word.eq.'ga') THEN
            CALL GraSoftCopy (Object, 'fil-gra-app.ps', Fehler)
         ELSE
            CALL GraChoice (Word(2:5), Object, Fehler)
            ENDIF

      ELSEIF (Word.eq.'p' .or. Word.eq.'a' .or. Word.eq.'pp') THEN
         CALL IdaPlot (Word, nJList, JList, Object, Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'c') THEN                                ! Curves
C  --------------------------------------------------------------------

      Print *, 'curves :'
      Print *, '   cc  = create'
      Print *, '   cf  = fit'
      Print *, '   cp  = parameters'
      Print *, '   cpa = parameters (alt)'
      Print *, '   ca  = auxiliary parameters'
      Print *, '   cnn = new number (change function definition)'
      Print *, '   cs  = setup for fitroutine'
      Print *, '   ci  = get parameters (-> integral file)'
      Print *, '   cg  = get data on a grid (-> full file)'

      ELSEIF (Word.eq.'cc') THEN
         CALL CuCreate (nJList, JList, .false., Fehler)
      ELSEIF (Word.eq.'cf') THEN
         CALL CuFitCall (nJList, JList, Fehler)
      ELSEIF (Word.eq.'cp') THEN
         CALL EditCPar (nJList, JList, .false., Fehler)
      ELSEIF (Word.eq.'cpa') THEN
         CALL CuSetPar (nJList, JList, .false., Fehler)
      ELSEIF (Word.eq.'cnn') THEN
         CALL CuRedef (nJList, JList, Fehler)
      ELSEIF (Word.eq.'ca') THEN
         CALL CuSetAux (nJList, JList, Fehler)
      ELSEIF (Word.eq.'cs') THEN
         CALL CuSetFit (Fehler)
      ELSEIF (Word.eq.'ci') THEN
         CALL CuGetPar (nJList, JList, Fehler)
      ELSEIF (Word.eq.'cg') THEN
         CALL GridCurve (nJList, JList, Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'m') THEN                   !  Data manipulations
C  --------------------------------------------------------------------

      Print *, 'manipulations :'         ! systematic arrangement 5nov91
      Print *, '   mcd = channels delete'
      Print *, '   mca =          add'
      Print *, '   mcaa=          add automatically'
      Print *, '   mco =          order (sort/sum)'
      Print *, '   mcg =          determine groups'
      Print *, '   mcm =          y -> f(y)'
      Print *, '   mch =          histogram binning'
      Print *, '   mcx =          x <-> y'
      Print *, '   mcs =          break into spectra'
      Print *, '   mgi = grid     interpolate'
      Print *, '   mge =          extrapolate'
      Print *, '   mgd =          delete'
      Print *, '   mga =          add'
      Print *, '   mgr =          redistribute'
      Print *, '   mgh =          points -> histogram'
      Print *, '   msd = spectra  delete'
      Print *, '   msa =          add'
      Print *, '   msaw=          add (weighted)'
      Print *, '   msj =          join'
      Print *, '   mso =          order'
      Print *, '   msx =          exchange <-> channels'
      Print *, '   mfj = files    join'
      Print *, '   mfs =          sum'
      Print *, '   mfx =          exchange <-> spectra'

      ELSEIF (Word.eq.'mca') THEN
         CALL OrgChSum (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mcaa') THEN
         CALL OrgChSAuto (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mcd') THEN
         CALL OrgChCut (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mco') THEN
         CALL OrgChSort (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mcg') THEN
         CALL OrgChGroup (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mcg') THEN
         CALL OrgChGroup (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mcm') THEN
         CALL OrgChSubs (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mch') THEN
         CALL OrgHistMake (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mcx') THEN
         CALL OrgChExch (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mcs') THEN
         CALL OrgChSpectra (nJList, JList, qOv, Fehler)

      ELSEIF (Word.eq.'mgi') THEN
         CALL GridIntExt (.false., nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mge') THEN
         CALL GridIntExt (.true., nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mga') THEN
         CALL GridSumCut (.true., nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mgd') THEN
         CALL GridSumCut (.false., nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mgr') THEN
         CALL GridRedis (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mgh') THEN
         CALL GridHist (nJList, JList, qOv, Fehler)

      ELSEIF (Word.eq.'msd') THEN
         CALL OrgSpectraCut (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'msa') THEN
         CALL OrgSpectraSum ('a', .false., nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'msaw') THEN
         CALL OrgSpectraSum ('a', .true.,  nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'msj') THEN
         CALL OrgSpectraSum ('j', .false., nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'mso') THEN
         CALL OrgSpectraSort (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'msx') THEN
         CALL OrgSpectraExch (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'msf') THEN
         CALL OrgSpectraBreak (nJList, JList, Fehler)

      ELSEIF (Word.eq.'mfj') THEN
         CALL OrgFileJoin (nJList, JList, qOv, Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'o') THEN                   !  Operations on data
C  --------------------------------------------------------------------

      Print *, 'operations :'
      Print *, '   oi  = integral properties'
      Print *, '   ox  = pointwise operation on x'
      Print *, '   oxs =   dito, selected subrange'
      Print *, '   oy  = pointwise operation on y'
      Print *, '   oys =   dito, selected subrange'
      Print *, '   oz  = pointwise operation on z (also oz1, oz2, ...)'
      Print *, '   of  = operate on y as function of x   '
      Print *, '   ot  = form tensor product'

      ELSEIF (Word.eq.'oi') THEN
         CALL OprIntegral (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'of') THEN
         CALL OprDifferential (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'ot') THEN
         CALL OprTensor (nJList, JList, qOv, Fehler)

      ELSEIF (Word.eq.'ox') THEN
         CALL OprPointwise ('x', nJList, JList, .false., qOv,
     *                           Object, Fehler)
      ELSEIF (Word.eq.'oxs') THEN
         CALL OprPointwise ('x', nJList, JList, .true.,  qOv,
     *                           Object, Fehler)
      ELSEIF (Word.eq.'oy') THEN
         CALL OprPointwise ('y', nJList, JList, .false., qOv,
     *                           Object, Fehler)
      ELSEIF (Word.eq.'oys') THEN
         CALL OprPointwise ('y', nJList, JList, .true.,  qOv,
     *                           Object, Fehler)
      ELSEIF (Word.eq.'oz') THEN
         CALL OprPointwise ('z1', nJList, JList, .false., qOv,
     *                            Object, Fehler)
      ELSEIF (Word.eq.'oz#') THEN
         CALL OprPointwise ('z'//ch1(iWord), nJList, JList, .false.,
     *                            qOv, Object, Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'r') THEN
C  --------------------------------------------------------------------

      Print *, 'read raw data'
      Print *, '   rf     = from FPI @ E13'
      Print *, '   rr     = from Raman @ E13'
      Print *, '   ro     = from Oke @ Lens'
      Print *, '   re     = from NRSE @ LLB'
      Print *, '   rec    =      NRSE/ correct (what for ?)'
      Print *, '   r04    = from IN 4'
      Print *, '   r05    = from IN 5'
      Print *, '   r06    = from IN 6'
      Print *, '   rtof   = from TOFTOF (FRM2)'
      Print *, '   rMI    = from Mibemol'
      Print *, '   rNE    = from NEAT (V3) @ HMI'
      Print *, '   rfcs   = from FCS @ NIST'
      Print *, '   rfoc   = from Focus @ Sinq'
      Print *, '   rfoco  = from Focus @ Sinq (not NEXUS old)'
      Print *, '   rdcs   = from DCS @ NIST'
      Print *, '   r10    = from IN10'
      Print *, '   r13    = from IN13'
      Print *, '   r16    = from IN16'
      Print *, '   r10e   = from IN10/ elastic scans'
      Print *, '   rhfbs  = from HFBS @ NIST'
      Print *, '   rhfbo  = from HFBS @ NIST (pre Aug03)'
      Print *, '   rfans  = from FANS @ NIST'
      Print *, '   rbt2   = from BT2 @ NIST (+DMC PSI)'
      Print *, '   r11    = from IN11'
      Print *, '   rrs    = from IRIS @ ISIS'
      Print *, 'read log files'
      Print *, '   rh     = history files from light scattering'

c      ELSEIF (Word.eq.'re') THEN
c         CALL RRawNRSE (Fehler)
c      ELSEIF (Word.eq.'rec') THEN
c         CALL DCorrNRSE (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'rf') THEN
         CALL RRawFPI (Fehler)
      ELSEIF (Word.eq.'rr') THEN
         CALL RRawRam (Fehler)
      ELSEIF (Word.eq.'ro') THEN
         CALL RRawOke (Fehler)
      ELSEIF (Word.eq.'rhfbs') THEN
         CALL RRT_In_Hfbs (Fehler)
      ELSEIF (Word.eq.'rhfbo') THEN
         CALL RRT_In_Hfbso (Fehler)
      ELSEIF (Word.eq.'rfans') THEN
         CALL RRT_In_Fans (Fehler)
      ELSEIF (Word.eq.'rbt#'.and.iWord.eq.2) THEN
         CALL RRT_In_BT2 (Fehler)
      ELSEIF (Word.eq.'r#' .and. iWord.eq.10) THEN
         CALL RRawBS ('IN10', Fehler)
      ELSEIF (Word.eq.'r#' .and. iWord.eq.13) THEN
         CALL RRawBS ('IN13', Fehler)
      ELSEIF (Word.eq.'r#' .and. iWord.eq.16) THEN
         CALL RRawBS ('IN16', Fehler)
      ELSEIF (Word.eq.'r#' .and. iWord.eq.4) THEN
         CALL RRawTOF ('IN4', Fehler)
      ELSEIF (Word.eq.'r#' .and. iWord.eq.5) THEN
         CALL RRawTOF ('IN5', Fehler)
      ELSEIF (Word.eq.'r#' .and. iWord.eq.6) THEN
         CALL RRawTOF ('IN6', Fehler)
      ELSEIF (Word.eq.'rtof') THEN
         CALL RRawTOF ('TOF', Fehler)
      ELSEIF (Word.eq.'rmi') THEN
         CALL RRawTOF ('MIB', Fehler)
      ELSEIF (Word.eq.'rfoc') THEN
         CALL RRawTOF ('FOCUS', Fehler)
      ELSEIF (Word.eq.'rfoco') THEN
         CALL RRawTOF ('FOCUSO',Fehler)
      ELSEIF (Word.eq.'rdcs') THEN
         CALL RRawTOF ('DCS',Fehler)
      ELSEIF (Word.eq.'rrs') THEN
C         CALL RRaw_IRS ('IRS',Fehler) !Artem I couldn’t find any implementation of GETPARR() or GETDAT() for RRaw_IRSIn
         print *,
     * "ERROR: i00.f 491: Does not work — I could not find"//
     * " any implementation of GETPARR() or GETDAT() for RRaw_IRSIn"
      ELSEIF (Word.eq.'rfcs') THEN
         CALL RRawTOF ('FCS', Fehler)
      ELSEIF (Word.eq.'rne') THEN
         CALL RRawTOF ('NEAT', Fehler)
      ELSEIF (Word.eq.'r#e' .and. iWord.eq.10) THEN
         CALL RRawIN10e (Fehler)
      ELSEIF (Word.eq.'r#' .and. iWord.eq.11) THEN
         CALL RRawIN11 (Fehler)
      ELSEIF (Word.eq.'rh') THEN
         CALL RRawHistory (Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'t') THEN                           !  Transforms
C  --------------------------------------------------------------------

      Print *, 'transforms :'
      Print *, '   tu  = conversion of units in x (and possibly in y)'
      Print *, '   tfc = Fourier cosine transform with Filon algorithm'
      Print *, '   tfs = Fourier sine transform with Filon algorithm'
      Print *, '   tff = Fast Fourier cosine transform'
      Print *, '   tfp = Fast Fourier complex transform'
      Print *, '   tfe = Fast Fourier transform of exponential'
      Print *, '   ts  = (anti)symmetrize'
      Print *, '   tm  = double by adding mirror image Y(-x)'
      Print *, '   tr  = representation of complex numbers'
      Print *, '   td  = deconvolution'
      Print *, '   tc  = convolution'

      ELSEIF (Word.eq.'tu') THEN
         CALL TraUnits (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'tfc') THEN
         CALL TraFilon (nJList, JList, qOv, .true., Fehler)
      ELSEIF (Word.eq.'tfs') THEN
         CALL TraFilon (nJList, JList, qOv, .false., Fehler)
      ELSEIF (Word.eq.'tff') THEN
         CALL TraFFTsingle (nJList, JList, Fehler)
      ELSEIF (Word.eq.'tfp') THEN
         CALL TraFFTpair   (nJList, JList, Fehler)
      ELSEIF (Word.eq.'ts') THEN
         CALL TraSymm (nJList, JList, Fehler)
      ELSEIF (Word.eq.'tm') THEN
         CALL TraDouble (nJList, JList, qOv, Fehler)
c      ELSEIF (Word.eq.'tr') THEN
c         CALL TraRepr (Fehler)
c      ELSEIF (Word.eq.'td') THEN
c         CALL TraDeconv (Fehler)
c      ELSEIF (Word.eq.'tc') THEN
c         CALL TraConv (Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.'_') THEN
C  --------------------------------------------------------------------

      Print *, 'incorporated data analysis programs :'
      Print *, '   _af  = frequency axis for FPI data'
      Print *, '   _coq = interpolation to constant q'
      Print *, '   _mph = multiphonon correction for DOS'
      Print *, '   _muc = DOS calculation on basis of MUPHCOR'
      Print *, '   _pm  = +- convention for energies'
      Print *, '   _sac = self absorption coefficients'
      Print *, '   _sg  = conversions between S,S~,g,G,...'
c     Print *, '   _sub = test subtraction D(x)-r*C(x-d)'
      Print *, '   _tx  = nonlinear transformation of x axis'
      Print *, '   _mcc = Monte-Carlo convolution'
      Print *, '   _mss = multiple scattering simulation'
      Print *, '   _fmm = full mode-coupling model'
      Print *, '   _of  = Optik FPI'
      Print *, '   _rmc = retain master curve'
      Print *, '   _msd = MSD file generation IRIS'
      Print *, '   __   = test'

      ELSEIF (Word.eq.'_af') THEN
         CALL NormFPI (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'_coq') THEN
         CALL CoQ (nJList, JList, Fehler)
      ELSEIF (Word.eq.'_mph') THEN
         CALL DOS (nJList, JList, Fehler)
      ELSEIF (Word.eq.'_muc') THEN
         CALL MUPHDOS (nJList, JList, Fehler)
      ELSEIF (Word.eq.'_pm') THEN
         CALL SEGconv (nJList, JList, Fehler)
      ELSEIF (Word.eq.'_sac') THEN
         CALL AbsCoeffs (nJList, JList, Fehler)
      ELSEIF (Word.eq.'_sg') THEN
         CALL DOSconv (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'_tx') THEN
         CALL TraAchseX (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'_mcc') THEN
         CALL MC_Conv (nJList, JList, Fehler)
C      ELSEIF (Word.eq.'_mss') THEN                      ! LINK ONLY WHEN NEEDED
C         CALL MScat (nJList, JList, Fehler)
       ELSEIF (Word.eq.'_fmm') THEN                      ! LINK ONLY WHEN NEEDED
         CALL FullMCT (nJList, JList, Fehler)
      ELSEIF (Word.eq.'_of') THEN
         CALL OptikFPI (nJList, JList, qOv, Fehler)
      ELSEIF (Word.eq.'_rmc') THEN
         CALL RetainMaster (nJList, JList, Fehler)
      ELSEIF (Word.eq.'_msd') THEN
C         CALL MSDCalc (Fehler) !Artem I couldn’t find any implementation of GETPARR() or GETDAT() for RRaw_IRSIn
          print *,
     * "ERROR: i00.f 589: Does not work — I could not find"//
     * " any implementation of GETPARR() or GETDAT() for RRaw_IRSIn"
      ELSEIF (Word.eq.'__') THEN
         CALL IdaTest (nJList, JList, Fehler)

C  --------------------------------------------------------------------
      ELSEIF (Word(1:3).eq.'qui') THEN                          !  Exit
C  --------------------------------------------------------------------

      GOTO 999

C  --------------------------------------------------------------------
      ELSEIF (Word.eq.' ') THEN                           !  Do nothing
C  --------------------------------------------------------------------

C  --------------------------------------------------------------------
      ELSE                                           !  End of commands
C  --------------------------------------------------------------------

         Fehler = 'Unknown command : '//Word
         ENDIF

C  --------------------------------------------------------------------
 99   CONTINUE                                         !  Error message
C  --------------------------------------------------------------------

      IF (Fehler.ne.'&ff') THEN
         CALL FehlerGong (Fehler, 3)
         Stream = ';!me'
         ENDIF

C  --------------------------------------------------------------------
                                             !  Renew default file list
C  --------------------------------------------------------------------

      CALL FileClean ()
      CALL LiGetDef (nJList, JList)

      Print *

C  --------------------------------------------------------------------
      GOTO 1                                          !  End of main loop
C  --------------------------------------------------------------------

C  --------------------------------------------------------------------
 999  CONTINUE                                                  !  Stop
C  --------------------------------------------------------------------

      CALL GraChoice ('-', Object, Fehler)

      END ! Ida_Main

C      Implementation history :
C         1990ff.  VAX / VMS                ILL Grenoble
C         1992ff.  SUN                      LLB Saclay
C         1992     VAX / BSD-Unix           CCNY New York
C         1992     Macintosh                CCNY New York
C         1993ff.  DEC-Station / Ultrix     TU Muenchen
C         1994ff.  Alpha / OSF              TU Muenchen
C         1995     VAX / VMS                ILL Grenoble
C         1995     Silicon Graphics         ILL Grenoble

C  ====================================================================
C     Auxiliary routines (for menu)
C  ====================================================================

      SUBROUTINE NextCommand (Stream, Word, Object, iWord,
     *                        nJList, JList, qOv)
C     ----------------------------------------------------
            ! separated (for better readability of ida$main) JWu 29dec99
         ! Decompose Stream and find next command

      IMPLICIT NONE
      INCLUDE 'i_dim.f'
      INCLUDE 'l_def.f'
      CHARACTER*(*) Stream, Word, Object
      CHARACTER*80  Task, ein, aus, CJList, Fehler
      INTEGER       nJList, JList(*), iWord, niWord, i, kJLC, jB, nF,
     *              MemBlockInq
      LOGICAL       qOv

 1    CONTINUE
      Fehler = '&ff'

C Extract Task from Stream :
         ! The Stream consists of several Tasks which are separated
         ! from each other by a ";" delimiter. The last Task in the
         ! Stream has always to be the "me" command.

      CALL StaTake (1, Stream, Task)
      IF (Task(1:1).ne.';') THEN
         aus = 'Stream beginnt nicht mit '';'' : Stream = '''
     *           // Stream(1:lenU(Stream)) // '''.'
         CALL Absturz ('IDA/main', aus)
         ENDIF
      CALL TakeVorKla (Stream, Task, ';', Fehler)
      IF (Fehler.ne.'&ff') GOTO 99
      IF (Task(1:1).eq.'!') THEN
         CALL DelVonBis (Task, 1, 1)
      ELSE
         Print *, Task ! echo
         ENDIF

C  Extract Word from Task :
         ! A Task consists of a Word and optionally an nJList, JList.
         ! Word and nJList, JList are separated by a blank " " delimiter.
         ! If the word contains an integer, its value is saved as iWord
         ! and its string representation then replaced by a "#" symbol.

      jB     = min0 (jPos1(Task, ' '), jPos1(Task, '=')+1)
      Word   = Task(1:jB-1)        ! what comes before ' ' or '='
      Object = Task(jB:len(Task))  ! what comes after
      CALL DelLeft (Object)

      IF    (jPos1('0123456789=*?/', Word(1:1)).le.14) THEN
         CJList = Word                  ! Word seems to contain list of files
         Task   = Object                ! therefore what comes later
         jB     = jPos1 (Task, ' ')     !    must be divided further
         Word   = Task(1:jB-1)
         IF (jB+1.le.len(Task)) THEN
            Object = Task(jB+1:len(Task))
            CALL DelLeft (Object)
         ELSE
            Object = ' '
            ENDIF

         IF     (CJList.eq.'?') THEN ! help demanded
            IF (Word.eq.' ') THEN
               Object = 'Ida'       ! general help
            ELSE ! help on one command (Word)
               Object = Word
               ENDIF
            Word = 'h'
            nJList = 0  ! still needed ?
         ELSEIF (CJList.eq.'=') THEN
            ! no, no file list was given: overwrite the default files
            qOv = .true.
         ELSE
            ! overwrite files or not ?
            kJLC = lenU (CJList)
            IF (CJList(kJLC:kJLC).eq.'=') THEN
               qOv = .true.
               CJList(kJLC:kJLC) = ' '
            ELSE
               qOv = .false.
               ENDIF
            ! decode list of files :
            nF = MemBlockInq ('nF')
            CALL DecJList (CJList, MF, nJList, JList, 1, nF, Fehler)
            IF (Fehler.ne.'&ff') GOTO 99

            ENDIF ! new file list

      ELSE ! default file-list
         qOv = .false.

         ENDIF

      CALL FindN (Word, 1, niWord, iWord) ! replace integer in Word 1 by #
      CALL Minuskeln (Word)

C  Execute control commands :
      IF     (Word.eq.'me') THEN ! menue
            ! default file-list added 1feb93.
 20      CONTINUE
         CALL EncJList (nJList, JList, CJList)
         IF (nJList.gt.0) THEN
            CALL Compose2 (aus, 'IDA ('//CJList, ') ')
         ELSE
            aus = 'IDA ()'
            ENDIF
         CALL FrageC (aus, ein)

         IF (ein.eq.' ') THEN
            CALL LiStackIO (-1, nJList, JList)
            GOTO 20
            ENDIF
         CALL DelLeft (ein)
         Stream = ';!' // ein(1:lenU(ein)) // ';!me' ! '!' means no echo
         Print *
         GOTO 1

      ELSEIF (Word(1:2).eq.'#*') THEN ! n times a command (JWu 5jul91)

         CALL DelVonBis (Task, 1, 2)
         IF (niWord.ne.1) CALL Absturz
     *                  ('Main', '#* : no integer prepared')
         DO i = 1, iWord
            CALL Insert (Stream, 1, ';'//Task(1:lenU(Task)))
            ENDDO
         niWord = 0
         GOTO 1

         ENDIF

      RETURN ! regular exit

 99   CONTINUE
      IF (Fehler.ne.'&ff') THEN
         CALL FehlerGong (Fehler, 3)
         Stream = ';!me'
         GOTO 1 ! try again
         ENDIF

      END ! NextCommand

      SUBROUTINE LiGetDef (nJList, JList)
C     -----------------------------------
         ! Get new default-file-list

      INCLUDE 'i_dim.f'
      INTEGER       nF, nFB, nFA !Artem: Add memory data because nF was undetermined.
      REAL*8        FMem
      COMMON / OLM / FMem(Mmem), nFA(0:MB*MF), nFB(0:MF), nF
      INTEGER JList(MF)
      nFold = nF
      nF = MemBlockInq ('nF')
      IF     (nF.gt.nFold) THEN
         ! save old default
         CALL LiStackIO (+1, nJList, JList)
         ! set new default
         nJList = max0 (0, nF-nFold)
         print *, nJList, nF,nFold
         DO lj = 1, nJList
            JList(lj) = nFold + lj
            ENDDO
      ELSEIF (nF.lt.nFold) THEN ! after delete
         CALL LiStackIO (-1, nJList, JList)
         ENDIF

      END ! LiGetDef

      SUBROUTINE LiStackIO (io, nJList, JList)
C     ----------------------------------------
            ! 11mar93
      INCLUDE   'i_dim.f'
      INTEGER    JList(MF), JLMem(MLM,MF), NJLMem(MLM)
      common /JLMEMORY/ JLMem, NJLMem !Artem: Add global variables for the stack.

      IF     (io.eq.-1) THEN
         ! get and remove a list from stack :
         nJList = NJLMem(1)
         DO lj = 1, NJLMem(1)
            JList(lj) = JLMem(1,lj)
            ENDDO
         DO iLM = 1, MLM-1
            NJLMem(iLM) = NJLMem(iLM+1)
            DO lj = 1, NJLMem(iLM)
               JLMem(iLM,lj) = JLMem(iLM+1,lj)
               ENDDO
            ENDDO

      ELSEIF (io.eq.0) THEN
         ! remove entries from all lists :
         DO iLM = 1, MLM
            DO lj = 1, nJList
               DO ljj = 1, NJLMem(iLM)
                  IF (JLMem(iLM,ljj).eq.JList(lj)) THEN
                     DO ljjj = ljj, NJLMem(iLM)-1
                        JLMem(iLM,ljjj) = JLMem(iLM,ljjj+1)
                        ENDDO
                     NJLMem(iLM)= NJLMem(iLM)-1
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO

      ELSEIF (io.eq.+1) THEN
         ! put onto stack :
         DO iLM = MLM, 2, -1
            NJLMem(iLM) = NJLMem(iLM-1)
            DO lj = 1, NJLMem(iLM)
               JLMem(iLM,lj) = JLMem(iLM-1,lj)
               ENDDO
            ENDDO
         NJLMem(1) = nJList
         DO lj = 1, NJLMem(1)
            JLMem(1,lj) = JList(lj)
            ENDDO

      ELSE
         CALL Absturz('LiStackIO', 'io o.o.r.')
         ENDIF

      END ! LiStackIO
