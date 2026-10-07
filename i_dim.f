C  ====================================================================
C
C      Library  IDA   :  Inelastic data treatment
C      Include  i_dim :     Array dimensions   (version for DEC-Alphas)
C
C  ====================================================================

          ! The Makefile.alp copies this file to i_dim.f
          ! which is then included into modules i?.f

      INTEGER       MF, MLM, MK, MZ, MBH, MB, MC, MmemSpe, MmemBlo,
     *              MDtof, MP, MFP, MRP, MTP, MWrk3dim, MK3V,MS3V,Mmem
      PARAMETER    (MF=120,           ! max # on-line-files
     *              MLM=8,            ! max # lists for on-line-files !Artem: Move MLM from a local parameter in i00.f (LiStackIO()) to a global definition.
     *              MK=1024,          ! max # spectra per file
     *              MZ=19,           ! max # z's
     *              MBH=6,           ! # header blocks
     *              MB=MBH+4*MK,     ! max # blocks per file
     *              MC=16384,         ! max # x-y-d data points per spectrum
     *              MmemSpe=3*MC+1,  ! max # entries per spectrum
     *              MmemBlo=MmemSpe, ! max # entries per mem-block
     *              MDtof=1, ! for call to IDOL MDtof=512*128, else MDtof=1
     *              MP=40,           ! in 92'formats, max # i/r/t params
     *              MFP=MP-20,       ! max # fitparams
     *              MRP=20,          ! max # r-params
     *              MTP=50,          ! max # doc-lines
     *              MWrk3dim=3,      ! 3rd dimension of (xy-z-*)-array
                       ! requirements: coq(3)
     *              MK3V=1,          ! work space for i77/FullMCT (typ =MK)
     *              MS3V=8,          ! work space for i77/FullMCT (typ 8000)
     *              Mmem=40000000)    ! max # entries in on-line-memory
                       ! verkraftet bis mindestens 3M6
