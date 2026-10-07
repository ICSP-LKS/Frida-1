C=======================================================================
C  napif_iface.inc  -- fixed-form include (F77 style)
C
C  Provides a generic interface NXGETDATA that accepts as DATA:
C    - INTEGER(C_INT8_T)  array  (typical BYTE / INTEGER*1)
C    - INTEGER(C_INT16_T) array  (typical INTEGER*2)
C    - INTEGER(C_INT)     array  (typical INTEGER*4)
C    - REAL(C_FLOAT)      array  (typical REAL*4)
C    - REAL(C_DOUBLE)     array  (typical REAL*8)
C    - TYPE(C_PTR)        array
C
C
C  Provides a generic interface NXGETATTR that accepts as DATA:
C    - INTEGER(C_INT)            (typical INTEGER*4)
C    - INTEGER(C_INT8_T)  array  (typical BYTE / INTEGER*1)
C    - INTEGER(C_INT16_T) array  (typical INTEGER*2)
C    - INTEGER(C_INT)     array  (typical INTEGER*4)
C    - REAL(C_FLOAT)             (typical REAL*4)
C    - REAL(C_DOUBLE)            (typical REAL*8)
C    - REAL(C_FLOAT)      array  (typical REAL*4)
C    - REAL(C_DOUBLE)     array  (typical REAL*8)
C    - TYPE(C_PTR)        array
C
C
C  Provides a generic interface NXPUTDATA that accepts as DATA:
C    - INTEGER(C_INT8_T)  array    (typical BYTE / INTEGER*1)
C    - INTEGER(C_INT16_T) array    (typical INTEGER*2)
C    - INTEGER(C_INT)     array    (typical INTEGER*4)
C    - INTEGER(C_INT)     array^2  (typical INTEGER*4)
C    - REAL(C_FLOAT)      array    (typical REAL*4)
C    - REAL(C_FLOAT)      array^2  (typical REAL*4)
C    - REAL(C_DOUBLE)     array    (typical REAL*8)
C    - REAL(C_DOUBLE)     array^2  (typical REAL*8)
C    - TYPE(C_PTR)        array
C
C
C  Provides a generic interface NXPUTATTR that accepts as DATA:
C    - INTEGER(C_INT8_T)         (typical BYTE / INTEGER*1)
C    - INTEGER(C_INT16_T)        (typical INTEGER*2)
C    - INTEGER(C_INT)            (typical INTEGER*4)
C    - REAL(C_FLOAT)             (typical REAL*4)
C    - REAL(C_DOUBLE)            (typical REAL*8)
C    - TYPE(C_PTR)        array
C
C
C  Provides a generic interface NXPUTSLAB that accepts as DATA:
C    - INTEGER(C_INT)              (typical INTEGER*4)
C    - INTEGER(C_INT)     array    (typical INTEGER*4)
C    - REAL(C_DOUBLE)              (typical REAL*8)
C    - REAL(C_DOUBLE)     array    (typical REAL*8)
C    - REAL(C_DOUBLE)     array^2  (typical REAL*8)
C    - TYPE(C_PTR)        array
C
C  Provides a generic interface NXMAKEDATA that accepts as DIM:
C    - INTEGER(C_INT)            (typical INTEGER*4)
C    - INTEGER(C_INT)     array  (typical INTEGER*4)
C
C
C  Provides a generic interface NXGETSLAB that accepts as DATA:
C    - INTEGER(C_INT)     array  (typical INTEGER*4)
C    - REAL(C_FLOAT)      array  (typical REAL*4)
C    - REAL(C_DOUBLE)     array  (typical REAL*8)
C    - TYPE(C_PTR)        array

C=======================================================================


C--- Rename ISO_C_BINDING kinds to avoid collisions with existing symbols
C      USE, INTRINSIC :: ISO_C_BINDING, ONLY :
C     *   C_INT     => C_INT,
C     *   FC_INT8_T  => C_INT8_T,
C     *   FC_INT16_T => C_INT16_T,
C     *   C_FLOAT      => C_FLOAT,
C     *   C_DOUBLE      => C_DOUBLE,
C     *   C_PTR     => C_PTR

      INTERFACE NXGETDATA

C----- Specific: DATA is INTEGER*4 rank-1 / C int buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_I4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*)
          INTEGER(C_INT) :: DATA
      END FUNCTION NXGETDATA_I4

C----- Specific: DATA is INTEGER*1 rank-1 / C 8-bit buffer -------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_NI1(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT8_T
          INTEGER(C_INT)    :: FILEID(*)
          INTEGER(C_INT8_T) :: DATA(*)
      END FUNCTION NXGETDATA_NI1

C----- Specific: DATA is INTEGER*2 rank-1 / C 16-bit buffer -------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_NI2(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT16_T
          INTEGER(C_INT)    :: FILEID(*)
          INTEGER(C_INT16_T) :: DATA(*)
      END FUNCTION NXGETDATA_NI2

C----- Specific: DATA is INTEGER*4 rank-1 / C int buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_NI4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*)
          INTEGER(C_INT) :: DATA(*)
      END FUNCTION NXGETDATA_NI4

C----- Specific: DATA is REAL*4 rank-1 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_R4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FLOAT
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_FLOAT) :: DATA
      END FUNCTION NXGETDATA_R4

C----- Specific: DATA is REAL*4 rank-1 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_NR4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FLOAT
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_FLOAT) :: DATA(*)
      END FUNCTION NXGETDATA_NR4

C----- Specific: DATA is REAL*8 rank-1 / C double buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_R8(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_DOUBLE) :: DATA
      END FUNCTION NXGETDATA_R8

C----- Specific: DATA is REAL*8 rank-1 / C double buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETDATA_NR8(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_DOUBLE) :: DATA(*)
      END FUNCTION NXGETDATA_NR8

C----- Specific: DATA is DATA is BYTE / C void buffer -------------------------
C      INTEGER(C_INT) FUNCTION NXGETDATA_PTR(FILEID, DATA)
C          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
C          INTEGER(C_INT) :: FILEID(*)
C          TYPE(C_PTR), VALUE :: DATA
C      END FUNCTION NXGETDATA_PTR

      END INTERFACE NXGETDATA

      INTERFACE NXGETATTR

C----- Specific: DATA is INTEGER*4 rank-0 / C int buffer ------------------------
      INTEGER(C_INT) FUNCTION NXGETATTR_I4(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*), DATA, DATALEN, TYPE
          CHARACTER*(*)   :: NAME
      END FUNCTION NXGETATTR_I4

C----- Specific: DATA is BYTE rank-1 / INTEGER*1 / C 8-bit buffer -----------------
      INTEGER(C_INT) FUNCTION NXGETATTR_NI1(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT8_T
          INTEGER(C_INT)    :: FILEID(*), DATALEN, TYPE
          CHARACTER*(*)      :: NAME
          INTEGER(C_INT8_T), contiguous :: DATA(:)
      END FUNCTION NXGETATTR_NI1

C----- Specific: DATA is BYTE rank-1 / INTEGER*2 / C 8-bit buffer -----------------
      INTEGER(C_INT) FUNCTION NXGETATTR_NI2(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT16_T
          INTEGER(C_INT)    :: FILEID(*), DATALEN, TYPE
          CHARACTER*(*)      :: NAME
          INTEGER(C_INT16_T), contiguous :: DATA(:)
      END FUNCTION NXGETATTR_NI2

C----- Specific: DATA is INTEGER*4 rank-1 / C int buffer ------------------------
      INTEGER(C_INT) FUNCTION NXGETATTR_NI4(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT8_T
          INTEGER(C_INT)    :: FILEID(*), DATA(*), DATALEN, TYPE
          CHARACTER*(*)      :: NAME
      END FUNCTION NXGETATTR_NI4

C----- Specific: DATA is REAL*4 rank-0 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETATTR_R4(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT,  C_FLOAT
          INTEGER(C_INT)    :: FILEID(*), DATALEN, TYPE
          CHARACTER*(*)      :: NAME
          REAL( C_FLOAT)        :: DATA
      END FUNCTION NXGETATTR_R4

C----- Specific: DATA is REAL*8 rank-0 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETATTR_R8(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT)    :: FILEID(*), DATALEN, TYPE
          CHARACTER*(*)     :: NAME
          REAL(C_DOUBLE)    :: DATA
      END FUNCTION NXGETATTR_R8

C----- Specific: DATA is REAL*4 rank-1 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETATTR_NR4(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FLOAT
          INTEGER(C_INT)            :: FILEID(*), DATALEN, TYPE
          CHARACTER*(*)             :: NAME
          REAL(C_FLOAT), contiguous :: DATA(:)
      END FUNCTION NXGETATTR_NR4

C----- Specific: DATA is REAL*8 rank-1 / C double buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXGETATTR_NR8(FILEID, NAME, DATA,
     *                                     DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT)             :: FILEID(*), DATALEN, TYPE
          CHARACTER*(*)              :: NAME
          REAL(C_DOUBLE), contiguous :: DATA(:)
      END FUNCTION NXGETATTR_NR8

C----- Specific: DATA is BYTE / C void buffer -------------------------
C      INTEGER(C_INT) FUNCTION NXGETATTR_PTR(FILEID, NAME, DATA,
C     *                                     DATALEN, TYPE)
C          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
C          INTEGER(C_INT)     :: FILEID(*), DATALEN, TYPE
C          CHARACTER*(*)      :: NAME
C          TYPE(C_PTR), VALUE :: DATA
C      END FUNCTION NXGETATTR_PTR

      END INTERFACE NXGETATTR


      INTERFACE NXPUTDATA

C----- Specific: DATA is INTEGER*1 rank-1 / C 8-bit buffer -------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NI1(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT8_T
          INTEGER(C_INT) :: FILEID(*)
          INTEGER(C_INT8_T), contiguous :: DATA(:)
      END FUNCTION NXPUTDATA_NI1

C----- Specific: DATA is INTEGER*2 rank-1 / C 16-bit buffer -------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NI2(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT16_T
          INTEGER(C_INT) :: FILEID(*)
          INTEGER(C_INT16_T), contiguous :: DATA(:)
      END FUNCTION NXPUTDATA_NI2

C----- Specific: DATA is INTEGER*4 rank-1 / C int buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NI4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*)
          INTEGER(C_INT), contiguous :: DATA(:)
      END FUNCTION NXPUTDATA_NI4

C----- Specific: DATA is INTEGER*4 rank-2 / C int buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NMI4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*)
          INTEGER(C_INT), contiguous :: DATA(:,:)
      END FUNCTION NXPUTDATA_NMI4

C----- Specific: DATA is REAL*4 rank-1 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NR4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FLOAT
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_FLOAT), contiguous :: DATA(:)
      END FUNCTION NXPUTDATA_NR4

C----- Specific: DATA is REAL*8 rank-1 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NR8(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_DOUBLE), contiguous :: DATA(:)
      END FUNCTION NXPUTDATA_NR8

C----- Specific: DATA is REAL*4 rank-2 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NMR4(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FLOAT
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_FLOAT), contiguous :: DATA(:,:)
      END FUNCTION NXPUTDATA_NMR4

C----- Specific: DATA is REAL*8 rank-2 / C float buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_NMR8(FILEID, DATA)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT) :: FILEID(*)
          REAL(C_DOUBLE), contiguous :: DATA(:,:)
      END FUNCTION NXPUTDATA_NMR8

C----- Specific: DATA is DATA is BYTE / C void buffer -------------------------
      INTEGER(C_INT) FUNCTION NXPUTDATA_PTR(FILEID, DATA)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
          INTEGER(C_INT) :: FILEID(*)
          TYPE(C_PTR), VALUE :: DATA
      END FUNCTION NXPUTDATA_PTR

      END INTERFACE NXPUTDATA


      INTERFACE NXPUTATTR

C----- Specific: DATA is INTEGER*1 rank-0 / C 8-bit buffer -------------------------
      INTEGER(C_INT) FUNCTION NXPUTATTR_I1(FILEID, NAME, DATA,
     * DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT8_T
          INTEGER(C_INT) FILEID(*), DATALEN, TYPE
          CHARACTER*(*) NAME
          INTEGER(C_INT8_T) :: DATA
      END FUNCTION NXPUTATTR_I1

C----- Specific: DATA is INTEGER*2 rank-0 / C 8-bit buffer -------------------------
      INTEGER(C_INT) FUNCTION NXPUTATTR_I2(FILEID, NAME, DATA,
     * DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_INT16_T
          INTEGER(C_INT) FILEID(*), DATALEN, TYPE
          CHARACTER*(*) NAME
          INTEGER(C_INT16_T) :: DATA
      END FUNCTION NXPUTATTR_I2

C----- Specific: DATA is INTEGER*4 rank-0 / C int buffer --------------------------------
      INTEGER(C_INT) FUNCTION NXPUTATTR_I4(FILEID, NAME, DATA,
     * DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) FILEID(*), DATALEN, TYPE
          CHARACTER*(*) NAME
          INTEGER(C_INT) :: DATA
      END FUNCTION NXPUTATTR_I4

C----- REAL*4 rank-0: DATA -----------------------------------------
      INTEGER(C_INT) FUNCTION NXPUTATTR_R4(FILEID, NAME, DATA,
     * DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FLOAT
          INTEGER(C_INT) FILEID(*), DATALEN, TYPE
          CHARACTER*(*) NAME
          REAL(C_FLOAT) :: DATA
      END FUNCTION NXPUTATTR_R4

C----- REAL*8 rank-0: DATA -----------------------------------------
      INTEGER(C_INT) FUNCTION NXPUTATTR_R8(FILEID, NAME, DATA,
     * DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT) FILEID(*), DATALEN, TYPE
          CHARACTER*(*) NAME
          REAL(C_DOUBLE) :: DATA
      END FUNCTION NXPUTATTR_R8

C----- Specific: BYTE / C void buffer -------------------------
      INTEGER(C_INT) FUNCTION NXPUTATTR_PTR(FILEID, NAME, DATA,
     * DATALEN, TYPE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
          INTEGER(C_INT) FILEID(*), DATALEN, TYPE
          CHARACTER*(*) NAME
          TYPE(C_PTR), VALUE :: DATA
      END FUNCTION NXPUTATTR_PTR

      END INTERFACE NXPUTATTR

      INTERFACE NXPUTSLAB

C----- Specific: DATA is INTEGER*4 rank-0 / C int buffer  -----------------
      INTEGER(C_INT) FUNCTION NXPUTSLAB_I4(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*), START(*), SIZE(*)
          INTEGER(C_INT) :: DATA
      END FUNCTION NXPUTSLAB_I4

C----- Specific: DATA is INTEGER*4 rank-1 / C int buffer  -----------------
      INTEGER(C_INT) FUNCTION NXPUTSLAB_NI4(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*), START(*), SIZE(*)
          INTEGER(C_INT), contiguous :: DATA(:)
      END FUNCTION NXPUTSLAB_NI4

C----- Specific: DATA is REAL*8 rank-0 / C double buffer  -----------------
      INTEGER(C_INT) FUNCTION NXPUTSLAB_R8(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT) :: FILEID(*), START(*), SIZE(*)
          REAL(C_DOUBLE) :: DATA
      END FUNCTION NXPUTSLAB_R8

C----- Specific: DATA is REAL*8 rank-1 / C double buffer  -----------------
      INTEGER(C_INT) FUNCTION NXPUTSLAB_NR8(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT)    :: FILEID(*), START(*), SIZE(*)
          REAL(C_DOUBLE), contiguous :: DATA(:)
      END FUNCTION NXPUTSLAB_NR8

C----- Specific: DATA is REAL*8 rank-2 / C double buffer  -----------------
      INTEGER(C_INT) FUNCTION NXPUTSLAB_NMR8(FILEID, DATA, START,
     * SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT)    :: FILEID(*), START(*), SIZE(*)
          REAL(C_DOUBLE), contiguous :: DATA(:,:)
      END FUNCTION NXPUTSLAB_NMR8

C----- Specific: DATA is Void* / C void* buffer  -----------------
      INTEGER(C_INT) FUNCTION NXPUTSLAB_PTR(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
          INTEGER(C_INT)    :: FILEID(*), START(*), SIZE(*)
          TYPE(C_PTR), VALUE :: DATA
      END FUNCTION NXPUTSLAB_PTR

      END INTERFACE NXPUTSLAB

      INTERFACE NXMAKEDATA

C----- Specific: DIM is INTEGER*4 rank-0 / C int buffer  -----------------
      INTEGER(C_INT) FUNCTION NXMAKEDATA_S(FILEID, LABEL, DATATYPE,
     * RANK, DIM)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*), DATATYPE, RANK, DIM
          CHARACTER*(*)   :: LABEL
      END FUNCTION NXMAKEDATA_S

C----- Specific: DIM is INTEGER*4 rank-0 / C int buffer  -----------------
      INTEGER(C_INT) FUNCTION NXMAKEDATA_A(FILEID, LABEL, DATATYPE,
     * RANK, DIM)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT) :: FILEID(*), DATATYPE, RANK, DIM(*)
          CHARACTER*(*)   :: LABEL
      END FUNCTION NXMAKEDATA_A

      END INTERFACE NXMAKEDATA

      INTERFACE NXGETSLAB

C----- Specific: DATA is INTEGER*4 rank-1 / C int buffer  -----------------
      INTEGER(C_INT) FUNCTION NXGETSLAB_NI4(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
          INTEGER(C_INT)    :: FILEID(*), START(*), SIZE(*)
          INTEGER(C_INT), contiguous :: DATA(:)
      END FUNCTION NXGETSLAB_NI4

C----- Specific: DATA is REAL*4 rank-1 / C double buffer  -----------------
      INTEGER(C_INT) FUNCTION NXGETSLAB_NR4(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FLOAT
          INTEGER(C_INT)    :: FILEID(*), START(*), SIZE(*)
          REAL(C_FLOAT), contiguous :: DATA(:)
      END FUNCTION NXGETSLAB_NR4

C----- Specific: DATA is REAL*8 rank-1 / C double buffer  -----------------
      INTEGER(C_INT) FUNCTION NXGETSLAB_NR8(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
          INTEGER(C_INT)    :: FILEID(*), START(*), SIZE(*)
          REAL(C_DOUBLE), contiguous :: DATA(:)
      END FUNCTION NXGETSLAB_NR8

C----- Specific: DATA is Void* / C void* buffer  -----------------
      INTEGER(C_INT) FUNCTION NXGETSLAB_PTR(FILEID, DATA, START, SIZE)
          USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
          INTEGER(C_INT)    :: FILEID(*), START(*), SIZE(*)
          TYPE(C_PTR), VALUE :: DATA
      END FUNCTION NXGETSLAB_PTR

      END INTERFACE NXGETSLAB
