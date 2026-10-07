      module nexus_api
  !--------------------------------------------------------------------
  ! NeXus (NAPI) Fortran interface helpers.
  !
  ! This module provides a thin, type-safe wrapper for the NeXus C API
  ! function NXIGETDATA(handle, void* data).
  !
  ! Key points:
  !   - NXhandle is defined in C as: typedef void* NXhandle;
  !     so the handle must be represented in Fortran as TYPE(C_PTR).
  !   - Both the handle and the data pointer must be passed BY VALUE to C,
  !     hence the VALUE attribute on the dummy arguments.
  !   - DATA should be passed as C_LOC(array(1)) from the caller.
  !   - Using a module provides an explicit interface automatically, which
  !     avoids ABI/calling-convention bugs and runtime crashes.
  !--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_long,
     * c_float, c_ptr, c_loc
      implicit none

      type, bind(C) :: NXlink
        integer(c_long) :: iTag
        integer(c_long) :: iRef
        character(c_char) :: targetPath(1024)
        integer(c_int) :: linkType
      end type NXlink

      contains

      integer(c_int) function NXGETINFO_PTR(handle, rank, dimension,
     *                                      datatype) result(status)
        !------------------------------------------------------------------
        ! Query metadata for the current NeXus dataset.
        !
        ! This is a Fortran wrapper around the NeXus C API call:
        !     NXstatus NXIGETINFO(NXhandle handle, int* rank,
        !                         int dimension[], int* datatype)
        !
        ! Arguments:
        !   handle    : NeXus object handle (NXhandle = void*) passed by value.
        !   rank      : (output) number of dimensions (C int*).
        !   dimension : (output) array of dimension sizes (C int[]).
        !               The caller must provide enough storage.
        !   datatype  : (output) NeXus datatype code (C int*).
        !
        ! Return:
        !   status    : NXstatus (typically 0 on success; nonzero on error).
        !
        ! Notes:
        !   - C dimension order is typically row-major (slowest index first),
        !     while Fortran is column-major (fastest index first). If the NeXus
        !     API returns dimensions in C order, reversing DIMENSION(1:RANK)
        !     makes it consistent with Fortran indexing expectations.
        !------------------------------------------------------------------
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none

        type(c_ptr),  value :: handle
        integer(c_int)      :: rank
        integer(c_int)      :: dimension(*)
        integer(c_int)      :: datatype
        integer             :: i
        integer(c_int)      :: tmp

        interface
          integer(c_int) function NXIGETINFO(h, r, dim, dt)
     * bind(C, name="nxigetinfo_")
            use, intrinsic :: iso_c_binding, only: c_int, c_ptr
            type(c_ptr),  value :: h
            integer(c_int)      :: r
            integer(c_int)      :: dim(*)
            integer(c_int)      :: dt
          end function NXIGETINFO
        end interface

        status = NXIGETINFO(handle, rank, dimension, datatype)

        ! Reverse DIMENSION(1:RANK) to convert C (row-major) ordering to
        ! Fortran (column-major) ordering, if required by your application.
        do i = 1, rank/2
          tmp = dimension(i)
          dimension(i) = dimension(rank-i+1)
          dimension(rank-i+1) = tmp
        end do

      end function NXGETINFO_PTR

      integer(c_int) function NXGETDATA_PTR(handle, data)
     * result(status)
        !------------------------------------------------------------------
        ! Read the current NeXus dataset into a user-supplied buffer.
        !
        ! This is a thin Fortran wrapper around the NeXus C API call:
        !     NXstatus NXIGETDATA(NXhandle handle, void* data)
        !
        ! Arguments:
        !   handle : NeXus object handle (NXhandle = void*) passed by value.
        !   data   : C pointer to the destination buffer (void*) passed by value.
        !            The caller is responsible for:
        !              - allocating/providing enough storage
        !              - using the correct element type and total size,
        !                consistent with NXgetinfo (rank/dimensions/datatype).
        !
        ! Return:
        !   status : NXstatus (typically 0 on success; nonzero on error).
        !
        ! Notes:
        !   - 'data' is an untyped (void*) pointer; no type conversion is done.
        !   - Typical usage is: NXGETDATA_PTR(handle, C_LOC(array(1))).
        !------------------------------------------------------------------
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none

        type(c_ptr), value :: handle   ! NXhandle (void*) passed by value
        type(c_ptr), value :: data     ! void* passed by value

        interface
          integer(c_int) function NXIGETDATA(h, p)
     * bind(C, name="nxigetdata_")
            use, intrinsic :: iso_c_binding, only: c_int, c_ptr
            type(c_ptr), value :: h
            type(c_ptr), value :: p
          end function NXIGETDATA
        end interface

        status = NXIGETDATA(handle, data)
      end function NXGETDATA_PTR

      integer(c_int) function NXGETSLAB_PTR(handle, data, start, size)
     * result(status)
        !------------------------------------------------------------------
        ! Read a hyperslab (subarray) of the current NeXus dataset into a
        ! user-supplied buffer.
        !
        ! This is a Fortran wrapper around the NeXus C API call:
        !     NXstatus NXIGETSLAB(NXhandle handle, void* data,
        !                         int start[], int size[])
        !
        ! Arguments:
        !   handle : NeXus object handle (NXhandle = void*) passed by value.
        !   data   : C pointer to destination buffer (void*) passed by value.
        !            The caller must provide enough storage for the slab.
        !   start  : Fortran-style slab start indices (1-based), one per dim.
        !            Expected length is at least RANK (see NXGETINFO_PTR).
        !   size   : Slab extents (number of elements) per dimension.
        !            Expected length is at least RANK.
        !
        ! Return:
        !   status : NXstatus (typically NX_OK on success; nonzero on error).
        !
        ! Notes:
        !   - The NeXus C API uses C ordering for dimension arrays and 0-based
        !     indices. Fortran is column-major and typically uses 1-based indices.
        !   - This wrapper:
        !       * queries RANK via NXGETINFO_PTR
        !       * reverses dimension order for the C call
        !       * converts START from 1-based (Fortran) to 0-based (C)
        !   - 'data' is a void*; no type conversion is done.
        !------------------------------------------------------------------
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none

        type(c_ptr), value :: handle
        type(c_ptr), value :: data
        integer(c_int)     :: start(*)
        integer(c_int)     :: size(*)

        integer(c_int), parameter :: NX_MAXRANK = 32
        integer(c_int), parameter :: NX_OK      = 1

        integer(c_int) :: cstart(NX_MAXRANK), csize(NX_MAXRANK)
        integer(c_int) :: rank, dimension(NX_MAXRANK), datatype
        integer        :: i

        interface
          integer(c_int) function NXIGETSLAB(h, p, s, n)
     * bind(C, name="nxigetslab_")
            use, intrinsic :: iso_c_binding, only: c_int, c_ptr
            type(c_ptr), value :: h
            type(c_ptr), value :: p
            integer(c_int)     :: s(*)
            integer(c_int)     :: n(*)
          end function NXIGETSLAB
        end interface

        ! Get dataset rank/dimensions (rank needed to interpret START/SIZE).
        status = NXGETINFO_PTR(handle, rank, dimension, datatype)
        if (status .ne. NX_OK) return

        if (rank .gt. NX_MAXRANK) then
          ! Rank exceeds local limits; return an error status (choose policy).
          status = -1
          return
        end if

        ! Convert Fortran inputs (column-major, 1-based) to C inputs
        ! (row-major ordering, 0-based start).
        do i = 1, rank
          cstart(i) = start(rank-i+1) - 1
          csize(i)  = size (rank-i+1)
        end do

        status = NXIGETSLAB(handle, data, cstart, csize)
      end function NXGETSLAB_PTR

      integer(c_int) function NXGETATTR_PTR(handle, name, data,
     * datalength, itype) result(status)
        !------------------------------------------------------------------
        ! Wrapper for:
        !   NXstatus NXigetattr(NXhandle handle, const char* name, void* data,
        !                       int* iDataLen, int* iType);
        !
        ! Uses EXTRACT_STRING to build a NUL-terminated C string in INAME(256)
        ! stored as bytes (INTEGER(C_INT8_T)).
        !------------------------------------------------------------------
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr,
     *   c_int8_t, c_loc
        implicit none

        type(c_ptr), value          :: handle
        character(*), intent(in)    :: name
        type(c_ptr), value          :: data
        integer(c_int), intent(out) :: datalength
        integer(c_int), intent(out) :: itype

        integer(c_int8_t) :: iname(256)

        interface
          integer(c_int) function NXIGETATTR(h, n, d, dl, t)
     * bind(C, name="nxigetattr_")
            use, intrinsic :: iso_c_binding, only: c_int, c_int8_t,
     * c_ptr
            type(c_ptr), value :: h
            integer(c_int8_t) :: n(*)     ! const char*
            type(c_ptr), value :: d     ! void*
            integer(c_int) :: dl        ! int*
            integer(c_int) :: t         ! int*
          end function NXIGETATTR

          subroutine EXTRACT_STRING(out, outlen, instr)
            use, intrinsic :: iso_c_binding, only: c_int8_t
            integer(c_int8_t), intent(out) :: out(*)
            integer,           intent(in)  :: outlen
            character(*),      intent(in)  :: instr
          end subroutine EXTRACT_STRING
        end interface

        ! Build NUL-terminated C string bytes in INAME
        call EXTRACT_STRING(iname, 256, name)

        ! Pass pointer to first byte as const char*
        status = NXIGETATTR(handle, iname, data,
     * datalength, itype)
      end function NXGETATTR_PTR

      integer(c_int) function NXPUTDATA_PTR(handle, data)
     * result(status)
        !------------------------------------------------------------------
        ! Write a user-supplied buffer into the current NeXus dataset.
        !
        ! This is a thin Fortran wrapper around the NeXus C API call:
        !     NXstatus NXIPUTDATA(NXhandle handle, void* data)
        ! (symbol name may be "nxiputdata_" depending on your library build)
        !
        ! Arguments:
        !   handle : NeXus object handle (NXhandle = void*) passed by value.
        !   data   : C pointer (void*) to the source buffer, passed by value.
        !            The caller is responsible for:
        !              - providing enough storage for the dataset
        !              - using the correct element type and layout, consistent
        !                with the dataset definition (rank/dimensions/datatype).
        !
        ! Return:
        !   status : NXstatus (typically 0 on success; nonzero on error).
        !
        ! Notes:
        !   - 'data' is an untyped (void*) pointer; no type conversion is done.
        !   - Typical usage is: NXPUTDATA_PTR(handle, C_LOC(array(1))).
        !------------------------------------------------------------------
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none

        type(c_ptr), value :: handle   ! NXhandle (void*) passed by value
        type(c_ptr), value :: data     ! void* passed by value

        interface
          integer(c_int) function NXIPUTDATA(h, p)
     * bind(C, name="nxiputdata_")
            use, intrinsic :: iso_c_binding, only: c_int, c_ptr
            type(c_ptr), value :: h
            type(c_ptr), value :: p
          end function NXIPUTDATA
        end interface

        status = NXIPUTDATA(handle, data)
      end function NXPUTDATA_PTR

      integer(c_int) function NXPUTSLAB_PTR(handle, data, start, size)
     * result(status)
        !------------------------------------------------------------------
        ! Write a hyperslab (subarray) from a user-supplied buffer into the
        ! current NeXus dataset.
        !
        ! This is a Fortran wrapper around the NeXus C API call:
        !     NXstatus NXIPUTSLAB(NXhandle handle, void* data,
        !                         int start[], int size[])
        !
        ! Arguments:
        !   handle : NeXus object handle (NXhandle = void*) passed by value.
        !   data   : C pointer (void*) to the source buffer, passed by value.
        !            The caller must provide enough storage for the slab.
        !   start  : Fortran-style slab start indices (1-based), one per dim.
        !            Expected length is at least RANK (see NXGETINFO_PTR).
        !   size   : Slab extents (number of elements) per dimension.
        !            Expected length is at least RANK.
        !
        ! Return:
        !   status : NXstatus (typically NX_OK on success; nonzero on error).
        !
        ! Notes:
        !   - The NeXus C API uses C ordering for dimension arrays and 0-based
        !     indices. Fortran is column-major and typically uses 1-based indices.
        !   - This wrapper:
        !       * queries RANK via NXGETINFO_PTR
        !       * reverses dimension order for the C call
        !       * converts START from 1-based (Fortran) to 0-based (C)
        !   - 'data' is a void*; no type conversion is done.
        !------------------------------------------------------------------
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none

        type(c_ptr), value :: handle
        type(c_ptr), value :: data
        integer(c_int)     :: start(*)
        integer(c_int)     :: size(*)

        integer(c_int), parameter :: NX_MAXRANK = 32
        integer(c_int), parameter :: NX_OK      = 1

        integer(c_int) :: cstart(NX_MAXRANK), csize(NX_MAXRANK)
        integer(c_int) :: rank, dimension(NX_MAXRANK), datatype
        integer        :: i

        interface
          integer(c_int) function NXIPUTSLAB(h, p, s, n)
     * bind(C, name="nxiputslab_")
            use, intrinsic :: iso_c_binding, only: c_int, c_ptr
            type(c_ptr), value :: h
            type(c_ptr), value :: p
            integer(c_int)     :: s(*)
            integer(c_int)     :: n(*)
          end function NXIPUTSLAB
        end interface

        ! Get dataset rank/dimensions (rank needed to interpret START/SIZE).
        status = NXGETINFO_PTR(handle, rank, dimension, datatype)
        if (status .ne. NX_OK) return

        if (rank .gt. NX_MAXRANK) then
          status = -1
          return
        end if

        ! Convert Fortran inputs (column-major, 1-based) to C inputs
        ! (row-major ordering, 0-based start).
        do i = 1, rank
          cstart(i) = start(rank-i+1) - 1
          csize(i)  = size (rank-i+1)
        end do

        status = NXIPUTSLAB(handle, data, cstart, csize)
      end function NXPUTSLAB_PTR

      integer(c_int) function NXPUTATTR_PTR(handle, name, data,
     * datalength, itype) result(status)
        !------------------------------------------------------------------
        ! Write (set) an attribute on the current NeXus object.
        !
        ! This is a Fortran wrapper around the NeXus C API call:
        !   NXstatus NXiputattr(NXhandle handle, const char* name, void* data,
        !                       int* iDataLen, int* iType)
        !
        ! Uses EXTRACT_STRING to build a NUL-terminated C string in INAME(256)
        ! stored as bytes (INTEGER(C_INT8_T)).
        !
        ! Arguments:
        !   handle     : NeXus handle (NXhandle = void*) passed by value.
        !   name       : Fortran attribute name; converted to C string bytes.
        !   data       : C pointer (void*) to the source value buffer, by value.
        !   datalength : (input) number of bytes in DATA to write (int* in C).
        !   itype      : (input) NeXus datatype code for the attribute (int* in C).
        !
        ! Return:
        !   status     : NXstatus (typically 0/NX_OK on success; nonzero on error).
        !
        ! Notes:
        !   - DATA is a void*; no type conversion is performed.
        !   - Typical usage: NXPUTATTR_PTR(h,'units',C_LOC(buf(1)),len,typ)
        !------------------------------------------------------------------
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr,
     * c_int8_t, c_loc
        implicit none

        type(c_ptr), value        :: handle
        character(*), intent(in)  :: name
        type(c_ptr), value        :: data
        integer(c_int), intent(in) :: datalength
        integer(c_int), intent(in) :: itype

        integer(c_int8_t), target :: iname(256)

        interface
          integer(c_int) function NXIPUTATTR(h, n, d, dl, t)
     * bind(C, name="nxiputattr_")
            use, intrinsic :: iso_c_binding, only: c_int,
     * c_ptr, c_int8_t
            type(c_ptr), value :: h
            type(c_ptr), value :: n     ! const char* (pointer to bytes)
            type(c_ptr), value :: d     ! void*
            integer(c_int) :: dl        ! int*
            integer(c_int) :: t         ! int*
          end function NXIPUTATTR

          subroutine EXTRACT_STRING(out, outlen, instr)
            use, intrinsic :: iso_c_binding, only: c_int8_t
            integer(c_int8_t), intent(out) :: out(*)
            integer,           intent(in)  :: outlen
            character(*),      intent(in)  :: instr
          end subroutine EXTRACT_STRING
        end interface

        ! Build NUL-terminated C string bytes in INAME (must include trailing 0 byte).
        call EXTRACT_STRING(iname, 256, name)

        ! Pass pointer to first byte as const char*
        status = NXIPUTATTR(handle, c_loc(iname(1)), data,
     * datalength, itype)
      end function NXPUTATTR_PTR

      end module nexus_api
