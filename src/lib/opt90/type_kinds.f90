!------------------------------------------------------------------------------
!m+
! name:
!       type_kinds
!
! purpose:
!       module to hold specification kinds for variable declaration.
!
! category:
!       ncep rtm
!
! calling sequence:
!       use type_kinds
!
! outputs:
!       byte     specification kind for byte (1-byte) integer variable
!       short    specification kind for short (2-byte) integer variable
!       long     specification kind for long (4-byte) integer variable
!       llong    specification kind for double long (8-byte) integer variable
!       single   specification kind for single precision (4-byte) real variable
!       double   specification kind for double precision (8-byte) real variable
!       quad     specification kind for quad precision (16-byte) real variable
!
!       ip_kind: generic specification kind for default integer
!       fp_kind: generic specification kind for default floating point
!
! modules:
!       none
!
! contains:
!       type_size:  public function to return the number of bytes used to
!                   represent the data type.
!
! side effects:
!       if the llong or quad type kinds are not available they default to the
!       long and double kind specifications.
!
! restrictions:
!       none
!
! example:
!       use type_kinds
!       integer( long ) :: i, j
!       real( single )  :: x, y
!
! creation history:
!       written by:     paul van delst, cimss@noaa/ncep 12-jun-2000
!                       pvandelst@ncep.noaa.gov
!
!
!  copyright (c) 2000 paul van delst
!
!  this program is free software; you can redistribute it and/or
!  modify it under the terms of the gnu general public license
!  as published by the free software foundation; either version 2
!  of the license, or (at your option) any later version.
!
!  this program is distributed in the hope that it will be useful,
!  but without any warranty; without even the implied warranty of
!  merchantability or fitness for a particular purpose.  see the
!  gnu general public license for more details.
!
!  you should have received a copy of the gnu general public license
!  along with this program; if not, write to the free software
!  foundation, inc., 59 temple place - suite 330, boston, ma  02111-1307, usa.
!m-
!------------------------------------------------------------------------------

module type_kinds


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------------
  ! default visibility
  ! ------------------

  private
  public :: type_size


  ! ----------
  ! intrinsics
  ! ----------

  intrinsic abs,                &
            kind,               &
            selected_int_kind,  &
            selected_real_kind, &
            size,               &
            transfer,           &
            trim


  ! -------------------
  ! integer definitions
  ! -------------------

  ! -- integer types
  integer, parameter, public  :: byte    = selected_int_kind(1)   ! byte  integer
  integer, parameter, public  :: short   = selected_int_kind(4)   ! short integer
  integer, parameter, public  :: long    = selected_int_kind(8)   ! long  integer
  integer, parameter, private :: llong_t = selected_int_kind(16)  ! llong integer
  integer, parameter, public  :: llong   = max( llong_t, long )
!  integer, parameter, public  :: llong   = ( ( ( 1 + sign( 1, llong_t ) ) / 2 ) * llong_t ) + &
!                                           ( ( ( 1 - sign( 1, llong_t ) ) / 2 ) * long    )

  ! -- expected 8-bit byte sizes of the integer kinds
  integer, parameter, public :: n_bytes_for_byte_kind  = 1
  integer, parameter, public :: n_bytes_for_short_kind = 2
  integer, parameter, public :: n_bytes_for_long_kind  = 4
  integer, parameter, public :: n_bytes_for_llong_kind = 8

  ! -- define arrays for default definition
  integer, parameter, private :: n_ip_kinds = 4
  integer, parameter, dimension( n_ip_kinds ), private :: ip_kind_types = (/ byte,  &
                                                                             short, &
                                                                             long,  &
                                                                             llong  /) 
  integer, parameter, dimension( n_ip_kinds ), private :: ip_byte_sizes = (/ n_bytes_for_byte_kind,  &
                                                                             n_bytes_for_short_kind, &
                                                                             n_bytes_for_long_kind,  &
                                                                             n_bytes_for_llong_kind  /)

  ! -- default values

  ! **** change the following to change the default integer type kind ***
  integer, parameter, private :: iip = 3  ! 1=byte, 2=short, 3=long, 4=llong

  integer, parameter, public  :: ip_kind             = ip_kind_types( iip )
  integer, parameter, public  :: n_bytes_for_ip_kind = ip_byte_sizes( iip )


  ! --------------------------
  ! floating point definitions
  ! --------------------------

  ! -- floating point types
  integer, parameter, public  :: single = selected_real_kind(6)  ! single precision
  integer, parameter, public  :: double = selected_real_kind(15) ! double precision
  integer, parameter, private :: quad_t = selected_real_kind(20) ! quad precision
  integer, parameter, public  :: quad   = max( quad_t, double )
!  integer, parameter, public  :: quad   = ( ( ( 1 + sign( 1, quad_t ) ) / 2 ) * quad_t ) + &
!                                          ( ( ( 1 - sign( 1, quad_t ) ) / 2 ) * double )

  ! -- expected 8-bit byte sizes of the floating point kinds
  integer, parameter, public :: n_bytes_for_single_kind = 4
  integer, parameter, public :: n_bytes_for_double_kind = 8
  integer, parameter, public :: n_bytes_for_quad_kind   = 16

  ! -- define arrays for default definition
  integer, parameter, private :: n_fp_kinds = 3
  integer, parameter, dimension( n_fp_kinds ), private :: fp_kind_types = (/ single, &
                                                                             double, &
                                                                             quad    /) 
  integer, parameter, dimension( n_fp_kinds ), private :: fp_byte_sizes = (/ n_bytes_for_single_kind, &
                                                                             n_bytes_for_double_kind, &
                                                                             n_bytes_for_quad_kind    /)

  ! -- default values

  ! **** change the following to change the default floating point kind ***
  integer, parameter, private :: ifp = 2  ! 1=single, 2=double, 3=quad

  integer, parameter, public  :: fp_kind             = fp_kind_types( ifp )
  integer, parameter, public  :: n_bytes_for_fp_kind = fp_byte_sizes( ifp )


contains


!--------------------------------------------------------------------------------
!s+
! name:
!       type_size
!
! purpose:
!       public function to determine the size (in bytes) of a particular data type.
!
! category:
!       general
!
! calling sequence:
!       result = type_size( type_kind, &
!                           expected_type_size = expected_type_size )
!
! input arguments:
!       type_kind:   string describing the definition of the data type. valid
!                    values are (case sensitive):
!                      "byte"
!                      "short"
!                      "long"
!                      "llong"
!                      "single"
!                      "double"
!                      "quad"
!                      "ip_kind"
!                      "fp_kind"
!
! optional input arguments:
!       none.
!
! output arguments:
!       none.
!
! optional output arguments:
!       expected_type_size:  integer argument containing the *expected* size
!                            of the kind type in bytes. useful if you will be
!                            using 8-byte integer (llong) or 16-byte floating
!                            point (quad) types and want to check if what they
!                            *should* be is what they actually *are*. not all
!                            compilers support the llong or quad types.
!
!                            if included in argument list and an invalid or
!                            unrecognised type kind is specified, the returned
!                            value is 0.
!
! function result:
!       the returned value is the number of bytes used to represent the data
!       type on the platform the source code was compiled.
!
!       if an invalid or unrecognised type_kind is specified, the returned
!       result is -1.
!
! calls:
!       none.
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       none known.
!
! restrictions:
!       assumes that a single character == 1 byte (8 bits).
!s-
!--------------------------------------------------------------------------------

  function type_size( type_kind,         &
                      expected_type_size )

    ! -- arguments
    character( * ), intent( in  )           :: type_kind
    integer,        intent( out ), optional :: expected_type_size

    ! -- function
    integer                            :: type_size

    ! -- local variables
    integer :: esize
    character( 1 ), dimension( 1 ) :: x

    select case ( trim( type_kind ) )

      case ( 'byte' )
        type_size = size( transfer( 0_byte, x ) )
        esize     = n_bytes_for_byte_kind

      case ( 'short' )
        type_size = size( transfer( 0_short, x ) )
        esize     = n_bytes_for_short_kind

      case ( 'long' )
        type_size = size( transfer( 0_long, x ) )
        esize     = n_bytes_for_long_kind

      case ( 'llong' )
        type_size = size( transfer( 0_llong, x ) )
        esize     = n_bytes_for_llong_kind

      case ( 'ip_kind' )
        type_size = size( transfer( 0_ip_kind, x ) )
        esize     = n_bytes_for_ip_kind

      case ( 'single' )
        type_size = size( transfer( 0.0_single, x ) )
        esize     = n_bytes_for_single_kind

      case ( 'double' )
        type_size = size( transfer( 0.0_double, x ) )
        esize     = n_bytes_for_double_kind

      case ( 'quad' )
        type_size = size( transfer( 0.0_quad, x ) )
        esize     = n_bytes_for_quad_kind

      case ( 'fp_kind' )
        type_size = size( transfer( 0.0_fp_kind, x ) )
        esize     = n_bytes_for_fp_kind

      case default
        write( *, '( /5x, "invalid type kind: ", a )' ) trim( type_kind )
        type_size = -1
        esize     = 0

    end select

    if ( present( expected_type_size ) ) expected_type_size = esize

  end function type_size

end module type_kinds

!-------------------------------------------------------------------------------
!                          -- modification history --
!-------------------------------------------------------------------------------
!
! $id$
!
! $date$
!
! $revision$
!
! $state$
!
! $log$
! revision 2.6  2001/08/31 20:47:00  paulv
! - updated definitions such that when the default type definition is changed
!   so is the assumed byte size of the result.
! - commented out correct definition of type for llong and quad. pgi compiler
!   has a bug in it that does not allow elemental intrinsic functions to be
!   used in parameter initialisation expressions.
! - function type_size altered to optional return the "expected" byte size
!   of specified type. this allows the user to check if the requested kind
!   type is supported by the compiler, i.e. if the actual and expected sizes
!   do not agree, then the kind type is unsupported.
!
! revision 2.5  2001/07/12 16:43:16  paulv
! - replaced possible llong (8-byte integer) kind definition from
!     llong   = ( ( abs( llong_t ) + llong_t ) * llong_t + &
!                 ( abs( llong_t ) - llong_t ) * long ) / &
!               ( 2 * abs( llong_t ) )
!   to
!     llong   = max( llong_t, long )
! - replaced possible quad (16-byte floating point) kind definition from
!     quad   = ( ( abs( quad_t ) + quad_t ) * quad_t + &
!                ( abs( quad_t ) - quad_t ) * double ) / &
!              ( 2 * abs( quad_t ) )
!   to
!     quad   = max( quad_t, double )
! - added commented definition for quad precision kind.
! - added comment in type_size function header explaining reliance on
!   1 character = 8 bits for function to return storage.
! - removed len intrinsic from character definitions.
!
! revision 2.4  2001/03/26 22:35:18  paulv
! - renamed ip_precision and fp_precision parameters with ip_kind and
!   fp_kind.
! - initialisation of ip_kind and fp_kind are now defined from the definitions
!   of the "fundamental" types. currently ip_kind = long and fp_kind = double.
!   this was changed from the default kind typing using kind( 0 ) and kind( 0.0 )
!   respectively to allow the user to easily and explicitly redefine some sort
!   of default integer and floating point kind.
!
! revision 2.3  2000/08/31 19:36:34  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 2.2  2000/08/31 15:55:34  paulv
! - added documentation delimiters.
! - added documentation for type_size function.
! - changed default module visibility from public to private.
! - added true default integer and floating point types.
!
! revision 2.1  2000/08/08 17:08:55  paulv
! - added definitions of 64-bit integers and reals. definitions default to
!   the largest integer and real available on a system if not available.
! - added type_size function to return the number of bytes used by a defined
!   data type - both integer and real.
!
! revision 1.1  2000/07/12 16:08:11  paulv
! initial checked in version
!
!
!
