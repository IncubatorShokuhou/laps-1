!------------------------------------------------------------------------------
!m+
! name:
!       file_utility
!
! purpose:
!       module containing generic file utility routines
!
! category:
!       ncep rtm
!
! calling sequence:
!       use file_utility
!
! outputs:
!       none
!
! modules:
!       none.
!
! contains:
!       get_lun:     public function to return a free logical unit number for
!                    file access.
!
!       file_exists: public function to determine if a named file exists.
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       none.
!
! restrictions:
!       none.
!
! creation history:
!       written by:     paul van delst, cimss@noaa/ncep 12-jul-2000
!                       pvandelst@ncep.noaa.gov
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

module file_utility


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------------
  ! default visibility
  ! ------------------

  private


  ! ----------------------------------
  ! explicit visibility of subprograms
  ! ----------------------------------

  public :: get_lun, &
            file_exists


contains


!------------------------------------------------------------------------------
!s+
! name:
!       get_lun
!
! purpose:
!       public function to obtain a free logical unit number for file access
!
! calling sequence:
!       result = get_lun()
!
! input arguments:
!       none.
!
! output arguments:
!       none.
!
! function result:
!       function returns a default integer that can be used as a logical unit
!       number to open and access a file.
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
!       none.
!
! procedure:
!       the search for a free logical unit number begins at 10. the logical
!       unit number if tested to see if it is connected to an open file. if
!       so, it is incremented by 1. this is repeated until a free logical
!       unit number is found.
!s-
!------------------------------------------------------------------------------

  function get_lun() result( lun )


    ! -----------------
    ! type declarations
    ! -----------------
 
    integer :: lun
    logical :: file_open


    ! --------------------------------------------
    ! initialise logical unit number and file_open
    ! --------------------------------------------

    lun = 9
    file_open = .true.


    ! ------------------------------
    ! start open loop for lun search
    ! ------------------------------

    lun_search: do

      ! -- increment logical unit number
      lun = lun + 1

      ! -- check if file is open
      inquire( lun, opened = file_open )

      ! -- is this lun available?
      if ( .not. file_open ) exit lun_search

    end do lun_search

  end function get_lun



!------------------------------------------------------------------------------
!s+
! name:
!       file_exists
!
! purpose:
!       public function to determine if a file exists.
!
! calling sequence:
!       result = file_exists( file_name )
!
! input arguments:
!       file_name:  name of the file the existence of which is to be determined.
!                   units:      none
!                   type:       character
!                   dimension:  scalar
!                   attributes: intent( in )
!
! output arguments:
!       none.
!
! function result:
!       function returns a logical result.
!
!       result = .true.  => file exists
!              = .false. => file does not exist
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
!       none.
!
! procedure:
!       the file name is inquired by file keyword. the result of the inquiry
!       is the function result.
!s-
!------------------------------------------------------------------------------

  function file_exists( file_name ) result ( existence )


    ! -----------------
    ! type declarations
    ! -----------------
 
    character( * ), intent( in ) :: file_name
    logical :: existence


    ! ---------------
    ! inquire by name
    ! ---------------

    inquire( file = file_name, exist = existence )

  end function file_exists

end module file_utility


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
! revision 1.3  2000/08/31 19:36:32  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.2  2000/08/24 15:33:42  paulv
! - in the get_lun subprogram, the loop to search for a free unit number
!   was changed from:
!
!     do while ( file_open )
!       ...search
!     end do
!
!   to
!
!     lun_search: do
!       ...search
!       if ( .not. file_open ) exit lun_search
!     end do lun_search
!
!   the earlier version is a deprecated use of the do with while.
!
! - the subprogram file_exists was added. note that the inquire statement
!   required the file =  keyword to work. simply using the file name in
!   the inquire returned an error (compiler assumed it was an inquire by
!   unit number?)
! - updated module and subprogram documentation.
!
! revision 1.1  2000/07/12 16:08:10  paulv
! initial checked in version
!
!
!

