!------------------------------------------------------------------------------
!m+
! name:
!       error_handler
!
! purpose:
!       module to define simple error codes and handle error conditions
!
! category:
!       ncep rtm
!
! calling sequence:
!       use error_handler
!
! outputs:
!       success:     code specifying successful completion.
!                    units:      none.
!                    type:       integer
!                    dimension:  scalar
!                    attributes: parameter, public
!
!       information: code specifying information output.
!                    units:      none.
!                    type:       integer
!                    dimension:  scalar
!                    attributes: parameter, public
!
!       warning:     code specifying warning state. execution can
!                    continue but results may be incorrect.
!                    units:      none.
!                    type:       integer
!                    dimension:  scalar
!                    attributes: parameter, public
!
!       failure:     code specifying severe error. execution cannot
!                    continue.
!                    units:      none.
!                    type:       integer
!                    dimension:  scalar
!                    attributes: parameter, public
!
!       undefined:   code specifying undefined completion status.
!                    units:      none.
!                    type:       integer
!                    dimension:  scalar
!                    attributes: parameter, public
!
!
! modules:
!       file_utility: module containing global file utility routines.
!                     only the get_lun() function is used in this
!                     module.
!
! contains:
!       display_message:  public subroutine to display error/status messages
!                         either to standard output (default) or to a log file.
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
! example:
!       use error_handler
!       error_status = calculate_widget_size()
!       if ( error_status /= success ) then
!         call display_message( routine_name, &
!                               'error calculating widget size', &
!                               error_status, &
!                               message_log = 'error_log.txt' )
!         return
!       end if
!
! creation history:
!       written by:     paul van delst, cimss@noaa/ncep 12-jun-2000
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

module error_handler


  ! ---------------------
  ! module use statements
  ! ---------------------
  ! fsl modification (bls, 26 nov 02):
  !  commented out the use file_utility statement here, because
  !  it causes problems when using pgf90.  it appears to pgf90
  !  that get_lun is defined in two different modules.  this
  !  use statement is instead placed into the subroutine below
  !  that actually calls "get_lun"
  !use file_utility, only: get_lun


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------------
  ! default visibility
  ! ------------------

  private


  ! ------------------------------------
  ! definitions of public parameter data
  ! ------------------------------------
 
  ! -- integer values that define the error state.
  ! -- note: these values are totally arbitrary. 
  integer, parameter, public :: success     = 0
  integer, parameter, public :: information = success + 1
  integer, parameter, public :: warning     = information + 1
  integer, parameter, public :: failure     = warning + 1
  integer, parameter, public :: undefined   = failure + 1


  ! -----------------------------------
  ! definitions of local parameter data
  ! -----------------------------------

  ! -- character descriptors of the error states
  integer,         parameter :: max_n_states = 5
  character( 11 ), parameter, dimension( 0:max_n_states-1 ) :: &
    state_descriptor = (/ 'success    ', &
                          'information', &
                          'warning    ', &
                          'failure    ', &
                          'undefined  ' /)


  ! ----------------------------------
  ! explicit visibility of subprograms
  ! ----------------------------------

  public :: display_message


contains



!------------------------------------------------------------------------------
!s+
! name:
!       display_message
!
! purpose:
!       recursive public routine to display messages.
!
! calling sequence:
!       call display_message( routine_name, &
!                             message,      &
!                             error_state,  &
!                             message_log  = message_log )
!
! input arguments:
!       routine_name: name of the routine in which the message originated.
!                     units:      none
!                     type:       character
!                     dimension:  scalar
!                     attributes: intent( in )
!
!       message:      message text
!                     units:      none
!                     type:       character
!                     dimension:  scalar
!                     attributes: intent( in )
!
!       error_state:  flag corresponding to one of the defined error states.
!                     if not, the error state is set to undefined.
!                     units:      none
!                     type:       integer
!                     dimension:  scalar
!                     attributes: intent( in )
!
! optional input arguments:
!       message_log:  character string specifying a filename in which any
!                     messages will be logged. if not specified, or if an
!                     error occurs opening the log file, the default action
!                     is to output messages to the screen.
!                     units:      none
!                     type:       character
!                     dimension:  scalar
!                     attributes: intent( in ), optional
!
! calls:
!      get_lun:   function to return a free logical unit number for
!                 file access.
!                 source: file_utility module
!
!      routine calls itself if the optional argument message_log is passed and
!      an error occurs opening the output log file.
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
!       output message format is:
!
!         "routine name"("state description") : "message"
!
!       for example, if an error occurs in this routine the output is:
!
!         "display_message(failure) : error opening message log file"
!s-
!------------------------------------------------------------------------------

  recursive subroutine display_message ( routine_name, &
                                         message,      &
                                         error_state,  &
                                         message_log   )


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#
    ! fsl modification (bls, 26 nov 02)
    !  added the use file_utility statement originally at top of
    !  module to this private subroutine.
    use file_utility, only: get_lun
    ! ---------
    ! arguments
    ! ---------

    character( * ), intent( in )           :: routine_name
    character( * ), intent( in )           :: message
    integer,        intent( in )           :: error_state
    character( * ), intent( in ), optional :: message_log


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: this_routine_name = 'display_message'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: error_state_to_use
    integer :: log_to_file
    integer :: file_id
    integer :: io_status

    character( 28 ) :: fmt_string


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic present, &
              trim


    !#--------------------------------------------------------------------------#
    !#                   -- check the input error state --                      #
    !#--------------------------------------------------------------------------#

    error_state_to_use = error_state
    if ( error_state < 0 .or. error_state > max_n_states ) then
      error_state_to_use = undefined
    end if



    !#--------------------------------------------------------------------------#
    !#      -- set the message log. if not specified, output to screen --       #
    !#--------------------------------------------------------------------------#

    if ( present( message_log ) ) then

      log_to_file = 1
      file_id     = get_lun()

      open( file_id, file     = message_log,  &
                     access   = 'sequential', &
                     form     = 'formatted',  &
                     status   = 'unknown',    &
                     position = 'append',     &
                     action   = 'readwrite',  & ! just read may cause probs on some
                                                ! systems using position = 'append'
                     iostat   = io_status )

      if ( io_status /= 0 ) then
        call display_message( this_routine_name, &
                              'error opening message log file', &
                              failure )
        log_to_file = 0
      end if

    else

      log_to_file = 0

    end if


    !#--------------------------------------------------------------------------#
    !#                         -- output the message --                         #
    !#--------------------------------------------------------------------------#

    fmt_string = '( 1x, a, "(", a, ") : ", a )'

    log_message: if ( log_to_file == 0 ) then
      write( *, fmt = fmt_string ) &
                trim( routine_name ), &
                trim( state_descriptor( error_state_to_use ) ), &
                trim( message )
    else
      write( file_id, fmt = fmt_string ) &
                      trim( routine_name ), &
                      trim( state_descriptor( error_state_to_use ) ), &
                      trim( message )
      close( file_id )
    end if log_message

  end subroutine display_message

end module error_handler


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
! revision 1.1  2002/11/15 15:21:31  birk
! added to cvs mainly to see how this compiles on other platforms, it currently
! seems to compile on the ibm
!
! revision 1.3  2000/08/31 19:36:32  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.2  2000/08/24 15:27:18  paulv
! - the display_message subprogram was made recursive so it can call itself
!   if an error occurs opening the message log file defined by the optional
!   input argument message_log.
! - the message log file is now closed after the message is written (as it
!   should have always been...oops).
! - updated module and subprogram documentation.
!
! revision 1.1  2000/07/12 16:08:10  paulv
! initial checked in version
!
!
!
