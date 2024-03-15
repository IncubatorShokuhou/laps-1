!------------------------------------------------------------------------------
!m+
! name:
!       coefficient_utility
!
! purpose:
!       module to hold utility data and routines for reading and writing rt
!       model coefficient data
!
! category:
!       ncep rtm
!
! calling sequence:
!       use coefficient_utility
!
! outputs:
!       none.
!
! modules:
!       error_handler:          module to define error codes and handle error
!                               conditions
!
!       file_utility:   .       module containing global file utility routines
!
! contains:
!       open_coefficient_file:  public function to open the sequential access
!                               coefficient files.
!
!       check_coefficient_file: private function to determine if the coefficient
!                               file is in the correct format, endian-wise.
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
! comments:
!       all of the array documentation lists the dimensions by a single letter.
!       throughout the rtm code these are:
!         i: array dimension is of i predictors (istd and iint are variants).
!         j: array dimension is of j absorbing species.
!         k: array dimension is of k atmospheric layers.
!         l: array dimension is of l spectral channels.
!         m: array dimension is of m profiles.
!       not all of these dimensions will appear in every module.
!
! creation history:
!       written by:     paul van delst, cimss@noaa/ncep 12-jun-2000
!                       paul.vandelst@ssec.wisc.edu
!                         or
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

module coefficient_utility


  ! ----------
  ! module use
  ! ----------

  use type_kinds
  use file_utility
  use error_handler


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private
  public :: open_coefficient_file


  ! -----------------
  ! module parameters
  ! -----------------

  integer( long ), parameter, private :: magic_number = 123456789_long


  ! ----------------
  ! module variables
  ! ----------------

  character( 128 ) :: message


  ! -----------------
  ! module intrinsics
  ! -----------------

  intrinsic present, &
            trim



contains


!------------------------------------------------------------------------------
!s+
! name:
!       open_coefficient_file
!
! purpose:
!       public function to open the sequential access coefficient files
!
! calling sequence:
!       result = open_coefficient_file( coefficient_file, &
!                                       file_id, &
!                                       message_log = message_log )
!
! input arguments:
!       coefficient_file: name of the file containing the coefficient data.
!                         units:      none
!                         type:       character
!                         dimension:  scalar
!                         attributes: intent( in )
!
! optional input arguments:
!       message_log:      character string specifying a filename in which any
!                         messages will be logged. if not specified, or if an
!                         error occurs opening the log file, the default action
!                         is to output messages to the screen.
!                         units:      none
!                         type:       character
!                         dimension:  scalar
!                         attributes: intent( in ), optional
!
! output arguments:
!       file_id:          file logical unit number.
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( out )
!
! optional output arguments:
!       none.
!
! function result:
!       result = success => file open was successful
!              = failure => error occurred during file open or
!                           error occurred during file check
!
! calls:
!      file_exists:             function to determine if a named file exists.
!                               source: file_utility module
!
!      get_lun:                 function to return a free logical unit number
!                               for file access.
!                               source: file_utility module
!
!      display_message:         subroutine to output messages
!                               source: error_handler module
!
!      check_coefficient_file:  function to check the file for the correct
!                               endian-ness.
!                               source: private module subprogram
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
!       the file is inquired to determine if it exists. if so, it is opened
!       and the endian-ness is checked.
!s-
!------------------------------------------------------------------------------

  function open_coefficient_file( coefficient_file, &
                                  file_id,          &
                                  for_output,       &
                                  no_check,         &
                                  message_log )     &
                                result( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    character( * ),           intent( in )  :: coefficient_file
    integer,                  intent( out ) :: file_id

    integer,        optional, intent( in )  :: for_output
    integer,        optional, intent( in )  :: no_check

    character( * ), optional, intent( in )  :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'open_coefficient_file'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: file_check
    integer :: file_output
    integer :: io_status

    character( 7 ) :: file_status
    character( 5 ) :: file_action



    !#--------------------------------------------------------------------------#
    !#                      -- check optional arguments --                      #
    !#--------------------------------------------------------------------------#

    if ( present( no_check ) ) then
      file_check = 0
    else
      file_check = 1
    end if

    if ( present( for_output ) ) then
      file_output = 1
      file_check  = 0
    else
      file_output = 0
    end if


    !#--------------------------------------------------------------------------#
    !#                      -- check data file existence --                     #
    !#--------------------------------------------------------------------------#

    if ( file_output == 0 ) then

      ! -- if data file does not exist, return an error
      if ( .not. file_exists( coefficient_file ) ) then
        error_status = failure
        write( message, '( "coefficient file, ", a, " not found." )' ) &
                        trim( coefficient_file )
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        return
      end if

      ! -- set open keywords for reading
      file_status = 'old   '
      file_action = 'read '

    else

      ! -- if data file does exist, output a warning message
      if ( file_exists( coefficient_file ) ) then
        write( message, '( "coefficient file, ", a, " will be overwritten." )' ) &
                        trim( coefficient_file )
        call display_message( routine_name, &
                              trim( message ), &
                              warning, &
                              message_log = message_log )
      end if

      ! -- set open keywords for writing
      file_status = 'replace'
      file_action = 'write'

    end if



    !#--------------------------------------------------------------------------#
    !#                        -- open the data file --                          #
    !#--------------------------------------------------------------------------#

    file_id = get_lun()
    open( file_id, file   = coefficient_file, &
                   status = trim( file_status ), &
                   action = trim( file_action ), &
                   access = 'sequential', &
                   form   = 'unformatted', &
                   iostat = io_status )

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error opening ", a, ". iostat = ", i5 )' ) &
                      coefficient_file, io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#        -- if file is opened for output, write the magic number --        #
    !#--------------------------------------------------------------------------#

    if ( file_output == 1 ) then

      write( file_id, iostat = io_status ) magic_number

      if ( io_status /= 0 ) then
        error_status = failure
        write( message, '( "error writing magic number to ", a, ". iostat = ", i5 )' ) &
                        trim( coefficient_file ), io_status
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        close( file_id )
        return
      end if

    end if


    
    !#--------------------------------------------------------------------------#
    !#             -- check the coefficient data file if required --            #
    !#--------------------------------------------------------------------------#

    if ( file_check == 0 ) then
      error_status = success
      return
    end if

    error_status = check_coefficient_file( file_id, &
                                           message_log = message_log )

    if ( error_status /= success ) then
      close( file_id )
      write( message, '( "error checking ", a, ". file closed." )' ) &
                      coefficient_file
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      return
    end if

  end function open_coefficient_file



!------------------------------------------------------------------------------
!p+
! name:
!       check_coefficient_file
!
! purpose:
!       private function to determine if the coefficient file is in the correct
!       format, endian-wise.
!
! calling sequence:
!       result = check_coefficient_file( file_id, &
!                                        message_log = message_log )
!
! input arguments:
!       file_id:          file logical unit number for the open file that is
!                         to be checked.
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in )
!
!
! optional input arguments:
!       message_log:      character string specifying a filename in which any
!                         messages will be logged. if not specified, or if an
!                         error occurs opening the log file, the default action
!                         is to output messages to the screen.
!                         units:      none
!                         type:       character
!                         dimension:  scalar
!                         attributes: intent( in ), optional
!
! output arguments:
!       none.
!
! optional output arguments:
!       none.
!
! function result:
!       result = success => file check was successful
!              = failure => error occurred reading a file record or
!                           8- and/or 32-bit integers not supported.
!
! calls:
!       display_message:  subroutine to output messages
!                         source: error_handler module
!
! modules:
!       type_kinds:  module containing data type kind definitions. only the
!                    byte and long kind types, and the definition of the number
!                    of bytes used for each type are used.
!
! contains:
!       swap_endian_long_integer:  function to byte-swap a long integer to
!                                  determine if input data byte-swapping may
!                                  be required.
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
!       the file is inquired to determine if it exists. if so, it is opened
!       and the endian-ness is checked.
!p-
!------------------------------------------------------------------------------

  function check_coefficient_file( file_id,  &
                                   message_log ) &
                                 result( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    integer,        intent( in )           :: file_id
    character( * ), intent( in ), optional :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ),  parameter :: routine_name = 'check_coefficient_file'


    ! ---------------
    ! local variables
    ! ---------------

    integer         :: io_status
    integer( long ) :: magic_number_read, dummy_long



    !#--------------------------------------------------------------------------#
    !#             -- check that the current compilation supports --            #
    !#             -- 1- and 4-byte integer types                 --            #
    !#--------------------------------------------------------------------------#

    if ( bit_size( 1_long ) /= 32 .or. &
         bit_size( 1_byte ) /=  8      ) then
      error_status = failure
      call display_message( routine_name, &
                            '8- and/or 32-bit integers not supported. '//&
                            'unable to determine endian-ness', &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#                        -- read the magic number --                       #
    !#--------------------------------------------------------------------------#

    read( file_id, iostat = io_status ) magic_number_read

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading file. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            message, &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#                      -- compare the magic numbers --                     #
    !#--------------------------------------------------------------------------#

    if ( magic_number_read /= magic_number ) then

      ! -- set the error status
      error_status = failure

      ! -- byte swap the file data.
      dummy_long = swap_endian_long_integer( magic_number_read )

      ! -- check the file data again
      if ( dummy_long /= magic_number ) then
        call display_message( routine_name, &
                              'unrecognised file format. invalid magic number.', &
                              error_status, &
                              message_log = message_log )
        return
      end if

      ! -- if we get here then the data does need to be byte-swapped.
      call display_message( routine_name, &
                            'data file needs to be byte-swapped.', &
                            error_status, &
                            message_log = message_log )
      return
    end if

    error_status = success


  contains


    !#--------------------------------------------------------------------------#
    !#          -- internal subprogram to byte-swap a long integer --           #
    !#--------------------------------------------------------------------------#

    function swap_endian_long_integer ( input ) &
                                      result ( output )


      ! -----------------
      ! type declarations
      ! -----------------

      ! -- argument and result
      integer( long ), intent( in ) :: input
      integer( long )               :: output


      ! -- local variables
      integer( byte ), dimension( n_bytes_for_long_kind ) :: binput, boutput
      integer( long )                                     :: linput, loutput
      integer :: i


      ! -------------------------------------------
      ! equivalence the byte array and long integer
      ! -------------------------------------------

      equivalence ( binput,  linput  ), &
                  ( boutput, loutput )


      ! ----------------------------------------------------------------
      ! loop over the number of bytes for swapping.
      !
      ! doing it this way is little bit faster (by about a factor of 4)
      ! than using the mvbits intrinsic ( on the systems tested; linux
      ! and aix):
      !
      !  call mvbits( input, 0,  8, output, 24 )  ! bits  0-7  --> 24-31
      !  call mvbits( input, 8,  8, output, 16 )  ! bits  8-15 --> 16-23
      !  call mvbits( input, 16, 8, output,  8 )  ! bits 16-23 -->  8-15
      !  call mvbits( input, 24, 8, output,  0 )  ! bits 24-31 -->  0-8
      !
      ! but only if the byte swap loop is inline (rather than by calling
      ! a generic byte swap routine with the number of bytes to swap is
      ! passed as an argument.)
      ! ----------------------------------------------------------------

      ! -- reassign the input argument. can't
      ! -- equivalence dummy arguments.
      linput = input

      ! -- loop over the bytes and swap
      do i = 1, n_bytes_for_long_kind
        boutput( i ) = binput( n_bytes_for_long_kind - ( i - 1 ) )
      end do

      ! -- assign the output argument. can't
      ! -- equivalence dummy arguments.
      output = loutput

    end function swap_endian_long_integer

  end function check_coefficient_file

end module coefficient_utility


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
! revision 1.8  2001/08/16 16:39:54  paulv
! - updated documentation
!
! revision 1.7  2001/08/09 20:38:15  paulv
! - changed magic number visibility attribute to private. ahhh....
! - moved all the spectral and transmittance coefficient data type and name
!   definitions into their respective modules. another ahhh.....
! - added optional for_output argument to open_coefficient_file function
!   so that the same function can be used to open coefficient files for
!   writing. it also means the magic number write can be done in this module
!   totally encapsulating that functionality in this module only. double ahhh....
!
! revision 1.6  2001/08/01 16:43:05  paulv
! - updated the definitions of data items and types in the transmittance
!   coefficient data file to reflect changes in code. the absorber space
!   levels are no longer calculated during model initialisation, but are
!   precalculated and stored in the transmittance coefficient data file.
!
! revision 1.5  2001/07/12 16:58:18  paulv
! - added use of type_kinds module at top of this module. previously it was
!   used only in the check_coefficient_file() function.
! - data file magic number definition now defined at top of module rather
!   than in the check_coefficient_file() function.
! - definitions for the number, type, and names of items in the transmittance
!   and spectral coefficient files moved from the transmittance_coefficients
!   and spectral coefficients module to this one. this was done to allow this
!   module to be used in both reading and writing/reformatting the coefficient
!   data files.
! - module-wide error message character string defined.
! - added no_check optional argument to the open_coefficient_file() function.
!   this was done to allow the function to be used to open the old format
!   coefficient files for reformatting by not performing a magic number check.
!
! revision 1.4  2000/08/31 19:36:31  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.3  2000/08/24 15:22:10  paulv
! - file access changed from direct to sequential. record length argument
!   no longer required by open_coefficient_file and check_coefficient_file
!   subprograms.
! - inquire statement in open_coefficient_file that checks for existence
!   of the file replaced by function file_exists in module file_utility.
! - check_coefficient_file used to return a warning status if either 8- or
!   32-bit integers were not supported. this condition now returns a
!   failure status as the magic number would not be read so any subsequent
!   attempt to read data would either fail or return junk.
! - the name of the swap_endian_fourbyte_integer subprogram was changed to
!   swap_endian_long_integer to remove any indication of how many bytes are
!   expected for this data type *apart* from the definition of
!   n_bytes_for_long_kind in the type_kinds module.
! - updated module and subprogram documentation.
!
! revision 1.2  2000/08/08 17:05:45  paulv
! cosmetic changes to highlight parameters in the source by making them
! uppercase.
!
! revision 1.1  2000/07/12 16:08:10  paulv
! initial checked in version
!
!
!

