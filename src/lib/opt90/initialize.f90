!------------------------------------------------------------------------------
!m+
! name:
!       initialize
!
! purpose:
!       module for rt model initialisation. the absorber space layer values
!       are calculated and the transmittance and spectral coefficient data
!       are read from their respective data files.
!
! category:
!       ncep rtm
!
! calling sequence:
!       use initialize
!
! outputs:
!       none.
!
! modules:
!       error_handler:               module to define error codes and handle
!                                    error conditions.
!
!       spectral_coefficients:       module containing the spectral coefficient
!                                    data and read methods.
!
!       transmittance_coefficients:  module containing the transmittance
!                                    coefficient data and read methods.
!
! contains:
!       initialize_rtm:   public function to initialise the rt model.
!
!       destroy_rtm:      public function to destroy the rt model space.
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       absorber space, transmittance coefficient, and spectral coefficient
!       data arrays are filled.
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
!       written by:     paul van delst, cimss@noaa/ncep 31-jul-2000
!                       paul.vandelst@ssec.wisc.edu
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

module initialize

  ! ----------
  ! module use
  ! ----------

  use error_handler
  use transmittance_coefficients
  use spectral_coefficients


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private
  public :: initialize_rtm
  public :: destroy_rtm


contains


!------------------------------------------------------------------------------
!s+
! name:
!       initialize_rtm
!
! purpose:
!       public function to initialise all the data arrays required by the
!       rt model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       result = initialize_rtm( tau_file      = tau_file,      &
!                                spectral_file = spectral_file, &
!                                path          = path,          &
!                                message_log   = message_log    )
!
! input arguments:
!       none.
!
! optional input arguments:
!       tau_file:      character string specifying a file name for the
!                      transmittance coefficient data file. if not
!                      specified, "transmittance_coefficients" is the
!                      default.
!                      units:      none
!                      type:       character
!                      dimension:  scalar
!                      attributes: intent( in ), optional
!
!       spectral_file: character string specifying a file name for the
!                      spectral coefficient data file. if not
!                      specified, "spectral_coefficients" is the
!                      default.
!                      units:      none
!                      type:       character
!                      dimension:  scalar
!                      attributes: intent( in ), optional
!
!       path:          character string specifying a file path for the
!                      input coefficient files. if not specified, the
!                      current directory is the default.
!                      units:      none
!                      type:       character
!                      dimension:  scalar
!                      attributes: intent( in ), optional
!
!       message_log:   character string specifying a filename in which any
!                      messages will be logged. if not specified, or if an
!                      error occurs opening the log file, the default action
!                      is to output messages to the screen.
!                      units:      none
!                      type:       character
!                      dimension:  scalar
!                      attributes: intent( in ), optional
!
! output arguments:
!       none.
!
! optional output arguments:
!       none.
!
! function result:
!       result = success => initialisation was successful
!              = failure => error occurred during initialisation
!
! calls:
!       display_message:             subroutine to output messages
!                                    source: error_handler module
!
!       read_tau_coefficients:       function to read the transmittance model
!                                    coefficients and predictor indices, and
!                                    absorber space definitions.
!                                    source: transmittance_coefficients module
!
!       read_spectral_coefficients:  function to read the spectral coefficients
!                                    for the various satellite/channels.
!                                    source: spectral_coefficients module
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       all public data arrays accessed by this module and its dependencies
!       are overwritten.
!
! restrictions:
!       if specified, the length of the combined path and filename strings
!       cannot exceed 255 characters.
!
!s-
!------------------------------------------------------------------------------

  function initialize_rtm( tau_file,      &
                           spectral_file, &
                           path,          &
                           message_log    ) &
                         result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    character( * ), optional, intent( in ) :: tau_file
    character( * ), optional, intent( in ) :: spectral_file
    character( * ), optional, intent( in ) :: path
    character( * ), optional, intent( in ) :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'initialize_rtm'


    ! ---------------
    ! local variables
    ! ---------------

    character( 255 ) :: transmittance_coefficient_file
    character( 255 ) :: spectral_coefficient_file


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic adjustl, &
              present, &
              trim



    !#--------------------------------------------------------------------------#
    !#           -- read the rt model transmittance coefficients --             #
    !#--------------------------------------------------------------------------#

    ! -- construct filename
    transmittance_coefficient_file = 'transmittance_coefficients'

    if ( present( tau_file ) ) &
      transmittance_coefficient_file = tau_file(1:len_trim(tau_file))

    if ( present( path ) ) &
      transmittance_coefficient_file = path(1:len_trim(path)) // &
                                       transmittance_coefficient_file


    ! -- read the data file
    error_status = read_tau_coefficients( &
         transmittance_coefficient_file(1:len_trim(transmittance_coefficient_file)) , &
         message_log = message_log        &
         )

    if ( error_status /= success ) then
      call display_message( routine_name, &
                            'error reading transmittance model coefficients', &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#               -- read the rt model spectral coefficients --              #
    !#--------------------------------------------------------------------------#

    ! -- construct filename
    spectral_coefficient_file = 'spectral_coefficients'

    if ( present( spectral_file ) ) &
      spectral_coefficient_file = spectral_file(1:len_trim(spectral_file))

    if ( present( path ) ) &
      spectral_coefficient_file = path(1:len_trim(path)) // &
                                  spectral_coefficient_file


    ! -- read the data file
    error_status = read_spectral_coefficients( &
          spectral_coefficient_file(1:len_trim(spectral_coefficient_file)), &
                                               message_log = message_log          )

    if ( error_status /= success ) then
      call display_message( routine_name, &
                            'error reading spectral coefficients', &
                            error_status, &
                            message_log = message_log )
      return
    end if


    !#--------------------------------------------------------------------------#
    !#                               -- done --                                 #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function initialize_rtm


!------------------------------------------------------------------------------
!s+
! name:
!       destroy_rtm
!
! purpose:
!       public function to destroy all the public allocated data arrays
!       creating during initialisation of the rt model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       result = destroy_rtm( message_log = message_log )
!
! input arguments:
!       none.
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
! output arguments:
!       none.
!
! optional output arguments:
!       none.
!
! function result:
!       result = success => deallocations were successful
!              = failure => error occurred during deallocations
!
! calls:
!       display_message:                subroutine to output messages
!                                       source: error_handler module
!
!       destroy_tau_coefficients:       function to deallocate the transmittance
!                                       model public arrays.
!                                       source: transmittance_coefficients module
!
!       destroy_spectral_coefficients:  function to deallocate the spectral
!                                       coefficient public arrays.
!                                       source: spectral_coefficients module
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       all public data arrays are deallocated.
!
! restrictions:
!       none.
!
!s-
!------------------------------------------------------------------------------

  function destroy_rtm( message_log ) result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    character( * ), optional, intent( in ) :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'destroy_rtm'



    !#--------------------------------------------------------------------------#
    !#      -- deallocate the rt model transmittance coefficient arrays --      #
    !#--------------------------------------------------------------------------#

    error_status = destroy_tau_coefficients( message_log = message_log )

    if ( error_status /= success ) then
      call display_message( routine_name, &
                            'error deallocating transmittance coefficients', &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#         -- deallocate the rt model spectral coefficient arrays --        #
    !#--------------------------------------------------------------------------#

    error_status = destroy_spectral_coefficients( message_log = message_log )

    if ( error_status /= success ) then
      call display_message( routine_name, &
                            'error deallocating spectral coefficients', &
                            error_status, &
                            message_log = message_log )
      return
    end if


    !#--------------------------------------------------------------------------#
    !#                               -- done --                                 #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function destroy_rtm

end module initialize


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
! revision 1.6  2001/08/31 21:17:25  paulv
! - add tau_file, spectral_file, and path optional arguments to initialize_rtm
!   to allow users to specify alternate file names and locations.
!
! revision 1.5  2001/08/16 16:40:56  paulv
! - updated documentation
!
! revision 1.4  2001/08/01 16:51:52  paulv
! - removed use of module absorber_space to reflect changes in code.
!   the absorber space levels are no longer calculated during model
!   initialisation, but are precalculated and stored in the transmittance
!   coefficient data file.
! - removed compute_absorber_space() function call for same reason.
!
! revision 1.3  2001/05/29 17:40:22  paulv
! - changed name of initialisation routine from rtm_initialize to initialize_rtm.
! - added destroy_rtm function to deallocate memory used in coefficient read.
!
! revision 1.2  2000/08/31 19:36:32  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.1  2000/08/08 16:39:35  paulv
! initial checkin
!
!
!
