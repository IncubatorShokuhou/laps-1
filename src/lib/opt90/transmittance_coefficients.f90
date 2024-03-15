!------------------------------------------------------------------------------
!m+
! name:
!       transmittance_coefficients
!
! purpose:
!       module to hold rt model transmittance coefficients 
!
! category:
!       ncep rtm
!
! calling sequence:
!       use transmittance_coefficients
!
! outputs:
!       absorber_space_levels:     array containing teh absorber space levels
!                                  used in generating the transmittance
!                                  coefficients.
!                                  units:      depends on absorber type
!                                  type:       double
!                                  dimension:  0:ka x j
!                                  attributes: save, public
!
!       predictor_index:           array containing the predictor indices used to
!                                  identify predictors in the transmittance model.
!                                  units:      none.
!                                  type:       integer
!                                  dimension:  0:i x l x j
!                                  attributes: save, public
!
!       tau_coefficients:          array containing the transmittance model
!                                  coefficients.
!                                  units:      absorber, predictor dependent
!                                  type:       double
!                                  dimension:  0:i x 0:ka x l x j
!                                  attributes: save, public
!
! modules:
!       type_kinds:                module containing data type kind definitions.
!
!       error_handler:             module to define error codes and handle error
!                                  conditions
!
!       parameters:                module containing parameter definitions for the
!                                  rt model.
!
!       coefficient_utility:       module containing coefficient file utility
!                                  subprograms.
!
! contains:
!       read_tau_coefficients:     public function to read the transmittance model
!                                  coefficients and predictor indices and fill the
!                                  public data arrays.
!
!       destroy_tau_coefficients:  public function to release the memory that
!                                  was allocated for the transmittance coefficient
!                                  data.
!
!       write_tau_coefficients:    public function to write the new format
!                                  transmittance coefficient data file.
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

module transmittance_coefficients


  ! ---------------------
  ! module use statements
  ! ---------------------

  use type_kinds
  use error_handler
  use parameters
  use coefficient_utility


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private
  public :: read_tau_coefficients
  public :: destroy_tau_coefficients
  public :: write_tau_coefficients



  ! -----------------
  ! module parameters
  ! -----------------

  ! -- transmittance coefficient file version numbers
  integer( long ), parameter, private :: min_valid_release = 1_long
  integer( long ), parameter, private :: max_valid_release = 1_long

  integer( long ), parameter, private :: min_valid_version = 4_long
  integer( long ), parameter, private :: max_valid_version = 4_long

  ! -- transmittance coefficient file descriptors:
  ! -- number of items in the transmittance coefficient file
  integer( long ), parameter, private :: n_transmittance_items = 4_long

  ! -- data types of the transmittance coefficient data
  !    5 = double (i.e. 8-byte float)
  !    4 = single (i.e. 4-byte float)
  !    3 = long   (i.e. 4-byte integer)
  integer( long ), parameter, private, &
                   dimension( n_transmittance_items ) :: transmittance_data_type = &
                                                         (/ 5_long, 5_long, 3_long, 5_long /)

  ! -- names of the data items (for error processing)
  character( * ), parameter, private, &
                  dimension( n_transmittance_items ) :: transmittance_data_name = &
                                                        (/ 'alpha                    ', &
                                                           'absorber_space_levels    ', &
                                                           'predictor_index          ', &
                                                           'transmittance_coefficient' /)

  ! ----------------
  ! module variables
  ! ----------------

  character( 128 ) :: message


  ! ---------------------------------------------------
  ! definitions of shared data.
  !
  ! note that the save attribute is specified to ensure
  ! that the data is retained even when this module is
  ! not being directly accessed.
  ! ---------------------------------------------------

  ! -- alpha value
  !
  ! dimension is: number of absorbers
  real( double ), save, &
                  public, &
                  allocatable, &
                  dimension( : ) :: alpha


  ! -- absorber space levels
  !
  ! eventual dimensions are:
  !   #1: 0->maximum number of absorber layers  (fixed by transmittance algorithm)
  !   #2: 1->number of absorbers                (fixed by transmittance algorithm)
  !
  ! although both dimensions are fixed by the transmittance
  ! algorithm, there is the possibility of two different sets
  ! of transmittance coefficients existing: when different
  ! profile sets are used to generate the coefficients. by
  ! including the data with the transmittance coefficients
  ! rather than calculating them every time the model is
  ! initialised using different coeficient sets is a simple
  ! matter of reading a different file.
  real( double ), save, &
                  public, &
                  allocatable, &
                  dimension( :, : ) :: absorber_space_levels



  ! -- absorber predictor index arrays
  !
  ! remember the 0'th index is an indicator of absorption for
  ! a particular channel. if, for any channel, l, 
  ! predictor_index( 0,l ) = 0 then there is no absorption due
  ! to the absorber in channel l.
  !
  ! eventual dimensions are:
  !   #1: 0->maximum number of predictors used - (fixed by transmittance algorithm)
  !   #2: 1->number of channels requested      - (user selectable)
  !   #3: 1->number of absorbers               - (fixed by transmittance algorithm)
  integer( long ), save, &
                   public, &
                   allocatable, &
                   dimension( :, :, : ) :: predictor_index



  ! -- coefficients for calculating transmittances
  !
  ! the 0'th predictor coefficient is the offset term used in:
  !
  !                                  n
  ! absorption_coefficient = b(0) + sum{ b(i)*x(i) }
  !                                 i=1
  !
  ! where n = n_predictors_used
  !
  ! eventual dimensions are:
  !   #1: 0->maximum number of predictors used - (fixed by transmittance algorithm)
  !   #2: 0->number of absorber layers         - (fixed by transmittance algorithm)
  !   #3: 1->number of channels requested      - (user selectable)
  !   #4: 1->number of absorbers               - (fixed by transmittance algorithm)
  real( double ), save, &
                  public, &
                  allocatable, &
                  dimension( :, :, :, : ) :: tau_coefficients



contains


!------------------------------------------------------------------------------
!s+
! name:
!       read_tau_coefficients
!
! purpose:
!       public function to read the transmittance model coefficients and
!       predictor indices and fill the public data arrays.
!
! calling sequence:
!       result = read_tau_coefficients( coefficient_file, &
!                                       message_log = message_log )
!
! input arguments:
!       coefficient_file: name of the file containing the coefficient data.
!                         units:      none
!                         type:       character
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
! optional outupt arguments:
!       none.
!
! function result:
!       result = success => coefficient read was successful
!              = failure => error occurred opening or accessing coefficient
!                           data file
!
! calls:
!      display_message:         subroutine to output messages
!                               source: error_handler module
!
!      open_coefficient_file:   function to open the transmittance coefficient
!                               data file.
!                               source: coefficient_utility module
!
!      get_max_n_channels:      routine to retrieve the value of the
!                               max_n_channels "pseudo-parameter".
!                               source: parameters module
!
!      set_max_n_channels:      routine to set the value of the
!                               max_n_channels "pseudo-parameter".
!                               source: parameters module
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
!s-
!------------------------------------------------------------------------------

  function read_tau_coefficients( coefficient_file, &
                                  message_log ) &
                                result ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    character( * ), intent( in )           :: coefficient_file
    character( * ), intent( in ), optional :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'read_tau_coefficients'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- these are dimensioned ( long ) as they are read in.
    integer( long ) :: file_release
    integer( long ) :: file_version

    integer( long ) :: n_predictors_to_use
    integer( long ) :: n_absorber_layers
    integer( long ) :: n_channels
    integer( long ) :: n_absorbers
    integer( long ) :: n_items
    integer, dimension( n_transmittance_items ) :: data_type

    integer :: i, j, k, l
    integer :: io_status
    integer :: file_id
    integer :: allocate_status

    logical :: min_relver
    logical :: max_relver

    ! -- maximum channels pseudo parameter
    integer :: max_n_channels
    logical :: is_set


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic trim



    !#--------------------------------------------------------------------------#
    !#           -- open the transmittance coefficient data file --             #
    !#--------------------------------------------------------------------------#

    error_status = open_coefficient_file( coefficient_file, &
                                          file_id, &
                                          message_log = message_log )

    if ( error_status /= success ) then
      call display_message( routine_name, &
                            'error opening '//trim( coefficient_file ), &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#                      -- read the "file header" --                        #
    !#--------------------------------------------------------------------------#

    ! ------------------------------------
    ! read the release/version information
    ! ------------------------------------

    read( file_id, iostat = io_status ) file_release, file_version

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading file release/version information. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if

    ! -- check that this file is not an old release/version
    if ( ( file_release + file_version ) < ( min_valid_release + min_valid_version ) ) then
      error_status = failure
      write( message, '( "need to update the coefficient file. ", &
                        &"file version is ", i1, ".", i2.2, &
                        &". oldest valid release is ",i1,".",i2.2,"." )' ) &
                      file_release,  file_version, &
                      min_valid_release, min_valid_version
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    ! -- check that this file is not a too new release/version
    ! -- i.e. update the software!
    if ( ( file_release + file_version ) > ( max_valid_release + max_valid_version ) ) then
      error_status = failure
      write( message, '( "need to update the coefficient read software. ", &
                        &"file version is ", i1, ".", i2.2, &
                        &". newest valid release is ",i1,".",i2.2,"." )' ) &
                      file_release,  file_version, &
                      max_valid_release, max_valid_version
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    ! -------------------
    ! read the dimensions
    ! -------------------

    read( file_id, iostat = io_status ) n_predictors_to_use, &
                                        n_absorber_layers, &
                                        n_channels, &
                                        n_absorbers

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading data dimensions. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    ! -----------------------------        
    ! read the number of data items
    ! -----------------------------

    read( file_id, iostat = io_status ) n_items

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading total number of data items. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if

    if ( n_items /= n_transmittance_items ) then
      error_status = failure
      call display_message( routine_name, &
                            'number of data items in '//&
                            trim( coefficient_file )//&
                            ' inconsistent with definition.', &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    ! -------------------
    ! read the data types
    ! -------------------

    read( file_id, iostat = io_status ) data_type
    
    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading data types. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    do l = 1, n_items

      if ( data_type( l ) /= transmittance_data_type( l ) ) then
        error_status = failure
        write( message, '( "invalid type for ", a, " data item." )' ) &
                        trim( transmittance_data_name( l ) )
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        close( file_id )
        return
      end if

    end do


    !#--------------------------------------------------------------------------#
    !#                     -- check the dimension values --                     #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------------------------------
    ! check the dimensions that are transmittance algorithm dependent:
    !   n_predictors_to_use, n_absorber_layers, n_absorbers
    ! ----------------------------------------------------------------

    ! -- number of predictors to use
    if ( n_predictors_to_use /= max_n_predictors_used ) then
      error_status = failure
      write( message, '( "predictor index dimension inconsistent. value is ", i3, &
                        &" but should be ", i1 )' ) &
                      n_predictors_to_use, max_n_predictors_used
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if

    ! -- number of absorber layers
    if ( n_absorber_layers /= max_n_absorber_layers ) then
      error_status = failure
      write( message, '( "absorber layer dimension inconsistent. value is ", i3, &
                        &" but should be ", i3 )' ) &
                      n_absorber_layers, max_n_absorber_layers
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if

    ! -- number of absorbers
    if ( n_absorbers /= max_n_absorbers ) then
      error_status = failure
      write( message, '( "absorber dimension inconsistent. value is ", i3, &
                        &" but should be ", i1 )' ) &
                      n_absorbers, max_n_absorbers
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    ! --------------------------------------------
    ! set the global definition for max_n_channels
    ! --------------------------------------------

    call get_max_n_channels( max_n_channels, is_set )

    if ( is_set  ) then
      if ( max_n_channels /= n_channels ) then
        error_status = warning
        write( message, '( "max_n_channels set to different value, ", i4, ", ", &
                          &"than defined in coefficient file, ", i4, ". overwriting." )' ) &
                        max_n_channels, n_channels
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        call set_max_n_channels( n_channels )
      end if
    else
      call set_max_n_channels( n_channels )
    end if



    !#--------------------------------------------------------------------------#
    !#          -- allocate arrays for transmittance coefficient data --        #
    !#--------------------------------------------------------------------------#

    ! -- check if arrays are already allocated
    if ( allocated( alpha                 ) .or. &
         allocated( absorber_space_levels ) .or. &
         allocated( predictor_index       ) .or. &
         allocated( tau_coefficients      )      ) then
      error_status = failure
      call display_message( routine_name, &
                            'transmittance coefficient data arrays already allocated.', &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if

    ! -- if not, allocate them
    allocate( alpha( n_absorbers ), &

              absorber_space_levels( 0:n_absorber_layers, &
                                     n_absorbers          ), &

              predictor_index( 0:n_predictors_to_use, &
                               n_channels,            &
                               n_absorbers            ), &

              tau_coefficients( 0:n_predictors_to_use, &
                                0:n_absorber_layers,   &
                                n_channels,            &
                                n_absorbers            ), &
              stat = allocate_status )

    if ( allocate_status /= 0 ) then
      error_status = failure
      call display_message( routine_name, &
                            'unable to allocate arrays for transmittance coefficient data.', &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if
   


    !#--------------------------------------------------------------------------#
    !#                 -- read the absorber space data --                       #
    !#--------------------------------------------------------------------------#

    ! -- alpha values
    read( file_id, iostat = io_status ) alpha

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading alpha values. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      return
    end if

    ! -- absorber space levels
    read( file_id, iostat = io_status ) absorber_space_levels

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading absorber space levels. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#                        -- loop over channels --                          #
    !#--------------------------------------------------------------------------#

    l_channel_loop: do l = 1, n_channels


      ! ----------------------
      ! read predictor indices
      ! ----------------------

      read( file_id, iostat = io_status ) predictor_index( :, l, : )

      if ( io_status /= 0 ) then
        error_status = failure
        write( message, '( "error reading channel ", i4, " predictor indices. iostat = ", i5 )' ) &
                        l, io_status
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        return
      end if


      ! -------------------------------
      ! read transmittance coefficients
      ! -------------------------------

      read( file_id, iostat = io_status ) tau_coefficients( :, :, l, : )

      if ( io_status /= 0 ) then
        error_status = failure
        write( message, '( "error reading channel ", i4, " coefficients. iostat = ", i5 )' ) &
                        l, io_status
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        return
      end if

    end do l_channel_loop



    !#--------------------------------------------------------------------------#
    !#                          -- close the file --                            #
    !#--------------------------------------------------------------------------#

    close( file_id, iostat = io_status )

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error closing ", a, ". iostat = ", i5 )' ) &
                      coefficient_file, io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      return
    end if


    !#--------------------------------------------------------------------------#
    !#                      -- successful completion --                         #
    !#--------------------------------------------------------------------------#

    ! -- output an info message
!    write( message, '( "file version: ", i1, ".", i2.2, 2x, &
!                      &"n_predictors_to_use=",i1,2x,&
!                      &"n_absorber_layers=",i3,2x,&
!                      &"n_channels=",i4,2x,&
!                      &"n_absorbers=",i1 )' ) &
!                    file_release, file_version, &
!                    n_predictors_to_use, &
!                    n_absorber_layers, &
!                    n_channels, &
!                    n_absorbers
!    call display_message( routine_name, &
!                          trim( message ), &
!                          information, &
!                          message_log = message_log )

    error_status = success

  end function read_tau_coefficients 



!------------------------------------------------------------------------------
!s+
! name:
!       destroy_tau_coefficients
!
! purpose:
!       public function to deallocate the alpha, absorber space level, 
!       predictor index and transmittance model coefficient public data arrays.
!
! calling sequence:
!       result = destroy_tau_coefficients( message_log = message_log )
!
! input arguments:
!       none.
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
! optional outupt arguments:
!       none.
!
! function result:
!       result = success => deallocation was successful
!              = failure => error occurred deallocating coefficient arrays.
!
! calls:
!      display_message:         subroutine to output messages
!                               source: error_handler module
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
!s-
!------------------------------------------------------------------------------

  function destroy_tau_coefficients( message_log ) &
                                   result ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    character( * ), intent( in ), optional :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'destroy_tau_coefficients'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: allocate_status



    !#--------------------------------------------------------------------------#
    !#         -- deallocate arrays for transmittance coefficient data --       #
    !#--------------------------------------------------------------------------#

    deallocate( alpha,                 &
                absorber_space_levels, &
                predictor_index,       &
                tau_coefficients,      &
                stat = allocate_status )

    if ( allocate_status /= 0 ) then
      error_status = failure
      call display_message( routine_name, &
                            'error occurred deallocating transmittance coefficient data arrays.', &
                            error_status, &
                            message_log = message_log )
      return
    end if
   


    !#--------------------------------------------------------------------------#
    !#           -- reset the max_n_channels "pseudo-parameter" --              #
    !#--------------------------------------------------------------------------#

    call reset_max_n_channels()


    !#--------------------------------------------------------------------------#
    !#                      -- successful completion --                         #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function destroy_tau_coefficients 





!------------------------------------------------------------------------------
!
! name:
!       write_tau_coefficients
!
! purpose:
!       function to write the transmittance coefficient data to file.
!
! category:
!       ncep rtm
!
! language:
!       fortran-90
!
! calling sequence:
!       result = write_tau_coefficients( coefficient_file,           &  ! input
!                                        alpha,                      &  ! input
!                                        absorber_space_levels,      &  ! input
!                                        predictor_index,            &  ! input, 0:iuse        x l x j
!                                        transmittance_coefficients, &  ! input, 0:iuse x 0:ka x l x j
!                                        message_log                 )  ! optional input
!
! input arguments:
!       coefficient_file:            name of the file to which the transmittance
!                                    coefficient data is to be written.
!                                    units:      none
!                                    type:       character
!                                    dimension:  scalar
!                                    attributes: intent( in )
!
!       alpha:                       vector containing the exponential parameter
!                                    for each absorber used to construct the
!                                    absorber space.
!                                    ** note: this variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!                                    units:      none
!                                    type:       double
!                                    dimension:  j
!                                    attributes: intent( in )
!
!       absorber_space_levels:       array containing the absorber space levels
!                                    used in the transmittance model.
!                                    ** note: this variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!                                    units:      varies with absorber
!                                    type:       double
!                                    dimension:  0:ka x j
!                                    attributes: intent( in )
!
!       predictor_index:             array containing the predictor indices
!                                    used to identify which predictors to use
!                                    in the transmittance model.
!                                    ** note: this variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!
!                                    remember the 0'th index is an indicator
!                                    of absorption for a particular channel.
!                                    if, for any channel, l,
!                                      predictor_index( 0,l ) = 0
!                                    then there is no absorption due to that
!                                    absorber in channel l. if the 0'th index
!                                    is non-zero, then it is the number of
!                                    valid predictor indices to follow (may
!                                    not always need all the predictors).
!                                    units:      none
!                                    type:       long
!                                    dimension:  0:iuse x l x j
!                                    attributes: intent( in )
!
!       transmittance_coefficients:  array containing the transmittance
!                                    model coefficients in absorber, and
!                                    predictor dependent units.
!                                    ** note: this variable is in a different scope
!                                             from the module data variable of the
!                                             same name.
!
!                                    the 0'th predictor coefficient is the
!                                    offset term used in:
!                                                         __ n
!                                                        \
!                                      abs_coeff = b(0) + >  b(i)*x(i)
!                                                        /__
!                                                           i=1
!
!                                    where n = n_predictors_used (or the 0'th
!                                    term of predictor_index if it is non-zero).
!                                    units:      none
!                                    type:       double
!                                    dimension:  0:iuse x 0:ka x l x j
!                                    attributes: intent( in )
!
! optional input arguments:
!       message_log:                 character string specifying a filename in which
!                                    any messages will be logged. if not specified,
!                                    or if an error occurs opening the log file,
!                                    the default action is to output messages to
!                                    the screen.
!                                    units:      none
!                                    type:       character
!                                    dimension:  scalar
!                                    attributes: intent( in ), optional
!
! output arguments:
!       none.
!
! optional output arguments:
!       none.
!
! function result:
!       result = success => data file write was successful
!              = failure => error occurred opening or writing to the
!                           data file
!
! calls:
!       display_message:        subroutine to output messages
!                               source: error_handler module
!
! contains:
!       none.
!
! include files:
!       none.
!
! externals:
!       none.
!
! common blocks:
!       none.
!
! files accessed:
!       none.
!
! side effects:
!       if a coefficient file with the specified name exists, it is
!       overwritten.
!
! restrictions:
!       none.
!
!------------------------------------------------------------------------------

  function write_tau_coefficients( coefficient_file,           &  ! input
                                   alpha,                      &  ! input
                                   absorber_space_levels,      &  ! input
                                   predictor_index,            &  ! input
                                   transmittance_coefficients, &  ! input

                                   release,                    &  ! optional input
                                   version,                    &  ! optional input
                                   message_log )               &  ! optional input
                                 result ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- input
    character( * ),                             intent( in ) :: coefficient_file            ! input
    real( double ),  dimension(            : ), intent( in ) :: alpha                       ! input, j
    real( double ),  dimension(     0:,    : ), intent( in ) :: absorber_space_levels       ! input, 0:ka x j
    integer( long ), dimension( 0:,     :, : ), intent( in ) :: predictor_index             ! input, 0:iuse x l x j
    real( double ),  dimension( 0:, 0:, :, : ), intent( in ) :: transmittance_coefficients  ! input, 0:iuse x 0:ka x l x j
    !                            ^   ^  ^  ^
    !                            i   k  l  j

    integer,        optional,                   intent( in ) :: release
    integer,        optional,                   intent( in ) :: version

    character( * ), optional,                   intent( in ) :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ),  parameter :: routine_name = 'write_tau_coefficients'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: n_alpha_absorbers 

    integer :: n_as_layers
    integer :: n_as_absorbers 

    integer :: n_pi_predictors
    integer :: n_pi_channels  
    integer :: n_pi_absorbers 

    ! -- these are typed ( long ) as they will be output
    integer( long ) :: n_tc_predictors
    integer( long ) :: n_tc_layers    
    integer( long ) :: n_tc_channels  
    integer( long ) :: n_tc_absorbers 

    integer( long ) :: file_release
    integer( long ) :: file_version

    integer :: file_id
    integer :: io_status

    integer :: i, j, k, l


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic size, &
              trim



    !#--------------------------------------------------------------------------#
    !#                           -- check input --                              #
    !#--------------------------------------------------------------------------#

    ! -------------------------------------------
    ! check the predictor index and transmittance
    ! coefficient input arrays
    ! -------------------------------------------

    n_alpha_absorbers = size( alpha )

    n_as_layers       = size( absorber_space_levels, dim = 1 ) - 1
    n_as_absorbers    = size( absorber_space_levels, dim = 2 )

    n_pi_predictors   = size( predictor_index, dim = 1 ) - 1
    n_pi_channels     = size( predictor_index, dim = 2 )
    n_pi_absorbers    = size( predictor_index, dim = 3 )

    n_tc_predictors   = size( transmittance_coefficients, dim = 1 ) - 1
    n_tc_layers       = size( transmittance_coefficients, dim = 2 ) - 1
    n_tc_channels     = size( transmittance_coefficients, dim = 3 )
    n_tc_absorbers    = size( transmittance_coefficients, dim = 4 )

    ! -- check number of predictors
    if ( n_pi_predictors /= n_tc_predictors ) then 
      error_status = failure
      call display_message( routine_name, &
                            'inconsistent predictor_index and '//&
                            'transmittance_coefficients predictor dimension.', &
                            error_status, &
                            message_log = message_log )
      return
    end if

    ! -- check number of layers
    if ( n_as_layers /= n_tc_layers ) then 
      error_status = failure
      call display_message( routine_name, &
                            'inconsistent absorber_space_levels and '//&
                            'transmittance_coefficients level dimension.', &
                            error_status, &
                            message_log = message_log )
      return
    end if

    ! -- check number of channels
    if ( n_pi_channels /= n_tc_channels ) then 
      error_status = failure
      call display_message( routine_name, &
                            'inconsistent predictor_index and '//&
                            'transmittance_coefficients channel dimension.', &
                            error_status, &
                            message_log = message_log )
      return
    end if

    ! -- check number of absorbers
    if ( n_alpha_absorbers /= n_tc_absorbers .or. &
         n_as_absorbers    /= n_tc_absorbers .or. &
         n_pi_absorbers    /= n_tc_absorbers      ) then 
      error_status = failure
      call display_message( routine_name, &
                            'inconsistent alpha/absorber_space_levels/predictor_index '//&
                            'and transmittance_coefficients absorber dimension.', &
                            error_status, &
                            message_log = message_log )
      return
    end if


    ! ----------------------------------
    ! assign release and version numbers
    ! ----------------------------------

    if ( present( release ) ) then
      file_release = int( release, long )
    else
      file_release = max_valid_release
    end if

    if ( present( version ) ) then
      file_version = int( version, long )
    else
      file_version = max_valid_version
    end if



    !#--------------------------------------------------------------------------#
    !#      -- open the transmittance coefficient data file for output --       #
    !#--------------------------------------------------------------------------#

    error_status = open_coefficient_file( coefficient_file,         &
                                          file_id,                  &
                                          for_output  = 1,          &
                                          message_log = message_log )

    if ( error_status /= success ) then
      call display_message( routine_name, &
                            'error opening '//trim( coefficient_file )//' for output.', &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#                     -- write the "file header" --                        #
    !#--------------------------------------------------------------------------#

    ! -------------------------------------
    ! write the release/version information
    ! -------------------------------------

    write( file_id, iostat = io_status ) file_release, file_version 

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error writing file release/version information to ", a, &
                        &". iostat = ", i5 )' ) &
                      trim( coefficient_file ), io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    ! --------------------
    ! write the dimensions
    ! --------------------

    write( file_id, iostat = io_status ) n_tc_predictors, &
                                         n_tc_layers,     &
                                         n_tc_channels,   &
                                         n_tc_absorbers

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error writing dimension values to ", a, ". iostat = ", i5 )' ) &
                      trim( coefficient_file ), io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    ! ----------------------------------------------
    ! write the number of data items per channel and
    ! the type of the data items.
    ! ----------------------------------------------

    write( file_id, iostat = io_status ) n_transmittance_items

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error writing data item count to ", a, ". iostat = ", i5 )' ) &
                      trim( coefficient_file ), io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    write( file_id, iostat = io_status ) transmittance_data_type

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error writing data item types to ", a, ". iostat = ", i5 )' ) &
                      trim( coefficient_file ), io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if


    !#--------------------------------------------------------------------------#
    !#                 -- write the absorber space data --                      #
    !#--------------------------------------------------------------------------#

    ! -- alpha values
    write( file_id, iostat = io_status ) alpha

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error writing alpha values. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if

    ! -- absorber space level data
    write( file_id, iostat = io_status ) absorber_space_levels

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error writing absorber space level data. iostat = ", i5 )' ) &
                      io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#                     -- write the data by channel --                      #
    !#--------------------------------------------------------------------------#

    do l = 1, n_tc_channels


      ! ---------------
      ! predictor index
      ! ---------------

      write( file_id, iostat = io_status ) predictor_index( :, l, : )

      if ( io_status /= 0 ) then
        error_status = failure
        write( message, '( "error writing channel ", i4, 1x, &
                          &"predictor indices. iostat = ", i5 )' ) &
                        l, io_status
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        close( file_id )
        return
      end if


      ! --------------------------
      ! transmittance coefficients
      ! --------------------------

      write( file_id, iostat = io_status ) transmittance_coefficients( :, :, l, : )

      if ( io_status /= 0 ) then
        error_status = failure
        write( message, '( "error writing channel ", i4, 1x, &
                          &"transmittance coefficients. iostat = ", i5 )' ) &
                        l, io_status
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        close( file_id )
        return
      end if


    end do



    !#--------------------------------------------------------------------------#
    !#                          -- close the file --                            #
    !#--------------------------------------------------------------------------#

    close( file_id, iostat = io_status )

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error closing ", a, ". iostat = ", i5 )' ) &
                      trim( coefficient_file ), io_status
      call display_message( routine_name, &
                            trim( message ), &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#                      -- successful completion --                         #
    !#--------------------------------------------------------------------------#

    ! -- output an info message
!    write( message, '( "file version: ", i1, ".", i2.2, 2x, &
!                      &"n_predictors=",i1,2x,&
!                      &"n_absorber_layers=",i3,2x,&
!                      &"n_channels=",i4,2x,&
!                      &"n_absorbers=",i1 )' ) &
!                    file_release, file_version, &
!                    n_tc_predictors, &
!                    n_tc_layers, &
!                    n_tc_channels, &
!                    n_tc_absorbers
!    call display_message( routine_name, &
!                          trim( message ), &
!                          information, &
!                          message_log = message_log )

    error_status = success

  end function write_tau_coefficients

end module transmittance_coefficients


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
! revision 1.1  2002/11/15 15:21:33  birk
! added to cvs mainly to see how this compiles on other platforms, it currently
! seems to compile on the ibm
!
! revision 1.11  2001/08/31 21:14:36  paulv
! - added min and max release/version parameters to allow for valid use of
!   data files within a specified range.
! - added the absorber space exponential parameter alpha to the transittance
!   file data. read_ and write_ functions updated as well as the valid version
!   number.
!
! revision 1.10  2001/08/16 17:17:54  paulv
! - updated documentation
! - the comparison of n_channels and max_n_channels is now done via the
!   max_n_channels methods in the parameters module.
!
! revision 1.9  2001/08/09 20:47:25  paulv
! - added the write_transmittance_coefficients function.
! - moved all the transmittance data type and name definitions from the
!   coefficient_utility module to this one. altered use statement of the
!   coefficient_utility module to reflect this change.
! - added valid_release and valid_version parameters for data file version
!   checking.
! - added data file release and version number read/write to the requisite
!   read/write function.
!
! revision 1.8  2001/08/01 17:07:02  paulv
! - the absorber space levels are no longer calculated during model
!   initialisation, but are precalculated and stored in the transmittance
!   coefficient data file. thus a shared data array, absorber_space_levels,
!   was added to this module as were array allocation/deallocation and
!   data read statements.
!
! revision 1.7  2001/07/12 17:49:04  paulv
! - removed definitions of the number, type, and name of the transmittance items
!   and moved them into the coefficient_utility module. they are now available
!   via:
!     use coefficient_utility, only: open_coefficient_file, &
!                                    n_transmittance_items,      &
!                                    transmittance_data_type,    &
!                                    transmittance_data_name
!   this was done to allow the coefficient_utility module to be used for
!   reformatting. however, this may change - now definitions for the contents
!   of the transmittance coefficient data file are distributed in two different
!   modules. i don't like that.
!
! revision 1.6  2001/05/29 17:50:31  paulv
! - added destroy_transmittance_coefficients function. this deallocates the
!   data arrays used to store the transmittance coefficient data.
!
! revision 1.5  2000/11/09 20:42:11  paulv
! - coefficient arrays are now allocatable.
! - input file format has changed to contain data dimension and type
!   information for file data checks and array allocation.
! - coefficient data is now read directly into the correct channel position
!   straight from the file. previosuly a dummy chunk was read and then placed
!   into the shared data array in a separate loop. the new method is faster
!   and a lot easier to follow.
!
! revision 1.4  2000/08/31 19:36:34  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.3  2000/08/24 16:17:02  paulv
! - dimension order of predictor_index changed from:
!
!     0:max_n_predictors_to_use x max_n_absorbers x max_n_channels
!   to:
!     0:max_n_predictors_to_use x max_n_channels  x max_n_absorbers
!
!   and for tau_coefficients from
!
!     0:max_n_predictors_to_use x max_n_absorbers x 0:max_n_absorber_layers x max_n_channels
!   to
!     0:max_n_predictors_to_use x 0:max_n_absorber_layers x max_n_channels x max_n_absorbers
!
!   this allowed for more efficent access to the various data on an absorber
!   by absorber basis (since if there is no significant channel absorption by
!   any particular absorber go to the next absorber).
! - removed references to the record length parameter. no longer needed as
!   file access is sequential rather than direct.
! - replaced error check after open_coefficient_file call with a simple
!   test for error_status /= success. the open function no longer returns
!   any error status other than success or failure (used to return warning
!   in some circumstances.)
! - "rec =" keyword removed from file read statement.
! - tau_coefficients array filled by looping over the number of absorber layers.
!   this is only marginally faster than using array syntax but it makes the code
!   a bit easier to understand (imo).
! - channel loop construct name changed from "channel_loop" to "l_channel_loop"
!   to indicate the loop counter variable is "l". this is not a big deal for
!   this situation but has proven useful in other modules with a high degree
!   of nested loops.
! - updated module and subprogram documentation.
!
! revision 1.2  2000/08/08 17:03:38  paulv
! module modified to:
! - read the transmittance coefficients correctly! and
! - to use the parameters module rather than the constants module.
!
! revision 1.1  2000/07/12 16:08:10  paulv
! initial checked in version
!
!
!
