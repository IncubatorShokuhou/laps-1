!------------------------------------------------------------------------------
!m+
! name:
!       spectral_coefficients
!
! purpose:
!       module to hold the rt model spectral coefficients and their access
!       routines. 
!
! category:
!       ncep rtm
!
! calling sequence:
!       use spectral_coefficients
!
! outputs:
!       frequency:                     channel frequency
!                                      units:      ghz
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       wavenumber:                    channel frequency
!                                      units:      cm^-1
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       planck_c1:                     first planck function coefficient
!                                      units:      mw/(m^2.sr.cm^-2)
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       planck_c2:                     second planck function coefficient
!                                      units:      k/cm^1
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       band_c1:                       polychromatic band correction offset.
!                                      units:      k
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       band_c2:                       polychromatic band correction slope.
!                                      units:      k/k
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       is_microwave_channel:          flag indicating if the particular satellite
!                                      channel is a microwave channel.
!                                      if = 0, ir channel
!                                         = 1, uw channel
!                                      units:      none
!                                      type:       long
!                                      dimension:  l
!                                      attributes: public, save
!
!       cosmic_background_temperature: cosmic background temperature in kelvin to
!                                      use in the radiative transfer. these values
!                                      are frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      rayleigh-jeans approximation, but the radiative
!                                      transfer using the full planck function
!                                      expression. ir channels use a fixed value of
!                                      2.736k.
!                                      units:      k
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       cosmic_background_radiance:    cosmic background radiance used to initialise
!                                      the radiative transfer. these values
!                                      are frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      rayleigh-jeans approximation, but the radiative
!                                      transfer using the full planck function
!                                      expression. ir channel values are all set to 0.0.
!                                      units:      mw/(m^2.sr.cm^-1)
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       solar_irradiance:              exterrestrial solar irradiance source values
!                                      (based on kurucz).
!                                      units:      mw/(m^2.cm^-1)
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       blackbody_irradiance:          exterrestrial solar irradiance source values
!                                      for a blackbody source at the equivalent solar
!                                      temperature of 5783k.
!                                      units:      mw/(m^2.cm^-1)
!                                      type:       double
!                                      dimension:  l
!                                      attributes: public, save
!
!       is_solar_channel:              array of integer flag values indicating if the
!                                      particular satellite channel possesses valid
!                                      solar irradiance values.
!                                      if = 0, no
!                                         = 1, yes
!                                      units:      none
!                                      type:       long
!                                      dimension:  l
!                                      attributes: public, save
!
! modules:
!       type_kinds:           module containing data type kind definitions.
!
!       error_handler:        module to define error codes and handle error
!                             conditions
!
!       parameters:           module containing parameter definitions for the
!                             rt model.
!
!       coefficient_utility:  module containing coefficient file utility subprograms.
!
! contains:
!       read_spectral_coefficients:     public function to read the spectral coefficients
!                                       for the satellite/channels used by the rt model
!                                       and fill the public data arrays.
!
!       destroy_spectral_coefficients:  public function to release the memory that
!                                       was allocated for the spectral coefficient data.
!
!       write_spectral_coefficients:    public function to write the new format spectral
!                                       coefficient data file.
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

module spectral_coefficients


  ! ---------------------
  ! module use statements
  ! ---------------------

  use type_kinds
  use file_utility
  use error_handler
  use parameters
  use coefficient_utility


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------------
  ! default visibilities
  ! ------------------

  private
  public :: read_spectral_coefficients
  public :: destroy_spectral_coefficients
  public :: write_spectral_coefficients



  ! -----------------
  ! module parameters
  ! -----------------

  ! -- spectral coefficient file version numbers
  integer( long ), parameter, private :: min_valid_release = 1_long
  integer( long ), parameter, private :: max_valid_release = 1_long

  integer( long ), parameter, private :: min_valid_version = 1_long
  integer( long ), parameter, private :: max_valid_version = 1_long

  ! -- spectral coefficient file descriptors:
  ! -- number of items in the spectral coefficient file
  integer( long ), parameter, private :: n_spectral_items = 12_long

  ! -- data types of the spectral coefficient data
  !    5 = double (i.e. 8-byte float)
  !    4 = single (i.e. 4-byte float)
  !    3 = long   (i.e. 4-byte integer)
  integer( long ), parameter, private, &
                   dimension( n_spectral_items ) :: spectral_data_type = &
                                                    (/ 5_long, 5_long, 5_long, 5_long, &
                                                       5_long, 5_long, 3_long, 5_long, &
                                                       5_long, 5_long, 5_long, 3_long /)

  ! -- names of the data items (for error processing)
  character( * ), parameter, private, &
                  dimension( n_spectral_items ) :: spectral_data_name = &
                                                     (/ 'frequency                    ', &
                                                        'wavenumber                   ', &
                                                        'planck_c1                    ', &
                                                        'planck_c2                    ', &
                                                        'band_c1                      ', &
                                                        'band_c2                      ', &
                                                        'is_microwave_channel         ', &
                                                        'cosmic_background_temperature', &
                                                        'cosmic_background_radiance   ', &
                                                        'solar_irradiance             ', &
                                                        'blackbody_irradiance         ', &
                                                        'is_solar_channel             ' /)

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

  ! -- channel frequency in ghz (mw only, 0.0 for ir)
  real( double ),  save, public, allocatable, dimension( : ) :: frequency

  ! -- channel frequency in cm^-1
  real( double ),  save, public, allocatable, dimension( : ) :: wavenumber

  ! -- first channel planck function coefficient in mw/(m^2.sr.cm^-2)
  real( double ),  save, public, allocatable, dimension( : ) :: planck_c1

  ! -- second channel planck function coefficient in k/cm^1
  real( double ),  save, public, allocatable, dimension( : ) :: planck_c2

  ! -- band correction (for polychromaticity) offset in k
  real( double ),  save, public, allocatable, dimension( : ) :: band_c1

  ! -- band correction (for polychromaticity) slope in k/k
  real( double ),  save, public, allocatable, dimension( : ) :: band_c2

  ! -- flag indicating whether channel is microwave or infrared.
  integer( long ), save, public, allocatable, dimension( : ) :: is_microwave_channel

  ! -- effective cosmic background temperature in kelvin.
  real( double ),  save, public, allocatable, dimension( : ) :: cosmic_background_temperature

  ! -- cosmic background radiance in mw/(m^2.sr.cm^-1).
  real( double ),  save, public, allocatable, dimension( : ) :: cosmic_background_radiance

  ! -- kurucz solar irradiance source function in mw/(m^2.cm^-1)
  real( double ),  save, public, allocatable, dimension( : ) :: solar_irradiance

  ! -- equivalent blackbody solar irradiance source function in mw/(m^2.cm^-1)
  real( double ),  save, public, allocatable, dimension( : ) :: blackbody_irradiance

  ! -- flag indicating whether channel is sensitive to solar contribution.
  integer( long ), save, public, allocatable, dimension( : ) :: is_solar_channel



contains


!------------------------------------------------------------------------------
!s+
! name:
!       read_spectral_coefficients
!
! purpose:
!       public function to read the spectral coefficients for the satellite/
!       channels used by the rt model and fill the public data arrays..
!
! calling sequence:
!       result = read_spectral_coefficients( coefficient_file, &
!                                            message_log  = message_log )
!
! input arguments:
!       coefficient_file: name of the file containing the spectral coefficient data.
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
! procedure:
!       all spectral coefficients are read in and stored by channel.
!s-
!------------------------------------------------------------------------------

  function read_spectral_coefficients( coefficient_file, &
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

    character( * ), parameter :: routine_name = 'read_spectral_coefficients'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- these are dimensioned ( long ) as they are read in
    integer( long ) :: file_release
    integer( long ) :: file_version

    integer( long ) :: n_channels
    integer( long ) :: n_items
    integer( long ), dimension( n_spectral_items ) :: data_type

    integer :: l
    integer :: io_status
    integer :: file_id
    integer :: allocate_status

    ! -- maximum channels pseudo parameter
    integer :: max_n_channels
    logical :: is_set


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic trim                          


    !#--------------------------------------------------------------------------#
    !#              -- open the spectral coefficient data file --               #
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


    ! ---------------------------
    ! read the number of channels
    ! ---------------------------

    read( file_id, iostat = io_status ) n_channels

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading total number of channels. iostat = ", i5 )' ) &
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

    if ( n_items /= n_spectral_items ) then
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

      if ( data_type( l ) /= spectral_data_type( l ) ) then
        error_status = failure
        write( message, '( "invalid type for ", a, " data item." )' ) &
                        trim( spectral_data_name( l ) )
        call display_message( routine_name, &
                              trim( message ), &
                              error_status, &
                              message_log = message_log )
        close( file_id )
        return
      end if

    end do


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
    !#          -- allocate arrays for spectral coefficient data --             #
    !#--------------------------------------------------------------------------#

    ! -- check if arrays are already allocated
    if ( allocated( frequency                     ) .or. &
         allocated( wavenumber                    ) .or. &
         allocated( planck_c1                     ) .or. &
         allocated( planck_c2                     ) .or. &
         allocated( band_c1                       ) .or. &
         allocated( band_c2                       ) .or. &
         allocated( is_microwave_channel          ) .or. &
         allocated( cosmic_background_temperature ) .or. &
         allocated( cosmic_background_radiance    ) .or. &
         allocated( solar_irradiance              ) .or. &
         allocated( blackbody_irradiance          ) .or. &
         allocated( is_solar_channel              )      ) then
      error_status = failure
      call display_message( routine_name, &
                            'spectral coefficient data arrays already allocated.', &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    endif

    ! -- if not, allocate them
    allocate( frequency( n_channels ),                     &
              wavenumber( n_channels ),                    &
              planck_c1( n_channels ),                     &
              planck_c2( n_channels ),                     &
              band_c1( n_channels ),                       &
              band_c2( n_channels ),                       &
              is_microwave_channel( n_channels ),          &
              cosmic_background_temperature( n_channels ), &
              cosmic_background_radiance( n_channels ),    &
              solar_irradiance( n_channels ),              &
              blackbody_irradiance( n_channels ),          &
              is_solar_channel( n_channels ),              &
              stat = allocate_status                       )

    if ( allocate_status /= 0 ) then
      error_status = failure
      call display_message( routine_name, &
                            'unable to allocate arrays for spectral coefficient data.', &
                            error_status, &
                            message_log = message_log )
      close( file_id )
      return
    end if
   


    !#--------------------------------------------------------------------------#
    !#                        -- loop over channels --                          #
    !#--------------------------------------------------------------------------#

    l_channel_loop: do l = 1, n_channels


      ! -------------------
      ! read channel record
      ! -------------------

      read( file_id, iostat = io_status )   &
        frequency( l ),                     &  ! frequency
        wavenumber( l ),                    &  ! wavenumber
        planck_c1( l ),                     &  ! first planck constant
        planck_c2( l ),                     &  ! second planck constant
        band_c1( l ),                       &  ! band correction offset
        band_c2( l ),                       &  ! band correction slope
        is_microwave_channel( l ),          &  ! microwave channel flag
        cosmic_background_temperature( l ), &  ! effective cosmic background temperature
        cosmic_background_radiance( l ),    &  ! cosmic background radiance
        solar_irradiance( l ),              &  ! solar irradiance source
        blackbody_irradiance( l ),          &  ! blackbody solar irradiance source
        is_solar_channel( l )                  ! solar use flag


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
!                      &"n_channels=",i4 )' ) &
!                    file_release, file_version, &
!                    n_channels
!    call display_message( routine_name, &
!                          trim( message ), &
!                          information, &
!                          message_log = message_log )
 
    error_status = success

  end function read_spectral_coefficients 




!------------------------------------------------------------------------------
!s+
! name:
!       destroy_spectral_coefficients
!
! purpose:
!       public function to deallocate the spectral coefficients public
!       data arrays..
!
! calling sequence:
!       result = destroy_spectral_coefficients( message_log = message_log )
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
!              = failure => error occurred deallocating coefficient data arrays
!
! calls:
!      display_message:         subroutine to output messages
!                               source: error_handler module
!
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

  function destroy_spectral_coefficients( message_log ) &
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

    character( * ), parameter :: routine_name = 'destroy_spectral_coefficients'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: allocate_status



    !#--------------------------------------------------------------------------#
    !#           -- deallocate arrays for spectral coefficient data --          #
    !#--------------------------------------------------------------------------#

    deallocate( frequency,                     &
                wavenumber,                    &
                planck_c1,                     &
                planck_c2,                     &
                band_c1,                       &
                band_c2,                       &
                is_microwave_channel,          &
                cosmic_background_temperature, &
                cosmic_background_radiance,    &
                solar_irradiance,              &
                blackbody_irradiance,          &
                is_solar_channel,              &
                stat = allocate_status         )

    if ( allocate_status /= 0 ) then
      error_status = failure
      call display_message( routine_name, &
                            'error occurred deallocating spectral coefficient data arrays.', &
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

  end function destroy_spectral_coefficients 





!------------------------------------------------------------------------------
!
! name:
!       write_spectral_coefficients
!
! purpose:
!       function to write the spectral coefficient data to file.
!
! category:
!       ncep rtm
!
! language:
!       fortran-90
!
! calling sequence:
!       result = write_spectral_coefficients( coefficient_file,              &  ! input
!                                             frequency,                     &  ! input, l
!                                             wavenumber,                    &  ! input, l
!                                             planck_c1,                     &  ! input, l
!                                             planck_c2,                     &  ! input, l
!                                             band_c1,                       &  ! input, l
!                                             band_c2,                       &  ! input, l
!                                             is_microwave_channel,          &  ! input, l
!                                             cosmic_background_temperature, &  ! input, l
!                                             cosmic_background_radiance     &  ! input, l
!                                             solar_irradiance,              &  ! input, l
!                                             blackbody_irradiance,          &  ! input, l
!                                             is_solar_channel,              &  ! input, l
!                                             message_log                    )  ! optional input
!
! input arguments:
!       coefficient_file:              name of the file to which the spectral
!                                      data is to be written.
!                                      units:      none
!                                      type:       character
!                                      dimension:  scalar
!                                      attributes: intent( in )
!
!       frequency:                     channel frequency
!                                      units:      ghz
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       wavenumber:                    channel frequency
!                                      units:      cm^-1
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       planck_c1:                     first planck function coefficient
!                                      units:      mw/(m^2.sr.cm^-2)
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       planck_c2:                     second planck function coefficient
!                                      units:      k/cm^1
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       band_c1:                       polychromatic band correction offset.
!                                      units:      k
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       band_c2:                       polychromatic band correction slope.
!                                      units:      k/k
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       is_microwave_channel:          flag indicating if the particular satellite
!                                      channel is a microwave channel.
!                                      if = 0, ir channel
!                                         = 1, uw channel
!                                      units:      none
!                                      type:       long integer
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       cosmic_background_temperature: cosmic background temperature in kelvin to
!                                      use in the radiative transfer. these values
!                                      may be frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      rayleigh-jeans approximation, but the radiative
!                                      transfer using the full planck function
!                                      expression. ir channels use a fixed value of
!                                      2.736k.
!                                      units:      k
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       cosmic_background_radiance:    cosmic background radiance used to initialise
!                                      the radiative transfer. these values
!                                      may be frequency dependent for microwave channels
!                                      to account for instrument calibration using the
!                                      rayleigh-jeans approximation, but the radiative
!                                      transfer using the full planck function
!                                      expression. ir channel values are all set to 0.0.
!                                      units:      mw/(m^2.sr.cm^-1)
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       solar_irradiance:              exterrestrial solar irradiance source values
!                                      (based on kurucz).
!                                      units:      mw/(m^2.cm^-1)
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       blackbody_irradiance:          exterrestrial solar irradiance source values
!                                      for a blackbody source at the equivalent solar
!                                      temperature of 5783k.
!                                      units:      mw/(m^2.cm^-1)
!                                      type:       double
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
!       is_solar_channel:              array of integer flag values indicating if the
!                                      particular satellite channel possesses valid
!                                      solar irradiance values.
!                                      if = 0, no
!                                         = 1, yes
!                                      units:      none
!                                      type:       long
!                                      dimension:  l, n_channels
!                                      attributes: intent( in )
!
! optional input arguments:
!       message_log:                   character string specifying a filename in which
!                                      any messages will be logged. if not specified,
!                                      or if an error occurs opening the log file,
!                                      the default action is to output messages to
!                                      the screen.
!                                      units:      none
!                                      type:       character
!                                      dimension:  scalar
!                                      attributes: intent( in ), optional
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

  function write_spectral_coefficients( coefficient_file,              &  ! input
                                        frequency,                     &  ! input
                                        wavenumber,                    &  ! input
                                        planck_c1,                     &  ! input
                                        planck_c2,                     &  ! input
                                        band_c1,                       &  ! input
                                        band_c2,                       &  ! input
                                        is_microwave_channel,          &  ! input
                                        cosmic_background_temperature, &  ! input
                                        cosmic_background_radiance,    &  ! input
                                        solar_irradiance,              &  ! input
                                        blackbody_irradiance,          &  ! input
                                        is_solar_channel,              &  ! input

                                        release,                       &  ! optional input
                                        version,                       &  ! optional input
                                        message_log )                  &  ! optional input
                                      result ( error_status )


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- input
    character( * ),                  intent( in ) :: coefficient_file              ! input
    real( double ),  dimension( : ), intent( in ) :: frequency                     ! input, l
    real( double ),  dimension( : ), intent( in ) :: wavenumber                    ! input, l
    real( double ),  dimension( : ), intent( in ) :: planck_c1                     ! input, l
    real( double ),  dimension( : ), intent( in ) :: planck_c2                     ! input, l
    real( double ),  dimension( : ), intent( in ) :: band_c1                       ! input, l
    real( double ),  dimension( : ), intent( in ) :: band_c2                       ! input, l
    integer( long ), dimension( : ), intent( in ) :: is_microwave_channel          ! input, l
    real( double ),  dimension( : ), intent( in ) :: cosmic_background_temperature ! input, l
    real( double ),  dimension( : ), intent( in ) :: cosmic_background_radiance    ! input, l
    real( double ),  dimension( : ), intent( in ) :: solar_irradiance              ! input, l
    real( double ),  dimension( : ), intent( in ) :: blackbody_irradiance          ! input, l
    integer( long ), dimension( : ), intent( in ) :: is_solar_channel              ! input, l


    integer,        optional,        intent( in ) :: release
    integer,        optional,        intent( in ) :: version

    character( * ), optional,        intent( in ) :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ),  parameter :: routine_name = 'write_spectral_coefficients'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- typed ( long ) for output
    integer( long ) :: n_channels

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

    n_channels = size( frequency )

    if ( size( wavenumber                    ) /= n_channels .or. &
         size( planck_c1                     ) /= n_channels .or. &
         size( planck_c2                     ) /= n_channels .or. &
         size( band_c1                       ) /= n_channels .or. &
         size( band_c1                       ) /= n_channels .or. &
         size( is_microwave_channel          ) /= n_channels .or. &
         size( cosmic_background_temperature ) /= n_channels .or. &
         size( cosmic_background_radiance    ) /= n_channels .or. &
         size( solar_irradiance              ) /= n_channels .or. &
         size( blackbody_irradiance          ) /= n_channels .or. &
         size( is_solar_channel              ) /= n_channels      ) then
      error_status = failure
      call display_message( routine_name, &
                            'inconsistent vector channel dimensions.', &
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
    !#        -- open the spectral coefficient data file for output --          #
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

    write( file_id, iostat = io_status ) n_channels

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

    write( file_id, iostat = io_status ) n_spectral_items

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


    write( file_id, iostat = io_status ) spectral_data_type

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
    !#                     -- write the data by channel --                      #
    !#--------------------------------------------------------------------------#

    do l = 1, n_channels

      write( file_id, iostat = io_status ) frequency( l ),                     &                    
                                           wavenumber( l ),                    &
                                           planck_c1( l ),                     &
                                           planck_c2( l ),                     &
                                           band_c1( l ),                       &
                                           band_c2( l ),                       &
                                           is_microwave_channel( l ),          &
                                           cosmic_background_temperature( l ), &
                                           cosmic_background_radiance( l ),    &
                                           solar_irradiance( l ),              &
                                           blackbody_irradiance( l ),          &
                                           is_solar_channel( l )

      if ( io_status /= 0 ) then
        error_status = failure
        write( message, '( "error writing channel ", i4, 1x, &
                          &"spectral data. iostat = ", i5 )' ) &
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
    ! -- output an info message
    write( message, '( "file version: ", i1, ".", i2.2, 2x, &
                      &"n_channels=",i4 )' ) &
                    file_release, file_version, &
                    n_channels
    call display_message( routine_name, &
                          trim( message ), &
                          information, &
                          message_log = message_log )

    error_status = success

  end function write_spectral_coefficients

end module spectral_coefficients


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
! revision 1.11  2001/08/31 21:11:41  paulv
! - added min and max release/version parameters to allow for valid use of
!   data files within a specified range.
!
! revision 1.10  2001/08/16 17:16:30  paulv
! - updated documentation
! - the comparison of n_channels and max_n_channels is now done via the
!   max_n_channels methods in the parameters module.
!
! revision 1.9  2001/08/09 20:45:33  paulv
! - added the write_spectral_coefficients function.
! - moved all the spectral data type and name definitions from the
!   coefficient_utility module to this one. altered use statement of the
!   coefficient_utility module to reflect this change.
! - added valid_release and valid_version parameters for data file version
!   checking.
! - added data file release and version number read/write to the requisite
!   read/write function.
!
! revision 1.8  2001/07/12 17:46:31  paulv
! - removed definitions of the number, type, and name of the spectral items
!   and moved them into the coefficient_utility module. they are now available
!   via:
!     use coefficient_utility, only: open_coefficient_file, &
!                                    n_spectral_items,      &
!                                    spectral_data_type,    &
!                                    spectral_data_name
!   this was done to allow the coefficient_utility module to be used for
!   reformatting. however, this may change - now definitions for the contents
!   of the spectral coefficient data file are distributed in two different
!   modules. i don't like that.
! - compressed the coefficient read statement.
!
! revision 1.7  2001/05/29 17:49:32  paulv
! - made all real valued spectral coefficients double precision. previously
!   the cosmic background and solar terms were single precision.
! - added precalculated cosmic background radiance term.
! - some cosmetic changes.
! - changed use_solar array name to is_solar_channel.
! - added destroy_spectral_coefficients function. this deallocates the
!   data arrays used to store the spectral coefficient data.
!
! revision 1.6  2001/01/09 21:23:27  paulv
! - added is_microwave_channel and cosmic_background_temperature arrays to
!   read.
! - updated module documentation.
!
! revision 1.5  2000/11/09 20:39:49  paulv
! - coefficient arrays are now allocatable.
! - input file format has changed to contain data dimension and type
!   information for file data checks and array allocation.
!
! revision 1.4  2000/08/31 19:36:33  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.3  2000/08/24 16:06:28  paulv
! - removed references to the record length parameter. no longer needed as
!   file access is sequential rather than direct.
! - replaced error check after open_coefficient_file call with a simple
!   test for error_status /= success. the open function no longer returns
!   any error status other than success or failure (used to return warning
!   in some circumstances.)
! - "rec =" keyword removed from file read statement.
! - channel loop construct name changed from "channel_loop" to "l_channel_loop"
!   to indicate the loop counter variable is "l". this is not a big deal for
!   this situation but has proven useful in other modules with a high degree
!   of nested loops.
! - updated module and subprogram documentation.
!
! revision 1.2  2000/08/08 17:04:02  paulv
! module modified to:
! - read the spectral coefficients correctly! and
! - to use the parameters module rather than the constants module.
!
! revision 1.1  2000/07/12 16:08:10  paulv
! initial checked in version
!
!
!
