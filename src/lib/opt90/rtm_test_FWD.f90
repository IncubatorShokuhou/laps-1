!------------------------------------------------------------------------------
!
! name:
!       rtm_test_fwd
!
! purpose:
!       program to test the fwd rtm routine
!
! category:
!       ncep rtm
!
! language:
!       fortran-90
!
! modules:
!       type_kinds:            module to hold specification kinds for variable
!                              declaration.
!
!       file_utility:          module containing generic file utility routines
!
!       error_handler:         module to define simple error codes and handle
!                              error conditions
!
!       parameters:            module to hold rt model parameter constants
!
!       initialize:            module for rt model initialisation.
!
!       forward_model:         module containing the rt forward model function.
!
!       spectral_coefficients: module to hold the rt model spectral coefficients
!                              and their access routines.
!                              note: this routine is used only to access a single
!                                    data entity for test purposes. does not need 
!                                    to be used to call the rt routines.
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
!       input:  binary, direct access profile data file.
!       output: ascii forward model output file.
!
! side effects:
!       output file(s) are overwritten.
!
! restrictions:
!       none.
!
! creation history:
!       written by:     paul van delst, cimss/ssec 20-aug-2001
!                       paul.vandelst@ssec.wisc.edu
!
!  copyright (c) 2001 paul van delst
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
!
!------------------------------------------------------------------------------


program rtm_test_fwd


  !#----------------------------------------------------------------------------#
  !#                             -- module usage --                             #
  !#----------------------------------------------------------------------------#

  use type_kinds
  use file_utility
  use error_handler
  use parameters
  use initialize
  use forward_model, only: compute_rtm

  ! -- for access to is_microwave_channel only.
  use spectral_coefficients



  !#----------------------------------------------------------------------------#
  !#                           -- type declarations --                          #
  !#----------------------------------------------------------------------------#

  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------------
  ! program parameters
  ! ------------------

  ! -- name
  character( * ),  parameter :: program_name = 'rtm_test_fwd'

  ! -- angle definitions
  integer,         parameter :: n_angles  = 11
  real( fp_kind ), parameter :: min_angle =  0.0_fp_kind
  real( fp_kind ), parameter :: max_angle = 60.0_fp_kind

  ! -- profile definitions
  character( * ),  parameter :: profile_file        = 'profiles_layer.bin'
  integer,         parameter :: n_profiles          = 20
  integer,         parameter :: n_layers            = 100

  ! -- other dimension parameters
  integer,         parameter :: n_absorbers  = max_n_absorbers
  integer,         parameter :: n_predictors = max_n_predictors

  ! -- emissivity parameters
  real( fp_kind ), parameter :: default_mw_emissivity = 0.6_fp_kind
  real( fp_kind ), parameter :: default_ir_emissivity = 0.96_fp_kind

  ! -- default solar angle secant (> 11.47 means no solar)
  real( fp_kind ), parameter :: default_secant_solar_angle = 12.0_fp_kind


  ! -----------------
  ! program variables
  ! -----------------

  ! -- error message
  character( 128 ) :: message 

  ! -- status variables
  integer :: error_status
  integer :: allocate_status
  integer :: io_status

  ! -- loop counters, dimensions and file lun
  integer :: i, j, k, l, m, il
  integer :: n_available_channels, l1, l2, begin_channel, end_channel, n_channels
  integer :: record_length, record_number
  integer :: lun
  integer :: lun_fwd

  ! -- secant angle arrays
  integer         :: i_angle
  real( fp_kind ) :: d_angle
  real( fp_kind ), dimension( n_angles ) :: view_angle
  real( fp_kind ), dimension( n_angles ) :: secant_view_angle
  real( fp_kind ), dimension( n_angles ) :: secant_solar_angle

  ! -- profile read array
  real( fp_kind ), dimension( n_layers ) :: level_pressure
  real( fp_kind ), dimension( n_layers ) :: layer_pressure
  real( fp_kind ), dimension( n_layers ) :: layer_temperature
  real( fp_kind ), dimension( n_layers ) :: layer_water_vapor
  real( fp_kind ), dimension( n_layers ) :: layer_ozone

  ! -- profile data arrays
  real( fp_kind ), dimension( n_layers, n_angles ) :: level_p
  real( fp_kind ), dimension( n_layers, n_angles ) :: layer_p
  real( fp_kind ), dimension( n_layers, n_angles ) :: layer_t
  real( fp_kind ), dimension( n_layers, n_angles ) :: layer_w
  real( fp_kind ), dimension( n_layers, n_angles ) :: layer_o

  ! -- number of channels processed for each profile
  integer, dimension( n_angles ) :: n_channels_per_profile

  ! -- allocatable arrays
  real( fp_kind ), dimension( : ),    allocatable :: surface_emissivity
  real( fp_kind ), dimension( : ),    allocatable :: surface_reflectivity
  integer,         dimension( : ),    allocatable :: channel_index

  real( fp_kind ), dimension( :, : ), allocatable :: tau
  real( fp_kind ), dimension( :, : ), allocatable :: flux_tau
  real( fp_kind ), dimension( :, : ), allocatable :: solar_tau

  real( fp_kind ), dimension( : ),    allocatable :: upwelling_radiance
  real( fp_kind ), dimension( : ),    allocatable :: brightness_temperature

  ! -- stuff for emissivities
  real( fp_kind ) :: angle_modifier
  real( fp_kind ), dimension( n_angles ) :: mw_emissivity
  real( fp_kind ), dimension( n_angles ) :: ir_emissivity


  ! ----------
  ! intrinsics
  ! ----------

  intrinsic cos, &
            min, max, &
            real, &
            trim



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                         -- initialize the rtm --                           #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  write( *, '( /5x, "initialising rtm...." )' )

  error_status = initialize_rtm( tau_file = 'test_transmittance_coefficients', &
                                 spectral_file = 'test_spectral_coefficients'  )

  if ( error_status /= success ) then
    call display_message( program_name, &
                          'error initialzing rtm', &
                          error_status )
    stop
  end if




  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#               -- setup data and arrays for calculations --                 #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  !#----------------------------------------------------------------------------#
  !#                    -- allocate arrays for rtm model --                     #
  !#----------------------------------------------------------------------------#

  ! ----------------------------------
  ! get the number of channels read in
  ! from the coefficient data files
  ! ----------------------------------

  call get_max_n_channels( n_channels )

  begin_channel = 1
  end_channel   = n_channels

!!!replace the above code with the following to process channel subsets
!!!
!!!  call get_max_n_channels( n_available_channels )
!!!
!!!  ! -- get user channel limits
!!!  write( *, '( /5x, "number of available channels: ", i4 )' ) n_available_channels
!!!  write( *, '( /5x, "enter begin,end channels to process: " )', &
!!!            advance = 'no' )
!!!  read( *, * ) l1, l2
!!!
!!!  ! -- only keep valid values
!!!  l1 = max( 1, min( l1, n_available_channels ) )
!!!  l2 = max( 1, min( l2, n_available_channels ) )
!!!  begin_channel = min( l1, l2 )
!!!  end_channel   = max( l1, l2 )
!!!  n_channels    = end_channel - begin_channel + 1

  write( *, '( /5x, "processing channels ", i4, " to ", i4, "...", / )' ) &
            begin_channel, end_channel


  ! ---------------------------------------------------
  ! allocate the forward model channel dependent arrays
  ! ---------------------------------------------------

  allocate( surface_emissivity(   n_channels * n_angles ),   &  ! input,  l*m
            surface_reflectivity( n_channels * n_angles ),   &  ! input,  l*m
            channel_index(        n_channels * n_angles ),   &  ! input,  l*m

            tau(       n_layers, n_channels * n_angles ),    &  ! output, k x l*m
            flux_tau(  n_layers, n_channels * n_angles ),    &  ! output, k x l*m
            solar_tau( n_layers, n_channels * n_angles ),    &  ! output, k x l*m

            upwelling_radiance(     n_channels * n_angles ), &  ! output, l*m
            brightness_temperature( n_channels * n_angles ), &  ! output, l*m

            stat = allocate_status )

  if ( allocate_status /= 0 ) then
    write( message, '( "error allocating forward model channel ", &
                      &"dependent arrays. stat = ", i5 )' ) &
                    allocate_status
    call display_message( program_name,    &
                          trim( message ), &
                          failure          )
    stop
  end if



  !#----------------------------------------------------------------------------#
  !#                 -- fill profile independent input arrays --                #
  !#----------------------------------------------------------------------------#

  ! -- compute delta angle value
  d_angle = ( max_angle - min_angle ) / real( n_angles-1, fp_kind )

  ! -- initialise angle x channel counter
  il = 0


  ! ----------------
  ! loop over angles
  ! ----------------

  do i_angle = 1, n_angles

    ! -- fill angle arrays
    view_angle( i_angle )         = min_angle + ( real( i_angle-1, fp_kind ) * d_angle )
    secant_view_angle( i_angle )  = one / cos( view_angle( i_angle ) * degrees_to_radians )
    secant_solar_angle( i_angle ) = default_secant_solar_angle

    ! -- fill channel count array, i.e. do them all
    n_channels_per_profile( i_angle ) = n_channels


    ! ------------------
    ! loop over channels
    ! ------------------

    do l = 1, n_channels

      ! -- increment angle x channel counter
      ! -- and set the channel index
      il = il + 1
      channel_index( il ) = begin_channel + l - 1

      ! -- assign pretend surface emissivity and reflectivity
      ! -- for uw, r = specular; for ir, r = isotropic
      ! -- the angle modifier is just something to provide
      ! -- a little bit of angular variation in the surface
      ! -- emissivities and reflectivities.
      angle_modifier = (cos( view_angle( i_angle ) * degrees_to_radians ))**(0.1_fp_kind)
      if ( is_microwave_channel( channel_index( il ) ) == 1 ) then
        surface_emissivity( il )   = default_mw_emissivity * angle_modifier
        surface_reflectivity( il ) = one - surface_emissivity( il )
        mw_emissivity( i_angle ) = surface_emissivity( il )
      else
        surface_emissivity( il )   = default_ir_emissivity * angle_modifier
        surface_reflectivity( il ) = ( one - surface_emissivity( il ) ) / pi
        ir_emissivity( i_angle ) = surface_emissivity( il )
      end if

    end do

  end do



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                   -- open the input profile data file --                   #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  ! ------------------------------
  ! set the record length in bytes
  ! ------------------------------

  record_length = n_layers * n_bytes_for_fp_kind


  ! ---------
  ! open file
  ! ---------

  lun = get_lun()
  open( lun, file   = profile_file,  &
             status = 'old',         &
             form   = 'unformatted', &
             access = 'direct',      &
             recl   = record_length, &
             iostat = io_status      )

  if ( io_status /= 0 ) then
    error_status = failure
    call display_message( program_name, &
                          'error occurred opening file '//profile_file, &
                          error_status )
    stop
  end if


  ! ---------------------------------
  ! read the interface pressure array
  ! ---------------------------------

  record_number = 1
  read( lun, rec    = record_number, &
             iostat = io_status      ) level_pressure

  if ( io_status /= 0 ) then
    error_status = failure
    call display_message( program_name, &
                          'error occurred reading interface pressure from '//&
                            profile_file, &
                          error_status )
    stop
  end if



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                     -- open the output data files --                       #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  lun_fwd = get_lun()

  open( lun_fwd, file   = 'rtm_test.fwd.output', &
                 status = 'replace',             &
                 form   = 'formatted',           &
                 access = 'sequential',          &
                 iostat = io_status              )

  if ( io_status /= 0 ) then
    error_status = failure
    call display_message( program_name, &
                          'error occurred opening file rtm_test.fwd.output', &
                          error_status )
    stop
  end if


  ! -- write some header information
  write( lun_fwd, fmt = 100 ) 'mw', mw_emissivity
  write( lun_fwd, fmt = 100 ) 'ir', ir_emissivity
  
  ! -- output the requisite dimensions
  write( lun_fwd, fmt = 110 ) n_profiles, n_angles, n_channels

  ! -- write a title
  write( lun_fwd, fmt = 120 )

  ! -- output format statements
  100 format( '! ', a, ' surface emissivity (fn. of angle) = ', 20( 2x, f5.3, : ) )
  110 format( 3( 1x, i5 ) )
  120 format( '  prof   ang    ch  angle  radiance     temp.', /, &
             &'---------------------------------------------' )



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                       -- begin loop over profiles --                       #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################


  m_profile_loop: do m = 1, n_profiles

    write( *, '( 5x, "processing profile: ", i2 )' ) m



    !#--------------------------------------------------------------------------#
    !#                           -- read a profile --                           #
    !#--------------------------------------------------------------------------#

    ! -- layer pressure
    record_number = record_number + 1
    read( lun, rec    = record_number, &
               iostat = io_status      ) layer_pressure

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading layer pressure for profile #", i2.2, &
                        &" from ", a, ". iostat = ", i5 )' ) &
                      m, profile_file, io_status
      call display_message( program_name,    &
                            trim( message ), &
                            error_status     )
      stop
    end if


    ! -- layer temperature
    record_number = record_number + 1
    read( lun, rec    = record_number, &
               iostat = io_status      ) layer_temperature

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading layer temperature for profile #", i2.2, &
                        &" from ", a, ". iostat = ", i5 )' ) &
                      m, profile_file, io_status
      call display_message( program_name,    &
                            trim( message ), &
                            error_status     )
      stop
    end if


    ! -- layer water vapor
    record_number = record_number + 1
    read( lun, rec    = record_number, &
               iostat = io_status      ) layer_water_vapor

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading layer water vapor for profile #", i2.2, &
                        &" from ", a, ". iostat = ", i5 )' ) &
                      m, profile_file, io_status
      call display_message( program_name,    &
                            trim( message ), &
                            error_status     )
      stop
    end if


    ! -- layer ozone
    record_number = record_number + 1
    read( lun, rec    = record_number, &
               iostat = io_status      ) layer_ozone

    if ( io_status /= 0 ) then
      error_status = failure
      write( message, '( "error reading layer ozone for profile #", i2.2, &
                        &" from ", a, ". iostat = ", i5 )' ) &
                      m, profile_file, io_status
      call display_message( program_name,    &
                            trim( message ), &
                            error_status     )
      stop
    end if


    ! ---------------------------------------------------------
    ! load rtm input profile arrays
    ! (using same profile for all angles - not efficient since
    !  all quantities are recalculated for each angle, but, eh,
    !  it's a test.)
    ! ---------------------------------------------------------

    do i_angle = 1, n_angles

      level_p( :, i_angle ) = level_pressure
      layer_p( :, i_angle ) = layer_pressure
      layer_t( :, i_angle ) = layer_temperature
      layer_w( :, i_angle ) = layer_water_vapor
      layer_o( :, i_angle ) = layer_ozone

    end do


    !#--------------------------------------------------------------------------#
    !#                            -- forward model --                           #
    !#--------------------------------------------------------------------------#

    error_status = compute_rtm( &
                                ! -- forward inputs
                                level_p,                &  ! input,  k x m
                                layer_p,                &  ! input,  k x m
                                layer_t,                &  ! input,  k x m
                                layer_w,                &  ! input,  k x m
                                layer_o,                &  ! input,  k x m

                                layer_t(n_layers,:),    &  ! input, m
                                surface_emissivity,     &  ! input, l*m
                                surface_reflectivity,   &  ! input, l*m

                                secant_view_angle,      &  ! input, m
                                secant_solar_angle,     &  ! input, m
                                n_channels_per_profile, &  ! input, m
                                channel_index,          &  ! input, l*m

                                ! -- forward outputs
                                tau,                    &  ! output, k x l*m
                                flux_tau,               &  ! output, k x l*m
                                solar_tau,              &  ! output, k x l*m

                                upwelling_radiance,     &  ! output, l*m
                                brightness_temperature  )  ! output, l*m

    if ( error_status /= success ) then
      call display_message( program_name, &
                            'error occured in compute_rtm', &
                            error_status )
      stop
    end if



    ! --------------------------------
    ! output the forward model results
    ! --------------------------------

    ! -- initialise the channel/profile counter
    il = 0

    ! -- loop over the number of angles
    do i_angle = 1, n_angles

      ! -- loop over the channels
      do l = 1, n_channels_per_profile( i_angle )

        ! -- increment the channel counter and output the
        ! -- channel radiance, and brightness temperature
        il = il + 1
        write( lun_fwd, '( 3( 1x, i5 ), 1x, f6.3, 2( 1x, f9.5 ) )' ) &
                        m, i_angle, l, &
                        view_angle( i_angle ), &
                        upwelling_radiance( il ), &
                        brightness_temperature( il )

      end do

    end do


  end do m_profile_loop

  close( lun_fwd )
  close( lun )



  !#----------------------------------------------------------------------------#
  !#                           -- destroy the rtm --                            #
  !#----------------------------------------------------------------------------#

  ! ---------------------------------------
  ! deallocate the channel dependent arrays
  ! ---------------------------------------

  ! -- forward model
  deallocate( surface_emissivity,     &  ! input,  l*m
              surface_reflectivity,   &  ! input,  l*m
              channel_index,          &  ! input,  l*m

              tau,                    &  ! output, k x l*m
              flux_tau,               &  ! output, k x l*m
              solar_tau,              &  ! output, k x l*m

              upwelling_radiance,     &  ! output, l*m
              brightness_temperature, &  ! output, l*m

              stat = allocate_status  )

  if ( allocate_status /= 0 ) then
    write( message, '( "error deallocating forward model channel ", &
                      &"dependent arrays. stat = ", i5 )' ) &
                    allocate_status
    call display_message( program_name,    &
                          trim( message ), &
                          warning          )
  end if


  ! ---------------------------------
  ! deallocate the coefficient arrays
  ! ---------------------------------

  write( *, '( /5x, "destroying rtm...." )' )

  error_status = destroy_rtm()

  if ( error_status /= success ) then
    call display_message( program_name, &
                          'error destroying rtm', &
                          error_status )
    stop
  end if

end program rtm_test_fwd


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
! revision 1.1  2001/09/13 22:08:13  paulv
! initial checkin.
!
!
!
!
