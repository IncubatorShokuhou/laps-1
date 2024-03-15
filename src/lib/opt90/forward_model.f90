!------------------------------------------------------------------------------
!m+
! name:
!       forward_model
!
! purpose:
!       module containing the ncep rt forward model function.
!
! category:
!       ncep rtm
!
! calling sequence:
!       use forward_model
!
! outputs:
!       none.
!
! modules:
!       type_kinds:            module to define kind types for variable declaration.
!
!       error_handler:         module to define error codes and handle error conditions
!
!       parameters:            module containing parameter definitions for the
!                              rt model.
!
!       spectral_coefficients: module containing the rt model spectral coefficients.
!
!       absorber_profile:      module containing routines for generating the absorber
!                              profiles.
!
!       predictors:            module containing routines for generating the predictor
!                              profiles.
!
!       transmittance:         module containing transmittance calculation routines.
!
!       radiance:              module containing radiance calculation routines.
!
! contains:
!       compute_rtm:           public function that calculates the forward model 
!                              top-of-atmosphere (toa) radiances and brightness 
!                              temperatures for an input atmospheric profile set and
!                              user specified satellites/channels.
!
!                              this function is simply a wrapper around the forward
!                              model so that the user doesn't have to declare the
!                              absorber/predictor/etc arrays in the calling routine.
!
!       forward_rtm:           public function that calculates top-of-atmosphere (toa)
!                              radiances and brightness temperatures for user specified
!                              profiles and satellites/channels.
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
!       written by:     paul van delst, cimss@noaa/ncep 15-july-2000
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

module forward_model


  ! ------------
  ! module usage
  ! ------------

  use type_kinds,            only : fp_kind
  use error_handler
  use parameters
  use spectral_coefficients, only : is_solar_channel, &
                                    is_microwave_channel
  use absorber_profile,      only : compute_absorber_amount, &
                                    find_absorber_layer_index
  use predictors,            only : compute_predictors
  use transmittance,         only : compute_transmittance
  use radiance,              only : compute_radiance


  ! -----------------------
  ! disable implicit typing
  ! -----------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private
  public :: compute_rtm
  public :: forward_rtm


  ! --------------------
  ! function overloading
  ! --------------------

  interface compute_rtm
    module procedure compute_rtm_rank1
    module procedure compute_rtm_rank2
  end interface ! compute_rtm

  interface forward_rtm
    module procedure forward_rtm_rank1
    module procedure forward_rtm_rank2
  end interface ! forward_rtm

contains


!--------------------------------------------------------------------------------
!s+
! name:
!       compute_rtm
!
! purpose:
!       public function that calculates the forward model top-of-atmosphere (toa)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
!       this function is simply a wrapper around the forward model so that the
!       user doesn't have to declare the absorber/predictor/etc arrays in the
!       calling routine.
!
! category:
!       ncep rtm
!
! calling sequence:
!       result = compute_rtm( &
!                             ! -- forward inputs
!                             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k x m
!
!                             surface_temperature,                         &  ! input, m
!                             surface_emissivity,                          &  ! input, l*m
!                             surface_reflectivity,                        &  ! input, l*m
!
!                             ! -- other inputs
!                             secant_view_angle,                           &  ! input, m
!                             secant_solar_angle,                          &  ! input, m
!                             n_channels_per_profile,                      &  ! input, m
!                             channel_index,                               &  ! input, l*m
!
!                             ! -- forward output
!                             tau,                                         &  ! input, k x l*m
!                             flux_tau,                                    &  ! input, k x l*m
!                             solar_tau,                                   &  ! input, k x l*m
!
!                             upwelling_radiance,                          &  ! input, l*m
!                             brightness_temperature,                      &  ! input, l*m
!
!                             ! optional inputs
!                             n_input_profiles = n_input_profiles,         &
!                             message_log      = message_log               )
!
! input arguments:
!
!       level_p:                   profile set layer interface pressure array. the toa
!                                  pressure is not included. toa pressure is parameterised
!                                  in the parameters module.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_p:                   profile set layer average pressure array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_t:                   profile set layer average temperature array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_w:      .            profile set layer average water vapor mixing ratio array
!                                  units:      g/kg
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_o:                   profile set layer average ozone mixing ratio array.
!                                  units:      ppmv
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       surface_temperature:       profile set surface temperature array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  m; m > or = 1 (i.e. scalar)
!                                  attributes: intent( in )
!
!       surface_emissivity:        profile set surface emissivity array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       surface_reflectivity:      profile set surface reflectivity array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       secant_view_angle:         secant of the satellite view angle measured
!                                  from nadir for each profile in the set.
!                                  units:      none
!                                  type:       real
!                                  dimension:  m; m > or = 1 (i.e. scalar)
!                                  attributes: intent( in )
!
!       secant_solar_angle:        secant of the solar zenith angle for each
!                                  profile in the set.
!                                  units:      none
!                                  type:       real
!                                  dimension:  m; m > or = 1 (i.e. scalar)
!                                  attributes: intent( in )
!
!       n_channels_per_profile:    the number of channels for each profile in the
!                                  set for which radiances are required.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  m; m > or = 1 (i.e. scalar)
!                                  attributes: intent( in )
!
!       channel_index:             channel index id array. each element is a unique
!                                  index to a (supported) sensor channel.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!
! optional input arguments:
!
!       n_input_profiles:          the number of profiles in the passed arrays to process.
!                                  if not specified, the default value is the second dimension
!                                  of the pressure array determined using the size intrinsic,
!                                  n_profiles. if n_input_profiles is specified and is < 1 or
!                                  greater than n_profiles, the default value is set to n_profiles.
!                                  this argument is ignored if the input profile arrays are
!                                  vectors, i.e. a single profile.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  scalar
!                                  attributes: intent( in ), optional
!
!       message_log:               character string specifying a filename in which any
!                                  messages will be logged. if not specified, or if an
!                                  error occurs opening the log file, the default action
!                                  is to output messages to the screen.
!                                  units:      none
!                                  type:       character
!                                  dimension:  scalar
!                                  attributes: intent( in ), optional
!
! output arguments:
!
!       tau:                       layer->toa transmittance for the satellite
!                                  view angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( out )
!
!       flux_tau:                  layer->sfc transmittance for the default
!                                  diffusivity angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( out )
!
!       solar_tau:                 layer->sfc transmittance for the solar
!                                  zenith angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( out )
!
!       upwelling_radiance:        toa radiances for each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
!       brightness_temperature:    toa brightness temperatures corresponding
!                                  to the toa radiances.
!                                  n.b.: set to zero upon output.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
! optional output arguments:
!       none.
!
! function result:
!       result = success => calculation was successful
!              = failure => error occurred
!
! calls:
!      display_message:            subroutine to output messages
!                                  source: error_handler module
!
!      get_max_n_channels:         routine to retrieve the value of the
!                                  max_n_channels "pseudo-parameter".
!                                  source: parameters module
!
!      forward_rtm:                function to construct the forward model and calculate
!                                  the transmittance profiles and toa radiance/temperatures.
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       all input adjoint quantities are set to zero on output.
!
! restrictions:
!       note the restrictions on the input array dimensions:
!         k == n_layers > 1
!         l == n_channels > 1
!         m == n_profiles > or = 1
!
! procedure:
!       see individual module function documentation.
!s-
!--------------------------------------------------------------------------------

  function compute_rtm_rank2( &
             ! -- forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k x m

             surface_temperature,                         &  ! input, m
             surface_emissivity,                          &  ! input, l*m
             surface_reflectivity,                        &  ! input, l*m

             ! -- other inputs
             secant_view_angle,                           &  ! input, m
             secant_solar_angle,                          &  ! input, m
             n_channels_per_profile,                      &  ! input, m
             channel_index,                               &  ! input, l*m

             ! -- forward output
             tau,                                         &  ! output, k x l*m
             flux_tau,                                    &  ! output, k x l*m
             solar_tau,                                   &  ! output, k x l*m

             upwelling_radiance,                          &  ! output, l*m
             brightness_temperature,                      &  ! output, l*m

             ! -- optional inputs
             n_input_profiles,                            &
             message_log )                                &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward inputs
    real( fp_kind ), dimension( :, : ), intent( in )  :: level_p                 ! k x m
    real( fp_kind ), dimension( :, : ), intent( in )  :: layer_p                 ! k x m
    real( fp_kind ), dimension( :, : ), intent( in )  :: layer_t                 ! k x m
    real( fp_kind ), dimension( :, : ), intent( in )  :: layer_w                 ! k x m
    real( fp_kind ), dimension( :, : ), intent( in )  :: layer_o                 ! k x m

    real( fp_kind ), dimension( : ),    intent( in )  :: surface_temperature     ! m
    real( fp_kind ), dimension( : ),    intent( in )  :: surface_emissivity      ! l*m
    real( fp_kind ), dimension( : ),    intent( in )  :: surface_reflectivity    ! l*m

    ! -- other inputs
    real( fp_kind ), dimension( : ),    intent( in )  :: secant_view_angle       ! m
    real( fp_kind ), dimension( : ),    intent( in )  :: secant_solar_angle      ! m
    integer,         dimension( : ),    intent( in )  :: n_channels_per_profile  ! m
    integer,         dimension( : ),    intent( in )  :: channel_index           ! l*m

    ! -- forward outputs
    real( fp_kind ), dimension( :, : ), intent( out ) :: tau                     ! k x l*m
    real( fp_kind ), dimension( :, : ), intent( out ) :: flux_tau                ! k x l*m
    real( fp_kind ), dimension( :, : ), intent( out ) :: solar_tau               ! k x l*m

    real( fp_kind ), dimension( : ),    intent( out ) :: upwelling_radiance      ! l*m
    real( fp_kind ), dimension( : ),    intent( out ) :: brightness_temperature  ! l*m

    ! -- optional input
    integer,        optional,           intent( in )  :: n_input_profiles        ! scalar
    character( * ), optional,           intent( in )  :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_rtm_rank2'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- array for integrated absorber amounts, 0:k x j x m
    real( fp_kind ), dimension( 0:size( layer_p, dim = 1 ), &
                                  max_n_absorbers,          &
                                  size( layer_p, dim = 2 )  ) :: absorber

    ! -- arrays for absorber space indexing, k x j x m
    integer,         dimension( size( layer_p, dim = 1 ), &
                                max_n_absorbers,          &
                                size( layer_p, dim = 2 )  ) :: tau_layer_index,      &
                                                               flux_tau_layer_index, &
                                                               solar_tau_layer_index

    ! -- arrays for predictors, imax x k x m
    real( fp_kind ), dimension( max_n_predictors,         &
                                size( layer_p, dim = 1 ), &
                                size( layer_p, dim = 2 )  ) :: tau_predictor,      &
                                                               flux_tau_predictor, &
                                                               solar_tau_predictor

    ! -- array for layer planck radiance term, k x l*m
    real( fp_kind ), dimension( size( layer_p, dim = 1 ),  &
                                size( upwelling_radiance ) ) :: layer_radiance
      

    ! -- array for downwelling radiance (flux + solar), l*m
    real( fp_kind ), dimension( size( upwelling_radiance ) ) :: downwelling_radiance



    ! ----------
    ! intrinsics
    ! ----------

    intrinsic adjustl, &
              maxval,  &
              size,    &
              trim



    !#--------------------------------------------------------------------------#
    !#         -- compute the forward radiances and temperatures --             #
    !#--------------------------------------------------------------------------#

    error_status = forward_rtm_rank2( &
                     ! -- forward inputs
                     level_p, layer_p, layer_t, layer_w, layer_o,  &  ! input,  k x m

                     surface_temperature,                          &  ! input,  m
                     surface_emissivity,                           &  ! input,  l*m
                     surface_reflectivity,                         &  ! input,  l*m

                     ! -- other inputs
                     secant_view_angle,                            &  ! input,  m
                     secant_solar_angle,                           &  ! input,  m
                     n_channels_per_profile,                       &  ! input,  m
                     channel_index,                                &  ! input,  l*m

                     ! -- outputs
                     absorber,                                     &  ! output, 0:k x j x m

                     tau_layer_index,                              &  ! output, k x j x m
                     flux_tau_layer_index,                         &  ! output, k x j x m
                     solar_tau_layer_index,                        &  ! output, k x j x m

                     tau_predictor,                                &  ! output, imax x k x m
                     flux_tau_predictor,                           &  ! output, imax x k x m
                     solar_tau_predictor,                          &  ! output, imax x k x m

                     tau,                                          &  ! output, k x l*m
                     flux_tau,                                     &  ! output, k x l*m
                     solar_tau,                                    &  ! output, k x l*m

                     layer_radiance,                               &  ! output, k x l*m
                     downwelling_radiance,                         &  ! output, l*m
                     upwelling_radiance,                           &  ! output, l*m

                     brightness_temperature,                       &  ! output, l*m

                     ! -- optional inputs
                     n_input_profiles = n_input_profiles,          &
                     message_log      = message_log                )


    ! -------------------------------
    ! check for successful completion
    ! -------------------------------

    if ( error_status /= success ) then

      call display_message( routine_name, &
                            'error occured in forward_rtm_rank2', &
                            error_status, &
                            message_log = message_log )
      return

    end if



    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success


  end function compute_rtm_rank2



  function compute_rtm_rank1( &
             ! -- forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k

             surface_temperature,                         &  ! input, scalar
             surface_emissivity,                          &  ! input, l
             surface_reflectivity,                        &  ! input, l

             ! -- other inputs
             secant_view_angle,                           &  ! input, scalar
             secant_solar_angle,                          &  ! input, scalar
             n_channels,                                  &  ! input, scalar
             channel_index,                               &  ! input, l

             ! -- forward output
             tau,                                         &  ! output, k x l
             flux_tau,                                    &  ! output, k x l
             solar_tau,                                   &  ! output, k x l

             upwelling_radiance,                          &  ! output, l
             brightness_temperature,                      &  ! output, l

             ! -- optional inputs
             n_input_profiles,                            &
             message_log )                                &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward inputs
    real( fp_kind ), dimension( : ),    intent( in )  :: level_p                 ! k
    real( fp_kind ), dimension( : ),    intent( in )  :: layer_p                 ! k
    real( fp_kind ), dimension( : ),    intent( in )  :: layer_t                 ! k
    real( fp_kind ), dimension( : ),    intent( in )  :: layer_w                 ! k
    real( fp_kind ), dimension( : ),    intent( in )  :: layer_o                 ! k

    real( fp_kind ),                    intent( in )  :: surface_temperature     ! scalar
    real( fp_kind ), dimension( : ),    intent( in )  :: surface_emissivity      ! l
    real( fp_kind ), dimension( : ),    intent( in )  :: surface_reflectivity    ! l

    ! -- other inputs
    real( fp_kind ),                    intent( in )  :: secant_view_angle       ! scalar
    real( fp_kind ),                    intent( in )  :: secant_solar_angle      ! scalar
    integer,                            intent( in )  :: n_channels              ! scalar
    integer,         dimension( : ),    intent( in )  :: channel_index           ! l

    ! -- forward outputs
    real( fp_kind ), dimension( :, : ), intent( out ) :: tau                     ! k x l
    real( fp_kind ), dimension( :, : ), intent( out ) :: flux_tau                ! k x l
    real( fp_kind ), dimension( :, : ), intent( out ) :: solar_tau               ! k x l

    real( fp_kind ), dimension( : ),    intent( out ) :: upwelling_radiance      ! l
    real( fp_kind ), dimension( : ),    intent( out ) :: brightness_temperature  ! l

    ! -- optional input. note that n_input_profiles is not used in this
    ! -- function. it is included here so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    integer,        optional,           intent( in )  :: n_input_profiles        ! scalar
    character( * ), optional,           intent( in )  :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_rtm_rank1'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- array for integrated absorber amounts, 0:k x j
    real( fp_kind ), dimension( 0:size( layer_p, dim = 1 ), &
                                  max_n_absorbers           ) :: absorber

    ! -- arrays for absorber space indexing, k x j
    integer,         dimension( size( layer_p, dim = 1 ), &
                                max_n_absorbers           ) :: tau_layer_index,      &
                                                               flux_tau_layer_index, &
                                                               solar_tau_layer_index

    ! -- arrays for predictors, imax x k
    real( fp_kind ), dimension( max_n_predictors,         &
                                size( layer_p, dim = 1 )  ) :: tau_predictor,      &
                                                               flux_tau_predictor, &
                                                               solar_tau_predictor

    ! -- array for layer planck radiance term, k x l
    real( fp_kind ), dimension( size( layer_p, dim = 1 ),  &
                                size( upwelling_radiance ) ) :: layer_radiance
      

    ! -- array for downwelling radiance (flux + solar), l
    real( fp_kind ), dimension( size( upwelling_radiance ) ) :: downwelling_radiance



    ! ----------
    ! intrinsics
    ! ----------

    intrinsic adjustl, &
              maxval,  &
              size,    &
              trim



    !#--------------------------------------------------------------------------#
    !#         -- compute the forward radiances and temperatures --             #
    !#--------------------------------------------------------------------------#

    error_status = forward_rtm_rank1( &
                     ! -- forward inputs
                     level_p, layer_p, layer_t, layer_w, layer_o,  &  ! input,  k

                     surface_temperature,                          &  ! input,  scalar
                     surface_emissivity,                           &  ! input,  l
                     surface_reflectivity,                         &  ! input,  l

                     ! -- other inputs
                     secant_view_angle,                            &  ! input,  scalar
                     secant_solar_angle,                           &  ! input,  scalar
                     n_channels,                                   &  ! input,  scalar
                     channel_index,                                &  ! input,  l

                     ! -- outputs
                     absorber,                                     &  ! output, 0:k x j

                     tau_layer_index,                              &  ! output, k x j
                     flux_tau_layer_index,                         &  ! output, k x j
                     solar_tau_layer_index,                        &  ! output, k x j

                     tau_predictor,                                &  ! output, imax x k
                     flux_tau_predictor,                           &  ! output, imax x k
                     solar_tau_predictor,                          &  ! output, imax x k

                     tau,                                          &  ! output, k x l
                     flux_tau,                                     &  ! output, k x l
                     solar_tau,                                    &  ! output, k x l

                     layer_radiance,                               &  ! output, k x l
                     downwelling_radiance,                         &  ! output, l
                     upwelling_radiance,                           &  ! output, l

                     brightness_temperature,                       &  ! output, l

                     message_log = message_log )


    ! -------------------------------
    ! check for successful completion
    ! -------------------------------

    if ( error_status /= success ) then

      call display_message( routine_name, &
                            'error occured in forward_rtm_rank1', &
                            error_status, &
                            message_log = message_log )
      return

    end if



    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success


  end function compute_rtm_rank1







!--------------------------------------------------------------------------------
!s+
! name:
!       forward_rtm
!
! purpose:
!       public function that calculates top-of-atmosphere (toa) radiances
!       and brightness temperatures for an input atmospheric profile or
!       profile set and user specified satellites/channels.
!
! category:
!       ncep rtm
!
! calling sequence:
!
!       result = forward_rtm( &
!                             ! -- inputs
!                             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k x m
!
!                             surface_temperature,                         &  ! input,  m
!                             surface_emissivity,                          &  ! input,  l*m
!                             surface_reflectivity,                        &  ! input,  l*m
!
!                             secant_view_angle,                           &  ! input,  m
!                             secant_solar_angle,                          &  ! input,  m
!                             n_channels_per_profile,                      &  ! input,  m
!                             channel_index,                               &  ! input,  l*m
!
!                             ! -- outputs
!                             absorber,                                    &  ! output, 0:k x j x m
!
!                             tau_layer_index,                             &  ! output, k x j x m
!                             flux_tau_layer_index,                        &  ! output, k x j x m
!                             solar_tau_layer_index,                       &  ! output, k x j x m
!
!                             tau_predictor,                               &  ! output, imax x k x m
!                             flux_tau_predictor,                          &  ! output, imax x k x m
!                             solar_tau_predictor,                         &  ! output, imax x k x m
!
!                             tau,                                         &  ! output, k x l*m
!                             flux_tau,                                    &  ! output, k x l*m
!                             solar_tau,                                   &  ! output, k x l*m
!
!                             layer_radiance,                              &  ! output, k x l*m
!                             downwelling_radiance,                        &  ! output, l*m
!                             upwelling_radiance,                          &  ! output, l*m
!
!                             brightness_temperature,                      &  ! output, l*m
!
!                             ! -- optional inputs
!                             n_input_profiles = n_input_profiles,         &
!                             message_log      = message_log               )
!                             
! input arguments:
!
!       level_p:                 profile set layer interface pressure array. the toa
!                                pressure is not included. toa pressure is parameterised
!                                in the parameters module.
!                                units:      hpa
!                                type:       real
!                                dimension:  k x m; k > 1, and m > or = 1
!                                attributes: intent( in )
!
!       layer_p:                 profile set layer average pressure array.
!                                units:      hpa
!                                type:       real
!                                dimension:  k x m; k > 1, and m > or = 1
!                                attributes: intent( in )
!
!       layer_t:                 profile set layer average temperature array.
!                                units:      kelvin
!                                type:       real
!                                dimension:  k x m; k > 1, and m > or = 1
!                                attributes: intent( in )
!
!       layer_w:      .          profile set layer average water vapor mixing ratio array
!                                units:      g/kg
!                                type:       real
!                                dimension:  k x m; k > 1, and m > or = 1
!                                attributes: intent( in )
!
!       layer_o:                 profile set layer average ozone mixing ratio array.
!                                units:      ppmv
!                                type:       real
!                                dimension:  k x m; k > 1, and m > or = 1
!                                attributes: intent( in )
!
!       surface_temperature:     profile set surface temperature array.
!                                units:      kelvin
!                                type:       real
!                                dimension:  m; m > or = 1 (i.e. scalar)
!                                attributes: intent( in )
!
!       surface_emissivity:      profile set surface emissivity array
!                                units:      none
!                                type:       real
!                                dimension:  l*m; l > 1, and m > or = 1
!                                            nb: this is a 1-d array.
!                                attributes: intent( in )
!
!       surface_reflectivity:    profile set surface reflectivity array
!                                units:      none
!                                type:       real
!                                dimension:  l*m; l > 1, and m > or = 1
!                                            nb: this is a 1-d array.
!                                attributes: intent( in )
!
!       secant_view_angle:       secant of the satellite view angle measured
!                                from nadir for each profile in the set.
!                                units:      none
!                                type:       real
!                                dimension:  m; m > or = 1 (i.e. scalar)
!                                attributes: intent( in )
!
!       secant_solar_angle:      secant of the solar zenith angle for each
!                                profile in the set.
!                                units:      none
!                                type:       real
!                                dimension:  m; m > or = 1 (i.e. scalar)
!                                attributes: intent( in )
!
!       n_channels_per_profile:  the number of channels for each profile in the
!                                set for which radiances are required.
!                                units:      none
!                                type:       integer
!                                dimension:  m; m > or = 1 (i.e. scalar)
!                                attributes: intent( in )
!
!       channel_index:           channel index id array. each element is a unique
!                                index to a (supported) sensor channel.
!                                units:      none
!                                type:       integer
!                                dimension:  l*m; l > 1, and m > or = 1
!                                            nb: this is a 1-d array.
!                                attributes: intent( in )
!
! optional input arguments:
!
!       n_input_profiles:        the number of profiles in the passed arrays to process.
!                                if not specified, the default value is the second dimension
!                                of the pressure array determined using the size intrinsic,
!                                n_profiles. if n_input_profiles is specified and is < 1 or
!                                greater than n_profiles, the default value is set to n_profiles.
!                                this argument is ignored if the input profile arrays are
!                                vectors, i.e. a single profile.
!                                units:      none
!                                type:       integer
!                                dimension:  scalar
!                                attributes: intent( in ), optional
!
!       message_log:             character string specifying a filename in which any
!                                messages will be logged. if not specified, or if an
!                                error occurs opening the log file, the default action
!                                is to output messages to the screen.
!                                units:      none
!                                type:       character
!                                dimension:  scalar
!                                attributes: intent( in ), optional
!
! output arguments:
!
!       absorber:                array of absorber amount for nadir view.
!                                units:      absorber dependent.
!                                type:       real
!                                dimension:  0:k x j x m
!                                attributes: intent( out )
!
!       tau_layer_index:         array of absorber space layer indices of the input
!                                absorber amounts at the satellite view angle.
!                                units:      none.
!                                type:       integer
!                                dimension:  k x j x m
!                                attributes: intent( out )
!
!       flux_tau_layer_index:    array of absorber space layer indices of the input
!                                absorber amounts at the default diffusivity angle.
!                                units:      none.
!                                type:       integer
!                                dimension:  k x j x m
!                                attributes: intent( out )
!
!       solar_tau_layer_index:   array of absorber space layer indices of the input
!                                absorber amounts at the solar zenith angle.
!                                units:      none.
!                                type:       integer
!                                dimension:  k x j x m
!                                attributes: intent( out )
!
!       tau_predictor:           predictor profiles for the layer->toa transmittance.
!                                units:      none.
!                                type:       real
!                                dimension:  i x k x m; k > 1, and m > or = 1
!                                attributes: intent( out )
!
!       flux_tau_predictor:      predictor profiles for the thermal flux transmittance.
!                                units:      none.
!                                type:       real
!                                dimension:  i x k x m; k > 1, and m > or = 1
!                                attributes: intent( out )
!
!       solar_tau_predictor:     predictor profiles for the solar transmittance.
!                                units:      none.
!                                type:       real
!                                dimension:  i x k x m; k > 1, and m > or = 1
!                                attributes: intent( out )
!
!       tau:                     layer->toa transmittance for the satellite
!                                view angle.
!                                units:      none.
!                                type:       real
!                                dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                attributes: intent( out )
!
!       flux_tau:                layer->sfc transmittance for the default
!                                diffusivity angle.
!                                units:      none.
!                                type:       real
!                                dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                attributes: intent( out )
!
!       solar_tau:               layer->sfc transmittance for the solar
!                                zenith angle.
!                                units:      none.
!                                type:       real
!                                dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                attributes: intent( out )
!
!       layer_radiance:          layer planck radiances at every layer for
!                                each channel/profile.
!                                units:      mw/(m^2.sr.cm^-1)
!                                type:       real
!                                dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                attributes: intent( out )
!
!       downwelling_radiance:    toa->sfc radiances for each channel/profile due
!                                to thermal flux and solar components.
!                                units:      mw/(m^2.sr.cm^-1)
!                                type:       real
!                                dimension:  l*m; l > 1, and m > or = 1
!                                            nb: this is a 1-d array.
!                                attributes: intent( out )
!
!       upwelling_radiance:      toa radiances for each channel/profile.
!                                units:      mw/(m^2.sr.cm^-1)
!                                type:       real
!                                dimension:  l*m; l > 1, and m > or = 1
!                                            nb: this is a 1-d array.
!                                attributes: intent( out )
!
!       brightness_temperature:  temperatures corresponding to the toa radiances
!                                for each channel/profile.
!                                units:      kelvin
!                                type:       real
!                                dimension:  l*m; l > 1, and m > or = 1
!                                            nb: this is a 1-d array.
!                                attributes: intent( out )
!
! optional output arguments:
!       none.
!
! function result:
!       result = success => calculation was successful
!              = failure => error occurred checking input arguments
!
! calls:
!      display_message:            subroutine to output messages
!                                  source: error_handler module
!
!      get_max_n_channels:         routine to retrieve the value of the
!                                  max_n_channels "pseudo-parameter".
!                                  source: parameters module
!
!      compute_absorber_amount:    subroutine to integrate the absorber profiles
!                                  source: absorber_profile module
!
!      find_absorber_layer_index:  subroutine to find the absorber space levels
!                                  that bracket the calculated absorber amounts
!                                  source: absorber_profile module
!                                  
!      compute_predictors:         subroutine to compute the transmittance predictor
!                                  profiles.
!                                  source: predictor module
!
!      compute_transmittance:      subroutine to compute the transmittance profiles.
!                                  source: transmittance module
!
!      compute_radiance:           subroutine to compute the toa radiances and
!                                  brightness temperatures.
!                                  source: radiance module
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
!       note the restrictions on the input array dimensions:
!         k == n_layers > 1
!         l == n_channels > 1
!         m == n_profiles > or = 1
!
! procedure:
!       see individual module function documentation.
!s-
!--------------------------------------------------------------------------------

  function forward_rtm_rank2( &

             ! -- inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k x m

             surface_temperature,                         &  ! input,  m
             surface_emissivity,                          &  ! input,  l*m
             surface_reflectivity,                        &  ! input,  l*m

             secant_view_angle,                           &  ! input,  m
             secant_solar_angle,                          &  ! input,  m
             n_channels_per_profile,                      &  ! input,  m
             channel_index,                               &  ! input,  l*m

             ! -- outputs
             absorber,                                    &  ! output, 0:k x j x m

             tau_layer_index,                             &  ! output, k x j x m
             flux_tau_layer_index,                        &  ! output, k x j x m
             solar_tau_layer_index,                       &  ! output, k x j x m

             tau_predictor,                               &  ! output, imax x k x m
             flux_tau_predictor,                          &  ! output, imax x k x m
             solar_tau_predictor,                         &  ! output, imax x k x m

             tau,                                         &  ! output, k x l*m
             flux_tau,                                    &  ! output, k x l*m
             solar_tau,                                   &  ! output, k x l*m

             layer_radiance,                              &  ! output, k x l*m
             downwelling_radiance,                        &  ! output, l*m
             upwelling_radiance,                          &  ! output, l*m

             brightness_temperature,                      &  ! output, l*m

             ! -- optional inputs
             n_input_profiles,                            &
             message_log )                                &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- inputs
    real( fp_kind ), dimension( :, : ),     intent( in )  :: level_p                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_p                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_t                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_w                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_o                 ! k x m

    real( fp_kind ), dimension( : ),        intent( in )  :: surface_temperature     ! m
    real( fp_kind ), dimension( : ),        intent( in )  :: surface_emissivity      ! l*m
    real( fp_kind ), dimension( : ),        intent( in )  :: surface_reflectivity    ! l*m

    real( fp_kind ), dimension( : ),        intent( in )  :: secant_view_angle       ! m
    real( fp_kind ), dimension( : ),        intent( in )  :: secant_solar_angle      ! m
    integer,         dimension( : ),        intent( in )  :: n_channels_per_profile  ! m
    integer,         dimension( : ),        intent( in )  :: channel_index           ! l*m

    ! -- outputs
    real( fp_kind ), dimension( 0:, :, : ), intent( out ) :: absorber                ! 0:k x j x m

    integer,         dimension(  :, :, : ), intent( out ) :: tau_layer_index         ! k x j x m
    integer,         dimension(  :, :, : ), intent( out ) :: flux_tau_layer_index    ! k x j x m
    integer,         dimension(  :, :, : ), intent( out ) :: solar_tau_layer_index   ! k x j x m

    real( fp_kind ), dimension( :, :, : ),  intent( out ) :: tau_predictor           ! imax x k x m
    real( fp_kind ), dimension( :, :, : ),  intent( out ) :: flux_tau_predictor      ! imax x k x m
    real( fp_kind ), dimension( :, :, : ),  intent( out ) :: solar_tau_predictor     ! imax x k x m

    real( fp_kind ), dimension( :, : ),     intent( out ) :: tau                     ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( out ) :: flux_tau                ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( out ) :: solar_tau               ! k x l*m

    real( fp_kind ), dimension( :, : ),     intent( out ) :: layer_radiance          ! k x l*m
    real( fp_kind ), dimension( : ),        intent( out ) :: downwelling_radiance    ! l*m
    real( fp_kind ), dimension( : ),        intent( out ) :: upwelling_radiance      ! l*m

    real( fp_kind ), dimension( : ),        intent( out ) :: brightness_temperature  ! l*m

    ! -- optional input
    integer,        optional,               intent( in )  :: n_input_profiles        ! scalar
    character( * ), optional,               intent( in )  :: message_log
    

    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'forward_rtm_rank2'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    character( 100 ) :: message
    character( 10 )  :: value_in, value_allowed

    integer :: m, n_profiles          ! profile loop variables
    integer :: l, l1, l2, n_channels  ! channel loop/index variables

    integer :: valid_solar

    ! -- maximum channels pseudo parameter
    integer :: max_n_channels
    logical :: is_set


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic adjustl, &
              any,     &
              cos,     &
              maxval,  &
              present, &
              size,    &
              trim



    !#--------------------------------------------------------------------------#
    !#          -- determine array dimensions and check the values --           #
    !#--------------------------------------------------------------------------#

    ! ----------------------------
    ! check the number of profiles
    ! ----------------------------

    ! -- number of atmospheric profiles. default size is the second
    ! -- dimension of the layer_p argument.
    n_profiles = size( layer_p, dim = 2 )
    if ( present( n_input_profiles ) ) then
      if ( n_input_profiles > 0 .and. n_input_profiles <= n_profiles ) then
        n_profiles = n_input_profiles
      else
        write( message, '( "invalid n_input_profiles value: ", i5, &
                          &". using pressure array dimension value of ", i5, "." )' ) &
                        n_input_profiles, n_profiles
        call display_message( routine_name,    &
                              trim( message ), &
                              warning,         &
                              message_log = message_log )
      end if
    end if

    ! -- check that the number of profiles is not greater than
    ! -- max_n_profiles. this is simply a limit to restrict the
    ! -- size of the input arrays so they're not too big.
    if ( n_profiles > max_n_profiles ) then

      error_status = failure
      write( value_in,      '( i5 )' ) n_profiles
      write( value_allowed, '( i5 )' ) max_n_profiles
      call display_message( routine_name, &
                            'number of passed profiles ('// &
                            trim( adjustl( value_in ) )// &
                            ') > maximum number of profiles allowed ('// &
                            trim( adjustl( value_allowed ) )//').', &
                            error_status, &
                            message_log = message_log )
      return

    end if


    ! ----------------------------
    ! check the number of channels
    ! ----------------------------

    ! -- check for a negative number of channels
    if ( any( n_channels_per_profile < 0 ) ) then
      error_status = failure
      call display_message( routine_name, &
                            'elements of n_channels_per_profile are negative.', &
                            error_status, &
                            message_log = message_log )
      return
    end if

    ! -- check that the maximum number of channels for any profile
    ! -- is not greater than than the number of channels with which
    ! -- the model was initialised
    n_channels = maxval( n_channels_per_profile )

    call get_max_n_channels( max_n_channels, is_set )

    if ( .not. is_set ) then
      error_status = failure
      call display_message( routine_name, &
                            'max_n_channels value not set. check that rtm is initialised.', &
                            error_status, &
                            message_log = message_log )
      return
    end if

    if ( n_channels > max_n_channels ) then
      error_status = failure
      write( value_in,      '( i5 )' ) n_channels
      write( value_allowed, '( i5 )' ) max_n_channels
      call display_message( routine_name, &
                            'number of requested channels ('// &
                            trim( adjustl( value_in ) )// &
                            ') > number of initialisation channels ('// &
                            trim( adjustl( value_allowed ) )//').', &
                            error_status, &
                            message_log = message_log )
      return
    end if




    !#--------------------------------------------------------------------------#
    !#                           -- profile loop --                             #
    !#--------------------------------------------------------------------------#

    ! -------------------------------
    ! initialise channel index finder
    ! -------------------------------

    l1 = 1


    ! ------------------
    ! begin profile loop
    ! ------------------

    m_profile_loop: do m = 1, n_profiles


      ! -----------------------------------------------
      ! check for the "no channel" case. if no channels
      ! required for this profile, go to the next one.
      ! -----------------------------------------------

      if ( n_channels_per_profile( m ) == 0 ) cycle m_profile_loop


      ! -------------------------------------------
      ! determine the end channel index index range
      ! -------------------------------------------

      l2 = l1 + n_channels_per_profile( m ) - 1


      ! ------------------------
      ! call the rank-1 function
      ! ------------------------

      error_status = forward_rtm_rank1( &

                       ! -- forward inputs
                       level_p( :, m ),                   &  ! input,  k
                       layer_p( :, m ),                   &  ! input,  k
                       layer_t( :, m ),                   &  ! input,  k
                       layer_w( :, m ),                   &  ! input,  k
                       layer_o( :, m ),                   &  ! input,  k

                       surface_temperature( m ),          &  ! input,  scalar
                       surface_emissivity( l1:l2 ),       &  ! input,  l
                       surface_reflectivity( l1:l2 ),     &  ! input,  l

                       ! -- other inputs
                       secant_view_angle( m ),            &  ! input,  scalar
                       secant_solar_angle( m ),           &  ! input,  scalar
                       n_channels_per_profile( m ),       &  ! input,  scalar
                       channel_index( l1:l2 ),            &  ! input,  l

                       ! -- outputs
                       absorber( 0:, :, m ),              &  ! output, 0:k x j

                       tau_layer_index( :, :, m ),        &  ! output, k x j
                       flux_tau_layer_index( :, :, m ),   &  ! output, k x j
                       solar_tau_layer_index( :, :, m ),  &  ! output, k x j

                       tau_predictor( :, :, m ),          &  ! output, imax x k
                       flux_tau_predictor( :, :, m ),     &  ! output, imax x k
                       solar_tau_predictor( :, :, m ),    &  ! output, imax x k

                       tau( :, l1:l2 ),                   &  ! output, k x l
                       flux_tau( :, l1:l2 ),              &  ! output, k x l
                       solar_tau( :, l1:l2 ),             &  ! output, k x l

                       layer_radiance( :, l1:l2 ),        &  ! output, k x l
                       downwelling_radiance( l1:l2 ),     &  ! output, l
                       upwelling_radiance( l1:l2 ),       &  ! output, l

                       brightness_temperature( l1:l2 ),   &  ! output, l

                       message_log = message_log )


      ! -------------------------------
      ! check for successful completion
      ! -------------------------------

      if ( error_status /= success ) then

        call display_message( routine_name, &
                              'error occured in forward_rtm_rank1', &
                              error_status, &
                              message_log = message_log )
        return

      end if


      ! ------------------------------
      ! update the channel begin index
      ! ------------------------------

      l1 = l2 + 1

    end do m_profile_loop



    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function forward_rtm_rank2


  function forward_rtm_rank1( &

             ! -- inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k

             surface_temperature,                         &  ! input,  scalar
             surface_emissivity,                          &  ! input,  l
             surface_reflectivity,                        &  ! input,  l

             secant_view_angle,                           &  ! input,  scalar
             secant_solar_angle,                          &  ! input,  scalar
             n_channels,                                  &  ! input,  scalar
             channel_index,                               &  ! input,  l

             ! -- outputs
             absorber,                                    &  ! output, 0:k x j

             tau_layer_index,                             &  ! output, k x j
             flux_tau_layer_index,                        &  ! output, k x j
             solar_tau_layer_index,                       &  ! output, k x j

             tau_predictor,                               &  ! output, imax x k
             flux_tau_predictor,                          &  ! output, imax x k
             solar_tau_predictor,                         &  ! output, imax x k

             tau,                                         &  ! output, k x l
             flux_tau,                                    &  ! output, k x l
             solar_tau,                                   &  ! output, k x l

             layer_radiance,                              &  ! output, k x l
             downwelling_radiance,                        &  ! output, l
             upwelling_radiance,                          &  ! output, l

             brightness_temperature,                      &  ! output, l

             ! -- optional inputs
             n_input_profiles,                            &
             message_log )                                &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- inputs
    real( fp_kind ), dimension( : ),     intent( in )  :: level_p                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_p                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_t                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_w                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_o                 ! k

    real( fp_kind ),                     intent( in )  :: surface_temperature     ! scalar
    real( fp_kind ), dimension( : ),     intent( in )  :: surface_emissivity      ! l
    real( fp_kind ), dimension( : ),     intent( in )  :: surface_reflectivity    ! l

    real( fp_kind ),                     intent( in )  :: secant_view_angle       ! scalar
    real( fp_kind ),                     intent( in )  :: secant_solar_angle      ! scalar
    integer,                             intent( in )  :: n_channels              ! scalar
    integer,         dimension( : ),     intent( in )  :: channel_index           ! l

    ! -- outputs
    real( fp_kind ), dimension( 0:, : ), intent( out ) :: absorber                ! 0:k x j

    integer,         dimension(  :, : ), intent( out ) :: tau_layer_index         ! k x j
    integer,         dimension(  :, : ), intent( out ) :: flux_tau_layer_index    ! k x j
    integer,         dimension(  :, : ), intent( out ) :: solar_tau_layer_index   ! k x j

    real( fp_kind ), dimension( :, : ),  intent( out ) :: tau_predictor           ! imax x k
    real( fp_kind ), dimension( :, : ),  intent( out ) :: flux_tau_predictor      ! imax x k
    real( fp_kind ), dimension( :, : ),  intent( out ) :: solar_tau_predictor     ! imax x k

    real( fp_kind ), dimension( :, : ),  intent( out ) :: tau                     ! k x l
    real( fp_kind ), dimension( :, : ),  intent( out ) :: flux_tau                ! k x l
    real( fp_kind ), dimension( :, : ),  intent( out ) :: solar_tau               ! k x l

    real( fp_kind ), dimension( :, : ),  intent( out ) :: layer_radiance          ! k x l
    real( fp_kind ), dimension( : ),     intent( out ) :: downwelling_radiance    ! l
    real( fp_kind ), dimension( : ),     intent( out ) :: upwelling_radiance      ! l

    real( fp_kind ), dimension( : ),     intent( out ) :: brightness_temperature  ! l

    ! -- optional input. note that n_input_profiles is not used in this
    ! -- function. it is specified so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    integer,        optional,            intent( in )  :: n_input_profiles        ! scalar
    character( * ), optional,            intent( in )  :: message_log
    

    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'forward_rtm_rank1'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    character( 100 ) :: message
    character( 10 )  :: value_in, value_allowed

    integer :: n_layers   ! layer dimension
    integer :: l          ! channel loop/index variables

    integer :: valid_solar

    ! -- maximum channels pseudo parameter
    integer :: max_n_channels
    logical :: is_set

    ! -- arrays for integrated absorber amounts.
    real( fp_kind ), dimension( 0:size( absorber, dim = 1 )-1, &
                                  size( absorber, dim = 2 )    ) :: tau_absorber,      &
                                                                    flux_tau_absorber, &
                                                                    solar_tau_absorber

    ! ----------
    ! intrinsics
    ! ----------

    intrinsic adjustl, &
              any,     &
              cos,     &
              maxval,  &
              size,    &
              trim



    !#--------------------------------------------------------------------------#
    !#           -- determine array dimensions and check input --               #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! check the number of channels - if zero
    ! then simply return.
    ! --------------------------------------

    if ( n_channels == 0 ) return


    ! ------------------
    ! get the dimensions
    ! ------------------

    n_layers = size( layer_p )


    ! ------------------------------------------------------------
    ! check that the number of layers is not greater than
    ! max_n_absorber_layers. if there are more input layers
    ! there may be non-unique assignation of absorber->atmospheric
    ! layers
    ! ------------------------------------------------------------

    if ( n_layers > max_n_absorber_layers ) then

      error_status = failure
      write( value_in,      '( i5 )' ) n_layers
      write( value_allowed, '( i5 )' ) max_n_absorber_layers
      call display_message( routine_name, &
                            'number of passed layers ('// &
                            trim( adjustl( value_in ) )// &
                            ') > maximum number of absorber layers ('// &
                            trim( adjustl( value_allowed ) )//').', &
                            error_status, &
                            message_log = message_log )
      return

    end if


    ! -----------------------------------
    ! perform a simple check on the input
    ! data for negative values
    ! -----------------------------------

    ! -- profile data
    if ( any( level_p < zero ) .or. &
         any( layer_p < zero ) .or. &
         any( layer_t < zero ) .or. &
         any( layer_w < zero ) .or. &
         any( layer_o < zero )      ) then

       write (6,*) 'level_p', level_p
       write (6,*) 'layer_p', layer_p
       write (6,*) 'layer_t', layer_t
       write (6,*) 'layer_w', layer_w
       write (6,*) 'layer_o', layer_o
       
      error_status = failure
      call display_message( routine_name, &
                            'negative values found in input profile data.', &
                            error_status, &
                            message_log = message_log )
      
    end if

    ! -- surface properties
    if (      surface_temperature  < zero   .or. &
         any( surface_emissivity   < zero ) .or. &
         any( surface_reflectivity < zero )      ) then
      error_status = failure
      call display_message( routine_name, &
                            'negative values found in surface properties (tsfc,esfc,rsfc).', &
                            error_status, &
                            message_log = message_log )
      return
    end if

    ! -- number of channels
    if ( n_channels < 0 ) then
      error_status = failure
      call display_message( routine_name, &
                            'negative number of channels passed..', &
                            error_status, &
                            message_log = message_log )
      return
    end if



    !#--------------------------------------------------------------------------#
    !#           -- calculate the profile generic absorber amounts --           #
    !#--------------------------------------------------------------------------#

    call compute_absorber_amount( level_p( : ),     &  ! input,  k
                                  layer_w( : ),     &  ! input,  k
                                  layer_o( : ),     &  ! input,  k

                                  absorber( 0:, : ) )  ! output, 0:k x j



    !#--------------------------------------------------------------------------#
    !#      -- calculate the predictors for the upwelling transmittance --      #
    !#--------------------------------------------------------------------------#

    ! ------------------------------------------------------
    ! modify absorber quantities by the angle secant
    ! could put a loop here but here's hoping the compiler
    ! recognises this as a group of loops over layer.
    ! ------------------------------------------------------

    tau_absorber( 0:, : ) = secant_view_angle * absorber( 0:, : )

    ! -- subtract the top pressure from the dry absorber profile
    tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - toa_pressure


    ! -----------------------------------------------------
    ! calculate the predictors for the satellite view angle
    ! -----------------------------------------------------

    call compute_predictors( layer_p( : ),          &  ! input,  k
                             layer_t( : ),          &  ! input,  k
                             layer_w( : ),          &  ! input,  k
                             tau_absorber( 0:, : ), &  ! input,  0:k x j

                             tau_predictor( :, : )  )  ! output, i x k


    ! ------------------------------------------------
    ! determine the absorber space levels that bracket
    ! the "average absorber" amounts at the view angle
    ! ------------------------------------------------

    call find_absorber_layer_index( tau_absorber( 0:, : ),  &  ! input, 0:k x j
                                    tau_layer_index( :, : ) )  ! output, k x j



    !#--------------------------------------------------------------------------#
    !#        -- calculate the predictors for the flux transmittance --         #
    !#--------------------------------------------------------------------------#

    ! ---------------------------------------------
    ! have any infrared channels been specified for
    ! the current profile? (microwave channels are
    ! flagged as == 1, so ir == 0).
    !
    ! for microwave channels the downwelling flux
    ! transmission is assumed == upwelling view
    ! angle transmission.
    ! ---------------------------------------------

    if ( any( is_microwave_channel( channel_index( 1:n_channels ) ) == 0 ) ) then


      ! ---------------------------------
      ! modify the nadir absorber amounts
      ! ---------------------------------

      flux_tau_absorber( 0:, : ) = secant_diffusivity_angle * absorber( 0:, : )
      
      ! -- subtract the top pressure from the dry absorber profile
      flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - toa_pressure


      ! --------------------------------
      ! calculate the predictors for the
      ! diffusivity angle
      ! --------------------------------

      ! -- calculate the integrated predictors only
      call compute_predictors( layer_p( : ),               &  ! input,  k
                               layer_t( : ),               &  ! input,  k
                               layer_w( : ),               &  ! input,  k
                               flux_tau_absorber( 0:, : ), &  ! input,  0:k x j

                               flux_tau_predictor( :, : ), &  ! output, i x k

                               no_standard = 1             )  ! optional input

      ! -- copy the angle independent (standard) predictors
      flux_tau_predictor( 1:max_n_standard_predictors, : ) = &
           tau_predictor( 1:max_n_standard_predictors, : ) 


      ! ------------------------------------------------
      ! determine the absorber space levels that bracket
      ! the absorber amounts at the diffusivity angle
      ! ------------------------------------------------

      call find_absorber_layer_index( flux_tau_absorber( 0:, : ),  &  ! input,  0:k x j
                                      flux_tau_layer_index( :, : ) )  ! output, k x j

    end if
     


    !#--------------------------------------------------------------------------#
    !#       -- calculate the predictors for the solar transmittance --         #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------------------
    ! have *any* solar sensitive channels been specified
    ! for the current profile (flagged as == 1)?
    !
    ! and
    !
    ! is the specified solar zenith angle valid?
    ! --------------------------------------------------

    if ( ( any( is_solar_channel( channel_index( 1:n_channels ) ) == 1 ) ) .and. &
         secant_solar_angle < max_secant_solar_angle ) then


      ! --------------------------------
      ! modify the nadir absorber amount
      ! --------------------------------

      solar_tau_absorber( 0:, : ) = secant_solar_angle * absorber( 0:, : )

      ! -- subtract the top pressure from the dry absorber profile
      solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - toa_pressure


      ! --------------------------------
      ! calculate the predictors for the
      ! solar zenith angle
      ! --------------------------------

      ! -- calculate the integrated predictors only
      call compute_predictors( layer_p( : ),                &  ! input,  k
                               layer_t( : ),                &  ! input,  k
                               layer_w( : ),                &  ! input,  k
                               solar_tau_absorber( 0:, : ), &  ! input,  0:k x j

                               solar_tau_predictor( :, : ), &  ! output, i x k

                               no_standard = 1              )  ! optional input

      ! -- copy the angle independent predictors
      solar_tau_predictor( 1:max_n_standard_predictors, : ) = &
            tau_predictor( 1:max_n_standard_predictors, : ) 


      ! ------------------------------------------------
      ! determine the absorber space levels that bracket
      ! the absorber amounts at the solar zenith angle
      ! ------------------------------------------------

      call find_absorber_layer_index( solar_tau_absorber( 0:, : ),  &  ! input, 0:k x j
                                      solar_tau_layer_index( :, : ) )  ! output, k x j

    end if
     


    !#--------------------------------------------------------------------------#
    !#                           -- channel loop --                             #
    !#--------------------------------------------------------------------------#

    l_channel_loop: do l = 1, n_channels


      ! --------------------------------------------------
      ! calculate the current channel layer transmittances
      ! for the satellite view angle
      ! --------------------------------------------------

      call compute_transmittance( tau_absorber( 0:, : ),   &   ! input, 0:k x j
                                  tau_predictor( :, : ),   &   ! input, i x k
                                  tau_layer_index( :, : ), &   ! input, k x j
                                  channel_index( l ),      &   ! input, scalar
                                  up,                      &   ! input, scalar

                                  tau( :, l )              )   ! output, k


      ! -----------------------------------------------------
      ! if the current channel is an infrared channel,
      ! then calculate the downwelling flux transmittance.
      !
      ! if the current channel is a microwave channel,
      ! then use the predictors for the upwelling
      ! transmittance calculations.
      !
      ! two things:
      ! - currently, the predictors and coefficients are
      !   the same for the up- and downwelling cases,
      !   hence the simple "upending" of the transmittances
      !   for the microwave (assumed specular) case.
      ! - the "upending", or assumption that the downwelling
      !   transmittances can be derived directly from the
      !   upwelling transmittances, is only valid for
      !   monochromatic transmittances. for broadband sensors
      !   this is an approximation - and depending on the
      !   spectral structure of the transmittances within
      !   a channel's response - not usually a good one.
      ! -----------------------------------------------------

      if ( is_microwave_channel( channel_index( l ) ) == 0 ) then

        ! -- ir channel
        call compute_transmittance( flux_tau_absorber( 0:, : ),   &   ! input, 0:k x j
                                    flux_tau_predictor( :, : ),   &   ! input, i x k
                                    flux_tau_layer_index( :, : ), &   ! input, k x j
                                    channel_index( l ),           &   ! input, scalar
                                    down,                         &   ! input, scalar

                                    flux_tau( :, l )              )   ! output, k

      else


        ! -- uw channel

!  this can be considered a "hook" for future versions where
!  downwelling will not be derived from the upwelling.
!
!        call compute_transmittance( tau_absorber( 0:, : ),   &   ! input, 0:k x j
!                                    tau_predictor( :, : ),   &   ! input, i x k
!                                    tau_layer_index( :, : ), &   ! input, k x j
!                                    channel_index( l ),      &   ! input, scalar
!                                    down,                    &   ! input, scalar
!
!                                    flux_tau( :, l )         )   ! output, k

        ! -- this gives the identical result as a separate call
        ! -- but without the extra computational burden.
        flux_tau( :, l ) = zero
        if ( tau( n_layers, l ) > tolerance ) then
          flux_tau( 1, l ) = tau( n_layers, l )
          flux_tau( 2:n_layers, l ) =     flux_tau( 1, l ) / &
          !                           ------------------------
                                       tau( 1:n_layers-1, l )
        end if

      end if



      ! ----------------------------------------------------
      ! if the current channel is a solar sensitive channel,
      !   and
      ! the solar angle is valid, then calculate the
      ! transmittance for direct solar.
      ! ----------------------------------------------------

      if ( is_solar_channel( channel_index( l ) ) == 1 .and. &
           secant_solar_angle < max_secant_solar_angle ) then

        valid_solar = 1

        call compute_transmittance( solar_tau_absorber( 0:, : ),   &   ! input, 0:k x j
                                    solar_tau_predictor( :, : ),   &   ! input, i x k
                                    solar_tau_layer_index( :, : ), &   ! input, k x j
                                    channel_index( l ),            &   ! input, scalar
                                    down,                          &   ! input, scalar

                                    solar_tau( :, l )              )   ! output, k

      else

        valid_solar = 0

        solar_tau( :, l ) = zero

      end if


      ! ---------------------------------------
      ! calculate the profile/channel radiances
      ! ---------------------------------------

      call compute_radiance( layer_t( : ),               &  ! input, k

                             surface_temperature,        &  ! input, scalar
                             surface_emissivity( l ),    &  ! input, scalar
                             surface_reflectivity( l ),  &  ! input, scalar

                             tau(      :, l ),           &  ! input, k
                             flux_tau( :, l ),           &  ! input, k
                             solar_tau( n_layers, l ),   &  ! input, scalar

                             secant_solar_angle,         &  ! input, scalar
                             valid_solar,                &  ! input, scalar
                             channel_index( l ),         &  ! input, scalar

                             layer_radiance( :, l ),     &  ! output, k
                             downwelling_radiance( l ),  &  ! output, scalar
                             upwelling_radiance( l ),    &  ! output, scalar

                             brightness_temperature( l ) )  ! output, scalar

    end do l_channel_loop



    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function forward_rtm_rank1

end module forward_model


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
! $name$
!
! $state$
!
! $log$
! revision 2.12  2001/11/07 15:03:15  paulv
! - added check for negative number of channels to forward_rtm_rank2() function.
! - added profile loop cycle statement to forward_rtm_rank2() function.
! - added check for negative number of channels to forward_rtm_rank1() function.
! - changed the logical if test for the number of input profiles in the
!   forward_rtm_rank2() function from
!     if ( n_input_profiles >= 1 .and. n_input_profiles <= n_profiles ) then
!   to
!     if ( n_input_profiles > 0 .and. n_input_profiles <= n_profiles ) then
!   the use of both the ">=" and "<=" realtional operators with the .and. i
!   found confusing.
!
! revision 2.11  2001/09/28 22:44:24  paulv
! - overloaded the compute_rtm and forward_rtm functions to accept both a
!   single or group of profiles. contained private functions are now
!     o compute_rtm_rank1 and forward_rtm_rank1 for single profile input
!     o compute_rtm_rank2 and forward_rtm_rank2 for multiple profile input
! - put n_input_profiles optional argument back in the compute_rtm and forward_rtm
!   argument lists. this allows more control over how many profiles are to be
!   processed rather than simply relying on the dimension of the input arrays.
!   now, as before,
!     n_profiles = size( layer_p, dim = 2 )
!   but also,
!     if ( present( n_input_profiles ) ) then
!       if ( n_input_profiles >= 1 .and. n_input_profiles <= n_profiles ) then
!         n_profiles = n_input_profiles
!     ....
! - changed surface_temperature argument check from
!     if ( any( surface_temperature < zero ) )
!   to
!     if (      surface_temperature < zero )
!   as the check is now done in the rank-1 forward_rtm function. this eliminates
!   the need to fully populate the input arrays with data when only certain
!   "chunks" may be processed. previously, the use of any() could generate
!   and error if the full surface_temperature array was not initialised.
! - added "name" to rcs keyword list.
!
! revision 2.10  2001/09/04 21:29:11  paulv
! - updated documentation.
!
! revision 2.9  2001/08/31 21:22:23  paulv
! - removed input data checks from compute_rtm. the same checks are performed
!   in the main routine, forward_rtm, so there was no need to replicate them.
! - added check for negative profile and surface data in forward_rtm.
! - maximum solar angle secant is no longer calculated in forward_rtm but
!   is declared as a parameter in the parameters module.
!
! revision 2.8  2001/08/16 16:36:29  paulv
! - updated documentation.
! - the channel dimension is now obtained by:
!     n_channels = maxval( n_channels_per_profile )
!   rather than
!     n_channels = size( channel_index ) / n_profiles
!   since the latter assumes the same number of channels will be processed
!   for each profile - which may not be the case. the new method determines
!   the largest number of channels to be processed for any particular
!   profile.
! - the comparison of n_channels and max_n_channels is now done via the
!   max_n_channels methods in the parameters module. and extra check to
!   see if max_n_channels has been set was added.
!
! revision 2.7  2001/08/01 16:47:59  paulv
! - removed use of absorber_space module.
! - only clauses added to other use statements so only those module members
!   required in this module are visible.
! - added compute_rtm function. this is a wrapper for the forward_rtm function.
! - updated input argument checking. now consistent with other model
!   components.
!
! revision 2.6  2001/07/12 18:41:37  paulv
! - commented out informational message output at start of function.
!
! revision 2.5  2001/05/29 18:21:02  paulv
! - now use type_kinds module parameter fp_kind to set the floating point
!   data type.
! - added absorber_space and predictors module use to provide access to
!   the absorber space definitions and the predictor calculation routines.
! - all angle parameters placed in the parameters module.
! - argument list increased to provide data for other component rtm calls
!   (e.g. tl or ad).
! - for full listing of differences, do a cvs diff between revisions 2.4
!   and 2.5 (this one).
!
! revision 2.4  2001/01/24 20:14:21  paulv
! - latest test versions.
!
! revision 2.3  2000/11/09 20:58:28  paulv
! - made radiance function call consistent with changes to that module.
!   (specifically regarding the solar, flux, and reflectivity terms).
!
! revision 2.2  2000/08/31 19:36:32  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 2.1  2000/08/24 18:08:28  paulv
! - many changes from initial version. too many to list here.
! - updated module and subprogram documentation.
!
!
!
!
!
