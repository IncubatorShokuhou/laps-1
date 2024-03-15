!------------------------------------------------------------------------------
!m+
! name:
!       tangent_linear_model
!
! purpose:
!       module containing the ncep rt tangent linear model function.
!
! category:
!       ncep rtm
!
! calling sequence:
!       use tangent_linear_model
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
!       tangent_linear_rtm:    public function that calculates top-of-atmosphere (toa)
!                              tangent-linear radiances and brightness temperatures for
!                              user specified profiles and satellites/channels. 
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

module tangent_linear_model


  ! ------------
  ! module usage
  ! ------------

  use type_kinds, only : fp_kind
  use error_handler
  use parameters
  use spectral_coefficients
  use absorber_profile
  use predictors
  use transmittance
  use radiance


  ! -----------------------
  ! disable implicit typing
  ! -----------------------

  implicit none


  ! --------------------
  ! default visibilities
  ! --------------------

  private
  public :: tangent_linear_rtm


  ! --------------------
  ! function overloading
  ! --------------------

  interface tangent_linear_rtm
    module procedure tangent_linear_rtm_rank1
    module procedure tangent_linear_rtm_rank2
  end interface ! tangent_linear_rtm


contains


!--------------------------------------------------------------------------------
!s+
! name:
!       tangent_linear_rtm
!
! purpose:
!       public function that calculates top-of-atmosphere (toa) tangent-linear
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
! category:
!       ncep rtm
!
! calling sequence:
!
!       result = tangent_linear_rtm( &
!                                    ! forward inputs
!                                    level_p, layer_p, layer_t, layer_w, layer_o,                &  ! input, k x m
!
!                                    surface_temperature,                                        &  ! input, m
!                                    surface_emissivity,                                         &  ! input, l*m
!                                    surface_reflectivity,                                       &  ! input, l*m
!
!                                    absorber,                                                   &  ! input, 0:k x j x m
!
!                                    tau_layer_index,                                            &  ! input, k x j x m
!                                    flux_tau_layer_index,                                       &  ! input, k x j x m
!                                    solar_tau_layer_index,                                      &  ! input, k x j x m
!
!                                    tau_predictor,                                              &  ! input, imax x k x m
!                                    flux_tau_predictor,                                         &  ! input, imax x k x m
!                                    solar_tau_predictor,                                        &  ! input, imax x k x m
!
!                                    tau,                                                        &  ! input, k x l*m
!                                    flux_tau,                                                   &  ! input, k x l*m
!                                    solar_tau,                                                  &  ! input, k x l*m
!
!                                    layer_radiance,                                             &  ! input, k x l*m
!                                    downwelling_radiance,                                       &  ! input, l*m
!                                    upwelling_radiance,                                         &  ! input, l*m
!
!                                    ! -- tangent-linear inputs
!                                    level_p_tl, layer_p_tl, layer_t_tl, layer_w_tl, layer_o_tl, &  ! input, k x m
!
!                                    surface_temperature_tl,                                     &  ! input, m
!                                    surface_emissivity_tl,                                      &  ! input, l*m
!                                    surface_reflectivity_tl,                                    &  ! input, l*m
!
!                                    ! -- other inputs
!                                    secant_view_angle,                                          &  ! input, m
!                                    secant_solar_angle,                                         &  ! input, m
!                                    n_channels_per_profile,                                     &  ! input, m
!                                    channel_index,                                              &  ! input, l*m
!
!                                    ! -- tangent-linear outputs
!                                    tau_tl,                                                     &  ! output, k x l*m
!                                    flux_tau_tl,                                                &  ! output, k x l*m
!                                    solar_tau_tl,                                               &  ! output, k x l*m
!
!                                    layer_radiance_tl,                                          &  ! output, k x l*m
!                                    downwelling_radiance_tl,                                    &  ! output, l*m
!                                    upwelling_radiance_tl,                                      &  ! output, l*m
!
!                                    brightness_temperature_tl,                                  &  ! output, l*m
!
!                                    ! optional inputs
!                                    message_log = message_log )
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
!       absorber:                  array of absorber amount for nadir view.
!                                  units:      absorber dependent.
!                                  type:       real
!                                  dimension:  0:k x j x m
!                                  attributes: intent( in )
!
!       tau_layer_index:           array of absorber space layer indices of the input
!                                  absorber amounts at the satellite view angle.
!                                  units:      none.
!                                  type:       integer
!                                  dimension:  k x j x m
!                                  attributes: intent( in )
!
!       flux_tau_layer_index:      array of absorber space layer indices of the input
!                                  absorber amounts at the default diffusivity angle.
!                                  units:      none.
!                                  type:       integer
!                                  dimension:  k x j x m
!                                  attributes: intent( in )
!
!       solar_tau_layer_index:     array of absorber space layer indices of the input
!                                  absorber amounts at the solar zenith angle.
!                                  units:      none.
!                                  type:       integer
!                                  dimension:  k x j x m
!                                  attributes: intent( in )
!
!       tau_predictor:             predictor profiles for the layer->toa transmittance.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  i x k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       flux_tau_predictor:        predictor profiles for the thermal flux transmittance.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  i x k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       solar_tau_predictor:       predictor profiles for the solar transmittance.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  i x k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       tau:                       layer->toa transmittance for the satellite
!                                  view angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       flux_tau:                  layer->sfc transmittance for the default
!                                  diffusivity angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       solar_tau:                 layer->sfc transmittance for the solar
!                                  zenith angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_radiance:            layer planck radiances at every layer for
!                                  each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       downwelling_radiance:      toa->sfc radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       upwelling_radiance:        toa radiances for each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       level_p_tl:                profile set layer interface pressure tangent-linear
!                                  array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_p_tl:                profile set layer average pressure tangent-linear
!                                  array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_t_tl:                profile set layer average temperature tangent-linear
!                                  array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_w_tl:      .         profile set layer average water vapor mixing ratio
!                                  tangent-linear array.
!                                  units:      g/kg
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       layer_o_tl:                profile set layer average ozone mixing ratio
!                                  tangent-linear array.
!                                  units:      ppmv
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in )
!
!       surface_temperature_tl:    profile set surface temperature tangent-linear array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  m; m > or = 1 (i.e. scalar)
!                                  attributes: intent( in )
!
!       surface_emissivity_tl:     profile set surface emissivity tangent-linear array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       surface_reflectivity_tl:   profile set surface reflectivity tangent-linear array
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
! optional input arguments:
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
!       tau_tl:                    layer->toa tangent-linear transmittance for the satellite
!                                  view angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( out )
!
!       flux_tau_tl:               layer->sfc tangent-linear transmittance for the default
!                                  diffusivity angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( out )
!
!       solar_tau_tl:              layer->sfc tangent-linear transmittance for the solar
!                                  zenith angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( out )
!
!       layer_radiance_tl:         layer planck tangent-linear radiances at every layer for
!                                  each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( out )
!
!       downwelling_radiance_tl:   toa->sfc tangent-linear radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
!       upwelling_radiance_tl:     toa tangent-linear radiances for each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
!       brightness_temperature_tl: tangent-linear temperatures corresponding to the
!                                  toa tangent-linear radiances.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
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
!      compute_absorber_amount_tl: subroutine to integrate the tangent-linear
!                                  absorber profiles
!                                  source: absorber_profile module
!
!      compute_predictors_tl:      subroutine to compute the tangent-linear 
!                                  transmittance predictor profiles.
!                                  source: predictor module
!
!      compute_transmittance_tl:   subroutine to compute the tangent-linear
!                                  transmittance profiles.
!                                  source: transmittance module
!
!      compute_radiance_tl:        subroutine to compute the toa tangent-linear 
!                                  radiances and brightness temperatures.
!                                  source: radiance module
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
!       note the restrictions on the input array dimensions:
!         k == n_layers > 1
!         l == n_channels > 1
!         m == n_profiles > or = 1
!
! procedure:
!       see individual module function documentation.
!s-
!--------------------------------------------------------------------------------

  function tangent_linear_rtm_rank2( &

             ! -- forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,                &  ! input, k x m

             surface_temperature,                                        &  ! input, m
             surface_emissivity,                                         &  ! input, l*m
             surface_reflectivity,                                       &  ! input, l*m

             absorber,                                                   &  ! input, 0:k x j x m

             tau_layer_index,                                            &  ! input, k x j x m
             flux_tau_layer_index,                                       &  ! input, k x j x m
             solar_tau_layer_index,                                      &  ! input, k x j x m

             tau_predictor,                                              &  ! input, imax x k x m
             flux_tau_predictor,                                         &  ! input, imax x k x m
             solar_tau_predictor,                                        &  ! input, imax x k x m

             tau,                                                        &  ! input, k x l*m
             flux_tau,                                                   &  ! input, k x l*m
             solar_tau,                                                  &  ! input, k x l*m

             layer_radiance,                                             &  ! input, k x l*m
             downwelling_radiance,                                       &  ! input, l*m
             upwelling_radiance,                                         &  ! input, l*m

             ! -- tangent-linear inputs
             level_p_tl, layer_p_tl, layer_t_tl, layer_w_tl, layer_o_tl, &  ! input, k x m

             surface_temperature_tl,                                     &  ! input, m
             surface_emissivity_tl,                                      &  ! input, l*m
             surface_reflectivity_tl,                                    &  ! input, l*m

             ! -- other inputs
             secant_view_angle,                                          &  ! input, m
             secant_solar_angle,                                         &  ! input, m
             n_channels_per_profile,                                     &  ! input, m
             channel_index,                                              &  ! input, l*m

             ! -- tangent-linear outputs
             tau_tl,                                                     &  ! output, k x l*m
             flux_tau_tl,                                                &  ! output, k x l*m
             solar_tau_tl,                                               &  ! output, k x l*m

             layer_radiance_tl,                                          &  ! output, k x l*m
             downwelling_radiance_tl,                                    &  ! output, l*m
             upwelling_radiance_tl,                                      &  ! output, l*m

             brightness_temperature_tl,                                  &  ! output, l*m

             ! -- optional inputs
             n_input_profiles,                                           &
             message_log )                                               &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward inputs
    real( fp_kind ), dimension( :, : ),     intent( in )  :: level_p                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_p                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_t                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_w                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_o                    ! k x m

    real( fp_kind ), dimension( : ),        intent( in )  :: surface_temperature        ! m
    real( fp_kind ), dimension( : ),        intent( in )  :: surface_emissivity         ! l*m
    real( fp_kind ), dimension( : ),        intent( in )  :: surface_reflectivity       ! l*m

    real( fp_kind ), dimension( 0:, :, : ), intent( in )  :: absorber                   ! 0:k x j x m

    integer,         dimension(  :, :, : ), intent( in )  :: tau_layer_index            ! k x j x m
    integer,         dimension(  :, :, : ), intent( in )  :: flux_tau_layer_index       ! k x j x m
    integer,         dimension(  :, :, : ), intent( in )  :: solar_tau_layer_index      ! k x j x m

    real( fp_kind ), dimension( :, :, : ),  intent( in )  :: tau_predictor              ! imax x k x m
    real( fp_kind ), dimension( :, :, : ),  intent( in )  :: flux_tau_predictor         ! imax x k x m
    real( fp_kind ), dimension( :, :, : ),  intent( in )  :: solar_tau_predictor        ! imax x k x m

    real( fp_kind ), dimension( :, : ),     intent( in )  :: tau                        ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: flux_tau                   ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: solar_tau                  ! k x l*m

    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_radiance             ! k x l*m
    real( fp_kind ), dimension( : ),        intent( in )  :: downwelling_radiance       ! l*m
    real( fp_kind ), dimension( : ),        intent( in )  :: upwelling_radiance         ! l*m

    ! -- tangent-linear inputs
    real( fp_kind ), dimension( :, : ),     intent( in )  :: level_p_tl                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_p_tl                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_t_tl                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_w_tl                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )  :: layer_o_tl                 ! k x m

    real( fp_kind ), dimension( : ),        intent( in )  :: surface_temperature_tl     ! m
    real( fp_kind ), dimension( : ),        intent( in )  :: surface_emissivity_tl      ! l*m
    real( fp_kind ), dimension( : ),        intent( in )  :: surface_reflectivity_tl    ! l*m

    ! -- other inputs
    real( fp_kind ), dimension( : ),        intent( in )  :: secant_view_angle          ! m
    real( fp_kind ), dimension( : ),        intent( in )  :: secant_solar_angle         ! m
    integer,         dimension( : ),        intent( in )  :: n_channels_per_profile     ! m
    integer,         dimension( : ),        intent( in )  :: channel_index              ! l*m

    ! -- tangent-linear outputs
    real( fp_kind ), dimension( :, : ),     intent( out ) :: tau_tl                     ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( out ) :: flux_tau_tl                ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( out ) :: solar_tau_tl               ! k x l*m

    real( fp_kind ), dimension( :, : ),     intent( out ) :: layer_radiance_tl          ! k x l*m
    real( fp_kind ), dimension( : ),        intent( out ) :: downwelling_radiance_tl    ! l*m
    real( fp_kind ), dimension( : ),        intent( out ) :: upwelling_radiance_tl      ! l*m

    real( fp_kind ), dimension( : ),        intent( out ) :: brightness_temperature_tl  ! l*m

    ! -- optional input
    integer,        optional,               intent( in )  :: n_input_profiles          ! scalar
    character( * ), optional,               intent( in )  :: message_log
    
    
    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'tangent_linear_rtm_rank2'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    character( 100 ) :: message
    character( 5 )   :: value_in, value_allowed

    integer :: m, n_profiles          ! profile loop/index variables
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

      error_status = tangent_linear_rtm_rank1( &
                       ! -- forward inputs
                       level_p( :, m ), layer_p( :, m ), layer_t( :, m ), layer_w( :, m ), layer_o( :, m ), &  ! input, k

                       surface_temperature( m ),           &  ! input,  scalar
                       surface_emissivity( l1:l2 ),        &  ! input,  l
                       surface_reflectivity( l1:l2 ),      &  ! input,  l

                       absorber( 0:, :, m ),               &  ! input, 0:k x j

                       tau_layer_index( :, :, m ),         &  ! input, k x j
                       flux_tau_layer_index( :, :, m ),    &  ! input, k x j
                       solar_tau_layer_index( :, :, m ),   &  ! input, k x j

                       tau_predictor( :, :, m ),           &  ! input, imax x k
                       flux_tau_predictor( :, :, m ),      &  ! input, imax x k
                       solar_tau_predictor( :, :, m ),     &  ! input, imax x k

                       tau( :, l1:l2 ),                    &  ! input, k x l
                       flux_tau( :, l1:l2 ),               &  ! input, k x l
                       solar_tau( :, l1:l2 ),              &  ! input, k x l

                       layer_radiance( :, l1:l2 ),         &  ! input, k x l
                       downwelling_radiance( l1:l2 ),      &  ! input, l
                       upwelling_radiance( l1:l2 ),        &  ! input, l

                       ! -- tangent-linear inputs
                       level_p_tl( :, m ), layer_p_tl( :, m ), layer_t_tl( :, m ), layer_w_tl( :, m ), layer_o_tl( :, m ), &  ! input, k

                       surface_temperature_tl,             &  ! input, scalar
                       surface_emissivity_tl( l1:l2 ),     &  ! input, l
                       surface_reflectivity_tl( l1:l2 ),   &  ! input, l

                       ! -- other inputs
                       secant_view_angle( m ),             &  ! input,  scalar
                       secant_solar_angle( m ),            &  ! input,  scalar
                       n_channels_per_profile( m ),        &  ! input,  scalar
                       channel_index( l1:l2 ),             &  ! input,  l

                       ! -- tangent-linear outputs
                       tau_tl( :, l1:l2 ),                 &  ! output, k x l*m
                       flux_tau_tl( :, l1:l2 ),            &  ! output, k x l*m
                       solar_tau_tl( :, l1:l2 ),           &  ! output, k x l*m

                       layer_radiance_tl( :, l1:l2 ),      &  ! output, k x l*m
                       downwelling_radiance_tl( l1:l2 ),   &  ! output, l*m
                       upwelling_radiance_tl( l1:l2 ),     &  ! output, l*m

                       brightness_temperature_tl( l1:l2 ), &  ! output, l*m

                       ! -- optional inputs
                       message_log = message_log )


      ! -------------------------------
      ! check for successful completion
      ! -------------------------------

      if ( error_status /= success ) then

        call display_message( routine_name, &
                              'error occured in tangent_linear_rtm_rank1', &
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

  end function tangent_linear_rtm_rank2



  function tangent_linear_rtm_rank1( &

             ! forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,                &  ! input, k

             surface_temperature,                                        &  ! input, scalar
             surface_emissivity,                                         &  ! input, l
             surface_reflectivity,                                       &  ! input, l

             absorber,                                                   &  ! input, 0:k x j

             tau_layer_index,                                            &  ! input, k x j
             flux_tau_layer_index,                                       &  ! input, k x j
             solar_tau_layer_index,                                      &  ! input, k x j

             tau_predictor,                                              &  ! input, imax x k
             flux_tau_predictor,                                         &  ! input, imax x k
             solar_tau_predictor,                                        &  ! input, imax x k

             tau,                                                        &  ! input, k x l
             flux_tau,                                                   &  ! input, k x l
             solar_tau,                                                  &  ! input, k x l

             layer_radiance,                                             &  ! input, k x l
             downwelling_radiance,                                       &  ! input, l
             upwelling_radiance,                                         &  ! input, l

             ! -- tangent-linear inputs
             level_p_tl, layer_p_tl, layer_t_tl, layer_w_tl, layer_o_tl, &  ! input, k

             surface_temperature_tl,                                     &  ! input, scalar
             surface_emissivity_tl,                                      &  ! input, l
             surface_reflectivity_tl,                                    &  ! input, l

             ! -- other inputs
             secant_view_angle,                                          &  ! input, scalar
             secant_solar_angle,                                         &  ! input, scalar
             n_channels,                                                 &  ! input, scalar
             channel_index,                                              &  ! input, l

             ! -- tangent-linear outputs
             tau_tl,                                                     &  ! output, k x l
             flux_tau_tl,                                                &  ! output, k x l
             solar_tau_tl,                                               &  ! output, k x l

             layer_radiance_tl,                                          &  ! output, k x l
             downwelling_radiance_tl,                                    &  ! output, l
             upwelling_radiance_tl,                                      &  ! output, l

             brightness_temperature_tl,                                  &  ! output, l

             ! -- optional inputs
             n_input_profiles,                                           &
             message_log )                                               &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward inputs
    real( fp_kind ), dimension( : ),     intent( in )  :: level_p                    ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_p                    ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_t                    ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_w                    ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_o                    ! k

    real( fp_kind ),                     intent( in )  :: surface_temperature        ! scalar
    real( fp_kind ), dimension( : ),     intent( in )  :: surface_emissivity         ! l
    real( fp_kind ), dimension( : ),     intent( in )  :: surface_reflectivity       ! l

    real( fp_kind ), dimension( 0:, : ), intent( in )  :: absorber                   ! 0:k x j

    integer,         dimension(  :, : ), intent( in )  :: tau_layer_index            ! k x j
    integer,         dimension(  :, : ), intent( in )  :: flux_tau_layer_index       ! k x j
    integer,         dimension(  :, : ), intent( in )  :: solar_tau_layer_index      ! k x j

    real( fp_kind ), dimension( :, : ),  intent( in )  :: tau_predictor              ! imax x k
    real( fp_kind ), dimension( :, : ),  intent( in )  :: flux_tau_predictor         ! imax x k
    real( fp_kind ), dimension( :, : ),  intent( in )  :: solar_tau_predictor        ! imax x k

    real( fp_kind ), dimension( :, : ),  intent( in )  :: tau                        ! k x l
    real( fp_kind ), dimension( :, : ),  intent( in )  :: flux_tau                   ! k x l
    real( fp_kind ), dimension( :, : ),  intent( in )  :: solar_tau                  ! k x l

    real( fp_kind ), dimension( :, : ),  intent( in )  :: layer_radiance             ! k x l
    real( fp_kind ), dimension( : ),     intent( in )  :: downwelling_radiance       ! l
    real( fp_kind ), dimension( : ),     intent( in )  :: upwelling_radiance         ! l

    ! -- tangent-linear inputs
    real( fp_kind ), dimension( : ),     intent( in )  :: level_p_tl                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_p_tl                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_t_tl                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_w_tl                 ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: layer_o_tl                 ! k

    real( fp_kind ),                     intent( in )  :: surface_temperature_tl     ! scalar
    real( fp_kind ), dimension( : ),     intent( in )  :: surface_emissivity_tl      ! l
    real( fp_kind ), dimension( : ),     intent( in )  :: surface_reflectivity_tl    ! l

    ! -- other inputs
    real( fp_kind ),                     intent( in )  :: secant_view_angle          ! scalar
    real( fp_kind ),                     intent( in )  :: secant_solar_angle         ! scalar
    integer,                             intent( in )  :: n_channels                 ! scalar
    integer,         dimension( : ),     intent( in )  :: channel_index              ! l

    ! -- tangent-linear outputs
    real( fp_kind ), dimension( :, : ),  intent( out ) :: tau_tl                     ! k x l
    real( fp_kind ), dimension( :, : ),  intent( out ) :: flux_tau_tl                ! k x l
    real( fp_kind ), dimension( :, : ),  intent( out ) :: solar_tau_tl               ! k x l

    real( fp_kind ), dimension( :, : ),  intent( out ) :: layer_radiance_tl          ! k x l
    real( fp_kind ), dimension( : ),     intent( out ) :: downwelling_radiance_tl    ! l
    real( fp_kind ), dimension( : ),     intent( out ) :: upwelling_radiance_tl      ! l

    real( fp_kind ), dimension( : ),     intent( out ) :: brightness_temperature_tl  ! l

    ! -- optional input. note that n_input_profiles is not used in this
    ! -- function. it is specified so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    integer,        optional,            intent( in )  :: n_input_profiles           ! scalar
    character( * ), optional,            intent( in )  :: message_log
    
    
    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'tangent_linear_rtm_rank1'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    character( 100 ) :: message
    character( 5 )   :: value_in, value_allowed

    integer :: k, n_layers  ! layer loop/index variables
    integer :: l            ! channel loop/index variable

    integer :: valid_solar

    ! -- maximum channels pseudo parameter
    integer :: max_n_channels
    logical :: is_set

    ! -- arrays for integrated absorber amounts, 0:k x j
    real( fp_kind ), dimension( 0:size( absorber, dim = 1 )-1, &
                                  size( absorber, dim = 2 )    ) :: tau_absorber,      &
                                                                    flux_tau_absorber, &
                                                                    solar_tau_absorber

    real( fp_kind ), dimension( 0:size( absorber, dim = 1 )-1, &
                                  size( absorber, dim = 2 )    ) :: absorber_tl, &
                                                                    tau_absorber_tl,      &
                                                                    flux_tau_absorber_tl, &
                                                                    solar_tau_absorber_tl

    ! -- arrays for tangent-linear predictors, imax x k
    real( fp_kind ), dimension( size( tau_predictor, dim = 1 ), &
                                size( tau_predictor, dim = 2 )  ) :: tau_predictor_tl,      &
                                                                     flux_tau_predictor_tl, &
                                                                     solar_tau_predictor_tl

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
    !#             -- determine array dimensions and check input --             #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! check the number of channels - if zero
    ! then simply return.
    ! --------------------------------------

    if ( n_channels == 0 ) return


    ! ------------------
    ! get the dimensions
    ! ------------------

    n_layers   = size( layer_p )


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
    ! profile data for negative values
    ! -----------------------------------

    ! -- profile data
    if ( any( level_p < zero ) .or. &
         any( layer_p < zero ) .or. &
         any( layer_t < zero ) .or. &
         any( layer_w < zero ) .or. &
         any( layer_o < zero )      ) then
      error_status = failure
      call display_message( routine_name, &
                            'negative values found in input profile data.', &
                            error_status, &
                            message_log = message_log )
      return
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
    !#   -- calculate the profile generic tangent-linear absorber amounts --    #
    !#--------------------------------------------------------------------------#

    call compute_absorber_amount_tl( &
                                     ! -- forward input
                                     level_p( : ),        &  ! input,  k
                                     layer_w( : ),        &  ! input,  k
                                     layer_o( : ),        &  ! input,  k

                                     ! -- tangent-linear input
                                     level_p_tl( : ),     &  ! input,  k
                                     layer_w_tl( : ),     &  ! input,  k
                                     layer_o_tl( : ),     &  ! input,  k

                                     ! -- tangent-linear output
                                     absorber_tl( 0:, : ) )  ! output, 0:k x j



    !#--------------------------------------------------------------------------#
    !#      -- calculate the predictors for the upwelling transmittance --      #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------------------
    ! modify absorber quantities by the angle secant
    ! could put a loop here but here's hoping the compiler
    ! recognises this as a group of loops over layer.
    ! ----------------------------------------------------

    tau_absorber(    0:, : ) = secant_view_angle * absorber(    0:, : )
    tau_absorber_tl( 0:, : ) = secant_view_angle * absorber_tl( 0:, : )

    ! -- subtract the top pressure from the dry absorber profile
    tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - toa_pressure


    ! ---------------------------------------
    ! calculate the tangent-linear predictors
    ! for the satellite view angle
    ! ---------------------------------------

    call compute_predictors_tl( &
                                ! -- forward input
                                layer_p( : ),              &  ! input,  k
                                layer_t( : ),              &  ! input,  k
                                layer_w( : ),              &  ! input,  k
                                tau_absorber( 0:, : ),     &  ! input,  0:k x j

                                ! -- tangent-linear input
                                layer_p_tl( : ),           &  ! input,  k
                                layer_t_tl( : ),           &  ! input,  k
                                layer_w_tl( : ),           &  ! input,  k
                                tau_absorber_tl( 0:, : ),  &  ! input,  0:k x j

                                ! -- tangent-linear output
                                tau_predictor_tl( :, : )   )  ! output, i x k
 


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

      flux_tau_absorber(    0:, : ) = secant_diffusivity_angle * absorber(    0:, : )
      flux_tau_absorber_tl( 0:, : ) = secant_diffusivity_angle * absorber_tl( 0:, : )
      
      ! -- subtract the top pressure from the dry absorber profile
      flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - toa_pressure


      ! ---------------------------------------
      ! calculate the tangent-linear predictors
      ! for the diffusivity angle
      ! ---------------------------------------

      ! -- calculate the integrated predictors only
      call compute_predictors_tl( &
                                  ! -- forward input
                                  layer_p( : ),                  &  ! input,  k
                                  layer_t( : ),                  &  ! input,  k
                                  layer_w( : ),                  &  ! input,  k
                                  flux_tau_absorber( 0:, : ),    &  ! input,  0:k x j

                                  ! -- tangent-linear input
                                  layer_p_tl( : ),               &  ! input,  k
                                  layer_t_tl( : ),               &  ! input,  k
                                  layer_w_tl( : ),               &  ! input,  k
                                  flux_tau_absorber_tl( 0:, : ), &  ! input,  0:k x j

                                  ! -- tangent-linear output
                                  flux_tau_predictor_tl( :, : ), &  ! output, i x k

                                  no_standard = 1                )  ! optional input

      ! -- copy the angle independent (standard) predictors
      flux_tau_predictor_tl( 1:max_n_standard_predictors, : ) = &
           tau_predictor_tl( 1:max_n_standard_predictors, : ) 

    end if



    !#--------------------------------------------------------------------------#
    !#        -- calculate the predictors for the solar transmittance --        #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------------------
    ! have *any* solar sensitive channels been specified
    ! for the current profile (flagged as == 1)?
    !
    ! and
    !
    ! is the specified solar zenith angle valid (<85)?
    ! --------------------------------------------------

    if ( ( any( is_solar_channel( channel_index( 1:n_channels ) ) == 1 ) ) .and. &
         secant_solar_angle < max_secant_solar_angle ) then


      ! --------------------------------
      ! modify the nadir absorber amount
      ! --------------------------------

      solar_tau_absorber(    0:, : ) = secant_solar_angle * absorber(    0:, : )
      solar_tau_absorber_tl( 0:, : ) = secant_solar_angle * absorber_tl( 0:, : )

      ! -- subtract the top pressure from the dry absorber profile
      solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - toa_pressure


      ! ---------------------------------------
      ! calculate the tangent-linear predictors
      ! for the solar zenith angle
      ! ---------------------------------------

      ! -- calculate the integrated predictors only
      call compute_predictors_tl( &
                                  ! -- forward input
                                  layer_p( : ),                   &  ! input,  k
                                  layer_t( : ),                   &  ! input,  k
                                  layer_w( : ),                   &  ! input,  k
                                  solar_tau_absorber( 0:, : ),    &  ! input,  0:k x j

                                  ! -- tangent-linear input
                                  layer_p_tl( : ),                &  ! input,  k
                                  layer_t_tl( : ),                &  ! input,  k
                                  layer_w_tl( : ),                &  ! input,  k
                                  solar_tau_absorber_tl( 0:, : ), &  ! input,  0:k x j

                                  ! -- tangent-linear output
                                  solar_tau_predictor_tl( :, : ), &  ! output, i x k

                                  no_standard = 1                 )  ! optional input

      ! -- copy the angle independent predictors
      solar_tau_predictor_tl( 1:max_n_standard_predictors, : ) = &
            tau_predictor_tl( 1:max_n_standard_predictors, : ) 

    end if
     

    !#--------------------------------------------------------------------------#
    !#                            -- channel loop --                            #
    !#--------------------------------------------------------------------------#

    l_channel_loop: do l = 1, n_channels


      ! --------------------------------------------
      ! calculate the current channel tangent-linear
      ! transmittances for the satellite view angle
      ! --------------------------------------------

      call compute_transmittance_tl( &
                                     ! -- forward input
                                     tau_absorber( 0:, : ),    &   ! input, 0:k x j
                                     tau_predictor( :, : ),    &   ! input, i x k
                                     tau( :, l ),              &   ! input, k

                                     ! -- tangent-linear input
                                     tau_absorber_tl( 0:, : ), &   ! input, 0:k x j
                                     tau_predictor_tl( :, : ), &   ! input, i x k

                                     ! -- other input
                                     tau_layer_index( :, : ),  &   ! input, k x j
                                     channel_index( l ),       &   ! input, scalar
                                     up,                       &   ! input, scalar

                                     ! -- tangent-linear output
                                     tau_tl( :, l )            )   ! output, k



      ! --------------------------------------------------
      ! if the current channel is an infrared channel,
      ! then calculate the downwelling flux transmittance.
      !
      ! if the current channel is a microwave channel,
      ! then use the predictors for the upwelling
      ! transmittance calculations.
      !
      ! this situation is a toss-up. one could return the
      ! layer optical depth rather than the transmittance
      ! and use that to derive the microwave flux
      ! transmittance but then there would be two sets
      ! of quantities to be carried through - the layer
      ! to boundary transmittances and the individual
      ! layer optical depths. if the speed penalty of
      ! calculating the microwave flux transmittances in
      ! this manner is too high, and memory is available,
      ! then this is an option.
      ! --------------------------------------------------

      if ( is_microwave_channel( channel_index( l ) ) == 0 ) then

        ! -- ir channel
        call compute_transmittance_tl( &
                                       ! -- forward input
                                       flux_tau_absorber( 0:, : ),    &   ! input, 0:k x j
                                       flux_tau_predictor( :, : ),    &   ! input, i x k
                                       flux_tau( :, l ),              &   ! input, k

                                       ! -- tangent-linear input
                                       flux_tau_absorber_tl( 0:, : ), &   ! input, 0:k x j
                                       flux_tau_predictor_tl( :, : ), &   ! input, i x k

                                       ! -- other input
                                       flux_tau_layer_index( :, : ),  &   ! input, k x j
                                       channel_index( l ),            &   ! input, scalar
                                       down,                          &   ! input, scalar

                                       ! -- tangent-linear output
                                       flux_tau_tl( :, l )            )   ! output, k

      else

        ! -- uw channel
!        call compute_transmittance_tl( &
!                                       ! -- forward input
!                                       tau_absorber( 0:, : ),    &   ! input, 0:k x j
!                                       tau_predictor( :, : ),    &   ! input, i x k
!                                       flux_tau( :, l ),         &   ! input, k
!
!                                       ! -- tangent-linear input
!                                       tau_absorber_tl( 0:, : ), &   ! input, 0:k x j
!                                       tau_predictor_tl( :, : ), &   ! input, i x k
!
!                                       ! -- other input
!                                       tau_layer_index( :, : ),  &   ! input, k x j
!                                       channel_index( l ),       &   ! input, scalar
!                                       down,                     &   ! input, scalar
!
!                                       ! -- tangent-linear output
!                                       flux_tau_tl( :, l )       )   ! output, k

        ! -- if total transmittance /= 0, calculate 
        ! -- tangent-linear flux transmittance
        flux_tau_tl( :, l ) = zero
        if ( tau( n_layers, l ) > tolerance ) then
          flux_tau_tl( 1, l ) = tau_tl( n_layers, l )
          do k = 2, n_layers

            flux_tau_tl( k, l ) = ( flux_tau_tl( 1, l ) / &
            !                       -------------------   
                                      tau( k-1, l )     ) - &

                                  ( flux_tau( 1, l ) * tau_tl( k-1, l ) / &
            !                       -----------------------------------
                                              tau( k-1, l )**2          )
          end do
        end if

      end if



      ! ----------------------------------------------------
      ! if the current channel is a solar sensitive channel,
      !   and
      ! the solar angle is valid, then calculate the
      ! transmittance for direct solar.
      ! ----------------------------------------------------

      if ( is_solar_channel( channel_index( l ) ) == 1 .and. &
           secant_solar_angle < max_secant_solar_angle       ) then

        valid_solar = 1

        call compute_transmittance_tl( &
                                       ! -- forward input
                                       solar_tau_absorber( 0:, : ),    &   ! input, 0:k x j
                                       solar_tau_predictor( :, : ),    &   ! input, i x k
                                       solar_tau( :, l ),              &   ! input, k

                                       ! -- tangent-linear input
                                       solar_tau_absorber_tl( 0:, : ), &   ! input, 0:k x j
                                       solar_tau_predictor_tl( :, : ), &   ! input, i x k

                                       ! -- other input
                                       solar_tau_layer_index( :, : ),  &   ! input, k x j
                                       channel_index( l ),             &   ! input, scalar
                                       down,                           &   ! input, scalar

                                       ! -- tangent-linear output
                                       solar_tau_tl( :, l )            )   ! output, k

      else

        valid_solar = 0
        solar_tau_tl( :, l ) = zero

      end if


      ! ------------------------------------------------------
      ! calculate the profile/channel tangent-linear radiances
      ! ------------------------------------------------------

      call compute_radiance_tl( &
                                ! -- forward input
                                layer_t( : ),                  &  ! input, k

                                surface_temperature,           &  ! input, scalar
                                surface_emissivity( l ),       &  ! input, scalar
                                surface_reflectivity( l ),     &  ! input, scalar

                                tau(      :, l ),              &  ! input, k
                                flux_tau( :, l ),              &  ! input, k
                                solar_tau( n_layers, l ),      &  ! input, scalar

                                layer_radiance( :, l ),        &  ! input, k
                                downwelling_radiance( l ),     &  ! input, scalar
                                upwelling_radiance( l ),       &  ! input, scalar

                                ! -- tangent-linear input
                                layer_t_tl( : ),               &  ! input, k

                                surface_temperature_tl,        &  ! input, scalar
                                surface_emissivity_tl( l ),    &  ! input, scalar
                                surface_reflectivity_tl( l ),  &  ! input, scalar

                                tau_tl(      :, l ),           &  ! input, k
                                flux_tau_tl( :, l ),           &  ! input, k
                                solar_tau_tl( n_layers, l ),   &  ! input, scalar

                                ! -- other input
                                secant_solar_angle,            &  ! input, scalar
                                valid_solar,                   &  ! input, scalar
                                channel_index( l ),            &  ! input, scalar

                                ! -- tangent-linear output
                                layer_radiance_tl( :, l ),     &  ! output, k
                                downwelling_radiance_tl( l ),  &  ! output, scalar
                                upwelling_radiance_tl( l ),    &  ! output, scalar

                                brightness_temperature_tl( l ) )  ! output, scalar


    end do l_channel_loop



    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function tangent_linear_rtm_rank1

end module tangent_linear_model


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
! revision 1.9  2001/11/07 15:08:07  paulv
! - changed the logical if test for the number of input profiles in the
!   [<tangent_linear><adjoint><k_matrix>]_rtm_rank2() functions from
!     if ( n_input_profiles >= 1 .and. n_input_profiles <= n_profiles ) then
!   to
!     if ( n_input_profiles > 0 .and. n_input_profiles <= n_profiles ) then
!   the use of both the ">=" and "<=" relational operators with the .and. i
!   found confusing.
!
! revision 1.8  2001/11/07 14:57:46  paulv
! - added check for negative number of channels to tangent_linear_rtm_rank2() function.
! - added profile loop cycle statement to tangent_linear_rtm_rank2() function.
! - added check for negative number of channels to tangent_linear_rtm_rank1() function.
!
! revision 1.7  2001/10/01 20:26:20  paulv
! - overloaded the tangent_linear_rtm function to accept both a
!   single or group of profiles. contained private functions are now
!     o tangent_linear_rtm_rank1 for single profile input
!     o tangent_linear_rtm_rank2 for multiple profile input
! - put n_input_profiles optional argument back in tangent_linear_rtm
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
!   as the check is now done in the rank-1 tangent_linear_rtm function. this eliminates
!   the need to fully populate the input arrays with data when only certain
!   "chunks" may be processed. previously, the use of any() could generate
!   and error if the full surface_temperature array was not initialised.
! - added "name" to rcs keyword list.
!
! revision 1.6  2001/09/04 21:29:11  paulv
! - updated documentation.
!
! revision 1.5  2001/08/31 21:33:04  paulv
! - added check for negative profile and surface data in tangent_linear_rtm.
! - maximum solar angle secant is no longer calculated in tangent_linear_rtm but
!   is declared as a parameter in the parameters module.
!
! revision 1.4  2001/08/16 16:37:52  paulv
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
! revision 1.3  2001/08/01 17:00:16  paulv
! - altered function declaration to avoid more than the standard 39
!   continuation lines allowed in fortran 90.
! - updated input argument checking. now consistent with other model
!   components.
!
! revision 1.2  2001/07/12 18:46:06  paulv
! - commented out informational message output at start of function.
! - corrected bug in the calculation of the thermal flux transmittance for
!   the infrared channels. the call to compute_transmittance_tl used the
!   absorber bracketing index array tau_layer_index (for the upwelling
!   transmittance) rather than flux_tau_layer_index.
!
! revision 1.1  2001/05/29 16:36:02  paulv
! initial checkin
!
!
!
