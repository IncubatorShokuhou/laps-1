!------------------------------------------------------------------------------
!m+
! name:
!       adjoint_model
!
! purpose:
!       module containing the ncep rt adjoint model function
!
! category:
!       ncep rtm
!
! calling sequence:
!       use adjoint_model
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
!       adjoint_rtm:           public function that calculates the adjoint of the
!                              top-of-atmosphere (toa) radiances and brightness
!                              temperatures for user specified profiles and
!                              satellites/channels.
!
!       compute_rtm_ad:        public function that calculates the adjoint of the
!                              top-of-atmosphere (toa) radiances and brightness
!                              temperatures for an input atmospheric profile
!                              set and user specified satellites/channels.
!
!                              this function is simply a wrapper around both the
!                              forward model (fortward_rtm) and the adjoint model
!                              (adjoint_rtm) so that the user doesn't have to declare
!                              the absorber/predictor/etc arrays in the calling routine.
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
!       written by:     paul van delst, cimss@noaa/ncep 05-feb-2001
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

module adjoint_model


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
  use forward_model, only: forward_rtm


  ! -----------------------
  ! disable implicit typing
  ! -----------------------

  implicit none


  ! --------------------
  ! default visibilities
  ! --------------------

  private
  public :: compute_rtm_ad
  public :: adjoint_rtm


  ! --------------------
  ! function overloading
  ! --------------------

  interface compute_rtm_ad
    module procedure compute_rtm_ad_rank1
    module procedure compute_rtm_ad_rank2
  end interface ! compute_rtm_ad

  interface adjoint_rtm
    module procedure adjoint_rtm_rank1
    module procedure adjoint_rtm_rank2
  end interface ! adjoint_rtm


contains


!--------------------------------------------------------------------------------
!s+
! name:
!       compute_rtm_ad
!
! purpose:
!       public function that calculates the adjoint of the top-of-atmosphere (toa)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
!       this function is simply a wrapper around both the forward model and the
!       adjoint model so that the user doesn't have to declare the absorber/
!       predictor/etc arrays in the calling routine.
!
! category:
!       ncep rtm
!
! calling sequence:
!       result = compute_rtm_ad( &
!                               ! -- forward inputs
!                               level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k x m
!
!                               surface_temperature,                         &  ! input, m
!                               surface_emissivity,                          &  ! input, l*m
!                               surface_reflectivity,                        &  ! input, l*m
!
!                               ! -- adjoint inputs
!                               tau_ad,                                      &  ! in/output, k x l*m
!                               flux_tau_ad,                                 &  ! in/output, k x l*m
!                               solar_tau_ad,                                &  ! in/output, k x l*m
!
!                               upwelling_radiance_ad,                       &  ! in/output, l*m
!                               brightness_temperature_ad,                   &  ! in/output, l*m
!
!                               ! -- other inputs
!                               secant_view_angle,                           &  ! input, m
!                               secant_solar_angle,                          &  ! input, m
!                               n_channels_per_profile,                      &  ! input, m
!                               channel_index,                               &  ! input, l*m
!
!                               ! -- forward output
!                               tau,                                         &  ! input, k x l*m
!                               flux_tau,                                    &  ! input, k x l*m
!                               solar_tau,                                   &  ! input, k x l*m
!
!                               upwelling_radiance,                          &  ! input, l*m
!                               brightness_temperature,                      &  ! input, l*m
!
!                               ! -- adjoint outputs
!                               level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output, k x m
!
!                               surface_temperature_ad,                      &  ! in/output, m
!                               surface_emissivity_ad,                       &  ! in/output, l*m
!                               surface_reflectivity_ad,                     &  ! in/output, l*m
!
!                               ! optional inputs
!                               message_log = message_log )
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
!       tau_ad:                    layer->toa adjoint transmittance for the satellite
!                                  view angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       flux_tau_ad:               layer->sfc adjoint transmittance for the default
!                                  diffusivity angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       solar_tau_ad:              layer->sfc adjoint transmittance for the solar
!                                  zenith angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       upwelling_radiance_ad:     toa adjoint radiances for each channel/profile.
!                                  units:      (m^2.sr.cm^-1)/mw
!                                  ** this argument is set to zero on output **.
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       brightness_temperature_ad: adjoint temperatures corresponding to the
!                                  toa adjoint radiances.
!                                  ** this argument is set to zero on output **.
!                                  units:      k^-1
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
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
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
!       level_p_ad:                profile set layer interface pressure adjoint array.
!                                  units:      hpa^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_p_ad:                profile set layer average pressure adjoint array.
!                                  units:      hpa^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_t_ad:                profile set layer average temperature adjoint array.
!                                  units:      k^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_w_ad:      .         profile set layer average water vapor mixing ratio
!                                  adjoint array.
!                                  units:      kg/g
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_o_ad:                profile set layer average ozone mixing ratio
!                                  adjoint array.
!                                  units:      ppmv^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       surface_temperature_ad:    profile set surface temperature adjoint
!                                  array.
!                                  units:      k^-1
!                                  type:       real
!                                  dimension:  m; m > or = 1 (i.e. scalar)
!                                  attributes: intent( in out )
!
!       surface_emissivity_ad:     profile set surface emissivity k-matrix adjoint
!                                  array.
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       surface_reflectivity_ad:    profile set surface reflectivity k-matrix adjoint
!                                  array.
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
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
!                                  source: forward_model module
!
!      adjoint_rtm:                function that calculates the adjoint of the (toa)
!                                  radiances/temperatures.
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       all input adjoint arguments are set to zero on output.
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

  function compute_rtm_ad_rank2( &
             ! -- forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k x m

             surface_temperature,                         &  ! input, m
             surface_emissivity,                          &  ! input, l*m
             surface_reflectivity,                        &  ! input, l*m

             ! -- k-matrix inputs
             tau_ad,                                      &  ! in/output, k x l*m
             flux_tau_ad,                                 &  ! in/output, k x l*m
             solar_tau_ad,                                &  ! in/output, k x l*m

             upwelling_radiance_ad,                       &  ! in/output, l*m
             brightness_temperature_ad,                   &  ! in/output, l*m

             ! -- other inputs
             secant_view_angle,                           &  ! input, m
             secant_solar_angle,                          &  ! input, m
             n_channels_per_profile,                      &  ! input, m
             channel_index,                               &  ! input, l*m

             ! -- forward output
             tau,                                         &  ! input, k x l*m
             flux_tau,                                    &  ! input, k x l*m
             solar_tau,                                   &  ! input, k x l*m

             upwelling_radiance,                          &  ! input, l*m
             brightness_temperature,                      &  ! input, l*m

             ! -- k-matrix outputs
             level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output, k x m

             surface_temperature_ad,                      &  ! in/output, m
             surface_emissivity_ad,                       &  ! in/output, l*m
             surface_reflectivity_ad,                     &  ! in/output, l*m

             ! -- optional inputs
             n_input_profiles,                            &
             message_log )                                &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations  --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward inputs
    real( fp_kind ), dimension( :, : ),     intent( in )     :: level_p                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_p                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_t                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_w                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_o                    ! k x m

    real( fp_kind ), dimension( : ),        intent( in )     :: surface_temperature        ! m
    real( fp_kind ), dimension( : ),        intent( in )     :: surface_emissivity         ! l*m
    real( fp_kind ), dimension( : ),        intent( in )     :: surface_reflectivity       ! l*m

    ! -- adjoint inputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: tau_ad                     ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: flux_tau_ad                ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: solar_tau_ad               ! k x l*m

    real( fp_kind ), dimension( : ),        intent( in out ) :: upwelling_radiance_ad      ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: brightness_temperature_ad  ! l*m

    ! -- other inputs
    real( fp_kind ), dimension( : ),        intent( in )     :: secant_view_angle          ! m
    real( fp_kind ), dimension( : ),        intent( in )     :: secant_solar_angle         ! m
    integer,         dimension( : ),        intent( in )     :: n_channels_per_profile     ! m
    integer,         dimension( : ),        intent( in )     :: channel_index              ! l*m

    ! -- forward outputs
    real( fp_kind ), dimension( :, : ),     intent( out )    :: tau                        ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( out )    :: flux_tau                   ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( out )    :: solar_tau                  ! k x l*m

    real( fp_kind ), dimension( : ),        intent( out )    :: upwelling_radiance         ! l*m
    real( fp_kind ), dimension( : ),        intent( out )    :: brightness_temperature     ! l*m

    ! -- adjoint outputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: level_p_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_p_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_t_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_w_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_o_ad                 ! k x m

    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_temperature_ad     ! m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_emissivity_ad      ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_reflectivity_ad    ! l*m

    ! -- optional input
    integer,        optional,               intent( in )     :: n_input_profiles           ! scalar
    character( * ), optional,               intent( in )     :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_rtm_ad_rank2'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    integer :: error_status_fwd
    integer :: error_status_ad

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

    ! -- array for forward and adjoint layer planck radiance term, k x l*m
    real( fp_kind ), dimension( size( layer_p, dim = 1 ),  &
                                size( upwelling_radiance ) ) :: layer_radiance,  &
                                                                layer_radiance_ad
      

    ! -- array for forward and adjoint downwelling radiance (flux + solar), l*m
    real( fp_kind ), dimension( size( upwelling_radiance ) ) :: downwelling_radiance,  &
                                                                downwelling_radiance_ad 



    ! ----------
    ! intrinsics
    ! ----------

    intrinsic adjustl, &
              maxval,  &
              size,    &
              trim



    !#--------------------------------------------------------------------------#
    !#                   -- compute the forward radiances --                    #
    !#--------------------------------------------------------------------------#

    error_status_fwd = forward_rtm( &
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

    if ( error_status_fwd /= success ) then

      error_status = failure
      call display_message( routine_name, &
                            'error occured in forward_rtm', &
                            error_status, &
                            message_log = message_log )
      return

    end if



    !#--------------------------------------------------------------------------#
    !#                     -- compute the adjoint rtm --                        #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! initialise all local adjoint variables
    ! --------------------------------------

    layer_radiance_ad( :, : )    = zero
    downwelling_radiance_ad( : ) = zero


    ! ----------------------
    ! call the adjoint model
    ! ----------------------

    error_status_ad = adjoint_rtm_rank2( &
                        ! -- forward inputs
                        level_p, layer_p, layer_t, layer_w, layer_o, &  ! input,  k x m

                        surface_temperature,                         &  ! input, m
                        surface_emissivity,                          &  ! input, l*m
                        surface_reflectivity,                        &  ! input, l*m

                        absorber,                                    &  ! input, 0:k x j x m

                        tau_layer_index,                             &  ! input, k x j x m
                        flux_tau_layer_index,                        &  ! input, k x j x m
                        solar_tau_layer_index,                       &  ! input, k x j x m

                        tau_predictor,                               &  ! input, imax x k x m
                        flux_tau_predictor,                          &  ! input, imax x k x m
                        solar_tau_predictor,                         &  ! input, imax x k x m

                        tau,                                         &  ! input, k x l*m
                        flux_tau,                                    &  ! input, k x l*m
                        solar_tau,                                   &  ! input, k x l*m

                        layer_radiance,                              &  ! input, k x l*m
                        downwelling_radiance,                        &  ! input, l*m
                        upwelling_radiance,                          &  ! input, l*m

                        ! -- adjoint inputs
                        tau_ad,                                      &  ! in/output, k x l*m
                        flux_tau_ad,                                 &  ! in/output, k x l*m
                        solar_tau_ad,                                &  ! in/output, k x l*m

                        layer_radiance_ad,                           &  ! in/output, k x l*m
                        downwelling_radiance_ad,                     &  ! in/output, l*m
                        upwelling_radiance_ad,                       &  ! in/output, l*m

                        brightness_temperature_ad,                   &  ! in/output, l*m

                        ! -- other inputs
                        secant_view_angle,                           &  ! input, m
                        secant_solar_angle,                          &  ! input, m
                        n_channels_per_profile,                      &  ! input, m
                        channel_index,                               &  ! input, l*m

                        ! -- adjoint outputs
                        level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output,  k x m

                        surface_temperature_ad,                      &  ! in/output, m
                        surface_emissivity_ad,                       &  ! in/output, l*m
                        surface_reflectivity_ad,                     &  ! in/output, l*m

                       ! -- optional inputs
                       n_input_profiles = n_input_profiles,          &
                       message_log      = message_log                )


    ! -------------------------------
    ! check for successful completion
    ! -------------------------------

    if ( error_status_ad /= success ) then

      error_status = failure
      call display_message( routine_name, &
                            'error occured in adjoint_rtm', &
                            error_status, &
                            message_log = message_log )
      return

    end if


    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success


  end function compute_rtm_ad_rank2





  function compute_rtm_ad_rank1( &
             ! -- forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o, &  ! input, k

             surface_temperature,                         &  ! input, scalar
             surface_emissivity,                          &  ! input, l
             surface_reflectivity,                        &  ! input, l

             ! -- k-matrix inputs
             tau_ad,                                      &  ! in/output, k x l
             flux_tau_ad,                                 &  ! in/output, k x l
             solar_tau_ad,                                &  ! in/output, k x l

             upwelling_radiance_ad,                       &  ! in/output, l
             brightness_temperature_ad,                   &  ! in/output, l

             ! -- other inputs
             secant_view_angle,                           &  ! input, scalar
             secant_solar_angle,                          &  ! input, scalar
             n_channels_per_profile,                      &  ! input, scalar
             channel_index,                               &  ! input, l

             ! -- forward output
             tau,                                         &  ! input, k x l
             flux_tau,                                    &  ! input, k x l
             solar_tau,                                   &  ! input, k x l

             upwelling_radiance,                          &  ! input, l
             brightness_temperature,                      &  ! input, l

             ! -- k-matrix outputs
             level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output, k

             surface_temperature_ad,                      &  ! in/output, scalar
             surface_emissivity_ad,                       &  ! in/output, l
             surface_reflectivity_ad,                     &  ! in/output, l

             ! -- optional inputs
             n_input_profiles,                            &
             message_log )                                &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations  --                         #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward inputs
    real( fp_kind ), dimension( : ),    intent( in )     :: level_p                    ! k
    real( fp_kind ), dimension( : ),    intent( in )     :: layer_p                    ! k
    real( fp_kind ), dimension( : ),    intent( in )     :: layer_t                    ! k
    real( fp_kind ), dimension( : ),    intent( in )     :: layer_w                    ! k
    real( fp_kind ), dimension( : ),    intent( in )     :: layer_o                    ! k

    real( fp_kind ),                    intent( in )     :: surface_temperature        ! scalar
    real( fp_kind ), dimension( : ),    intent( in )     :: surface_emissivity         ! l
    real( fp_kind ), dimension( : ),    intent( in )     :: surface_reflectivity       ! l

    ! -- adjoint inputs
    real( fp_kind ), dimension( :, : ), intent( in out ) :: tau_ad                     ! k x l
    real( fp_kind ), dimension( :, : ), intent( in out ) :: flux_tau_ad                ! k x l
    real( fp_kind ), dimension( :, : ), intent( in out ) :: solar_tau_ad               ! k x l

    real( fp_kind ), dimension( : ),    intent( in out ) :: upwelling_radiance_ad      ! l
    real( fp_kind ), dimension( : ),    intent( in out ) :: brightness_temperature_ad  ! l

    ! -- other inputs
    real( fp_kind ),                    intent( in )     :: secant_view_angle          ! scalar
    real( fp_kind ),                    intent( in )     :: secant_solar_angle         ! scalar
    integer,                            intent( in )     :: n_channels_per_profile     ! scalar
    integer,         dimension( : ),    intent( in )     :: channel_index              ! l

    ! -- forward outputs
    real( fp_kind ), dimension( :, : ), intent( out )    :: tau                        ! k x l
    real( fp_kind ), dimension( :, : ), intent( out )    :: flux_tau                   ! k x l
    real( fp_kind ), dimension( :, : ), intent( out )    :: solar_tau                  ! k x l

    real( fp_kind ), dimension( : ),    intent( out )    :: upwelling_radiance         ! l
    real( fp_kind ), dimension( : ),    intent( out )    :: brightness_temperature     ! l

    ! -- adjoint outputs
    real( fp_kind ), dimension( : ),    intent( in out ) :: level_p_ad                 ! k
    real( fp_kind ), dimension( : ),    intent( in out ) :: layer_p_ad                 ! k
    real( fp_kind ), dimension( : ),    intent( in out ) :: layer_t_ad                 ! k
    real( fp_kind ), dimension( : ),    intent( in out ) :: layer_w_ad                 ! k
    real( fp_kind ), dimension( : ),    intent( in out ) :: layer_o_ad                 ! k

    real( fp_kind ),                    intent( in out ) :: surface_temperature_ad     ! scalar
    real( fp_kind ), dimension( : ),    intent( in out ) :: surface_emissivity_ad      ! l
    real( fp_kind ), dimension( : ),    intent( in out ) :: surface_reflectivity_ad    ! l

    ! -- optional input. note that n_input_profiles is not used in this
    ! -- function. it is included here so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    integer,        optional,           intent( in )     :: n_input_profiles           ! scalar
    character( * ), optional,           intent( in )     :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_rtm_ad_rank1'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    integer :: error_status_fwd
    integer :: error_status_ad

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

    ! -- array for forward and adjoint layer planck radiance term, k x l
    real( fp_kind ), dimension( size( layer_p, dim = 1 ),  &
                                size( upwelling_radiance ) ) :: layer_radiance,  &
                                                                layer_radiance_ad
      

    ! -- array for forward and adjoint downwelling radiance (flux + solar), l
    real( fp_kind ), dimension( size( upwelling_radiance ) ) :: downwelling_radiance,  &
                                                                downwelling_radiance_ad 



    ! ----------
    ! intrinsics
    ! ----------

    intrinsic adjustl, &
              maxval,  &
              size,    &
              trim



    !#--------------------------------------------------------------------------#
    !#                   -- compute the forward radiances --                    #
    !#--------------------------------------------------------------------------#

    error_status_fwd = forward_rtm( &
                                    ! -- forward inputs
                                    level_p, layer_p, layer_t, layer_w, layer_o,  &  ! input,  k

                                    surface_temperature,                          &  ! input,  scalar
                                    surface_emissivity,                           &  ! input,  l
                                    surface_reflectivity,                         &  ! input,  l

                                    ! -- other inputs
                                    secant_view_angle,                            &  ! input,  scalar
                                    secant_solar_angle,                           &  ! input,  scalar
                                    n_channels_per_profile,                       &  ! input,  scalar
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

    if ( error_status_fwd /= success ) then

      error_status = failure
      call display_message( routine_name, &
                            'error occured in forward_rtm', &
                            error_status, &
                            message_log = message_log )
      return

    end if



    !#--------------------------------------------------------------------------#
    !#                     -- compute the adjoint rtm --                        #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! initialise all local adjoint variables
    ! --------------------------------------

    layer_radiance_ad( :, : )    = zero
    downwelling_radiance_ad( : ) = zero


    ! ----------------------
    ! call the adjoint model
    ! ----------------------

    error_status_ad = adjoint_rtm_rank1( &
                        ! -- forward inputs
                        level_p, layer_p, layer_t, layer_w, layer_o, &  ! input,  k

                        surface_temperature,                         &  ! input, scalar
                        surface_emissivity,                          &  ! input, l
                        surface_reflectivity,                        &  ! input, l

                        absorber,                                    &  ! input, 0:k x j

                        tau_layer_index,                             &  ! input, k x j
                        flux_tau_layer_index,                        &  ! input, k x j
                        solar_tau_layer_index,                       &  ! input, k x j

                        tau_predictor,                               &  ! input, imax x k
                        flux_tau_predictor,                          &  ! input, imax x k
                        solar_tau_predictor,                         &  ! input, imax x k

                        tau,                                         &  ! input, k x l
                        flux_tau,                                    &  ! input, k x l
                        solar_tau,                                   &  ! input, k x l

                        layer_radiance,                              &  ! input, k x l
                        downwelling_radiance,                        &  ! input, l
                        upwelling_radiance,                          &  ! input, l

                        ! -- adjoint inputs
                        tau_ad,                                      &  ! in/output, k x l
                        flux_tau_ad,                                 &  ! in/output, k x l
                        solar_tau_ad,                                &  ! in/output, k x l

                        layer_radiance_ad,                           &  ! in/output, k x l
                        downwelling_radiance_ad,                     &  ! in/output, l
                        upwelling_radiance_ad,                       &  ! in/output, l

                        brightness_temperature_ad,                   &  ! in/output, l

                        ! -- other inputs
                        secant_view_angle,                           &  ! input, scalar
                        secant_solar_angle,                          &  ! input, scalar
                        n_channels_per_profile,                      &  ! input, scalar
                        channel_index,                               &  ! input, l

                        ! -- adjoint outputs
                        level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output,  k

                        surface_temperature_ad,                      &  ! in/output, scalar
                        surface_emissivity_ad,                       &  ! in/output, l
                        surface_reflectivity_ad,                     &  ! in/output, l

                        message_log = message_log )


    ! -------------------------------
    ! check for successful completion
    ! -------------------------------

    if ( error_status_ad /= success ) then

      error_status = failure
      call display_message( routine_name, &
                            'error occured in adjoint_rtm_rank1', &
                            error_status, &
                            message_log = message_log )
      return

    end if


    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success


  end function compute_rtm_ad_rank1






!--------------------------------------------------------------------------------
!s+
! name:
!       adjoint_rtm
!
! purpose:
!       public function that calculates the adjoint of the top-of-atmosphere (toa)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
! category:
!       ncep rtm
!
! calling sequence:
!       result = adjoint_rtm( &
!                             ! -- forward inputs
!                             level_p, layer_p, layer_t, layer_w, layer_o,                &  ! input, k x m
!
!                             surface_temperature,                                        &  ! input, m
!                             surface_emissivity,                                         &  ! input, l*m
!                             surface_reflectivity,                                       &  ! input, l*m
!
!                             absorber,                                                   &  ! input, 0:k x j x m
!
!                             tau_layer_index,                                            &  ! input, k x j x m
!                             flux_tau_layer_index,                                       &  ! input, k x j x m
!                             solar_tau_layer_index,                                      &  ! input, k x j x m
!
!                             tau_predictor,                                              &  ! input, imax x k x m
!                             flux_tau_predictor,                                         &  ! input, imax x k x m
!                             solar_tau_predictor,                                        &  ! input, imax x k x m
!
!                             tau,                                                        &  ! input, k x l*m
!                             flux_tau,                                                   &  ! input, k x l*m
!                             solar_tau,                                                  &  ! input, k x l*m
!
!                             layer_radiance,                                             &  ! input, k x l*m
!                             downwelling_radiance,                                       &  ! input, l*m
!                             upwelling_radiance,                                         &  ! input, l*m
!
!                             ! -- adjoint inputs
!                             tau_ad,                                                     &  ! in/output, k x l*m
!                             flux_tau_ad,                                                &  ! in/output, k x l*m
!                             solar_tau_ad,                                               &  ! in/output, k x l*m
!
!                             layer_radiance_ad,                                          &  ! in/output, k x l*m
!                             downwelling_radiance_ad,                                    &  ! in/output, l*m
!                             upwelling_radiance_ad,                                      &  ! in/output, l*m
!
!                             brightness_temperature_ad,                                  &  ! in/output, l*m
!
!                             ! -- other inputs
!                             secant_view_angle,                                          &  ! input, m
!                             secant_solar_angle,                                         &  ! input, m
!                             n_channels_per_profile,                                     &  ! input, m
!                             channel_index,                                              &  ! input, l*m
!
!                             ! -- adjoint outputs
!                             level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output, k x m
!
!                             surface_temperature_ad,                                     &  ! in/output, m
!                             surface_emissivity_ad,                                      &  ! in/output, l*m
!                             surface_reflectivity_ad,                                    &  ! in/output, l*m
!
!                             ! optional inputs
!                             message_log = message_log )
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
!       tau_ad:                    layer->toa adjoint transmittance for the satellite
!                                  view angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       flux_tau_ad:               layer->sfc adjoint transmittance for the default
!                                  diffusivity angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       solar_tau_ad:              layer->sfc adjoint transmittance for the solar
!                                  zenith angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_radiance_ad:         layer planck adjoint radiances at every layer for
!                                  each channel/profile.
!                                  ** this argument is set to zero on output **.
!                                  units:      (m^2.sr.cm^-1)/mw
!                                  type:       real
!                                  dimension:  k x l*m; k > 1, l > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       downwelling_radiance_ad:   toa->sfc adjoint radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  ** this argument is set to zero on output **.
!                                  units:      (m^2.sr.cm^-1)/mw
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       upwelling_radiance_ad:     toa adjoint radiances for each channel/profile.
!                                  units:      (m^2.sr.cm^-1)/mw
!                                  ** this argument is set to zero on output **.
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       brightness_temperature_ad: adjoint temperatures corresponding to the
!                                  toa adjoint radiances.
!                                  ** this argument is set to zero on output **.
!                                  units:      k^-1
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
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
!       level_p_ad:                profile set layer interface pressure adjoint
!                                  array.
!                                  units:      hpa^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_p_ad:                profile set layer average pressure adjoint
!                                  array.
!                                  units:      hpa^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_t_ad:                profile set layer average temperature adjoint
!                                  array.
!                                  units:      k^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_w_ad:      .         profile set layer average water vapor mixing ratio
!                                  adjoint array.
!                                  units:      kg/g
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       layer_o_ad:                profile set layer average ozone mixing ratio
!                                  adjoint array.
!                                  units:      ppmv^-1
!                                  type:       real
!                                  dimension:  k x m; k > 1, and m > or = 1
!                                  attributes: intent( in out )
!
!       surface_temperature_ad:    profile set surface temperature adjoint array.
!                                  units:      k^-1
!                                  type:       real
!                                  dimension:  m; m > or = 1 (i.e. scalar)
!                                  attributes: intent( in out )
!
!       surface_emissivity_ad:     profile set surface emissivity adjoint array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       surface_reflectivity_ad:   profile set surface reflectivity adjoint array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m; l > 1, and m > or = 1
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
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
!      compute_absorber_amount_ad: subroutine to calculate the adjoint
!                                  absorber profiles
!                                  source: absorber_profile module
!
!      compute_predictors_ad:      subroutine to compute the adjoint 
!                                  transmittance predictor profiles.
!                                  source: predictor module
!
!      compute_transmittance_ad:   subroutine to compute the adjoint
!                                  transmittance profiles.
!                                  source: transmittance module
!
!      compute_radiance_ad:        subroutine to compute the toa adjoint 
!                                  radiances and brightness temperatures.
!                                  source: radiance module
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       all input adjoint arguments are set to zero on output.
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

  function adjoint_rtm_rank2( &

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

             ! -- adjoint inputs
             tau_ad,                                                     &  ! in/output, k x l*m
             flux_tau_ad,                                                &  ! in/output, k x l*m
             solar_tau_ad,                                               &  ! in/output, k x l*m

             layer_radiance_ad,                                          &  ! in/output, k x l*m
             downwelling_radiance_ad,                                    &  ! in/output, l*m
             upwelling_radiance_ad,                                      &  ! in/output, l*m

             brightness_temperature_ad,                                  &  ! in/output, l*m

             ! -- other inputs
             secant_view_angle,                                          &  ! input, m
             secant_solar_angle,                                         &  ! input, m
             n_channels_per_profile,                                     &  ! input, m
             channel_index,                                              &  ! input, l*m

             ! -- adjoint outputs
             level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output, k x m

             surface_temperature_ad,                                     &  ! in/output, m
             surface_emissivity_ad,                                      &  ! in/output, l*m
             surface_reflectivity_ad,                                    &  ! in/output, l*m

             ! -- optional inputs
             n_input_profiles,                                           &
             message_log )                                               &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations  --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward inputs
    real( fp_kind ), dimension( :, : ),     intent( in )     :: level_p                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_p                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_t                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_w                    ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_o                    ! k x m

    real( fp_kind ), dimension( : ),        intent( in )     :: surface_temperature        ! m
    real( fp_kind ), dimension( : ),        intent( in )     :: surface_emissivity         ! l*m
    real( fp_kind ), dimension( : ),        intent( in )     :: surface_reflectivity       ! l*m

    real( fp_kind ), dimension( 0:, :, : ), intent( in )     :: absorber                   ! 0:k x j x m

    integer,         dimension(  :, :, : ), intent( in )     :: tau_layer_index            ! k x j x m
    integer,         dimension(  :, :, : ), intent( in )     :: flux_tau_layer_index       ! k x j x m
    integer,         dimension(  :, :, : ), intent( in )     :: solar_tau_layer_index      ! k x j x m

    real( fp_kind ), dimension( :, :, : ),  intent( in )     :: tau_predictor              ! imax x k x m
    real( fp_kind ), dimension( :, :, : ),  intent( in )     :: flux_tau_predictor         ! imax x k x m
    real( fp_kind ), dimension( :, :, : ),  intent( in )     :: solar_tau_predictor        ! imax x k x m

    real( fp_kind ), dimension( :, : ),     intent( in )     :: tau                        ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: flux_tau                   ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in )     :: solar_tau                  ! k x l*m

    real( fp_kind ), dimension( :, : ),     intent( in )     :: layer_radiance             ! k x l*m
    real( fp_kind ), dimension( : ),        intent( in )     :: downwelling_radiance       ! l*m
    real( fp_kind ), dimension( : ),        intent( in )     :: upwelling_radiance         ! l*m

    ! -- adjoint inputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: tau_ad                     ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: flux_tau_ad                ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: solar_tau_ad               ! k x l*m

    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_radiance_ad          ! k x l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: downwelling_radiance_ad    ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: upwelling_radiance_ad      ! l*m

    real( fp_kind ), dimension( : ),        intent( in out ) :: brightness_temperature_ad  ! l*m

    ! -- other inputs
    real( fp_kind ), dimension( : ),        intent( in )     :: secant_view_angle          ! m
    real( fp_kind ), dimension( : ),        intent( in )     :: secant_solar_angle         ! m
    integer,         dimension( : ),        intent( in )     :: n_channels_per_profile     ! m
    integer,         dimension( : ),        intent( in )     :: channel_index              ! l*m

    ! -- adjoint outputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: level_p_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_p_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_t_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_w_ad                 ! k x m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_o_ad                 ! k x m

    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_temperature_ad     ! m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_emissivity_ad      ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_reflectivity_ad    ! l*m

    ! -- optional input
    integer,        optional,               intent( in )     :: n_input_profiles           ! scalar
    character( * ), optional,               intent( in )     :: message_log
    
    
    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'adjoint_rtm_rank2'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    character( 100 ) :: message
    character( 5 )   :: value_in, value_allowed

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

      error_status = adjoint_rtm_rank1( &
                       ! -- forward inputs
                       level_p( :, m ), layer_p( :, m ), layer_t( :, m ), layer_w( :, m ), layer_o( :, m ), &  ! input, k

                       surface_temperature) m ),          &  ! input, scalar
                       surface_emissivity( l1:l2),        &  ! input, l
                       surface_reflectivity( l1:l2),      &  ! input, l

                       absorber( 0:, :, m ),              &  ! input, 0:k x j

                       tau_layer_index( :, :, m ),        &  ! input, k x j
                       flux_tau_layer_index( :, :, m ),   &  ! input, k x j
                       solar_tau_layer_index( :, :, m ),  &  ! input, k x j

                       tau_predictor( :, :, m ),          &  ! input, imax x k
                       flux_tau_predictor( :, :, m ),     &  ! input, imax x k
                       solar_tau_predictor( :, :, m ),    &  ! input, imax x k

                       tau( :, l1:l2 ),                   &  ! input, k x l
                       flux_tau( :, l1:l2 ),              &  ! input, k x l
                       solar_tau( :, l1:l2 ),             &  ! input, k x l

                       layer_radiance( :, l1:l2 ),        &  ! input, k x l
                       downwelling_radiance( l1:l2),      &  ! input, l
                       upwelling_radiance( l1:l2),        &  ! input, l

                       ! -- adjoint inputs
                       tau_ad( :, l1:l2 ),                &  ! in/output, k x l
                       flux_tau_ad( :, l1:l2 ),           &  ! in/output, k x l
                       solar_tau_ad( :, l1:l2 ),          &  ! in/output, k x l

                       layer_radiance_ad( :, l1:l2 ),     &  ! in/output, k x l
                       downwelling_radiance_ad( l1:l2),   &  ! in/output, l
                       upwelling_radiance_ad( l1:l2),     &  ! in/output, l

                       brightness_temperature_ad( l1:l2), &  ! in/output, l

                       ! -- other inputs
                       secant_view_angle( m ),            &  ! input, scalar
                       secant_solar_angle( m ),           &  ! input, scalar
                       n_channels_per_profile( m ),       &  ! input, scalar
                       channel_index( l1:l2),             &  ! input, l

                       ! -- adjoint outputs
                       level_p_ad( :, m ), layer_p_ad( :, m ), layer_t_ad( :, m ), layer_w_ad( :, m ), layer_o_ad( :, m ), &  ! in/output, k

                       surface_temperature_ad( m ),       &  ! in/output, scalar
                       surface_emissivity_ad( l1:l2),     &  ! in/output, l
                       surface_reflectivity_ad( l1:l2),   &  ! in/output, l

                       ! -- optional inputs
                       message_log = message_log )


      ! -------------------------------
      ! check for successful completion
      ! -------------------------------

      if ( error_status /= success ) then

        call display_message( routine_name, &
                              'error occured in adjoint_rtm_rank1', &
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

  end function adjoint_rtm_rank2



  function adjoint_rtm_rank1( &

             ! -- forward inputs
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

             ! -- adjoint inputs
             tau_ad,                                                     &  ! in/output, k x l
             flux_tau_ad,                                                &  ! in/output, k x l
             solar_tau_ad,                                               &  ! in/output, k x l

             layer_radiance_ad,                                          &  ! in/output, k x l
             downwelling_radiance_ad,                                    &  ! in/output, l
             upwelling_radiance_ad,                                      &  ! in/output, l

             brightness_temperature_ad,                                  &  ! in/output, l

             ! -- other inputs
             secant_view_angle,                                          &  ! input, scalar
             secant_solar_angle,                                         &  ! input, scalar
             n_channels,                                                 &  ! input, scalar
             channel_index,                                              &  ! input, l

             ! -- adjoint outputs
             level_p_ad, layer_p_ad, layer_t_ad, layer_w_ad, layer_o_ad, &  ! in/output, k

             surface_temperature_ad,                                     &  ! in/output, scalar
             surface_emissivity_ad,                                      &  ! in/output, l
             surface_reflectivity_ad,                                    &  ! in/output, l

             ! optional inputs
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
    real( fp_kind ), dimension( : ),     intent( in )     :: level_p                    ! k
    real( fp_kind ), dimension( : ),     intent( in )     :: layer_p                    ! k
    real( fp_kind ), dimension( : ),     intent( in )     :: layer_t                    ! k
    real( fp_kind ), dimension( : ),     intent( in )     :: layer_w                    ! k
    real( fp_kind ), dimension( : ),     intent( in )     :: layer_o                    ! k

    real( fp_kind ),                     intent( in )     :: surface_temperature        ! scalar
    real( fp_kind ), dimension( : ),     intent( in )     :: surface_emissivity         ! l
    real( fp_kind ), dimension( : ),     intent( in )     :: surface_reflectivity       ! l

    real( fp_kind ), dimension( 0:, : ), intent( in )     :: absorber                   ! 0:k x j

    integer,         dimension(  :, : ), intent( in )     :: tau_layer_index            ! k x j
    integer,         dimension(  :, : ), intent( in )     :: flux_tau_layer_index       ! k x j
    integer,         dimension(  :, : ), intent( in )     :: solar_tau_layer_index      ! k x j

    real( fp_kind ), dimension( :, : ),  intent( in )     :: tau_predictor              ! imax x k 
    real( fp_kind ), dimension( :, : ),  intent( in )     :: flux_tau_predictor         ! imax x k 
    real( fp_kind ), dimension( :, : ),  intent( in )     :: solar_tau_predictor        ! imax x k 

    real( fp_kind ), dimension( :, : ),  intent( in )     :: tau                        ! k x l
    real( fp_kind ), dimension( :, : ),  intent( in )     :: flux_tau                   ! k x l
    real( fp_kind ), dimension( :, : ),  intent( in )     :: solar_tau                  ! k x l

    real( fp_kind ), dimension( :, : ),  intent( in )     :: layer_radiance             ! k x l
    real( fp_kind ), dimension( : ),     intent( in )     :: downwelling_radiance       ! l
    real( fp_kind ), dimension( : ),     intent( in )     :: upwelling_radiance         ! l

    ! -- adjoint inputs
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: tau_ad                     ! k x l
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: flux_tau_ad                ! k x l
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: solar_tau_ad               ! k x l

    real( fp_kind ), dimension( :, : ),  intent( in out ) :: layer_radiance_ad          ! k x l
    real( fp_kind ), dimension( : ),     intent( in out ) :: downwelling_radiance_ad    ! l
    real( fp_kind ), dimension( : ),     intent( in out ) :: upwelling_radiance_ad      ! l

    real( fp_kind ), dimension( : ),     intent( in out ) :: brightness_temperature_ad  ! l

    ! -- other inputs
    real( fp_kind ),                     intent( in )     :: secant_view_angle          ! scalar
    real( fp_kind ),                     intent( in )     :: secant_solar_angle         ! scalar
    integer,                             intent( in )     :: n_channels                 ! scalar
    integer,         dimension( : ),     intent( in )     :: channel_index              ! l

    ! -- adjoint outputs
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: level_p_ad                 ! k
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: layer_p_ad                 ! k
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: layer_t_ad                 ! k
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: layer_w_ad                 ! k
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: layer_o_ad                 ! k

    real( fp_kind ),                     intent( in out ) :: surface_temperature_ad     ! scalar
    real( fp_kind ), dimension( : ),     intent( in out ) :: surface_emissivity_ad      ! l
    real( fp_kind ), dimension( : ),     intent( in out ) :: surface_reflectivity_ad    ! l

    ! -- optional input. note that n_input_profiles is not used in this
    ! -- function. it is specified so if a user specifies it by mistake
    ! -- for rank-1 profile input the code won't (hopefully) fall in a heap.
    integer,        optional,            intent( in )  :: n_input_profiles              ! scalar
    character( * ), optional,            intent( in )  :: message_log
    
    
    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'adjoint_rtm_rank1'


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

    ! -- arrays for adjoints of integrated absorber amounts, 0:k x j
    real( fp_kind ), dimension( 0:size( absorber, dim = 1 )-1, &
                                  size( absorber, dim = 2 )    ) :: absorber_ad, &
                                                                    tau_absorber_ad,      &
                                                                    flux_tau_absorber_ad, &
                                                                    solar_tau_absorber_ad

    ! -- arrays for adjoint predictors, imax x k
    real( fp_kind ), dimension( size( tau_predictor, dim = 1 ), &
                                size( tau_predictor, dim = 2 )  ) :: tau_predictor_ad,      &
                                                                     flux_tau_predictor_ad, &
                                                                     solar_tau_predictor_ad

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
    !#           -- scale absorber amounts and initialise adjoints --           #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------------------
    ! modify absorber quantities by the angle secant
    ! could put a loop here but here's hoping the compiler
    ! recognises this as a group of loops over layer.
    ! ----------------------------------------------------

    ! -- upwelling transmittance
    tau_absorber( 0:, : ) = secant_view_angle * absorber( 0:, : )
    tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - toa_pressure

    ! -- flux transmittance
    if ( any( is_microwave_channel( channel_index( 1:n_channels ) ) == 0 ) ) then
      flux_tau_absorber( 0:, : ) = secant_diffusivity_angle * absorber( 0:, : )
      flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - toa_pressure
    end if

    ! -- solar transmittance
    if ( ( any( is_solar_channel( channel_index( 1:n_channels ) ) == 1 ) ) .and. &
         secant_solar_angle < max_secant_solar_angle ) then
      solar_tau_absorber( 0:, : ) = secant_solar_angle( m ) * absorber( 0:, : )
      solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - toa_pressure
    end if


    ! ----------------------------------
    ! initialise local adjoint variables
    ! ----------------------------------

    ! -- absorber arrays, 0:k x j
    absorber_ad( 0:, : )           = zero
    tau_absorber_ad( 0:, : )       = zero
    flux_tau_absorber_ad( 0:, : )  = zero
    solar_tau_absorber_ad( 0:, : ) = zero


    ! -- predictor arrays, imax x k
    tau_predictor_ad( :, : )       = zero
    flux_tau_predictor_ad( :, : )  = zero
    solar_tau_predictor_ad( :, : ) = zero 



    !#--------------------------------------------------------------------------#
    !#                            -- channel loop --                            #
    !#--------------------------------------------------------------------------#

    l_channel_loop: do l = 1, n_channels


      ! ----------------------------------------------------
      ! set the "this is a channel influenced by solar" flag
      ! ----------------------------------------------------

      valid_solar = 0

      if ( is_solar_channel( channel_index( l ) ) == 1 .and. &
           secant_solar_angle < max_secant_solar_angle       ) valid_solar = 1



      ! ------------------------------------------------
      ! calculate the adjoint of the current channel toa
      ! radiance or brightness temperature
      ! ------------------------------------------------

      call compute_radiance_ad( &
                                ! -- forward input
                                layer_t( : ),                   &  ! input, k

                                surface_temperature,            &  ! input, scalar
                                surface_emissivity( l ),        &  ! input, scalar
                                surface_reflectivity( l ),      &  ! input, scalar

                                tau(      :, l ),               &  ! input, k
                                flux_tau( :, l ),               &  ! input, k
                                solar_tau( n_layers, l ),       &  ! input, scalar

                                layer_radiance( :, l ),         &  ! input, k
                                downwelling_radiance( l ),      &  ! input, scalar
                                upwelling_radiance( l ),        &  ! input, scalar

                                ! -- adjoint input
                                layer_radiance_ad( :, l ),      &  ! in/output, k
                                downwelling_radiance_ad( l ),   &  ! in/output, scalar
                                upwelling_radiance_ad( l ),     &  ! in/output, scalar

                                brightness_temperature_ad( l ), &  ! in/output, scalar

                                ! -- other input
                                secant_solar_angle,             &  ! input, scalar
                                valid_solar,                    &  ! input, scalar
                                channel_index( l ),             &  ! input, scalar

                                ! -- adjoint output
                                layer_t_ad( : ),                &  ! in/output, k

                                surface_temperature_ad,         &  ! in/output, scalar
                                surface_emissivity_ad( l ),     &  ! in/output, scalar
                                surface_reflectivity_ad( l ),   &  ! in/output, scalar

                                tau_ad( :, l ),                 &  ! in/output, k
                                flux_tau_ad( :, l ),            &  ! in/output, k
                                solar_tau_ad( n_layers, l )     )  ! in/output, scalar


      ! ----------------------------------------------------
      ! if the current channel is a solar sensitive channel,
      !   and
      ! the solar angle is valid, then calculate the adjoint
      ! of the transmittance for direct solar.
      ! ----------------------------------------------------

      if ( valid_solar == 1 ) then

        call compute_transmittance_ad( &
                                       ! -- forward input
                                       solar_tau_absorber( 0:, : ),    &   ! input, 0:k x j
                                       solar_tau_predictor( :, : ),    &   ! input, i x k
                                       solar_tau( :, l ),              &   ! input, k

                                       ! -- adjoint input
                                       solar_tau_ad( :, l ),           &   ! in/output, k

                                       ! -- other input
                                       solar_tau_layer_index( :, : ),  &   ! input, k x j
                                       channel_index( l ),             &   ! input, scalar
                                       down,                           &   ! input, scalar

                                       ! -- adjoint output
                                       solar_tau_absorber_ad( 0:, : ), &   ! in/output, 0:k x j
                                       solar_tau_predictor_ad( :, : )  )   ! in/output, i x k

      end if




      ! --------------------------------------------------
      ! if the current channel is an infrared channel,
      ! then calculate the adjoint of the downwelling flux
      ! transmittance using the flux absorber amounts.
      !
      ! if the current channel is a microwave channel,
      ! then calculate the adjoint of the flux transmittance
      ! using the upwelling absorber amounts.
      ! --------------------------------------------------

      if ( is_microwave_channel( channel_index( l ) ) == 0 ) then

        ! -- ir channel
        call compute_transmittance_ad( &
                                       ! -- forward input
                                       flux_tau_absorber( 0:, : ),    &   ! input, 0:k x j
                                       flux_tau_predictor( :, : ),    &   ! input, i x k
                                       flux_tau( :, l ),              &   ! input, k

                                       ! -- adjoint input
                                       flux_tau_ad( :, l ),           &   ! in/output, k

                                       ! -- other input
                                       flux_tau_layer_index( :, : ),  &   ! input, k x j
                                       channel_index( l ),            &   ! input, scalar
                                       down,                          &   ! input, scalar

                                       ! -- adjoint output
                                       flux_tau_absorber_ad( 0:, : ), &   ! in/output, 0:k x j
                                       flux_tau_predictor_ad( :, : )  )   ! in/output, i x k

      else

        ! -- uw channel. use upwelling transmittance predictors
!        call compute_transmittance_ad( &
!                                       ! -- forward input
!                                       tau_absorber( 0:, : ),    &   ! input, 0:k x j
!                                       tau_predictor( :, : ),    &   ! input, i x k
!                                       flux_tau( :, l ),         &   ! input, k
!
!                                       ! -- adjoint input
!                                       flux_tau_ad( :, l ),      &   ! in/output, k
!
!                                       ! -- other input
!                                       tau_layer_index( :, : ),  &   ! input, k x j
!                                       channel_index( l ),       &   ! input, scalar
!                                       down,                     &   ! input, scalar
!
!                                       ! -- adjoint output
!                                       tau_absorber_ad( 0:, : ), &   ! in/output, 0:k x j
!                                       tau_predictor_ad( :, : )  )   ! in/output, i x k

        ! -- if total transmittance /= 0, calculate 
        ! -- adjoint flux transmittance
        if ( tau( n_layers, l ) > tolerance ) then
          do k = n_layers, 2, -1
            flux_tau_ad( 1, l ) = flux_tau_ad( 1, l ) + ( flux_tau_ad( k, l ) / &
            !                                             -------------------
                                                              tau( k-1, l )   )

            tau_ad( k-1, l ) = tau_ad( k-1, l ) - ( flux_tau_ad( k, l ) * flux_tau( 1, l ) / &
            !                                       --------------------------------------
                                                                 tau( k-1, l )**2          )
            flux_tau_ad( k, l ) = zero
          end do
          tau_ad( n_layers, l ) = tau_ad( n_layers, l ) + flux_tau_ad( 1, l )
          flux_tau_ad( 1, l ) = zero
        else
          flux_tau_ad( :, l ) = zero
        end if

      end if


      ! -----------------------------------------------------
      ! calculate the adjoint of the upwelling transmittances
      ! for the satellite view angle
      ! -----------------------------------------------------

      call compute_transmittance_ad( &
                                     ! -- forward input
                                     tau_absorber( 0:, : ),    &   ! input, 0:k x j
                                     tau_predictor( :, : ),    &   ! input, i x k
                                     tau( :, l ),              &   ! input, k

                                     ! -- adjoint input
                                     tau_ad( :, l ),           &   ! in/output, k

                                     ! -- other input
                                     tau_layer_index( :, : ),  &   ! input, k x j
                                     channel_index( l ),       &   ! input, scalar
                                     up,                       &   ! input, scalar

                                     ! -- adjoint output
                                     tau_absorber_ad( 0:, : ), &   ! in/output, 0:k x j
                                     tau_predictor_ad( :, : )  )   ! in/output, i x k

    end do l_channel_loop



    !#--------------------------------------------------------------------------#
    !#    -- calculate the predictor adjoints for the solar transmittance --    #
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


      ! --------------------------------------
      ! adjoint of the standard predictor copy
      ! --------------------------------------

      tau_predictor_ad( 1:max_n_standard_predictors, : ) = &
              tau_predictor_ad( 1:max_n_standard_predictors, : ) + &
        solar_tau_predictor_ad( 1:max_n_standard_predictors, : )

      solar_tau_predictor_ad( 1:max_n_standard_predictors, : ) = zero


      ! -------------------------------------
      ! adjoint of the integrateed predictors
      ! -------------------------------------

      call compute_predictors_ad( &
                                  ! -- forward input
                                  layer_p( : ),                   &  ! input,  k
                                  layer_t( : ),                   &  ! input,  k
                                  layer_w( : ),                   &  ! input,  k
                                  solar_tau_absorber( 0:, : ),    &  ! input,  0:k x j

                                  ! -- adjoint input
                                  solar_tau_predictor_ad( :, : ), &  ! in/output, i x k

                                  ! -- adjoint output
                                  layer_p_ad( : ),                &  ! in/output,  k
                                  layer_t_ad( : ),                &  ! in/output,  k
                                  layer_w_ad( : ),                &  ! in/output,  k
                                  solar_tau_absorber_ad( 0:, : ), &  ! in/output,  0:k x j

                                  no_standard = 1                 )  ! optional input


      ! -------------------------------------------------
      ! adjoint of the nadir absorber amount modification
      ! -------------------------------------------------

      absorber_ad( 0:, : ) = absorber_ad( 0:, : ) + &
                             ( secant_solar_angle * solar_tau_absorber_ad( 0:, : ) )

    end if



    !#--------------------------------------------------------------------------#
    !#     -- calculate the predictor adjoints for the flux transmittance --    #
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


      ! --------------------------------------
      ! adjoint of the standard predictor copy
      ! --------------------------------------

      tau_predictor_ad( 1:max_n_standard_predictors, : ) = &
             tau_predictor_ad( 1:max_n_standard_predictors, : ) + &
        flux_tau_predictor_ad( 1:max_n_standard_predictors, : )

      flux_tau_predictor_ad( 1:max_n_standard_predictors, : ) = zero


      ! -------------------------------------
      ! adjoint of the integrateed predictors
      ! -------------------------------------

      call compute_predictors_ad( &
                                  ! -- forward input
                                  layer_p( : ),                  &  ! input,  k
                                  layer_t( : ),                  &  ! input,  k
                                  layer_w( : ),                  &  ! input,  k
                                  flux_tau_absorber( 0:, : ),    &  ! input,  0:k x j

                                  ! -- adjoint input
                                  flux_tau_predictor_ad( :, : ), &  ! in/output, i x k

                                  ! -- adjoint output
                                  layer_p_ad( : ),               &  ! in/output,  k
                                  layer_t_ad( : ),               &  ! in/output,  k
                                  layer_w_ad( : ),               &  ! in/output,  k
                                  flux_tau_absorber_ad( 0:, : ), &  ! in/output,  0:k x j

                                  no_standard = 1                )  ! optional input


      ! -------------------------------------------------
      ! adjoint of the nadir absorber amount modification
      ! -------------------------------------------------

      absorber_ad( 0:, : ) = absorber_ad( 0:, : ) + &
                             ( secant_diffusivity_angle * flux_tau_absorber_ad( 0:, : ) )


    end if



    !#--------------------------------------------------------------------------#
    !#  -- calculate the predictor adjoints for the upwelling transmittance --  #
    !#--------------------------------------------------------------------------#

    ! -----------------------------
    ! adjoint of all the predictors
    ! -----------------------------

    call compute_predictors_ad( &
                                ! -- forward input
                                layer_p( : ),             &  ! input,  k
                                layer_t( : ),             &  ! input,  k
                                layer_w( : ),             &  ! input,  k
                                tau_absorber( 0:, : ),    &  ! input,  0:k x j

                                ! -- adjoint input
                                tau_predictor_ad( :, : ), &  ! in/output, i x k

                                ! -- adjoint output
                                layer_p_ad( : ),          &  ! in/output,  k
                                layer_t_ad( : ),          &  ! in/output,  k
                                layer_w_ad( : ),          &  ! in/output,  k
                                tau_absorber_ad( 0:, : )  )  ! in/output,  0:k x j


    ! -------------------------------------------------
    ! adjoint of the nadir absorber amount modification
    ! -------------------------------------------------

    absorber_ad( 0:, : ) = absorber_ad( 0:, : ) + &
                           ( secant_view_angle * tau_absorber_ad( 0:, : ) )



    !#--------------------------------------------------------------------------#
    !#   -- calculate the absorber adjoints for the upwelling transmittance --  #
    !#--------------------------------------------------------------------------#

    call compute_absorber_amount_ad( &
                                     ! -- forward input
                                     level_p( : ),         &  ! input,  k
                                     layer_w( : ),         &  ! input,  k
                                     layer_o( : ),         &  ! input,  k

                                     ! -- adjoint input
                                     absorber_ad( 0:, : ), &  ! in/output, 0:k x j

                                     ! -- adjoint output
                                     level_p_ad( : ),      &  ! in/ouput,  k
                                     layer_w_ad( : ),      &  ! in/ouput,  k
                                     layer_o_ad( : )       )  ! in/ouput,  k



    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function adjoint_rtm_rank1

end module adjoint_model


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
! revision 1.9  2001/11/07 15:08:08  paulv
! - changed the logical if test for the number of input profiles in the
!   [<tangent_linear><adjoint><k_matrix>]_rtm_rank2() functions from
!     if ( n_input_profiles >= 1 .and. n_input_profiles <= n_profiles ) then
!   to
!     if ( n_input_profiles > 0 .and. n_input_profiles <= n_profiles ) then
!   the use of both the ">=" and "<=" relational operators with the .and. i
!   found confusing.
!
! revision 1.8  2001/11/06 20:59:40  paulv
! - corrected adjoint variable units documentation.
! - added check for negative number of channels to adjoint_rtm_rank2() function.
! - added profile loop cycle statement to adjoint_rtm_rank2() function.
! - added check for negative number of channels to adjoint_rtm_rank1() function.
!
! revision 1.7  2001/10/01 20:13:28  paulv
! - overloaded the compute_rtm_k and adjoint_rtm functions to accept both a
!   single or group of profiles. contained private functions are now
!     o compute_rtm_ad_rank1 and adjoint_rtm_rank1 for single profile input
!     o compute_rtm_ad_rank2 and adjoint_rtm_rank2 for multiple profile input
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
!   the check for n_input_profiles is only performed in the adjoint_rtm_rank2
!   function.
! - changed surface_temperature argument check from
!     if ( any( surface_temperature < zero ) )
!   to
!     if (      surface_temperature < zero )
!   as the check is now done in the rank-1 adjoint_rtm function. this eliminates
!   the need to fully populate the input arrays with data when only certain
!   "chunks" may be processed. previously, the use of any() could generate
!   and error if the full surface_temperature array was not initialised.
! - added "name" to rcs keyword list.
!
! revision 1.6  2001/09/28 22:56:40  paulv
! - overloaded the adjoint_rtm function to accept both a
!   single or group of profiles. contained private functions are now
!     o adjoint_rtm_rank1 for single profile input
!     o adjoint_rtm_rank2 for multiple profile input
! - put n_input_profiles optional argument back in the adjoint_rtm
!   argument list. this allows more control over how many profiles are to be
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
!   as the check is now done in the rank-1 adjoint_rtm function. this eliminates
!   the need to fully populate the input arrays with data when only certain
!   "chunks" may be processed. previously, the use of any() could generate
!   and error if the full surface_temperature array was not initialised.
! - added "name" to rcs keyword list.
!
! revision 1.5  2001/08/31 21:31:33  paulv
! - added check for negative profile and surface data in adjoint_rtm.
! - maximum solar angle secant is no longer calculated in adjoint_rtm but
!   is declared as a parameter in the parameters module.
! - added compute_rtm_ad function. this is a wrapper for the adjoint_rtm function.
!
! revision 1.4  2001/08/16 16:35:16  paulv
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
! revision 1.3  2001/08/01 16:40:12  paulv
! - added documentation.
! - altered function declaration to avoid more than the standard 39
!   continuation lines allowed in fortran 90.
! - updated input argument checking. now consistent with other model
!   components.
!
! revision 1.2  2001/07/12 19:01:03  paulv
! - added absorber quantity modification for the various zenith angles (view,
!   diffusivity, and solar.)
! - added local adjoint variable initialisation.
! - corrected bug in main transmittance adjoint calculation. the
!   direction argument was set to down when it should be up.
! - corrected bug in calculation of the predictor adjoint for the ir flux
!   transmittance. the forward input absorber amount argument was specified
!   as solar_tau_absorber rather than flux_tau_absorber.
!
! revision 1.1  2001/05/29 16:33:55  paulv
! initial checkin
!
!
!
