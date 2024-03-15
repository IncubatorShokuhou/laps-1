!------------------------------------------------------------------------------
!m+
! name:
!       k_matrix_model
!
! purpose:
!       module containing the ncep rt k-matrix model functions
!
! category:
!       ncep rtm
!
! calling sequence:
!       use k_matrix_model
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
!       forward_model:         module containing the forward model function.
!
! contains:
!       compute_rtm_k:         public function that calculates the k-matrix of the 
!                              top-of-atmosphere (toa) radiances and brightness 
!                              temperatures for an input atmospheric profile set and
!                              user specified satellites/channels.
!
!                              this function is simply a wrapper around both the forward
!                              model and the k-matrix model so that the user doesn't have
!                              to declare the absorber/predictor/etc. arrays in the calling
!                              routine.
!
!       k_matrix_rtm:          public function that calculates the k-matrix of the
!                              top-of-atmosphere (toa) radiances and brightness
!                              temperatures for user specified profiles and
!                              satellite/channel transmittance profiles, radiances
!                              and brightness temperatures.
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
!       written by:     paul van delst, cimss@noaa/ncep 18-jul-2001
!                       pvandelst@ncep.noaa.gov
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
!m-
!------------------------------------------------------------------------------

module k_matrix_model


  ! ------------
  ! module usage
  ! ------------

  use type_kinds,            only : fp_kind
  use error_handler
  use parameters
  use spectral_coefficients, only : is_solar_channel, &
                                    is_microwave_channel
  use absorber_profile,      only : compute_absorber_amount_ad
  use predictors,            only : compute_predictors_ad
  use transmittance,         only : compute_transmittance_ad
  use radiance,              only : compute_radiance_ad
  use forward_model,         only : forward_rtm


  ! -----------------------
  ! disable implicit typing
  ! -----------------------

  implicit none


  ! --------------------
  ! default visibilities
  ! --------------------

  private
  public :: compute_rtm_k
  public :: k_matrix_rtm


contains





!--------------------------------------------------------------------------------
!s+
! name:
!       compute_rtm_k
!
! purpose:
!       public function that calculates the k-matrix of the top-of-atmosphere (toa)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
!       this function is simply a wrapper around both the forward model and the
!       k-matrix model so that the user doesn't have to declare the absorber/
!       predictor/etc arrays in the calling routine.
!
! category:
!       ncep rtm
!
! calling sequence:
!       result = compute_rtm_k( &
!                               ! -- forward inputs
!                               level_p, layer_p, layer_t, layer_w, layer_o,           &  ! input, k x m
!
!                               surface_temperature,                                   &  ! input, m
!                               surface_emissivity,                                    &  ! input, l*m
!                               surface_reflectivity,                                  &  ! input, l*m
!
!                               ! -- k-matrix inputs
!                               tau_k,                                                 &  ! in/output, k x l*m
!                               flux_tau_k,                                            &  ! in/output, k x l*m
!                               solar_tau_k,                                           &  ! in/output, k x l*m
!
!                               upwelling_radiance_k,                                  &  ! in/output, l*m
!                               brightness_temperature_k,                              &  ! in/output, l*m
!
!                               ! -- other inputs
!                               secant_view_angle,                                     &  ! input, m
!                               secant_solar_angle,                                    &  ! input, m
!                               n_channels_per_profile,                                &  ! input, m
!                               channel_index,                                         &  ! input, l*m
!
!                               ! -- forward output
!                               tau,                                                   &  ! input, k x l*m
!                               flux_tau,                                              &  ! input, k x l*m
!                               solar_tau,                                             &  ! input, k x l*m
!
!                               upwelling_radiance,                                    &  ! input, l*m
!                               brightness_temperature,                                &  ! input, l*m
!
!                               ! -- k-matrix outputs
!                               level_p_k, layer_p_k, layer_t_k, layer_w_k, layer_o_k, &  ! in/output, k x l*m
!
!                               surface_temperature_k,                                 &  ! in/output, l*m
!                               surface_emissivity_k,                                  &  ! in/output, l*m
!                               surface_reflectivity_k,                                &  ! in/output, l*m
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
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_p:                   profile set layer average pressure array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_t:                   profile set layer average temperature array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_w:      .            profile set layer average water vapor mixing ratio array
!                                  units:      g/kg
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_o:                   profile set layer average ozone mixing ratio array.
!                                  units:      ppmv
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       surface_temperature:       profile set surface temperature array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       surface_emissivity:        profile set surface emissivity array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       surface_reflectivity:      profile set surface reflectivity array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       tau_k:                     layer->toa adjoint transmittance for the satellite
!                                  view angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in out )
!
!       flux_tau_k:                layer->sfc adjoint transmittance for the default
!                                  diffusivity angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in out )
!
!       solar_tau_k:               layer->sfc adjoint transmittance for the solar
!                                  zenith angle.
!                                  ** this argument is set to zero on output **.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in out )
!
!       upwelling_radiance_k:      toa adjoint radiances for each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  ** this argument is set to zero on output **.
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       brightness_temperature_k:  adjoint temperatures corresponding to the
!                                  toa adjoint radiances.
!                                  ** this argument is set to zero on output **.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       secant_view_angle:         secant of the satellite view angle measured
!                                  from nadir for each profile in the set.
!                                  units:      none
!                                  type:       real
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       secant_solar_angle:        secant of the solar zenith angle for each
!                                  profile in the set.
!                                  units:      none
!                                  type:       real
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       n_channels_per_profile:    the number of channels for each profile in the
!                                  set for which radiances are required.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       channel_index:             channel index id array. each element is a unique
!                                  index to a (supported) sensor channel.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  l*m
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
!                                  dimension:  k x l*m
!                                  attributes: intent( out )
!
!       flux_tau:                  layer->sfc transmittance for the default
!                                  diffusivity angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( out )
!
!       solar_tau:                 layer->sfc transmittance for the solar
!                                  zenith angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( out )
!
!       upwelling_radiance:        toa radiances for each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
!       brightness_temperature:    toa brightness temperatures corresponding
!                                  to the toa radiances.
!                                  n.b.: set to zero upon output.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( out )
!
!       level_p_k:                 profile set layer interface pressure k-matrix
!                                  adjoint array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_p_k:                 profile set layer average pressure k-matrix
!                                  adjoint array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_t_k:                 profile set layer average temperature k-matrix
!                                  adjoint array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_w_k:      .          profile set layer average water vapor mixing ratio
!                                  k-matrix adjoint array.
!                                  units:      g/kg
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_o_k:                 profile set layer average ozone mixing ratio
!                                  k-matrix adjoint array.
!                                  units:      ppmv
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       surface_temperature_k:     profile set surface temperature k-matrix adjoint
!                                  array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       surface_emissivity_k:      profile set surface emissivity k-matrix adjoint
!                                  array.
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       surface_reflectivity_k:    profile set surface reflectivity k-matrix adjoint
!                                  array.
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
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
!      k_matrix_rtm:               function that calculates the k-matrix of the (toa)
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
!       the code has is not overloaded for scalar input so the input
!       arguments must be dimensioned accordingly, even if only one
!       profile or channel is being passed.
!
! procedure:
!       see individual module function documentation.
!s-
!--------------------------------------------------------------------------------

  function compute_rtm_k( &
             ! -- forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,           &  ! input, k x m

             surface_temperature,                                   &  ! input, m
             surface_emissivity,                                    &  ! input, l*m
             surface_reflectivity,                                  &  ! input, l*m

             ! -- k-matrix inputs
             tau_k,                                                 &  ! in/output, k x l*m
             flux_tau_k,                                            &  ! in/output, k x l*m
             solar_tau_k,                                           &  ! in/output, k x l*m

             upwelling_radiance_k,                                  &  ! in/output, l*m
             brightness_temperature_k,                              &  ! in/output, l*m

             ! -- other inputs
             secant_view_angle,                                     &  ! input, m
             secant_solar_angle,                                    &  ! input, m
             n_channels_per_profile,                                &  ! input, m
             channel_index,                                         &  ! input, l*m

             ! -- forward output
             tau,                                                   &  ! input, k x l*m
             flux_tau,                                              &  ! input, k x l*m
             solar_tau,                                             &  ! input, k x l*m

             upwelling_radiance,                                    &  ! input, l*m
             brightness_temperature,                                &  ! input, l*m

             ! -- k-matrix outputs
             level_p_k, layer_p_k, layer_t_k, layer_w_k, layer_o_k, &  ! in/output, k x l*m

             surface_temperature_k,                                 &  ! in/output, l*m
             surface_emissivity_k,                                  &  ! in/output, l*m
             surface_reflectivity_k,                                &  ! in/output, l*m

             ! optional inputs
             message_log )                                          &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
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

    ! -- k-matrix inputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: tau_k                      ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: flux_tau_k                 ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: solar_tau_k                ! k x l*m

    real( fp_kind ), dimension( : ),        intent( in out ) :: upwelling_radiance_k       ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: brightness_temperature_k   ! l*m

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

    ! -- k-matrix outputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: level_p_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_p_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_t_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_w_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_o_k                  ! k x l*m

    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_temperature_k      ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_emissivity_k       ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_reflectivity_k     ! l*m

    ! -- optional input
    character( * ), optional,               intent( in )     :: message_log


    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_rtm_k'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    character( 100 ) :: message
    character( 5 )   :: value_in, value_allowed

    integer :: m, n_profiles          ! profile loop variables
    integer :: l, l1, l2, n_channels  ! channel loop/index variables
    integer :: k, n_layers            ! layer loop variables

    integer :: error_status_fwd
    integer :: error_status_k

    ! -- maximum channels pseudo parameter
    integer :: max_n_channels
    logical :: is_set

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

    ! -- array for forward and k-matrix layer planck radiance term, k x l*m
    real( fp_kind ), dimension( size( layer_p, dim = 1 ),  &
                                size( upwelling_radiance ) ) :: layer_radiance,  &
                                                                layer_radiance_k
      

    ! -- array for forward and k-matrix downwelling radiance (flux + solar), l*m
    real( fp_kind ), dimension( size( upwelling_radiance ) ) :: downwelling_radiance,  &
                                                                downwelling_radiance_k 



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
    !#                   -- compute the k-matrix profiles --                    #
    !#--------------------------------------------------------------------------#

    ! -----------------------------------------------
    ! initialise all local adjoint/k-matrix variables
    ! -----------------------------------------------

    layer_radiance_k( :, : )    = zero
    downwelling_radiance_k( : ) = zero


    ! -----------------------
    ! call the k-matrix model
    ! -----------------------

    error_status_k = k_matrix_rtm( &
                                   ! -- forward inputs
                                   level_p, layer_p, layer_t, layer_w, layer_o,           &  ! input,  k x m

                                   surface_temperature,                                   &  ! input, m
                                   surface_emissivity,                                    &  ! input, l*m
                                   surface_reflectivity,                                  &  ! input, l*m

                                   absorber,                                              &  ! input, 0:k x j x m

                                   tau_layer_index,                                       &  ! input, k x j x m
                                   flux_tau_layer_index,                                  &  ! input, k x j x m
                                   solar_tau_layer_index,                                 &  ! input, k x j x m

                                   tau_predictor,                                         &  ! input, imax x k x m
                                   flux_tau_predictor,                                    &  ! input, imax x k x m
                                   solar_tau_predictor,                                   &  ! input, imax x k x m

                                   tau,                                                   &  ! input, k x l*m
                                   flux_tau,                                              &  ! input, k x l*m
                                   solar_tau,                                             &  ! input, k x l*m

                                   layer_radiance,                                        &  ! input, k x l*m
                                   downwelling_radiance,                                  &  ! input, l*m
                                   upwelling_radiance,                                    &  ! input, l*m

                                   ! -- k-matrix inputs
                                   tau_k,                                                 &  ! in/output, k x l*m
                                   flux_tau_k,                                            &  ! in/output, k x l*m
                                   solar_tau_k,                                           &  ! in/output, k x l*m

                                   layer_radiance_k,                                      &  ! in/output, k x l*m
                                   downwelling_radiance_k,                                &  ! in/output, l*m
                                   upwelling_radiance_k,                                  &  ! in/output, l*m

                                   brightness_temperature_k,                              &  ! in/output, l*m

                                   ! -- other inputs
                                   secant_view_angle,                                     &  ! input, m
                                   secant_solar_angle,                                    &  ! input, m
                                   n_channels_per_profile,                                &  ! input, m
                                   channel_index,                                         &  ! input, l*m

                                   ! -- k-matrix outputs
                                   level_p_k, layer_p_k, layer_t_k, layer_w_k, layer_o_k, &  ! in/output,  k x l*m

                                   surface_temperature_k,                                 &  ! in/output, l*m
                                   surface_emissivity_k,                                  &  ! in/output, l*m
                                   surface_reflectivity_k,                                &  ! in/output, l*m

                                   message_log = message_log )


    ! -------------------------------
    ! check for successful completion
    ! -------------------------------

    if ( error_status_k /= success ) then

      error_status = failure
      call display_message( routine_name, &
                            'error occured in k_matrix_rtm', &
                            error_status, &
                            message_log = message_log )
      return

    end if


    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success


  end function compute_rtm_k





!--------------------------------------------------------------------------------
!s+
! name:
!       k_matrix_rtm
!
! purpose:
!       public function that calculates the k-matrix of the top-of-atmosphere (toa)
!       radiances and brightness temperatures for an input atmospheric profile
!       set and user specified satellites/channels.
!
! category:
!       ncep rtm
!
! calling sequence:
!       result = k_matrix_rtm( &
!                              ! -- forward inputs
!                              level_p, layer_p, layer_t, layer_w, layer_o,           &  ! input, k x m
!
!                              surface_temperature,                                   &  ! input, m
!                              surface_emissivity,                                    &  ! input, l*m
!                              surface_reflectivity,                                  &  ! input, l*m
!
!                              absorber,                                              &  ! input, 0:k x j x m
!
!                              tau_layer_index,                                       &  ! input, k x j x m
!                              flux_tau_layer_index,                                  &  ! input, k x j x m
!                              solar_tau_layer_index,                                 &  ! input, k x j x m
!
!                              tau_predictor,                                         &  ! input, imax x k x m
!                              flux_tau_predictor,                                    &  ! input, imax x k x m
!                              solar_tau_predictor,                                   &  ! input, imax x k x m
!
!                              tau,                                                   &  ! input, k x l*m
!                              flux_tau,                                              &  ! input, k x l*m
!                              solar_tau,                                             &  ! input, k x l*m
!
!                              layer_radiance,                                        &  ! input, k x l*m
!                              downwelling_radiance,                                  &  ! input, l*m
!                              upwelling_radiance,                                    &  ! input, l*m
!
!                              ! -- k-matrix inputs
!                              tau_k,                                                 &  ! in/output, k x l*m
!                              flux_tau_k,                                            &  ! in/output, k x l*m
!                              solar_tau_k,                                           &  ! in/output, k x l*m
!
!                              layer_radiance_k,                                      &  ! in/output, k x l*m
!                              downwelling_radiance_k,                                &  ! in/output, l*m
!                              upwelling_radiance_k,                                  &  ! in/output, l*m
!
!                              brightness_temperature_k,                              &  ! in/output, l*m
!
!                              ! -- other inputs
!                              secant_view_angle,                                     &  ! input, m
!                              secant_solar_angle,                                    &  ! input, m
!                              n_channels_per_profile,                                &  ! input, m
!                              channel_index,                                         &  ! input, l*m
!
!                              ! -- k-matrix outputs
!                              level_p_k, layer_p_k, layer_t_k, layer_w_k, layer_o_k, &  ! in/output, k x l*m
!
!                              surface_temperature_k,                                 &  ! in/output, l*m
!                              surface_emissivity_k,                                  &  ! in/output, l*m
!                              surface_reflectivity_k,                                &  ! in/output, l*m
!
!                              ! optional inputs
!                              message_log = message_log )
!
! input arguments:
!
!       level_p:                   profile set layer interface pressure array. the toa
!                                  pressure is not included. toa pressure is parameterised
!                                  in the parameters module.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_p:                   profile set layer average pressure array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_t:                   profile set layer average temperature array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_w:      .            profile set layer average water vapor mixing ratio array
!                                  units:      g/kg
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       layer_o:                   profile set layer average ozone mixing ratio array.
!                                  units:      ppmv
!                                  type:       real
!                                  dimension:  k x m
!                                  attributes: intent( in )
!
!       surface_temperature:       profile set surface temperature array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       surface_emissivity:        profile set surface emissivity array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       surface_reflectivity:      profile set surface reflectivity array
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
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
!                                  dimension:  i x k x m
!                                  attributes: intent( in )
!
!       flux_tau_predictor:        predictor profiles for the thermal flux transmittance.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  i x k x m
!                                  attributes: intent( in )
!
!       solar_tau_predictor:       predictor profiles for the solar transmittance.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  i x k x m
!                                  attributes: intent( in )
!
!       tau:                       layer->toa transmittance for the satellite
!                                  view angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in )
!
!       flux_tau:                  layer->sfc transmittance for the default
!                                  diffusivity angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in )
!
!       solar_tau:                 layer->sfc transmittance for the solar
!                                  zenith angle.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in )
!
!       layer_radiance:            layer planck radiances at every layer for
!                                  each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in )
!
!       downwelling_radiance:      toa->sfc radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       upwelling_radiance:        toa radiances for each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in )
!
!       tau_k:                     layer->toa adjoint transmittance for the satellite
!                                  view angle.
!                                  n.b.: set to zero upon output.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in out )
!
!       flux_tau_k:                layer->sfc adjoint transmittance for the default
!                                  diffusivity angle.
!                                  n.b.: set to zero upon output.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in out )
!
!       solar_tau_k:               layer->sfc adjoint transmittance for the solar
!                                  zenith angle.
!                                  n.b.: set to zero upon output.
!                                  units:      none.
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in out )
!
!       layer_radiance_k:          layer planck adjoint radiances at every layer for
!                                  each channel/profile.
!                                  n.b.: set to zero upon output.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k x l*m
!                                  attributes: intent( in out )
!
!       downwelling_radiance_k:    toa->sfc adjoint radiances for each channel/profile due
!                                  to thermal flux and solar components.
!                                  n.b.: set to zero upon output.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       upwelling_radiance_k:      toa adjoint radiances for each channel/profile.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  n.b.: set to zero upon output.
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       brightness_temperature_k:  adjoint temperatures corresponding to the
!                                  toa adjoint radiances.
!                                  n.b.: set to zero upon output.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       secant_view_angle:         secant of the satellite view angle measured
!                                  from nadir for each profile in the set.
!                                  units:      none
!                                  type:       real
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       secant_solar_angle:        secant of the solar zenith angle for each
!                                  profile in the set.
!                                  units:      none
!                                  type:       real
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       n_channels_per_profile:    the number of channels for each profile in the
!                                  set for which radiances are required.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  m
!                                  attributes: intent( in )
!
!       channel_index:             channel index id array. each element is a unique
!                                  index to a (supported) sensor channel.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  l*m
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
!       level_p_k:                 profile set layer interface pressure k-matrix
!                                  adjoint array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_p_k:                 profile set layer average pressure k-matrix
!                                  adjoint array.
!                                  units:      hpa
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_t_k:                 profile set layer average temperature k-matrix
!                                  adjoint array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_w_k:      .          profile set layer average water vapor mixing ratio
!                                  k-matrix adjoint array.
!                                  units:      g/kg
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       layer_o_k:                 profile set layer average ozone mixing ratio
!                                  k-matrix adjoint array.
!                                  units:      ppmv
!                                  type:       real
!                                  dimension:  k x l*m
!                                              nb: this is a 2-d array.
!                                  attributes: intent( in out )
!
!       surface_temperature_k:     profile set surface temperature k-matrix adjoint
!                                  array.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       surface_emissivity_k:      profile set surface emissivity k-matrix adjoint
!                                  array.
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
!                                              nb: this is a 1-d array.
!                                  attributes: intent( in out )
!
!       surface_reflectivity_k:    profile set surface reflectivity k-matrix adjoint
!                                  array.
!                                  units:      none
!                                  type:       real
!                                  dimension:  l*m
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
!       all input adjoint quantities are set to zero on output.
!
! restrictions:
!       the code has is not overloaded for scalar input so the input
!       arguments must be dimensioned accordingly, even if only one
!       profile or channel is being passed.
!
! procedure:
!       see individual module function documentation.
!s-
!--------------------------------------------------------------------------------

  function k_matrix_rtm( &

             ! -- forward inputs
             level_p, layer_p, layer_t, layer_w, layer_o,           &  ! input, k x m

             surface_temperature,                                   &  ! input, m
             surface_emissivity,                                    &  ! input, l*m
             surface_reflectivity,                                  &  ! input, l*m

             absorber,                                              &  ! input, 0:k x j x m

             tau_layer_index,                                       &  ! input, k x j x m
             flux_tau_layer_index,                                  &  ! input, k x j x m
             solar_tau_layer_index,                                 &  ! input, k x j x m

             tau_predictor,                                         &  ! input, imax x k x m
             flux_tau_predictor,                                    &  ! input, imax x k x m
             solar_tau_predictor,                                   &  ! input, imax x k x m

             tau,                                                   &  ! input, k x l*m
             flux_tau,                                              &  ! input, k x l*m
             solar_tau,                                             &  ! input, k x l*m

             layer_radiance,                                        &  ! input, k x l*m
             downwelling_radiance,                                  &  ! input, l*m
             upwelling_radiance,                                    &  ! input, l*m

             ! -- k-matrix inputs
             tau_k,                                                 &  ! in/output, k x l*m
             flux_tau_k,                                            &  ! in/output, k x l*m
             solar_tau_k,                                           &  ! in/output, k x l*m

             layer_radiance_k,                                      &  ! in/output, k x l*m
             downwelling_radiance_k,                                &  ! in/output, l*m
             upwelling_radiance_k,                                  &  ! in/output, l*m

             brightness_temperature_k,                              &  ! in/output, l*m

             ! -- other inputs
             secant_view_angle,                                     &  ! input, m
             secant_solar_angle,                                    &  ! input, m
             n_channels_per_profile,                                &  ! input, m
             channel_index,                                         &  ! input, l*m

             ! -- k-matrix outputs
             level_p_k, layer_p_k, layer_t_k, layer_w_k, layer_o_k, &  ! in/output, k x l*m

             surface_temperature_k,                                 &  ! in/output, l*m
             surface_emissivity_k,                                  &  ! in/output, l*m
             surface_reflectivity_k,                                &  ! in/output, l*m

             ! optional inputs
             message_log )                                          &

           result ( error_status )



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
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

    ! -- k-matrix inputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: tau_k                      ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: flux_tau_k                 ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: solar_tau_k                ! k x l*m

    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_radiance_k           ! k x l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: downwelling_radiance_k     ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: upwelling_radiance_k       ! l*m

    real( fp_kind ), dimension( : ),        intent( in out ) :: brightness_temperature_k   ! l*m

    ! -- other inputs
    real( fp_kind ), dimension( : ),        intent( in )     :: secant_view_angle          ! m
    real( fp_kind ), dimension( : ),        intent( in )     :: secant_solar_angle         ! m
    integer,         dimension( : ),        intent( in )     :: n_channels_per_profile     ! m
    integer,         dimension( : ),        intent( in )     :: channel_index              ! l*m

    ! -- k-matrix outputs
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: level_p_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_p_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_t_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_w_k                  ! k x l*m
    real( fp_kind ), dimension( :, : ),     intent( in out ) :: layer_o_k                  ! k x l*m

    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_temperature_k      ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_emissivity_k       ! l*m
    real( fp_kind ), dimension( : ),        intent( in out ) :: surface_reflectivity_k     ! l*m

    ! -- optional input
    character( * ), optional,               intent( in )     :: message_log
    
    
    ! ---------------
    ! function result
    ! ---------------

    integer :: error_status


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'k_matrix_rtm'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- scalars
    character( 100 ) :: message
    character( 5 )   :: value_in, value_allowed

    integer :: m, n_profiles          ! profile loop variables
    integer :: l, l1, l2, n_channels  ! channel loop/index variables
    integer :: k, n_layers            ! layer loop variables

    integer :: valid_solar

    ! -- maximum channels pseudo parameter
    integer :: max_n_channels
    logical :: is_set

    ! -- arrays for integrated absorber amounts, 0:k x j
    real( fp_kind ), dimension( 0:size( absorber, dim = 1 )-1, &
                                  size( absorber, dim = 2 )    ) :: tau_absorber,      &
                                                                    flux_tau_absorber, &
                                                                    solar_tau_absorber

    ! -- arrays for k-matrix adjoints of integrated absorber amounts, 0:k x j
    real( fp_kind ), dimension( 0:size( absorber, dim = 1 )-1, &
                                  size( absorber, dim = 2 )    ) :: absorber_k, &
                                                                    tau_absorber_k,      &
                                                                    flux_tau_absorber_k, &
                                                                    solar_tau_absorber_k

    ! -- arrays for k-matrix adjoint predictors, imax x k
    real( fp_kind ), dimension( size( tau_predictor, dim = 1 ), &
                                size( tau_predictor, dim = 2 )  ) :: tau_predictor_k,      &
                                                                     flux_tau_predictor_k, &
                                                                     solar_tau_predictor_k

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
    !#                   -- determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! ------------------
    ! get the dimensions
    ! ------------------

    n_layers   = size( layer_p, dim = 1 )
    n_profiles = size( layer_p, dim = 2 )
    n_channels = maxval( n_channels_per_profile ) ! result is scalar for rank-1 array


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


    ! ------------------------------------------------------
    ! check that the number of profiles is not greater than
    ! max_n_profiles. this is simply a limit to restrict the
    ! size of the input arrays so they're not too big.
    ! ------------------------------------------------------

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


    ! -----------------------------------------------------
    ! check that the number of channels is not greater than
    ! than the number of channels with which the model was
    ! initialised
    ! -----------------------------------------------------

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
    if ( any( surface_temperature  < zero ) .or. &
         any( surface_emissivity   < zero ) .or. &
         any( surface_reflectivity < zero )      ) then
      error_status = failure
      call display_message( routine_name, &
                            'negative values found in surface properties (tsfc,esfc,rsfc).', &
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


      ! -------------------------------------------
      ! determine the end channel index index range
      ! -------------------------------------------

      l2 = l1 + n_channels_per_profile( m ) - 1


      ! ----------------------------------------------------
      ! modify absorber quantities by the angle secant
      ! could put a loop here but here's hoping the compiler
      ! recognises this as a group of loops over layer.
      ! ----------------------------------------------------

      ! -- upwelling transmittance
      tau_absorber( 0:, : ) = secant_view_angle( m ) * absorber( 0:, :, m )
      tau_absorber( 1:, 2 ) = tau_absorber( 1:, 2 ) - toa_pressure

      ! -- flux transmittance
      if ( any( is_microwave_channel( channel_index( l1:l2 ) ) == 0 ) ) then
        flux_tau_absorber( 0:, : ) = secant_diffusivity_angle * absorber( 0:, :, m )
        flux_tau_absorber( 0:, : ) = secant_view_angle( m ) * absorber( 0:, :, m )
        flux_tau_absorber( 1:, 2 ) = flux_tau_absorber( 1:, 2 ) - toa_pressure
      end if

      ! -- solar transmittance
      if ( ( any( is_solar_channel( channel_index( l1:l2 ) ) == 1 ) ) .and. &
           secant_solar_angle( m ) < max_secant_solar_angle ) then
        solar_tau_absorber( 0:, : ) = secant_solar_angle( m ) * absorber( 0:, :, m )
        solar_tau_absorber( 1:, 2 ) = solar_tau_absorber( 1:, 2 ) - toa_pressure
      end if


      !#------------------------------------------------------------------------#
      !#                           -- channel loop --                           #
      !#------------------------------------------------------------------------#

      l_channel_loop: do l = l1, l2


        ! ---------------------------------------------------------
        ! initialise local, channel dependent, adjoint variables
        !
        ! this initialisation may slow things down a bit if 
        ! a) the compiler isn't too clever about optimising and
        !    generates code that loops over each array separately
        ! b) there are a large number of channels.
        !
        ! the solution to (a) above is put in explicit loops, which
        ! may be necessary.
        !
        ! the alternative to (b) above is to define the arrays
        ! with a channel dimension which could be problematical
        ! memorywise for lots of channels.
        !
        ! these will eventually be moved outside both the channel
        ! and profile loop as once they are initialised, they are
        ! set to zero in the relevant adjoint routine on output.
        !----------------------------------------------------------

        ! -- absorber arrays, 0:k x j
        absorber_k( 0:, : )           = zero
        tau_absorber_k( 0:, : )       = zero
        flux_tau_absorber_k( 0:, : )  = zero
        solar_tau_absorber_k( 0:, : ) = zero


        ! -- predictor arrays, imax x k
        tau_predictor_k( :, : )       = zero
        flux_tau_predictor_k( :, : )  = zero
        solar_tau_predictor_k( :, : ) = zero


        ! ----------------------------------------------------
        ! set the "this is a channel influenced by solar" flag
        ! ----------------------------------------------------

        valid_solar = 0

        if ( is_solar_channel( channel_index( l ) ) == 1 .and. &
             secant_solar_angle( m ) < max_secant_solar_angle        ) valid_solar = 1



        ! ------------------------------------------------
        ! calculate the adjoint of the current channel toa
        ! radiance or brightness temperature
        ! ------------------------------------------------

        call compute_radiance_ad( &
                                  ! -- forward input
                                  layer_t( :, m ),               &  ! input, k

                                  surface_temperature( m ),      &  ! input, scalar
                                  surface_emissivity( l ),       &  ! input, scalar
                                  surface_reflectivity( l ),     &  ! input, scalar

                                  tau(      :, l ),              &  ! input, k
                                  flux_tau( :, l ),              &  ! input, k
                                  solar_tau( n_layers, l ),      &  ! input, scalar

                                  layer_radiance( :, l ),        &  ! input, k
                                  downwelling_radiance( l ),     &  ! input, scalar
                                  upwelling_radiance( l ),       &  ! input, scalar

                                  ! -- k-matrix input
                                  layer_radiance_k( :, l ),      &  ! in/output, k
                                  downwelling_radiance_k( l ),   &  ! in/output, scalar
                                  upwelling_radiance_k( l ),     &  ! in/output, scalar

                                  brightness_temperature_k( l ), &  ! in/output, scalar

                                  ! -- other input
                                  secant_solar_angle( m ),       &  ! input, scalar
                                  valid_solar,                   &  ! input, scalar
                                  channel_index( l ),            &  ! input, scalar

                                  ! -- k-matrix output
                                  layer_t_k( :, l ),             &  ! in/output, k

                                  surface_temperature_k( l ),    &  ! in/output, scalar
                                  surface_emissivity_k( l ),     &  ! in/output, scalar
                                  surface_reflectivity_k( l ),   &  ! in/output, scalar

                                  tau_k( :, l ),                 &  ! in/output, k
                                  flux_tau_k( :, l ),            &  ! in/output, k
                                  solar_tau_k( n_layers, l )     )  ! in/output, scalar



        !#----------------------------------------------------------------------#
        !# --    calculate the k-matrix adjoint result for the solar term    -- #
        !#----------------------------------------------------------------------#

        ! ----------------------------------------------------
        ! if the current channel is a solar sensitive channel,
        !   and
        ! the solar angle is valid, then calculate the adjoint
        ! of the transmittance for direct solar.
        ! ----------------------------------------------------

        solar_if_block: if ( valid_solar == 1 ) then


          ! -----------------------------------------
          ! adjoint of the direct solar transmittance
          ! -----------------------------------------

          call compute_transmittance_ad( &
                                         ! -- forward input
                                         solar_tau_absorber( 0:, : ),      &   ! input, 0:k x j
                                         solar_tau_predictor( :, :, m ),   &   ! input, i x k
                                         solar_tau( :, l ),                &   ! input, k

                                         ! -- k-matrix input
                                         solar_tau_k( :, l ),              &   ! in/output, k

                                         ! -- other input
                                         solar_tau_layer_index( :, :, m ), &   ! input, k x j
                                         channel_index( l ),               &   ! input, scalar
                                         down,                             &   ! input, scalar

                                         ! -- k-matrix output
                                         solar_tau_absorber_k( 0:, : ),    &   ! in/output, 0:k x j
                                         solar_tau_predictor_k( :, : )     )   ! in/output, i x k


          ! -----------------------------------------------
          ! k-matrix adjoint of the standard predictor copy
          ! -----------------------------------------------

          tau_predictor_k( 1:max_n_standard_predictors, : ) = &
                  tau_predictor_k( 1:max_n_standard_predictors, : ) + &
            solar_tau_predictor_k( 1:max_n_standard_predictors, : )

          solar_tau_predictor_k( 1:max_n_standard_predictors, : ) = zero


          ! ----------------------------------------------
          ! k-matrix adjoint of the integrateed predictors
          ! ----------------------------------------------

          call compute_predictors_ad( &
                                      ! -- forward input
                                      layer_p( :, m ),               &  ! input,  k
                                      layer_t( :, m ),               &  ! input,  k
                                      layer_w( :, m ),               &  ! input,  k
                                      solar_tau_absorber( 0:, : ),   &  ! input,  0:k x j

                                      ! -- k-matrix input
                                      solar_tau_predictor_k( :, : ), &  ! in/output, i x k

                                      ! -- k-matrix output
                                      layer_p_k( :, l ),             &  ! in/output,  k
                                      layer_t_k( :, l ),             &  ! in/output,  k
                                      layer_w_k( :, l ),             &  ! in/output,  k
                                      solar_tau_absorber_k( 0:, : ), &  ! in/output,  0:k x j

                                      no_standard = 1                )  ! optional input


          ! ----------------------------------------------------------
          ! k-matrix adjoint of the nadir absorber amount modification
          ! ----------------------------------------------------------

          absorber_k( 0:, : ) = absorber_k( 0:, : ) + &
                                  ( secant_solar_angle( m ) * solar_tau_absorber_k( 0:, : ) )

        end if solar_if_block



        !#----------------------------------------------------------------------#
        !# --    calculate the k-matrix adjoint result for the flux term     -- #
        !#----------------------------------------------------------------------#

        ! --------------------------------------------------
        ! if the current channel is an infrared channel,
        ! then calculate the adjoint of the downwelling flux
        ! transmittance using the flux absorber amounts.
        !
        ! if the current channel is a microwave channel,
        ! then calculate the adjoint of the flux transmittance
        ! using the upwelling absorber amounts.
        ! --------------------------------------------------

        flux_if_block: if ( is_microwave_channel( channel_index( l ) ) == 0 ) then


          ! ------------------------------------
          ! adjoint of the ir flux transmittance
          ! ------------------------------------

          call compute_transmittance_ad( &
                                         ! -- forward input
                                         flux_tau_absorber( 0:, : ),      &   ! input, 0:k x j
                                         flux_tau_predictor( :, :, m ),   &   ! input, i x k
                                         flux_tau( :, l ),                &   ! input, k

                                         ! -- k-matrix input
                                         flux_tau_k( :, l ),              &   ! in/output, k

                                         ! -- other input
                                         flux_tau_layer_index( :, :, m ), &   ! input, k x j
                                         channel_index( l ),              &   ! input, scalar
                                         down,                            &   ! input, scalar

                                         ! -- k-matrix output
                                         flux_tau_absorber_k( 0:, : ),    &   ! in/output, 0:k x j
                                         flux_tau_predictor_k( :, : )     )   ! in/output, i x k


          ! -----------------------------------------------
          ! k-matrix adjoint of the standard predictor copy
          ! -----------------------------------------------

          tau_predictor_k( 1:max_n_standard_predictors, : ) = &
                 tau_predictor_k( 1:max_n_standard_predictors, : ) + &
            flux_tau_predictor_k( 1:max_n_standard_predictors, : )

          flux_tau_predictor_k( 1:max_n_standard_predictors, : ) = zero


          ! ----------------------------------------------
          ! k-matrix adjoint of the integrateed predictors
          ! ----------------------------------------------

          call compute_predictors_ad( &
                                      ! -- forward input
                                      layer_p( :, m ),              &  ! input,  k
                                      layer_t( :, m ),              &  ! input,  k
                                      layer_w( :, m ),              &  ! input,  k
                                      flux_tau_absorber( 0:, : ),   &  ! input,  0:k x j

                                      ! -- k-matrix input
                                      flux_tau_predictor_k( :, : ), &  ! in/output, i x k

                                      ! -- k-matrix output
                                      layer_p_k( :, l ),            &  ! in/output,  k
                                      layer_t_k( :, l ),            &  ! in/output,  k
                                      layer_w_k( :, l ),            &  ! in/output,  k
                                      flux_tau_absorber_k( 0:, : ), &  ! in/output,  0:k x j

                                      no_standard = 1               )  ! optional input


          ! ----------------------------------------------------------
          ! k-matrix adjoint of the nadir absorber amount modification
          ! ----------------------------------------------------------

          absorber_k( 0:, : ) = absorber_k( 0:, : ) + &
                                   ( secant_diffusivity_angle * flux_tau_absorber_k( 0:, : ) )


        else  ! we have a microwave channel....

          ! ----------------------------------------------
          ! if total transmittance /= 0, calculate adjoint
          ! flux transmittance for the microwave
          ! ----------------------------------------------
        
          if ( tau( n_layers, l ) > tolerance ) then
            do k = n_layers, 2, -1
              flux_tau_k( 1, l ) = flux_tau_k( 1, l ) + ( flux_tau_k( k, l ) / &
              !                                           ------------------
                                                             tau( k-1, l )   )

              tau_k( k-1, l ) = tau_k( k-1, l ) - ( flux_tau_k( k, l ) * flux_tau( 1, l ) / &
              !                                     -------------------------------------
                                                                tau( k-1, l )**2          )
              flux_tau_k( k, l ) = zero
            end do
            tau_k( n_layers, l ) = tau_k( n_layers, l ) + flux_tau_k( 1, l )
            flux_tau_k( 1, l ) = zero
          else
            flux_tau_k( :, l ) = zero
          end if


        end if flux_if_block



        !#----------------------------------------------------------------------#
        !# --   calculate the k-matrix adjoint result for the transmittance  -- #
        !#----------------------------------------------------------------------#

        ! -----------------------------------------------------
        ! calculate the adjoint of the upwelling transmittances
        ! for the satellite view angle
        ! -----------------------------------------------------

        call compute_transmittance_ad( &
                                       ! -- forward input
                                       tau_absorber( 0:, : ),      &   ! input, 0:k x j
                                       tau_predictor( :, :, m ),   &   ! input, i x k
                                       tau( :, l ),                &   ! input, k

                                       ! -- k-matrix input
                                       tau_k( :, l ),              &   ! in/output, k

                                       ! -- other input
                                       tau_layer_index( :, :, m ), &   ! input, k x j
                                       channel_index( l ),         &   ! input, scalar
                                       up,                         &   ! input, scalar

                                       ! -- k-matrix output
                                       tau_absorber_k( 0:, : ),    &   ! in/output, 0:k x j
                                       tau_predictor_k( :, : )     )   ! in/output, i x k



        ! --------------------------------------
        ! k-matrix adjoint of all the predictors
        ! --------------------------------------

        call compute_predictors_ad( &
                                    ! -- forward input
                                    layer_p( :, m ),         &  ! input,  k
                                    layer_t( :, m ),         &  ! input,  k
                                    layer_w( :, m ),         &  ! input,  k
                                    tau_absorber( 0:, : ),   &  ! input,  0:k x j

                                    ! -- k-matrix input
                                    tau_predictor_k( :, : ), &  ! in/output, i x k

                                    ! -- k-matrix output
                                    layer_p_k( :, l ),       &  ! in/output,  k
                                    layer_t_k( :, l ),       &  ! in/output,  k
                                    layer_w_k( :, l ),       &  ! in/output,  k
                                    tau_absorber_k( 0:, : )  )  ! in/output,  0:k x j


        ! ----------------------------------------------------------
        ! k-matrix adjoint of the nadir absorber amount modification
        ! ----------------------------------------------------------

        absorber_k( 0:, : ) = absorber_k( 0:, : ) + &
                                 ( secant_view_angle( m ) * tau_absorber_k( 0:, : ) )



        !#----------------------------------------------------------------------#
        !#            -- calculate the absorber k-matrix adjoints --            #
        !#----------------------------------------------------------------------#

        call compute_absorber_amount_ad( &
                                         ! -- forward input
                                         level_p( :, m ),     &  ! input,  k
                                         layer_w( :, m ),     &  ! input,  k
                                         layer_o( :, m ),     &  ! input,  k

                                         ! -- k-matrix input
                                         absorber_k( 0:, : ), &  ! in/output, 0:k x j

                                         ! -- k-matrix output
                                         level_p_k( :, l ),   &  ! in/ouput,  k
                                         layer_w_k( :, l ),   &  ! in/ouput,  k
                                         layer_o_k( :, l )    )  ! in/ouput,  k


      end do l_channel_loop


      ! ------------------------------
      ! update the channel begin index
      ! ------------------------------

      l1 = l2 + 1


    end do m_profile_loop



    !#--------------------------------------------------------------------------#
    !#                              -- done --                                  #
    !#--------------------------------------------------------------------------#

    error_status = success

  end function k_matrix_rtm

end module k_matrix_model


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
! revision 1.4  2001/09/04 21:29:11  paulv
! - updated documentation.
!
! revision 1.3  2001/08/31 21:28:52  paulv
! - removed input data checks from compute_rtm_k. the same checks are performed
!   in the main routine, k_matrix_rtm, so there was no need to replicate them.
! - added check for negative profile and surface data in k_matrix_rtm.
! - maximum solar angle secant is no longer calculated in k_matrix_rtm but
!   is declared as a parameter in the parameters module.
!
! revision 1.2  2001/08/16 16:38:57  paulv
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
! revision 1.1  2001/08/01 16:34:02  paulv
! initial checkin.
!
!
!
!
