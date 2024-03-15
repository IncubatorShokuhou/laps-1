!------------------------------------------------------------------------------
!m+
! name:
!       radiance
!
! purpose:
!       rt model radiance module
!
! category:
!       ncep rtm
!
! calling sequence:
!       use radiance
!
! outputs:
!       none.
!
! modules:
!       type_kinds:              module containing data type kind definitions.
!
!       parameters:              module containing parameter definitions for
!                                the rt model.
!
!       spectral_coefficients:   module containing the rt model spectral
!                                coefficients.
!
!       sensor_planck_routines:  module containing all the forward, tangent-
!                                linear, and adjoint planck radiance and
!                                temperature subroutines. 
!
! contains:
!       compute_radiance:        public subroutine to calculate the channel toa
!                                radiance and brightness temperature.
!
!       compute_radiance_tl:     public subroutine to calculate the tangent-
!                                linear toa radiance and brightness temperature.
!
!       compute_radiance_ad:     public subroutine to calculate the adjoint of
!                                the toa radiance and brightness temperature.
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
!       written by:     paul van delst, cimss/ssec 11-jul-2000
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

module radiance


  ! ---------------------
  ! module use statements
  ! ---------------------

  use type_kinds, only : fp_kind
  use parameters
  use spectral_coefficients
  use sensor_planck_routines


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private
  public  :: compute_radiance
  public  :: compute_radiance_tl
  public  :: compute_radiance_ad


contains


!--------------------------------------------------------------------------------
!s+
! name:
!       compute_radiance
!
! purpose:
!       public subroutine to calculate the toa radiance and brightness temperature.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_radiance( temperature,           &  ! input, k
!
!                              surface_temperature,   &  ! input, scalar
!                              surface_emissivity,    &  ! input, scalar
!                              surface_reflectivity,  &  ! input, scalar
!
!                              tau,                   &  ! input, k
!                              flux_tau,              &  ! input, k
!                              solar_tau,             &  ! input, scalar
!
!                              secant_solar_angle,    &  ! input, scalar
!                              valid_solar,           &  ! input, scalar
!                              channel_index,         &  ! input, scalar
!
!                              layer_radiance,        &  ! output, k
!                              downwelling_radiance,  &  ! output, scalar
!                              upwelling_radiance,    &  ! output, scalar
!
!                              brightness_temperature )  ! output, scalar
!
!
! input arguments:
!       temperature:             profile layer average temperature array.
!                                units:      kelvin
!                                type:       float
!                                dimension:  k
!                                attributes: intent( in )
!
!       surface_temperature:     surface boundary temperature.
!                                units:      kelvin
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( in )
!
!       surface_emissivity:      surface boundary emissivity
!                                units:      none
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( in )
!
!       surface_reflectivity:    surface boundary reflectivity
!                                units:      none
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( in )
!
!       tau:                     layer-to-space transmittance profile for
!                                a particular satellite view angle.
!                                units:      none
!                                type:       real
!                                dimension:  k
!                                attributes: intent( in )
!
!       flux_tau:                layer-to-surface transmittance profile for
!                                either the diffuse approximation angle (ir)
!                                or the satellite view angle (mw). the latter
!                                assumes specular reflectivity.
!                                units:      none
!                                type:       real
!                                dimension:  k
!                                attributes: intent( in )
!
!       solar_tau:               total space-to-surface transmittance at the
!                                solar zenith angle.
!                                units:      none
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( in )
!
!       secant_solar_angle:      secant of the solar zenith angle corresponding
!                                to that used in calculating the total solar
!                                transmittance.
!                                units:      none
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( in )
!
!       valid_solar:             flag indicating if the solar component should
!                                be included.
!                                if = 0, no solar (if sensor channel frequency
!                                        is less than a preset cutoff or if solar
!                                        zenith angle is greater than its preset
!                                         cutoff.)
!                                   = 1, include solar
!                                units:      none.
!                                type:       integer
!                                dimension:  scalar
!                                attributes: intent( in )
!
!       channel_index:           channel index id. this is a unique index
!                                to a (supported) sensor channel.
!                                units:      none
!                                type:       integer
!                                dimension:  scalar
!                                attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       layer_radiance:          channel planck radiance for every input layer.
!                                units:      mw/(m^2.sr.cm^-1)
!                                type:       real
!                                dimension:  k.
!                                attributes: intent( out )
!
!       downwelling_radiance:    channel radiance at surface due to downwelling
!                                flux and solar components.
!                                units:      mw/(m^2.sr.cm^-1)
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( out )
!
!       upwelling_radiance:      channel toa radiance simulating the satellite
!                                sensor measurement. this is composed of the
!                                reflected downwelling propagated through the
!                                atmosphere as well as the upwelling only component.
!                                units:      mw/(m^2.sr.cm^-1)
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( out )
!
!       brightness_temperature:  channel toa brightness temperature.
!                                units:      kelvin
!                                type:       real
!                                dimension:  scalar
!                                attributes: intent( out )
!
! optional output arguments:
!       none.
!
! calls:
!       sensor_planck_radiance:    function to compute the planck radiance
!                                  for a specified channel given the temperature.
!                                  source: sensor_planck_routines module
!
!       sensor_planck_temperature: function to compute the planck temperature
!                                  for a specified channel given the radiance.
!                                  source: sensor_planck_routines module
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
!       the downwelling radiance is first initialised to the space emissisio
!       boundary term using precalculated cosmic background radiances,
!
!         r_down = cbr * flux_tau(1)
!
!       where the emissivity of space is implicitly assumed to be 1.0 and
!       flux_tau(1) is the space-to-ground transmittance.
!
!       the contributions of all the layers except the surface layer to the
!       downwelling flux is accumulated,
!
!                            __k-1
!                           \
!         r_down = r_down +  >  b(k) * dflux_tau(k)
!                           /__
!                              k=1
!
!       the surface layer contribution is then added explicitly,
!
!         r_down = r_down + ( b(k) * ( 1 - flux_tau(k) ) )
!
!       to avoid exceeding the arrays bounds of flux_tau 
!       (i.e. flux_tau(k+1) == 1.0 ) or requiring flux_tau to be
!       dimensioned 0:k.
!
!       the solar term is then added if required,
!
!         r_down = r_down + ( solar_irradiance * solar_tau * cos(solar_theta) )
!
!       the downwelling radiance is then reflected off the surface, added
!       to the surface emission term, propagated upwards through the atmosphere
!       and used to initialise the upwelling radiance term,
!
!         r_up = ( ( e_sfc * b_sfc ) + ( r_sfc * r_down ) ) * tau(k)
!
!       the contributions of all the layers except the top layer to the
!       upwelling radiance is accumulated,
!
!                        __ 2
!                       \
!         r_up = r_up +  >  b(k) * dtau(k)
!                       /__
!                          k=k
!
!       the top layer contribution is then added explicitly,
!
!         r_up = r_up + ( b(1) * ( 1 - tau(1) ) )
!
!       to avoid exceeding the arrays bounds of tau (i.e. tau(0) == 1.0 )
!       or requiring tau to be dimensioned 0:k.
!
!       the final upwelling radiance is then converted to a brightness
!       temperature.
!
!s-      
!--------------------------------------------------------------------------------


  subroutine compute_radiance( temperature,           &  ! input, k      

                               surface_temperature,   &  ! input, scalar 
                               surface_emissivity,    &  ! input, scalar 
                               surface_reflectivity,  &  ! input, scalar 

                               tau,                   &  ! input, k      
                               flux_tau,              &  ! input, k      
                               solar_tau,             &  ! input, scalar 

                               secant_solar_angle,    &  ! input, scalar 
                               valid_solar,           &  ! input, scalar 
                               channel_index,         &  ! input, scalar 

                               layer_radiance,        &  ! output, k     
                               downwelling_radiance,  &  ! output, scalar
                               upwelling_radiance,    &  ! output, scalar

                               brightness_temperature )  ! output, scalar



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ), intent( in )  :: temperature

    real( fp_kind ),                 intent( in )  :: surface_temperature
    real( fp_kind ),                 intent( in )  :: surface_emissivity
    real( fp_kind ),                 intent( in )  :: surface_reflectivity

    real( fp_kind ), dimension( : ), intent( in )  :: tau
    real( fp_kind ), dimension( : ), intent( in )  :: flux_tau
    real( fp_kind ),                 intent( in )  :: solar_tau

    real( fp_kind ),                 intent( in )  :: secant_solar_angle
    integer,                         intent( in )  :: valid_solar
    integer,                         intent( in )  :: channel_index

    real( fp_kind ), dimension( : ), intent( out ) :: layer_radiance
    real( fp_kind ),                 intent( out ) :: downwelling_radiance
    real( fp_kind ),                 intent( out ) :: upwelling_radiance

    real( fp_kind ),                 intent( out ) :: brightness_temperature

 

    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_radiance'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: k, n_layers
    integer :: l

    real( fp_kind ) :: surface_b


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic size



    !#--------------------------------------------------------------------------#
    !#                   -- determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    n_layers = size( temperature )



    !#--------------------------------------------------------------------------#
    !#              -- calculate the downwelling thermal flux --                #
    !#--------------------------------------------------------------------------#

    ! -- assign the channel index to a short name
    l = channel_index


    ! --------------------------------------------
    ! initialise the downwelling radiance to the
    ! space emission boundary term reaching the
    ! surface. thhe cosmic background radiance is
    ! zero for infrared channels and precalculated
    ! for microwave channels. the emissivity of
    ! space is assumed to be 1.0.
    !
    ! cosmic background data from the
    ! spectral_coefficients module
    ! --------------------------------------------

    downwelling_radiance = cosmic_background_radiance( l ) * flux_tau( 1 )


    ! --------------------------------
    ! loop over layers from toa->sfc-1
    ! --------------------------------

    k_down_layer_loop: do k = 1, n_layers - 1

      ! -- calculate the planck layer radiance
      call sensor_planck_radiance( l,                  &  ! input
                                   temperature( k ),   &  ! input
                                   layer_radiance( k ) )  ! output

      ! -- accumulate absorption and emission for current layer.
      ! -- lte assumed.
      downwelling_radiance = downwelling_radiance + &
                             ( layer_radiance( k ) * ( flux_tau( k+1 ) - flux_tau( k ) ) )

    end do k_down_layer_loop


    ! ----------------------------------
    ! flux bottom layer (closest to sfc)
    ! ----------------------------------

    ! -- lowest layer planck radiance
    call sensor_planck_radiance( l,                         &  ! input
                                 temperature( n_layers ),   &  ! input
                                 layer_radiance( n_layers ) )  ! output

    ! -- contribution of lowest layer. note that at the
    ! -- surface, a transmittance of 1.0 is used.
    downwelling_radiance = downwelling_radiance + &
                           ( layer_radiance( n_layers ) * ( one - flux_tau( n_layers ) ) )



    !#--------------------------------------------------------------------------#
    !#                -- calculate the downwelling solar terms --               #
    !#--------------------------------------------------------------------------#

    solar_term: if ( valid_solar == 1 ) then

      downwelling_radiance = downwelling_radiance + &

                             ( solar_irradiance( l ) * solar_tau / &
      !                        ---------------------------------
                                      secant_solar_angle         )

    end if solar_term



    !#--------------------------------------------------------------------------#
    !#   -- reflect the downwelling radiance and add the surface emission --    #
    !#   -- and use it to initialise the upwelling radiance               --    #
    !#--------------------------------------------------------------------------#

    ! -- calculate the surface term
    call sensor_planck_radiance( l,                   &  ! input
                                 surface_temperature, &  ! input
                                 surface_b            )  ! output

    ! -- initialise upwelling radiance
    upwelling_radiance = ( ( surface_emissivity   * surface_b            ) + &
                           ( surface_reflectivity * downwelling_radiance )   ) * tau( n_layers )



    !#--------------------------------------------------------------------------#
    !#                  -- calculate the upwelling radiance --                  #
    !#--------------------------------------------------------------------------#

    ! --------------------------------
    ! loop over layers from sfc->toa-1
    ! --------------------------------

    k_up_layer_loop: do k = n_layers, 2, -1

      upwelling_radiance = upwelling_radiance + &
                           ( layer_radiance( k ) * ( tau( k-1 ) - tau( k ) ) )

    end do k_up_layer_loop


    ! --------------------------
    ! top layer (closest to toa)
    ! --------------------------

    upwelling_radiance = upwelling_radiance + &
                         ( layer_radiance( 1 ) * ( one - tau( 1 ) ) )



    !#--------------------------------------------------------------------------#
    !#           -- convert the radiances to brightness temperatures --         #
    !#--------------------------------------------------------------------------#

    call sensor_planck_temperature( l,                     &  ! input
                                    upwelling_radiance,    &  ! input
                                    brightness_temperature )  ! output

  end subroutine compute_radiance



!--------------------------------------------------------------------------------
!s+
! name:
!       compute_radiance_tl
!
! purpose:
!       public subroutine to calculate the tangent-linear toa radiance and
!       brightness temperature.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_radiance_tl( &
!
!                                 ! -- forward inputs
!                                 temperature,              &  ! input, k
!
!                                 surface_temperature,      &  ! input, scalar
!                                 surface_emissivity,       &  ! input, scalar
!                                 surface_reflectivity,     &  ! input, scalar
!
!                                 tau,                      &  ! input, k
!                                 flux_tau,                 &  ! input, k
!                                 solar_tau,                &  ! input, scalar
!
!                                 layer_radiance,           &  ! input, k
!                                 downwelling_radiance,     &  ! input, scalar
!                                 upwelling_radiance,       &  ! input, scalar
!
!                                 ! -- tangent-linear inputs
!                                 temperature_tl,           &  ! input, k
!
!                                 surface_temperature_tl,   &  ! input, scalar
!                                 surface_emissivity_tl,    &  ! input, scalar
!                                 surface_reflectivity_tl,  &  ! input, scalar
!
!                                 tau_tl,                   &  ! input, k
!                                 flux_tau_tl,              &  ! input, k
!                                 solar_tau_tl,             &  ! input, scalar
!
!                                 ! -- other inputs
!                                 secant_solar_angle,       &  ! input, scalar
!                                 valid_solar,              &  ! input, scalar
!                                 channel_index,            &  ! input, scalar
!
!                                 ! -- tangent-linear outputs
!                                 layer_radiance_tl,        &  ! output, k
!                                 downwelling_radiance_tl,  &  ! output, scalar
!                                 upwelling_radiance_tl,    &  ! output, scalar
!
!                                 brightness_temperature_tl )  ! output, scalar
!
!
! input arguments:
!       temperature:               profile layer average temperature array.
!                                  units:      kelvin
!                                  type:       float
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       surface_temperature:       surface boundary temperature.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_emissivity:        surface boundary emissivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_reflectivity:      surface boundary reflectivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       tau:                       layer-to-space transmittance profile for
!                                  a particular satellite view angle.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k.
!                                  attributes: intent( in )
!
!       flux_tau:                  layer-to-surface transmittance profile for
!                                  either the diffuse approximation angle (ir)
!                                  or the satellite view angle (mw). the latter
!                                  assumes specular reflectivity.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k.
!                                  attributes: intent( in )
!
!       solar_tau:                 total space-to-surface transmittance at the
!                                  solar zenith angle.
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       layer_radiance:            channel planck radiance for every layer.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k.
!                                  attributes: intent( in )
!
!       downwelling_radiance:      channel radiance at surface due to downwelling
!                                  flux and solar components.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       upwelling_radiance:        channel toa radiance simulating the satellite
!                                  sensor measurement. this is composed of the
!                                  reflected downwelling propagated through the
!                                  atmosphere as well as the upwelling only component.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       temperature_tl:            tangent-linear temperature profile.
!                                  units:      kelvin
!                                  type:       float
!                                  dimension:  k.
!                                  attributes: intent( in )
!
!       surface_temperature_tl:    tangent-linear surface boundary temperature.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_emissivity_tl:     tangent-linear surface boundary emissivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_reflectivity_tl:   tangent-linear surface boundary reflectivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       tau_tl:                    tangent-linear layer-to-space transmittance
!                                  profile.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k.
!                                  attributes: intent( in )
!
!       flux_tau_tl:               tangent-linear layer-to-surface flux transmittance
!                                  profile.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k.
!                                  attributes: intent( in )
!
!       solar_tau_tl:              tangent-linear total space-to-surface solar
!                                  transmittance.
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       secant_solar_angle:        secant of the solar zenith angle corresponding
!                                  to that used in calculating the total solar
!                                  transmittance.
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       valid_solar:               flag indicating if the solar component should
!                                  be included.
!                                  if = 0, no solar (if sensor channel frequency
!                                          is less than a preset cutoff or if solar
!                                          zenith angle is greater than its preset
!                                           cutoff.)
!                                     = 1, include solar
!                                  units:      none.
!                                  type:       integer
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       channel_index:             channel index id. this is a unique index
!                                  to a (supported) sensor channel.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       layer_radiance_tl:         tangent-linear channel planck radiance for
!                                  every layer.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k.
!                                  attributes: intent( out )
!
!       downwelling_radiance_tl:   tangent-linear channel radiance at surface
!                                  due to downwelling flux and solar components.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( out )
!
!       upwelling_radiance_tl:     tangent-linear channel toa radiance.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( out )
!
!       brightness_temperature_tl: tangent-linear channel toa brightness
!                                  temperature.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( out )
!
! optional output arguments:
!       none.
!
! calls:
!       sensor_planck_radiance_tl:    function to compute the tangent-linear
!                                     planck radiance.
!                                     source: sensor_planck_routines module
!
!       sensor_planck_temperature_tl: function to compute the tangent-linear
!                                     planck temperature.
!                                     source: sensor_planck_routines module
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
!       the downwelling radiance is first initialised to the space emission
!       boundary term using precalculated cosmic background radiances,
!
!         r_down_tl = cbr * flux_tau_tl(1)
!
!       and the emissivity of space is implicitly assumed to be 1.0 and
!       flux_tau_tl(1) is the tangent-linear form of the space-to-ground
!       transmittance.
!
!       the contributions of all the layers except the surface layer to the
!       downwelling flux is accumulated,
!
!                                  __k-1
!                                 \
!         r_down_tl = r_down_tl +  >  ( b(k) * dflux_tau_tl(k) ) + ( b_tl(k) * dflux_tau(k) )
!                                 /__
!                                    k=1
!
!       the surface layer contribution is then added explicitly,
!
!         r_down_tl = r_down_tl + ( b(k) * ( -flux_tau_tl(k) ) ) + ( b_tl(k) * ( 1 - flux_tau(k) ) )
!
!       to avoid exceeding the arrays bounds of flux_tau and flux_tau_tl
!       (i.e. flux_tau(k+1) == 1.0 ) or requiring them to be
!       dimensioned 0:k.
!
!       the solar term is then added if required,
!
!         r_down_tl = r_down_tl + ( solar_irradiance * solar_tau_tl * cos(solar_theta) )
!
!       the downwelling radiance is then reflected off the surface, added
!       to the surface emission term, propagated upwards through the atmosphere
!       and used to initialise the upwelling radiance term,
!
!         r_up_tl = ( e_sfc    * b_sfc_tl  * tau(k) ) + &
!                   ( e_sfc_tl * b_sfc     * tau(k) ) + &
!                   ( r_sfc    * r_down_tl * tau(k) ) + &
!                   ( r_sfc_tl * r_down    * tau(k) ) + &
!                   ( ( ( e_sfc * b_sfc ) + ( r_sfc * r_down ) ) * tau_tl(k) )
!
!       the contributions of all the layers except the top layer to the
!       upwelling radiance is accumulated,
!
!                              __ 2
!                             \
!         r_up_tl = r_up_tl +  >  ( b(k) * dtau_tl(k) ) + ( b_tl(k) * dtau(k) )
!                             /__
!                                k=k
!
!       the top layer contribution is then added explicitly,
!
!         r_up_tl = r_up_tl + ( b(1) * ( -tau_tl(1) ) ) + ( b_tl(1) * ( 1 - tau(1) ) )
!
!       to avoid exceeding the arrays bounds of tau (i.e. tau(0) == 1.0 )
!       or tau_tl or requiring them to be dimensioned 0:k.
!
!       the final tangent-linear upwelling radiance is then converted to a
!       tangent-linear  brightness temperature.
!
!s-      
!--------------------------------------------------------------------------------


  subroutine compute_radiance_tl( &
                                  ! -- forward inputs
                                  temperature,              &  ! input, k

                                  surface_temperature,      &  ! input, scalar
                                  surface_emissivity,       &  ! input, scalar
                                  surface_reflectivity,     &  ! input, scalar

                                  tau,                      &  ! input, k
                                  flux_tau,                 &  ! input, k
                                  solar_tau,                &  ! input, scalar

                                  layer_radiance,           &  ! input, k
                                  downwelling_radiance,     &  ! input, scalar
                                  upwelling_radiance,       &  ! input, scalar

                                  ! -- tangent-linear inputs
                                  temperature_tl,           &  ! input, k

                                  surface_temperature_tl,   &  ! input, scalar
                                  surface_emissivity_tl,    &  ! input, scalar
                                  surface_reflectivity_tl,  &  ! input, scalar

                                  tau_tl,                   &  ! input, k
                                  flux_tau_tl,              &  ! input, k
                                  solar_tau_tl,             &  ! input, scalar

                                  ! -- other inputs
                                  secant_solar_angle,       &  ! input, scalar
                                  valid_solar,              &  ! input, scalar
                                  channel_index,            &  ! input, scalar

                                  ! -- tangent-linear outputs
                                  layer_radiance_tl,        &  ! output, k
                                  downwelling_radiance_tl,  &  ! output, scalar
                                  upwelling_radiance_tl,    &  ! output, scalar

                                  brightness_temperature_tl )  ! output, scalar



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( : ), intent( in )  :: temperature

    real( fp_kind ),                 intent( in )  :: surface_temperature
    real( fp_kind ),                 intent( in )  :: surface_emissivity
    real( fp_kind ),                 intent( in )  :: surface_reflectivity

    real( fp_kind ), dimension( : ), intent( in )  :: tau
    real( fp_kind ), dimension( : ), intent( in )  :: flux_tau
    real( fp_kind ),                 intent( in )  :: solar_tau

    real( fp_kind ), dimension( : ), intent( in )  :: layer_radiance
    real( fp_kind ),                 intent( in )  :: downwelling_radiance
    real( fp_kind ),                 intent( in )  :: upwelling_radiance

    ! -- tangent-linear input
    real( fp_kind ), dimension( : ), intent( in )  :: temperature_tl

    real( fp_kind ),                 intent( in )  :: surface_temperature_tl
    real( fp_kind ),                 intent( in )  :: surface_emissivity_tl
    real( fp_kind ),                 intent( in )  :: surface_reflectivity_tl

    real( fp_kind ), dimension( : ), intent( in )  :: tau_tl
    real( fp_kind ), dimension( : ), intent( in )  :: flux_tau_tl
    real( fp_kind ),                 intent( in )  :: solar_tau_tl

    ! -- other input
    real( fp_kind ),                 intent( in )  :: secant_solar_angle
    integer,                         intent( in )  :: valid_solar
    integer,                         intent( in )  :: channel_index

    ! -- tangent-linear output
    real( fp_kind ), dimension( : ), intent( out ) :: layer_radiance_tl
    real( fp_kind ),                 intent( out ) :: downwelling_radiance_tl
    real( fp_kind ),                 intent( out ) :: upwelling_radiance_tl

    real( fp_kind ),                 intent( out ) :: brightness_temperature_tl

 

    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_radiance_tl'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: k, n_layers
    integer :: l

    real( fp_kind ) :: surface_b
    real( fp_kind ) :: surface_b_tl


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic size



    !#--------------------------------------------------------------------------#
    !#                   -- determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    n_layers = size( temperature )



    !#--------------------------------------------------------------------------#
    !#       -- calculate the tangent-linear downwelling thermal flux --        #
    !#--------------------------------------------------------------------------#

    ! -- assign the channel index to a short name
    l = channel_index


    ! ---------------------------------------------
    ! initialise the tangent-linear downwelling
    ! radiance to the space emission boundary term
    ! reaching the surface. the cosmic background
    ! radiance is zero for infrared channels and
    ! precalculated for microwave channels. the
    ! emissivity of space is assumed to be 1.0.
    !
    ! cosmic background data from the
    ! spectral_coefficients module
    ! ---------------------------------------------

    downwelling_radiance_tl = cosmic_background_radiance( l ) * flux_tau_tl( 1 )


    ! --------------------------------
    ! loop over layers from toa->sfc-1
    ! --------------------------------

    k_down_layer_loop: do k = 1, n_layers - 1

      ! -- calculate the tangent-linear layer planck radiance
      call sensor_planck_radiance_tl( l,                     &  ! input
                                      temperature( k ),      &  ! input
                                      temperature_tl( k ),   &  ! input
                                      layer_radiance_tl( k ) )  ! output

      ! -- accumulate tangent-linear absorption and emission for current layer.
      ! -- lte assumed.
      downwelling_radiance_tl = downwelling_radiance_tl + &
                                ( layer_radiance(k)    * ( flux_tau_tl(k+1) - flux_tau_tl(k) ) ) + &
                                ( layer_radiance_tl(k) * ( flux_tau(k+1)    - flux_tau(k)    ) )

    end do k_down_layer_loop


    ! ----------------------------------
    ! flux bottom layer (closest to sfc)
    ! ----------------------------------

    ! -- lowest layer tangent-linear planck radiance
    call sensor_planck_radiance_tl( l,                            &  ! input
                                    temperature( n_layers ),      &  ! input
                                    temperature_tl( n_layers ),   &  ! input
                                    layer_radiance_tl( n_layers ) )  ! output

    ! -- contribution of lowest layer. note that at the
    ! -- surface, a transmittance of 1.0 is used.
    downwelling_radiance_tl = downwelling_radiance_tl + &
                              ( layer_radiance(n_layers)    * (      -flux_tau_tl(n_layers) ) ) + &
                              ( layer_radiance_tl(n_layers) * ( one - flux_tau(n_layers)    ) )



    !#--------------------------------------------------------------------------#
    !#               -- calculate the tangent-linear solar term --              #
    !#--------------------------------------------------------------------------#

    solar_term: if ( valid_solar == 1 ) then

      downwelling_radiance_tl = downwelling_radiance_tl + &

                                ( solar_irradiance(l) * solar_tau_tl / &
      !                           ------------------------------------
                                           secant_solar_angle          )

    end if solar_term



    !#--------------------------------------------------------------------------#
    !#      -- reflect the tangent-linear downwelling radiance, add the  --     #
    !#      -- surface emission and use it to initialise the tangent-    --     #
    !#      -- linear upwelling radiance                                 --     #
    !#--------------------------------------------------------------------------#

    ! -- calculate the surface terms
    call sensor_planck_radiance( l,                   &  ! input
                                 surface_temperature, &  ! input
                                 surface_b            )  ! output

    call sensor_planck_radiance_tl( l,                      &  ! input
                                    surface_temperature,    &  ! input
                                    surface_temperature_tl, &  ! input
                                    surface_b_tl            )  ! output

    ! -- initialise the tangent-linear upwelling radiance
    upwelling_radiance_tl = ( tau(n_layers) * surface_emissivity   * surface_b_tl            ) + &
                            ( tau(n_layers) * surface_reflectivity * downwelling_radiance_tl ) + &

                            ( ((surface_emissivity   * surface_b           ) + &
                               (surface_reflectivity * downwelling_radiance)) * tau_tl(n_layers) ) + &

                            ( tau(n_layers) * surface_b            * surface_emissivity_tl   ) + &
                            ( tau(n_layers) * downwelling_radiance * surface_reflectivity_tl )



    !#--------------------------------------------------------------------------#
    !#           -- calculate the tangent-linear upwelling radiance --          #
    !#--------------------------------------------------------------------------#

    ! --------------------------------
    ! loop over layers from sfc->toa-1
    ! --------------------------------

    k_up_layer_loop: do k = n_layers, 2, -1

      upwelling_radiance_tl = upwelling_radiance_tl + &
                              ( layer_radiance(k)    * ( tau_tl(k-1) - tau_tl(k) ) ) + &
                              ( layer_radiance_tl(k) * ( tau(k-1)    - tau(k)    ) )

    end do k_up_layer_loop


    ! --------------------------
    ! top layer (closest to toa)
    ! --------------------------

    upwelling_radiance_tl = upwelling_radiance_tl + &
                            ( layer_radiance(1)    * (      -tau_tl(1) ) ) + &
                            ( layer_radiance_tl(1) * ( one - tau(1)    ) )



    !#--------------------------------------------------------------------------#
    !#    -- convert the tangent-linear radiance to brightness temperature --   #
    !#--------------------------------------------------------------------------#

    call sensor_planck_temperature_tl( l,                        &  ! input
                                       upwelling_radiance,       &  ! input
                                       upwelling_radiance_tl,    &  ! input
                                       brightness_temperature_tl )  ! output

  end subroutine compute_radiance_tl






!--------------------------------------------------------------------------------
!s+
! name:
!       compute_radiance_ad
!
! purpose:
!       public subroutine to calculate the adjoint of the toa radiance and
!       brightness temperature.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_radiance_ad( &
!                                 ! -- forward inputs
!                                 temperature,               &  ! input, k
!
!                                 surface_temperature,       &  ! input, scalar
!                                 surface_emissivity,        &  ! input, scalar
!                                 surface_reflectivity,      &  ! input, scalar
!
!                                 tau,                       &  ! input, k
!                                 flux_tau,                  &  ! input, k
!                                 solar_tau,                 &  ! input, scalar
!
!                                 layer_radiance,            &  ! input, k
!                                 downwelling_radiance,      &  ! input, scalar
!                                 upwelling_radiance,        &  ! input, scalar
!
!                                 ! -- adjoint inputs
!                                 layer_radiance_ad,         &  ! in/output, k
!                                 downwelling_radiance_ad,   &  ! in/output, scalar
!                                 upwelling_radiance_ad,     &  ! in/output, scalar
!
!                                 brightness_temperature_ad, &  ! in/output, scalar
!
!                                 ! -- other inputs
!                                 secant_solar_angle,        &  ! input, scalar
!                                 valid_solar,               &  ! input, scalar
!                                 channel_index,             &  ! input, scalar
!
!                                 ! -- adjoint outputs
!                                 temperature_ad,            &  ! in/output, k
!
!                                 surface_temperature_ad,    &  ! in/output, scalar
!                                 surface_emissivity_ad,     &  ! in/output, scalar
!                                 surface_reflectivity_ad,   &  ! in/output, scalar
!
!                                 tau_ad,                    &  ! in/output, k
!                                 flux_tau_ad,               &  ! in/output, k
!                                 solar_tau_ad               )  ! in/output, scalar
!
!
! input arguments:
!       temperature:               profile layer average temperature array.
!                                  units:      kelvin
!                                  type:       float
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       surface_temperature:       surface boundary temperature.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_emissivity:        surface boundary emissivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_reflectivity:      surface boundary reflectivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       tau:                       layer-to-space transmittance profile for
!                                  a particular satellite view angle.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       flux_tau:                  layer-to-surface transmittance profile for
!                                  either the diffuse approximation angle (ir)
!                                  or the satellite view angle (mw). the latter
!                                  assumes specular reflectivity.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       solar_tau:                 total space-to-surface transmittance at the
!                                  solar zenith angle.
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       layer_radiance:            channel planck radiance for every layer.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       downwelling_radiance:      channel radiance at surface due to downwelling
!                                  flux and solar components.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       upwelling_radiance:        channel toa radiance simulating the satellite
!                                  sensor measurement. this is composed of the
!                                  reflected downwelling propagated through the
!                                  atmosphere as well as the upwelling only component.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       layer_radiance_ad:         adjoint of the channel planck radiance for every
!                                  layer.
!                                  ** this argument is set to zero on output **.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       downwelling_radiance_ad:   adjoint of the channel radiance at surface due
!                                  to downwelling flux and solar components.
!                                  ** this argument is set to zero on output **.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       upwelling_radiance_ad:     adjoint of the channel toa radiance simulating
!                                  the satellite sensor measurement.
!                                  ** this argument is set to zero on output **.
!                                  units:      mw/(m^2.sr.cm^-1)
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       brightness_temperature_ad: adjoint of the channel toa brightness temperature.
!                                  ** this argument is set to zero on output **.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       secant_solar_angle:        secant of the solar zenith angle corresponding
!                                  to that used in calculating the total solar
!                                  transmittance.
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       valid_solar:               flag indicating if the solar component should
!                                  be included.
!                                  if = 0, no solar (if sensor channel frequency
!                                          is less than a preset cutoff or if solar
!                                          zenith angle is greater than its preset
!                                           cutoff.)
!                                     = 1, include solar
!                                  units:      none.
!                                  type:       integer
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       channel_index:             channel index id. this is a unique index
!                                  to a (supported) sensor channel.
!                                  units:      none
!                                  type:       integer
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!
! optional input arguments:
!       none.
!
! output arguments:
!       temperature_ad:            adjoint of the profile layer average temperature.
!                                  units:      kelvin
!                                  type:       float
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       surface_temperature_ad:    adjoint of the surface boundary temperature.
!                                  units:      kelvin
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_emissivity_ad:     adjoint of the surface boundary emissivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       surface_reflectivity_ad:   adjoint of the surface boundary reflectivity
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
!       tau_ad:                    adjoint of the layer-to-space transmittance
!                                  profile for a particular satellite view angle.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       flux_tau_ad:               adjoint of the layer-to-surface transmittance
!                                  profile for either the diffuse approximation
!                                  angle (ir) or the satellite view angle (mw).
!                                  the latter assumes specular reflectivity.
!                                  units:      none
!                                  type:       real
!                                  dimension:  k
!                                  attributes: intent( in )
!
!       solar_tau_ad:              adjoint of the total space-to-surface 
!                                  transmittance at the solar zenith angle.
!                                  units:      none
!                                  type:       real
!                                  dimension:  scalar
!                                  attributes: intent( in )
!
! optional output arguments:
!       none.
!
! calls:
!       sensor_planck_radiance_ad:    function to compute the adjoint
!                                     planck radiance.
!                                     source: sensor_planck_routines module
!
!       sensor_planck_temperature_ad: function to compute the adjoint
!                                     planck temperature.
!                                     source: sensor_planck_routines module
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
!       none.
!
! procedure:
!       this adjoint code is derived directly from the transpose of the
!       tangent-linear model. starting with the final statement of the
!       tangent-linear model,
!
!         r_up_tl = r_up_tl + ( b(1) * ( -tau_tl(1) ) ) + ( b_tl(1) * ( 1 - tau(1) ) )
!
!       the adjoint of the components is given by,
!
!         tau_ad(1) = -b(1) * r_up_ad
!         b_ad(1)   = ( 1 - tau(1) ) * r_up_ad
!
!       r_up_ad is not set to zero at this point as it is being summed. the 
!       surface to space tangent-linear summation is given by,
!
!                              __ 2
!                             \
!         r_up_tl = r_up_tl +  >  ( b(k) * dtau_tl(k) ) + ( b_tl(k) * dtau(k) )
!                             /__
!                                k=k
!
!       the adjoint of the components are,
!
!         b_ad(k)     = dtau(k) * r_up_ad
!         tau_ad(k)   = tau_ad(k)   - ( b(k) * r_up_ad )
!         tau_ad(k-1) = tau_ad(k-1) + ( b(k) * r_up_ad )
!         t_ad(k)     = planck_ad(t(k),b_ad(k))
!         b_ad(k)     = 0.0
!
!       next comes the tangent-linea surface term,
!
!         r_up_tl = ( e_sfc    * b_sfc_tl  * tau(k) ) + &
!                   ( e_sfc_tl * b_sfc     * tau(k) ) + &
!                   ( r_sfc    * r_down_tl * tau(k) ) + &
!                   ( r_sfc_tl * r_down    * tau(k) ) + &
!                   ( ( ( e_sfc * b_sfc ) + ( r_sfc * r_down ) ) * tau_tl(k) )
!
!       with the adjoints,
!
!         r_ad      = tau(k) * r_down * r_up_ad
!         e_ad      = tau(k) * b(sfc) * r_up_ad
!         tau_ad(k) = tau_ad(k) + ( e.b(sfc) + r.r_down ) * r_up_ad
!         r_down_ad = r_down_ad + ( tau(k) * r * r_up_ad )
!         b_ad(sfc) =               tau(k) * e * r_up_ad
!         r_up_ad   = 0.0
!         t_ad(sfc) = planck_ad(t(sfc),b_ad(sfc))
!
!       the tangent-linear solar term, if used, is
!
!         r_down_tl = r_down_tl + ( solar_irradiance * solar_tau_tl * cos(solar_theta) )
!
!       with its adjoint being,
!
!         solar_tau_ad = solar_tau_ad + ( solar_irradiance * cos(solar_theta) * r_down_ad )
!
!       the tangent-linear surface layer contribution to the downwelling flux
!       is calculated separately,
!
!         r_down_tl = r_down_tl + ( b(k) * ( -flux_tau_tl(k) ) ) + ( b_tl(k) * ( 1 - flux_tau(k) ) )
!         
!       and its component adjoints are,
!
!         b_ad(k)        = b_ad(k) + ( 1 - flux_tau(k) ) * r_down_ad
!         flux_tau_ad(k) = flux_tau_ad(k) - ( b(k) * r_down_ad )
!         t_ad(k)        = planck_ad(t(k),b_ad(k))
!         b_ad(k)        = 0.0
!
!       as with the upwelling tangent-linear summation, the downwelling is
!       given by,
!
!                                  __k-1
!                                 \
!         r_down_tl = r_down_tl +  >  ( b(k) * dflux_tau_tl(k) ) + ( b_tl(k) * dflux_tau(k) )
!                                 /__
!                                    k=1
!
!       the adjoint of the components are,
!
!         b_ad(k)          = b_ad(k ) + ( dflux_tau(k) * r_down_ad )
!         flux_tau_ad(k)   = flux_tau_ad(k)   - ( b(k) * r_down_ad )
!         flux_tau_ad(k-1) = flux_tau_ad(k-1) + ( b(k) * r_down_ad )
!         t_ad(k)          = planck_ad(t(k),b_ad(k))
!         b_ad(k)          = 0.0
!
!       the final step is to determine the adjoints of the cosmic background
!       tangent-linear term,
!
!         r_down_tl = cbr * flux_tau_tl(1)
!
!       which is
!
!         flux_tau_ad(1) = flux_tau_ad(1) + ( cbr * r_down_ad )
!         r_down_ad      = 0.0
!
!s-      
!--------------------------------------------------------------------------------


  subroutine compute_radiance_ad( &
                                  ! -- forward inputs
                                  temperature,               &  ! input, k

                                  surface_temperature,       &  ! input, scalar
                                  surface_emissivity,        &  ! input, scalar
                                  surface_reflectivity,      &  ! input, scalar

                                  tau,                       &  ! input, k
                                  flux_tau,                  &  ! input, k
                                  solar_tau,                 &  ! input, scalar

                                  layer_radiance,            &  ! input, k
                                  downwelling_radiance,      &  ! input, scalar
                                  upwelling_radiance,        &  ! input, scalar

                                  ! -- adjoint inputs
                                  layer_radiance_ad,         &  ! in/output, k
                                  downwelling_radiance_ad,   &  ! in/output, scalar
                                  upwelling_radiance_ad,     &  ! in/output, scalar

                                  brightness_temperature_ad, &  ! in/output, scalar

                                  ! -- other inputs
                                  secant_solar_angle,        &  ! input, scalar
                                  valid_solar,               &  ! input, scalar
                                  channel_index,             &  ! input, scalar

                                  ! -- adjoint outputs
                                  temperature_ad,            &  ! in/output, k

                                  surface_temperature_ad,    &  ! in/output, scalar
                                  surface_emissivity_ad,     &  ! in/output, scalar
                                  surface_reflectivity_ad,   &  ! in/output, scalar

                                  tau_ad,                    &  ! in/output, k
                                  flux_tau_ad,               &  ! in/output, k
                                  solar_tau_ad               )  ! in/output, scalar




    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( : ), intent( in )     :: temperature

    real( fp_kind ),                 intent( in )     :: surface_temperature
    real( fp_kind ),                 intent( in )     :: surface_emissivity
    real( fp_kind ),                 intent( in )     :: surface_reflectivity

    real( fp_kind ), dimension( : ), intent( in )     :: tau
    real( fp_kind ), dimension( : ), intent( in )     :: flux_tau
    real( fp_kind ),                 intent( in )     :: solar_tau

    real( fp_kind ), dimension( : ), intent( in )     :: layer_radiance
    real( fp_kind ),                 intent( in )     :: downwelling_radiance
    real( fp_kind ),                 intent( in )     :: upwelling_radiance

    ! -- adjoint inputs
    real( fp_kind ), dimension( : ), intent( in out ) :: layer_radiance_ad
    real( fp_kind ),                 intent( in out ) :: downwelling_radiance_ad
    real( fp_kind ),                 intent( in out ) :: upwelling_radiance_ad

    real( fp_kind ),                 intent( in out ) :: brightness_temperature_ad

    ! -- other input
    real( fp_kind ),                 intent( in )     :: secant_solar_angle
    integer,                         intent( in )     :: valid_solar
    integer,                         intent( in )     :: channel_index

    ! -- adjoint outputs
    real( fp_kind ), dimension( : ), intent( in out ) :: temperature_ad

    real( fp_kind ),                 intent( in out ) :: surface_temperature_ad
    real( fp_kind ),                 intent( in out ) :: surface_emissivity_ad
    real( fp_kind ),                 intent( in out ) :: surface_reflectivity_ad

    real( fp_kind ), dimension( : ), intent( in out ) :: tau_ad
    real( fp_kind ), dimension( : ), intent( in out ) :: flux_tau_ad
    real( fp_kind ),                 intent( in out ) :: solar_tau_ad


 

    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_radiance_ad'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: k, n_layers
    integer :: l

    real( fp_kind ) :: surface_b
    real( fp_kind ) :: surface_b_ad

    real( fp_kind ) :: b_ur_ad
    real( fp_kind ) :: b_dr_ad


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic size



    !#--------------------------------------------------------------------------#
    !#                   -- determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    n_layers = size( temperature )



    !#--------------------------------------------------------------------------#
    !#         -- calculate the adjoint of the brightness temperature --        #
    !#--------------------------------------------------------------------------#

    ! -- assign the channel index to a short name
    l = channel_index

    ! -- upwelling radiance adjoint
    call sensor_planck_temperature_ad( l,                         &  ! input
                                       upwelling_radiance,        &  ! input
                                       brightness_temperature_ad, &  ! input
                                       upwelling_radiance_ad      )  ! in/output
    brightness_temperature_ad = zero



    !#--------------------------------------------------------------------------#
    !#      -- calculate the adjoints of the upwelling radiance term --         #
    !#--------------------------------------------------------------------------#

    ! --------------------------
    ! top layer (closest to toa)
    ! --------------------------

    ! -- adjoint of top layer radiance
    layer_radiance_ad( 1 ) = layer_radiance_ad( 1 ) + ( ( one - tau( 1 ) ) * upwelling_radiance_ad )

    ! -- adjoint of top layer dtau
    tau_ad( 1 ) = tau_ad( 1 ) + ( -layer_radiance( 1 ) * upwelling_radiance_ad )
    ! note: no upwelling_radiance_ad = 0 here since
    !       upwelling_radiance_tl = upwelling_radiance_tl + (...)

    ! -- adjoint of top layer temperature
    call sensor_planck_radiance_ad( l,                    &  ! input
                                    temperature(1),       &  ! input
                                    layer_radiance_ad(1), &  ! input
                                    temperature_ad(1)     )  ! in/output
    layer_radiance_ad( 1 ) = zero


    ! --------------------------------
    ! loop over layers from toa-1->sfc
    ! --------------------------------

    k_down_layer_loop: do k = 2, n_layers

      ! -- adjoint of layer radiance
      layer_radiance_ad( k ) = layer_radiance_ad( k ) + &
                               ( ( tau( k-1 ) - tau( k ) ) * upwelling_radiance_ad )

      ! -- adjoint of dtau
      b_ur_ad = layer_radiance( k ) * upwelling_radiance_ad
      tau_ad( k )   = tau_ad(k)     - b_ur_ad
      tau_ad( k-1 ) = tau_ad( k-1 ) + b_ur_ad
      ! note: no upwelling_radiance_ad = 0 here since
      !       upwelling_radiance_tl = upwelling_radiance_tl + (...)

      ! -- adjoint of layer temperature
      call sensor_planck_radiance_ad( l,                      &  ! input
                                      temperature( k ),       &  ! input
                                      layer_radiance_ad( k ), &  ! input
                                      temperature_ad( k )     )  ! in/output
      layer_radiance_ad( k ) = zero
 
    end do k_down_layer_loop



    !#--------------------------------------------------------------------------#
    !#     -- calculate the adjoints of the surface tangent-linear terms --     #
    !#--------------------------------------------------------------------------#

    ! -- recalculate surface planck radiance
    call sensor_planck_radiance( l,                   &  ! input
                                 surface_temperature, &  ! input
                                 surface_b            )  ! output


    ! -- surface property adjoints
    surface_reflectivity_ad = surface_reflectivity_ad + &
                              ( tau( n_layers ) * downwelling_radiance * upwelling_radiance_ad )
 
    surface_emissivity_ad   = surface_emissivity_ad + &
                              ( tau( n_layers ) * surface_b * upwelling_radiance_ad )
 
    ! -- total transmittance adjoint
    tau_ad( n_layers ) = tau_ad( n_layers ) + &
                         ( ( ( surface_emissivity   * surface_b            ) + &
                             ( surface_reflectivity * downwelling_radiance ) ) * upwelling_radiance_ad )
 
    ! -- downweling radiance adjoint
    downwelling_radiance_ad = downwelling_radiance_ad + &
                              ( tau( n_layers ) * surface_reflectivity * upwelling_radiance_ad )
 
    ! -- surface emission adjoint. this quantity
    ! -- is initialised each call.
    surface_b_ad = ( tau( n_layers ) * surface_emissivity * upwelling_radiance_ad )

    ! -- set upwelling adjoint to zero.
    ! --  no more impact on gradient vector
    upwelling_radiance_ad = zero

    ! -- surface temperature adjoint
    call sensor_planck_radiance_ad( l,                     & ! input
                                    surface_temperature,   & ! input
                                    surface_b_ad,          & ! input
                                    surface_temperature_ad ) ! in/output
    ! -- no need to zero surface_b_ad as it is initialised each call.
  


    !#--------------------------------------------------------------------------#
    !#      -- calculate the adjoints of the solar transmittance term --        #
    !#--------------------------------------------------------------------------#

    solar_term: if ( valid_solar == 1 ) then

      solar_tau_ad = solar_tau_ad + &
                     ( ( solar_irradiance( l ) / secant_solar_angle ) * downwelling_radiance_ad )

      ! note: no downwelling_radiance_ad = 0 here since
      !       downwelling_radiance_tl = downwelling_radiance_tl + (tl solar term)

    end if solar_term




    !#--------------------------------------------------------------------------#
    !#        -- calculate the adjoints of the downwelling flux term --         #
    !#--------------------------------------------------------------------------#

    ! -----------------------------
    ! bottom layer (closest to sfc)
    ! -----------------------------

    ! -- adjoint of bottom layer radiance
    layer_radiance_ad( n_layers ) = layer_radiance_ad( n_layers ) + &
                                    ( ( one - flux_tau( n_layers ) ) * downwelling_radiance_ad )

    ! -- adjoint of flux transmittance
    flux_tau_ad( n_layers ) = flux_tau_ad( n_layers ) - &
                              ( layer_radiance( n_layers ) * downwelling_radiance_ad )
    ! note: no downwelling_radiance_ad = 0 here since
    !       downwelling_radiance_tl = downwelling_radiance_tl + (tl flux term)
 
    ! -- adjoint of layer temperature
    call sensor_planck_radiance_ad( l,                             &  ! input
                                    temperature( n_layers ),       &  ! input
                                    layer_radiance_ad( n_layers ), &  ! input
                                    temperature_ad( n_layers )     )  ! in/output
    layer_radiance_ad( n_layers ) = zero


    ! --------------------------------
    ! loop over layers from sfc-1->toa
    ! --------------------------------

    k_up_layer_loop: do k = n_layers - 1, 1, -1

      ! -- adjoint of layer radiance
      layer_radiance_ad( k ) = layer_radiance_ad( k ) + &
                               ( ( flux_tau( k+1 ) - flux_tau( k ) ) * downwelling_radiance_ad )

      ! -- adjoint of dtau
      b_dr_ad = layer_radiance( k ) * downwelling_radiance_ad
      flux_tau_ad( k )   = flux_tau_ad( k )   - b_dr_ad
      flux_tau_ad( k+1 ) = flux_tau_ad( k+1 ) + b_dr_ad
      ! note: no downwelling_radiance_ad = 0 here since
      !       downwelling_radiance_tl = downwelling_radiance_tl + (...)

      ! -- adjoint of layer temperature
      call sensor_planck_radiance_ad( l,                 &  ! input
                                 temperature( k ),       &  ! input
                                 layer_radiance_ad( k ), &  ! input
                                 temperature_ad( k )     )  ! in/output
      layer_radiance_ad( k ) = zero

    end do k_up_layer_loop


    ! --------------------------------------------
    ! background term. note that the emissivity of
    ! space is implicitly assumed to 1.0.
    ! --------------------------------------------

    flux_tau_ad( 1 ) = flux_tau_ad( 1 ) + &
                       ( cosmic_background_radiance( l ) * downwelling_radiance_ad )

    downwelling_radiance_ad = zero

  end subroutine compute_radiance_ad

end module radiance


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
! revision 2.4  2001/08/16 17:11:57  paulv
! - updated documentation
!
! revision 2.3  2001/08/08 20:02:12  paulv
! - removed sensor planck function routines and placed them in their own
!   module, sensor_planck_routines. some routines were required for other
!   uses so their private subprogram status wasn't amenable to code-sharing.
! - updated header documentation.
!
! revision 2.2  2001/08/01 16:58:32  paulv
! - corrected bug in compute_radiance() function. the initialisation
!   statement of the downwelling radiance was,
!      downwelling_radiance = cosmic_background_radiance( l )
!   and has been changed to,
!      downwelling_radiance = cosmic_background_radiance( l ) * flux_tau( 1 )
!   i.e. the transmission of the space emission term to the surface. this
!   was a holdover from earlier versions of the functions when the transmittances
!   were calculated and passed as *layer* rather than layer-to-surface
!   transmittances.
! - removed initialisation and zeroing of the adjoint of the surface emissioni
!   term surface_b_ad. this is used in only one place so there is no need to
!   do,
!     surface_b_ad = zero
!     surface_b_ad = surface_b_ad + &
!                    ( tau( n_layers ) * surface_emissivity * upwelling_radiance_ad )
!     ....use surface_b_ad...
!     surface_b_ad = zero
!   when,
!     surface_b_ad = ( tau( n_layers ) * surface_emissivity * upwelling_radiance_ad )
!   will do.
! - updated documentation.
!
! revision 2.1  2001/05/29 18:05:29  paulv
! - all tangent-linear and adjoint routines included.
! - no more optional arguments of downwelling flux and surface reflectivity -
!   they are expected.
! - altered the method of calculating the layer contributions. changed code
!   from:
!     layer_radiance(k) = (1-tau)*b(t) + tau*layer_radiance(k-1)
!   to:
!     layer_radiance(k) = b(t) * dtau
!
! revision 1.5  2001/01/24 20:14:21  paulv
! - latest test versions.
!
! revision 1.4  2000/11/09 20:46:07  paulv
! - added solar term.
! - downwelling flux transmittance term is now an optional argument.
!   if not specified, the layer radiances are calculated during the
!   upwelling radiance integration.
! - surface reflectivity is an optional argument. if not specified the
!   surface emissivity is used to generate an isotropic reflectivity.
!
! revision 1.3  2000/08/31 19:36:33  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.2  2000/08/24 15:48:34  paulv
! - replaced "regular" reflectivity for reflected downwelling thermal with
!   the isotropic reflectivity in the compute_radiance subprogram.
! - updated module and subprogram documentation.
!
! revision 1.1  2000/08/21 20:59:34  paulv
! initial checkin.
!
!
!
!
