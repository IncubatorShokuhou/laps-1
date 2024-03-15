!------------------------------------------------------------------------------
!m+
! name:
!       transmittance
!
! purpose:
!       rt model transmittance module
!
! category:
!       ncep rtm
!
! calling sequence:
!       use transmittance
!
! outputs:
!       none.
!
! modules:
!       type_kinds:                 module containing data type kind definitions.
!
!       parameters:                 module containing parameter definitions for the
!                                   rt model.
!
!       absorber_space:             module containing the rt model absorber space
!                                   definitions.
!
!       transmitance_coefficients:  module containing the rt model absorber space
!                                   definitions and transmittance coefficients.
!
!       predictors:                 module containing the predictor calculation
!                                   routines and data.
!
! contains:
!       compute_transmittance:      public subroutine to calculate the layer
!                                   transmittances of an input atmospheric profile.
!
!       compute_transmittance_tl:   public subroutine to calculate the tangent-linear
!                                   layer transmittances of an input atmospheric profile.
!
!       compute_transmittance_ad:   public subroutine to calculate the layer transmittance
!                                   adjoints of an input atmospheric profile.
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
!       adapted from code written by: thomas j.kleespies
!                                     noaa/nesdis/ora
!                                     thomas.j.kleespies@noaa.gov
!
!                                     and
!
!                                     john derber
!                                     noaa/ncep/emc
!                                     john.derber@noaa.gov
!
!
!  copyright (c) 2000 thomas kleespies, john derber, paul van delst
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

module transmittance

  ! ------------
  ! module usage
  ! ------------

  use type_kinds, only : fp_kind
  use parameters
  use transmittance_coefficients


  ! -----------------------
  ! disable implicit typing
  ! -----------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private
  public :: compute_transmittance
  public :: compute_transmittance_tl
  public :: compute_transmittance_ad


contains


!------------------------------------------------------------------------------
!s+
! name:
!       compute_transmittance
!
! purpose:
!       public subroutine to calculate the layer transmittances given an
!       input atmospheric profile.
!
! calling sequence:
!       call compute_transmittance( absorber,      &   ! input, 0:k x j
!                                   predictor,     &   ! input, i x k
!                                   layer_index,   &   ! input, k x j
!                                   channel_index, &   ! input, scalar
!                                   direction,     &   ! input, scalar
!
!                                   tau            )   ! output, k
!
! input arguments:
!       absorber:         profile level integrated absorber amount array.
!                         units:      varies with absorber.
!                         type:       real( fp_kind )
!                         dimension:  0:k x j
!                         attributes: intent( in )
!
!       predictor:        profile layer predictors array.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  i x k
!                         attributes: intent( in )
!
!       layer_index:      index array array associating the input absorber
!                         layer amount to the bracketing absorber space levels.
!                         units:      none
!                         type:       integer
!                         dimension:  k x j
!                         attributes: intent( in )
!
!       channel_index:    channel index id. this is a unique index associated
!                         with a (supported) sensor channel.
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in )
!
!       direction:        direction identifier.
!                         if = 0, calculate layer->surface transmittances (i.e. down)
!                            = 1, calculate layer->space   transmittances (i.e. up)
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in )
!
! optional input arguments:
!        none.
!
! output arguments:
!        tau:             layer to boundary transmittances for the input atmosphere
!                         and channel.
!                         units:      none
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( out )
!
! optional output arguments:
!       none.
!
! calls:
!       none.
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
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!       the layer absorption coefficient is calculated using,
!
!                                   __ i
!                                  \   _
!         abs_coeff(k) = b(0,k,l) + >  b(i,k,l).x(k,l)
!                                  /__
!                                     i=1
!             _
!       where b == interpolated coefficients
!             x == predictors
!             i == predictor index
!             k == layer index
!             l == channel index
!
!       the input coefficients are linearly interpolated in absorber space
!       from the precalculated absorber levels to that defined by the
!       input profile.
!
!       the absorber layer optical depth is then calculated using,
!
!         optical_depth(k) = abs_coeff(k) * da(k)
!
!       where da == layer absorber difference.
!
!       the layer transmittance is then calculated using,
!
!         tau(k) = exp( -optical_depth(k) )
!s-
!------------------------------------------------------------------------------

  subroutine compute_transmittance( absorber,      &   ! input, 0:k x j
                                    predictor,     &   ! input, i x k
                                    layer_index,   &   ! input, k x j
                                    channel_index, &   ! input, scalar
                                    direction,     &   ! input, scalar, 0==down, 1==up

                                    tau            )   ! output, k



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( 0:, : ), intent( in )  :: absorber         ! input, 0:k x j
    real( fp_kind ), dimension( :, : ),  intent( in )  :: predictor        ! input, i x k
    integer,         dimension( :, : ),  intent( in )  :: layer_index      ! input, k x j
    integer,                             intent( in )  :: channel_index    ! input, scalar
    integer,                             intent( in )  :: direction        ! input, scalar, 0==down, 1==up

    real( fp_kind ), dimension( : ),     intent( out ) :: tau              ! output, k


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_transmittance'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: l
    integer :: k, k1, k2, dk, n_layers
    integer :: j, n_absorbers
    integer :: i, ip, n_predictors

    real( fp_kind ) :: ave_absorber, d_absorber
    real( fp_kind ) :: gradient, b, b1, b2
    real( fp_kind ) :: absorption_coefficient
    real( fp_kind ) :: total_od, od_tolerance

    real( fp_kind ), dimension( size( tau ) ) :: optical_depth


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic exp, &
              present, &
              size



    !#--------------------------------------------------------------------------#
    !#                   -- determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! -- number of atmospheric layers and absorbers. the "-1"
    ! -- for the layer assign is because absorber is based on
    ! -- levels.
    n_layers    = size( absorber, dim = 1 ) - 1
    n_absorbers = size( absorber, dim = 2 )

    ! -- number of predictors *to use*. currently this is
    ! -- pretty much hardwired to 5 but it could change
    ! -- in the future. this value is explicitly
    ! -- checked as the tau_coefficients array contained
    ! -- in the transmittance_coefficients module is
    ! -- allocated dynamically when initialising the 
    ! -- radiative transfer model.
 
    ! -- the "-1" is used because the first element of the
    ! -- predictor index array is used to determine if there
    ! -- is *any* absorption for a particular absorber species.
    n_predictors = size( tau_coefficients, dim = 1 ) - 1



    !#--------------------------------------------------------------------------#
    !#                -- calculate the layer optical depths --                  #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------
    ! assign the channel index to a short name
    ! ----------------------------------------

    l = channel_index


    ! ---------------------------
    ! initilise the optical depth
    ! ---------------------------

    optical_depth( : ) = zero


    ! -----------------------------------------------------
    ! loop over each absorber for optical depth calculation
    ! -----------------------------------------------------

    j_absorber_loop: do j = 1, n_absorbers


      ! -----------------------------------------
      ! check if there is any absorption for this
      ! absorber/channel combination.
      !
      ! this check is the reason why all channels
      ! cannot be processed at once and why the
      ! layer loop is within the absorber loop.
      ! -----------------------------------------

      if ( predictor_index( 0, l, j ) == 0 ) cycle j_absorber_loop



      !#------------------------------------------------------------------------#
      !#                    -- begin loop over layers --                        #
      !#------------------------------------------------------------------------#

      k_layer_od_loop: do k = 1, n_layers


        ! -----------------------------------
        ! calculate the current layer average
        ! absorber amount and difference
        ! -----------------------------------

        ave_absorber = point_5 * ( absorber( k, j ) + absorber( k-1, j ) )
        d_absorber   = absorber( k, j ) - absorber( k-1, j )


        ! -----------------------------------------------------------
        ! to linearly interpolate the tau_coeffs to the actual user
        ! space absorber amount, need the gradient across the layer
        ! -----------------------------------------------------------

        k2 = layer_index( k, j )
        k1 = k2 - 1

        gradient = ( absorber_space_levels( k2, j ) - ave_absorber                   ) / &
        !          -------------------------------------------------------------------
                   ( absorber_space_levels( k2, j ) - absorber_space_levels( k1, j ) )


        ! -------------------------------
        ! calculate absorption coeficient
        ! -------------------------------

        ! -- offset term
        b1 = tau_coefficients( 0, k1, l, j )
        b2 = tau_coefficients( 0, k2, l, j )
        absorption_coefficient = b2 + ( gradient * ( b1 - b2 ) )


        i_predictor_loop: do i = 1, n_predictors

          ! -- current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )
          b  = b2 + ( gradient * ( b1 - b2 ) )

          !  - sum current predictor's contribution
          ip = predictor_index( i, l, j )
          absorption_coefficient = absorption_coefficient + ( b * predictor( ip, k ) ) 

        end do i_predictor_loop


        ! -- **** it would be nice to not need this at all! ****
        absorption_coefficient = max( absorption_coefficient, zero )


        ! -----------------------
        ! calculate optical_depth
        ! -----------------------

        optical_depth( k ) = optical_depth( k ) + &
                             ( absorption_coefficient * d_absorber )


      end do k_layer_od_loop

    end do j_absorber_loop



    !#--------------------------------------------------------------------------#
    !#           -- calculate the layer->boundary transmittances --             #
    !#                                                                          #
    !# this step involves another loop over layers. one *could* reverse the     #
    !# order of the j absorber loop and k layer od loop and calculate the tau   #
    !# values outside the absorber loop. however, this would involve an if test #
    !# for every layer even if predictor_index( 0, l ) == 0. i figured it is    #
    !# better to skip the layer loop altogether if there is no absorption for   #
    !# a particular absorber, despite the need for an extra layer loop here to  #
    !# calculate the transmittance.                                             #
    !#--------------------------------------------------------------------------#

    ! ----------------------
    ! determine loop indices
    ! ----------------------

    if ( direction == up ) then
      ! -- transmittance going up, e.g. upwelling. layer 1 => toa
      k1 = 1
      k2 = n_layers
      dk = 1
    else
      ! -- transmittance going down, e.g. downwelling flux, solar. layer n_layers => sfc
      k1 =  n_layers
      k2 =  1
      dk = -1
    end if


    ! ------------------------------------------
    ! initialise the total optical depth and set
    ! its tolerance value. the tolerance value
    ! should prevent floating point underflows.
    ! ------------------------------------------

    total_od     = zero
    od_tolerance = abs( log( tolerance ) )


    ! ----------------------------------------------
    ! loop over layers for transmittance calculation
    ! ----------------------------------------------

    k_layer_tau_loop: do k = k1, k2, dk

      ! -- update total optical depth
      total_od = total_od + optical_depth( k )

      ! -- if optical depth is < than tolerance, calculate transmittance
      if ( total_od < od_tolerance ) then

        tau( k ) = exp( -total_od )

      else

        ! -- the following *could* be slow if dk is negative
        ! -- if so, then tau( k ) = zero would work.
        tau( k:k2:dk ) = zero
        exit k_layer_tau_loop

      end if

    end do k_layer_tau_loop

  end subroutine compute_transmittance




!------------------------------------------------------------------------------
!s+
! name:
!       compute_transmittance_tl
!
! purpose:
!       public subroutine to calculate the tangent-linear layer transmittances
!       of an input atmospheric profile.
!
! calling sequence:
!        call compute_transmittance_tl( &
!                                       ! -- forward input
!                                       absorber,      &   ! input, 0:k x j
!                                       predictor,     &   ! input, i x k
!                                       tau,           &   ! input, k
!
!                                       ! -- tangent-liner input
!                                       absorber_tl,   &   ! input, 0:k x j
!                                       predictor_tl,  &   ! input, i x k
!
!                                       ! -- other input
!                                       layer_index,   &   ! input, k x j
!                                       channel_index, &   ! input, scalar
!                                       direction,     &   ! input, scalar, 0==down, 1==up
!
!                                       ! -- tangent-liner output
!                                       tau_tl         )   ! output, k
!
! input arguments:
!       absorber:         profile level integrated absorber amount array.
!                         units:      varies with absorber.
!                         type:       real( fp_kind )
!                         dimension:  0:k x j
!                         attributes: intent( in )
!
!       predictor:        profile layer predictors array.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  i x k
!                         attributes: intent( in )
!
!       absorber_tl:      profile level tangent-linear integrated absorber
!                         amount array.
!                         units:      varies with absorber.
!                         type:       real( fp_kind )
!                         dimension:  0:k x j
!                         attributes: intent( in )
!
!       predictor_tl:     profile layer tangent-linear predictors array.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  i x k
!                         attributes: intent( in )
!
!       layer_index:      index array array associating the input absorber
!                         layer amount to the bracketing absorber space levels.
!                         units:      none
!                         type:       integer
!                         dimension:  k x j
!                         attributes: intent( in )
!
!       channel_index:    channel index id. this is a unique index associated
!                         with a (supported) sensor channel.
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in )
!
!       direction:        direction identifier.
!                         if = 0, calculate layer->surface transmittances (i.e. down)
!                            = 1, calculate layer->space   transmittances (i.e. up)
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in )
!
! optional input arguments:
!        none.
!
! output arguments:
!        tau_tl:          layer to boundary tangent-linear transmittances for the
!                         input atmosphere and channel.
!                         units:      none
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( out )
!
! optional output arguments:
!       none.
!
! calls:
!       none.
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
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!       the layer absorption coefficient is calculated using,
!
!                                   __ i
!                        _         \   _
!         abs_coeff(k) = b(0,k,l) + >  b(i,k,l).x(k,l)
!                                  /__
!                                     i=1
!
!       and the tangent-linear absorption coefficient from,
!
!                                         __ i
!                           _            \   _                    _
!         abs_coeff_tl(k) = b_tl(0,k,l) + >  b_tl(i,k,l).x(k,l) + b(i,k,l).x_tl(k,l)
!                                        /__
!                                           i=1
!             _  _
!       where b, b_tl == interpolated coefficients
!             x, x_tl == predictors
!             i       == predictor index
!             k       == layer index
!             l       == channel index
!
!       the input coefficients are linearly interpolated in absorber space
!       from the precalculated absorber levels to that defined by the
!       input profile.
!
!       the absorber layer optical depth is then calculated using,
!
!         optical_depth(k) = abs_coeff(k).da(k)
!
!       and the tangent-linear term is determined from,
!
!         optical_depth_tl(k) = abs_coeff_tl(k).da(k) + abs_coeff(k).da_tl(k)
!
!       where da, da_tl == layer absorber difference.
!
!       the tangent-linear transmittance is then calculated using,
!
!                               __ k
!                              \
!         tau_tl(k) = -tau(k) . >  optical_depth_tl(i)
!                              /__
!                                 i=k
!
!s-
!------------------------------------------------------------------------------

  subroutine compute_transmittance_tl( &
                                       ! -- forward input
                                       absorber,      &   ! input, 0:k x j
                                       predictor,     &   ! input, i x k
                                       tau,           &   ! input, k

                                       ! -- tangent-liner input
                                       absorber_tl,   &   ! input, 0:k x j
                                       predictor_tl,  &   ! input, i x k

                                       ! -- other input
                                       layer_index,   &   ! input, k x j
                                       channel_index, &   ! input, scalar
                                       direction,     &   ! input, scalar, 0==down, 1==up

                                       ! -- tangent-liner output
                                       tau_tl         )   ! output, k



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( 0:, : ), intent( in )  :: absorber         ! 0:k x j
    real( fp_kind ), dimension( :, : ),  intent( in )  :: predictor        ! i x k
    real( fp_kind ), dimension( : ),     intent( in )  :: tau              ! k

    ! -- tangent_linear input
    real( fp_kind ), dimension( 0:, : ), intent( in )  :: absorber_tl      ! 0:k x j
    real( fp_kind ), dimension( :, : ),  intent( in )  :: predictor_tl     ! i x k

    ! -- other input
    integer,         dimension( :, : ),  intent( in )  :: layer_index      ! k x j
    integer,                             intent( in )  :: channel_index    ! scalar
    integer,                             intent( in )  :: direction        ! scalar, 0==down, 1==up

    ! -- tangent_linear output
    real( fp_kind ), dimension( : ),     intent( out ) :: tau_tl           ! k


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_transmittance_tl'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: l
    integer :: k, k1, k2, dk, n_layers
    integer :: j, n_absorbers
    integer :: i, ip, n_predictors

    real( fp_kind ) :: ave_absorber,    d_absorber
    real( fp_kind ) :: ave_absorber_tl, d_absorber_tl

    real( fp_kind ) :: b1, b2
    real( fp_kind ) :: gradient,    b
    real( fp_kind ) :: gradient_tl, b_tl

    real( fp_kind ) :: absorption_coefficient
    real( fp_kind ) :: absorption_coefficient_tl

    real( fp_kind ) :: total_od_tl

    real( fp_kind ), dimension( size( tau ) ) :: optical_depth_tl


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic exp, &
              present, &
              size



    !#--------------------------------------------------------------------------#
    !#                   -- determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! -- number of atmospheric layers and absorbers. the "-1"
    ! -- for the layer assign is because absorber is based on
    ! -- levels.
    n_layers    = size( absorber, dim = 1 ) - 1
    n_absorbers = size( absorber, dim = 2 )

    ! -- number of predictors *to use*. currently this is
    ! -- pretty much hardwired to 5 but it could change
    ! -- in the future. this value is explicitly
    ! -- checked as the tau_coefficients array contained
    ! -- in the transmittance_coefficients module is
    ! -- allocated dynamically when initialising the 
    ! -- radiative transfer model.
 
    ! -- the "-1" is used because the first element of the
    ! -- predictor index array is used to determine if there
    ! -- is *any* absorption for a particular absorber species.
    n_predictors = size( tau_coefficients, dim = 1 ) - 1



    !#--------------------------------------------------------------------------#
    !#                -- calculate the layer optical depths --                  #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------
    ! assign the channel index to a short name
    ! ----------------------------------------

    l = channel_index


    ! ------------------------------------------
    ! initilise the tangent-linear optical depth
    ! ------------------------------------------

    optical_depth_tl( : ) = zero


    ! -----------------------------------------------------
    ! loop over each absorber for optical depth calculation
    ! -----------------------------------------------------

    j_absorber_loop: do j = 1, n_absorbers


      ! -----------------------------------------
      ! check if there is any absorption for this
      ! absorber/channel combination.
      !
      ! this check is the reason why all channels
      ! cannot be processed at once and why the
      ! layer loop is within the absorber loop.
      ! -----------------------------------------

      if ( predictor_index( 0, l, j ) == 0 ) cycle j_absorber_loop



      !#------------------------------------------------------------------------#
      !#                    -- begin loop over layers --                        #
      !#------------------------------------------------------------------------#

      k_layer_od_loop: do k = 1, n_layers


        ! -----------------------------------
        ! calculate the current layer average
        ! absorber amounts and differences
        ! -----------------------------------

        ave_absorber    = point_5 * ( absorber(    k, j ) + absorber(    k-1, j ) )
        ave_absorber_tl = point_5 * ( absorber_tl( k, j ) + absorber_tl( k-1, j ) )

        d_absorber    = absorber(    k, j ) - absorber(    k-1, j )
        d_absorber_tl = absorber_tl( k, j ) - absorber_tl( k-1, j )


        ! -----------------------------------------------------------
        ! to linearly interpolate the tau_coeffs to the actual user
        ! space absorber amount, need the gradient across the layer
        ! -----------------------------------------------------------

        k2 = layer_index( k, j )
        k1 = k2 - 1

        gradient = ( absorber_space_levels( k2, j ) - ave_absorber                   ) / &
        !          -------------------------------------------------------------------
                   ( absorber_space_levels( k2, j ) - absorber_space_levels( k1, j ) )

        gradient_tl =                        ( -ave_absorber_tl )                         / &
        !             -------------------------------------------------------------------
                      ( absorber_space_levels( k2, j ) - absorber_space_levels( k1, j ) )


        ! --------------------------------
        ! calculate absorption coeficients
        ! --------------------------------

        ! -- offset term
        b1 = tau_coefficients( 0, k1, l, j )
        b2 = tau_coefficients( 0, k2, l, j )

        absorption_coefficient    = b2 + ( gradient    * ( b1 - b2 ) )
        absorption_coefficient_tl =        gradient_tl * ( b1 - b2 )


        i_predictor_loop: do i = 1, n_predictors

          ! -- current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )

          b     = b2 + ( gradient    * ( b1 - b2 ) )
          b_tl  =        gradient_tl * ( b1 - b2 )

          !  - sum current predictor's contribution
          ip = predictor_index( i, l, j )

          absorption_coefficient    = absorption_coefficient + ( b * predictor( ip, k ) ) 
          absorption_coefficient_tl = absorption_coefficient_tl        + &
                                      ( b_tl * predictor(    ip, k ) ) + &
                                      ( b    * predictor_tl( ip, k ) )

        end do i_predictor_loop


        ! -- **** it would be nice to not need this at all! ****
        absorption_coefficient = max( absorption_coefficient, zero )


        ! --------------------------------------
        ! calculate tangent-linear optical depth
        ! --------------------------------------

        optical_depth_tl( k ) = optical_depth_tl( k ) + &
                                ( absorption_coefficient_tl * d_absorber    ) + &
                                ( absorption_coefficient    * d_absorber_tl )

      end do k_layer_od_loop

    end do j_absorber_loop



    !#--------------------------------------------------------------------------#
    !#    -- calculate the layer->boundary tangent-linear transmittances --     #
    !#                                                                          #
    !# this step involves another loop over layers. one *could* reverse the     #
    !# order of the j absorber loop and k layer od loop and calculate the tau   #
    !# values outside the absorber loop. however, this would involve an if test #
    !# for every layer even if predictor_index( 0, l ) == 0. i figured it is    #
    !# better to skip the layer loop altogether if there is no absorption for   #
    !# a particular absorber, despite the need for an extra layer loop here to  #
    !# calculate the transmittance.                                             #
    !#--------------------------------------------------------------------------#

    ! ----------------------
    ! determine loop indices
    ! ----------------------

    if ( direction == up ) then
      ! -- transmittance going up, e.g. upwelling. layer 1 => toa
      k1 = 1
      k2 = n_layers
      dk = 1
    else
      ! -- transmittance going down, e.g. downwelling flux, solar. layer n_layers => sfc
      k1 =  n_layers
      k2 =  1
      dk = -1
    end if


    ! ----------------------------------------------
    ! loop over layers for transmittance calculation
    ! ----------------------------------------------

    total_od_tl = zero

    k_layer_tau_loop: do k = k1, k2, dk

      total_od_tl = total_od_tl + optical_depth_tl( k )
      tau_tl( k ) = -tau( k ) * total_od_tl

    end do k_layer_tau_loop

  end subroutine compute_transmittance_tl





!------------------------------------------------------------------------------
!s+
! name:
!       compute_transmittance_ad
!
! purpose:
!       public subroutine to calculate the adjoint of the layer transmittances
!       of an input atmospheric profile.
!
! calling sequence:
!       call compute_transmittance_ad( &
!                                      ! -- forward input
!                                      absorber,      &   ! input, 0:k x j
!                                      predictor,     &   ! input, i x k
!                                      tau,           &   ! input, k
!
!                                      ! -- adjoint input
!                                      tau_ad,        &   ! in/output, k
!
!                                      ! -- other input
!                                      layer_index,   &   ! input, k x j
!                                      channel_index, &   ! input, scalar
!                                      direction,     &   ! input, scalar, 0==down, 1==up
!
!                                      ! -- adjoint output
!                                      absorber_ad,   &   ! in/output, 0:k x j
!                                      predictor_ad   )   ! in/output, i x k
!
! input arguments:
!       absorber:         profile level integrated absorber amount array.
!                         units:      varies with absorber.
!                         type:       real( fp_kind )
!                         dimension:  0:k x j
!                         attributes: intent( in )
!
!       predictor:        profile layer predictors array.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  i x k
!                         attributes: intent( in )
!
!       tau:              layer to boundary transmittances for the input
!                         atmosphere and channel.
!                         units:      none
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( out )
!
!       tau_ad:           adjoint of the layer to boundary transmittances
!                         for the input atmosphere and channel.
!                         **this argument is set to zero on output.**
!                         units:      none
!                         type:       real
!                         dimension:  k
!                         attributes: intent( in out )
!
!       layer_index:      index array array associating the input absorber
!                         layer amount to the bracketing absorber space levels.
!                         units:      none
!                         type:       integer
!                         dimension:  k x j
!                         attributes: intent( in )
!
!       channel_index:    channel index id. this is a unique index associated
!                         with a (supported) sensor channel.
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in )
!
!       direction:        direction identifier.
!                         if = 0, calculate layer->surface transmittances (i.e. down)
!                            = 1, calculate layer->space   transmittances (i.e. up)
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       absorber_ad:      adjoint of the profile level integrated absorber
!                         amount array.
!                         units:      varies with absorber.
!                         type:       real( fp_kind )
!                         dimension:  0:k x j
!                         attributes: intent( in out )
!
!       predictor_ad:     adjoint of the profile layer predictors.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  i x k
!                         attributes: intent( in out )
!
! optional output arguments:
!       none.
!
! calls:
!       none.
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       the input argument tau_ad is set to zero upon output.
!
! restrictions:
!       none.
!
! procedure:
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!s-
!------------------------------------------------------------------------------

  subroutine compute_transmittance_ad( &
                                       ! -- forward input
                                       absorber,      &   ! input, 0:k x j
                                       predictor,     &   ! input, i x k
                                       tau,           &   ! input, k

                                       ! -- adjoint input
                                       tau_ad,        &   ! in/output, k

                                       ! -- other input
                                       layer_index,   &   ! input, k x j
                                       channel_index, &   ! input, scalar
                                       direction,     &   ! input, scalar, 0==down, 1==up

                                       ! -- adjoint output
                                       absorber_ad,   &   ! in/output, 0:k x j
                                       predictor_ad   )   ! in/output, i x k



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( 0:, : ), intent( in )     :: absorber       ! 0:k x j
    real( fp_kind ), dimension( :, : ),  intent( in )     :: predictor      ! i x k
    real( fp_kind ), dimension( : ),     intent( in )     :: tau            ! k

    ! -- adjoint input
    real( fp_kind ), dimension( : ),     intent( in out ) :: tau_ad         ! k

    ! -- other input
    integer,         dimension( :, : ),  intent( in )     :: layer_index    ! k x j
    integer,                             intent( in )     :: channel_index  ! scalar
    integer,                             intent( in )     :: direction      ! scalar, 0==down, 1==up

    ! -- adjoint output
    real( fp_kind ), dimension( 0:, : ), intent( in out ) :: absorber_ad    ! 0:k x j
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: predictor_ad   ! i x k


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_transmittance_ad'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: l
    integer :: k, k1, k2, dk, n_layers
    integer :: j, n_absorbers
    integer :: i, ip, n_predictors

    real( fp_kind ) :: ave_absorber,    d_absorber
    real( fp_kind ) :: ave_absorber_ad, d_absorber_ad

    real( fp_kind ) :: d_absorber_space

    real( fp_kind ) :: b1o, b2o
    real( fp_kind ) :: b1,  b2
    real( fp_kind ) :: gradient,    b
    real( fp_kind ) :: gradient_ad, b_ad

    real( fp_kind ) :: absorption_coefficient
    real( fp_kind ) :: absorption_coefficient_ad

    real( fp_kind ) :: total_od_ad

    real( fp_kind ), dimension( size( tau ) ) :: optical_depth_ad


    ! ----------
    ! intrinsics
    ! ----------

    intrinsic size



    !#--------------------------------------------------------------------------#
    !#                   -- determine array dimensions --                       #
    !#--------------------------------------------------------------------------#

    ! -- number of atmospheric layers and absorbers. the "-1"
    ! -- for the layer assign is because absorber is based on
    ! -- levels.
    n_layers    = size( absorber, dim = 1 ) - 1
    n_absorbers = size( absorber, dim = 2 )

    ! -- number of predictors *to use*. currently this is
    ! -- pretty much hardwired to 5 but it could change
    ! -- in the future. this value is explicitly
    ! -- checked as the tau_coefficients array contained
    ! -- in the transmittance_coefficients module is
    ! -- allocated dynamically when initialising the 
    ! -- radiative transfer model.
 
    ! -- the "-1" is used because the first element of the
    ! -- predictor index array is used to determine if there
    ! -- is *any* absorption for a particular absorber species.
    n_predictors = size( tau_coefficients, dim = 1 ) - 1



    !#--------------------------------------------------------------------------#
    !#                       -- some initialisations --                         #
    !#--------------------------------------------------------------------------#

    ! ----------------------------------------
    ! assign the channel index to a short name
    ! ----------------------------------------

    l = channel_index


    ! -------------------------------------
    ! initilise the local adjoint variables
    ! -------------------------------------

    optical_depth_ad( : ) = zero

    total_od_ad = zero
    gradient_ad = zero



    !#--------------------------------------------------------------------------#
    !#       -- calculate the layer->boundary transmittance adjoint --          #
    !#                                                                          #
    !# this step involves another loop over layers. one *could* reverse the     #
    !# order of the j absorber loop and k layer od loop and calculate the tau   #
    !# values outside the absorber loop. however, this would involve an if test #
    !# for every layer even if predictor_index( 0, l ) == 0. i figured it is    #
    !# better to skip the layer loop altogether if there is no absorption for   #
    !# a particular absorber, despite the need for an extra layer loop here to  #
    !# calculate the transmittance.                                             #
    !#--------------------------------------------------------------------------#

    ! --------------------------------------
    ! determine loop indices. same direction
    ! as forward and tangent-linear models
    ! --------------------------------------

    if ( direction == up ) then
      ! -- transmittance going up, e.g. upwelling. layer 1 => toa
      k1 = 1
      k2 = n_layers
      dk = 1
    else
      ! -- transmittance going down, e.g. downwelling flux, solar. layer n_layers => sfc
      k1 =  n_layers
      k2 =  1
      dk = -1
    end if


    ! -----------------------------------------------
    ! loop over layers for transmittance calculation.
    ! note the loop index reversal.
    ! -----------------------------------------------

    k_layer_tau_loop: do k = k2, k1, -dk

      total_od_ad = total_od_ad - ( tau( k ) * tau_ad( k ) )
      tau_ad( k ) = zero

      optical_depth_ad( k ) = optical_depth_ad( k ) + total_od_ad
      ! note: no total_od_ad = zero here because
      !       total_od_tl = total_od_tl + (....)

    end do k_layer_tau_loop


    !#--------------------------------------------------------------------------#
    !#                        -- loop over absorbers --                         #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: do j = 1, n_absorbers


      ! -----------------------------------------
      ! check if there is any absorption for this
      ! absorber/channel combination.
      !
      ! this check is the reason why all channels
      ! cannot be processed at once and why the
      ! layer loop is within the absorber loop.
      ! -----------------------------------------

      if ( predictor_index( 0, l, j ) == 0 ) cycle j_absorber_loop



      !#------------------------------------------------------------------------#
      !#                        -- loop over layers --                          #
      !#------------------------------------------------------------------------#

      k_layer_od_loop: do k = n_layers, 1, -1


        ! -----------------------------------
        ! calculate the current layer average
        ! absorber amounts and differences
        ! -----------------------------------

        ave_absorber = point_5 * ( absorber( k, j ) + absorber( k-1, j ) )
        d_absorber   = absorber( k, j ) - absorber( k-1, j )


        ! -----------------------------------
        ! assign absorber space layer indices
        ! -----------------------------------

        k2 = layer_index( k, j )
        k1 = k2 - 1



        !#----------------------------------------------------------------------#
        !#           -- here repeat the forward calculation of the   --         #
        !#           -- absorption coefficient for the current layer --         #
        !#----------------------------------------------------------------------#

        ! ---------------------------------------------------------
        ! to linearly interpolate the tau_coeffs to the actual user
        ! space absorber amount, need the gradient across the layer
        ! ---------------------------------------------------------

        d_absorber_space = absorber_space_levels( k2, j ) - absorber_space_levels( k1, j )

        gradient = ( absorber_space_levels( k2, j ) - ave_absorber ) / &
        !          -------------------------------------------------
                                   d_absorber_space


        ! --------------------------------
        ! calculate absorption coeficients
        ! --------------------------------

        ! -- offset term (save these values for later)
        b1o = tau_coefficients( 0, k1, l, j )
        b2o = tau_coefficients( 0, k2, l, j )

        absorption_coefficient = b2o + ( gradient * ( b1o - b2o ) )


        i_predictor_loop_fwd: do i = 1, n_predictors

          ! -- current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )

          b = b2 + ( gradient * ( b1 - b2 ) )

          !  - sum current predictor's contribution
          ip = predictor_index( i, l, j )
          absorption_coefficient = absorption_coefficient + ( b * predictor( ip, k ) ) 

        end do i_predictor_loop_fwd


        ! -- **** it would be nice to not need this at all! ****
        absorption_coefficient = max( absorption_coefficient, zero )



        !#----------------------------------------------------------------------#
        !#                  -- begin adjoint calculations --                    #
        !#----------------------------------------------------------------------#

        ! ---------------------------------------------------------------
        ! adjoints of the optical depth.
        !
        ! these quantities are local to the k_layer_od_loop
        ! and are equal to zero at this point so a straight
        ! initialisation is used, i.e. there is no
        !   d_adbsorber_ad            = d_absorber_ad + (...)
        !   absorption_coefficient_ad = absorption_coefficient_ad + (...)
        ! this also eliminates the need to zero out the two
        ! quanitities later in the loop once they no longer
        ! have an impact on the gradient vector result.
        !
        ! also not that there is no
        !   optical_depth_ad( k ) = zero
        ! because
        !   optical_depth_tl( k ) = optical_depth_tl( k ) + (....)
        ! ---------------------------------------------------------------

        d_absorber_ad = absorption_coefficient * optical_depth_ad( k )   ! .... (1)
        absorption_coefficient_ad = d_absorber * optical_depth_ad( k )


        ! ---------------------------------
        ! adjoint of absorption coefficient
        ! ---------------------------------

        ! -- loop over predictors
        ! -- note that b_ad is local to this loop only
        ! -- hence rather than
        ! --   b_ad = b_ad + (....)
        ! -- and, at end of loop,
        ! --   b_ad = zero
        ! -- it is simply initialised each iteration
        i_predictor_loop_ad: do i = n_predictors, 1, -1

          ! -- current predictor interpolated coefficient
          b1 = tau_coefficients( i, k1, l, j )
          b2 = tau_coefficients( i, k2, l, j )

          b = b2 + ( gradient * ( b1 - b2 ) )

          ! -- adjoints of the absorption coefficient
          ip = predictor_index( i, l, j )

          predictor_ad( ip, k ) = predictor_ad( ip, k ) + ( b * absorption_coefficient_ad )
          b_ad = predictor( ip, k ) * absorption_coefficient_ad
          ! note: no absorption_coefficient_ad = zero here because
          !       absorption_coefficient_tl = absorption_coefficient_tl + (....)

          ! -- coefficient adjoint
          gradient_ad = gradient_ad + ( ( b1 - b2 ) * b_ad )

        end do i_predictor_loop_ad


        ! -- offset term.
        ! -- note that absorption_coefficient_ad is not set to zero
        ! -- after this statement due to its initialisation each layer.
        gradient_ad = gradient_ad + ( ( b1o - b2o ) * absorption_coefficient_ad )


        ! -------------------------------------------------
        ! adjoint of coefficient layer gradient.
        !
        ! the ave_absorber_ad quantity is also local to the
        ! k_layer_od_loop so a straight initialisation can
        ! be used rather than,
        !   ave_absorber_ad = ave_absorber_ad - (....)
        ! -------------------------------------------------

        ave_absorber_ad = -( gradient_ad / d_absorber_space )   ! .... (2)
        gradient_ad = zero


        ! ---------------------------------------------------
        ! adjoints of the current layer average
        ! absorber amount and difference.
        !
        ! neither d_absorber_ad nor ave_absorber_ad need
        ! to be set to zero after this as they are explicitly
        ! reassigned each layer iteration at (1) and (2) above
        ! respectively.
        ! ---------------------------------------------------

        absorber_ad( k-1, j ) = absorber_ad( k-1, j ) - d_absorber_ad + ( point_5 * ave_absorber_ad )
        absorber_ad( k,   j ) = absorber_ad( k,   j ) + d_absorber_ad + ( point_5 * ave_absorber_ad )

      end do k_layer_od_loop

    end do j_absorber_loop

  end subroutine compute_transmittance_ad

end module transmittance


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
! revision 1.8  2001/08/16 17:19:11  paulv
! - updated documentation
!
! revision 1.7  2001/08/01 17:04:00  paulv
! - the absorber space levels are no longer calculated during model
!   initialisation, but are precalculated and stored in the transmittance
!   coefficient data file. this means that,
!     use absorber_space, only : absorber_space_levels
!   was deleted as the absorber space level array is now available from
!   the transmittance_coefficients module.
!
! revision 1.6  2001/07/12 18:38:28  paulv
! - use of absorber_space module now includes an only clause so that the
!   absorber_space_levels is all that is available.
! - direction specification changed from
!     if ( direction == 0 ) then
!       ...do downwelling stuff...
!     else
!       ...do upwelling stuff...
!     end if
!   to
!     if ( direction == up ) then
!       ...do upwelling stuff...
!     else
!       ...do downwelling stuff...
!     end if
!   since the upwelling case is required for every call, but the downwelling
!   may not be. also, the parameter up is now used in the if rather than
!   an actual number (0 in this case).
! - changed
!      od_tolerance = abs( alog( tolerance ) )
!   to
!      od_tolerance = abs( log( tolerance ) )
! - corrected bug in the forward calculation of the absorption coefficient
!   in transmittance_ad. the offset coefficients are defined as
!     b1o = tau_coefficients( 0, k1, l, j )
!     b2o = tau_coefficients( 0, k2, l, j )
!   and the offset term was initialised as
!     absorption_coefficient = b2 + ( gradient * ( b1o - b2o ) )
!   instead of
!     absorption_coefficient = b2o + ( gradient * ( b1o - b2o ) )
!   where in the former, b2 was specified rather than b2o
!
! revision 1.5  2001/05/29 18:00:08  paulv
! - added adjoint form of the transmittance calculation.
! - removed the find_absorber_space_layer  routine. now resides in the
!   absorber_profile module. the absorber space bracket layer indices are
!   now passed as arguments from the calling routine.
! - the predictor indices and transmittance coefficients are no longer passed
!   as arguments but read from the transmittance_coefficients module.
!
! revision 1.4  2000/11/14 18:42:32  paulv
! - merged branch incorporating tangent-linear code into main truck.
!   optical depth debug code still present.
!
! revision 1.3.1.2  2000/11/14 18:34:56  paulv
! - finished adding tangent-linear code. optical depth debug code still
!   present - output sent to unit numbers 51 and 61.
!
! revision 1.3.1.1  2000/11/09 20:49:35  paulv
! - adding tangent linear forms of the optical depth and transmittance
!   computation. in progress and incomplete.
! - removed code that finds the absorber space bracket layers into its own
!   subroutine. both the forward and tangent linear routines use the same
!   search method.
!
! revision 1.3  2000/08/31 19:36:33  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.2  2000/08/24 16:55:33  paulv
! - added optional no_standard input argument to the compute_transmittance
!   subprogram to prevent the (angle independent) standard predictors from being
!   recalculated when only the path angle has changed in the calling procedure.
! - the  profile data integration has been removed from the compute_transmittance
!   subprogram and is now performed outside of this module in the absorber_profile
!   module. this has a number of consequences:
!   o the view_angle input argument was removed and replaced with the path-angle
!     scaled absorber_amounts argument.
!   o the interface pressure and ozone profile data are no longer required
!     and have been removed from the input argument list.
! - all profile integration and predictor calculation code has been removed
!   from the compute_transmittance subprogram.
! - updated module and subprogram documentation.
!
! revision 1.1  2000/08/22 15:57:26  paulv
! initial checkin.
!
!
!
!
