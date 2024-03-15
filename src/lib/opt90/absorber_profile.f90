!------------------------------------------------------------------------------
!m+
! name:
!       absorber_profile
!
! purpose:
!       module containing routines to compute and assemble the integrated
!       absorber profiles.
!
! category:
!       ncep rtm
!
! calling sequence:
!       use absorber_profile
!
! outputs:
!       none.
!
! modules:
!       parameters:  module containing parameter definitions for the
!                    rt model.
!
! contains:
!       compute_absorber_amount:     public subroutine to compute the integrated
!                                    absorber profiles. currently the absorbers
!                                    are:
!                                      - water vapor
!                                      - dry/fixed gases
!                                      - ozone
!
!       compute_absorber_amount_tl:  public subroutine to compute the tangent-
!                                    linear form of the integrated absorber 
!                                    profiles.
!
!       compute_absorber_amount_ad:  public subroutine to compute the adjoint of
!                                    the integrated absorber profiles.
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
!       written by:     paul van delst, cimss@noaa/ncep 01-aug-2000
!                       pvandelst@ncep.noaa.gov
!
!       adapted from code written by: thomas j.kleespies
!                                     noaa/nesdis/ora
!                                     tkleespies@nesdis.noaa.gov
!
!  copyright (c) 2000 thomas kleespies, paul van delst
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

module absorber_profile


  ! ---------------------
  ! module use statements
  ! ---------------------

  use type_kinds, only : fp_kind
  use parameters
  use transmittance_coefficients, only : alpha, &
                                         absorber_space_levels


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------------
  ! default visibility
  ! ------------------

  private


  ! ----------------------------------
  ! explicit visibility of subprograms
  ! ----------------------------------

  public :: compute_absorber_amount
  public :: compute_absorber_amount_tl
  public :: compute_absorber_amount_ad
  public :: find_absorber_layer_index


contains



!--------------------------------------------------------------------------------
!s+
! name:
!       compute_absorber_amount
!
! purpose:
!       public subroutine to compute the integrated profiles for all the
!       absorbers. currently the number of absorbers are:
!         - water vapor
!         - dry/fixed gases (pressure == absorber amount)
!         - ozone
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_absorber_amount( pressure,    &  ! input,  k
!                                     water_vapor, &  ! input,  k
!                                     ozone,       &  ! input,  k
!                                     absorber     )  ! output, 0:k x j
!
! input arguments:
!       pressure:     profile level pressure array.
!                     units:      hpa
!                     type:       real( fp_kind )
!                     dimension:  k
!                     attributes: intent( in )
!
!       water_vapor:  profile layer water vapor mixing ratio array.
!                     units:      g/kg
!                     type:       real( fp_kind )
!                     dimension:  k
!                     attributes: intent( in )
!
!       ozone:        profile layer ozone mixing ratio array.
!                     units:      ppmv
!                     type:       real( fp_kind )
!                     dimension:  k
!                     attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       absorber:     profile level integrated absorber amount array.
!                     units:      varies with absorber
!                     type:       real( fp_kind )
!                     dimension:  0:k x j
!                     attributes: intent( out )
!
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
!       the function calculates and accumulates the integrated path length
!       through k atmospheric layers,
!
!                     __ k
!                 1  \
!         u(k) = ---  >  q(i).dp(i)
!                 g  /__
!                       i=1
!
!       with the units of u dependent on those of the input mixing ratio, q.
!
!       the exception is the dry gas absorber amount which is represented
!       by the pressure (which is already an integrated quantity).
!
!       also,
!
!         u(0) = 0.0          for water vapor and ozone
!              = toa_pressure for dry/fixed gases
!
!       the routine loop over layers, k, with each absorber amount calculated
!       independently (like looping over absorber, j). this is opposite
!       to what is recommended (since the arrays are dimensioned [k,j]) but 
!       later use of the arrays require the [k,j] ordering where the absorber
!       loop over j is the outside loop (see compute_transmittance()).
!s-
!--------------------------------------------------------------------------------

  subroutine compute_absorber_amount( pressure,    &  ! input,  k
                                      water_vapor, &  ! input,  k
                                      ozone,       &  ! input,  k
                                      absorber     )  ! output, 0:k x j


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ),     intent( in )  :: pressure     ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: water_vapor  ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: ozone        ! input,  k

    real( fp_kind ), dimension( 0:, : ), intent( out ) :: absorber     ! output, 0:k x j


    ! ---------------
    ! local variables
    ! ---------------

    integer :: k
    real( fp_kind ) :: dp



    !#--------------------------------------------------------------------------#
    !#              -- initialise 0'th level absorber amounts --                #
    !#                                                                          #
    !# this is done so that layer differences and averages can be calculated    #
    !# simply in the predictor and transmittance routines.                      #
    !#--------------------------------------------------------------------------#

    absorber( 0, : ) = zero



    !#--------------------------------------------------------------------------#
    !#             -- assemble the nadir level absorber profiles --             #
    !#--------------------------------------------------------------------------#

    ! ---------------------------------------------
    ! toa layer. pressure is dimensioned as k hence
    ! the top input layer is treated separately
    ! ---------------------------------------------

    dp = pressure( 1 )

    ! -- integrated absorber amount
    absorber( 1, 1 ) = reciprocal_gravity * dp * water_vapor( 1 )
    absorber( 1, 2 ) = dp
    absorber( 1, 3 ) = reciprocal_gravity * dp * ozone( 1 )


    ! --------------------------------
    ! loop over layers, toa - 1 -> sfc
    ! --------------------------------

    k_layer_loop: do k = 2, size( pressure )

      ! -- layer pressure difference
      dp = pressure( k ) - pressure( k-1 )

      ! -- integrated absorber amounts
      absorber( k, 1 ) = absorber( k-1, 1 ) + ( reciprocal_gravity * dp * water_vapor( k ) )
      absorber( k, 2 ) = pressure( k )
      absorber( k, 3 ) = absorber( k-1, 3 ) + ( reciprocal_gravity * dp * ozone( k ) )

    end do k_layer_loop

  end subroutine compute_absorber_amount





!--------------------------------------------------------------------------------
!s+
! name:
!       compute_absorber_amount_tl
!
! purpose:
!       public subroutine to compute the tangent linear form of the integrated
!       profiles for all the absorbers. currently the number of absorbers
!       are:
!         - water vapor
!         - dry/fixed gases (pressure == absorber amount)
!         - ozone
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_absorber_amount_tl( &
!                                        ! -- forward input
!                                        pressure,       &  ! input, k
!                                        water_vapor,    &  ! input, k
!                                        ozone,          &  ! input, k
!
!                                        ! -- tangent-linear input
!                                        pressure_tl,    &  ! input, k
!                                        water_vapor_tl, &  ! input, k
!                                        ozone_tl,       &  ! input, k
!
!                                        ! -- tangent-linear output
!                                        absorber_tl     )  ! output, 0:k x j
!
! input arguments:
!       pressure:        profile level pressure array.
!                        units:      hpa
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
!       water_vapor:     profile layer water vapor mixing ratio array.
!                        units:      g/kg
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
!       ozone:           profile layer ozone mixing ratio array.
!                        units:      ppmv
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
!       pressure_tl:     profile level tangent-linear pressure array,
!                        i.e. the pressure perturbation.
!                        units:      hpa
!                        type:       real( fp_kind )
!                        dimension:  k, number of levels - 1
!                        attributes: intent( in )
!
!       water_vapor_tl:  profile layer tangent-linear water vapor mixing
!                        ratio array, i.e. the water vapor mixing
!                        ratio perturbation.
!                        units:      g/kg
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
!       ozone_tl:        profile layer tangent-linear ozone mixing ratio
!                        array i.e. the ozone mixing ratio perturbation.
!                        units:      ppmv
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       absorber_tl:     profile level tangent-linear average integrated 
!                        absorber amount array.
!                        units:      varies with absorber.
!                        type:       real( fp_kind )
!                        dimension:  0:k x j
!                        attributes: intent( out )
!
!
! optional output arguments:
!       none.
!
! calls:
!       none
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
!       the function calculates and accumulates the tangent-linear of the
!       integrated path length through k atmospheric layers,
!                        __ k
!                    1  \
!         u_tl(k) = ---  >  [ q_tl(i).dp(i) + q(i).dp_tl(i) ]
!                    g  /__
!                          i=1
!
!       with the units of u dependent on those of the input mixing ratio, q.
!
!       the exception is the dry gas absorber amount which is represented
!       by the tangent linear of the pressure (i.e. a perturbation).
!
!       also,
!
!         u(0)    = 0.0          for water vapor and ozone
!                 = toa_pressure for dry/fixed gases
!
!         u_tl(0) = 0.0 for all absorbers
!
!       the routine loop over layers, k, with each absorber amount calculated
!       independently (like looping over absorber, j). this is opposite
!       to what is recommended (since the arrays are dimensioned [k,j]) but 
!       later use of the arrays require the [k,j] ordering where the absorber
!       loop over j is the outside loop (see compute_transmittance_tl()).
!s-
!--------------------------------------------------------------------------------

  subroutine compute_absorber_amount_tl( &
                                         ! -- forward input
                                         pressure,       &  ! input,  k
                                         water_vapor,    &  ! input,  k
                                         ozone,          &  ! input,  k

                                         ! -- tangent-linear input
                                         pressure_tl,    &  ! input,  k
                                         water_vapor_tl, &  ! input,  k
                                         ozone_tl,       &  ! input,  k

                                         ! -- tangent-linear output
                                         absorber_tl     )  ! output, 0:k x j


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( : ),     intent( in )  :: pressure        ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: water_vapor     ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: ozone           ! input,  k

    ! -- tangent-linear input
    real( fp_kind ), dimension( : ),     intent( in )  :: pressure_tl     ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: water_vapor_tl  ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: ozone_tl        ! input,  k

    ! -- tangent-linear output
    real( fp_kind ), dimension( 0:, : ), intent( out ) :: absorber_tl     ! output, 0:k x j



    ! ---------------
    ! local variables
    ! ---------------

    integer :: k

    real( fp_kind ) :: dp
    real( fp_kind ) :: dp_tl



    !#--------------------------------------------------------------------------#
    !#            -- initialise 0'th level tl absorber amounts --               #
    !#                                                                          #
    !# this is done so that layer differences and averages can be calculated    #
    !# simply in the predictor and transmittance routines.                      #
    !#--------------------------------------------------------------------------#

    absorber_tl( 0, : ) = zero



    !#--------------------------------------------------------------------------#
    !#        -- assemble the nadir tangent-linear absorber profiles --         #
    !#--------------------------------------------------------------------------#

    ! -----------------------------------------------
    ! toa layer. profile inputs are dimensioned as k
    ! hence the top input layer is treated separately
    ! -----------------------------------------------

    absorber_tl( 1, 1 ) = reciprocal_gravity * (( water_vapor_tl( 1 ) * pressure( 1 )    ) + &
                                                ( water_vapor( 1 )    * pressure_tl( 1 ) ))
    absorber_tl( 1, 2 ) = pressure_tl( 1 )
    absorber_tl( 1, 3 ) = reciprocal_gravity * (( ozone_tl( 1 ) * pressure( 1 )    ) + &
                                                ( ozone( 1 )    * pressure_tl( 1 ) ))


    ! --------------------------------
    ! loop over layers, toa - 1 -> sfc
    ! --------------------------------

    k_layer_loop: do k = 2, size( pressure )

      ! -- layer pressure differences
      dp    = pressure( k )    - pressure( k-1 )
      dp_tl = pressure_tl( k ) - pressure_tl( k-1 )

      ! -- integrated tl absorber amounts
      absorber_tl( k, 1 ) = absorber_tl( k-1, 1 ) + &
                            ( reciprocal_gravity * (( water_vapor_tl( k ) * dp    ) + &
                                                    ( water_vapor( k )    * dp_tl )) )
      absorber_tl( k, 2 ) = pressure_tl( k )

      absorber_tl( k, 3 ) = absorber_tl( k-1, 3 ) + &
                            ( reciprocal_gravity * (( ozone_tl( k ) * dp    ) + &
                                                    ( ozone( k )    * dp_tl )) )

    end do k_layer_loop

  end subroutine compute_absorber_amount_tl





!--------------------------------------------------------------------------------
!s+
! name:
!       compute_absorber_amount_ad
!
! purpose:
!       public subroutine to compute the adjoint of the integrated
!       profiles for all the absorbers. currently the number of absorbers
!       are:
!         - water vapor
!         - dry/fixed gases (pressure == absorber amount)
!         - ozone
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_absorber_amount_ad( &
!                                        ! -- forward input
!                                        pressure,       &  ! input, k
!                                        water_vapor,    &  ! input, k
!                                        ozone,          &  ! input, k
!
!                                        ! -- adjoint input
!                                        absorber_ad,    &  ! in/output, 0:k x j
!
!                                        ! -- adjoint output
!                                        pressure_ad,    &  ! in/output, k
!                                        water_vapor_ad, &  ! in/output, k
!                                        ozone_ad        )  ! in/output, k
!
! input arguments:
!       pressure:        profile level pressure array.
!                        units:      hpa
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
!       water_vapor:     profile layer water vapor mixing ratio array.
!                        units:      g/kg
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
!       ozone:           profile layer ozone mixing ratio array.
!                        units:      ppmv
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in )
!
!       absorber_ad:     profile level adjoint average integrated 
!                        absorber amount array.
!                        ** this argument is set to zero on output **.
!                        units:      varies with absorber
!                        type:       real( fp_kind )
!                        dimension:  0:k x j
!                        attributes: intent( out )
!
! optional input arguments:
!       none.
!
! output arguments:
!       pressure_ad:     profile layer adjoint pressure array.
!                        units:      hpa
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in out )
!
!       water_vapor_ad:  profile layer adjoint water vapor array.
!                        units:      g/kg
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in out )
!
!       ozone_ad:        profile layer adjoint ozone water vapor array.
!                        units:      ppmv
!                        type:       real( fp_kind )
!                        dimension:  k
!                        attributes: intent( in out )
!
!
! optional output arguments:
!       none.
!
! calls:
!       none
!
! externals:
!       none
!
! common blocks:
!       none.
!
! side effects:
!       the input argument absorber_ad is set to zero on output.
!
! restrictions:
!       none.
!
! procedure:
!       the function calculates the adjoints of the integrated path length
!       through k atmospheric layers. given the tangent-linear expression,
!
!                        __ k
!                    1  \
!         u_tl(k) = ---  >  [ q_tl(i).dp(i) + q(i).dp_tl(i) ]
!                    g  /__
!                          i=1
!
!       for the profile inputs (water vapor and ozone), the adjoints are
!       given by,
!
!                              1
!         dp_ad     = dp_ad + --- q(k).u_ad(k)
!                              g
!
!                                1
!         q_ad(k)   = q_ad(k) + --- dp.u_ad(k)
!                                g
!
!
!         u_ad(k-1) = u_ad(k-1) + u_ad(k)
!
!         u_ad(k)   = 0.0
!
!       looping over the layer index, k, from the surface (max.) to the
!       top of the atmosphere. for the dry absorber component, which uses
!       pressure as the absorber amount proxy, the tangent-linear expression
!       is simply,
!
!         u_tl(k) = p_tl(k)
!
!       the adjoint is thus,
!
!         p_ad(k) = p_ad(k) + u_ad(k)
!
!         u_ad(k) = 0.0
!
!       the actual execution of the above expressions in the code are done to
!       prevent recalculation of the the same quantities in individual loops.
!
!       typically, the adjoint quantities that correspond to inputs in the
!       forward model (pressure_ad, water_vapor_ad, and ozone_ad) are set to
!       0.0 before entering this routine, and the adjoint quantities 
!       corresponding to outputs in the forward model (absorber_ad) are set to
!       1.0 before entering this routine. this will return the gradient vector
!       of the output variable with respect to the input variable (partial
!       derivative.) e.g. if the forward problem is defined as,
!
!
!         u = f(q,p)
!
!       then the tl form would be,
!
!                 df          df
!         u_tl = ---- q_tl + ---- p_tl   (note: the df/dx is partial)
!                 dq          dp
!
!       if u_ad = 1.0 and q_ad,p_ad = 0.0 on input, then on output,
!
!                 df            df
!         q_ad = ---- , p_ad = ---- , and u_ad = 0.0
!                 dq            dp 
! 
!s-
!--------------------------------------------------------------------------------

  subroutine  compute_absorber_amount_ad( &
                                          ! -- forward input
                                          pressure,       &  ! input, k
                                          water_vapor,    &  ! input, k
                                          ozone,          &  ! input, k

                                          ! -- adjoint input
                                          absorber_ad,    &  ! in/output, 0:k x j

                                          ! -- adjoint output
                                          pressure_ad,    &  ! in/output, k
                                          water_vapor_ad, &  ! in/output, k
                                          ozone_ad        )  ! in/output, k



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( : ),     intent( in )     :: pressure        ! input, k
    real( fp_kind ), dimension( : ),     intent( in )     :: water_vapor     ! input, k
    real( fp_kind ), dimension( : ),     intent( in )     :: ozone           ! input, k

    ! -- adjoint input
    real( fp_kind ), dimension( 0:, : ), intent( in out ) :: absorber_ad     ! in/output, 0:k x j

    ! -- adjoint output
    real( fp_kind ), dimension( : ),     intent( in out ) :: pressure_ad     ! in/output, k
    real( fp_kind ), dimension( : ),     intent( in out ) :: water_vapor_ad  ! in/output, k
    real( fp_kind ), dimension( : ),     intent( in out ) :: ozone_ad        ! in/output, k


    ! ---------------
    ! local variables
    ! ---------------

    integer :: k

    real( fp_kind ) :: dp
    real( fp_kind ) :: dp_ad



    !#--------------------------------------------------------------------------#
    !#           -- assemble the nadir adjoint absorber profiles --             #
    !#--------------------------------------------------------------------------#

    ! --------------------------------
    ! loop over layers, sfc -> toa - 1
    ! --------------------------------

    k_layer_loop: do k = size( pressure ), 2, -1

      ! -- layer pressure differences
      dp = pressure( k ) - pressure( k-1 )

      ! -- ozone absorber adjoint
      ozone_ad( k )       = ozone_ad( k ) + &
                            ( reciprocal_gravity * dp * absorber_ad( k, 3 ) )

      ! -- pressure absorber adjoint
      pressure_ad( k )    = pressure_ad( k ) + absorber_ad( k, 2 )

      ! -- water vapor absorber adjoint
      water_vapor_ad( k ) = water_vapor_ad( k ) + &
                            ( reciprocal_gravity * dp * absorber_ad( k, 1 ) )

      ! -- layer pressure difference adjoint
      dp_ad = reciprocal_gravity * ( ( water_vapor( k ) * absorber_ad( k, 1 ) ) + &
                                     ( ozone( k )       * absorber_ad( k, 3 ) )   )
      pressure_ad( k )   = pressure_ad( k )   + dp_ad
      pressure_ad( k-1 ) = pressure_ad( k-1 ) - dp_ad
      ! note: no dp_ad = zero here as it is reinitialised every layer

      ! -- previous layer absorber amounts
      absorber_ad( k-1, 3 ) = absorber_ad( k-1, 3 ) + absorber_ad( k, 3 )
      absorber_ad( k,   3 ) = zero

      absorber_ad( k,   2 ) = zero

      absorber_ad( k-1, 1 ) = absorber_ad( k-1, 1 ) + absorber_ad( k, 1 )
      absorber_ad( k,   1 ) = zero

    end do k_layer_loop


    ! -----------------------------------------------
    ! toa layer. profile inputs are dimensioned as k
    ! hence the top input layer is treated separately
    ! -----------------------------------------------

    ! -- ozone
    ozone_ad( 1 )       = ozone_ad( 1 ) + &
                          ( reciprocal_gravity * pressure( 1 ) * absorber_ad( 1, 3 ) )

    ! -- pressure
    pressure_ad( 1 )    = pressure_ad( 1 ) + &
                          ( reciprocal_gravity * ( ( ozone( 1 )       * absorber_ad( 1, 3 ) ) + &
                                                   ( water_vapor( 1 ) * absorber_ad( 1, 1 ) )   ) ) + &
                          absorber_ad( 1, 2 )
    ! -- water vapor
    water_vapor_ad( 1 ) = water_vapor_ad( 1 ) + &
                          ( reciprocal_gravity * pressure( 1 ) * absorber_ad( 1, 1 ) )

    ! -- zero absorber amount adjoints
    absorber_ad( 0:1, : ) = zero

  end subroutine compute_absorber_amount_ad






  subroutine find_absorber_layer_index( absorber,   &  ! input,  0:k x j
                                        layer_index )  ! output, k x j

    real( fp_kind ), dimension( 0:, : ), intent( in )  :: absorber
    integer,         dimension(  :, : ), intent( out ) :: layer_index

    integer :: j,  n_absorbers
    integer :: k,  n_layers

    real( fp_kind ) :: ave_absorber


    intrinsic abs, &
              size



    !#--------------------------------------------------------------------------#
    !#                  -- determine the array dimension --                     #
    !#--------------------------------------------------------------------------#

    n_layers    = size( absorber, dim = 1 ) - 1
    n_absorbers = size( absorber, dim = 2 )



    !#--------------------------------------------------------------------------#
    !#                       -- loop over absorbers --                          #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: do j = 1, n_absorbers


      ! ------------------------------
      ! loop over the user layer space
      ! ------------------------------

      k_user_space_loop: do k = 1, n_layers


        ! ---------------------------------------
        ! calculate layer average absorber amount
        ! ---------------------------------------

        ave_absorber = point_5 * ( absorber( k, j ) + absorber( k-1, j ) )


        ! -------------------------------------------------------
        ! calculate the absorber space bracket layer, ka, for:
        !
        !   abs_spc(ka-1) < absorber < abs_spc(ka)
        !
        ! this is done using the exponential parameter, alpha,
        ! used in determining the absorber space layers to 
        ! begin with:
        !
        !                   exp( ka.alpha ) - 1
        !   a(ka) = a(1) . ---------------------
        !                    exp( alpha ) - 1
        !
        ! so given an absorber amount, a(k), the required ka
        ! value is found using:
        !
        !          1        (  a(k)                        )
        !   k = ------- . ln( ------( exp(alpha)-1 )  +  1 )
        !        alpha      (  a(1)                        )
        !
        ! the use of the ceiling intrinsic function provides the
        ! integer value equal to or greater than the argument.
        ! then use of max( min( ka, max_n_absorbers_layers ), 1 )
        ! is to ensure that the calculated value never exceeds the 
        ! maximum number of absorber space layers - which could
        ! happen if an input absorber profile contained amounts
        ! which exceeded the maximum absorber space amount - and
        ! is never less than 1 - which could happen if the input
        ! absorber profile has adjacent layers with zero amounts.
        ! -------------------------------------------------------

        layer_index( k, j ) = max( min( int( ceiling( ( one / alpha(j) ) * &
                                                 log( ( ave_absorber / absorber_space_levels(1,j) ) * &
                                                      ( exp( alpha(j) ) - one ) + one ) ) ), &
                                   max_n_absorber_layers ), 1 )


      end do k_user_space_loop

    end do j_absorber_loop

  end subroutine find_absorber_layer_index

end module absorber_profile


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
! revision 1.11  2001/09/25 15:51:29  paulv
! - changed the calculation of the bracketing absorber space layer in
!   sbroutine find_absorber_layer_index from
!     min( ka, max_n_absorbers_layers )
!   to
!     max( min( ka, max_n_absorbers_layers ), 1 )
!   so as to avoid the result ever being zero - which could happen before if
!   adjacent layers of the input absorber profile were zero.
!
! revision 1.10  2001/08/31 20:41:18  paulv
! - altered method of searching for bracketing absorber space layers in
!   find_absorber_layer_index. previosuly a trickle down search was performed.
!   now the actual corresponding layer is calculated using the exponential
!   factor used in generating the absorber space.
!
! revision 1.9  2001/08/16 16:30:38  paulv
! - updated documentation.
!
! revision 1.8  2001/08/01 16:36:34  paulv
! - removed use of module absorber_space and replaced it with
!     use transmittance_coefficients, only : absorber_space_levels
!   to reflect changes in code. the absorber space levels are no longer
!   calculated during model initialisation, but are precalculated and stored
!   in the transmittance coefficient data file.
!
! revision 1.7  2001/06/05 21:18:10  paulv
! - changed adjoint routine slightly to make adjoint calcs a bit clearer
!   when looking at the tangent-linear code.
! - corrected bug in toa layer pressure_ad calculation.
!
! revision 1.6  2001/05/29 17:32:51  paulv
! - some cosmetic changes
! - removed subtraction of the toa_pressure parameter from the dry absorber
!   calculation. this was causing the upper level channels to produce
!   spurious numbers in the forward calculation.
! - added the  find_absorber_layer_index routine. removed it from the forward_model
!   module. it seemed more appropriate in this one.
! - using pressure array data directly in first layer calcs rather than
!   dp variable.
!
! revision 1.5  2001/03/26 18:45:59  paulv
! - now use type_kinds module parameter fp_kind to set the floating point
!   data type.
! - module parameter reciprocal_gravity moved to parameters module.
! - only clause used in use parameters statement. only parameters available
!   in absorber_profile module are zero, toa_pressure, and reciprocal_gravity.
! - output absorber argument is now dimensioned as 0:k. this eliminates the
!   need for using an absorber_km1 variable in other routines that use the
!   absorber array variable where the layer loop always goes from 1 -> n_layers.
! - removed output arguments of ave_absorber and delta_absorber due to the
!   absorber dimension change to 0:k. calculating the average and layer
!   difference absorber amounts in other routines can now be done simply
!   by averaging or differencing abosrber(k) and absorber(k-1) even for
!   layer #1.
! - integration of absorber amount for the toa layer is done outside of the
!   layer loop. this avoids the need for a pressure_km1 variable since
!   pressure is dimensioned 1:k.
! - layer loop, thus, goes from 2 -> n_layers.
! - changed order of argument list in compute_predictors_tl and its
!   associated routines. all forward arguments are listed followed by
!   the tangent-linear arguments rather than interspersing them as before.
! - added adjoint routine compute_absorber_amount_ad.
!
! revision 1.4  2001/03/26 18:30:54  paulv
! - removed integrate_absorber_profile and integrate_absorber_profile_tl
!   functions. integration is now done in-line in the main routines.
!
! revision 1.3  2000/11/09 20:29:40  paulv
! - added tangent linear form of the absorber profile routines.
!
! revision 1.2  2000/08/31 19:36:31  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.1  2000/08/24 13:11:27  paulv
! initial checkin.
!
!
!
!
