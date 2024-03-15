!------------------------------------------------------------------------------
!m+
! name:
!       predictors
!
! purpose:
!       rt model predictor module
!
! category:
!       ncep rtm
!
! calling sequence:
!       use predictors
!
! outputs:
!       none.
!
! modules:
!       type_kinds:     module containing data type kind definitions.
!
!       parameters:     module containing parameter definitions for the
!                       rt model.
!
!       error_handler:  module to define error codes and handle error
!                       conditions
!
! contains:
!       compute_predictors:             public subroutine to calculate all
!                                       the predictors. this routine calls the
!                                       private functions that follow.
!
!       compute_std_predictors:         private function to compute the standard,
!                                       absorber independent predictor set.
!
!       compute_int_predictors:         private function to compute the integrated,
!                                       absorber dependent predictor set.
!
!       compute_predictors_tl:          public subroutine to calculate all
!                                       tangent-linear predictors. this routine
!                                       calls the private functions that follow.
!
!       compute_std_predictors_tl:      private function to compute the standard,
!                                       absorber independent tangent-linear 
!                                       predictor set.
!
!       compute_int_predictors_tl:      private function to compute the integrated,
!                                       absorber dependent tangent-linear predictor
!                                       set.
!
!       compute_predictors_ad:          public subroutine to calculate the
!                                       predictor adjoints. this routine calls
!                                       the functions that follow.
!
!       compute_std_predictors_ad:      private function to compute the adjoint of
!                                       the standard, absorber independent predictor
!                                       set.
!
!       compute_int_predictors_ad:      private function to compute the adjoint of
!                                       integrated, absorber dependent predictor set.
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
!       written by:     paul van delst, cimss@noaa/ncep 11-jul-2000
!                       pvandelst@ncep.noaa.gov
!
!       adapted from code written by: thomas j.kleespies
!                                     noaa/nesdis/ora
!                                     tkleespies@nesdis.noaa.gov
!
!                                     and modified by:
!
!                                     john derber
!                                     noaa/ncep/emc
!                                     jderber@ncep.noaa.gov
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

module predictors


  ! ---------------------
  ! module use statements
  ! ---------------------

  use type_kinds, only : fp_kind
  use parameters
  use error_handler


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private
  public :: compute_predictors
  public :: compute_predictors_tl
  public :: compute_predictors_ad


contains



!--------------------------------------------------------------------------------
!s+
! name:
!       compute_predictors
!
! purpose:
!       public routine to calculate the forward transmittance model predictors.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_predictors ( pressure,    &  ! input,  k
!                                 temperature, &  ! input,  k
!                                 water_vapor, &  ! input,  k
!                                 absorber,    &  ! input,  0:k x j
!
!                                 predictor,   &  ! output, i x k
!
!                                 no_standard  )  ! optional input
!
! input arguments:
!       pressure:         profile layer average pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer average temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       water_vapor:      profile layer average water vapor mixing ratio array.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       absorber:         profile level integrated absorber amount array.
!                         units:      varies with absorber.
!                         type:       real( fp_kind )
!                         dimension:  0:k x j
!                         attributes: intent( in )
!
! optional input arguments:
!       no_standard:      if present, the standard predictors are not calculated.
!                         this prevents recalculation of the standard predictors
!                         is only the view angle has changed - which only affects
!                         the integrated predictors.
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in ), optional
!
! output arguments:
!       predictor:        profile layer predictors array.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  i x k
!                         attributes: intent( out ), target
!
! optional output arguments:
!       none.
!
! calls:
!       compute_std_predictors:   private function to compute the standard
!                                 (absorber independent) predictor set.
!
!       compute_int_predictors:   private function to compute the absorber
!                                 integrated predictor set.
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
!       the predictors used in the ncep transmittance model are organised in
!       the following manner:
!
!       --------------------------------------------
!       | 1 | 2 | 3 | ... | 9 | 10 | 11 | ... | 27 |
!       --------------------------------------------
!
!       \                    / \                   /
!        \                  /   \                 /
!         ------------------     -----------------
!                  |                      |
!                  v                      v
!
!              standard               integrated
!             predictors            predictors for
!                                   each absorber
!
!       pointers are used to reference the module scope predictor data array
!       before the calls are made to the standard and integrated predictor
!       calculation functions. this eliminates array copying.
!s-
!--------------------------------------------------------------------------------

  subroutine compute_predictors ( pressure,    &  ! input,  k
                                  temperature, &  ! input,  k
                                  water_vapor, &  ! input,  k
                                  absorber,    &  ! input,  0:k x j

                                  predictor,   &  ! output, i x k

                                  no_standard  )  ! optional input



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ),     intent( in )           :: pressure      ! k
    real( fp_kind ), dimension( : ),     intent( in )           :: temperature   ! k
    real( fp_kind ), dimension( : ),     intent( in )           :: water_vapor   ! k
    real( fp_kind ), dimension( 0:, : ), intent( in )           :: absorber      ! 0:k x j

    real( fp_kind ), dimension( :, : ),  intent( out ), target  :: predictor     ! i x k

    integer,                             intent( in ), optional :: no_standard


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_predictors'


    ! ---------------
    ! local variables
    ! ---------------

    character( 80 ) :: message

    integer :: i1, i2, j

    real( fp_kind ), pointer, dimension( :, : ) :: std_predictors   ! istd x k
    real( fp_kind ), pointer, dimension( :, : ) :: int_predictors   ! iint x k



    !#--------------------------------------------------------------------------#
    !#            -- calculate the standard predictors if needed --             #
    !#--------------------------------------------------------------------------#

    if ( .not. present( no_standard ) ) then

      ! -- alias the input predictor array
      std_predictors => predictor( 1:max_n_standard_predictors, : )

      call compute_std_predictors( pressure,      &
                                   temperature,   &
                                   water_vapor,   &
                                   std_predictors )

    end if


                                                   
    !#--------------------------------------------------------------------------#
    !#                -- calculate the integrated predictors --                 #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: do j = 1, size( absorber, dim = 2 )

      ! -- determine indices of current absorber predictors
      i1 = max_n_standard_predictors + ( ( j - 1 ) * max_n_integrated_predictors ) + 1
      i2 = i1 + max_n_integrated_predictors - 1

      ! -- alias the input predictor array
      int_predictors => predictor( i1:i2, : )

      ! -- compute the predictors for the current absorber
      call compute_int_predictors( pressure,          &
                                   temperature,       &
                                   absorber( 0:, j ), &
                                   int_predictors     )

    end do j_absorber_loop

  end subroutine compute_predictors



!--------------------------------------------------------------------------------
!p+
! name:
!       compute_std_predictors
!
! purpose:
!       private function to compute the standard, absorber independent
!       predictor set for the forward transmittance model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_std_predictors( pressure,    &
!                                    temperature, &
!                                    water_vapor, &
!                                    predictor    )
!
! input arguments:
!       pressure:         profile layer average pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer average temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       water_vapor:      profile layer average water vapor mixing ratio array.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
! optional input arguments:
!       none
!
! output arguments:
!       predictor:        array containing the calculated standard predictor set.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  istd x k
!                         attributes: intent( out )
!
! optional output arguments:
!       none
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
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!       the standard predictors are the following:
!
!         1) temperature, t
!         2) pressure, p
!         3) t^2
!         4) p^2
!         5) t.p
!         6) t^2.p
!         7) t.p^2
!         8) t^2.p^2
!         9) water vapor mixing ratio, q
!p-
!--------------------------------------------------------------------------------

  subroutine compute_std_predictors( p,        &  ! input,  k
                                     t,        &  ! input,  k
                                     w,        &  ! input,  k
                                     predictor )  ! output, istd x k



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ),     intent( in )  :: p          ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: t          ! k
    real( fp_kind ), dimension( : ),     intent( in )  :: w          ! k

    real( fp_kind ), dimension( :, : ),  intent( out ) :: predictor  ! istd x k


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_std_predictors'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: k

    real( fp_kind ) :: p2
    real( fp_kind ) :: t2



    !#--------------------------------------------------------------------------#
    !#                -- calculate the standard predictor set --                #
    !#--------------------------------------------------------------------------#

    k_layer_loop: do k = 1, size( p )


      ! -- precalculate the squared terms
      p2 = p( k ) * p( k )
      t2 = t( k ) * t( k )

      ! -- calculate and assign the absorber independent predictors
      predictor( 1, k ) = t( k )
      predictor( 2, k ) = p( k )
      predictor( 3, k ) = t2
      predictor( 4, k ) = p2
      predictor( 5, k ) = t( k ) * p( k )
      predictor( 6, k ) = t2     * p( k )
      predictor( 7, k ) = t( k ) * p2
      predictor( 8, k ) = t2     * p2
      predictor( 9, k ) = w( k )

    end do k_layer_loop

  end subroutine compute_std_predictors






!--------------------------------------------------------------------------------
!p+
! name:
!       compute_int_predictors
!
! purpose:
!       private function to compute the integrated, absorber dependent
!       predictor set for the forward transmittance model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_int_predictors( pressure,    &
!                                    temperature, &
!                                    absorber,    &
!                                    predictor    )
!
! input arguments:
!       pressure:         profile layer average pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer average temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       absorber:         profile level integrated absorber amount array.
!                         units:      varies with absorber
!                         type:       real( fp_kind )
!                         dimension:  0:k
!                         attributes: intent( in )
!
! optional input arguments:
!       none
!
! output arguments:
!       predictor:        array containing the calculated integrated predictor set
!                         for the passed absorber.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  iint x k
!                         attributes: intent( out )
!
! optional output arguments:
!       none
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
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!       the integrated predictors consist of six that are repeated for every
!       absorber:
!
!         1) t*
!         2) p*
!         3) t**
!         4) p**
!         5) t***
!         6) p***
!
!       the reference above details the predictor value for input level values. in
!       the equations below, the "x" quantities can be replaced with temperature, t, 
!       of pressure, p, to determine the required predictor:
!
!                           __ k
!                    1     \ 
!         x*(k) = --------  >  ( x(i) + x(i-1) ) ( a(i) - a(i-1) )
!                  c.a(k)  /__
!                             i=1
!
!                            __ k
!                     1     \
!        x**(k) = ---------  >  ( x(i)a(i) + x(i-1)a(i-1) ) ( a(i) - a(i-1) )
!                  c.a^2(k) /__
!                              i=1
!
!                            __ k
!                     1     \ 
!       x***(k) = ---------  >  ( x(i)a^2(i) + x(i-1)a^2(i-1) ) ( a(i) - a(i-1) )
!                  c.a^3(k) /__
!                              i=1
!
!       to accomodate input layer values (i.e. layer average values) from the
!       ncep gdas, the predictor formulations were modified to:
!
!                   __ k                        __k-1
!                  \                           \ 
!                   >  x(i)( a(i)-a(i-1) )      >  x(i)( a(i)-a(i-1) )
!                  /__                         /__ 
!                     i=1                         i=1
!         x*(k) = ------------------------- + -------------------------
!                            2.a(k)                    2.a(k-1)
!
!                   __ k
!                  \ 
!                   >  x(i)( a(i) + a(i-1) ) ( a(i) - a(i-1) )
!                  /__
!                     i=1
!        x**(k) = --------------------------------------------- +
!                                       2.a^2(k) 
!                   
!                   __k-1
!                  \ 
!                   >  x(i)( a(i) + a(i-1) ) ( a(i) - a(i-1) )
!                  /__
!                     i=1
!                 ---------------------------------------------
!                                       2.a^2(k-1)
!
!                      __ k
!                     \ 
!                 3 .  >  x(i)( a^2(i) + a^2(i-1) ) ( a(i) - a(i-1) )
!                     /__
!                        i=1
!       x***(k) = ---------------------------------------------------- +
!                                       4.a^3(k)                    
!
!                      __k-1
!                     \ 
!                 3 .  >  x(i)( a^2(i) + a^2(i-1) ) ( a(i) - a(i-1) )
!                     /__
!                        i=1
!                 ----------------------------------------------------
!                                       4.a^3(k-1)
!
!       thus the transmittance model coefficients calculated using the level
!       predictor formulation are used with the predictors constructed above 
!       with layer values.
!p-
!--------------------------------------------------------------------------------

  subroutine compute_int_predictors( pressure,    &  ! input,  k
                                     temperature, &  ! input,  k
                                     absorber,    &  ! input,  0:k
                                     predictor    )  ! output, iint x k



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ),    intent( in )  :: pressure     ! k
    real( fp_kind ), dimension( : ),    intent( in )  :: temperature  ! k
    real( fp_kind ), dimension( 0: ),   intent( in )  :: absorber     ! 0:k

    real( fp_kind ), dimension( :, : ), intent( out ) :: predictor    ! iint x k


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_int_predictors'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: i, k
    integer :: n_layers
    integer :: n_predictors

    real( fp_kind ) :: d_absorber
    real( fp_kind ) :: factor_1
    real( fp_kind ) :: factor_2
    real( fp_kind ) :: inverse_1
    real( fp_kind ) :: inverse_2
    real( fp_kind ) :: inverse_3
    real( fp_kind ) :: absorber_3


    ! -- square of the absorber amount. 0:k
    real( fp_kind ), dimension( 0:size( pressure ) ) :: absorber_2

    ! -- intermediate summation array. iint
    real( fp_kind ), dimension( size( predictor, dim=1 ) ) :: s

    ! -- level predictor, iint x 0:k
    real( fp_kind ), dimension( size( predictor, dim=1 ), 0:size( pressure ) ) :: x



    !#--------------------------------------------------------------------------#
    !#                  -- determine the number predictors --                   #
    !#--------------------------------------------------------------------------#

    n_predictors = size( predictor, dim = 1 )



    !#--------------------------------------------------------------------------#
    !#                         -- initialise values --                          #
    !#--------------------------------------------------------------------------#

    absorber_2( 0 ) = zero

    s( : )    = zero
    x( :, 0 ) = zero



    !#--------------------------------------------------------------------------#
    !#               -- calculate the integrated predictor set --               #
    !#--------------------------------------------------------------------------#

    k_layer_loop: do k = 1, size( pressure )


      ! -----------------------------------------
      ! calculate absorber multiplicative factors
      ! -----------------------------------------

      absorber_2( k ) = absorber( k ) * absorber( k )

      d_absorber = absorber( k ) - absorber( k-1 )                      ! for * terms
      factor_1   = ( absorber( k )   + absorber( k-1 )   ) * d_absorber ! for ** terms
      factor_2   = ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber ! for *** terms


      ! -------------------------------
      ! calculate the intermediate sums
      ! -------------------------------

      s( 1 ) = s( 1 ) + ( temperature( k ) * d_absorber )  ! t*
      s( 2 ) = s( 2 ) + ( pressure( k )    * d_absorber )  ! p*

      s( 3 ) = s( 3 ) + ( temperature( k ) * factor_1 )    ! t**
      s( 4 ) = s( 4 ) + ( pressure( k )    * factor_1 )    ! p**

      s( 5 ) = s( 5 ) + ( temperature( k ) * factor_2 )    ! t***
      s( 6 ) = s( 6 ) + ( pressure( k )    * factor_2 )    ! p***


      ! -------------------------------------------------------
      ! calculate the normalising factors for the integrated
      ! predictors. note that the checks below, the if tests to
      ! determine if the absorber products are represenatble
      ! are to minimise the number of calcs. i.e if inverse_1
      ! is toast because absorber(k) is too small there's no
      ! need to check any further.
      ! -------------------------------------------------------

      ! -- is inverse_1 representable?
      inverse_1_check: if ( absorber( k ) > tolerance ) then

        inverse_1 = one / absorber( k )

        ! -- is inverse_2 representable
        inverse_2_check: if ( absorber_2( k ) > tolerance ) then

          inverse_2  = inverse_1 * inverse_1
          absorber_3 = absorber( k ) * absorber_2( k )
         
          ! -- is inverse_3 representable
          inverse_3_check: if ( absorber_3 > tolerance ) then

            inverse_3 = inverse_2 * inverse_1

          else

            inverse_3 = zero

          end if inverse_3_check

        else

          inverse_2 = zero
          inverse_3 = zero

        end if inverse_2_check

      else

        inverse_1 = zero
        inverse_2 = zero
        inverse_3 = zero

      end if inverse_1_check


      ! ---------------------------------------------
      ! scale and normalise the integrated predictors
      ! ---------------------------------------------

      x( 1, k ) = point_5  * s( 1 ) * inverse_1  ! t*
      x( 2, k ) = point_5  * s( 2 ) * inverse_1  ! p*

      x( 3, k ) = point_5  * s( 3 ) * inverse_2  ! t**
      x( 4, k ) = point_5  * s( 4 ) * inverse_2  ! p**

      x( 5, k ) = point_75 * s( 5 ) * inverse_3  ! t***
      x( 6, k ) = point_75 * s( 6 ) * inverse_3  ! p***


      ! ----------------------------
      ! sum predictors across layers
      ! ----------------------------

      do i = 1, n_predictors

        predictor( i, k ) = x( i, k ) + x( i, k-1 )

      end do

    end do k_layer_loop

  end subroutine compute_int_predictors






!--------------------------------------------------------------------------------
!s+
! name:
!       compute_predictors_tl
!
! purpose:
!       public routine to calculate the tangent-linear transmittance model
!       predictors.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_predictors_tl ( pressure,       &  ! input,  k
!                                    temperature,    &  ! input,  k
!                                    water_vapor,    &  ! input,  k
!                                    absorber,       &  ! input,  0:k x j
!
!                                    pressure_tl,    &  ! input,  k
!                                    temperature_tl, &  ! input,  k
!                                    water_vapor_tl, &  ! input,  k
!                                    absorber_tl,    &  ! input,  0:k x j
!
!                                    predictor_tl,   &  ! output, i x k
!
!                                    no_standard     )  ! optional input
!
! input arguments:
!       pressure:            profile layer average pressure array.
!                            units:      hpa
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
!       temperature:         profile layer average temperature array.
!                            units:      kelvin
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
!       water_vapor:         profile layer average water vapor mixing ratio array.
!                            units:      g/kg
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
!       absorber:            profile level integrated absorber amount array.
!                            units:      varies with absorber
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
!       pressure_tl:         profile tangent-linear pressure array, i.e. the
!                            pressure perturbation.
!                            units:      hpa
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
!       temperature_tl:      profile tangent-linear layer average temperature array
!                            i.e. the temperature perturbation.
!                            units:      kelvin
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
!       water_vapor_tl:      profile tangent-linear water vapor mixing
!                            ratio array, i.e. the water vapor mixing
!                            ratio perturbation.
!                            units:      g/kg
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
!       absorber_tl:         profile level integrated tangent-linear 
!                            absorber amount array.
!                            units:      varies with absorber
!                            type:       real( fp_kind )
!                            dimension:  k
!                            attributes: intent( in )
!
! optional input arguments:
!       no_standard:         if present, the standard predictors are not calculated.
!                            this prevents recalculation of the standard predictors
!                            is only the view angle has changed - which only affects
!                            the integrated predictors.
!                            units:      none
!                            type:       integer
!                            dimension:  scalar
!                            attributes: intent( in ), optional
!
! output arguments:
!       predictor_tl:        profile layer tangent-linear predictors array.
!                            units:      varies with predictor type.
!                            type:       real( fp_kind )
!                            dimension:  iint x k
!                            attributes: intent( out ), target
!
!
! optional output arguments:
!       none.
!
! calls:
!       compute_std_predictors_tl:   private function to compute the tangent-
!                                    linear form of the standard (absorber
!                                    independent) predictor set.
!
!       compute_int_predictors_tl:   private function to compute the tangent-
!                                    linear form of the absorber integrated
!                                    predictor set.
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
!       the predictors used in the ncep transmittance model are organised in
!       the following manner:
!
!       --------------------------------------------
!       | 1 | 2 | 3 | ... | 9 | 10 | 11 | ... | 27 |
!       --------------------------------------------
!
!       \                    / \                   /
!        \                  /   \                 /
!         ------------------     -----------------
!                  |                      |
!                  v                      v
!
!              standard               integrated
!             predictors            predictors for
!                                   each absorber
!
!       pointers are used to reference the module scope predictor data array
!       before the calls are made to the standard and integrated predictor
!       calculation functions. this eliminates array copying.
!s-
!--------------------------------------------------------------------------------

  subroutine compute_predictors_tl ( pressure,       &  ! input,  k
                                     temperature,    &  ! input,  k
                                     water_vapor,    &  ! input,  k
                                     absorber,       &  ! input,  0:k x j

                                     pressure_tl,    &  ! input,  k
                                     temperature_tl, &  ! input,  k
                                     water_vapor_tl, &  ! input,  k
                                     absorber_tl,    &  ! input,  0:k x j

                                     predictor_tl,   &  ! output, i x k

                                     no_standard     )  ! optional input



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ),     intent( in )           :: pressure        ! k
    real( fp_kind ), dimension( : ),     intent( in )           :: temperature     ! k
    real( fp_kind ), dimension( : ),     intent( in )           :: water_vapor     ! k
    real( fp_kind ), dimension( 0:, : ), intent( in )           :: absorber        ! 0:k x j

    real( fp_kind ), dimension( : ),     intent( in )           :: pressure_tl     ! k
    real( fp_kind ), dimension( : ),     intent( in )           :: temperature_tl  ! k
    real( fp_kind ), dimension( : ),     intent( in )           :: water_vapor_tl  ! k
    real( fp_kind ), dimension( 0:, : ), intent( in )           :: absorber_tl     ! 0:k x j

    real( fp_kind ), dimension( :, : ),  intent( out ), target  :: predictor_tl    ! i x k

    integer,                             intent( in ), optional :: no_standard


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_predictors_tl'


    ! ---------------
    ! local variables
    ! ---------------

    character( 80 ) :: message

    integer :: i1, i2, j

    real( fp_kind ), pointer, dimension( :, : ) :: std_predictors_tl   ! istd x k
    real( fp_kind ), pointer, dimension( :, : ) :: int_predictors_tl   ! iint x k




    !#--------------------------------------------------------------------------#
    !#     -- calculate the tangent-linear standard predictors if needed --     #
    !#--------------------------------------------------------------------------#

    if ( .not. present( no_standard ) ) then

      ! -- alias the input predictor array
      std_predictors_tl => predictor_tl( 1:max_n_standard_predictors, : )

      call compute_std_predictors_tl( pressure,         &
                                      temperature,      &
                                      water_vapor,      &

                                      pressure_tl,      &
                                      temperature_tl,   &
                                      water_vapor_tl,   &

                                      std_predictors_tl )

    end if


                                                   
    !#--------------------------------------------------------------------------#
    !#                -- calculate the integrated predictors --                 #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: do j = 1, size( absorber_tl, dim = 2 )

      ! -- determine indices of current absorber predictors
      i1 = max_n_standard_predictors + ( ( j - 1 ) * max_n_integrated_predictors ) + 1
      i2 = i1 + max_n_integrated_predictors - 1

      ! -- alias the input predictor array
      int_predictors_tl => predictor_tl( i1:i2, : )

      ! -- calculate tangent-linear predictors for current absorber
      call compute_int_predictors_tl( pressure,             &
                                      temperature,          &
                                      absorber( 0:, j ),    &

                                      pressure_tl,          &
                                      temperature_tl,       &
                                      absorber_tl( 0:, j ), &

                                      int_predictors_tl     )

    end do j_absorber_loop

  end subroutine compute_predictors_tl





!--------------------------------------------------------------------------------
!p+
! name:
!       compute_std_predictors_tl
!
! purpose:
!       private function to compute the standard, absorber independent
!       predictor set for the tangent-linear transmittance model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_std_predictors_tl( pressure,       &  ! input,  k
!                                       temperature,    &  ! input,  k
!                                       water_vapor,    &  ! input,  k
!
!                                       pressure_tl,    &  ! input,  k
!                                       temperature_tl, &  ! input,  k
!                                       water_vapor_tl, &  ! input,  k
!
!                                       predictors_tl   )  ! output, istd x k
!
! input arguments:
!       pressure:         profile layer average pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer average temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       water_vapor:      profile layer average water vapor mixing ratio array.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       pressure_tl:      profile layer tangent-linear pressure array, i.e. the
!                         pressure perturbation.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature_tl:   profile layer tangent-linear layer average temperature array
!                         i.e. the temperature perturbation.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       water_vapor_tl:   profile layer tangent-linear water vapor mixing
!                         ratio array, i.e. the water vapor mixing
!                         ratio perturbation.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
! optional input arguments:
!       none
!
! output arguments:
!       predictor_tl:     array containing the calculated tangent-linear
!                         standard predictor set.
!                         units:      varies depending on predictor index.
!                         type:       real( fp_kind )
!                         dimension:  istd x k
!                         attributes: intent( out )
!
! optional output arguments:
!       none
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
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!       the standard predictors are the following:
!
!         1) temperature, t
!         2) pressure, p
!         3) t^2
!         4) p^2
!         5) t.p
!         6) t^2.p
!         7) t.p^2
!         8) t^2.p^2
!         9) water vapor mixing ratio, q
!
!       thus the tangent-linear form of these are
!
!         1) dt
!         2) dp
!         3) 2t.dt
!         4) 2p.dp
!         5) p.dt + t.dp
!         6) 2tp.dt + t^2.dp
!         7) 2tp.dp + p^2.dt
!         8) 2t(p^2).dt + 2(t^2)p.dp
!         9) dq
!
!p-
!--------------------------------------------------------------------------------

  subroutine compute_std_predictors_tl( p,           &  ! input,  k
                                        t,           &  ! input,  k
                                        w,           &  ! input,  k

                                        p_tl,        &  ! input,  k
                                        t_tl,        &  ! input,  k
                                        w_tl,        &  ! input,  k

                                        predictor_tl  )  ! output, istd x k



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ),     intent( in )  :: p              ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: t             ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: w             ! input,  k

    real( fp_kind ), dimension( : ),     intent( in )  :: p_tl          ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: t_tl          ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )  :: w_tl          ! input,  k

    real( fp_kind ), dimension( :, : ),  intent( out ) :: predictor_tl  ! output, istd x k


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_std_predictors_tl'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: k

    real( fp_kind ) :: p2, p2_tl
    real( fp_kind ) :: t2, t2_tl



    !#--------------------------------------------------------------------------#
    !#        -- calculate the tangent-linear standard predictor set --         #
    !#--------------------------------------------------------------------------#

    k_layer_loop: do k = 1, size( p )


      ! -- precalculate the squared terms
      p2 = p( k ) * p( k )
      t2 = t( k ) * t( k )

      ! -- tangent-linear of squared terms
      p2_tl = two * p( k ) * p_tl( k )
      t2_tl = two * t( k ) * t_tl( k )
      
      ! -- calculate and assign the absorber independent predictors
      predictor_tl( 1, k ) = t_tl( k )
      predictor_tl( 2, k ) = p_tl( k )
      predictor_tl( 3, k ) = t2_tl
      predictor_tl( 4, k ) = p2_tl
      predictor_tl( 5, k ) = ( t( k ) * p_tl( k ) ) + ( p( k ) * t_tl( k ) )
      predictor_tl( 6, k ) = ( p( k ) * t2_tl     ) + ( t2     * p_tl( k ) )
      predictor_tl( 7, k ) = ( t( k ) * p2_tl     ) + ( p2     * t_tl( k ) )
      predictor_tl( 8, k ) = ( t2     * p2_tl     ) + ( p2     * t2_tl     )
      predictor_tl( 9, k ) = w_tl( k )

    end do k_layer_loop

  end subroutine compute_std_predictors_tl





!--------------------------------------------------------------------------------
!p+
! name:
!       compute_int_predictors_tl
!
! purpose:
!       private function to compute the integrated, absorber dependent
!       predictor set for the tangent-linear transmittance model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_int_predictors_tl( pressure,       &  ! input, k
!                                       temperature,    &  ! input, k
!                                       absorber,       &  ! input, 0:k
!
!                                       pressure_tl,    &  ! input, k
!                                       temperature_tl, &  ! input, k
!                                       absorber_tl,    &  ! input, 0:k
!
!                                       predictor_tl    )  ! input, iint x k
!
! input arguments:
!       pressure:         profile layer pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       absorber_tl:      profile level integrated absorber amount array.
!                         units:      varies with absorber
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       pressure_tl:      profile layer tangent-linear pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature_tl:   profile layer tangent-linear temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       absorber_tl:      profile level tangent-linear integrated absorber
!                         amount array.
!                         units:      varies with absorber
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
! optional input arguments:
!       none
!
! output arguments:
!       predictor_tl:     array containing the calculated tangent-linear
!                         integrated predictor set for every absorber.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  iint x k
!                         attributes: intent( out )
!
! optional output arguments:
!       none
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
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!       the integrated predictors consist of six that are repeated for every
!       absorber:
!
!         1) t*
!         2) p*
!         3) t**
!         4) p**
!         5) t***
!         6) p***
!
!       the reference above details the predictor value for input level values. in
!       the equations below, the "x" quantities can be replaced with temperature, t, 
!       of pressure, p, to determine the required predictor, and a is the level
!       absorber amount:
!
!                           __ k
!                    1     \ 
!         x*(k) = --------  >  ( x(i) + x(i-1) ) ( a(i) - a(i-1) )
!                  c.a(k)  /__
!                             i=1
!
!                            __ k
!                     1     \
!        x**(k) = ---------  >  ( x(i)a(i) + x(i-1)a(i-1) ) ( a(i) - a(i-1) )
!                  c.a^2(k) /__
!                              i=1
!
!                            __ k
!                     1     \ 
!       x***(k) = ---------  >  ( x(i)a^2(i) + x(i-1)a^2(i-1) ) ( a(i) - a(i-1) )
!                  c.a^3(k) /__
!                              i=1
!
!       to accomodate input layer values (i.e. layer average values) from the
!       ncep gdas, the predictor formulations were modified to:
!
!                   __ k                        __k-1
!                  \                           \ 
!                   >  x(i)( a(i)-a(i-1) )      >  x(i)( a(i)-a(i-1) )
!                  /__                         /__ 
!                     i=1                         i=1
!         x*(k) = ------------------------- + -------------------------
!                            2.a(k)                    2.a(k-1)
!
!                   __ k
!                  \ 
!                   >  x(i)( a(i) + a(i-1) ) ( a(i) - a(i-1) )
!                  /__
!                     i=1
!        x**(k) = --------------------------------------------- +
!                                       2.a^2(k) 
!                   
!                   __k-1
!                  \ 
!                   >  x(i)( a(i) + a(i-1) ) ( a(i) - a(i-1) )
!                  /__
!                     i=1
!                 ---------------------------------------------
!                                       2.a^2(k-1)
!
!                      __ k
!                     \ 
!                 3 .  >  x(i)( a^2(i) + a^2(i-1) ) ( a(i) - a(i-1) )
!                     /__
!                        i=1
!       x***(k) = ---------------------------------------------------- +
!                                       4.a^3(k)                    
!
!                      __k-1
!                     \ 
!                 3 .  >  x(i)( a^2(i) + a^2(i-1) ) ( a(i) - a(i-1) )
!                     /__
!                        i=1
!                 ----------------------------------------------------
!                                       4.a^3(k-1)
!
!
!       the tangent-linear form of these layer predictor formulations are:
!
!                    __ k                                             
!                   \                                                 
!                    >  dx(i)( a(i)-a(i-1) ) + x(i)( da(i)-da(i-1) )  
!                   /__                                               
!                      i=1                                            
!         dx*(k) = -------------------------------------------------- -
!                                       2.a(k)                        
!
!                           __ k
!                          \
!                   da(k) . >  x(i)( a(i)-a(i-1) )
!                          /__
!                             i=1
!                  -------------------------------- +
!                               2.a^2(k)
!                 
!                    __k-1                                              
!                   \                                                   
!                    >  dx(i)( a(i)-a(i-1) ) + x(i)( da(i)-da(i-1) )    
!                   /__                                                 
!                      i=1                                              
!                  -------------------------------------------------- - 
!                                      2.a(k-1)                         
!
!                             __k-1
!                            \
!                   da(k-1) . >  x(i)( a(i)-a(i-1) )
!                            /__
!                               i=1
!                  ----------------------------------
!                              2.a^2(k-1)
!
!                 
!                 
!                    __ k                                                                                                       
!                   \                                                                                                           
!                    >  [dx(i)( a(i)+a(i-1) ) + x(i)( da(i)+da(i-1) )]( a(i)-a(i-1) ) + x(i)( a(i)+a(i-1) )( da(i)-da(i-1) )    
!                   /__                                                                                                         
!                      i=1                                                                                                      
!        dx**(k) = ---------------------------------------------------------------------------------------------------------- - 
!                                                                  2.a^2(k)                                                     
!
!                           __ k
!                          \
!                   da(k) . >  x(i)( a(i)+a(i-1) )( a(i)-a(i-1) )
!                          /__
!                             i=1
!                  ----------------------------------------------- +
!                                      a^3(k)
!
!                    __k-1                                                                                                      
!                   \                                                                                                           
!                    >  [dx(i)( a(i)+a(i-1) ) + x(i)( da(i)+da(i-1) )]( a(i)-a(i-1) ) + x(i)( a(i)+a(i-1) )( da(i)-da(i-1) )    
!                   /__                                                                                                         
!                      i=1                                                                                                      
!                  ---------------------------------------------------------------------------------------------------------- - 
!                                                                 2.a^2(k-1)                                                    
!
!                             __k-1
!                            \
!                   da(k-1) . >  x(i)( a(i)+a(i-1) )( a(i)-a(i-1) )
!                            /__
!                               i=1
!                  -------------------------------------------------
!                                      a^3(k-1)
!
!
!                 
!                      __ k                                                                                                       
!                     \            2    2                 2     2                                2    2                           
!                  3 . >  [dx(i)( a(i)+a(i-1) ) + x(i)( da(i)+da(i-1) )]( a(i)-a(i-1) ) + x(i)( a(i)+a(i-1) )( da(i)-da(i-1) )    
!                     /__                                                                                                         
!                        i=1                                                                                                      
!       dx***(k) = ------------------------------------------------------------------------------------------------------------ - 
!                                                                  4.a^3(k)                                                       
!
!                            __ k
!                           \          2    2
!                  9.da(k) . >  x(i)( a(i)+a(i-1) )( a(i)-a(i-1) )
!                           /__
!                              i=1
!                  ------------------------------------------------ +
!                                      4.a^4(k)
!
!                      __k-1                                                                                                      
!                     \            2    2                 2     2                                2    2                           
!                  3 . >  [dx(i)( a(i)+a(i-1) ) + x(i)( da(i)+da(i-1) )]( a(i)-a(i-1) ) + x(i)( a(i)+a(i-1) )( da(i)-da(i-1) )    
!                     /__                                                                                                         
!                        i=1                                                                                                      
!                  ------------------------------------------------------------------------------------------------------------ - 
!                                                                 4.a^3(k-1)                                                      
!
!                            __k-1
!                           \          2    2
!                  9.da(k) . >  x(i)( a(i)+a(i-1) )( a(i)-a(i-1) )
!                           /__
!                              i=1
!                  ------------------------------------------------
!                                      4.a^4(k-1)
!
!
!       thus the transmittance model coefficients calculated using the level
!       predictor formulation are used with the tangent-linear predictors
!       constructed above with layer values.
!
!p-
!--------------------------------------------------------------------------------

  subroutine compute_int_predictors_tl( pressure,       &  ! input,  k
                                        temperature,    &  ! input,  k
                                        absorber,       &  ! input,  0:k

                                        pressure_tl,    &  ! input,  k
                                        temperature_tl, &  ! input,  k
                                        absorber_tl,    &  ! input,  0:k

                                        predictor_tl    )  ! output, iint x k
 


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    real( fp_kind ), dimension( : ),    intent( in )  :: pressure        ! k
    real( fp_kind ), dimension( : ),    intent( in )  :: temperature     ! k
    real( fp_kind ), dimension( 0: ),   intent( in )  :: absorber        ! 0:k

    real( fp_kind ), dimension( : ),    intent( in )  :: pressure_tl     ! k
    real( fp_kind ), dimension( : ),    intent( in )  :: temperature_tl  ! k
    real( fp_kind ), dimension( 0: ),   intent( in )  :: absorber_tl     ! 0:k

    real( fp_kind ), dimension( :, : ), intent( out ) :: predictor_tl    ! iint x k


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_int_predictors_tl'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: i, k
    integer :: n_predictors

    real( fp_kind ) :: d_absorber
    real( fp_kind ) :: d_absorber_tl
    real( fp_kind ) :: factor_1
    real( fp_kind ) :: factor_1_tl
    real( fp_kind ) :: factor_2
    real( fp_kind ) :: factor_2_tl
    real( fp_kind ) :: inverse_1
    real( fp_kind ) :: inverse_2
    real( fp_kind ) :: inverse_3
    real( fp_kind ) :: inverse_4
    real( fp_kind ) :: absorber_3
    real( fp_kind ) :: absorber_4
    real( fp_kind ) :: inverse_1_tl
    real( fp_kind ) :: inverse_2_tl
    real( fp_kind ) :: inverse_3_tl


    ! -- square of the absorber amount. 0:k
    real( fp_kind ), dimension( 0:size( pressure ) ) :: absorber_2

    ! -- intermediate summation arrays. iint
    real( fp_kind ), dimension( size( predictor_tl, dim=1 ) ) :: s
    real( fp_kind ), dimension( size( predictor_tl, dim=1 ) ) :: s_tl

    ! -- level predictor, iint x 0:k
    real( fp_kind ), dimension( size( predictor_tl, dim=1 ), 0:size( pressure ) ) :: x_tl



    !#--------------------------------------------------------------------------#
    !#          -- determine the number of layers and predictors --             #
    !#--------------------------------------------------------------------------#

    n_predictors = size( predictor_tl, dim = 1 )



    !#--------------------------------------------------------------------------#
    !#                         -- initialise values --                          #
    !#--------------------------------------------------------------------------#

    absorber_2( 0 ) = zero

    s( : )       = zero
    s_tl( : )    = zero
    x_tl( :, 0 ) = zero



    !#--------------------------------------------------------------------------#
    !#               -- calculate the integrated predictor set --               #
    !#--------------------------------------------------------------------------#

    k_layer_loop: do k = 1, size( pressure )


      ! --------------------------------
      ! calculate multiplicative factors
      ! --------------------------------

      absorber_2( k ) = absorber( k ) * absorber( k )

      ! -- for the * terms
      d_absorber    = absorber( k )    - absorber( k-1 )
      d_absorber_tl = absorber_tl( k ) - absorber_tl( k-1 )

      ! -- for the ** terms
      factor_1    = ( absorber( k ) + absorber( k-1 ) ) * d_absorber
      factor_1_tl = ( ( absorber( k )    + absorber( k-1 )    ) * d_absorber_tl ) + &
                    ( ( absorber_tl( k ) + absorber_tl( k-1 ) ) * d_absorber    )

      ! -- for the *** terms       
      factor_2    = ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber
      factor_2_tl = ( ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber_tl ) + &
                    ( two * ( ( absorber( k )   * absorber_tl( k )   ) + &
                              ( absorber( k-1 ) * absorber_tl( k-1 ) ) ) * d_absorber )


      ! -------------------------------
      ! calculate the intermediate sums
      ! -------------------------------

      ! -- t*
      s( 1 )    = s( 1 )    + ( temperature( k )    * d_absorber    )     ! forward predictor
      s_tl( 1 ) = s_tl( 1 ) + ( temperature_tl( k ) * d_absorber    ) + &
                              ( temperature( k )    * d_absorber_tl )

      ! -- p*
      s( 2 )    = s( 2 )    + ( pressure( k )       * d_absorber    )     ! forward predictor
      s_tl( 2 ) = s_tl( 2 ) + ( pressure_tl( k )    * d_absorber    ) + &
                              ( pressure( k )       * d_absorber_tl )

      ! -- t**
      s( 3 )    = s( 3 )    + ( temperature( k )    * factor_1    )       ! forward predictor
      s_tl( 3 ) = s_tl( 3 ) + ( temperature_tl( k ) * factor_1    ) + &
                              ( temperature( k )    * factor_1_tl )

      ! -- p**
      s( 4 )    = s( 4 )    + ( pressure( k )       * factor_1    )       ! forward predictor
      s_tl( 4 ) = s_tl( 4 ) + ( pressure_tl( k )    * factor_1    ) + &
                              ( pressure( k )       * factor_1_tl )

      ! -- t***
      s( 5 )    = s( 5 )    + ( temperature( k )    * factor_2    )       ! forward predictor
      s_tl( 5 ) = s_tl( 5 ) + ( temperature_tl( k ) * factor_2    ) + &
                              ( temperature( k )    * factor_2_tl )

      ! -- p***
      s( 6 )    = s( 6 )    + ( pressure( k )       * factor_2    )       ! forward predictor
      s_tl( 6 ) = s_tl( 6 ) + ( pressure_tl( k )    * factor_2    ) + &
                              ( pressure( k )       * factor_2_tl )


      ! ------------------------------------------------------
      ! calculate the normalising factors for the integrated
      ! tangent-linear predictors. note that the checks below,
      ! the if tests to determine if the absorber products are
      ! represenatble are to minimise the number of calcs. i.e
      ! if inverse_1 is toast because absorber(k) is too small
      ! there's no need to check any further.
      ! ------------------------------------------------------

      ! -- is inverse_1 representable?
      inverse_1_check: if ( absorber( k ) > tolerance ) then

        inverse_1 = one / absorber( k )

        ! -- is inverse_2 representable
        inverse_2_check: if ( absorber_2( k ) > tolerance ) then

          inverse_2    =  inverse_1 * inverse_1
          inverse_1_tl = -inverse_2 * absorber_tl( k )
          absorber_3   =  absorber( k ) * absorber_2( k )
         
          ! -- is inverse_3 representable
          inverse_3_check: if ( absorber_3 > tolerance ) then

            inverse_3    =  inverse_2 * inverse_1
            inverse_2_tl = -inverse_3 * absorber_tl( k ) * two
            absorber_4   =  absorber( k ) * absorber_3

            ! -- is inverse_4 represenatble?
            inverse_4_check: if ( absorber_4 > tolerance ) then

              inverse_4    =  inverse_3 * inverse_1
              inverse_3_tl = -inverse_4 * absorber_tl( k ) * three

            else

              inverse_3_tl = zero

            end if inverse_4_check

          else

            inverse_3 = zero

            inverse_2_tl = zero
            inverse_3_tl = zero

          end if inverse_3_check

        else

          inverse_2 = zero
          inverse_3 = zero

          inverse_1_tl = zero
          inverse_2_tl = zero
          inverse_3_tl = zero

        end if inverse_2_check

      else

        inverse_1 = zero
        inverse_2 = zero
        inverse_3 = zero

        inverse_1_tl = zero
        inverse_2_tl = zero
        inverse_3_tl = zero

      end if inverse_1_check


      ! ------------------------------------------------------------
      ! scale and normalise the tangent-linear integrated predictors
      ! ------------------------------------------------------------

      ! -- t*
      x_tl( 1, k ) = point_5  * ( ( s_tl( 1 ) * inverse_1    ) + &
                                  ( s( 1 )    * inverse_1_tl ) )

      ! -- p*
      x_tl( 2, k ) = point_5  * ( ( s_tl( 2 ) * inverse_1    ) + &
                                  ( s( 2 )    * inverse_1_tl ) )

      ! -- t**
      x_tl( 3, k ) = point_5  * ( ( s_tl( 3 ) * inverse_2    ) + &
                                  ( s( 3 )    * inverse_2_tl ) )

      ! -- p**
      x_tl( 4, k ) = point_5  * ( ( s_tl( 4 ) * inverse_2    ) + &
                                  ( s( 4 )    * inverse_2_tl ) )

      ! -- t***
      x_tl( 5, k ) = point_75 * ( ( s_tl( 5 ) * inverse_3    ) + &
                                  ( s( 5 )    * inverse_3_tl ) )

      ! -- p***
      x_tl( 6, k ) = point_75 * ( ( s_tl( 6 ) * inverse_3    ) + &
                                  ( s( 6 )    * inverse_3_tl ) )


      ! ----------------------------
      ! sum predictors across layers
      ! ----------------------------

      do i = 1, n_predictors

        predictor_tl( i, k ) = x_tl( i, k ) + x_tl( i, k - 1 )

      end do

    end do k_layer_loop

  end subroutine compute_int_predictors_tl




!--------------------------------------------------------------------------------
!s+
! name:
!       compute_predictors_ad
!
! purpose:
!       public routine to calculate the adjoint transmittance model
!       predictors.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_predictors_ad ( &
!                                    ! -- forward input
!                                    pressure,       &  ! input,  k
!                                    temperature,    &  ! input,  k
!                                    water_vapor,    &  ! input,  k
!                                    absorber,       &  ! input,  0:k x j
!
!                                    ! -- adjoint input
!                                    predictor_ad,   &  ! in/output, i x k
!
!                                    ! -- adjoint output
!                                    pressure_ad,    &  ! in/output,  k
!                                    temperature_ad, &  ! in/output,  k
!                                    water_vapor_ad, &  ! in/output,  k
!                                    absorber_ad,    &  ! in/output,  0:k x j
!
!                                    ! -- optional input
!                                    no_standard     )  ! optional input
!
! input arguments:
!       pressure:         profile layer average pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer average temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       water_vapor:      profile layer average water vapor mixing ratio array.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       absorber:         profile level integrated absorber amount array.
!                         units:      varies with absorber
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       predictor_ad:     profile layer predictor adjoint array.
!                         ** this argument is set to zero on output **.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  i x k
!                         attributes: intent( in out ), target
!
!
! optional input arguments:
!       no_standard:      if present, the standard predictors are not calculated.
!                         this prevents recalculation of the standard predictors
!                         is only the view angle has changed - which only affects
!                         the integrated predictors.
!                         units:      none
!                         type:       integer
!                         dimension:  scalar
!                         attributes: intent( in ), optional
!
! output arguments:
!       pressure_ad:      profile layer adjoint pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
!       temperature_ad:   profile layer adjoint temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
!       water_vapor_ad:   profile layer adjoint water vapor array.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
!       absorber_ad:      profile level adjoint integrated absorber
!                         amount array.
!                         units:      varies with absorber
!                         type:       real( fp_kind )
!                         dimension:  0:k
!                         attributes: intent( in out )
!
! optional output arguments:
!       none.
!
! calls:
!       compute_std_predictors_tl:   private function to compute the tangent-
!                                    linear form of the standard (absorber
!                                    independent) predictor set.
!
!       compute_int_predictors_tl:   private function to compute the tangent-
!                                    linear form of the absorber integrated
!                                    predictor set.
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
!       the predictors used in the ncep transmittance model are organised in
!       the following manner:
!
!       --------------------------------------------
!       | 1 | 2 | 3 | ... | 9 | 10 | 11 | ... | 27 |
!       --------------------------------------------
!
!       \                    / \                   /
!        \                  /   \                 /
!         ------------------     -----------------
!                  |                      |
!                  v                      v
!
!              standard               integrated
!             predictors            predictors for
!                                   each absorber
!
!       pointers are used to reference the module scope predictor data array
!       before the calls are made to the standard and integrated predictor
!       calculation functions. this eliminates array copying.
!s-
!--------------------------------------------------------------------------------

  subroutine compute_predictors_ad ( &
                                     ! -- forward input
                                     pressure,       &  ! input,  k
                                     temperature,    &  ! input,  k
                                     water_vapor,    &  ! input,  k
                                     absorber,       &  ! input,  0:k x j

                                     ! -- adjoint input
                                     predictor_ad,   &  ! in/output, i x k

                                     ! -- adjoint output
                                     pressure_ad,    &  ! in/output,  k
                                     temperature_ad, &  ! in/output,  k
                                     water_vapor_ad, &  ! in/output,  k
                                     absorber_ad,    &  ! in/output,  0:k x j

                                     ! -- optional input
                                     no_standard     )  ! optional input



    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( : ),     intent( in )             :: pressure        ! k
    real( fp_kind ), dimension( : ),     intent( in )             :: temperature     ! k
    real( fp_kind ), dimension( : ),     intent( in )             :: water_vapor     ! k
    real( fp_kind ), dimension( 0:, : ), intent( in )             :: absorber        ! 0:k x j

    ! -- adjoint input
    real( fp_kind ), dimension( :, : ),  intent( in out ), target :: predictor_ad    ! i x k

    ! -- adjoint output
    real( fp_kind ), dimension( : ),     intent( in out )         :: pressure_ad     ! k
    real( fp_kind ), dimension( : ),     intent( in out )         :: temperature_ad  ! k
    real( fp_kind ), dimension( : ),     intent( in out )         :: water_vapor_ad  ! k
    real( fp_kind ), dimension( 0:, : ), intent( in out )         :: absorber_ad     ! 0:k x j

    ! -- optional input
    integer,                             intent( in ), optional   :: no_standard


    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_predictors_ad'


    ! ---------------
    ! local variables
    ! ---------------

    character( 80 ) :: message

    integer :: i1, i2, j

    real( fp_kind ), pointer, dimension( :, : ) :: std_predictors_ad   ! istd x k
    real( fp_kind ), pointer, dimension( :, : ) :: int_predictors_ad   ! iint x k



    !#--------------------------------------------------------------------------#
    !#         -- calculate the adjoint of the integrated predictors --         #
    !#--------------------------------------------------------------------------#

    j_absorber_loop: do j = 1, size( absorber, dim = 2 )

      ! -- determine indices of current absorber predictors
      i1 = max_n_standard_predictors + ( ( j - 1 ) * max_n_integrated_predictors ) + 1
      i2 = i1 + max_n_integrated_predictors - 1

      ! -- alias the input predictor array
      int_predictors_ad => predictor_ad( i1:i2, : )

      ! -- compute the predictor adjoints for the current absorber
      call compute_int_predictors_ad( &
                                      ! -- forward input
                                      pressure,            &  ! input,  k
                                      temperature,         &  ! input,  k
                                      absorber( 0:, j ),   &  ! input,  0:k

                                      ! -- adjoint input
                                      int_predictors_ad,   &  ! in/output, iint x k
                                                              
                                      ! -- adjoint output
                                      pressure_ad,         &  ! in/output,  k
                                      temperature_ad,      &  ! in/output,  k
                                      absorber_ad( 0:, j ) )  ! in/output,  0:k
                                                              
    end do j_absorber_loop                                    




    !#--------------------------------------------------------------------------#
    !#     -- calculate the adjoint of the standard predictors if needed --     #
    !#--------------------------------------------------------------------------#

    if ( .not. present( no_standard ) ) then

      ! -- alias the input predictor array
      std_predictors_ad => predictor_ad( 1:max_n_standard_predictors, : )

      ! -- compute the predictor adjoints
      call compute_std_predictors_ad( &
                                      ! -- forward input
                                      pressure,          &  ! input,  k
                                      temperature,       &  ! input,  k
                                      water_vapor,       &  ! input,  k

                                      ! -- adjoint input
                                      std_predictors_ad, &  ! in/output, istd x k

                                      ! -- adjoint output
                                      pressure_ad,       &  ! in/output,  k
                                      temperature_ad,    &  ! in/output,  k
                                      water_vapor_ad     )  ! in/output,  k
    else

      predictor_ad( 1:max_n_standard_predictors, : ) = zero

    end if

  end subroutine compute_predictors_ad





!--------------------------------------------------------------------------------
!p+
! name:
!       compute_std_predictors_ad
!
! purpose:
!       private function to compute the standard, absorber independent
!       predictor set for the adjoint transmittance model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_std_predictors_ad( &
!                                       ! -- forward input
!                                       pressure,       &  ! input,  k
!                                       temperature,    &  ! input,  k
!                                       water_vapor,    &  ! input,  k
!
!                                       ! -- adjoint input
!                                       predictor_ad,   &  ! in/output, istd x k
!
!                                       ! -- adjoint output
!                                       pressure_ad,    &  ! in/output,  k
!                                       temperature_ad, &  ! in/output,  k
!                                       water_vapor_ad  )  ! in/output,  k
!
! input arguments:
!       pressure:         profile layer average pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer average temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       water_vapor:      profile layer average water vapor mixing ratio array.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       predictor_ad:     adjoint of the layer predictor arrays.
!                         ** this argument is set to zero on output **.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  iint x k
!                         attributes: intent( in out )
!
! optional input arguments:
!       none
!
! output arguments:
!       pressure_ad:      profile layer adjoint pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
!       temperature_ad:   profile layer adjoint temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
!       water_vapor_ad:   profile layer adjoint water vapor array.
!                         units:      g/kg
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
! optional output arguments:
!       none
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
!       mcmillin, l.m., l.j. crone, m.d. goldberg, and t.j. kleespies,
!         "atmospheric transmittance of an absorbing gas. 4. optran: a
!          computationally fast and accurate transmittance model for absorbing
!          with fixed and with variable mixing ratios at variable viewing
!          angles.", applied optics, 1995, v34, pp6269-6274.
!
!       the standard predictors are the following:
!
!         pred(1) = t, temperature
!         pred(2) = p, pressure
!         pred(3) = t^2
!         pred(4) = p^2
!         pred(5) = t.p
!         pred(6) = t^2.p
!         pred(7) = t.p^2
!         pred(8) = t^2.p^2
!         pred(9) = w, water vapor mixing ratio
!
!       the tangent-linear form of these are
!
!         dpred(1) = dt
!         dpred(2) = dp
!         dpred(3) = 2t.dt
!         dpred(4) = 2p.dp
!         dpred(5) = p.dt + t.dp
!         dpred(6) = 2tp.dt + t^2.dp
!         dpred(7) = 2tp.dp + p^2.dt
!         dpred(8) = 2t(p^2).dt + 2(t^2)p.dp
!         dpred(9) = dw
!
!       thus, the adjoint of the predictors are
!
!         d#t = 2t( p^2.d#pred(8) + p.d#pred(6) + d#pred(3) ) + 
!                   p^2.d#pred(7) + p.d#pred(5) + d#pred(1)
!
!         d#p = 2p( t^2.d#pred(8) + t.d#pred(7) + d#pred(4) ) + 
!                   t^2.d#pred(6) + t.d#pred(5) + d#pred(2)
!
!         d#w = d#pred(9)
!
!       where the "#"ed terms represent the adjoints.
!
!p-
!--------------------------------------------------------------------------------

  subroutine compute_std_predictors_ad( &
                                        ! -- forward input
                                        p,            &  ! input,  k
                                        t,            &  ! input,  k
                                        w,            &  ! input,  k

                                        ! -- adjoint input
                                        predictor_ad, &  ! in/output, istd x k

                                        ! -- adjoint output
                                        p_ad,         &  ! in/output,  k
                                        t_ad,         &  ! in/output,  k
                                        w_ad          )  ! in/output,  k




    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( : ),     intent( in )     :: p             ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )     :: t             ! input,  k
    real( fp_kind ), dimension( : ),     intent( in )     :: w             ! input,  k

    ! -- adjoint input
    real( fp_kind ), dimension( :, : ),  intent( in out ) :: predictor_ad  ! in/output, istd x k

    ! -- adjoint output
    real( fp_kind ), dimension( : ),     intent( in out ) :: p_ad          ! in/output,  k
    real( fp_kind ), dimension( : ),     intent( in out ) :: t_ad          ! in/output,  k
    real( fp_kind ), dimension( : ),     intent( in out ) :: w_ad          ! in/output,  k



    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_std_predictors_ad'


    ! ---------------
    ! local variables
    ! ---------------

    integer :: i, n_predictors
    integer :: k

    real( fp_kind ) :: p2, p2_ad
    real( fp_kind ) :: t2, t2_ad



    !#--------------------------------------------------------------------------#
    !#              -- determine the number of predictors --                    #
    !#--------------------------------------------------------------------------#

    n_predictors = size( predictor_ad, dim = 1 )



    !#--------------------------------------------------------------------------#
    !#        -- calculate the adjoints of the standard predictor set --        #
    !#                                                                          #
    !# don't have to loop backwards here as this is a parallel loop.            #
    !#                                                                          #
    !# pressure and temperature squared adjoint terms are not zeroed out every  #
    !# loop iteration as they are local to each iteration and can be simply     #
    !# re-assigned.                                                             #
    !#--------------------------------------------------------------------------#

    k_layer_loop: do k = 1, size( p )


      ! -- precalculate the squared terms
      p2 = p( k ) * p( k )
      t2 = t( k ) * t( k )

      ! -- pressure squared adjoint
      p2_ad = ( t2   * predictor_ad( 8, k ) ) + &   ! predictor #8, t^2.p^2
              ( t(k) * predictor_ad( 7, k ) ) + &   ! predictor #7, t.p^2
                       predictor_ad( 4, k )         ! predictor #4, p^2

      ! -- temperature squared adjoint
      t2_ad = ( p2   * predictor_ad( 8, k ) ) + &   ! predictor #8, t^2.p^2
              ( p(k) * predictor_ad( 6, k ) ) + &   ! predictor #6, t^2.p
                       predictor_ad( 3, k )         ! predictor #3, t^2

      ! -- water vapor adjoint
      w_ad( k ) = w_ad( k ) + predictor_ad( 9, k )  ! predictor #9, w

      ! -- temperature adjoint
      t_ad( k ) = t_ad( k ) + &
                  ( p2   * predictor_ad( 7, k ) ) + &   ! predictor #7, t.p^2
                  ( p(k) * predictor_ad( 5, k ) ) + &   ! predictor #5, t.p
                           predictor_ad( 1, k )   + &   ! predictor #1, t
                  ( two * t(k) * t2_ad )                ! t^2 term

      ! -- pressure adjoint
      p_ad( k ) = p_ad( k ) + &
                  ( t2   * predictor_ad( 6, k ) ) + &   ! predictor #6, t^2.p
                  ( t(k) * predictor_ad( 5, k ) ) + &   ! predictor #5, t.p
                           predictor_ad( 2, k )   + &   ! predictor #2, p
                  ( two * p(k) * p2_ad )                ! p^2 term

      ! -- zero output adjoint
      predictor_ad( :, k ) = zero

    end do k_layer_loop

  end subroutine compute_std_predictors_ad





!--------------------------------------------------------------------------------
!p+
! name:
!       compute_int_predictors_ad
!
! purpose:
!       private function to compute the integrated, absorber dependent predictor
!       set for the adjoint transmittance model.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call compute_int_predictors_ad( &
!                                       ! -- forward input
!                                       pressure,       &  ! input,  k
!                                       temperature,    &  ! input,  k
!                                       absorber,       &  ! input,  0:k
!
!                                       ! -- adjoint input
!                                       predictor_ad,   &  ! in/output, iint x k
!
!                                       ! -- adjoint output
!                                       pressure_ad,    &  ! in/output,  k
!                                       temperature_ad, &  ! in/output,  k
!                                       absorber_ad     )  ! in/output,  0:k
!
! input arguments:
!       pressure:         profile layer pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       temperature:      profile layer temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in )
!
!       absorber :        profile level integrated absorber amount array.
!                         units:      varies with absorber.
!                         type:       real( fp_kind )
!                         dimension:  0:k
!                         attributes: intent( in )
!
!       predictor_ad:     adjoint of the layer predictor arrays.
!                         ** this argument is set to zero on output **.
!                         units:      varies with predictor type.
!                         type:       real( fp_kind )
!                         dimension:  iint x k
!                         attributes: intent( in out )
!
! optional input arguments:
!       none
!
! output arguments:
!       pressure_ad:      profile layer adjoint pressure array.
!                         units:      hpa
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
!       temperature_ad:   profile layer adjoint temperature array.
!                         units:      kelvin
!                         type:       real( fp_kind )
!                         dimension:  k
!                         attributes: intent( in out )
!
!       absorber_ad:      profile level adjoint integrated absorber
!                         amount array.
!                         units:      varies with absorber
!                         type:       real( fp_kind )
!                         dimension:  0:k
!                         attributes: intent( in out )
!
!
! optional output arguments:
!       none
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
!       the input argument predictor_ad is set to zero on output.
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
!       the integrated predictors consist of six that are repeated for every
!       absorber:
!
!         1) t*
!         2) p*
!         3) t**
!         4) p**
!         5) t***
!         6) p***
!
!p-
!--------------------------------------------------------------------------------

  subroutine compute_int_predictors_ad( &
                                        ! -- forward input
                                        pressure,       &  ! input,  k
                                        temperature,    &  ! input,  k
                                        absorber,       &  ! input,  0:k

                                        ! -- adjoint input
                                        predictor_ad,   &  ! in/output, iint x k

                                        ! -- adjoint output
                                        pressure_ad,    &  ! in/output,  k
                                        temperature_ad, &  ! in/output,  k
                                        absorber_ad     )  ! in/output,  0:k

 


    !#--------------------------------------------------------------------------#
    !#                         -- type declarations --                          #
    !#--------------------------------------------------------------------------#

    ! ---------
    ! arguments
    ! ---------

    ! -- forward input
    real( fp_kind ), dimension( : ),    intent( in )     :: pressure        ! k
    real( fp_kind ), dimension( : ),    intent( in )     :: temperature     ! k
    real( fp_kind ), dimension( 0: ),   intent( in )     :: absorber        ! 0:k

    ! -- adjoint input
    real( fp_kind ), dimension( :, : ), intent( in out ) :: predictor_ad    ! iint x k

    ! -- adjoint output
    real( fp_kind ), dimension( : ),    intent( in out ) :: pressure_ad     ! k
    real( fp_kind ), dimension( : ),    intent( in out ) :: temperature_ad  ! k
    real( fp_kind ), dimension( 0: ),   intent( in out ) :: absorber_ad     ! 0:k



    ! ----------------
    ! local parameters
    ! ----------------

    character( * ), parameter :: routine_name = 'compute_int_predictors_ad'


    ! ---------------
    ! local variables
    ! ---------------

    ! -- square of the absorber amount. 0:k
    real( fp_kind ), dimension( 0:size( pressure ) ) :: absorber_2

    ! -- multiplicative factors, k
    real( fp_kind ), dimension( size( pressure ) ) :: d_absorber
    real( fp_kind ), dimension( size( pressure ) ) :: factor_1
    real( fp_kind ), dimension( size( pressure ) ) :: factor_2

    ! -- intermediate summation arrays, iint x 0:k and iint
    real( fp_kind ), dimension( size( predictor_ad, dim=1 ), 0:size( pressure ) ) :: s
    real( fp_kind ), dimension( size( predictor_ad, dim=1 ) ) :: s_ad

    ! -- level predictor, iint x 0:k
    real( fp_kind ), dimension( size( predictor_ad, dim=1 ), 0:size( pressure ) ) :: x_ad

    ! -- scalars
    integer :: i, n_predictors
    integer :: k, n_layers

    real( fp_kind ) :: d_absorber_ad
    real( fp_kind ) :: factor_1_ad
    real( fp_kind ) :: factor_2_ad

    real( fp_kind ) :: inverse_1
    real( fp_kind ) :: inverse_2
    real( fp_kind ) :: inverse_3
    real( fp_kind ) :: inverse_4
    real( fp_kind ) :: absorber_3
    real( fp_kind ) :: absorber_4
    real( fp_kind ) :: inverse_1_ad
    real( fp_kind ) :: inverse_2_ad
    real( fp_kind ) :: inverse_3_ad



    !#--------------------------------------------------------------------------#
    !#                       -- assign the dimensions --                        #
    !#--------------------------------------------------------------------------#

    n_predictors = size( predictor_ad, dim=1 )
    n_layers     = size( pressure )

!    if ( n_predictors /= max_n_integrated_predictors ) then....



    !#--------------------------------------------------------------------------#
    !#          -- recalculate the intermediate forward model sums --           #
    !#--------------------------------------------------------------------------#

    ! -----------------
    ! initialise arrays
    ! -----------------

    absorber_2( 0 ) = zero
    s( :, 0: )      = zero


    ! ----------------
    ! loop over layers
    ! ----------------

    k_layer_loop_forward: do k = 1, n_layers


      ! -----------------------------------------
      ! calculate absorber multiplicative factors
      ! and save for adjoint calculation.
      ! -----------------------------------------

      absorber_2( k ) = absorber( k ) * absorber( k )

      d_absorber( k ) = absorber( k ) - absorber( k-1 )                           ! for * terms
      factor_1( k )   = ( absorber( k )   + absorber( k-1 )   ) * d_absorber( k ) ! for ** terms
      factor_2( k )   = ( absorber_2( k ) + absorber_2( k-1 ) ) * d_absorber( k ) ! for *** terms


      ! ----------------------------------------
      ! calculate and save the intermediate sums
      ! ----------------------------------------

      s( 1, k ) = s( 1, k-1 ) + ( temperature( k ) * d_absorber( k ) )  ! t*
      s( 2, k ) = s( 2, k-1 ) + ( pressure( k )    * d_absorber( k ) )  ! p*

      s( 3, k ) = s( 3, k-1 ) + ( temperature( k ) * factor_1( k ) )    ! t**
      s( 4, k ) = s( 4, k-1 ) + ( pressure( k )    * factor_1( k ) )    ! p**

      s( 5, k ) = s( 5, k-1 ) + ( temperature( k ) * factor_2( k ) )    ! t***
      s( 6, k ) = s( 6, k-1 ) + ( pressure( k )    * factor_2( k ) )    ! p***

    end do k_layer_loop_forward



    !#--------------------------------------------------------------------------#
    !#                -- initialise local adjoint variables --                  #
    !#--------------------------------------------------------------------------#

    x_ad( :, 0: ) = zero
    s_ad( : )     = zero

    absorber_ad( 0 ) = zero



    !#--------------------------------------------------------------------------#
    !#            -- calculate the integrated predictor adjoints --             #
    !#--------------------------------------------------------------------------#


    ! -----------------------------------
    ! here loop order does matter as this
    ! is a sequential loop
    ! -----------------------------------

    k_layer_loop_adjoint: do k = n_layers, 1, -1


      ! -------------------------------------------------------
      ! calculate the normalising factors for the integrated
      ! predictors. note that the checks below, the if tests to
      ! determine if the absorber products are represenatble
      ! are to minimise the number of calcs. i.e if inverse_1
      ! is toast because absorber(k) is too small there's no
      ! need to check any further.
      ! -------------------------------------------------------

      ! -- is inverse_1 representable?
      inverse_1_check: if ( absorber( k ) > tolerance ) then

        inverse_1 = one / absorber( k )

        ! -- is inverse_2 representable
        inverse_2_check: if ( absorber_2( k ) > tolerance ) then

          inverse_2  =  inverse_1 * inverse_1
          absorber_3 =  absorber( k ) * absorber_2( k )
         
          ! -- is inverse_3 representable
          inverse_3_check: if ( absorber_3 > tolerance ) then

            inverse_3  =  inverse_2 * inverse_1
            absorber_4 =  absorber( k ) * absorber_3

            ! -- is inverse_4 represenatble?
            inverse_4_check: if ( absorber_4 > tolerance ) then

              inverse_4 =  inverse_3 * inverse_1

            else

              inverse_4 = zero

            end if inverse_4_check

          else

            inverse_3 = zero
            inverse_4 = zero

          end if inverse_3_check

        else

          inverse_2 = zero
          inverse_3 = zero
          inverse_4 = zero

        end if inverse_2_check

      else

        inverse_1 = zero
        inverse_2 = zero
        inverse_3 = zero
        inverse_4 = zero

      end if inverse_1_check


      ! --------------------------------------------
      ! adjoint of predictor summation across layers
      ! --------------------------------------------

      do i = 1, n_predictors

        x_ad( i, k )   = x_ad( i, k )   + predictor_ad( i, k )
        x_ad( i, k-1 ) = x_ad( i, k-1 ) + predictor_ad( i, k )
        predictor_ad( i, k ) = zero

      end do


      ! ----------------------------------------------------------------
      ! adjoint of the scaled and normalised level integrated predictors
      !
      ! note that the adjoint variables inverse_x_ad are local to this
      ! loop iteration so they are simply assigned when they are first
      ! used.
      ! ----------------------------------------------------------------

      ! -- p*** and t***, predictor indices #6 and 5
      s_ad( 6 )    = s_ad( 6 ) + ( point_75 * inverse_3 * x_ad( 6, k ) )
      s_ad( 5 )    = s_ad( 5 ) + ( point_75 * inverse_3 * x_ad( 5, k ) )
      inverse_3_ad = point_75 * ( ( s( 6, k ) * x_ad( 6, k ) ) + &
                                  ( s( 5, k ) * x_ad( 5, k ) ) )

      ! -- p** and t**, predictor indices #4 and 3
      s_ad( 4 )    = s_ad( 4 ) + ( point_5 * inverse_2 * x_ad( 4, k ) )
      s_ad( 3 )    = s_ad( 3 ) + ( point_5 * inverse_2 * x_ad( 3, k ) )
      inverse_2_ad = point_5 * ( ( s( 4, k ) * x_ad( 4, k ) ) + &
                                 ( s( 3, k ) * x_ad( 3, k ) ) )

      ! -- p* and t*, predictor indices #2 and 1
      ! -- simply assign a value for inverse_1_ad
      s_ad( 2 )    = s_ad( 2 ) + ( point_5 * inverse_1 * x_ad( 2, k ) )
      s_ad( 1 )    = s_ad( 1 ) + ( point_5 * inverse_1 * x_ad( 1, k ) )
      inverse_1_ad = point_5 * ( ( s( 2, k ) * x_ad( 2, k ) ) + &
                                 ( s( 1, k ) * x_ad( 1, k ) ) )

      ! -- zero out the level predictor adjoint array for current k.
      x_ad( :, k ) = zero

      ! -- adjoint of inverse terms. note that the inverse_x_ad
      ! -- terms are *not* zeroed out as they are re-assigned values
      ! -- each loop iteration above.
      absorber_ad( k ) = absorber_ad( k ) - (         inverse_2 * inverse_1_ad ) - &
                                            ( two *   inverse_3 * inverse_2_ad ) - &
                                            ( three * inverse_4 * inverse_3_ad )


      ! ----------------------------------------------------------
      ! pressure and temperature adjoints of the intermediate sums
      ! ----------------------------------------------------------

      ! -- pressure
      pressure_ad( k ) = pressure_ad( k ) + ( factor_2( k )   * s_ad( 6 ) ) + &  ! p***
                                            ( factor_1( k )   * s_ad( 4 ) ) + &  ! p**
                                            ( d_absorber( k ) * s_ad( 2 ) )      ! p*


      ! -- temperature
      temperature_ad( k ) = temperature_ad( k ) + ( factor_2( k )   * s_ad( 5 ) ) + &  ! t***
                                                  ( factor_1( k )   * s_ad( 3 ) ) + &  ! t**
                                                  ( d_absorber( k ) * s_ad( 1 ) )      ! t*


      ! -------------------------------------
      ! adjoint of the multiplicative factors
      ! -------------------------------------

      ! --------------------------------------------------
      ! note that the adjoint variables factor_x_ad and
      ! d_absorber_ad are local to this loop iteration
      ! so they are simply assigned when they are first
      ! used.
      !
      ! note there are no
      !   s_ad() = 0
      ! because all the tangent-linear forms are
      !   s_tl() = s_tl() + (...)
      ! summing from the previous layer.
      !
      ! note that the factor_x_ad and d_absorber_ad
      ! terms are *not* zeroed out as they are re-assigned
      ! values each loop iteration.
      ! --------------------------------------------------

      ! -- multiplicative factors
      factor_2_ad = ( pressure( k )    * s_ad( 6 ) ) + &
                    ( temperature( k ) * s_ad( 5 ) )

      factor_1_ad = ( pressure( k )    * s_ad( 4 ) ) + &
                    ( temperature( k ) * s_ad( 3 ) )

      d_absorber_ad = ( pressure( k )    * s_ad( 2 ) ) + &
                      ( temperature( k ) * s_ad( 1 ) )

      ! -- adjoint of factor_2
      absorber_ad( k-1 ) = absorber_ad( k-1 ) + ( two * absorber( k-1 ) * d_absorber( k ) * factor_2_ad )
      absorber_ad(  k  ) = absorber_ad(  k  ) + ( two * absorber(  k  ) * d_absorber( k ) * factor_2_ad )
      d_absorber_ad      = d_absorber_ad      + ( ( absorber_2( k ) + absorber_2( k-1 ) ) * factor_2_ad )

      ! -- adjoint of factor_1
      absorber_ad( k-1 ) = absorber_ad( k-1 ) + ( d_absorber( k ) * factor_1_ad )
      absorber_ad(  k  ) = absorber_ad(  k  ) + ( d_absorber( k ) * factor_1_ad )
      d_absorber_ad      = d_absorber_ad      + ( ( absorber( k ) + absorber( k-1 ) ) * factor_1_ad )

      ! -- adjoint of d_absorber
      absorber_ad( k-1 ) = absorber_ad( k-1 ) - d_absorber_ad
      absorber_ad(  k  ) = absorber_ad(  k  ) + d_absorber_ad

    end do k_layer_loop_adjoint

  end subroutine compute_int_predictors_ad

end module predictors


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
! revision 2.10  2001/08/16 16:45:05  paulv
! - updated documentation
!
! revision 2.9  2001/07/12 18:22:33  paulv
! - added more robust code for calculating the powers of the inverse absorber
!   amount. the squares, cubes and fourth powers of absorber amounts were
!   causing floating point underflows which, when used in a denominator were
!   greatly inflating any precision errors. so, the solution i adopted was,
!   prior to each inverse calculation, to check the value of the absorber
!   quantity (e.g. absorber**2, absorber**3, or absorber**4). if they are less
!   than a tolerance value (defined using the epsilon intrinsic) then their
!   inverse is set to zero - as are any higher power inverses. this prevents
!   too small values from being used.
!   this is mostly a problem for water vapor near at the top of the defined
!   atmosphere, although some low ozone values can be found (rarely).
! - changed all loop initialisation to vector expressions, e.g. from
!     do i = 1, n_predictors
!       s( i )       = zero
!       s_tl( i )    = zero
!       x_tl( i, 0 ) = zero
!     end do
!   to
!     s( : )       = zero
!     s_tl( : )    = zero
!     x_tl( :, 0 ) = zero
! - corrected bug in definition of the intermediate summation array in
!   compute_int_predictors_ad from
!     real( fp_kind ), dimension( size( predictor_ad, dim=1 ), size( pressure ) ) :: s
!   to
!     real( fp_kind ), dimension( size( predictor_ad, dim=1 ), 0:size( pressure ) ) :: s
! - added initialisation of absorber_2( 0 )
!     absorber_2( 0 ) = zero
!   in compute_int_predictors_ad.
! - corrected bug in intermediate summation loop. the sums were being calculated
!   like:
!     do k = 1, n_layers
!       s( i, k ) = s( i, k ) + .......
!   rather than:
!     do k = 1, n_layers
!       s( i, k ) = s( i, k-1 ) + .......
!   hence the need for the redefintion of s(k) to s(0:k).
!
! revision 2.8  2001/05/29 17:51:34  paulv
! - corrected some documentation errors.
!
! revision 2.7  2001/05/04 14:39:43  paulv
! - removed shared predictor arrays from module.
! - now use type_kinds module parameter fp_kind to set the floating point
!   data type.
! - added adjoint form of routines to module. use of no_standard optional
!   keyword in compute_predictors_ad has not yet been looked into.
! - changed names of private integrated predictor routines from
!   compute_stard_predictors to compute_int_predictors. the difference between
!   "stard" and "std" is small enough to cause mild confusion when scanning
!   the code.
! - changed references to parameter maximums to input array sizes, e.g.
!     j_absorber_loop: do j = 1, max_n_absorbers
!   becomes
!     j_absorber_loop: do j = 1, size( absorber, dim = 2 )
! - shortened variable names in standard predictor routines, e.g. p instead
!   of pressure. these calcs are simple enough that short names are clear.
! - absorber arrays are now dimensioned as 0:k. this eliminates the
!   need for using an absorber_km1 variable in computing the absorber layer
!   difference, d_absorber, and average, ave_absorber where the layer loop
!   always goes from 1 -> n_layers.
! - simplified compute_int_predictors_tl. this is a combination of the
!   change in the absorber array dimensioning and just thinking about it
!   for a while.
! - updated header documentation. adjoint routines not yet fully documented.
!
! revision 2.6  2001/04/03 20:02:56  paulv
! - commented out shared predictor data arrays. predictor arrays are now
!   passed arguments.
! - removed reference to profile number. calls to routines are now a simgle
!   profile passed per call.
! - removed planned allocation of predictor arrays.
! - correted bug in 1st and 2nd order predictor calculation.
!
! revision 2.5  2001/01/24 20:14:21  paulv
! - latest test versions.
!
! revision 2.4  2000/11/09 20:36:07  paulv
! - added tangent linear form of routines in this module.
! - changed names of private routines used to compute predictors. the
!   appearance of "standard" and "integrated" in routine names have been
!   replaced with "std" and "stard" respectively.
! - initial setup in place for dynamic allocation of predictor arrays. these
!   arrays are still statically allocated with parameterised dimensions.
!
! revision 2.3  2000/08/31 19:36:33  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 2.2  2000/08/24 16:43:29  paulv
! - added optional no_standard argument to compute_predictors subprogram.
!   using this argument prevents the standard predictors (which are angle
!   independent) from being recalculated when only the path angle has changed
!   in the calling procedure.
! - updated module and subprogram documentation.
!
! revision 2.1  2000/08/21 21:03:16  paulv
! - standard and integrated predictor sets calculated in separate
!   functions.
! - predictor values saved in linear store. simplifies application of
!   predictors when calculating absorption coefficients.
! - wrapper function "compute_predictors" added to simplify predictor
!   calculation.
!
! revision 1.1  2000/08/08 16:34:17  paulv
! initial checkin
!
!
!
