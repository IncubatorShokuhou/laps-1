!------------------------------------------------------------------------------
!m+
! name:
!       sensor_planck_routines
!
! purpose:
!       module containing the sensor planck function routines.
!
! category:
!       ncep rtm
!
! calling sequence:
!       use sensor_planck_routines
!
! outputs:
!       none.
!
! modules:
!       type_kinds:                    module containing data type kind definitions.
!
!       parameters:                    module containing parameter definitions for
!                                      the rt model.
!
!       spectral_coefficients:         module containing the rt model spectral
!                                      coefficients
!
! contains:
!       sensor_planck_radiance:        public subroutine to calculate the instrument
!                                      channel radiance.
!
!       sensor_planck_radiance_tl:     public subroutine to calculate the tangent-linear
!                                      instrument channel radiance.
!
!       sensor_planck_radiance_ad:     public subroutine to calculate the adjoint
!                                      instrument channel radiance.
!
!       sensor_planck_temperature:     public subroutine to calculate the instrument
!                                      channel brightness temperature.
!
!       sensor_planck_temperature_tl:  public subroutine to calculate the tangent-linear
!                                      instrument channel brightness temperature.
!
!       sensor_planck_temperature_ad:  public subroutine to calculate the adjoint
!                                      instrument channel brightness temperature.
!
! externals:
!       none
!
! common blocks:
!       none.
!
! restrictions:
!       these functions are called frequently so no input checking is
!       performed.
!
! creation history:
!       written by:     paul van delst, cimss/ssec 08-aug-2001
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
!m-
!------------------------------------------------------------------------------

module sensor_planck_routines


  ! ---------------------
  ! module use statements
  ! ---------------------

  use type_kinds, only : fp_kind
  use parameters
  use spectral_coefficients


  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none


  ! ------------
  ! visibilities
  ! ------------

  private

  public  :: sensor_planck_radiance
  public  :: sensor_planck_radiance_tl
  public  :: sensor_planck_radiance_ad

  public  :: sensor_planck_temperature
  public  :: sensor_planck_temperature_tl
  public  :: sensor_planck_temperature_ad


contains


!--------------------------------------------------------------------------------
!p+
! name:
!       sensor_planck_radiance
!
! purpose:
!       subroutine to calculate the instrument channel radiance.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call sensor_planck_radiance( channel,     &  ! input
!                                    temperature, &  ! input
!                                    radiance     )  ! output
!
! input arguments:
!       channel:     channel index id. this is a unique index
!                    to a (supported) sensor channel.
!                    units:      none
!                    type:       integer
!                    dimension:  scalar
!                    attributes: intent( in )
!
!       temperature: temperature for which the planck radiance is
!                    to be calculated.
!                    units:      kelvin, k
!                    type:       real
!                    dimension:  scalar
!                    attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       radiance:    channel planck radiance.
!                    units:      mw/(m^2.sr.cm^-1)
!                    type:       real
!                    dimension:  scalar
!                    attributes: intent( out )
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
!       none.
!
! restrictions:
!       spectral coefficients are obtained from the spectral_coefficients module
!       so only radiances for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! procedure:
!       first a polychromatic correction is applied to give an effective
!       temperature,
!
!         t_eff = bc1 + ( bc2 * t )
!
!       the sensor radiance is then calculated using the effective temperature:
!
!                       pc1
!         r = ------------------------
!              exp( pc2 / t_eff ) - 1
!
!       the bc1, bc2, pc1, and pc2 values are obtained from the 
!       spectral_coefficients module which is filled during the initialisation
!       phase.
!p-
!--------------------------------------------------------------------------------

  subroutine sensor_planck_radiance( channel,     &  ! input
                                     temperature, &  ! input
                                     radiance     )  ! output

    ! -- arguments
    integer,         intent( in )  :: channel
    real( fp_kind ), intent( in )  :: temperature
    real( fp_kind ), intent( out ) :: radiance

    ! -- local
    real( fp_kind ) :: effective_temperature

    intrinsic exp


    ! -------------------------------------
    ! apply the polychromaticity correction
    ! to obtain an effective temperature
    ! -------------------------------------

    effective_temperature = band_c1( channel ) + ( band_c2( channel ) * temperature )


    ! -----------------------------
    ! calculate the planck radiance
    ! -----------------------------

    radiance =                  planck_c1( channel )  / &
    !          -------------------------------------------------------------
               ( exp( planck_c2( channel ) / effective_temperature ) - one )

  end subroutine sensor_planck_radiance




!--------------------------------------------------------------------------------
!p+
! name:
!       sensor_planck_radiance_tl
!
! purpose:
!       subroutine to calculate the tangent-linear instrument channel radiance.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call sensor_planck_radiance_tl( channel,        &  ! input
!                                       temperature,    &  ! input
!                                       temperature_tl, &  ! input
!                                       radiance_tl     )  ! output
!
! input arguments:
!       channel:        channel index id. this is a unique index
!                       to a (supported) sensor channel.
!                       units:      none
!                       type:       integer
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       temperature:    temperature for which the tangent-linear planck radiance
!                       is to be calculated.
!                       units:      kelvin, k
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       temperature_tl: tangent-linear temperature for which the tangent-linear
!                       planck radiance is required.
!                       units:      kelvin, k
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       radiance_tl:    tangent-linear planck radiance.
!                       units:      mw/(m^2.sr.cm^-1)
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( out )
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
!       none.
!
! restrictions:
!       spectral coefficients are obtained from the spectral_coefficients module
!       so only radiances for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! procedure:
!       first a polychromatic correction is applied to give an effective
!       temperature,
!
!         t_eff = bc1 + ( bc2 . t )
!
!       the sensor tangent-linear radiance is then calculated by first computing
!       the exponent term,
!
!          exponent = exp( pc2 / t_eff )
!
!       and then the actual operator,
!
!                 pc1 . pc2 . bc1 . exponent
!         f = ---------------------------------
!              ( t_eff . ( exponent - 1 ) )^2
!
!       which is the derivate of the planck equation wrt temperature. the
!       tangent-linear radiance is then determined by,
!
!         dr = f . dt
!
!       where dt is the input tangent-linear temeprature.
!
!       the bc1, bc2, pc1, and pc2 values are obtained from the 
!       spectral_coefficients module which is filled during the initialisation
!       phase.
!p-
!--------------------------------------------------------------------------------

  subroutine sensor_planck_radiance_tl( channel,        &  ! input
                                        temperature,    &  ! input
                                        temperature_tl, &  ! input
                                        radiance_tl     )  ! output

    ! -- arguments
    integer,         intent( in )  :: channel
    real( fp_kind ), intent( in )  :: temperature
    real( fp_kind ), intent( in )  :: temperature_tl
    real( fp_kind ), intent( out ) :: radiance_tl

    ! -- local
    real( fp_kind ) :: effective_temperature
    real( fp_kind ) :: exponent
    real( fp_kind ) :: f

    intrinsic exp


    ! -------------------------------------
    ! apply the polychromaticity correction
    ! -------------------------------------

    effective_temperature = band_c1( channel ) + ( band_c2( channel ) * temperature )


    ! --------------------------------------
    ! calculate the planck function operator
    ! --------------------------------------

    ! -- the exponent term
    exponent = exp( planck_c2( channel ) / effective_temperature )

    ! -- the operator, call it f
    f =  planck_c1( channel ) * planck_c2( channel ) * exponent * band_c2( channel ) / &
    !   -----------------------------------------------------------------------------
                      ( effective_temperature * ( exponent - one ) )**2


    ! -------------------------------------
    ! calculate the tangent-linear radiance
    ! -------------------------------------

    radiance_tl = f * temperature_tl

  end subroutine sensor_planck_radiance_tl




!--------------------------------------------------------------------------------
!p+
! name:
!       sensor_planck_radiance_ad
!
! purpose:
!       subroutine to calculate the adjoint instrument channel radiance.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call sensor_planck_radiance_ad( channel,       &  ! input
!                                       temperature,   &  ! input
!                                       radiance_ad,   &  ! input
!                                       temperature_ad )  ! in/output
!
! input arguments:
!       channel:        channel index id. this is a unique index
!                       to a (supported) sensor channel.
!                       units:      none
!                       type:       integer
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       temperature:    temperature for which the tangent-linear planck radiance
!                       is to be calculated.
!                       units:      kelvin
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       radiance_ad:    adjoint planck radiance.
!                       units:      mw/(m2.sr.cm-1)
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       temperature_ad: adjoint planck temperature
!                       units:      kelvin
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in out )
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
!       the input adjoint radiance argument, radiance_ad, is not set to zero
!       before returning to the calling routine.
!
! restrictions:
!       spectral coefficients are obtained from the spectral_coefficients module
!       so only radiances for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! procedure:
!       first a polychromatic correction is applied to give an effective
!       temperature,
!
!         t_eff = bc1 + ( bc2 . t )
!
!       the sensor tangent-linear radiance is then calculated by first computing
!       the exponent term,
!
!          exponent = exp( pc2 / t_eff )
!
!       and then the actual operator,
!
!                 pc1 . pc2 . bc1 . exponent
!         f = ---------------------------------
!              ( t_eff . ( exponent - 1 ) )^2
!
!       which is the derivate of the planck equation wrt temperature. the
!       adjoint temperature is then determined from,
!
!         t_ad = t_ad + ( f . r_ad )
!
!       where t_ad and r_ad on the lhs are the input adjoint temperature and
!       radiance respectively.
!
!       the bc1, bc2, pc1, and pc2 values are obtained from the 
!       spectral_coefficients module which is filled during the initialisation
!       phase.
!p-
!--------------------------------------------------------------------------------

  subroutine sensor_planck_radiance_ad( channel,       &  ! input
                                        temperature,   &  ! input
                                        radiance_ad,   &  ! input
                                        temperature_ad )  ! in/output

    ! -- arguments
    integer,         intent( in )     :: channel
    real( fp_kind ), intent( in )     :: temperature
    real( fp_kind ), intent( in )     :: radiance_ad
    real( fp_kind ), intent( in out ) :: temperature_ad

    ! -- local
    real( fp_kind ) :: effective_temperature
    real( fp_kind ) :: exponent
    real( fp_kind ) :: f

    intrinsic exp


    ! -------------------------------------
    ! apply the polychromaticity correction
    ! -------------------------------------

    effective_temperature = band_c1( channel ) + ( band_c2( channel ) * temperature )


    ! --------------------------------------
    ! calculate the planck function operator
    ! --------------------------------------

    ! -- the exponent term
    exponent = exp( planck_c2( channel ) / effective_temperature )

    ! -- the operator, call it f
    f =  planck_c1( channel ) * planck_c2( channel ) * exponent * band_c2( channel ) / &
    !   -----------------------------------------------------------------------------
                      ( effective_temperature * ( exponent - one ) )**2


    ! ---------------------------------
    ! calculate the adjoint temperature
    ! ---------------------------------

    temperature_ad = temperature_ad + ( f * radiance_ad )

  end subroutine sensor_planck_radiance_ad




!--------------------------------------------------------------------------------
!p+
! name:
!       sensor_planck_temperature
!
! purpose:
!       subroutine to calculate the instrument channel brightness temperature.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call sensor_planck_temperature( channel,    &  ! input
!                                       radiance,   &  ! input
!                                       temperature )  ! output
!
! input arguments:
!       channel:     channel index id. this is a unique index
!                    to a (supported) sensor channel.
!                    units:      none
!                    type:       integer
!                    dimension:  scalar
!                    attributes: intent( in )
!
!       radiance:    radiance for which the planck temperature is desired.
!                    units:      mw/(m^2.sr.cm^-1)
!                    type:       real
!                    dimension:  scalar
!                    attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       temperature: planck temperature.
!                    units:      kelvin, k
!                    type:       real
!                    dimension:  scalar
!                    attributes: intent( in )
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
!       none.
!
! restrictions:
!       spectral coefficients are obtained from the spectral_coefficients module
!       so only temperatures for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! procedure:
!       first the effective temperature is calculated from the inverse planck function,
!
!                        pc2
!         t_eff = ------------------
!                  log( pc1/r + 1 )
!
!       the polychromatic correction is then removed to provide the brightness
!       temperature,
!
!              t_eff - bc1
!         t = -------------
!                  bc2
!
!       the bc1, bc2, pc1, and pc2 values are obtained from the 
!       spectral_coefficients module which is filled during the initialisation
!       phase.
!p-
!--------------------------------------------------------------------------------

  subroutine sensor_planck_temperature( channel,    &  ! input
                                        radiance,   &  ! input
                                        temperature )  ! output

    ! -- arguments
    integer,         intent( in )  :: channel
    real( fp_kind ), intent( in )  :: radiance
    real( fp_kind ), intent( out ) :: temperature

    ! -- local
    real( fp_kind ) :: effective_temperature

    intrinsic log


    ! -----------------------------------
    ! calculate the effective temperature
    ! -----------------------------------

    effective_temperature =              planck_c2( channel )  / &
    !                       ------------------------------------------------
                            log( ( planck_c1( channel ) / radiance ) + one )

    ! -------------------------------------
    ! apply the polychromatic correction to 
    ! obtain the true temperature
    ! -------------------------------------

    temperature = ( effective_temperature - band_c1( channel ) ) / &
    !             ----------------------------------------------
                                band_c2( channel )

  end subroutine sensor_planck_temperature




!--------------------------------------------------------------------------------
!p+
! name:
!       sensor_planck_temperature_tl
!
! purpose:
!       subroutine to calculate the tangent-linear instrument channel
!       brightness temperature.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call sensor_planck_temperature_tl( channel,       &  ! input
!                                          radiance,      &  ! input
!                                          radiance_tl,   &  ! input
!                                          temperature_tl )  ! output
!
! input arguments:
!       channel:        channel index id. this is a unique index
!                       to a (supported) sensor channel.
!                       units:      none
!                       type:       integer
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       radiance:       radiance at which the tangent-linear planck temperature
!                       is desired.
!                       units:      mw/(m^2.sr.cm^-1)
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       radiance_tl:    tangent-linear radiance for which the tangent-linear
!                       planck temperature is desired.
!                       units:      mw/(m^2.sr.cm^-1)
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       temperature_tl: tangent-linear planck temperature.
!                       units:      kelvin, k
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
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
!       none.
!
! restrictions:
!       spectral coefficients are obtained from the spectral_coefficients module
!       so only temperatures for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! procedure:
!       first the logarithm argument is calculated,
!
!         a = pc1/r + 1
!
!       the inverse planck function operator is then calculated,
!
!                      pc1 . pc2
!         f = ------------------------------
!              bc2 . a . ( r . log( a ) )^2
!
!       and the tangent-linear temperature is then given by,
!
!         dt = f . dr
!
!       the bc1, bc2, pc1, and pc2 values are obtained from the 
!       spectral_coefficients module which is filled during the initialisation
!       phase.
!p-
!--------------------------------------------------------------------------------


  subroutine sensor_planck_temperature_tl( channel,       &  ! input
                                           radiance,      &  ! input
                                           radiance_tl,   &  ! input
                                           temperature_tl )  ! output

    ! -- arguments
    integer,         intent( in )  :: channel
    real( fp_kind ), intent( in )  :: radiance
    real( fp_kind ), intent( in )  :: radiance_tl
    real( fp_kind ), intent( out ) :: temperature_tl

    ! -- local
    real( fp_kind ) :: argument
    real( fp_kind ) :: f

    intrinsic log


    ! --------------------------------------
    ! calculate the planck function operator
    ! --------------------------------------

    ! -- the logarithm argument
    argument = ( planck_c1( channel ) / radiance ) + one

    ! -- the operator, call it f
    f =             planck_c1( channel ) * planck_c2( channel ) / &
    !   -----------------------------------------------------------------------
         ( band_c2( channel ) * argument * ( radiance * log( argument ) )**2 )


    ! ----------------------------------------
    ! calculate the tangent-linear temperature
    ! ----------------------------------------

    temperature_tl = f * radiance_tl

  end subroutine sensor_planck_temperature_tl



!--------------------------------------------------------------------------------
!p+
! name:
!       sensor_planck_temperature_ad
!
! purpose:
!       subroutine to calculate the adjoint instrument channel
!       brightness temperature.
!
! category:
!       ncep rtm
!
! calling sequence:
!       call sensor_planck_temperature_ad( channel,        &  ! input
!                                          radiance,       &  ! input
!                                          temperature_ad, &  ! input
!                                          radiance_ad     )  ! in/output
!
! input arguments:
!       channel:        channel index id. this is a unique index
!                       to a (supported) sensor channel.
!                       units:      none
!                       type:       integer
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       radiance:       radiance at which the adjoint radiance is desired.
!                       units:      mw/(m^2.sr.cm^-1)
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
!       temperature_ad: adjoint planck temperature.
!                       units:      kelvin, k
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in )
!
! optional input arguments:
!       none.
!
! output arguments:
!       radiance_ad:    adjoint radiance.
!                       units:      k.m^2.sr.cm^-1/mw
!                       type:       real
!                       dimension:  scalar
!                       attributes: intent( in out )
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
!       the input adjoint temperature argument, temperature_ad, is not set to zero
!       before returning to the calling routine.
!
! restrictions:
!       spectral coefficients are obtained from the spectral_coefficients module
!       so only temperatures for those sensors which are included in the spectral
!       coefficient data file can be calculated.
!
! procedure:
!       first the logarithm argument is calculated,
!
!         a = pc1/r + 1
!
!       the inverse planck function operator is then calculated,
!
!                      pc1 . pc2
!         f = ------------------------------
!              bc2 . a . ( r . log( a ) )^2
!
!       which is the derivate of the planck temperature wrt radiance. the
!       adjoint radiance is then determined from,
!
!         r_ad = r_ad + ( f . t_ad )
!
!       where r_ad and t_ad on the lhs are the input adjoint radiance and
!       temperature respectively.
!
!       the bc1, bc2, pc1, and pc2 values are obtained from the 
!       spectral_coefficients module which is filled during the initialisation
!       phase.
!p-
!--------------------------------------------------------------------------------


  subroutine sensor_planck_temperature_ad( channel,        &  ! input
                                           radiance,       &  ! input
                                           temperature_ad, &  ! input
                                           radiance_ad     )  ! in/output

    ! -- arguments
    integer,         intent( in )     :: channel
    real( fp_kind ), intent( in )     :: radiance
    real( fp_kind ), intent( in )     :: temperature_ad
    real( fp_kind ), intent( in out ) :: radiance_ad

    ! -- local
    real( fp_kind ) :: argument
    real( fp_kind ) :: f

    intrinsic log


    ! --------------------------------------
    ! calculate the planck function operator
    ! --------------------------------------

    ! -- the logarithm argument
    argument = ( planck_c1( channel ) / radiance ) + one

    ! -- the operator, call it f
    f =             planck_c1( channel ) * planck_c2( channel ) / &
    !   -----------------------------------------------------------------------
         ( band_c2( channel ) * argument * ( radiance * log( argument ) )**2 )


    ! ------------------------------
    ! calculate the adjoint radiance
    ! ------------------------------

    radiance_ad = radiance_ad + ( f * temperature_ad )

  end subroutine sensor_planck_temperature_ad

end module sensor_planck_routines


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
! revision 1.1  2001/08/08 20:04:03  paulv
! initial checkin.
! - routines were extracted from the radiance module and placed in their
!   own module to facilitate code-sharing.
!
!
!
