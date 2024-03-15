!------------------------------------------------------------------------------
!m+
! name:
!       parameters
!
! purpose:
!       module to hold rt model parameter constants
!
! category:
!       ncep rtm
!
! calling sequence:
!       use parameters
!
! outputs:
!       parameters
!       ----------
!       max_n_absorbers:             integer parameter defining the maximum
!                                    number of absorbing species.
!
!       max_n_predictors_used:       integer parameter defining the maximum
!                                    number of predictors used in the absorption
!                                    coefficient calculation.
!
!       max_n_standard_predictors:   integer parameter defining the number of
!                                    standard (i.e. absorber independent) predictors.
!
!       max_n_integrated_predictors: integer parameter defining the number of
!                                    integrated (i.e. absorber dependent) predictors.
!
!       max_n_predictors:            integer parameter defining the total number
!                                    of predictors for all absorbers,
!                                      max_n_standard_predictors + &
!                                      ( max_n_absorbers * max_n_integrated_predictors )
!
!       max_n_layers:                integer parameter defining the maximum
!                                    number of atmospheric layers allowed for
!                                    input.
!
!       max_n_absorber_layers:       integer parameter defining the maximum
!                                    number of
!                                    absorber space layers.
!
!       the following are defined numerical parameters. parameter definitions
!       for non-integer constants are used throughout the code to facilitate
!       changes to the default floating point precision if required:
!
!       zero      = 0.0
!       one       = 1.0
!       two       = 2.0
!       three     = 3.0
!       point_5   = 0.5
!       point_75  = 0.75
!
!       pi                 = 3.14159265
!       degrees_to_radians = pi / 180.0
!
!       toa_pressure = 0.005
!
!
!       pseudo-parameters
!       -----------------
!
!       these values are not parameters in the fortran sense in that they are
!       defined at run-time based on user inputs but once defined, they are
!       (or should be) invariant.
!
!       max_n_channels:              integer defining the maximum number of
!                                    instrument channels. this defines the
!                                    valid channels for the valid satellites
!                                    that the user has selected.
!                                    upon rtm initialisation and destruction
!                                    the value is set to -1.
!                                    the value of max_n_channels can only be
!                                    accessed through its own methods:
!                                      - set_max_n_channels()
!                                      - reset_max_n_channels()
!                                      - get_max_n_channels()
!       
!
! modules:
!       none.
!
! contains:
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
!       none.
!
! creation history:
!       written by:     paul van delst, cimss@noaa/ncep 31-jul-2000
!                       pvandelst@ncep.noaa.gov
!
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

module parameters

  ! ---------------------
  ! module use statements
  ! ---------------------

  use type_kinds, only : fp_kind


  ! ------------------
  ! default visibility
  ! ------------------

  private


  ! ---------------------------
  ! current number of absorbers
  ! ---------------------------

  integer, public, parameter :: max_n_absorbers = 3


  ! --------------------
  ! number of predictors.
  ! --------------------

  integer, public, parameter :: max_n_predictors_used = 5

  integer, public, parameter :: max_n_standard_predictors   = 9
  integer, public, parameter :: max_n_integrated_predictors = 6

  integer, public, parameter :: max_n_predictors = max_n_standard_predictors + &
                                                   ( max_n_absorbers * max_n_integrated_predictors )



  ! --------------------------------------
  ! number of absorber layers in algorithm
  ! --------------------------------------

  integer, public, parameter :: max_n_absorber_layers = 300


  ! ----------------------------------------------------------
  ! number of channels (for all satellites - really the number
  ! of satellites x number of channels used per satellite)
  !
  ! this is also the number of lines in the satellite 
  ! information file.
  !
  ! eventually the max_n_profiles and max_n_layers values
  ! will be dynamic, i.e they will be defined by user inputs.
  ! for now, however, they're hardwired.
  ! --------------------------------------------------------

  ! -- accessed via set_<name>, reset_<name>, and get_<name> routines
  integer, private, parameter :: reset_value = -1
  integer, private, save      :: max_n_channels = reset_value

  integer, public, parameter :: max_n_profiles = 128
  integer, public, parameter :: max_n_layers   = 100
!!!  integer, private :: max_n_profiles = reset_value
!!!  integer, private :: max_n_layers   = reset_value


  ! -----
  ! flags
  ! -----

  ! -- direction flags for transmittance calculation
  integer, public, parameter :: down = 0
  integer, public, parameter :: up   = 1


  ! --------------------
  ! numerical parameters
  ! --------------------

  ! -- numbers
  real( fp_kind ), public, parameter :: zero      = 0.0_fp_kind
  real( fp_kind ), public, parameter :: one       = 1.0_fp_kind
  real( fp_kind ), public, parameter :: two       = 2.0_fp_kind
  real( fp_kind ), public, parameter :: three     = 3.0_fp_kind
  real( fp_kind ), public, parameter :: point_5   = 0.5_fp_kind
  real( fp_kind ), public, parameter :: point_75  = 0.75_fp_kind

  ! -- precision/tolerance
  real( fp_kind ), public, parameter :: tolerance = epsilon( one )

  ! -- constant to allow degrees->radians conversion
  real( fp_kind ), public, parameter :: pi = 3.14159265358979323_fp_kind
  real( fp_kind ), public, parameter :: degrees_to_radians = pi / 180.0_fp_kind

  ! -- top-of-atmosphere pressure in hpa
  real( fp_kind ), public, parameter :: toa_pressure = 0.005_fp_kind

  ! -- reciprocal gravity (scaled by 100 for use with pressure in hpa)
  real( fp_kind ), public, parameter :: reciprocal_gravity = one / 980.665_fp_kind

  ! -- diffusivity angle secant = acos( 3/5 ) in degrees (~53.13)
  real( fp_kind ), public, parameter :: secant_diffusivity_angle = 5.0_fp_kind / three

  ! -- maximum solar zenith angle secant definition. should be determined
  ! -- by the maximum angle secant used in generating the transmittance
  ! -- model coefficients, i.e. a secant of 2.25 => 63.6deg. users have
  ! -- requested the value be 85deg => secant of ~11.47.
  real( fp_kind ), public, parameter :: max_solar_angle = 85.0_fp_kind
  real( fp_kind ), public, parameter :: max_secant_solar_angle = 11.473711738554476_fp_kind

!!!! the following is preferred but not allowed in fortran 90          !!!!
!!!! use of non-integer result intrinsics allowed in fortran 95 though !!!!
!  real( fp_kind ), private, parameter :: max_solar_angle = 85.0_fp_kind
!  real( fp_kind ), public,  parameter :: max_secant_solar_angle = one / cos( degrees_to_radians * max_solar_angle )
                                                                 

  ! ---------------------
  ! subprogram visibility
  ! ---------------------

  public :: set_max_n_channels, reset_max_n_channels, get_max_n_channels


contains


  ! -------------------------------------------
  ! subroutines to set and get the value of the 
  ! "pseudo-parameter" max_n_channels
  ! -------------------------------------------

  ! -- set the value
  subroutine set_max_n_channels( value )
    integer, intent( in ) :: value
    max_n_channels = value
  end subroutine set_max_n_channels

  ! -- reset the value
  subroutine reset_max_n_channels()
    max_n_channels = reset_value
  end subroutine reset_max_n_channels

  ! -- get the value and test if it's been set
  subroutine get_max_n_channels( value, is_set )
    integer, intent( out )           :: value
    logical, intent( out ), optional :: is_set
    value = max_n_channels
    if ( present( is_set ) ) then
      is_set = .false.
      if ( value /= reset_value ) is_set = .true.
    end if
  end subroutine get_max_n_channels

end module parameters


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
! revision 1.8  2001/08/31 21:09:45  paulv
! - added the secant of the maximum solar angle as a parameter. the secant
!   value was calculated and expressed as an explicit number since fortran 90
!   does not allow intrinsics that have other than integer results in a
!   parameter initialisation expression.
! - changed  max_solar_zenith_angle name to max_solar_angle.
!
! revision 1.7  2001/08/16 16:44:14  paulv
! - updated documentation.
! - changed max_n_channels attributes from public to private, save. the value
!   of max_n_channels is now accessed via its public methods:
!     set_max_n_channels
!     reset_max_n_channels()
!    get_max_n_channels()
! - added reset_value parameter.
! - removed point_333 parameter.
!
! revision 1.6  2001/07/12 16:46:12  paulv
! - added private statement to prevent definitions in module type_kinds
!   being available outside the scope of this module.
!
! revision 1.5  2001/05/29 17:42:55  paulv
! - now use type_kinds module parameter fp_kind to set the floating point
!   data type. all real declarations are now typed with fp_kind.
! - added direction flags for transmittance calculation.
!
! revision 1.4  2000/11/09 20:32:11  paulv
! - removed max_n_channels as a parameter. it is now a "pseudo" parameter
!   in that it is determined by the number of channels for which coefficients
!   are defined.
! - added some more numerical parameters.
!
! revision 1.3  2000/08/31 19:36:33  paulv
! - added documentation delimiters.
! - updated documentation headers.
!
! revision 1.2  2000/08/24 15:39:45  paulv
! - changed the parameter name that references how many predictors of the
!   total set to use from max_n_predictors_to_use to max_n_predictors_used.
!   i felt this would clarify (for me at least) that while the maximum
!   number of predictors is set, the number that is actually used can be
!   less than that.
! - current maximum number of layers is 100. this is a temporary limit for
!   testing purposes.
! - the parameter reciprocal_gravity was removed from this module and placed
!   in the absorber_profile module where it is used.
!
! revision 1.1  2000/08/08 16:57:21  paulv
! initial checkin
!
!
!
!
