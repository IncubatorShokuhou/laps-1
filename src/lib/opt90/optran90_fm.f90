!------------------------------------------------------------------------------
!
! name:
!       optran90_fm ! initial interface to f90 code to the laps system
!                     built on optran90 test code
!                     and paul van delst (acknowledged below)
!                     forecast systems laboratory
!
!  copywrite (c) 2003 daniel birkenheuer
!
!
!
! creation history:
!       written by:     paul van delst, cimss/ssec 20-aug-2001
!                       paul.vandelst@ssec.wisc.edu
!                       modified by dan birkenheuer
!                       birk@fsl.noaa.gov for interface to laps
!                       used in laps per permission of paul van delst
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
!
!------------------------------------------------------------------------------


  subroutine optran90_fm ( kk, &
          lev_p, lay_p, lay_t, lay_w, lay_o, &
          laps_sfc_t, &
          sfc_emis, &
          sfc_refl, &
          sec_za, &
          sec_solar, &
          mchan, &
          tau90, &
          flx_tau, &
          sol_tau, &
          up_radiance, &
          bright_temp, & 
          sndr_coeff,&
          sndr_trans, &
          sndr_coeff_len, &
          sndr_trans_len &
          )


  !#----------------------------------------------------------------------------#
  !#                             -- module usage --                             #
  !#----------------------------------------------------------------------------#

  use type_kinds
  use file_utility
  use error_handler
  use parameters
  use initialize
  use forward_model, only: compute_rtm

  ! -- for access to is_microwave_channel only.
  use spectral_coefficients



  !#----------------------------------------------------------------------------#
  !#                           -- type declarations --                          #
  !#----------------------------------------------------------------------------#

  ! ---------------------------
  ! disable all implicit typing
  ! ---------------------------

  implicit none

  integer kk
  integer mchan

  real tau90 (kk,mchan)
  real sec_solar
  real flx_tau (kk,mchan)
  real bright_temp (mchan)
  real up_radiance (mchan)
  real sol_tau (kk,mchan)
  real sec_za
  real lay_w (kk)
  real lay_t (kk)
  real lev_p (kk+1)
  real lay_p (kk)
  real lay_o (kk)
  real sfc_refl
  real sfc_emis
  real laps_sfc_t

! optran 90 coefficient data information

  character*256 sndr_coeff, sndr_trans
  integer sndr_coeff_len, sndr_trans_len



  ! ------------------
  ! program parameters
  ! ------------------

  ! path for laps

  character*256 fname
  integer len
  integer channels_used


  ! -- name
  character( * ),  parameter :: program_name = 'rtm_test_fwd'

  ! -- angle definitions
  integer,         parameter :: n_angles  = 1
  real( fp_kind ), parameter :: min_angle =  0.0_fp_kind
  real( fp_kind ), parameter :: max_angle = 0.0_fp_kind

  ! -- profile definitions
  character( * ),  parameter :: profile_file        = 'profiles_layer.bin'
  integer,         parameter :: n_profiles          = 1
  integer n_layers    ! arbitrary level number to be kk


  ! -- other dimension parameters
  integer,         parameter :: n_absorbers  = max_n_absorbers
  integer,         parameter :: n_predictors = max_n_predictors

  ! -- emissivity parameters
  real( fp_kind ), parameter :: default_mw_emissivity = 0.6_fp_kind
  real( fp_kind ), parameter :: default_ir_emissivity = 0.96_fp_kind

  ! -- default solar angle secant (> 11.47 means no solar)
  real( fp_kind ), parameter :: default_secant_solar_angle = 12.0_fp_kind


  ! -----------------
  ! program variables
  ! -----------------

  ! -- error message
  character( 128 ) :: message 

  ! -- status variables
  integer :: error_status
  integer :: allocate_status
  integer :: io_status

  ! -- loop counters, dimensions and file lun
  integer :: i, j, k, l, m, il
  integer :: n_available_channels, l1, l2, begin_channel, end_channel, n_channels
  integer :: record_length, record_number
  integer :: lun
  integer :: lun_fwd

  ! -- secant angle arrays
  integer         :: i_angle
  real( fp_kind ) :: d_angle
  real( fp_kind ), dimension( n_angles ) :: view_angle
  real( fp_kind ), dimension( n_angles ) :: secant_view_angle
  real( fp_kind ), dimension( n_angles ) :: secant_solar_angle
  real( fp_kind ), dimension( n_angles ) :: surface_temperature

  ! -- profile read array
  real( fp_kind ), dimension( : ), allocatable  :: level_pressure
  real( fp_kind ), dimension( : ), allocatable  :: layer_pressure
  real( fp_kind ), dimension( : ), allocatable  :: layer_temperature
  real( fp_kind ), dimension( : ), allocatable  :: layer_water_vapor
  real( fp_kind ), dimension( : ), allocatable  :: layer_ozone



  ! -- profile data arrays
  real( fp_kind ), dimension( :, : ), allocatable :: level_p
  real( fp_kind ), dimension( :, : ), allocatable :: layer_p
  real( fp_kind ), dimension( :, : ), allocatable :: layer_t
  real( fp_kind ), dimension( :, : ), allocatable :: layer_w
  real( fp_kind ), dimension( :, : ), allocatable :: layer_o
 
  ! -- number of channels processed for each profile
  integer, dimension( n_angles ) :: n_channels_per_profile

  ! -- allocatable arrays
  real( fp_kind ), dimension( : ),    allocatable :: surface_emissivity
  real( fp_kind ), dimension( : ),    allocatable :: surface_reflectivity
  integer,         dimension( : ),    allocatable :: channel_index

  real( fp_kind ), dimension( :, : ), allocatable :: tau
  real( fp_kind ), dimension( :, : ), allocatable :: flux_tau
  real( fp_kind ), dimension( :, : ), allocatable :: solar_tau

  real( fp_kind ), dimension( : ),    allocatable :: upwelling_radiance
  real( fp_kind ), dimension( : ),    allocatable :: brightness_temperature

  ! -- stuff for emissivities
  real( fp_kind ) :: angle_modifier
  real( fp_kind ), dimension( n_angles ) :: mw_emissivity
  real( fp_kind ), dimension( n_angles ) :: ir_emissivity
  integer ::  first_time = 1


  ! ----------
  ! intrinsics
  ! ----------

  intrinsic cos, &
            min, max, &
            real, &
            trim



  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#                         -- initialize the rtm --                           #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  n_layers = kk  ! profile dependent length


  if (first_time == 1) then ! first time it is called
     first_time = 0

     call get_directory ('static',fname,len)
     fname = fname(1:len)//'optranlib/'
     len = len_trim(fname)
     
     error_status = initialize_rtm( tau_file = sndr_trans(1:sndr_trans_len),&
                   path = fname(1:len),  &
               spectral_file = sndr_coeff(1:sndr_coeff_len)  )

     if ( error_status /= success ) then
        call display_message( program_name, &
             'error initialzing rtm', &
             error_status )
        stop
     end if

  endif ! first_time called


  !##############################################################################
  !##############################################################################
  !##############################################################################
  !#                                                                            #
  !#               -- setup data and arrays for calculations --                 #
  !#                                                                            #
  !##############################################################################
  !##############################################################################
  !##############################################################################

  !#----------------------------------------------------------------------------#
  !#                    -- allocate arrays for rtm model --                     #
  !#----------------------------------------------------------------------------#

  ! ----------------------------------
  ! get the number of channels read in
  ! from the coefficient data files
  ! ----------------------------------

! allocate arrays for variable profile contents
!

  allocate (      level_p                  ( n_layers,n_angles), &
                  layer_p                  ( n_layers,n_angles), &    
                  layer_t                  ( n_layers,n_angles), & 
                  layer_w                  ( n_layers,n_angles), & 
                  layer_o                  ( n_layers,n_angles), &
                  level_pressure           ( n_layers), &
                  layer_pressure           ( n_layers), & 
                  layer_temperature        ( n_layers), & 
                  layer_water_vapor        ( n_layers), & 
                  layer_ozone              ( n_layers)   )

  call get_max_n_channels( n_channels )

  begin_channel = 1
  end_channel   = n_channels


  ! ---------------------------------------------------
  ! allocate the forward model channel dependent arrays
  ! ---------------------------------------------------

  allocate( surface_emissivity(   n_channels * n_angles ),   &  ! input,  l*m
            surface_reflectivity( n_channels * n_angles ),   &  ! input,  l*m
            channel_index(        n_channels * n_angles ),   &  ! input,  l*m

            tau(       n_layers, n_channels * n_angles ),    &  ! output, k x l*m
            flux_tau(  n_layers, n_channels * n_angles ),    &  ! output, k x l*m
            solar_tau( n_layers, n_channels * n_angles ),    &  ! output, k x l*m

            upwelling_radiance(     n_channels * n_angles ), &  ! output, l*m
            brightness_temperature( n_channels * n_angles ), &  ! output, l*m

            stat = allocate_status )

  if ( allocate_status /= 0 ) then
    write( message, '( "error allocating forward model channel ", &
                      &"dependent arrays. stat = ", i5 )' ) &
                    allocate_status
    call display_message( program_name,    &
                          trim( message ), &
                          failure          )
    stop
  end if





    !#--------------------------------------------------------------------------#
    !#                           -- assign variable from laps input             #
    !#--------------------------------------------------------------------------#

      do k = 1, kk
         level_p(k, 1 ) = lev_p(k+1)
      enddo
      layer_p(1:kk, 1 ) = lay_p(1:kk)
      layer_t(1:kk, 1 ) = lay_t(1:kk)
      layer_w(1:kk, 1 ) = lay_w(1:kk)
      layer_o(1:kk, 1 ) = lay_o(1:kk)
      surface_temperature = laps_sfc_t
      surface_emissivity = sfc_emis    !array assignment to constant
      surface_reflectivity = sfc_refl  !array assignment to constant
      secant_view_angle = sec_za
      secant_solar_angle = sec_solar
      n_channels_per_profile  = n_channels  ! don't understand

      do k = 1,n_channels
         channel_index(k) = k  ! don't understand
      enddo

      ! note the (in,out) cautions on upwelling radiance and brightness_temperature
      ! this is preset here.

      upwelling_radiance = 0.
      brightness_temperature = 1.0  ! recommended preset


   

    !#--------------------------------------------------------------------------#
    !#                            -- forward model --                           #
    !#--------------------------------------------------------------------------#


    error_status = compute_rtm( &
                                ! -- forward inputs
                                level_p,                &  ! input,  k x m
                                layer_p,                &  ! input,  k x m
                                layer_t,                &  ! input,  k x m
                                layer_w,                &  ! input,  k x m
                                layer_o,                &  ! input,  k x m

!                                layer_t(n_layers,:),    &  ! input, m
                                surface_temperature,    &  ! input, m
                                surface_emissivity,     &  ! input, l*m
                                surface_reflectivity,   &  ! input, l*m

                                secant_view_angle,      &  ! input, m
                                secant_solar_angle,     &  ! input, m
                                n_channels_per_profile, &  ! input, m
                                channel_index,          &  ! input, l*m

                                ! -- forward outputs
                                tau,                    &  ! output, k x l*m
                                flux_tau,               &  ! output, k x l*m
                                solar_tau,              &  ! output, k x l*m

                                upwelling_radiance,     &  ! output, l*m
                                brightness_temperature  )  ! output, l*m

    if ( error_status /= success ) then
      call display_message( program_name, &
                            'error occured in compute_rtm', &
                            error_status )
      stop
    end if

! here the returned arrays (upwelling_radiance and brightness_temperature)
! have dimensions of the insturment channels.  the laps array however is a 
! fixed one of 18.  this takes and only uses the part of the laps array needed

    
    channels_used = size (upwelling_radiance)

    do i = 1,channels_used
       up_radiance(i) = upwelling_radiance(i)
       bright_temp(i) = brightness_temperature(i)
    enddo




  !#----------------------------------------------------------------------------#
  !#                           -- destroy the rtm --                            #
  !#----------------------------------------------------------------------------#

  ! ---------------------------------------
  ! deallocate the channel dependent arrays
  ! ---------------------------------------

  ! -- forward model

  deallocate (  layer_p,  layer_t  ,layer_w  ,layer_o ,  layer_pressure,  layer_temperature , layer_water_vapor ,  layer_ozone,  &
                level_p, level_pressure  )

  deallocate( surface_emissivity,     &  ! input,  l*m
              surface_reflectivity,   &  ! input,  l*m
              channel_index,          &  ! input,  l*m

              tau,                    &  ! output, k x l*m
              flux_tau,               &  ! output, k x l*m
              solar_tau,              &  ! output, k x l*m

              upwelling_radiance,     &  ! output, l*m
              brightness_temperature, &  ! output, l*m

              stat = allocate_status  )

  if ( allocate_status /= 0 ) then
    write( message, '( "error deallocating forward model channel ", &
                      &"dependent arrays. stat = ", i5 )' ) &
                    allocate_status
    call display_message( program_name,    &
                          trim( message ), &
                          warning          )
  end if


  ! ---------------------------------
  ! deallocate the coefficient arrays
  ! ---------------------------------

 ! this utility is now handled by the new module below.  called from 
 ! variational.f after all calls to the forward model are complete.

end subroutine  optran90_fm

subroutine optran_deallocate (istatus)

!  use type_kinds
!  use file_utility
!  use error_handler
!  use parameters

  use initialize
  integer :: istatus
  
  istatus = destroy_rtm()
  
!  if ( error_status /= success ) then
!     call display_message( program_name, &
!          'error destroying rtm', &
!          error_status )
!     stop
!  end if
  

end subroutine optran_deallocate


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
! revision 1.3  2003/07/10 17:09:24  birk
! ready for laps optran90
!
! revision 1.2  2002/11/18 20:01:39  birk
! changes made to avoid compilation errors on jet, statement order specifics.
!
! revision 1.1  2002/11/15 15:21:32  birk
! added to cvs mainly to see how this compiles on other platforms, it currently
! seems to compile on the ibm
!
! revision 1.1  2001/09/13 22:08:13  paulv
! initial checkin.
!
!
!
!
