!====================================================================
!>
!! crtm amsu-b module for stmas assimilation
!!
!! \history \n
!! creation: yuanfu xie 2013
!
!====================================================================
module crtm_kmatrix

  use crtm_module
  use prmtrs_stmas_cloud, only : satellite_obs

  implicit none

  character(*), parameter :: program_name = 'lvd_da'

  ! parameters which could be ported to sat_da.nl later.
  integer, parameter :: num_lvd_sensor = 1
  integer, parameter :: num_lvd_absorber = 2
  integer, parameter :: num_lvd_cloud = 0
  integer, parameter :: num_lvd_aerosol = 0
  real(fp),    parameter :: lvd_zenith_angle=-34.04_fp, lvd_scan_angle=-30.20_fp

  integer :: num_lvd_channels

  type(crtm_channelinfo_type), allocatable:: lvd_channelinfo(:)
  type(crtm_geometry_type),    allocatable:: lvd_geometryinfo(:)
  type(crtm_atmosphere_type),  allocatable:: lvd_atm(:)     ! profiles
  type(crtm_atmosphere_type),  allocatable:: lvd_atm_k(:,:) ! channel x profiles
  type(crtm_surface_type),     allocatable:: lvd_sfc(:)     ! profiles
  type(crtm_surface_type),     allocatable:: lvd_sfc_k(:,:) ! channel x profile
  ! rtsolution and rtsolution_k: channels x profiles
  type(crtm_rtsolution_type),  allocatable:: lvd_rtsolution(:,:), &
                                             lvd_rtsolution_k(:,:)
contains

!====================================================================
!> 
!! crtm configuration for goes.
!!
!! \history \n
!! creation: yuanfu xie 2013, s albers (modified for goes)
!
!====================================================================

subroutine conf4lvd(lvd,numgrid,sensor_id,coe_path,moist_unit,&
                      badflag)

  implicit none

  character(*), intent(in) :: sensor_id,coe_path,moist_unit
  integer, intent(in) :: numgrid(4)
  real,    intent(in) :: badflag
  type(satellite_obs), intent(in) :: lvd
  real(fp) :: angle(4,lvd%numpro)

  ! local variables:
  character(256) :: message
  integer :: istatus,np

  ! allocate memory for amsu-b assimilation:
  allocate(lvd_channelinfo(num_lvd_sensor), &
           lvd_geometryinfo(lvd%numpro), &
           stat=istatus)

  istatus = crtm_init( (/sensor_id/), lvd_channelinfo, &
                        emiscoeff_file='wu-smith.emiscoeff.bin', &
                        file_path=trim(coe_path)//'coefficient_data/', &
                        quiet = .true.)  
  if ( istatus /= success ) then
    message = 'error initializing crtm'
    call display_message( program_name, message, failure )
    stop
  else
    message = 'crtm initialized'
    call display_message( program_name, message, success )
  endif
  num_lvd_channels = sum(lvd_channelinfo%n_channels)
  if (num_lvd_channels .ne. lvd%numchs) then
    print*,'numbers of amsub channels do not match!'
    stop
  endif

  ! allocate memory for crtm arrays:
  allocate(lvd_rtsolution(lvd%numchs,lvd%numpro), &
           lvd_rtsolution_k(lvd%numchs,lvd%numpro), &
           lvd_atm(lvd%numpro), &
           lvd_sfc(lvd%numpro), &
           lvd_atm_k(lvd%numchs,lvd%numpro), &
           lvd_sfc_k(lvd%numchs,lvd%numpro), &
           stat=istatus)

  ! create crtm:
  call crtm_atmosphere_create(lvd_atm,numgrid(3),num_lvd_absorber, &
                              num_lvd_cloud,num_lvd_aerosol)
  if ( any(.not. crtm_atmosphere_associated(lvd_atm)) ) then
     message = 'error allocating crtm atmosphere structures'
     call display_message( program_name, message, failure )
     stop
  endif
  call crtm_atmosphere_create(lvd_atm_k,numgrid(3),num_lvd_absorber, &
                              num_lvd_cloud,num_lvd_aerosol)
  if ( any(.not. crtm_atmosphere_associated(lvd_atm)) ) then
     message = 'error allocating crtm atmosphere_k structures'
     call display_message( program_name, message, failure )
     stop
  endif

  ! data type conversion: real to real(fp) of crtm:
  angle = lvd%angles

  ! geometry setting:
  call crtm_geometry_setvalue( lvd_geometryinfo, &
                               ifov = lvd%numpro, &
                               sensor_zenith_angle  = angle(1,1:lvd%numpro), & 
                               sensor_azimuth_angle = angle(2,1:lvd%numpro), &
                               source_zenith_angle  = angle(3,1:lvd%numpro), &
                               source_azimuth_angle = angle(4,1:lvd%numpro) ) 

  ! pass the non-control profile values to crtm:
  do np=1,lvd%numpro

    if (lvd%lndcvr(np) .gt. 0) then
      lvd_sfc(np)%land_coverage  = 1.0_fp
      lvd_sfc(np)%water_coverage  = 0.0_fp
    else
      lvd_sfc(np)%land_coverage  = 0.0_fp
      lvd_sfc(np)%water_coverage  = 1.0_fp
    endif
    lvd_sfc(np)%land_type      = 1 ! compacted soil, according crtm 2.1.1
    lvd_sfc(np)%water_type     = 1 ! sea_water
    lvd_atm(np)%absorber_id    = (/ h2o_id , o3_id /)
    lvd_atm(np)%absorber(:,2)  = 2.53e-3  ! fake value as mhs does not use o3

    ! currently, use mass mixing ratio only crtm user guide 2.1.1 section 4.6:
    if (moist_unit .eq. 'mass_mixing_ratio_units') then
      lvd_atm(np)%absorber_units = (/ mass_mixing_ratio_units,  &
                                        volume_mixing_ratio_units /)      
    else
      print*,'not implemented, check the namelist, sat_da.nl'
      stop
    endif
    
  enddo

end subroutine conf4lvd

!====================================================================
!>
!! crtm forward and k-matrix calculation for amsu-b.
!!
!! \history \n
!! creation: yuanfu xie 2013
!
!====================================================================

subroutine lvd_costgrad(numgrid,u_sfc,v_sfc,pres,sh,temp,lvd,tf, &
                          badflag,cost,grad,lvd_channels_used, &
                          iu,iv,it,iq)

  implicit none

  type(satellite_obs), intent(in) :: lvd

  integer, intent(in) :: numgrid(4),lvd_channels_used(*)
  integer, intent(in) :: tf,iu,iv,it,iq  ! time frame and indices for control variables
  real,    intent(in) :: pres(numgrid(3))
  real,    intent(in) ::   sh(lvd%numpro,numgrid(3))
  real,    intent(in) :: temp(lvd%numpro,numgrid(3))
  real,    intent(in) :: u_sfc(lvd%numpro)
  real,    intent(in) :: v_sfc(lvd%numpro)
  real,    intent(in) :: badflag
  real,   intent(out) :: cost,grad(numgrid(1),numgrid(2),numgrid(3),numgrid(4),*)

  ! local variables:
  real, parameter :: pi=3.1415926
  character(256) :: message
  double precision :: f
  integer :: np,nl,nc,istatus
  real    :: bt,spd,weight

  print*,'amsub cost function and gradient calculation...'

  do np=1,lvd%numpro

    ! 2d fields:
    lvd_sfc(np)%land_temperature = temp(np,1)
    lvd_sfc(np)%wind_speed = sqrt(u_sfc(np)**2+v_sfc(np)**2)
    lvd_sfc(np)%wind_direction = (2.0*pi+asin(v_sfc(np)/lvd_sfc(np)%wind_speed))*180.0/pi

    ! 3d fields:
    ! note: future improvement would be using terrain-following coordinate:
    ! currently, we do not treat terrain!!!
    lvd_atm(np)%level_pressure(0) = 1.5*pres(numgrid(3))-0.5*pres(numgrid(3)-1)
    do nl=1,numgrid(3)-1
      lvd_atm(np)%level_pressure(nl) = 0.5*(pres(numgrid(3)-nl)+pres(numgrid(3)-nl+1))
      lvd_atm(np)%pressure(nl)       = pres(numgrid(3)-nl+1)
    enddo
    lvd_atm(np)%level_pressure(numgrid(3)) = 1.5*pres(1)-0.5*pres(2)
    lvd_atm(np)%pressure(numgrid(3)) = pres(1)

    ! assign meteorological profiles to crtm structured data:

    ! convert specific humidity to mass mixing ratio:
    ! from the book of "fundamentals of atmospheric modeling by mark z. jacobson,
    ! mass mixing ratio (mmr) = rho_v/rho_d (equation 2.26) and 
    ! specific humidity (sh)  = rho_v/(rho_d+rho_v) (equation 2.27). thus
    ! mmr = sh/(1-sh).
    do nl=1,numgrid(3)
      ! crtm mass mixing ratio is in g/kg: table 4.7 crtm user guide 2.1.1:
      lvd_atm(np)%absorber(nl,1) = 1.0e3*sh(np,numgrid(3)-nl+1)/(1.0e3-sh(np,numgrid(3)-nl+1))
      if (lvd_atm(np)%absorber(nl,1) .lt. 0.0) then
        print*,'negative vapor? ',lvd_atm(np)%absorber(nl,1),nl,np
        lvd_atm(np)%absorber(nl,1) = 0.0
      endif

      lvd_atm(np)%temperature(nl) = temp(np,numgrid(3)-nl+1)
    enddo
    
  enddo

  ! clean up arrays for k-matrix:
  call crtm_atmosphere_zero( lvd_atm_k )
  call crtm_surface_zero( lvd_sfc_k )
  lvd_rtsolution_k%brightness_temperature = one
  lvd_rtsolution_k%radiance = zero

  print*,'calling k_matrix...'

  ! k-matrix:
  istatus = crtm_k_matrix( lvd_atm         , &
                           lvd_sfc         , &
                           lvd_rtsolution_k, &
                           lvd_geometryinfo, &
                           lvd_channelinfo , &
                           lvd_atm_k       , &
                           lvd_sfc_k       , &
                           lvd_rtsolution  )
  if ( istatus /= success ) then
    message = 'error in crtm k_matrix model'
    call display_message( program_name, message, failure )
    stop
  endif

  ! cost and grad:
  f = 0.0d0
  ! break function and gradient for calculating weight:
  weight = 0.0
  do np=1,lvd%numpro

    do nc=1,lvd%numchs
      if (lvd_channels_used(nc) .ne. 1) cycle

      bt = lvd_rtsolution(nc,np)%brightness_temperature

      ! cost:
      f = f+(bt-lvd%satval(nc,np))**2
      weight = weight+0.5
    enddo

  enddo
  ! pass the double precision cost function value to the output variable:
  ! weight = sqrt(weight)

  cost = 1.0*f/weight

  do np=1,lvd%numpro

    do nc=1,lvd%numchs
      if (lvd_channels_used(nc) .ne. 1) cycle

      bt = lvd_rtsolution(nc,np)%brightness_temperature

      ! grad:
      ! surface:
      spd = sqrt(u_sfc(np)**2+v_sfc(np)**2)
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,it) = &
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,it)+(bt-lvd%satval(nc,np))/weight* &
        lvd_sfc_k(nc,np)%land_temperature
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iu) = &
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iu)+(bt-lvd%satval(nc,np))/weight* &
        (lvd_sfc_k(nc,np)%wind_speed*     u_sfc(np)/spd - &
         lvd_sfc_k(nc,np)%wind_direction*sign(1.0,u_sfc(np))*v_sfc(np)/(spd*spd))*180.0/pi
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iu) = &
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iv)+(bt-lvd%satval(nc,np))/weight* &
        (lvd_sfc_k(nc,np)%wind_speed*     v_sfc(np)/spd + &
         lvd_sfc_k(nc,np)%wind_direction*abs(u_sfc(np))/(spd*spd))*180.0/pi

      ! 3d profiles:
      do nl=1,numgrid(3)
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,it) = &
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,it) + &
          (bt-lvd%satval(nc,np))/weight*lvd_atm_k(nc,np)%temperature(nl)
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,iq) = &
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,iq) + &
          (bt-lvd%satval(nc,np))/weight*lvd_atm_k(nc,np)%absorber(nl,1)* &
          1.0e6/(1.0e3-sh(np,numgrid(3)-nl+1))**2
      enddo
    enddo

  enddo

end subroutine lvd_costgrad

!====================================================================
!>
!! destroy the crtm for amsu-b.
!!
!! \history \n
!! creation: yuanfu xie 2013
!
!====================================================================

subroutine dstrylvd

  ! local variables:
  integer :: istatus

  istatus = crtm_destroy(lvd_channelinfo)

end subroutine dstrylvd

end module crtm_kmatrix

