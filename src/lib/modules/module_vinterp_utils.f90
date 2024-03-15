!dis   
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis    
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis   
!dis 

module vinterp_utils

  ! module that contains arrays and routines for linear and logarithmic
  ! interpolation 

  
  implicit none
  
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compute_lin_weights(p_mid_mb,p_lower_mb,p_upper_mb, &
                                  weight_bot, weight_top)

    ! computes the weighting coefficient for the top bounding level
    ! of two pressure levels to interpolate to a level in-between using
    ! linear interpolation.
   
    implicit none
    real, intent(in)             :: p_mid_mb  ! desired pressure level
    real, intent(in)             :: p_lower_mb ! lower bounding pressure
    real, intent(in)             :: p_upper_mb ! upper bounding pressure
    real, intent(out)            :: weight_bot ! weight given to bottom level
    real, intent(out)            :: weight_top ! weight given to top level

    weight_bot = (p_mid_mb - p_upper_mb) / (p_lower_mb - p_upper_mb) 
    weight_top = 1.0 - weight_bot
    
  end subroutine compute_lin_weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine compute_log_weights(p_mid_mb,p_lower_mb,p_upper_mb, &
                                 weight_bot, weight_top)
    
    ! computes weighting coefficient for upper pressure level that bounds
    ! a desired pressure level for logarithmic interpolation.
    
    ! note: pressures must be in mb!!!
    
    implicit none
    real, intent(in)             :: p_mid_mb  ! desired pressure level
    real, intent(in)             :: p_lower_mb ! lower bounding pressure
    real, intent(in)             :: p_upper_mb ! upper bounding pressure
    real, intent(out)            :: weight_bot ! weight given to bottom level
    real, intent(out)            :: weight_top ! weight given to top level
    
    weight_bot = ( alog(p_mid_mb) - alog(p_upper_mb) ) / &
                 ( alog(p_lower_mb) - alog(p_upper_mb) ) 
    weight_top = 1.0 - weight_bot

  end subroutine compute_log_weights
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine vinterp_3d(var_in, trap_bot_ind, trap_top_ind, &
                        weight_top, below_ground_value, var_out, &
                        nx, ny, nzin, nzout)

    ! interpolates var_in to var_out using array of trapping indices
    ! and top weight values.  allows the user to set a default
    ! 2d array to use for "below ground" values

    implicit none
    integer, intent(in)            :: nx
    integer, intent(in)            :: ny
    integer, intent(in)            :: nzin
    integer, intent(in)            :: nzout
    real, intent(in)               :: var_in (nx,ny,nzin)
    integer, intent(in)            :: trap_bot_ind(nx,ny,nzout)
    integer, intent(in)            :: trap_top_ind(nx,ny,nzout)
    real, intent(in)               :: weight_top(nx,ny,nzout)
    real, intent(in)               :: below_ground_value(nx,ny)
    real, intent(out)              :: var_out(nx,ny,nzout)
    integer                        :: i,j,k

    do j = 1,ny 
      do i = 1,nx
        do k= 1,nzout

          ! is level below ground?  if so, zero it out
          if (trap_bot_ind(i,j,k).lt. 1) then
            var_out(i,j,k) = below_ground_value(i,j)
          
          ! is it above model top? if so, replicate top value
          else if (trap_bot_ind(i,j,k).eq.nzin) then
            var_out(i,j,k) = var_in(i,j,nzin)

          else
            var_out(i,j,k) = weight_top(i,j,k) * &
                             var_in(i,j,trap_top_ind(i,j,k)) + &
                             (1.-weight_top(i,j,k)) * &
                             var_in(i,j,trap_bot_ind(i,j,k))
          endif
        enddo
      enddo
    enddo
 
    return
  end subroutine vinterp_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module vinterp_utils  
