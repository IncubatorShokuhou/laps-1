cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 


      subroutine adjust_heights(temp_3d,heights_3d,ht_500_fg
     1                         ,ni,nj,nk,k_ref,istatus)

!     this routine changes the reference for the hydrostatic integration
!     of temperature into height to the model 500mb ht field.
!     the heights normally are referenced using unreduced absolute pressure
!     at the surface (laps var = 'ps').

!     this is currently being used for experimental purposes only.

!     dec 1994    steve albers

      real temp_3d(ni,nj,nk)    ! input
      real heights_3d(ni,nj,nk) ! input/output
      real ht_500_fg(ni,nj)     ! local

      write(6,*)' subroutine adjust_heights'

      resid_max = 0.

      do i = 1,ni
      do j = 1,nj

!         calculate residual = model 500 ht - first guess height
          residual = ht_500_fg(i,j) - heights_3d(i,j,k_ref)

!         apply residual to column of height values.
          do k = 1,nk
              heights_3d(i,j,k) = heights_3d(i,j,k) + residual
          enddo ! k

          resid_max = max(resid_max,abs(residual))

      enddo ! j
      enddo ! i

      write(6,*)'         maximum height adjustment (m) = ',resid_max

      istatus = 1
      return
      end
