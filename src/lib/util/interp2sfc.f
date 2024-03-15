cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis
      subroutine interp_to_sfc(sfc_2d,field_3d,heights_3d,ni,nj,nk,
     &                         badflag,interp_2d)
c
c=================================================================
c
c     routine to interpolate a 3-d field to a 2-d surface.  the
c     2-d surface can be the actual "surface" (the topography),
c     or any other 2-d surface within the grid.
c
c     based on code in the laps wind analysis, steve albers, fsl.
c
c     original: 04-09-99  peter a. stamus, noaa/fsl
c     changes:
c
c=================================================================
c
      real sfc_2d(ni,nj), field_3d(ni,nj,nk), heights_3d(ni,nj,nk)
      real interp_2d(ni,nj)
      
      logical ltest_vertical_grid
c
      write(6,*)' interpolate 3-d field to 2-d surface (interp_to_sfc)'

c
c..... interpolate from the 3-d grid to the 2-d surface at each point.
c
      do j=1,nj
      do i=1,ni
c
         if(.not. ltest_vertical_grid('sigma_ht'))then
            zlow = height_to_zcoord2(sfc_2d(i,j),heights_3d,ni,nj,nk,       
     &                                                  i,j,istatus)
            if(istatus .ne. 1)then
               write(6,*) 
     &             ' error in height_to_zcoord2 in interp_to_sfc',
     &                 istatus
               write(6,*) i,j,zlow,sfc_2d(i,j),(heights_3d(i,j,k)
     1                                          ,k=1,nk)
               return
            endif

         elseif(.true.)then ! sigma_ht grid
            zlow = 1.       ! surface defined as lowest sigma_ht level

         else               ! .false. 
            zlow = rlevel_of_field(sfc_2d(i,j),heights_3d,ni,nj,nk,       
     &                                                  i,j,istatus)
            if(istatus .ne. 1)then
               write(6,*) 
     &             ' error in rlevel_of_field in interp_to_sfc',
     &                 istatus
               write(6,*) i,j,zlow,sfc_2d(i,j),(heights_3d(i,j,k)
     1                                          ,k=1,nk)
               return
            endif

         endif
c
         klow = max(zlow, 1.)
         khigh = klow + 1
         fraclow = float(khigh) - zlow
         frachigh = 1.0 - fraclow
c
         if( field_3d(i,j,klow)  .eq. badflag .or.
     &       field_3d(i,j,khigh) .eq. badflag) then

            write(6,3333)i,j
 3333       format(' warning: cannot interpolate to sfc at ',2i5)
            interp_2d(i,j) = badflag

         else

            interp_2d(i,j) = field_3d(i,j,klow ) * fraclow  +
     &                       field_3d(i,j,khigh) * frachigh

         endif
c
      enddo !i
      enddo !j
c
c..... that's all.
c
      return
      end
