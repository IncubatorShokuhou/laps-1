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


       subroutine contour_settings(a,ni,nj,clow,chigh,cint,zoom,density
     1                            ,scale)       

       real a(ni,nj)

       write(6,*)' subroutine contour_settings...'

       scale_pcp = 1. / ((100./2.54))
       if(abs(scale-scale_pcp) .le. .0001) then ! denom = (in/m) for precip
           clow = 0.0
           chigh = 10.
           cint = -0.01
           write(6,*)' contour_settings for precip: ',cint,scale
           return
       endif

       call get_r_missing_data(r_missing_data,istatus)

       call array_range(a,ni,nj,rmin,rmax,r_missing_data)

       write(6,*)' rmax/rmin prior to scaling = ',rmax,rmin

       rmax = rmax / scale
       rmin = rmin / scale

       write(6,*)' rmax/rmin after scaling = ',rmax,rmin

       range = (rmax-rmin) / (sqrt(zoom) * density)
       call range2cint(range,cint)

       range_eff = rmax-rmin
       call range2cint(range_eff,cint_eff)

       cint_2 = cint * 2.

!      if(cint .ge. 1.0)then
!          if(rmin .ge. 0.)then
!              clow  = int(rmin)/int(cint_2) * int(cint_2)
!              chigh = int(rmax)/int(cint) * int(cint) + int(cint)
!          else ! rmin < 0. (perhaps more generic)

!              this is using an effective cint value that is independent of 
!              zoom/density
               clow =  nint(rmin/cint_eff - 0.5) * cint_eff 
               chigh = nint(rmax/cint_eff - 0.5) * cint_eff + cint_eff

!          endif
!      else
!          clow =  int(rmin/cint) * cint 
!          chigh = int(rmax/cint) * cint + cint
!      endif
 
       write(6,*)' range/zoom/clow/chigh/cint/scale = '
     1            ,range,zoom,clow,chigh,cint,scale      

       return
       end


       subroutine range2cint(range,cint)

       if(range .gt. 20000)then
           cint = 5000.
       elseif(range .gt. 8000)then
           cint = 1000.
       elseif(range .gt. 2000)then
           cint = 400.
       elseif(range .gt. 600)then
           cint = 100.
       elseif(range .gt. 250)then
           cint = 50.
       elseif(range .gt. 200)then
           cint = 20.
       elseif(range .gt. 60)then ! from 60-200, cint = 10 (6  - 20 contours)
           cint = 10.
       elseif(range .gt. 30)then ! from  30-60, cint = 5  (6  - 12 contours)
           cint = 5.
       elseif(range .gt. 6)then  ! from   6-30, cint = 2  (3  - 15 contours)
           cint = 2.
       elseif(range .gt. 3)then  ! from   3- 6, cint = 1  (3  -  6 contours)
           cint = 1.
       elseif(range .gt. 2)then  ! from   2- 3, cint = 1  (4  -  6 contours)
           cint = 0.5
       elseif(range .gt. 0.6)then ! from .6- 2, cint =.25 (2  -  8 contours)
           cint = 0.25
       else ! range < 1          
           cint = 0.1
       endif

       return
       end


        subroutine condition_precip(nx_l,ny_l,c_type,accum_2d,scale
     1                             ,pthr)           

        real accum_2d(nx_l,ny_l), image_cutoff
        character*(*)c_type

        if(.true.)return

        image_cutoff = 0.2

!       eliminate "minor" maxima
        do j = 1,ny_l
        do i = 1,nx_l
            if(c_type(1:2) .eq. 'pa')then
                if(accum_2d(i,j)/scale .lt. 0.005)then
                    accum_2d(i,j) = 0.0
                elseif(accum_2d(i,j)/scale .gt. pthr .and. 
     1                 accum_2d(i,j)/scale .lt. image_cutoff .and.
     1                 c_type .eq. 'pai'                   )then
                    accum_2d(i,j) = image_cutoff*scale
                endif
            else
                if(accum_2d(i,j)/scale .lt. pthr)then
                    accum_2d(i,j) = 0.0
                endif
!               if(accum_2d(i,j)/scale .le. 0.0001)then
!                   accum_2d(i,j) = -1e-6
!               endif
            endif
        enddo ! i
        enddo ! j

        return
        end

        subroutine condition_cape(nx_l,ny_l,c_type,r_missing_data
     1                           ,field_2d)

        real field_2d(nx_l,ny_l), image_cutoff
        character*(*)c_type

        image_cutoff = 300.

        if(c_type .eq. 'pe')then ! contour case
!           change flag value 
            do i = 1,nx_l
            do j = 1,ny_l
                if(field_2d(i,j) .eq. r_missing_data)
     1              field_2d(i,j) = -1.           
            enddo ! j
            enddo ! i

        else                     ! image case
!           change flag value 
            do i = 1,nx_l
            do j = 1,ny_l
                if(field_2d(i,j) .eq. -1.0)then
                    field_2d(i,j) = r_missing_data    ! or could be 0.0
                endif

                if(field_2d(i,j) .gt. 0.0 .and. 
     1             field_2d(i,j) .lt. image_cutoff)then
                    field_2d(i,j) = image_cutoff
                endif

            enddo ! j
            enddo ! i

        endif

        return
        end
