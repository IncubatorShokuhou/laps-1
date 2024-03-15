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
        subroutine plot_winds_2d(u,v,interval,size_in,zoom
     1          ,imax,jmax,lat,lon,r_missing_data,namelist_parms)

        include 'trigd.inc'

        include 'lapsplot.inc'

        real u(imax,jmax),v(imax,jmax)
        real lat(imax,jmax),lon(imax,jmax)
        real mspkt

        logical l_barbs

        data mspkt/.518/

        l_barbs = .true.

!       this variable keeps the barbs away from the boundary
        isize = 0 ! interval + 1

        relsize = size_in

        write(6,*)
        write(6,*) ' plot_winds_2d: interval/size=',interval,relsize
        write(6,*)
        write(6,*) ' winds are assumed to be grid north at this point'       

        if(namelist_parms%l_sphere)then
            write(6,*) ' aspect ratio is non-unity for spherical proj'       
        endif

        do j = 1+isize,jmax-isize,interval

          aspect = 1.0

!         adjust barb spacing for spherical projection (by powers of two)
!         it is assumed a 'latlon' grid is being used

          if(namelist_parms%l_sphere)then
            arg = cosd(lat(1,j))
            ratio_log = nint(log(arg) / log(0.5))
            projfrac = 0.5**ratio_log

            interval_i = nint(float(interval) / projfrac)

            if(arg .gt. 0.)then
                aspect = 1.0 / arg
            endif

            write(6,*)'j/intvl/aspect/lat=',j,interval_i,aspect,lat(1,j)       

          else
            interval_i = interval

          endif

          do i = 1+isize,imax-isize,interval_i

            alat = lat(i,j)
            alon = lon(i,j)

            if( u(i,j) .ne. r_missing_data
     1    .and. v(i,j) .ne. r_missing_data
     1    .and. abs(u(i,j)) .lt. 1e6               ! deals with old data
     1    .and. abs(v(i,j)) .lt. 1e6
     1    .and. aspect      .le. 10.               ! cap on aspect ratio
     1                                                  )then

                call         uv_to_disp(u(i,j),
     1                          v(i,j),
     1                          dir,
     1                          speed)
                spd_kt = speed / mspkt
                call latlon_to_rlapsgrid(alat,alon,lat,lon,imax,jmax
     1                                                  ,ri,rj,istatus)

                if(l_barbs)then 
                    call plot_windob(dir,spd_kt,ri,rj,lat,lon,imax,jmax
     1                              ,relsize,aspect,'grid')

                else ! plot wind arrows
!                   call plot_windarrow()

                endif

            endif


          enddo ! i
        enddo ! j

        return
        end
