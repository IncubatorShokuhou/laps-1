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
c
c
        subroutine vlog2vis(lvis,vis,ni,nj)
c
c       routine to convert gridded visibility from log( miles ) to
c       miles.  if the analyzed log (vis) is less that -2, set the
c       visibility to zero.
c
c       01-14-92        p. stamus
c
        real vis(ni,nj), lvis(ni,nj)

        vis_scale = 7.
 
        do j=1,nj
        do i=1,ni
          if(lvis(i,j) .le. -1.5 * vis_scale) then
            vis(i,j) = 0.
          else
            vis(i,j) = 10. ** ( lvis(i,j) / vis_scale )
          endif
        enddo !i
        enddo !j
c
        return
        end



        subroutine viss2log(vis_s,mxstn,numsfc,smsng)
c
c       routine to convert visibility from miles to log( miles ) for
c       the individual observations.  if the observed visibility is
c       zero, the set the log( vis ) to -10 for the analysis (so that
c       its a very small number).
c
c       01-14-92        p. stamus
c
        real vis_s(mxstn)
 
        vis_scale = 7.

        do i=1,numsfc
          if(vis_s(i) .lt. 0.) then
            vis_s(i) = smsng
          elseif(vis_s(i) .eq. 0.) then
            vis_s(i) = -10.
          else
            vis_s(i) = log10( vis_s(i) ) * vis_scale
          endif
        enddo !i
c
        return
        end



        subroutine visg2log(vis,ni,nj,smsng)
c
c       routine to convert visibility from miles to log( miles ) for
c       the ni x nj laps grid.  if the visibility is zero, then set
c       the log( vis ) to -10 for the analysis.
c
c       01-14-92        p. stamus
c
        real vis(ni,nj)

        vis_scale = 7.

        do j=1,nj
        do i=1,ni
          if(vis(i,j) .lt. 0.) then
            vis(i,j) = smsng
          elseif(vis(i,j) .eq. 0.) then
            vis(i,j) = -10.
          else
            vis(i,j) = log10( vis(i,j) ) * vis_scale
          endif
        enddo !i
        enddo !j
c
        return
        end


c
        subroutine enhance_vis(i4time,vis,hum,topo,ni,nj,kcloud)
c
c==============================================================================
c
c       routine to call other routines to adjust the visibility analysis
c       based on other data (radar, cloud, etc.).
c            ** may want to put the spline call in here someday...
c
c       original:  ??-??-93  peter a. stamus  noaa/fsl
c       changes:   02-03-94  rewritten
c                  08-26-97  changes for dynamic laps.
c       
c       notes:
c          1.  the variables here are:
c                  i4time = time for this analysis.
c                  vis    = visibility (units are not changed) 
c                  hum    = relative humidity (0 to 100 percent)
c                  topo   = laps topography (meters)
c                  ni,nj  = laps grid dimensions
c                  kcloud = cloud grid dimension in vertical
c          2.  units of visibility (miles, meters) are not changed
c              in this routine.
c
c==============================================================================
c
        real vis(ni,nj), hum(ni,nj), topo(ni,nj)
        real vismod(ni,nj)  !work array
c
        print *,' in enhance_vis routine...'
c
c..... radar adjustment.            ! still disabled...2-3-94 pas (no data)
c
c       call constant(vis,-10.,ni,nj)
c       call get_radar_visibility(i4time,vis,istatus)
c
c
c..... low cloud/humidity adjustment.
c.....................................
c..... get the modification array using the cloud data from lc3 and the
c..... surface relative humidity.  the multiply the visibilities by the
c..... modification factor to get the adjusted visibility.
c
        call get_vismods(i4time,hum,topo,vismod,ni,nj,kcloud)
c
        do j=1,nj
        do i=1,ni
           vis(i,j) = vis(i,j) * vismod(i,j)
        enddo !i
        enddo !j
c
c..... that's it...
c
        print *,' enhance_vis routine...done.'
        return
        end
c
c
        subroutine get_vismods(i4time,hum,topo,vismod,ni,nj,kcloud)
c
c==============================================================================
c
c       routine to set a visibility adjustment based on low cloud data from
c       the laps 3-d cloud analysis and the laps surface humidity analysis.
c       the adjustment is the percentage (0-1) that the visibility is 
c       reduced for given cloud amounts/humidities.  the adjustment is put
c       into the vismod array, and is multiplied by the vis array in the
c       calling routine (enhance_vis).
c
c       original:  ??-??-93  peter a. stamus
c       changes:   02-03-94  rewritten
c                  08-26-97  changes for dynamic laps
c	           08-19-98  initialize k_hold array.
c
c==============================================================================
c
        real hum(ni,nj), vismod(ni,nj), topo(ni,nj)
	real cld_hts(kcloud), cld_pres(kcloud)
        real clouds_3d(ni,nj,kcloud)
c
        integer k_hold(ni,nj), lvl(kcloud)
c
        character ext*31, var(kcloud)*3, comment(kcloud)*125
        character units(kcloud)*10, lvl_coord(kcloud)*4, dir*256
c
c..... start by setting up the default values for vismod (1.0=no adjustment)
c
        call constant(vismod,1.0,ni,nj)
c
c..... get the low cloud data from the nearest lc3 file (timewise).
c
        icnt = 0
        i4time_c = i4time
        do k=1,kcloud
           lvl(k) = k
           var(k) = 'lc3'
        enddo !k
        ext = 'lc3'
	call get_directory('lc3', dir, len)
 500    call read_laps_data(i4time_c,dir,ext,ni,nj,kcloud,kcloud,
     &       var,lvl,lvl_coord,units,comment,clouds_3d,istatus)
c
        if(istatus .ne. 1) then
           if(istatus .eq. 0) then  !no data
              if(icnt .lt. 1) then  !just try 1-hr for now.
           print *,' no data for given i4time...trying 1-hr earlier.'
               i4time_c = i4time_c - 3600
               icnt = icnt + 1
               go to 500
              else
               print *,' lc3 data too old.'
               print *,' no visibility modification done.'
               return
              endif 
           else
              print *,' bad return from lc3 read: istatus = ',istatus
              print *,' no visibility modification done.'
              return
           endif
        endif
c
c..... check time difference.  don't use if cloud analysis is too old.
c
c        if((i4time - i4time_nearest) .gt. 5400) then  !1.5 hours
c           print *,' lc3 data too old.'
c           print *,' no visibility modification done.'
c           return
c        endif
c
c..... decode the cloud heights and pressures.
c
        print *,' got cloud data.'
        do k=1,kcloud
           read(comment(k),100,err=999) cld_hts(k), cld_pres(k)
100        format(2e20.7)
        enddo !k
c
c..... find the vertical level from 'cld_hts' just below the surface.  will
c..... then start cloud checks at the next level up.
c
        do j=1,nj
        do i=1,ni
	   k_hold(i,j) = 0
           do k=1,kcloud
              if(cld_hts(k) .lt. topo(i,j)) k_hold(i,j) = k
           enddo !k
        enddo !i
        enddo !j
c
c..... now loop over the grid and check the lowest 3 levels above the 
c..... surface.  check for fog first, then check for low clouds.
c
        do j=1,nj
        do i=1,ni
           k_start = k_hold(i,j) + 1     ! 1st level above the surface
           k_end = k_start + 3           ! this is usually within 300 m
c
c.....     check for fog in the layer just above the surface.
c
           if(clouds_3d(i,j,k_start) .gt. 0.65) then
              if(hum(i,j) .gt. 70.) vismod(i,j) = 0.90
              if(hum(i,j) .gt. 80.) vismod(i,j) = 0.75
              if(hum(i,j) .gt. 90.) vismod(i,j) = 0.55
              if(hum(i,j) .gt. 95.) vismod(i,j) = 0.35
              go to 200
           endif
c
c.....     if no fog in lowest layer, find the maximum value in the 3 levels
c.....     above the surface.  then adjust the vismod based on humidity.
c
           amax_layer = 0.
           do k=k_start,k_end
              if(clouds_3d(i,j,k) .gt. amax_layer) then
                 amax_layer = clouds_3d(i,j,k)
              endif
           enddo !k
c
           if(amax_layer .gt. 0.65) then
              if(hum(i,j) .gt. 70.) vismod(i,j) = 0.95
              if(hum(i,j) .gt. 80.) vismod(i,j) = 0.80
              if(hum(i,j) .gt. 90.) vismod(i,j) = 0.60
              if(hum(i,j) .gt. 95.) vismod(i,j) = 0.40
           endif
c
200     continue
        enddo !i
        enddo !j
c
c..... that is all.
c
        return ! normal return
c
 999    print *,' error reading comment field from lc3.'
        print *,' no visibility modification done.'
        return
        end
