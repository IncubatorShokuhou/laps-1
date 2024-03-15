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

       subroutine nowradwsi_to_laps(ctype,
     &                    filename,
     &                    nlines,nelems,
     &                    imax,jmax,
     &                    lat,lon,
     &                    validtime,
     &                    rdbz,
     &                    maxradars,
     &                    radar_dist_min,
     &                    istatus)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c routine reads netcdf nowrad (wsi) data using subroutine read_wsi_cdf.
c nowrad data is then remapped to laps domain given lat/lon of domain.
c routine automatically moves the boundary for any domain with the nowrad
c confines.

      integer maxradars
      integer nsites_present
      integer nsites_absent
      real    present_site_loc_i(maxradars)
      real    present_site_loc_j(maxradars)
      real    absent_site_loc_i(maxradars)
      real    absent_site_loc_j(maxradars)
      real    present_site_lat(maxradars)
      real    present_site_lon(maxradars)
      real    absent_site_lat(maxradars)
      real    absent_site_lon(maxradars)
      real    radar_lat(maxradars)
      real    radar_lon(maxradars)
      real    radar_elev(maxradars)
      real    radar_ri(maxradars)
      real    radar_rj(maxradars)

      real lat(imax,jmax)
      real lon(imax,jmax)
      real height_grid(imax,jmax)
      real rdbz(imax,jmax)
      real radar_dist_min(imax,jmax)
      real ri(imax,jmax)
      real rj(imax,jmax)
      real distmin,dist
      real rdistmax
      real xgrddismx
      real ygrddismx
      real grid_spacing_m
      real dlat,dlon
      real lov, latin, la1,la2,lo1,lo2
      real centerlon,toplat,dx,dy
      real nw(2),se(2)
      real pi,rad2dg,rii,rjj
      real r_missing_data
      real heigth_mean_grid
      real height_mean
      real height_sum
      real ridis,rjdis

      integer image_to_dbz(0:15)
      integer validtime
      integer istatus
      integer jstatus
      integer status
      integer i_base_value
      integer increment

      character filename*200
      character ctype*3

      integer image(nelems,nlines)
      integer indxofthosein(maxradars)

      common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc

c ----------------------------------------------------------------
c read data and get navigation info to build lat/lon to i/j table
c ----------------------------------------------------------------
      istatus=1

      nradars_in = 0 ! initialize this

      call read_nowrad_cdf(ctype,filename,nlines,nelems,
     + dlat,dlon,la1,lo1,la2,lo2,centerlon,toplat,validtime,
     + dx,dy,lov, latin, image, istatus)
c
c this routine converts the raw bytes to integers scaled 0 - 64
c it also determines the wsi grid i/j locations for absent
c or presently reporting radars.

      call cvt_wsi_nowrad(ctype,nelems,nlines,image
     +,nsites_present,present_site_loc_i,present_site_loc_j
     +,nsites_absent,absent_site_loc_i,absent_site_loc_j
     +,maxradars,istatus)
      if(istatus .ne. 1)then
          write(6,*)'bad status returned from cvt_wsi_nowrad'
          return
      endif

c     call cvt_wsi_nowrad(ctype,nelems,nlines,image,istatus)

      pi=acos(-1.)
      rad2dg=180.0/pi

      if(ctype.eq.'wfo')then

         lo1=-lo1
         lo2=-lo2
         if(lov.lt.-180.0)lov=lov+360.
         if(lo1.lt.-180.0)lo1=lo1+360.
         if(lo2.lt.-180.0)lo2=lo2+360.

         if(lov.gt.180.0)lov=lov-360.
         if(lo1.gt.180.0)lo1=lo1-360.
         if(lo2.gt.180.0)lo2=lo2-360.

         call gen_rirj_lam(imax,jmax,lat,lon,nelems,nlines
     &        ,la1,lo1,la2,lo2,latin,lov,dx,dy,ri,rj
     &        ,istatus)
c
 
      elseif(ctype.eq.'wsi')then
c
c load common block for cyclindrical equidistant
c
           nx=nelems
           ny=nlines
           nz=1
           nw(1)=la1
           nw(2)=lo1
           se(1)=la2
           se(2)=lo2
           rlatc=(toplat - (dlat*(nlines-1)*0.5))*rad2dg
           rlonc=centerlon*rad2dg

c--------------------------------------------------
c compute the wsi grid i/j values for domain lat/lon
c--------------------------------------------------           
           call latlon_2_ceij(imax*jmax,lat,lon,ri,rj)

           call check_domain_vrc(imax,jmax,ri,rj,nx,ny)

c--------------------------------------------------
c compute radar site lat/lon locations both present
c and missing for the wsi mosaic.
c--------------------------------------------------
           call ceij_2_latlon(nsites_present,present_site_loc_i
     .,present_site_loc_j,present_site_lat,present_site_lon)

           call ceij_2_latlon(nsites_absent,absent_site_loc_i
     .,absent_site_loc_j,absent_site_lat,absent_site_lon)

      endif

      write(*,*)' nowrad data prepared for requested time '
      write(*,*)
      write(*,*)'valid time ',validtime

c ---------------------------------------------------
c    build radar count level to dbz look up
c ---------------------------------------------------
      increment=5
      i_base_value=7
      do i = 1,15
         image_to_dbz(i)=(i-1)*increment+i_base_value
      end do
      image_to_dbz(0)=0
c ---------------------------------------------------
c compute grid ratio 
c ---------------------------------------------------
      call get_grid_spacing(grid_spacing_m,istatus)
      r_grid_ratio=sqrt(dx*dx + dy*dy)/grid_spacing_m
c ---------------------------------------------------
c  remap nowrad-wsi data to laps grid.
c ---------------------------------------------------
      call process_nowrad_z(imax,jmax,
     &                  r_grid_ratio,
     &                  image_to_dbz,
     &                  image,
     &                  ri,
     &                  rj,
     &                  nlines,nelems, ! input array dimensions
     &                  rdbz,
     &                  istatus)

c----------------------------------------------------
c determine which radars in wsi grid are in this domain
c and compute/save the minimum distance to radar.
c----------------------------------------------------
      call get_r_missing_data(r_missing_data,istatus)
      call  read_static_grid(imax,jmax,'avg',height_grid,istatus)

      if(ctype.eq.'wsi')then

         height_sum=0
         do j=1,jmax
         do i=1,imax
            height_sum=height_sum+height_grid(i,j)
         enddo
         enddo
         height_mean_grid=height_sum/(imax*jmax)
         nradars_tot=0
         rdistmax=480000.
         xgrddismx=imax+rdistmax/grid_spacing_m
         ygrddismx=jmax+rdistmax/grid_spacing_m
         xgrddismn=-rdistmax/grid_spacing_m
         ygrddismn=-rdistmax/grid_spacing_m
         height_sum=0.0

         print*,'mx num grid points for search x: ',xgrddismx
         print*,'mx num grid points for search y: ',ygrddismx
         print*,'mn num grid points for search x: ',xgrddismn
         print*,'mn num grid points for search y: ',ygrddismn

         do i=1,nsites_present
c
c fudge factor on site lat/lon due to unknown systematic
c offset.
           present_site_lat(i)=present_site_lat(i)-0.07
           present_site_lon(i)=present_site_lon(i)+0.02

           call latlon_to_rlapsgrid(present_site_lat(i),
     &                              present_site_lon(i),
     &                              lat,lon,            !laps lat/lon arrays
     &                              imax,jmax,          !laps horiz domain
     &                              rii,rjj,            !output: real i,j, scalars
     &                              jstatus)

           ridis=abs(rii)
           rjdis=abs(rjj)

c          if(ridis.le.xgrddismx .and. rjdis.le.ygrddismx)then

           if( (rii.le.xgrddismx .and. rii.ge.xgrddismn) .and.
     &         (rjj.le.ygrddismx .and. rjj.ge.ygrddismn)  )then

              nradars_tot = nradars_tot + 1
              radar_lat(nradars_tot)=present_site_lat(i)
              radar_lon(nradars_tot)=present_site_lon(i)
              radar_ri(nradars_tot)=rii
              radar_rj(nradars_tot)=rjj
              radar_elev(nradars_tot)=r_missing_data

              if( (rii.gt.0.0  .and.  rii.lt.imax) .and.
     &            (rjj.gt.0.0  .and.  rjj.lt.jmax)  )then

                 call bilinear_interp_extrap(rii,rjj,imax,jmax
     &                    ,height_grid,result,istatus)

                 height_sum=height_sum+result
                 radar_elev(nradars_tot)=result
                 nradars_in=nradars_in+1
                 indxofthosein(nradars_in)=nradars_tot
              endif

           endif

        enddo

        if(nradars_in.gt.0)then
           print*
           print*,'found ',nradars_in,' in this domain'
           height_mean=height_sum/nradars_in
        else
           print*
           print*,'interesting! no radars in domain'
           print*,'set mean radar height to mean grid elev'
           print*,' = ',height_mean_grid
           height_mean = height_mean_grid
        endif
        print*
        print*,'found ',nradars_tot,' to consider in min dist array'
        print*,'mean height of radars = ',height_mean
        print*

        do i=1,nradars_tot
        do j=1,nradars_in
         ii=indxofthosein(j)
         if(ii.eq.i)then
          print*,'#/ri/rj/lat/lon/elev: ',ii,radar_ri(ii),radar_rj(ii)
     .,radar_lat(ii),radar_lon(ii),radar_elev(ii)
         endif
        enddo
        enddo

        where(radar_elev.eq.r_missing_data)radar_elev=height_mean_grid
        
        print*,'determine nearest radar distance array'
        print*

        do j=1,jmax
        do i=1,imax

           distmin=r_missing_data
           do k=1,nradars_tot
              if(radar_elev(k).lt.10000.)then
                call latlon_to_radar(lat(i,j),lon(i,j),height_grid(i,j)
     1,azimuth,slant_range,elev,radar_lat(k),radar_lon(k),radar_elev(k))     
                if(slant_range.lt.distmin.and.slant_range.lt.rdistmax)
     &             distmin=slant_range
              endif
           enddo
           radar_dist_min(i,j)=distmin
       enddo
       enddo 

      endif

      goto 16
19    print*,'error in nowrad_2_laps, terminating'
      goto 16
14    print*,'nowrad data not found for given time or'
      print*,'data could be bad. no vrc produced'

16    print*,'finished in nowrad_2_laps'
      return
      end
