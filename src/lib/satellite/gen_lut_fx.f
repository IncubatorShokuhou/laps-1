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
       subroutine gen_lut_fx(isat,jtype,kchl,
     &nx_l,ny_l,lat,lon,ri_laps,rj_laps,jstatus)
c
c
c     returns satellite pixel locations for each point on model (laps) domain
c     initially implemented for himawari data with 1-d lat/lon metadata arrays
c
      implicit none

      integer   nx_l,ny_l

      real    lat(nx_l,ny_l)
      real    lon(nx_l,ny_l)   ! laps lat/lon data  -  input

      real    xlat(nx_l+2,ny_l+2)
      real    xlon(nx_l+2,ny_l+2) !expanded domain lats/lons
 
      real    pi

      real    r_missing_data
      real    rls,rle,res,ree

      real    ri(nx_l+2,ny_l+2)
      real    rj(nx_l+2,ny_l+2)
      real    rel_ri(nx_l+2,ny_l+2)
      real    rel_rj(nx_l+2,ny_l+2)
      real    ri_laps(nx_l,ny_l)         ! output
      real    rj_laps(nx_l,ny_l) ! output

      real emission_angle_d(nx_l+2,ny_l+2)
      real phase_angle_d(nx_l+2,ny_l+2)
      real specular_ref_angle_d(nx_l+2,ny_l+2)
      real azimuth_d(nx_l+2,ny_l+2)
      real range_m
      real sublat_d_a(nx_l+2,ny_l+2), sublon_d_a(nx_l+2,ny_l+2)

      integer isat,jtype,kchl
c     integer i,j,n,nc
      integer i,j,nc,nt
      integer ii,jj
      integer indx
      integer n1,np
      integer i1,j1
      integer nx,ny
      integer istatus
      integer jstatus
      integer linestart,lineend
      integer elemstart,elemend
      integer nijout
      integer isector /2/       ! 0=dynamic, 1=fulldisk, 2=conus
      integer i4time_latest, ifstat, i4time_sys
      integer i4_elapsed, ishow_timer

      logical lpoint

      character*200 table_path
      character*255 path
      character*255 fname
      character*200 cname
      character*9   cftime9
c     character*6   csatid  !satellite data identifier {gmssat only for now}
      character*3   cdtype !satellite data type       {'hko' only one for mercator currently}
      character*3   ct     !                         {'vis', 'ir', or 'wv'}
c
c ===========================================================
c the following include file contains static navigation parameters for
c all possible satellites.  it contains fortran logic that uses csatid,
c and csattyp to set the appropriate navigation parameters.
c
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'

      logical     lwrite
      data lwrite /.true./

      integer nxfx,nyfx
      real*8  dx,dy,xmin,ymin
      real*8  scale_x,scale_y,offset_x,offset_y
      real*8  sub_lon_degrees
      real    rlatc,rlonc,swlat,swlon,ne(2),lat0,lon0,geos_dist_m
      character*1 cgrddef /'n'/ ! latitude ordered n to s 

      common /fxgrid/offset_x,offset_y,dx,dy,xmin,ymin,sub_lon_degrees

      jstatus = -1
      cdtype=c_sat_types(jtype,isat)

      ct=c_channel_types(kchl,jtype,isat)
      nc=index(ct,' ')-1
      if(nc.le.0)nc=3
      call lvd_file_specifier(ct,indx,istatus)
c
c retrieve the latest nav info from a file if it exists
c
      print*,'setting fixed grid navigation information ',cdtype,indx

      if(cdtype .ne. 'gr2')then ! goes16 fixed grid data, sector pab
         if(indx.eq.1)then ! vis
            nxfx = 2048
            nyfx = 2048
            rlatc=(55. + (-5.)) / 2.
            rlonc=(68. + 148. ) / 2.
            xmin = 2048.
            ymin = 0.
            offset_x = -101353. * 1d-6 ! radians
            offset_y = +128323. * 1d-6 ! radians
            dx = 14. * 1d-6        ! radians
            dy = 14. * 1d-6        ! radians
            sub_lon_degrees = -75.0
         else              ! ir
            nxfx = 512
            nyfx = 512
            rlatc=(55. + (-5.)) / 2.
            rlonc=(68. + 148. ) / 2.
            xmin = 512.
            ymin = 0.
            offset_x = -101353. * 1d-6 ! radians
            offset_y = +128323. * 1d-6 ! radians
            dx = 56. * 1d-6        ! radians
            dy = 56. * 1d-6        ! radians
            sub_lon_degrees = -75.0
         endif

      else ! level 2 'gr2' type
         if(isector .eq. 0)then ! dynamic sector

!           call s_len(path_to_raw_sat(kchl,jtype,isat),np)
!           fname=path_to_raw_sat(kchl,jtype,isat)(1:np)
            fname=fname(1:np)//'*_'//ct
            call get_latest_file_time(fname,i4time_latest)
            call make_fnam_lp(i4time_latest,cftime9,ifstat)
!           fname=path_to_raw_sat(kchl,jtype,isat)(1:np)
            fname=fname(1:np)//cftime9//'_'//ct
         
            call rdgr2head(fname, nx, ny, xmin, ymin,
     +         offset_x, offset_y,
     +         dx, dy, sub_lon_degrees)

            continue ! read netcdf header (under construction)

         elseif(isector .eq. 1)then ! full disk
            if(indx.eq.1)then ! vis
               nxfx = 21696
               nyfx = 21696
               xmin = 0. ! float(nxfx)
               ymin = 0.
               offset_x = -0.151844    ! radians
               offset_y = +0.151844    ! radians
               dx = 14. * 1d-6         ! radians
               dy = 14. * 1d-6         ! radians
               sub_lon_degrees = -75.0 ! depends on satellite
            else              ! ir
               nxfx = 5424
               nyfx = 5424
               xmin = 0. ! float(nxfx)
               ymin = 0.
               offset_x = -0.151844    ! radians
               offset_y = +0.151844    ! radians
               dx = 56. * 1d-6         ! radians
               dy = 56. * 1d-6         ! radians
               sub_lon_degrees = -75.0 ! depends on satellite
            endif

         else ! isector = 2 (conus)
            if(indx.eq.1)then ! vis
               nxfx = 6000
               nyfx = 10000
               xmin = 0. ! float(nxfx)
               ymin = 0.
               offset_x = -0.101353    ! radians
               offset_y = +0.128233    ! radians
               dx = 14. * 1d-6         ! radians
               dy = 14. * 1d-6         ! radians
               sub_lon_degrees = -75.0 ! depends on satellite
            else              ! ir
               nxfx = 1500
               nyfx = 2500
               xmin = 0. ! float(nxfx)
               ymin = 0.
               offset_x = -0.101332    ! radians
               offset_y = +0.128212    ! radians
               dx = 56. * 1d-6         ! radians
               dy = 56. * 1d-6         ! radians
               sub_lon_degrees = -75.0 ! depends on satellite
            endif
         endif

      endif

      print*,' '
      print*,'gen_lut_fx: sat nav latlon parameters'
      print*,'-------------------------------------'
      print*,'nxfx:  ',nxfx
      print*,'nyfx:  ',nyfx
      print*,'rlatc: ',rlatc
      print*,'rlonc: ',rlonc
      print*,'sw: ',swlat,swlon
      print*,'ne: ',ne(1),ne(2)
      print*,'dx (radians):',dx
      print*,'dy (radians):',dy
      print*,'offset_x (radians):',offset_x
      print*,'offset_y (radians):',offset_y
      print*,'xmin:  ',xmin
      print*,'ymin:  ',ymin
      print*,'sub_lon_degrees:',sub_lon_degrees
c
c expand domain lats/lons: extra row and column "laps domain" used
c to build the domain relative look up table
c
      nx=nx_l+2
      ny=ny_l+2
      pi = acos(-1.0)
 
      call expand_domain(nx_l,ny_l,lat,lon,nx,ny,xlat,xlon,
     &istatus)
c
c laps domain as specified in laps lat/lon arrays.
c
      call latlon_2_fxij(nx*ny,xlat,xlon,ri,rj)

      i4_elapsed = ishow_timer()
      
      write(6,*)' filter out points beyond the limb'
      call get_systime_i4(i4time_sys,istatus)

      range_m = 42155680.00
      sublat_d_a(:,:) = 0.              ! only scalar might really be needed
      sublon_d_a(:,:) = sub_lon_degrees ! only scalar might really be needed
      call satgeom(i4time_sys,xlat,xlon,nx,ny
     1   ,sublat_d_a,sublon_d_a,range_m,r_missing_data
     1   ,phase_angle_d,specular_ref_angle_d,emission_angle_d
     1   ,azimuth_d,istatus)

      where (emission_angle_d(:,:) .le. 5.0)
          ri(:,:) = r_missing_data
          rj(:,:) = r_missing_data
      endwhere

      i4_elapsed = ishow_timer()

      write(6,*)'sat ri/rj for expanded model domain corners'
      write(6,*)'ri1/rj1 (sw) ',ri(1,1),rj(1,1)
      write(6,*)'ri2/rj2 (se) ',ri(nx,1),rj(nx,1)
      write(6,*)'ri3/rj3 (nw) ',ri(1,ny),rj(1,ny)
      write(6,*)'ri4/rj4 (ne) ',ri(nx,ny),rj(nx,ny)
      write(6,*)
c
c get new i/j start/end values for this domain
c
      call get_sat_boundary(nx,ny,nx,ny,0,nxfx,nyfx,
     &ri,rj,linestart,lineend,elemstart,elemend,
     &rls,rle,res,ree,istatus)

      if(istatus.ne.1)then
         write(6,*)'warning: laps domain outside sat data cover!'
      endif
c
c compute ri, rj relative look up table for the block of data surrounding
c the laps domain.
c
      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         print*,'error getting r_missing_data'
         return
      endif

      do j = 1,ny
      do i = 1,nx
       if(ri(i,j).ne.r_missing_data.and.rj(i,j).ne.r_missing_data)then
          rel_ri(i,j) = ri(i,j) - res
          rel_rj(i,j) = rj(i,j) - rls 
          if(lpoint)then
             i1=i
             j1=j
             lpoint=.false.
          endif
       else
          rel_ri(i,j) = r_missing_data
          rel_rj(i,j) = r_missing_data
          nijout=nijout+1
       endif

      enddo
      enddo

      if(nijout.gt.0)then
         print*,'found ',nijout,' points outside domain'
      endif
c
c put the expanded domain ri/rj's into the original laps domain
c
      jj = 0
      do j = 2,ny-1
      jj = jj+1
      ii = 0
      do i = 2,nx-1
         ii = ii+1
         ri_laps(ii,jj) = rel_ri(i,j)
         rj_laps(ii,jj) = rel_rj(i,j)
      enddo
      enddo

      if(lwrite)then
        do i = 1,nx_l,100
        do j = 1,ny_l,100
          print*,'i,j,ri,rj: ',i,j,ri_laps(i,j),rj_laps(i,j)
        enddo
        enddo
      endif

      call get_directory('static',path,n1)
      path=path(1:n1)//'/lvd/'
      n1=index(path,' ')-1

      nt=3
      if(indx.eq.3.or.indx.eq.4.or.indx.eq.5)then
         ct='ir'
         nt=2
         if(indx.eq.3)ct='wv'
      endif
      cname=path(1:n1)//c_sat_id(isat)//'-fxij-'//ct(1:nt)
      n1=index(cname,' ')-1
      table_path = cname(1:n1)//'-'//cdtype//'.lut'

      n1=index(table_path,' ')
      write(6,*)'not writing lat/lon to i/j look up table'
      write(6,*)table_path(1:n1)

c     call write_table (table_path,nx_l,ny_l,xlat,xlon,
c    &ri_laps,rj_laps,istatus)
c     if(istatus .ne. 1)then
c        write(6,*)'error writing look-up table'
c        goto 1000
c     endif

      if(elemstart.le.0)elemstart=1
      if(elemend.gt.nxfx)elemend=nxfx
      if(linestart.le.0)linestart=1
      if(lineend.gt.nyfx)lineend=nyfx

      geos_dist_m = 40000e3

      if(indx.eq.1)then

         i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         r_la1(jtype,isat) = rlatc
         r_lo1(jtype,isat) = rlonc
         r_resolution_x_vis(jtype,isat) = dx * geos_dist_m ! m
         r_resolution_y_vis(jtype,isat) = dy * geos_dist_m ! m
         n_pixels_vis(jtype,isat) = nxfx
         n_lines_vis(jtype,isat)  = nyfx

      elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

         i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         r_la1(jtype,isat) = rlatc
         r_lo1(jtype,isat) = rlonc
         r_resolution_x_ir(jtype,isat) = dx * geos_dist_m ! m
         r_resolution_y_ir(jtype,isat) = dy * geos_dist_m ! m
         n_pixels_ir(jtype,isat) = nxfx
         n_lines_ir(jtype,isat)  = nyfx

      elseif(indx.eq.3)then

         i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         r_la1(jtype,isat) = rlatc
         r_lo1(jtype,isat) = rlonc
         r_resolution_x_wv(jtype,isat) = dx * geos_dist_m ! m
         r_resolution_y_wv(jtype,isat) = dy * geos_dist_m ! m
         n_pixels_wv(jtype,isat) = nxfx
         n_lines_wv(jtype,isat)  = nyfx

      endif

      jstatus = 1

1000  write(6,*)' returning from gen_lut_fx'

      return
      end
