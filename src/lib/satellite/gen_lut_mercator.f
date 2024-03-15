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
       subroutine gen_lut_mercator(isat,jtype,kchl,
     &nx_l,ny_l,lat,lon,ri_laps,rj_laps,jstatus)
c
c
c     returns satellite pixel locations for each point on model (laps) domain
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
      real    rj_laps(nx_l,ny_l)         ! output

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

      real    nxmc,nymc,nzmc
      real    dx,dy
      real    rlatc,rlonc,sw(2),ne(2)
      common /mcgrid/rlonc,rlatc,nxmc,nymc,sw,ne,dx,dy

      jstatus = -1
      cdtype=c_sat_types(jtype,isat)

      ct=c_channel_types(kchl,jtype,isat)
      nc=index(ct,' ')-1
      if(nc.le.0)nc=3
      call lvd_file_specifier(ct,indx,istatus)
c
c retrieve the latest nav info from a file if it exists
c

      if(cdtype.eq.'hko')then  !hong kong gms satellite data

         if(indx.eq.1.or.indx.eq.3.or.
     .      indx.eq.4.or.indx.eq.5)then   

            print*,'hko navigation information - 3-20-03'

            nxmc = 1024
            nymc = 1024
            rlatc=21.4543
            rlonc=+114.0000
            sw(1)=13.1233
            sw(2)=+105.0
            ne(1)=29.7921
            ne(2)=+122.9894
            dx = 1.9492
            dy = 1.8069

         endif

      endif

      print*,'sat nav parameters - mercator'
      print*,'-----------------------------'
      print*,'nxmc:  ',nxmc
      print*,'nymc:  ',nymc
      print*,'rlatc: ',rlatc
      print*,'rlonc: ',rlonc
      print*,'sw: ',sw(1),sw(2)
      print*,'ne: ',ne(1),ne(2)
      print*,'dx:    ',dx
      print*,'dy:    ',dy
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
      if(cdtype.eq.'hko')then           !only other type is taiwan gms

         call  latlon_2_mcij(nx*ny,xlat,xlon,ri,rj)

      endif

      write(6,*)'sat ri/rj corners for domain'
      write(6,*)'ri1/rj1 (sw) ',ri(1,1),rj(1,1)
      write(6,*)'ri2/rj2 (se) ',ri(nx,1),rj(nx,1)
      write(6,*)'ri3/rj3 (nw) ',ri(1,ny),rj(1,ny)
      write(6,*)'ri4/rj4 (ne) ',ri(nx,ny),rj(nx,ny)
      write(6,*)
c
c get new i/j start/end values for this domain
c
      call get_sat_boundary(nx,ny,nx,ny,0,nxmc,nymc,
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
        do i = 1,nx_l,10
        do j = 1,ny_l,10
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
      cname=path(1:n1)//c_sat_id(isat)//'-llij-'//ct(1:nt)
      n1=index(cname,' ')-1
      table_path = cname(1:n1)//'-'//cdtype//'.lut'

      n1=index(table_path,' ')
      write(6,*)'write lat/lon to i/j look up table'
      write(6,*)table_path(1:n1)

c     call write_table (table_path,nx_l,ny_l,xlat,xlon,
c    &ri_laps,rj_laps,istatus)
c     if(istatus .ne. 1)then
c        write(6,*)'error writing look-up table'
c        goto 1000
c     endif

      if(elemstart.le.0)elemstart=1
      if(elemend.gt.nxmc)elemend=nxmc
      if(linestart.le.0)linestart=1
      if(lineend.gt.nymc)lineend=nymc

      if(indx.eq.1)then

         i_start_vis(jtype,isat)= elemstart
         i_end_vis(jtype,isat) = elemend
         j_start_vis(jtype,isat) = linestart
         j_end_vis(jtype,isat) = lineend
         r_la1(jtype,isat) = rlatc
         r_lo1(jtype,isat) = rlonc
         r_resolution_x_vis(jtype,isat) = dx  !*1000.
         r_resolution_y_vis(jtype,isat) = dy  !*1000.
         n_pixels_vis(jtype,isat) = nxmc
         n_lines_vis(jtype,isat)  = nymc

      elseif(indx.eq.2.or.indx.eq.4.or.indx.eq.5)then

         i_start_ir(jtype,isat) = elemstart
         i_end_ir(jtype,isat) = elemend
         j_start_ir(jtype,isat) = linestart
         j_end_ir(jtype,isat) = lineend
         r_la1(jtype,isat) = rlatc
         r_lo1(jtype,isat) = rlonc
         r_resolution_x_ir(jtype,isat) = dx  !*1000.
         r_resolution_y_ir(jtype,isat) = dy  !*1000.
         n_pixels_ir(jtype,isat) = nxmc
         n_lines_ir(jtype,isat)  = nymc

      elseif(indx.eq.3)then

         i_start_wv(jtype,isat) = elemstart
         i_end_wv(jtype,isat) = elemend
         j_start_wv(jtype,isat) = linestart
         j_end_wv(jtype,isat) = lineend
         r_la1(jtype,isat) = rlatc
         r_lo1(jtype,isat) = rlonc
         r_resolution_x_wv(jtype,isat) = dx  !*1000.
         r_resolution_y_wv(jtype,isat) = dy  !*1000.
         n_pixels_wv(jtype,isat) = nxmc
         n_lines_wv(jtype,isat)  = nymc

      endif

      jstatus = 1

1000  return
      end
