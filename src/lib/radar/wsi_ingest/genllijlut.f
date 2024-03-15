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
      subroutine gen_rirj_ce(imax,jmax,lat,lon,ri,rj
     +    ,istatus)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
c routine reads netcdf nowrad (wsi) data using subroutine read_wsi_cdf.
c nowrad data is then remapped to laps domain given lat/lon of domain.
c routine automatically moves the boundary for any domain with the nowrad
c confines.
c
       implicit none

       integer imax,jmax
       real lat(imax,jmax)
       real lon(imax,jmax)
       real ri(imax,jmax)
       real rj(imax,jmax)

       real dx,dy
       real rla1,rlo1
       real rla2,rlo2
       real rlat,rlon
       real rlatc,rlonc
       real dlat,dlon
       real nw(2),se(2)
       real r_missing_data

       integer i,j
       integer n,n1
       integer nx,ny,nz
       integer istatus
       integer nlines
       integer nelems
       integer i_out
       integer j_out

       character table_path*255
       character path*100
       character file*255

       common /cegrid/nx,ny,nz,nw,se,rlatc,rlonc
c
c ***************************************************************************
c
      istatus=1
c
c need to get latest filetime from appropriate data subdirectory
c
      write(6,*)'nav parameters'
      write(6,*)'dx     ',dx
      write(6,*)'dy     ',dy
      write(6,*)'nelems ',nelems
      write(6,*)'nlines ',nlines
      write(6,*)'rla1   ',rla1
      write(6,*)'rlo1   ',rlo1
      write(6,*)'rla2   ',rla2
      write(6,*)'rlo2   ',rlo2
      write(6,*)'rlon   ',rlon
      write(6,*)'rlat   ',rlat
c
      rlatc = rlat
      rlonc = rlon
      nx = nelems
      ny = nlines
      nw(1)=rla1
      nw(2)=rlo1
      se(1)=rla2
      se(2)=rlo2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  build ri/rj look up for laps domain
c
      call latlon_2_ceij(imax*jmax,lat,lon,ri,rj)

      write(6,*)'ri/rj corners for domain: '
      write(6,*)'--------------------------'
      write(6,*)'(sw) ',ri(1,1),rj(1,1)
      write(6,*)'(se) ',ri(imax,1),rj(imax,1)
      write(6,*)'(nw) ',ri(1,jmax),rj(1,jmax)
      write(6,*)'(ne) ',ri(imax,jmax),rj(imax,jmax)
      write(6,*)
      write(6,*)'check for domain coverage'
c
c check for laps domain outside data domain
c
       call check_domain_vrc(imax,jmax,ri,rj,nx,ny)

c
c output
c
      write(6,*)'ri/rj corners for domain after check: '
      write(*,*)'--------------------------------------'
      write(6,*)'(sw) ',ri(1,1),rj(1,1)
      write(6,*)'(se) ',ri(imax,1),rj(imax,1)
      write(6,*)'(nw) ',ri(1,jmax),rj(1,jmax)
      write(6,*)'(ne) ',ri(imax,jmax),rj(imax,jmax)
      write(6,*)

      if(.false.)then
        do i = 1,imax,10
        do j = 1,jmax,10
           write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)
        enddo
        enddo
      endif

      goto 16

901   write(6,*)'error reading parm file ',file(1:n)

16    write(6,*)'finished in get_llij_lut_wsi'
      return
      end
c
c *******************************************************************
c *******************************************************************
c
      subroutine gen_rirj_lam(imax,jmax,lat,lon,nelems,nlines,
     &        rla1,rlo1,rla2,rlo2,rlatin,rlov,dx,dy,
     &        ri,rj,istatus)
c
      implicit none

      integer imax,jmax

      real    lat(imax,jmax)
      real    lon(imax,jmax)
      real    rla1,rlo1
      real    rla2,rlo2
      real    dx,dy
      real    du,dv
      real    u_orig,v_orig
      real    rlatin,rlov
      real    pi
      real    u,v

      real    ri1,rj1
      real    ri2,rj2
      real    ri3,rj3
      real    ri4,rj4

      real    lapterm
      real    lovterm
      real    latterm
      real    lonterm
      real    dxterm
      real    dyterm

      real    ri(imax,jmax)
      real    rj(imax,jmax)

      integer   i,j,n,nn
      integer   n1,n2
      integer   nlines,nelems
      integer   istatus

      character*200 table_path
      character*200 file
      character*255 path
      character*200 cname
      character*3   ctype
c
c =================================
      pi = acos(-1.0)
c
      write(6,*)'nav parameters'
      write(6,*)'dx     ',dx
      write(6,*)'dy     ',dy
      write(6,*)'nelems ',nelems
      write(6,*)'nlines ',nlines
      write(6,*)'la1   ',rla1
      write(6,*)'lo1   ',rlo1
c     write(6,*)'la2   ',rla2
c     write(6,*)'lo2   ',rlo2
      write(6,*)'lov   ',rlov
      write(6,*)'latin ',rlatin

c compute lat lons  fsl conus grid (or read fm disk)
c determine ri/rj pair for the four laps domain corners.
c compute i,j  in lambert grid for desired points.
c get delta u and delta v in the lambert grid

      dxterm = dx/1000.
      dyterm = dy/1000.
      call getdudv_lam(rlov,rlatin,dxterm,dyterm,
     &       rla1,rlo1,du,dv,u_orig,v_orig)
      lapterm=(rlatin*pi)/180.
      lovterm=(rlov*pi)/180.

      do j = 1, jmax
      do i = 1, imax
         latterm=(lat(i,j)*pi)/180.
         lonterm=(lon(i,j)*pi)/180.
         call getuv_lam (lapterm,lovterm,
     &       latterm,lonterm,u,v)
         call uv_ij (nlines,u_orig,v_orig,du,dv,
     &       u,v,ri(i,j),rj(i,j))
      enddo
      enddo

      call check_domain_vrc(imax,jmax,ri,rj,nelems,nlines)
c
      if(.false.)then
         do i = 1,imax,10
         do j = 1,jmax,10
            write(6,*)'i,j,ri,rj: ',i,j,ri(i,j),rj(i,j)
         enddo
         enddo
      endif

c
      write(6,*)'lambert-conus (wfo) corners for domain: '
      write(6,*)'----------------------------------------'
      write(6,*)'(sw) ',ri(1,1),rj(1,1)
      write(6,*)'(se) ',ri(imax,1),rj(imax,1)
      write(6,*)'(nw) ',ri(1,jmax),rj(1,jmax)
      write(6,*)'(ne) ',ri(imax,jmax),rj(imax,jmax)
      write(6,*)

      goto 900

901   write(6,*)'error openning file ',file(1:n)

900   return
      end

      subroutine check_domain_vrc(imax,jmax,ri,rj,nx,ny)
c
c check for laps domain outside data domain
c
       implicit none

       integer imax,jmax
       integer i_out,j_out
       integer i,j
       integer nx,ny
       integer istatus

       real  ri(imax,jmax),rj(imax,jmax)
       real  r_missing_data

       call get_r_missing_data(r_missing_data,istatus)
       i_out=0
       j_out=0
       do j=1,jmax
         do i=1,imax
           if(ri(i,j).gt.nx)then
              i_out=i_out+1
              ri(i,j)=r_missing_data
           endif
           if(ri(i,j).le.0.)then
              i_out=i_out+1
              ri(i,j)=r_missing_data
           endif
           if(rj(i,j).gt.ny)then
              j_out=j_out+1
              rj(i,j)=r_missing_data
           endif
           if(rj(i,j).le.0.0)then
              j_out=j_out+1
              rj(i,j)=r_missing_data
           endif
         enddo
       enddo

       if(i_out.gt.0)then
          print*,'warning! laps domain outside of data domain'
          print*,'i out =',i_out
       endif
       if(j_out.gt.0)then
          print*,'warning! laps domain outside of data domain'
          print*,'j out =',j_out
       endif
       if(i_out.gt.0.or.j_out.gt.0)then
          print*,'check yor laps domain (nest7grid.parms) settings'
          print*,'if you do not want partial coverage of domain'
       endif

       return
       end

