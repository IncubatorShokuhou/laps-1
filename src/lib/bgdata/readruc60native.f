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
      subroutine read_ruc60_native(path,fname,af,nx,ny,nz,
     .                   pr,ht,tp,sh,uw,vw,gproj,lon0_ps,istatus)
c
c *** subroutine to read 60 km ruc data on the native polar stereographic,
c        hybrid-b grid.
c *** code modified from b. schwartz auto netcdf generator.
c
      implicit none
c

      include 'netcdf.inc'
      integer ncid, ntp, nvdim, lenstr, nvs, ndsize
c
      integer nx,ny,nz,rcode
      real   cp,g,cpog
      parameter (cp=1004.686,g=9.80665,cpog=cp/g)
c
c *** ruc arrays.
c
      real   pr(nx,ny,nz),       !output ruc pressure (mb)
     .       ht(nx,ny,nz),       !output ruc height (m)
     .       tp(nx,ny,nz),       !output ruc temperature (k)
     .       sh(nx,ny,nz),       !output ruc specific humidity (kg/kg)
     .       uw(nx,ny,nz),       !output ruc u-wind (m/s)
     .       vw(nx,ny,nz),       !output ruc v-wind (m/s)
     .       th(nx,ny,nz),       !ruc virtual potential temperature
     .       pc(nx,ny,nz)        !ruc condensation pressure
c
      real   psi(nx,ny),psj(nx,ny),
     .       lat(nx,ny),lon(nx,ny),
     .       mr,tv
c
      integer start(10),count(10)
      integer vdims(10)
      integer ndims,nvars,ngatts,recdim,nrecs
      character*31 dummy
c
      integer i,j,k,l,istatus
c
      character*(*) path
      character*9   fname
      character*4   af
      character*255 cdfname
      character*2   gproj
      
c
c *** common block variables for polar-stereographic grid.
c
      integer nx_ps,ny_ps,nz_ps  !no. of ps domain grid points
      real   lat0,lon0,rota,       !pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !sw lat, lon, ne lat, lon
      common /psgrid/nx_ps,ny_ps,nz_ps,lat0,lon0,rota,sw,ne
      real   lon0_ps             !returned for wind rotations
c_______________________________________________________________________________
c      
c *** open the netcdf file.
c
c      l=index(path,' ')-1
      call s_len(path,l)
      cdfname=path(1:l)//'/'//fname//af
      print *,'reading - ',cdfname(1:l+14)
      rcode=nf_open(cdfname,nf_nowrite,ncid)
      call ncinq(ncid,ndims,nvars,ngatts,recdim,rcode)
      call ncdinq(ncid,recdim,dummy,nrecs,rcode)
      if (nrecs .lt. 1) then
         print *,'not enough records in netcdf file.'
         istatus=0
         return
      endif
c
c *** read netcdf data.
c *** statements to fill uw.
c
      call ncvinq(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,1,start,count,uw)
c
c *** statements to fill vw.
c
      call ncvinq(ncid,2,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,2,start,count,vw)
c
c *** statements to fill ht (mont. stream fucntion / g).
c
      call ncvinq(ncid,3,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,3,start,count,ht)
c
c *** statements to fill pc.
c
      call ncvinq(ncid,4,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,4,start,count,pc)
c
c *** statements to fill pr.
c
      call ncvinq(ncid,5,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,5,start,count,pr)
c
c *** statements to fill th.
c
      call ncvinq(ncid,6,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,6,start,count,th)
c
c *** close netcdf file.
c
      rcode= nf_close(ncid)
c
c *** convert pascals to mb.
c *** compute temp and sh from thetav and pc.
c *** compute height from msf.
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=pr(i,j,k)*0.01
         call thvpc2tq(th(i,j,k),pc(i,j,k),pr(i,j,k),
     .                 tp(i,j,k),sh(i,j,k))
         mr=sh(i,j,k)/(1.-sh(i,j,k))
         tv=tp(i,j,k)*(1.+0.61*mr)
         ht(i,j,k)=ht(i,j,k)-cpog*tv
      enddo
      enddo
      enddo
c
c *** fill the polar-stereographic common block variables.
c
      gproj='ps'
      nx_ps=nx
      ny_ps=ny
      nz_ps=nz
      rota=0.
      lat0=90.0
      lon0=-105.0
      lon0_ps=lon0
      sw(1)=22.83730698
      sw(2)=-120.4905014
      ne(1)=45.98867416
      ne(2)=-60.82944107
c
c *** convert ruc winds from grid north to true north.
c commented js 11-09-00.
c     do j=1,ny
c     do i=1,nx
c        psi(i,j)=float(i)
c        psj(i,j)=float(j)
c     enddo
c     enddo
c     call psij_2_latlon(nx*ny,psi,psj,lat,lon)
c
c     call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz)
c
      istatus=1
      return
      end
