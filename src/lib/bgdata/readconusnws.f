      subroutine get_conus_dims(fname,nx,ny,nz)
c
      implicit none
c
      include 'netcdf.inc'
c
      integer nx,ny,nz, j, ncid
      integer vdims(10),start(10),count(10)
      integer ntp,nvdim,nvs,lenstr,ndsize,rcode
c
      character*(*) fname
      character*31  dummy
c_______________________________________________________________________________
c
c *** open the netcdf file.
c
      rcode=nf_open(fname,nf_nowrite,ncid)
c
c *** statements to fill nx.
c
      call ncvinq(ncid,16,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_int(ncid,16,start,count,nx)
c
c *** statements to fill ny.
c
      call ncvinq(ncid,17,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_int(ncid,17,start,count,ny)
c
c *** statements to fill nz.
c
      call ncvinq(ncid,18,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_int(ncid,18,start,count,nz)
c
      rcode=nf_close(ncid)
c

      return 
      end
c
c===============================================================================
c
      subroutine read_conus_nws(path,fname,af,nx,ny,nz
     .          ,pr,ht,tp,sh,uw,vw
     .          ,gproj,lon0_lc,lat1_lc,lat2_lc,istatus)
c
      implicit none
c
      include 'netcdf.inc'
c
      integer nx,ny,nz,len
      integer nvs, ndsize, ncid, ntp, nvdim, lenstr
c
      integer rcode
c
c *** output arrays.
c
      real   pr(nx,ny,nz)
     .      ,ht(nx,ny,nz)
     .      ,tp(nx,ny,nz)
     .      ,sh(nx,ny,nz)
     .      ,uw(nx,ny,nz)
     .      ,vw(nx,ny,nz)
     .      ,prn(nz)
c
      real   lci(nx,ny),lcj(nx,ny),
     .       lat(nx,ny),lon(nx,ny)
c
      integer start(10),count(10)
      integer vdims(10) 
      character*31 dummy
c
      integer i,j,k,istatus
c
      character*(*) path
      character*256 bgname
      character*9   fname
      character*4   af
      character*2   gproj
c
      real   msgflg
c
c *** common block variables for lambert-conformal grid.
c
      integer nx_lc,ny_lc,nz_lc  !no. of lc domain grid points
      real   lat1,lat2,lon0,     !lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)         !sw lat, lon, ne lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
      real   lon0_lc
      real   lat1_lc,lat2_lc
c_______________________________________________________________________________
c
      istatus = 1

      msgflg=1.e30
c
c *** open the netcdf file.
c
      call s_len(path,len)
      bgname=path(1:len)//'/'//fname//af
      rcode=nf_open(bgname,nf_nowrite,ncid)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** read netcdf data.
c *** statements to fill prn.
c
      call ncvinq(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            return
         endif
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,1,start,count,prn)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** statements to fill ht.
c
      call ncvinq(ncid,2,dummy,ntp,nvdim,vdims,nvs,rcode)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            return
         endif
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,2,start,count,ht)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** statements to fill tp.
c
      call ncvinq(ncid,3,dummy,ntp,nvdim,vdims,nvs,rcode)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            return
         endif
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,3,start,count,tp)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** statements to fill sh.
c
      call ncvinq(ncid,4,dummy,ntp,nvdim,vdims,nvs,rcode)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            return
         endif
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,4,start,count,sh)
      rcode=nf_get_vara_real(ncid,3,start,count,tp)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
      rcode=nf_get_vara_real(ncid,3,start,count,tp)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** statements to fill uw.
c
      call ncvinq(ncid,5,dummy,ntp,nvdim,vdims,nvs,rcode)
      rcode=nf_get_vara_real(ncid,3,start,count,tp)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            return
         endif
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,5,start,count,uw)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** statements to fill vw.
c
      call ncvinq(ncid,6,dummy,ntp,nvdim,vdims,nvs,rcode)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         if(rcode.ne.nf_noerr) then
            print *, nf_strerror(rcode)
            return
         endif
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=nf_get_vara_real(ncid,6,start,count,vw)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** close netcdf file.
c
      rcode=nf_close(ncid)
      if(rcode.ne.nf_noerr) then
        print *, nf_strerror(rcode)
        return
      endif
c
c *** fill missing value flag (in netcdf file, missing = -99999.)
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=prn(k)
         if (ht(i,j,k) .lt. -10000.) ht(i,j,k)=msgflg
         if (tp(i,j,k) .lt. -10000.) tp(i,j,k)=msgflg
         if (sh(i,j,k) .lt. -10000.) sh(i,j,k)=msgflg
         if (uw(i,j,k) .lt. -10000.) uw(i,j,k)=msgflg
         if (vw(i,j,k) .lt. -10000.) vw(i,j,k)=msgflg
      enddo
      enddo
      enddo
c         
c *** fill lambert-conformal common block variables.
c
      gproj='lc'
      nx_lc=nx
      ny_lc=ny
      lat1=25.0
      lat1_lc=lat1
      lat2=25.0
      lat2_lc=lat2
      lon0=-95.0
      lon0_lc=lon0
      sw(1)=12.19
      sw(2)=-133.459
      ne(1)=57.29
      ne(2)=-49.3849
c
c *** convert ruc winds from grid north to true north.
c
cc      do j=1,ny
cc      do i=1,nx
cc         lci(i,j)=float(i)
cc         lcj(i,j)=float(j)
cc      enddo
cc      enddo
cc      call lcij_2_latlon(nx*ny,lci,lcj,lat,lon)
c
cc      call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz)
c
      istatus=0

      return
      end
