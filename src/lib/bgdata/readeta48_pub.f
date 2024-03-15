      subroutine get_eta48_dims(filename,nx,ny,nz
     &,stdlat1,stdlat2,lon0,la1,lo1,la2,lo2,istatus)
c
c updated to get all nav info - js 4-01
c
      implicit none
      include 'netcdf.inc'

      integer nz, nrecs, nx, ny, nf_fid, nf_vid
      integer nf_status,istatus
      real    stdlat1,stdlat2
      real    lon0
      real    la1,lo1
      real    la2,lo2

      character filename*200
c
c  open netcdf file for reading
c
      print*,'filename =',filename

      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', filename
      endif
c
c  fill all dimension values
c
c
c  get size of isolevel (returned as nz)
c
      nf_status = nf_inq_dimid(nf_fid,'isolevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim isolevel'
      endif

      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nz)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim isolevel'
      endif

c
c get size of record
c
      nf_status = nf_inq_dimid(nf_fid,'record',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim record'
      endif

      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nrecs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim record'
      endif

      if(nrecs.ne.1) then
         print *, 'error in eta48 format'
         stop
      endif
c
c get size of x (returned as nx)
c
      nf_status = nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
      endif

      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
      endif

c
c get size of y (returned as ny)
c
      nf_status = nf_inq_dimid(nf_fid,'y',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
      endif

      nf_status = nf_inq_dimlen(nf_fid,nf_vid,ny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
      endif
c
c     variable        netcdf long name
c      intlat1      "1st latitude of intersect earth-cone"
c                    returned as stdlat1
c
        nf_status = nf_inq_varid(nf_fid,'intlat1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var intlat1'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stdlat1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var intlat1'
      endif
c
c     variable        netcdf long name
c      intlat2      "2nd latitude of intersect earth-cone"
c                    returned as stdlat2
c
        nf_status = nf_inq_varid(nf_fid,'intlat2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var intlat2'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stdlat2)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var intlat2'
      endif
c
c     variable        netcdf long name
c      la1          "first latitude"
c
        nf_status = nf_inq_varid(nf_fid,'la1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la1'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,la1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la1'
      endif
c
c     variable        netcdf long name
c      la2          "last latitude"
c
        nf_status = nf_inq_varid(nf_fid,'la2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la2'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,la2)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la2'
      endif
c
c     variable        netcdf long name
c      lo1          "first longitude"
c
        nf_status = nf_inq_varid(nf_fid,'lo1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo1'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lo1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo1'
      endif
c
c     variable        netcdf long name
c      lo2          "last longitude"
c
        nf_status = nf_inq_varid(nf_fid,'lo2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo2'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lo2)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo2'
      endif
c
c     variable        netcdf long name
c      lon0         "true longitude"
c
        nf_status = nf_inq_varid(nf_fid,'lon0',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lon0'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lon0)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lon0'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      if(nf_status.ne.0)then
         print*,'set istatus to error condition'
         istatus = 0
         return
      endif

      istatus = 1
      return
      end
c
c
c  subroutine to read the file "eta 48 km awips regional conus " 
c
      subroutine read_eta_conusc(fname, nx,ny,nz, ht,p,t,uw,vw,
     +     rh, pvv, ht_sfc, p_sfc, rh_sfc, td_sfc, t_sfc, uw_sfc, vw_sfc
     +     ,mslp,istatus)

c kml: changes made april 2004
c td_sfc (model 2m dew point) is now being read in
c td_sfc will ultimately be used in subroutine sfcbkgd
c kml: end

      include 'netcdf.inc'
      integer nx,ny,nz, nf_fid, nf_status, k
      character*(*) fname
      real mslp( nx,  ny),
     +    ht( nx, ny, nz), ht_sfc( nx,  ny),
     +     p( nx, ny, nz),  p_sfc( nx,  ny), 
     +    rh( nx,  ny,  nz), rh_sfc( nx,  ny), 
     +                       td_sfc( nx,  ny),
     +     t( nx,  ny,  nz),  t_sfc( nx,  ny), 
     +    uw( nx,  ny,  nz), uw_sfc( nx,  ny), 
     +    vw( nx,  ny,  nz), vw_sfc( nx,  ny),
     +    pvv(nx,  ny,  nz),
     +    tmp(nz)
      real qsfc
      real make_ssh, make_td
      integer nxny,nxnynz
      character c8_project*8
      logical reverse_fields
      data reverse_fields/.false./

      istatus = 1

      nf_status = nf_open(fname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',fname
        return
      endif

      nxny = nx*ny
      nxnynz=nx*ny*nz

      call get_c8_project(c8_project,istatus)
      if(istatus.ne. 1)then
         print*,'error: returned from get_c8_project'
         print*,'error: current routine: read_eta_conusc'
         return
      endif
c
c     variable        netcdf long name
c      p            "isobaric levels" 
c
      call read_netcdf_real(nf_fid,'isolevel',nz,tmp,0,0,nf_status)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim isolevel'
        return
      endif

c
c if p is going from top to bottom, reverse it and set a flag to 
c also reverse the other fields
c 
      if(tmp(1).lt.tmp(nz)) then

         do k=1,nz
         do j=1,ny
         do i=1,nx
            p(i,j,nz-k+1) = tmp(k)
         enddo
         enddo
         enddo
         reverse_fields=.true.
      else
         do k=1,nz
         do j=1,ny
         do i=1,nx
            p(i,j,k) = tmp(k)
         enddo
         enddo 
         enddo
      endif


c
c     variable        netcdf long name
c      ht           "geopotential height" 
c
      call read_netcdf_real(nf_fid,'gh',nxnynz,ht,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the height field is missing'
         print*, 'aborting file processing for eta file ',fname
         return
      endif
         
      if(reverse_fields) call swap_array(nxny,nz,ht)
c
c     variable        netcdf long name
c      rh           "relative humidity" 
c
      call read_netcdf_real(nf_fid,'rh',nxnynz,rh,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the rh field is missing'
         print*, 'aborting file processing for eta file ',fname
         return
      endif

      if(reverse_fields) call swap_array(nxny,nz,rh)

c
c     variable        netcdf long name
c      t            "temperature" 
c
      call read_netcdf_real(nf_fid,'t',nxnynz,t,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the temp field is missing'
         print*, 'aborting file processing for eta file ',fname
         return
      endif


      if(reverse_fields) call swap_array(nxny,nz,t)
c
c     variable        netcdf long name
c      uw           "u-component of wind" 
c
      call read_netcdf_real(nf_fid,'uw',nxnynz,uw,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the u-wind field is missing'
         print*, 'aborting file processing for eta file ',fname
         return
      endif


      if(reverse_fields) call swap_array(nxny,nz,uw)
c
c     variable        netcdf long name
c      vw           "v-component of wind" 
c
      call read_netcdf_real(nf_fid,'vw',nxnynz,vw,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the v-wind field is missing'
         print*, 'aborting file processing for eta file ',fname
         return
      endif


      if(reverse_fields) call swap_array(nxny,nz,vw)
c
c     variable        netcdf long name
c      pvv           "pressure vertical velocity"
c
      call read_netcdf_real(nf_fid,'pvv',nxnynz,pvv,0,0,nf_status)
      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the w-wind field is missing'
         print*, 'aborting file processing for eta file ',fname
         return
      endif

      if(reverse_fields) call swap_array(nxny,nz,pvv)

c
c     variable        netcdf long name
c      mslp         "eta mean sea level pressure" 
c

      call read_netcdf_real(nf_fid,'emspmsl',nxny,mslp,0,0,nf_status)

c
c     variable        netcdf long name
c      ht_sfc       "surface geopotential height" 
c
      call read_netcdf_real(nf_fid,'gh_sfc',nxny,ht_sfc,0,0,nf_status)

c
c     variable        netcdf long name
c      p_sfc        "surface pressure" 
c
      call read_netcdf_real(nf_fid,'p_sfc',nxny,p_sfc,0,0,nf_status)

c fsl netcdf file only has sfc rh
c ---------------------------------------
      if(c8_project.eq.'nimbus')then
c
c     variable        netcdf long name
c      rh_sfc       "relative humidity 2m fixed height abv ground" 
c
      print*,'fsl-nimbus: read sfc rh'
      call read_netcdf_real(nf_fid,'rh_2mfh',nxny,rh_sfc,0,0,nf_status)
c
c     variable        netcdf long name
c      td_sfc        "dew point temp 2m fixed height abv ground"
c
      else
      print*,'read td_2mfh: kml upgrade'
      call read_netcdf_real(nf_fid,'td_2mfh',nxny,td_sfc,0,0,nf_status)
      endif
c ---------------------------------------
c
c     variable        netcdf long name
c      t_sfc        "temperature 2m fixed height abv ground" 
c
      call read_netcdf_real(nf_fid,'t_2mfh',nxny,t_sfc,0,0,nf_status)
c
c     variable        netcdf long name
c      uw_sfc       "u wind component 10m fixed height abv ground" 
c
      call read_netcdf_real(nf_fid,'uw_10mfh',nxny,uw_sfc,0,0,nf_status)
c
c     variable        netcdf long name
c      vw_sfc       "v wind component 10m fixed height abv ground" 
c
      call read_netcdf_real(nf_fid,'vw_10mfh',nxny,vw_sfc,0,0,nf_status)

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif
c
c routines in file lib/bgdata/sfcbkgd.f require td; thus, for fsl netcdf
c we must derive this from rh.
c
      if(c8_project.eq.'nimbus')then
       do j=1,ny
       do i=1,nx
         qsfc=make_ssh(p_sfc(i,j)/100.,t_sfc(i,j)-273.15,
     &rh_sfc(i,j)/100.,-132.)
         td_sfc(i,j)=make_td(p_sfc(i,j)/100.,t_sfc(i,j)-273.15,
     & qsfc,-132.)+273.15
       enddo
       enddo
      endif

      istatus = 0

      return
      end


      subroutine swap_array(n1,n2,a1)
      integer i,j,n1,n2
      real a1(n1,n2),a2(n1,n2)
      do j=1,n2
         do i=1,n1
            a2(i,j)=a1(i,n2-j+1)
         enddo
      enddo
      do j=1,n2
         do i=1,n1
            a1(i,j)=a2(i,j)
         enddo
      enddo
      return
      end
