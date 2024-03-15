      subroutine read_radar_v01_cdf(filename,nxv01,nyv01,nzv01,
     &cref_com,cvel_com,cnyq_com,refd,veld,nyqd,istatus)
c
c note: nxv01= lon
c       nyv01= lat
c       nzv01= level
c
      real*4      refd(nxv01,nyv01,nzv01)
      real*4      veld(nxv01,nyv01,nzv01)
      real*4      nyqd(nxv01,nyv01,nzv01)

      character       cvel_com(nzv01)*126
      character       cref_com(nzv01)*126
      character       cnyq_com(nzv01)*126
      character*(*)   filename

      include 'netcdf.inc'
c
c  open netcdf file for reading
c
      n=index(filename,' ')-1
      print*,'open netcdf file',filename(1:n)
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',filename(1:n)
      endif
c
c  fill all dimension values
c
c
c get size of lat
c
      nf_status = nf_inq_dimid(nf_fid,'lat',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lat'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,lat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lat'
      endif
c
c get size of level
c
      nf_status = nf_inq_dimid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,level)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim level'
      endif
c
c get size of lon
c
      nf_status = nf_inq_dimid(nf_fid,'lon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lon'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,lon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lon'
      endif
c
c get size of n_2d_grids
c
      nf_status = nf_inq_dimid(nf_fid,'n_2d_grids',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim n_2d_grids'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,n_2d_grids)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim n_2d_grids'
      endif
      if(lon.ne.nxv01)then
         print *,'lon ne nxv01'
         stop
      endif
      if(lat.ne.nyv01)then
         print *,'lat ne nyv01'
         stop
      endif
      if(level.ne.nzv01)then
         print *,'level ne nzv01'
         stop
      endif

      call main_sub(nf_fid , nyv01, nzv01, nxv01, n_2d_grids,
     +cref_com,cvel_com,cnyq_com,refd,veld,nyqd)

      end
c
c
      subroutine main_sub(nf_fid, lat, level, lon, n_2d_grids,
     +refd_comment,veld_comment,nyqd_comment,refd,veld,nyqd)
      include 'netcdf.inc'
      integer lat, level, lon, n_2d_grids, nf_fid, nf_vid, nf_status
      character*18 asctime
      character*12 laps_domain_file
      character*132 model
      character*126 nyqd_comment(level)
      character*132 origin
      character*126 refd_comment(level)
      character*126 veld_comment(level)
      integer fctimes_var, imax, jmax, kdim, kmax, level_var(level),
     +     lvl(n_2d_grids), num_variables, nyqd_fcinv(level),
     +     refd_fcinv(level), veld_fcinv(level), version

      real nyqd( lon,  lat, level), refd( lon,  lat, level), veld(
     +     lon,  lat, level)

      call read_netcdf(nf_fid , lat, level, lon, n_2d_grids, asctime,
     +     fctimes_var, imax, jmax, kdim, kmax, laps_domain_file,
     +     level_var, lvl, model, num_variables, nyqd, nyqd_comment,
     +     nyqd_fcinv, origin, refd, refd_comment, refd_fcinv, veld,
     +     veld_comment, veld_fcinv, version)

c
c the netcdf variables are filled - your code goes here
c
      return
      end
      subroutine read_netcdf(nf_fid , lat, level, lon, n_2d_grids,
     +     asctime, fctimes_var, imax, jmax, kdim, kmax,
     +     laps_domain_file, level_var, lvl, model, num_variables,
     +     nyqd, nyqd_comment, nyqd_fcinv, origin, refd,
     +     refd_comment, refd_fcinv, veld, veld_comment, veld_fcinv,
     +     version)

      include 'netcdf.inc'
      integer lat, level, lon, n_2d_grids, nf_fid, nf_vid, nf_status

      character*18 asctime
      character*12 laps_domain_file
      character*132 model
      character*126 nyqd_comment(level)
      character*132 origin
      character*126 refd_comment(level)
      character*126 veld_comment(level)
      integer fctimes_var, imax, jmax, kdim, kmax, level_var(level),
     +     lvl(n_2d_grids), num_variables, nyqd_fcinv(level),
     +     refd_fcinv(level), veld_fcinv(level), version

      real nyqd( lon,  lat, level), refd( lon,  lat, level), veld(
     +     lon,  lat, level)

c
c     variable        netcdf long name
c      asctime
        nf_status = nf_inq_varid(nf_fid,'asctime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var asctime'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,asctime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ asctime '
      endif
c
c     variable        netcdf long name
c      laps_domain_file 
        nf_status = nf_inq_varid(nf_fid,'laps_domain_file',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var laps_domain_file'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,laps_domain_file)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ laps_domain_file '
      endif
c
c     variable        netcdf long name
c      model         
        nf_status = nf_inq_varid(nf_fid,'model',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var model'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,model)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ model '
      endif
c
c     variable        netcdf long name
c      nyqd_comment
       nf_status = nf_inq_varid(nf_fid,'nyqd_comment',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyqd_comment'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,nyqd_comment)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ nyqd_comment '
      endif
c
c     variable        netcdf long name
c      origin    
        nf_status = nf_inq_varid(nf_fid,'origin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var origin'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,origin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ origin '
      endif
c
c     variable        netcdf long name
c      refd_comment 
        nf_status = nf_inq_varid(nf_fid,'refd_comment',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var refd_comment'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,refd_comment)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ refd_comment '
      endif
c
c     variable        netcdf long name
c      veld_comment    
        nf_status = nf_inq_varid(nf_fid,'veld_comment',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var veld_comment'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,veld_comment)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ veld_comment '
      endif
c
c     variable        netcdf long name
c      fctimes_var  "forecast times" 
c
        nf_status = nf_inq_varid(nf_fid,'fctimes',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var fctimes'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,fctimes_var)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ fctimes_var '
      endif
c
c     variable        netcdf long name
c      imax
        nf_status = nf_inq_varid(nf_fid,'imax',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var imax'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,imax)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ imax '
      endif
c
c     variable        netcdf long name
c      jmax
        nf_status = nf_inq_varid(nf_fid,'jmax',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var jmax'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,jmax)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ jmax '
      endif
c
c     variable        netcdf long name
c      kdim
        nf_status = nf_inq_varid(nf_fid,'kdim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var kdim'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,kdim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ kdim '
      endif
c
c     variable        netcdf long name
c      kmax
        nf_status = nf_inq_varid(nf_fid,'kmax',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var kmax'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,kmax)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ kmax '
      endif
c
c     variable        netcdf long name
c      level_var    "level" 
c
        nf_status = nf_inq_varid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var level'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,level_var)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ level_var '
      endif
c
c     variable        netcdf long name
c      lvl
        nf_status = nf_inq_varid(nf_fid,'lvl',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lvl'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,lvl)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ lvl '
      endif
c
c     variable        netcdf long name
c      num_variables
        nf_status = nf_inq_varid(nf_fid,'num_variables',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var num_variables'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,num_variables)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ num_variables '
      endif
c
c     variable        netcdf long name
c      nyqd_fcinv
        nf_status = nf_inq_varid(nf_fid,'nyqd_fcinv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyqd_fcinv'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,nyqd_fcinv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ nyqd_fcinv '
      endif
c
c     variable        netcdf long name
c      refd_fcinv
        nf_status = nf_inq_varid(nf_fid,'refd_fcinv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var refd_fcinv'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,refd_fcinv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ refd_fcinv '
      endif
c
c     variable        netcdf long name
c      veld_fcinv
        nf_status = nf_inq_varid(nf_fid,'veld_fcinv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var veld_fcinv'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,veld_fcinv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ veld_fcinv '
      endif
c
c     variable        netcdf long name
c      version
        nf_status = nf_inq_varid(nf_fid,'version',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var version'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,version)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ version '
      endif
c
c     variable        netcdf long name
c      nyqd         "3d radar nyquist velocity" 
c
      print*,'read nyqd'
      nf_status = nf_inq_varid(nf_fid,'nyqd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyqd'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,nyqd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ nyqd '
      endif
c
c     variable        netcdf long name
c      refd         "3d radar" 
c
      print*,'read refd'
      nf_status = nf_inq_varid(nf_fid,'refd',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var refd'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,refd)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ refd '
      endif
c
c     variable        netcdf long name
c      veld         "3d radar" 
c
      print*,'read veld'
      nf_status = nf_inq_varid(nf_fid,'veld',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var veld'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,veld)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ veld '
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
