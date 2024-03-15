      subroutine rdgr2head(fname, nx, ny, xmin, ymin,
     +     offset_x, offset_y,
     +     dx, dy, sub_lon_degrees)

      integer nf_fid, nf_vid, nf_status
      integer nx, ny 
      real dx, dy, xmin, ymin, offset_x, offset_y
      character*(*) fname
      include 'netcdf.inc'
c
c  open netcdf file for reading
c
      call s_len(fname,nf)
      print*,'open and read file ',fname(1:nf)
      nf_status = nf_open(fname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',fname(1:nf)
      endif


c   variables of type real
c
c     variable        netcdf long name
c      dx           "x grid increment"
c
        nf_status = nf_inq_varid(nf_fid,'dx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,dx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
      endif
c
c     variable        netcdf long name
c      dy           "y grid increment"
c
      nf_status = nf_inq_varid(nf_fid,'dy',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,dy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
      endif
c
c     variable        netcdf long name
c      nominal_satellite_subpoint_lon 
c
      nf_status = nf_inq_varid(nf_fid,'nominal_satellite_subpoint_lon'
     +                        ,nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nominal_satellite_subpoint_lon'
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,sub_lon_degrees)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nominal_satellite_subpoint_lon'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      nx           "number of x point"
c
        nf_status = nf_inq_varid(nf_fid,'nx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nx'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,nx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nx'
      endif
c
c     variable        netcdf long name
c      ny           "number of y points"
c
        nf_status = nf_inq_varid(nf_fid,'ny',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ny'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,ny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ny'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      reftime      "reference time"
c
        nf_status = nf_inq_varid(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reftime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,reftime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reftime'
      endif
c
c     variable        netcdf long name
c      valtime      "valid time"
c
        nf_status = nf_inq_varid(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var valtime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,valtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var valtime'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      earth_shape  
c
        nf_status = nf_inq_varid(nf_fid,'earth_shape',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var earth_shape'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,earth_shape)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var earth_shape'
      endif
c
c     variable        netcdf long name
c      grid_name    
c
        nf_status = nf_inq_varid(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_name'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,grid_name)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_name'
      endif
c
c     variable        netcdf long name
c      grid_type    "grib-1 grid type"
c
        nf_status = nf_inq_varid(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_type'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,grid_type)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_type'
      endif
c
c     variable        netcdf long name
c      origin_name  
c
        nf_status = nf_inq_varid(nf_fid,'origin_name',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var origin_name'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,origin_name)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var origin_name'
      endif
c
c     variable        netcdf long name
c      process_name 
c
        nf_status = nf_inq_varid(nf_fid,'process_name',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var process_name'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,process_name)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var process_name'
      endif

      return
      end
