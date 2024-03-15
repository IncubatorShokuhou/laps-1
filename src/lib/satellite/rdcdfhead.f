      subroutine rdcdfhead(fname, nx, ny, center_id, 
     +     process_id, wmo_sat_id, dx, dy, la1, latin, lo1, 
     +     lov, reftime, valtime, earth_shape, grid_name, grid_type, 
     +     origin_name, process_name, wavelength, x_dim, y_dim)
c
      integer nf_fid, nf_vid, nf_status
      integer nx, ny, center_id, process_id,
     +     wmo_sat_id
      real dx, dy, la1, latin, lo1, lov
      double precision reftime, valtime
      character*132 origin_name
      character*132 x_dim
      character*132 y_dim
      character*132 earth_shape
      character*132 wavelength
      character*132 grid_name
      character*132 process_name
      character*132 grid_type
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
c      latin        "conic tangent latitude"
c
        nf_status = nf_inq_varid(nf_fid,'latin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,latin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin'
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
c      lov          "orientation of grid"
c
        nf_status = nf_inq_varid(nf_fid,'lov',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lov'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,lov)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lov'
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
c
c     variable        netcdf long name
c      center_id    "center id"
c
        nf_status = nf_inq_varid(nf_fid,'center_id',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var center_id'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,center_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var center_id'
      endif
c
c     variable        netcdf long name
c      process_id   "process id"
c
        nf_status = nf_inq_varid(nf_fid,'process_id',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var process_id'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,process_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var process_id'
      endif
c
c variable not in some netcdf files. commented on 9-6-06: jrs.
c ------------------------------------------------------------
c     variable        netcdf long name
c      wmo_sat_id   "wmo satellite id number"
 
      nf_status = nf_inq_varid(nf_fid,'wmo_sat_id',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wmo_sat_id'
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,wmo_sat_id)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var wmo_sat_id'
        endif
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
c
c     variable        netcdf long name
c      wavelength   "wavelength of satellite data"
c
        nf_status = nf_inq_varid(nf_fid,'wavelength',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wavelength'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,wavelength)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wavelength'
      endif
c
c     variable        netcdf long name
c      x_dim        "longitude dimension"
c
        nf_status = nf_inq_varid(nf_fid,'x_dim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var x_dim'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,x_dim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var x_dim'
      endif
c
c     variable        netcdf long name
c      y_dim        "latitude dimension"
c
        nf_status = nf_inq_varid(nf_fid,'y_dim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var y_dim'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,y_dim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var y_dim'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
