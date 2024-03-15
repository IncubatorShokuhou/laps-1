      subroutine rd_static_attr_sub(staticfile, nx, ny
     .,la1, latin1, latin2, lo1, lov, dx, dy
     .,c8_maproj,istatus)
c
      include 'netcdf.inc'
      integer nf_fid, nf_vid, nf_status
      integer nx, ny
      integer istatus
      real dx, dy, la1, latin1, latin2, lo1, lov

      character*132 grid_type
      character*132 grid_name
      character*132 earth_shape
      character*8   c8_maproj
      character*(*) staticfile
c
c  open netcdf file for reading
c
      istatus = 0
 
      nf_status = nf_open(trim(staticfile),nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',trim(staticfile)
        return
      endif
c
c   variables of type real
c
c     variable        netcdf long name
c      dx           "x grid increment"
c
      nf_status = nf_inq_varid(nf_fid,'dx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,dx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dx'
        return
      endif
c
c     variable        netcdf long name
c      dy           "y grid increment"
c
      nf_status = nf_inq_varid(nf_fid,'dy',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,dy)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var dy'
        return
      endif
c
c     variable        netcdf long name
c      la1          "first latitude"
c
      nf_status = nf_inq_varid(nf_fid,'la1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la1'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,la1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la1'
        return
      endif
c
c     variable        netcdf long name
c      latin1       "orientation of grid"
c
      nf_status = nf_inq_varid(nf_fid,'latin1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin1'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,latin1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin1'
        return
      endif
c
c     variable        netcdf long name
c      latin2       "orientation of grid"
c
      nf_status = nf_inq_varid(nf_fid,'latin2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin2'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,latin2)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin2'
        return
      endif
c
c     variable        netcdf long name
c      lo1          "first longitude"
c
      nf_status = nf_inq_varid(nf_fid,'lo1',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo1'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,lo1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo1'
        return
      endif
c
c     variable        netcdf long name
c      lov          "orientation of grid"
c
      nf_status = nf_inq_varid(nf_fid,'lov',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lov'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,lov)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lov'
        return
      endif
c
c     variable        netcdf long name
c      grid_spacing 
c
      nf_status = nf_inq_varid(nf_fid,'grid_spacing',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_spacing'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,grid_spacing)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_spacing'
        return
      endif
c
c   variables of type int
c
c
c     variable        netcdf long name
c      nx           "number of x points"
c
      nf_status = nf_inq_varid(nf_fid,'nx',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nx'
        return
      endif
      nf_status = nf_get_var_int(nf_fid,nf_vid,nx)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nx'
        return
      endif
c
c     variable        netcdf long name
c      ny           "number of y points"
c
      nf_status = nf_inq_varid(nf_fid,'ny',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ny'
        return
      endif
      nf_status = nf_get_var_int(nf_fid,nf_vid,ny)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ny'
        return
      endif
c
c      earth_shape  
c
      nf_status = nf_inq_varid(nf_fid,'earth_shape',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var earth_shape'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,earth_shape)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var earth_shape'
        return 
      endif
c
c     variable        netcdf long name
c      grid_name    
c
      nf_status = nf_inq_varid(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_name'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,grid_name)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_name'
        return
      endif
c
c     variable        netcdf long name
c      grid_type    "grib-1 grid type"
c
      nf_status = nf_inq_varid(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_type'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,grid_type)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_type'
        return
      endif

      call downcase(grid_type,grid_type)
      call s_len(grid_type,lgt)
      if(grid_type(1:lgt).eq.'tangential'.or.
     &grid_type(1:lgt).eq.'secant')then
         c8_maproj='lambert'
      elseif(grid_type(1:lgt).eq.'mercator')then
         c8_maproj='mercator'
      else
         c8_maproj='polar'
      endif
      call s_len(grid_type,lgt)

      if(trim(grid_type).eq.'tangential lambert conformal'
     &.or.trim(grid_type).eq.'secant lambert conformal')then
         c8_maproj='lambert'
      elseif(trim(grid_type).eq.'mercator')then
         c8_maproj='mercator'
      elseif(trim(grid_type).eq.'polar')then
         c8_maproj='polar'
      endif
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 1
      return
      end
