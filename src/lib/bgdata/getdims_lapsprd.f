      subroutine getdims_lapsprd(fullname,x,y,z,istatus)
      implicit none
      include 'netcdf.inc'
      integer record, x, y, z
      integer nf_fid, nf_vid, nf_status, ln
      integer istatus
      character*(*) fullname
c
c  open netcdf file for reading
c
      istatus = 0

      call s_len(fullname,ln)
      nf_status = nf_open(fullname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',fullname(1:ln)
        return
      endif
c
c get size of x
c
      nf_status = nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,x)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
        return
      endif
c
c get size of y
c
      nf_status = nf_inq_dimid(nf_fid,'y',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,y)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
        return
      endif
c
c get size of z
c
      nf_status = nf_inq_dimid(nf_fid,'z',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,z)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z'
        return
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 1

      return
      end
c
c  subroutine to read the file "laps fua/fsf files
c  - forecast model attributes" 
c
      subroutine read_lapsprd_attr(fullname,
     +     dx, dy, la1, lo1, latin1, latin2, lov, 
     +     grid_type, la2,lo2, istatus)
c
      implicit none
      include 'netcdf.inc'
      integer nf_fid, nf_vid, nf_status
      integer nx, ny      !<-- not returned as these are x/y from dims call
      integer ln
      integer istatus
      real dx, dy, la1, la2, lo1, lo2, lov
      real latin1, latin2
      character*132 grid_type_internal
      character*30 grid_type
      character*(*) fullname
c
c  open netcdf file for reading
c
      istatus = 0
      call s_len(fullname,ln)
      nf_status = nf_open(fullname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',fullname(1:ln)
        return
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
c      la2          "last latitude"
c
      nf_status = nf_inq_varid(nf_fid,'la2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la1'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,la2)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la2'
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
c      lo2          "last longitude"
c
      nf_status = nf_inq_varid(nf_fid,'lo2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo2'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,lo2)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lo2'
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

c   variables of type char
c
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
      nf_status = nf_get_var_text(nf_fid,nf_vid,grid_type_internal)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_type'
        return 
      endif
      grid_type=grid_type_internal

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return 
      endif

      istatus = 1
      return
      end
