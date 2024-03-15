      subroutine rd_wsi_3dradar_cdf(cfname,lines,elems,
     &radsperline,radsperelem,la1,lo1,la2,lo2,centerlon,
     &validtime,datalevel,data_levels,num_levels,level_prefix,
     &image,istatus)

      include 'netcdf.inc'
c     integer datalevel, elems, lines, nf_fid, nf_vid, nf_status
      integer datalevel, elems, lines
      integer image(elems,lines)
      integer istatus
      integer validtime

      character*(*) cfname
      character  level_prefix(datalevel)*50

      real dx,dy,la1,lo1,centerlon
      real radsperelem,radsperline
c    +   , toplat

      istatus = 1

      call rd_wsi3d_sub(cfname,datalevel,elems,lines,
     .   dx,dy,la1,lo1,la2,lo2,centerlon,
     .   radsperline,radsperelem,num_levels,data_levels,
     .   level_prefix,image,validtime,istatus)
      if(istatus.ne.0)then
         print*,'error reading wsi 3d radar file'
         return
      endif

      istatus = 0

      return
      end
c
      subroutine rd_wsi3d_sub(cfname,datalevel,elems,lines,
     .   dx,dy,la1,lo1,la2,lo2,centerlon,radsperline,
     .   radsperelem,num_levels,data_levels,level_prefix,image,
     .   validtime,istatus)


      include 'netcdf.inc'
      integer datalevel, elems, lines
      character*50 grid_name
      character*50 grid_type
      character*50 level_prefix(datalevel)
      character*20 product_units
      character*50 x_dim
      character*50 y_dim
      character*(*) cfname
      integer nx, ny, data_levels(datalevel),image(elems,lines), 
     +   num_levels, validtime
      real dx,dy,la1,la2,lo1,lo2,centerlon,difflon,
     +   radsperelem,  radsperline, toplat

      call read_cdf_wsi3d_head(cfname,datalevel, dx, dy,
     +    la1, la2, lo1, lo2, nx, ny, centerlon,
     +    data_levels,
     +    difflon, grid_name, grid_type, level_prefix,
     +    num_levels, product_units, radsperelem, radsperline,
     +    toplat, validtime, x_dim, y_dim, istatus)
c

      call read_cdf_wsi3d_image(cfname,elems,lines,image,istatus)

      return
      end
c
c  subroutine to read the file header "wsi nexrad data" 
c
      subroutine read_cdf_wsi3d_head(cfname,datalevel, dx, dy,
     +    la1, la2, lo1, lo2, nx, ny,  centerlon,
     +    data_levels,difflon,grid_name,grid_type,level_prefix,
     +    num_levels, product_units, radsperelem, radsperline,
     +    toplat, validtime, x_dim, y_dim,istatus)
      include 'netcdf.inc'
c     integer datalevel, elems, lines, nf_fid, nf_vid, nf_status
      integer datalevel, nf_fid, nf_vid, nf_status

      character*50 grid_name
      character*50 grid_type
      character*50 level_prefix(datalevel)
      character*20 product_units
      character*50 x_dim
      character*50 y_dim
      character*(*) cfname
      integer nx, ny, data_levels(datalevel), 
     +   num_levels, validtime
      real dx, dy, la1, la2, lo1, lo2, centerlon,
     +   difflon, radsperelem,  radsperline, toplat
c
c  open netcdf file for reading
c
      nf_status = nf_open(cfname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        call s_len(cfname,nf)
        print *, nf_strerror(nf_status)
        print *,'nf_open ',cfname(1:nf)
        return
      endif
c
c     variable        netcdf long name
c      grid_name    "grid name" 
c
      istatus=1
      nf_status = nf_inq_varid(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid_name'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,grid_name)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ grid_name '
        return
      endif
c
c     variable        netcdf long name
c      grid_type    "grid type" 
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
        print *,'in nf_get_var_ grid_type '
        return
      endif
c
c     variable        netcdf long name
c      level_prefix "level prefix" 
c
      nf_status = nf_inq_varid(nf_fid,'level_prefix',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var level_prefix'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,level_prefix)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ level_prefix '
        return
      endif
c
c     variable        netcdf long name
c      product_units"product units" 
c
      nf_status = nf_inq_varid(nf_fid,'product_units',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var product_units'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,product_units)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ product_units '
        return
      endif
c
c     variable        netcdf long name
c      x_dim        "longitude dimension" 
c
      nf_status = nf_inq_varid(nf_fid,'x_dim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var x_dim'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,x_dim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ x_dim '
        return
      endif
c
c     variable        netcdf long name
c      y_dim        "latitude dimension" 
c
      nf_status = nf_inq_varid(nf_fid,'y_dim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var y_dim'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,y_dim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ y_dim '
        return
      endif
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
        print *,'in nf_get_var_ nx '
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
        print *,'in nf_get_var_ ny '
        return
      endif
c
c     variable        netcdf long name
c      data_levels  "data levels" 
c
      nf_status = nf_inq_varid(nf_fid,'data_levels',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var data_levels'
        return
      endif
      nf_status = nf_get_var_int(nf_fid,nf_vid,data_levels)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ data_levels '
        return
      endif
c
c     variable        netcdf long name
c      num_levels   "number of levels" 
c
      nf_status = nf_inq_varid(nf_fid,'num_levels',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var num_levels'
        return
      endif
      nf_status = nf_get_var_int(nf_fid,nf_vid,num_levels)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ num_levels '
        return
      endif
c
c     variable        netcdf long name
c      validtime    "valid time" 
c
      nf_status = nf_inq_varid(nf_fid,'validtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var validtime'
        return
      endif
      nf_status = nf_get_var_int(nf_fid,nf_vid,validtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ validtime '
        return
      endif
c
c     variable        netcdf long name
c      dx           
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
        print *,'in nf_get_var_ dx '
        return
      endif
c
c     variable        netcdf long name
c      dy           
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
        print *,'in nf_get_var_ dy '
        return
      endif
c
c     variable        netcdf long name
c      la1          
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
        print *,'in nf_get_var_ la1 '
        return
      endif
c
c     variable        netcdf long name
c      la2          
c
      nf_status = nf_inq_varid(nf_fid,'la2',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var la2'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,la2)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ la2 '
        return
      endif
c
c     variable        netcdf long name
c      lo1          
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
        print *,'in nf_get_var_ lo1 '
        return
      endif
c
c     variable        netcdf long name
c      lo2          
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
        print *,'in nf_get_var_ lo2 '
        return
      endif
c
c     variable        netcdf long name
c      centerlon    "center longitude" 
c
      nf_status = nf_inq_varid(nf_fid,'centerlon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var centerlon'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,centerlon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ centerlon '
        return
      endif
c
c     variable        netcdf long name
c      difflon      "difference longitude" 
c
      nf_status = nf_inq_varid(nf_fid,'difflon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var difflon'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,difflon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ difflon '
        return
      endif
c
c     variable        netcdf long name
c      radsperelem  "radians per element" 
c
      nf_status = nf_inq_varid(nf_fid,'radsperelem',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radsperelem'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,radsperelem)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ radsperelem '
        return
      endif
c
c     variable        netcdf long name
c      radsperline  "radians per line" 
c
      nf_status = nf_inq_varid(nf_fid,'radsperline',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radsperline'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,radsperline)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ radsperline '
        return
      endif
c
c     variable        netcdf long name
c      toplat       "top latitude" 
c
      nf_status = nf_inq_varid(nf_fid,'toplat',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var toplat'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,toplat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ toplat '
        return
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0

      return
      end

c
c --------------------------------------------------------
c
      subroutine get_wsi_3drad_dims(cfname,datalevel,elems,lines,
     &                              istatus)

      include 'netcdf.inc'
      integer datalevel, elems, lines, nf_fid, nf_vid, nf_status
      integer istatus
      character*(*) cfname

      istatus=1  !return error status 
c
c  open netcdf file for reading
c
      nf_status = nf_open(cfname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        call s_len(cfname,nf)
        print *, nf_strerror(nf_status)
        print *,'nf_open',cfname(1:nf)
        return
      endif
c
c  fill all dimension values
c
c
c get size of datalevel
c
      nf_status = nf_inq_dimid(nf_fid,'datalevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim datalevel'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,datalevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim datalevel'
        return
      endif
c
c get size of elems
c
      nf_status = nf_inq_dimid(nf_fid,'elems',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim elems'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,elems)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim elems'
        return
      endif
c
c get size of lines
c
      nf_status = nf_inq_dimid(nf_fid,'lines',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lines'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,lines)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim lines'
        return
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0
      return
      end
c
c ----------------------------------------------------
c
      subroutine read_cdf_wsi3d_image(cfname,
     &              elems,lines,image,istatus)
c
c     variable        netcdf long name
c      image        "image pixel values"
c
      include 'netcdf.inc'
      integer elems,lines
      integer image(elems,lines)
      character*(*) cfname
c
c  open netcdf file for reading
c
      nf_status = nf_open(cfname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        call s_len(cfname,nf)
        print *, nf_strerror(nf_status)
        print *,'nf_open ',cfname(1:nf)
        return
      endif

      nf_status = nf_inq_varid(nf_fid,'image',nf_vid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in var image'
        return
      endif
      nf_status = nf_get_var_int(nf_fid,nf_vid,image)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         print *,'in nf_get_var_ image '
        return
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      return
      end
