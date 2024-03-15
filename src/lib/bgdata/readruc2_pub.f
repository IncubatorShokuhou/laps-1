      subroutine get_ruc2_dims(filename,cmodel,nx,ny,nz
     &,stdlat1,stdlat2,lon0,la1,lo1,la2,lo2,istatus)
c
c updated to get all nav info - js 4-01
c
      implicit none
      include 'netcdf.inc'

      real    stdlat1,stdlat2
      real    lon0
      real    la1,lo1
      real    la2,lo2

      integer  nx, ny, nz, nf_fid, nf_vid, nf_status
      character filename*200
      character cmodel*(*)
      integer istatus

      istatus = 0
c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', filename
        return
      endif
c
c  fill all dimension values
c
c
c get size of record
c
c      nf_status = nf_inq_dimid(nf_fid,'record',nf_vid)
c      if(nf_status.ne.nf_noerr) then
c        print *, nf_strerror(nf_status)
c        print *,'dim record'
c      endif
c      nf_status = nf_inq_dimlen(nf_fid,nf_vid,record)
c      if(nf_status.ne.nf_noerr) then
c        print *, nf_strerror(nf_status)
c        print *,'dim record'
c      endif
c
c get size of x
c
      nf_status = nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
        return
      endif

      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nx)
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

      nf_status = nf_inq_dimlen(nf_fid,nf_vid,ny)
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

      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nz)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z'
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
      if(cmodel.eq.'ruc40_native')then
        nf_status = nf_inq_varid(nf_fid,'latin',nf_vid)
      else
        nf_status = nf_inq_varid(nf_fid,'latin1',nf_vid)
      endif
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stdlat1)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var latin - stdlat'
      endif

c not in netcdf file
      stdlat2=stdlat1
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
        nf_status = nf_get_var_real(nf_fid,nf_vid,lon0)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var lov - lon0'
      endif

c these not in netcdf file!
      la2=55.4818
      lo2=-57.3794

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close ruc2'
        return
      endif

      istatus=1 
      return
      end
c
c
c
c  subroutine to read the file "hybrid-b 40km rapid update cycle" 
c
      subroutine read_ruc2_hybb(filename, nx, ny, nz, mmsp
     +     ,hgt, p, qv, u, v, vpt, w,istatus)
      implicit none
      include 'netcdf.inc'
c     integer nx, ny, nz, nf_fid, nf_vid, nf_status,istatus
      integer nx, ny, nz, nf_fid, nf_status,istatus
      character*200 filename
      character*4   cvar
      integer nxny,nxnynz
      real mmsp(nx,ny), hgt( nx,  ny,  nz), 
     +     p( nx,  ny,  nz), qv( nx,  ny,  nz), 
     +     u( nx,  ny,  nz), v( nx,  ny,  nz), 
     +     vpt( nx,  ny,  nz), w( nx,  ny,  nz)

      istatus = 1
      nxny=nx*ny
      nxnynz=nx*ny*nz

c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ', filename
        return
      endif

c
c     variable        netcdf long name
c      mmsp         "maps mean sea level pressure" 
c
      cvar='mmsp'
      call read_netcdf_real(nf_fid,cvar,nxny,mmsp,0,0,nf_status)

c
c     variable        netcdf long name
c      hgt          "geopotential height" 
c
      cvar='hgt '
      call read_netcdf_real(nf_fid,cvar,nxnynz,hgt,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the height field is missing'
         print*, 'aborting file processing for ruc2 file ',filename
         return
      endif


c
c     variable        netcdf long name
c      p            "pressure" 
c
      cvar='p   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,p,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
      print*, 'a substantial portion of the pressure field is missing'
         print*, 'aborting file processing for ruc2 file ',filename
         return
      endif
c
c     variable        netcdf long name
c      qv           "water vapor mixing ratio" 
c
      cvar='qv  '
      call read_netcdf_real(nf_fid,cvar,nxnynz,qv,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
      print*,'a substantial portion of the water vapor field is missing'
         print*, 'aborting file processing for ruc2 file ',filename
         return
      endif

c
c     variable        netcdf long name
c      u            "u-component of wind" 
c
      cvar='u   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,u,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the u-wind field is missing'
         print*, 'aborting file processing for ruc2 file ',filename
         return
      endif

c
c     variable        netcdf long name
c      v            "v-component of wind" 
c
      cvar='v   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,v,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the v-wind field is missing'
         print*, 'aborting file processing for ruc2 file ',filename
         return
      endif

c
c     variable        netcdf long name
c      vpt          "virtual potential temperature" 
c
      cvar='vpt '
      call read_netcdf_real(nf_fid,cvar,nxnynz,vpt,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'a substantial portion of the vpt field is missing'
         print*, 'aborting file processing for ruc2 file ',filename
         return
      endif

c
c     variable        netcdf long name
c      w            "vertical velocity" 
c
      cvar='w   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,w,0,0,nf_status)


      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0
      return
      end
