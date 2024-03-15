      subroutine readavnpublicdims(fname,x,y,numisolevel,record,
     +istatus)

      include 'netcdf.inc'
      integer       numisolevel
      integer       record, x, y
      integer       nf_fid, nf_vid, nf_status
      integer       istatus
      character*200 fname
c
      istatus=0
c
c  open netcdf file for reading
c
      nf_status = nf_open(fname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        call s_len(fname,len)
        print *, nf_strerror(nf_status)
        print *,'nf_open: ',fname(1:len)
        return
      endif
c
c  fill all dimension values
c
c
c get size of numisolevel
c
      nf_status = nf_inq_dimid(nf_fid,'numisolevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        nf_status = nf_inq_dimid(nf_fid,'isolevel',nf_vid)
        if(nf_status.ne.nf_noerr) then
           print *, nf_strerror(nf_status)
           print *,'dim numisolevel'
           return
        endif
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,numisolevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim numisolevel'
        return
      endif
c
c get size of record
c
      nf_status = nf_inq_dimid(nf_fid,'record',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim record'
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,record)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim record'
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

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
        return
      endif

      istatus=1

c     call main_sub(nf_fid, numisolevel, record, x, y)

      return
      end
c
c *********************************************************************
c  subroutine to read the file "global 26 layer spectral aviation model" 
c
      subroutine read_avn_netcdf(fname, numisolevel, record, x, y, 
     +     version, accs, gh, gh_s, p, pmsl, pw, rh, t, t_2m, t_s,
     +     uw, vw, ww, rh_2m, isolevel, reftime, valtime, grid, model,
     +     nav, origin, istatus)
c
      include 'netcdf.inc'
      integer numisolevel, record, x, y,nf_fid, nf_vid, nf_status
      integer version
      real accs( x,  y, record), gh( x,  y,  numisolevel, record),
     +     gh_s( x,  y, record), p( x,  y, record), pmsl( x,  y,
     +     record), pw( x,  y, record), rh( x,  y,  numisolevel,
     +     record), t( x,  y,  numisolevel, record), t_s( x,  y,
     +     record), t_2m( x,  y, record), uw( x,  y,  numisolevel,
     +     record), vw( x,  y, numisolevel, record), rh_2m( x, y,
     +     record), ww( x,  y, numisolevel, record)
      double precision isolevel(numisolevel), reftime(record),
     +     valtime(record)
      character*132 origin
      character*132 model
      character*132 nav
      character*132 grid
      character*255 fname
c
c
      istatus = 1
c
c  open netcdf file for reading
c
      nf_status = nf_open(fname,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        call s_len(fname,len)
        print *, nf_strerror(nf_status)
        print *,'nf_open: ',fname(1:len)
        return
      endif
c
c     variables of type real
c
c     variable        netcdf long name
c      accs         "accumulated snow"
c
      nf_status = nf_inq_varid(nf_fid,'accs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var accs'
c       return
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,accs)
        if(nf_status.ne.nf_noerr) then
           print *, nf_strerror(nf_status)
           print *,'in var accs'
           return
        endif
      endif
c
c     variable        netcdf long name
c      gh           "geopotential height"
c
      nf_status = nf_inq_varid(nf_fid,'gh',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gh'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,gh)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gh'
        return
      endif
c
c     variable        netcdf long name
c      gh_s         "geopotential height of surface"
c
      nf_status = nf_inq_varid(nf_fid,'gh_s',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gh_s'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,gh_s)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gh_s'
        return
      endif
c
c     variable        netcdf long name
c      p            "pressure at surface"
c
      nf_status = nf_inq_varid(nf_fid,'p',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var p'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,p)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var p'
        return
      endif
c
c     variable        netcdf long name
c      pmsl         "pressure at mean sea level"
c
      nf_status = nf_inq_varid(nf_fid,'pmsl',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pmsl'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,pmsl)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pmsl'
        return
      endif
c
c     variable        netcdf long name
c      pw           "precipitable water"
c
      nf_status = nf_inq_varid(nf_fid,'pw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pw'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,pw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pw'
        return
      endif
c
c     variable        netcdf long name
c      rh           "relative humidity"
c
      nf_status = nf_inq_varid(nf_fid,'rh',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rh'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,rh)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rh'
        return
      endif
c
c     variable        netcdf long name
c      t            "temperature"
c
      nf_status = nf_inq_varid(nf_fid,'t',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,t)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t'
        return
      endif
c
c     variable        netcdf long name
c      t_2m          "temperature at 2 meters above surface"
c
      nf_status = nf_inq_varid(nf_fid,'t_2m',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t_2m'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,t_2m)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t_2m'
        return
      endif
c
c     variable        netcdf long name
c      t_s          "temperature at surface"
c
      nf_status = nf_inq_varid(nf_fid,'t_s',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t_s'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,t_s)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var t_s'
        return
      endif
c
c     variable        netcdf long name
c      uw           "eastward wind"
c
      nf_status = nf_inq_varid(nf_fid,'uw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var uw'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,uw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var uw'
        return
      endif
c
c     variable        netcdf long name
c      vw           "northward wind"

      nf_status = nf_inq_varid(nf_fid,'vw',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vw'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,vw)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vw'
        return
      endif
c
c     variable        netcdf long name
c      rh_2m          "relative humidity at 2 meters above surface"
c
      nf_status = nf_inq_varid(nf_fid,'rh_2m',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rh_2m'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,rh_2m)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rh_2m'
        return
      endif
c
c     variable        netcdf long name
c      ww           "pressure vertical velocity"

      nf_status = nf_inq_varid(nf_fid,'pvv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pvv (pres vert velocity)'
        return
      endif
      nf_status = nf_get_var_real(nf_fid,nf_vid,ww)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var ww'
        return
      endif
c
c   variables of type int
c
c
c     variable        netcdf long name
c      version      
c
      nf_status = nf_inq_varid(nf_fid,'version',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var version'
        return
      endif
      nf_status = nf_get_var_int(nf_fid,nf_vid,version)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var version'
        return
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      isolevel     "isobaric levels"
c
      nf_status = nf_inq_varid(nf_fid,'isolevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var isolevel'
        return
      endif
      nf_status = nf_get_var_double(nf_fid,nf_vid,isolevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var isolevel'
        return
      endif
c
c     variable        netcdf long name
c      reftime      "reference time"
c
      nf_status = nf_inq_varid(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reftime'
        return
      endif
      nf_status = nf_get_var_double(nf_fid,nf_vid,reftime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reftime'
        return
      endif
c
c     variable        netcdf long name
c      valtime      "valid time"
c
      nf_status = nf_inq_varid(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var valtime'
        return
      endif
      nf_status = nf_get_var_double(nf_fid,nf_vid,valtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var valtime'
        return
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      grid_type 
c
      nf_status = nf_inq_varid(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,grid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var grid'
        return
      endif
c
c     variable        netcdf long name
c      model        
c
      nf_status = nf_inq_varid(nf_fid,'model',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var model'
        return
      endif
      nf_status = nf_get_var_text(nf_fid,nf_vid,model)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var model'
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
