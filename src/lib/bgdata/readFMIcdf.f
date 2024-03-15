      subroutine get_bgdata_data
     +                   (i4time_sys,ilaps_cycle_time,nx_l,ny_l
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer numisolevel, record, x, y,nf_fid, nf_vid, nf_status
c
c  open netcdf file for reading
c
      nf_status=nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),filename
        istatus=0
        return
      endif
c
c  fill all dimension values
c
c
c get size of numisolevel
c
      nf_status=nf_inq_dimid(nf_fid,'numisolevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim numisolevel'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,numisolevel)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim numisolevel'
      endif
c
c get size of record
c
      nf_status=nf_inq_dimid(nf_fid,'record',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim record'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,record)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim record'
      endif
c
c get size of x
c
      nf_status=nf_inq_dimid(nf_fid,'x',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,x)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim x'
      endif
c
c get size of y
c
      nf_status=nf_inq_dimid(nf_fid,'y',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
      endif
      nf_status=nf_inq_dimlen(nf_fid,nf_vid,y)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim y'
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
c
c
c
      subroutine read_fmi_netcdf(filename, numisolevel, record, x, y, 
     +     isolevel, dpt_2m, gh, gh_s, msl, p, pvv, sh, t, t_2m, t_s,
     +     uw, uw_10m, vw, vw_10m, model, origin, reftime, valtime,
     +     istatus)
c
      include 'netcdf.inc'
      integer numisolevel, record, x, y,nf_fid, nf_vid, nf_status
      integer isolevel(numisolevel)
      real dpt_2m( x,  y, record), gh( x,  y,  numisolevel, record),
     +     gh_s( x,  y, record), msl( x,  y, record), p( x,  y,
     +     record), pvv( x,  y,  numisolevel, record), sh( x,  y, 
     +     numisolevel, record), t( x,  y,  numisolevel, record),
     +     t_2m( x,  y, record), uw( x,  y,  numisolevel, record),
     +     uw_10m( x,  y, record), vw( x,  y,  numisolevel, record),
     +     vw_10m( x,  y, record), t_s(x,  y,  record),
     +     sst(x,  y,  record), skt(x,  y,  record), lsm(x, y, record)
      double precision reftime(record), valtime(record)
      character*132 origin
      character*132 model
      character*255 filename

      istatus =0
c
c  open netcdf file for reading
c
      nf_status=nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),trim(filename)
        istatus=1
        return
      endif

c   variables of type real
c
c     variable        netcdf long name
c     dpt_2m        "2 metre dewpoint temperature"
c
      nf_status=nf_inq_varid(nf_fid,'dpt_2m',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for dpt_2m'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,dpt_2m)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for dpt_2m'
       endif
      endif
c
c     variable        netcdf long name
c     gh            "geopotential"
c
      nf_status=nf_inq_varid(nf_fid,'gh',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for gh'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,gh)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for gh'
       endif
      endif
c convert to geometric height
c
      gh=gh/9.81
c
c     variable        netcdf long name
c     gh_s          "geopotential at surface"
c
      nf_status=nf_inq_varid(nf_fid,'gh_s',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for gh_s'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,gh_s)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for gh_s'
       endif
      endif
c
c     variable        netcdf long name
c     msl           "pressure at mean sea level"
c
      nf_status=nf_inq_varid(nf_fid,'msl',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for msl'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,msl)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for msl'
       endif
      endif
c
c     variable        netcdf long name
c     p             "pressure at surface"
c
      nf_status=nf_inq_varid(nf_fid,'p',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for p'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,p)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for p'
       endif
      endif
c
c     variable        netcdf long name
c     pvv           "vertical velocity"
c
      nf_status=nf_inq_varid(nf_fid,'pvv',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for pvv'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,pvv)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for pvv'
       endif
      endif
c
c     variable        netcdf long name
c     sh            "specific humidity"
c
      nf_status=nf_inq_varid(nf_fid,'sh',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sh'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,sh)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sh'
       endif
      endif
c
c     variable        netcdf long name
c     t             "temperature"
c
      nf_status=nf_inq_varid(nf_fid,'t',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for t'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,t)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for t'
       endif
      endif
c
c     variable        netcdf long name
c     t_2m          "2 metre temperature"
c
      nf_status=nf_inq_varid(nf_fid,'t_2m',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for t_2m'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,t_2m)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for t_2m'
       endif
      endif
c
c     variable        netcdf long name
c     lsm          "land/sea mask"
c
      nf_status=nf_inq_varid(nf_fid,'lsm',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for lsm'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,lsm)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for lsm'
       endif
      endif
c
c     variable        netcdf long name
c     sst          "sea surface temperature"
c
      nf_status=nf_inq_varid(nf_fid,'sst',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for sst'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,sst)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for sst'
       endif
      endif
c
c     variable        netcdf long name
c     skt          "skin temperature"
c
      nf_status=nf_inq_varid(nf_fid,'skt',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for skt'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,skt)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for skt'
       endif
      endif
c
c prepare surface temperature array for return. sst for water
c points and skin temps for land points.
c
      do j=1,y
       do i=1,x
        if(lsm(i,j,1).eq.1)then
         t_s(i,j,1)=sst(i,j,1)
        else
         t_s(i,j,1)=skt(i,j,1)
        endif
       enddo
      enddo
c
c     variable        netcdf long name
c     uw            "eastward wind"
c
      nf_status=nf_inq_varid(nf_fid,'uw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for uw'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,uw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for uw'
       endif
      endif
c
c     variable        netcdf long name
c     uw_10m        "eastward wind"
c
      nf_status=nf_inq_varid(nf_fid,'uw_10m',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for uw_10m'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,uw_10m)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for uw_10m'
       endif
      endif
c
c     variable        netcdf long name
c     vw            "northward wind"
c
      nf_status=nf_inq_varid(nf_fid,'vw',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for vw'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,vw)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for vw'
       endif
      endif
c
c     variable        netcdf long name
c     vw_10m        "northward wind"
c
      nf_status=nf_inq_varid(nf_fid,'vw_10m',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for vw_10m'
      else
       nf_status=nf_get_var_real(nf_fid,nf_vid,vw_10m)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for vw_10m'
       endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c     isolevel      "isobaric levels"
c
      nf_status=nf_inq_varid(nf_fid,'isolevel',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for isolevel'
      else
       nf_status=nf_get_var_int(nf_fid,nf_vid,isolevel)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for isolevel'
       endif
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c     reftime       "reference time"
c
      nf_status=nf_inq_varid(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for reftime'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,reftime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for reftime'
       endif
      endif
c
c     variable        netcdf long name
c     valtime       "valid time"
c
      nf_status=nf_inq_varid(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for valtime'
      else
       nf_status=nf_get_var_double(nf_fid,nf_vid,valtime)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for valtime'
       endif
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c     model         "model name"
c
      nf_status=nf_inq_varid(nf_fid,'model',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for model'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,model)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for model'
       endif
      endif
c
c     variable        netcdf long name
c     origin        "origin name"
c
      nf_status=nf_inq_varid(nf_fid,'origin',nf_vid)
      if(nf_status.ne.nf_noerr) then
       print *, nf_strerror(nf_status),' for origin'
      else
       nf_status=nf_get_var_text(nf_fid,nf_vid,origin)
       if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status),' for origin'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
