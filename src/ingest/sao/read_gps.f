cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis
c
c  subroutine to read the file 
c
      subroutine read_gps(nf_fid, recnum, 
     +     pressure, staelev, stalat, stalon, temperature, 
     +     relativehumidity, timeobs, stanam)
c
      include 'netcdf.inc'
      integer recnum,nf_fid, nf_vid, nf_status

      real drydelay(recnum), formalerror(recnum), pressure(recnum),
     +     staelev(recnum), stalat(recnum), stalon(recnum),
     +     temperature(recnum), relativehumidity(recnum), 
     +     totaldelay(recnum),
     +     watervapor(recnum), wetdelay(recnum)
      double precision timeobs(recnum)
      character*5 stanam(recnum)
      character*80 stalongnam(recnum)


c   variables of type real
c
c     variable        netcdf long name
c      drydelay     "dry component gps signal delay"
c
        nf_status = nf_inq_varid(nf_fid,'drydelay',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drydelay'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,drydelay)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var drydelay'
      endif
c
c     variable        netcdf long name
c      formalerror  "formal error"
c
        nf_status = nf_inq_varid(nf_fid,'formalerror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var formalerror'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,formalerror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var formalerror'
      endif
c
c     variable        netcdf long name
c      pressure     "pressure used for pwv calculation"
c
        nf_status = nf_inq_varid(nf_fid,'pressure',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pressure'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,pressure)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var pressure'
      endif
c
c     variable        netcdf long name
c      staelev      "station elevation (above msl)"
c
        nf_status = nf_inq_varid(nf_fid,'staelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staelev'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,staelev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staelev'
      endif
c
c     variable        netcdf long name
c      stalat       "station latitude"
c
        nf_status = nf_inq_varid(nf_fid,'stalat',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalat'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stalat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalat'
      endif
c
c     variable        netcdf long name
c      stalon       "station longitude"
c
        nf_status = nf_inq_varid(nf_fid,'stalon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalon'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stalon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalon'
      endif
c
c     variable        netcdf long name
c      temperature  "temperature used for calculation"
c
        nf_status = nf_inq_varid(nf_fid,'temperature',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,temperature)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var temperature'
      endif
c
c     variable        netcdf long name
c      rh  "percent"
c
        nf_status = nf_inq_varid(nf_fid,'relativehumidity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relativehumidity'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,relativehumidity)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var relativehumidity'
      endif
c
c     variable        netcdf long name
c      totaldelay   "total gps signal delay"
c
        nf_status = nf_inq_varid(nf_fid,'totaldelay',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var totaldelay'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,totaldelay)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var totaldelay'
      endif
c
c     variable        netcdf long name
c      watervapor   "water vapor"
c
        nf_status = nf_inq_varid(nf_fid,'watervapor',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var watervapor'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,watervapor)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var watervapor'
      endif
c
c     variable        netcdf long name
c      wetdelay     "wet component gps signal delay"
c
        nf_status = nf_inq_varid(nf_fid,'wetdelay',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wetdelay'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wetdelay)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wetdelay'
      endif

c   variables of type int
c

c   variables of type double
c
c
c     variable        netcdf long name
c      timeobs      "time of observation"
c
        nf_status = nf_inq_varid(nf_fid,'timeobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,timeobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      stalongnam   "station location"
c
        nf_status = nf_inq_varid(nf_fid,'stalongnam',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalongnam'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,stalongnam)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalongnam'
      endif
c
c     variable        netcdf long name
c      stanam       "alphanumeric station name"
c
        nf_status = nf_inq_varid(nf_fid,'stanam',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stanam'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,stanam)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stanam'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
