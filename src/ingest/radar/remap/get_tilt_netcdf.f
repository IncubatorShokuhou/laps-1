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

      subroutine get_tilt_netcdf_hdr(filename,nf_fid
     1                               ,radarname
     1                               ,sitelat                        
     1                               ,sitelon                        
     1                               ,sitealt                        
     1                               ,elevationangle
     1                               ,numradials 
     1                               ,elevationnumber
     1                               ,vcp
     1                               ,radialazim
     1                               ,resolutionv
     1                               ,gatesizev,gatesizez
     1                               ,firstgaterangev,firstgaterangez
     1                               ,v_bin_max, z_bin_max, radial_max ! i
     1                               ,v_bin,     z_bin,     radial     ! o
     1                               ,istatus)

!     argument list
      character*(*) filename
      character*5  radarname
      integer v_bin_max, z_bin_max, radial_max
      real radialazim(radial_max)

!     local
      real radialelev(radial_max)
      character*132 sitename
      double precision esendtime, esstarttime, radialtime(radial_max)

!.............................................................................

      include 'netcdf.inc'
      integer v_bin, z_bin, radial,nf_fid, nf_vid, nf_status

!     include 'remap_constants.dat' ! for debugging only
!     include 'remap.cmn' ! for debugging only

      write(6,*)' get_tilt_netcdf_hdr: reading ',filename
c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',filename
        istatus = 0
        return
      endif
c
c  fill all dimension values
c
c
c get size of v_bin
c
      nf_status = nf_inq_dimid(nf_fid,'v_bin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim v_bin'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,v_bin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim v_bin'
      endif
c
c get size of z_bin
c
      nf_status = nf_inq_dimid(nf_fid,'z_bin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z_bin'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,z_bin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z_bin'
      endif
c
c get size of radial
c
      nf_status = nf_inq_dimid(nf_fid,'radial',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radial'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,radial)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radial'
      endif

!.....test whether dimensions of netcdf file are within bounds...............

      if(v_bin .gt. v_bin_max)then
          write(6,*)' v_bin > permitted dimensions ',v_bin,v_bin_max
          stop
      endif

      if(z_bin .gt. z_bin_max)then
          write(6,*)' z_bin > permitted dimensions ',z_bin,z_bin_max     
          stop
      endif

      if(radial .gt. radial_max)then
          write(6,*)' radial > permitted dimensions ',radial,radial_max       
          stop
      endif

      istatus = 1
      return

      end

      subroutine get_tilt_netcdf_data(filename,nf_fid
     1                               ,radarname
     1                               ,sitelat                        
     1                               ,sitelon                        
     1                               ,sitealt                        
     1                               ,elevationangle
     1                               ,numradials 
     1                               ,elevationnumber
     1                               ,vcp
     1                               ,nyquist
     1                               ,radialazim
     1                               ,z  
     1                               ,v
     1                               ,resolutionv
     1                               ,gatesizev,gatesizez
     1                               ,firstgaterangev,firstgaterangez
     1                               ,z_scale, z_offset
     1                               ,v_scale, v_offset
     1                               ,v_bin_in, z_bin_in, radial_in    ! i
     1                               ,istatus)

!     argument list
      character*(*) filename
      character*5  radarname
      integer v_bin_in, z_bin_in, radial_in
      integer v(v_bin_in,radial_in), z(z_bin_in,radial_in)
      real radialazim(radial_in)

!     local
      real radialelev(radial_in)
      character*132 sitename
      double precision esendtime, esstarttime, radialtime(radial_in)

!.............................................................................

      include 'netcdf.inc'
      integer v_bin, z_bin, radial,nf_fid, nf_vid, nf_status

!     include 'remap_constants.dat' ! for debugging only
!     include 'remap.cmn' ! for debugging only

      write(6,*)' get_tilt_netcdf_data: reading ',filename
      write(6,*)'                       v/z/rad = '
     1         ,v_bin_in,z_bin_in,radial_in
c
c  open netcdf file for reading
c
!     nf_status = nf_open(filename,nf_nowrite,nf_fid)
!     if(nf_status.ne.nf_noerr) then
!       print *, nf_strerror(nf_status)
!       print *,'nf_open ',filename
!       istatus = 0
!       return
!     endif
c
c  fill all dimension values
c
c
c get size of v_bin
c
      nf_status = nf_inq_dimid(nf_fid,'v_bin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim v_bin'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,v_bin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim v_bin'
      endif
c
c get size of z_bin
c
      nf_status = nf_inq_dimid(nf_fid,'z_bin',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z_bin'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,z_bin)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim z_bin'
      endif
c
c get size of radial
c
      nf_status = nf_inq_dimid(nf_fid,'radial',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radial'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,radial)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radial'
      endif

!.....test whether dimensions of netcdf file are within bounds...............

      if(v_bin .ne. v_bin_in)then
          write(6,*)' v_bin != permitted dimensions ',v_bin,v_bin_in
          stop
      endif

      if(z_bin .ne. z_bin_in)then
          write(6,*)' z_bin != permitted dimensions ',z_bin,z_bin_in     
          stop
      endif

      if(radial .ne. radial_in)then
          write(6,*)' radial != permitted dimensions ',radial,radial_in       
          stop
      endif

      call read_netcdf(nf_fid, v_bin_in, z_bin_in, radial_in,        ! i
!............................................................................
     +     v, vcp, 
     +     z, elevationnumber, numgatesv, numgatesz, numradials, 
     +     atmosattenfactor, calibconst, elevationangle, 
     +     firstgaterangev, firstgaterangez, gatesizev, gatesizez, 
     +     z_scale, z_offset, v_scale, v_offset,
     +     nyquist, powdiffthreshold, radialazim, radialelev, 
     +     resolutionv, sitealt, sitelat, sitelon, unambigrange, 
     +     esendtime, esstarttime, radialtime, radarname, sitename)

      if(sitealt .eq. 0. .or. sitelat .eq. 0. .or. sitelon .eq. 0.)then       
          write(6,*)' warning, no site info in get_tilt_netcdf_data'       
          istatus = 0
          return
      else
          write(6,*)' site info:',sitealt, sitelat, sitelon       
          istatus = 1
          return
      endif

      end
c
c
c
c  subroutine to read the file "wsr-88d wideband data" 
c
      subroutine read_netcdf(nf_fid, v_bin, z_bin, radial, v, vcp, 
     +     z, elevationnumber, numgatesv, numgatesz, numradials, 
     +     atmosattenfactor, calibconst, elevationangle, 
     +     firstgaterangev, firstgaterangez, gatesizev, gatesizez, 
     +     z_scale, z_offset, v_scale, v_offset,
     +     nyquist, powdiffthreshold, radialazim, radialelev, 
     +     resolutionv, sitealt, sitelat, sitelon, unambigrange, 
     +     esendtime, esstarttime, radialtime, radarname, sitename)
c
      include 'netcdf.inc'
!     include 'remap_constants.dat' ! for debugging only
!     include 'remap.cmn' ! for debugging only
      integer v_bin, z_bin, radial,nf_fid, nf_vid, nf_status
      integer v( v_bin, radial), vcp, z( z_bin,
     +     radial), elevationnumber, numgatesv, numgatesz, numradials
      real atmosattenfactor, calibconst, elevationangle,
     +     firstgaterangev, firstgaterangez, gatesizev, gatesizez,
     +     nyquist, powdiffthreshold, radialazim(radial),
     +     radialelev(radial), resolutionv, sitealt, sitelat,
     +     sitelon, unambigrange
      double precision esendtime, esstarttime, radialtime(radial)
      character*5 radarname
      character*132 sitename


c   variables of type real
c
c     variable        netcdf long name
c      atmosattenfactor"atmospheric attenuation factor"
c
        nf_status = nf_inq_varid(nf_fid,'atmosattenfactor',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var atmosattenfactor'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,atmosattenfactor)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var atmosattenfactor'
      endif
c
c     variable        netcdf long name
c      calibconst   "system gain calibration constant"
c
        nf_status = nf_inq_varid(nf_fid,'calibconst',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var calibconst'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,calibconst)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var calibconst'
      endif
c
c     variable        netcdf long name
c      elevationangle"elevation angle"
c
        nf_status = nf_inq_varid(nf_fid,'elevationangle',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationangle'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,elevationangle)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationangle'
      endif
c
c     variable        netcdf long name
c      firstgaterangev"range to 1st doppler gate"
c
        nf_status = nf_inq_varid(nf_fid,'firstgaterangev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstgaterangev'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,firstgaterangev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstgaterangev'
      endif
c
c     variable        netcdf long name
c      firstgaterangez"range to 1st reflectivity gate"
c
        nf_status = nf_inq_varid(nf_fid,'firstgaterangez',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstgaterangez'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,firstgaterangez)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var firstgaterangez'
      endif
c
c     variable        netcdf long name
c      gatesizev    "doppler gate spacing"
c
        nf_status = nf_inq_varid(nf_fid,'gatesizev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gatesizev'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,gatesizev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gatesizev'
      endif
c
c     variable        netcdf long name
c      gatesizez    "reflectivity gate spacing"
c
        nf_status = nf_inq_varid(nf_fid,'gatesizez',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gatesizez'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,gatesizez)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var gatesizez'
      endif
c
c     variable        netcdf long name
c      z_scale  "reflectivity scale value"
c
      nf_status = nf_inq_varid(nf_fid,'z_scale',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var z_scale'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,z_scale)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var z_scale'
        endif
      endif
c
c     variable        netcdf long name
c      z_offset  "reflectivity offset value"
c
      nf_status = nf_inq_varid(nf_fid,'z_offset',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var z_offset'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,z_offset)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var z_offset'
        endif
      endif
c
c     variable        netcdf long name
c      nyquist      "nyquist velocity"
c
        nf_status = nf_inq_varid(nf_fid,'nyquist',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyquist'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,nyquist)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyquist'
      endif
c
c     variable        netcdf long name
c      powdiffthreshold"range de-aliasing threshold"
c
        nf_status = nf_inq_varid(nf_fid,'powdiffthreshold',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var powdiffthreshold'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,powdiffthreshold)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var powdiffthreshold'
      endif
c
c     variable        netcdf long name
c      radialazim   "radial azimuth angle"
c
        nf_status = nf_inq_varid(nf_fid,'radialazim',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialazim'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,radialazim)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialazim'
      endif
c
c     variable        netcdf long name
c      radialelev   "radial elevation angle"
c
        nf_status = nf_inq_varid(nf_fid,'radialelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialelev'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,radialelev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialelev'
      endif
c
c     variable        netcdf long name
c      resolutionv  "doppler velocity resolution"
c
        nf_status = nf_inq_varid(nf_fid,'resolutionv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var resolutionv'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,resolutionv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var resolutionv'
      endif
c
c     variable        netcdf long name
c      v_scale  "velocity scale value"
c
      nf_status = nf_inq_varid(nf_fid,'v_scale',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var v_scale'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,v_scale)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var v_scale'
        endif
      endif
c
c     variable        netcdf long name
c      v_offset  "velocity offset value"
c
      nf_status = nf_inq_varid(nf_fid,'v_offset',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var v_offset'
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,v_offset)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var v_offset'
        endif
      endif
c
c     variable        netcdf long name
c      sitealt      "altitude of site above mean sea level"
c
        nf_status = nf_inq_varid(nf_fid,'sitealt',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitealt'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,sitealt)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitealt'
      endif
c
c     variable        netcdf long name
c      sitelat      "latitude of site"
c
        nf_status = nf_inq_varid(nf_fid,'sitelat',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitelat'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,sitelat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitelat'
      endif
c
c     variable        netcdf long name
c      sitelon      "longitude of site"
c
        nf_status = nf_inq_varid(nf_fid,'sitelon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitelon'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,sitelon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitelon'
      endif
c
c     variable        netcdf long name
c      unambigrange "unambiguous range"
c
        nf_status = nf_inq_varid(nf_fid,'unambigrange',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var unambigrange'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,unambigrange)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var unambigrange'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      v            "velocity"
c
        nf_status = nf_inq_varid(nf_fid,'v',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var v'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,v)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var v'
      endif
c
c     variable        netcdf long name
c      vcp          "volume coverage pattern"
c
        nf_status = nf_inq_varid(nf_fid,'vcp',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vcp'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,vcp)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var vcp'
      endif
c
c     variable        netcdf long name
c      w            "spectrum width"
c
!        nf_status = nf_inq_varid(nf_fid,'w',nf_vid)
!      if(nf_status.ne.nf_noerr) then
!        print *, nf_strerror(nf_status)
!        print *,'in var w'
!      endif
!        nf_status = nf_get_var_int(nf_fid,nf_vid,w)
!      if(nf_status.ne.nf_noerr) then
!        print *, nf_strerror(nf_status)
!        print *,'in var w'
!      endif
c
c     variable        netcdf long name
c      z            "reflectivity"
c
        nf_status = nf_inq_varid(nf_fid,'z',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var z'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,z)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var z'
      endif
c
c     variable        netcdf long name
c      elevationnumber"elevation number"
c
        nf_status = nf_inq_varid(nf_fid,'elevationnumber',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationnumber'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,elevationnumber)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationnumber'
      endif
c
c     variable        netcdf long name
c      numgatesv    "number of doppler gates"
c
        nf_status = nf_inq_varid(nf_fid,'numgatesv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesv'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,numgatesv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesv'
      endif
c
c     variable        netcdf long name
c      numgatesz    "number of reflectivity gates"
c
        nf_status = nf_inq_varid(nf_fid,'numgatesz',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesz'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,numgatesz)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesz'
      endif
c
c     variable        netcdf long name
c      numradials   "number of radials"
c
        nf_status = nf_inq_varid(nf_fid,'numradials',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numradials'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,numradials)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numradials'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      esendtime    "end time of elevation scan"
c
        nf_status = nf_inq_varid(nf_fid,'esendtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var esendtime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,esendtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var esendtime'
      endif
c
c     variable        netcdf long name
c      esstarttime  "start time of elevation scan"
c
        nf_status = nf_inq_varid(nf_fid,'esstarttime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var esstarttime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,esstarttime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var esstarttime'
      endif
c
c     variable        netcdf long name
c      radialtime   "time of radial"
c
        nf_status = nf_inq_varid(nf_fid,'radialtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialtime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,radialtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialtime'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      radarname    "official name of the radar"
c
        nf_status = nf_inq_varid(nf_fid,'radarname',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radarname'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,radarname)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radarname'
      endif
c
c     variable        netcdf long name
c      sitename     "long name of the radar site"
c
        nf_status = nf_inq_varid(nf_fid,'sitename',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitename'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,sitename)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var sitename'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
