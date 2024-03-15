
      subroutine get_vol_netcdf_hdr(filename,
     +        gater, gater_hi, gatev, gatev_hi, radialr, radialr_hi,
     +        radialv, radialv_hi, scanr, scanr_hi, scanv,
     +        scanv_hi,nf_fid, nf_vid, nf_status)

!     argument list
      character*(*) filename

      include 'netcdf.inc'
      integer gater, gater_hi, gatev, gatev_hi, radialr, radialr_hi,
     +     radialv, radialv_hi, scanr, scanr_hi, scanv,
     +     scanv_hi,nf_fid, nf_vid, nf_status
c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open ',filename
      endif
c
c  fill all dimension values
c
c
c get size of gater
c
      nf_status = nf_inq_dimid(nf_fid,'gater',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gater'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,gater)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gater'
        gater = 0
      endif
c
c get size of gater_hi
c
      nf_status = nf_inq_dimid(nf_fid,'gater_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gater_hi'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,gater_hi)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gater_hi'
        gater_hi = 0
      endif
c
c get size of gatev
c
      nf_status = nf_inq_dimid(nf_fid,'gatev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gatev'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,gatev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gatev'
        gatev = 0
      endif
c
c get size of gatev_hi
c
      nf_status = nf_inq_dimid(nf_fid,'gatev_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gatev_hi'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,gatev_hi)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim gatev_hi'
        gatev_hi = 0
      endif
c
c get size of radialr
c
      nf_status = nf_inq_dimid(nf_fid,'radialr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialr'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,radialr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialr'
        radialr = 0
      endif
c
c get size of radialr_hi
c
      nf_status = nf_inq_dimid(nf_fid,'radialr_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialr_hi'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,radialr_hi)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialr_hi'
        radialr_hi = 0
      endif
c
c get size of radialv
c
      nf_status = nf_inq_dimid(nf_fid,'radialv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialv'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,radialv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialv'
        radialv = 0
      endif
c
c get size of radialv_hi
c
      nf_status = nf_inq_dimid(nf_fid,'radialv_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialv_hi'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,radialv_hi)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim radialv_hi'
        radialv_hi = 0
      endif
c
c get size of scanr
c
      nf_status = nf_inq_dimid(nf_fid,'scanr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanr'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,scanr)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanr'
        scanr = 0
      endif
c
c get size of scanr_hi
c
      nf_status = nf_inq_dimid(nf_fid,'scanr_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanr_hi'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,scanr_hi)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanr_hi'
        scanr_hi = 0
      endif
c
c get size of scanv
c
      nf_status = nf_inq_dimid(nf_fid,'scanv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanv'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,scanv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanv'
        scanv = 0
      endif
c
c get size of scanv_hi
c
      nf_status = nf_inq_dimid(nf_fid,'scanv_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanv_hi'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,scanv_hi)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim scanv_hi'
        scanv_hi = 0
      endif
!     call get_vol_netcdf_data(nf_fid, gater, gater_hi, gatev, gatev_hi, radialr,
!    +     radialr_hi, radialv, radialv_hi, scanr, scanr_hi, scanv,
!    +     scanv_hi)

      return
      end
c
c
      subroutine get_vol_netcdf_data(nf_fid, gater, gater_hi, gatev, 
     +     gatev_hi,
     +     radialr, radialr_hi, radialv, radialv_hi, scanr, scanr_hi,
     +     scanv, scanv_hi,
     +     reflectivity, reflectivity_hi,
     +     radialvelocity, radialvelocity_hi, 
     +     elevationr, elevationr_hi,
     +     elevationv, elevationv_hi,
     +     azimuthr, azimuthr_hi,
     +     azimuthv, azimuthv_hi,
     +     distancer, distancer_hi,
     +     distancev, distancev_hi,
     +     nyquistvelocityv, nyquistvelocityv_hi)

      include 'netcdf.inc'
      integer gater, gater_hi, gatev, gatev_hi, radialr, radialr_hi,
     +     radialv, radialv_hi, scanr, scanr_hi, scanv,
     +     scanv_hi,nf_fid, nf_vid, nf_status
      integer radialvelocity( gatev,  radialv, scanv),
     +     radialvelocity_hi( gatev_hi,  radialv_hi, scanv_hi),
     +     reflectivity( gater,  radialr, scanr), reflectivity_hi(
     +     gater_hi,  radialr_hi, scanr_hi), spectrumwidth( gatev, 
     +     radialv, scanv), spectrumwidth_hi( gatev_hi,  radialv_hi,
     +     scanv_hi), numgatesr(scanr), numgatesr_hi(scanr_hi),
     +     numgatesv(scanv), numgatesv_hi(scanv_hi),
     +     numradialsr(scanr), numradialsr_hi(scanr_hi),
     +     numradialsv(scanv), numradialsv_hi(scanv_hi), timer(
     +     radialr, scanr), timer_hi( radialr_hi, scanr_hi), timev(
     +     radialv, scanv), timev_hi( radialv_hi, scanv_hi)
      real azimuthr( radialr, scanr), azimuthr_hi( radialr_hi,
     +     scanr_hi), azimuthv( radialv, scanv), azimuthv_hi(
     +     radialv_hi, scanv_hi), distancer(gater),
     +     distancer_hi(gater_hi), distancev(gatev),
     +     distancev_hi(gatev_hi), elevationr( radialr, scanr),
     +     elevationr_hi( radialr_hi, scanr_hi), elevationv( radialv,
     +     scanv), elevationv_hi( radialv_hi, scanv_hi),
     +     nyquistvelocityv(scanv), nyquistvelocityv_hi(scanv_hi)


      call read_netcdf_vol(nf_fid, gater, gater_hi, gatev, gatev_hi, 
     +     radialr, radialr_hi, radialv, radialv_hi, scanr, scanr_hi, 
     +     scanv, scanv_hi, radialvelocity, radialvelocity_hi, 
     +     reflectivity, reflectivity_hi, spectrumwidth, 
     +     spectrumwidth_hi, numgatesr, numgatesr_hi, numgatesv, 
     +     numgatesv_hi, numradialsr, numradialsr_hi, numradialsv, 
     +     numradialsv_hi, timer, timer_hi, timev, timev_hi, 
     +     azimuthr, azimuthr_hi, azimuthv, azimuthv_hi, distancer, 
     +     distancer_hi, distancev, distancev_hi, elevationr, 
     +     elevationr_hi, elevationv, elevationv_hi,
     +     nyquistvelocityv, nyquistvelocityv_hi)
c
c the netcdf variables are filled - your code goes here
c
      return
      end
c
c  subroutine to read the file 
c
      subroutine read_netcdf_vol(nf_fid, gater, gater_hi, gatev, 
     +     gatev_hi, 
     +     radialr, radialr_hi, radialv, radialv_hi, scanr, scanr_hi, 
     +     scanv, scanv_hi, radialvelocity, radialvelocity_hi, 
     +     reflectivity, reflectivity_hi, spectrumwidth, 
     +     spectrumwidth_hi, numgatesr, numgatesr_hi, numgatesv, 
     +     numgatesv_hi, numradialsr, numradialsr_hi, numradialsv, 
     +     numradialsv_hi, timer, timer_hi, timev, timev_hi, 
     +     azimuthr, azimuthr_hi, azimuthv, azimuthv_hi, distancer, 
     +     distancer_hi, distancev, distancev_hi, elevationr, 
     +     elevationr_hi, elevationv, elevationv_hi,
     +     nyquistvelocityv, nyquistvelocityv_hi)
c
   
      use mem_namelist, only: r_missing_data

      include 'netcdf.inc'
      integer gater, gater_hi, gatev, gatev_hi, radialr, radialr_hi, 
     +     radialv, radialv_hi, scanr, scanr_hi, scanv, 
     +     scanv_hi,nf_fid, nf_vid, nf_status
      integer radialvelocity( gatev,  radialv, scanv),
     +     radialvelocity_hi( gatev_hi,  radialv_hi, scanv_hi),
     +     reflectivity( gater,  radialr, scanr), reflectivity_hi(
     +     gater_hi,  radialr_hi, scanr_hi), spectrumwidth( gatev, 
     +     radialv, scanv), spectrumwidth_hi( gatev_hi,  radialv_hi,
     +     scanv_hi), numgatesr(scanr), numgatesr_hi(scanr_hi),
     +     numgatesv(scanv), numgatesv_hi(scanv_hi),
     +     numradialsr(scanr), numradialsr_hi(scanr_hi),
     +     numradialsv(scanv), numradialsv_hi(scanv_hi), timer(
     +     radialr, scanr), timer_hi( radialr_hi, scanr_hi), timev(
     +     radialv, scanv), timev_hi( radialv_hi, scanv_hi)
      real azimuthr( radialr, scanr), azimuthr_hi( radialr_hi,
     +     scanr_hi), azimuthv( radialv, scanv), azimuthv_hi(
     +     radialv_hi, scanv_hi), distancer(gater),
     +     distancer_hi(gater_hi), distancev(gatev),
     +     distancev_hi(gatev_hi), elevationr( radialr, scanr),
     +     elevationr_hi( radialr_hi, scanr_hi), elevationv( radialv,
     +     scanv), elevationv_hi( radialv_hi, scanv_hi),
     +     nyquistvelocityr(scanr), nyquistvelocityr_hi(scanr_hi),
     +     nyquistvelocityv(scanv), nyquistvelocityv_hi(scanv_hi)


      i_missing_data = 255

c   variables of type real
c
c     variable        netcdf long name
c      azimuthr     "azimuth angle in degrees: 0 = true north, 90 = east"
c
      nf_status = nf_inq_varid(nf_fid,'azimuthr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var azimuthr'
        azimuthr = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,azimuthr)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var azimuthr'
          azimuthr = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      azimuthr_hi  "azimuth angle in degrees: 0 = true north, 90 = east"
c
      nf_status = nf_inq_varid(nf_fid,'azimuthr_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var azimuthr_hi'
        azimuthr_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,azimuthr_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var azimuthr_hi'
          azimuthr_hi = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      azimuthv     "azimuth angle in degrees: 0 = true north, 90 = east"
c
      nf_status = nf_inq_varid(nf_fid,'azimuthv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var azimuthv'
        azimuthv = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,azimuthv)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var azimuthv'
          azimuthv = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      azimuthv_hi  "azimuth angle in degrees: 0 = true north, 90 = east"
c
      nf_status = nf_inq_varid(nf_fid,'azimuthv_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var azimuthv_hi'
        azimuthv_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,azimuthv_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var azimuthv_hi'
          azimuthv_hi = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      distancer    "radial distance to start of gate"
c
      nf_status = nf_inq_varid(nf_fid,'distancer',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var distancer'
        distancer = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,distancer)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var distancer'
          distancer = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      distancer_hi "radial distance to start of gate"
c
      nf_status = nf_inq_varid(nf_fid,'distancer_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var distancer_hi'
        distancer_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,distancer_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var distancer_hi'
          distancer_hi = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      distancev    "radial distance to start of gate"
c
      nf_status = nf_inq_varid(nf_fid,'distancev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var distancev'
        distancev = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,distancev)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var distancev'
          distancev = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      distancev_hi "radial distance to start of gate"
c
      nf_status = nf_inq_varid(nf_fid,'distancev_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var distancev_hi'
        distancev_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,distancev_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var distancev_hi'
          distancev_hi = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      elevationr   "elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
c
      nf_status = nf_inq_varid(nf_fid,'elevationr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationr'
        elevationr = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,elevationr)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var elevationr'
          elevationr = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      elevationr_hi"elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
c
      nf_status = nf_inq_varid(nf_fid,'elevationr_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationr_hi'
        elevationr_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,elevationr_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var elevationr_hi'
          elevationr_hi = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      elevationv   "elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
c
      nf_status = nf_inq_varid(nf_fid,'elevationv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationv'
        elevationv = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,elevationv)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var elevationv'
          elevationv = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      elevationv_hi"elevation angle in degres: 0 = parallel to pedestal base, 90 = perpendicular"
c
      nf_status = nf_inq_varid(nf_fid,'elevationv_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var elevationv_hi'
        elevationv_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,elevationv_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var elevationv_hi'
          elevationv_hi = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      nyquistvelocityr"nyquist velocity"
c
      nf_status = nf_inq_varid(nf_fid,'nyquistvelocityr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyquistvelocityr'
        nyquistvelocityr = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,nyquistvelocityr)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var nyquistvelocityr'
          nyquistvelocityr = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      nyquistvelocityr_hi"nyquist velocity"
c
      nf_status = nf_inq_varid(nf_fid,'nyquistvelocityr_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyquistvelocityr_hi'
        nyquistvelocityr_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,nyquistvelocityr_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var nyquistvelocityr_hi'
          nyquistvelocityr_hi = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      nyquistvelocityv"nyquist velocity"
c
      nf_status = nf_inq_varid(nf_fid,'nyquistvelocityv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyquistvelocityv'
        nyquistvelocityv = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,nyquistvelocityv)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var nyquistvelocityv'
          nyquistvelocityv = r_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      nyquistvelocityv_hi"nyquist velocity"
c
      nf_status = nf_inq_varid(nf_fid,'nyquistvelocityv_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var nyquistvelocityv_hi'
        nyquistvelocityv_hi = r_missing_data
      else
        nf_status = nf_get_var_real(nf_fid,nf_vid,nyquistvelocityv_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var nyquistvelocityv_hi'
          nyquistvelocityv_hi = r_missing_data
        endif
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      radialvelocity"radial velocity"
c
      nf_status = nf_inq_varid(nf_fid,'radialvelocity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialvelocity'
        radialvelocity = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,radialvelocity)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var radialvelocity'
          radialvelocity = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      radialvelocity_hi"radial velocity_hi"
c
      nf_status = nf_inq_varid(nf_fid,'radialvelocity_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var radialvelocity_hi'
        radialvelocity_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,radialvelocity_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var radialvelocity_hi'
          radialvelocity_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      reflectivity "reflectivity"
c
      nf_status = nf_inq_varid(nf_fid,'reflectivity',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reflectivity'
        reflectivity = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,reflectivity)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var reflectivity'
          reflectivity = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      reflectivity_hi"reflectivity_hi"
c
      nf_status = nf_inq_varid(nf_fid,'reflectivity_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var reflectivity_hi'
        reflectivity_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,reflectivity_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var reflectivity_hi'
          reflectivity_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      spectrumwidth"radial spectrum"
c
      nf_status = nf_inq_varid(nf_fid,'spectrumwidth',nf_vid)
      if(nf_status.ne.nf_noerr .or. .true.) then
        print *, nf_strerror(nf_status)
        print *,'in var spectrumwidth'
        spectrumwidth = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,spectrumwidth)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var spectrumwidth'
          spectrumwidth = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      spectrumwidth_hi"radial spectrum_hi"
c
      nf_status = nf_inq_varid(nf_fid,'spectrumwidth_hi',nf_vid)
      if(nf_status.ne.nf_noerr .or. .true.) then
        print *, nf_strerror(nf_status)
        print *,'in var spectrumwidth_hi'
        spectrumwidth_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,spectrumwidth_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var spectrumwidth_hi'
          spectrumwidth_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numgatesr    "number of valid gates in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numgatesr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesr'
        numgatesr = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numgatesr)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numgatesr'
          numgatesr = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numgatesr_hi "number of valid gates in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numgatesr_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesr_hi'
        numgatesr_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numgatesr_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numgatesr_hi'
          numgatesr_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numgatesv    "number of valid gates in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numgatesv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesv'
        numgatesv = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numgatesv)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numgatesv'
          numgatesv = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numgatesv_hi "number of valid gates in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numgatesv_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numgatesv_hi'
        numgatesv_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numgatesv_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numgatesv_hi'
          numgatesv_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numradialsr  "number of valid radials in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numradialsr',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numradialsr'
        numradialsr = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numradialsr)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numradialsr'
          numradialsr = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numradialsr_hi"number of valid radials in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numradialsr_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numradialsr_hi'
        numradialsr_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numradialsr_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numradialsr_hi'
          numradialsr_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numradialsv  "number of valid radials in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numradialsv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numradialsv'
        numradialsv = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numradialsv)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numradialsv'
          numradialsv = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      numradialsv_hi"number of valid radials in this scan"
c
      nf_status = nf_inq_varid(nf_fid,'numradialsv_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numradialsv_hi'
        numradialsv_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,numradialsv_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var numradialsv_hi'
          numradialsv_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      timer        "time since base date"
c
      nf_status = nf_inq_varid(nf_fid,'timer',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timer'
        timer = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,timer)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var timer'
          timer = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      timer_hi     "time since base date"
c
      nf_status = nf_inq_varid(nf_fid,'timer_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timer_hi'
        timer_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,timer_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var timer_hi'
          timer_hi = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      timev        "time since base date"
c
      nf_status = nf_inq_varid(nf_fid,'timev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timev'
        timev = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,timev)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var timev'
          timev = i_missing_data
        endif
      endif
c
c     variable        netcdf long name
c      timev_hi     "time since base date"
c
      nf_status = nf_inq_varid(nf_fid,'timev_hi',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timev_hi'
        timev_hi = i_missing_data
      else
        nf_status = nf_get_var_int(nf_fid,nf_vid,timev_hi)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'in var timev_hi'
          timev_hi = i_missing_data
        endif
      endif

c   variables of type double
c


c   variables of type char
c

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
