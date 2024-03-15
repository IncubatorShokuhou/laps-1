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
        subroutine  get_scandata(
     :          i_tilt,             ! input  (integer)
     :          r_missing_data,     ! input  (real)
     :          max_rays_in,        ! input  (integer)
     :          n_rays,             ! output (integer)
     :          n_gates,            ! output (integer)
     :          gate_spacing_m_ret, ! output (real)
     :          elevation_deg,      ! output (real)
     :          i_scan_mode,        ! output (integer)
     :          v_nyquist_tilt,     ! output (real)
     :          v_nyquist_ray,      ! output (real array)
     :          istatus )           ! output (integer)

c
c     purpose:
c
c       get info for radar data.
c
      implicit none
c
c     input variables
c
      integer i_tilt
      real r_missing_data
      integer max_rays_in
c
c     output variables
c
      integer n_rays
      integer n_gates
      real gate_spacing_m_ret
      real elevation_deg
      integer i_scan_mode
      real v_nyquist_tilt
      real v_nyquist_ray(max_rays_in)
      integer istatus
c
c     include file
c
      include 'remap_dims.inc'
      include 'remap_constants.dat'
      include 'remap_buffer.cmn'
c
c     misc internal variables
c
      integer i
c
      n_rays = n_rays_cmn
      n_gates = max_gates
      gate_spacing_m_ret = gate_spacing_m
      elevation_deg = elev_cmn
      i_scan_mode = 1
c
c     transfer nyquist info from common to output array
c
      do 100 i = 1,n_rays
        if (v_nyquist_ray_a_cmn(i) .eq. 0.) then
          v_nyquist_ray(i) = r_missing_data
        else
          v_nyquist_ray(i) = v_nyquist_ray_a_cmn(i)
        end if
  100 continue
c
c     check all nyquists in this tilt.
c     if they are all the same, report that as the tilt nyquist,
c     else report r_missing
c
      v_nyquist_tilt = v_nyquist_ray(1)
      do 200 i = 1,n_rays
        if (v_nyquist_ray(i) .eq. r_missing_data) then
          v_nyquist_tilt = r_missing_data
          go to 999
        else if (v_nyquist_ray(i) .ne. v_nyquist_tilt) then
          v_nyquist_tilt = r_missing_data
          go to 999
        end if
  200 continue

  999 continue

      istatus = 1
      return
      end
