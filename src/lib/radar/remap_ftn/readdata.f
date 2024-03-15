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
      subroutine read_data_88d(
     :               i_tilt,           ! input
     :               vel_thr_rtau,
     :               r_missing_data,
     :               namelist_parms,
     :               gate_spacing_m_ret,   ! output
     :               i_scan_mode,
     :               num_sweeps,
     :               elevation_deg,
     :               n_rays,
     :               n_gates,
     :               slant_ranges_m,
     :               velocity,
     :               reflect,
     :               az_array,
     :               v_nyquist_tilt,
     :               istatus)
c
c     subroutine read_data_xm
c
c       steve albers    15-mar-1988
c       sa                 dec-1988   modified more for rt purposes
c       sa                     1994   wsr-88d remapper
c       keith brewster     jun-1994   clean-up
c
      implicit none
c
c     include files
c
      include         'remap_constants.dat'
      include         'remap_dims.inc'
      include         'remap_buffer.cmn'
c
c     input variables
c
      integer i_tilt
      real vel_thr_rtau
      real r_missing_data
c
c     output variables
c
      real gate_spacing_m_ret
      integer num_sweeps
      real elevation_deg
      integer n_rays
      integer n_gates
      real slant_ranges_m(max_gates)
      real velocity(max_gates,max_ray_tilt)
      real reflect(max_gates,max_ray_tilt)
      real az_array(max_ray_tilt)
      real v_nyquist_tilt
      integer istatus
c
c     misc local variables
c
      integer i,i_scan_mode,i_ray
      real v_nyquist_ray(max_ray_tilt)
c
c     get housekeeping info for radar scan.
c
      call get_scandata(
     :          i_tilt,
     :          r_missing_data,
     :          max_ray_tilt,
     :          n_rays,
     :          n_gates,
     :          gate_spacing_m_ret,
     :          elevation_deg,
     :          i_scan_mode,
     :          v_nyquist_tilt,
     :          v_nyquist_ray,
     :          istatus)
      if ( istatus .ne. 1 ) go to 998
c
      if (gate_spacing_m_ret .ne. gate_spacing_m) then
        write(6,805) gate_spacing_m_ret,gate_spacing_m
  805   format(' error, returned gate spacing is different from',/,
     :         ' parameter value',2f12.2)
        istatus = 0
        go to 998
      end if

      call get_azimuths_deg ( i_tilt,
     :                  max_ray_tilt,
     :                  az_array,
     :                  istatus )
      if( istatus .ne. 1 )  go to 998

      do 50 i = 1,n_gates
        slant_ranges_m(i) = float(i) * gate_spacing_m_ret
   50 continue

      call lgate_lut_gen(gsp_ref_m_cmn,gsp_vel_m_cmn
     1                  ,n_ref_gates_cmn,n_vel_gates_cmn) 

      write(6,815)
  815 format(' read_data_88d > loading arrays')

      do 100 i_ray = 1,n_rays
        call ld_ray(
     :              i_ray,
     :              n_gates,
     :              r_missing_data,
!    :              i_missing_data,
     :              velocity(1,i_ray),
     :              reflect(1,i_ray),
     :              istatus)
        if( istatus .ne. 1 )  go to 998
  100 continue

      if(namelist_parms%l_line_ref_qc)then
          write(6,*)' read_data_88d: calling rayqckz...'
!         iscale = 2
!         miss = nint(b_missing_data)
!         call rayqckz(max_ref_gates,n_rays,iscale
!    1                ,miss,reflect,az_array)
      endif

!     this has been disabled since it apparently uses up the memory when
!     we run on the ibm
      if(.false.)then
!     if(namelist_parms%l_unfold)then
          write(6,*)' read_data_88d: v_nyquist_tilt = ',v_nyquist_tilt       
          if(v_nyquist_tilt .ne. r_missing_data .and. 
     1       v_nyquist_tilt .ne. 0.                   )then
              write(6,*)' read_data_88d: calling unfold...'       
!             call unfold(max_vel_gates,n_rays,velocity
!    1                   ,az_array,v_nyquist_tilt,r_missing_data)       
              v_nyquist_ray = r_missing_data ! prevent further dealiasing
              v_nyquist_tilt = r_missing_data ! prevent further dealiasing
          else
              write(6,*)' nyquist is missing - skipping unfold call'
          endif
      endif

      return

  998 write(6,825)
  825 format(' bad status in readdata')
      return
      end
