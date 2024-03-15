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
      subroutine ofm (kk,laps_p,laps_t,laps_q,laps_sfc_t,
     1     psfc, jday, lat, za,
     1                                tbest, 
     1     radest,
     1     sec_za,              !optran 90 variable 
     1     sfc_emis,            !optran 90 variable
     1     sfc_refl,            !optran 90 variable
     1     sec_solar            !optran 90 variable
     1     )


c   this routine interfaces goes 8/10 satellite broadcast network data (and
c   local gvar data) to the laps moisture analysis.  in 1999, this routine
c   was modified from an earlier version that used the university of
c   wisconsin -- madison's forward model to a newer model developed at
c   nesdis.  optran (optical transmittance) forward model was developed by
c   thomas kleespies (nesdis) and questions about this model should be
c   directed to him.  forecast systems laboratory does not in any way
c   guarantee the validity of optran and distributes this software on an
c   as-is basis.  moreover, fsl has permission to distribute optran as part
c   of laps to federal agencies.  non-federal entities need to inquire with
c   nesdis to establish their rights and obligations with regard to optran.
c   
c   the version of optran with which this software is used, has been
c   modified by fsl to include both sounder and imager channels for a
c   particular satellite in one call to the routine.  thus a user only need
c   to setup optran for a particular satellite.  after doing such, either
c   the imager or sounding instrument can be used with the software without
c   further recompilation.  
      
 

c      use module_sfc_structure

      implicit none

      save
      
      include '../include/constants_optran.inc'
      
      include '../include/trans.inc'
      
      integer mchan
      parameter (mchan=18) ! adjusted for optran 90.
      
      integer kk                ! kk the laps profile dim, 
      integer nk                ! the size of the composite vectors
      integer start_level       !lowest level of climo used
      real laps_p(kk),laps_t(kk),laps_q(kk), lat, laps_sfc_t, psfc
      real za                   !zenith angle (degrees)
      real sec_za               !secant of za (for optran 90)
      real sfc_emis             !surface emissivity for optran 90
      real sfc_refl             !surface reflectivity for optran 90
      real sec_solar            !local secant of solar angle for optran 90
      real, dimension (:,:), allocatable :: tau90,flux_tau,
     1     solar_tau            !tau for optran 90 return variable
      real, dimension (:), allocatable ::  upwelling_radiance, 
     1     brightness_temperature !upwelling radiance 
      integer compute_rtm      ! optran 90 wrapper function is type int
      integer jday, istatus
      logical first_call
      data first_call /.true./

c     optran 90 variables

      real, dimension (0:nlevel,1)::  level_p,level_t, level_w, level_o
      real, dimension (nlevel,1)  ::  layer_p,layer_t, layer_w,layer_o
      integer nk_90 ! optran 90 reduced nk counter
      
c     climo variables
      real c_p (40)
      real c_t (40)
      real c_m (40)
      
      real conversion
      real mdf
      
c     satellite number

      integer goes_number
      common /sat_id/ goes_number
      
      real o(nlevel)
      
      real tbest(nchan)         ! brightness temps from estimated tau
      real radest(nchan)
      
      integer ichan
      
      integer channels(mchan) 
      data channels/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/
      
      integer i

c     optran 90 coefficient data information

      character*256 sndr_coeff, sndr_trans
      integer sndr_coeff_len, sndr_trans_len
      common /optn_coef/ sndr_coeff, sndr_trans, sndr_coeff_len,
     1     sndr_trans_len

      
c     ----------------------do ozone first call ----

      if (first_call) then
         first_call = .false.
         
c     determine climo start level
         
         call climate_sm(lat, jday, c_p, c_t, c_m, istatus)
         
         start_level = 40       ! (high c_p)
         do while (  c_p(start_level) .ge. laps_p(kk) )
            start_level = start_level - 1
         enddo
         
c     fill pressure level
         
         nk = 0
         
         do i = 1,start_level
            nk = nk +1
            p(nk)   = c_p(i)
         enddo
         
         do i = kk,1,-1
            nk      = nk+1
            p(nk)   = laps_p(i)
         enddo
         
         write(6,*) 'ozone levels used = ',nk
         
c     grab ozone
         
         call oh_zone (p,o,nk,1,istatus)
         
c     convert (g/g) to (g/kg)
         
         conversion    = 1.0e+3
         
c     ozone called only one time
         
      endif
      
c------------------------------------------------
      
c     00000000000000000000000000000
c     assemble vectors module 1
c000000000000000000000000000000
      
c     build vectors, count total levels (not including surface)
      
c     determine climo start level
      
      start_level = 40          ! (high c_p)
      do while (  c_p(start_level) .ge. laps_p(kk) )
         start_level = start_level - 1
      enddo

      nk = 0
      
      do i = 1,start_level
         nk = nk +1
         p(nk)   = c_p(i)
         t(nk,1) = c_t(i)
         q(nk,1) = c_m(i)       ! note this is g/kg mixing ratio
      enddo
      
      do i = kk,1, -1
         if(laps_p(i) .lt. psfc) then ! accept
            nk      = nk+1
            p(nk)   = laps_p(i)
            t(nk,1) = laps_t(i)
            q(nk,1) = laps_q(i) * conversion ! note this is g/kg mixing ratio
         endif
      enddo
      
c     now at sfc or cloud top
      
      nk = nk + 1
      p(nk) = psfc
      t(nk,1) = laps_sfc_t
      q(nk,1) = q(nk-1,1)       ! approximation for now
      
      if(nk.gt.nlevel) then
         write(6,*) 'array dimension error'
         write(6,*) 'module ofm.f'
         write(6,*) 'parameter nlevel too small'
      endif
      
c0000000000000000000000
c     end assembling vectors
c00000000000000000000000

cc insert here new code for optran 90 interface  db  2002
c upper level layer computations.

      do i = 1,nk
         level_p (i-1,1) = p(i)
         level_t (i-1,1) = t(i,1)
         level_w (i-1,1) = q(i,1)
         level_o (i-1,1) = o(i)
      enddo
      nk_90 = nk - 1            ! (start counter at zero for layers)

c     generate layers

      call make_optran_layers ( level_t, level_w, level_o, 
     1     layer_t, layer_w, layer_o, level_p,layer_p,
     1     nk_90, istatus)

c surface parameters

      
c     call to optran 90



c      if (nk .eq. 10000) then

         allocate (tau90 (1:nk_90, 1:mchan))
         allocate (flux_tau (1:nk_90,1:mchan))
         allocate (solar_tau (1:nk_90,1:mchan))
         allocate (upwelling_radiance (1:mchan))
         allocate (brightness_temperature (1:mchan))
         
         tau90 = 0.0
         flux_tau = 0.0
         solar_tau = 0.0
         upwelling_radiance = 0.0
         brightness_temperature = 0.0
         
         call optran90_fm (     ! call to optran 90 forward model
     1        nk_90,level_p, layer_p, layer_t, layer_w, layer_o,
     1        laps_sfc_t,
     1        sfc_emis,
     1        sfc_refl,
     1        sec_za,
     1        sec_solar,
     1        mchan,
     1        tau90,
     1        flux_tau,
     1        solar_tau,
     1        upwelling_radiance,
     1        brightness_temperature,
     1        sndr_coeff,
     1        sndr_trans,
     1        sndr_coeff_len,
     1        sndr_trans_len
     1        )
         
         
         
         do ichan = 1 , mchan
            
            tbest (ichan) = brightness_temperature (ichan) 
            radest (ichan) = upwelling_radiance(ichan)
            
         enddo
         

        
         
         deallocate (tau90,flux_tau)
         deallocate (upwelling_radiance, brightness_temperature)
         deallocate (solar_tau)
         


      return
      
      end
