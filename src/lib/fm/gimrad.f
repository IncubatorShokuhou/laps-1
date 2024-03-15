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

      function gimrad(tau,temp,tskin,kchan,lsfc,psfc,emiss)
c  12 dec 94  correct error in sfc tau extrapolation  
c $ function gimrad(tau,temp,tskin,kchan,lsfc,psfc,emiss)   (hmw)
c $ gimrad - compute goes radiance from temperature profile
c $ input:
c $     tau = (r) array of goes transmittances
c $     temp = (r) array of atmospheric temperatures
c $     tskin = (r) temperature of surface
c $     kchan = (i) channel number
c $     lsfc = (i) subscript of level below or at surface
c $     psfc = (i) pressure of surface ( or cloud)
c $     emiss = (r) surface emissivity
c $ output description:
c $     goes radiance
c $$ gimrad = sounder, goes, convert
      dimension tau(*),temp(*)
      common/atmos/p(40),t(40),w(40),ozo(40)
      common /atmrad/atrd,taus,refl
c     prepare surface correction
      dp = alog(psfc/p(lsfc))/alog(p(lsfc)/p(lsfc-1))
c     dp is negative to interpolate, positive to exrapolate (psfc>nl)
      dtau = tau(lsfc-1)-tau(lsfc)
      taus = tau(lsfc)-dtau*dp
      refl = 0.
      t1 = temp(1)
      b1 = plango(t1,kchan)
      tau1 = tau(1)
      rad = 0.
c     integrate down to level at or below surface
      do 110 i = 2,lsfc
      t2 = temp(i)
      b2 = plango(t2,kchan)
      tau2 = tau(i)
      dtau = tau1-tau2
      dr =.5*(b1+b2)*dtau
      if (taus .gt. 0.1.and.emiss .lt. 1.00)then
c     do not add reflected for last level unless psfc .gt. 1000
      if (i .eq. lsfc.and.psfc .le. 1000.)go to 105
      taub = 0.5*(tau1+tau2)
      taufac = taus/taub
      refl = refl+taufac*taufac*dr
      endif
  105 rad = rad+dr
      b1 = b2
  110 tau1 = tau2
c     add (subtract) increment of atmos radiance to reach surface
c     dp will be negative if psfc  < 1000 mb
c     dr falls out as the delta radiance of layer
      rad = rad+dr*dp
c     add increment of reflected radiance for layer down to surface
      if (taus .gt. 0.1.and.emiss .lt. 1.00)then
         if (psfc .lt. 1000.)then
         taub = 0.5*(tau(lsfc-1)+taus)
c     change dp to increment rather than decrement
         dp = 1.+dp
         else
         taub = 0.5*(tau(lsfc)+taus)
         endif
      taufac = taus/taub
      refl = refl+taufac*taufac*dr*dp
      endif
      atrd = rad
      rad = rad+(1.-emiss)*refl
      bs = plango(tskin,kchan)
      rad = rad+emiss*bs*taus
      rad = amax1(rad,.001)
      gimrad = rad
      return
      end
