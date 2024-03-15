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
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine  int_layerpw(x,data,kstart,qs,ps,p_1d,
     1     p1,p2,p3,lpw1,lpw2,lpw3,kk,mdf)

c     subroutine to integrate the layer precipitable water from the
c     specific humidity data. the colum is weighted by function x
c     prior to integration for purposes of meshing with func_o.f
      
c     birkenheuer  01/09/2001
      
      implicit none
      
c     parameter variables
      
      integer kk    
      real data(kk)
      real x(3)
      integer kstart
      real qs
      real ps
      real p_1d (kk) ! laps pressures  
      real lpw1,lpw2,lpw3
      real mdf
      real p1,p2,p3 ! layer tops
      
      
c     variables requireing dynamic allocation
      
c     none
      
c     internal variables
      
      integer i,j,k
      real tem_data(kk)         ! data modified by multiplier x(i)


c     generate tem_data
      do k = 1,kk
         if(data(k).ne.mdf .or. data(k).ge.0.0)then
            if(p_1d(k).ge.700. )  tem_data(k) = data(k)*x(1)
            if(p_1d(k).ge.500. .and. p_1d(k) .lt.700. )
     1           tem_data(k) = data(k)*x(2)
            if(p_1d(k).ge.100. .and. p_1d(k) .lt.500. ) 
     1           tem_data(k) = data(k)*x(3)
         else
            tem_data(k) = mdf
         endif
      enddo



c     integrate q for the pressure layers from tem data
      lpw1 = 0.0
      lpw2 = 0.0
      lpw3 = 0.0

      do k = kstart,kk-1
         if (p_1d(k) .ge. p1 .and. lpw1 .ne. mdf )then
            lpw1 = lpw1 +(tem_data(k)+tem_data(k+1))/2. *
     1           (p_1d(k)-p_1d(k+1))
            if(tem_data(k) .eq. mdf) lpw1 = mdf
         elseif (p_1d(k) .ge. p2 .and. lpw2 .ne. mdf )then
            lpw2 = lpw2 +(tem_data(k)+tem_data(k+1))/2. *
     1           (p_1d(k)-p_1d(k+1))
            if(tem_data(k) .eq. mdf ) lpw2 = mdf
         elseif (p_1d(k) .ge. p3 .and. lpw3 .ne. mdf )then
            lpw3 = lpw3 +(tem_data(k)+tem_data(k+1))/2. *
     1           (p_1d(k)-p_1d(k+1))
            if(tem_data(k) .eq. mdf) lpw3 = mdf
         endif
      enddo                     !k



c     change units of integrated(q) to  g/kg     
      if (lpw1 .ne. mdf) lpw1 = lpw1 *1.e3 
      if (lpw2 .ne. mdf) lpw2 = lpw2 *1.e3 
      if (lpw3 .ne. mdf) lpw3 = lpw3 *1.e3 
c     add surface layer qs (already in g/kg) to first level
c      if (lpw1 .ne. mdf) lpw1 = lpw1 
c     1     + ( qs + data(kstart) *1.e3 ) /2.
c     1     * (ps-p_1d(kstart) )
c     comvert g/kg to mm <-- millimeters, not cm!
      if (lpw1 .ne. mdf) lpw1 = lpw1 / 10. / 9.8
      if (lpw2 .ne. mdf) lpw2 = lpw2 / 10. / 9.8
      if (lpw3 .ne. mdf) lpw3 = lpw3 / 10. / 9.8

c     check bounds

      if (p1.gt.ps)  lpw1 = mdf
      if (p2.gt.ps)  lpw2 = mdf
      if (p3.gt.ps)  lpw3 = mdf

c     check zeros

      if (lpw1.eq. 0.0)  lpw1 = mdf
      if (lpw2.eq. 0.0)  lpw2 = mdf
      if (lpw3.eq. 0.0)  lpw3 = mdf
   

      
      return
      
      end
