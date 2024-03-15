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

      subroutine laps_slwc_revb(cb_pa,cb_k,gt_pa,gt_k,ct_k,
     1           adiabatic_lwc,adjusted_lwc,adjusted_slwc,
     2           i_status1,i_status2)
c
c.......................history.............................
c
c        written: ca. 1982 by w. a. cooper in hp fortran 4
c
c....... calculates temperature t and liquid water content from
c..      cloud base pressure p0 and temperature t0, for adiabatic
c..      ascent to the pressure p.
c..     ->  input:  cloud base pressure p0 and temperature t0
c..                 pressure at observation level p
c..     ->  output: "adiabatic" temperature t and liquid water content
c
c        modified: november 1989 by paul lawson for laps/wisp.  routine
c                  now calculates adiabatic liquid water content
c                  (adiabatic_lwc) using cloud base pressure and grid-top
c                  temperature and pressure.  also calculated are adjusted_lwc,
c                  which adjusts adiabatic_lwc using an empirical cloud
c                  water depletion algorithm, and adjusted_slwc, which is
c                  adiabatic_lwc in regions where t < 0 c adjusted
c                  using an empirical algorithm by marcia politovich.
c
c                  subroutine is now hardwired for stratiform cloud only.
c                  can be modified to include cu with input from laps main.
c
c                  revb: ca 12/89 calculate adiabatic lwc by going from cloud
c                        base to laps grid level instead to cloud top, thus
c                        helping to better calculate in layer clouds.
c                        add tg (grid temperature) to calcualtion.
c
c                  revc: 2/27/90 correct error in code.  zero-out slwc when grid
c                        temperature (gt) > 0.
c
c                  revd: 1/27/99 correct apparent error in test for i_status1
c                        (steve albers)
c
c
c        outputs:  adiabatic_lwc
c                  adjusted_lwc
c                  adjusted_slwc
c                  i_status1 - 1 when -20 < cld_top_temp < 0 for stratus
c                              0 otherwise
c                  i_status2 - 1 when valid input data provided from main
c
      data eps/0.622/,cpd/1.0042e3/,cw/4.218e3/,rd/287.05/,alhv/2.501e6/
      integer cty
      integer i_status1, i_status2
      i_status1=1
      i_status2=1
c   2 print *,'enter: p-base(mb), t-base(c), p-top, t-top, cld type'
c     read(5,*) p0, t0, p, ctt, cty
c     if(cty.ne.0.and.cty.ne.1) go to 2
c
c     hardwire cloud type (cty) for stratus for now
c
      cty=0
c
c.....convert pa to mb and kelvin to celcius
c
      p0 = cb_pa/100.
      p  = gt_pa/100.
      t0 = cb_k - 273.15
      tg = gt_k - 273.15
      ctt= ct_k - 273.15
c     print *, 'ctt in sub = ', ctt
c
c     check for valid input data...
c
        if(p0.gt.1013..or.p0.lt.50.) then
          i_status2=0
          return
        else
        endif
c
c
        if(t0.gt.50..or.t0.lt.-70.) then
          i_status2=0
          return
        else
        endif
c
c
        if(p.gt.1013..or.p.lt.50.) then
          i_status2=0
          return
        else
        endif
c
c     set i_status1 = f if 0 < cld top < -20 c (for stratus).
c
      if(ctt.ge.0..or.ctt.le.-20.) i_status1=0
c
      tk=t0+273.15
      e=vapor(t0)
      r=eps*e/(p0-e)
      cpt=cpd+r*cw
      thetaq=tk*(1000./(p0-e))**(rd/cpt)*exp(alhv*r/(cpt*tk))
c 1st approx
      t1=tk
      e=vapor(t1-273.15)
      rv=eps*e/(p-e)
      t1=thetaq/((1000./(p-e))**(rd/cpt)*exp(alhv*rv/(cpt*t1)))
c successive approximations
      do 1 i=1,10
      e=vapor(t1-273.15)
      rv=eps*e/(p-e)
      t1=(thetaq/((1000./(p-e))**(rd/cpt)*exp(alhv*rv/(cpt*t1)))
     $   +t1)/2.
      t=t1-273.15
c     print *, p0,t0,p,t,e,rv,thetaq
1     continue
c get lwc
      e=vapor(t)
      rv=eps*e/(p-e)
      tw=r-rv
      adiabatic_lwc=tw*p*28.9644/(8.314e7*t1)*1.e9
      if(adiabatic_lwc.lt.0.) adiabatic_lwc=0.
c     print *, 'adiabtic lwc = ', adiabatic_lwc
      if(tg.ge.0.) then
c
      adjusted_slwc=0.                                          ! added 2/27/90
c

         if(cty.eq.0.) then
           if(ctt.lt.-20.) then
             adjusted_lwc=0.
           elseif(ctt.lt.-15..and.ctt.ge.-20.) then
             adjusted_lwc=adiabatic_lwc/8.
           elseif(ctt.lt.-10..and.ctt.ge.-15.) then
             adjusted_lwc=adiabatic_lwc/4.
           else
             adjusted_lwc=adiabatic_lwc/2.
           endif
         else
           if(ctt.lt.-25.) then
             adjusted_lwc=0.
           elseif(ctt.lt.-15..and.ctt.ge.-25.) then
             adjusted_lwc=adiabatic_lwc/8.
           elseif(ctt.lt.-10..and.ctt.ge.-15.) then
             adjusted_lwc=adiabatic_lwc/4.
           else
             adjusted_lwc=adiabatic_lwc/2.
           endif
         endif
      else
         if(cty.eq.0.) then
           if(ctt.lt.-20.) then
             adjusted_lwc=0.
             adjusted_slwc=0.
           elseif(ctt.lt.-15..and.ctt.ge.-20.) then
             adjusted_lwc=adiabatic_lwc/8.
             adjusted_slwc=adiabatic_lwc/8.
           elseif(ctt.lt.-10..and.ctt.ge.-15.) then
             adjusted_lwc=adiabatic_lwc/4.
             adjusted_slwc=adiabatic_lwc/4.
           else
             adjusted_lwc=adiabatic_lwc/2.
             adjusted_slwc=adiabatic_lwc/2.
           endif
         else
           if(ctt.lt.-25.) then
             adjusted_lwc=0.
             adjusted_slwc=0.
           elseif(ctt.lt.-15..and.ctt.ge.-25.) then
             adjusted_lwc=adiabatic_lwc/8.
             adjusted_slwc=adiabatic_lwc/8.
           elseif(ctt.lt.-10..and.ctt.ge.-15.) then
             adjusted_lwc=adiabatic_lwc/4.
             adjusted_slwc=adiabatic_lwc/4.
           else
             adjusted_lwc=adiabatic_lwc/2.
             adjusted_slwc=adiabatic_lwc/2.
           endif
         endif
      endif
c     print *,'adjusted lwc = ', adjusted_lwc
c     print *,'adjusted slwc = ', adjusted_slwc
      end
c  function to calculate vapor pressure:
c
      function vapor(tfp)
c input is in degrees c.  if gt 0, assumed to be dew point.  if
c less than 0, assumed to be frost point.
c routine codes goff-gratch formula
      tvap=273.16+tfp
      if(tfp.gt.0.) go to 1
c this is ice saturation vapor pressure
      if(tvap.le.0) tvap=1e-20
      e=-9.09718*(273.16/tvap-1.)-3.56654*alog10(273.16/tvap)
     $  +0.876793*(1.-tvap/273.16)
      vapor=6.1071*10.**e
      return
 1    continue
c this is water saturation vapor pressure
      if(tvap.le.0) tvap=1e-20
      e=-7.90298*(373.16/tvap-1.)+5.02808*alog10(373.16/tvap)
     $  -1.3816e-7*(10.**(11.344*(1.-tvap/373.16))-1.)
     $  +8.1328e-3*(10.**(3.49149*(1-373.16/tvap))-1)
      vapor=1013.246*10.**e
      return
      end

