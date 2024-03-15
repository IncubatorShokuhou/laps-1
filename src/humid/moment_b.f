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
cdis cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine moment_b(data,n,ave,adev,sdev,var,skew,curt,istatus)

c alternative moment program adapted from numerical recipies not to comflict
c with the other routine modifed and installed by jim edwards (moment.f)
c found under the satellite library area.


      integer n
      real adev,ave,curt,sdev,skew,var,data(n)
      integer j, istatus
      real p,s,ep
      istatus = 0 ! bad
      if(n.le.1) then
         write (6,*) 'must be at least 2 in moment comp'
         return
      endif
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      if(var.lt.0.) then
         write(6,*) 'negative var detected, taking absolute value'
         var = abs(var)
      endif
      sdev=sqrt(var)
      if(var.gt.1.e-10) then               
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      else
        skew = 0.0
        curt = 0.0
      write (6,*) 'skew and curt returned as zero, variance is small'
        return
      endif
      istatus = 1 !good
      return
      end
