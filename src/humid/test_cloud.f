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
c fortran 90
      subroutine test_cloud (t1,t2,q1,q2,p,qsat, qadjust)

      implicit none

      real :: tbar, ubar, qbar, delu, delq, delt
      real :: terror, qerror
      real :: t1,t2, q1,q2, p, tref = -132.0
      real :: qsat
      real :: qadjust,uadj
      real :: temp_uncertainty = 0.005 ! 0.5% precision error in temp
      real :: q_uncertainty = 0.01 ! 1% precision error in sh

      real, external :: ssh2    !function for saturation specific 
                                ! humidity

c     tref = temperature for ice phase reference
c     p = pressure in mb (something close is all that is needed)
c     level 1 is low, level 2 is higher, level 1 is (k)
c     bar values are mean for the range
c     del values are allowed deltas
c     u = relative humidity (fraction)
c     t = temperature degk
c     terror, qerror are the approximate error in the analyzed values of
                                !these components, the natural error
                                !in these variables (precision) will 
                                !have a direct bearing on the probabily 
                                !of cloud.   
c     q = specific humidity 
c     qbarsat = the computed saturated q based on tbar
c     qadjust = the min amount that qbar can be increased to make a cloud
c     routine based on climo scheme developed by peixoto and oort 1997
c     refer to noaa tech note for reference

      tbar = (t1+t2)/2.         ! average temperature in k
      qbar = (q1+q2)/2.         ! average specific humidity

      terror = temp_uncertainty* tbar !approximate temperature error
      
      qerror = q_uncertainty* qbar !approximate sh error
      
      delq = abs(q2-q1)         !box moisture variabiliy
      delt = abs(t2-t1)         !box  temperature variabilty
      
      ubar = q1/ qsat           !use gridpoint value here 
      
      delu = ubar * ( (delq+qerror)/qbar + 20.*(delt+terror)/tbar )
      
c     this is the accepted variation of relhumidity in the "box"
c     immediately adjacent to a a level  we now apply this variation
c     to that level (k) being evaluated.
     
c     if this range extends rh to saturation, then the box is
c     assumed to support cloud.  if clouds are associated with this
c     box, then no adjustment to moisture is needed if the cloud
c     analysis spcifies clouds

c     test the box for cloudy conditions.

      if (ubar+delu/2. >= 1.) then ! clouds are seen as a possibility
         qadjust = 0.0  ! no need to enhance q if clouds are present
      else                      !clouds are not a possiblity compute
                                !the minimum amount of extra q needed 
                                ! to make a cloud somewhere in the box

         uadj = 1. - (ubar+delu/2.)

         qadjust = uadj * qsat

         if (qadjust < 0 ) qadjust = 0.0 !slight roundoff correction

c     this is the approximate increment in q needed to achieve saturation

      endif

      return

      end subroutine test_cloud
