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

      function eslo(t)

      implicit none

      real eslo
      real a0,a1,a2,a3,a4,a5,a6
      real t
      real esice ! function

c     baker, schlatter  17-may-1982	  original version.
c     added functionality to extend to temps below -47c (by approx)
c     d. birkenheuer 5/19/99.  approx based on formula via pielke's book
c     tuned by dan to match lowe's function 0>t>-47c.

c     this function returns the saturation vapor pressure over liquid
c     water eslo (millibars) given the temperature t (celsius). the
c     formula is due to lowe, paul r.,1977: an approximating polynomial
c     for the computation of saturation vapor pressure, journal of applied
c     meteorology, vol 16, no. 1 (january), pp. 100-103.
c     the polynomial coefficients are a0 through a6.

      data a0,a1,a2,a3,a4,a5,a6
     1     /6.107799961,     4.436518521e-01, 1.428945805e-02,
     2     2.650648471e-04, 3.031240396e-06, 2.034080948e-08,
     3     6.136820929e-11/


      if (t.gt.-46.) then ! lowe's function
         eslo = a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t)))))
         if (eslo.lt.0.) eslo = 0.
      elseif (t.gt.-132.) then ! pielke's tuned function (liquid data unknown)
         eslo = 6.10813510 * exp (t*17.66795997 / 
     1        (t+273.15-29.96392251) )
      else ! adjust to ice, liquid data unknown
         eslo = esice(t)
      endif
      return
      end
