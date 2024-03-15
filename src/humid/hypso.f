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
c
c
c
      subroutine hypso (
     1     z,                   !height meters n levels     
     1     t,                   !temp kelvin, nlevels
     1     rmd,                 !missing data flag for real numbers
     1     nn,                  !levels nn
     1     p,                   !pressure hpa, n, levels
     1     abort                !abort flag, 0= abort, 1= good
     1     )

c     notes:
c
c     this routine assumes dry air but could use virtual temperature
c     as input if so desired.
c
c     even though this routine is tailored for radiometer use, 
c     the 2km cutoff in height (according to stick ware) is not
c     performed here but in the routine above this one.
c
c     author: dan birkenheuer 12/14/2010
c     purpose: add radiometrics (brand) of radiometer data to laps
c


      implicit none

c     parameters for radiometer hypsometric eqn to determine pressure

      integer nn                !levels passed
      real z(nn)                !height (m)
      real t(nn)                !temperature (c)
      real rmd                  !missing data flag
      real p(nn)                !pressure (hpa) (variable to be "descovered"
      integer abort             !0=fail, 1= success

c     internal variables and constants

      real r                    !gas constant
      real g                    !gravity
      real c2k                  !conversion to kelvin
      real p2p                  !converson pa to hpa
      real delz                 !delta z (m)
      real tbar                 !average temperature (k)
      integer n                 !index 

c     code

c     assign constants
      r = 287.04                !gas constant for dry air
      g = 9.80665
      c2k = 273.15
      p2p = 0.01 

c     run hypsometric equation to determine p from sfc values an hts

      if ( (nn .eq. 1) .or. (z(1) .eq. rmd) ) then ! bad p value abt
         abort = 0 
         return
      endif


      do n = 1, nn-1            !stop shy of top level

         delz = z(n+1) - z(n)   !z2-z1 (m)
         tbar = (t(n+1)+t(n)+2*c2k)/2.0 ! (k)
         p(n+1) = p(n)*p2p/(exp (g/(r*tbar)*(delz))) ! in pa
         p(n+1) = p(n+1)/p2p    !in hpa (as returned desire)


c     this is coded for clarity.  improvement would be to remove the
c     p2p operationa in mult and later divide to improve speed.  some
c     complilers might do this.

      enddo

      abort = 1

      return
      end


