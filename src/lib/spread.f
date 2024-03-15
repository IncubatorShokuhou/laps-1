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

      subroutine spread(a,imax,jmax,kmax,i,j,k,b)
c this subroutine currently inserts the observation into the
c three dimensional array at one point location.
      dimension a(imax,jmax,kmax)

      if(.true.)return

      a(i,j,k) = b
c     if(b .ne. 1.00)write(6,101)i,j,k,b
101   format(' spread ',3i3,f8.2)

      return
      end

      subroutine spread2(a_array,b_array,i_array,j_array,n,max_array,kma
     1x
     1                          ,i,j,k,a,b)

c     this subroutine inserts the cloud sounding into the analysis arrays
c     at one point location.

      real a_array(max_array,kmax)
      real b_array(max_array,kmax)
      integer i_array(max_array)
      integer j_array(max_array)

      a_array(n,k) = a
      b_array(n,k) = b
      i_array(n) = i
      j_array(n) = j

      return
      end
