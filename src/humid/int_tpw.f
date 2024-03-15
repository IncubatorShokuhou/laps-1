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
      subroutine  int_tpw(data,kstart,qs,ps,p_3d,tpw,mdf,ii,jj,kk)

c     subroutine to integrate the total precipitable water from the
c     specific humidity data.  consolodating this process, since it is
c     done numerous times in the code
      
c     birkenheuer  feb 9 1993
      
      implicit none
      
c     include 'lapsparms.for'
c     include 'parmtrs.inc'
      
c     parameter variables
      
      integer ::  ii,jj,kk
      real, dimension(ii,jj,kk) :: data,p_3d
      integer kstart (ii,jj)
      real, dimension(ii,jj) :: qs, ps, tpw
      real :: mdf               !laps missing data flag
      real, dimension (500) :: p_vert,data_vert
      real, dimension(3) :: xdum = (/1.,1.,1./) !dummy variable for call
      
      
c     variables requiring dynamic allocation
      
c     none
      
c     internal variables
      
      integer i,j,k
      
c     integrate the tpw field
      
      do j = 1,jj
         do i = 1,ii
            
            tpw(i,j) = 0.0

c     fill vertical specific arrays
            do k = 1,kk
               p_vert(k) = p_3d(i,j,k)
               data_vert(k) = data(i,j,k)
            enddo               !end k

            call int_ipw (xdum,p_vert,data_vert,kstart(i,j),
     1           qs(i,j),ps(i,j),tpw(i,j),mdf,kk)

         enddo                  !i
      enddo                     !j

      return
      
      end
