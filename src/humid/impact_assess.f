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
      subroutine impact_assess (data_back, data_final, tpw1,tpw2, 
     1     ii,jj,kk, 
     1     gps_data,gps_w, p, mdf)

c     this routine compares the original and final analysis gridpoints to the
c     gps value for tpw at gps sites.  

c     the purpose of this routine is to assess the change in tpw over gps areas
c     for an assessment of improvement to the analysis.

c     parameter variables

      implicit none

      integer :: ii,jj,kk       !passed in dimensions
      real, dimension (ii,jj,kk) :: data_back, data_final
      real :: p                 !location pressure
      real, dimension (ii,jj) :: gps_data, gps_w,tpw1,tpw2
      real :: mdf               !laps missing data flag

      integer :: i,j,k

c     compute the ipw data for each place where we have gps data

      do i = 1,ii
         do j = 1,jj

c     call to ipw routine at gps point

            if(gps_w(i,j) == 1.0) then
               write(6,*) 
     1              'tassess',tpw1(i,j),tpw2(i,j),gps_data(i,j),i,j
            endif

c     ipw_data_back = f(i,j)
c     call int_ipw (x,p,data_back(i,j),ipw_data_back(i,j),mdf,kk)
c     call int_ipw (x,p,dadta_final(i,j),ipw_data_final(i,j),mdf,kk)
c     ipw_data_final = f(i,j)

         enddo                  !j
      enddo                     !i

c     print out comparisons of co-located points only
 

c     do n = 1,nn
c     write(6,*) ipw_data_back(point(2,n),point(3,n), ipw.....

      return

      end subroutine impact_assess
 
