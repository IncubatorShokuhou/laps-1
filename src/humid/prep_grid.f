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


      subroutine prep_grid(m,n,data,pn_max,points,pn,istatus)

c     $log: prep_grid.for,v $
c     revision 1.1  1996/08/30  20:48:53  birk
c     initial revision
c

      implicit none

      integer m,n,pn,pn_max,istatus

      real points(pn_max,3),data(m,n),weight_t,dist,weight


      integer i,j,k

      istatus = 0

c     i perimeter set

      do i = 1,m,m-1
         do j = 1,n


            weight_t = 0.0
            data(i,j) = 0.0


            do k = 1,pn

               if (points(k,2).eq.i  .and. points(k,3) .eq.j) then
                  data(i,j) = points(k,1)
                  go to 22
               elseif (points(k,2) .ne. 0) then
                  
                  dist = sqrt( (i-points(k,2))**2+(j-points(k,3))**2)
                  if (dist > 1.e-9) then
                     weight = 1./dist
                     data(i,j) = data(i,j) + weight*points(k,1)
                     weight_t = weight_t + weight
                  else
                     weight = 1.e9
                     data(i,j) = data(i,j) + weight*points(k,1)
                     weight_t = weight_t + weight
                  endif
               endif
               
            enddo               !k
            
            if(weight_t > 1.e-9) then
               data(i,j) = data(i,j) / weight_t
            else
               continue
            endif
            
 22         continue
            
         enddo                  !j
      enddo                     !k
      
      
      if(weight_t.eq.0.) then 
         istatus = 0
         write(6,*) 'returning from prep_grid.f on istatus 0'
         return
      endif

c     j perimeter set
      

      do i = 1,m
         do j = 1,n,n-1
            
            
            weight_t = 0.0
            data(i,j) = 0.0
            
            
            do k = 1,pn
               
               if (points(k,2).eq.i  .and. points(k,3) .eq.j) then
                  data(i,j) = points(k,1)
                  go to 23
               elseif (points(k,2).ne.0) then
                  
                  dist = sqrt( (i-points(k,2))**2+(j-points(k,3))**2)
                  if (dist > 1.e-9) then
                     weight = 1./dist
                     data(i,j) = data(i,j) + weight*points(k,1)
                     weight_t = weight_t + weight
                  else
                     weight = 1.e9
                     data(i,j) = data(i,j) + weight*points(k,1)
                     weight_t = weight_t + weight
                  endif
               endif
               
            enddo               !k
            
            if(weight_t > 1.e-9) then
               data(i,j) = data(i,j) / weight_t
            else
               continue
            endif
            
 23         continue
              
         enddo                  !j
      enddo                     !k
      
      
      if(weight_t.eq.0.) then 
         istatus = 0
         write(6,*) 'returning from prep_grid.f on istatus 0'
         return
      endif

      istatus = 1

      return
      end
      





