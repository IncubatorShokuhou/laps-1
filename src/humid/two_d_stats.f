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
      subroutine two_d_stats (ii,jj,x,excluded_value)

      implicit none
      integer ii,jj
      real x(ii,jj)
      real excluded_value


c     internal variables

      integer i,j
      real maxa,mina
      integer maxi,maxj, mini,minj
      integer counter
      real x_linear (ii*jj)
      integer istatus
      real ave,adev,sdev,var,skew,curt

c     compute the max and min values of the field

      maxa = -1.e32
      mina = 1.e32

      do i = 1,ii
         do j = 1,jj
            if(x(i,j) .ne. excluded_value) then
               mina = min (mina, x(i,j))
               maxa = max (maxa, x(i,j))
               if(mina.eq.x(i,j)) then
                  mini = i
                  minj = j
               endif
               if(maxa.eq.x(i,j)) then
                  maxi = i
                  maxj = j
               endif
            endif
         enddo
      enddo

      write (6,*)
      write (6,*) 'computed 2-d statistics'
      write (6,*) 'max and min values in the field'
      write (6,*) 'max,i,j', maxa,maxi,maxj
      write (6,*) 'min,i,j', mina,mini,minj
      write (6,*)

c     begin computation of moment statistics.

      counter = 0

      do i = 1,ii
         do j = 1,jj
            if (x(i,j) .ne. excluded_value) then
               counter = counter + 1
               x_linear(counter) = x(i,j)
            endif
         enddo
      enddo

      call moment_b (x,ii*jj,ave,adev,sdev,var,skew,curt,istatus)

      write (6,*) 'moment output'

      write (6,*) 'field average ',ave, '+/-',sdev

      return
      end



      
