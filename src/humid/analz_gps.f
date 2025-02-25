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
      subroutine analz_gps (lat,lon,wt,nn,glat,glon,
     1                       data_out,data_weights,
     1     gps_points,idotj,gps_count,ii,jj,istatus)

      implicit none

c     wt = gps total water cm
c     nn = total number gps obs (even ones with missing data)
c     glat, glon are laps gridpoints, ii,jj domain dimensions of laps grid
c     data_out = output analyzed gps data on the laps grid
c     data_weights = weight(distance from ob) at each laps gridpoint
c     
c     internal variables of special interest
c
c     r50 = radius(meters) from ob where weight drops to 0.5 (50%)

      integer ii,jj,nn,istatus
      real lat(nn),lon(nn),wt(nn)
      real data_out(ii,jj)
      real data_weights(ii,jj)
      real glat(ii,jj)
      real glon(ii,jj)
      integer idotj
      real gps_points(idotj,3)
      integer gps_count

c     volitile array points
      real points (nn,3)
      integer mask (ii,jj)
      real r50
      real ri, rj


      integer i,j,n,ncount
      real sum

c
      istatus =  0

c     initialaize data_out on first call  (-1 used as missing)
c     since zero water is valid number

      do i = 1,ii
         do j = 1,jj
            data_out(i,j) = -1.
            mask (i,j) = 0
            data_weights(i,j) = 0.0
         enddo
      enddo



c     foreach n element of wt, determine its location in ii,jj space

      ncount = 0

      do n = 1,nn
c         write(6,*) 'gpstemp latloncheck', lat(n),lon(n)
         if (abs(lat(n)) <= 90.000) then
          if(abs(lon(n)) <= 180.0 ) then
            call  latlon_to_rlapsgrid(lat(n),lon(n),glat,glon,ii,jj,
     1        ri,rj, istatus)
          endif
         else
            istatus = 0 ! floating point error check for jet
         endif

         if (istatus == 1 .and. wt(n) >  0.0) then
            ncount = ncount + 1
            i = nint (ri)
            j = nint (rj)
            data_out (i,j) = wt(n)
            mask(i,j) = 1
            data_weights (i,j) = 1.
            points(ncount,1) = wt(n)
            points(ncount,2) = i
            points(ncount,3) = j

            write(6,*)lat(n),lon(n), glat(i,j), glon(i,j), wt(n)
         else
            continue
         endif

      enddo

c     now that data_out is as full as it is going to get with the water
c     process the remainder.

c     fill outbound variables for gps comparison stuff for seth gutman
      gps_count = ncount
      do i = 1, ncount
         gps_points(i,1) = points(i,1)
         gps_points(i,2) = points(i,2)
         gps_points(i,3) = points(i,3)
      enddo

c     compute the fraction of data_out that is empty

      sum = 0.

      do i = 1,ii
         do j = 1,jj
            if (data_out (i,j) ==  -1.) then
               sum = sum + 1
               data_out(i,j) = points(1,1) ! helps converge
            endif
         enddo
      enddo

      if(ncount.eq.0) then      !abort gps here
         write(6,*) 'no gps data avail to process'
         istatus = 0
         return
      endif

      write (6,*) ncount, ' out of ', n, ' total avial data used'
      write (6,*) 'gps field fraction empty = ', sum/ii/jj

c      if (sum/ii/jj .gt. 0.75) return ! return if not enough data

      
c     now have fairly full data array.  now analyze

      call prep_grid (ii,jj,data_out,nn,points,ncount,istatus)
      if (istatus == 0) then
         write (6,*) 'error in prep_grid'
         return
      endif
      call slv_laplc (data_out,mask,ii,jj)

c     prep the weighting array for the above analyzed sheet

      r50 = 750.e+3 ! m radius of influence (here 750 km)    

      call weight_field (data_weights, mask,  ii,jj,r50 , istatus)

      if (istatus .ne. 1) then! test weight_field
         write (6,*) 'failure in weight_field from analz_gps'
         return
      endif

c      call slv_laplc (data_weights, mask, ii,jj)

c     test nan values coming out of data_out and data_weights

      call check_nan2 (data_out,ii,jj,istatus)
      if (istatus.ne.1) then 
         write (6,*) 'nan detected in data_out in analz_gps.f'
         return
      endif
      call check_nan2 (data_weights,ii,jj,istatus)
      if (istatus.ne.1) then 
         write (6,*) 'nan detected in data_weights in analz_gps.f'
         return
      endif

      istatus = 1

      return

      end


