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
      subroutine  process_gps (ii,jj,gps_data_out,gps_data_weights,
     1     tpw,glat,glon,time_diff,gps_points,idotj,gps_count,
     1     path,filetime,gps_switch,istatus)

      implicit none

c     input variables


      character*9 filetime
      integer ii,jj,istatus,gps_switch
      integer time_diff         !time allowed for latency (sec)
      real gps_data_out(ii,jj),tpw(ii,jj)
      real glat(ii,jj), glon(ii,jj)
      integer i4time
      character*256 path
      real gps_data_weights(ii,jj)
      integer idotj
      integer gps_count
      real gps_points(idotj,3)
      
      integer gps_n, gps_num
      parameter (gps_n = 20000)
      real gps_tpw(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)

      integer i,j

c     call to normal public read left intact (gps_switch 1)

      if (gps_switch .eq. 1 ) then ! routine data read
         
         
         call read_gps (path, filetime, time_diff,
     1        gps_tpw, gps_error, gps_lat,
     1        gps_lon, gps_num, gps_n,
     1        istatus)
         
         if (
     1        istatus .ne. 1
     1        .or.
     1        gps_num .eq. 0
     1        ) then            ! failure
            
            write(6,*) 'failure to acquire gvap data'
            istatus = 0
         return                 !istatus = fail
         
      else
         
         write(6,*) gps_num, ' number of stations read in file'
         
      endif
      
      elseif (gps_switch .eq. 2) then ! madis read

         call read_madis_gps (path, filetime, time_diff,
     1        gps_tpw, gps_error, gps_lat,
     1        gps_lon, gps_num, gps_n,
     1        istatus)
         
         if (
     1        istatus .ne. 1
     1        .or.
     1        gps_num .eq. 0
     1        ) then            ! failure
            
            write(6,*) 'failure to acquire gvap data'
            istatus = 0
            return              !istatus = fail
            
         else
            
            write(6,*) gps_num, ' number of stations read in file'
            
         endif
         
      endif !  end figuring out which gps to read
      
      


c     carry on as before with gps processing

      call analz_gps (gps_lat,gps_lon,gps_tpw,gps_num,glat,
     1     glon,gps_data_out,
     1     gps_data_weights,
     1     gps_points,idotj,gps_count,ii,jj,istatus)
      if(istatus.ne.1) then ! failure to get data
         return
      endif

      call check_nan2(gps_data_out,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'data_out corrupted in processing gps data '
         write(6,*) 'var:gps_data_out     routine:process_gps.v'
         return
      endif
      call check_nan2(gps_data_weights,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'data_out corrupted in processing gps data '
         write(6,*) 'var:gps_data_weights     routine:process_gps.v'
         return
      endif

      istatus = 1

      return

      end
