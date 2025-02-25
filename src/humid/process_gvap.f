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
      subroutine  process_gvap (ii,jj,sfc_data,data_out,data_weights,
     1     gw1,gw2,gw3,gww1,gww2,gww3,gvap_p,mdf,
     1     glat,glon,time_diff,ihop_flag,
     1     path_to_gvap12,path_to_gvap10,filetime,print_switch,istatus)

      use module_sfc_structure

      implicit none

      

c     input variables

      type (lbsi), dimension(ii,jj) :: sfc_data

      character*9 filename,filetime
      integer ii,jj,istatus,print_switch
      integer time_diff         !time allowed for latency (sec)
      real data_out(ii,jj)
      real gw1(ii,jj),gww1(ii,jj)
      real gw2(ii,jj),gww2(ii,jj)
      real gw3(ii,jj),gww3(ii,jj)
      real glat(ii,jj), glon(ii,jj)
      real gvap_p (ii,jj)
      real mdf
      integer i4time
      integer ihop_flag
      character*256 path_to_gvap12,path_to_gvap10
      real data_weights(ii,jj)

      integer nstations,nn
      parameter (nstations = 100000)
      real lat(nstations)
      real lon(nstations)
      real wt(nstations)
      real w1(nstations)
      real w2(nstations)
      real w3(nstations)
      real gvap_pres(nstations)    !surface pressure for sigma coord



      integer i,j


c     code

      filename = filetime

      do j = 1,jj
         do i = 1,ii
            gww1(i,j) = mdf
            gww2(i,j) = mdf
            gww3(i,j) = mdf
            gw1(i,j) = mdf
            gw2(i,j) = mdf
            gw3(i,j) = mdf
         enddo
      enddo


      call read_gvap (filename, nstations, path_to_gvap12,
     1     path_to_gvap10,
     1     time_diff, ihop_flag, lat,lon, wt,w1,w2,w3, gvap_pres, nn,
     1     istatus)

      if (
     1     istatus .ne. 1
     1     .or.
     1     nn .eq. 0
     1     ) then               ! failure

         write(6,*) 'failure to acquire gvap data'
         istatus = 0
         return                 !istatus = fail

      else


         write(6,*) nn, ' number of stations read in file'
         if (nn.gt.nstations) then ! exceeded dimension
            istatus  = 0
            write (6,*) 'nstations exceeded (parameter dimension)'
            write (6,*) 'readjust and recompile code'
            write (6,*) 'warning only, not fatal'
            write (6,*) 'gvap data not used'
            return
         endif
c     correct longitute to negative for west
         do i = 1,nn
            lon(i) = lon(i) * (-1.0)
         enddo

      endif

      call analz_gvap (lat,lon,wt,w1,w2,w3,gvap_pres,nn,
     1     glat,glon,sfc_data,data_out,
     1     gw1,gw2,gw3,gww1,gww2,gww3,gvap_p,
     1     data_weights,ii,jj,print_switch,istatus)

c     temporary assignment made until gvap_p read okay from data bases
c      gvap_p = mdf

      if(istatus.ne.1) then ! failure to get data
         return
      endif

      call check_nan2(data_out,ii,jj,istatus)
      if(istatus.ne.1) then
         write(6,*) 'nan detected var:data_out  routine:process_gvap.f'
         return
      endif

      istatus = 1

      return

      end
