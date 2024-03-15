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

      subroutine read_madis_gps (path, filename, time_diff,
     1     gps_tpw, gps_error, gps_lat, 
     1     gps_lon,gps_num,
     1     gps_n, istatus)


      implicit none

      include 'netcdf.inc'

c     parameter list variables
      integer gps_n, gps_num
      real gps_tpw(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)
      integer time_diff
      character*256 path 
      character*9 filename
      character*13 filefound, cvt_i4time_wfo_fname13
      

c     internal
      integer istatus, ptg_index,i4time
      integer file_name_length

      integer recnum, nf_fid, nf_vid, nf_status
      character*120 extension
      integer extension_index
      character*120 desired_ext
      integer de_index


c     prep code
      call s_len(path, ptg_index)
c     create i4time locally from input filename (time) reference variable
      call i4time_fname_lp (filename, i4time, istatus)
c     adjust time to read prior hour (madis data of current hour will not
c     contain any data
      i4time = i4time - 3600 ! 3600 = 1 hour, subtraction, one hour earlier
c     create filefound (wfo mode name) from local i4time just generated
      filefound = cvt_i4time_wfo_fname13(i4time)
c
c  open desired netcdf file for reading
c
      nf_status = nf_open(path(1:ptg_index)//'/'//
     1     filefound,
     1     nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         istatus = 0
         write(6,*) 'failure getting madis gps data'
         return
      else
         istatus = 1
   
      endif

c
c  fill all dimension values
c
c
c get size of recnum
c
      nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
        istatus = 0
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,recnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
        istatus = 0
      endif
      call read_madis_gps_data (nf_fid , recnum, gps_tpw, gps_error, 
     1     gps_lat, gps_lon,
     1     gps_n,gps_num)

      

      return
      end

 
c
c
      subroutine read_madis_gps_data (nf_fid,recnum,gps_tpw,gps_error, 
     1     gps_lat, gps_lon,
     1     gps_n,gps_num)
      include 'netcdf.inc'
c     parameter list variables
      integer gps_n, gps_num
      real gps_tpw(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)
c
      integer recnum, nf_fid, nf_vid, nf_status,i,j
      character*80 stalongnam(recnum)
      character*5 stanam(recnum)
      real formalerror(recnum), stalat(recnum), stalon(recnum), 
     +   watervapor(recnum)
      double precision observationtime(recnum),latesttime

      call read_gps_madis_basics (nf_fid , recnum, formalerror, stalat,
     +    stalon,  watervapor, observationtime)
c
c the netcdf variables are filled - your code goes here
c
c     okay for the madis problem we are almost home.  we now have to 
c     make a subtle mod to the old reader to make it look like the old reader
c     in the old reader, the recnum was the total number of data, now
c     that is not the case.  since madis uses recnum differently.

c     but in the old code we did introduce gps_num and now we will actually 
c     use it as intended.

      latesttime = 0.

      do j = 1, recnum
         if(watervapor(j) .lt. 10000.) then ! good data
            if (latesttime .le. observationtime(j)) then
               latesttime = observationtime(j)
            endif
         endif
      enddo

c     at this point latesttime is the desired time to trap

      i = 0
      do j = 1, recnum
         
         if (watervapor(j) .lt. 10000. .and.
     +        latesttime .eq. observationtime(j)) then ! presume good data
            i = i +1
            gps_tpw(i) = watervapor(j)
            gps_error(i) = formalerror(j)
            gps_lat(i) = stalat(j)
            gps_lon(i) = stalon(j)
         endif
      enddo

      gps_num = i
      
      return
      end



      subroutine read_gps_madis_basics (nf_fid , recnum, formalerror, 
     +     latitude, longitude, totalcolumnpwv, observationtime)

      include 'netcdf.inc'
      integer recnum, nf_fid, nf_vid, nf_status



      real formalerror(recnum), totalcolumnpwv(recnum),
     +     latitude(recnum), longitude(recnum)
      double precision observationtime(recnum)
c


c
c     variable        netcdf long name
c      formalerror  "formal error" 
c

c     gather observation time

        nf_status = nf_inq_varid(nf_fid,'observationtime',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var observationtime'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,observationtime)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ observationtime '
      endif


c     get formal error


        nf_status = nf_inq_varid(nf_fid,'formalerror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var formalerror'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,formalerror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ formalerror '
      endif
c
c     variable        netcdf long name
c      stalat       "station latitude" 
c
        nf_status = nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalat'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,latitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stalat '
      endif
c
c     variable        netcdf long name
c      stalon       "station longitude" 
c
        nf_status = nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalon'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,longitude)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stalon '
      endif
c
c     variable        netcdf long name
c      watervapor   "water vapor" 
c
        nf_status = nf_inq_varid(nf_fid,'totalcolumnpwv',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var watervapor'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,totalcolumnpwv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ watervapor '
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end
