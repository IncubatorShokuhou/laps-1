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

      subroutine read_gps (path, filename, time_diff,
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
      character*9 filefound
      

c     internal
      integer istatus, ptg_index
      integer file_name_length

      integer recnum, nf_fid, nf_vid, nf_status
      character*120 extension
      integer extension_index
      character*120 desired_ext
      integer de_index


c     prep code
      call s_len(path, ptg_index)

      file_name_length = 14
      desired_ext = 'nc'


!gfortran modifications begin
      extension_index = 0
!gfortran modifications end
      call get_newest_file (filename, time_diff,file_name_length,
     1     path,ptg_index,filefound,desired_ext, de_index,
     1     extension, extension_index, istatus)

!gfortran modifications begin
      if (istatus.ne.1) then    !failure, return with failure code
         write(6,*) 'no gps data'
         return
      endif
!gfortran modifications end
c     prep filename to fit


c
c  open netcdf file for reading
c
      nf_status = nf_open(path(1:ptg_index)//'/'//
     1     filefound//'0030o.'//extension(1:extension_index),
     1     nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
         print *, nf_strerror(nf_status)
         istatus = 0
         write(6,*) 'failure getting gps data'
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
      call read_gps_data (nf_fid , recnum, gps_tpw, gps_error, 
     1     gps_lat, gps_lon,
     1     gps_n,gps_num)

      return
      end

 
c
c
      subroutine read_gps_data (nf_fid, recnum, gps_tpw, gps_error, 
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
      integer recnum, nf_fid, nf_vid, nf_status,i
      character*80 stalongnam(recnum)
      character*5 stanam(recnum)
      real formalerror(recnum), stalat(recnum), stalon(recnum), 
     +   watervapor(recnum)
      call read_gps_basics (nf_fid , recnum, formalerror, stalat,
     +    stalon, stalongnam, stanam, watervapor)
c
c the netcdf variables are filled - your code goes here
c

      gps_num = recnum
      do i = 1, recnum
         gps_tpw(i) = watervapor(i)
         gps_error(i) = formalerror(i)
         gps_lat(i) = stalat(i)
         gps_lon(i) = stalon(i)
      enddo

      return
      end
      subroutine read_gps_basics (nf_fid , recnum, formalerror, stalat,
     +    stalon, stalongnam, stanam, watervapor)
      include 'netcdf.inc'
      integer recnum, nf_fid, nf_vid, nf_status

      character*80 stalongnam(recnum)
      character*5 stanam(recnum)
      real formalerror(recnum), stalat(recnum), stalon(recnum), 
     +   watervapor(recnum)
c
c     variable        netcdf long name
c      stalongnam   "station location" 
c
        nf_status = nf_inq_varid(nf_fid,'stalongnam',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalongnam'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,stalongnam)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stalongnam '
      endif
c
c     variable        netcdf long name
c      stanam       "alphanumeric station name" 
c
        nf_status = nf_inq_varid(nf_fid,'stanam',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stanam'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,stanam)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stanam '
      endif
c
c     variable        netcdf long name
c      formalerror  "formal error" 
c
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
        nf_status = nf_inq_varid(nf_fid,'stalat',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalat'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stalat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stalat '
      endif
c
c     variable        netcdf long name
c      stalon       "station longitude" 
c
        nf_status = nf_inq_varid(nf_fid,'stalon',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var stalon'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,stalon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ stalon '
      endif
c
c     variable        netcdf long name
c      watervapor   "water vapor" 
c
        nf_status = nf_inq_varid(nf_fid,'watervapor',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var watervapor'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,watervapor)
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
