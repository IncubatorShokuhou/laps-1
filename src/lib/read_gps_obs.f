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

      subroutine read_gps_obs (lun, path, i4beg, i4end,
     1     imax, jmax, latgrid, longrid, bad_sfc,
     1     gps_tpw, gps_wet, gps_error, gps_xy, gps_elv, 
     1     gps_tim, gps_indomain, gps_n, istatus)
c     this code is originated from dan birkenheuer.

c     on july 16, 2009, yuanfu xie modified it to read an additional 
c     variable, wet delays, station elevation and gps obs time.

c     on august 2, 2009, yuanfu xie modified it to read all gps data
c     files between i4begin and i4end.

c     reworked a bit by steve albers in 2011

      implicit none

      include 'netcdf.inc'

c     parameter list variables
      integer max_files
      parameter (max_files=3000)
      integer lun, gps_n, gps_num, gps_indomain, verbose
      integer i4times(max_files),nfiles,gps_i,ifile,i,istat_latlon
      integer imax, jmax
      real latgrid(imax,jmax)
      real longrid(imax,jmax)
      real gps_tpw(gps_n)
      real gps_wet(gps_n)
      real gps_error(gps_n)
      real gps_xy(2,gps_n)
      real gps_elv(gps_n)
      real gps_tim(gps_n)

      real gps_lat(gps_n)
      real gps_lon(gps_n)
      integer i4beg, i4end
      character*256 path,filespec,c_filenames(max_files) 
      
c     internal
      integer istatus, ptg_index
      integer file_name_length

      integer recnum, nf_fid, nf_vid, nf_status
      real x, y, bad_sfc

      verbose = 0

c     set file specs for get_file_times:
      call s_len(path, ptg_index)
      filespec(1:ptg_index) = path(1:ptg_index)
      filespec(ptg_index+1:ptg_index+9) = '*0030o.nc'  ! hardcode for now by yuanfu

c     get all filenames under the path:
      call get_file_times(filespec(1:ptg_index+9),max_files,c_filenames,
     1                    i4times,nfiles,istatus)

      gps_i = 0    ! total gps data read
c     loop through all available files: yuanfu
      do ifile=1,nfiles
      
c
c       open netcdf file for reading
c
        ! read data file between i4beg and i4end:
        if (i4times(ifile) .ge. i4beg .and. 
     1      i4times(ifile) .le. i4end) then
        nf_status = nf_open(c_filenames(ifile),nf_nowrite,nf_fid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          istatus = 0
          write(6,*) 'failure getting gps data'
          nf_status = nf_close(nf_fid)
          return
        else
          istatus = 1
          write(6,*)' opened gps file: ',trim(c_filenames(ifile))
        endif
c
c       fill all dimension values
c
c
c       get size of recnum
c
        nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'dim recnum'
          istatus = 0
          goto 900  
        endif
        nf_status = nf_inq_dimlen(nf_fid,nf_vid,recnum)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'dim recnum'
          istatus = 0
          goto 900  
        endif

        ! check if recnum is larger than space allocated:
        if (gps_i+recnum .gt. gps_n) then
          print *,' too many gps obs, increase your gps_n and rerun!'
     1	,gps_i+recnum,' > ',gps_n
          stop
        endif

        call read_gps_data (nf_fid , recnum, gps_tpw(gps_i+1), 
     1       gps_wet(gps_i+1), gps_error(gps_i+1), gps_lat(gps_i+1), 
     1       gps_lon(gps_i+1), gps_elv(gps_i+1), gps_tim(gps_i+1), 
     1       gps_n,gps_num)

c       accumulate all:
        gps_i = gps_i+gps_num
        write(6,*)' gps_num / gps_i = ',gps_num,gps_i

        ! finish read this qualified file:
        endif

 900  enddo
      ! check if there is any gps obs:
      if (gps_i .ge. 1) istatus = 1
      if (gps_i .le. 0) istatus = 0

      ! write the gps wetdelay information into laps hmg file:
      gps_indomain = 0
      do i=1,gps_i ! loop through all obs from the files
        call latlon_to_rlapsgrid(gps_lat(i),gps_lon(i),latgrid,
     1             longrid, imax, jmax, x, y, istat_latlon)

        if(verbose .ge. 1 .or. i .le. 10)then
          write(6,101)i,gps_lat(i),gps_lon(i),x,y,gps_wet(i)
 101      format(i7,4f9.2,e16.6)
        endif

        if (x .ge. 1 .and. x .le. imax .and.
     1      y .ge. 1 .and. y .le. jmax .and. 
     1    gps_wet(i) .ne. bad_sfc) then
          gps_indomain = gps_indomain + 1
          gps_tpw(gps_indomain) = gps_tpw(i)
          gps_wet(gps_indomain) = gps_wet(i)
          gps_error(gps_indomain) = gps_error(i)
          gps_elv(gps_indomain) = gps_elv(i)
          gps_tim(gps_indomain) = gps_tim(i)
          gps_xy(1,gps_indomain) = x
          gps_xy(2,gps_indomain) = y

          if(lun .gt. 0)then! write wet delays and tpw into hmg file:
            write(lun, *) x, y,0, gps_wet(i), 'gpswet'
            write(lun, *) x, y,0, gps_tpw(i), 'gpstpw'
          endif

          if(verbose .ge. 1)then
            write(6, *) x, y,0, gps_wet(i), 'gpswet'
            write(6, *) x, y,0, gps_tpw(i), 'gpstpw'
          endif

        endif
      enddo

      return

      end

c
c
      subroutine read_gps_data (nf_fid, recnum, gps_tpw, 
     1     gps_wet, gps_error, gps_lat, gps_lon, gps_elv,
     1     gps_tim, gps_n,gps_num)
      include 'netcdf.inc'

c     this code is originated from dan birkenheuer.
c     on july 16, 2009, yuanfu xie modified it to read an additional
c     variable, wet delays, station elevation and gps obs time.

c     parameter list variables
      integer gps_n, gps_num
      real gps_tpw(gps_n)
      real gps_wet(gps_n)
      real gps_error(gps_n)
      real gps_lat(gps_n)
      real gps_lon(gps_n)
      real gps_elv(gps_n)
      real gps_tim(gps_n)
c
      integer recnum, nf_fid, nf_vid, nf_status,i
      character*80 stalongnam(recnum)
      character*5 stanam(recnum)
      real formalerror(recnum), stalat(recnum), stalon(recnum), 
     +     staelv(recnum), timobs(recnum), watervapor(recnum), 
     +     wetdelay(recnum)
      call read_gps_basics (nf_fid , recnum, formalerror, stalat,
     +    stalon, staelv, timobs, stalongnam, stanam, watervapor, 
     +    wetdelay) 
      ! yuanfu: wetdelay, staelv -- station elevation, and gps obs time
c
c the netcdf variables are filled - your code goes here
c

      gps_num = recnum
      do i = 1, recnum
         gps_tpw(i) = watervapor(i)
         gps_wet(i) = wetdelay(i)
         gps_error(i) = formalerror(i)
         gps_lat(i) = stalat(i)
         gps_lon(i) = stalon(i)
         gps_elv(i) = staelv(i)
         gps_tim(i) = timobs(i)
      enddo

      return
      end
      subroutine read_gps_basics (nf_fid , recnum, formalerror, 
     +           stalat,stalon, staelv, timobs, stalongnam, 
     +           stanam, watervapor,wetdelay)	
      ! yuanfu: add wetdelay and staelv -- station elevation

c     this code is originated from dan birkenheuer.
c     on july 16, 2009, yuanfu xie modified it to read an additional
c     variable, wet delays and station elevation.

      include 'netcdf.inc'
      integer recnum, nf_fid, nf_vid, nf_status

      character*80 stalongnam(recnum)
      character*5 stanam(recnum)
      real formalerror(recnum), stalat(recnum), stalon(recnum), 
     +     staelv(recnum), timobs(recnum), watervapor(recnum), 
     +     wetdelay(recnum)	
      ! yuanfu: add wetdelay and staelv
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
c      staelv       "station elevation" 
c
        nf_status = nf_inq_varid(nf_fid,'staelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staelev'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,staelv)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ staelev '
      endif
c
c     variable        netcdf long name
c      timobs       "time of observation" 
c
        nf_status = nf_inq_varid(nf_fid,'timeobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,timobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ timeobs '
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
c
c     variable        netcdf long name
c      wetdelay   "wet component gps signal delay" 
c
        nf_status = nf_inq_varid(nf_fid,'wetdelay',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var wetdelay'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,wetdelay)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in nf_get_var_ wetdelay '
      endif


c     close file:
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end


      subroutine get_gps_path(path_to_gps_out,istatus)

      use mem_namelist, only: path_to_gps

      character*256 path_to_gps_out

      logical l_exist

      path_to_gps_out = path_to_gps 

      inquire(file=trim(path_to_gps_out),exist=l_exist)
      if(l_exist .eqv. .true.)then
          write(6,*)' path to gps from namelist exists'         
      else
          write(6,*)
     1              ' path to gps non-existent, trying /public'
          path_to_gps_out = '/public/data/gpsmet/netcdf/'
      endif

      istatus = 1

      return
      end
