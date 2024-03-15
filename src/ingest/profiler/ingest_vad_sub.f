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

      subroutine get_vad_data(i4time_sys,ilaps_cycle_time
     1                                    ,nx_l,ny_l
     1                                    ,filename,istatus)

      character*170 filename

!.............................................................................

      include 'netcdf.inc'
      integer maxlevels, recnum,nf_fid, nf_vid, nf_status
c
c  open netcdf file for reading
c
      nf_status = nf_open(filename,nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_open',filename
      endif
c
c  fill all dimension values
c
c
c get size of maxlevels
c
      nf_status = nf_inq_dimid(nf_fid,'maxlevels',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxlevels'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,maxlevels)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim maxlevels'
      endif
c
c get size of recnum
c
      nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,recnum)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'dim recnum'
      endif
      call main_sub(nf_fid, maxlevels, recnum,
!.............................................................................
     1              i4time_sys,ilaps_cycle_time,nx_l,ny_l,istatus)
!.............................................................................

      return
      end
c
c
      subroutine main_sub(nf_fid, maxlevels, recnum,
!.............................................................................
     1              i4time_sys,ilaps_cycle_time,nx_l,ny_l,istatus)

      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

      character*9 a9time_ob
      character*4 c4_staname

!.............................................................................

      include 'netcdf.inc'
      integer maxlevels, recnum,nf_fid
      integer numlevels(recnum), winddir( maxlevels, recnum)
      real obheight( maxlevels, recnum), staelev(recnum),
     +     stalat(recnum), stalon(recnum), windspeed( maxlevels,
     +     recnum)
      double precision timeobs(recnum)
      character*4 staname(recnum)
      character*350 rawmsg(recnum)
      character rmserror( maxlevels, recnum)

!............................................................................

      call read_netcdf(nf_fid, maxlevels, recnum, numlevels, winddir, 
     +     obheight, staelev, stalat, stalon, windspeed, timeobs, 
     +     rawmsg, rmserror, staname)
c
c the netcdf variables are filled - your code goes here
c
!............................................................................

      call get_latlon_perimeter(nx_l,ny_l,1.0
     1                           ,lat_a,lon_a,topo_a
     1                           ,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error reading laps perimeter'
          return
      endif

      num_vad = recnum
      write(6,*)' # of vad = ',recnum

      do ista = 1,num_vad

          write(6,*)
          write(6,*)' vad #',ista

!         write(6,*)' location = '
!    1             ,stalat(ista),stalon(ista),staelev(ista)

          if(stalat(ista) .le. rnorth .and. stalat(ista) .ge. south 
     1                                .and.
     1       stalon(ista) .ge. west   .and. stalon(ista) .le. east      
     1                                                             )then
              continue
          else ! outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,stalat(ista),stalon(ista)
              goto 900
          endif


          if(abs(timeobs(ista))      .lt. 3d9)then
              call c_time2fname(nint(timeobs(ista)),a9time_ob)
          else
              write(6,*)' bad observation time - reject',timeobs(ista)
              goto 900
          endif

          call cv_asc_i4time(a9time_ob,i4time_ob)
          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. (ilaps_cycle_time / 2) )then ! outside time window
              write(6,*)' time - reject '
     1           ,a9time_ob,i4_resid,ilaps_cycle_time / 2
              goto 900        
          endif

          write(6,1)a9time_ob
 1        format(' ob time:',1x,a9)

          id_num = 0

          c4_staname = 'k'//staname(ista)(1:3)

          if(abs(staelev(ista)) .gt. 999999.)then
              staelev_out = 999999.
          else
              staelev_out = staelev(ista)
          endif

          write(6,401)id_num,numlevels(ista),
     1                stalat(ista),stalon(ista),staelev_out,
     1                c4_staname,a9time_ob,'vad     '
          write(1,401)id_num,numlevels(ista),
     1                stalat(ista),stalon(ista),staelev_out,
     1                c4_staname,a9time_ob,'vad     '
401       format(i12,i12,f11.3,f15.3,f15.0,5x,a4,5x,a9,1x,a8)

          do ilvl = 1,numlevels(ista)
!             set the rms error for the level
              rms = 10.
              if(rmserror(ilvl,ista)(1:1) .eq. 'a')rms = 1.03
              if(rmserror(ilvl,ista)(1:1) .eq. 'b')rms = 2.06
              if(rmserror(ilvl,ista)(1:1) .eq. 'c')rms = 3.09
              if(rmserror(ilvl,ista)(1:1) .eq. 'd')rms = 4.12
              if(rmserror(ilvl,ista)(1:1) .eq. 'e')rms = 5.14
              if(rmserror(ilvl,ista)(1:1) .eq. 'f')rms = 6.17
              if(rmserror(ilvl,ista)(1:1) .eq. 'g')rms = 7.20

              write(1,301,err=303)obheight(ilvl,ista),
     1                            float(winddir(ilvl,ista)), 
     1                            windspeed(ilvl,ista),rms       
              write(6,301,err=303)obheight(ilvl,ista),
     1                            float(winddir(ilvl,ista)), 
     1                            windspeed(ilvl,ista),rms       
301           format(1x,f6.0,f6.0,2f6.1)
303           continue
          enddo ! ilvl

 900  enddo ! i

!............................................................................

      return
      end
c
c  subroutine to read the file "vad wind data : selected by receipt time : time range from 883602000 to 883603800" 
c
      subroutine read_netcdf(nf_fid, maxlevels, recnum, numlevels, 
     +     winddir, obheight, staelev, stalat, stalon, windspeed, 
     +     timeobs, rawmsg, rmserror, staname)
c
      include 'netcdf.inc'
      integer maxlevels, recnum,nf_fid, nf_vid, nf_status
      integer numlevels(recnum), winddir( maxlevels, recnum)
      real obheight( maxlevels, recnum), staelev(recnum),
     +     stalat(recnum), stalon(recnum), windspeed( maxlevels,
     +     recnum)
      double precision timeobs(recnum)
      character*4 staname(recnum)
      character*350 rawmsg(recnum)
      character rmserror( maxlevels, recnum)


c   variables of type real
c
c     variable        netcdf long name
c      obheight     "observation height above msl"
c
        nf_status = nf_inq_varid(nf_fid,'obheight',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var obheight'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,obheight)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var obheight'
      endif
c
c     variable        netcdf long name
c      staelev      "elevation of station above msl"
c
        nf_status = nf_inq_varid(nf_fid,'staelev',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staelev'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,staelev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staelev'
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
        print *,'in var stalat'
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
        print *,'in var stalon'
      endif
c
c     variable        netcdf long name
c      windspeed    "wind speed"
c
        nf_status = nf_inq_varid(nf_fid,'windspeed',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif
        nf_status = nf_get_var_real(nf_fid,nf_vid,windspeed)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var windspeed'
      endif

c   variables of type int
c
c
c     variable        netcdf long name
c      numlevels    "number of observations for station"
c
        nf_status = nf_inq_varid(nf_fid,'numlevels',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numlevels'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,numlevels)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var numlevels'
      endif
c
c     variable        netcdf long name
c      winddir      "wind direction"
c
        nf_status = nf_inq_varid(nf_fid,'winddir',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif
        nf_status = nf_get_var_int(nf_fid,nf_vid,winddir)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var winddir'
      endif

c   variables of type double
c
c
c     variable        netcdf long name
c      timeobs      "time of observation"
c
        nf_status = nf_inq_varid(nf_fid,'timeobs',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif
        nf_status = nf_get_var_double(nf_fid,nf_vid,timeobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var timeobs'
      endif


c   variables of type char
c
c
c     variable        netcdf long name
c      rawmsg       "receipt format vad winds message"
c
        nf_status = nf_inq_varid(nf_fid,'rawmsg',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rawmsg'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,rawmsg)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rawmsg'
      endif
c
c     variable        netcdf long name
c      rmserror     "rms vector wind error"
c
        nf_status = nf_inq_varid(nf_fid,'rmserror',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rmserror'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,rmserror)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var rmserror'
      endif
c
c     variable        netcdf long name
c      staname      "three letter wmo site identifier"
c
        nf_status = nf_inq_varid(nf_fid,'staname',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staname'
      endif
        nf_status = nf_get_var_text(nf_fid,nf_vid,staname)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'in var staname'
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'nf_close'
      endif

      return
      end


