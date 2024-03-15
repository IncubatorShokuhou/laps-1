c
        subroutine get_local_towerobs(maxsta,maxlvls,                    ! i
     &                 i4time_sys,lun_out,
     &                 path_to_local_data,tower_format,ext,
     &                 itime_before,itime_after,
     &                 lat,lon,ni,nj,                                    ! i
     &                 nsta,                                             ! o
     &                 stalat_s,stalon_s,staelev_s,                      ! o
     &                 stname_s,                                         ! o
     &                 soilmoist_p,                                      ! o
     &                 istatus)

c       read tower data from either rsa or nimbus netcdf cdls
c
c.....  input variables/arrays
c
        integer maxlvls ! raw/processed stations for snd file
        integer maxobs ! raw stations in netcdf files
        integer maxsta ! processed stations for snd file 

        parameter (maxobs=10000) ! raw stations in netcdf files

        character*(*) path_to_local_data, tower_format, ext

        real    lat(ni,nj), lon(ni,nj)
c
c.....  local variables/arrays
c
        double precision d_timeobs

!       obs arrays (raw files)
        integer     nobs,nlvl(maxobs)
        real      lvls_m(maxlvls,maxobs)
        real      fill_stationp, fill_lvls, stationp
        integer   i4time
        character*51  stationname
        character*6   c_staid
	real      dd(maxlvls,maxobs), ff(maxlvls,maxobs)

        integer  i4time_ob_a(maxobs), before, after
c
c.....  variables returned from 'read_local_tower'
c
        integer    nsnd_all ! combined # of obs over multiple files
	integer    wmoid(maxobs)
        real     stalat(maxobs),stalon(maxobs)
        real     staelev(maxobs)
        real     soilmoisture(maxobs)
        character  c5_staid(maxobs)*5, a9time_ob(maxobs)*9
        character  a9time*9
!       character  c8_obstype(maxobs)*8
        real     height_m(maxobs,maxlvls), pressure_mb(maxobs,maxlvls)
        real     temp_c(maxobs,maxlvls), dewpoint_c(maxobs,maxlvls)
        real     dir_deg(maxobs,maxlvls),spd_mps(maxobs,maxlvls)
	character  stname(maxobs)*6
        integer    tempqcflag(maxobs,maxlvls)
        integer    prsqcflag(maxobs,maxlvls)
        integer    rhqcflag(maxobs,maxlvls)
        integer    wsqcflag(maxobs,maxlvls)
        integer    wdqcflag(maxobs,maxlvls)
        integer    smqcflag(maxobs)

        logical    l_closest_time(maxobs), l_closest_time_i
c
c.....  output arrays used by 'write_snd' 
c
        integer    nsnd_all_s ! combined # of obs over multiple files
        integer    nlvl_s(maxsta)
	integer  wmoid_s(maxsta)
        real     stalat_s(maxsta,maxlvls),stalon_s(maxsta,maxlvls)
        real     staelev_s(maxsta)
        character  c5_staid_s(maxsta)*5, a9time_ob_s(maxsta,maxlvls)*9
        character  c8_obstype_s(maxsta)*8
        real     height_m_s(maxsta,maxlvls)
        real     pressure_mb_s(maxsta,maxlvls)
        real     temp_c_s(maxsta,maxlvls)
        real     dewpoint_c_s(maxsta,maxlvls)      
        real     dir_deg_s(maxsta,maxlvls),spd_mps_s(maxsta,maxlvls)       
        real     soilmoist_p(maxsta)       
	character  stname_s(maxsta)*5

c.....  unknown vars.
	character  save_stn(maxsta)*6

	integer    rtime
	integer    recnum, nf_fid, nf_vid, nf_status
	character  timech*9, time*4
        character*13 filename13, cvt_i4time_wfo_fname13
        character*150 data_file 
c
c.....  start.
c
c       get r_missing_data
        call get_r_missing_data(r_missing_data,istatus)
c       
c.....	set istatus flag for the local data to bad until we find otherwise.
c
	istatus = -1

        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

!       call get_box_size(box_size,istatus)
!       if(istatus .ne. 1)return

        box_size = 1.0

c
c.....  figure out the size of the "box" in gridpoints.  user defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c

        call get_grid_spacing_cen(grid_spacing_m,istatus)
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing_m / 1000.) !in km
c
c.....	zero out the counters.
c
        nobs = 0	        ! # of local obs in the laps grid
c
c.....  get the data from the netcdf file(s).  first, open the file(s).
c.....  if none are there, return to tower_driver.
c
        ix = 1
c
c.....  set up the time window.
c
	i4time_before = i4time_sys - itime_before
	i4time_after  = i4time_sys + itime_after

!       ob times contained in each file
        i4_contains_early = 0 
        i4_contains_late = 3599

        call get_filetime_range(i4time_before,i4time_after                
     1                         ,i4_contains_early,i4_contains_late       
     1                         ,3600                                     
     1                         ,i4time_file_before,i4time_file_after)

        write(6,*)' i4time_file_before,i4time_file_after:'
     1             ,i4time_file_before,i4time_file_after

        write(6,*)' tower format = ',tower_format

        do i4time_file = i4time_file_before, i4time_file_after, +3600       

            call s_len(path_to_local_data,len_path)

            if(tower_format .eq. 'wfo' .or. tower_format .eq. 'rsa')then
                filename13= cvt_i4time_wfo_fname13(i4time_file)
                if(len_path .lt. 1)goto590
 	        data_file = path_to_local_data(1:len_path)//'/'
     1                                                    //filename13
            else
                call make_fnam_lp(i4time_file,a9time,istatus)
                if(len_path .lt. 1)goto590
 	        data_file = path_to_local_data(1:len_path)//'/'
     1                                       //a9time//'0100o'       
            endif

            write(6,*)' ldad tower file = ',trim(data_file)
            call s_len(data_file,lenf)
c
c.....  call the read routine.
c
	    call read_local_tower(data_file,lenf,                 ! i 
     &         tower_format,                                      ! i
     &         maxobs, maxlvls,                                   ! i
     &         r_missing_data,                                    ! i
     &         nsnd_file, nlvl(ix), lvls_m(1,ix),                 ! o
     &         staelev(ix), stalat(ix), stalon(ix),               ! o
     &         soilmoisture(ix),                                  ! o
     &         temp_c(ix,1), dewpoint_c(ix,1),                    ! o
     &         height_m(ix,1),                                    ! o
     &         dir_deg(ix,1), spd_mps(ix,1),                      ! o
c    &         pressure_mb(ix,1), prsqcflag(ix,1),                ! o
     &         tempqcflag(ix,1), rhqcflag(ix,1),                  ! o
     &         wsqcflag(ix,1), wdqcflag(ix,1),                    ! o
     &         smqcflag(ix),                                      ! o
     &         a9time_ob(ix), stname(ix), wmoid(ix),              ! o
     &         istatus)                                           ! o
 
	    if(istatus .ne. 1)then
                write(6,*)
     1          '     warning: bad status return from read_local_tower'       
                nsnd_file = 0
            else
                write(6,*)'     nsnd_file = ',nsnd_file
            endif

            ix = ix + nsnd_file

590     enddo ! i4time_file

        nsnd_all = ix - 1
        write(6,*)' nsnd_all = ',nsnd_all

        if(nsnd_all .gt. maxobs)then
            write(6,*)' error: nsnd_all > maxobs ',nsnd_all,maxobs
            nsta = 0
            istatus = 0
            return
        endif
c
c.....  post process the soundings.....
c
        print *
	print *,'  appending snd file, # of obs (in grid) = ',nsnd_all       

        pressure_mb_s = r_missing_data

!       flag those reports that are the closest times for the station
        do i = 1,nsnd_all
            l_closest_time(i) = l_closest_time_i(wmoid,a9time_ob
     1                                          ,nsnd_all,i,i4time_sys       
     1                                          ,istatus)       
        enddo ! i

!       transfer arrays (with various qc steps)

        nsta = 0
        do i = 1,nsnd_all
            if(nlvl(i) .gt. 0 .and. nlvl(i) .le. maxlvls 
     1                        .and. l_closest_time(i)     )then 

!               valid sounding - use for output
                nsta = nsta + 1
                write(6,601)nsta,i,stname(i),wmoid(i),a9time_ob(i)
 601            format(' valid sounding - transferring to output arrays'      
     1                ,2i5,2x,a7,2x,i10,2x,a10)
                stalat_s(nsta,:) = stalat(i)           
                stalon_s(nsta,:) = stalon(i)           
                staelev_s(nsta) = staelev(i)           
                stname_s(nsta) = stname(i)(1:5)           

                wmoid_s(nsta) = wmoid(i)           
                a9time_ob_s(nsta,:) = a9time_ob(i)           
                c8_obstype_s(nsta) = 'tower   '

!               single level data
                soilmoist_p(nsta) = soilmoisture(i)

                nlvl_s(nsta) = nlvl(i)

!               qc the levels 
                do il = nlvl(i),1,-1
                    if(     lvls_m(il,i) .ge. 1e10 
!    1                 .or. lvls_m(il,i) .le. 0.          
     1                                             )then
                        write(6,*)' error: invalid lvls_m',i,wmoid(i)
     1                           ,il,lvls_m(il,i)
                        nlvl_s(nsta) = il-1
                    else
                        if(il .ge. 2)then
                            if(lvls_m(il,i) .le. lvls_m(il-1,i))then
                                write(6,*)' error: levels out of order'
     1                                   ,i,lvls_m(il-1,i),lvls_m(il,i)      
                                go to 1500
                            endif 
                        endif
                    endif
                enddo ! il

                do il = 1,nlvl_s(nsta)
                    height_m_s(nsta,il) = lvls_m(il,i) + staelev_s(nsta)       

                    if(iqc_rsa(tempqcflag(nsta,il)) .ne. -1)then
                        temp_c_s(nsta,il) = temp_c(i,il)           
                    else
                        temp_c_s(nsta,il) = r_missing_data
                    endif

                    if(iqc_rsa(rhqcflag(nsta,il)) .ne. -1)then
                        dewpoint_c_s(nsta,il) = dewpoint_c(i,il)           
                    else
                        dewpoint_c_s(nsta,il) = r_missing_data
                    endif

                    if(iqc_rsa(wdqcflag(nsta,il)) .ne. -1 .and.
     1                 iqc_rsa(wsqcflag(nsta,il)) .ne. -1      )then
                        dir_deg_s(nsta,il) = dir_deg(i,il)           
                        spd_mps_s(nsta,il) = spd_mps(i,il)           
                    else
                        dir_deg_s(nsta,il) = r_missing_data
                        spd_mps_s(nsta,il) = r_missing_data
                    endif
                enddo

                go to 1600 

 1500           write(6,*)' sounding rejected: ' ,nsta

                nsta = nsta - 1

 1600           continue ! normal status

            endif
        enddo ! i
c
c.....  call the routine to write the snd file.
c
        if(ext(1:3) .eq. 'snd')then

            if(nsta .gt. 0)then
                call open_ext(lun_out,i4time_sys,'snd',istatus)
            endif

            call write_snd(lun_out                               ! i
     1                    ,maxsta,maxlvl,nsta                    ! i
     1                    ,wmoid_s                               ! i
     1                    ,stalat_s,stalon_s,staelev_s           ! i
     1                    ,stname_s,a9time_ob_s,c8_obstype_s     ! i
     1                    ,nlvl_s                                ! i
     1                    ,height_m_s                            ! i
     1                    ,pressure_mb_s                         ! i
     1                    ,temp_c_s                              ! i
     1                    ,dewpoint_c_s                          ! i
     1                    ,dir_deg_s                             ! i
     1                    ,spd_mps_s                             ! i
     1                    ,istatus)                              ! o

            if(istatus .ne. 1)then
                write(6,*)
     1       ' get_local_towerobs: bad status returned from write_snd'       
            endif

            return

        elseif(ext(1:3) .eq. 'lso')then
            write(6,*)
     1       ' lso option, pass data back instead of writing data out'       

        else
            write(6,*)' error - unknown ext in get_local_towerobs',ext
            nsta = 0
            istatus = 0
            return

        endif
c
 990    continue               ! no data available
        istatus = 0
        print *,' no data available from get_local_towerobs'
        return
c
        end

         subroutine read_local_tower(filename,fn_len,             ! i 
     &         tower_format,                                      ! i
     &         maxobs, maxlvls,                                   ! i
     &         r_missing_data,                                    ! i
     &         nobs, nlvl, lvls_m,                                ! o
     &         staelev, stalat, stalon,                           ! o
     &         soilmoisture,                                      ! o
     &         temp_c, dewpoint_c, height_m,                      ! o
     &         dir_deg, spd_mps,                                  ! o
c    &         pressure_pa, prsqcflag,                            ! o 
     &         tempqcflag, rhqcflag, wsqcflag, wdqcflag,          ! o
     &         smqcflag,                                          ! o
     &         a9time_ob, stname, wmoid,                          ! o
     &         istatus)                                           ! o

      include 'netcdf.inc'

      character*(*) filename, tower_format 
      integer       maxobs ! raw stations for snd file
      integer       maxlvls ! raw/processed stations for snd file
      real        r_missing_data 
      integer       nobs,nlvls,lev_set
      real        lvls_m(maxlvls,maxobs)
      real        staelev(maxobs)
      real        stalat(maxobs),stalon(maxobs)
      real        soilmoisture(maxobs)
      real        dd(maxlvls,maxobs), ff(maxlvls,maxobs)
      real        temp_k, rh_pct,stationp,ws,wd
      integer       tempqcflag(maxobs,maxlvls)
      integer       prsqcflag(maxobs,maxlvls)
      integer       rhqcflag(maxobs,maxlvls)
      integer       wsqcflag(maxobs,maxlvls)
      integer       wdqcflag(maxobs,maxlvls)
      integer       smqcflag(maxobs)
      integer       tempqf, prsqf, rhqf, wsqf, wdqf
      real        height_m(maxobs,maxlvls)
      real        pressure_pa(maxobs,maxlvls)      
      real        temp_c(maxobs,maxlvls), dewpoint_c(maxobs,maxlvls)
      real        dir_deg(maxobs,maxlvls),spd_mps(maxobs,maxlvls)
      real        sp_fill,levels_fill
      real        fill_t, fill_ws, fill_rh
      integer       fill_wd
      integer       wmoid(maxobs),nlvl(maxobs),sp_id
      integer       istatus
      integer       sta_id,sn_id,wd_id,ws_id,rh_id,temp_id,lev_id,ot_id
      integer       tempqc_id,prsqc_id,rhqc_id,wsqc_id,wdqc_id
      integer       index_1(1), index_2(2),start(2),count(2)
      integer       start1(1),count1(1)
      integer       pi_len, sn_len, fn_len, obno
      character     stname(maxobs)*6, stationname*51, c_staid*6
      character     a9time_ob(maxobs)*9, rh_var*30
      double precision d_timeobs

c     open data_file
      nf_status = nf_open(filename(1:fn_len),nf_nowrite,nf_fid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error on nf_open - open file: ',filename(1:fn_len)
      endif

c     read dim recnum -> nobs
      nf_status = nf_inq_dimid(nf_fid,'recnum',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error finding dim recnum' 
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nobs)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error reading dim recnum'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif

c     verify nobs .lt. maxobs

      print *, 'number of records in file: ',nobs
      if (nobs .gt. maxobs) then
        print *,
     1'warning: nobs is greater than maxobs: ',
     1'nobs / maxobs',nobs,' / ',maxobs
        nobs = maxobs
        print *, 'warning: reading only the first ',
     1nobs,' records'
      endif
      
c     read dim level -> nlvls
      nf_status = nf_inq_dimid(nf_fid,'level',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error finding dim level' 
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,nlvls)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error reading dim level'
      endif

c lw
      print *, 'nlvls / nobs ',nlvls,' / ',nobs

c     verify nlvls .lt. maxlvls
      if (nlvls .gt. maxlvls) then
        print *,'nlvls is greater than maxlvls: ',nlvls,' ',maxlvls
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif

c     read dim level stanamlen
      nf_status = nf_inq_dimid(nf_fid,'stanamlen',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error finding dim stanamlen' 
        print *, 'set sn_len = 51'
        sn_len = 51
      else
        nf_status = nf_inq_dimlen(nf_fid,nf_vid,sn_len)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'error reading dim stanamlen'
          print *, 'set sn_len = 51'
          sn_len = 51
        endif
      endif

c     read dim level provideridlen
      nf_status = nf_inq_dimid(nf_fid,'provideridlen',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error finding dim provideridlen' 
        print *, 'set pi_len = 6'
        pi_len = 6
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,pi_len)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error reading dim provideridlen'
        print *, 'set pi_len = 6'
        pi_len = 6
      endif

c     setup start and count
      start1(1) = 1
      count1(1) = nobs

c     read var elevation(recnum) -> staelev(maxobs)
      nf_status = nf_inq_varid(nf_fid,'elevation',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var elevation'
      endif
      nf_status = nf_get_vara_real(nf_fid,nf_vid,
     1 start1,count1,staelev)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading var elevation'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif 

c     read var latitude(recnum -> stalat(maxobs)
      nf_status = nf_inq_varid(nf_fid,'latitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var latitude'
      endif
      nf_status = nf_get_vara_real(nf_fid,nf_vid,
     1 start1,count1,stalat)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading var latitude'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif 

c     read var longitude(recnum) -> stalon(maxobs)
      nf_status = nf_inq_varid(nf_fid,'longitude',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var longitude'
      endif
      nf_status = nf_get_vara_real(nf_fid,nf_vid,
     1 start1,count1,stalon)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading var longitude'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif 

c     read var soilmoisture(recnum) -> soilmoisture(maxobs)
      nf_status = nf_inq_varid(nf_fid,'soilmoisture',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'warning: could not find var soilmoisture'
      else
        nf_status = nf_get_vara_real(nf_fid,nf_vid,
     1                               start1,count1,soilmoisture)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading var soilmoisture'
          print *, 'aborting read'
          nf_status = nf_close(nf_fid)
          istatus = 0
          return
        endif 
      endif
      
c     read dim smqcflag -> smqcflag
      nf_status = nf_inq_dimid(nf_fid,'smqcflag',nf_vid)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'could not find dim smqcflag' 
      endif
      nf_status = nf_inq_dimlen(nf_fid,nf_vid,smqcflag)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'could not read dim smqcflag'
      endif

c lw
      print *, 'b4 read levels: nlvls / nobs ',nlvls,' / ',nobs

c     read var levels(recnum,level) -> lvls_m(maxlvls,maxobs)
      start(1) = 1
      start(2) = 1
      count(1) = nlvls
      count(2) = nobs
      nf_status = nf_inq_varid(nf_fid,'levels',lev_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var levels'
      endif
      nf_status = nf_get_vara_real(nf_fid,lev_id,
     1 start,count,lvls_m)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading var levels'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif 

c     read _fillvalue for levels
      nf_status = nf_inq_attlen(nf_fid, lev_id,'_fillvalue',
     1                          levels_fill)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading levels _fillvalue attribute'
        levels_fill = 1.0e+38
      endif 

c     lw _fillvalue variable returns 1.4012985e-45...hard wire for now
      levels_fill = 9.9999997e+37

c     get varids for stationid, stationname, temperature,stationpressure
c       relhumidity, windspeed, winddir, observationtime

      nf_status = nf_inq_varid(nf_fid,'stationid',sta_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var stationid'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif

      nf_status = nf_inq_varid(nf_fid,'stationname',sn_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var stationname'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif

      nf_status = nf_inq_varid(nf_fid,'temperature',temp_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var temperature'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif

      nf_status = nf_inq_varid(nf_fid,'tempqcflag',tempqc_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'warning: could not find var tempqcflag'
!       nf_status = nf_close(nf_fid)
!       istatus = 0
!       return
        istat_tempqcflag = 0
      else
        istat_tempqcflag = 1
      endif

c     read _fillvalue for temperature
      nf_status = nf_inq_attlen(nf_fid, temp_id,'_fillvalue',
     1                          t_fill)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'use guess for filling temperature _fillvalue attribute'
        t_fill = 1.e+38
      endif 

      nf_status = nf_inq_varid(nf_fid,'stationpressure',sp_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'warning: could not find var stationpressure'
!       nf_status = nf_close(nf_fid)
!       istatus = 0
!       return
        istat_stationpressure = 0
      else
        istat_stationpressure = 1
      endif

      if(istat_stationpressure .eq. 1)then
        nf_status = nf_inq_varid(nf_fid,'prsqcflag',spqc_id)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'finding var spqcflag'
          print *, 'aborting read'
          nf_status = nf_close(nf_fid)
          istatus = 0
          return
        endif

c       read _fillvalue for stationpressure
        nf_status = nf_inq_attlen(nf_fid, sp_id,'_fillvalue',
     1                            sp_fill)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading stationpressure _fillvalue attribute'
          sp_fill = 3.4028e+38
        endif 
      endif

      if(.true.)then ! nimbus
          rh_var = 'relativehumidity'
      else ! rsa
          rh_var = 'relhumidity'
      endif

      nf_status = nf_inq_varid(nf_fid,trim(rh_var),rh_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'warning: could not find var ',rh_var     
!       nf_status = nf_close(nf_fid)
!       istatus = 0
!       return
        istat_relhumidity = 0
      else
        istat_relhumidity = 1
      endif

      if(istat_relhumidity .eq. 1)then
        nf_status = nf_inq_varid(nf_fid,'rhqcflag',rhqc_id)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'finding var rhqcflag'
!         print *, 'aborting read'
!         nf_status = nf_close(nf_fid)
!         istatus = 0
!         return
        endif

c       read _fillvalue for relhumidity
        nf_status = nf_inq_attlen(nf_fid, rh_id,'_fillvalue',
     1                            rh_fill)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading relhumidity _fillvalue attribute'
          rh_fill = 3.4028e+38
        endif 
      endif

      nf_status = nf_inq_varid(nf_fid,'windspeed',ws_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'error, could not find var windspeed',nf_noerr,nf_status
        print *, 'aborting read ',nf_fid
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif

      nf_status = nf_inq_varid(nf_fid,'wsqcflag',wsqc_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'could not find var wsqcflag'
!       print *, 'aborting read'
!       nf_status = nf_close(nf_fid)
!       istatus = 0
!       return
      endif

c     read _fillvalue for windspeed
      nf_status = nf_inq_attlen(nf_fid, ws_id,'_fillvalue',
     1                          ws_fill)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading windspeed _fillvalue attribute'
        ws_fill = 1.e+38
      endif 

      nf_status = nf_inq_varid(nf_fid,'winddir',wd_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'finding var winddir'
        print *, 'aborting read'
        nf_status = nf_close(nf_fid)
        istatus = 0
        return
      endif

      nf_status = nf_inq_varid(nf_fid,'wdqcflag',wdqc_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'could not find var wdqcflag'
!       print *, 'aborting read'
!       nf_status = nf_close(nf_fid)
!       istatus = 0
!       return
      endif

c     read _fillvalue for winddir
      nf_status = nf_inq_attlen(nf_fid, wd_id,'_fillvalue',
     1                          wd_fill)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'reading winddir _fillvalue attribute'
        wd_fill = 2147483647
      endif 

      nf_status = nf_inq_varid(nf_fid,'observationtime',ot_id)
      if(nf_status.ne.nf_noerr) then
        print *, nf_strerror(nf_status)
        print *,'could not find observationtime'
        print *, 'try for timeobs '

        nf_status = nf_inq_varid(nf_fid,'timeobs',ot_id)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'could not find timeobs'
          print *, 'aborting routine '
          nf_status = nf_close(nf_fid)
          istatus = 0
          return
        endif

      endif

      print *, 'lw nobs = ',nobs
      lev_set = 0
      do obno = 1, nobs

        if(obno .le. 50
     1     .or. (obno .eq. (obno/100 * 100) )     
     1     .or. (obno .ge. 1700 .and. obno .le. 1800)      
     1                                                 )then
            id = 1
        else
            id = 0
        endif

        if(id.eq.1)write(6,*)
        if(id.eq.1)write(6,*)' sa: obno,stalat,stalon,staelev'
     1                 ,obno,stalat(obno),stalon(obno),staelev(obno)       

        lev_set = 0

        index_1(1) = obno
        index_2(2) = obno
        start(2) = obno
        start(1) = 1
        count(2) = 1

        do lno = 1, 10
          if(id.eq.1)print *, 'lw level ',lno,'= ',lvls_m(lno,obno)
        enddo

c       read var observationtime(recnum) -> d_timeobs
        nf_status = nf_get_var1_double(nf_fid,ot_id,index_1,d_timeobs)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading var levels'
          print *, 'aborting read'
          nf_status = nf_close(nf_fid)
          istatus = 0
          return
        endif 

        i4_tim = int(d_timeobs + 315619200)
        call make_fnam_lp(i4_tim,a9time_ob(obno),istatus)

c       read var stationname(obno,stanamlen) -> stationname
        count(1) = sn_len 
        nf_status = nf_get_vara_text(nf_fid,sn_id,start,count,
     1                               stationname)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading var levels'
          print *, 'aborting read'
          nf_status = nf_close(nf_fid)
          istatus = 0
          return
        endif

c       truncate stname(obno)  
        call left_justify(stationname)
        call remove_blanks(stationname)
        stname(obno) = stationname(1:6)

c       read var stationid(recnum,provideridlen) -> c_staid 
        count(1) = pi_len 
        nf_status = nf_get_vara_text(nf_fid,sta_id,start,count,
     1                               c_staid)
        if(nf_status.ne.nf_noerr) then
          print *, nf_strerror(nf_status)
          print *,'reading var stationid'
          print *, 'aborting read'
          nf_status = nf_close(nf_fid)
          istatus = 0
          return
        endif
        call s_len(c_staid,len_sta)
        stname(obno) = '      '
        stname(obno)(1:len_sta) = c_staid(1:len_sta)          

c       needs doing
c       convert string to iwmostanum(maxobs) (cvt s to 0 and n to 1)
        call s_len(c_staid,lensta)
        do ic = 1, lensta 
          chr = ichar(c_staid(ic:ic))
          if ((chr.ge.48).and.(chr.le.57)) then 
            ichr = chr - 48
          else
            if (chr.eq.83) then
              ichr = 0
            else
              if (chr.eq.78) then
                ichr = 1
              else
                ichr = 0
              endif
            endif 
          endif
          if (ic.eq.1) ichr1 = ichr
          if (ic.eq.2) ichr2 = ichr
          if (ic.eq.3) ichr3 = ichr
          if (ic.eq.4) ichr4 = ichr
          if (ic.eq.5) ichr5 = ichr
          if (ic.eq.6) ichr6 = ichr
        enddo
        if (lensta .eq.1) wmoid(obno) = ichr1
        if (lensta .eq.2) 
     1    wmoid(obno) = ichr1*10 + ichr2
        if (lensta .eq.3) 
     1    wmoid(obno) = ichr1*100 + ichr2*10 + ichr3
        if (lensta .eq.4) 
     1    wmoid(obno) = ichr1*1000 + ichr2*100 + ichr3*10 + 
     1                  ichr4
        if (lensta .eq.5) 
     1    wmoid(obno) = ichr1*10000 + ichr2*1000 + ichr3*100 + 
     1                  ichr4*10 +ichr5
        if (lensta .eq.6) 
     1    wmoid(obno) = ichr1*100000 + ichr2*10000 + ichr3*1000 + 
     1                  ichr4*100 +ichr5*10 + ichr6

        if(id.eq.1)write(6,*) 'lw c_staid wmoid >',c_staid(1:lensta),
     1'<  >',wmoid(obno),'<'
        
        if(id.eq.1)print *,'lw obno nobs nlvls a9time staname'
     1            ,obno,nobs,nlvls,a9time_ob(obno),' ',stationname(1:20)       
     1            ,' ',stname(obno)

        lvl = 1
        do while ((lvl.le.nlvls).and.(lev_set.eq.0))

          if(id.eq.1)print *, 'lw levels_fill lvls_m ',
     1               levels_fill, lvls_m(lvl,obno)

          if(lvls_m(lvl,obno) .ne. levels_fill)then

            height_m(obno,lvl) = lvls_m(lvl,obno)
            index_2(1) = lvl
c           read stationpressure
            if(istat_stationpressure .eq. 1)then
              nf_status = nf_get_var1_real(nf_fid,sp_id,index_1
     1                                                 ,stationp)
              if(nf_status.ne.nf_noerr) then
                print *, nf_strerror(nf_status)
                print *,'error reading var stationpressure'
              endif 

c             read prsqcflag
              nf_status=nf_get_var1_real(nf_fid,spqc_id,index_1,prsqf)    
              if(nf_status.ne.nf_noerr) then
                print *, nf_strerror(nf_status)
                print *,'reading var prsqcflag'
              endif 

c             write prsqcflag
              prsqcflag(obno,lvl) = prsqf

              if(id.eq.1)write(6,*) 'lw o l stationp ',obno,lvl,stationp       

c             check stationpressure for _fillvalue
              if (stationp .eq. sp_fill) stationp = r_missing_data
              if (lvl.eq.1) then
                pressure_pa(obno,lvl) = stationp
              else
                pressure_pa(obno,lvl) = r_missing_data
              endif
            else
                pressure_pa(obno,lvl) = r_missing_data
            endif

c           read var temperature(recnum,lvl) -> temp
            nf_status = nf_get_var1_real(nf_fid,temp_id,index_2,temp_k)
            if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'reading var temperature'
            endif 

c           read var tempqcflag(recnum,lvl) -> tempqcflag
            if(istat_tempqcflag .eq. 1)then
                nf_status = nf_get_var1_real(nf_fid,tempqc_id,index_2,
     1                                       tempqf)
                if(nf_status.ne.nf_noerr) then
                  print *, nf_strerror(nf_status)
                  print *,'reading var tempqcflag'
                endif 

c               write tempqcflag
                tempqcflag(obno,lvl) = tempqf
            endif

c           convert temp_k to temp_c
            if ((temp_k .eq. -9999.).or.(temp_k .eq. t_fill)) then
              temp_c(obno,lvl) = r_missing_data
            else
              if ((temp_k .ge. 227.50).and.(temp_k .le. 328.15)) then
                temp_c(obno,lvl) = temp_k - 273.15
              else
                temp_c(obno,lvl) = r_missing_data
              endif
            endif

            if(id.eq.1)
     1      write(6,*) 'lw temp_k temp_c',temp_k,'   ',temp_c(obno,lvl)

c           read var relhumidity(recnum,level) -> rh
            nf_status = nf_get_var1_real(nf_fid,rh_id,index_2,rh_pct)
            if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'reading var relhumidity'
            endif 

c           read var rhqcflag(recnum,level) -> rhqcflag
            nf_status=nf_get_var1_real(nf_fid,rhqc_id,index_2,rhqf)
            if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'reading var rhqcflag'
            endif 

c           write rhqcflag
            rhqcflag(obno,lvl) = rhqf

c           convert rh to dewpoint
            if ((rh_pct .eq. -9999.).or.(rh_pct .eq. rh_fill)) then
              dewpoint_c(obno,lvl) = r_missing_data
            else
              if ((rh_pct .ge. 0).and.(rh_pct .le. 100)
     1            .and.(temp_c(obno,lvl).ne.r_missing_data)) then
                dewpoint_c(obno,lvl)=dwpt(temp_c(obno,lvl), rh_pct) ! celsius
              else
                dewpoint_c(obno,lvl) = r_missing_data
              endif
            endif

            if(id.eq.1)write(6,*)'lw rh_pct dpt_c ',rh_pct,'   '
     1                ,dewpoint_c(obno,lvl)       

c           read var windspeed(recnum,level) -> ws
            nf_status = nf_get_var1_real(nf_fid,ws_id,index_2,ws)
            if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'reading var windspeed'
            endif 

c           read var wsqcflag(recnum,level) -> wsqcflag
            nf_status=nf_get_var1_real(nf_fid,wsqc_id,index_2,wsqf)
            if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'reading var wsqcflag'
            endif 

c           write wsqcflag
            wsqcflag(obno,lvl) = wsqf

            if ((ws .eq. -9999.).or.(ws .eq. ws_fill)) then
              spd_mps(obno,lvl) = r_missing_data
            else
              if ((ws .ge. 0.).and.(ws .le. 150.)) then
                spd_mps(obno,lvl) = ws
              else
                spd_mps(obno,lvl) = r_missing_data
              endif
            endif

            if(id.eq.1)write(6,*) 'lw spd_mps      ',spd_mps(obno,lvl)

c           read var winddir(recnum,level) -> dd(lvl,obno)
            nf_status = nf_get_var1_real(nf_fid,wd_id,index_2,wd)
            if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'reading var winddir'
            endif 

c           read var wdqcflag(recnum,level) -> wdqcflag(lvl,obno)
            nf_status=nf_get_var1_real(nf_fid,wdqc_id,index_2,wdqf)
            if(nf_status.ne.nf_noerr) then
              print *, nf_strerror(nf_status)
              print *,'reading var wdqcflag'
            endif 

c           write wdqcflag
            wdqcflag(obno,lvl) = wdqf

            if ((wd .eq. -9999.).or.(wd .eq. wd_fill)) then
              dir_deg(obno,lvl) = r_missing_data
            else
              if ((ws .ge. 0).and.(ws .le. 360)) then
                dir_deg(obno,lvl) = wd
              else
                dir_deg(obno,lvl) = r_missing_data
              endif
            endif

            if(id.eq.1)write(6,*) 'lw dir_deg      ',dir_deg(obno,lvl)

          else
            nlvl(obno) = lvl - 1
            lev_set = 1
            if(id.eq.1)print *, 'lw lvl nlvl(obno) ',lvl, nlvl(obno)
          endif
          lvl = lvl + 1
        enddo
      enddo

      write(6,*) 'end of read_local_tower'

      return
      end

      subroutine remove_blanks(string)

      character*(*)string

      len1 = len(string)
      call s_len2(string,len2)

      len_loop = min(len2,len1-1)

      do i = 1,len_loop
          if(string(i:i) .eq. ' ')then ! blank found
!             ip1 = i+1
!             len1m1 = len1-1
              string(i:len1-1) = string(i+1:len1)
              string(len1:len1) = ' '
          endif
      enddo ! i

      return
      end


