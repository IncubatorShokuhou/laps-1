           
      subroutine get_satsnd_afwa(i4time_sys,i4_satsnd_window
     1                          ,nx_l,ny_l
     1                          ,lun_in,filename,lun_out,istatus)

!     steve albers fsl    may 1999

      character*(*) filename

!.............................................................................

!     character*6 c6_a1acid
      character*9 a9_timeobs,a9_recpttime,a9time_ob 
      character*1000 c_line
      character*9 c_read
      character*5 c5_staid
      character*8 c8_obstype

      integer max_levels

      parameter (max_levels = 50)
      real rheight(max_levels)
      real rh(max_levels)
      real pressure(max_levels)
      real temp_k(max_levels)

      real lat_a(nx_l,ny_l)
      real lon_a(nx_l,ny_l)
      real topo_a(nx_l,ny_l)

!............................................................................

      open(lun_in,file=filename,status='old')

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_r_missing_data'
          return
      endif

      call get_domain_perimeter(nx_l,ny_l,'nest7grid',lat_a,lon_a, 
     1            topo_a,1.0,rnorth,south,east,west,istatus)
      if(istatus .ne. 1)then
          write(6,*)' error in get_laps_perimeter'
          return
      endif
  
      i = 0

      do while (.true.)

          read(lun_in,51,err=890,end=999)c_line
 51       format(a)

          read(c_line,101,err=890)     !    name             units & factor
     1         i_a1cycc,                       
     1         i_a1gwc,
     1         i_a1jul,                ! julian hour          hr since 673650000
     1         i_a1lat,                ! latitude             deg * 100
     1         i_a1lon,                ! longitude            deg * -100 
     1         i_a1nlvl

!    1         i_a1type,
!    1         i_a1min,                ! time-report-minutes
!    1         i_a1kind,
!    1         i_a1pla,
!    1         i_a1dval,
!    1         i_a1holm,
!    1         i_a1fltp,               ! temperature          kelvins * 10
!    1         i_a1wd,                 ! wind-direction       deg
!    1         i_a1wfls,               ! wind-speed           m/s * 10
!    1         c6_a1acid 

 101      format(6(i9,2x))

!         nlvls = i_a1nlvl
          nlvls = 16          ! hardwired as per afwa documentation

          iblk = 1
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)i_height
 102          format(i9)
              rheight(lvl) = i_height
          enddo ! lvl

          iblk = 2
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)i_rh
              rh(lvl) = i_rh
          enddo ! lvl

          iblk = 3
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)i_pressure
              pressure(lvl) = float(i_pressure) / 10.
          enddo ! lvl

          iblk = 4
          do lvl = 1,nlvls
              istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
              iend = istart+8
              c_read = c_line(istart:iend) 
              read(c_read,102,err=890)i_temp
              temp_k(lvl) = float(i_temp) / 10.
          enddo ! lvl

!         read hours & minutes
          iblk = 5
          lvl = 2
          istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
          iend = istart+8
          c_read = c_line(istart:iend) 
          read(c_read,112,err=890)i_hr,i_min
 112      format(5x,2i2)

!         read satellite id
          iblk = 5
          lvl = 3
          istart = 66 + (iblk-1)*nlvls*11 + (lvl-1)*11 + 1
          iend = istart+8
          c_read = c_line(istart:iend) 
          c5_staid = c_read(1:2)//c_read(4:6)

          i = i + 1

          write(6,*)
          write(6,*)' satsnd #',i

          stalat =  float(i_a1lat)/100.
          stalon = +float(i_a1lon)/100.

          write(6,2)stalat,stalon
 2        format(' lat, lon '/f8.3,f10.3)  

          if(stalat .le. rnorth .and. stalat .ge. south .and.
     1       stalon .ge. west   .and. stalon .le. east        )then       
              continue
          else ! outside lat/lon perimeter - reject
              write(6,*)' lat/lon - reject'       
!    1                 ,stalat,stalon
              goto 900
          endif

          i_hr_jul = i_a1jul - (i_a1jul/24)*24

          write(6,*)' i_hr, i_hr_jul, i_min ', i_hr, i_hr_jul, i_min       

          if(i_hr .ne. i_hr_jul)then
              write(6,*)' warning: i_hr discrepancy - reject '
     1                 ,i_hr,i_hr_jul
              goto900
          endif

          call afwa_julhr_i4time(i_a1jul,i_min,i4time_ob)

          call make_fnam_lp(i4time_ob,a9time_ob,istatus)
          if(istatus .ne. 1)goto900

          a9_recpttime = '         '

          i4_resid = abs(i4time_ob - i4time_sys)
          if(i4_resid .gt. i4_satsnd_window)then ! outside time window
              write(6,*)' time - reject '
     1           ,a9time_ob,i4_resid,i4_satsnd_window
              goto 900        
          endif

          write(6,1)a9time_ob
 1        format(' time:'/1x,a9) 

          staelev = 0.

          c8_obstype = 'satsnd'

          iwmostanum = 0

          write(6,511,err=990)
     1             iwmostanum,nlvls,stalat
     1            ,stalon,staelev,c5_staid,a9time_ob,c8_obstype
          write(lun_out,511,err=990)
     1             iwmostanum,nlvls,stalat
     1            ,stalon,staelev,c5_staid,a9time_ob,c8_obstype

  511     format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

          dir = r_missing_data
          spd = r_missing_data

          do lvl = 1,nlvls
              if(abs(rheight(lvl)) .gt. 100000.)then
                  rheight(lvl) = r_missing_data
              endif

              if(abs(pressure(lvl)) .gt. 100000.)then
                  pressure(lvl) = r_missing_data
              endif

              if(temp_k(lvl) .ge. 150. .and. temp_k(lvl) .le. 400.)then       
                  temp_c = temp_k(lvl) - 273.15
              else
                  temp_c = r_missing_data
              endif

!             convert rh(lvl) to dewpoint
              if(rh(lvl) .gt. 0. .and. rh(lvl) .le. 100.
     1                           .and. temp_c .ne. r_missing_data)then       
                  dewpoint_c = dwpt(temp_c,rh(lvl))
              else
                  dewpoint_c = r_missing_data
              endif

              if(i .le. 100)write(6,*)rheight(lvl),pressure(lvl)
     1              ,temp_c
     1              ,dewpoint_c
     1              ,lvl

              write(lun_out,*)rheight(lvl),pressure(lvl)
     1              ,temp_c
     1              ,dewpoint_c
     1              ,dir,spd

          enddo ! lvl

          go to 900

 890      write(6,*)' warning (get_satsnd_afwa): read/write error'

 900  enddo ! read line of afwa file

!............................................................................

 990  write(6,*)' error in ingest_satsnd_af'
      istatus=0
      return
  
 999  write(6,*)' end of afwa file detected'
      istatus = 1
      return
      end

