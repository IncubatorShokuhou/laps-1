           
      subroutine write_snd(lun_out                         ! i
     1                    ,maxsnd,maxlvl,nsnd              ! i
     1                    ,iwmostanum                      ! i
     1                    ,stalat,stalon,staelev           ! i
     1                    ,c5_staid,a9time_ob,c8_obstype   ! i
     1                    ,nlvl                            ! i
     1                    ,height_m                        ! i
     1                    ,pressure_mb                     ! i
     1                    ,temp_c                          ! i
     1                    ,dewpoint_c                      ! i
     1                    ,dir_deg                         ! i
     1                    ,spd_mps                         ! i
     1                    ,istatus)                        ! o

!     steve albers fsl    2001

!     write routine for 'snd' file

!     for missing data values, 'r_missing_data' should be passed in 

!.............................................................................

      integer iwmostanum(maxsnd),nlvl(maxsnd)
      real stalat(maxsnd,maxlvl),stalon(maxsnd,maxlvl),staelev(maxsnd)       
      character c5_staid(maxsnd)*5,a9time_ob(maxsnd,maxlvl)*9
     1         ,c8_obstype(maxsnd)*8,c_line*200

      character*9 a9_time

      character*5 c5_sta

      real height_m(maxsnd,maxlvl)
      real pressure_mb(maxsnd,maxlvl)
      real temp_c(maxsnd,maxlvl)
      real dewpoint_c(maxsnd,maxlvl)
      real dir_deg(maxsnd,maxlvl)
      real spd_mps(maxsnd,maxlvl)

      integer nsnd_total
      data nsnd_total/0/
      save nsnd_total,i4time_sys

      integer max_lvls
      parameter (max_lvls=200)
      real liquid_a(max_lvls)

      common /write_snd_data/ cloud_base_temp,cloud_integrated_liquid      
     1                       ,liquid_a

!............................................................................

!     get raob time window
      call get_windob_time_window('raob',i4_wind_ob,istatus)
      if(istatus .ne. 1)goto 990

      call get_tempob_time_window('raob',i4_temp_ob,istatus)
      if(istatus .ne. 1)goto 990

      i4_raob_window = max(i4_wind_ob,i4_temp_ob)

      if(nsnd_total .eq. 0)then ! initialize variables
!         get systime
          call getenv('laps_a9time',a9_time)
          call s_len(a9_time,ilen)

          if(ilen .eq. 9)then
!           write(6,*)' systime (from env) = ',a9_time
            call i4time_fname_lp(a9_time,i4time_sys,istatus)
          else
            call get_systime(i4time_sys,a9_time,istatus)
            if(istatus .ne. 1)go to 990
!           write(6,*)' systime = ',a9_time
          endif
      endif

      do isnd = 1,nsnd

        nsnd_total = nsnd_total + 1

!       reject observation times outside time window        
        call i4time_fname_lp(a9time_ob(isnd,1),i4time_raob,status)
        if(abs(i4time_raob - i4time_sys) .gt. i4_raob_window)then
            write(6,*)a9time_ob(isnd,1),
     1                ' is outside time window with sounding # ',isnd
            goto 900
        endif

!       default qc for 'iwmostanum' and 'c5_staid'
        call check_nan(iwmostanum(isnd),istatus)
        if(istatus .ne. 1)then
            write(6,*)' nan detected for wmoid, set to sounding # ',isnd
            iwmo_out = isnd
        else
            iwmo_out = iwmostanum(isnd)
        endif

        if(iwmostanum(isnd) .lt. 0 .or. 
     1     iwmostanum(isnd) .gt. 99999)then
            write(6,*)' wmoid out of bounds, set to sounding # ',isnd
            iwmo_out = isnd
        endif

        call s_len(c5_staid(isnd),len_sta)
        if(len_sta .gt. 0)then
            c5_sta = c5_staid(isnd)

        else 
            write(6,*)' filling in blank staid with wmoid # '
     1               ,isnd,iwmo_out

            write(c5_sta,101)iwmo_out
 101        format(i5.5)

        endif

!       write sounding header

        if(c8_obstype(isnd) .ne. 'radiomtr')then
            write(6,511,err=990)
     1             iwmo_out,nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_sta,a9time_ob(isnd,1),c8_obstype(isnd)

            write(lun_out,511,err=990)
     1             iwmo_out,nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_sta,a9time_ob(isnd,1),c8_obstype(isnd)

        else
            write(6,511,err=990)
     1             iwmo_out,nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_sta,a9time_ob(isnd,1),c8_obstype(isnd)
     1            ,cloud_base_temp,cloud_integrated_liquid

            write(lun_out,511,err=990)
     1             iwmo_out,nlvl(isnd)
     1            ,stalat(isnd,1),stalon(isnd,1),staelev(isnd)
     1            ,c5_sta,a9time_ob(isnd,1),c8_obstype(isnd)
     1            ,cloud_base_temp,cloud_integrated_liquid

        endif

  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8,f8.1,f8.5)

        do lvl = 1,nlvl(isnd)

!         write sounding level (the character array helps keep everything
!                               in one line when using free format)

          if(c8_obstype(isnd) .ne. 'radiomtr')then
            write(c_line,*)height_m(isnd,lvl)," "
     1              ,pressure_mb(isnd,lvl)," "
     1              ,temp_c(isnd,lvl)," "
     1              ,dewpoint_c(isnd,lvl)," "
     1              ,dir_deg(isnd,lvl)," "
     1              ,spd_mps(isnd,lvl)," "
     1              ,a9time_ob(isnd,lvl)," "
     1              ,stalat(isnd,lvl)," "
     1              ,stalon(isnd,lvl)
          else
            write(c_line,*)height_m(isnd,lvl)," "
     1              ,pressure_mb(isnd,lvl)," "
     1              ,temp_c(isnd,lvl)," "
     1              ,dewpoint_c(isnd,lvl)," "
     1              ,dir_deg(isnd,lvl)," "
     1              ,spd_mps(isnd,lvl)," "
     1              ,a9time_ob(isnd,lvl)," "
     1              ,stalat(isnd,lvl)," "
     1              ,stalon(isnd,lvl)," "
     1              ,liquid_a(lvl)
          endif
          call s_len2(c_line,len_line)
          write(lun_out,521)c_line(1:len_line)
  521     format(a)

          if(nsnd_total .le. 100)then
              write(6,521)c_line(1:len_line)
          endif

        enddo ! lvl

 900    continue

      enddo ! isnd

      go to 999

 990  write(6,*)' error in write_snd'
      istatus=0
      return

 999  istatus = 1
      return
      end

