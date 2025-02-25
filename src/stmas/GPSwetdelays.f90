subroutine gpswdelay

   use prmtrs_stmas

   implicit none

   character*200 :: fname, line
   character*9   :: ftime
   integer       :: length, istatus, i, j, k, m, i4
   integer, parameter :: lun = 33

   ! namelist data -- copied from dan's humid code

   integer covar_switch
   integer print_switch
   integer raob_switch
   integer raob_lookback
   real raob_radius
   integer endian  ! 1 = big, 0 = little , big default
   integer goes_switch
   integer cloud_switch
   integer cloud_d
   integer tiros_switch
   integer sounder_switch
   integer sat_skip
   integer gvap_switch
   integer ihop_flag
   integer time_diff         !time allowed for latency (sec)
   integer sfc_mix
   integer mod_4dda_1
   real mod_4dda_factor
   real t_ref, x, y
   integer gps_switch
   character*256 path_to_gvap12, path_to_gvap10, path_to_gps, path2covar

   namelist /moisture_switch_nl/ covar_switch, print_switch &
      , raob_switch, raob_lookback, endian &
      , raob_radius, goes_switch, cloud_switch, cloud_d &
      , tiros_switch, sounder_switch, sat_skip &
      , gvap_switch, ihop_flag, time_diff, gps_switch &
      , sfc_mix, mod_4dda_1, mod_4dda_factor &
      , t_ref, path_to_gvap12, path_to_gvap10, path_to_gps &
      , path2covar

!     fill namelist variables with default values -------------------------------------
!
!     set namelist parameters to defaults
   covar_switch = 0
   print_switch = 0
   cloud_switch = 1
   cloud_d = 1
   raob_switch = 0
   raob_lookback = 0
   endian = 1 ! big endian is default
   raob_radius = 45000.0  ! meters (45km)
   goes_switch = 0
   sounder_switch = 0
   tiros_switch = 0
   sat_skip = 0
   gvap_switch = 1
   ihop_flag = 0
   time_diff = 0
   gps_switch = 1
   sfc_mix = 0
   mod_4dda_1 = 0
   mod_4dda_factor = 0.02
   t_ref = -132.0
   path_to_gvap12 = ' '
   path_to_gvap10 = ' '
   path_to_gps = ' '
   path2covar = ' '

   call get_directory('static', fname, length)
   open (23, file=fname(1:length)//'moisture_switch.nl', status='old', err=6)
   read (23, moisture_switch_nl, end=6)
   close (23)

   ! generate filename:
   call make_fnam_lp(lapsi4t, ftime, istatus)

   ! open hmg file for recording gps data info and move to end of file: by yuanfu
   call open_lapsprd_file(lun, lapsi4t, 'hmg', istatus)
   if (istatus .eq. 1) then
      ! existing file:
1     read (lun, *, end=2) line
      goto 1
2     continue
   end if

   ! read gps data:
   call read_gps_obs(lun, path_to_gps, lapsi4t - icycle/2, lapsi4t + icycle/2, &
                     fcstgrd(1), fcstgrd(2), latitude, longitud, badsfcdt, &
                     gps_tpw, gps_wet, gps_err, gps_xyt, gps_elv, &
                     gps_tim, num_gps, max_gps, istatus)
   ! close hmg file:
   close (lun)

   j = 0
   k = 0
   m = 0
   do i = 1, num_gps
      if ((gps_tpw(i) .lt. 0.0) .and. (gps_wet(i) .ge. 0.0)) then
         print *, 'good wet delay without tpw'
         j = j + 1
      else if ((gps_tpw(i) .ge. 0.0) .and. (gps_wet(i) .lt. 0.0)) then
         print *, 'bad wet with good tpw'
         k = k + 1
      end if

      ! gps observation i4time from gps_tim that is seconds from 1/1/1970.
      ! i4time is seconds from 1/1/1960. we need to add gps_i4t
      i4 = gps_tim(i) + i4t_gps

      ! check gps data on the analysis time domain: by yuanfu
      if (i4 .ge. itime2(1) .and. i4 .le. itime2(2)) then
         m = m + 1
         gps_tpw(m) = gps_tpw(i)
         gps_wet(m) = gps_wet(i)
         gps_err(m) = gps_err(i)
         gps_xyt(1, m) = gps_xyt(1, i)
         gps_xyt(2, m) = gps_xyt(2, i)
         gps_elv(m) = gps_elv(i)
         gps_tim(m) = i4
      end if
   end do

   ! number of valid gps data:
   num_gps = m

   print *, 'total good wet bad tpw: ', j
   print *, 'total good tpw bad wet: ', k

   print *, 'actual available: ', m, ' out of total gps data: ', num_gps

   return

6  continue

end subroutine gpswdelay
