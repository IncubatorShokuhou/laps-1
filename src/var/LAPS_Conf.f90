!dis    forecast systems laboratory
!dis    noaa/oar/erl/fsl
!dis    325 broadway
!dis    boulder, co     80303
!dis
!dis    forecast research division
!dis    local analysis and prediction branch
!dis    laps
!dis
!dis    this software and its documentation are in the public domain and
!dis    are furnished "as is."  the united states government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  they assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  all modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  if significant modifications or enhancements
!dis    are made to this software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

subroutine laps_conf

!==========================================================
!  this routine initializes laps configuration parameters:
!  i4time, dimensions and so on.
!
!  history:
!         creation: yuanfu xie        3-2006
!==========================================================

   use laps_parm
   use mem_namelist

   implicit none

   ! local variables:
   character :: a9time*9, fnm*9, hr*2, mins*2, jday*5
   integer :: status, i4time_sys
!  integer ::thresh_2_radarobs_lvl_unfltrd, &
!              thresh_4_radarobs_lvl_unfltrd, &
!              thresh_9_radarobs_lvl_unfltrd, i4time_sys
!  real :: weight_bkg_const_wind,weight_radar,rms_thresh_wind
!  integer :: max_obs

   ! wind parameters:
   character*150 :: static_dir, filename
   integer       :: len_dir

!  namelist /wind_nl/ l_raob, l_cdw, l_radial, &
!                     thresh_2_radarobs_lvl_unfltrd, &
!                     thresh_4_radarobs_lvl_unfltrd, &
!                     thresh_9_radarobs_lvl_unfltrd, &
!                     weight_bkg_const_wind, &
!                     weight_radar, &
!                     rms_thresh_wind, &
!                     max_pr,max_pr_lvls,max_obs

   ! spatial dimensions:
   call get_grid_dim_xy(n(1), n(2), status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error in horizontal dimensions'
      stop
   end if
   call get_laps_dimensions(n(3), status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error in vertical dimension'
      stop
   end if
   print *, 'n = ', n

   ! system time:
   call get_systime(i4time, a9time, status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error in system times'
      stop
   else
      call get_directory('log', filename, len_dir)
      filename = filename(1:len_dir)//'i4time.txt'
      ! open(10,file='i4time.txt')
      open (10, file=filename(1:len_dir + 10))
      write (10, *) i4time
      close (10)
   end if
   call get_systime_all(i4time_sys, fnm, hr, mins, asctime, jday, status)
   if (i4time .ne. i4time_sys) then
      print *, 'laps_conf: error in reading background at wrong time'
      stop
   end if
   call get_laps_cycle_time(timelen, status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error in laps cycle time'
      stop
   end if
   print *, 'time: ', i4time, a9time, timelen

   ! missing data:
   call get_r_missing_data(rmissing, status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error obtaining real missing_data'
      stop
   end if
   call get_i2_missing_data(imissing, status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error obtaining integer missing_data'
      stop
   end if

   ! get wind parameters:
   call get_directory('static', static_dir, len_dir)
   filename = static_dir(1:len_dir)//'wind.nl'
   call read_namelist_laps('wind', filename)
   ! get temperature parameters:
   filename = static_dir(1:len_dir)//'temp.nl'
   call read_namelist_laps('temp_anal', filename)

   ! radar info:
   call get_max_radars(max_radars, status)
   if (status .ne. 1) then
      write (6, *) 'error obtaining max_radars'
      stop
   end if

   ! meso, sao and pirep:
   call get_meso_sao_pirep(n_meso, n_sao, n_pirep, status)
   if (status .ne. 1) then
      write (6, *) 'error obtaining n_meso, n_sao, and n_pirep'
      stop
   end if

   ! allocate laps dynamic arrays:
   call laps_allc

   ! read lat/lon and topo data:
   call read_static_grid(n(1), n(2), 'lat', lat, status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error get laps lat'
      stop
   end if
   call read_static_grid(n(1), n(2), 'lon', lon, status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error get laps lon'
      stop
   end if
   call read_static_grid(n(1), n(2), 'avg', topo, status)
   if (status .ne. 1) then
      write (6, *) 'laps_conf: error get laps topography'
      stop
   end if
   print *, 'lat/lon/avg: ', lat(1, 1), lon(1, 1), topo(1, 1)

   ! grid spacing:
   call get_grid_spacing_actual(lat(n(1)/2 + 1, n(2)/2 + 1), &
                                lon(n(1)/2 + 1, n(2)/2 + 1), &
                                dxy, status)
   print *, 'grid spacing: ', dxy

end subroutine laps_conf
