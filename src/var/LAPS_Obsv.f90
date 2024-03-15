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

subroutine laps_obsv

!==========================================================
!  this routine retrieves observation through laps ingest.
!
!  history:
!         creation: yuanfu xie        3-2006        wind obs;
!        modified: yuanfu xie        5-2006        add temp obs;
!==========================================================

   use laps_parm
   use mem_namelist

   implicit none

   ! local variables:
   integer :: useable_radar, istatus_remap_pro, status
   integer :: init_timer, ntmin, ntmax
   real :: u(n(1), n(2), n(3)), v(n(1), n(2), n(3))
   real :: ulaps(n(1), n(2), n(3)), vlaps(n(1), n(2), n(3))
   real :: wt(n(1), n(2), n(3)), weight_prof

   integer :: max_snd_grd, max_snd_lvl
   logical :: l_adj_hgt
   real :: bg_weight

   ! profiler data:
   ntmin = -1
   ntmax = +1

   status = init_timer()

   call get_wind_3d_obs(n(1), n(2), n(3), rmissing, imissing, &
                        i4time, height3d, height1d, max_pr, max_pr_lvls, &
                        weight_prof, l_raob, l_cdw, n_sao, n_pirep, lat, lon, &
                        ntmin, ntmax, u, v, ulaps, vlaps, wt, maxxobs, obs_point, &
                        nobs_point, rlat_radar, rlon_radar, rhgt_radar, &
                        useable_radar, n_grid_vel, istatus_remap_pro, status)

   ! temperature parameters:
   ! call get_temp_parms(l_raob,l_use_raob,l_adj_hgt,       &
   !                    bg_weight,rms_thresh,              &
   !                    pres_mix_thresh,max_snd_grd,max_snd_lvl,maxtobs,   &
   !                    status)

   ! assign temperature parameters using the new temperature module:
   l_raob = l_read_raob_t
   maxtobs = max_obs
   call get_temp_3d_obs(max_snd_grid, max_snd_levels, max_obs)

end subroutine laps_obsv

subroutine get_temp_3d_obs(max_snd_grd, max_snd_lvl, max_obs)

!==========================================================
!  this routine retrieves 3d temperature observation data.
!
!  history:
!        creation: yuanfu xie        5-2006.
!==========================================================

   use laps_parm

   implicit none

   integer, intent(in) :: max_snd_grd, max_snd_lvl, max_obs

   integer :: i_obstype
   integer :: status

   include 'tempobs.inc'

   ! sonde temperature data:
   call get_temp_snd(max_snd_grd, max_snd_lvl, temp_obs, max_obs, status)

   ! acar temperature data:
   call get_temp_acar(temp_obs, max_obs, status)

   ! copy to the output array:
   obs_temp(1:max_obs, 1:12) = temp_obs(1:max_obs, 1:12)

end subroutine get_temp_3d_obs

subroutine get_temp_snd(max_snd_grd, max_snd_lvl, temp_obs, maxaobs, error)

!==========================================================
!  this routine retrieves temperature obs from sonde data.
!
!  history:
!        adapted from insert_tobs: yuanfu xie        5-2006.
!==========================================================

   use laps_parm

   implicit none

   integer, intent(in) :: max_snd_grd, max_snd_lvl, maxaobs
   integer, intent(out) :: error
   real, intent(out) :: temp_obs(maxaobs, 12)        ! tempobs.inc

   integer :: i4_window_raob
   integer :: status, n_rass, n_snde, n_tsnd
   real :: lattsnd(max_snd_grd), lontsnd(max_snd_grd)
   real :: tsnd(max_snd_grd, n(3))
   real :: inst_err_tsnd(max_snd_grd)
   real :: bias_tsnd(max_snd_grd, n(3)), bias_htlow(max_snd_grd)
   character*5 c5name(max_snd_grd)
   character*8 c8obstype(max_snd_grd)
   logical :: l_struct

   integer :: k, isnd, n_qc_snd
   integer :: igrid(max_snd_grd), jgrid(max_snd_grd)
   real :: ri, rj, p_pa, sh, tvir, tamb, devirt_sh
   logical :: l_string_contains, l_qc

   i4_window_raob = 0                 ! same as insert_tobs
   l_struct = .true.

   ! read sonde temperature:
   call read_tsnd(i4time, height3d, temptr3d, sphumd3d, &
                  pressr3d, lattsnd, lontsnd, lat, lon, &
                  max_snd_grd, max_snd_lvl, tsnd, inst_err_tsnd, c5name, &
                  c8obstype, l_raob, l_struct, &
                  i4_window_raob, bias_htlow, &
                  n_rass, n_snde, n_tsnd, timelen, n(1), &
                  n(2), n(3), rmissing, status)

   ! qc temperature observation data:
   n_qc_snd = 0
   n_tobs = 0
   do isnd = 1, n_tsnd

      ! location on the grid:
      call latlon_to_rlapsgrid(lattsnd(isnd), lontsnd(isnd), &
                               lat, lon, n(1), n(2), ri, rj, status)

      ! count obs within the domain:
      if (status .eq. 1) then
         igrid(isnd) = nint(ri)
         jgrid(isnd) = nint(rj)

         ! find the sonde data:
         do k = 1, n(3)
            if (tsnd(isnd, k) .ne. rmissing) then
               ! rass observes the virtue temperature:
               if (l_string_contains(c8obstype(isnd), 'rass', status)) then
                  ! convert from virtual temperature to temperature
                  tvir = tsnd(isnd, k)
                  sh = sphumd3d(igrid(isnd), jgrid(isnd), k)
                  p_pa = pressr3d(igrid(isnd), jgrid(isnd), k)
                  tamb = devirt_sh(tvir, sh, p_pa)
               else
                  sh = 0.
                  tamb = tsnd(isnd, k)
               end if
               ! save bias:
               bias_tsnd(isnd, k) = tamb - &
                                    temptr3d(igrid(isnd), jgrid(isnd), k)

               ! hard qc:
               if (abs(bias_tsnd(isnd, k)) .gt. 10.) then
                  l_qc = .true.
                  write (6, *) ' abs(temp - first guess) > 10., temp not used'
               end if
            end if
         end do

         ! good sonde obs:
         if (.not. l_qc) then
            n_qc_snd = n_qc_snd + 1
            ! count the observations:
            do k = 1, n(3)
               if (bias_tsnd(isnd, k) .ne. rmissing) then
                  n_tobs = n_tobs + 1

                  if (n_tobs .gt. maxaobs) then
                     write (6, *) 'too many temperature sonde obs'
                     error = 0
                     return
                  end if

                  ! pass the obs to temperature obs array:
                  ! note: referring to tempobs.inc under include
                  temp_obs(n_tobs, 1) = igrid(isnd)
                  temp_obs(n_tobs, 2) = jgrid(isnd)
                  temp_obs(n_tobs, 3) = k
                  temp_obs(n_tobs, 4) = igrid(isnd)
                  temp_obs(n_tobs, 5) = jgrid(isnd)
                  temp_obs(n_tobs, 6) = k
                  temp_obs(n_tobs, 8) = &
                     temptr3d(igrid(isnd), jgrid(isnd), k) &
                     + bias_tsnd(isnd, k)
                  temp_obs(n_tobs, 9) = bias_tsnd(isnd, k)
                  temp_obs(n_tobs, 10) = &
                     1.0/inst_err_tsnd(isnd)**2
                  temp_obs(n_tobs, 11) = inst_err_tsnd(isnd)
               end if
            end do
         end if
      end if

   end do

   write (6, *) '# of sonde temperature obs: ', n_tobs, n_qc_snd

end subroutine get_temp_snd

subroutine get_temp_acar(temp_obs, maxaobs, status)

!==========================================================
!  this routine reads in the temperature obs from acar.
!
!  history:
!        adapted from insert_tobs: yuanfu xie        5-2006.
!==========================================================

   use laps_parm
   use mem_namelist

   implicit none

   integer, intent(in) :: maxaobs
   integer, intent(out) :: status
   real, intent(inout) :: temp_obs(maxaobs, 12)

   integer :: n_good_acars

   call rd_acars_t(i4time, height3d, temptr3d &
                   , n_pirep, n_good_acars, 'pin' &
                   , n(1), n(2), n(3), lat, lon, rmissing &
                   , temp_obs, maxaobs, n_tobs, status)

end subroutine get_temp_acar

