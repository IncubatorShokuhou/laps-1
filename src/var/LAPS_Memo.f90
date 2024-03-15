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

subroutine laps_allc

!==========================================================
!  this routine allocates memory for laps dynamic arrays.
!
!  history:
!         creation: yuanfu xie        3-2006
!==========================================================

   use laps_parm
   use mem_namelist

   implicit none

   integer :: status

   ! one dimension arrays:
   allocate (pressr1d(n(3)), height1d(n(3)), stat=status)
   if (status .ne. 0) then
      write (6, *) 'laps_allc: error to allocate memory for one-d arrays'

      ! deallocate all previously allocated arrays:
      call laps_remv
      stop
   end if
   allocate (rlat_radar(max_radars), rlon_radar(max_radars), &
             rhgt_radar(max_radars), n_grid_vel(max_radars))

   allocate (grid_radar_vel(n(1), n(2), n(3), max_radars), stat=status)
   if (status .ne. 0) then
      write (6, *) 'laps_allc: error to allocate memory for 4-d arrays'

      ! deallocate all previously allocated arrays:
      call laps_remv
      stop
   end if

   ! lat/lon and topography:
   allocate (lat(n(1), n(2)), lon(n(1), n(2)), topo(n(1), n(2)), &
             stat=status)
   if (status .ne. 0) then
      write (6, *) 'laps_allc: error to allocate memory for lat/lon/topo'

      ! deallocate all previously allocated arrays:
      call laps_remv
      stop
   end if

   ! 3d variables:
   allocate (height3d(n(1), n(2), n(3)), pressr3d(n(1), n(2), n(3)), &
             temptr3d(n(1), n(2), n(3)), sphumd3d(n(1), n(2), n(3)), &
             u_wind3d(n(1), n(2), n(3)), v_wind3d(n(1), n(2), n(3)), &
             stat=status)
   if (status .ne. 0) then
      write (6, *) 'laps_allc: error to allocate memory for 3d variables'

      ! deallocate all previously allocated arrays:
      call laps_remv
      stop
   end if

end subroutine laps_allc

subroutine laps_remv

!==========================================================
!  this routine deallocates memory for laps dynamic arrays.
!
!  history:
!         creation: yuanfu xie        3-2006
!==========================================================

   use laps_parm

   implicit none

   integer :: status

   ! one dimension arrays:
   deallocate (pressr1d, height1d, stat=status)

   deallocate (rlat_radar, rlon_radar, rhgt_radar, n_grid_vel, &
               grid_radar_vel, stat=status)

   ! lat/lon and topography:
   deallocate (lat, lon, topo, stat=status)

   ! 3d variables:
   deallocate (height3d, pressr3d, temptr3d, sphumd3d, &
               u_wind3d, v_wind3d, stat=status)

end subroutine laps_remv
