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

module laps_parm

!==========================================================
!  this module defines all necessary constant and variables
!  for laps parameters.
!
!  history:
!         creation: yuanfu xie        3-2006
!==========================================================

   implicit none

   integer, parameter :: maxxobs = 100000        ! maximum obs limit

   include 'barnesob.inc'

   character*16 :: asctime

   integer :: n(3)                ! spatial dimensions
   integer :: timelen                ! time interval
   integer :: i4time                ! system time
   integer :: imissing                ! integer missing data
   integer :: nobs                ! number of observations
   integer :: n_tobs                ! number of temp observations
   integer :: nobs_point                ! number of point observations
   integer :: n_meso                ! number of mesonet data
   integer :: n_sao                ! number of sfc aviation obs
   ! integer :: n_pirep                ! number of pilot report obs
   ! integer :: max_pr                ! max profiles
   integer :: max_pr_lvls        ! max levels
   ! integer :: max_radars                ! maximum number of radars
   integer :: maxtobs                ! maximum number of temp obs
   integer :: istat_radar_vel

   logical :: l_raob, l_cdw, l_radial
   ! logical :: l_use_raob

   real :: rmissing                ! real missing data
   real :: dxy

   ! dynamic arrays:
   real, allocatable, dimension(:) :: pressr1d, height1d
   real, allocatable, dimension(:) :: rlat_radar, rlon_radar
   real, allocatable, dimension(:) :: rhgt_radar, n_grid_vel
   real, allocatable, dimension(:, :) :: lat, lon, topo
   real, allocatable, dimension(:, :, :) :: height3d, pressr3d
   real, allocatable, dimension(:, :, :) :: temptr3d, sphumd3d
   real, allocatable, dimension(:, :, :) :: u_wind3d, v_wind3d
   real, allocatable, dimension(:, :, :, :) :: grid_radar_vel

   real :: obs_temp(maxxobs, 12)
   type(barnesob) :: obs_point(maxxobs)

end module laps_parm
