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

subroutine laps_bkgd

!==========================================================
!  this routine retrieves background fields by laps ingest.
!  the 3 d background fields are:
!        height, pressure, temperature, specific humidity,
!       u wind, v wind.
!
!  history:
!         creation: yuanfu xie        3-2006
!==========================================================

   use laps_parm

   implicit none

   ! local variables:
   character*3 :: varname
   integer :: status

   ! height:
   varname = 'ht'
   call get_modelfg_3d(i4time, varname, n(1), n(2), n(3), &
                       height3d, status)
   if (status .ne. 1) then
      write (6, *) 'laps_pars: error retrieving height'
      call laps_remv
      stop
   end if

   ! pressure:
   call get_pres_3d(i4time, n(1), n(2), n(3), pressr3d, status)
   if (status .ne. 1) then
      write (6, *) 'laps_pars: error retrieving pressure'
      call laps_remv
      stop
   end if
   call get_pres_1d(i4time, n(3), pressr1d, status)

   ! temperature:
   varname = 't3'
   call get_modelfg_3d(i4time, varname, n(1), n(2), n(3), &
                       temptr3d, status)
   if (status .ne. 1) then
      write (6, *) 'laps_pars: error retrieving temperature'
      call laps_remv
      stop
   end if

   ! specific:
   ! varname = 'rh'
   ! call get_modelfg_3d(i4time,varname,n(1),n(2),n(3), &
   !                sphumd3d,status)
   ! if (status .ne. 1) then
   !  write(6,*) 'laps_pars: error retrieving relative humidity'
   !  call laps_remv
   !  stop
   !endif
   varname = 'sh'
   call get_modelfg_3d(i4time, varname, n(1), n(2), n(3), &
                       sphumd3d, status)
   if (status .ne. 1) then
      write (6, *) 'laps_pars: error retrieving specific humidity'
      call laps_remv
      stop
   end if

   ! u wind:
   varname = 'u3'
   call get_modelfg_3d(i4time, varname, n(1), n(2), n(3), &
                       u_wind3d, status)
   if (status .ne. 1) then
      write (6, *) 'laps_pars: error retrieving u wind'
      call laps_remv
      stop
   end if

   ! v wind:
   varname = 'v3'
   call get_modelfg_3d(i4time, varname, n(1), n(2), n(3), &
                       v_wind3d, status)
   if (status .ne. 1) then
      write (6, *) 'laps_pars: error retrieving v wind'
      call laps_remv
      stop
   end if

end subroutine laps_bkgd
