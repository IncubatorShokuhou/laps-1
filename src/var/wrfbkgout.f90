!dis
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis
!dis
!dis
!dis
!dis

subroutine wrfbkgout(times, imax, jmax, kmax, ptop, znu, znw, dxy, &
                     mapfac, lat, lon, dam, pdam, t, geo, sh, u, v, topo)

!==========================================================
!  this routine writes the background fields into wrf_inout
!  in netcdf format to be used by gsi.
!
!  history: mar. 2006 by yuanfu xie.
!           sep. 2006 by yuanfu xie: compute mub following
!                                     wrf convention (see
!                                    arw description pp.36)
!                oct. 2007 by yuanfu xie: add true ph and phb.
!==========================================================

   implicit none

   include 'netcdf.inc'

   character*19, intent(in) :: times
   integer, intent(in) :: imax, jmax, kmax          ! 3d array dimensions
   real, intent(in) :: lat(imax, jmax)        ! latitude
   real, intent(in) :: lon(imax, jmax)        ! longitude
   real, intent(in) :: ptop                ! pressure top
   real, intent(in) :: znu(kmax - 1), znw(kmax)        ! staggered eta,eta
   real, intent(in) :: dxy                        ! grid spacing
   real, intent(in) :: mapfac(imax, jmax)        ! map factor projection
   real, intent(in) :: dam(imax, jmax)        ! dry air mass (column)
   real, intent(in) :: pdam(imax, jmax)        ! perturbation
   real, intent(in) :: t(imax, jmax, kmax)   ! temperature
   real, intent(in) :: geo(imax, jmax, kmax)   ! temperature
   real, intent(in) :: sh(imax, jmax, kmax)        ! specific humidity
   real, intent(in) :: u(imax, jmax, kmax)   ! u
   real, intent(in) :: v(imax, jmax, kmax)   ! v
   real, intent(in) :: topo(imax, jmax)                ! topography

   ! local variables:
   character*10 :: empty
   character*125 :: filename                        ! nc filename
   integer :: namelen                                ! nc filename length
   integer :: ncid                                ! file id
   integer :: tmid, uid, vid, tid, muid, mubid        ! var ids
   integer :: qid, mapid, ptid, znuid, znwid        ! var ids
   integer :: latid, lonid, rxid, ryid                ! var ids
   integer :: phid, phbid, lmkid, iceid, sstid, vgid        ! var ids
   integer :: slid, vfid, snwid, u10id, v10id        ! var ids
   integer :: smsid, tslbid, tskid                ! var ids
   integer :: start(4), count(4)                ! netcdf start/count
   integer :: i, j, k, ierr                                ! error flag
   integer :: time, dlen, ndim(4), ndm1(4), btop, we, sn, nd(4), scal, nsol
   integer :: itnc(imax - 1, jmax - 1, kmax - 1)
   real :: unc(imax, jmax - 1, kmax - 1), tmp(imax, jmax, kmax)
   real :: vnc(imax - 1, jmax, kmax - 1), tnc(imax - 1, jmax - 1, kmax - 1)
   real :: zmp(kmax), tmp1(imax*jmax*kmax)
!! ltsung add new variables for subroutine staggerxy
   real :: u_out(imax, jmax - 1, kmax), v_out(imax - 1, jmax, kmax), &
           t_out(imax - 1, jmax - 1, kmax), mu_out(imax - 1, jmax - 1, 1), &
           mub_out(imax - 1, jmax - 1, 1), qvapor_out(imax - 1, jmax - 1, kmax), &
           mapfac_m_out(imax - 1, jmax - 1, 1), xlat_out(imax - 1, jmax - 1, 1), &
           xlong_out(imax - 1, jmax - 1, 1)
   real :: unc1(imax, jmax - 1, kmax - 1), vnc1(imax - 1, jmax, kmax - 1), &
           tnc1(imax - 1, jmax - 1, kmax - 1), qnc1(imax - 1, jmax - 1, kmax - 1)
   real :: a0, p0, t0, terrain

   empty = ' '

   ! create the netcdf file:
   call get_directory('log', filename, namelen)        ! yuanfu: use laps;
   filename = filename(1:namelen)//'wrf_inout.nc'
   ! check if the directory path is too long:
   if (namelen + 12 .gt. 125) then
      print *, 'error: need longer name for filename'
      stop
   end if
   ! ncid = nccre('wrf_inout.nc',ncnoclob,ierr)
   ncid = nccre(filename(1:namelen + 12), ncnoclob, ierr)

   ! global attributes:
   call ncaptc(ncid, ncglobal, 'title', ncchar, 22, &
               'laps background ingest', ierr)
   call ncaptc(ncid, ncglobal, 'start_date', ncchar, &
               len(times), times, ierr)
   call ncapt(ncid, ncglobal, 'west-east_grid_dimension', &
              nclong, 1, imax, ierr)
   call ncapt(ncid, ncglobal, 'south-north_grid_dimension', &
              nclong, 1, jmax, ierr)
   call ncapt(ncid, ncglobal, 'bottom-top_grid_dimension', &
              nclong, 1, kmax, ierr)
   call ncaptc(ncid, ncglobal, 'gridtype', ncchar, 1, 'c', ierr)
   call ncapt(ncid, ncglobal, 'west-east_patch_start_unstag', &
              nclong, 1, 1, ierr)
   call ncapt(ncid, ncglobal, 'west-east_patch_end_unstag', &
              nclong, 1, imax - 1, ierr)
   call ncapt(ncid, ncglobal, 'west-east_patch_start_stag', &
              nclong, 1, 1, ierr)
   call ncapt(ncid, ncglobal, 'west-east_patch_end_stag', &
              nclong, 1, imax, ierr)
   call ncapt(ncid, ncglobal, 'south-north_patch_start_unstag', &
              nclong, 1, 1, ierr)
   call ncapt(ncid, ncglobal, 'south-north_patch_end_unstag', &
              nclong, 1, jmax - 1, ierr)
   call ncapt(ncid, ncglobal, 'south-north_patch_start_stag', &
              nclong, 1, 1, ierr)
   call ncapt(ncid, ncglobal, 'south-north_patch_end_stag', &
              nclong, 1, jmax, ierr)
   call ncapt(ncid, ncglobal, 'bottom-top_patch_start_unstag', &
              nclong, 1, 1, ierr)
   call ncapt(ncid, ncglobal, 'bottom-top_patch_end_unstag', &
              nclong, 1, kmax - 1, ierr)
   call ncapt(ncid, ncglobal, 'bottom-top_patch_start_stag', &
              nclong, 1, 1, ierr)
   call ncapt(ncid, ncglobal, 'bottom-top_patch_end_stag', &
              nclong, 1, kmax, ierr)
   call ncapt(ncid, ncglobal, 'dx', ncfloat, 1, dxy, ierr)
   call ncapt(ncid, ncglobal, 'dy', ncfloat, 1, dxy, ierr)

   ! create dimensions:
   time = ncddef(ncid, 'time', ncunlim, ierr)
   dlen = ncddef(ncid, 'datestrlen', 19, ierr)
   ndim(1) = ncddef(ncid, 'west_east_stag', imax, ierr)
   ndim(2) = ncddef(ncid, 'south_north_stag', jmax, ierr)
   ndim(3) = ncddef(ncid, 'bottom_top_stag', kmax, ierr)
   ndm1(1) = ncddef(ncid, 'west_east', imax - 1, ierr)
   ndm1(2) = ncddef(ncid, 'south_north', jmax - 1, ierr)
   ndm1(3) = ncddef(ncid, 'bottom_top', kmax - 1, ierr)
   nsol = ncddef(ncid, 'soil_layers_stag', 4, ierr)
   scal = ncddef(ncid, 'ext_scalar', 1, ierr)

   ! create variables:
   ! times:
   nd(1) = dlen
   nd(2) = time
   tmid = ncvdef(ncid, 'times', ncchar, 2, nd, ierr)
   ! u:
   nd(1) = ndim(1)
   nd(2:3) = ndm1(2:3)
   nd(4) = time
   uid = ncvdef(ncid, 'u', ncfloat, 4, nd, ierr)
   call ncapt(ncid, uid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, uid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, uid, 'description', ncchar, 16, &
               'x-wind component', ierr)
   call ncaptc(ncid, uid, 'units', ncchar, 7, 'm s{-1}', ierr)
   call ncaptc(ncid, uid, 'stagger', ncchar, 1, 'x', ierr)
   ! v:
   nd(1) = ndm1(1)
   nd(2) = ndim(2)
   nd(3) = ndm1(3)
   nd(4) = time
   vid = ncvdef(ncid, 'v', ncfloat, 4, nd, ierr)
   call ncapt(ncid, vid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, vid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, vid, 'description', ncchar, 16, &
               'y-wind component', ierr)
   call ncaptc(ncid, vid, 'units', ncchar, 7, 'm s{-1}', ierr)
   call ncaptc(ncid, vid, 'stagger', ncchar, 1, 'y', ierr)
   ! t:
   nd(1:3) = ndm1(1:3)
   nd(4) = time
   tid = ncvdef(ncid, 't', ncfloat, 4, nd, ierr)
   call ncapt(ncid, tid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, tid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, tid, 'description', ncchar, 45, &
               'perturbation potential temperature (theta-t0)', ierr)
   call ncaptc(ncid, tid, 'units', ncchar, 1, 'k', ierr)
   call ncaptc(ncid, tid, 'stagger', ncchar, 0, empty, ierr)
   ! mu:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   muid = ncvdef(ncid, 'mu', ncfloat, 3, nd, ierr)
   call ncapt(ncid, muid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, muid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, muid, 'description', ncchar, 35, &
               'perturbation dry air mass in column', ierr)
   call ncaptc(ncid, muid, 'units', ncchar, 7, 'pascals', ierr)
   call ncaptc(ncid, muid, 'stagger', ncchar, 0, empty, ierr)
   ! mub:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   mubid = ncvdef(ncid, 'mub', ncfloat, 3, nd, ierr)
   call ncapt(ncid, mubid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, mubid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, mubid, 'description', ncchar, 32, &
               'base state dry air mass in column', ierr)
   call ncaptc(ncid, mubid, 'units', ncchar, 7, 'pascals', ierr)
   call ncaptc(ncid, mubid, 'stagger', ncchar, 0, empty, ierr)
   ! qvapor:
   nd(1:3) = ndm1(1:3)
   nd(4) = time
   qid = ncvdef(ncid, 'qvapor', ncfloat, 4, nd, ierr)
   call ncapt(ncid, qid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, qid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, qid, 'description', ncchar, 1, '-', ierr)
   call ncaptc(ncid, qid, 'units', ncchar, 1, '-', ierr)
   call ncaptc(ncid, qid, 'stagger', ncchar, 0, empty, ierr)
   ! mapfac_m:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   mapid = ncvdef(ncid, 'mapfac_m', ncfloat, 3, nd, ierr)
   call ncapt(ncid, mapid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, mapid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, mapid, 'description', ncchar, 29, &
               'map scale factor on mass grid', ierr)
   call ncaptc(ncid, mapid, 'units', ncchar, 13, 'dimensionless', ierr)
   call ncaptc(ncid, mapid, 'stagger', ncchar, 0, empty, ierr)
   ! p_top:
   nd(1) = scal
   nd(2) = time
   ptid = ncvdef(ncid, 'p_top', ncfloat, 2, nd, ierr)
   call ncapt(ncid, ptid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, ptid, 'memoryorder', ncchar, 3, '0  ', ierr)
   call ncaptc(ncid, ptid, 'description', ncchar, 0, empty, ierr)
   call ncaptc(ncid, ptid, 'units', ncchar, 1, '-', ierr)
   call ncaptc(ncid, ptid, 'stagger', ncchar, 0, empty, ierr)
   ! znu:
   nd(1) = ndm1(3)
   nd(2) = time
   znuid = ncvdef(ncid, 'znu', ncfloat, 2, nd, ierr)
   call ncapt(ncid, znuid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, znuid, 'memoryorder', ncchar, 3, 'z  ', ierr)
   call ncaptc(ncid, znuid, 'description', ncchar, 32, &
               'eta values on half (mass) levels', ierr)
   call ncaptc(ncid, znuid, 'units', ncchar, 13, 'dimensionless', ierr)
   call ncaptc(ncid, znuid, 'stagger', ncchar, 0, empty, ierr)
   ! znw:
   nd(1) = ndim(3)
   nd(2) = time
   znwid = ncvdef(ncid, 'znw', ncfloat, 2, nd, ierr)
   call ncapt(ncid, znwid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, znwid, 'memoryorder', ncchar, 3, 'z  ', ierr)
   call ncaptc(ncid, znwid, 'description', ncchar, 29, &
               'eta values on full (w) levels', ierr)
   call ncaptc(ncid, znwid, 'units', ncchar, 13, 'dimensionless', ierr)
   call ncaptc(ncid, znwid, 'stagger', ncchar, 1, 'z', ierr)
   ! xlat:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   latid = ncvdef(ncid, 'xlat', ncfloat, 3, nd, ierr)
   call ncapt(ncid, latid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, latid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, latid, 'description', ncchar, 27, &
               'latitude, south is negative', ierr)
   call ncaptc(ncid, latid, 'units', ncchar, 6, 'degree', ierr)
   call ncaptc(ncid, latid, 'stagger', ncchar, 0, empty, ierr)
   ! xlong:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   lonid = ncvdef(ncid, 'xlong', ncfloat, 3, nd, ierr)
   call ncapt(ncid, lonid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, lonid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, lonid, 'description', ncchar, 28, &
               'longitude, west is negative', ierr)
   call ncaptc(ncid, lonid, 'units', ncchar, 6, 'degree', ierr)
   call ncaptc(ncid, lonid, 'stagger', ncchar, 0, empty, ierr)
   ! rdx:
   nd(1) = scal
   nd(2) = time
   rxid = ncvdef(ncid, 'rdx', ncfloat, 2, nd, ierr)
   call ncapt(ncid, rxid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, rxid, 'memoryorder', ncchar, 3, '0  ', ierr)
   call ncaptc(ncid, rxid, 'description', ncchar, 21, &
               'inverse x grid length', ierr)
   call ncaptc(ncid, rxid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, rxid, 'stagger', ncchar, 0, empty, ierr)
   ! rdy:
   nd(1) = scal
   nd(2) = time
   ryid = ncvdef(ncid, 'rdy', ncfloat, 2, nd, ierr)
   call ncapt(ncid, ryid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, ryid, 'memoryorder', ncchar, 3, '0  ', ierr)
   call ncaptc(ncid, ryid, 'description', ncchar, 21, &
               'inverse y grid length', ierr)
   call ncaptc(ncid, ryid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, ryid, 'stagger', ncchar, 0, empty, ierr)

   ! extra variables requested by gsi:
   ! ph:
   nd(1:2) = ndm1(1:2)
   nd(3) = ndim(3)
   nd(4) = time
   phid = ncvdef(ncid, 'ph', ncfloat, 4, nd, ierr)
   call ncapt(ncid, phbid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, phbid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, phbid, 'description', ncchar, 23, &
               'perturbation geopotential', ierr)
   call ncaptc(ncid, phbid, 'units', ncchar, 10, 'm{2} s{-2}', ierr)
   call ncaptc(ncid, phbid, 'stagger', ncchar, 1, 'z', ierr)
   ! phb:
   nd(1:2) = ndm1(1:2)
   nd(3) = ndim(3)
   nd(4) = time
   phbid = ncvdef(ncid, 'phb', ncfloat, 4, nd, ierr)
   call ncapt(ncid, phbid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, phbid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, phbid, 'description', ncchar, 23, &
               'base-state geopotential', ierr)
   call ncaptc(ncid, phbid, 'units', ncchar, 10, 'm{2} s{-2}', ierr)
   call ncaptc(ncid, phbid, 'stagger', ncchar, 1, 'z', ierr)
   ! landmask:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   lmkid = ncvdef(ncid, 'landmask', ncfloat, 3, nd, ierr)
   call ncapt(ncid, lmkid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, lmkid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, lmkid, 'description', ncchar, 9, &
               'land mask', ierr)
   call ncaptc(ncid, lmkid, 'units', ncchar, 4, 'flag', ierr)
   call ncaptc(ncid, lmkid, 'stagger', ncchar, 0, empty, ierr)
   ! xice:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   iceid = ncvdef(ncid, 'xice', ncfloat, 3, nd, ierr)
   call ncapt(ncid, iceid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, iceid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, iceid, 'description', ncchar, 7, &
               'sea ice', ierr)
   call ncaptc(ncid, iceid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, iceid, 'stagger', ncchar, 0, empty, ierr)
   ! sst:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   sstid = ncvdef(ncid, 'sst', ncfloat, 3, nd, ierr)
   call ncapt(ncid, sstid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, sstid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, sstid, 'description', ncchar, 23, &
               'sea surface temperature', ierr)
   call ncaptc(ncid, sstid, 'units', ncchar, 1, 'k', ierr)
   call ncaptc(ncid, sstid, 'stagger', ncchar, 0, empty, ierr)
   ! ivgtyp:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   vgid = ncvdef(ncid, 'ivgtyp', nclong, 3, nd, ierr)
   call ncapt(ncid, vgid, 'fieldtype', nclong, 1, 106, ierr)
   call ncaptc(ncid, vgid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, vgid, 'description', ncchar, 15, &
               'vegetation type', ierr)
   call ncaptc(ncid, vgid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, vgid, 'stagger', ncchar, 0, empty, ierr)
   ! isltyp:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   slid = ncvdef(ncid, 'isltyp', nclong, 3, nd, ierr)
   call ncapt(ncid, slid, 'fieldtype', nclong, 1, 106, ierr)
   call ncaptc(ncid, slid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, slid, 'description', ncchar, 9, &
               'soil type', ierr)
   call ncaptc(ncid, slid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, slid, 'stagger', ncchar, 0, empty, ierr)
   ! vegfra:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   vfid = ncvdef(ncid, 'vegfra', ncfloat, 3, nd, ierr)
   call ncapt(ncid, vfid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, vfid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, vfid, 'description', ncchar, 19, &
               'vegetation fraction', ierr)
   call ncaptc(ncid, vfid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, vfid, 'stagger', ncchar, 0, empty, ierr)
   ! snow:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   snwid = ncvdef(ncid, 'snow', ncfloat, 3, nd, ierr)
   call ncapt(ncid, snwid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, snwid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, snwid, 'description', ncchar, 21, &
               'snow water equivalent', ierr)
   call ncaptc(ncid, snwid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, snwid, 'stagger', ncchar, 0, empty, ierr)
   ! u10:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   u10id = ncvdef(ncid, 'u10', ncfloat, 3, nd, ierr)
   call ncapt(ncid, u10id, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, u10id, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, u10id, 'description', ncchar, 9, &
               'u at 10 m', ierr)
   call ncaptc(ncid, u10id, 'units', ncchar, 3, 'm/s', ierr)
   call ncaptc(ncid, u10id, 'stagger', ncchar, 0, empty, ierr)
   ! v10:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   v10id = ncvdef(ncid, 'v10', ncfloat, 3, nd, ierr)
   call ncapt(ncid, v10id, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, v10id, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, v10id, 'description', ncchar, 9, &
               'v at 10 m', ierr)
   call ncaptc(ncid, v10id, 'units', ncchar, 3, 'm/s', ierr)
   call ncaptc(ncid, v10id, 'stagger', ncchar, 0, empty, ierr)
   ! smois:
   nd(1:2) = ndm1(1:2)
   nd(3) = nsol
   nd(4) = time
   smsid = ncvdef(ncid, 'smois', ncfloat, 4, nd, ierr)
   call ncapt(ncid, smsid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, smsid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, smsid, 'description', ncchar, 13, &
               'soil moisture', ierr)
   call ncaptc(ncid, smsid, 'units', ncchar, 0, empty, ierr)
   call ncaptc(ncid, smsid, 'stagger', ncchar, 1, 'z', ierr)
   ! tslb:
   nd(1:2) = ndm1(1:2)
   nd(3) = nsol
   nd(4) = time
   tslbid = ncvdef(ncid, 'tslb', ncfloat, 4, nd, ierr)
   call ncapt(ncid, tslbid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, tslbid, 'memoryorder', ncchar, 3, 'xyz', ierr)
   call ncaptc(ncid, tslbid, 'description', ncchar, 15, &
               'soil tempeature', ierr)
   call ncaptc(ncid, tslbid, 'units', ncchar, 1, 'k', ierr)
   call ncaptc(ncid, tslbid, 'stagger', ncchar, 1, 'z', ierr)
   ! tsk:
   nd(1:2) = ndm1(1:2)
   nd(3) = time
   tskid = ncvdef(ncid, 'tsk', ncfloat, 3, nd, ierr)
   call ncapt(ncid, tskid, 'fieldtype', nclong, 1, 104, ierr)
   call ncaptc(ncid, tskid, 'memoryorder', ncchar, 3, 'xy ', ierr)
   call ncaptc(ncid, tskid, 'description', ncchar, 24, &
               'surface skin temperature', ierr)
   call ncaptc(ncid, tskid, 'units', ncchar, 1, 'k', ierr)
   call ncaptc(ncid, tskid, 'stagger', ncchar, 0, empty, ierr)

   ! end defining mode:
   call ncendf(ncid, ierr)

   !++++++++++++++++++++++++++++
   ! assign values to variables:
   !++++++++++++++++++++++++++++

   ! 1. times:
   start = (/1, 1, 1, 1/)
   count(1) = 19
   count(2:4) = 1
   call ncvptc(ncid, tmid, start, count, times, 19, ierr)

   ! 2. u:
   ! stagger: y and z:
   call staggerxy_3d(u, imax, jmax, kmax, imax, jmax - 1, kmax, u_out)
   call staggerz(u_out, imax, jmax - 1, kmax, ptop, dam, znw, 4, 0, unc1)
   count(1) = imax
   count(2) = jmax - 1
   count(3) = kmax - 1
   call ncvpt(ncid, uid, start, count, unc1, ierr)

   ! 3. v:
   ! stagger: x and z:
   call staggerxy_3d(v, imax, jmax, kmax, imax - 1, jmax, kmax, v_out)
   call staggerz(v_out, imax - 1, jmax, kmax, ptop, dam, znw, 4, 0, vnc1)
   count(1) = imax - 1
   count(2) = jmax
   count(3) = kmax - 1
   call ncvpt(ncid, vid, start, count, vnc1, ierr)

   ! 4. t:
   ! stagger: x, y and z:
   call staggerxy_3d(t, imax, jmax, kmax, imax - 1, jmax - 1, kmax, t_out)
   call staggerz(t_out, imax - 1, jmax - 1, kmax, ptop, dam, znw, 4, 0, tnc1)
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = kmax - 1
   call ncvpt(ncid, tid, start, count, tnc1, ierr)

   ! 16. ph:
   ! stagger: z:
   t_out = 0.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = kmax
   call ncvpt(ncid, phid, start, count, t_out, ierr)
   ! 17. phb:
   ! stagger: z:
   call staggerxy_3d(geo, imax, jmax, kmax, imax - 1, jmax - 1, kmax, t_out)
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = kmax
   call ncvpt(ncid, phbid, start, count, t_out, ierr)

   ! 5. mub:
   ! stagger: x, and y:
   ! see arw description page 36 for reference pressure:
   p0 = 1.0e5
   t0 = 300.0 !273.15 testing and will be added to laps_parameter
   a0 = 50.0

   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1

   do j = 1, count(2)
      do i = 1, count(1)
         terrain = 0.25*(topo(i, j) + topo(i + 1, j) &
                         + topo(i, j + 1) + topo(i + 1, j + 1))
         mub_out(i, j, 1) = p0*exp(-t0/a0 + sqrt((t0/a0)**2 - &
                                                 2.0*terrain*9.806/a0/287.0))
      end do
   end do

   call ncvpt(ncid, mubid, start, count, mub_out, ierr)

   ! 6. mu:
   ! stagger: x, and y:
   ! according to arw description, this should be dry pressure-reference pressure:
   call staggerxy_3d(dam, imax, jmax, kmax, imax - 1, jmax - 1, 1, mu_out)
   mu_out = mu_out - mub_out

   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, muid, start, count, mu_out, ierr)

   ! 7. qvapor:
   ! stagger: x, y and z:
   call staggerxy_3d(sh, imax, jmax, kmax, imax - 1, jmax - 1, kmax, qvapor_out)
   call staggerz(qvapor_out, imax - 1, jmax - 1, kmax, ptop, dam, znw, 4, 0, qnc1)
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = kmax - 1
   call ncvpt(ncid, qid, start, count, qnc1, ierr)

   ! 8. mapfac_m:
   ! stagger: x, and y:
   call staggerxy_3d(mapfac, imax, jmax, kmax, imax - 1, jmax - 1, 1, mapfac_m_out)
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, mapid, start, count, mapfac_m_out, ierr)

   ! 9. p_top:
   count(1) = 1
   count(2) = 1
   call ncvpt(ncid, ptid, start, count, ptop, ierr)

   ! 10. znu:
   ! stagger: z:
   count(1) = kmax - 1
   count(2) = 1
   call ncvpt(ncid, znuid, start, count, znu, ierr)

   ! 11. znw:
   ! unstagger: z:
   count(1) = kmax
   count(2) = 1
   call ncvpt(ncid, znwid, start, count, znw, ierr)

   ! 12. xlat:
   ! stagger: x, and y:
   call staggerxy_3d(lat, imax, jmax, kmax, imax - 1, jmax - 1, 1, xlat_out)
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, latid, start, count, xlat_out, ierr)

   ! 13. xlong:
   ! stagger: x, and y:
   call staggerxy_3d(lon, imax, jmax, kmax, imax - 1, jmax - 1, 1, xlong_out)
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, lonid, start, count, xlong_out, ierr)

   ! 14. rdx:
   count(1:2) = 1
   call ncvpt(ncid, rxid, start, count, 1.0/dxy, ierr)

   ! 15. rdy:
   count(1:2) = 1
   call ncvpt(ncid, ryid, start, count, 1.0/dxy, ierr)

   !+++++++++++++++++++++++++++++++++++++++++++
   ! assign fake values to the extra variables:
   !+++++++++++++++++++++++++++++++++++++++++++
   ! 18. landmask:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1) = 1.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, lmkid, start, count, tnc, ierr)

   ! 19. xice:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1) = 0.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, iceid, start, count, tnc, ierr)

   ! 20. sst:
   ! stagger: x, and y:
   ! arw description states: sea level temperature for mub
   ! not sure if gsi uses sst for that.
   tnc(1:imax - 1, 1:jmax - 1, 1) = t0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, sstid, start, count, tnc, ierr)

   ! 21. ivgtyp:
   ! stagger: x, and y:
   itnc(1:imax - 1, 1:jmax - 1, 1) = 2
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, vgid, start, count, itnc, ierr)

   ! 22. isltyp:
   ! stagger: x, and y:
   itnc(1:imax - 1, 1:jmax - 1, 1) = 1
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, slid, start, count, itnc, ierr)

   ! 23. vegfra:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1) = 0.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, vfid, start, count, tnc, ierr)

   ! 24. snow:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1) = 0.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, snwid, start, count, tnc, ierr)

   ! 25. u10:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1) = 0.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, u10id, start, count, tnc, ierr)

   ! 26. v10:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1) = 0.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, v10id, start, count, tnc, ierr)

   ! 27. smois:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1:nsol) = 0.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 4
   count(4) = 1
   call ncvpt(ncid, smsid, start, count, tnc, ierr)

   ! 28. tslb:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1:nsol) = 280.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 4
   count(4) = 1
   call ncvpt(ncid, tslbid, start, count, tnc, ierr)
   ! 29. tsk:
   ! stagger: x, and y:
   tnc(1:imax - 1, 1:jmax - 1, 1) = 280.0
   count(1) = imax - 1
   count(2) = jmax - 1
   count(3) = 1
   call ncvpt(ncid, tskid, start, count, tnc, ierr)

   ! close the netcdf file:
   call ncclos(ncid, ierr)

   return
end
