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

module laps_static

   implicit none
   integer                     :: x              ! x-dimension of grid
   integer                     :: y              ! y-dimension of grid
   integer                     :: z2             ! z-dimension of 2-d grids
   integer                     :: z3             ! z-dimension of 3-d press
   real                        :: la1            ! sw corner lat
   real                        :: lo1            ! se corner lon
   real                        :: la2            ! ne corner lat
   real                        :: lo2            ! ne corner lon
   real                        :: dx             ! x-direction grid spacing
   real                        :: dy             ! y-direction grid spacing
   real                        :: lov            ! orientation longitude
   real                        :: latin1         ! standard lat 1
   real                        :: latin2         ! standard lat 2
   real, allocatable           :: topo(:, :)      ! laps topographic height
   real, allocatable           :: lats(:, :)      ! laps latitude array
   real, allocatable           :: lons(:, :)      ! laps longitude array
   character(len=132)         :: grid_type      ! map projection type

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_horiz_grid_spec(laps_data_root)

      implicit none
      character(len=*), intent(in) :: laps_data_root
      character(len=256)           :: static_file
      character(len=132)           :: dum
      integer :: cdfid, rcode, xid, yid, vid
      integer, dimension(2) :: startc, countc
      integer, dimension(4) :: start, count
      include "netcdf.inc"

      static_file = trim(laps_data_root)//'/static/static.nest7grid'

      ! open the static file which is a netcdf file

      cdfid = ncopn(trim(static_file), ncnowrit, rcode)

      ! get x/y dimensions.

      xid = ncdid(cdfid, 'x', rcode)
      yid = ncdid(cdfid, 'y', rcode)
      call ncdinq(cdfid, xid, dum, x, rcode)
      call ncdinq(cdfid, yid, dum, y, rcode)

      ! get the grid projection type

      startc = (/1, 1/)
      countc = (/132, 1/)
      vid = ncvid(cdfid, 'grid_type', rcode)
      call ncvgtc(cdfid, vid, startc, countc, grid_type, 132, rcode)

      ! get the corner points by reading lat/lon array
      allocate (lats(x, y))
      allocate (lons(x, y))
      vid = ncvid(cdfid, 'lat', rcode)
      start = (/1, 1, 1, 1/)
      count = (/x, y, 1, 1/)
      call ncvgt(cdfid, vid, start, count, lats, rcode)
      vid = ncvid(cdfid, 'lon', rcode)
      start = (/1, 1, 1, 1/)
      count = (/x, y, 1, 1/)
      call ncvgt(cdfid, vid, start, count, lons, rcode)
      la1 = lats(1, 1)
      lo1 = lons(1, 1)
      la2 = lats(x, y)
      lo2 = lons(x, y)
      if (lo1 .gt. 180) lo1 = lo1 - 360.
      if (lo2 .gt. 180) lo2 = lo2 - 360.
      if (la1 .gt. 270) la1 = la1 - 360.
      if (la2 .gt. 270) la2 = la2 - 360.
      print '(a,2f10.2)', 'sw corner lat/lon = ', la1, lo1
      print '(a,2f10.2)', 'ne corner lat/lon = ', la2, lo2

      ! get other projection parameters

      vid = ncvid(cdfid, 'dx', rcode)
      call ncvgt1(cdfid, vid, 1, dx, rcode)
      dx = dx/1000

      vid = ncvid(cdfid, 'dy', rcode)
      call ncvgt1(cdfid, vid, 1, dy, rcode)
      dy = dy/1000

      vid = ncvid(cdfid, 'lov', rcode)
      call ncvgt1(cdfid, vid, 1, lov, rcode)
      if (lov .gt. 180) then
         lov = lov - 360
      end if

      vid = ncvid(cdfid, 'latin1', rcode)
      call ncvgt1(cdfid, vid, 1, latin1, rcode)

      vid = ncvid(cdfid, 'latin2', rcode)
      call ncvgt1(cdfid, vid, 1, latin2, rcode)

      ! get the topography

      allocate (topo(x, y))
      vid = ncvid(cdfid, 'avg', rcode)
      start = (/1, 1, 1, 1/)
      count = (/x, y, 1, 1/)
      call ncvgt(cdfid, vid, start, count, topo, rcode)

   end subroutine get_horiz_grid_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module laps_static
