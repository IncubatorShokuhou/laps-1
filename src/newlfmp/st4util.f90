module st4util

   implicit none

   save

   integer :: ncid

end module

!===============================================================================

subroutine get_st4_dims(fname, nx, ny, nz)

   use st4util

   implicit none

   include 'netcdf.inc'

   integer :: nx, ny, nz, nid, icode
   character(len=*) :: fname
   logical :: there

! open stage iv file, and leave open for future use.

   inquire (file=trim(fname), exist=there)
   if (there) then
      icode = nf_open(trim(fname), nf_nowrite, ncid)
      if (ncid <= 0) then
         print *, 'could not open st4 file: ', trim(fname)
         stop
      end if
   else
      print *, 'could not find nmm file: ', trim(fname)
      stop
   end if

! read stage iv grid dimensions.

   icode = nf_inq_dimid(ncid, 'x', nid)
   icode = nf_inq_dimlen(ncid, nid, nx)
   icode = nf_inq_dimid(ncid, 'y', nid)
   icode = nf_inq_dimlen(ncid, nid, ny)
   icode = nf_inq_dimid(ncid, 'record', nid)
   icode = nf_inq_dimlen(ncid, nid, nz)

   return
end

!===============================================================================

subroutine fill_st4_grid

   use lfmgrid
   use st4util
   use constants

   implicit none

   include 'netcdf.inc'

   integer :: nid, icode, i, j
   real, parameter :: lat1 = 23.117, lon1 = 240.977

! stage iv data uses the ncep 255 grid, which is a user-defined grid.
! grid parameter settings as of 11 may 2004 are given at the following urls.
! http://www.emc.ncep.noaa.gov/mmb/ylin/pcpanl/qanda/gridchange/gridchange.htm
! http://wwwt.emc.ncep.noaa.gov/mmb/ylin/pcpanl/qanda/#gridinfo
! polar stereographic true at 60n
! y-axis is parallel to 105w.
! grid point spacing is 4.7625 at 60n
! lat1=23.117, lon1=240.977
! lat1, lon1 are the latitude and longitude (degrees east) of gridpoint (1,1)

   nprojection = 'polar stereographic'
   ntruelat1 = 60.0
   ntruelat2 = ntruelat1
   nstdlon = -105.0
   ngrid_spacingx = 4.7625
   ngrid_spacingy = ngrid_spacingx
!icode=nf_get_att_real(ncid,nf_global,'dx',ngrid_spacingx)
!icode=nf_get_att_real(ncid,nf_global,'dy',ngrid_spacingy)

! caculate latitude/longitude locations of stage iv grid points
   do j = 1, ny
   do i = 1, nx
      call w3fb07(float(i), float(j), lat1, lon1, ngrid_spacingx*1000., &
                  nstdlon + 360., nlat(i, j), nlon(i, j))
      nlon(i, j) = nlon(i, j) - 360.
   end do
   end do

! read accumulated precip.
! and convert from mm to m.

   icode = nf_inq_varid(ncid, 'apcp', nid)
   if (icode .eq. 0) then
      icode = nf_get_var_real(ncid, nid, npcp_tot)
   end if
   npcp_tot = npcp_tot/1000.

   return
end

!===============================================================================
subroutine w3fb07(xi, xj, alat1, alon1, dx, alonv, alat, alon)
!$$$   subprogram  documentation  block
!
! subprogram:  w3fb07        grid coords to lat/lon for grib
!   prgmmr: stackpole        org: nmc42       date:88-04-05
!
! abstract: converts the coordinates of a location on earth given in a
!   grid coordinate system overlaid on a polar stereographic map pro-
!   jection true at 60 degrees n or s latitude to the
!   natural coordinate system of latitude/longitude
!   w3fb07 is the reverse of w3fb06.
!   uses grib specification of the location of the grid
!
! program history log:
!   88-01-01  original author:  stackpole, w/nmc42
!   90-04-12  r.e.jones   convert to cray cft77 fortran
!
! usage:  call w3fb07(xi,xj,alat1,alon1,dx,alonv,alat,alon)
!   input argument list:
!     xi       - i coordinate of the point  real*4
!     xj       - j coordinate of the point  real*4
!     alat1    - latitude  of lower left point of grid (point 1,1)
!                latitude <0 for southern hemisphere; real*4
!     alon1    - longitude of lower left point of grid (point 1,1)
!                  east longitude used throughout; real*4
!     dx       - mesh length of grid in meters at 60 deg lat
!                 must be set negative if using
!                 southern hemisphere projection; real*4
!                   190500.0 lfm grid,
!                   381000.0 nh pe grid, -381000.0 sh pe grid, etc.
!     alonv    - the orientation of the grid.  i.e.,
!                the east longitude value of the vertical meridian
!                which is parallel to the y-axis (or columns of
!                the grid) along which latitude increases as
!                the y-coordinate increases.  real*4
!                   for example:
!                   255.0 for lfm grid,
!                   280.0 nh pe grid, 100.0 sh pe grid, etc.
!
!   output argument list:
!     alat     - latitude in degrees (negative in southern hemi.)
!     alon     - east longitude in degrees, real*4
!
!   remarks: formulae and notation loosely based on hoke, hayes,
!     and renninger's "map projections and grid systems...", march 1981
!     afgwc/tn-79/003
!
! attributes:
!   language: cray cft77 fortran
!   machine:  cray y-mp8/832
!
!$$$
!
   data rerth/6.3712e+6/, pi/3.1416/
   data ss60/1.86603/
!
!        preliminary variables and redifinitions
!
!        h = 1 for northern hemisphere; = -1 for southern
!
!        reflon is longitude upon which the positive x-coordinate
!        drawn through the pole and to the right lies
!        rotated around from orientation (y-coordinate) longitude
!        differently in each hemisphere
!
   if (dx .lt. 0) then
      h = -1.0
      dxl = -dx
      reflon = alonv - 90.0
   else
      h = 1.0
      dxl = dx
      reflon = alonv - 270.0
   end if
!
   radpd = pi/180.0
   degprd = 180.0/pi
   rebydx = rerth/dxl
!
!        radius to lower left hand (ll) corner
!
   ala1 = alat1*radpd
   rmll = rebydx*cos(ala1)*ss60/(1.+h*sin(ala1))
!
!        use ll point info to locate pole point
!
   alo1 = (alon1 - reflon)*radpd
   polei = 1.-rmll*cos(alo1)
   polej = 1.-h*rmll*sin(alo1)
!
!        radius to the i,j point (in grid units)
!
   xx = xi - polei
   yy = (xj - polej)*h
   r2 = xx**2 + yy**2
!
!        now the magic formulae
!
   if (r2 .eq. 0) then
      alat = h*90.
      alon = reflon
   else
      gi2 = (rebydx*ss60)**2
      alat = degprd*h*asin((gi2 - r2)/(gi2 + r2))
      arccos = acos(xx/sqrt(r2))
      if (yy .gt. 0) then
         alon = reflon + degprd*arccos
      else
         alon = reflon - degprd*arccos
      end if
   end if
   if (alon .lt. 0) alon = alon + 360.
!
   return
end
