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

module map_utils

! module that defines constants, data structures, and
! subroutines used to convert grid indices to lat/lon
! and vice versa.
!
! supported projections
! ---------------------
! cylindrical lat/lon (code = proj_latlon)
! mercator (code = proj_merc)
! lambert conformal (code = proj_lc)
! polar stereographic (code = proj_ps)
!
! remarks
! -------
! the routines contained within were adapted from routines
! obtained from the ncep w3 library.  the original ncep routines were less
! flexible (e.g., polar-stereo routines only supported truelat of 60n/60s)
! than what we needed, so modifications based on equations in hoke, hayes, and
! renninger (afgwc/tn/79-003) were added to improve the flexibility.
! additionally, coding was improved to f90 standards and the routines were
! combined into this module.
!
! assumptions
! -----------
!  grid definition:
!    for mercator, lambert conformal, and polar-stereographic projections,
!    the routines within assume the following:
!
!       1.  grid is dimensioned (i,j) where i is the east-west direction,
!           positive toward the east, and j is the north-south direction,
!           positive toward the north.
!       2.  origin is at (1,1) and is located at the southwest corner,
!           regardless of hemispere.
!       3.  grid spacing (dx) is always positive.
!       4.  values of true latitudes must be positive for nh domains
!           and negative for sh domains.
!
!     for the latlon projection, the grid origin may be at any of the
!     corners, and the deltalat and deltalon values can be signed to
!     account for this using the following convention:
!       origin location        deltalat sign      deltalon sign
!       ---------------        -------------      -------------
!        sw corner                  +                   +
!        ne corner                  -                   -
!        nw corner                  -                   +
!        se corner                  +                   -
!
!  data definitions:
!       1. any arguments that are a latitude value are expressed in
!          degrees north with a valid range of -90 -> 90
!       2. any arguments that are a longitude value are expressed in
!          degrees east with a valid range of -180 -> 180.
!       3. distances are in meters and are always positive.
!       4. the standard longitude (stdlon) is defined as the longitude
!          line which is parallel to the y-axis (j-direction), along
!          which latitude increases (not the absolute value of latitude, but
!          the actual latitude, such that latitude increases continuously
!          from the south pole to the north pole) as j increases.
!       5. one true latitude value is required for polar-stereographic and
!          mercator projections, and defines at which latitude the
!          grid spacing is true.  for lambert conformal, two true latitude
!          values must be specified, but may be set equal to each other to
!          specify a tangent projection instead of a secant projection.
!
! usage
! -----
! to use the routines in this module, the calling routines must have the
! following statement at the beginning of its declaration block:
!   use map_utils
!
! the use of the module not only provides access to the necessary routines,
! but also defines a structure of type (proj_info) that can be used
! to declare a variable of the same type to hold your map projection
! information.  it also defines some integer parameters that contain
! the projection codes so one only has to use those variable names rather
! than remembering the acutal code when using them.  the basic steps are
! as follows:
!
!   1.  ensure the "use map_utils" is in your declarations.
!   2.  declare the projection information structure as type(proj_info):
!         type(proj_info) :: proj
!   3.  populate your structure by calling the map_set routine:
!         call map_set(code,lat1,lon1,dx,stdlon,truelat1,truelat2,nx,ny,proj)
!       where:
!         code (input) = one of proj_latlon, proj_merc, proj_lc, or proj_ps
!         lat1 (input) = latitude of grid origin point (i,j)=(1,1)
!                         (see assumptions!)
!         lon1 (input) = longitude of grid origin
!         dx (input) = grid spacing in meters (ignored for latlon projections)
!         stdlon (input) = standard longitude for proj_ps and proj_lc,
!               deltalon (see assumptions) for proj_latlon,
!               ignored for proj_merc
!         truelat1 (input) = 1st true latitude for proj_ps, proj_lc, and
!                proj_merc, deltalat (see assumptions) for proj_latlon
!         truelat2 (input) = 2nd true latitude for proj_lc,
!                ignored for all others.
!         nx = number of points in east-west direction
!         ny = number of points in north-south direction
!         proj (output) = the structure of type (proj_info) that will be fully
!                populated after this call
!
!   4.  now that the proj structure is populated, you may call any
!       of the following routines:
!
!       latlon_to_ij(proj, lat, lon, i, j)
!       ij_to_latlon(proj, i, j, lat, lon)
!       truewind_to_gridwind(lon, proj, ugrid, vgrid, utrue, vtrue)
!       gridwind_to_truewind(lon, proj, utrue, vtrue, ugrid, vgrid)
!       compare_projections(proj1, proj2, same_proj)
!
!       it is incumbent upon the calling routine to determine whether or
!       not the values returned are within your domain bounds.  all values
!       of i, j, lat, and lon are real values.
!
!
! references
! ----------
!  hoke, hayes, and renninger, "map preojections and grid systems for
!       meteorological applications." afgwc/tn-79/003(rev), air weather
!       service, 1985.
!
!  ncar mm5v3 modeling system, regridder program, module_first_guess_map.f
!  ncep routines w3fb06, w3fb07, w3fb08, w3fb09, w3fb11, w3fb12
!
! history
! -------
! 27 mar 2001 - original version
!               brent l. shaw, noaa/fsl (csu/cira)
! 02 apr 2001 - added routines to rotate winds from true to grid
!               and vice versa.
!               brent l. shaw, noaa/fsl (csu/cira)
! 09 apr 2001 - added compare_projections routine to compare two
!               sets of projection parameters.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   ! define some private constants
   real, private, parameter     :: pi = 3.1415927
   real, private, parameter    :: deg_per_rad = 180./pi
   real, private, parameter    :: rad_per_deg = pi/180.

   ! mean earth radius in m.  the value below is consistent
   ! with nceps routines and grids.
   real, public, parameter     :: earth_radius_m = 6371200.

   ! define public parameters

   ! projection codes for proj_info structure:
   integer, public, parameter  :: proj_latlon = 0
   integer, public, parameter  :: proj_merc = 1
   integer, public, parameter  :: proj_lc = 3
   integer, public, parameter  :: proj_ps = 5

   ! define data structures to define various projections

   type proj_info

      integer          :: code     ! integer code for projection type
      real             :: lat1    ! sw latitude (1,1) in degrees (-90->90n)
      real             :: lon1    ! sw longitude (1,1) in degrees (-180->180e)
      real             :: dx       ! grid spacing in meters at truelats, used
      ! only for ps, lc, and merc projections
      real             :: dlat     ! lat increment for lat/lon grids
      real             :: dlon     ! lon increment for lat/lon grids
      real             :: stdlon   ! longitude parallel to y-axis (-180->180e)
      real             :: truelat1 ! first true latitude (all projections)
      real             :: truelat2 ! second true lat (lc only)
      real             :: hemi     ! 1 for nh, -1 for sh
      real             :: cone     ! cone factor for lc projections
      real             :: polei    ! computed i-location of pole point
      real             :: polej    ! computed j-location of pole point
      real             :: rsw      ! computed radius to sw corner
      real             :: rebydx   ! earth radius divided by dx
      logical          :: init     ! flag to indicate if this struct is
      ! ready for use
      integer          :: nx
      integer          :: ny
   end type proj_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine map_init(proj)
      ! initializes the map projection structure to missing values

      implicit none
      type(proj_info), intent(inout)  :: proj

      proj%lat1 = -999.9
      proj%lon1 = -999.9
      proj%dx = -999.9
      proj%stdlon = -999.9
      proj%truelat1 = -999.9
      proj%truelat2 = -999.9
      proj%hemi = 0.0
      proj%cone = -999.9
      proj%polei = -999.9
      proj%polej = -999.9
      proj%rsw = -999.9
      proj%init = .false.
      proj%nx = -99
      proj%ny = -99
   end subroutine map_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine map_set(proj_code, lat1, lon1, dx, stdlon, truelat1, truelat2, &
                      idim, jdim, proj)
      ! given a partially filled proj_info structure, this routine computes
      ! polei, polej, rsw, and cone (if lc projection) to complete the
      ! structure.  this allows us to eliminate redundant calculations when
      ! calling the coordinate conversion routines multiple times for the
      ! same map.
      ! this will generally be the first routine called when a user wants
      ! to be able to use the coordinate conversion routines, and it
      ! will call the appropriate subroutines based on the
      ! proj%code which indicates which projection type  this is.

      implicit none

      ! declare arguments
      integer, intent(in)               :: proj_code
      real, intent(in)                  :: lat1
      real, intent(in)                  :: lon1
      real, intent(in)                  :: dx
      real, intent(in)                  :: stdlon
      real, intent(in)                  :: truelat1
      real, intent(in)                  :: truelat2
      integer, intent(in)               :: idim
      integer, intent(in)               :: jdim
      type(proj_info), intent(out)      :: proj

      ! local variables

      ! executable code

      ! first, check for validity of mandatory variables in proj
      if (abs(lat1) .gt. 90.001) then
         print '(a)', 'latitude of origin corner required as follows:'
         print '(a)', '    -90n <= lat1 < = 90.n'
         stop 'map_init'
      end if
      if (abs(lon1) .gt. 180.) then
         print '(a)', 'longitude of origin required as follows:'
         print '(a)', '   -180e <= lon1 <= 180w'
         stop 'map_init'
      end if
      if ((dx .le. 0.) .and. (proj_code .ne. proj_latlon)) then
         print '(a)', 'require grid spacing (dx) in meters be positive!'
         stop 'map_init'
      end if
      if ((abs(stdlon) .gt. 180.) .and. (proj_code .ne. proj_merc)) then
         print '(a)', 'need orientation longitude (stdlon) as: '
         print '(a)', '   -180e <= lon1 <= 180w'
         stop 'map_init'
      end if
      if (abs(truelat1) .gt. 90.) then
         print '(a)', 'set true latitude 1 for all projections!'
         stop 'map_init'
      end if

      call map_init(proj)
      proj%code = proj_code
      proj%lat1 = lat1
      proj%lon1 = lon1
      proj%dx = dx
      proj%stdlon = stdlon
      proj%truelat1 = truelat1
      proj%truelat2 = truelat2
      proj%nx = idim
      proj%ny = jdim
      if (proj%code .ne. proj_latlon) then
         proj%dx = dx
         if (truelat1 .lt. 0.) then
            proj%hemi = -1.0
         else
            proj%hemi = 1.0
         end if
         proj%rebydx = earth_radius_m/dx
      end if
      pick_proj:select case(proj%code)

      case (proj_ps)
      !print '(a)', 'setting up polar stereographic map...'
      call set_ps(proj)

      case (proj_lc)
      !print '(a)', 'setting up lambert conformal map...'
      if (abs(proj%truelat2) .gt. 90.) then
         print '(a)', 'second true latitude not set, assuming a tangent'
         print '(a,f10.3)', 'projection at truelat1: ', proj%truelat1
         proj%truelat2 = proj%truelat1
      else
         ! ensure truelat1 < truelat2
         proj%truelat1 = min(truelat1, truelat2)
         proj%truelat2 = max(truelat1, truelat2)
      end if
      call set_lc(proj)

      case (proj_merc)
      !print '(a)', 'setting up mercator map...'
      call set_merc(proj)

      case (proj_latlon)
      !print '(a)', 'setting up cylindrical equidistant latlon map...'
      ! convert lon1 to 0->360 notation
      if (proj%lon1 .lt. 0.) proj%lon1 = proj%lon1 + 360.
      proj%dlat = truelat1
      proj%dlon = stdlon

      case default
      print '(a,i2)', 'unknown projection code: ', proj%code
      stop 'map_init'

      end select pick_proj
      proj%init = .true.
      return
   end subroutine map_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine latlon_to_ij(proj, lat, lon, i, j)
      ! converts input lat/lon values to the cartesian (i,j) value
      ! for the given projection.

      implicit none
      type(proj_info), intent(in)          :: proj
      real, intent(in)                     :: lat
      real, intent(in)                     :: lon
      real, intent(out)                    :: i
      real, intent(out)                    :: j

      if (.not. proj%init) then
         print '(a)', 'you have not called map_set for this projection!'
         stop 'latlon_to_ij'
      end if

      select case (proj%code)

      case (proj_latlon)
         call llij_latlon(lat, lon, proj, i, j)

      case (proj_merc)
         call llij_merc(lat, lon, proj, i, j)

      case (proj_ps)
         call llij_ps(lat, lon, proj, i, j)

      case (proj_lc)
         call llij_lc(lat, lon, proj, i, j)

      case default
         print '(a,i2)', 'unrecognized map projection code: ', proj%code
         stop 'latlon_to_ij'

      end select
      return
   end subroutine latlon_to_ij
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ij_to_latlon(proj, i, j, lat, lon)
      ! computes geographical latitude and longitude for a given (i,j) point
      ! in a grid with a projection of proj

      implicit none
      type(proj_info), intent(in)          :: proj
      real, intent(in)                    :: i
      real, intent(in)                    :: j
      real, intent(out)                   :: lat
      real, intent(out)                   :: lon

      if (.not. proj%init) then
         print '(a)', 'you have not called map_set for this projection!'
         stop 'ij_to_latlon'
      end if
      select case (proj%code)

      case (proj_latlon)
         call ijll_latlon(i, j, proj, lat, lon)

      case (proj_merc)
         call ijll_merc(i, j, proj, lat, lon)

      case (proj_ps)
         call ijll_ps(i, j, proj, lat, lon)

      case (proj_lc)
         call ijll_lc(i, j, proj, lat, lon)

      case default
         print '(a,i2)', 'unrecognized map projection code: ', proj%code
         stop 'ij_to_latlon'

      end select
      return
   end subroutine ij_to_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_ps(proj)
      ! initializes a polar-stereographic map projection from the partially
      ! filled proj structure. this routine computes the radius to the
      ! southwest corner and computes the i/j location of the pole for use
      ! in llij_ps and ijll_ps.
      implicit none

      ! declare args
      type(proj_info), intent(inout)    :: proj

      ! local vars
      real                              :: ala1
      real                              :: alo1
      real                              :: reflon
      real                              :: scale_top

      ! executable code
      reflon = proj%stdlon + 90.

      ! cone factor
      proj%cone = 1.0

      ! compute numerator term of map scale factor
      scale_top = 1.+proj%hemi*sin(proj%truelat1*rad_per_deg)

      ! compute radius to lower-left (sw) corner
      ala1 = proj%lat1*rad_per_deg
      proj%rsw = proj%rebydx*cos(ala1)*scale_top/(1.+proj%hemi*sin(ala1))

      ! find the pole point
      alo1 = (proj%lon1 - reflon)*rad_per_deg
      proj%polei = 1.-proj%rsw*cos(alo1)
      proj%polej = 1.-proj%hemi*proj%rsw*sin(alo1)
      print '(a,2f10.1)', 'computed (i,j) of pole point: ', proj%polei, proj%polej
      return
   end subroutine set_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine llij_ps(lat, lon, proj, i, j)
      ! given latitude (-90 to 90), longitude (-180 to 180), and the
      ! standard polar-stereographic projection information via the
      ! public proj structure, this routine returns the i/j indices which
      ! if within the domain range from 1->nx and 1->ny, respectively.

      implicit none

      ! delcare input arguments
      real, intent(in)               :: lat
      real, intent(in)               :: lon
      type(proj_info), intent(in)     :: proj

      ! declare output arguments
      real, intent(out)              :: i !(x-index)
      real, intent(out)              :: j !(y-index)

      ! declare local variables

      real                           :: reflon
      real                           :: scale_top
      real                           :: ala
      real                           :: alo
      real                           :: rm

      ! begin code

      reflon = proj%stdlon + 90.

      ! compute numerator term of map scale factor

      scale_top = 1.+proj%hemi*sin(proj%truelat1*rad_per_deg)

      ! find radius to desired point
      ala = lat*rad_per_deg
      rm = proj%rebydx*cos(ala)*scale_top/(1.+proj%hemi*sin(ala))
      alo = (lon - reflon)*rad_per_deg
      i = proj%polei + rm*cos(alo)
      j = proj%polej + proj%hemi*rm*sin(alo)

      return
   end subroutine llij_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ijll_ps(i, j, proj, lat, lon)

      ! this is the inverse subroutine of llij_ps.  it returns the
      ! latitude and longitude of an i/j point given the projection info
      ! structure.

      implicit none

      ! declare input arguments
      real, intent(in)                    :: i    ! column
      real, intent(in)                    :: j    ! row
      type(proj_info), intent(in)        :: proj

      ! declare output arguments
      real, intent(out)                   :: lat     ! -90 -> 90 north
      real, intent(out)                   :: lon     ! -180 -> 180 east

      ! local variables
      real                                :: reflon
      real                                :: scale_top
      real                                :: xx, yy
      real                                :: gi2, r2
      real                                :: arccos

      ! begin code

      ! compute the reference longitude by rotating 90 degrees to the east
      ! to find the longitude line parallel to the positive x-axis.
      reflon = proj%stdlon + 90.

      ! compute numerator term of map scale factor
      scale_top = 1.+proj%hemi*sin(proj%truelat1*rad_per_deg)

      ! compute radius to point of interest
      xx = i - proj%polei
      yy = (j - proj%polej)*proj%hemi
      r2 = xx**2 + yy**2

      ! now the magic code
      if (r2 .eq. 0.) then
         lat = proj%hemi*90.
         lon = reflon
      else
         gi2 = (proj%rebydx*scale_top)**2.
         lat = deg_per_rad*proj%hemi*asin((gi2 - r2)/(gi2 + r2))
         arccos = acos(xx/sqrt(r2))
         if (yy .gt. 0) then
            lon = reflon + deg_per_rad*arccos
         else
            lon = reflon - deg_per_rad*arccos
         end if
      end if

      ! convert to a -180 -> 180 east convention
      if (lon .gt. 180.) lon = lon - 360.
      if (lon .lt. -180.) lon = lon + 360.
      return

   end subroutine ijll_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_lc(proj)
      ! initialize the remaining items in the proj structure for a
      ! lambert conformal grid.

      implicit none

      type(proj_info), intent(inout)     :: proj

      real                               :: arg
      real                               :: deltalon1
      real                               :: tl1r
      real                               :: ctl1r

      ! compute cone factor
      call lc_cone(proj%truelat1, proj%truelat2, proj%cone)
      ! print '(a,f8.6)', 'computed cone factor: ', proj%cone
      ! compute longitude differences and ensure we stay out of the
      ! forbidden "cut zone"
      deltalon1 = proj%lon1 - proj%stdlon
      if (deltalon1 .gt. +180.) deltalon1 = deltalon1 - 360.
      if (deltalon1 .lt. -180.) deltalon1 = deltalon1 + 360.

      ! convert truelat1 to radian and compute cos for later use
      tl1r = proj%truelat1*rad_per_deg
      ctl1r = cos(tl1r)

      ! compute the radius to our known lower-left (sw) corner
      proj%rsw = proj%rebydx*ctl1r/proj%cone* &
                 (tan((90.*proj%hemi - proj%lat1)*rad_per_deg/2.)/ &
                  tan((90.*proj%hemi - proj%truelat1)*rad_per_deg/2.))**proj%cone

      ! find pole point
      arg = proj%cone*(deltalon1*rad_per_deg)
      proj%polei = 1.-proj%hemi*proj%rsw*sin(arg)
      proj%polej = 1.+proj%rsw*cos(arg)
      !print '(a,2f10.3)', 'computed pole i/j = ', proj%polei, proj%polej
      return
   end subroutine set_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine lc_cone(truelat1, truelat2, cone)

      ! subroutine to compute the cone factor of a lambert conformal projection

      implicit none

      ! input args
      real, intent(in)             :: truelat1  ! (-90 -> 90 degrees n)
      real, intent(in)             :: truelat2  !   "   "  "   "     "

      ! output args
      real, intent(out)            :: cone

      ! locals

      ! begin code

      ! first, see if this is a secant or tangent projection.  for tangent
      ! projections, truelat1 = truelat2 and the cone is tangent to the
      ! earth surface at this latitude.  for secant projections, the cone
      ! intersects the earth surface at each of the distinctly different
      ! latitudes
      if (abs(truelat1 - truelat2) .gt. 0.1) then

         ! compute cone factor following:
         cone = (alog(cos(truelat1*rad_per_deg)) - alog(cos(truelat2*rad_per_deg)))/ &
                (alog(tan((90.-abs(truelat1))*rad_per_deg*0.5)) - &
                 alog(tan((90.-abs(truelat2))*rad_per_deg*0.5)))
      else
         cone = sin(abs(truelat1)*rad_per_deg)
      end if
      return
   end subroutine lc_cone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ijll_lc(i, j, proj, lat, lon)

      ! subroutine to convert from the (i,j) cartesian coordinate to the
      ! geographical latitude and longitude for a lambert conformal projection.

      ! history:
      ! 25 jul 01: corrected by b. shaw, noaa/fsl
      !
      implicit none

      ! input args
      real, intent(in)              :: i        ! cartesian x coordinate
      real, intent(in)              :: j        ! cartesian y coordinate
      type(proj_info), intent(in)    :: proj     ! projection info structure

      ! output args
      real, intent(out)             :: lat      ! latitude (-90->90 deg n)
      real, intent(out)             :: lon      ! longitude (-180->180 e)

      ! locals
      real                          :: inew
      real                          :: jnew
      real                          :: r
      real                          :: chi, chi1, chi2
      real                          :: r2
      real                          :: xx
      real                          :: yy

      ! begin code

      chi1 = (90.-proj%hemi*proj%truelat1)*rad_per_deg
      chi2 = (90.-proj%hemi*proj%truelat2)*rad_per_deg

      ! see if we are in the southern hemispere and flip the indices
      ! if we are.
      if (proj%hemi .eq. -1.) then
         inew = -i + 2.
         jnew = -j + 2.
      else
         inew = i
         jnew = j
      end if

      ! compute radius**2 to i/j location
      xx = inew - proj%polei
      yy = proj%polej - jnew
      r2 = (xx*xx + yy*yy)
      r = sqrt(r2)/proj%rebydx

      ! convert to lat/lon
      if (r2 .eq. 0.) then
         lat = proj%hemi*90.
         lon = proj%stdlon
      else

         ! longitude
         lon = proj%stdlon + deg_per_rad*atan2(proj%hemi*xx, yy)/proj%cone
         lon = amod(lon + 360., 360.)

         ! latitude.  latitude determined by solving an equation adapted
         ! from:
         !  maling, d.h., 1973: coordinate systems and map projections
         ! equations #20 in appendix i.

         if (chi1 .eq. chi2) then
            chi = 2.0*atan((r/tan(chi1))**(1./proj%cone)*tan(chi1*0.5))
         else
            chi = 2.0*atan((r*proj%cone/sin(chi1))**(1./proj%cone)*tan(chi1*0.5))
         end if
         lat = (90.0 - chi*deg_per_rad)*proj%hemi

      end if

      if (lon .gt. +180.) lon = lon - 360.
      if (lon .lt. -180.) lon = lon + 360.
      return
   end subroutine ijll_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine llij_lc(lat, lon, proj, i, j)

      ! subroutine to compute the geographical latitude and longitude values
      ! to the cartesian x/y on a lambert conformal projection.

      implicit none

      ! input args
      real, intent(in)              :: lat      ! latitude (-90->90 deg n)
      real, intent(in)              :: lon      ! longitude (-180->180 e)
      type(proj_info), intent(in)      :: proj     ! projection info structure

      ! output args
      real, intent(out)             :: i        ! cartesian x coordinate
      real, intent(out)             :: j        ! cartesian y coordinate

      ! locals
      real                          :: arg
      real                          :: deltalon
      real                          :: tl1r
      real                          :: rm
      real                          :: ctl1r

      ! begin code

      ! compute deltalon between known longitude and standard lon and ensure
      ! it is not in the cut zone
      deltalon = lon - proj%stdlon
      if (deltalon .gt. +180.) deltalon = deltalon - 360.
      if (deltalon .lt. -180.) deltalon = deltalon + 360.

      ! convert truelat1 to radian and compute cos for later use
      tl1r = proj%truelat1*rad_per_deg
      ctl1r = cos(tl1r)

      ! radius to desired point
      rm = proj%rebydx*ctl1r/proj%cone* &
           (tan((90.*proj%hemi - lat)*rad_per_deg/2.)/ &
            tan((90.*proj%hemi - proj%truelat1)*rad_per_deg/2.))**proj%cone

      arg = proj%cone*(deltalon*rad_per_deg)
      i = proj%polei + proj%hemi*rm*sin(arg)
      j = proj%polej - rm*cos(arg)

      ! finally, if we are in the southern hemisphere, flip the i/j
      ! values to a coordinate system where (1,1) is the sw corner
      ! (what we assume) which is different than the original ncep
      ! algorithms which used the ne corner as the origin in the
      ! southern hemisphere (left-hand vs. right-hand coordinate?)
      if (proj%hemi .eq. -1.) then
         i = 2.-i
         j = 2.-j
      end if
      return
   end subroutine llij_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine set_merc(proj)

      ! sets up the remaining basic elements for the mercator projection

      implicit none
      type(proj_info), intent(inout)       :: proj
      real                                 :: clain

      !  preliminary variables

      clain = cos(rad_per_deg*proj%truelat1)
      proj%dlon = proj%dx/(earth_radius_m*clain)

      ! compute distance from equator to origin, and store in the
      ! proj%rsw tag.

      proj%rsw = 0.
      if (proj%lat1 .ne. 0.) then
         proj%rsw = (alog(tan(0.5*((proj%lat1 + 90.)*rad_per_deg))))/proj%dlon
      end if
      return
   end subroutine set_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine llij_merc(lat, lon, proj, i, j)

      ! compute i/j coordinate from lat lon for mercator projection

      implicit none
      real, intent(in)              :: lat
      real, intent(in)              :: lon
      type(proj_info), intent(in)    :: proj
      real, intent(out)              :: i
      real, intent(out)              :: j
      real                          :: deltalon

      deltalon = lon - proj%lon1
      if (deltalon .lt. -180.) deltalon = deltalon + 360.
      if (deltalon .gt. 180.) deltalon = deltalon - 360.
      i = 1.+(deltalon/(proj%dlon*deg_per_rad))
      j = 1.+(alog(tan(0.5*((lat + 90.)*rad_per_deg))))/ &
          proj%dlon - proj%rsw

      return
   end subroutine llij_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ijll_merc(i, j, proj, lat, lon)

      ! compute the lat/lon from i/j for mercator projection

      implicit none
      real, intent(in)               :: i
      real, intent(in)               :: j
      type(proj_info), intent(in)    :: proj
      real, intent(out)             :: lat
      real, intent(out)             :: lon

      lat = 2.0*atan(exp(proj%dlon*(proj%rsw + j - 1.)))*deg_per_rad - 90.
      lon = (i - 1.)*proj%dlon*deg_per_rad + proj%lon1
      if (lon .gt. 180.) lon = lon - 360.
      if (lon .lt. -180.) lon = lon + 360.
      return
   end subroutine ijll_merc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine llij_latlon(lat, lon, proj, i, j)

      ! compute the i/j location of a lat/lon on a latlon grid.
      implicit none
      real, intent(in)             :: lat
      real, intent(in)             :: lon
      type(proj_info), intent(in)  :: proj
      real, intent(out)            :: i
      real, intent(out)            :: j

      real                         :: deltalat
      real                         :: deltalon
      real                         :: lon360
      real                         :: latinc
      real                         :: loninc

      ! compute deltalat and deltalon as the difference between the input
      ! lat/lon and the origin lat/lon

      deltalat = lat - proj%lat1

      ! to account for issues around the dateline, convert the incoming
      ! longitudes to be 0->360.
      if (lon .lt. 0) then
         lon360 = lon + 360.
      else
         lon360 = lon
      end if
      deltalon = lon360 - proj%lon1
      if (deltalon .lt. 0) deltalon = deltalon + 360.

      ! compute i/j
      i = deltalon/proj%dlon + 1.
      j = deltalat/proj%dlat + 1.
      return
   end subroutine llij_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ijll_latlon(i, j, proj, lat, lon)

      ! compute the lat/lon location of an i/j on a latlon grid.
      implicit none
      real, intent(in)             :: i
      real, intent(in)             :: j
      type(proj_info), intent(in)  :: proj
      real, intent(out)            :: lat
      real, intent(out)            :: lon

      real                         :: deltalat
      real                         :: deltalon
      real                         :: lon360

      ! compute deltalat and deltalon

      deltalat = (j - 1.)*proj%dlat
      deltalon = (i - 1.)*proj%dlon
      lat = proj%lat1 + deltalat
      lon = proj%lon1 + deltalon

      if ((abs(lat) .gt. 90.) .or. (abs(deltalon) .gt. 360.)) then
         ! off the earth for this grid
         lat = -999.
         lon = -999.
      else
         lon = lon + 360.
         lon = amod(lon, 360.)
         if (lon .gt. 180.) lon = lon - 360.
      end if

      return
   end subroutine ijll_latlon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine gridwind_to_truewind(lon, proj, ugrid, vgrid, utrue, vtrue)

      ! subroutine to convert a wind from grid north to true north.

      implicit none

      ! arguments
      real, intent(in)          :: lon     ! longitude of point in degrees
      type(proj_info), intent(in):: proj    ! projection info structure
      real, intent(in)          :: ugrid   ! u-component, grid-relative
      real, intent(in)          :: vgrid   ! v-component, grid-relative
      real, intent(out)         :: utrue   ! u-component, earth-relative
      real, intent(out)         :: vtrue   ! v-component, earth-relative

      ! locals
      real                      :: alpha
      real                      :: diff

      if ((proj%code .eq. proj_ps) .or. (proj%code .eq. proj_lc)) then
         diff = lon - proj%stdlon
         if (diff .gt. 180.) diff = diff - 360.
         if (diff .lt. -180.) diff = diff + 360.
         alpha = diff*proj%cone*rad_per_deg*sign(1., proj%truelat1)
         utrue = vgrid*sin(alpha) + ugrid*cos(alpha)
         vtrue = vgrid*cos(alpha) - ugrid*sin(alpha)
      elseif ((proj%code .eq. proj_merc) .or. (proj%code .eq. proj_latlon)) then
         utrue = ugrid
         vtrue = vgrid
      else
         print '(a)', 'unrecognized projection.'
         stop 'gridwind_to_truewind'
      end if

      return
   end subroutine gridwind_to_truewind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine truewind_to_gridwind(lon, proj, utrue, vtrue, ugrid, vgrid)

      ! subroutine to compute grid-relative u/v wind components from the earth-
      ! relative values for a given projection.

      implicit none

      ! arguments
      real, intent(in)                 :: lon
      type(proj_info), intent(in)       :: proj
      real, intent(in)                 :: utrue
      real, intent(in)                 :: vtrue
      real, intent(out)                :: ugrid
      real, intent(out)                :: vgrid

      ! locals
      real                             :: alpha
      real                             :: diff

      if ((proj%code .eq. proj_ps) .or. (proj%code .eq. proj_lc)) then

         diff = proj%stdlon - lon
         if (diff .gt. 180.) diff = diff - 360.
         if (diff .lt. -180.) diff = diff + 360.
         alpha = diff*proj%cone*rad_per_deg*sign(1., proj%truelat1)
         ugrid = vtrue*sin(alpha) + utrue*cos(alpha)
         vgrid = vtrue*cos(alpha) - utrue*sin(alpha)

      elseif ((proj%code .eq. proj_merc) .or. (proj%code .eq. proj_latlon)) then
         ugrid = utrue
         vgrid = vtrue
      else
         print '(a)', 'unrecognized map projection.'
         stop 'truewind_to_gridwind'
      end if
      return
   end subroutine truewind_to_gridwind
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compare_projections(proj1, proj2, same_proj)

      ! subroutine to compare two proj_info structures to determine if the
      ! maps are the same.

      implicit none

      ! arguments
      type(proj_info), intent(in)   :: proj1
      type(proj_info), intent(in)   :: proj2
      logical, intent(out)          :: same_proj

      ! locals

      same_proj = .false.

      ! make sure both structures are initialized

      if ((.not. proj1%init) .or. (.not. proj2%init)) then
         print '(a)', 'compare_projections: map_set not called yet!'
         same_proj = .false.
      else
         same_proj = .true.
      end if

      ! check projection type
      if (same_proj) then
         if (proj1%code .ne. proj2%code) then
            print '(a)', 'compare_projections: different projection type.'
         else
            same_proj = .true.
         end if
      end if

      ! check corner lat/lon
      if (same_proj) then
         if ((abs(proj1%lat1 - proj2%lat1) .gt. 0.001) .or. &
             (abs(proj1%lon1 - proj2%lon1) .gt. 0.001)) then
            print '(a)', 'compare_projections: different lat1/lon1'
            same_proj = .false.
         end if
      end if

      ! compare dx
      if ((same_proj) .and. (proj1%code .ne. proj_latlon)) then
         if (abs(proj1%dx - proj2%dx) .gt. 1.) then
            print '(a)', 'compare_projections: different dx'
            same_proj = .false.
         end if
      end if
      ! compare dimensions
      if ((same_proj) .and. ((proj1%nx .ne. proj2%nx) .or. &
                             (proj1%ny .ne. proj2%ny))) then
         print '(a)', 'compare_projections: different dimensions'
         same_proj = .false.
      end if

      if ((same_proj) .and. (proj1%code .eq. proj_latlon)) then
         if ((proj1%dlat .ne. proj2%dlat) .or. &
             (proj1%dlon .ne. proj2%dlon)) then
            print '(a)', 'compare_projections: different dlat/dlon'
            same_proj = .false.
         end if
      end if

      ! compare stdlon for polar and lc projections
      if ((same_proj) .and. ((proj1%code .eq. proj_ps) .or. &
                             (proj1%code .eq. proj_lc))) then
         if (proj1%stdlon .ne. proj2%stdlon) then
            print '(a)', 'compare_projections: different stdlon.'
            same_proj = .false.
         end if
      end if
      ! compare true latitude1 if polar stereo, mercator, or lambert
      if ((same_proj) .and. ((proj1%code .eq. proj_ps) .or. &
                             (proj1%code .eq. proj_lc) .or. &
                             (proj1%code .eq. proj_merc))) then
         if (proj1%truelat1 .ne. proj2%truelat1) then
            print '(a)', 'compare_projections: different truelat1'
            same_proj = .false.
         end if
      end if

      ! compare true latitude2 if lc
      if ((same_proj) .and. (proj1%code .eq. proj_lc)) then
         if (proj1%truelat2 .ne. proj2%truelat2) then
            print '(a)', 'compare_projections: different truelat2'
            same_proj = .false.
         end if
      end if

      return
   end subroutine compare_projections
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compute_msf_lc(lat, truelat1, truelat2, msf)

      ! computes the map scale factor for a lambert conformal grid at a given
      ! latitude.

      implicit none
      real, intent(in)            :: lat  ! latitude where msf is requested
      real, intent(in)            :: truelat1
      real, intent(in)            :: truelat2
      real, intent(out)           :: msf

      real                        :: cone
      real                        :: psi1, psix, pole

      call lc_cone(truelat1, truelat2, cone)
      if (truelat1 .ge. 0.) then
         psi1 = (90.-truelat1)*rad_per_deg
         pole = 90.
      else
         psi1 = (90.+truelat1)*rad_per_deg
         pole = -90.
      end if
      psix = (pole - lat)*rad_per_deg
      msf = (sin(psi1)/sin(psix)) &
            *((tan(psix*0.5)/tan(psi1*0.5))**cone)
      return
   end subroutine compute_msf_lc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine compute_msf_ps(lat, truelat1, msf)

      ! computes the map scale factor for a polar stereographic grid at a given
      ! latitude.

      implicit none
      real, intent(in)            :: lat  ! latitude where msf is requested
      real, intent(in)            :: truelat1
      real, intent(out)           :: msf

      real                        :: psi1, psix, pole

      if (truelat1 .ge. 0.) then
         psi1 = (90.-truelat1)*rad_per_deg
         pole = 90.
      else
         psi1 = (90.+truelat1)*rad_per_deg
         pole = -90.
      end if
      psix = (pole - lat)*rad_per_deg
      msf = ((1.+cos(psi1))/(1.0 + cos(psix)))
      return
   end subroutine compute_msf_ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module map_utils
