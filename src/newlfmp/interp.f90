subroutine lfm_hinterp

   use lfmgrid
   use map_utils

   implicit none

   integer :: i, j, k
   real, allocatable, dimension(:, :) :: ri, rj
   logical l_bilinear

   l_bilinear = .true.

   allocate (ri(lx, ly), rj(lx, ly))

! compute i,j locations in native grid of output grid points.

   if (trim(mtype) == 'nmm') then
      do j = 1, ly
      do i = 1, lx
         call latlon_to_ij_nmm(llat(i, j), llon(i, j), ri(i, j), rj(i, j))
      end do
      end do
   else if (trim(mtype) == 'st4') then
      do j = 1, ly
      do i = 1, lx
! grid parameters for stage iv valid as of may 11, 2004
! see st4util.f90 for details
         call w3fb06(llat(i, j), llon(i, j) + 360., 23.117, 240.977, 4762.5, 255., &
                     ri(i, j), rj(i, j))
      end do
      end do
   else
      do j = 1, ly
      do i = 1, lx
         call latlon_to_ij(proj, llat(i, j), llon(i, j), ri(i, j), rj(i, j))
         ri(i, j) = max(1., min(float(nx), ri(i, j)))
         rj(i, j) = max(1., min(float(ny), rj(i, j)))
      end do
      end do
   end if

! horizontally interpolate all grids.

   do k = 1, nvar2d + nvar3d*nz
      if (ngrid(1, 1, k) < rmsg) then
         if (k <= nvar2d) then ! sfc grids
            if (trim(mtype) == 'nmm') then
               call hinterp_nmm(nx, ny, lx, ly, ri, rj, ngrid(1, 1, k), sgrid(1, 1, k))
            else
               do j = 1, ly
               do i = 1, lx
                  call gdtost_lfm(ngrid(1, 1, k), nx, ny, ri(i, j), rj(i, j), sgrid(i, j, k))
                  if (trim(mtype) == 'st4') then
                     if (sgrid(i, j, k) < 0) sgrid(i, j, k) = 0
                  end if
               end do
               end do
            end if
         else ! 3d grids
            if (trim(mtype) == 'nmm') then
               call hinterp_nmm(nx, ny, lx, ly, ri, rj, ngrid(1, 1, k), hgrid(1, 1, k - nvar2d))
            elseif (l_bilinear .and. k .ge. k_micro) then
               if (k .eq. k_micro) then
                  write (6, *) ' using bilinear interpolation for microphysical variables starting at k index ', k_micro
               end if
               call bilinear_laps_2d(ri, rj, nx, ny, lx, ly, ngrid(1, 1, k), hgrid(1, 1, k - nvar2d))
            else
               do j = 1, ly
               do i = 1, lx
                  call gdtost_lfm(ngrid(1, 1, k), nx, ny, ri(i, j), rj(i, j), hgrid(i, j, k - nvar2d))
               end do
               end do
            end if
         end if
      end if
   end do

   deallocate (ri, rj)

   return
end

!===============================================================================

subroutine lfm_vinterp(icall)

   use lfmgrid

   implicit none

   real :: lapse, gor, rog
   real, allocatable, dimension(:, :, :) :: logp
   parameter(lapse=0.0065, gor=9.8/287.04, rog=287.04/9.8)

   integer :: i, j, k, kk, icall
   integer :: i4_elapsed, ishow_timer, kkguess

   real :: pla, dz, plo, phi, slope

! generate a table of log pressure for future use.

   allocate (logp(lx, ly, nz))
   logp = alog(hpsig)

   if (.not. large_pgrid .and. icall .eq. 1) then

! interpolate 3d horizontally interpolated data to isobaric surfaces.
! assume that height and temp are always available, but check for
!  all missing in other fields.  no need to interpolate if all missing.

      do k = 1, lz
         pla = lprsl(k)
         kkguess = 1
         do j = 1, ly
         do i = 1, lx
            if (hpsig(i, j, 1) <= lprs(k)) then
               if (abs(hpsig(i, j, 1) - lprs(k)) < .01) then
                  dz = 0.
               else
                  dz = htsig(i, j, 1)/(gor/(pla - alog(hpsig(i, j, 1))) - lapse*0.5)
               end if
               if (.not. large_ngrid) zprs(i, j, k) = hzsig(i, j, 1) - dz
               tprs(i, j, k) = htsig(i, j, 1) + lapse*dz
            elseif (hpsig(i, j, nz) >= lprs(k)) then
               tprs(i, j, k) = htsig(i, j, nz)
               if (.not. large_ngrid) zprs(i, j, k) = hzsig(i, j, nz) + rog*htsig(i, j, nz)*alog(hpsig(i, j, nz)/lprs(k))
            else
               do kk = kkguess, kkguess
                  if (hpsig(i, j, kk) >= lprs(k) .and. hpsig(i, j, kk + 1) <= lprs(k)) then

                     plo = logp(i, j, kk)
                     phi = logp(i, j, kk + 1)
                     slope = (plo - pla)/(plo - phi)

                     if (.not. large_ngrid) zprs(i, j, k) = hzsig(i, j, kk) - slope*(hzsig(i, j, kk) - hzsig(i, j, kk + 1))
                     tprs(i, j, k) = htsig(i, j, kk) - slope*(htsig(i, j, kk) - htsig(i, j, kk + 1))
                     goto 10
                  end if
               end do
               do kk = 1, nz - 1
                  if (hpsig(i, j, kk) >= lprs(k) .and. hpsig(i, j, kk + 1) <= lprs(k)) then

                     plo = logp(i, j, kk)
                     phi = logp(i, j, kk + 1)
                     slope = (plo - pla)/(plo - phi)

                     if (.not. large_ngrid) zprs(i, j, k) = hzsig(i, j, kk) - slope*(hzsig(i, j, kk) - hzsig(i, j, kk + 1))
                     tprs(i, j, k) = htsig(i, j, kk) - slope*(htsig(i, j, kk) - htsig(i, j, kk + 1))
                     goto 10
                  end if
               end do
            end if
10          continue
         end do ! i
         end do ! j
      end do

      write (6, *) ' completed t,z vinterp'

      i4_elapsed = ishow_timer()

      if (minval(hmrsig) < rmsg) call vinterp_single(logp, hmrsig, shprs)

      if (.not. large_pgrid) then
         if (.true.) then ! linear
            if (minval(husig) < rmsg) call vinterp_single(logp, husig, uprs)
            if (minval(hvsig) < rmsg) call vinterp_single(logp, hvsig, vprs)
         else           ! cubic spline
            if (minval(husig) < rmsg) call vinterp_single_cubic(logp, husig, uprs)
            if (minval(hvsig) < rmsg) call vinterp_single_cubic(logp, hvsig, vprs)
         end if
         if (.not. large_ngrid) then
            if (minval(hwsig) < rmsg) call vinterp_single(logp, hwsig, wprs)
!     if (minval(htkesig) < rmsg) call vinterp_single(logp,htkesig,tkeprs)
         end if
      end if

   end if ! large_pgrid

   if (make_micro) then
      if (.not. large_pgrid) then
         if (icall .eq. 1) then
            if (minval(hcldliqmr_sig) < rmsg) call vinterp_single(logp, hcldliqmr_sig, cldliqmr_prs)
            if (minval(hcldicemr_sig) < rmsg) call vinterp_single(logp, hcldicemr_sig, cldicemr_prs)
            if (minval(hrainmr_sig) < rmsg) call vinterp_single(logp, hrainmr_sig, rainmr_prs)
            if (minval(hsnowmr_sig) < rmsg) call vinterp_single(logp, hsnowmr_sig, snowmr_prs)
            if (minval(hgraupelmr_sig) < rmsg) call vinterp_single(logp, hgraupelmr_sig, graupelmr_prs)
         elseif (icall .eq. 2) then
            if (minval(hzdr_sig) < rmsg) call vinterp_single(logp, hzdr_sig, zdr_prs)
            if (minval(hldr_sig) < rmsg) call vinterp_single(logp, hldr_sig, ldr_prs)
         end if
      end if
      if (icall .eq. 2) then
         if (minval(hrefl_sig) < rmsg) call vinterp_single(logp, hrefl_sig, refl_prs)
      end if
   end if

   deallocate (logp)

   return
end

!===============================================================================

subroutine vinterp_single(logp, sig, prs)

   use lfmgrid

   implicit none

   integer :: i, j, k, kk
   integer :: i4_elapsed, ishow_timer, kkguess
   real :: logp(lx, ly, nz), sig(lx, ly, nz), prs(lx, ly, lz), pla, plo, phi, slope

   do k = 1, lz
      pla = lprsl(k)
      kkguess = 1
      do j = 1, ly
      do i = 1, lx
         if (hpsig(i, j, 1) <= lprs(k)) then
            prs(i, j, k) = sig(i, j, 1)
         elseif (hpsig(i, j, nz) >= lprs(k)) then
            prs(i, j, k) = sig(i, j, nz)
         else
            do kk = kkguess, kkguess
               if (hpsig(i, j, kk) >= lprs(k) .and. hpsig(i, j, kk + 1) <= lprs(k)) then
                  plo = logp(i, j, kk)
                  phi = logp(i, j, kk + 1)
                  slope = (plo - pla)/(plo - phi)
                  prs(i, j, k) = sig(i, j, kk) - slope*(sig(i, j, kk) - sig(i, j, kk + 1))
                  goto 10
               end if
            end do
            do kk = 1, nz - 1
               if (hpsig(i, j, kk) >= lprs(k) .and. hpsig(i, j, kk + 1) <= lprs(k)) then
                  plo = logp(i, j, kk)
                  phi = logp(i, j, kk + 1)
                  slope = (plo - pla)/(plo - phi)
                  prs(i, j, k) = sig(i, j, kk) - slope*(sig(i, j, kk) - sig(i, j, kk + 1))
                  kkguess = kk
                  goto 10
               end if
            end do
         end if
10       continue
      end do ! i
      end do ! j
   end do

   write (6, *) ' completed vinterp_single'

   i4_elapsed = ishow_timer()

   return
end

!===============================================================================

subroutine vinterp_single_cubic(logp, sig, field_prs_grid)

   use lfmgrid

   implicit none

   integer :: i, j, k, kk, left, idebug
   real :: logp(lx, ly, nz), sig(lx, ly, nz), field_prs_grid(lx, ly, lz), pla, plo, phi, slope
   real :: ypp(nz), ypval, yppval
   real :: mlogp(nz), mlprsl(lz)
   logical :: l_sigma = .false. ! used for viewing wrf sigma coords for debugging

   write (6, *) ' running vinterp_single_cubic'
   write (6, *) ' nz,lprsl = ', nz, lprsl

   mlprsl(:) = -lprsl(:) ! laps grid

   do j = 1, ly
   do i = 1, lx
      idebug = (i*j) - 1
      mlogp(:) = -logp(i, j, :) ! native grid
      if (idebug .eq. 0) write (6, *) ' i,j,mlogp = ', i, j, lprs(:)
      if (idebug .eq. 0) write (6, *) ' i,j,mlogp = ', i, j, mlogp
      if (idebug .eq. 0) write (6, *) ' i,j,hpsig = ', i, j, hpsig(i, j, :)
      if (idebug .eq. 0) write (6, *) ' i,j,sig = ', i, j, sig(i, j, :)
      if (idebug .eq. 0) write (6, *) ' calling spline_cubic_set'
      call spline_cubic_set(nz, mlogp(:), sig(i, j, :), 0, 0., 0, 0., ypp)
      left = 1

      do k = 1, lz
         if (idebug .eq. 0) write (6, *) i, j, k, lprsl(k), mlprsl(k), logp(i, j, 1), mlogp(1), sig(i, j, 1)
         if (l_sigma) then
            if (idebug .eq. 0) write (6, *) ' setting to sigma value for debugging'
            kk = max((nz - (lz - k)), 1)
            field_prs_grid(i, j, k) = sig(i, j, kk)
         elseif (mlprsl(k) .ge. mlogp(1)) then ! laps pressure value at or above the native surface location
            if (idebug .eq. 0) write (6, *) ' calling spline_cubic_val2'
            call spline_cubic_val2(nz, mlogp(:), sig(i, j, :), ypp, left, mlprsl(k), field_prs_grid(i, j, k), ypval, yppval)
         else                              ! below the surface
            if (idebug .eq. 0) write (6, *) ' setting to surface value'
            field_prs_grid(i, j, k) = sig(i, j, 1)
         end if
         if (idebug .eq. 0) write (6, *) ' interpolated value = ', field_prs_grid(i, j, k)
      end do ! k

   end do ! i
   end do ! j

   return
end

!===============================================================================

subroutine gdtost_lfm(a, ix, iy, stax, stay, staval)

! subroutine to return stations back-interpolated values(staval)
!    from uniform grid points using overlapping-quadratics.
! gridded values of input array a dimensioned a(ix,iy),where
!    ix=grid points in x, iy = grid points in y.
! station location given in terms of grid relative station x (stax)
!    and station column.
! values greater than 1.0e30 indicate missing data.

   dimension a(ix, iy), r(4), scr(4)

   iy1 = int(stay) - 1
   iy2 = iy1 + 3
   ix1 = int(stax) - 1
   ix2 = ix1 + 3
   staval = 1e30
   fiym2 = float(iy1) - 1
   fixm2 = float(ix1) - 1
   ii = 0
   do i = ix1, ix2
      ii = ii + 1
      if (i >= 1 .and. i <= ix) then
         jj = 0
         do j = iy1, iy2
            jj = jj + 1
            if (j >= 1 .and. j <= iy) then
               r(jj) = a(i, j)
            else
               r(jj) = 1e30
            end if
         end do
         yy = stay - fiym2
         if (yy == 2.0) then
            scr(ii) = r(2)
         else
            call binom_lfm(1., 2., 3., 4., r(1), r(2), r(3), r(4), yy, scr(ii))
         end if
      else
         scr(ii) = 1e30
      end if
   end do
   xx = stax - fixm2
   if (xx == 2.0) then
      staval = scr(2)
   else
      call binom_lfm(1., 2., 3., 4., scr(1), scr(2), scr(3), scr(4), xx, staval)
   end if

   return
end

!===============================================================================

subroutine binom_lfm(x1, x2, x3, x4, y1, y2, y3, y4, xxx, yyy)

   yyy = 1e30
   if (x2 > 1.e19 .or. x3 > 1.e19 .or. &
       y2 > 1.e19 .or. y3 > 1.e19) return

   wt1 = (xxx - x3)/(x2 - x3)
   wt2 = 1.0 - wt1

   if (y4 < 1.e19 .and. x4 < 1.e19) then
      yz22 = (xxx - x3)*(xxx - x4)/((x2 - x3)*(x2 - x4))
      yz22 = wt1*(xxx - x4)/(x2 - x4)
      yz23 = (xxx - x2)*(xxx - x4)/((x3 - x2)*(x3 - x4))
      yz23 = wt2*(xxx - x4)/(x3 - x4)
      yz24 = (xxx - x2)*(xxx - x3)/((x4 - x2)*(x4 - x3))
   else
      yz22 = wt1
      yz23 = wt2
      yz24 = 0.0
   end if

   if (y1 < 1.e19 .and. x1 < 1.e19) then
      yz11 = (xxx - x2)*(xxx - x3)/((x1 - x2)*(x1 - x3))
      yz12 = (xxx - x1)*(xxx - x3)/((x2 - x1)*(x2 - x3))
      yz12 = wt1*(xxx - x1)/(x2 - x1)
      yz13 = (xxx - x1)*(xxx - x2)/((x3 - x1)*(x3 - x2))
      yz13 = wt2*(xxx - x1)/(x3 - x1)
   else
      yz11 = 0.0
      yz12 = wt1
      yz13 = wt2
   end if

   if (yz11 == 0. .and. yz24 == 0.) then
      yyy = wt1*y2 + wt2*y3
   else
      yyy = wt1*(yz11*y1 + yz12*y2 + yz13*y3) + wt2*(yz22*y2 + yz23*y3 + yz24*y4)
   end if

   return
end

!===============================================================================

subroutine latlon_to_ij_nmm(glatd, glond, ri, rj)

   use lfmgrid
   use constants

   implicit none

   real :: glatd, glond, ri, rj
   real :: dlm, dph, glat, glon, tlat, tlon, tph0, tlm0, x, y, z

   glat = glatd*deg2rad
   glon = glond*deg2rad
   dph = ngrid_spacingy*deg2rad
   dlm = ngrid_spacingx*deg2rad
   tph0 = ntruelat1*deg2rad
   tlm0 = nstdlon*deg2rad

   x = cos(tph0)*cos(glat)*cos(glon - tlm0) + sin(tph0)*sin(glat)
   y = -cos(glat)*sin(glon - tlm0)
   z = cos(tph0)*sin(glat) - sin(tph0)*cos(glat)*cos(glon - tlm0)
   tlat = rad2deg*atan(z/sqrt(x*x + y*y))
   tlon = rad2deg*atan(y/x)

   ri = float((nx + 1)/2) - (tlon/ngrid_spacingx/2.)
   rj = float((ny + 1)/2 + 1) + (tlat/ngrid_spacingy)

   return
end

!===============================================================================

subroutine w3fb06(alat, alon, alat1, alon1, dx, alonv, xi, xj)
!$$$   subprogram  documentation  block
!
! subprogram:  w3fb06        lat/lon to pola (i,j) for grib
!   prgmmr: stackpole        org: nmc42       date:88-04-05
!
! abstract: converts the coordinates of a location on earth given in
!   the natural coordinate system of latitude/longitude to a grid
!   coordinate system overlaid on a polar stereographic map pro-
!   jection true at 60 degrees n or s latitude. w3fb06 is the reverse
!   of w3fb07. uses grib specification of the location of the grid
!
! program history log:
!   88-01-01  original author:  stackpole, w/nmc42
!   90-04-12  r.e.jones   convert to cray cft77 fortran
!
! usage:  call w3fb06 (alat,alon,alat1,alon1,dx,alonv,xi,xj)
!   input argument list:
!     alat     - latitude in degrees (negative in southern hemis)
!     alon     - east longitude in degrees, real*4
!     alat1    - latitude  of lower left point of grid (point (1,1))
!     alon1    - longitude of lower left point of grid (point (1,1))
!                all real*4
!     dx       - mesh length of grid in meters at 60 deg lat
!                 must be set negative if using
!                 southern hemisphere projection.
!                   190500.0 lfm grid,
!                   381000.0 nh pe grid, -381000.0 sh pe grid, etc.
!     alonv    - the orientation of the grid.  i.e.,
!                the east longitude value of the vertical meridian
!                which is parallel to the y-axis (or columns of
!                of the grid)along which latitude increases as
!                the y-coordinate increases.  real*4
!                   for example:
!                   255.0 for lfm grid,
!                   280.0 nh pe grid, 100.0 sh pe grid, etc.
!
!   output argument list:
!     xi       - i coordinate of the point specified by alat, alon
!     xj       - j coordinate of the point; both real*4
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
!        radius to desired point and the i j too
!
   ala = alat*radpd
   rm = rebydx*cos(ala)*ss60/(1.+h*sin(ala))
!
   alo = (alon - reflon)*radpd
   xi = polei + rm*cos(alo)
   xj = polej + h*rm*sin(alo)
!
   return
end

!===============================================================================

subroutine hinterp_nmm(nx, ny, lx, ly, ri, rj, in, out)

   implicit none

   integer :: nx, ny, lx, ly, i, j, k
   real :: rih, riv, rjh
   real :: d1, d2, vbot, vtop
   integer, dimension(lx, ly) :: h1, h2, v1, v2, j1, j2
   real, dimension(lx, ly) :: ri, rj &
                              , hfact1, hfact2 &
                              , vfact1, vfact2 &
                              , jfact1, jfact2
   real, dimension(nx, ny) :: in
   real, dimension(lx, ly) :: out

   do j = 1, ly
   do i = 1, lx
      rih = max(1., min(float(nx), ri(i, j) + 0.5))
      riv = max(1., min(float(nx), ri(i, j)))
      rjh = max(1., min(float(ny), rj(i, j)))
      if (rih .eq. float(nx) .and. (rjh .gt. 1. .and. rjh .lt. float(ny))) then
         rih = ri(i, j) + 0.5
         riv = ri(i, j)
         rjh = rj(i, j)
         h1(i, j) = nx - 1
         h2(i, j) = h1(i, j) + 1
         j1(i, j) = min(ny - 1, int(rjh))
         j2(i, j) = j1(i, j) + 1
         d1 = rih - float(h2(i, j))
         jfact1(i, j) = float(j2(i, j)) - rjh
         jfact2(i, j) = 1.-jfact1(i, j)
         vbot = in(h2(i, j), j1(i, j)) &
                + d1*(in(h2(i, j), j1(i, j)) - in(h1(i, j), j1(i, j)))
         vtop = in(h2(i, j), j2(i, j)) &
                + d1*(in(h2(i, j), j2(i, j)) - in(h1(i, j), j2(i, j)))
         out(i, j) = jfact1(i, j)*vbot + jfact2(i, j)*vtop
!   else if (rih .eq. 1.0 .and. (rjh .gt. 1. .and. rjh .lt. float(ny)) ) then
!   rih=ri(i,j)+0.5
!    riv=ri(i,j)
!    rjh=rj(i,j)
!    h1(i,j)=2
!    h2(i,j)=1
!    j1(i,j)=min(ny-1,int(rjh))
!    j2(i,j)=j1(i,j)+1
!    d1 = rih-float(h2(i,j))
!    jfact1(i,j)=float(j2(i,j))-rjh
!    jfact2(i,j)=1.-jfact1(i,j)
!    vbot = in(h2(i,j),j1(i,j)) &
!       + d1*( in(h2(i,j),j1(i,j))-in(h1(i,j),j1(i,j)) )
!    vtop = in(h2(i,j),j2(i,j)) &
!       + d1*( in(h2(i,j),j2(i,j))-in(h1(i,j),j2(i,j)) )
!    out(i,j) = jfact1(i,j) * vbot + jfact2(i,j) * vtop
!   else if (rjh .eq. float(ny) .and. (rih .gt. 1. .and. rih .lt. float(nx)) ) then
!    rih=ri(i,j)+0.5
!    riv=ri(i,j)
!    rjh=rj(i,j)
!    h1(i,j)=int(rih)
!    h2(i,j)=h1(i,j)+1
!    j1(i,j)=ny-1
!    j2(i,j)=j1(i,j)+1
!    d1 = rjh-float(j2(i,j))
!    hfact1(i,j)=float(h2(i,j))-rih
!    hfact2(i,j)=1.-hfact1(i,j)
!    vbot = in(h2(i,j),j2(i,j)) &
!       + d1*( in(h2(i,j),j2(i,j))-in(h1(i,j),j2(i,j)) )
!    vtop = in(h1(i,j),j2(i,j)) &
!       + d1*( in(h1(i,j),j2(i,j))-in(h1(i,j),j1(i,j)) )
!    out(i,j) = hfact1(i,j) * vbot + hfact2(i,j) * vtop
!   else if (rjh .eq. 1.0 .and. (rih .gt. 1. .and. rih .lt. float(nx)) ) then
!    rih=ri(i,j)+0.5
!    riv=ri(i,j)
!    rjh=rj(i,j)
!    h1(i,j)=int(rih)
!    h2(i,j)=h1(i,j)+1
!    j1(i,j)=1
!    j2(i,j)=2
!    d1 = rjh-float(j2(i,j))
!    hfact1(i,j)=float(h2(i,j))-rih
!    hfact2(i,j)=1.-hfact1(i,j)
!    vbot = in(h2(i,j),j2(i,j)) &
!       + d1*( in(h2(i,j),j2(i,j))-in(h1(i,j),j2(i,j)) )
!    vtop = in(h1(i,j),j2(i,j)) &
!       + d1*( in(h1(i,j),j2(i,j))-in(h1(i,j),j1(i,j)) )
!    out(i,j) = hfact1(i,j) * vbot + hfact2(i,j) * vtop
      else
         h1(i, j) = min(nx - 1, int(rih))
         h2(i, j) = h1(i, j) + 1
         v1(i, j) = min(nx - 1, int(riv))
         v2(i, j) = v1(i, j) + 1
         j1(i, j) = min(ny - 1, int(rjh))
         j2(i, j) = j1(i, j) + 1
         hfact1(i, j) = float(h2(i, j)) - rih
         hfact2(i, j) = 1.-hfact1(i, j)
         vfact1(i, j) = float(v2(i, j)) - riv
         vfact2(i, j) = 1.-vfact1(i, j)
         jfact1(i, j) = float(j2(i, j)) - rjh
         jfact2(i, j) = 1.-jfact1(i, j)
         vbot = hfact1(i, j)*in(h1(i, j), j1(i, j)) &
                + hfact2(i, j)*in(h2(i, j), j1(i, j))
         vtop = hfact1(i, j)*in(h1(i, j), j2(i, j)) &
                + hfact2(i, j)*in(h2(i, j), j2(i, j))
         out(i, j) = jfact1(i, j)*vbot + jfact2(i, j)*vtop
      end if
   end do
   end do

   return
end
