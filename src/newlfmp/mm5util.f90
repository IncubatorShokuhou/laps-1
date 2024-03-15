module mm5util

   implicit none

   save

   integer :: bhi(50, 20)
   real :: bhr(20, 20)

end module

!===============================================================================

subroutine open_mm5(fname)

! open an mm5v3 file.

   implicit none

   character(len=*) :: fname
   open (1, file=trim(fname), form='unformatted', status='old', err=900)

   return

900 continue
   print *, 'could not open mm5 file: ', trim(fname)
   stop

end subroutine

!===============================================================================

subroutine get_mm5_dims(fname, nx, ny, nz)

   use mm5util

   implicit none

   integer :: nx, ny, nz, hdr_flag
   character(len=80) :: bhic(50, 20), bhrc(20, 20)
   character(len=*) :: fname
   logical :: opened

! open mm5 file, if not already opened.

   inquire (1, opened=opened)
   if (.not. opened) call open_mm5(fname)

! read big header information and fill grid dimensions.

   rewind (1)
   read (1, end=900) hdr_flag
   if (hdr_flag == 0) then
      read (1, end=900) bhi, bhr, bhic, bhrc
   else
      print *, 'get_mm5_dims: mm5 big header not found...'
      stop
   end if

! mm5 grid dimensions are specified for wind (dot) points,
!   but all fields will be specified at mass (cross) points.
!   so horizontal grid dimensions are set to one less.

   nx = bhi(17, 1) - 1
   ny = bhi(16, 1) - 1
   nz = bhi(12, 11)

   return

900 continue
   print *, 'get_mm5_dims: end of file reached...'
   stop

end subroutine

!===============================================================================

subroutine read_short_header(ndim, nx, ny, nz, name)

   implicit none

   integer :: ndim, nx, ny, nz, i1, i2, i3, i4, i5

   real :: xtime

   character(len=46) :: description
   character(len=25) :: units
   character(len=24) :: current_date
   character(len=9)  :: name
   character(len=4)  :: staggering, ordering

   read (1) ndim, i1, i2, i3, i4, nx, ny, nz, i5, xtime, staggering, ordering &
      , current_date, name, units, description
!xtime=xtime*60.

   return
end subroutine

!===============================================================================

subroutine fill_mm5_grid

   use lfmgrid
   use mm5util

   implicit none

   integer :: flag, ndim, fnx, fny, fnz, i, j, k
   real :: ptop
   real, allocatable, dimension(:) :: fld1d, sigh
   real, allocatable, dimension(:, :) :: fld2d, pstar, ncon_pcp_tot
   real, allocatable, dimension(:, :, :) :: fld3d, pp
   character(len=9) :: name
   logical :: there

! allocate local variables.

   allocate (sigh(nz), pstar(nx, ny), pp(nx, ny, nz), ncon_pcp_tot(nx, ny))
   ncon_pcp_tot = 0.
   npcp_tot = 0.

! read model data.

   read (1) flag
   do while (flag == 1)

      call read_short_header(ndim, fnx, fny, fnz, name)
      if (ndim == 3) then
         allocate (fld3d(fnx, fny, fnz))
         read (1) (((fld3d(i, j, k), i=1, fnx), j=1, fny), k=1, fnz)

         if (trim(name) == 'u') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call dot2crs(fnx, fny, fnz, fld3d)
            call fill_3dfield(nusig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'v') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call dot2crs(fnx, fny, fnz, fld3d)
            call fill_3dfield(nvsig, fld3d, nx, ny, nz)

         elseif (trim(name) == 't') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(ntsig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'q') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(nmrsig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'clw') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(ncldliqmr_sig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'rnw') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(nrainmr_sig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'ice') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(ncldicemr_sig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'snow') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(nsnowmr_sig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'graupel') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(ngraupelmr_sig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'w') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz + 1)
            call fill_2dfield(nwsfc, fld3d(1, 1, nz + 1), nx, ny)
            do k = 1, nz
            do j = 1, fny
            do i = 1, fnx
               fld3d(i, j, k) = (fld3d(i, j, k) + fld3d(i, j, k + 1))*0.5
            end do
            end do
            end do
            call fill_3dfield(nwsig, fld3d, nx, ny, nz)

         elseif (trim(name) == 'pp') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_3dfield(pp, fld3d, nx, ny, nz)

         end if
         deallocate (fld3d)

      elseif (ndim == 2) then
         allocate (fld2d(fnx, fny))
         read (1) ((fld2d(i, j), i=1, fnx), j=1, fny)

         if (trim(name) == 'pstarcrs') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(pstar, fld2d, nx, ny)

         elseif (trim(name) == 'ground t') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nground_t, fld2d, nx, ny)

         elseif (trim(name) == 't2') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            if (maxval(fld2d) <= 0.) then
               ntsfc = rmsg
            else
               call fill_2dfield(ntsfc, fld2d, nx, ny)
            end if

         elseif (trim(name) == 'q2') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            if (maxval(fld2d) <= 0.) then
               nmrsfc = rmsg
            else
               call fill_2dfield(nmrsfc, fld2d, nx, ny)
            end if

         elseif (trim(name) == 'u10') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            if (maxval(fld2d) <= 0.) then
               nusfc = rmsg
            else
               call fill_2dfield(nusfc, fld2d, nx, ny)
            end if

         elseif (trim(name) == 'v10') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            if (maxval(fld2d) <= 0.) then
               nvsfc = rmsg
            else
               call fill_2dfield(nvsfc, fld2d, nx, ny)
            end if

         elseif (trim(name) == 'rain non') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(npcp_tot, fld2d, nx, ny)

         elseif (trim(name) == 'rain con') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(ncon_pcp_tot, fld2d, nx, ny)

         elseif (trim(name) == 'terrain') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nzsfc, fld2d, nx, ny)

         elseif (trim(name) == 'latitcrs') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nlat, fld2d, nx, ny)

         elseif (trim(name) == 'longicrs') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nlon, fld2d, nx, ny)

         elseif (trim(name) == 'pbl hgt') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            if (maxval(fld2d) <= 1.) then
               npblhgt = rmsg
            else
               call fill_2dfield(npblhgt, fld2d, nx, ny)
            end if

         elseif (trim(name) == 'lwout') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nlwout, fld2d, nx, ny)

         elseif (trim(name) == 'swout') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nswout, fld2d, nx, ny)

         elseif (trim(name) == 'lwdown') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nlwdown, fld2d, nx, ny)

         elseif (trim(name) == 'swdown') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nswdown, fld2d, nx, ny)

         elseif (trim(name) == 'shflux') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nshflux, fld2d, nx, ny)

         elseif (trim(name) == 'lhflux') then
            call grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)
            call fill_2dfield(nlhflux, fld2d, nx, ny)

         end if
         deallocate (fld2d)

      elseif (ndim == 1) then
         allocate (fld1d(fnx))
         read (1) (fld1d(i), i=1, fnx)

         if (trim(name) == 'sigmah') then
            do k = 1, nz
               sigh(nz + 1 - k) = fld1d(k)
            end do

         end if

         deallocate (fld1d)
      else
         print *, 'unexpected ndim =', ndim
         stop
      end if
      read (1) flag
   end do
   close (1)

! generate 3d height, 3d pressure, and surface pressure.

   ptop = bhr(2, 2)*0.001
   call height_mm5(ntsig, nmrsig, nzsfc, pstar, pp, sigh, ptop, nx, ny, nz, nzsig)

   call pressure_mm5(nx, ny, nz, ptop, sigh, pstar, pp, npsfc, npsig)

! fill total precip and convert from cm to m.

   npcp_tot = (npcp_tot + ncon_pcp_tot)*0.01

   deallocate (pstar, pp, sigh, ncon_pcp_tot)

   if (fcsttime == 0.) then
      nlwout = rmsg
      nswout = rmsg
      nlwdown = rmsg
      nswdown = rmsg
      nshflux = rmsg
      nlhflux = rmsg
   end if

! fill native map projection settings.

   ngrid_spacingx = bhr(9, 1)
   ntruelat1 = bhr(5, 1)
   ntruelat2 = bhr(6, 1)
   nstdlon = bhr(3, 1)
   select case (bhi(7, 1))
   case (1)
      nprojection = 'lambert conformal'
   case (2)
      nprojection = 'polar stereographic'
   case (3)
      nprojection = 'mercator'
   end select

   return
end subroutine

!===============================================================================

subroutine grid_size_check(name, ndim, fnx, fny, fnz, nx, ny, nz)

   implicit none

   integer :: ndim, fnx, fny, fnz, nx, ny, nz

   character(len=9) :: name

   if (ndim == 1) return
   if (ndim == 2) then
      if (fnx /= ny + 1 .or. fny /= nx + 1) then
         print *, 'mm5 grid size failure - ', name
         print *, '   mm5 field size:', fny, fnx
         print *, '   mm5 post size :', nx, ny
         stop
      end if
   elseif (ndim == 3) then
      if (fnx /= ny + 1 .or. fny /= nx + 1 .or. fnz /= nz) then
         print *, 'mm5 grid size failure for - ', name
         print *, '   mm5 field size:', fny, fnx, fnz
         print *, '   mm5 post size :', nx, ny, nz
         stop
      end if
   end if

   return
end subroutine

!===============================================================================

subroutine fill_3dfield(var, fld, nx, ny, nz)

   implicit none

   integer :: nx, ny, nz, i, j, k, kk

   real :: var(nx, ny, nz), fld(ny + 1, nx + 1, nz)

   do k = 1, nz
      kk = nz + 1 - k
      do j = 1, ny
      do i = 1, nx
         var(i, j, k) = fld(j, i, kk)
      end do
      end do
   end do

   return
end subroutine

!===============================================================================

subroutine fill_2dfield(var, fld, nx, ny)

   implicit none

   integer :: nx, ny, i, j

   real :: var(nx, ny), fld(ny + 1, nx + 1)

   do j = 1, ny
   do i = 1, nx
      var(i, j) = fld(j, i)
   end do
   end do

   return
end subroutine

!===============================================================================

subroutine dot2crs(nx, ny, nz, field)

   implicit none

   integer :: nx, ny, nz, i, j, k
   real, dimension(nx, ny, nz) :: field

   do k = 1, nz
   do j = 1, ny - 1
   do i = 1, nx - 1
      field(i, j, k) = (field(i, j, k) + &
                        field(i + 1, j, k) + &
                        field(i, j + 1, k) + &
                        field(i + 1, j + 1, k))*0.25
   end do
   end do
   end do

   return
end subroutine

!===============================================================================

subroutine height_mm5(tp, mr, ter, ps, pph, sigh, ptop, nx, ny, nz, phih)

! code from mm5 to compute heights.

! input:   tp    temperature                cross    3d
!          mr    mixing ratio               cross    3d
!          ter   terrain                    cross    2d
!          ps    p* = psurf - ptop          cross    2d
!          sigh  sigma on half levels                1d
!          ptop  pressure at model lid
!          nx    dot point dimension w-e
!          ny    dot point dimension s-n

! stack:   phif  height on full levels      cross    3d
!          alnp  slab of pressure diff      cross    2d

! output:  phih  height on half levels      cross    3d

   implicit none

   real :: r, g, gi
   parameter(r=287.04, g=9.8, gi=1./g)

   integer :: nx, ny, nz, i, j, k

   real :: tp(nx, ny, nz), mr(nx, ny, nz), ter(nx, ny), ps(nx, ny), sigh(nz), ptop &
           , phif(nx, ny, nz + 1), alnp(nx, ny), sigf(nz + 1), phih(nx, ny, nz) &
           , pph(nx, ny, nz), ppf(nx, ny, nz + 1), denom, ptop_pa

   ptop_pa = ptop*1000.

! surface values of height on full levels.

   do j = 1, ny
   do i = 1, nx
      phif(i, j, 1) = ter(i, j)
   end do
   end do

   sigf(nz + 1) = 0.
   do k = nz, 1, -1
      sigf(k) = 2.*sigh(k) - sigf(k + 1)
   end do

! interpolate pressure perturbation to full levels.

   do k = 1, nz + 1
      if (k == 1) then
         do j = 1, ny
         do i = 1, nx
            ppf(i, j, 1) = pph(i, j, 1)
         end do
         end do
      elseif (k == nz + 1) then
         do j = 1, ny
         do i = 1, nx
            ppf(i, j, k) = pph(i, j, nz)
         end do
         end do
      else
         do j = 1, ny
         do i = 1, nx
            ppf(i, j, k) = 0.5*(pph(i, j, k - 1) + pph(i, j, k))
         end do
         end do
      end if
   end do

! integrate hydrostatic eqn.

   do k = 2, nz + 1
      do j = 1, ny
      do i = 1, nx
         alnp(i, j) = alog((sigf(k - 1)*ps(i, j) + ppf(i, j, k - 1) + ptop_pa) &
                           /(sigf(k)*ps(i, j) + ppf(i, j, k) + ptop_pa))*r
      end do
      end do

      do j = 1, ny
      do i = 1, nx
         phif(i, j, k) = phif(i, j, k - 1) + (tp(i, j, k - 1)*(1.+0.61*mr(i, j, k - 1))*alnp(i, j))*gi
      end do
      end do
   end do

! logarithmic interpolation to half-levels.

   do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         denom = 1./(sigf(k)*ps(i, j) + ppf(i, j, k) + ptop_pa)
         alnp(i, j) = alog((sigh(k)*ps(i, j) + pph(i, j, k) + ptop_pa)*denom) &
                      /alog((sigf(k + 1)*ps(i, j) + ppf(i, j, k + 1) + ptop_pa)*denom)
      end do
      end do

      do j = 1, ny
      do i = 1, nx
         phih(i, j, k) = phif(i, j, k) + (phif(i, j, k + 1) - phif(i, j, k))*alnp(i, j)
      end do
      end do
   end do

   return
end subroutine

!===============================================================================

subroutine pressure_mm5(nx, ny, nz, ptop, sigmah, ps, pp, psfc, p3)

! compute mm5 surface pressure and 3d pressure fields.

   implicit none

   integer :: nx, ny, nz, i, j, k

   real :: ptop &  !mm5 top pressure (kpa)          -input
           , sigmah(nz) &  !half sigma levels               -input
           , ps(nx, ny) &  !p star (pa)                     -input
           , pp(nx, ny, nz) &  !pressure perturbation (pa)      -input
           , psfc(nx, ny) &  !surface pressure (pa)           -output
           , p3(nx, ny, nz) &  !3-d pressure (pa)               -output
           , ptop_pa

   ptop_pa = ptop*1000.

! compute surface pressure.

   do j = 1, ny
   do i = 1, nx
      psfc(i, j) = ps(i, j) + pp(i, j, 1) + ptop_pa
   end do
   end do

! compute 3-d pressure at half sigma levels.

   do k = 1, nz
   do j = 1, ny
   do i = 1, nx
      p3(i, j, k) = sigmah(k)*ps(i, j) + pp(i, j, k) + ptop_pa
   end do
   end do
   end do

   return
end subroutine
