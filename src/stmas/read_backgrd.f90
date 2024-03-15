module read_backgrd
!*************************************************
! read in background informations
! history: january 2008, divorced from input_bg_obs module by zhongjie he.
!*************************************************

   use prmtrs_stmas
   use drawcountour, only: drcontour, drcontour_2d

   public tpg, bk0, c00, d00, x00, y00, t00, p00, z00, dx0, dy0, dz0, dt0, heightu, heightl
   public rdlapsbkg, rdbckgrnd, alloctbkg, dealctbkg, rdbkgtest

   real                  :: dx0, dy0
   real, allocatable :: x00(:, :), y00(:, :), t00(:), p00(:), z00(:), dz0(:, :, :)
   real, allocatable :: tpg(:, :), bk0(:, :, :, :, :), c00(:, :), d00(:, :)
   real, allocatable :: heightu(:, :), heightl(:, :)

! ===========added by zhongjie he just for draw pictures.
   public ori_lon, ori_lat, end_lon, end_lat
   real                  :: ori_lon, ori_lat, end_lon, end_lat
! =======================================================

!***************************************************
!!comment:
!   this module is mainly used by the module of input_bg_obs, the aim is to read background informations from models or some initial fields.
!   subroutines:
!      alloctbkg: memory allocate for x00, y00, p00, z00, tpg, bk0, c00, d00.
!      dealctbkg: memory deallocate for x00, y00, p00, z00, tpg, bk0, c00, d00.
!      rdlapsbkg: read in background from laps system.
!      rdbckgrnd: read in background from data files of 'fort.11'.
!      rdbkgtest: just a temporary subroutine to construction a background for a test case.
!   arrays:
!      tpg: topography
!      bk0: background fields on the original grid points (such as models or laps systems).
!      c00: coriolis force on the original gird.
!      d00: rotate angle's of the coordinate to the east-north coordinate. (guessed by zhongjie he)
!      x00: longitude of the original grid points.
!      y00: latitude of the original grid points.
!      p00: presures at each level
!      z00: height of each level
!      dx0: grid spacing in x direction, unit is meter
!      dy0: grid spacing in y direction, unit is meter
!   variables:
!      ori_lon: longitude of the east boundary of the area
!      ori_lat: latitude of the south boundary of the area
!      end_lon: longitude of the west boundary of the area
!      end_lat: latitude of the north boundary of teh area
contains

   subroutine alloctbkg
!*************************************************
! memory allocate for x00, y00, p00, z00, tpg, bk0, c00, d00
! history: september 2007, coded by wei li.
!          october   2007, by yuanfu xie (use laps)
!          january   2008, divorced frome memoryalc subroutine by zhongjie he
!*************************************************
      implicit none
! --------------------
      integer  :: er
! --------------------
      allocate (bk0(fcstgrd(1), fcstgrd(2), fcstgrd(3), fcstgrd(4), numstat + 3), stat=er)
      if (er .ne. 0) stop 'bk0 allocate wrong'
      allocate (x00(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) stop 'x00 allocate wrong'
      allocate (y00(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) stop 'y00 allocate wrong'
      allocate (tpg(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) stop 'tpg allocate wrong'
      allocate (c00(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) stop 'c00 allocate wrong'
      allocate (d00(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) stop 'd00 allocate wrong'
      allocate (p00(fcstgrd(3)), stat=er)
      if (er .ne. 0) stop 'p00 allocate wrong'
      allocate (z00(fcstgrd(3)), stat=er)
      if (er .ne. 0) stop 'z00 allocate wrong'
      allocate (dz0(fcstgrd(1), fcstgrd(2), fcstgrd(3)), stat=er)
      if (er .ne. 0) stop 'dz0 allocate wrong'
      allocate (heightu(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) stop 'heightu allocate wrong'
      allocate (heightl(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) stop 'heightl allocate wrong'
      allocate (t00(fcstgrd(4)), stat=er)
      if (er .ne. 0) stop 't00 allocate wrong'

      ! allocate space for laps lat/lon/topography grids (yuanfu):
      allocate (latitude(fcstgrd(1), fcstgrd(2)), longitud(fcstgrd(1), fcstgrd(2)), &
                topogrph(fcstgrd(1), fcstgrd(2)), stat=er)
      if (er .ne. 0) then
         print *, 'alloctbkg: cannot allocate memory for laps lat/lon/topography'
         stop
      end if

      return
   end subroutine alloctbkg

   subroutine dealctbkg
!*************************************************
! memory deallocate for x00, y00, p00, z00, tpg, bk0, c00, d00
! history: september 2007, coded by wei li.
!          october   2007, by yuanfu xie (use laps)
!          january   2008, separated frome memoryalc subroutine by zhongjie he
!*************************************************
      implicit none
! --------------------
      ! deallocate(bk0)                for keeping the background to the end by yuanfu
      deallocate (x00)
      deallocate (y00)
      deallocate (tpg)
      deallocate (c00)
      deallocate (d00)
      deallocate (p00)
      deallocate (z00)
      deallocate (dz0)
      deallocate (heightu)
      deallocate (heightl)
      deallocate (t00)

      ! deallocate laps variables:
      if (if_test .ne. 0) deallocate (latitude, longitud, topogrph)

      return
   end subroutine dealctbkg

   subroutine rdlapsbkg
!*************************************************
! read in background from laps system
! history: september 2007, coded by wei li.
!          october   2007, by yuanfu xie (use laps)
!
!          modified dec. 2013 by yuanfu reading in
!          surface pressure.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, ln
      character(len=256) :: dr
      character*4  :: lvl_coord
      character*10 :: unit
      character*125 :: comment
      real  :: ox, oy, ex, ey

      character*4  :: vn            ! varname added by yuanfu
      character*9  :: tm            ! ascii time added by yuanfu
      integer      :: st            ! status of laps calls added by yuanfu
      integer*4    :: i4            ! system i4 time added by yuanfu
      real         :: p1(fcstgrd(3)), ds(2)
      real :: d2r

      ! include a statement function for converting sh to 'rh' = sh/s2r(p) by yuanfu xie:
      include 'sh2rh.inc'

      d2r = 3.14159/180.0

!        'la' and 'lo' are respectively the latitude and longitude of each gridpoint
!        'p1' is the pressure of each level
! --------------------

      if (ifpcdnt .ne. 1) then              !   by zhongjiehe
         print *, 'error! laps is pressure coordinate! please set ifpcdnt to 1!'
         stop
      end if                              !   end of modified by zhongjie he

      ! laps cycle/system times:
      ! call get_systime(itime2(2),tm,st)        ! move to laps_config
      ! call get_laps_cycle_time(icycle,st)! move to laps_config

      ! stmas time window temporarily uses (-0.5cycle, 0.5*cycle):
      !itime2(1) = itime2(2)-(fcstgrd(4)-1)*icycle
      itime2(1) = lapsi4t - icycle/2
      itime2(2) = lapsi4t + icycle/2

      ! use laps ingest to read background fields: added by yuanfu

      ! currently only working for fcstgrd(4)=3:
      if (fcstgrd(4) .ne. 3) then
         print *, 'currently, stmas 3d has been tested only for 3 time frames!'
         print *, 'change the rdlapsbkg code for an analysis with different time frames!'
         stop
      end if

      ! background time:
      do t = 1, fcstgrd(4)
         t00(t) = (t - 1)/float(fcstgrd(4) - 1)*(itime2(2) - itime2(1))
      end do

      ! pressure levels:
      call get_pres_1d(i4, fcstgrd(3), p1, st)

      do t = 1, 2                        ! read in previous two time frames
         ! system time:
         i4 = lapsi4t + (t - 2)*icycle        ! use half laps time frames by yuanfu

         ! height field:
         vn = 'ht'
         call get_modelfg_3d(i4, vn, fcstgrd(1), fcstgrd(2), fcstgrd(3), bk0(1, 1, 1, t, pressure), st)
         if (st .ne. 1) then
            print *, 'rdlapsbkg: no ht background available at time frame: ', t
            stop
         end if
         ! temperature field:
         vn = 't3'
         call get_modelfg_3d(i4, vn, fcstgrd(1), fcstgrd(2), fcstgrd(3), bk0(1, 1, 1, t, temprtur), st)
         if (st .ne. 1) then
            print *, 'rdlapsbkg: no t3 background available at time frame: ', t
            stop
         end if
         ! wind u:
         vn = 'u3'
         call get_modelfg_3d(i4, vn, fcstgrd(1), fcstgrd(2), fcstgrd(3), bk0(1, 1, 1, t, u_cmpnnt), st)
         if (st .ne. 1) then
            print *, 'rdlapsbkg: no u3 background available at time frame: ', t
            stop
         end if
         ! wind v:
         vn = 'v3'
         call get_modelfg_3d(i4, vn, fcstgrd(1), fcstgrd(2), fcstgrd(3), bk0(1, 1, 1, t, v_cmpnnt), st)
         if (st .ne. 1) then
            print *, 'rdlapsbkg: no v3 background available at time frame: ', t
            stop
         end if
         ! specific humidity:
         vn = 'sh'
         call get_modelfg_3d(i4, vn, fcstgrd(1), fcstgrd(2), fcstgrd(3), bk0(1, 1, 1, t, humidity), st)
         if (st .ne. 1) then
            print *, 'rdlapsbkg: no sh background available at time frame: ', t
            if (numstat .gt. 4) stop
         end if
         print *, 'minvalue of laps sh: ', &
            minval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), t, humidity))

         ! convert sh to 'rh' = sh/s2r(p) by yuanfu xie:
         do i = 1, fcstgrd(1)
         do j = 1, fcstgrd(2)
         do k = 1, fcstgrd(3)
            bk0(i, j, k, t, humidity) = bk0(i, j, k, t, humidity)*1000.0/s2r(p1(k)/100.0) ! lga humidity: kg/kg
            if (p1(k) .le. 10000.0 .and. bk0(i, j, k, t, humidity) .gt. 1.0) &
               bk0(i, j, k, t, humidity) = 1.0  ! lga may have constant sh above 100hpa.
         end do
         end do
         end do

         ! surface pressure:
         if (t .eq. 2) then
            vn = 'ps'
            call get_directory('lsx', dr, ln)
            call read_laps(i4, i4, dr, 'lsx', fcstgrd(1), fcstgrd(2), 1, 1, vn, 0, &
                           lvl_coord, unit, comment, p_sfc_f, st)
         end if
      end do

      ! interpolate the background at previous and after time frames:
      ! interpolation:        previous time frame
      do i = 1, fcstgrd(1)
      do j = 1, fcstgrd(2)
      do k = 1, fcstgrd(3)
         bk0(i, j, k, 1, pressure) = 0.5*(bk0(i, j, k, 1, pressure) + bk0(i, j, k, 2, pressure))
         bk0(i, j, k, 1, temprtur) = 0.5*(bk0(i, j, k, 1, temprtur) + bk0(i, j, k, 2, temprtur))
         bk0(i, j, k, 1, u_cmpnnt) = 0.5*(bk0(i, j, k, 1, u_cmpnnt) + bk0(i, j, k, 2, u_cmpnnt))
         bk0(i, j, k, 1, v_cmpnnt) = 0.5*(bk0(i, j, k, 1, v_cmpnnt) + bk0(i, j, k, 2, v_cmpnnt))
         bk0(i, j, k, 1, humidity) = 0.5*(bk0(i, j, k, 1, humidity) + bk0(i, j, k, 2, humidity))
      end do
      end do
      end do
      ! extrapolation: after time frame
      do i = 1, fcstgrd(1)
      do j = 1, fcstgrd(2)
      do k = 1, fcstgrd(3)
         bk0(i, j, k, 3, pressure) = 1.5*bk0(i, j, k, 2, pressure) - 0.5*bk0(i, j, k, 1, pressure)
         bk0(i, j, k, 3, temprtur) = 1.5*bk0(i, j, k, 2, temprtur) - 0.5*bk0(i, j, k, 1, temprtur)
         bk0(i, j, k, 3, u_cmpnnt) = 1.5*bk0(i, j, k, 2, u_cmpnnt) - 0.5*bk0(i, j, k, 1, u_cmpnnt)
         bk0(i, j, k, 3, v_cmpnnt) = 1.5*bk0(i, j, k, 2, v_cmpnnt) - 0.5*bk0(i, j, k, 1, v_cmpnnt)
         bk0(i, j, k, 3, humidity) = amax1(0.0, 1.5*bk0(i, j, k, 2, humidity) - 0.5*bk0(i, j, k, 1, humidity))
      end do
      end do
      end do

      ! grid spacing:
      call get_grid_spacing_actual(latitude((fcstgrd(1) - 1)/2 + 1, (fcstgrd(2) - 1)/2 + 1), &
                                   longitud((fcstgrd(1) - 1)/2 + 1, (fcstgrd(2) - 1)/2 + 1), ds, st)

      dx0 = ds(1)
      dy0 = dx0 !ds(2)

      ! convert to real 8:
      y00 = latitude
      x00 = longitud
      tpg = topogrph
      p00 = p1

      ! end of yuanfu's modification

      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         d00(i, j) = 0.0d0
         c00(i, j) = 2.0*7.29e-5*sin(d2r*y00(i, j))
      end do
      end do

      ! modified by yuanfu for checking if rain and snow is analyzed:
      if (numstat .gt. 5) then
      do t = 1, fcstgrd(4)
      do k = 1, fcstgrd(3)
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         !added by shuyuan for ref  20100811
         bk0(i, j, k, t, rour_cmpnnt) = 0.0
         bk0(i, j, k, t, rous_cmpnnt) = 0.0

      end do
      end do
      end do
      end do
      end if

!===============
      do k = 1, fcstgrd(3)
         z_fcstgd(k) = p00(k)
      end do
      ! added by yuanfu for checking if vertical grid is uniform:
      uniform = 1    ! default uniform
      do k = 2, fcstgrd(3) - 1
         ! note: pressure coordinate is in pascal and 1 pascal is used as a threshold:
         if (abs(z_fcstgd(k) - z_fcstgd(k + 1) - z_fcstgd(1) + z_fcstgd(2)) .ge. 1.0) uniform = 0
      end do
!===============
      ox = x00(1, 1)
      ex = x00(fcstgrd(1), fcstgrd(2))
      oy = y00(1, 1)
      ey = y00(fcstgrd(1), fcstgrd(2))

      ori_lon = ox
      ori_lat = oy
      end_lon = ex
      end_lat = ey

      return
   end subroutine rdlapsbkg

   subroutine rdbkgtest
!*************************************************
! read in test data of background
! history: february 2008, coded by zhongjie he.
!*************************************************
      implicit none

      ! removed as it is not used.

   end subroutine rdbkgtest

end module read_backgrd
