module output_anals
!*************************************************
! output the analysis result
! history: september 2007, coded by wei li.
!          february 2008,  zhongjie he.
!
!          december 2013, yuanfu xie:
!          a) change ana dimension for saving memory
!          b) modify the tpw calculation
!*************************************************

   use prmtrs_stmas
   use generaltools

   use read_backgrd, only: ori_lon, ori_lat, end_lon, end_lat, bk0
   use drawcountour, only: drcontour, drcontour_2d
   use readobserves, only: x_radar, y_radar

   public outptlaps, outputana, tmpoutput
   private bkgmemrls, drcontour_0, radialwnd
!***************************************************
!!comment:
!   this module is used by the main program, to output the analized fields.
!   subroutines:
!      outptlaps : output the analized fields to laps.
!      outputana : output the analized fileds to some files and draw some pictures to check.
!      tmpoutput : output the analized fileds to some files and draw some pictures to check.
!      bkgmemrls : release the memorys.
!      drcontour_0: write some fields into the files with the format of surfer( a software to draw pictures).
!      radialwnd : translate the u, v fields to radial wind field, used to draw picture.

contains

   subroutine outptlaps
!*************************************************
! get background for analysis
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s, ln
      character(len=200) :: dr
      real     :: xb(maxgrid(1)), yb(maxgrid(2)), zb(maxgrid(3)), tb(maxgrid(4))
      real     :: xf(fcstgrd(1)), yf(fcstgrd(2)), zf(fcstgrd(3)), tf(fcstgrd(4))
      real, allocatable :: ana(:, :, :)
      integer  :: fg(maxdims), mg(maxdims)
      real     :: z1(1), t1(1), z2(1), t2(1), dx, dy

      ! parameters for converting q, p and t to rh (copied from lib/degrib/rrpr.f90:
      real, parameter :: svp1 = 611.2
      real, parameter :: svp2 = 17.67
      real, parameter :: svp3 = 29.65
      real, parameter :: svpt0 = 273.15
      real, parameter :: eps = 0.622
      real, parameter :: r_dry = 287.0
      real            :: tmp, sph, density

      ! variables needed for laps: yuanfu
      character*9   :: a9
      character*3   :: wn(3) = (/'u3', 'v3', 'om'/)        ! wind names
      character*3   :: sn(2) = (/'su', 'sv'/)                ! wind names
      character*4   :: wu(3) = (/'m/s ', 'm/s ', 'pa/s'/) ! wind units
      ! added by shuyuan 20100722
      character*3   :: qw(2) = (/'rai', 'sno'/)       ! rain water content   snow water content
      character*3   :: re(1) = (/'ref'/)                ! reflectivity  dbz
      character*125 :: qwc(2) = (/'rour', 'rous'/) ! qw comments
      character*125 :: rc(1) = (/'reflectivity'/) ! reflectivity comments
      character*10   :: units_3d(2) = (/'kg/m**3', 'kg/m**3'/)
      !----------------------------------------------------------
      character*125 :: wc(3) = (/'3dwind', '3dwind', '3dwind'/) ! wind comments
      character*125 :: sc(2) = (/'sfcwind', 'sfcwind'/) ! sfc wind comments
      integer       :: i4, st, iframe, ilv(fcstgrd(3))  ! pressure in mb needed using humid write routines
      ! real          :: la(fcstgrd(1),fcstgrd(2)),lo(fcstgrd(1),fcstgrd(2))
      ! real          :: tp(fcstgrd(1),fcstgrd(2)),ds,lv(fcstgrd(3))
      real          :: ds, lv(fcstgrd(3))

      ! yuanfu: make all large arrays into allocatable ones from automatic:
      real, allocatable :: ht(:, :, :), rh(:, :, :), t3(:, :, :), w3(:, :, :, :), sf(:, :, :), sp(:, :), tpw(:, :)

      real          :: height_to_zcoord3, ssh2, make_rh, make_ssh, rm, a, dlnp
!added by shuyuan 20100722 for reflectivity
      real          :: ref_out(fcstgrd(1), fcstgrd(2), fcstgrd(3))
      integer       :: istatus, n_3d_fields
!yuanfu test of cloud ice and liquid for temperature:
      character :: ext*31, unit*10, comment*30
      integer   :: i4_tol, i4_ret, iqc, nref, n2d, n3d
      integer   :: istatus2d(fcstgrd(1), fcstgrd(2)), istatus3d(fcstgrd(1), fcstgrd(2))
      integer   :: cloud_base, cloud_top, k400
      real :: closest_radar(fcstgrd(1), fcstgrd(2)), at
      real ::  rlat, rlon, rhgt
      real :: fraction

      ! yuanfu add integer arrays for interpolation dimensions:
      integer :: nanal(4), nbkgd(4)

      ! include a statement function for converting sh to 'rh' = sh/s2r(p) by yuanfu xie:
      include 'sh2rh.inc'

      i4_tol = 900
      i4_ret = 0
! --------------------

      ! allocate memory:
      allocate (ht(fcstgrd(1), fcstgrd(2), fcstgrd(3)), rh(fcstgrd(1), fcstgrd(2), fcstgrd(3)), &
                t3(fcstgrd(1), fcstgrd(2), fcstgrd(3)), w3(fcstgrd(1), fcstgrd(2), fcstgrd(3), 2), &
                sf(fcstgrd(1), fcstgrd(2), 2), sp(fcstgrd(1), fcstgrd(2)), tpw(fcstgrd(1), fcstgrd(2)), &
                stat=st)

      do i = 1, fcstgrd(1)
         xf(i) = (i - 1)*1.0d0
      end do
      do j = 1, fcstgrd(2)
         yf(j) = (j - 1)*1.0d0
      end do
      dx = ((fcstgrd(1) - 1)*1.0d0)/((maxgrid(1) - 1)*1.0d0)
      do i = 1, maxgrid(1)
         xb(i) = xf(1) + (i - 1)*dx
      end do
      dy = ((fcstgrd(2) - 1)*1.0d0)/((maxgrid(2) - 1)*1.0d0)
      do j = 1, maxgrid(2)
         yb(j) = yf(1) + (j - 1)*dy
      end do

!===========
      do k = 1, fcstgrd(3)
         zf(k) = z_fcstgd(k)
      end do
      do k = 1, maxgrid(3)
         zb(k) = z_maxgid(k)
      end do
!=================
!  open(2,file=dr(1:ln)//'p_level.dat',status='old',action='read')
!  do k=1,maxgrid(3)
!    read(2,*)zb(k)
!  enddo
!  close(2)
!==================
      do t = 1, fcstgrd(4)
         tf(t) = (itime2(2) - itime2(1))*(t - 1)/(fcstgrd(4) - 1)
      end do
      do t = 1, maxgrid(4)
         tb(t) = (itime2(2) - itime2(1))*(t - 1)/(maxgrid(4) - 1)
      end do

      ! output time frame: 2 current by yuanfu xie 2014 feb.
      iframe = 2
      nanal = maxgrid
      nanal(4) = 1 ! only output time frame
      nbkgd = fcstgrd
      nbkgd(4) = 1 ! only output time frame

      ! height from hydrostatic:
      ! write(15,*) maxgrid,grdbkgd0(1:maxgrid(1),1:maxgrid(2),1:maxgrid(3),2,3)
      ! do k=2,maxgrid(3)
      !   do j=1,maxgrid(2)
      !    do i=1,maxgrid(1)
      !      grdbkgd0(i,j,k,1:maxgrid(4),3) = grdbkgd0(i,j,k-1,1:maxgrid(4),3)-287.0/9.8*0.5* &
!           (grdbkgd0(i,j,k,1:maxgrid(4),4)+grdbkgd0(i,j,k-1,1:maxgrid(4),4))* &
      !        (log(ppm(k))-log(ppm(k-1)))
      !    enddo
      !  enddo
      ! enddo
      ! write(15,*) maxgrid,grdbkgd0(1:maxgrid(1),1:maxgrid(2),1:maxgrid(3),2,3)
      print *, 'increment 1: ', maxval(abs(grdbkgd0(1:maxgrid(1), 1:maxgrid(2), 1:maxgrid(3), iframe, 1)))
      print *, 'increment 2: ', maxval(abs(grdbkgd0(1:maxgrid(1), 1:maxgrid(2), 1:maxgrid(3), iframe, 2)))
      print *, 'increment 3: ', maxval(abs(grdbkgd0(1:maxgrid(1), 1:maxgrid(2), 1:maxgrid(3), iframe, 3)))
      print *, 'increment 4: ', maxval(abs(grdbkgd0(1:maxgrid(1), 1:maxgrid(2), 1:maxgrid(3), iframe, 4)))
      print *, 'increment 5: ', maxval(abs(grdbkgd0(1:maxgrid(1), 1:maxgrid(2), 1:maxgrid(3), iframe, 5)))
      if (numstat .gt. 5) print *, 'increment 6: ', maxval(abs(grdbkgd0(1:maxgrid(1), 1:maxgrid(2), maxgrid(3), iframe, raincont)))

      ! yuanfu: change ana to a 3 dimensional array to save space:
      ! allocate memory:
      allocate (ana(fcstgrd(1), fcstgrd(2), fcstgrd(3)), stat=st)

      ! pressure levels:
      call get_pres_1d(lapsi4t, fcstgrd(3), lv, st)
      ilv = lv/100.0   ! integer pressures needed in humid lh3 output

      ! loop through all analysis variables:
      do s = 1, numstat
         print *, '------------------------------------------------------'
         ana = 0.0
         if (uniform .eq. 0) then
            call bkgtofine(1, nanal, xb, yb, zb, tb(iframe), &
                           nbkgd, xf, yf, zf, tf(iframe), grdbkgd0(:, :, :, iframe, s), ana)
         else
            call uniform_interpolation3(maxgrid, fcstgrd, grdbkgd0(:, :, :, iframe, s), ana)
         end if
         ! call bkgtofine(1,maxgrid,xb,yb,zb,tb,fcstgrd,xf,yf,zf,tf,grdbkgd0(:,:,:,:,s),ana)

         ! make sure interpolated sh analysis positive:
         if (s .eq. humidity) then
            do k = 1, fcstgrd(3)
            do j = 1, fcstgrd(2)
            do i = 1, fcstgrd(1)

               ! convert 'rh' = sh/s2r(p) back to sh by yuanfu xie:
               ana(i, j, k) = &
                  max(-bk0(i, j, k, iframe, humidity), ana(i, j, k))*s2r(lv(k)/100.0)
               bk0(i, j, k, iframe, humidity) = bk0(i, j, k, iframe, humidity)*s2r(lv(k)/100.0)

            end do
            end do
            end do
         end if

         ! add increment ana to bk0:
         do k = 1, fcstgrd(3)
         do j = 1, fcstgrd(2)
         do i = 1, fcstgrd(1)
            bk0(i, j, k, iframe, s) = bk0(i, j, k, iframe, s) + ana(i, j, k)
         end do
         end do
         end do
         print *, 'bko_max=', maxval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), iframe, s)), s
         print *, '        ', maxval(ana)
         print *, 'bko_min=', minval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), iframe, s)), s
         print *, '        ', minval(ana)
      end do

      call get_grid_spacing_actual(latitude((fcstgrd(1) - 1)/2 + 1, (fcstgrd(2) - 1)/2 + 1), &
                                   longitud((fcstgrd(1) - 1)/2 + 1, (fcstgrd(2) - 1)/2 + 1), ds, st)

      call get_r_missing_data(rm, st)

      ! output: wind by yuanfu --

      ! interpolate surface wind:
      ! ht = bk0(:,:,:,1,3)
      ! w3 = bk0(:,:,:,1,1:2)
      ! do j=1,fcstgrd(2)
      ! do i=1,fcstgrd(1)
      !   sp(i,j) = height_to_zcoord3(tp(i,j),ht,lv,fcstgrd(1),fcstgrd(2),fcstgrd(3),i,j,st)
      !   k = int(sp(i,j))
      !   a = sp(i,j)-k
      !   sf(i,j,1:2) = (1.0-a)*w3(i,j,k,1:2)+a*w3(i,j,k+1,1:2)
      ! enddo
      ! enddo

      ! write surface wind (this is also reqired to plot 3d wind on-the-fly):
      ! call put_laps_multi_2d(i4,'lwm',sn,wu,sc,sf,fcstgrd(1),fcstgrd(2),2,st)
      ! call wind_post_process(i4,'lw3',wn,wu,wc,w3(1,1,1,1),w3(1,1,1,2), &
      !                        fcstgrd(1),fcstgrd(2),fcstgrd(3),3,sf(1,1,1),sf(1,1,2),tp, &
      !                        la,lo,ds,sp,rm,.true.,st)

      ! write temperature: w3 array serves as working array:
      ! t3(:,:,:) = bk0(:,:,:,1,4)
      ! call write_temp_anal(i4,fcstgrd(1),fcstgrd(2),fcstgrd(3),t3,ht,w3,a9,st)

      ! output time frame:
      iframe = 2        ! temporary testing by yuanfu xie

      ! use laps write balanced field:
      ! vertical velocity:
      t3 = 0.0 ! also test for x y boundaries: extrapolate later yuanfu
      ! center finite differnce require 2*delta x (or y):
      ds = 2.0*ds
      do k = 2, fcstgrd(3)
         do j = 2, fcstgrd(2) - 1
            do i = 2, fcstgrd(1) - 1
               ! v3(k)=v3(k-1)-dz*(ux+vy):
               t3(i, j, k) = t3(i, j, k - 1) - 0.5* &
                             (lv(k) - lv(k - 1))*( &
                             (bk0(i + 1, j, k, iframe, 1) - bk0(i - 1, j, k, iframe, 1))/ds + &
                             (bk0(i, j + 1, k, iframe, 2) - bk0(i, j - 1, k, iframe, 2))/ds)
            end do
         end do
      end do

      ! according a discussion with dan, specific humidity is adjusted by q_r (rain content) using ssh2 routine of laps:
      goto 111
      do k = 1, fcstgrd(3)
         do j = 1, fcstgrd(2)
            do i = 1, fcstgrd(1)
               ! assume saturation: td = t:
               if (bk0(i, j, k, iframe, raincont) .gt. 0.0) &
             bk0(i, j, k, iframe, 5) = ssh2(lv(k)/100.0, bk0(i, j, k, iframe, 4) - 273.15, bk0(i, j, k, iframe, 4) - 273.15, -132.0)
            end do
         end do
      end do
111   continue ! skip sh adjustment according to q_r

      ! adjust sh by bounds:
      goto 222
      do k = 1, fcstgrd(3)
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         bk0(i, j, k, iframe, 5) = max(bk0(i, j, k, iframe, 5), bk0(i, j, k, iframe, raincont))
      end do
      end do
      end do
222   continue

      ! converted from p, q t: note: q is in g/kg
      do k = 1, fcstgrd(3)
         do j = 1, fcstgrd(2)
            do i = 1, fcstgrd(1)
               ! sph = bk0(i,j,k,iframe,5)*0.001
               ! tmp = bk0(i,j,k,iframe,4)
               ! rh(i,j,k) = 1.e2 * (lv(k)*sph/(sph*(1.-eps) + eps))/(svp1*exp(svp2*(tmp-svpt0)/(tmp-svp3)))
               if (bk0(i, j, k, iframe, 5) .ge. 0.0) then
                  rh(i, j, k) = make_rh(lv(k)/100.0, bk0(i, j, k, iframe, 4) - 273.15, bk0(i, j, k, iframe, 5), -132.0)

                  ! ensure no rh greater than 1.0:
                  if (rh(i, j, k) .gt. 1.0) then
                     rh(i, j, k) = 1.0
                     bk0(i, j, k, iframe, 5) = &
                        make_ssh(lv(k)/100.0, bk0(i, j, k, iframe, 4) - 273.15, 1.0, -132.0)
                  end if

                  ! percentage:
                  rh(i, j, k) = rh(i, j, k)*100.0
               else
                  rh(i, j, k) = 0.0
               end if
            end do
         end do
      end do
      print *, 'max/min rh: ', maxval(rh), minval(rh)

      ! height from hydrostatic:
      ! bk0(1:fcstgrd(1),1:fcstgrd(2),1,iframe,3) = 0.0
      goto 1
      do k = 4, fcstgrd(3)
         do j = 1, fcstgrd(2)
            do i = 1, fcstgrd(1)
               bk0(i, j, k, iframe, 3) = bk0(i, j, k - 1, iframe, 3) - 287.0/9.8*0.5* &
                                         (bk0(i, j, k, iframe, 4) + bk0(i, j, k - 1, iframe, 4))* &
                                         (log(lv(k)) - log(lv(k - 1)))
            end do
         end do
      end do
1     continue

      ! temperature adjustment according to rain and snow: testing by yuanfu
      ! skip for now as cwb domain for morakot, this causes very large temp increment
      goto 11
      do k = 1, fcstgrd(3)
         do j = 1, fcstgrd(2)
            do i = 1, fcstgrd(1)
               ! dry density:
               density = lv(k)/r_dry/bk0(i, j, k, iframe, 4)
               ! adjusted temperature:
               bk0(i, j, k, iframe, 4) = lv(k)/r_dry/ &
                                         (density + bk0(i, j, k, iframe, raincont)*0.001 + bk0(i, j, k, iframe, snowcont)*0.001)
            end do
         end do
      end do
11    continue
      ! temperature adjustment according to cloud ice and liquid: test by yuanfu
      ext = "vrz"
      bk0(:, :, :, 1, numstat + 1) = 0.0      ! reflectivity
      call read_multiradar_3dref(lapsi4t, &
                                 900, 0, & ! 900 tolerate
                                 .true., -10.0, & ! apply_map: true; ref missing data value: 10
                                 fcstgrd(1), fcstgrd(2), fcstgrd(3), ext, latitude, longitud, topogrph, &
                                 .false., .false., & ! l_low_fill: false; l_high_fill: false
                                bk0(1, 1, 1, iframe, 3), bk0(1, 1, 1, 1, numstat + 1), rlat, rlon, rhgt, unit, iqc, closest_radar, &
                                 nref, n2d, n3d, istatus2d, istatus3d)
      ext = "lwc"
      bk0(:, :, :, iframe, numstat + 1) = 0.0 ! cloud liquid
      bk0(:, :, :, iframe, numstat + 2) = 0.0 ! cloud ice
      call get_laps_3dgrid(lapsi4t, i4_tol, i4_ret, &
                           fcstgrd(1), fcstgrd(2), fcstgrd(3), ext, "lwc", &
                           unit, comment, bk0(1, 1, 1, iframe, numstat + 1), st)
      call get_laps_3dgrid(lapsi4t, i4_tol, i4_ret, &
                           fcstgrd(1), fcstgrd(2), fcstgrd(3), ext, "ice", &
                           unit, comment, bk0(1, 1, 1, iframe, numstat + 2), st)
      print *, 'max cloud liquid: ', maxval(bk0(:, :, :, iframe, numstat + 1))
      print *, 'max cloud ice   : ', maxval(bk0(:, :, :, iframe, numstat + 2))
      goto 333
      ! adjust temperature according to cloud:
      do j = 1, fcstgrd(2)
         do i = 1, fcstgrd(1)
            at = 0.0
            cloud_base = fcstgrd(3) + 1
            cloud_top = 0
            do k = 1, fcstgrd(3)
               if (lv(k) .eq. 40000) then
                  if (bk0(i, j, k, 1, numstat + 1) .gt. 5 .and. &
                      bk0(i, j, k, 1, numstat + 1) .lt. 100.0) &
                     at = bk0(i, j, k, 1, numstat + 1)/10.0
                  k400 = k
               end if
               if (bk0(i, j, k, iframe, numstat + 1) .gt. 0.0 .or. &
                   bk0(i, j, k, iframe, numstat + 2) .gt. 0.0) then
                  cloud_base = min0(k, cloud_base)
                  cloud_top = max0(k, cloud_top)
               end if
            end do

            ! only adjust when cloud top is above 400mb
            if (cloud_top .lt. k400 .or. at .le. 0.0) cycle

            if (cloud_base .le. k400) then
               do k = cloud_base, k400
                  bk0(i, j, k, iframe, 4) = bk0(i, j, k, iframe, 4) + &
                                            at*(k - cloud_base)/float(max(k400 - cloud_base, 1))
               end do
            end if
            do k = k400 + 1, cloud_top
               bk0(i, j, k, iframe, 4) = bk0(i, j, k, iframe, 4) + &
                                         at*(cloud_top - k)/float(max(cloud_top - k400, 1))
            end do
         end do
      end do
333   continue

      ! total precipitable water:
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         tpw(i, j) = 0.0
         do k = 1, fcstgrd(3) - 1
            if (lv(k) .le. p_sfc_f(i, j)) then
               ! above topography: summed up
               tpw(i, j) = tpw(i, j) + 0.5* &
                           (bk0(i, j, k, iframe, 5) + bk0(i, j, k + 1, iframe, 5))* &
                           (lv(k) - lv(k + 1))/100.0 ! pressure in mb

            elseif (lv(k + 1) .le. p_sfc_f(i, j)) then
               tpw(i, j) = tpw(i, j) + bk0(i, j, k + 1, iframe, 5)* &
                           (p_sfc_f(i, j) - lv(k + 1))/100.0
            end if
         end do
         ! from g/kg to cm:
         tpw(i, j) = tpw(i, j)/100.0/9.8 ! following dan's int_ipw.f routine
      end do
      end do
      print *, 'tpw max/min: ', maxval(tpw), minval(tpw)

      ! kg/kg for balance netcdf:
      bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), iframe, 5) = 0.001* &
                                                                 bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), iframe, 5)
      call write_bal_laps(lapsi4t, bk0(1, 1, 1, iframe, 3), bk0(1, 1, 1, iframe, 1), &
                          bk0(1, 1, 1, iframe, 2), bk0(1, 1, 1, iframe, 4), &
                          t3, rh, bk0(1, 1, 1, iframe, 5), fcstgrd(1), fcstgrd(2), &
                          fcstgrd(3), lv, st)

      ! write temperature and wind into lapsprd:
      call write_temp_anal(lapsi4t, fcstgrd(1), fcstgrd(2), fcstgrd(3), &
                           bk0(1, 1, 1, iframe, 4), bk0(1, 1, 1, iframe, 3), dr, st)

      ! write wind by components:
      call put_laps_multi_3d_jacket(lapsi4t, 'lw3', wn, wu, wc, &
                                    bk0(1, 1, 1, iframe, 1), bk0(1, 1, 1, iframe, 2), t3, &
                                    fcstgrd(1), fcstgrd(2), fcstgrd(3), 3, st)

      ! write sp to file following dan's lq3_drive1a.f:
      write (*, *) 'write-file...'
      call writefile(lapsi4t, 'stmas sh analysis', ilv, bk0(1, 1, 1, iframe, 5), &
                     fcstgrd(1), fcstgrd(2), fcstgrd(3), st)

      ! write total precipitable water to lh4 file:
      write (*, *) 'write-lh4...'
      call write_lh4(lapsi4t, tpw, 1.0, fcstgrd(1), fcstgrd(2), st)

      ! write rh to lh3:
      write (*, *) 'write-lh3...'
      call lh3_compress(bk0(1, 1, 1, iframe, 5), bk0(1, 1, 1, iframe, 4), lapsi4t, ilv, &
                        -132.0, fcstgrd(1), fcstgrd(2), fcstgrd(3), 0, st)

      ! end of output to laps by yuanfu

! --------added by shuyuan 201007---------------------------------
! calculate and output reflectivity
      ! check if rain and snow is analyzed:
      if (numstat .le. 5) goto 555 ! by yuanfu
      do k = 1, fcstgrd(3)
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         ref_out(i, j, k) = 0.
         if (bk0(i, j, k, iframe, raincont) .gt. 0.0) then
            ref_out(i, j, k) = ref_out(i, j, k) + 43.1 + 17.5*alog10(bk0(i, j, k, iframe, raincont))
         else
            ! laps uses -10 as base value for reflectivity:
            ref_out(i, j, k) = -10.
         end if

         if (ref_out(i, j, k) .lt. 0.) then
            ! laps uses -10 as base value for reflectivity:
            ref_out(i, j, k) = -10.
         end if
      end do
      end do
      end do
      print *, 'ref_max=', maxval(ref_out(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3)))
      print *, 'ref_min=', minval(ref_out(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3)))
      call put_laps_3d(lapsi4t, 'lps', re, 'dbz', rc, ref_out(1, 1, 1), fcstgrd(1), fcstgrd(2), fcstgrd(3))

      !rain content(air density* rain water mixing ratio(kg/m3),
      !snow content(air density * snow water mixing ratio)
      n_3d_fields = 2   ! variable number
      istatus = 0
      !convert rc sc to rour rous
      do k = 1, fcstgrd(3)
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         if (bk0(i, j, k, iframe, raincont) .ne. 0) then
            bk0(i, j, k, iframe, raincont) = bk0(i, j, k, iframe, raincont)/1000.  ! kg/m3
         end if
         if (bk0(i, j, k, iframe, snowcont) .ne. 0) then
            bk0(i, j, k, iframe, snowcont) = bk0(i, j, k, iframe, snowcont)/1000.   !20100907
         end if
      end do
      end do
      end do
      print *, 'bko6_max=', maxval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), iframe, raincont))
      print *, 'bko7_max=', maxval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), iframe, snowcont))
      ! out put snow content(sno) and rai
      call put_laps_3d_multi_r(lapsi4t, 'lwc', qw, units_3d, qwc, &
                               bk0(1, 1, 1, iframe, raincont), bk0(1, 1, 1, iframe, snowcont), &
                               fcstgrd(1), fcstgrd(2), fcstgrd(3), &
                               fcstgrd(1), fcstgrd(2), fcstgrd(3), n_3d_fields, istatus)
!--end output ref rai sno  by shuyuan---------------------

      ! skip output rain and snow if not analyzed:
555   continue

      call bkgmemrls
      deallocate (z_fcstgd)
      ! finally release the background arrays: by yuanfu
      deallocate (bk0)
      ! finally release memory for lat/lon/topo:
      deallocate (latitude, longitud, topogrph)
      return
   end subroutine outptlaps

   subroutine outputana
!*************************************************
! get background for analysis
! history: september 2007, coded by wei li.
!          february 2008, by zhongjie he
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, l, t, s, ln
      character(len=200) :: dr
      real     :: xb(maxgrid(1)), yb(maxgrid(2)), zb(maxgrid(3)), tb(maxgrid(4))
      real     :: xf(fcstgrd(1)), yf(fcstgrd(2)), zf(fcstgrd(3)), tf(fcstgrd(4))
      real     :: bk1(fcstgrd(1), fcstgrd(2), fcstgrd(3), fcstgrd(4), numstat)
      integer  :: fg(maxdims), mg(maxdims)
      real     :: z1(1), t1(1), z2(1), t2(1), dx, dy, dt

! added by zhongjie he==========
      real  :: uu(fcstgrd(1), fcstgrd(2), fcstgrd(3))
      real  :: vv(fcstgrd(1), fcstgrd(2), fcstgrd(3))
      real  :: ww(fcstgrd(1), fcstgrd(2), fcstgrd(3))
      real  :: uv(fcstgrd(1), fcstgrd(2), fcstgrd(3))
      real  :: w0(fcstgrd(1), fcstgrd(2), fcstgrd(3), fcstgrd(4))

      real  :: d1, d2, dd
!===============================

! --------------------
      do i = 1, fcstgrd(1)
         xf(i) = (i - 1)*1.0d0
      end do
      do j = 1, fcstgrd(2)
         yf(j) = (j - 1)*1.0d0
      end do
      dx = ((fcstgrd(1) - 1)*1.0d0)/((maxgrid(1) - 1)*1.0d0)
      do i = 1, maxgrid(1)
         xb(i) = xf(1) + (i - 1)*dx
      end do
      dy = ((fcstgrd(2) - 1)*1.0d0)/((maxgrid(2) - 1)*1.0d0)
      do j = 1, maxgrid(2)
         yb(j) = yf(1) + (j - 1)*dy
      end do
      call get_directory('static', dr, ln)
!================
      do k = 1, fcstgrd(3)
         zf(k) = z_fcstgd(k)
      end do
      do k = 1, maxgrid(3)
         zb(k) = z_maxgid(k)
      end do
!=============== just for test data
!  open(2,file=dr(1:ln)//'p_level.dat',status='old',action='read')
!  do k=1,maxgrid(3)
!    read(2,*)zb(k)
!  enddo
!  close(2)

!  do k=1,maxgrid(3)
!    zb(k)=zf(1)+(k-1)*(zf(fcstgrd(3))-zf(1))/(maxgrid(3)-1.0)
!  enddo
! ===============

      do t = 1, fcstgrd(4)
         tf(t) = (t - 1)*1.0
      end do
      if (maxgrid(4) .ge. 2) dt = ((fcstgrd(4) - 1)*1.0)/(maxgrid(4) - 1.0)
      do t = 1, maxgrid(4)
         tb(t) = tf(1) + (t - 1)*dt
      end do

! added by zhongjie he==========

!  call fcst2bkgd(numdims,ngptobs,numstat,maxgrid,xb,yb,zb,tb,fcstgrd,xf,yf,zf,tf,difftout,bk1)
      call bkgtofine(numstat, maxgrid, xb, yb, zb, tb, fcstgrd, xf, yf, zf, tf, grdbkgd0, bk1)
!  call fcst2bkgd(numdims,ngptobs,1,maxgrid,xb,yb,zb,tb,fcstgrd,xf,yf,zf,tf,wwwout,w0)
      call bkgtofine(1, maxgrid, xb, yb, zb, tb, fcstgrd, xf, yf, zf, tf, wwwout, w0)

      deallocate (difftout)
      deallocate (wwwout)

      do k = 1, fcstgrd(3)
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         uu(i, j, k) = bk1(i, j, k, 1, u_cmpnnt)
         vv(i, j, k) = bk1(i, j, k, 1, v_cmpnnt)
!    ww(i,j,k)=w0(i,j,k,1)
!    d1=ori_lon+(end_lon-ori_lon)/float(fcstgrd(1)-1)*(i-1)
!    d1=d1-x_radar
!    d2=ori_lat+(end_lat-ori_lat)/float(fcstgrd(2)-1)*(j-1)
!    d2=d2-y_radar
!    dd=sqrt(d1*d1+d2*d2)
!    uv(i,j,k)=uu(i,j,k)*d1/dd+vv(i,j,k)*d2/dd
      end do
      end do
      end do

      call drcontour(ori_lon, end_lon, ori_lat, end_lat, numdims, fcstgrd, uu, 'ud.grd')
      call drcontour(ori_lon, end_lon, ori_lat, end_lat, numdims, fcstgrd, vv, 'vd.grd')
!  call drcontour(ori_lon,end_lon,ori_lat,end_lat,numdims,fcstgrd,ww,'wd.grd')
!  call drcontour(ori_lon,end_lon,ori_lat,end_lat,numdims,fcstgrd,uv,'rd.grd')

!===============================

!  call fcst2bkgd(numdims,ngptobs,numstat,maxgrid,xb,yb,zb,tb,fcstgrd,xf,yf,zf,tf,grdbkgd0,bk1)
      call bkgtofine(numstat, maxgrid, xb, yb, zb, tb, fcstgrd, xf, yf, zf, tf, grdbkgd0, bk1)

! added by zhongjie he==========
      do k = 1, fcstgrd(3)
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         uu(i, j, k) = bk1(i, j, k, 1, u_cmpnnt)
         vv(i, j, k) = bk1(i, j, k, 1, v_cmpnnt)
         ww(i, j, k) = w0(i, j, k, 1)
!    d1=ori_lon+(end_lon-ori_lon)/float(fcstgrd(1)-1)*(i-1)
!    d1=d1-x_radar
!    d2=ori_lat+(end_lat-ori_lat)/float(fcstgrd(2)-1)*(j-1)
!    d2=d2-y_radar
!    dd=sqrt(d1*d1+d2*d2)
!    uv(i,j,k)=uu(i,j,k)*d1/dd+vv(i,j,k)*d2/dd
      end do
      end do
      end do

      call drcontour(ori_lon, end_lon, ori_lat, end_lat, numdims, fcstgrd, uu, 'ua.grd')
      call drcontour(ori_lon, end_lon, ori_lat, end_lat, numdims, fcstgrd, vv, 'va.grd')
      call drcontour(ori_lon, end_lon, ori_lat, end_lat, numdims, fcstgrd, ww, 'wa.grd')
!  call drcontour(ori_lon,end_lon,ori_lat,end_lat,numdims,fcstgrd,uv,'ra.grd')
!===============================

      call bkgmemrls
      deallocate (z_fcstgd)
      deallocate (z_maxgid)

!============= to check the analysis field======
      open (1, file='stmas_anal.dat', action='write')
      write (1, *) fcstgrd(1:4)
      do l = 1, fcstgrd(4)
      do k = 1, fcstgrd(3)
      do j = 1, fcstgrd(2)
      do i = 1, fcstgrd(1)
         write (1, *) bk1(i, j, k, l, u_cmpnnt), bk1(i, j, k, l, v_cmpnnt), bk1(i, j, k, l, temprtur), w0(i, j, k, l)
      end do
      end do
      end do
      end do
      close (1)
!===========================================

      return
   end subroutine outputana

   subroutine tmpoutput
!*************************************************
! output analysis
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, s
      integer  :: ic, jc
! --------------------
      open (2, file='analysis.dat')
      do s = 1, numstat
         do t = 1, 1      ! numgrid(4)
         do k = 1, 1      ! numgrid(3)
         do j = 1, numgrid(2)
         do i = 1, numgrid(1)
            write (2, *) grdbkgd0(i, j, k, t, s)
         end do
         end do
         end do
         end do
      end do
      close (2)
!  call drcontour_0(pressure,'z.grd')
      call drcontour_0(u_cmpnnt, 'u.grd')
      call drcontour_0(v_cmpnnt, 'v.grd')
!  call drcontour_0(4,'t.grd')
      ic = (numgrid(1) + 1)/2
      jc = (numgrid(2) + 1)/2
      call radialwnd(ic, jc)
      call bkgmemrls
      return
   end subroutine tmpoutput

   subroutine bkgmemrls
!*************************************************
! release memory for background array
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
      deallocate (grdbkgd0)
      deallocate (xx0)
      deallocate (yy0)
      deallocate (cr0)
      deallocate (dg0)
      deallocate (dn0)
      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then
         deallocate (zz0)
         deallocate (zzb)
      elseif (ifpcdnt .eq. 1) then
         deallocate (pp0, ppm)
      end if
      return
   end subroutine bkgmemrls

   subroutine drcontour_0(s, fn)
!*************************************************
! draw the analysis (affiliate)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer   :: i, j, k, t, s, iu
      real      :: mx, mn
      character(len=5) :: fn
! --------------------
      k = 2 !(numgrid(3)+1)/2 !numgrid(3)
      t = (numgrid(4) + 1)/2 !numgrid(4) !1
      iu = 2
      open (iu, file=fn, status='unknown')
      mx = -1000.0
      mn = 1000.0
      do i = 1, numgrid(1)
      do j = 1, numgrid(2)
         if (grdbkgd0(i, j, k, t, s) .lt. 9000.0) then
            if (grdbkgd0(i, j, k, t, s) .gt. mx) mx = grdbkgd0(i, j, k, t, s)
            if (grdbkgd0(i, j, k, t, s) .lt. mn) mn = grdbkgd0(i, j, k, t, s)
         end if
      end do
      end do
      write (iu, '(a4)') 'dsaa'
      write (iu, '(2i4)') numgrid(1), numgrid(2)
      write (iu, *) oripstn(1), oripstn(1) + (numgrid(1) - 1)*grdspac(1)
      write (iu, *) oripstn(2), oripstn(2) + (numgrid(2) - 1)*grdspac(2)
      write (iu, *) mn, mx
      do i = 1, numgrid(1)
      do j = 1, numgrid(2)
         if (grdbkgd0(i, j, k, t, s) .gt. 9000.0) grdbkgd0(i, j, k, t, s) = 2.e38
      end do
      end do
      do j = 1, numgrid(2)
         write (iu, *) (grdbkgd0(i, j, k, t, s), i=1, numgrid(1))
      end do
      close (iu)
      return
   end subroutine drcontour_0

   subroutine radialwnd(ic, jc)
!*************************************************
! draw radial wind (affiliate)
! history: august 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: iu, ic, jc, i, j, k, t
      real     :: dx, dy, dd, mn, mx
      real     :: uv(numgrid(1), numgrid(2), numgrid(3), numgrid(4))
! --------------------
      do t = 1, numgrid(4)
      do k = 1, numgrid(3)
      do j = 1, numgrid(2)
      do i = 1, numgrid(1)
         if (i .eq. ic .and. j .eq. jc) then
            uv(i, j, k, t) = 2.e38
            cycle
         end if
         dx = 1.0d0*(i - ic)
         dy = 1.0d0*(j - jc)
         dd = sqrt(dx*dx + dy*dy)
         uv(i, j, k, t) = grdbkgd0(i, j, k, t, u_cmpnnt)*dx/dd + grdbkgd0(i, j, k, t, v_cmpnnt)*dy/dd
      end do
      end do
      end do
      end do
      k = (numgrid(3) + 1)/2 !numgrid(3)
      t = (numgrid(4) + 1)/2 !numgrid(4) !1
      iu = 2
      open (iu, file='radial.grd', status='unknown')
      mx = -1000.0
      mn = 1000.0
      do i = 1, numgrid(1)
      do j = 1, numgrid(2)
         if (uv(i, j, k, t) .lt. 9000.0) then
            if (uv(i, j, k, t) .gt. mx) mx = uv(i, j, k, t)
            if (uv(i, j, k, t) .lt. mn) mn = uv(i, j, k, t)
         end if
      end do
      end do
      write (iu, '(a4)') 'dsaa'
      write (iu, '(2i4)') numgrid(1), numgrid(2)
      write (iu, *) oripstn(1), oripstn(1) + (numgrid(1) - 1)*grdspac(1)
      write (iu, *) oripstn(2), oripstn(2) + (numgrid(2) - 1)*grdspac(2)
      write (iu, *) mn, mx
      do i = 1, numgrid(1)
      do j = 1, numgrid(2)
         if (uv(i, j, k, t) .gt. 9000.0) uv(i, j, k, t) = 2.e38
      end do
      end do
      do j = 1, numgrid(2)
         write (iu, *) (uv(i, j, k, t), i=1, numgrid(1))
      end do
      close (iu)
      return
   end subroutine radialwnd

!!added by shuyuan 20100729!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine put_laps_3d_multi_r(i4time, ext, var_3d &  ! in
                                  , units_3d, comment_3d &            ! in
                                  , array1, array2 &           ! in
                                  , nx_l_1, ny_l_1, nz_l_1 &             ! in
                                  , nx_l_2, ny_l_2, nz_l_2 &            ! in
                                  , n_3d_fields, istatus)              !in out

      real array1(nx_l_1, ny_l_1, nz_l_1) !variable for output
      real array2(nx_l_2, ny_l_2, nz_l_2) !variable for output

      character*125 comment_3d(n_3d_fields) !variable comment
      character*10 units_3d(n_3d_fields)    !variable unit
      character*3 var_3d(n_3d_fields)      !variable name
      character*(*) ext                    ! file name suffix and file derictory
      integer :: l

      write (6, *) ' subroutine put_laps_3d_multi_r...'

      l = 1
      call put_laps_multi_3d(i4time, ext, var_3d(l), units_3d(l), &
                             comment_3d(l), array1, nx_l_1, ny_l_1, nz_l_1, 1, istatus)
      if (istatus .ne. 1) return
      if (l .eq. n_3d_fields) return

      l = 2
      call put_laps_multi_3d_append(i4time, ext, var_3d(l), units_3d(l), &
                                    comment_3d(l), array2, nx_l_2, ny_l_2, nz_l_2, 1, istatus)
      if (istatus .ne. 1) return
      if (l .eq. n_3d_fields) return

      write (6, *) ' error: n_3d_fields exceeds limit ', n_3d_fields

      return
   end subroutine put_laps_3d_multi_r

!!!!from cloud_deriv_subs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine put_laps_multi_3d_append(i4time, ext, var_2d, units_2d, &
                                       comment_2d, field_3d, ni, nj, nk, nf, istatus)

      logical ltest_vertical_grid

      character*150 directory
      character*(*) ext

      character*125 comment_3d(nk*nf), comment_2d(nf)
      character*10 units_3d(nk*nf), units_2d(nf)
      character*3 var_3d(nk*nf), var_2d(nf)
      integer lvl_3d(nk*nf)
      character*4 lvl_coord_3d(nk*nf)

      real field_3d(ni, nj, nk, nf)

      istatus = 0

      call get_directory(ext, directory, len_dir)

      do l = 1, nf
         write (6, 11) directory, ext(1:5), var_2d(l)
11       format(' writing 3d ', a50, 1x, a5, 1x, a3)
      end do ! l

      do l = 1, nf
         do k = 1, nk

            iscript_3d = (l - 1)*nk + k

            units_3d(iscript_3d) = units_2d(l)
            comment_3d(iscript_3d) = comment_2d(l)
            if (ltest_vertical_grid('height')) then
               lvl_3d(iscript_3d) = zcoord_of_level(k)/10
               lvl_coord_3d(iscript_3d) = 'msl'
            elseif (ltest_vertical_grid('pressure')) then

               lvl_3d(iscript_3d) = nint(zcoord_of_level(k))/100
               lvl_coord_3d(iscript_3d) = 'hpa'
            else
               write (6, *) ' error, vertical grid not supported,' &
                  , ' this routine supports pressure or height'
               istatus = 0
               return
            end if

            var_3d(iscript_3d) = var_2d(l)

         end do ! k
      end do ! l

      call write_laps_multi(i4time, directory, ext, ni, nj, &
                            nk*nf, nk*nf, var_3d, lvl_3d, lvl_coord_3d, units_3d, &
                            comment_3d, field_3d, istatus)

      if (istatus .ne. 1) return

      istatus = 1

      return
   end subroutine put_laps_multi_3d_append
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module output_anals
