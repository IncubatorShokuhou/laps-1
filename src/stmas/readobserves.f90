module readobserves
!*************************************************
! read in radar observations
! history: january 2008, coded by zhongjie he.
!*************************************************
   use prmtrs_stmas
   use generaltools, only: vrtclpstn, vrtclpstn8, getobdate, interpltn, interpltn_xie, zpconvert, obstogrid
   use read_backgrd, only: bk0, x00, y00, p00, z00, t00, heightu, heightl

   private handleobs, ht_to_prs, lurao, lundx, lunew
   public rdradrobs, rdbufrobs, rdbufrobs_xie, allocatob, dealoctob, rdobstest
   public nobsmax, obp, obs, obe, nst, oba      !  , nallobs

   public x_radar, y_radar                      ! just for test

   integer, parameter :: lurao = 11, lundx = 21, lunew = 50, nobsmax = 10000000
   real, allocatable :: obp(:, :, :), obs(:, :), obe(:, :), oba(:, :)
   integer, allocatable :: nst(:)

!**************************************************
!comment:
!   this module is mainly used by input_bg_and_obs.f90.
!   subroutines:
!      allocatob: memory allocate for obp, obs, obe, oba and stt.
!      dealoctob: memory deallocate for obp, obs, obe, and stt.
!      rdradrobs: read in data of radial wind from laps.
!      rdbufrobs: read in conventional data from laps.
!      handleobs: determine the position of the observation in the background field.
!      ht_to_prs: translate the height of the observation to presure according to the background, used for the case of presure coordinate.
!      rdobstest: just a test subroutine for reading some ideal datas for a test case.
!
!   arrays:
!      obp: position of each observation in the background.
!      obs: value of observation.
!      obe: error of observation.
!      nst: number of each variable observation datas.
!      oba: azimuth and tilt angles of observation.
!
!   variable:
!       x_radar: the longitude of radar, just used by output_analysis.f90 to draw pictures in the test case
!       y_radar: the latitude of radar, just used by output_analysis.f90 to draw pictures in the test case
!**************************************************
contains

   subroutine allocatob
!*************************************************
! memory allocate for obp, obs, obe, oba and stt
! history: september 2007, coded by wei li
! history: january 2008, departed form the module 'input_bg_obs' by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer  :: s, er
! --------------------

      allocate (obp(numdims, nobsmax, numstat + 3), stat=er)
      if (er .ne. 0) stop 'obp allocate wrong'
      allocate (obs(nobsmax, numstat + 3), stat=er)
      if (er .ne. 0) stop 'obs allocate wrong'
      allocate (obe(nobsmax, numstat + 3), stat=er)
      if (er .ne. 0) stop 'obe allocate wrong'
      allocate (nst(numstat + 3), stat=er)
      if (er .ne. 0) stop 'nst allocate wrong'
      do s = 1, numstat + 3
         nst(s) = 0
      end do
      allocate (oba(nobsmax, 3), stat=er)
      if (er .ne. 0) stop 'oba allocate wrong'

      ! allocate grid space for laps radial wind:
      ! allocate(idxradwn(4,fcstgrd(1)*fcstgrd(2)*fcstgrd(3)*fcstgrd(4)), &
      !            lapsradw(fcstgrd(1)*fcstgrd(2)*fcstgrd(3)*fcstgrd(4)), stat=er)
      ! if (er .ne. 0) then
      !   print*,'size too large: ',fcstgrd(1)*fcstgrd(2)*fcstgrd(3)*fcstgrd(4)
      !   print*,'size too large: ',fcstgrd(1),fcstgrd(2),fcstgrd(3),fcstgrd(4)
      !   print*,'allocatob: cannot allocate memory for laps radial wind'
      !   stop
      ! endif

      return
   end subroutine allocatob

   subroutine dealoctob
!*************************************************
! memory deallocate for obp, obs, obe, and stt
! history: september 2007, coded by wei li
! history: january 2008, departed form the module 'input_bg_obs' by zhongjie he.
!*************************************************
      implicit none
! --------------------

      deallocate (obp)
      deallocate (obs)
      deallocate (obe)
      deallocate (nst)
      deallocate (oba)

      ! deallocate laps radial wind:
      ! deallocate(idxradwn, lapsradw)

      return
   end subroutine dealoctob

   subroutine rdlapsrdr
!*************************************************
!  this routine reads gridded radar data from laps
!  remapped dataset using laps:
!               get_multiradar_vel and
!               qc_radar_obs
!       assuming the background has been read in.
!
!  history: march 2008, coded by yuanfu xie.
!*************************************************

      use mem_namelist              ! laps wind parameter module

      implicit none

      ! local variables:
      character*31 :: radext(max_radars)    ! possible radar name extensions
      character*4  :: radnam(max_radars)    ! radar station names
      character*8  :: sid                        ! radar station name in 8 byte
      integer      :: nradar                ! number of radar available
      integer      :: sttrad, sttnqy         ! radar and its nyquist status
      integer      :: ngrdrd(max_radars)    ! number of gridpoints with measurable vel
      integer      :: radtim(max_radars)    ! radar observation time
      integer      :: radids(max_radars)    ! radar idx
      integer      :: ioffset(max_radars), joffset(max_radars) ! offset for the new get_multiradar_vel
      integer      :: i, j, k, l, m, ix0, iy0, ix1, iy1        ! grid indices, number of radar, time frame
      logical      :: cluttr                ! .true. -- remove 3d radar clutter
      real         :: radvel(fcstgrd(1), fcstgrd(2), fcstgrd(3), max(3, max_radars))  ! 3 to allow cloud liquid read in
      ! radar 4d velocity grid
      real         :: radnqy(fcstgrd(1), fcstgrd(2), fcstgrd(3), max(3, max_radars))  ! 3 to allow cloud ice read in
      ! radar 4d nyquist velocity
      real         :: uvzero(fcstgrd(1), fcstgrd(2), fcstgrd(3), 3)    ! increased to 3 sharing with cloud read
      ! zero uv grids used for calling laps qc_radar_obs
      real         :: volnqy(max_radars)            ! volume nyquist velocity
      real         :: radlat(max_radars), radlon(max_radars), radhgt(max_radars)
      real         :: uvgrid(2)

      integer      :: toltim, inc        ! radar data tolerate window, time increment
      real         :: xradar, yradar, zradar  ! radar station location in analysis grid
      real         :: height_to_zcoord2     ! laps routine for converting height to grid
      real         :: xspace, yspace ! spacings

      ! variables for using handle obs routine:
      integer      :: ip
      real         :: op(4), prslvl(fcstgrd(3)), az, ea, ob, oe
      real         :: height_of_level, oh, sr                ! oh obs height, sr radar slant range

      !====== used by the modification of zhongjie he
      integer      :: i0, j0        ! index of the grid point at the west-south corner to radar station
      real         :: xrdr, yrdr    ! longitude and latitude of radar station used to calculate az and ea
      real         :: re           ! radius of earth
      !====== end modification of zhongjie he
!jhui
      integer :: tt, tt1, inc1, nn
      real    :: refscl
      real    :: lat(fcstgrd(1), fcstgrd(2))
      real    :: lon(fcstgrd(1), fcstgrd(2))
      real    :: topo(fcstgrd(1), fcstgrd(2))
      real    :: rlaps_land_frac(fcstgrd(1), fcstgrd(2))
      real    :: grid_spacing_cen_m
      integer :: istatus, i4_tol, i4_ret, iqc_2dref
      character :: units*10, comment*125, radar_name*4, iext*31, c_filespec*255
      real :: heights_3d(fcstgrd(1), fcstgrd(2), fcstgrd(3))
      real :: radar_ref_3d(fcstgrd(1), fcstgrd(2), fcstgrd(3), fcstgrd(4))
      real :: closest_radar(fcstgrd(1), fcstgrd(2))
      real ::  rlat_radar(5), rlon_radar(5), rheight_radar(5)
      integer :: n_ref_grids, n_2dref, n_3dref
!  integer :: istat_radar_2dref_a,istat_radar_3dref_a
!!variable's defination are not same with fountion read_multiradar_3dref ,modified by shuyuan 20100525
      integer :: istat_radar_2dref_a(fcstgrd(1), fcstgrd(2))
      integer :: istat_radar_3dref_a(fcstgrd(1), fcstgrd(2))
!! liu 20100525
      character*40 c_vars_req
      character*180 c_values_req
      integer :: i4time_radar
      real     :: tempref, make_ssh, make_td, tw, make_rh

! add the following definition since build on 5/23/2011 failed. hj 5/23/2011
      integer nx_r, ny_r

      include 'main_sub.inc'
      include 'laps_cloud.inc'
      real :: cld_prs(kcloud)
      real :: cloud3d(fcstgrd(1), fcstgrd(2), kcloud)

      if (ifpcdnt .ne. 1) then              !   by zhongjiehe
         print *, 'error! laps is pressure coordinate! please set ifpcdnt to 1!'
         stop
      end if                              !   end of modified by zhongjie he

      ! cluttr = .true.               ! true. -- remove 3d radar clutter for old get_multiradar_vel
      toltim = 600                ! default 10 minute tolerate time window
      re = 6365.0e3                   ! added by zhongjie he
      if (fcstgrd(4) .gt. 1) toltim = t00(fcstgrd(4)) - t00(1)
      ! get pressure levels:
      call get_pres_1d(lapsi4t, fcstgrd(3), prslvl, sttrad)

      ! read radar data close to each time frame:
      inc = toltim/(fcstgrd(4) - 1)

      do m = 1, fcstgrd(4)
         ! get unfolded radar through laps:
         call get_multiradar_vel(itime2(1) + (m - 1)*inc, inc/2, radtim, max_radars, &
                                 nradar, radext, rmissing, &
                                 fcstgrd(1), fcstgrd(2), fcstgrd(3), &
                                 latitude, longitud, &
                                 fcstgrd(1), fcstgrd(2), 0, &
                                 radvel, radnqy, radids, volnqy, &
                                 ioffset, joffset, cluttr, ngrdrd, &
                                 radlat, radlon, radhgt, radnam, sttrad, sttnqy)

         ! set uv zero grids:
         uvzero = 0.0

         ! nyquist unfolding:
         do l = 1, nradar
            ! qc and unfolding radar nyquist:
            nx_r = fcstgrd(1)
            ny_r = fcstgrd(2)
            ioffset = 0
            joffset = 0
            call qc_radar_obs(fcstgrd(1), fcstgrd(2), fcstgrd(3), &
                              rmissing, nx_r, ny_r, ioffset, joffset, radvel(1, 1, 1, l), radnqy(1, 1, 1, l), ngrdrd(l), &
                              latitude, longitud, radlat(l), radlon(l), radhgt(l), &
                              uvzero(1, 1, 1, 1), uvzero(1, 1, 1, 2), bk0(1, 1, 1, 1, 1), bk0(1, 1, 1, 1, 2), &
                              volnqy(l), l_correct_unfolding, l_grid_north, sttrad)
            print *, 'status of qc_radar_obs: ', sttrad, ' for time frame: ', m

            ! assign azimuth and elevation angles: bk0 stores uvzt
            call latlon_to_rlapsgrid(radlat(l), radlon(l), latitude, longitud, &
                                     fcstgrd(1), fcstgrd(2), xradar, yradar, sttrad)

            ! yuanfu found zradar is not being used: turn off for now 07/2009
            ! for radar stations outside domain, this caused problem as
            ! height_to_zcoord2 uses height(xradar,yradar...):
            !zradar = height_to_zcoord2(radhgt(l),bk0(1,1,1,1,3), &
            !                           fcstgrd(1),fcstgrd(2),fcstgrd(3), &
            !                           nint(xradar),nint(yradar),sttrad)

            ! laps grid spacing:
            call get_grid_spacing_actual_xy(radlat(l), radlon(l), xspace, yspace, sttrad)

            do k = 1, fcstgrd(3)
               ! note laps remap routine currently uses standard atmosphere
               ! for obs height. when it changes, the height fed to
               ! latlon_to_radar has to changed:
               oh = height_of_level(k)
               do j = 1, fcstgrd(2)
                  do i = 1, fcstgrd(1)

                     if (radvel(i, j, k, l) .ne. rmissing) then

!               print*,'rdlapsrdr: --radial wind: ', &
!                 radvel(i,j,k,l),radnqy(i,j,k,l),i,j,k,l,xradar,yradar,ngrdrd(l),volnqy(l)

                        ! compute azimuth and elevation angles using laps routine
                        ! latlon_to_radar.
                        call latlon_to_radar(latitude(i, j), longitud(i, j), oh, az, sr, ea, &
                                             radlat(l), radlon(l), radhgt(l))

                        op(2) = latitude(i, j)      ! op(1): longitude; op(2): latitude
                        op(1) = longitud(i, j)
                        op(3) = prslvl(k)          ! in pascal

                        ! for given tolerated time window, assume available radar is
                        ! at the time frame:
                        op(4) = radtim(l) - itime2(1)        ! actual radar time
                        ! op(4) = t00(m)

                        ip = 1
                        ob = radvel(i, j, k, l)                ! positive wind is toward the station by yuanfu
                        oe = 0.5
!jhui
!               oe=0.3
                        sid(1:4) = radnam(l)

                        call handleobs_sigma(op, ob, oe, numstat + 1, nallobs, ip, az, ea, sid)

                     end if
                  end do        ! i
               end do                ! j
            end do                ! k
         end do                ! l -- radar levels

      end do                ! m -- time frames
      print *, 'number of radar radial wind taken as obs: ', nst(numstat + 1)

      call get_laps_domain_95(fcstgrd(1), fcstgrd(2), lat, lon, topo &
                              , rlaps_land_frac, grid_spacing_cen_m, istatus)
      if (istatus .ne. 1) then
         write (6, *) ' error getting laps domain'
         return
      end if
      write (6, *) ' actual grid spacing in domain center = ', grid_spacing_cen_m

!=======reflectivity==for time cycle ,read multitime data file *.vrz======
!=========added by shuyuan liu 20100830==================

      i4_tol = 900
      i4_ret = 0
      ref_base = -10
      nn = 0
      refscl = 0.0
      iext = "vrz"
      bk0(:, :, :, :, numstat + 1) = 0.0
      radar_ref_3d = 0.0
      ! get radar at current and previous time:
      do l = 1, 2  !for l   time
         call get_filespec(iext(1:3), 2, c_filespec, istatus)
         call get_file_time(c_filespec, lapsi4t + (l - 2)*icycle, i4time_radar)

         call read_multiradar_3dref(lapsi4t + (l - 2)*icycle, i4_tol, i4_ret, &!i
                                    .true., ref_base, &                              ! i
                                    fcstgrd(1), fcstgrd(2), fcstgrd(3), iext, &   ! i
                                    lat, lon, topo, .false., .false., heights_3d, &  ! i
                                    radar_ref_3d(1, 1, 1, l), &                                      ! o
                                    rlat_radar, rlon_radar, rheight_radar, radar_name, &  ! o
                                    iqc_2dref, closest_radar, &                          ! o
                                    n_ref_grids, n_2dref, n_3dref, istat_radar_2dref_a, &  ! o
                                    istat_radar_3dref_a)                              ! o
         do k = 1, fcstgrd(3)
            do j = 1, fcstgrd(2)
               do i = 1, fcstgrd(1)
               if (radar_ref_3d(i, j, k, l) .gt. 5. .and. radar_ref_3d(i, j, k, l) .lt. 100) then
                  ! modified shuyuan 20100719

                  tempref = (radar_ref_3d(i, j, k, l) - 43.1)/17.5
                  tempref = (10.**tempref)       !g/m3
                  ! shuyuan 20100719
                  refscl = refscl + tempref**2

                  ! check if rain and snow is analyzed:
                  if (numstat .le. 5) goto 555
                  nn = nn + 1
                  op(2) = latitude(i, j)      ! op(1): longitude; op(2): latitude
                  op(1) = longitud(i, j)
                  op(3) = prslvl(k)          ! in pascal
!           op(4) = radtim(l)-itime2(1)    ! actual radar time
                  op(4) = t00(l)
                  ob = tempref
                  oe = 0.01  ! shuyuan   test 0.1 0.01 1
                  sid(1:3) = "vrz"
                  call handleobs_sigma(op, ob, oe, numstat + 3, nallobs, ip, az, ea, sid)
                  ! skip reflectivity:
555               continue

                  ! modified by yuanfu xie nov. 2011 for adding radar reflectivity generated sh obs:
                  !nst(humidity) = nst(humidity)+1
                  !obp(1,nst(humidity),humidity) = i-1
                  !obp(2,nst(humidity),humidity) = j-1
                  !obp(3,nst(humidity),humidity) = k-1
                  !obp(4,nst(humidity),humidity) = l-1
                  !obs(nst(humidity),humidity) = make_ssh(prslvl(k)/100.0,bk0(i,j,k,l,temprtur)-273.15,0.75,-132.0)
                  !obe(nst(humidity),humidity) = 1.0
                  !nallobs = nallobs+1
               end if
               end do
            end do
         end do
      end do  ! for l
      ! interpolation to the three time frames of stmas analysis:
      radar_ref_3d(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 3) = &
         1.5*radar_ref_3d(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) - &
         0.5*radar_ref_3d(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)
      radar_ref_3d(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) = &
         0.5*radar_ref_3d(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) + &
         0.5*radar_ref_3d(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)

      ! calculate the sh bounds from radar reflectivity:
      do l = 1, fcstgrd(4)
         do k = 1, fcstgrd(3)
            do j = 1, fcstgrd(2)
               do i = 1, fcstgrd(1)
                  if (radar_ref_3d(i, j, k, l) .gt. 5. .and. radar_ref_3d(i, j, k, l) .lt. 100) then

                     ! assume 75% satured rh where reflectivity present as sh lower bounds:
                     if (bk0(i, j, k, l, temprtur) .gt. 273.15) then
                        if (radar_ref_3d(i, j, k, l) .gt. 45.0) then
                           bk0(i, j, k, l, numstat + 1) = &
                              make_ssh(prslvl(k)/100.0, bk0(i, j, k, l, temprtur) - 273.15, 1.0, 0.0)
                        else
                           bk0(i, j, k, l, numstat + 1) = &
                         make_ssh(prslvl(k)/100.0, bk0(i, j, k, l, temprtur) - 273.15, 0.1 + 0.9*radar_ref_3d(i, j, k, l)/45.0, 0.0)
                        end if
                     else
                        if (radar_ref_3d(i, j, k, l) .gt. 30.0) then
                           bk0(i, j, k, l, numstat + 1) = &
                              make_ssh(prslvl(k)/100.0, bk0(i, j, k, l, temprtur) - 273.15, 1.0, 0.0)
                        else
                           bk0(i, j, k, l, numstat + 1) = &
                         make_ssh(prslvl(k)/100.0, bk0(i, j, k, l, temprtur) - 273.15, 0.2 + 0.8*radar_ref_3d(i, j, k, l)/30.0, 0.0)
                        end if
                     end if

                  end if
               end do
            end do
         end do
      end do
      print *, 'range of reflectivity derived bound: ', &
         maxval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2, numstat + 1)), &
         minval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2, numstat + 1))

      ! interpolate radref to grdbkgd0 as sh lower bound:
      if (maxgrid(3) .ne. fcstgrd(3) .or. maxgrid(4) .ne. fcstgrd(4)) then
         print *, 'the analysis grid does not match multigrid in z or t, please check!'
         print *, 'the sh lower bound calculation assumes they are the same'
         stop
      end if

      ! read in laps cloud liquid and ice:
      iext = "lwc"
      do l = 1, 2
         call get_laps_3dgrid(lapsi4t + (l - 2)*icycle, i4_tol, i4_ret, &
                              fcstgrd(1), fcstgrd(2), fcstgrd(3), iext, iext, &
                              units, comment, radvel(1, 1, 1, l), sttrad)
      end do
      ! interpolation to the three time frames of stmas analysis:
      radvel(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 3) = &
         1.5*radvel(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) - &
         0.5*radvel(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)
      radvel(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) = &
         0.5*radvel(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) + &
         0.5*radvel(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)
      print *, 'laps cloud liquid: ', minval(radvel(:, :, :, 2))*1000.0, maxval(radvel(:, :, :, 2))*1000.0
      do l = 1, 2
         call get_laps_3dgrid(lapsi4t + (l - 2)*icycle, i4_tol, i4_ret, &
                              fcstgrd(1), fcstgrd(2), fcstgrd(3), iext, "ice", &
                              units, comment, radnqy(1, 1, 1, l), sttrad)
      end do
      ! interpolation to the three time frames of stmas analysis:
      radnqy(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 3) = &
         1.5*radnqy(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) - &
         0.5*radnqy(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)
      radnqy(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) = &
         0.5*radnqy(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) + &
         0.5*radnqy(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)
      print *, 'laps cloud ice: ', minval(radnqy(:, :, :, 2))*1000.0, maxval(radnqy(:, :, :, 2))*1000.0

      ! read in laps cloud analysis for bounds of sh: use of array uzero
      uvzero = 0.0
      iext = "lc3"
      do l = 1, 2
         call get_clouds_3dgrid(lapsi4t + (l - 2)*icycle, j, fcstgrd(1), fcstgrd(2), &
                                kcloud, iext, cloud3d, cld_hts, cld_prs, i)
         if (j .ne. lapsi4t + (l - 2)*icycle) then
            print *, 'lc3 time does not match the analysis, skip'
            cycle
         else
            print *, 'readobservs: found file -- ', j
         end if

         do k = 1, fcstgrd(3)
            do m = 1, kcloud - 1
               if (prslvl(k) .le. cld_prs(m) .and. (prslvl(k) .ge. cld_prs(m + 1))) then
                  xspace = (log(cld_prs(m)) - log(prslvl(k)))/(log(cld_prs(m)) - log(cld_prs(m + 1)))
                  yspace = 1.0 - xspace
                  do j = 1, fcstgrd(2)
                  do i = 1, fcstgrd(1)
                     uvzero(i, j, k, l) = yspace*cloud3d(i, j, m) + xspace*cloud3d(i, j, m + 1)
                  end do
                  end do
               end if
            end do
         end do
      end do
      ! interpolation to the three time frames of stmas analysis:
      uvzero(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 3) = &
         1.5*uvzero(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) - &
         0.5*uvzero(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)
      uvzero(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) = &
         0.5*uvzero(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) + &
         0.5*uvzero(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1)

      ! convert cloud to sh bounds:
      do l = 1, 3
         do k = 1, fcstgrd(3)
         do j = 1, fcstgrd(2)
         do i = 1, fcstgrd(1)
            ! xspace = 0.0
            ! if (uvzero(i,j,k,l) .gt. 0.1) then
            !   xspace = make_ssh(prslvl(k)/100.0,bk0(i,j,k,l,temprtur)-273.15, &
            !                     0.8*(uvzero(i,j,k,l)**0.2),0.0)
            ! endif
            ! if (xspace .gt. bk0(i,j,k,l,numstat+1)) bk0(i,j,k,l,numstat+1)=xspace

            ! adjust sh bounds based on cloud fraction:
            if ((uvzero(i, j, k, l) .lt. 0.1) .and. (radar_ref_3d(i, j, k, l) .gt. 5.0)) then
               ! xspace = make_td(prslvl(k)/100.0,bk0(i,j,k,l,temprtur)-273.15, & ! dew
               !                  bk0(i,j,k,l,numstat+1),0.0)
               ! yspace = tw(bk0(i,j,k,l,temprtur)-273.15,xspace,prslvl(k)/100.0) ! wet bulb t
               ! xspace = make_rh(prslvl(k)/100.0,bk0(i,j,k,l,temprtur)-273.15,& ! rh
               !                  bk0(i,j,k,l,numstat+1),0.0)
               bk0(i, j, k, l, numstat + 1) = 0.0 !make_ssh(prslvl(k)/100.0,yspace,xspace,0.0)
            end if
         end do
         end do
         end do
      end do
      print *, 'range of reflectivity and cloud derived bound: ', &
         maxval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2, numstat + 1)), &
         minval(bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2, numstat + 1))

      ! assign reflectivity derived bounds to grdbkgd0:
      do i = 1, maxgrid(1)
         ix0 = float(i - 1)/float(maxgrid(1) - 1)*(fcstgrd(1) - 1) + 1
         ix1 = min(ix0 + 1, fcstgrd(1))
         do j = 1, maxgrid(2)
            iy0 = float(j - 1)/float(maxgrid(2) - 1)*(fcstgrd(2) - 1) + 1
            iy1 = min(iy0 + 1, fcstgrd(2))
            do k = 1, maxgrid(3)
               do l = 1, maxgrid(4)
                  ! simple shift instead of interpolation:
                  grdbkgd0(i, j, k, l, numstat + 1) = 0.25*(bk0(ix0, iy0, k, l, numstat + 1) + &
                             bk0(ix1, iy0, k, l, numstat + 1) + bk0(ix0, iy1, k, l, numstat + 1) + bk0(ix1, iy1, k, l, numstat + 1))
                  grdbkgd0(i, j, k, l, numstat + 2) = 1000.0 !grdbkgd0(i,j,k,l,numstat+1) ! test upper bound
               end do
            end do
         end do
      end do
      print *, 'max dbz over finest grid: ', maxval(grdbkgd0(:, :, :, :, numstat + 1)), &
         minval(grdbkgd0(:, :, :, :, numstat + 1))
   end subroutine rdlapsrdr

   subroutine rdradrobs
!*************************************************
! read in radar observations
! history: january 2008, coded by zhongjie he.
!*************************************************
      implicit none
! --------------------
      integer, parameter :: nh = 12, nr = 7, ml = 1500
      character(len=8)       :: ss
      integer         :: l, nl, ir, dt, nc, nw, o, os, is, ip
      real            :: hd(nh)
      real            :: ra(nr, ml)
      real            :: x, y, p, t, uv, zz, tt, ue, ms, ea, az, ze
      real            :: op(numdims), ob, oe

      real            :: x0, y0, t0, sh, aa, id, td
      integer         :: yr, mn, dy, hr
      integer         :: ii

      integer, external :: ireadsb, ireadmg, i4dy

!  character*150          :: od
!  character*9            :: a9
!  integer                :: ld,i4,n4,st

      character*80            :: hdrstr, datstr

      integer          :: fg

      data hdrstr/'sstn clat clon selv anaz anel year mnth days hour minu mgpt'/
      data datstr/'stdm suplat suplon heit rwnd rwaz rstd'/

!  call get_systime(i4,a9,st)
!  call get_filespec('bufr',2,od,st)
!  call get_file_time(od,i4,n4)
!  call s_len(od,ld)
!  call make_fnam_lp(n4,a9,st)
!  od(ld+4:ld+9) = od(ld-4:ld)
!  od(ld-5:ld+3) = a9
!  ld = ld+9

      print *, 'reading radar bufr file...........'
      open (unit=lurao, file='./radarii.bufr', form='unformatted' &
            , status='old', action='read')
      call openbf(lurao, 'in', lurao)

      fg = 0
      nc = 0
      nw = 0
!  call get_config(is)
!  if(is.ne.1)stop 'laps parameters are wrong!!!'
      o = nallobs
      do while (ireadmg(lurao, ss, dt) .eq. 0)

!    write(6,*) 'read_radar: bufr file date is ',dt

         do while (ireadsb(lurao) .eq. 0)

            if (fg .eq. 1) exit             ! just for test, by zhongjie he

!     read header. extract station information
            call ufbint(lurao, hd, nh, 1, ir, hdrstr)

!      iyr = hd(7)
!      imo = hd(8)
!      idy = hd(9)
!      ihr = hd(10)
!      imn = hd(11)
!      isc = izero
            t0 = hd(10)*60 + hd(11) ! base time   (hours)
            x0 = hd(3)            ! station lon (degrees)
            y0 = hd(2)            ! station lat (degrees)
            sh = hd(4)            ! station elevation
            aa = hd(5)            ! azimuth of radia (degrees)
            ea = hd(6)            ! elevation angle    (degrees)
            if (x0 .lt. 0.0) x0 = x0 + 360.0
            if (x0 .ge. 360.0) x0 = x0 - 360.0

!     go through the data levels
            call ufbint(lurao, ra, nr, ml, nl, datstr)

! ===============just for test by zhongjie he
            print *, 'lon=', x0, 'lat=', y0
            if (x0 - 150 .ge. 118 .and. x0 - 150 .le. 124 .and. y0 - 20 .ge. 21. .and. y0 - 20 .le. 25 .and. nl .ge. 100) then
               fg = 1
            else
               cycle
            end if

            x_radar = x0 - 150
            y_radar = y0 - 150
!===============

            if (nl > ml) then
               write (6, *) 'read_radar:  ***error*** increase read_radar bufr size since', &
                  'number of levs=', nl, ' > maxlevs=', ml
               stop
            end if

            do l = 1, nl
               nw = nw + 1

! for t location
               t = (t0 + ra(1, l))/60.         ! unit is hours
! for x location
               x = ra(3, l)

               x = x - 150                    ! just for a test
!        print*, 'x=x-150 is just for test ----------------------'

               if (x .lt. 0.0) x = x + 360.0
               if (x .ge. 360.0) x = x - 360.0
! for y location
               y = ra(2, l)

               y = y - 20                     ! just for a test
!        print*, 'y=y-20 is just for test ---------------------'

! for zz observation
               zz = ra(4, l)
! for azimuth of wind
               az = ra(6, l)

! transform
               op(1) = x
               op(2) = y
               ip = 1
               call ht_to_prs(op(1), op(2), zz, p, is)
               if (is .ne. 1) cycle
               op(3) = p
! for radial wind observation
               uv = ra(5, l)
               ue = ra(7, l)
! output
               ze = 2.0
               ob = uv
               oe = ue
               os = numstat + 1
               call handleobs(op, ob, oe, os, o, ip, az, ea)

               nc = nc + 1

            end do
         end do
      end do
      nallobs = o
      print *, 'the number of radar data is', nw, 'and', nc, 'available'
      call closbf(lurao)
      close (lurao)
!  close(lundx)
!  close(lunew)
      print *, 'nallobs=', nallobs
      return
   end subroutine rdradrobs

   subroutine rdbufrobs_xie
!*************************************************
! read convensional obs from a bufr file, replace
! the old rdbufrobs developed by wei li for more
! efficiency.
!
! history: jan. 2009 by yuanfu xie
!          apr. 2009 by yuanfu xie for using laps
!          routine for height computation.
!*************************************************

      implicit none
! --------------------
      integer, parameter  :: nh = 6, nr = 16, ln = 255
      character(len=8)    :: ss, sid, typ
      integer             :: l, nl, ir, dt, nc, nw, o, os, is, ip, i, j, k, kk
      integer             :: idx(2, numdims)        ! interpolation indices
      real                :: coe(2, numdims)        ! interpolation coefficents
      real                :: bd, p1
      real(kind=8)        :: ra(nr, ln), hd(nh)
      real                :: x, y, z, p, t, del, uu, vv, tt, zz, ze, op(numdims), di, sp
      integer, external :: ireadsb, ireadmg, i4dy
      real                :: height_to_pressure, rlevel_of_field        ! laps functions

      real                :: az, ea            ! added by zhongjie he

      ! variables for laps obs: by yuanfu
      character*150       :: od
      character*9         :: a9, wfo_fname13_to_fname9
      character*13              :: yyyymmdd_hhmm
      integer             :: ld, i4, n4, st
!jhui
      integer :: data_acar
      integer :: tw_u, tw_t, tw_sh, tw_p
      equivalence(sid, hd(5))
!jhui

      ! include a statement function for converting sh to 'rh' = sh/s2r(p) by yuanfu xie:
      include 'sh2rh.inc'

      data_acar = 0
      tw_u = 0
      tw_t = 0
      tw_sh = 0
      tw_p = 0

      call get_filespec('bufr', 2, od, st)        ! get path to bufr directory
      call get_file_time(od, lapsi4t, n4)        ! get nearest i4time into n4
      call s_len(od, ld)
      call make_fnam_lp(n4, a9, st)
      od(ld + 4:ld + 9) = od(ld - 4:ld)
      od(ld - 5:ld + 3) = a9
      ld = ld + 9

      ! open the bufr file:
      open (unit=lurao, file=od(1:ld), form='unformatted', status='old', action='read')

      call openbf(lurao, 'in', lurao)        ! open bufr in channel

      nc = 0                ! number of valid obs
      nw = 0                ! number of obs in bufr

      o = 0
      do while (ireadmg(lurao, ss, dt) .eq. 0)                ! open a bufr message
         do while (ireadsb(lurao) .eq. 0)                ! open a subset

            call ufbint(lurao, hd, nh, 1, ir, 'xob yob elv dhr sid typ')        ! read header
!jhui
            if (sid == "acar    ") then
               data_acar = data_acar + 1
            end if
            ! read obs:
            call ufbint(lurao, ra, nr, ln, nl, &
                        'xdr ydr prss pob poe hrdr uob vob woe zob zoe tob toe qob qoe pmo')

            do l = 1, nl                ! for all levels
               nw = nw + 1

               ! for x longitude
               op(1) = hd(1)
               if (ra(1, l) .lt. bufrmiss) op(1) = ra(1, l)        ! balloon drift x
               if (op(1) .ge. bufrmiss) cycle                        ! invalid obs longitude
               if (op(1) .lt. 0.0) op(1) = op(1) + 360.0        ! eastward longitude
               ! for y latitude
               op(2) = hd(2)
               if (ra(2, l) .lt. bufrmiss) op(2) = ra(2, l)        ! balloon drift y
               if (op(2) .ge. bufrmiss) cycle                        ! invalid obs latitude

               ! check whether the obs in the horizontal domain:
               if (if_test .ne. 1) then
                  call latlon_to_rlapsgrid(op(2), op(1), y00, x00, fcstgrd(1), fcstgrd(2), x, y, is)
               else
                  call obstogrid(op(2), op(1), y00, x00, fcstgrd(1), fcstgrd(2), x, y, is)
               end if
               if (x .lt. 1.0 .or. y .lt. 1.0 .or. &
                   x .gt. fcstgrd(1) .or. y .gt. fcstgrd(2) .or. &
                   is .ne. 1) cycle                                ! outside horizontal domain

               ! for t location
               ! fit into laps time frames by yuanfu:
               i4 = mod(dt, 100)*100
               write (yyyymmdd_hhmm, 1) 20000000 + dt/100, '_', i4
1              format(i8, a1, i4)
               a9 = wfo_fname13_to_fname9(yyyymmdd_hhmm)
               call cv_asc_i4time(a9, i4, st)
               ! if (i4 .ne. lapsi4t) then
               !   write(*,*) 'analysis time does not match: ',i4,itime2(2)
               !   stop
               ! endif
               t = i4 - itime2(1) + 3600*hd(4)
               if (t .lt. t00(1) .or. t .gt. t00(fcstgrd(4))) cycle        ! out of time window

               ! for p location
               p = bufrmiss
               ! use observed pressure in pascal as bufr saves in mb:
               if (ra(4, l) .lt. bufrmiss) p = ra(4, l)*100.0d0
               ! for zz observation
               zz = bufrmiss                                        ! initial
               if (ra(10, l) .lt. bufrmiss) zz = ra(10, l)         ! use observed height
               ze = ra(11, l)                                        ! obs error in height
               if (p .ge. bufrmiss .and. zz .ge. bufrmiss) cycle ! invalid vertical obs

               ! interpolation: background to obs:
               idx(1, 1) = int(x)
               idx(1, 2) = int(y)
               idx(2, 1:2) = min0(idx(1, 1:2) + 1, fcstgrd(1:2))
               coe(1, 1) = idx(2, 1) - x
               coe(1, 2) = idx(2, 2) - y
               coe(2, 1:2) = 1.0 - coe(1, 1:2)        ! assume horizontal grid distance is 1
               if (fcstgrd(4) .eq. 1) then
                  print *, 'rdbufrobs -- error: temporal grid has one gridpoint!'
                  stop
               end if
               del = (itime2(2) - itime2(1))/(fcstgrd(4) - 1)        ! delta time
               idx(1, 4) = int(t/del) + 1
               idx(2, 4) = min0(idx(1, 4) + 1, fcstgrd(4))
               coe(2, 4) = t/del + 1 - idx(1, 4)
               coe(1, 4) = 1.0 - coe(2, 4)

               ! for missing pressure, convert from background height:
               ip = 0                ! direct pressure obs; 1 becomes hydrostatic obs

               ! fill mslp if available:
               if (zz .le. 10.0 .and. ra(16, l) .ne. bufrmiss) p = ra(16, l)*100.0d0

               if (p .ge. bufrmiss) then        ! use hydrostatic relation
                  ip = 1                ! hydrostatic obs
                  p = height_to_pressure(zz, bk0(1, 1, 1, idx(1, 4), pressure), p00, &
                                         fcstgrd(1), fcstgrd(2), fcstgrd(3), idx(1, 1), idx(1, 2))
               end if

               ! find the vertical level:
               kk = rlevel_of_field(p, p00, 1, 1, fcstgrd(3), 1, 1, is)
               if (is .ne. 1) cycle        ! out of vertical domain

               ! vertical interpolation:
               if (kk .eq. fcstgrd(3)) then
                  idx(1, 3) = kk - 1
                  coe(1, 3) = 0.0
               else
                  idx(1, 3) = kk
                  coe(1, 3) = (p - p00(kk + 1))/(p00(kk) - p00(kk + 1))
               end if
               idx(2, 3) = idx(1, 3) + 1
               coe(2, 3) = 1.0 - coe(1, 3)
               ! vertical grid positin:
               z = idx(1, 3) + coe(2, 3)

               ! for wind observation
               if (ra(9, l) .lt. bufrmiss .and. &
                   ra(7, l) .lt. bufrmiss .and. ra(8, l) .lt. bufrmiss) then

                  ! save obs into data arrays:
                  o = o + 2        ! total single obs count
                  nst(u_cmpnnt) = nst(u_cmpnnt) + 1
                  nst(v_cmpnnt) = nst(v_cmpnnt) + 1        ! wind obs counts
                  obp(1, nst(u_cmpnnt), u_cmpnnt) = x - 1
                  obp(2, nst(u_cmpnnt), u_cmpnnt) = y - 1
                  obp(3, nst(u_cmpnnt), u_cmpnnt) = z - 1
                  obp(4, nst(u_cmpnnt), u_cmpnnt) = t/(t00(2) - t00(1)) ! wind obs locations
                  obp(1, nst(v_cmpnnt), v_cmpnnt) = x - 1
                  obp(2, nst(v_cmpnnt), v_cmpnnt) = y - 1
                  obp(3, nst(v_cmpnnt), v_cmpnnt) = z - 1
                  obp(4, nst(v_cmpnnt), v_cmpnnt) = t/(t00(2) - t00(1))
!jhui
                  if (obp(1, nst(u_cmpnnt), u_cmpnnt) .gt. 38 .and. &
                      obp(1, nst(u_cmpnnt), u_cmpnnt) .lt. 112 .and. &
                      obp(2, nst(u_cmpnnt), u_cmpnnt) .gt. 33 .and. &
                      obp(2, nst(u_cmpnnt), u_cmpnnt) .lt. 107) then
                     tw_u = tw_u + 1
                  end if

                  ! save increment:
                  uu = 0.0
                  vv = 0.0
                  do kk = 1, 2
                     do k = 1, 2
                        do j = 1, 2
                           do i = 1, 2
                              uu = uu + coe(i, 1)*coe(j, 2)*coe(k, 3)*coe(kk, 4)* &
                                   bk0(idx(i, 1), idx(j, 2), idx(k, 3), idx(kk, 4), u_cmpnnt)
                              vv = vv + coe(i, 1)*coe(j, 2)*coe(k, 3)*coe(kk, 4)* &
                                   bk0(idx(i, 1), idx(j, 2), idx(k, 3), idx(kk, 4), v_cmpnnt)
                           end do
                        end do
                     end do
                  end do
                  obs(nst(u_cmpnnt), u_cmpnnt) = ra(7, l) - uu
                  obs(nst(v_cmpnnt), v_cmpnnt) = ra(8, l) - vv
                  obe(nst(u_cmpnnt), u_cmpnnt) = 0.5        !ra(9,l)        ! wind error
                  obe(nst(v_cmpnnt), v_cmpnnt) = 0.5        !ra(9,l)
!jhui
!          obe(nst(u_cmpnnt),u_cmpnnt) =0.3        !ra(9,l)        ! wind error
!          obe(nst(v_cmpnnt),v_cmpnnt) =0.3        !ra(9,l)
                  ! save obs into pig file of laps as wind direction and speed:
                  ! uu = ra(7,l)
                  ! vv = ra(8,l)
                  ! call uv_to_disp(uu,vv,di,sp)
                  ! write(pigobs_channel,*) x-1,y-1,z-1,di,sp,sid

                  ! add a threshold check: yuanfu june 2010
                  if ((abs(obs(nst(u_cmpnnt), u_cmpnnt)) .gt. 20.0) .or. &
                      (abs(obs(nst(v_cmpnnt), v_cmpnnt)) .gt. 20.0)) then
                     ! over the threshold, remove this data:
                     o = o - 2
                     nst(u_cmpnnt) = nst(u_cmpnnt) - 1
                     nst(v_cmpnnt) = nst(v_cmpnnt) - 1
                  end if
               end if

               ! for temperature observation
               if (ra(12, l) .lt. bufrmiss .and. ra(13, l) .lt. bufrmiss) then
                  ! save obs into data arrays:
                  o = o + 1        ! total single obs count
                  nst(temprtur) = nst(temprtur) + 1        ! temperature obs counts
                  obp(1, nst(temprtur), temprtur) = x - 1
                  obp(2, nst(temprtur), temprtur) = y - 1
                  obp(3, nst(temprtur), temprtur) = z - 1
                  obp(4, nst(temprtur), temprtur) = t/(t00(2) - t00(1)) ! temperatuer obs location
                  obs(nst(temprtur), temprtur) = ra(12, l) + 273.15d0        ! in kelvin
                  obe(nst(temprtur), temprtur) = 0.5        !ra(13,l)                        ! obs error
!jhui
!          obe(nst(temprtur),temprtur) =0.3        !ra(13,l)                        ! obs error
!jhui
                  if (obp(1, nst(temprtur), temprtur) .gt. 38 .and. &
                      obp(1, nst(temprtur), temprtur) .lt. 112 .and. &
                      obp(2, nst(temprtur), temprtur) .gt. 33 .and. &
                      obp(2, nst(temprtur), temprtur) .lt. 107) then
                     tw_t = tw_t + 1
                  end if

                  ! save obs innovation:
                  tt = 0.0
                  do kk = 1, 2
                     do k = 1, 2
                        do j = 1, 2
                           do i = 1, 2
                              tt = tt + coe(i, 1)*coe(j, 2)*coe(k, 3)*coe(kk, 4)* &
                                   bk0(idx(i, 1), idx(j, 2), idx(k, 3), idx(kk, 4), temprtur)
                           end do
                        end do
                     end do
                  end do

                  obs(nst(temprtur), temprtur) = obs(nst(temprtur), temprtur) - tt

                  ! add a threshold check: yuanfu june 2010
                  if (abs(obs(nst(temprtur), temprtur)) .gt. 10.0) then
                     ! over the threshold, remove this data:
                     o = o - 1
                     nst(temprtur) = nst(temprtur) - 1
                  end if

               end if

               ! for specific humidity observation
               if (ra(14, l) .lt. bufrmiss .and. ra(15, l) .lt. bufrmiss) then
                  ! save obs into data arrays:
                  o = o + 1        ! total single obs count
                  nst(humidity) = nst(humidity) + 1        ! humidity obs counts
                  obp(1, nst(humidity), humidity) = x - 1
                  obp(2, nst(humidity), humidity) = y - 1
                  obp(3, nst(humidity), humidity) = z - 1
                  obp(4, nst(humidity), humidity) = t/(t00(2) - t00(1)) ! humidity obs location

!jhui
                  if (obp(1, nst(humidity), humidity) .gt. 38 .and. &
                      obp(1, nst(humidity), humidity) .lt. 112 .and. &
                      obp(2, nst(humidity), humidity) .gt. 33 .and. &
                      obp(2, nst(humidity), humidity) .lt. 107) then
                     tw_sh = tw_sh + 1
                  end if
                  ! save obs innovation:
                  tt = 0.0
                  do kk = 1, 2
                     do k = 1, 2
                        do j = 1, 2
                           do i = 1, 2
                              tt = tt + coe(i, 1)*coe(j, 2)*coe(k, 3)*coe(kk, 4)* &
                                   bk0(idx(i, 1), idx(j, 2), idx(k, 3), idx(kk, 4), humidity)
                           end do
                        end do
                     end do
                  end do
                  ! convert sh obs to 'rh' = sh/s2r(p) by yuanfu xie: dec. 2013
                  obs(nst(humidity), humidity) = ra(14, l)/1000.0/s2r(p/100.0) - tt        ! humidity obs: bufr uses mg/kg
                  obe(nst(humidity), humidity) = 1.0        !ra(15,l)        ! obs error

                  ! add a threshold check: yuanfu june 2010
                  if (abs(obs(nst(humidity), humidity)) .gt. 10.0) then
                     ! over the threshold, remove this data:
                     o = o - 1
                     nst(humidity) = nst(humidity) - 1
                  end if

               end if
               ! for height observation: if pressure is not derived from height:
               if (ip .eq. 0 .and. zz .lt. bufrmiss .and. ze .lt. bufrmiss) then

                  ! remove all high level height obs:
                  !if (l .gt. 3) cycle

                  ! save obs into data arrays:
                  o = o + 1        ! total single obs count
                  nst(pressure) = nst(pressure) + 1        ! height obs counts
                  obp(1, nst(pressure), pressure) = x - 1
                  obp(2, nst(pressure), pressure) = y - 1
                  obp(3, nst(pressure), pressure) = z - 1
                  obp(4, nst(pressure), pressure) = t/(t00(2) - t00(1)) ! height obs location

!jhui
                  if (obp(1, nst(pressure), pressure) .gt. 38 .and. &
                      obp(1, nst(pressure), pressure) .lt. 112 .and. &
                      obp(2, nst(pressure), pressure) .gt. 33 .and. &
                      obp(2, nst(pressure), pressure) .lt. 107) then
                     tw_p = tw_p + 1
                  end if

                  ! save obs innovation:
                  tt = 0.0
                  do kk = 1, 2
                     do k = 1, 2
                        do j = 1, 2
                           do i = 1, 2
                              tt = tt + coe(i, 1)*coe(j, 2)*coe(k, 3)*coe(kk, 4)* &
                                   bk0(idx(i, 1), idx(j, 2), idx(k, 3), idx(kk, 4), pressure)
                           end do
                        end do
                     end do
                  end do
                  obs(nst(pressure), pressure) = zz - tt                ! height obs
                  obe(nst(pressure), pressure) = 2.0        !ze        ! obs error

                  ! add a threshold check: yuanfu june 2010
                  if (abs(zz - tt) .gt. 50) then ! use 50m as the threshold value now.
                     ! over the threshold, remove this data:
                     o = o - 1
                     nst(pressure) = nst(pressure) - 1
                  end if

               end if

               nc = nc + 1                        ! number of valid obs counts

            end do
         end do
      end do
      nallobs = o
      print *, 'the number of location is', nw, 'and', nc, 'available'
      print *, 'rdbufrobs: number of all obs = ', nallobs
      call closbf(lurao)
      close (lurao)
!  close(lundx)
!  close(lunew)

      return
   end subroutine rdbufrobs_xie

   subroutine handleobs(op, ob, oe, os, o, ip, az, ea)
!*************************************************
! calculate the difference between observation and background
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: i, j, k, t, m, n, np(maxdims), nn(maxdims), is, o, os, ip
      real  :: x, y, p
      real  :: di, sp        ! yuanfu for output pig file
      real  :: ac(numdims, ngptobs), oc(numdims), co(ngptobs), ht
      real  :: op(numdims), ob, oe

      integer  :: uu, vv
      real     :: az, tu, tv, ea               ! added by zhongjie he
      real     :: d2r                       ! converter from deg to radian

      d2r = 3.14159/180.0

! --------------------
      uu = u_cmpnnt
      vv = v_cmpnnt
      call latlon_to_rlapsgrid(op(2), op(1), y00, x00, fcstgrd(1), fcstgrd(2), x, y, is)
!  call obstogrid(op(2),op(1),y00,x00,fcstgrd(1),fcstgrd(2),x,y,is)

      if (x .lt. 1.0 .or. y .lt. 1.0 .or. x .gt. fcstgrd(1) .or. y .gt. fcstgrd(2) .or. is .ne. 1) return
      call vrtclpstn(fcstgrd(3), limit_3, p00, op(3), p, is)

      if (is .ne. 1) return
      do n = 1, maxdims
         np(n) = 1
      end do
      np(1) = int(x)
      np(2) = int(y)
      np(3) = int(p)
      do n = 1, maxdims
         if (np(n) .eq. fcstgrd(n) .and. fcstgrd(n) .ne. 1) np(n) = fcstgrd(n) - 1
      end do
      oc(1) = x
      oc(2) = y
      oc(3) = p
      m = 0
!==============================================================
!  do t=np(4),min0(np(4)+1,fcstgrd(4))
!  do k=np(3),min0(np(3)+1,fcstgrd(3))
!  do j=np(2),min0(np(2)+1,fcstgrd(2))
!  do i=np(1),min0(np(1)+1,fcstgrd(1))
!    nn(1)=i
!    nn(2)=j
!    nn(3)=k
!    nn(4)=t
!    m=m+1
!    do n=1,numdims
!      ac(n,m)=nn(n)*1.0d0
!    enddo
!  enddo
!  enddo
!  enddo
!  enddo
!======================= modified by zhongjie he ==============
      do t = np(4), np(4) + 1
      do k = np(3), np(3) + 1
      do j = np(2), np(2) + 1
      do i = np(1), np(1) + 1
         nn(1) = min0(i, fcstgrd(1))
         nn(2) = min0(j, fcstgrd(2))
         nn(3) = min0(k, fcstgrd(3))
         nn(4) = min0(t, fcstgrd(4))
         m = m + 1
         do n = 1, numdims
            ac(n, m) = nn(n)*1.0d0
         end do
      end do
      end do
      end do
      end do
!==============================================================
      call interpltn(numdims, ngptobs, co, ac, oc)
      !call interpltn_xie(numdims,ngptobs,co,ac,oc,3,numgrid(3),ppm)
      ht = 0.0
      m = 0
      tu = 0.0
      tv = 0.0
      if (os .le. numstat) then
         do t = np(4), min0(np(4) + 1, fcstgrd(4))
         do k = np(3), min0(np(3) + 1, fcstgrd(3))
         do j = np(2), min0(np(2) + 1, fcstgrd(2))
         do i = np(1), min0(np(1) + 1, fcstgrd(1))
            m = m + 1
            ht = ht + co(m)*bk0(i, j, k, t, os)
         end do
         end do
         end do
         end do
      elseif (os .eq. numstat + 1) then        ! for radar data , by zhongjie he
         do t = np(4), min0(np(4) + 1, fcstgrd(4))
         do k = np(3), min0(np(3) + 1, fcstgrd(3))
         do j = np(2), min0(np(2) + 1, fcstgrd(2))
         do i = np(1), min0(np(1) + 1, fcstgrd(1))
            m = m + 1
            tu = tu + co(m)*bk0(i, j, k, t, uu)
            tv = tv + co(m)*bk0(i, j, k, t, vv)
         end do
         end do
         end do
         end do
         ht = (tu*sin(d2r*az) + tv*cos(d2r*az))*cos(d2r*ea)
      elseif (os .eq. numstat + 2) then        ! for sfmr data , by zhongjie he
         do t = np(4), min0(np(4) + 1, fcstgrd(4))
         do k = np(3), min0(np(3) + 1, fcstgrd(3))
         do j = np(2), min0(np(2) + 1, fcstgrd(2))
         do i = np(1), min0(np(1) + 1, fcstgrd(1))
            m = m + 1
            tu = tu + co(m)*bk0(i, j, k, t, uu)
            tv = tv + co(m)*bk0(i, j, k, t, vv)
         end do
         end do
         end do
         end do
         ht = sqrt(tu*tu + tv*tv)
      end if

      if (os .eq. 4 .and. abs(ob - ht) .ge. 8.0) return
      if (os .eq. 3 .and. abs(ob - ht) .ge. 50.0) return
      if (ip .eq. 1 .and. os .eq. 3) return
      ob = ob - ht

      call vrtclpstn8(maxgrid(3), limit_3, pp0, op(3), p, is)

      if (is .ne. 1) return

      oc(3) = p
      o = o + 1
      nst(os) = nst(os) + 1
      if (o .gt. nobsmax) stop 'number of observations exceeded'
      do n = 1, numdims
         obp(n, nst(os), os) = oc(n) - 1.0d0
      end do
      obs(nst(os), os) = ob
      obe(nst(os), os) = oe
      if (os .eq. numstat + 1) then
         oba(nst(os), 1) = az
         oba(nst(os), 2) = ea
      end if

      return
   end subroutine handleobs

   subroutine ht_to_prs(x0, y0, z, p, is)
!*************************************************
! convert height to pressure to use as many data as possible
! history: september 2007, coded by wei li.
!*************************************************
      implicit none
! --------------------
      integer  :: is, i, j, k, t, n, m, np(2), zz
      integer  :: lm(2)
      real  :: x, y, x0, y0
      real  :: z, p, ac(2, 4), oc(2), co(4), ht(fcstgrd(3))
! --------------------
      zz = pressure                                                        !  by zhongjie he
      call latlon_to_rlapsgrid(y0, x0, y00, x00, fcstgrd(1), fcstgrd(2), x, y, is)
!   call obstogrid(y0,x0,y00,x00,fcstgrd(1),fcstgrd(2),x,y,is)        !  by zhongjie he

      if (x .lt. 1.0 .or. y .lt. 1.0 .or. x .gt. fcstgrd(1) .or. y .gt. fcstgrd(2) .or. is .ne. 1) then
         is = 0
         return
      end if
      oc(1) = x
      oc(2) = y
      np(1) = int(x)
      np(2) = int(y)
      do n = 1, 2
         if (np(n) .eq. fcstgrd(n) .and. fcstgrd(n) .ne. 1) np(n) = fcstgrd(n) - 1
      end do
      m = 0
!==============================================
!  do j=np(2),min0(np(2)+1,fcstgrd(2))
!  do i=np(1),min0(np(1)+1,fcstgrd(1))
!    m=m+1
!    ac(1,m)=i*1.0d0
!    ac(2,m)=j*1.0d0
!  enddo
!  enddo
!============== modified by zhongjie he =======
      do n = 1, 2
         lm(n) = np(n) + 1
         if (numdims .lt. n) lm(n) = np(n)
      end do
      do j = np(2), lm(2)
      do i = np(1), lm(1)
         m = m + 1
         ac(1, m) = min0(i, fcstgrd(1))*1.0
         ac(2, m) = min0(j, fcstgrd(2))*1.0
      end do
      end do
!==============================================
      call interpltn(2, 4, co, ac, oc)
      !call interpltn_xie(2,4,co,ac,oc,3,numgrid(3),ppm)
      t = 1
      do k = 1, fcstgrd(3)
         ht(k) = 0.0
         m = 0
         do j = np(2), min0(np(2) + 1, fcstgrd(2))
         do i = np(1), min0(np(1) + 1, fcstgrd(1))
            m = m + 1
            ht(k) = ht(k) + co(m)*bk0(i, j, k, t, zz)
         end do
         end do
      end do
      call zpconvert(fcstgrd(3), ht, p00, z, p, is)
      return
   end subroutine ht_to_prs

   subroutine rdobstest
!*************************************************
! read in test data
! history: february 2008, coded by zhongjie he.
!*************************************************
      implicit none
! --------------------
      character*8     :: sid
      integer         :: l, o, os, is, ip
      real            :: x, y, p, uv, zz, ue, ea, az, ze, t
      real            :: op(numdims), ob, oe
      integer         :: nc, nw, lmax
      real            :: mf                 ! number flag of missed data

      mf = -9999.0
      nc = 0
      nw = 0
      sid = 'test obs'
      o = nallobs
      open (11, file='obs_radar_cnvtn.dat', action='read')
!  open(11,file='obs_radar.dat',action='read')
      x_radar = 350.
      y_radar = 150
      read (11, *) lmax
      do l = 1, lmax
         nw = nw + 1
         read (11, '(8f10.3,i3,f12.3)') x, y, zz, t, az, ea, uv, ue, os, p
!      read(11,*) x,y,zz,t,az,ea,uv,ue,os,p
!      read(11,'(7f10.3,i3,f12.3)') x,y,zz,az,ea,uv,ue,os,p      !  for presure coordinate and p is contented in the file of obs
         if (abs(x - mf) .le. 1e-5) cycle
         if (abs(y - mf) .le. 1e-5) cycle
         if (abs(az - mf) .le. 1e-5) cycle
         if (abs(ea - mf) .le. 1e-5) cycle
         if (abs(uv - mf) .le. 1e-5) cycle

!================= for test ===============
!      if(mod(l,11).ne.1) cycle
!      if(os.eq.numstat+1 ) cycle
!================= for test ===============

! transform
         op(1) = x
         op(2) = y
         ip = 0
         if (os .eq. numstat + 1) ip = 1

         if (ifpcdnt .eq. 1) then
            if (abs(p - mf) .le. 1e-5) then
               call ht_to_prs(op(1), op(2), zz, p, is)
               ip = 1
               if (is .ne. 1) cycle
            end if
            op(3) = p
            if (numdims .ge. 4) op(4) = t
            ze = 2.0
            ob = uv
            oe = ue
            call handleobs_sigma(op, ob, oe, os, o, ip, az, ea, sid)
            nc = nc + 1

            ob = zz
            oe = ze
            os = pressure
            call handleobs_sigma(op, ob, oe, os, o, ip, az, ea, sid)
            nc = nc + 1
         else
            op(3) = zz
            if (numdims .ge. 4) op(4) = t
            ze = 2.0
            ob = uv
            oe = ue
            call handleobs_sigma(op, ob, oe, os, o, ip, az, ea, sid)
            nc = nc + 1
         end if

      end do
      nallobs = o
      close (11)

      do os = 1, numstat + 2
         if (os .eq. u_cmpnnt) print *, 'the number of observed u data is:', nst(os)
         if (os .eq. v_cmpnnt) print *, 'the number of observed v data is:', nst(os)
         if (os .eq. w_cmpnnt) print *, 'the number of observed w data is:', nst(os)
         if (os .eq. pressure) print *, 'the number of observed pressure data is:', nst(os)
         if (os .eq. temprtur) print *, 'the number of observed temprtur data is:', nst(os)
         if (os .eq. numstat + 1 .and. nst(numstat + 1) .ge. 1) print *, 'the number of radar data is:', nst(os)
         if (os .eq. numstat + 2 .and. nst(numstat + 2) .ge. 1) print *, 'the number of sfmr data is:', nst(os)
      end do

      return
   end subroutine rdobstest

   subroutine handleobs_sigma(op, ob, oe, os, o, ip, az, ea, sid)
!*************************************************
! calculate the difference between observation and background
! history: march 2007, coded by zhongjie he.
!*************************************************
      implicit none
! --------------------
      character*8, intent(in) :: sid        ! station name
      integer  :: i, j, k, t, m, n, np(maxdims), nn(maxdims), is, o, os, ip
      integer  :: lm(maxdims)
      real  :: x, y, p, tm
      real  :: di, sp                ! yuanfu for output pig file
      real  :: ac(numdims, ngptobs), oc(numdims), co(ngptobs), ht
      real  :: op(numdims), ob, oe

      integer  :: uu, vv, ww
      real     :: az, tu, tv, ea, hs, he, tw
      real     :: d2r
! --------------------

      d2r = 3.14159/180.0

      uu = u_cmpnnt
      vv = v_cmpnnt
      ww = w_cmpnnt

      if (ip .eq. 1 .and. os .eq. 3) return

      if (if_test .ne. 1) then
         call latlon_to_rlapsgrid(op(2), op(1), y00, x00, fcstgrd(1), fcstgrd(2), x, y, is)
      else
         call obstogrid(op(2), op(1), y00, x00, fcstgrd(1), fcstgrd(2), x, y, is)
      end if

      if (x .lt. 1.0 .or. y .lt. 1.0 .or. x .gt. fcstgrd(1) .or. y .gt. fcstgrd(2) .or. is .ne. 1) return

      if (ifpcdnt .eq. 1) then                      ! for pressure coordinate
         call vrtclpstn(fcstgrd(3), limit_3, p00, op(3), p, is)
      elseif (ifpcdnt .eq. 2) then                  ! for height coordinate
         call vrtclpstn(fcstgrd(3), limit_3, z00, op(3), p, is)
      elseif (ifpcdnt .eq. 0) then                  ! for sigma coordinate
         do n = 1, 2
            np(n) = 1
         end do
         np(1) = int(x)
         np(2) = int(y)
         np(3) = 1
         np(4) = 1
         do n = 1, maxdims
            if (np(n) .eq. fcstgrd(n) .and. fcstgrd(n) .ne. 1) np(n) = fcstgrd(n) - 1
         end do
         oc(1) = x
         oc(2) = y
         oc(3) = 1
         oc(4) = 1
         m = 0
!============================================================
!    do t=np(4),min0(np(4)+1,fcstgrd(4))
!    do k=np(3),min0(np(3)+1,fcstgrd(3))
!    do j=np(2),min0(np(2)+1,fcstgrd(2))
!    do i=np(1),min0(np(1)+1,fcstgrd(1))
!      nn(1)=i
!      nn(2)=j
!      nn(3)=k
!      nn(4)=t
!      m=m+1
!      do n=1,numdims
!        ac(n,m)=nn(n)*1.0d0
!      enddo
!    enddo
!    enddo
!    enddo
!    enddo
!=====================modified by zhongjie he ==============
         do n = 1, maxdims
            lm(n) = np(n) + 1
            if (numdims .lt. n) lm(n) = np(n)
         end do
         do t = np(4), lm(4)
         do k = np(3), lm(3)
         do j = np(2), lm(2)
         do i = np(1), lm(1)
            nn(1) = min0(i, fcstgrd(1))
            nn(2) = min0(j, fcstgrd(2))
            nn(3) = min0(k, fcstgrd(3))
            nn(4) = min0(t, fcstgrd(4))
            m = m + 1
            do n = 1, numdims
               ac(n, m) = nn(n)*1.0d0
            end do
         end do
         end do
         end do
         end do
!============================================================
         ! call interpltn_xie(numdims,ngptobs,co,ac,oc,3,numgrid(3),ppm)
         call interpltn(numdims, ngptobs, co, ac, oc)
         m = 0
         hs = 0.0
         he = 0.0
         do t = np(4), min0(np(4) + 1, fcstgrd(4))
         do k = np(3), min0(np(3) + 1, fcstgrd(3))
         do j = np(2), min0(np(2) + 1, fcstgrd(2))
         do i = np(1), min0(np(1) + 1, fcstgrd(1))
            m = m + 1
            hs = hs + co(m)*heightl(i, j)
            he = he + co(m)*heightu(i, j)
         end do
         end do
         end do
         end do
         op(3) = (op(3) - hs)/(he - hs)
         call vrtclpstn(fcstgrd(3), limit_3, z00, op(3), p, is)
      end if
!jhui
!  for radar ingest data that outside 0-3600
      if (os .gt. numstat) then
         if (op(4) .gt. 3600) op(4) = 3600
         if (op(4) .lt. 0) op(4) = 0
         if (numdims .ge. 4) call vrtclpstn(fcstgrd(4), limit_4, t00, op(4), tm, is)
!  else
!  if(numdims.ge.4) call vrtclpstn(fcstgrd(4),limit_4,t00,op(4),tm,is)
      end if

      if (is .ne. 1) return
      do n = 1, maxdims
         np(n) = 1
      end do
      np(1) = int(x)
      np(2) = int(y)
      np(3) = int(p)
      if (numdims .ge. 4) np(4) = int(tm)
      do n = 1, maxdims
         if (np(n) .eq. fcstgrd(n) .and. fcstgrd(n) .ne. 1) np(n) = fcstgrd(n) - 1
      end do
      oc(1) = x
      oc(2) = y
      oc(3) = p
      if (numdims .ge. 4) oc(4) = tm
      m = 0
!===================================================
!  do t=np(4),min0(np(4)+1,fcstgrd(4))
!  do k=np(3),min0(np(3)+1,fcstgrd(3))
!  do j=np(2),min0(np(2)+1,fcstgrd(2))
!  do i=np(1),min0(np(1)+1,fcstgrd(1))
!    nn(1)=i
!    nn(2)=j
!    nn(3)=k
!    nn(4)=t
!    m=m+1
!    do n=1,numdims
!      ac(n,m)=nn(n)*1.0d0
!    enddo
!  enddo
!  enddo
!  enddo
!  enddo
!=================== modified by zhongjie he========
      do n = 1, maxdims
         lm(n) = np(n) + 1
         if (numdims .lt. n) lm(n) = np(n)
      end do
      do t = np(4), lm(4)
      do k = np(3), lm(3)
      do j = np(2), lm(2)
      do i = np(1), lm(1)
         nn(1) = min0(i, fcstgrd(1))
         nn(2) = min0(j, fcstgrd(2))
         nn(3) = min0(k, fcstgrd(3))
         nn(4) = min0(t, fcstgrd(4))
         m = m + 1
         do n = 1, numdims
            ac(n, m) = nn(n)*1.0d0
         end do
      end do
      end do
      end do
      end do
!===================================================
      ! call interpltn_xie(numdims,ngptobs,co,ac,oc,3,numgrid(3),ppm)
      call interpltn(numdims, ngptobs, co, ac, oc)
      ht = 0.0
      m = 0
      tu = 0.0
      tv = 0.0
      tw = 0.0
      if (os .le. numstat) then
         do t = np(4), min0(np(4) + 1, fcstgrd(4))
         do k = np(3), min0(np(3) + 1, fcstgrd(3))
         do j = np(2), min0(np(2) + 1, fcstgrd(2))
         do i = np(1), min0(np(1) + 1, fcstgrd(1))
            m = m + 1
            ht = ht + co(m)*bk0(i, j, k, t, os)
         end do
         end do
         end do
         end do
      elseif (os .eq. numstat + 1) then        ! for radar data , by zhongjie he
         do t = np(4), min0(np(4) + 1, fcstgrd(4))
         do k = np(3), min0(np(3) + 1, fcstgrd(3))
         do j = np(2), min0(np(2) + 1, fcstgrd(2))
         do i = np(1), min0(np(1) + 1, fcstgrd(1))
            m = m + 1
            tu = tu + co(m)*bk0(i, j, k, t, uu)
            tv = tv + co(m)*bk0(i, j, k, t, vv)
            if (ww .ne. 0) tw = tw + co(m)*bk0(i, j, k, t, ww)
         end do
         end do
         end do
         end do
         ht = (tu*sin(d2r*az) + tv*cos(d2r*az))*cos(d2r*ea) + tw*sin(d2r*ea)
      elseif (os .eq. numstat + 2) then        ! for sfmr data , by zhongjie he
         ! for sfmr data, no substract the background as its operator is nonlinear
         ht = 0.0
      elseif (os .eq. numstat + 3) then
         ht = 0.0
      end if

      if (os .eq. temprtur .and. abs(ob - ht) .ge. 8.0) return
      if (os .eq. pressure .and. abs(ob - ht) .ge. 50.0) return
      ob = ob - ht

      if (ifpcdnt .eq. 0 .or. ifpcdnt .eq. 2) then              ! for sigma and height coordinate
         call vrtclpstn8(maxgrid(3), limit_3, zzb, op(3), p, is)
      else                               ! for presure coordinate
         call vrtclpstn8(maxgrid(3), limit_3, pp0, op(3), p, is)
      end if

      if (is .ne. 1) return

      oc(3) = p
      o = o + 1
      nst(os) = nst(os) + 1
      if (o .gt. nobsmax) stop 'number of observations exceeded'
      do n = 1, numdims
         obp(n, nst(os), os) = oc(n) - 1.0d0
      end do
      obs(nst(os), os) = ob
      obe(nst(os), os) = oe
      if (os .eq. numstat + 1) then
         oba(nst(os), 1) = az
         oba(nst(os), 2) = ea
      end if

      return
   end subroutine handleobs_sigma

end module readobserves
