subroutine read_laps_radar
!==================================================================
!  this routine reads in radar through laps ingest codes.
!
!  note:
!  1. the original rdlapsrdr uses laps get_multiradar_vel to ingest
!     which requires the number of gridpoints times the number of
!     radar and causes problems for large domain data. the routine
!     avoids those arrays for efficiency in computer memory by
!     repeating what getradar.f under lib does.
!  2. this routine assumes to be used as a member of module of the
!     readobserves in readobserves.f90 and it uses all variables
!     defined in stmas.
!
!  history: august 2012
!  author : yuanfu xie at esrl/gsd/fab
!==================================================================

   use prmtrs_stmas

   implicit none

   ! get radar radial wind obs:
   ! print*,'reading radial wind... due to radar format change: no radial wind ingested!'
   if (fcstgrd(1) .lt. 10000) call radialwind

   ! get radar reflectivity obs:
   print *, 'reading reflectivity...'
   call reflectivity

end subroutine read_laps_radar

subroutine radialwind
!==================================================================
!  this routine reads in radial wind data through laps ingest.
!
!  history: august 2012
!  author : yuanfu xie at esrl/gsd/fab
!==================================================================

   use prmtrs_stmas
   use read_backgrd
   use readobserves
   use mem_namelist

   implicit none

   include 'main_sub.inc'

   ! local variables:
   character :: ext*31, rname*4, sid*8, dir*300
   integer   :: i, j, k, ll, iradar, len_dir, istatus, istat, i4radar, num
   logical   :: apply_map ! laps does not use this for radial wind
   real      :: volnyqt, rlat, rlon, rhgt, op(4), az, sr, ea

   ! for radial wind data on maxgrid:
   integer   :: istart, iend, jstart, jend
   real      :: ri, rj, grid_spacing_actual_mx, grid_spacing_actual_my
   real      :: xradius, yradius

   ! functions:
   integer :: init_timer, ishow_timer
   real :: height_of_level

   ! allocatable arrays:
   real, allocatable :: velo(:, :, :), nyqt(:, :, :), zero(:, :, :)

   ! allocate memory of local arrays:
   allocate (velo(fcstgrd(1), fcstgrd(2), fcstgrd(3)), &
             nyqt(fcstgrd(1), fcstgrd(2), fcstgrd(3)), &
             zero(fcstgrd(1), fcstgrd(2), fcstgrd(3)), &
             stat=istatus)

   ! ingest:
   print *, 'starting ingest radial wind: ', init_timer()
   do iradar = 1, max_radars ! look through all radar

      write (6, *) ' looking for potential radar # ', iradar

      ! get file:
      if (iradar .lt. 10) write (ext, 1) iradar
      if (iradar .ge. 10 .and. iradar .lt. 100) write (ext, 2) iradar
      if (iradar .ge. 100 .and. iradar .lt. 1000) write (ext, 3) iradar
      if (iradar .ge. 1000) then
         print *, 'too many radars: modify this code and rerun!'
      end if
1     format('v0', i1, ' ')
2     format('v', i2, ' ')
3     format('v', i3, ' ')

      ! get that directory:
      call get_directory(ext, dir, len_dir)
      dir = dir(1:len_dir)//'*.'//ext

      ! frequency of reading radial wind data:
      do ll = 2, 2 ! read data at current only !and previous cycle

         call get_file_time(dir, lapsi4t + (ll - 2)*int(0.5*icycle), i4radar)

         ! check file availability:
         if (abs(i4radar - (lapsi4t + (ll - 2)*0.5*icycle)) .ge. 0.5*icycle) then
            write (6, 10) lapsi4t + (ll - 2)*icycle, dir(1:len_dir)
10          format('no radial wind avail at ', i10, ' under ', &
                   /, ' in directory: ', a50)
         else ! find good file:
            num = 0 ! number of grid points where radial wind has value
            call read_radar_vel(i4radar, apply_map, &
                                fcstgrd(1), fcstgrd(2), fcstgrd(3), ext, &
                                velo, nyqt, volnyqt, rlat, rlon, rhgt, rname, num, istatus, istat)
            if (num .gt. 0) then
               ! radar qc for unfolding radial wind:
               zero = 0.0
               call qc_radar_obs(fcstgrd(1), fcstgrd(2), fcstgrd(3), &
                                 rmissing, fcstgrd(1), fcstgrd(2), 0, 0, &
                                 velo, nyqt, num, latitude, longitud, &
                                 rlat, rlon, rhgt, zero, zero, &
                                 bk0(1, 1, 1, ll, 1), bk0(1, 1, 1, ll, 2), volnyqt, &
                                 l_correct_unfolding, l_grid_north, istatus)

               ! estimate this radar coverage:
               call latlon_to_rlapsgrid(rlat, rlon, latitude, longitud, &
                                        fcstgrd(1), fcstgrd(2), ri, rj, istatus)

               ! grid spacing:
               call get_grid_spacing_actual_xy(rlat, rlon &
                                               , grid_spacing_actual_mx &
                                               , grid_spacing_actual_my &
                                               , istatus)

               ! radius on fcstgrid:
               xradius = 300000.0/grid_spacing_actual_mx  ! radius of radial wind in x
               yradius = 300000.0/grid_spacing_actual_my  ! radius of radial wind in y

               ! start and end points on fcstgrd:
               istart = max(nint(ri - xradius), 1)
               jstart = max(nint(rj - yradius), 1)
               iend = min(nint(ri + xradius), fcstgrd(1))
               jend = min(nint(rj + yradius), fcstgrd(2))

               ! pass radial wind to stmas data structure:
               do k = 1, fcstgrd(3)
               do j = jstart, jend
               do i = istart, iend

                  if (velo(i, j, k) .ne. rmissing) then

                     ! calculate azimuth/elevation/slant range:
                     call latlon_to_radar(latitude(i, j), longitud(i, j), &
                                          bk0(i, j, k, ll, 3), az, sr, ea, rlat, rlon, rhgt)

                     op(1) = longitud(i, j)
                     op(2) = latitude(i, j)
                     op(3) = z_fcstgd(k)
                     ! treat data as observed at analysis time frames:
                     op(4) = 0.5*icycle ! i4radar-itime2(1)

                     call handleobs_sigma(op, velo(i, j, k), 0.5, & ! radial obs error: 0.5
                                          numstat + 1, nallobs, 1, az, ea, sid)
                  end if
               end do
               end do
               end do
            end if
         end if

      end do

   end do
   print *, 'ending ingest radial wind: ', ishow_timer()

   ! release local allocatables:
   deallocate (velo, nyqt, zero, stat=istatus)

end subroutine radialwind

subroutine reflectivity
!==================================================================
!  this routine reads in reflectivity data through laps ingest for
!  deriving humidity (sh) lower bounds.
!  note: cloud fraction is also read in for detecting non-cloudness.
!
!  history: august 2012
!  author : yuanfu xie at esrl/gsd/fab
!==================================================================

   use prmtrs_stmas
   use read_backgrd
   use readobserves
   use mem_namelist

   implicit none

   ! local variables:
   character :: ext*31, rname*4, unit*10, comment*150
   integer :: i, j, k, ll, iqc, nref, n2d, n3d, istatus, ix0, ix1, iy0, iy1
   integer :: istatus2d(fcstgrd(1), fcstgrd(2)), &
              istatus3d(fcstgrd(1), fcstgrd(2))
   logical :: reflectivity_bound
   real    :: rhc, tref ! laps uses these scalars
   real    :: rmax(2), rlow(2)

   ! functions:
   real :: make_ssh

   ! include laps cloud parameters:
   include 'laps_cloud.inc'

   ! allocatables: reflectivity, cloud fraction on p-coor and on height-coord:
   real, allocatable :: refl(:, :, :, :), cldf(:, :, :, :), cldh(:, :, :, :)
   real :: cld_hgt(kcloud), cld_prs(kcloud)

   ! include a statement function for converting sh to 'rh' = sh/s2r(p) by yuanfu xie:
   include 'sh2rh.inc'

   ! allocate memory:
   allocate (refl(fcstgrd(1), fcstgrd(2), fcstgrd(3), 3), &
             cldf(fcstgrd(1), fcstgrd(2), fcstgrd(3), 3), &
             cldh(fcstgrd(1), fcstgrd(2), kcloud, 3), &
             stat=istatus)

   ! reference temperaure for converison between sh and rh:
   tref = 0.0 !-132.0

   ! get reflectivity at pressure levels:
   ext = 'lps '
   refl = 0.0
   do ll = 1, 2
      call get_laps_3dgrid(lapsi4t + (ll - 2)*icycle, icycle/2, i, &
                           fcstgrd(1), fcstgrd(2), fcstgrd(3), &
                           ext, 'ref', unit, comment, refl(1, 1, 1, ll), istatus)
      if (i .ne. lapsi4t + (ll - 2)*icycle) then
         print *, 'does not find a lps at: ', lapsi4t + (ll - 2)*icycle
         cycle
      end if
   end do
   ! interpolate:
   refl(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) = &
      0.5*(refl(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) + &
           refl(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2))
   refl(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 3) = &
      2.0*refl(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) - &
      refl(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) ! extrapolation

   ! reflectivity derived bounds:
   bk0(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1:fcstgrd(4), numstat + 1) = 0.0
   rmax(1) = 45.0
   rlow(1) = 0.5
   rmax(2) = 30.0
   rlow(2) = 0.6
   do ll = 1, fcstgrd(4)
   do k = 1, fcstgrd(3)
   do j = 1, fcstgrd(2)
   do i = 1, fcstgrd(1)
      if (refl(i, j, k, ll) .gt. 5.0 .and. refl(i, j, k, ll) .le. 100) then

         ! iqc: 1 rain; 2 snow:
         iqc = 2
         if (bk0(i, j, k, ll, temprtur) .gt. 273.15) iqc = 1

         if (refl(i, j, k, ll) .gt. rmax(iqc)) then
            bk0(i, j, k, ll, numstat + 1) = &
               make_ssh(z_fcstgd(k)/100.0, bk0(i, j, k, ll, temprtur) - 273.15, 1.0, tref)
         else
            bk0(i, j, k, ll, numstat + 1) = &
               make_ssh(z_fcstgd(k)/100.0, bk0(i, j, k, ll, temprtur) - 273.15, &
                        1.0, tref)
            !          rlow(iqc)+(1.0-rlow(iqc))*refl(i,j,k,ll)/rmax(iqc),tref)
         end if
      end if
   end do
   end do
   end do
   end do

   ! read laps cloud fraction on laps cld_hgt coordinate and then interpolate to pressure levels:
   ext = 'lc3 '
   cldf = 0.0
   do ll = 1, 2 ! lc3 at current and previous time frames
      call get_clouds_3dgrid(lapsi4t + (ll - 2)*icycle, i, &
                             fcstgrd(1), fcstgrd(2), kcloud, &
                             ext, cldh(1, 1, 1, ll), cld_hgt, cld_prs, istatus)
      if (i .ne. lapsi4t + (ll - 2)*icycle) then
         print *, 'does not find a lc3 at: ', lapsi4t + (ll - 2)*icycle
         cycle
      end if

      ! interpolate from laps_height to analysis pressure coordinate:
      print *, 'interpolating lc3 to lcp...'
      call interp_height_pres_fast(fcstgrd(1), fcstgrd(2), fcstgrd(3), kcloud, cldf(1, 1, 1, ll), &
                                   cldh(1, 1, 1, ll), bk0(1, 1, 1, ll, pressure), cld_hgt, istatus)

   end do
   print *, 'cloud fraction before interpolation: ', maxval(cldf(:, :, :, 1)), minval(cldf(:, :, :, 1))

   ! interpolate cloud fraction to all time frames:
   cldf(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) = &
      0.5*(cldf(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) + &
           cldf(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2))
   cldf(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 3) = &
      2.0*cldf(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 2) - &
      cldf(1:fcstgrd(1), 1:fcstgrd(2), 1:fcstgrd(3), 1) ! extrapolation

   ! adjust bound based on cloud:
   reflectivity_bound = .false.
   do ll = 1, fcstgrd(4)
   do k = 1, fcstgrd(3)
   do j = 1, fcstgrd(2)
   do i = 1, fcstgrd(1)
      ! temporarily remove sh lower bounds where there is no cloud:
      if (cldf(i, j, k, ll) .lt. 0.1 .and. refl(i, j, k, ll) .gt. 5.0) &
         bk0(i, j, k, ll, numstat + 1) = 0.0

      ! covering the following 5 lines: testing no bound for cloudy area without reflectivity:
      !if (cldf(i,j,k,ll) .ge. 0.1) then
      !  rhc = make_ssh(z_fcstgd(k)/100.0,bk0(i,j,k,ll,temprtur)-273.15, &
      !                 cldf(i,j,k,ll)**0.2,0.0)
      !  bk0(i,j,k,ll,numstat+1) = amax1(bk0(i,j,k,ll,numstat+1),rhc)
      !endif

      ! rh =100% if both cloud and reflectivity occur:
      if (cldf(i, j, k, ll) .ge. 0.1 .and. refl(i, j, k, ll) .ge. 5.0) &
         bk0(i, j, k, ll, numstat + 1) = &
         make_ssh(z_fcstgd(k)/100.0, bk0(i, j, k, ll, temprtur) - 273.15, 1.0, tref)
      if (ll .eq. 2 .and. cldf(i, j, k, ll) .ge. 0.1 .and. refl(i, j, k, ll) .ge. 5.0) &
         reflectivity_bound = .true.

   end do
   end do
   end do
   end do
   print *, 'reflectivity bound derived: ', reflectivity_bound

   ! convert derived bounds for sh to 'rh' = sh/s2r(p) by yuanfu xie dec. 2013:
   do ll = 1, fcstgrd(4)
   do k = 1, fcstgrd(3)
   do j = 1, fcstgrd(2)
   do i = 1, fcstgrd(1)
      bk0(i, j, k, ll, numstat + 1) = bk0(i, j, k, ll, numstat + 1)/s2r(z_fcstgd(k)/100.0)
   end do
   end do
   end do
   end do

   ! assign reflectivity derived bounds to grdbkgd0:
   do j = 1, maxgrid(2)
      iy0 = float(j - 1)/float(maxgrid(2) - 1)*(fcstgrd(2) - 1) + 1
      do i = 1, maxgrid(1)
         ix0 = float(i - 1)/float(maxgrid(1) - 1)*(fcstgrd(1) - 1) + 1
         do ll = 1, maxgrid(4)
         do k = 1, maxgrid(3)
            grdbkgd0(i, j, k, ll, numstat + 1) = bk0(ix0, iy0, k, ll, numstat + 1)
            grdbkgd0(i, j, k, ll, numstat + 2) = 1000.0 ! no uppper bound now
         end do
         end do
      end do
   end do
   print *, 'max dbz over finest grid: ', maxval(grdbkgd0(:, :, :, :, numstat + 1)), &
      minval(grdbkgd0(:, :, :, :, numstat + 1))
   ! deallocate:

   deallocate (refl, cldf, cldh, stat=istatus)

end subroutine reflectivity

