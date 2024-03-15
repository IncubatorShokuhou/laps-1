program gsi2laps

!!       history:
!!       creation: lungtsung cheng    8-2006

   implicit none

!--- decalare variables ---------------------------------------

!//common variables

   integer :: imax, jmax, kmax
   integer :: imax1, jmax1, kmax1
   integer :: i, j, k, ii, jj, kk, istatus
   integer :: i4time, namelen
   character :: filename*125

!// for subroutine read_gsi_output

   integer, parameter :: ii1 = 8
   real, allocatable, dimension(:, :, :) :: uu1, vv1, tt1, qq1
   real, allocatable, dimension(:, :) :: psfc
   real, allocatable, dimension(:) :: eta, znu
   real :: ptop

!// for subroutine get_systime

   character*9 :: a9time

!// for subroutine get_pres_1d

   real, allocatable, dimension(:) :: pres_1d

!// for subroutine get_r_missing_data

   real :: r_missing_data

!// for subroutine get_pres_3d

   real, allocatable, dimension(:, :, :) :: pres_3d

!// for subroutine get_modelfg_3d

   real, allocatable, dimension(:, :, :) :: heights_3d

!// for subroutine get_modelfg_3d

   real, allocatable, dimension(:, :, :) :: t_laps_bkg

!// for subroutine dryairmass

   real, allocatable, dimension(:, :) :: dam, pdam

!// for subroutine unstagger

   real, allocatable, dimension(:, :, :) :: uu2, vv2, qq2

!// for subroutine mass2laps

   real, allocatable, dimension(:, :, :) :: tt2, uvar1, vvar1, &
                                            wvar1, tvar1, qvar1, tvar
   real, parameter :: cp = 1004.0, rc = 287.0, t0 = 300.0 !273.15

!// for subroutine read_static_grid

   real, allocatable, dimension(:, :) :: lat, lon, topo

!// for subroutine get_grid_spacing_actual

   real :: grid_spacing_m
   integer :: icen, jcen

!// for subroutine wind_post_process

   real, allocatable, dimension(:, :) :: rk_terrain, &
                                         uanl_sfcitrp, vanl_sfcitrp
   character*3 :: var_3d(3) = (/'u3', 'v3', 'om'/)
   character*4 :: units_3d(3) = (/'m/s ', 'm/s ', 'pa/s'/)
   character*125 :: comment_3d(3) = (/'3dwind', &
                                      '3dwind', '3dwind'/)
   logical, parameter :: l_grid_north_out = .true.
   real :: height_to_zcoord2, fraclow, frachigh
   integer :: klow, khigh

!// for subroutine make_fnam_lp

   character*9 :: asc9_tim

!// for subroutine write_temp_anal

   real, allocatable, dimension(:, :, :, :) ::  output_4d
   character*10 :: comment_2d

!// for subroutine writefile

   integer, allocatable, dimension(:) :: lvl
   character*125, parameter :: &
      commentline = 'maps with clouds and surface effects only'

!// for subroutine lh3_compress

   real :: t_ref

!//variables for lwm: yuanfu xie

   real, allocatable, dimension(:, :, :) :: out_sfc_2d
   character*3 :: ext_xie, var_a_xie(2), units_a_xie(2)
   character*7 :: comment_a_xie(2)

!-----> main program <-----------------------------------------------------

! to read path of wrf_inout.

   call get_directory('log', filename, namelen)        ! yuanfu: use laps;
   filename = filename(1:namelen - 4)//'tmpdirejetlaps/wrf_inout'

! to read gsi dimension
! input :: filename
! output:: imax,jmax,kmax

   call read_gsi_dim(imax, jmax, kmax, filename, istatus)
   if (istatus .ne. 0) then
      print *, 'cannot read dimensions of wrf_inout.'
      print *, 'filenam: ', filename
      call exit(1)
   end if
   imax1 = imax - 1
   jmax1 = jmax - 1
   kmax1 = kmax - 1

! to allocate variables

   allocate (output_4d(imax, jmax, kmax, 2))
   allocate (uu1(imax, jmax1, kmax1))
   allocate (uu2(imax, jmax, kmax))
   allocate (vv1(imax1, jmax, kmax1))
   allocate (vv2(imax, jmax, kmax))
   allocate (tt1(imax1, jmax1, kmax1))
   allocate (tt2(imax, jmax, kmax))
   allocate (uvar1(imax, jmax, kmax))
   allocate (vvar1(imax, jmax, kmax))
   allocate (wvar1(imax, jmax, kmax))
   allocate (tvar(imax, jmax, kmax))
   allocate (tvar1(imax, jmax, kmax))
   allocate (qq1(imax1, jmax1, kmax1))
   allocate (qq2(imax, jmax, kmax))
   allocate (qvar1(imax, jmax, kmax))
   allocate (heights_3d(imax, jmax, kmax))
   allocate (out_sfc_2d(imax, jmax, 2))
   allocate (pres_3d(imax, jmax, kmax))
   allocate (t_laps_bkg(imax, jmax, kmax))
   allocate (psfc(imax1, jmax1))
   allocate (uanl_sfcitrp(imax, jmax))
   allocate (vanl_sfcitrp(imax, jmax))
   allocate (topo(imax, jmax))
   allocate (lat(imax, jmax))
   allocate (lon(imax, jmax))
   allocate (rk_terrain(imax, jmax))
   allocate (dam(imax, jmax))
   allocate (pdam(imax, jmax))
   allocate (eta(kmax))
   allocate (znu(kmax1))
   allocate (pres_1d(kmax))
   allocate (lvl(kmax))

! to write out variables of wrf_inout(after gsi analysis)
! input :: ii1,imax,jmax,kmax,imax1,jmax1,kmax1,filename
! output:: uu1(u wind),vv1(v wind),tt1(temperature),
!          psfc(surface pressure),ptop(top pressure),
!          eta(on constant p levels),znu(eta on mass levels)),&
!          qq1(specific humidity)

   call read_gsi_output(ii1, imax, jmax, kmax, imax1, jmax1, kmax1, &
                        filename, uu1, vv1, tt1, psfc, ptop, eta, znu, qq1, istatus)
   if (istatus .ne. 0) then
      print *, 'cannot read in data of wrf_inout.'
      print *, 'filenam: ', filename
      call exit(1)
   end if

! to write out default system time
!      call get_systime(i4time,a9time,istatus)
!      if( istatus.ne.1 )then
!          print*,'cannot get system time.'
!          call exit(1)
!      endif
   call get_directory('log', filename, namelen)
   filename = filename(1:namelen)//'i4time.txt'
   ! open(10,file='i4time.txt')
   open (10, file=filename(1:namelen + 10))
   read (10, *) i4time
   close (10)

! to write out 1_d pressure

   call get_pres_1d(i4time, kmax, pres_1d, istatus)
   if (istatus .ne. 1) then
      print *, 'cannot get vertical coordinates of laps.'
      call exit(1)
   end if

! to write out background data

   call get_r_missing_data(r_missing_data, istatus)

   call get_modelfg_3d(i4time, 'ht ', imax, jmax, kmax, heights_3d, istatus)

   call get_modelfg_3d(i4time, 't3 ', imax, jmax, kmax, t_laps_bkg, istatus)

   call get_pres_3d(i4time, imax, jmax, kmax, pres_3d, istatus)

   call dryairmass(dam, pdam, imax, jmax, kmax, pres_1d, heights_3d, pres_3d, &
                   t_laps_bkg)

! to unstagger output data  of wrf_inout
! input :: imax,jmax,kmax,imax1,jmax1,kmax1,tt1,uu1,vv1,qq1
!          ptop,psfc,znu
! output(unstagger):: tt2(temperature),uu2(u wind),
!                     vv2(v wind),qq2(specific humidity)

   call unstagger(imax, jmax, kmax, imax1, jmax1, kmax1, tt1, &
                  uu1, vv1, qq1, tt2, uu2, vv2, qq2, ptop, dam, znu)

! xie: mass to laps transfer:

   call mass2laps(tvar, imax, jmax, kmax, pres_1d, dam, eta, 4, 0, tt2)
   do k = 1, kmax
      tvar1(1:imax, 1:jmax, k) = (tvar(1:imax, 1:jmax, k) + t0)* &
                                 (pres_1d(k)/100000.0)**(rc/cp)
   end do
   call mass2laps(uvar1, imax, jmax, kmax, pres_1d, dam, eta, 4, 0, uu2)
   call mass2laps(vvar1, imax, jmax, kmax, pres_1d, dam, eta, 4, 0, vv2)
   call mass2laps(qvar1, imax, jmax, kmax, pres_1d, dam, eta, 2, 0, qq2)

! for background
   call read_static_grid(imax, jmax, 'lat', lat, istatus)

   call read_static_grid(imax, jmax, 'lon', lon, istatus)

   call read_static_grid(imax, jmax, 'avg', topo, istatus)

   icen = imax/2 + 1
   jcen = jmax/2 + 1

   call get_grid_spacing_actual(lat(icen, jcen), lon(icen, jcen), &
                                grid_spacing_m, istatus)

   do j = 1, jmax
      do i = 1, imax
         rk_terrain(i, j) = height_to_zcoord2(topo(i, j), &
                                              heights_3d, imax, jmax, kmax, i, j, istatus)
         klow = max(rk_terrain(i, j), 1.)
         khigh = klow + 1
         fraclow = float(khigh) - rk_terrain(i, j)
         frachigh = 1.0 - fraclow
         if (uvar1(i, j, klow) .eq. r_missing_data .or. &
             vvar1(i, j, klow) .eq. r_missing_data .or. &
             uvar1(i, j, khigh) .eq. r_missing_data .or. &
             vvar1(i, j, khigh) .eq. r_missing_data) then
            uanl_sfcitrp(i, j) = r_missing_data
            vanl_sfcitrp(i, j) = r_missing_data
         else
            uanl_sfcitrp(i, j) = uvar1(i, j, klow)*fraclow + &
                                 uvar1(i, j, khigh)*frachigh
            vanl_sfcitrp(i, j) = vvar1(i, j, klow)*fraclow + &
                                 vvar1(i, j, khigh)*frachigh
         end if
      end do
   end do

! to write wind to lw3

! old-version-wind_post_process
!      call wind_post_process(i4time,'lw3',var_3d,units_3d,&
!           comment_3d,uvar1,vvar1,imax,jmax,kmax,3,uanl_sfcitrp,&
!                 vanl_sfcitrp,topo,lat,lon,grid_spacing_m,rk_terrain,&
!                 r_missing_data,l_grid_north_out,istatus)

! wind_post_process for new version laps
! modify : lungtsung cheng  12-2007

   call wind_post_process(i4time &   ! i
                          , uvar1, vvar1 &   ! i
                          , wvar1 &   ! o
                          , imax, jmax, kmax, 3 &   ! i
                          , heights_3d &   ! i
                          , uanl_sfcitrp, vanl_sfcitrp &   ! i
                          , topo, lat, lon, grid_spacing_m &   ! i
                          , r_missing_data, l_grid_north_out &   ! i
                          , istatus)

   print *, "wvar1=", wvar1(1, 1, 1:10)

   call write_wind_output(i4time, 'lw3', var_3d &   ! i
                          , uvar1, vvar1 &   ! i
                          , wvar1 &   ! i
                          , uanl_sfcitrp, vanl_sfcitrp &   ! i
                          , imax, jmax, kmax, 3 &   ! i
                          , r_missing_data &   ! i
                          , istatus)

   !----------------------------------------------------------
! write lwm for display lw3: yuanfu xie

!       write out derived winds file (sfc wind)

   ext_xie = 'lwm'
   var_a_xie(1) = 'su'
   var_a_xie(2) = 'sv'
   do i = 1, 2
      units_a_xie(i) = 'm/s'
      comment_a_xie(i) = 'sfcwind'
   end do

   call move(uanl_sfcitrp, out_sfc_2d(1, 1, 1), imax, jmax)
   call move(vanl_sfcitrp, out_sfc_2d(1, 1, 2), imax, jmax)

   call put_laps_multi_2d(i4time, ext_xie, var_a_xie &
                          , units_a_xie, comment_a_xie, out_sfc_2d, imax, jmax, 2, istatus)
   !----------------------------------------------------------

   call make_fnam_lp(i4time, asc9_tim, istatus)
   comment_2d(1:9) = asc9_tim

! to write temperature to lt1

   call write_temp_anal(i4time, imax, jmax, kmax, tvar1, heights_3d, &
                        output_4d, comment_2d, istatus)

!
   do k = 1, kmax
      lvl(k) = pres_1d(k)*.01
   end do

! to write specific humidity to lq3

   call writefile(i4time, commentline, lvl, qvar1, imax, jmax, kmax, istatus)

!
   t_ref = 0.0

! to write relative humidity to lh3

   call lh3_compress(qvar1, tvar1, i4time, lvl, t_ref, imax, jmax, kmax, 1, istatus)

end program

!------------------> subroutine <------------------------------------------

subroutine read_gsi_dim(imax, jmax, kmax, filename, istatus)

   implicit none
   include 'netcdf.inc'
   character*125, intent(in) :: filename
   integer, intent(out) :: imax, jmax, kmax
   integer :: kk, n(3), ncid
   integer :: istatus, vid, vtype, vn, &
              vdims(maxvdims), vnatt
   character(maxncnam) :: vname

   ncid = ncopn(filename, ncnowrit, istatus)
   if (istatus .ne. 0) return

   vid = ncvid(ncid, 'u', istatus)
   if (istatus .ne. 0) return

   call ncvinq(ncid, vid, vname, vtype, vn, vdims, vnatt, istatus)
   if (istatus .ne. 0) return

   do kk = 1, 3
      call ncdinq(ncid, vdims(kk), vname, n(kk), istatus)
      if (istatus .ne. 0) return
   end do

   imax = n(1)
   jmax = n(2) + 1
   kmax = n(3) + 1

   call ncclos(ncid, istatus)

   return
end

subroutine read_gsi_output(ii1, imax, jmax, kmax, imax1, jmax1, kmax1, &
                           filename, uu1, vv1, tt1, psfc, ptop, eta, znu, qq1, istatus)

   implicit none
   include 'netcdf.inc'
   integer, intent(in) :: ii1, imax, jmax, kmax, imax1, jmax1, kmax1
   character*125, intent(in) :: filename
   real, intent(out) :: uu1(imax, jmax1, kmax1), &
                        vv1(imax1, jmax, kmax1), &
                        tt1(imax1, jmax1, kmax1), &
                        psfc(imax1, jmax1), ptop, &
                        eta(kmax), znu(kmax1), &
                        qq1(imax1, jmax1, kmax1)
   integer :: ncid
   integer :: i, j, k, ii, jj, kk
   integer :: istatus, vid, vtype, vn, &
              vdims(maxvdims), vnatt
   integer, dimension(2) :: start2, count2
   integer, dimension(3) :: start1, count1, n
   integer, dimension(4) :: start, count
   character(maxncnam) :: vname
   character*6 :: vargsi(8) = (/'u     ', 'v     ', 't     ', 'mub   ', &
                                'p_top ', 'znw   ', 'znu   ', 'qvapor'/)

   ncid = ncopn(filename, ncnowrit, istatus)
   if (istatus .ne. 0) return

   gsi_vars: do ii = 1, ii1

      vid = ncvid(ncid, trim(vargsi(ii)), istatus)
      if (istatus .ne. 0) return

      call ncvinq(ncid, vid, vname, vtype, vn, vdims, vnatt, istatus)
      if (istatus .ne. 0) return

      select case (ii)

      case (1)
         do kk = 1, 3
            call ncdinq(ncid, vdims(kk), vname, n(kk), istatus)
            if (istatus .ne. 0) return
         end do
         start = 1
         count = (/n(1), n(2), n(3), 1/)
         call ncvgt(ncid, vid, start, count, uu1, istatus)
         if (istatus .ne. 0) return

      case (2)
         n(1:3) = (/imax1, jmax, kmax1/)
         start = 1
         count = (/n(1), n(2), n(3), 1/)
         call ncvgt(ncid, vid, start, count, vv1, istatus)
         if (istatus .ne. 0) return

      case (3)
         do kk = 1, 3
            call ncdinq(ncid, vdims(kk), vname, n(kk), istatus)
            if (istatus .ne. 0) return
         end do
         start = 1
         count = (/n(1), n(2), n(3), 1/)
         call ncvgt(ncid, vid, start, count, tt1, istatus)
         if (istatus .ne. 0) return

      case (4)
         do kk = 1, 2
            call ncdinq(ncid, vdims(kk), vname, n(kk), istatus)
            if (istatus .ne. 0) return
         end do
         start1 = 1
         count1 = (/n(1), n(2), 1/)
         call ncvgt(ncid, vid, start1, count1, psfc, istatus)
         if (istatus .ne. 0) return

      case (5)
         call ncdinq(ncid, vdims, vname, n(1), istatus)
         if (istatus .ne. 0) return
         start2 = 1
         count2 = (/n(1), 1/)
         call ncvgt(ncid, vid, start2, count2, ptop, istatus)
         if (istatus .ne. 0) return

      case (6)
         start2 = 1
         count2 = (/kmax, 1/)
         call ncvgt(ncid, vid, start2, count2, eta, istatus)
         if (istatus .ne. 0) return

      case (7)
         start2 = 1
         count2 = (/kmax1, 1/)
         call ncvgt(ncid, vid, start2, count2, znu, istatus)
         if (istatus .ne. 0) return

      case (8)
         do kk = 1, 3
            call ncdinq(ncid, vdims(kk), vname, n(kk), istatus)
            if (istatus .ne. 0) return
         end do
         start = 1
         count = (/n(1), n(2), n(3), 1/)
         call ncvgt(ncid, vid, start, count, qq1, istatus)
         if (istatus .ne. 0) return

      end select

   end do gsi_vars

   call ncclos(ncid, istatus)
   if (istatus .ne. 0) return

   return
end

subroutine unstagger(imax, jmax, kmax, imax1, jmax1, kmax1, tt1, &
                     uu1, vv1, qq1, tt2, uu2, vv2, qq2, ptop, psfc, znu)

   implicit none
   integer :: i, j, k
   integer, intent(in) :: imax, jmax, kmax, imax1, jmax1, kmax1
   real, intent(in) ::   tt1(imax1, jmax1, kmax1), &
                       uu1(imax, jmax1, kmax1), &
                       vv1(imax1, jmax, kmax1), &
                       qq1(imax1, jmax1, kmax1), &
                       ptop, psfc(imax1, jmax1), znu(kmax1)
   real, intent(out) :: tt2(imax, jmax, kmax), &
                        uu2(imax, jmax, kmax), &
                        vv2(imax, jmax, kmax), &
                        qq2(imax, jmax, kmax)
   real :: sz(imax1, jmax1, kmax), uz(imax, jmax1, kmax), &
           vz(imax1, jmax, kmax), qz(imax1, jmax1, kmax), &
           qout(imax, jmax, kmax)

! get unstagger temperature variable tt2

   call unstaggerz(tt1, imax1, jmax1, kmax, ptop, psfc, znu, 4, 0, sz)
   call untaggerxy_3d(sz, imax1, jmax1, kmax, imax, jmax, tt2)

! get unstagger wind u component uu2

   call unstaggerz(uu1, imax, jmax1, kmax, ptop, psfc, znu, 4, 0, uz)
   call untaggerxy_3d(uz, imax, jmax1, kmax, imax, jmax, uu2)

! get unstagger wind v component vv2

   call unstaggerz(vv1, imax1, jmax, kmax, ptop, psfc, znu, 4, 0, vz)
   call untaggerxy_3d(vz, imax1, jmax, kmax, imax, jmax, vv2)

! get unstagger specific humidity qq2

   call unstaggerz(qq1, imax1, jmax1, kmax, ptop, psfc, znu, 4, 0, qz)
   call untaggerxy_3d(qz, imax1, jmax1, kmax, imax, jmax, qout)
   qq2 = abs(qout)

   return
end

