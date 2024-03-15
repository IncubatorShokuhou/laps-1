module lfmgrid

   use map_utils

   implicit none

   save

   integer, parameter :: maxdomains = 10
   real, parameter :: rmsg = 1.e37

   integer :: nx, ny, nz, lx, ly, lz, nvar2d, nvar3d, nvar3dh, nvar2dout, nvar3dout

   type(proj_info) :: proj

   logical :: out_cdf = .true., out_grib = .false., out_v5d = .true.
   logical :: make_micro = .true., make_firewx = .true., make_points(maxdomains) = .false.
   logical :: verbose = .true.
   logical :: realtime = .true.
   logical :: write_to_lapsdir = .true.
   logical :: make_donefile = .true.
   logical :: large_pgrid = .false.
   logical :: large_ngrid = .false.
   logical :: l_process_p = .false.
   logical :: l_process_t = .false.
   logical :: l_process_mr = .false.
   logical :: l_process_z = .false.
   logical :: l_process_uv = .false.
   logical :: l_process_w = .false.

   integer :: domnum, fcsttime, laps_reftime, laps_valtime, precip_dt = 3600  ! in seconds
   integer :: k_micro = 999
   integer :: n3d_pts_thr = 33000000
   integer :: p3d_pts_thr = 22000000
   integer :: advection_time = 0
   integer :: i4_adv_cld = -1
   integer :: i4_adv_pcp = -1

   character(len=256) :: filename, filename0
   character(len=32)  :: mtype
   character(len=20)  :: c_m2z = 'rams'
   character(len=3)   :: domnum_fstr
   character(len=4)   :: wrf_version = '3'

! point forecast variables.

   integer :: point_tz_utcoffset = 0
   character(len=5) :: point_vent_units = 'm^2/s'
   character(len=3) :: point_tz_label = 'utc'
   character(len=3) :: point_windspd_units = 'kts'
   character(len=1) :: point_temp_units = 'f'

! native grid map specifications.

   real :: ngrid_spacingx, ngrid_spacingy, nstdlon, ntruelat1, ntruelat2
   character(len=32) :: nprojection

! laps grid map specifications.

   real :: grid_spacing, stdlon, truelat1, truelat2
   character(len=32) :: projection

!real, target, dimension(:,:,:) :: ngrid   ! native model grid
!real, target, dimension(:,:,:) :: hgrid   ! horizontally interpolated 3d grid
!real, target, dimension(:,:,:) :: pgrid   ! isobaric output grid
!real, target, dimension(:,:,:) :: sgrid   ! horizontally interpolated sfc output grid
   real, target, allocatable, dimension(:, :, :) :: ngrid   ! native model grid
   real, target, allocatable, dimension(:, :, :) :: hgrid   ! horizontally interpolated 3d grid
   real, target, allocatable, dimension(:, :, :) :: pgrid   ! isobaric output grid
   real, target, allocatable, dimension(:, :, :) :: sgrid   ! horizontally interpolated sfc output grid

! native grid variables.

   real, allocatable, dimension(:, :) :: nlat, nlon
   real, pointer, dimension(:, :) :: &
      nzsfc, npsfc, ntsfc, nmrsfc, nusfc &
      , nvsfc, nwsfc, nground_t, npblhgt, npcp_tot &
      , nlwout, nswout, nlwdown, nswdown, nshflux &
      , nlhflux, nsnowc
   real, pointer, dimension(:, :, :) :: &
      npsig, nzsig, ntsig, nmrsig, nusig &
      , nvsig, nwsig, ntkesig &
      , ncldliqmr_sig, ncldicemr_sig, nrainmr_sig, nsnowmr_sig, ngraupelmr_sig &
      , nrefl_sig

! horizontally interpolated grid variables.

   real, pointer, dimension(:, :, :) :: &
      hpsig, hzsig, htsig, hmrsig, husig &
      , hvsig, hwsig, htkesig &
      , hcldliqmr_sig, hcldicemr_sig, hrainmr_sig, hsnowmr_sig, hgraupelmr_sig &
      , hrefl_sig, hzdr_sig, hldr_sig

! laps isobaric grid variables.

   real, allocatable, dimension(:) :: lprs, lprsl
   real, pointer, dimension(:, :, :) :: &
      zprs, tprs, shprs, uprs &
      , vprs, wprs, omprs, cldliqmr_prs &
      , cldicemr_prs, rainmr_prs, snowmr_prs, graupelmr_prs, refl_prs &
      , zdr_prs, ldr_prs &
      , pcptype_prs

   integer, allocatable, dimension(:) :: lvls3d
   character(len=3), allocatable, dimension(:) :: name3d
   character(len=10), allocatable, dimension(:) :: units3d
   character(len=4), allocatable, dimension(:) :: lvltype3d
   character(len=132), allocatable, dimension(:) :: com3d

! laps surface grid variables.

   real, allocatable, dimension(:, :) :: llat, llon
   real, pointer, dimension(:, :) :: &
      zsfc, psfc, tsfc, mrsfc, usfc &
      , vsfc, wsfc, ground_t, pblhgt, rhsfc &
      , pcp_tot, thetasfc, ztw0, ztw1 &
      , thetaesfc, tdsfc, redp, pmsl, upbl &
      , vpbl, u80, v80, clwmrsfc, icemrsfc, snowmrsfc, rainmrsfc &
      , graupmrsfc, cldbase, cldtop, cldamt, ceiling &
      , intliqwater, intcldice, totpcpwater, max_refl, echo_tops, refl_sfc &
      , pcptype_sfc, pcp_inc, snow_inc, snow_tot, snow_cvr &
      , srhel, uhel, bshr, scp, cape, cin, liftedind &
      , visibility, heatind, lwout, swout, lwdown &
      , swdown, shflux, lhflux, vnt_index, ham_index &
      , hah_index, fwi_index, fwx_index, upflux, bt11u &
      , cldalb, simvis
   integer, allocatable, dimension(:) :: lvls2d
   character(len=3), allocatable, dimension(:) :: name2d
   character(len=10), allocatable, dimension(:) :: units2d
   character(len=4), allocatable, dimension(:) :: lvltype2d
   character(len=132), allocatable, dimension(:) :: com2d

! grib parameters.

   integer :: table_version = 2, center_id = 59, subcenter_id = 2, process_id = 255 &
              , funit, startbyte, nbytes, igds(18)
   integer, allocatable, dimension(:) :: param, leveltype, level1, level2 &
                                         , timerange, scalep10 &
                                         , paramua, scalep10ua
   character(len=256) :: gribfile
   logical, allocatable, dimension(:) :: gribit, gribitua

contains

!===============================================================================

   subroutine alloc_native_grid

      implicit none

      integer :: ct

      if (trim(mtype) /= 'st4') then
         nvar2d = 17
         if (.not. large_ngrid) then
            nvar3d = 7 ! add 1 for tkesig if needed
            l_process_uv = .true.
            l_process_w = .true.
         else ! large_ngrid
            l_process_w = .false.
            if (.not. large_pgrid) then
               nvar3d = 6
               l_process_uv = .true.
            else
               nvar3d = 4
               l_process_uv = .false.
            end if
         end if
         if (make_micro) then
            nvar3d = nvar3d + 5
            if (c_m2z == 'wrf') then
               nvar3d = nvar3d + 1
            end if
         end if

         allocate (nlat(nx, ny), nlon(nx, ny))

         print *, 'native grid allocation 2d/3d/2dpts/3dpts = ', nvar2d, nvar3d, nx*ny*nvar2d, nx*ny*nz*nvar3d

         allocate (ngrid(nx, ny, nvar2d + nvar3d*nz))

         ngrid = rmsg

         ct = 1

         nzsfc => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         npsfc => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         ntsfc => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nmrsfc => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nusfc => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nvsfc => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nwsfc => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nground_t => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         npblhgt => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         npcp_tot => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nlwout => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nswout => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nlwdown => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nswdown => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nshflux => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nlhflux => ngrid(1:nx, 1:ny, ct); ct = ct + 1
         nsnowc => ngrid(1:nx, 1:ny, ct); ct = ct + 1

         npsig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz; l_process_p = .true.
         nzsig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz; l_process_z = .true.
         ntsig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz; l_process_t = .true.
         nmrsig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz; l_process_mr = .true.

         if (l_process_uv) then
            nusig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
            nvsig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
         end if

         if (l_process_w) then
            nwsig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
         end if

         if (make_micro) then
            k_micro = ct
            ncldliqmr_sig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
            ncldicemr_sig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
            nrainmr_sig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
            nsnowmr_sig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
            ngraupelmr_sig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
            if (c_m2z == 'wrf') then
               nrefl_sig => ngrid(1:nx, 1:ny, ct:ct + nz - 1); ct = ct + nz
            end if
         end if
      else
         nvar2d = 1
         nvar3d = 0

         allocate (nlat(nx, ny), nlon(nx, ny))

         allocate (ngrid(nx, ny, nvar2d + nvar3d*nz))

         ngrid = rmsg

         ct = 1

         npcp_tot => ngrid(1:nx, 1:ny, ct); ct = ct + 1
      end if

      return
   end subroutine

!===============================================================================

   subroutine alloc_hinterp_grid

      implicit none

      integer :: ct

      if (.not. large_ngrid) then ! state variables
         nvar3dh = 7
      else
         if (l_process_uv) then ! process u,v
            nvar3dh = 6
         else
            nvar3dh = 4
         end if
      end if

      if (make_micro) then
         nvar3dh = nvar3dh + 8
      end if

      print *, 'hinterp grid allocation 3d/3dpts = ', nvar3dh, lx*ly*nz*nvar3dh
      print *, 'ctmax (predicted) = ', nz*nvar3dh
      print *, 'nvar3dh / nvar3d = ', nvar3dh, nvar3d

      allocate (hgrid(lx, ly, nvar3dh*nz))

      hgrid = rmsg

      ct = 1

      hpsig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz

      hzsig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz

      htsig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
      hmrsig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz

      if (l_process_uv) then
         husig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hvsig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
!  htkesig  =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
      end if

      if (l_process_w) then
         hwsig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
      end if

      if (make_micro) then
         hcldliqmr_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hcldicemr_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hrainmr_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hsnowmr_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hgraupelmr_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hrefl_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hzdr_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
         hldr_sig => hgrid(1:lx, 1:ly, ct:ct + nz - 1); ct = ct + nz
      end if

      print *, 'ctmax (actual) = ', ct - 1
      if (nz*nvar3dh .ne. ct - 1) then
         write (6, *) ' error: ctmax actual is different from predicted'
         stop
      end if

      return
   end subroutine

!===============================================================================

   subroutine alloc_isobaric_grid

      implicit none

      integer :: ct

      if (.not. large_pgrid) then ! state variables
         nvar3dout = 7
      else
         nvar3dout = 0
      end if

      if (make_micro) then
         if (.not. large_pgrid) then
            nvar3dout = nvar3dout + 9
         else
            nvar3dout = nvar3dout + 1
         end if
      end if

      allocate (pgrid(lx, ly, nvar3dout*lz), name3d(nvar3dout*lz), units3d(nvar3dout*lz) &
                , lvltype3d(nvar3dout*lz), com3d(nvar3dout*lz), lvls3d(nvar3dout*lz))

      if (out_grib) then
         allocate (gribitua(nvar3dout), paramua(nvar3dout), scalep10ua(nvar3dout))
      end if

      pgrid = rmsg

      ct = 1

      if (.not. large_pgrid) then
         zprs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'ht '; com3d(ct:ct + lz - 1) = 'geopotential height'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
!  rhprs =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='rh3'; com3d(ct:ct+lz-1)='relative humidity'         ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
         tprs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 't3 '; com3d(ct:ct + lz - 1) = 'temperature'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
         shprs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'sh '; com3d(ct:ct + lz - 1) = 'specific humidity'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
         uprs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'u3 '; com3d(ct:ct + lz - 1) = 'u-component wind'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
         vprs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'v3 '; com3d(ct:ct + lz - 1) = 'v-component wind'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
!  tkeprs=>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='tke'; com3d(ct:ct+lz-1)='turbulent kinetic energy'  ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
         wprs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'w3 '; com3d(ct:ct + lz - 1) = 'vertical velocity'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
         omprs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'om '; com3d(ct:ct + lz - 1) = 'pressure vertical velocity'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
      end if

      if (make_micro) then
         if (.not. large_pgrid) then
            cldliqmr_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'lwc'; com3d(ct:ct + lz - 1) = 'cloud liq.'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            cldicemr_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'ice'; com3d(ct:ct + lz - 1) = 'cloud ice'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            rainmr_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'rai'; com3d(ct:ct + lz - 1) = 'rain conc.'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            snowmr_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'sno'; com3d(ct:ct + lz - 1) = 'snow conc.'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            graupelmr_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'pic'; com3d(ct:ct + lz - 1) = 'graupel conc.'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            refl_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'ref'; com3d(ct:ct + lz - 1) = 'radar ref.'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            zdr_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'zdr'; com3d(ct:ct + lz - 1) = 'radar zdr'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            ldr_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'ldr'; com3d(ct:ct + lz - 1) = 'radar ldr'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
            pcptype_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'pty'; com3d(ct:ct + lz - 1) = 'precip. type'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
         else ! allocate just reflectivity
            refl_prs => pgrid(1:lx, 1:ly, ct:ct + lz - 1); name3d(ct:ct + lz - 1) = 'ref'; com3d(ct:ct + lz - 1) = 'radar ref.'; lvls3d(ct:ct + lz - 1) = nint(lprs); ct = ct + lz
         end if
      end if

      if (out_grib) then
         ct = 1
         gribitua(ct) = .true.; paramua(ct) = 7; scalep10(ct) = 0; ct = ct + 1  ! zprs
!  gribitua(ct)=.false.;                                ; ct=ct+1  ! rhprs
         gribitua(ct) = .true.; paramua(ct) = 11; scalep10(ct) = 2; ct = ct + 1  ! tprs
         gribitua(ct) = .true.; paramua(ct) = 51; scalep10(ct) = 8; ct = ct + 1  ! shprs
         gribitua(ct) = .true.; paramua(ct) = 33; scalep10(ct) = 1; ct = ct + 1  ! uprs
         gribitua(ct) = .true.; paramua(ct) = 34; scalep10(ct) = 1; ct = ct + 1  ! vprs
!   gribitua(ct)=.true.; paramua(ct)=158; scalep10(ct)=3 ; ct=ct+1  ! tkeprs
         gribitua(ct) = .true.; paramua(ct) = 40; scalep10(ct) = 3; ct = ct + 1  ! wprs
         gribitua(ct) = .true.; paramua(ct) = 39; scalep10(ct) = 3; ct = ct + 1  ! omprs
         if (make_micro) then
            gribitua(ct) = .true.; paramua(ct) = 153; scalep10(ct) = 6; ct = ct + 1  ! cldliqmr_prs
            gribitua(ct) = .true.; paramua(ct) = 178; scalep10(ct) = 6; ct = ct + 1  ! cldicemr_prs
            gribitua(ct) = .true.; paramua(ct) = 170; scalep10(ct) = 6; ct = ct + 1  ! rainmr_prs
            gribitua(ct) = .true.; paramua(ct) = 171; scalep10(ct) = 6; ct = ct + 1  ! snowmr_prs
            gribitua(ct) = .true.; paramua(ct) = 179; scalep10(ct) = 6; ct = ct + 1  ! graupelmr_prs
            gribitua(ct) = .true.; paramua(ct) = 128; scalep10(ct) = 0; ct = ct + 1  ! refl_prs
!      gribitua(ct)=.true.; paramua(ct)=128; scalep10(ct)=0 ; ct=ct+1  ! zdr_prs
!      gribitua(ct)=.true.; paramua(ct)=128; scalep10(ct)=0 ; ct=ct+1  ! ldr_prs
            gribitua(ct) = .true.; paramua(ct) = 136; scalep10(ct) = 0; ct = ct + 1  ! pcptype_prs
         end if
      end if

      return
   end subroutine

!===============================================================================

   subroutine alloc_surface_grid

      implicit none

      integer :: ct

      if (trim(mtype) /= 'st4') then
         nvar2dout = 55
         if (make_micro) nvar2dout = nvar2dout + 5
         if (make_firewx) nvar2dout = nvar2dout + 7

         allocate (sgrid(lx, ly, nvar2dout), name2d(nvar2dout), units2d(nvar2dout) &
                   , lvltype2d(nvar2dout), com2d(nvar2dout), lvls2d(nvar2dout))

         if (out_grib) then
            allocate (gribit(nvar2dout), param(nvar2dout), leveltype(nvar2dout) &
                      , level1(nvar2dout), level2(nvar2dout), timerange(nvar2dout), scalep10(nvar2dout))
         end if

         sgrid = rmsg

         ct = 1

! the order of the first set of variables needs to match the order of
!   the surface variables in the native grid.

         zsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'ter'; com2d(ct) = 'model terrain'; ct = ct + 1
         psfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'psf'; com2d(ct) = 'surface pressure'; ct = ct + 1
         tsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'tsf'; com2d(ct) = 'sfc temperature'; ct = ct + 1
         mrsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'msf'; com2d(ct) = 'sfc mixing ratio'; ct = ct + 1
         usfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'usf'; com2d(ct) = 'sfc u-component wind'; ct = ct + 1
         vsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'vsf'; com2d(ct) = 'sfc v-component wind'; ct = ct + 1
         wsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'wsf'; com2d(ct) = 'sfc vertical velocity'; ct = ct + 1
         ground_t => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'tgd'; com2d(ct) = 'ground temperature'; ct = ct + 1
         pblhgt => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'blh'; com2d(ct) = 'boundary layer depth'; ct = ct + 1
         pcp_tot => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'rto'; com2d(ct) = 'run-total liq. precip accum'; ct = ct + 1
         lwout => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lwo'; com2d(ct) = 'outgoing lw radiation'; ct = ct + 1
         swout => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'swo'; com2d(ct) = 'outgoing sw radiation'; ct = ct + 1
         lwdown => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lwi'; com2d(ct) = 'incoming lw radiation'; ct = ct + 1
         swdown => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'swi'; com2d(ct) = 'incoming sw radiation'; ct = ct + 1
         shflux => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'shf'; com2d(ct) = 'sensible heat flux'; ct = ct + 1
         lhflux => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lhf'; com2d(ct) = 'latent heat flux'; ct = ct + 1
         upflux => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'umf'; com2d(ct) = 'upslope moisture flux'; ct = ct + 1
         bt11u => sgrid(1:lx, 1:ly, ct); name2d(ct) = 's8a'; com2d(ct) = '11u brightness temperature'; ct = ct + 1
         simvis => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'smv'; com2d(ct) = 'albedo (vis satellite)'; ct = ct + 1

         thetasfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'th '; com2d(ct) = 'sfc potential temperature'; ct = ct + 1
         thetaesfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'the'; com2d(ct) = 'sfc equiv. potential temperature'; ct = ct + 1
         rhsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'rh '; com2d(ct) = 'sfc relative humidity'; ct = ct + 1
         tdsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'dsf'; com2d(ct) = 'sfc dewpoint temperature'; ct = ct + 1
         redp => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'p  '; com2d(ct) = 'reduced pressure'; ct = ct + 1
         pmsl => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'slp'; com2d(ct) = 'sea-level pressure'; ct = ct + 1
         ztw0 => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'tw0'; com2d(ct) = 'height of wet-bulb = zero'; ct = ct + 1
         ztw1 => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'tw1'; com2d(ct) = 'height of wet-bulb = 1.3'; ct = ct + 1
         cldbase => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lcb'; com2d(ct) = 'cloud base asl'; ct = ct + 1
         cldtop => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lct'; com2d(ct) = 'cloud top asl'; ct = ct + 1
         cldamt => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lcv'; com2d(ct) = 'cloud opacity'; ct = ct + 1
         cldalb => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'cla'; com2d(ct) = 'cloud albedo'; ct = ct + 1
         ceiling => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'cce'; com2d(ct) = 'cloud ceiling agl'; ct = ct + 1
         intliqwater => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lil'; com2d(ct) = 'integrated cloud liquid'; ct = ct + 1
         intcldice => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lic'; com2d(ct) = 'integrated cloud ice'; ct = ct + 1
! intgraupel =>sgrid(1:lx,1:ly,ct); name2d(ct)='lig'; com2d(ct)='integrated graupel'              ; ct=ct+1
         totpcpwater => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'tpw'; com2d(ct) = 'total precipitable water'; ct = ct + 1
         max_refl => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lmr'; com2d(ct) = 'composite reflectivity'; ct = ct + 1
         echo_tops => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lmt'; com2d(ct) = 'radar echo tops'; ct = ct + 1
         refl_sfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'llr'; com2d(ct) = '1km agl reflectivity'; ct = ct + 1
         pcptype_sfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'spt'; com2d(ct) = 'sfc precip. type'; ct = ct + 1
         pcp_inc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'r01'; com2d(ct) = 'incremental tot. liq. precip'; ct = ct + 1
         snow_inc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 's01'; com2d(ct) = 'incremental snow depth'; ct = ct + 1
         snow_tot => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'sto'; com2d(ct) = 'run-total snow accum'; ct = ct + 1
         snow_cvr => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'sc'; com2d(ct) = 'snow cover fraction'; ct = ct + 1
         srhel => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'lhe'; com2d(ct) = 'storm relative helicity'; ct = ct + 1
         uhel => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'uhe'; com2d(ct) = 'updraft helicity'; ct = ct + 1
         cape => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'pbe'; com2d(ct) = 'cape'; ct = ct + 1
         cin => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'nbe'; com2d(ct) = 'cin'; ct = ct + 1
         bshr => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'bsh'; com2d(ct) = 'bulk shear'; ct = ct + 1
         scp => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'scp'; com2d(ct) = 'supercell composite parameter'; ct = ct + 1
         liftedind => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'li '; com2d(ct) = 'lifted index'; ct = ct + 1
         visibility => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'vis'; com2d(ct) = 'sfc. visibility'; ct = ct + 1
         heatind => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'hi '; com2d(ct) = 'heat index'; ct = ct + 1
         u80 => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'u80'; com2d(ct) = 'u-component wind at 80m'; ct = ct + 1
         v80 => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'v80'; com2d(ct) = 'v-component wind at 80m'; ct = ct + 1

         if (make_micro) then
            clwmrsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'scl'; com2d(ct) = 'sfc cloud liq water mr'; ct = ct + 1
            icemrsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'sic'; com2d(ct) = 'sfc cloud ice mr'; ct = ct + 1
            snowmrsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'ssn'; com2d(ct) = 'sfc prec. snow mr'; ct = ct + 1
            rainmrsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'srn'; com2d(ct) = 'sfc prec. rain mr'; ct = ct + 1
            graupmrsfc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'sgr'; com2d(ct) = 'sfc prec. graupel mr'; ct = ct + 1
         end if

         if (make_firewx) then
            upbl => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'upb'; com2d(ct) = 'u-component wind in pbl'; ct = ct + 1
            vpbl => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'vpb'; com2d(ct) = 'v-component wind in pbl'; ct = ct + 1
            vnt_index => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'vnt'; com2d(ct) = 'ventilation index'; ct = ct + 1
            ham_index => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'ham'; com2d(ct) = 'mid-level haines index'; ct = ct + 1
            hah_index => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'hah'; com2d(ct) = 'high-level haines index'; ct = ct + 1
            fwi_index => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'fwi'; com2d(ct) = 'fosberg fire wx index'; ct = ct + 1
            fwx_index => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'fwx'; com2d(ct) = 'laps/kelsch fire wx index'; ct = ct + 1
         end if

         if (out_grib) then
            ct = 1
            gribit(ct) = .true.; param(ct) = 7; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! zsfc
            gribit(ct) = .true.; param(ct) = 1; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! psfc
            gribit(ct) = .true.; param(ct) = 11; leveltype(ct) = 105; level1(ct) = 2; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! tsfc
            gribit(ct) = .false.; ct = ct + 1  ! mrsfc
            gribit(ct) = .true.; param(ct) = 33; leveltype(ct) = 105; level1(ct) = 10; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! usfc
            gribit(ct) = .true.; param(ct) = 34; leveltype(ct) = 105; level1(ct) = 10; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! vsfc
            gribit(ct) = .true.; param(ct) = 40; leveltype(ct) = 105; level1(ct) = 10; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 4; ct = ct + 1  ! wsfc
            gribit(ct) = .true.; param(ct) = 11; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! ground_t
            gribit(ct) = .true.; param(ct) = 221; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! pblhgt
            gribit(ct) = .true.; param(ct) = 142; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 4; scalep10(ct) = 1; ct = ct + 1  ! pcp_tot
            gribit(ct) = .true.; param(ct) = 212; leveltype(ct) = 8; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! lwout
            gribit(ct) = .true.; param(ct) = 211; leveltype(ct) = 8; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! swout
            gribit(ct) = .true.; param(ct) = 112; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! lwdown
            gribit(ct) = .true.; param(ct) = 111; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! swdown
            gribit(ct) = .true.; param(ct) = 122; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! shflux
            gribit(ct) = .true.; param(ct) = 121; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! lhflux
!     gribit(ct)=.true.; param(ct)=121; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! upflux

            gribit(ct) = .true.; param(ct) = 13; leveltype(ct) = 105; level1(ct) = 2; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! thetasfc
            gribit(ct) = .true.; param(ct) = 14; leveltype(ct) = 105; level1(ct) = 2; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! thetaesfc
            gribit(ct) = .true.; param(ct) = 52; leveltype(ct) = 105; level1(ct) = 2; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! rhsfc
            gribit(ct) = .true.; param(ct) = 17; leveltype(ct) = 105; level1(ct) = 2; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! tdsfc
            gribit(ct) = .false.; ct = ct + 1  ! redp
            gribit(ct) = .true.; param(ct) = 2; leveltype(ct) = 102; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! pmsl
            gribit(ct) = .false.; ct = ct + 1  ! ztw0
            gribit(ct) = .false.; ct = ct + 1  ! ztw1
            gribit(ct) = .true.; param(ct) = 138; leveltype(ct) = 102; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! cldbase
            gribit(ct) = .true.; param(ct) = 139; leveltype(ct) = 102; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! cldtop
            gribit(ct) = .true.; param(ct) = 71; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! cldamt
            gribit(ct) = .true.; param(ct) = 137; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! ceiling
            gribit(ct) = .false.; ct = ct + 1  ! intliqwater
            gribit(ct) = .true.; param(ct) = 54; leveltype(ct) = 200; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! totpcpwater
            gribit(ct) = .true.; param(ct) = 129; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! max_refl
            gribit(ct) = .true.; param(ct) = 130; leveltype(ct) = 102; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! echo_tops
            gribit(ct) = .true.; param(ct) = 128; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! refl_sfc
            gribit(ct) = .true.; param(ct) = 136; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! pcptype_sfc
            gribit(ct) = .true.; param(ct) = 61; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 4; scalep10(ct) = 1; ct = ct + 1  ! pcp_inc
            gribit(ct) = .true.; param(ct) = 66; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 4; scalep10(ct) = 4; ct = ct + 1  ! snow_inc
            gribit(ct) = .true.; param(ct) = 141; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 4; scalep10(ct) = 4; ct = ct + 1  ! snow_tot
            gribit(ct) = .true.; param(ct) = 190; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! srhel
            gribit(ct) = .true.; param(ct) = 157; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! cape
            gribit(ct) = .true.; param(ct) = 156; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 0; ct = ct + 1  ! cin
            gribit(ct) = .true.; param(ct) = 131; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 1; ct = ct + 1  ! liftedind
            gribit(ct) = .true.; param(ct) = 20; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = -3; ct = ct + 1  ! visibility
            gribit(ct) = .false.; ct = ct + 1  ! heatind

            if (make_micro) then
               gribit(ct) = .true.; param(ct) = 153; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 4; ct = ct + 1  ! clwmrsfc
               gribit(ct) = .true.; param(ct) = 178; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 4; ct = ct + 1  ! icemrsfc
               gribit(ct) = .true.; param(ct) = 171; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 4; ct = ct + 1  ! snowmrsfc
               gribit(ct) = .true.; param(ct) = 170; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 4; ct = ct + 1  ! rainmrsfc
               gribit(ct) = .true.; param(ct) = 179; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 0; scalep10(ct) = 4; ct = ct + 1  ! graupmrsfc
            end if

            if (make_firewx) then
               gribit(ct) = .false.; ct = ct + 1  ! upbl
               gribit(ct) = .false.; ct = ct + 1  ! vpbl
               gribit(ct) = .false.; ct = ct + 1  ! vnt_index
               gribit(ct) = .false.; ct = ct + 1  ! ham_index
               gribit(ct) = .false.; ct = ct + 1  ! hah_index
               gribit(ct) = .false.; ct = ct + 1  ! fwi_index
               gribit(ct) = .false.; ct = ct + 1  ! fwx_index
            end if
         end if
      else
         nvar2dout = 1

         allocate (sgrid(lx, ly, nvar2dout), name2d(nvar2dout), units2d(nvar2dout) &
                   , lvltype2d(nvar2dout), com2d(nvar2dout), lvls2d(nvar2dout))

         if (out_grib) then
            allocate (gribit(nvar2dout), param(nvar2dout), leveltype(nvar2dout) &
                      , level1(nvar2dout), level2(nvar2dout), timerange(nvar2dout), scalep10(nvar2dout))
         end if

         sgrid = rmsg

         ct = 1
         pcp_inc => sgrid(1:lx, 1:ly, ct); name2d(ct) = 'r01'; com2d(ct) = 'incremental tot. liq. precip'; ct = ct + 1
         if (out_grib) then
            ct = 1
            gribit(ct) = .true.; param(ct) = 61; leveltype(ct) = 1; level1(ct) = 0; level2(ct) = 0; timerange(ct) = 4; scalep10(ct) = 1; ct = ct + 1  ! pcp_inc
         end if
      end if

      return
   end subroutine

!===============================================================================

   subroutine dealloc_grid(gtype)

      implicit none

      character(len=*) :: gtype

      select case (trim(gtype))
      case ('native')
         deallocate (ngrid, nlat, nlon)
      case ('horiz')
         deallocate (hgrid)
      case ('isobaric')
         deallocate (pgrid, lprs, lprsl, name3d, units3d, lvltype3d, com3d, lvls3d)
      case ('surface')
         deallocate (sgrid, llat, llon, name2d, units2d, lvltype2d, com2d, lvls2d)
         if (out_grib) deallocate (gribit, param, leveltype &
                                   , level1, level2, timerange, scalep10)
      end select

      return
   end subroutine

end module
