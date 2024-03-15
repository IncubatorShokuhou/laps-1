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

module postproc_lfm

   ! contains routines necessary to post-process mm5v3 data file for laps.

   ! brent shaw, noaa/oar/fsl/frd,  dec 2000

   use mm5v3_io
   use setup
   use constants
   use vinterp_utils
   use fire

   implicit none

   private

   integer                    :: current_lun
   character(len=24)          :: time_to_proc
   integer                    :: i, j, k

   ! some variables needed on sigma for various derivations

   real, allocatable          :: psig(:, :, :)
   real, allocatable          :: tsig(:, :, :)
   real, allocatable          :: tvsig(:, :, :)
   real, allocatable          :: thetasig(:, :, :)
   real, allocatable          :: mrsig(:, :, :)
   real, allocatable          :: zsig(:, :, :)
   real, allocatable          :: rhsig(:, :, :)
   real, allocatable          :: rhodrysig(:, :, :)
   real, allocatable          :: rhomoistsig(:, :, :)
   real, allocatable          :: tvprs(:, :, :)
   real, allocatable          :: mrprs(:, :, :)
   ! some variables to hold interpolation coefficients
   integer, allocatable       :: trap_bot_ind(:, :, :)
   integer, allocatable       :: trap_top_ind(:, :, :)
   real, allocatable          :: weight_top_lin(:, :, :)
   real, allocatable          :: weight_top_log(:, :, :)

   ! output variables are on pressure or are at the surface
   real, allocatable, public  :: tprs(:, :, :)
   real, allocatable, public  :: tdprs(:, :, :)
   real, allocatable, public  :: thetaprs(:, :, :)
   real, allocatable, public  :: uprs(:, :, :)
   real, allocatable, public  :: vprs(:, :, :)
   real, allocatable, public  :: wprs(:, :, :)
   real, allocatable, public  :: omprs(:, :, :)
   real, allocatable, public  :: shprs(:, :, :)
   real, allocatable, public  :: rhprs(:, :, :)
   real, allocatable, public  :: zprs(:, :, :)
   real, allocatable, public  :: cldliqmr_prs(:, :, :)
   real, allocatable, public  :: cldicemr_prs(:, :, :)
   real, allocatable, public  :: snowmr_prs(:, :, :)
   real, allocatable, public  :: rainmr_prs(:, :, :)
   real, allocatable, public  :: graupelmr_prs(:, :, :)
   real, allocatable, public  :: refl_prs(:, :, :)
   real, allocatable, public  :: pcptype_prs(:, :, :)
   real, allocatable, public  :: abs_vort(:, :, :)
   real, allocatable, public  :: tkeprs(:, :, :)
   real, allocatable, public  :: psfc(:, :)
   real, allocatable, public  :: pmsl(:, :)
   real, allocatable, public  :: redp(:, :)
   real, allocatable, public  :: tsfc(:, :)
   real, allocatable, public  :: tdsfc(:, :)
   real, allocatable, public  :: usfc(:, :)
   real, allocatable, public  :: vsfc(:, :)
   real, allocatable, public  :: upbl(:, :)
   real, allocatable, public  :: vpbl(:, :)
   real, allocatable, public  :: wsfc(:, :)
   real, allocatable, public  :: rhsfc(:, :)
   real, allocatable, public  :: cldbase(:, :)
   real, allocatable, public  :: cldtop(:, :)
   real, allocatable, public  :: cldamt(:, :)
   real, allocatable, public  :: ceiling(:, :)
   real, allocatable, public  :: heatind(:, :)
   real, allocatable, public  :: intliqwater(:, :)
   real, allocatable, public  :: totpcpwater(:, :)
   real, allocatable, public  :: pcp_inc(:, :)
   real, allocatable, public  :: pcp_init(:, :)
   real, allocatable, public  :: con_pcp_inc(:, :)
   real, allocatable, public  :: con_pcp_init(:, :)
   real, allocatable, public  :: snow_inc(:, :)
   real, allocatable, public  :: snow_init(:, :)
   real, allocatable, public  :: pcp_tot(:, :)
   real, allocatable, public  :: con_pcp_tot(:, :)
   real, allocatable, public  :: snow_tot(:, :)
   real, allocatable, public  :: pcptype_sfc(:, :)
   real, allocatable, public  :: thetasfc(:, :)
   real, allocatable, public  :: thetaesfc(:, :)
   real, allocatable, public  :: cape(:, :)
   real, allocatable, public  :: cin(:, :)
   real, allocatable, public  :: liftedind(:, :)
   real, allocatable, public  :: srhel(:, :)
   real, allocatable, public  :: max_refl(:, :)
   real, allocatable, public  :: echo_tops(:, :)
   real, allocatable, public  :: refl_sfc(:, :)
   real, allocatable, public  :: visibility(:, :)
   real, allocatable, public  :: thick_10_5(:, :)
   real, allocatable, public  :: snowcover(:, :)
   real, allocatable, public  :: lwout(:, :)
   real, allocatable, public  :: swout(:, :)
   real, allocatable, public  :: lwdown(:, :)
   real, allocatable, public  :: swdown(:, :)
   real, allocatable, public  :: albedo(:, :)
   real, allocatable, public  :: shflux(:, :)
   real, allocatable, public  :: lhflux(:, :)
   real, allocatable, public  :: pblhgt(:, :)
   real, allocatable, public  :: ground_t(:, :)
   real, allocatable, public  :: clwmrsfc(:, :)
   real, allocatable, public  :: icemrsfc(:, :)
   real, allocatable, public  :: rainmrsfc(:, :)
   real, allocatable, public  :: snowmrsfc(:, :)
   real, allocatable, public  :: graupmrsfc(:, :)
   real, allocatable, public  :: vnt_index(:, :)
   real, allocatable, public  :: ham_index(:, :)
   real, allocatable, public  :: hah_index(:, :)
   real, allocatable, public  :: fwi_index(:, :)

   real, parameter            :: nonecode = 0.
   real, parameter            :: raincode = 1.
   real, parameter            :: snowcode = 2.
   real, parameter            :: zraincode = 3.
   real, parameter            :: sleetcode = 4.
   real, parameter            :: hailcode = 5.
   real, parameter            :: drizzlecode = 6.
   real, parameter            :: zdrizzlecode = 7.
   real, parameter            :: rainsnowcode = 8.
   real, parameter            :: rainicecode = 9.

   integer                    :: k300, k500, k700, k850, k1000

   real, external             :: dewpt
   real, external             :: relhum
   real, external             :: potential_temp
   real, external             :: eq_potential_temp
   logical                    :: initialize
   real                       :: smth
   public process_one_lfm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine process_one_lfm(lun, time)

      ! main driver subroutine for processing one time of mm5v3
      implicit none
      integer, intent(in)           :: lun
      character(len=24), intent(in)  :: time

      current_lun = lun

      ! depending on time step used in the model, there can be
      ! some rounding errors in the computation of the time
      ! string in the subheaders of the model output (i.e., it can
      ! be off by 1 second so it does not exactly match.  since the
      ! i/o routine can be called with a time string less than
      ! 1960-01-01 to tell it to use whatever it finds, we can
      ! excercise this option when using split output, since we
      ! know by the file number which output time is contained
      ! within.
      if (split_output) then
         time_to_proc = '0000-00-00_00:00:00.0000'
      else
         time_to_proc = time
      end if

      ! find some key indices of pressure levels for use in some routines
      do k = 1, kprs
         if (prslvl(k) .eq. 100000) k1000 = k
         if (prslvl(k) .eq. 85000) k850 = k
         if (prslvl(k) .eq. 70000) k700 = k
         if (prslvl(k) .eq. 50000) k500 = k
         if (prslvl(k) .eq. 30000) k300 = k
      end do

      print '(a)', 'process_one_lfm: allocating internal arrays...'
      call allocate_internal
      print '(a)', 'process_one_lfm: getting state variables on sigma...'
      call get_sigma_vars
      print '(a)', 'process_one_lfm: interp thermo variables to pressure...'
      call interp_thermo_3d
      print '(a)', 'process_one_lfm: interp momentum variables to pressure...'
      call interp_winds
      print '(a)', 'process_one_lfm: processing clouds and reflectivity...'
      call get_clouds_reflectivity
      print '(a)', 'process_one_lfm: getting precipitation amounts...'
      call get_precip
      print '(a)', 'process_one_lfm: computing stability indices...'
      call make_stability
      print '(a)', 'process_one_lfm: miscellaneous derivations...'
      call make_misc
      print '(a)', 'process_one_lfm: deallocating internal arrays...'
      call deallocate_internal

   end subroutine process_one_lfm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine allocate_internal

      ! this routine allocates all of the arrays internal to this module that
      ! may be shared among more than one routine in this module.  it also
      ! allocates the public (output) arrays if they are not already allocated
      implicit none

      allocate (psig(nx, ny, ksigh))
      allocate (tsig(nx, ny, ksigh))
      allocate (tvsig(nx, ny, ksigh))
      allocate (thetasig(nx, ny, ksigh))
      allocate (mrsig(nx, ny, ksigh))
      allocate (rhsig(nx, ny, ksigh))
      allocate (rhodrysig(nx, ny, ksigh))
      allocate (rhomoistsig(nx, ny, ksigh))
      allocate (mrprs(nx, ny, kprs))
      allocate (tvprs(nx, ny, kprs))

      if (.not. allocated(zsig)) then

         ! the non-hydrostatic version of mm5 (only version now
         ! supported since version 3) has constant height values
         ! for each sigma based on the reference state.  thus, they
         ! only need to be allocated/computed the first time through.

         initialize = .true.
         allocate (zsig(nx, ny, ksigh))
      else
         initialize = .false.
      end if
      allocate (trap_top_ind(nx, ny, kprs))
      allocate (trap_bot_ind(nx, ny, kprs))
      allocate (weight_top_lin(nx, ny, kprs))
      allocate (weight_top_log(nx, ny, kprs))

      ! allocation of public variables
      if (.not. allocated(psfc)) allocate (psfc(nx, ny))
      if (.not. allocated(tsfc)) allocate (tsfc(nx, ny))
      if (.not. allocated(thetasfc)) allocate (thetasfc(nx, ny))
      if (.not. allocated(thetaesfc)) allocate (thetaesfc(nx, ny))
      if (.not. allocated(rhsfc)) allocate (rhsfc(nx, ny))
      if (.not. allocated(tdsfc)) allocate (tdsfc(nx, ny))
      if (.not. allocated(thetaprs)) allocate (thetaprs(nx, ny, kprs))
      if (.not. allocated(zprs)) allocate (zprs(nx, ny, kprs))
      if (.not. allocated(rhprs)) allocate (rhprs(nx, ny, kprs))
      if (.not. allocated(tprs)) allocate (tprs(nx, ny, kprs))
      if (.not. allocated(tdprs)) allocate (tdprs(nx, ny, kprs))
      if (.not. allocated(shprs)) allocate (shprs(nx, ny, kprs))
      if (.not. allocated(redp)) allocate (redp(nx, ny))
      if (.not. allocated(pmsl)) allocate (pmsl(nx, ny))
      if (.not. allocated(usfc)) allocate (usfc(nx, ny))
      if (.not. allocated(uprs)) allocate (uprs(nx, ny, kprs))
      if (.not. allocated(vsfc)) allocate (vsfc(nx, ny))
      if (.not. allocated(vprs)) allocate (vprs(nx, ny, kprs))
      if (.not. allocated(wsfc)) allocate (wsfc(nx, ny))
      if (.not. allocated(wprs)) allocate (wprs(nx, ny, kprs))
      if (.not. allocated(omprs)) allocate (omprs(nx, ny, kprs))
      if (.not. allocated(cldbase)) allocate (cldbase(nx, ny))
      if (.not. allocated(cldtop)) allocate (cldtop(nx, ny))
      if (.not. allocated(cldamt)) allocate (cldamt(nx, ny))
      if (.not. allocated(ceiling)) allocate (ceiling(nx, ny))
      if (.not. allocated(intliqwater)) allocate (intliqwater(nx, ny))
      if (.not. allocated(totpcpwater)) allocate (totpcpwater(nx, ny))
      if (.not. allocated(max_refl)) allocate (max_refl(nx, ny))
      if (.not. allocated(echo_tops)) allocate (echo_tops(nx, ny))
      if (.not. allocated(cldliqmr_prs)) allocate (cldliqmr_prs(nx, ny, kprs))
      if (.not. allocated(cldicemr_prs)) allocate (cldicemr_prs(nx, ny, kprs))
      if (.not. allocated(rainmr_prs)) allocate (rainmr_prs(nx, ny, kprs))
      if (.not. allocated(snowmr_prs)) allocate (snowmr_prs(nx, ny, kprs))
      if (.not. allocated(graupelmr_prs)) allocate (graupelmr_prs(nx, ny, kprs))
      if (.not. allocated(refl_prs)) allocate (refl_prs(nx, ny, kprs))
      if (.not. allocated(refl_sfc)) allocate (refl_sfc(nx, ny))
      if (.not. allocated(pcptype_sfc)) allocate (pcptype_sfc(nx, ny))
      if (.not. allocated(pcptype_prs)) allocate (pcptype_prs(nx, ny, kprs))
      if (.not. allocated(pcp_init)) allocate (pcp_init(nx, ny))
      if (.not. allocated(pcp_inc)) allocate (pcp_inc(nx, ny))
      if (.not. allocated(pcp_tot)) allocate (pcp_tot(nx, ny))
      if (.not. allocated(con_pcp_init)) allocate (pcp_init(nx, ny))
      if (.not. allocated(con_pcp_inc)) allocate (pcp_inc(nx, ny))
      if (.not. allocated(con_pcp_tot)) allocate (pcp_tot(nx, ny))
      if (.not. allocated(snow_init)) allocate (snow_init(nx, ny))
      if (.not. allocated(snow_inc)) allocate (snow_inc(nx, ny))
      if (.not. allocated(snow_tot)) allocate (snow_tot(nx, ny))
      if (.not. allocated(srhel)) allocate (srhel(nx, ny))
      if (.not. allocated(cape)) allocate (cape(nx, ny))
      if (.not. allocated(cin)) allocate (cin(nx, ny))
      if (.not. allocated(liftedind)) allocate (liftedind(nx, ny))
      if (.not. allocated(visibility)) allocate (visibility(nx, ny))
      if (.not. allocated(heatind)) allocate (heatind(nx, ny))
      if (.not. allocated(lwout)) allocate (lwout(nx, ny))
      if (.not. allocated(swout)) allocate (swout(nx, ny))
      if (.not. allocated(lwdown)) allocate (lwdown(nx, ny))
      if (.not. allocated(swdown)) allocate (swdown(nx, ny))
      if (.not. allocated(albedo)) allocate (albedo(nx, ny))
      if (.not. allocated(shflux)) allocate (shflux(nx, ny))
      if (.not. allocated(lhflux)) allocate (lhflux(nx, ny))
      if (.not. allocated(pblhgt)) allocate (pblhgt(nx, ny))
      if (.not. allocated(ground_t)) allocate (ground_t(nx, ny))
      if (.not. allocated(clwmrsfc)) allocate (clwmrsfc(nx, ny))
      if (.not. allocated(icemrsfc)) allocate (icemrsfc(nx, ny))
      if (.not. allocated(rainmrsfc)) allocate (rainmrsfc(nx, ny))
      if (.not. allocated(snowmrsfc)) allocate (snowmrsfc(nx, ny))
      if (.not. allocated(graupmrsfc)) allocate (graupmrsfc(nx, ny))
      if ((.not. allocated(abs_vort)) .and. (make_v5d(domain_num))) &
         allocate (abs_vort(nx, ny, kprs))
      if ((.not. allocated(thick_10_5)) .and. (make_v5d(domain_num))) &
         allocate (thick_10_5(nx, ny))
      if ((.not. allocated(snowcover)) .and. (make_v5d(domain_num))) &
         allocate (snowcover(nx, ny))
      return

   end subroutine allocate_internal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_sigma_vars

      ! this subroutine populates some of the state variables needed on
      ! sigma by reading them in and/or deriving them

      implicit none
      real, allocatable             :: pstar(:, :)
      real, allocatable             :: ppsig(:, :, :)
      real, allocatable             :: pbase(:, :, :)
      real, allocatable             :: phb(:, :, :)
      real, allocatable             :: ph(:, :, :)
      real, allocatable             :: zsigf(:, :, :)
      real, allocatable             :: mu(:, :)
      real, allocatable             :: mub(:, :)
      integer                       :: status
      real                          :: tvbar
      real, external                :: mixsat
      real                          :: dz, dtdz, pbot, tvbot, zbot
      real, allocatable             :: mrsfc(:, :)
      real, external                :: tcvp
      logical                       :: made_pbl

      made_pbl = .false.
      ! get the 3d pressure

      if (mtype .eq. 'mm5') then
         ! compute pressure on sigma and sfc from the perturbation pressure and
         ! pstar

         allocate (pstar(nx, ny))
         allocate (ppsig(nx, ny, ksigh))

         call get_mm5_3d(current_lun, 'pp       ', time_to_proc, ppsig, &
                         'd    ', status)
         if (status .ne. 0) then
            print '(a)', 'problem getting peturbation pressure!  aborting...'
            call abort
         end if
         call get_mm5_2d(current_lun, 'pstarcrs ', time_to_proc, pstar, &
                         'd   ', status)
         if (status .ne. 0) then
            print '(a)', 'problem getting pstar!  aborting...'
            call abort
         end if

         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(ppsig, nx, ny, ksigh, smth)
               call smooth(pstar, nx, ny, 1, smth)
            end do
         end if

         do k = 1, ksigh
            psig(:, :, k) = pstar(:, :)*sigmah(k) + ptop + ppsig(:, :, k)
         end do
         psfc = pstar + ppsig(:, :, 1) + ptop

         ! free up memory by deallocating local var ppsig and pstar

         deallocate (ppsig)
         deallocate (pstar)

         ! get the temperature data
         call get_mm5_3d(current_lun, 't        ', time_to_proc, tsig, &
                         'd    ', status)
         if (status .ne. 0) then
            print '(a)', 'problem getting temperature!  aborting...'
            call abort
         end if

         ! get the mixing ratio data
         call get_mm5_3d(current_lun, 'q        ', time_to_proc, mrsig, &
                         'd    ', status)

         if (status .ne. 0) then
            print '(a)', 'problem getting mixing ratio!  aborting...'
            call abort
         end if

         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(tsig, nx, ny, ksigh, smth)
               call smooth(mrsig, nx, ny, ksigh, smth)
            end do
         end if
         where (mrsig .le. 0.) mrsig = 0.000001

         ! compute tv on sigma
         print *, 'computing virtual temperature on sigma layers...'
         tvsig = tsig*(1.+0.61*mrsig)
         print '(a,2f7.1)', '   min/max = ', minval(tvsig), maxval(tvsig)
         ! compute theta on sigma
         print *, 'computing theta on sigma layers...'

         thetasig = tsig*(100000./psig)**kappa
         print '(a,2f7.1)', '   min/max = ', minval(thetasig), maxval(thetasig)
         rhodrysig = psig/(r*tsig)
         rhomoistsig = psig/(r*tvsig)

         ! heights are computed using the hypsometric equation, starting with
         ! the surface pressure and terrain height and integrating upward
         ! one level at a time.  the earlier versions of this program computed
         ! the heights as a static field based on the mm5 base state per the
         ! documentation.  however, sometimes the heights at the first sigma
         ! level were slightly below ground when computed this way.  this new
         ! way gives nearly identical results, but the heights are never below
         ! ground on any sigma surface.

         print *, 'computing heights on sigma layers...'
         do j = 1, ny
            do i = 1, nx
               ! initialize the values for the bottom of the layer we
               ! are going to apply the hypsometric equation to.
               zbot = terdot(i, j)
               pbot = psfc(i, j)
               tvbot = tvsig(i, j, 1)
               do k = 1, ksigh
                  ! compute mean virtual temp for this layer
                  tvbar = 0.5*(tvsig(i, j, k) + tvbot)
                  ! use hypsometric equation to compute thickness
                  dz = r*tvbar*alog(pbot/psig(i, j, k))/grav
                  ! add thickness to height of bottom of layer
                  zsig(i, j, k) = zbot + dz

                  ! set the bottom level to the current level for next go-around.
                  zbot = zsig(i, j, k)
                  pbot = psig(i, j, k)
                  tvbot = tvsig(i, j, k)
               end do
            end do
         end do

      elseif (mtype(1:3) .eq. 'wrf') then
         ! get pressure on sigma...we do this by getting the base state
         ! and perturbation pressures and adding them together
         allocate (pbase(nx, ny, ksigh))
         allocate (ppsig(nx, ny, ksigh))
         call get_wrfnc_3d(current_lun, "pb", "a", nx, ny, ksigh, 1, pbase, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf base state pressure.'
            call abort
         end if
         call get_wrfnc_3d(current_lun, "p", "a", nx, ny, ksigh, 1, ppsig, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf perturbation pressure.'
            call abort
         end if
         psig = pbase + ppsig
         deallocate (pbase)
         deallocate (ppsig)
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(psig, nx, ny, ksigh, smth)
            end do
         end if
         ! get theta on sigma
         call get_wrfnc_3d(current_lun, "t", "a", nx, ny, ksigh, 1, thetasig, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf perturbation theta.'
            call abort
         end if
         thetasig = thetasig + 300.
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(thetasig, nx, ny, ksigh, smth)
            end do
         end if

         ! get q on sigma
         call get_wrfnc_3d(current_lun, "qvapor", "a", nx, ny, ksigh, 1, mrsig, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf mixing ratio.'
            call abort
         end if
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(mrsig, nx, ny, ksigh, smth)
            end do
         end if
         ! compute temperature on sigma
         tsig = thetasig/((100000./psig)**kappa)

         ! compute dry density
         rhodrysig = psig/(r*tsig)

         ! compute virtual temperature on sigma
         tvsig = tsig*(1.+0.61*mrsig)

         ! get gph and destagger vertically
         allocate (ph(nx, ny, ksigf))
         allocate (phb(nx, ny, ksigf))
         allocate (zsigf(nx, ny, ksigf))
         call get_wrfnc_3d(current_lun, "ph", "a", nx, ny, ksigf, 1, ph, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf perturbation geopotential.'
            call abort
         end if
         call get_wrfnc_3d(current_lun, "phb", "a", nx, ny, ksigf, 1, phb, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf base-state geopotential.'
            call abort
         end if

         zsigf = (ph + phb)/grav
         do k = 1, ksigh
            zsig(:, :, k) = 0.5*(zsigf(:, :, k) + zsigf(:, :, k + 1))
         end do

         !  check lowest level heights
         do j = 1, ny
            do i = 1, nx
               if (zsig(i, j, 1) .le. terdot(i, j)) then
                  print *, 'zsig1 < ter', i, j, zsig(i, j, 1), terdot(i, j)
                  stop
               end if
            end do
         end do
         deallocate (ph)
         deallocate (phb)
         deallocate (zsigf)
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(zsig, nx, ny, ksigh, smth)
            end do
         end if

         ! compute moist density on sigma
         rhomoistsig = psig/(r*tvsig)

         ! compute psfc
         ! get sfc dry pressure (mu+mub+ptop)
         allocate (mu(nx, ny))
         allocate (mub(nx, ny))
         call get_wrfnc_2d(current_lun, "mub", "a", nx, ny, 1, mub, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf base-state mu.'
            call abort
         end if
         call get_wrfnc_2d(current_lun, "mu", "a", nx, ny, 1, mu, status)
         if (status .ne. 0) then
            print *, 'could not properly obtain wrf perturbation mu.'
            call abort
         end if
         psfc = mu + mub
         psfc = psfc + ptop
         deallocate (mu)
         deallocate (mub)
         ! compute integrated column vapor pressure
         do j = 1, ny
            do i = 1, nx
               psfc(i, j) = psfc(i, j) + tcvp(psig(i, j, :), &
                                              mrsig(i, j, :), zsig(i, j, :), &
                                              rhomoistsig(i, j, :), ksigh)
            end do
         end do
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(psfc, nx, ny, 1, smth)
            end do
         end if

      end if

      ! compute rh on sigma levels, as it will be better to interpolate
      ! rh vertically than mixing ratio

      print *, 'computing rh on sigma...'
      do k = 1, ksigh
         do j = 1, ny
            do i = 1, nx
               rhsig(i, j, k) = relhum(tsig(i, j, k), mrsig(i, j, k), psig(i, j, k))*100.
            end do
         end do
      end do
      print *, '    min/max = ', minval(rhsig), maxval(rhsig)
      where (rhsig .gt. 100.) rhsig = 100.
      where (rhsig .lt. 1.) rhsig = 1.
      ! print diagnostics from center column
      print '(a)', 'diagnostics from domain center:'
      print '(a,f7.0)', 'terrain height at center = ', terdot(nx/2, ny/2)
      print '(a)', &
         '----------------------------------------------------------------------'
      print '(a)', &
         'level   sigma  pressure(pa)  height     t      tv     theta   rh   qv'
      print '(a)', &
         '----------------------------------------------------------------------'
      do k = 1, ksigh
         print '(i5,2x,f6.4,2x,f12.1,2x,f6.0,2x,f6.2,2x,f6.2,2x,f6.2,2x,f4.0,2x,f8.6)', &
            k, sigmah(k), psig(nx/2, ny/2, k), zsig(nx/2, ny/2, k), &
            tsig(nx/2, ny/2, k), tvsig(nx/2, ny/2, k), thetasig(nx/2, ny/2, k), &
            rhsig(nx/2, ny/2, k), mrsig(nx/2, ny/2, k)

      end do

      ! get surface temperature and moisture.  mm5v3 may have t2 and q2
      ! variables from similarity theory.  wrfv1 has th2 and q2 variables.

      allocate (mrsfc(nx, ny))
      tsfc(:, :) = 0.
      mrsfc(:, :) = 0.

      ! 2m temperature
      if (mtype .eq. 'mm5') then
         call get_mm5_2d(current_lun, 't2       ', time_to_proc, tsfc, &
                         'd   ', status)

      elseif (mtype(1:3) .eq. 'wrf') then

         ! get th2 and convert to temperature if non-zero
         call get_wrfnc_2d(current_lun, 'th2', 'a', nx, ny, 1, tsfc, status)
         thetasfc = tsfc
         if ((minval(tsfc) .gt. 100.) .and. (maxval(tsfc) .lt. 1000.)) then
            tsfc = tsfc/((100000./psfc)**kappa)
         else
            tsfc = 0.
            thetasfc = 0.
         end if
      end if

      ! 2m qvapor
      if (mtype .eq. 'mm5') then
         call get_mm5_2d(current_lun, 'q2       ', time_to_proc, mrsfc, &
                         'd   ', status)

      elseif (mtype(1:3) .eq. 'wrf') then

         ! get q2
         call get_wrfnc_2d(current_lun, 'q2', 'a', nx, ny, 1, mrsfc, status)

      end if

      ! make sure we had q2.  if not, use lowest sigma value

      if (maxval(mrsfc) .lt. 0.000001) then
         print *, 'using lowest model level mixing ratio for surface.'
         mrsfc(:, :) = mrsig(:, :, 1)
      else
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(mrsfc, nx, ny, 1, smth)
            end do
         end if
      end if

!diagnostic
      dz = zsig(nx/2, ny/2, 1) - (terdot(nx/2, ny/2) + 2.)
      dtdz = (tsig(nx/2, ny/2, 2) - tsig(nx/2, ny/2, 1))/ &
             (zsig(nx/2, ny/2, 2) - zsig(nx/2, ny/2, 1))
      print '(a,4f6.1,f10.5)', 'sfctemptest:t1 tsim texp dz dtdz =', tsig(nx/2, ny/2, 1), &
         tsfc(nx/2, ny/2), tsig(nx/2, ny/2, 1) - dtdz*dz, dz, dtdz

      ! make sure we have tsfc
      if (maxval(tsfc) .lt. 150.) then
         print '(a)', 't2 not available...will extrapolate from lowest level'

         ! assume a constant mixing ratio
         ! loop over all horizontal points for extrapolation

         do j = 1, ny
            do i = 1, nx

               ! compute dz between lowest sigma layer and 2m level
               ! added parenthesis around 2nd term -- bls 27 sep 01
               dz = zsig(i, j, 1) - (terdot(i, j) + 2.)
               ! compute the lapse rate of temp for the two lowest
               ! sigma levels.  use of the hypsometric equation
               ! did not work well for these thin layers.

               dtdz = (tsig(i, j, 2) - tsig(i, j, 1))/ &
                      (zsig(i, j, 2) - zsig(i, j, 1))
               tsfc(i, j) = tsig(i, j, 1) - dtdz*dz
               if (mtype(1:3) .eq. 'wrf') then
                  thetasfc(i, j) = potential_temp(tsfc(i, j), psfc(i, j))
               end if
            end do
         end do
      else
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(tsfc, nx, ny, 1, smth)
            end do
         end if
      end if
      ! compute some things that are derived from t and q
      print *, 'computing rh...'
      do j = 1, ny

         do i = 1, nx

            ! compute sfc relative humidity
            rhsfc(i, j) = min(relhum(tsfc(i, j), mrsfc(i, j), psfc(i, j)), 1.)*100.

            ! compute sfc dewpoint
            tdsfc(i, j) = dewpt(tsfc(i, j), rhsfc(i, j)*0.01)

            if (mtype(1:3) .ne. 'wrf') then
               ! compute theta at the surface
               thetasfc(i, j) = potential_temp(tsfc(i, j), psfc(i, j))
            end if
            ! compute thetae at the surface
            thetaesfc(i, j) = eq_potential_temp(tsfc(i, j), psfc(i, j), mrsfc(i, j), &
                                                rhsfc(i, j)*0.01)

         end do
      end do
      print *, '   min/max tsfc      = ', minval(tsfc), maxval(tsfc)
      print *, '   min/max rhsfc     = ', minval(rhsfc), maxval(rhsfc)
      print *, '   min/max tdsfc     = ', minval(tdsfc), maxval(tdsfc)
      print *, '   min/max thetasfc  = ', minval(thetasfc), maxval(thetasfc)
      print *, '   min/max thetaesfc = ', minval(thetaesfc), maxval(thetaesfc)

      print *, 'diagnostic from lowest sigma layer at domain center:'
      print '("p:",f8.1," t:",f7.1," mr:",f6.4," rh:",f5.1)', &
         psig(nx/2, ny/2, 1), tsig(nx/2, ny/2, 1), mrsig(nx/2, ny/2, 1), rhsig(nx/2, ny/2, 1)
      print *, 'diagnostic from surface at domain center:'
      print '("p:",f8.1," t:",f7.1," mr:",f6.4," rh:",f5.1," td:",f7.1)', &
         psfc(nx/2, ny/2), tsfc(nx/2, ny/2), mrsfc(nx/2, ny/2), rhsfc(nx/2, ny/2), tdsfc(nx/2, ny/2)
      deallocate (mrsfc)

      ! get a few other miscellaneous variables that were added for diagnostics

      if (mtype .eq. 'mm5') then
         call get_mm5_2d(current_lun, 'lwout    ', time_to_proc, lwout, &
                         'd   ', status)
         if (status .ne. 0) lwout(:, :) = 1.e37

         call get_mm5_2d(current_lun, 'swout    ', time_to_proc, swout, &
                         'd   ', status)
         if (status .ne. 0) swout(:, :) = 1.e37

         call get_mm5_2d(current_lun, 'lwdown   ', time_to_proc, lwdown, &
                         'd   ', status)
         if (status .ne. 0) lwdown(:, :) = 1.e37

         call get_mm5_2d(current_lun, 'swdown   ', time_to_proc, swdown, &
                         'd   ', status)
         if (status .ne. 0) swdown(:, :) = 1.e37

         call get_mm5_2d(current_lun, 'shflux   ', time_to_proc, shflux, &
                         'd   ', status)
         if (status .ne. 0) shflux(:, :) = 1.e37

         call get_mm5_2d(current_lun, 'lhflux   ', time_to_proc, lhflux, &
                         'd   ', status)
         if (status .ne. 0) lhflux(:, :) = 1.e37

         if (use_model_pbl) then
            print *, 'trying to obtain mm5 pblhgt'
            call get_mm5_2d(current_lun, 'pbl hgt  ', time_to_proc, pblhgt, &
                            'd   ', status)
            if ((status .ne. 0) .or. (maxval(pblhgt) .le. 0)) then
               print *, '  mm5 pblhgt not available, generating pblhgt with laps algorithm'
               call model_pblhgt(thetasig, thetasfc, psig, zsig, terdot, nx, ny, ksigh, pblhgt)
               made_pbl = .true.
            else
               print *, ' pblhgt found and used'
               made_pbl = .false.
            end if
         else
            print *, 'using internally generated pblhgt based on lfmpost.nl settings'
            call model_pblhgt(thetasig, thetasfc, psig, zsig, terdot, nx, ny, ksigh, pblhgt)
            made_pbl = .true.
         end if
         if (minval(pblhgt) .le. 0.) then
            print *, 'correcting pbl returned'
            do j = 1, ny
               do i = 1, nx
                  if (pblhgt(i, j) .le. 0) then
                     print *, 'pbl < 0 at i/j/val', i, j, pblhgt(i, j)
                     pblhgt(i, j) = zsig(i, j, 1)
                  end if
               end do
            end do
         end if
         print *, 'min/max pbl height: ', minval(pblhgt), maxval(pblhgt)
         print *, 'min/max zsig(:,:,1): ', minval(zsig(:, :, 1)), maxval(zsig(:, :, 1))
         call get_mm5_2d(current_lun, 'ground t ', time_to_proc, ground_t, &
                         'd   ', status)
         if (status .ne. 0) ground_t(:, :) = 1.e37

      elseif (mtype(1:3) .eq. 'wrf') then
         lwout(:, :) = 1.e37
         swout(:, :) = 1.e37
         lwdown(:, :) = 1.e37
         swdown(:, :) = 1.e37
         albedo(:, :) = 1.e37
         call get_wrfnc_2d(current_lun, 'hfx', 'a', nx, ny, 1, shflux, status)
         call get_wrfnc_2d(current_lun, 'qfx', 'a', nx, ny, 1, lhflux, status)
         call get_wrfnc_2d(current_lun, 'gsw', 'a', nx, ny, 1, swdown, status)

         ! unlike mm5, wrf downward shortwave is net...i.e., albedo
         ! effect has been removed.  lets get the wrf albedo (required
         ! making albedo an output variable in the wrf registry) and
         ! divide it back out.
         call get_wrfnc_2d(current_lun, 'albedo', 'a', nx, ny, 1, albedo, status)
         if (status .eq. 0) then
            ! make sure destaggering did not create values outside the range
            where (albedo .lt. 0.) albedo = 0.
            where (albedo .gt. 1.0) albedo = 1.0
            print '(a)', '  swdown being changed from net to total'
            print '(a,2f10.1)', 'min/max net swdown:', minval(swdown), maxval(swdown)
            swdown = swdown/(1.-albedo)
            print '(a,2f10.1)', 'min/max tot swdown:', minval(swdown), maxval(swdown)
         else
            print '(a)', 'swdown is still net value...albedo not obtained!'
         end if

         call get_wrfnc_2d(current_lun, 'glw', 'a', nx, ny, 1, lwdown, status)
         call model_pblhgt(thetasig, thetasfc, psig, zsig, terdot, nx, ny, ksigh, pblhgt)
         made_pbl = .true.
         call get_wrfnc_2d(current_lun, 'tsk', 'a', nx, ny, 1, ground_t, status)
      end if

!  commented out smoothing of these fields. 01/21/2004 bls
!  if (do_smoothing) then
!    do smth = 0.5, -0.5, -1
!       call smooth(lwout,nx,ny,1,smth)
!       call smooth(swout,nx,ny,1,smth)
!       call smooth(lwdown,nx,ny,1,smth)
!       call smooth(swdown,nx,ny,1,smth)
!       call smooth(shflux,nx,ny,1,smth)
!       call smooth(lhflux,nx,ny,1,smth)
!       if (.not. made_pbl) call smooth(pblhgt,nx,ny,1,smth)
!       call smooth(ground_t,nx,ny,1,smth)
!    enddo
!  endif

   end subroutine get_sigma_vars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine interp_thermo_3d

      ! subroutine to compute two 3-d arrays of weighting coefficients
      ! and an integer array of trapping indices that are used for
      ! the vertical interpolation from sigma to pressure

      implicit none

      real                    :: deltalnp, deltap
      real                    :: dtvdlnps, dqdp
      integer                 :: status
      integer                 :: ks, kp
      real                    :: p_lower_mb, p_upper_mb, p_mid_mb
      logical                 :: found_trap
      real                    :: tvbar, tvtop, tvbot, dz
      real, external          :: mixsat
      real                    :: weight_bot

      ! interpolate some of the basic variables to pressure levels.  while
      ! doing this lets save the trapping indices and the weight assigned
      ! to the top trapping value for use in later interpolations

      trap_bot_ind = -1
      trap_top_ind = -1
      weight_top_lin = 0.0
      weight_top_log = 0.0
      pressure_loop: do kp = 1, kprs
         ns_loop: do j = 1, ny
            ew_loop: do i = 1, nx

               ! case 1: is the pressure level below ground?
               if (prslvl(kp) .gt. psig(i, j, 1)) then
                  trap_bot_ind(i, j, kp) = 0
                  trap_top_ind(i, j, kp) = 1
                  weight_top_lin(i, j, kp) = 1.0
                  weight_top_log(i, j, kp) = 1.0

                  ! estimate the height of this level that is below
                  ! ground.  first, estimate mean virtual temperature by
                  ! using the mm5 base lapse rate to extrapolate the lowest
                  ! virtual temp (above the boundary layer) downward.

                  deltalnp = alog(prslvl(kp)) - alog(psig(i, j, 1))
                  tvbot = tvsig(i, j, 1) + deltalnp*dtdlnpbase
                  tvbar = (tvsig(i, j, 1) + tvbot)*0.5
                  tprs(i, j, kp) = tvbot/(1.+0.61*mrsig(i, j, 1))

                  ! derive relative humidity assuming constant mixing ratio
                  mrprs(i, j, kp) = mrsig(i, j, 1)
                  rhprs(i, j, kp) = min(relhum(tprs(i, j, kp), mrprs(i, j, kp), &
                                               prslvl(kp))*100., 100.)
                  thetaprs(i, j, kp) = potential_temp(tprs(i, j, kp), prslvl(kp))
                  if (abs(prslvl(kp) - psig(i, j, 1)) .lt. 0.1) then
                     zprs(i, j, kp) = zsig(i, j, 1)
                     dz = 0.0
                  else
                     dz = tvbar*rog*alog(prslvl(kp)/psig(i, j, 1))
                     zprs(i, j, kp) = zsig(i, j, 1) - dz
                  end if

                  ! case 2: is the pressure level above the model top?
               else if (prslvl(kp) .lt. psig(i, j, ksigh)) then
                  trap_bot_ind(i, j, kp) = ksigh
                  trap_top_ind(i, j, kp) = 0
                  weight_top_lin(i, j, kp) = 0.0
                  weight_top_log(i, j, kp) = 0.0

                  ! now, we simply assume we are high enough up that the
                  ! atmosphere is isothermal (above the tropopause).  for
                  ! mixing ratio, we will use half the value of either the
                  ! model top sigma layer or the next lowest pressure level,
                  ! whichever is physically higher in the atmosphere.  the
                  ! idea is to reduce the moisture toward extremely dry as
                  ! you approach the outer edges of the atmosphere.

                  tprs(i, j, kp) = tsig(i, j, ksigh)
                  if (prslvl(kp - 1) .gt. psig(i, j, ksigh)) then
                     mrprs(i, j, kp) = mrsig(i, j, ksigh)*0.5
                  else
                     mrprs(i, j, kp) = mrprs(i, j, kp - 1)*0.5
                  end if
                  tvtop = tprs(i, j, kp)*(1.+0.61*mrprs(i, j, kp))
                  tvbar = (tvsig(i, j, ksigh) + tvtop)*0.5

                  ! derive relative humidity from the new mixing ratio
                  rhprs(i, j, kp) = max(relhum(tprs(i, j, kp), mrprs(i, j, kp), &
                                               prslvl(kp))*100., 1.)
                  rhprs(i, j, kp) = min(rhprs(i, j, kp), 100.)
                  thetaprs(i, j, kp) = potential_temp(tprs(i, j, kp), prslvl(kp))
                  if (abs(prslvl(kp) - psig(i, j, ksigh)) .lt. 0.1) then
                     zprs(i, j, kp) = zsig(i, j, ksigh)
                     dz = 0.0
                  else
                     dz = tvbar*rog*alog(psig(i, j, ksigh)/prslvl(kp))
                     zprs(i, j, kp) = zsig(i, j, ksigh) + dz
                  end if

                  ! case 3: we can trap this level between 2 valid sigma levels
               else
                  sigma_loop: do ks = 1, ksigh - 1
                     if ((prslvl(kp) .le. psig(i, j, ks)) .and. &
                         (prslvl(kp) .ge. psig(i, j, ks + 1))) then

                        p_lower_mb = psig(i, j, ks)/100.
                        p_upper_mb = psig(i, j, ks + 1)/100.
                        p_mid_mb = prslvl(kp)/100.
                        call compute_lin_weights(p_mid_mb, p_lower_mb, p_upper_mb, &
                                                 weight_bot, weight_top_lin(i, j, kp))
                        call compute_log_weights(p_mid_mb, p_lower_mb, p_upper_mb, &
                                                 weight_bot, weight_top_log(i, j, kp))
                        trap_bot_ind(i, j, kp) = ks
                        trap_top_ind(i, j, kp) = ks + 1
                        tprs(i, j, kp) = weight_top_log(i, j, kp)*tsig(i, j, ks + 1) + &
                                         (1.0 - weight_top_log(i, j, kp))*tsig(i, j, ks)
                        rhprs(i, j, kp) = weight_top_lin(i, j, kp)*rhsig(i, j, ks + 1) + &
                                          (1.-weight_top_lin(i, j, kp))*rhsig(i, j, ks)
                        zprs(i, j, kp) = weight_top_log(i, j, kp)*zsig(i, j, ks + 1) + &
                                         (1.-weight_top_log(i, j, kp))*zsig(i, j, ks)

                        ! diagnose theta and mixing ratio
                        thetaprs(i, j, kp) = potential_temp(tprs(i, j, kp), prslvl(kp))
                        mrprs(i, j, kp) = mixsat(tprs(i, j, kp), prslvl(kp))* &
                                          rhprs(i, j, kp)*0.01

                        exit sigma_loop
                     end if
                  end do sigma_loop
               end if
            end do ew_loop
         end do ns_loop
      end do pressure_loop

      ! compute 1000-500mb thickness if make_v5d is set.

      if (make_v5d(domain_num)) thick_10_5 = zprs(:, :, k500) - zprs(:, :, k1000)

      ! convert theta into temperature
      !print *, 'converting interpolated theta to temp..'
      !do k = 1,kprs
      !  tprs(:,:,k) = thetaprs(:,:,k) / (p0/prslvl(k))**kappa
      !enddo
      print *, 'computing dewpoint on pressure...'
      tdprs = tprs/((-rvolv*alog(rhprs*.01)*tprs) + 1.0)

      ! compute specific humidity and tv pressure
      print *, 'computing specific humidity on pressure...'
      do k = 1, kprs
         do j = 1, ny
            do i = 1, nx
               shprs(i, j, k) = mrprs(i, j, k)/(1.+mrprs(i, j, k))
               tvprs(i, j, k) = tprs(i, j, k)*(1.+0.61*mrprs(i, j, k))
            end do
         end do
         print *, 'min/max sh at ', prslvl(k), ' = ', &
            minval(shprs(:, :, k)), maxval(shprs(:, :, k))
      end do

      ! get reduced pressure for laps usage
      print *, 'reducing pressure to ', redp_lvl, ' meters'
      call interp_press_to_z(prslvl, zprs, redp_lvl, redp, nx, ny, kprs)

      ! use same routine to interpolate sea-level pressure.  this method
      ! is used in lieu of original reduction routine, because it will
      ! keep the msl field consistent with the height field, which has been
      ! properly reduced.  it also produces a bit smoother field over the mountains.

      print *, 'reducing pressure to msl...'
      call interp_press_to_z(prslvl, zprs, 0., pmsl, nx, ny, kprs)

      print *, 'diagnostic print of interpolated values at center'
      print '(a)', &
         '--------------------------------------------------------------------'
      print '(a)', &
         'level  pressure(pa)  height  theta   rh  temp   dewpt  mr      tv'
      print '(a)', &
         '--------------------------------------------------------------------'
      do k = 1, kprs
         print '(i5,2x,f12.1,2x,f6.0,2x,f6.2,2x,f4.0,2x,f5.1,2x,f5.1,2x,f8.6,2x,f5.1)', &
            k, prslvl(k), zprs(nx/2, ny/2, k), &
            thetaprs(nx/2, ny/2, k), &
            rhprs(nx/2, ny/2, k), tprs(nx/2, ny/2, k), tdprs(nx/2, ny/2, k), &
            mrprs(nx/2, ny/2, k), tvprs(nx/2, ny/2, k)

      end do

      print *, '      min/max sea level press (pa): ', minval(pmsl), maxval(pmsl)
      print *, '      min/max reduced press (pa):     ', minval(redp), &
         maxval(redp)
      return

   end subroutine interp_thermo_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine interp_winds

      ! interpolates the momentum variables to pressure and surface

      implicit none

      real, allocatable       :: usig(:, :, :)
      real, allocatable       :: vsig(:, :, :)
      real, allocatable       :: wsig(:, :, :)
      real, allocatable       :: wsigf(:, :, :)
      real, allocatable       :: below_ground(:, :)
      real, allocatable       :: tkesig(:, :, :)
      real, allocatable       :: tkesigcomp(:, :, :)
      integer                 :: status
      real                    :: recipdx, dudy, dvdx

      allocate (below_ground(nx, ny))
      ! get u-component of wind
      allocate (usig(nx, ny, ksigh))
      if (mtype .eq. 'mm5') then
         call get_mm5_3d(current_lun, 'u        ', time_to_proc, usig, &
                         'd    ', status)
      elseif (mtype(1:3) .eq. 'wrf') then
         ! get/destagger wrf u wind
         call get_wrfnc_3d(current_lun, 'u', 'a', nx, ny, ksigh, 1, usig, status)

      end if

      if (status .ne. 0) then
         print '(a)', 'problem getting u-wind!  aborting...'
         call abort
      end if
      if (do_smoothing) then
         do smth = 0.5, -0.5, -1
            call smooth(usig, nx, ny, ksigh, smth)
         end do
      end if
      ! see if u10 is available (should be if using mrf pbl
      ! under mm5v3 r4)
      usfc(:, :) = 0.
      if (mtype .eq. 'mm5') then
         call get_mm5_2d(current_lun, 'u10      ', time_to_proc, usfc, &
                         'd   ', status)
      elseif (mtype(1:3) .eq. 'wrf') then
         ! get wrf u10
         call get_wrfnc_2d(current_lun, 'u10', 'a', nx, ny, 1, usfc, status)
      end if
      if ((status .ne. 0.) .or. &
          ((maxval(usfc) .eq. 0.) .and. (minval(usfc) .eq. 0))) then
         print *, 'u10 not available, using lowest sigma layer for usfc'
         usfc(:, :) = usig(:, :, 1)
      else
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(usfc, nx, ny, 1, smth)
            end do
         end if
      end if
      below_ground = usfc
      call vinterp_3d(usig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, uprs, &
                      nx, ny, ksigh, kprs)

      ! get v-component of wind
      allocate (vsig(nx, ny, ksigh))

      if (mtype .eq. 'mm5') then
         call get_mm5_3d(current_lun, 'v        ', time_to_proc, vsig, &
                         'd    ', status)
      elseif (mtype(1:3) .eq. 'wrf') then

         ! get wrf v
         call get_wrfnc_3d(current_lun, 'v', 'a', nx, ny, ksigh, 1, vsig, status)

      end if
      if (status .ne. 0) then
         print '(a)', 'problem getting v-wind!  aborting...'
         call abort
      end if
      if (do_smoothing) then
         do smth = 0.5, -0.5, -1
            call smooth(vsig, nx, ny, ksigh, smth)
         end do
      end if

      ! try to use v10 for surface v if available
      vsfc(:, :) = 0.
      if (mtype .eq. 'mm5') then
         call get_mm5_2d(current_lun, 'v10      ', time_to_proc, vsfc, &
                         'd   ', status)
      elseif (mtype(1:3) .eq. 'wrf') then
         ! get wrf v10
         call get_wrfnc_2d(current_lun, 'v10', 'a', nx, ny, 1, vsfc, status)
      end if
      if ((status .ne. 0.) .or. &
          ((maxval(vsfc) .eq. 0.) .and. (minval(vsfc) .eq. 0))) then
         print *, 'v10 not available, using lowest sigma layer for vsfc'
         vsfc(:, :) = vsig(:, :, 1)
      else
         if (do_smoothing) then
            do smth = 0.5, -0.5, -1
               call smooth(vsfc, nx, ny, 1, smth)
            end do
         end if
      end if
      below_ground = vsfc
      call vinterp_3d(vsig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, vprs, &
                      nx, ny, ksigh, kprs)

      ! get w-component of wind
      allocate (wsig(nx, ny, ksigh))
      allocate (wsigf(nx, ny, ksigf))

      if (mtype .eq. 'mm5') then
         call get_mm5_3d(current_lun, 'w        ', time_to_proc, wsigf, &
                         'd    ', status)
      elseif (mtype(1:3) .eq. 'wrf') then
         ! get wrf w on full layers
         call get_wrfnc_3d(current_lun, 'w', 'a', nx, ny, ksigf, 1, wsigf, status)
      end if
      if (status .ne. 0) then
         print '(a)', 'problem getting w-wind!  aborting...'
         call abort
      end if
      do k = 1, ksigh
         wsig(:, :, k) = 0.5*(wsigf(:, :, k) + wsigf(:, :, k + 1))
      end do
      deallocate (wsigf)
      wsfc(:, :) = wsig(:, :, 1)
      if (do_smoothing) then
         do smth = 0.5, -0.5, -1
            call smooth(wsig, nx, ny, ksigh, smth)
            call smooth(wsfc, nx, ny, 1, smth)
         end do
      end if
      below_ground(:, :) = 0.0
      call vinterp_3d(wsig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, wprs, &
                      nx, ny, ksigh, kprs)
      deallocate (wsig)

      ! compute omega for laps
      do k = 1, kprs
         omprs(:, :, k) = -(prslvl(k)*wprs(:, :, k)*grav)/(r*tvprs(:, :, k))
      end do

      ! compute the storm relative helicity

      call helicity(usig, vsig, zsig, terdot, nx, ny, ksigh, srhel)
      print *, 'min/max helicity: ', minval(srhel), maxval(srhel)

      ! compute ventilation index
      call ventilation(usig, vsig, zsig, pblhgt, terdot, nx, ny, ksigh, upbl, vpbl, vnt_index)
      print *, 'min/max ventilation : ', minval(vnt_index), maxval(vnt_index)

      ! compute tke
      allocate (tkesigcomp(nx, ny, ksigh))
      call compute_tke(psig, tsig, usig, vsig, zsig, terdot, nx, ny, ksigh, tkesigcomp)
      print '(a,2f10.3)', 'min/max tke computed from dtf3: ', &
         minval(tkesigcomp), maxval(tkesigcomp)
      deallocate (usig)
      deallocate (vsig)
      ! get tke
      print *, 'getting tke....'
      allocate (tkesig(nx, ny, ksigh))
      below_ground = 0.0
      if (mtype .eq. 'mm5') then
         call get_mm5_3d(current_lun, 'tke      ', time_to_proc, tkesig, &
                         'd    ', status)
      elseif (mtype(1:3) .eq. 'wrf') then
         call get_wrfnc_3d(current_lun, 'tke_myj', 'a', nx, ny, ksigf, 1, tkesig, status)
      end if
      if (status .ne. 0) then
         tkesig = tkesigcomp
      else
         print '(a,2f10.3)', 'min/max tke from model output: ', &
            minval(tkesig), maxval(tkesig)
      end if
      deallocate (tkesigcomp)
      print *, 'vertically interpolating tke'
      call vinterp_3d(tkesig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, tkeprs, &
                      nx, ny, ksigh, kprs)
      deallocate (tkesig)
      deallocate (below_ground)

      print *, '      min/max u wind (m/s):   ', minval(uprs), maxval(uprs)
      print *, '      min/max v wind (m/s):   ', minval(vprs), maxval(vprs)
      print *, '      min/max w wind (m/s):   ', minval(wprs), maxval(wprs)
      print *, '      min/max omega (pa/s):   ', minval(omprs), maxval(omprs)
      print *, '      min/max sr helicity:    ', minval(srhel), maxval(srhel)
      print *, '      min/max tke:            ', minval(tkeprs), maxval(tkeprs)
      ! if making vis5d output, compute vorticity
      if (make_v5d(domain_num)) then
         recipdx = 1.0/(2.0*grid_spacing)
         do k = 1, kprs
            do j = 2, ny - 1
               do i = 2, nx - 1
                  dvdx = (vprs(i + 1, j, k) - vprs(i - 1, j, k))*recipdx
                  dudy = (uprs(i, j + 1, k) - uprs(i, j - 1, k))*recipdx
                  abs_vort(i, j, k) = (mapfac_d(i, j)*(dvdx - dudy)) + coriolis(i, j)
               end do
            end do
         end do
      end if
      return

   end subroutine interp_winds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_clouds_reflectivity
      ! gets the cloud parameters on sigma surfaces, computes various
      ! parameters, then interpolates to sigma

      implicit none
      real, allocatable             :: below_ground(:, :)
      real, allocatable             :: cldliqmr_sig(:, :, :)
      real, allocatable             :: cldicemr_sig(:, :, :)
      real, allocatable             :: condmr_sig(:, :, :)
      real, allocatable             :: rainmr_sig(:, :, :)
      real, allocatable             :: snowmr_sig(:, :, :)
      real, allocatable             :: graupelmr_sig(:, :, :)
      real, allocatable             :: refl_sig(:, :, :)
      real                          :: cldliqthresh
      real                          :: cldicethresh
      real                          :: snowthresh
      integer                       :: status
      real                          :: zero_thresh

      ! set up the zero_thresh (threshold below which arrays
      ! are zeroed out for microphysics species) for mm5
      ! and wrf separately

      if (mtype .eq. 'mm5') then
         zero_thresh = 1.e-10
      else
         zero_thresh = 1.e-6
      end if
      ! get each of the cloud species mixing ratio arrays
      allocate (cldliqmr_sig(nx, ny, ksigh))
      allocate (cldicemr_sig(nx, ny, ksigh))
      allocate (rainmr_sig(nx, ny, ksigh))
      allocate (snowmr_sig(nx, ny, ksigh))
      allocate (graupelmr_sig(nx, ny, ksigh))

      ! cloud liquid and rain water
      if (clwflag) then
         if (mtype(1:3) .eq. 'mm5') then
            call get_mm5_3d(current_lun, 'clw      ', time_to_proc, cldliqmr_sig, &
                            'd    ', status)
         elseif (mtype(1:3) .eq. 'wrf') then
            call get_wrfnc_3d(current_lun, 'qcloud', 'a', nx, ny, ksigh, 1, cldliqmr_sig, &
                              status)
         end if
         if (status .ne. 0) then
            print '(a)', 'problem getting clw!  setting to 0.'
            cldliqmr_sig = 0.0
         end if
         where (cldliqmr_sig .lt. zero_thresh) cldliqmr_sig = 0.0
         clwmrsfc(:, :) = cldliqmr_sig(:, :, 1)

         if (mtype .eq. 'mm5') then
            call get_mm5_3d(current_lun, 'rnw      ', time_to_proc, rainmr_sig, &
                            'd    ', status)
         elseif (mtype(1:3) .eq. 'wrf') then
            call get_wrfnc_3d(current_lun, 'qrain', 'a', nx, ny, ksigh, 1, rainmr_sig, &
                              status)
         end if
         if (status .ne. 0) then
            print '(a)', 'problem getting rnw!  setting to 0.'
            rainmr_sig = 0.0
         end if
         where (rainmr_sig .lt. zero_thresh) rainmr_sig = 0.0
         rainmrsfc(:, :) = rainmr_sig(:, :, 1)

      else
         cldliqmr_sig = 0.0
         rainmr_sig = 0.0
         clwmrsfc = 0.0
         rainmrsfc = 0.0
      end if

      ! snow and ice
      if (iceflag) then
         if (mtype .eq. 'mm5') then
            call get_mm5_3d(current_lun, 'ice      ', time_to_proc, cldicemr_sig, &
                            'd    ', status)
         elseif (mtype(1:3) .eq. 'wrf') then
            call get_wrfnc_3d(current_lun, 'qice', 'a', nx, ny, ksigh, 1, cldicemr_sig, &
                              status)
         end if
         if (status .ne. 0) then
            print '(a)', 'problem getting ice!  setting to 0.'
            cldicemr_sig = 0.0
         end if
         where (cldicemr_sig .lt. zero_thresh) cldicemr_sig = 0.0
         icemrsfc(:, :) = cldicemr_sig(:, :, 1)

         if (mtype .eq. 'mm5') then
            call get_mm5_3d(current_lun, 'snow     ', time_to_proc, snowmr_sig, &
                            'd    ', status)
         elseif (mtype(1:3) .eq. 'wrf') then
            call get_wrfnc_3d(current_lun, 'qsnow', 'a', nx, ny, ksigh, 1, snowmr_sig, &
                              status)
         end if

         if (status .ne. 0) then
            print '(a)', 'problem getting snow!  setting to 0.'
            snowmr_sig = 0.0
         end if
         where (snowmr_sig .lt. zero_thresh) snowmr_sig = 0.0
         snowmrsfc(:, :) = snowmr_sig(:, :, 1)
      else
         cldicemr_sig = 0.0
         snowmr_sig = 0.0
         icemrsfc = 0.0
         snowmrsfc = 0.0
      end if

      ! graupel
      if (graupelflag) then
         if (mtype .eq. 'mm5') then
            call get_mm5_3d(current_lun, 'graupel  ', time_to_proc, graupelmr_sig, &
                            'd    ', status)
         elseif (mtype(1:3) .eq. 'wrf') then
            call get_wrfnc_3d(current_lun, 'qgraup', 'a', nx, ny, ksigh, 1, &
                              graupelmr_sig, status)
         end if
         if (status .ne. 0) then
            print '(a)', 'problem getting graupel!  setting to 0.'
            graupelmr_sig = 0.0
         end if
         where (graupelmr_sig .lt. zero_thresh) graupelmr_sig = 0.0
         graupmrsfc(:, :) = graupelmr_sig(:, :, 1)
      else
         graupelmr_sig = 0.0
         graupmrsfc = 0.0
      end if

      ! now compute the bases, tops, coverage, and ceiling

      call clouds(nx, ny, ksigh, grid_spacing, cldliqmr_sig, cldicemr_sig, &
                  snowmr_sig, zsig, terdot, cldbase, cldtop, ceiling, cldamt)

      ! compute total condensate so we can compute integrated liquid water
      allocate (condmr_sig(nx, ny, ksigh))
      condmr_sig = cldliqmr_sig + cldicemr_sig + snowmr_sig + rainmr_sig &
                   + graupelmr_sig
      print *, 'calling integrated_liquid'
      print *, 'min/max condmr_sig = ', minval(condmr_sig), maxval(condmr_sig)
      print *, 'min/max mrsig = ', minval(mrsig), maxval(mrsig)
      print *, 'min/max rhodrysig = ', minval(rhodrysig), maxval(rhodrysig)
      print *, 'min/max zsig = ', minval(zsig), maxval(zsig)
      call integrated_liquid(nx, ny, ksigh, condmr_sig, mrsig, rhodrysig, zsig, terdot, &
                             intliqwater, totpcpwater)
      deallocate (condmr_sig)

      ! compute reflectivities
      allocate (refl_sig(nx, ny, ksigh))
      call reflectivity(nx, ny, ksigh, rhomoistsig, zsig, &
                        rainmr_sig, cldicemr_sig, snowmr_sig, graupelmr_sig, &
                        refl_sig, max_refl, echo_tops)

      ! interpolate our 3d cloud mixing ratios and reflectivities to
      ! the pressure surfaces
      allocate (below_ground(nx, ny))

      below_ground(:, :) = 0.0
      call vinterp_3d(cldliqmr_sig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, cldliqmr_prs, &
                      nx, ny, ksigh, kprs)

      call vinterp_3d(cldicemr_sig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, cldicemr_prs, &
                      nx, ny, ksigh, kprs)

      call vinterp_3d(rainmr_sig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, rainmr_prs, &
                      nx, ny, ksigh, kprs)

      call vinterp_3d(snowmr_sig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, snowmr_prs, &
                      nx, ny, ksigh, kprs)

      call vinterp_3d(graupelmr_sig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, graupelmr_prs, &
                      nx, ny, ksigh, kprs)

      call vinterp_3d(refl_sig, trap_bot_ind, trap_top_ind, &
                      weight_top_lin, below_ground, refl_prs, &
                      nx, ny, ksigh, kprs)

      deallocate (below_ground)

      ! diagnostic print of top 4 layers

      print *, 'mean condensate values:'
      print *, 'level   species   min     max    mean'
      print *, '(mb)             (g/kg)  (g/kg) (g/kg)'
      print *, '------ -------- ------- ------- -------'
      do k = kprs - 3, kprs
         print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
            prslvl(k)*0.01, 'ice     ', minval(cldicemr_prs(:, :, k))*1000., &
            maxval(cldicemr_prs(:, :, k))*1000., &
            sum(cldicemr_prs(:, :, k))*1000./(nx*ny)
         print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
            prslvl(k)*0.01, 'snow    ', minval(snowmr_prs(:, :, k))*1000., &
            maxval(snowmr_prs(:, :, k))*1000., &
            sum(snowmr_prs(:, :, k))*1000./(nx*ny)
         print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
            prslvl(k)*0.01, 'clw     ', minval(cldliqmr_prs(:, :, k))*1000., &
            maxval(cldliqmr_prs(:, :, k))*1000., &
            sum(cldliqmr_prs(:, :, k))*1000./(nx*ny)
         print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
            prslvl(k)*0.01, 'rain    ', minval(rainmr_prs(:, :, k))*1000., &
            maxval(rainmr_prs(:, :, k))*1000., &
            sum(rainmr_prs(:, :, k))*1000./(nx*ny)
         print '(f6.1,1x,a8,1x,f7.4,1x,f7.4,1x,f7.4)', &
            prslvl(k)*0.01, 'graupel ', minval(graupelmr_prs(:, :, k))*1000., &
            maxval(graupelmr_prs(:, :, k))*1000., &
            sum(graupelmr_prs(:, :, k))*1000./(nx*ny)
      end do

      ! make  precip codes
      pcptype_sfc = nonecode
      pcptype_prs = nonecode
      do j = 1, ny
         do i = 1, nx
            do k = 1, kprs
               if (rainmr_prs(i, j, k) .gt. 0) then
                  if (snowmr_prs(i, j, k) .gt. 0.) then
                     if (graupelmr_prs(i, j, k) .gt. snowmr_prs(i, j, k)) then
                        pcptype_prs(i, j, k) = rainicecode
                     else
                        pcptype_prs(i, j, k) = rainsnowcode
                     end if
                  else
                     if (graupelmr_prs(i, j, k) .gt. 0) then
                        pcptype_prs(i, j, k) = rainicecode
                     else
                        pcptype_prs(i, j, k) = raincode
                     end if
                  end if
               else
                  if (snowmr_prs(i, j, k) .gt. 0) then
                     if (graupelmr_prs(i, j, k) .gt. snowmr_prs(i, j, k)) then
                        pcptype_prs(i, j, k) = sleetcode
                     else
                        pcptype_prs(i, j, k) = snowcode
                     end if
                  else
                     if (graupelmr_prs(i, j, k) .gt. 0) then
                        pcptype_prs(i, j, k) = sleetcode
                     else
                        pcptype_prs(i, j, k) = nonecode
                     end if
                  end if
               end if
            end do
            if (rainmr_sig(i, j, 1) .gt. 0) then
               if (snowmr_sig(i, j, 1) .gt. 0.) then
                  if (graupelmr_sig(i, j, 1) .gt. snowmr_sig(i, j, 1)) then
                     pcptype_sfc(i, j) = rainicecode
                  else
                     pcptype_sfc(i, j) = rainsnowcode
                  end if
               else
                  if (graupelmr_sig(i, j, 1) .gt. 0) then
                     pcptype_sfc(i, j) = rainicecode
                  else
                     pcptype_sfc(i, j) = raincode
                  end if
               end if
            else
               if (snowmr_sig(i, j, 1) .gt. 0) then
                  if (graupelmr_sig(i, j, 1) .gt. snowmr_sig(i, j, 1)) then
                     pcptype_sfc(i, j) = sleetcode
                  else
                     pcptype_sfc(i, j) = snowcode
                  end if
               else
                  if (graupelmr_sig(i, j, 1) .gt. 0) then
                     pcptype_sfc(i, j) = sleetcode
                  else
                     pcptype_sfc(i, j) = nonecode
                  end if
               end if
            end if
         end do
      end do
      ! set surface reflectivity to be the maximum value in the lowest
      ! third of the atmosphere.
      !
      ! refl_sfc(:,:) = maxval(refl_sig(:,:,1:ksigh/3),dim=3)

      ! no...go back to using only the lowest level for consistency with
      ! the precip type icons.

      refl_sfc(:, :) = refl_sig(:, :, 1)

      deallocate (refl_sig)
      deallocate (cldliqmr_sig)
      deallocate (cldicemr_sig)
      deallocate (rainmr_sig)
      deallocate (snowmr_sig)
      deallocate (graupelmr_sig)
      print *, '      min/max sfc reflectivity (dbz): ', minval(refl_sfc), &
         maxval(refl_sfc)
      print *, '      min/max 3d reflectivity  (dbz): ', minval(refl_prs), &
         maxval(refl_prs)
      print *, '      cloud base (m) at domain center: ', cldbase(nx/2, ny/2)
      print *, '      cloud top (m) at domain center:  ', cldtop(nx/2, ny/2)
      print *, '      cloud cover (fraction) at center:', cldamt(nx/2, ny/2)
      print *, '      ceiling (alg m) at center:       ', ceiling(nx/2, ny/2)
      print *, '      precip type codes at center:'
      print *, '      level(mb)       code'
      print'(6x,"    sfc  ",9x,f2.0)', pcptype_sfc(nx/2, ny/2)
      do k = 1, kprs
         print '(6x,f9.0,9x,f2.0)', prslvl(k)*0.01, pcptype_prs(nx/2, ny/2, k)
      end do
      print *, '      min/max int liquid water (kg/m2):', minval(intliqwater), &
         maxval(intliqwater)
      print *, '      min/max tot precip water (kg/m2):', minval(totpcpwater), &
         maxval(totpcpwater)

      return

   end subroutine get_clouds_reflectivity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_precip
      ! gets the total precipitation fields from the run and subtracts off
      ! intial and current total values to produce incremental amounts as well
      ! as total amounts.  also derives snowfall.

      implicit none

      real, allocatable                 :: raincon(:, :)
      real, allocatable                 :: rainnon(:, :)
      integer                           :: status
      integer, allocatable              :: fallen_precip_type(:, :)

      allocate (raincon(nx, ny))
      allocate (rainnon(nx, ny))
      allocate (fallen_precip_type(nx, ny))

      if (mtype .eq. 'mm5') then
         call get_mm5_2d(current_lun, 'rain con ', time_to_proc, raincon, &
                         'd   ', status)
      elseif (mtype(1:3) .eq. 'wrf') then
         call get_wrfnc_2d(current_lun, 'rainc', 'a', nx, ny, 1, raincon, status)
         ! convert to cm
         raincon = raincon*0.1
      end if
      if (status .ne. 0) then
         print '(a)', 'problem rain con'
         raincon = 0.0
      end if
      where (raincon .lt. .00001) raincon = 0.

      if (mtype .eq. 'mm5') then
         call get_mm5_2d(current_lun, 'rain non ', time_to_proc, rainnon, &
                         'd   ', status)
      elseif (mtype(1:3) .eq. 'wrf') then
         call get_wrfnc_2d(current_lun, 'rainnc', 'a', nx, ny, 1, rainnon, status)
         ! convert to cm
         rainnon = rainnon*0.1
      end if

      if (status .ne. 0) then
         print '(a)', 'problem rain non'
         rainnon = 0.0
      end if
      where (rainnon .lt. .00001) rainnon = 0.

      if (make_v5d(domain_num)) then
         if (mtype .eq. 'mm5') then
            call get_mm5_2d(current_lun, 'snowcovr ', time_to_proc, snowcover, &
                            'd   ', status)
         elseif (mtype(1:3) .eq. 'wrf') then
            call get_wrfnc_2d(current_lun, 'acsnow', 'a', nx, ny, 1, snowcover, status)
         end if
         if (status .ne. 0) then
            print '(a)', 'problem snowcovr'
            snowcover = 0.0
         end if
      end if

      ! convert from cm to m
      rainnon = rainnon*0.01
      raincon = raincon*0.01

      if (initialize) then
         pcp_init = raincon + rainnon
         pcp_inc = 0.0
         pcp_tot = 0.0
         con_pcp_init = raincon
         con_pcp_inc = 0.0
         con_pcp_tot = 0.0
         snow_inc = 0.0
         snow_tot = 0.0
      else
         pcp_inc = 0.0
         con_pcp_inc = 0.0
         snow_inc = 0.0
         con_pcp_inc = raincon - con_pcp_init - con_pcp_tot
         pcp_inc = raincon + rainnon - pcp_init - pcp_tot
         where (pcp_inc .lt. 0) pcp_inc = 0.
         where (con_pcp_inc .lt. 0) con_pcp_inc = 0.
         pcp_tot = pcp_tot + pcp_inc
         con_pcp_tot = con_pcp_tot + con_pcp_inc
         call wintprec(tsig, zsig, zprs, psfc, tsfc, terdot, pcp_inc, nx, ny, &
                       ksigh, kprs, k700, k850, k1000, raincon, &
                       fallen_precip_type)
         call snowfall(tsfc, pcp_inc, fallen_precip_type, nx, ny, &
                       snow_inc, snow_tot)
      end if
      deallocate (raincon)
      deallocate (rainnon)
      deallocate (fallen_precip_type)
      print *, '      min/max inc. precip (m): ', minval(pcp_inc), maxval(pcp_inc)
      print *, '      min/max conv. precip   : ', minval(con_pcp_inc), &
         maxval(con_pcp_inc)
      print *, '      min/max acc. precip (m): ', minval(pcp_tot), maxval(pcp_tot)
      print *, '      min/max conv. precip   : ', minval(con_pcp_tot), &
         maxval(con_pcp_tot)
      print *, '      min/max inc. snow   (m): ', minval(snow_inc), maxval(snow_inc)
      print *, '      min/max acc. snow   (m): ', minval(snow_tot), maxval(snow_tot)

      return
   end subroutine get_precip
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine make_stability

      ! make severe weather stability indices.  currently just creates cape, cin,
      ! and lifted index

      implicit none
      real, allocatable                     :: thetaesig(:, :, :)
      real, allocatable                     :: zsigf(:, :, :)
      integer                               :: itest, jtest
      allocate (thetaesig(nx, ny, ksigh))
      allocate (zsigf(nx, ny, ksigf))

      itest = nx/2
      jtest = ny/2

      ! compute thetae on sigma
      print *, ' --- computing thetae on sigma --- '
      do k = 1, ksigh
         do j = 1, ny
            do i = 1, nx
               thetaesig(i, j, k) = eq_potential_temp(tsig(i, j, k), psig(i, j, k), &
                                                      mrsig(i, j, k), rhsig(i, j, k)*0.01)
            end do
         end do
      end do
      ! compute heights on full sigma surfaces
      print *, ' --- computing heights on full sigmas ---'
      do k = 2, ksigh
         zsigf(:, :, k) = (zsig(:, :, k - 1) + zsig(:, :, k))*0.5
      end do
      zsigf(:, :, 1) = terdot
      zsigf(:, :, ksigf) = 2*zsig(:, :, ksigh) - zsigf(:, :, ksigh)

      ! call routine to compute cape, cin, and li

      print *, ' --- computing cape/cin/li ---'
      call capecin(psig*0.01, tsig, thetaesig, thetasig, rhsig*0.01, &
                   zsigf, tprs, liftedind, cape, cin, k500, nx, ny, ksigh, kprs)
      cin = -cin
      print *, 'stability stuff'
      print *, 'k   pressure  t      theta  thetae rh   z'
      do k = 1, ksigh
         print '(i4,f10.1,3f7.1,f5.0,f8.0)', &
            k, psig(itest, jtest, k), tsig(itest, jtest, k), thetasig(itest, jtest, k), &
            thetaesig(itest, jtest, k), rhsig(itest, jtest, k), zsigf(itest, jtest, k)
      end do
      print *, 'cape =', cape(itest, jtest)
      print *, 'cin = ', cin(itest, jtest)
      print *, 'li = ', liftedind(itest, jtest)
      deallocate (thetaesig)
      deallocate (zsigf)

      print *, '      min/max cape (j/kg):  ', minval(cape), maxval(cape)
      print *, '      min/max cin  (j/kg):  ', minval(cin), maxval(cin)
      print *, '      min/max li   (k):     ', minval(liftedind), maxval(liftedind)
      return
   end subroutine make_stability
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine make_misc
      ! makes miscellaneous derived parameters

      implicit none
      real, external       :: heatindex
      real, allocatable    :: tdsig(:, :, :)
      ! compute surface visibility from relative humidity and temp/dewpoint

      visibility = 6000.0*(tsfc - tdsfc)/(rhsfc**1.75)

      ! convert to meters from km
      visibility = visibility*1000.
      where (visibility .gt. 99990.) visibility = 99990.
      do j = 1, ny
         do i = 1, nx
            ! compute heat index if temp is above 80f (300k)
            if (tsfc(i, j) .ge. 300.) then
               heatind(i, j) = heatindex(tsfc(i, j), rhsfc(i, j))
            else
               heatind(i, j) = tsfc(i, j)
            end if

         end do
      end do

      ! compute fire indices

      ! we need dewpoint on sigma for the haines indices
      allocate (tdsig(nx, ny, ksigh))
      tdsig = tsig/((-rvolv*alog(rhsig*0.01)*tsig) + 1.0)

      ! mid-level haines index
      call haines_layer(psig*0.01, tsig, tdsig, ham_index, nx, ny, ksigh, &
                        850., 700.)

      ! high-level haines index
      call haines_layer(psig*0.01, tsig, tdsig, hah_index, nx, ny, ksigh, &
                        700., 500.)

      deallocate (tdsig)

      ! fosberg fwi

      call fosberg_fwi(tsfc, rhsfc, psfc*0.01, usfc, vsfc, nx, ny, fwi_index)
      print *, '      min/max surface visibility (m): ', minval(visibility), &
         maxval(visibility)
      print *, '      min/max surface heat index (k): ', minval(heatind), &
         maxval(heatind)
      return
   end subroutine make_misc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine deallocate_internal

      ! this routine deallocates all of the arrays internal to this module that
      ! may be shared among more than one routine in this module.  it does not
      ! deallocate the zsig array, however, which gets reused with each
      ! time step!!

      implicit none

      deallocate (psig)
      deallocate (tsig)
      deallocate (tvsig)
      deallocate (thetasig)
      deallocate (mrsig)
      deallocate (rhsig)
      deallocate (rhodrysig)
      deallocate (rhomoistsig)
      deallocate (trap_top_ind)
      deallocate (trap_bot_ind)
      deallocate (weight_top_lin)
      deallocate (weight_top_log)
      deallocate (tvprs)
      deallocate (mrprs)
      return

   end subroutine deallocate_internal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module postproc_lfm
