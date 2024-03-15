module wrf_lga

! a module containing necessary items to convert a wrf forecast
! file to lga or to swim for the current laps domain

  use map_utils
  use wrf_netcdf
  use time_utils
  use constants
  use horiz_interp
  use mem_namelist, only: r_missing_data
  use mem_allsky
  implicit none

  private
  integer                :: icentw, jcentw
  integer                :: icentl, jcentl
  real                   :: rmissingflag
  ! laps pressure levels 
  real, allocatable      :: pr_laps(:)

  ! lga variables
  real, allocatable      :: ht(:,:,:)
  real, allocatable      :: t3(:,:,:)
  real, allocatable      :: sh(:,:,:)
  real, allocatable      :: u3(:,:,:)
  real, allocatable      :: v3(:,:,:)
  real, allocatable      :: om(:,:,:)
  ! lgb variables
  real, allocatable      :: usf(:,:)
  real, allocatable      :: vsf(:,:)
  real, allocatable      :: tsf(:,:)
  real, allocatable      :: tsk(:,:) ! surface skin temp. (added by wei-ting 130312)
  real, allocatable      :: dsf(:,:)
  real, allocatable      :: slp(:,:)
  real, allocatable      :: psf(:,:)
  real, allocatable      :: rsf(:,:)
  real, allocatable      :: p(:,:)
  real, allocatable      :: pcp(:,:) ! rainnc+rainc (added by wei-ting 130312)
  ! laps static variables
  real, allocatable      :: topo_laps(:,:)
  real, allocatable      :: lat(:,:)
  real, allocatable      :: lon(:,:)
! integer                :: nxl, nyl, nzl
  character(len=200)     :: laps_data_root
  character(len=10)      :: laps_domain_name  
  real                   :: redp_lvl
  ! wrf on pressure levels
  real, allocatable      :: ht_wrfp(:,:,:)
  real, allocatable      :: t3_wrfp(:,:,:)
  real, allocatable      :: sh_wrfp(:,:,:)
  real, allocatable      :: u3_wrfp(:,:,:)
  real, allocatable      :: v3_wrfp(:,:,:)
  real, allocatable      :: om_wrfp(:,:,:)
  real, allocatable      :: qc_wrfp(:,:,:)
  real, allocatable      :: qi_wrfp(:,:,:)
  real, allocatable      :: qr_wrfp(:,:,:)
  real, allocatable      :: qs_wrfp(:,:,:)
  real, allocatable      :: aod_wrfp(:,:,:) ! extinction coefficient
  real, allocatable      :: usf_wrf(:,:)
  real, allocatable      :: vsf_wrf(:,:)
  real, allocatable      :: tsf_wrf(:,:)
  real, allocatable      :: tsk_wrf(:,:) ! surface skin temp. (added by wei-ting 130312)
  real, allocatable      :: dsf_wrf(:,:)
  real, allocatable      :: slp_wrf(:,:)
  real, allocatable      :: psf_wrf(:,:)
  real, allocatable      :: rsf_wrf(:,:)
  real, allocatable      :: lmk_wrf(:,:)
  real, allocatable      :: snc_wrf(:,:)
  real, allocatable      :: sna_wrf(:,:)
  real, allocatable      :: p_wrf(:,:)
  real, allocatable      :: pcp_wrf(:,:) ! rainnc+rainc (added by wei-ting 130312)
  real, allocatable      :: tvb_wrf(:,:)  ! mean virtual temperature in lowest 60mb
  ! wrf on native variables
  real, allocatable      :: pr_wrfs(:,:,:)
  real, allocatable      :: ht_wrfs(:,:,:)
  real, allocatable      :: dz_wrfs(:,:,:) ! if we want layer thicknesses
  real, allocatable      :: aod_wrfs(:,:,:)
  real, allocatable      :: t3_wrfs(:,:,:)
  real, allocatable      :: sh_wrfs(:,:,:)
  real, allocatable      :: u3_wrfs(:,:,:)
  real, allocatable      :: v3_wrfs(:,:,:)
  real, allocatable      :: om_wrfs(:,:,:)
  real, allocatable      :: rho_wrfs(:,:,:) ! density
  real, allocatable      :: mr_wrfs(:,:,:) ! mixing ratio
  real, allocatable      :: qc_wrfs(:,:,:) ! cloud liquid mixing ratio
  real, allocatable      :: qi_wrfs(:,:,:) ! cloud ice mixing ratio
  real, allocatable      :: qr_wrfs(:,:,:) ! rain  mixing ratio
  real, allocatable      :: qs_wrfs(:,:,:) ! snow mixing ratio
  ! wrf static variables 
  integer                :: cdf,cdp ! added cdp by wei-ting (130312)
  type(proj_info)        :: wrfgrid
  real, allocatable      :: topo_wrf(:,:)
  character(len=19)      :: reftime
  integer                :: tau_hr, tau_min,tau_sec
  integer                :: itimestep,projcode,istat_aod
  integer                :: nxw,nyw,nzw
  real                   :: dx_m, dy_m,dt
  real                   :: lat1_wrf, lon1_wrf
  real                   :: truelat1_wrf, truelat2_wrf, stdlon_wrf

  public wrf2lga, wrf2swim
contains

  subroutine wrf2swim(wrffile_in,i4time,itype_aod,nxl,nyl,nzl,latl,lonl,pres_1d,land_frac,snow_cover,snow_albedo_max,istatus)

     implicit none

     character(len=150)           :: wrffile_in 
     character(len=256)           :: wrffile
     integer, intent(in)          :: i4time
     integer                      :: i4reftime
     character(len=13)            :: reftime13
     integer, intent(out)         :: istatus
     integer                      :: k,k1000,bg_valid,icaller,itype_aod
     integer,external             :: cvt_wfo_fname13_i4time
     integer                      :: nxl,nyl,nzl
     real                         :: ricen,rjcen,latwcen,lonwcen
     real                         :: i_ll, j_ll, i_ul, j_ul, i_ur, j_ur, i_lr, j_lr 
     real                         :: latl(nxl,nyl),lonl(nxl,nyl),land_frac(nxl,nyl),pres_1d(nzl)
     real                         :: snow_cover(nxl,nyl),snow_albedo_max(nxl,nyl)
     logical                      :: need_hinterp
      istatus = 1

      wrffile = wrffile_in

     
     ! get some laps setup stuff
     rmissingflag = r_missing_data

     print *, "subroutine wrf2swim: itype_aod = ",itype_aod

     call find_domain_name(laps_data_root,laps_domain_name,istatus)
!    print *, "laps_data_root = ", trim(laps_data_root)
!    print *, "domain name = ", trim(laps_domain_name)
     print *, "dims:   ", nxl,nyl,nzl
     print *, "rmissingflag:   ", rmissingflag

    ! allocate static fields

     allocate(pr_laps(nzl))
     allocate(lat(nxl,nyl))
     allocate(lon(nxl,nyl))
     allocate(topo_laps(nxl,nyl))
     call get_laps_domain(nxl,nyl,laps_domain_name,lat,lon,topo_laps,istatus)
     if (istatus .ne. 1) then
       print *, "error reading laps static info."
       return
     endif
     pr_laps = pres_1d
     find_k1000: do k = 1,nzl
      if (nint(pres_1d(k)) .eq. 100000) then
        k1000 = k
        exit find_k1000
      endif
     enddo find_k1000

     ! print some config stuff
!    print *, "laps_data_root = ", trim(laps_data_root)
!    print *, "domain name = ", trim(laps_domain_name)
     print *, "dims:   ", nxl,nyl,nzl

     ! get the wrf config
     call open_wrfnc(wrffile,cdf,istatus) 
     call get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr,tau_min,tau_sec,istatus)
     reftime13 = reftime(1:4) // reftime(6:7) // reftime(9:13) // reftime(15:16)
     i4reftime = cvt_wfo_fname13_i4time(reftime13)
     bg_valid = i4reftime + tau_hr * 3600 + tau_min * 60 + tau_sec
     call get_wrf2_map(cdf,'t',projcode,lat1_wrf,lon1_wrf,stdlon_wrf, &
             truelat1_wrf,truelat2_wrf,dx_m,dy_m,nxw,nyw,nzw,istatus)
     call map_set(projcode,lat1_wrf,lon1_wrf,dx_m,stdlon_wrf,truelat1_wrf,truelat2_wrf, &
                  nxw,nyw,wrfgrid)
 
     ! make sure laps domain covers the wrf domain, and see if it is an exact match
     call latlon_to_ij(wrfgrid,lat(1,1),lon(1,1),i_ll,j_ll)
     call latlon_to_ij(wrfgrid,lat(1,nyl),lon(1,nyl),i_ul,j_ul)
     call latlon_to_ij(wrfgrid,lat(nxl,nyl),lon(nxl,nyl),i_ur,j_ur)
     call latlon_to_ij(wrfgrid,lat(nxl,1),lon(nxl,1),i_lr,j_lr)

     ricen = (float(nxw) + 1.) / 2.
     rjcen = (float(nyw) + 1.) / 2.
     call ij_to_latlon(wrfgrid,ricen,rjcen,latwcen,lonwcen)
     print *, "lat / lon of wrf domain center:"
     print *, latwcen,lonwcen

     print *, "location of laps corners in wrf domain (ll,ul,ur,lr):"
     print *, i_ll,j_ll
     print *, i_ul,j_ul
     print *, i_ur,j_ur
     print *, i_lr,j_lr
     if (nint(i_ll) .lt. 1 .or. nint(j_ll) .lt. 1  .or.  &
         nint(i_ul) .lt. 1 .or. nint(j_ul) .gt. nyw .or. &
         nint(i_ur) .gt. nxw .or. nint(j_ur) .gt. nyw .or. &
         nint(i_lr) .gt. nxw .or. nint(j_lr) .lt. 1) then
       print *, "error: laps domain exceeds bounds of wrf background!"
       print *, "consider reducing size of laps domain"
       istatus = 0 
       return
     else
       need_hinterp = .true.
       if (nint(i_ll) .eq. 1 .and. nint(j_ll) .eq. 1 .and. &
           nint(i_ul) .eq. 1 .and. nint(j_ul) .eq. nyw .and. &
           nint(i_ur) .eq. nxw .and. nint(j_ur) .eq. nyw .and. &
           nint(i_lr) .eq. nxw .and. nint(j_lr) .eq. 1 ) then
         print *, "exact match between laps and background.  no hinterp needed!"
         need_hinterp = .false.
       endif
     endif

     icentl = nxl/2
     jcentl = nyl/2
     icentw = nxw/2
     jcentw = nyw/2 
     ! get wrf on sigma
     allocate (pr_wrfs(nxw,nyw,nzw)) 
     allocate (ht_wrfs(nxw,nyw,nzw))
     allocate (dz_wrfs(nxw,nyw,nzw))
     allocate (aod_wrfs(nxw,nyw,nzw))
     allocate (t3_wrfs(nxw,nyw,nzw))
     allocate (sh_wrfs(nxw,nyw,nzw)) 
     allocate (mr_wrfs(nxw,nyw,nzw))
     allocate (rho_wrfs(nxw,nyw,nzw))
     allocate (qc_wrfs(nxw,nyw,nzw))
     allocate (qi_wrfs(nxw,nyw,nzw))
     allocate (qr_wrfs(nxw,nyw,nzw))
     allocate (qs_wrfs(nxw,nyw,nzw))
     allocate (usf_wrf(nxw,nyw))
     allocate (vsf_wrf(nxw,nyw))
     allocate (tsf_wrf(nxw,nyw))
     allocate (tsk_wrf(nxw,nyw)) ! surface skin temp. (added by wei-ting 130312)
     allocate (rsf_wrf(nxw,nyw))
     allocate (lmk_wrf(nxw,nyw))
     allocate (snc_wrf(nxw,nyw))
     allocate (sna_wrf(nxw,nyw))
     allocate (dsf_wrf(nxw,nyw))
     allocate (slp_wrf(nxw,nyw))
     allocate (psf_wrf(nxw,nyw))
     allocate (p_wrf(nxw,nyw))
     allocate (pcp_wrf(nxw,nyw)) ! rainnc+rainc (added by wei-ting 130312)
     allocate (topo_wrf(nxw,nyw))
     allocate (tvb_wrf(nxw,nyw))

     icaller = 1
!    itype_aod = 2 ! 1 is taod5503d, 2 is tracer_1a
     call fill_wrfs(icaller,itype_aod,istatus)
     if (istatus .ne. 1) then
       print *, "problem getting wrf data"
       return
     endif

     ! vertically interpolate to pressure levels
     print *, "allocating arrays for wrf on pressure"
       ! allocate wrfp
     allocate (ht_wrfp(nxw,nyw,nzl))
     allocate (t3_wrfp(nxw,nyw,nzl))
     allocate (sh_wrfp(nxw,nyw,nzl))
     allocate (qc_wrfp(nxw,nyw,nzl))
     allocate (qi_wrfp(nxw,nyw,nzl))
     allocate (qr_wrfp(nxw,nyw,nzl))
     allocate (qs_wrfp(nxw,nyw,nzl))
     allocate (aod_wrfp(nxw,nyw,nzl))
     ht_wrfp = rmissingflag
     t3_wrfp = rmissingflag
     sh_wrfp = rmissingflag
     qc_wrfp = rmissingflag
     qi_wrfp = rmissingflag
     qr_wrfp = rmissingflag
     qs_wrfp = rmissingflag
     aod_wrfp = rmissingflag
  
     ! vertically interpolate
     print *, "calling vinterp_wrfarw2p ",nzl,pr_laps(1)
     call vinterp_wrfarw2p(icaller,nzl,istatus)
     if (istatus .ne. 1) then
       print *, "problem vertically interpolating wrf data"
       return
     endif

     ! dealloc wrfs
     print *, "deallocating wrf sigma var"
     deallocate(pr_wrfs,ht_wrfs,dz_wrfs,aod_wrfs,sh_wrfs,t3_wrfs,mr_wrfs, &  
        rho_wrfs,qc_wrfs,qi_wrfs,qr_wrfs,qs_wrfs)

     ! horizontally interpolate to laps grid
       ! allocate lga/lgb

     print *, "allocating lga variables"
     allocate(ht(nxl,nyl,nzl))
     allocate(sh(nxl,nyl,nzl))
     allocate(t3(nxl,nyl,nzl))
     ! allocate lgb (2d) variables
     print *, "allocating lgb variables"
     allocate (usf(nxl,nyl))
     allocate (vsf(nxl,nyl))
     allocate (tsf(nxl,nyl))
     allocate (tsk(nxl,nyl)) ! surface skin temp. (added by wei-ting 130312)
     allocate (rsf(nxl,nyl))
     allocate (dsf(nxl,nyl))
     allocate (slp(nxl,nyl))
     allocate (psf(nxl,nyl))
     allocate (p  (nxl,nyl))
     allocate (pcp(nxl,nyl)) ! precitation (added by wei-ting 130312)
     usf = rmissingflag
     vsf = rmissingflag
     tsf = rmissingflag
     tsk = rmissingflag ! surface skin temp. (added by wei-ting 130312)
     dsf = rmissingflag
     slp = rmissingflag
     psf = rmissingflag
     p   = rmissingflag
     pcp = rmissingflag ! rainnc+rainc (added by wei-ting 130312)

     if (need_hinterp) then
       print *, "problem in wrf2swim: hinterp indicated as needed" 
       istatus = 0
       return
     else
       ht = ht_wrfp
       t3 = t3_wrfp
       sh = sh_wrfp
       clwc_3d = qc_wrfp
       cice_3d = qi_wrfp
       rain_3d = qr_wrfp
       snow_3d = qs_wrfp
       if(istat_aod .eq. 1 .and. mode_aero_cld .eq. 3)then
         write(6,*)' transferring aod-3d to all-sky array'
         aod_3d(:,:,:) = aod_wrfp(:,:,:)
       endif
       psf = psf_wrf
       tsf = tsf_wrf
       tsk = tsk_wrf ! surface skin temp. (added by wei-ting 130312)
       dsf = dsf_wrf
       rsf = rsf_wrf
       usf = usf_wrf
       vsf = vsf_wrf
       pcp = pcp_wrf ! rainnc+rainc (added by wei-ting 130312)
       land_frac = lmk_wrf ! land mask
       snow_cover = snc_wrf
       snow_albedo_max = sna_wrf
     endif 
!     pcp = 0 ! since pcp hasn't be used for now, assume that the value is 0

     print *, "deallocating wrf press vars"
     deallocate (ht_wrfp,t3_wrfp,sh_wrfp,qc_wrfp,qi_wrfp,qr_wrfp,qs_wrfp,aod_wrfp)

     print *, "deallocating wrf sfc vars"
     deallocate (usf_wrf,vsf_wrf,tsf_wrf,tsk_wrf,rsf_wrf,lmk_wrf,snc_wrf,sna_wrf,dsf_wrf,slp_wrf,&
         psf_wrf,p_wrf,pcp_wrf,tvb_wrf) ! added tsk_wrf & pcp_wrf by wei-ting (130312)
 
     ! create mslp and reduced pressure
     print *, "topo/slp/psf/p/tsf/tsk/dsf/rsf/usf/vsf",topo_laps(icentl,jcentl), &
       slp(icentl,jcentl),psf(icentl,jcentl),p(icentl,jcentl),tsf(icentl,jcentl), &
       tsk(icentl,jcentl),dsf(icentl,jcentl),rsf(icentl,jcentl),usf(icentl,jcentl), &
       vsf(icentl,jcentl) ! added tsk by wei-ting (130312)



     print *, "deallocating laps static vars"
     deallocate (lat,lon,topo_laps,topo_wrf)

     ! deallocate memory
     print *, "deallocating lga vars"
     deallocate(ht,t3,sh,pr_laps)
     print *, "deallocating lgb vars"
     deallocate(usf,vsf,tsf,tsk,rsf,dsf,slp,psf,p,pcp) ! added tsk & pcp by wei-ting (130312)
     print *, "successful processing of ", trim(wrffile) ! modified by wei-ting
     print *, "return from wrf2swim"
     return 
  end subroutine wrf2swim

  subroutine wrf2lga(wrffile,i4time,cmodel,istatus )

     implicit none

     character(len=256), intent(in) :: wrffile(2) ! add 1 dim. and (len=255 -> len=256) by wei-ting (130312) to contain previous time
     character(len=12), intent(in) :: cmodel
     integer, intent(in)          :: i4time
     integer                      :: i4reftime
     integer                      :: nxl,nyl,nzl
     character(len=13)            :: reftime13
     integer, intent(out)         :: istatus
     integer                      :: k,k1000,bg_valid,icaller,itype_aod
     integer,external             :: cvt_wfo_fname13_i4time
     real                         :: i_ll, j_ll, i_ul, j_ul, i_ur, j_ur, i_lr, j_lr 
     logical                      :: need_hinterp
      istatus = 1

     
     ! get some laps setup stuff
     call find_domain_name(laps_data_root,laps_domain_name,istatus)
     print *, "laps_data_root = ", trim(laps_data_root)
     print *, "domain name = ", trim(laps_domain_name)
    
     call get_grid_dim_xy(nxl,nyl,istatus)
     if (istatus .ne. 1) then
       print *, "could not get laps xy dims"
       return
     endif
     call get_laps_dimensions(nzl,istatus)
     if (istatus .ne. 1) then
       print *, "could not get laps z dim"
       return
     endif
     print *, "dims:   ", nxl,nyl,nzl
     call get_r_missing_data(rmissingflag,istatus) 
     call get_laps_redp(redp_lvl,istatus)
    ! allocate static fields

     allocate(pr_laps(nzl))
     allocate(lat(nxl,nyl))
     allocate(lon(nxl,nyl))
     allocate(topo_laps(nxl,nyl))
     call get_laps_domain(nxl,nyl,laps_domain_name,lat,lon,topo_laps,istatus)
     if (istatus .ne. 1) then
       print *, "error reading laps static info."
       return
     endif
     call get_pres_1d(i4time,nzl,pr_laps,istatus)
     if (istatus .ne. 1) then
       print *, "error reading laps pressure levels"
       return
     endif
     find_k1000: do k = 1,nzl
      if (nint(pr_laps(k)) .eq. 100000) then
        k1000 = k
        exit find_k1000
      endif
     enddo find_k1000
     ! print some config stuff
     print *, "laps_data_root = ", trim(laps_data_root)
     print *, "domain name = ", trim(laps_domain_name)
     print *, "dims:   ", nxl,nyl,nzl


     ! get the wrf config
     call open_wrfnc(wrffile(1),cdf,istatus) ! modified by wei-ting to use right time ( wrffile -> wrffile(1) )
     call open_wrfnc(wrffile(2),cdp,istatus) ! added by wei-ting to use previous time (just for calculate pcp)
     call get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr,tau_min,tau_sec,istatus)
     reftime13 = reftime(1:4) // reftime(6:7) // reftime(9:13) // reftime(15:16)
     i4reftime = cvt_wfo_fname13_i4time(reftime13)
     bg_valid = i4reftime + tau_hr * 3600 + tau_min * 60 + tau_sec
     call get_wrf2_map(cdf,'t',projcode,lat1_wrf,lon1_wrf,stdlon_wrf, &
             truelat1_wrf,truelat2_wrf,dx_m,dy_m,nxw,nyw,nzw,istatus)
     call map_set(projcode,lat1_wrf,lon1_wrf,dx_m,stdlon_wrf,truelat1_wrf,truelat2_wrf, &
                  nxw,nyw,wrfgrid)
 
     ! make sure laps domain covers the wrf domain, and see if it is an exact match
     call latlon_to_ij(wrfgrid,lat(1,1),lon(1,1),i_ll,j_ll)
     call latlon_to_ij(wrfgrid,lat(1,nyl),lon(1,nyl),i_ul,j_ul)
     call latlon_to_ij(wrfgrid,lat(nxl,nyl),lon(nxl,nyl),i_ur,j_ur)
     call latlon_to_ij(wrfgrid,lat(nxl,1),lon(nxl,1),i_lr,j_lr)
     print *, "location of laps corners in wrf domain (ll,ul,ur,lr):"
     print *, i_ll,j_ll
     print *, i_ul,j_ul
     print *, i_ur,j_ur
     print *, i_lr,j_lr
     if (nint(i_ll) .lt. 1 .or. nint(j_ll) .lt. 1  .or.  &
         nint(i_ul) .lt. 1 .or. nint(j_ul) .gt. nyw .or. &
         nint(i_ur) .gt. nxw .or. nint(j_ur) .gt. nyw .or. &
         nint(i_lr) .gt. nxw .or. nint(j_lr) .lt. 1) then
       print *, "laps domain exceeds bounds of wrf background!"
       print *, "consider reducing size of laps domain"
       istatus = 0 
       return
     else
       need_hinterp = .true.
       if (nint(i_ll) .eq. 1 .and. nint(j_ll) .eq. 1 .and. &
           nint(i_ul) .eq. 1 .and. nint(j_ul) .eq. nyw .and. &
           nint(i_ur) .eq. nxw .and. nint(j_ur) .eq. nyw .and. &
           nint(i_lr) .eq. nxw .and. nint(j_lr) .eq. 1 ) then
         print *, "exact match between laps and background.  no hinterp needed!"
         need_hinterp = .false.
       endif
     endif

     icentl = nxl/2
     jcentl = nyl/2
     icentw = nxw/2
     jcentw = nyw/2 
     ! get wrf on sigma
     allocate (pr_wrfs(nxw,nyw,nzw)) 
     allocate (ht_wrfs(nxw,nyw,nzw))
     allocate (t3_wrfs(nxw,nyw,nzw))
     allocate (sh_wrfs(nxw,nyw,nzw)) 
     allocate (u3_wrfs(nxw,nyw,nzw))
     allocate (v3_wrfs(nxw,nyw,nzw))
     allocate (om_wrfs(nxw,nyw,nzw)) 
     allocate (mr_wrfs(nxw,nyw,nzw))
     allocate (rho_wrfs(nxw,nyw,nzw))
     allocate (usf_wrf(nxw,nyw))
     allocate (vsf_wrf(nxw,nyw))
     allocate (tsf_wrf(nxw,nyw))
     allocate (tsk_wrf(nxw,nyw)) ! surface skin temp. (added by wei-ting 130312)
     allocate (rsf_wrf(nxw,nyw))
     allocate (dsf_wrf(nxw,nyw))
     allocate (slp_wrf(nxw,nyw))
     allocate (psf_wrf(nxw,nyw))
     allocate (p_wrf(nxw,nyw))
     allocate (pcp_wrf(nxw,nyw)) ! rainnc+rainc (added by wei-ting 130312)
     allocate (topo_wrf(nxw,nyw))
     allocate (tvb_wrf(nxw,nyw))

     icaller = 2
     itype_aod = 2 ! 1 is taod5503d, 2 is tracer_1a (probably a dummy value)
     call fill_wrfs(icaller,itype_aod,istatus)
     if (istatus .ne. 1) then
       print *, "problem getting wrf data"
       return
     endif

     ! vertically interpolate to pressure levels
     print *, "allocating arrays for wrf on pressure"
       ! allocate wrfp
     allocate (ht_wrfp(nxw,nyw,nzl))
     allocate (t3_wrfp(nxw,nyw,nzl))
     allocate (sh_wrfp(nxw,nyw,nzl))
     allocate (u3_wrfp(nxw,nyw,nzl))
     allocate (v3_wrfp(nxw,nyw,nzl))
     allocate (om_wrfp(nxw,nyw,nzl))
     ht_wrfp = rmissingflag
     t3_wrfp = rmissingflag
     u3_wrfp = rmissingflag
     v3_wrfp = rmissingflag
     om_wrfp = rmissingflag
  
     ! vertically interpolate
     print *, "calling vinterp_wrfarw2p ",nzl,pr_laps(1)
     call vinterp_wrfarw2p(icaller,nzl,istatus)
     if (istatus .ne. 1) then
       print *, "problem vertically interpolating wrf data"
       return
     endif

     ! dealloc wrfs
     print *, "deallocating wrf sigma var"
     deallocate(pr_wrfs,ht_wrfs,sh_wrfs,t3_wrfs,u3_wrfs,v3_wrfs,om_wrfs,mr_wrfs, &
        rho_wrfs)

     ! horizontally interpolate to laps grid
       ! allocate lga/lgb

     print *, "allocating lga variables"
     allocate(ht(nxl,nyl,nzl))
     allocate(sh(nxl,nyl,nzl))
     allocate(t3(nxl,nyl,nzl))
     allocate(u3(nxl,nyl,nzl))
     allocate(v3(nxl,nyl,nzl))
     allocate(om(nxl,nyl,nzl))
     ! allocate lgb (2d) variables
     print *, "allocating lgb variables"
     allocate (usf(nxl,nyl))
     allocate (vsf(nxl,nyl))
     allocate (tsf(nxl,nyl))
     allocate (tsk(nxl,nyl)) ! surface skin temp. (added by wei-ting 130312)
     allocate (rsf(nxl,nyl))
     allocate (dsf(nxl,nyl))
     allocate (slp(nxl,nyl))
     allocate (psf(nxl,nyl))
     allocate (p  (nxl,nyl))
     allocate (pcp(nxl,nyl)) ! precitation (added by wei-ting 130312)
     usf = rmissingflag
     vsf = rmissingflag
     tsf = rmissingflag
     tsk = rmissingflag ! surface skin temp. (added by wei-ting 130312)
     dsf = rmissingflag
     slp = rmissingflag
     psf = rmissingflag
     p   = rmissingflag
     pcp = rmissingflag ! rainnc+rainc (added by wei-ting 130312)

     if (need_hinterp) then
       ht = rmissingflag
       sh = rmissingflag
       t3 = rmissingflag
       u3 = rmissingflag
       v3 = rmissingflag
       om = rmissingflag
       usf = rmissingflag
       vsf = rmissingflag
       tsf = rmissingflag
       tsk = rmissingflag ! surface skin temp. (added by wei-ting 130312)
       dsf = rmissingflag
       slp = rmissingflag
       psf = rmissingflag
       p   = rmissingflag
       pcp = rmissingflag ! rainnc+rainc (added by wei-ting 130312)
                                                                                                                            

       ! call hinterplga
       print *, "calling hinterp_wrf2lga"
       call hinterp_wrf2lga(nxl,nyl,nzl,istatus)
       if (istatus .ne. 1) then 
         print *, "problem in hinterp_wrf2lga" 
         istatus = 0
         return
       endif 
     else
       ht = ht_wrfp
       t3 = t3_wrfp
       sh = sh_wrfp
       u3 = u3_wrfp
       v3 = v3_wrfp
       om = om_wrfp
       psf = psf_wrf
       tsf = tsf_wrf
       tsk = tsk_wrf ! surface skin temp. (added by wei-ting 130312)
       dsf = dsf_wrf
       rsf = rsf_wrf
       usf = usf_wrf
       vsf = vsf_wrf
       pcp = pcp_wrf ! rainnc+rainc (added by wei-ting 130312)
     endif 
!     pcp = 0 ! since pcp hasn't be used for now, assume that the value is 0

     print *, "deallocating wrf press vars"
     deallocate (ht_wrfp,t3_wrfp,sh_wrfp,u3_wrfp,v3_wrfp,om_wrfp)

     print *, "deallocating wrf sfc vars"
     deallocate (usf_wrf,vsf_wrf,tsf_wrf,tsk_wrf,rsf_wrf,dsf_wrf,slp_wrf,&
         psf_wrf,p_wrf,pcp_wrf,tvb_wrf) ! added tsk_wrf & pcp_wrf by wei-ting (130312)
 
     ! create mslp and reduced pressure
     print *, "creating reduced pressure arrays"
     call make_derived_pressures(nxl,nyl,nzl)
     print *, "topo/slp/psf/p/tsf/tsk/dsf/rsf/usf/vsf",topo_laps(icentl,jcentl), &
       slp(icentl,jcentl),psf(icentl,jcentl),p(icentl,jcentl),tsf(icentl,jcentl), &
       tsk(icentl,jcentl),dsf(icentl,jcentl),rsf(icentl,jcentl),usf(icentl,jcentl), &
       vsf(icentl,jcentl) ! added tsk by wei-ting (130312)



     print *, "deallocating laps static vars"
     deallocate (lat,lon,topo_laps,topo_wrf)

     ! write the data out
     print *, "writing lga"
     call write_lga(nxl,nyl,nzl,i4reftime,bg_valid,cmodel,rmissingflag, &
           pr_laps*0.01,ht,t3,sh,u3,v3,om,istatus)
     if(istatus .ne. 1) then
       print *, "error writing lga"
       return
     endif

     print *, "writing lgb" ! modified by wei-ting (lga -> lgb)
     call write_lgb(nxl,nyl,i4reftime,bg_valid,cmodel,rmissingflag, &
           usf,vsf,tsf,tsk,rsf,psf,slp,dsf,p,pcp,istatus) ! added tsk & pcp by wei-ting (130312)
     if(istatus .ne. 1) then
       print *, "error writing lgb"
       return
     endif


     ! deallocate memory
     print *, "deallocating lga vars"
     deallocate(ht,t3,sh,u3,v3,om,pr_laps)
     print *, "deallocating lgb vars"
     deallocate(usf,vsf,tsf,tsk,rsf,dsf,slp,psf,p,pcp) ! added tsk & pcp by wei-ting (130312)
     print *, "successful processing of ", trim(wrffile(1)) ! modified by wei-ting
     return 
  end subroutine wrf2lga
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fill_wrfs(icaller,itype_aod,istatus)
 
    implicit none
    integer, intent(out)  :: istatus
    integer               :: status,i,j,k,icaller,itype_aod
    real, allocatable     :: dum3d(:,:,:)
    real, allocatable     :: dum3df(:,:,:)
    real, allocatable     :: dum3df2(:,:,:)
    real, allocatable     :: dum2dt1(:,:) ! added by wei-ting (130312) to get rainnc
    real, allocatable     :: dum2dt2(:,:) ! added by wei-ting (130312) to get rainnc
    real, external        :: mixsat, relhum, dewpt2
    real                  :: rh, tr_conv 
    real                  :: tvbar, tvbar_nlevs
    real, parameter       :: tvbar_thick = 6000.
    ! varialbles have already been allocated by our driver routine
    ! so just start getting them
    istatus = 1
    print *, " allocating arrays ",icaller
    allocate(dum3d(nxw,nyw,nzw))
    allocate(dum3df(nxw,nyw,nzw+1))
    allocate(dum3df2(nxw,nyw,nzw+1))
    allocate(dum2dt1(nxw,nyw))
    allocate(dum2dt2(nxw,nyw))

    ! get 3d pressure array
    ! get pressures
    print *, "getting pb"
    call get_wrfnc_3d(cdf,"pb","t",nxw,nyw,nzw,1,dum3d,status)
    if (status .gt. 0) then
      print *, "+++critical:  could not get base pressure!"
      istatus = 0
      return
    endif
    pr_wrfs = dum3d
                        
    print *, "getting p"                                                                    
    call get_wrfnc_3d(cdf,"p","t",nxw,nyw,nzw,1,dum3d,status)
    if (status .gt. 0) then
      print *, "+++critical:  could not get pert pressure!"
      istatus = 0
      return
    endif
    pr_wrfs = pr_wrfs + dum3d
    print *, "min/max wrf 3d pressure: ",minval(pr_wrfs),maxval(pr_wrfs)
  
    ! get heights
    print *, "getting phb" 
    call get_wrfnc_3d(cdf,"phb","t",nxw,nyw,nzw+1,1,dum3df,status)
    if (status.ne.0) then
      print *, 'could not properly obtain wrf base-state geopotential.'
      istatus = 0
      return
    endif
    dum3df2 = dum3df
  
    print *, "getting ph"
    call get_wrfnc_3d(cdf,"ph","t",nxw,nyw,nzw+1,1,dum3df,status)
    if (status.ne.0) then
      print *, 'could not properly obtain wrf geopotential.'
      istatus = 0
      return
    endif
    print *,"destaggering (vertically) heights"
    dum3df2 = (dum3df2 + dum3df) / grav
    do k = 1,nzw
      ht_wrfs(:,:,k) = 0.5 * (dum3df2(:,:,k) + dum3df2(:,:,k+1))
      if(icaller .eq. 1)then
        dz_wrfs(:,:,k) = dum3df2(:,:,k+1) - dum3df2(:,:,k)
      endif
    enddo
  
    print *, "getting aod"
    if(itype_aod .eq. 1)then
      call get_wrfnc_3d(cdf,"taod5503d","t",nxw,nyw,nzw,1,dum3df,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf taod5503d - set to missing'
        aod_wrfs(:,:,:) = r_missing_data
        istat_aod = 0
      else
        print *, 'success reading wrf taod5503d - divide by dz'
        aod_wrfs(:,:,:) = dum3df(:,:,:) / dz_wrfs(:,:,:)
        print *, "min/max wrf 3d aod: ",minval(aod_wrfs),maxval(aod_wrfs)
        istat_aod = 1
      endif
    elseif(itype_aod .eq. 2)then
      call get_wrfnc_3d(cdf,"tracer_1a","t",nxw,nyw,nzw,1,dum3df,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf tracer_1a - set to missing'
        aod_wrfs(:,:,:) = r_missing_data
        istat_aod = 0
      else
        print *, 'success reading wrf tracer_1a - ,multiply by tr_conv'
        tr_conv = 1e-9             ! ug/kg to kg/kg (dimensionless)
        tr_conv = tr_conv * 3.0    ! smoke mee m^2/g
        tr_conv = tr_conv * 1000.0 ! smoke mee m^2/kg
        print *, "tracer_1a conversion factor (sans rho) is ",tr_conv
        aod_wrfs(:,:,:) = dum3df(:,:,:) * tr_conv
        print *, "min/max wrf 3d aod: ",minval(aod_wrfs),maxval(aod_wrfs)
        istat_aod = 1
      endif
    elseif(itype_aod .eq. 3)then
      print *, 'reading 3-d aerosols from total_ext_file'
      open(65,file='/users/albers/data/projects/muri/wrfchem_raw/total_ext_file',status='old',form='unformatted')
!     read(65)dum3df
      read(65,err=91)dum3df
      goto 92
      close(65)
91    print *, 'read error - setting dum3df to 0.'
      dum3df = 0.
92    print *, 'range of dum3df is',minval(dum3df),maxval(dum3df)
    else
      print *, 'option for no wrf aerosols to be read in'
    endif

    ! get theta and convert to temperature
    print *, "getting theta"
    call get_wrfnc_3d(cdf, "t","t",nxw,nyw,nzw,1,dum3d,status)
    if (status.ne.0) then
      print *, 'could not properly obtain wrf perturbation theta.'
      istatus = 0
      return
    endif

    print *, "computing temp"
    dum3d = dum3d + 300.
    do k = 1, nzw
      do j = 1, nyw
        do i = 1,nxw
          t3_wrfs(i,j,k) = dum3d(i,j,k)/ ((100000./pr_wrfs(i,j,k))**kappa)
        enddo
      enddo
    enddo

    ! get q on sigma
    print *, "getting q"
    call get_wrfnc_3d(cdf, "qvapor","t",nxw,nyw,nzw,1,mr_wrfs,status)
    if (status.ne.0) then
      print *, 'could not properly obtain wrf mixing ratio.'
      istatus = 0
      return
    endif

    ! derive specific humidity
    print *, "computing sh"
    sh_wrfs = mr_wrfs / (1. + mr_wrfs)

    if(icaller .eq. 1)then
      print *, "getting qc"
      call get_wrfnc_3d(cdf, "qcloud","t",nxw,nyw,nzw,1,qc_wrfs,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf qc mixing ratio.'
        istatus = 0
        return
      endif

      print *, "getting qi"
      call get_wrfnc_3d(cdf, "qice","t",nxw,nyw,nzw,1,qi_wrfs,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf qi mixing ratio.'
        istatus = 0
        return
      endif

      print *, "getting qr"
      call get_wrfnc_3d(cdf, "qrain","t",nxw,nyw,nzw,1,qr_wrfs,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf qr mixing ratio.'
        istatus = 0
        return
      endif

      print *, "getting qs"
      call get_wrfnc_3d(cdf, "qsnow","t",nxw,nyw,nzw,1,qs_wrfs,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf qs mixing ratio.'
        istatus = 0
        return
      endif

    endif

    if(icaller .eq. 2)then
    ! get u on sigma   
      print *, "getting u"
      call get_wrfnc_3d(cdf, "u","t",nxw,nyw,nzw,1,u3_wrfs,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf u-comp.'
        istatus = 0
        return
      endif

    ! get v on sigma
      print *, "getting v"
      call get_wrfnc_3d(cdf, "v","t",nxw,nyw,nzw,1,v3_wrfs,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf v-comp.'
        istatus = 0
        return
      endif

    ! get w on sigma
      print *, "getting w"
      call get_wrfnc_3d(cdf, "w","t",nxw,nyw,nzw+1,1,dum3df,status)
      if (status.ne.0) then
        print *, 'could not properly obtain wrf w-comp'
        istatus = 0
        return 
      endif
      print*, "destaggering (vertically) w"
      do k = 1,nzw
        om_wrfs(:,:,k) = 0.5*(dum3df(:,:,k)+dum3df(:,:,k+1))
      enddo
    endif 

    ! now, derive density, virtual potential temp, and omega
    print *, "computing rho, theta-v, and omega"
    do k = 1,nzw
      do j=1,nyw
        do i = 1,nxw
          rho_wrfs(i,j,k) = pr_wrfs(i,j,k) / ( r * t3_wrfs(i,j,k)*(1.+0.61*sh_wrfs(i,j,k)))
          if(icaller .eq. 2)then
            om_wrfs(i,j,k) = -1. * rho_wrfs(i,j,k) * grav * om_wrfs(i,j,k)
          endif
        enddo
      enddo
    enddo

    if(icaller .eq. 1)then ! convert mixing ratio to density
      print *, "convert qc,qi,qr,qs to density"
      qc_wrfs(:,:,:) = qc_wrfs(:,:,:) * rho_wrfs(:,:,:)
      qi_wrfs(:,:,:) = qi_wrfs(:,:,:) * rho_wrfs(:,:,:)
      qr_wrfs(:,:,:) = qr_wrfs(:,:,:) * rho_wrfs(:,:,:)
      qs_wrfs(:,:,:) = qs_wrfs(:,:,:) * rho_wrfs(:,:,:)
      if(itype_aod .eq. 2)then
        print *, "apply density conversion to aod to complete conversion to extinction coefficient (m^-1)"
        print *, "representative value of density is ",rho_wrfs(1,1,1)
        aod_wrfs(:,:,:) = aod_wrfs(:,:,:) * rho_wrfs(:,:,:)
        print *, "revised min/max wrf 3d aod: ",minval(aod_wrfs),maxval(aod_wrfs)
      endif
    endif
   
    ! get surface fields
    
    print *, "getting wrf topo"
    call get_wrfnc_2d(cdf, "hgt","t",nxw,nyw,1,topo_wrf,status)
    if (status .ne. 0) then
      print *, "could not get topo from wrf"
      istatus = 0 
      return
    endif
    
    print *, "getting wrf landmask"
    call get_wrfnc_2d(cdf, "landmask","t",nxw,nyw,1,lmk_wrf,status)
    if (status .ne. 0) then
      print *, "could not get land mask from wrf"
      istatus = 0 
      return
    endif
    
    print *, "getting wrf snowc"
    call get_wrfnc_2d(cdf, "snowc","t",nxw,nyw,1,snc_wrf,status)
    if (status .ne. 0) then
      print *, "could not get snow cover from wrf"
      istatus = 0 
      return
    endif
    
    print *, "getting wrf snoalb"
    call get_wrfnc_2d(cdf, "snoalb","t",nxw,nyw,1,sna_wrf,status)
    if (status .ne. 0) then
      print *, "could not get snow albedo from wrf"
      istatus = 0 
      return
    endif

    if(icaller .eq. 2)then
      print *, "getting usf"
      usf_wrf = u3_wrfs(:,:,1)
      print *, "getting vsf"
      vsf_wrf = v3_wrfs(:,:,1)
    endif

    print *, "getting psf"
    call get_wrfnc_2d(cdf, "psfc","t",nxw,nyw,1,psf_wrf,status)
    if ((status .ne. 0).or.(maxval(psf_wrf) .lt. 10000.))then
      print *, "could not get psfc, using lowest sigma level"
      psf_wrf = pr_wrfs(:,:,1)
    endif
 
    print *, "getting t2" 
    call get_wrfnc_2d(cdf, "t2","t",nxw,nyw,1,tsf_wrf,status)
    if ((status .ne. 0).or.(maxval(tsf_wrf) .lt. 100.))then
      print *, "could not get t2, using lowest sigma level"
      tsf_wrf = t3_wrfs(:,:,1)
    endif
    
    ! added tsk(skin temp.) by wei-ting (130312)
    print *, "getting tsk" 
    call get_wrfnc_2d(cdf, "tsk","t",nxw,nyw,1,tsk_wrf,status)
    if ((status .ne. 0).or.(maxval(tsk_wrf) .lt. 100.))then
      print *, "could not get tsk, using lowest sigma level"
      tsk_wrf = t3_wrfs(:,:,1)
    endif
    
    ! added pcp (rainnc+rainc) by wei-ting (130312) & modified (130326)
    print *, "getting precipitation"
    print *, "!!!!! this precipitaion is an accumulation per n hours. !!!!!"
    print *, "!!!!! n depends on the time difference of each wrfout.  !!!!!"
    print *, "   getting rainnc(t)"
    call get_wrfnc_2d(cdf,"rainnc","t",nxw,nyw,1,dum2dt2,status)
    if (status .ne. 0) then
      print *, "   could not get rainnc(t), setting the value = 0"
      dum2dt2 = 0
    endif
    pcp_wrf = dum2dt2
    print *, "   getting rainc(t)"
    call get_wrfnc_2d(cdf,"rainc","t",nxw,nyw,1,dum2dt2,status)
    if (status .ne. 0) then
      print *, "   could not get rainc(t), setting the value = 0"
      dum2dt2 = 0
    endif
    pcp_wrf = pcp_wrf+dum2dt2

    print *, "   getting rainnc(t-1)"
    if (cdp .lt. 0 .and. cdf .gt. 0) then
      print *, "   could not get rainnc(t-1), maybe result from t = initial time!"
      print *, "   set rainnc(t-1) = 0"
      dum2dt1 = 0
    else
      call get_wrfnc_2d(cdp,"rainnc","t",nxw,nyw,1,dum2dt1,status)
      if (status .ne. 0) then
         print *, "   could not get rainnc(t-1), setting the value = 0"
         dum2dt1 = 0
      endif
    endif
    pcp_wrf = pcp_wrf-dum2dt1
    print *, "   getting rainc(t-1)"
    if (cdp .lt. 0 .and. cdf .gt. 0) then
      print *, "   could not get rainc(t-1), maybe result from t = initial time!"
      print *, "   set rainc(t-1) = 0"
      dum2dt1 = 0
    else
      call get_wrfnc_2d(cdp,"rainc","t",nxw,nyw,1,dum2dt1,status)
      if (status .ne. 0) then
         print *, "   could not get rainc(t-1), setting the value = 0"
         dum2dt1 = 0
      endif
    endif
    pcp_wrf = pcp_wrf-dum2dt1
    where ( pcp_wrf < 0 ) ; pcp_wrf = 0 ; endwhere ! keep pcp >= 0
    print *, "min/max wrf precipitation : ",minval(pcp_wrf),maxval(pcp_wrf)
    ! end of reading rainnc+rainc

    ! qvapor at 2m
    print *, "getting q2"
    call get_wrfnc_2d(cdf, "q2","t",nxw,nyw,1,rsf_wrf,status)
    if ((status .ne. 0).or.(maxval(rsf_wrf) .lt. 0.0001))then
      print *, "could not get q2, using lowest sigma level"
      rsf_wrf = sh_wrfs(:,:,1)
    else 

      ! because 2m qv and t are derived from the pbl scheme and
      ! the wrf is apparently not checking for saturation, clean this
      ! up now
      print *, "checking q2 for supersaturation"
      do j = 1, nyw
        do i= 1, nxw
          rsf_wrf(i,j) = min(rsf_wrf(i,j),mixsat(rsf_wrf(i,j),psf_wrf(i,j)))
          ! compute dewpoint
          rh = relhum(tsf_wrf(i,j),rsf_wrf(i,j),psf_wrf(i,j))
          dsf_wrf(i,j) = dewpt2(tsf_wrf(i,j),rh)
          ! compute tvbar
          tvbar = 0.
          tvbar_nlevs = 0.
          comptvb:  do k = 1,nzw
            if ((psf_wrf(i,j)-pr_wrfs(i,j,k)) .le. tvbar_thick) then
               tvbar = tvbar + (t3_wrfs(i,j,k)*(1.+0.61*sh_wrfs(i,j,k)) )
               tvbar_nlevs = tvbar_nlevs + 1.
            else
              exit comptvb
            endif
          enddo comptvb
          tvb_wrf(i,j) = tvbar / tvbar_nlevs
        enddo
      enddo
    endif 
    ! smooth the tvb_wrf field
    call smooth2(nxw,nyw,4,tvb_wrf)
    ! convert mr to sh
    rsf_wrf(:,:) = rsf_wrf(:,:)/(1. + rsf_wrf(:,:))

    ! diagnostics
    print *, "wrf sigma data from center of wrf domain"
    if(icaller .eq. 1)then ! called from 'wrf2swim'
      print *, "k   press     height   dz     temp   sh         qc        qi       qr       qs      aod"
      print *, "--- --------  -------  -----  -----  -------    -------   ------   -----    -----   ------"
    else                   ! called from 'wrf2lga'
      print *, "k   press     height   temp   sh       u       v      om"
      print *, "--- --------  -------  -----  -------  ------  ------ -----------"
    endif
    do k = 1,nzw
      if(icaller .eq. 1)then ! called from 'wrf2swim'
        print ('(i3,1x,f8.1,2x,f7.0,2x,f6.1,2x,f5.1,2x,f7.5,2x,4f9.5,2x,f8.6)'), &
          k,pr_wrfs(icentw,jcentw,k),ht_wrfs(icentw,jcentw,k),dz_wrfs(icentw,jcentw,k), t3_wrfs(icentw,jcentw,k), &
          sh_wrfs(icentw,jcentw,k),qc_wrfs(icentw,jcentw,k),qi_wrfs(icentw,jcentw,k), &
          qr_wrfs(icentw,jcentw,k),qs_wrfs(icentw,jcentw,k),aod_wrfs(icentw,jcentw,k)
      else                   ! called from 'wrf2lga'
        print ('(i3,1x,f8.1,2x,f7.0,2x,f5.1,2x,f7.5,2x,f6.1,2x,f6.1,2x,f11.8)'), &
          k,pr_wrfs(icentw,jcentw,k),ht_wrfs(icentw,jcentw,k), t3_wrfs(icentw,jcentw,k), &
          sh_wrfs(icentw,jcentw,k),u3_wrfs(icentw,jcentw,k),v3_wrfs(icentw,jcentw,k),&
          om_wrfs(icentw,jcentw,k)
      endif
    enddo
    print *, "sfc:"
    print ('(f8.1,2x,f7.0,2x,f5.1,2x,f5.1,2x,f5.1,2x,f5.1,2x,f7.5,2x,f6.1,2x,f6.1)'), &
        psf_wrf(icentw,jcentw),topo_wrf(icentw,jcentw),tsf_wrf(icentw,jcentw), &
        tsk_wrf(icentw,jcentw),dsf_wrf(icentw,jcentw),tvb_wrf(icentw,jcentw), &
        rsf_wrf(icentw,jcentw), usf_wrf(icentw,jcentw), vsf_wrf(icentw,jcentw)
        ! added tsk_wrf by wei-ting (130312)
    
    print *, "deallocating arrays"
    deallocate(dum3d,dum3df,dum3df2,dum2dt1,dum2dt2)
    return
  end subroutine fill_wrfs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine vinterp_wrfarw2p(icaller,nzl,istatus)

    implicit none

    integer, intent(out)  :: istatus
    integer  :: i,j,k,ks,kp, ksb,kst,icaller,nzl
    real     :: lpb,lpt,lp, wgtb, wgtt
    real,parameter     :: dtdlnpbase = 50.0
    real               :: tvbot, tvbar, deltalnp
    real               :: dz
    istatus = 1

    write(6,*)' subroutine vinterp_wrfarw2p: ',icaller,nzl,pr_laps(1)

    ! loop over horizontal domain, interpolating vertically for each
    ! column
    do j = 1, nyw
      do i = 1, nxw   

        pressloop: do kp = 1,nzl
  
          ! initialize kst and ksb, which will hold the vertical sigma
          ! index values of the top and bottom bounding layers
          kst = 0
          ksb = 0 

          ! find bounding levels in raw data
          sigmaloop: do ks = 1,nzw
            if (pr_wrfs(i,j,ks) .le. pr_laps(kp)) then   

              kst = ks
              ksb = kst - 1
              exit sigmaloop
            endif
          enddo sigmaloop

          if (kst .gt. 1) then ! interpolate between two bounding points
            lp = alog(pr_laps(kp))
            lpt = alog(pr_wrfs(i,j,kst))
            lpb = alog(pr_wrfs(i,j,ksb))
            wgtb = (lpt - lp) / (lpt - lpb)
            wgtt = 1.0 - wgtb

            ! height
            ht_wrfp(i,j,kp) = wgtb * ht_wrfs(i,j,ksb) + &
                              wgtt * ht_wrfs(i,j,kst)

            ! temp
            t3_wrfp(i,j,kp) = wgtb * t3_wrfs(i,j,ksb) + &
                              wgtt * t3_wrfs(i,j,kst)

            ! sh
            sh_wrfp(i,j,kp) = wgtb * sh_wrfs(i,j,ksb) + &
                              wgtt * sh_wrfs(i,j,kst)

            if(icaller .eq. 1)then ! called from 'wrf2swim'
              qc_wrfp(i,j,kp) = wgtb * qc_wrfs(i,j,ksb) + &
                                wgtt * qc_wrfs(i,j,kst)

              qi_wrfp(i,j,kp) = wgtb * qi_wrfs(i,j,ksb) + &
                                wgtt * qi_wrfs(i,j,kst)

              qr_wrfp(i,j,kp) = wgtb * qr_wrfs(i,j,ksb) + &
                                wgtt * qr_wrfs(i,j,kst)

              qs_wrfp(i,j,kp) = wgtb * qs_wrfs(i,j,ksb) + &
                                wgtt * qs_wrfs(i,j,kst)

              if(istat_aod .eq. 1)then
                aod_wrfp(i,j,kp) = wgtb * aod_wrfs(i,j,ksb) + &
                                   wgtt * aod_wrfs(i,j,kst)
              endif

            else

            ! u3
              u3_wrfp(i,j,kp) = wgtb * u3_wrfs(i,j,ksb) + &
                                wgtt * u3_wrfs(i,j,kst)

            ! v3
              v3_wrfp(i,j,kp) = wgtb * v3_wrfs(i,j,ksb) + &
                                wgtt * v3_wrfs(i,j,kst)
 
            ! om
              om_wrfp(i,j,kp) = wgtb * om_wrfs(i,j,ksb) + &
                                wgtt * om_wrfs(i,j,kst)
            endif


          elseif (kst .eq. 1) then ! extrapolate downward
            lpt = alog(pr_wrfs(i,j,kst))
            lpb = alog(pr_laps(kp))
            deltalnp = lpb - lpt
            tvbot = tvb_wrf(i,j) + deltalnp*dtdlnpbase
            tvbar = 0.5*(tvb_wrf(i,j) + tvbot)
            ! height
 
            if ((pr_laps(kp) - pr_wrfs(i,j,1)).lt. 500.) then
              ! very small difference in pressures, so
              ! assume 10 m per 100 pa, because
              ! hypsometric eq breaks down in these cases
              dz =  0.1 * (pr_laps(kp) - pr_wrfs(i,j,1))
            else
              dz = tvbar * rog * alog(pr_laps(kp)/pr_wrfs(i,j,1))
            endif
            ht_wrfp(i,j,kp) = ht_wrfs(i,j,1) - dz
            ! temp
            t3_wrfp(i,j,kp) = tvbot/(1.+0.61*mr_wrfs(i,j,1))
            ! sh
            sh_wrfp(i,j,kp) = sh_wrfs(i,j,1)

            if(icaller .eq. 1)then ! called from 'wrf2swim'
              qc_wrfp(i,j,kp) = qc_wrfs(i,j,1)
              qi_wrfp(i,j,kp) = qi_wrfs(i,j,1)
              qr_wrfp(i,j,kp) = qr_wrfs(i,j,1)
              qs_wrfp(i,j,kp) = qs_wrfs(i,j,1)
              if(istat_aod .eq. 1)then
                aod_wrfp(i,j,kp) = aod_wrfs(i,j,1)
              endif
              
            else
            
            ! u3
              u3_wrfp(i,j,kp) = u3_wrfs(i,j,1)

            ! v3
              v3_wrfp(i,j,kp) = v3_wrfs(i,j,1)

            ! om
              om_wrfp(i,j,kp) = 0.

            endif

          else ! kst never got set .. extrapolate upward
            ! assume isothermal (above tropopause)
            t3_wrfp(i,j,kp) = t3_wrfs(i,j,nzw)

            if(icaller .eq. 1)then ! called from 'wrf2swim'
              qc_wrfp(i,j,kp) = qc_wrfs(i,j,nzw)
              qi_wrfp(i,j,kp) = qi_wrfs(i,j,nzw)
              qr_wrfp(i,j,kp) = qr_wrfs(i,j,nzw)
              qs_wrfp(i,j,kp) = qs_wrfs(i,j,nzw)
              if(istat_aod .eq. 1)then
                aod_wrfp(i,j,kp) = aod_wrfs(i,j,nzw)
              endif
            else
              u3_wrfp(i,j,kp) = u3_wrfs(i,j,nzw)
              v3_wrfp(i,j,kp) = v3_wrfs(i,j,nzw)
              om_wrfp(i,j,kp) = 0.
            endif

            dz = t3_wrfs(i,j,nzw) * rog * alog(pr_wrfs(i,j,nzw)/pr_laps(kp))
            ht_wrfp(i,j,kp) = ht_wrfs(i,j,nzw) + dz 
            ! reduce moisture toward zero
            if ( pr_laps(kp-1) .gt. pr_wrfs(i,j,nzw) ) then
              sh_wrfp(i,j,kp) = 0.5*sh_wrfs(i,j,nzw)
            else
              sh_wrfp(i,j,kp) = 0.5*sh_wrfp(i,j,kp-1)
            endif
          endif

        enddo pressloop
      enddo
    enddo

    ! if we do smoothing, do it here
    
    ! print some diagnostics
    print *, "wrf press data from center of wrf domain ",nzl
    if(icaller .eq. 1)then ! called from 'wrf2swim'
      print *, "kp  press     height   temp   sh        qc       qi      qr       qs"
      print *, "--- --------  -------  -----  -------   ------   ------  -----    -----"
    else
      print *, "kp  press     height   temp   sh       u       v      om"
      print *, "--- --------  -------  -----  -------  ------  ------ -----------"
    endif
    do k = 1,nzl
      if(icaller .eq. 1)then ! called from 'wrf2swim'
        print ('(i3,1x,f8.1,2x,f7.0,2x,f5.1,2x,f7.5,2x,4f8.5,2x,f6.1,2x,f11.8)'), &
          k,pr_laps(k),ht_wrfp(icentw,jcentw,k), t3_wrfp(icentw,jcentw,k), &
          sh_wrfp(icentw,jcentw,k),qc_wrfp(icentw,jcentw,k),qi_wrfp(icentw,jcentw,k), &
          qr_wrfp(icentw,jcentw,k),qs_wrfp(icentw,jcentw,k)
      else
        print ('(i3,1x,f8.1,2x,f7.0,2x,f5.1,2x,f7.5,2x,f6.1,2x,f6.1,2x,f11.8)'), &
          k,pr_laps(k),ht_wrfp(icentw,jcentw,k), t3_wrfp(icentw,jcentw,k), &
          sh_wrfp(icentw,jcentw,k),u3_wrfp(icentw,jcentw,k),v3_wrfp(icentw,jcentw,k),&
          om_wrfp(icentw,jcentw,k)
      endif
    enddo
    
    return
  end subroutine vinterp_wrfarw2p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hinterp_wrf2lga(nxl,nyl,nzl,istatus)

  implicit none
  integer, intent(out)  :: nxl,nyl,nzl,istatus
  integer               :: i,j,k
  real, allocatable     :: xloc(:,:), yloc(:,:),dum2d(:,:)
  real, allocatable     :: topo_wrfl(:,:)  ! wrf topo interpolated to laps grid
  real                  :: ri,rj,dtopo,dz, wgt1, wgt2,lp
  istatus = 1 

  allocate(dum2d(nxl,nyl))
  ! first, generate xloc/yloc locations
  allocate(xloc(nxl,nyl))
  allocate(yloc(nxl,nyl))
  do j = 1, nyl
    do i = 1, nxl
      call latlon_to_ij(wrfgrid,lat(i,j),lon(i,j),ri,rj)
      xloc(i,j) = ri
      yloc(i,j) = rj
    enddo
  enddo

  ! now, horizontally interpolate level by level
  do k=1,nzl
    ! height
    call interpolate_standard(nxw,nyw,ht_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,method_linear, &
                              dum2d)
    ht(:,:,k) = dum2d


    ! temp
    call interpolate_standard(nxw,nyw,t3_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,method_linear, &
                              dum2d)
    t3(:,:,k) = dum2d

    ! sh
    call interpolate_standard(nxw,nyw,sh_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,method_linear, &
                              dum2d)
    sh(:,:,k) = dum2d

    ! u3
    call interpolate_standard(nxw,nyw,u3_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,method_linear, &
                              dum2d)
    u3(:,:,k) = dum2d

    ! v3
    call interpolate_standard(nxw,nyw,v3_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,method_linear, &
                              dum2d)
    v3(:,:,k) = dum2d

    ! om
    call interpolate_standard(nxw,nyw,om_wrfp(:,:,k), &
                              nxl,nyl,xloc,yloc,method_linear, &
                              dum2d)
    om(:,:,k) = dum2d

  enddo

  print *, "lga press data from center of laps domain"
  print *, "kp  press     height   temp   sh       u       v      om"
  print *, "--- --------  -------  -----  -------  ------  ------ -----------"
  do k = 1,nzl
    print ('(i3,1x,f8.1,2x,f7.0,2x,f5.1,2x,f7.5,2x,f6.1,2x,f6.1,2x,f11.8)'), &
          k,pr_laps(k),ht(icentl,jcentl,k), t3(icentl,jcentl,k), &
          sh(icentl,jcentl,k),u3(icentl,jcentl,k),v3(icentl,jcentl,k),&
          om(icentl,jcentl,k)
  enddo

  deallocate(dum2d)

  allocate(topo_wrfl(nxl,nyl))
  ! do the surface variables
  call interpolate_standard(nxw,nyw,usf_wrf, nxl,nyl,xloc,yloc, &
      method_linear,usf)
  call interpolate_standard(nxw,nyw,vsf_wrf, nxl,nyl,xloc,yloc, &
      method_linear,vsf)
  call interpolate_standard(nxw,nyw,tsf_wrf, nxl,nyl,xloc,yloc, &
      method_linear,tsf)
  call interpolate_standard(nxw,nyw,tsk_wrf, nxl,nyl,xloc,yloc, &
      method_linear,tsk) ! added tsk by wei-ting (130312)
  call interpolate_standard(nxw,nyw,dsf_wrf, nxl,nyl,xloc,yloc, &
      method_linear,dsf)
  call interpolate_standard(nxw,nyw,rsf_wrf, nxl,nyl,xloc,yloc, &
      method_linear,rsf)
  call interpolate_standard(nxw,nyw,psf_wrf, nxl,nyl,xloc,yloc, &
      method_linear,psf)
  call interpolate_standard(nxw,nyw,topo_wrf, nxl,nyl,xloc,yloc, &
      method_linear,topo_wrfl)
  call interpolate_standard(nxw,nyw,pcp_wrf, nxl,nyl,xloc,yloc, &
      method_linear,pcp) ! added pcp by wei-ting (130312)
  where ( pcp < 0 ) ; pcp = 0 ; endwhere ! keep pcp >= 0 added by wei-ting (130312)

  ! adjust for terrain differences between wrf and laps
  print *, "adjusting wrf surface to laps surface"
  do j = 1 , nyl
    do i = 1 , nxl
      dtopo = topo_wrfl(i,j) - topo_laps(i,j)
      if (abs(dtopo) .gt. 10.) then
        if (dtopo .gt. 0) then ! move downward to laps level
          downloop:  do k = nzl , 1, -1
            if (ht(i,j,k) .le. topo_laps(i,j)) then
              dz = topo_wrfl(i,j) - ht(i,j,k)
              wgt1 = (topo_wrfl(i,j) - topo_laps(i,j)) / dz
              wgt2 = 1.0 - wgt1
              lp = wgt1*alog(pr_laps(k)) + wgt2*alog(psf(i,j))
              psf(i,j) = exp(lp)
              tsf(i,j) = wgt1*t3(i,j,k)  + wgt2*tsf(i,j)
              tsk(i,j) = wgt1*t3(i,j,k)  + wgt2*tsk(i,j) ! added tsk by wei-ting (130312)
              rsf(i,j) = wgt1*sh(i,j,k)  + wgt2*rsf(i,j)
              usf(i,j) = wgt1*u3(i,j,k)  + wgt2*usf(i,j)
              vsf(i,j) = wgt1*v3(i,j,k)  + wgt2*vsf(i,j)

              exit downloop
            endif
          enddo downloop
        else  ! move upward to laps level
          uploop:  do k = 1, nzl
            if (ht(i,j,k) .ge. topo_laps(i,j)) then
              dz = ht(i,j,k) - topo_wrfl(i,j)
              wgt2 = (topo_laps(i,j) - topo_wrfl(i,j)) / dz
              wgt1 = 1.0 - wgt2
              psf(i,j) = wgt1 * psf(i,j) + wgt2*pr_laps(k)
              tsf(i,j) = wgt1 * tsf(i,j) + wgt2*t3(i,j,k)
              tsk(i,j) = wgt1 * tsk(i,j) + wgt2*t3(i,j,k) ! added tsk by wei-ting (130312)
              rsf(i,j) = wgt1 * rsf(i,j) + wgt2*sh(i,j,k)
              usf(i,j) = wgt1 * usf(i,j) + wgt2*u3(i,j,k)
              vsf(i,j) = wgt1 * vsf(i,j) + wgt2*v3(i,j,k) 
 
              exit uploop
            endif
          enddo uploop
        endif
      endif
   
    enddo
  enddo
  deallocate(xloc,yloc,topo_wrfl)

  end subroutine hinterp_wrf2lga

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_derived_pressures(nxl,nyl,nzl)

    implicit none
    integer :: i,j,k,nxl,nyl,nzl
    real    :: wgt1, wgt2,lp1,lp2,dz,lp
    do j = 1, nyl 
      do i =  1, nxl

        ! sea-level pressure
        slploop: do k = 1, nzl
          if (ht(i,j,k) .ge. 0) then
            lp1 = alog(pr_laps(k-1))
            lp2 = alog(pr_laps(k))
            dz = ht(i,j,k) - ht(i,j,k-1)
            wgt1 = ht(i,j,k)/dz
            wgt2 = 1. - wgt1
            lp = wgt1 * lp1 + wgt2 * lp2
            slp(i,j) = exp(lp) 
            exit slploop
          endif
        enddo slploop

        ! laps reduced pressure
        redploop: do k = 1, nzl
          if (ht(i,j,k) .ge. redp_lvl) then
            lp1 = alog(pr_laps(k-1))
            lp2 = alog(pr_laps(k))
            dz = ht(i,j,k) - ht(i,j,k-1)
            wgt1 = (ht(i,j,k)-redp_lvl)/dz
            wgt2 = 1. - wgt1
            lp = wgt1 * lp1 + wgt2 * lp2
            p(i,j) = exp(lp)
            exit redploop
          endif
        enddo redploop
      enddo
    enddo
    return
  end subroutine make_derived_pressures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module wrf_lga
