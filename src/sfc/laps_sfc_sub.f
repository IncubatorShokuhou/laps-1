c


	subroutine laps_sfc_sub(i4time)
c
c
c*****************************************************************************
c
c	driver program for the laps variational surface analysis.  this 
c	program gets the correct time and passes it to the routines (former
c	programs) that get the mesonet, sao, and vas data; grid and quality
c	control it; and perform the variational analysis.
c
c
c	history:
c		p. stamus	03-09-90  original version.
c				03-22-90  check for batch or interactive.
c				03-29-90  added check for old vas data.
c				04-06-90  pass del,gam,ak to mdat for hdr.
c				04-11-90  pass lat/lon & topo to routines
c				12-10-90  change to rt_dev for lat/lon, topo.
c				11-08-91  changes for new grids.
c				10-15-92  add steve's get_laps_config call.
c				01-06-93  new version...new lso and lvd stuff.
c				04-19-93  common for fnorm (barnes routine)
c				08-05-93  replace fnorm with fnorm2. 
c                               12-08-93  change to format in 'systime.dat'.
c                               03-04-94  change static file read (read to rd)
c                               07-20-94  add include file.
c                               02-03-95  background reads here...then pass.
c                               07-20-95  bag fnorm2 for new barnes wt method.
c                               08-08-95  add code for verify routine.
c                               03-22-96  fixes for 30 min cycle.
c                               04-10-96  more 30 min...bkgwts file.
c                               11-07-96  use 3d sfc wind for bkg, adj wts.
c                               12-13-96  more porting changes...common for
c                                         sfc data, lgs grids. bag stations.in
c                               03-26-97  add ability to do interactive runs.
c                                         (removes need for 'laps_sfci')
c                                         remove equivs.
c                               09-11-97  changes for dynamic laps.
c                               01-20-98  move wt calcs to correct place.
c	                        07-11-98  new background calls.  still temp fix
c                                           until can remove weight code (new 
c                                           spline doesn't use).
c                               09-24-98  carry background flags for each var.
c                                           rm ceil qc check.
c                               09-30-98  housekeeping.
c                               12-02-98  remove status check for lt1.
c                               07-07-99  general upgrades and cleanup.  change
c                                           read_surface_obs to read_surface_data.
c                                           rm *4 from all declarations.
c                               09-19-99  check t/td bkgs until lgb can do it.
c                               12-01-99  rotate bkg winds to grid north.
c                               12-17-99  add option to use either lso or the
c                                           kalman estimate lso_qc.
c                               01-23-00  make that option readable via namelist.
c                               01-26-00  skip qcdata call if using lso_qc.
c
c
c       notes:
c
c       1. running 'laps_sfc.x' makes the analysis use the time from 
c          '../sched/systime.dat'.  running 'laps_sfc.x -i' allows the
c          user to enter the desired start time.
c
c       2. the background weights for the station pressure (wt_stnp) are
c          not currently used.  they are included because they might be in
c          the future.
c
c*****************************************************************************

!       global parameters
        use mem_namelist, only: ni=>nx_l,nj=>ny_l,nk=>nk_laps
     1                         ,mxstn=>maxstns
     1                         ,laps_cycle_time
     1                         ,grid_spacing=>grid_spacing_m
     1                         ,max_snd_grid,max_snd_levels

!       surface parameters
        use mem_namelist, only: use_lso_qc,skip_internal_qc 
     1                         ,itheta, redp_lvl, del, gam, ak
     1                         ,l_require_lso
     1                         ,bad_t,bad_td,bad_u,bad_v,bad_p
     1                         ,bad_mp,bad_th,bad_the
     1                         ,bad_vis,bad_tb8
     1                         ,thresh_t,thresh_td,thresh_mslp
     1                         ,rms_wind,rms_temp,rms_dewpoint
     1                         ,bad_tgd_land,bad_tgd_water

!       static fields
        use mem_grid, only: lat,lon,topo,ldf

	include 'laps_sfc.inc'
c
	real x1a(ni), x2a(nj), y2a(ni,nj)
c
	integer i4time
	integer jstatus(20)		! 20 is standard for prodgen drivers
	integer narg, iargc
c
	character atime*24, filename*9, filename_last*9
	character infile_last*256
	character units*10,comment*125,var_s*3
	character use*6
c
c.....  stuff for backgrounds.
c
	real u_bk(ni,nj), v_bk(ni,nj), t_bk(ni,nj), td_bk(ni,nj) ! kt & deg f
	real wt_u(ni,nj), wt_v(ni,nj), wt_t(ni,nj), wt_td(ni,nj)
	real rp_bk(ni,nj), mslp_bk(ni,nj), stnp_bk(ni,nj)
	real wt_rp(ni,nj), wt_mslp(ni,nj), wt_stnp(ni,nj)
	real vis_bk(ni,nj), wt_vis(ni,nj)
	real tgd_bk_f(ni,nj), wt_tgd(ni,nj)
	real wt(ni,nj)
        real dum_2d(ni,nj)
c
	integer back_t, back_td, back_rp, back_uv, back_vis, back_sp
	integer back_mp, back_tgd
	character var_req*4, ext_bk*31, back*9
c
c..... stuff for the sfc data and other station info (lso +)
c
        include 'sfcob.inc'
        type (sfcob) obs(mxstn)

	real lat_s(mxstn), lon_s(mxstn), elev_s(mxstn)
	real t_s(mxstn), t_ea(mxstn), max24t(mxstn), min24t(mxstn)
        real td_s(mxstn), td_ea(mxstn), rh_s(mxstn), rh_ea(mxstn)

        real dd_s(mxstn), ddg_s(mxstn), dd_ea(mxstn)
        real ff_s(mxstn), ffg_s(mxstn), ff_ea(mxstn)

        real alt_s(mxstn), alt_ea(mxstn), delp(mxstn)
	real pstn_s(mxstn), pmsl_s(mxstn), p_ea(mxstn), pred_s(mxstn)

	real store_hgt(mxstn,5) 

        real vis_s(mxstn), vis_ea(mxstn)
        real solar_s(mxstn), solar_ea(mxstn)

        real sfct(mxstn), sfct_ea(mxstn)
        real sfcm(mxstn), sfcm_ea(mxstn)
        real pcp1(mxstn), pcp3(mxstn), pcp6(mxstn), pcp24(mxstn)
        real snow(mxstn), snow_ea(mxstn), pcp_ea(mxstn)

	real rii(mxstn), rjj(mxstn)
c
	integer kloud_s(mxstn), obstime(mxstn)
        integer wmoid(mxstn), delpch(mxstn)
	integer ii(mxstn), jj(mxstn)
c
	character atime_s*24
	character store_amt(mxstn,5)*4
        character stations(mxstn)*20, provider(mxstn)*11,stn20*20      
        character reptype(mxstn)*6, autostntype(mxstn)*6
        character wx_s(mxstn)*25 
c
c.....  work arrays for the qc routine.
c
        integer rely(26,mxstn)
	character stn3(mxstn)*3
c
c       stuff for surface verification history
        integer mo,mf,mt
        parameter (mo=100000)
        parameter (mf=10)
        parameter (mt=24)

        character*5 stn_a(mo)
        real bkg_a(mo,mf,mt)
        real obs_a(mo,mf,mt)
        real diff_a(mo,mf,mt)
        real bias_a(mo,mf) 
        real obs_mean(mo,mf)      
        real obs_std(mo,mf) 

        real lapse_t, lapse_td
c
c.....  stuff for intermediate grids (old lgs file)
c
	real u1(ni,nj), v1(ni,nj)
	real t1(ni,nj), td1(ni,nj), tb81(ni,nj)
	real rp1(ni,nj), sp1(ni,nj), mslp1(ni,nj)
	real vis1(ni,nj), elev1(ni,nj)
c
c
c*************************************************************
c.....	start here.  first see if this is an interactive run.
c*************************************************************
c
	call tagit('laps_sfc',20000126)
	narg = iargc()
cc	print *,' narg = ', narg
c
c.....  now get the time tendencies
c
        call make_fnam_lp(i4time,filename,istatus)

	if(narg .eq. 0) then
c
	   ihours = 1	! default # of hrs back for time-tendencies
c
	else
c
 974	   write(6,975)
 975	   format(
     &    ' hours to go back for time-tendencies (1..6)[1]: ',$)
	   read(5,976) ihours
 976	   format(i1)
	   if(ihours .eq. 0) then
	      print *, ' will skip backgrounds.'
	   elseif(ihours.lt.1 .or. ihours.gt.6) then
	      print *, ' error: hrs out of bounds.  try again.'
	      go to 974
	   endif
c
	endif
c
c
500	dt = ihours * laps_cycle_time    ! multiples of time to get bkgs.
	i4time_last = i4time - ifix( dt )
	call cv_i4tim_asc_lp(i4time, atime, status)	   ! get the atime
	call make_fnam_lp(i4time_last,filename_last,istatus)  ! make earlier filename	
c
c
!       fill surface namelist data structure
        sfc_nl_parms%rms_wind = rms_wind
        sfc_nl_parms%rms_temp = rms_temp
        sfc_nl_parms%rms_dewpoint = rms_dewpoint
        sfc_nl_parms%bad_tgd_land  = bad_tgd_land
        sfc_nl_parms%bad_tgd_water = bad_tgd_water

        istat = init_timer()
c
c.....  continue...get static stuff.
c
 501	continue

!       compute stats on topo data
        call stats(topo,ni,nj)
        write(6,*)' center terrain height is ',topo(ni/2,nj/2)
c
c.....  read in the obs and calculate a weight based on distance to each
c.....  station.
c
c.....	read in the surface obs:  dd/ff in deg/kt, t and td in f, elev in m,
c.....	                          and the pressure variable. cld hts are msl.
c
c
	if(use_lso_qc .eq. 1) then
	   use = 'lso_qc'
	else
	   use = 'lso   '
	endif
c
	write(6,305) filename(1:9), use
 305	format(' getting surface data at: ',a9,' from ',a6)
c
c zero out values
        n_obs_g = 0
        n_obs_b = 0
	if(use_lso_qc .ne. 1) then
	   call read_surface_data(i4time,atime_s,n_obs_g,n_obs_b,      !regular lso
     &       obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &       lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &       alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &       pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &       td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &       sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,mxstn,
     &       istatus)
	else
	   call read_surface_dataqc(i4time,atime_s,n_obs_g,n_obs_b,    !qc'd lso
     &       obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &       lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &       alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &       pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &       td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &       sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,mxstn,
     &       istatus)
	endif
c
        if( (istatus .ne. +1 .and. istatus .ne. -1) 
     1                       .or.
     1      (istatus .eq. -1 .and. l_require_lso)     )then
            write(6,*)' bad istatus from read_surface_data..',istatus       
            stop
        endif

	print *,' '
	write(6,320) use,atime_s,n_obs_g,n_obs_b
320	format(1x,a6,' data valid time: ',a24,' num obs: ',2i6)
c

        if(.true.)then
	    call read_sfc_snd(i4time,atime_s,n_obs_g,n_obs_b, ! regular snd
     &       obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &       lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &       alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &       pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &       td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &       sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,mxstn,
     &       lat,lon,ni,nj,nk,                                           ! i
     &       max_snd_grid,max_snd_levels,                                ! i
     &       topo,                                                       ! i
     &       istatus)
        endif

	print *,' '
c
c.....  copy 3 characters of the station name to another array for use
c.....  in the qc_data routine.  this can be removed once the qc is changed
c.....  to the stand-alone kalman.  blank out the array first.
c
	do i=1,mxstn
	   stn3(i)(1:3) = '   '
	enddo !i
c
        i_rh_convert = 0
	do i=1,n_obs_b ! preprocess the obs
!           rearrange station strings
            stn20 = stations(i)
            call right_justify(stn20)
            stn3(i) = stn20(18:20)
            call left_justify(stations(i))

!           convert rh to dewpoint if dewpoint is missing 
            if(t_s(i) .ne. badflag .and. td_s(i) .eq. badflag 
     1                             .and. rh_s(i) .ne. badflag)then
                t_c = f_to_c(t_s(i))
                dwpt_c = dwpt(t_c,rh_s(i))
                td_s(i) = c_to_f(dwpt_c)
                i_rh_convert = i_rh_convert + 1
            endif

!           put obs into data structure 
            obs(i)%stn3      = stn3(i)
            obs(i)%stn20     = stations(i)
            obs(i)%t_f       = t_s(i)
            obs(i)%t_ea_f    = t_ea(i)
            obs(i)%td_f      = td_s(i)
            obs(i)%td_ea_f   = td_ea(i)
            obs(i)%sfct_f    = sfct(i)
            obs(i)%dd_deg    = dd_s(i)
            obs(i)%dd_ea_deg = dd_ea(i)
            obs(i)%ff_kt     = ff_s(i)
            obs(i)%ff_ea_kt  = ff_ea(i)
        enddo ! i

        if(i_rh_convert .gt. 0)then
            write(6,*)'# of dewpoints converted from rh = ',i_rh_convert       
        endif
c
c.....	find the i,j location of each station, then calculate the
c.....  background weights (based on station density).
c
	call find_ij(lat_s,lon_s,lat,lon,n_obs_b,mxstn,
     &               ni,nj,ii,jj,rii,rjj)
c
        call get_r_missing_data(r_missing_data,istatus)

        do ista=1,n_obs_b ! write station info including elevation
                          ! insert elev/ldf info into obs data structure
           obs(ista)%ri = rii(ista)
           obs(ista)%rj = rjj(ista)
           obs(ista)%i = nint(rii(ista)) ! non-staggered "major" grid
           obs(ista)%j = nint(rjj(ista)) ! non-staggered "major" grid

           if(       obs(ista)%i .ge. 1 .and. obs(ista)%i .le. ni
     1         .and. obs(ista)%j .ge. 1 .and. obs(ista)%j .le. nj )then ! in domain

               elev_diff = elev_s(ista) - topo(obs(ista)%i,obs(ista)%j)
               obs(ista)%elev_diff = elev_diff

               rland_frac_grid = ldf(obs(ista)%i,obs(ista)%j) ! for testing

               if(.true.)then ! nearest gridpoint to determine land fraction
                   obs(ista)%ldf = ldf(obs(ista)%i,obs(ista)%j)

               else
                   call bilinear_laps(obs(ista)%ri,obs(ista)%rj,ni,nj
     1                               ,ldf,obs(ista)%ldf)
               endif

               write(6,999)ista,stations(ista)(1:5), reptype(ista)(1:6)
     &                    ,autostntype(ista)(1:6), rii(ista), rjj(ista)       
     &                    ,elev_s(ista),topo(obs(ista)%i,obs(ista)%j)
     &                    ,obs(ista)%elev_diff,obs(ista)%ldf
     &                    ,rland_frac_grid

           else
               write(6,999)ista,stations(ista)(1:5), reptype(ista)(1:6)
     &                    ,autostntype(ista)(1:6), rii(ista), rjj(ista)

               obs(ista)%elev_diff = r_missing_data
               obs(ista)%ldf       = r_missing_data

           endif

        enddo !ista
 999    format(1x,i5,': ',a5,2x,a6,2x,a6,' is at i,j: ',f6.1,',',f6.1
     1        ,2x,3f8.0,2x,2f8.3)
c
        i4_elapsed = ishow_timer()

        call zero(wt, ni,nj)
        if(n_obs_b .gt. 0.)then
	    call bkgwts(lat,lon,topo,n_obs_b,lat_s,lon_s,elev_s,
     &                  rii,rjj,wt,ni,nj,mxstn,istatus)
        endif
c
c.....  zero the weights then get the background data.  try for an fsf
c.....  forecast first, except for surface winds from the 3d analysis, then
c.....  a previous laps analysis. if both laps and fsf are missing, print a
c.....  warning. ground temp is read from 'get_modelfg_2d', that looks for
c.....  an fsf forecast according to the fdda parameter, otherwise lgb.
c
        i4_elapsed = ishow_timer()

        call zero(wt_u, ni,nj)
        call zero(wt_v, ni,nj)
        call zero(wt_t, ni,nj)
        call zero(wt_td, ni,nj)
        call zero(wt_rp, ni,nj)
        call zero(wt_mslp, ni,nj)
        call zero(wt_stnp, ni,nj)
        call zero(wt_vis, ni,nj)
        call zero(wt_tgd, ni,nj)
c
	back_t = 0
	back_td = 0
	back_rp = 0
	back_sp = 0
	back_mp = 0
	back_uv = 0
	back_vis = 0
	back_tgd = 0
c
	if(ihours .eq. 0) then     ! skip the backgrounds
	   ilaps_bk = 0
	   print *,' skipping backgrounds...'
	   go to 600
	endif
c
c.....  get the backgrounds.  convert units while we're here.
c
	call get_bkgwind_sfc(i4time,ext_bk,ibkg_time,u_bk,v_bk,
     &            laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_u, ni,nj)
	   call move(wt, wt_v, ni,nj)
c
	   call conv_ms2kt(u_bk,u_bk,ni,nj)
	   call conv_ms2kt(v_bk,v_bk,ni,nj)
	   back_uv = 1
	else
	   print *,'     no background available'
	   call zero(u_bk,ni,nj)
	   call zero(v_bk,ni,nj)
	endif

        i4_elapsed = ishow_timer()
c
	print *,' '
	print *,' getting temperature background....'
	var_req = 'temp'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,dum_2d,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   ilaps_bk = 1
	   call make_fnam_lp(ibkg_time,back,istatus)
!          write(6,951) ext_bk(1:6), back
	   call conv_k2f(dum_2d,t_bk,ni,nj) ! conv k to deg f
	   call move(wt, wt_t, ni,nj)
	   back_t = 1
	else
	   print *,'     no background available for ',var_req
	   call zero(t_bk,ni,nj)
	endif
c
	print *,' '
	print *,' getting msl pressure background....'
	var_req = 'mslp'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,mslp_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call multcon(mslp_bk,0.01,ni,nj) ! conv pa to mb
	   call move(wt, wt_mslp, ni,nj)
	   back_mp = 1
	else
	   print *,'     no background available for ',var_req
	   call zero(mslp_bk,ni,nj)
	endif
c
	print *,' '
	print *,' getting station pressure background....'
	var_req = 'sfcp'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,stnp_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!          write(6,951) ext_bk(1:6), back
	   call move(wt, wt_stnp, ni,nj)
	   call multcon(stnp_bk,0.01,ni,nj) ! conv pa to mb
	   back_sp = 1
	else
	   print *,'     no background available for ',var_req
	   call zero(stnp_bk,ni,nj)
	endif
c
	print *,' '
	print *,' getting visibility background....'
	var_req = 'visb'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,vis_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_vis, ni,nj)
	   call conv_m2miles(vis_bk,vis_bk,ni,nj)
	   call visg2log(vis_bk,ni,nj,badflag) ! conv miles to log(miles)
	   back_vis = 1
	else
	   print *,'     no background available for ',var_req
	   vis_bk = 1.0 ! equivalent to 10 mile visibility
	endif
c
	print *,' '
	print *,' getting dew point temperature background....'
	var_req = 'dewp'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,dum_2d,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_td, ni,nj)
	   call conv_k2f(dum_2d,td_bk,ni,nj) ! conv k to deg f
	   back_td = 1
	else
	   print *,'     no background available for ',var_req
	   call zero(td_bk,ni,nj)
	endif
	print *,' '
c
	print *,' '
	print *,' getting reduced pressure background....'
	var_req = 'redp'
	call get_background_sfc(i4time,var_req,ext_bk,ibkg_time,rp_bk,
     &          laps_cycle_time,ni,nj,istatus)
	if(istatus .eq. 1) then
	   call make_fnam_lp(ibkg_time,back,istatus)
!	   write(6,951) ext_bk(1:6), back
	   call move(wt, wt_rp, ni,nj)
	   call multcon(rp_bk,0.01,ni,nj) ! conv pa to mb
	   back_rp = 1
	else
	   print *,'     no background available for ',var_req
	   call zero(rp_bk, ni,nj)
	endif
c
	print *,' '
	print *,' getting ground temperature background....'
	var_req = 'tgd'
        call get_modelfg_2d(i4time,var_req,ni,nj,dum_2d,istatus)       
!       istatus = 0 ! temporary statement until 'get_modelfg_2d' is fixed
	if(istatus .eq. 1) then
	   call make_fnam_lp(i4time,back,istatus)
	   write(6,*) 'read successful for ',var_req
	   call move(wt, wt_tgd, ni,nj)
           call conv_k2f(dum_2d,tgd_bk_f,ni,nj) ! conv k to deg f
	   back_tgd = 1
	elseif(back_t .eq. 1 .and. back_td .eq. 1)then
           write(6,*)' using the temp & dewpoint bkg as a substitute'    
           do j = 1,nj
           do i = 1,ni
               tgd_bk_f(i,j) = t_bk(i,j)  *       ldf(i,j) 
     1                       + td_bk(i,j) * (1. - ldf(i,j))       
           enddo ! i
           enddo ! j
	else
	   print *,'     no background available for ',var_req
	   call zero(tgd_bk_f, ni,nj)
	endif
c
 951	format(3x,'using background from ',a6,' at ',a9)
c
c.....  adjust wts for winds; if above 2500 m, increase wt by order of mag
c.....  so background winds will have more influence.
c
	do j=1,nj
	do i=1,ni
	   if(topo(i,j) .gt. 2500.) then
	      wt_u(i,j) = wt_u(i,j) * 10.
	      wt_v(i,j) = wt_v(i,j) * 10.
	   endif
	enddo !i
	enddo !j
c
 600	continue

        i4_elapsed = ishow_timer()
c
c.....	qc the surface data.
c
	if(use_lso_qc .eq. 1) then        !data already qc'd
	   print *,' ** using pre-qcd lso_qc...skipping lsx qc step.'       
	   go to 521
	endif
	if(skip_internal_qc .eq. 1) then  !skip qc routine
	   print *, ' **  skipping lsx internal qc routine **  '
	   go to 521
	endif
c
	call get_directory('lso', infile_last, len)
	infile_last = infile_last(1:len) // filename_last // '.lso'
c
	call qcdata(filename,infile_last,rely,mxstn,
     &     t_s, td_s, dd_s, ff_s, ddg_s, ffg_s, pstn_s, pmsl_s, alt_s, 
     &     vis_s, stn3, rii, rjj, ii, jj, n_obs_b, n_sao_b, n_sao_g,
     &     ni,nj,t_bk,td_bk,mslp_bk,
     &     thresh_t,thresh_td,thresh_mslp,
     &     back_t,back_td,back_mp,
     &     istatus)
c
	if(istatus .eq. 1) then
	   jstatus(2) = 1
	elseif(istatus .eq. 0) then
	   jstatus(2) = 0
	   print *, ' +++ no data for qc routine. +++'
	   go to 521
	else
	   print *,
     &        ' +++ error.  problem in qc routine. +++'
	   jstatus(2) = -2
	   go to 521
	endif
c
c.....	check each of the primary analysis variables.
c
	do mm=1,n_obs_b
	  if(rely(7,mm) .lt. 0) then	! temperature
	     print *, 'qc: bad t at ',stations(mm),' with rely/value '       
     1              ,rely(7,mm),t_s(mm)
	     t_s(mm) = badflag
             obs(mm)%t_f = badflag
	  endif
 121	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(8,mm) .lt. 0) then	! dewpt
	     print *, 'qc: bad td at ',stations(mm)
     1              ,' with rely/value ',rely(8,mm),td_s(mm)
	     td_s(mm) = badflag
             obs(mm)%td_f = badflag
	  endif
 122	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(9,mm) .lt. 0) then	! wind direction
	     print *, 'qc: bad dir at ',stations(mm),' with rely/value '       
     1              ,rely(9,mm),dd_s(mm)
	     dd_s(mm) = badflag
             obs(mm)%dd_deg = badflag
	  endif
 123	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(10,mm) .lt. 0) then	! wind speed
             print *, 'qc: bad spd at ',stations(mm)
     1               ,' with rely/value ',rely(10,mm),ff_s(mm)
	     ff_s(mm) = badflag
             obs(mm)%ff_kt = badflag
	  endif
 124	enddo  !mm

	do mm=1,n_obs_b
          if(rely(14,mm) .lt. 0) then	! mslp
	      print *, 'qc: bad mslp at ',stations(mm)
     1               ,' with rely/value ',rely(14,mm),pmsl_s(mm)
	      pmsl_s(mm) = badflag
	  endif
 125	enddo  !mm

	do mm=1,n_obs_b
	  if(rely(15,mm) .lt. 0) then	! altimeter 
	      print *, 'qc: bad alt at ',stations(mm)
     1               ,' with rely/value ',rely(15,mm),alt_s(mm)
	      alt_s(mm) = badflag
	  endif
 126	enddo  !mm

	do mm=1,n_obs_b
          if(rely(25,mm) .lt. 0) then	! visibility 
	      print *, 'qc: bad vis at ',stations(mm)
     1               ,' with rely/value ',rely(25,mm),vis_s(mm)
	      vis_s(mm) = badflag
	  endif
 128	enddo  !mm

 521	continue                          

        i4_elapsed = ishow_timer()

!       just experimental for now
        call read_sfc_verif_history(51,r_missing_data
     1                             ,mo,mf,mt
     1                             ,stn_a,bkg_a,obs_a,diff_a,nsta
     1                             ,istatus)

        i4_elapsed = ishow_timer()

 	call sfc_verif_qc(r_missing_data
     1                   ,mo,mf,mt
     1                   ,stn_a,bkg_a,obs_a,diff_a,nsta
     1                   ,bias_a,obs_mean,obs_std)

        i4_elapsed = ishow_timer()

        if(.true.)then
            call apply_qc_info(n_obs_b,r_missing_data
     1                        ,mo,mf,mt
     1                        ,stn_a,bkg_a,obs_a,diff_a,nsta
     1                        ,bias_a,obs_mean,obs_std
     1                        ,mxstn,obs,stations
     1                        ,t_s,t_ea,td_s,td_ea
     1                        ,dd_s,dd_ea,ff_s,ff_ea)

        else

!           apply qc information                            
            iv_t   = 1
            iv_td  = 2
            iv_spd = 5
            iv_dir = 9

            lapse_t = -.01167
            lapse_td = -.007

            do mm = 1,n_obs_b
              do ista = 1,nsta
                if(stations(mm)(1:5) .eq. stn_a(ista))then
                   goto 130
                endif
              enddo 
              goto 150
130           continue
              if(mm .le. 50)then
                write(6,*)' found a qc match',mm,ista,stn_a(ista)
              endif

!             test biases and flag ob using expected accuracy info

!             temperature bias
              if(bias_a(ista,iv_t) .ne. r_missing_data)then
                biast_corr  = bias_a(ista,iv_t)  
     1                      - lapse_t *obs(mm)%elev_diff

                if(abs(biast_corr) .gt. 8.0)then
                  if(t_s(mm) .ne. badflag)then
                    obs(mm)%t_ea_f = abs(biast_corr)
                    t_ea(mm)       = abs(biast_corr)
                    write(6,135)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),biast_corr
     1                       ,bias_a(ista,iv_t),obs(mm)%elev_diff
     1                       ,obs(mm)%ldf
135		    format(' setting t_ea based on bias '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,2f8.2,f8.1,f7.2)
                  endif                
                endif
              endif

!             stuck temperature
              if(obs_std(ista,iv_t) .ne. r_missing_data)then
                if(obs_std(ista,iv_t)  .lt. 0.1)then
                  if(t_s(mm) .ne. badflag)then
                    obs(mm)%t_ea_f = 50.                  
                    t_ea(mm)       = 50.                 
                    write(6,136)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),obs_std(ista,iv_t) 
136		    format(' stuck temperature for   '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,f8.1,'deg')
                    do it = 1,mt              
                      write(6,*)' stuck temp time ',it 
     1                         ,obs_a(ista,iv_t,it)
     1                         ,obs_a(ista,iv_t,it)
                    enddo
                  endif
                endif
              endif

!             dewpoint bias
              if(bias_a(ista,iv_td) .ne. r_missing_data)then
                biastd_corr = bias_a(ista,iv_td) 
     1                      - lapse_td*obs(mm)%elev_diff
      
                if(abs(biastd_corr) .gt. 8.0)then
                  if(td_s(mm) .ne. badflag)then
                    obs(mm)%td_ea_f = abs(biastd_corr)
                    td_ea(mm)       = abs(biastd_corr)
                    write(6,137)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),biastd_corr
     1                       ,bias_a(ista,iv_td),obs(mm)%elev_diff
137		    format(' setting td_ea based on bias'
     1                    ,i6,i6,4x,2i5,1x,a5,1x,2f8.2,f8.1)
                  endif                
!                 if(abs(biastd_corr) .gt. 1000.)then
!                   do it = 1,mt              
!                     write(6,*)' td processing error',it 
!    1                         ,bkg_a(ista,iv_td,it)
!    1                         ,obs_a(ista,iv_td,it)
!    1                         ,diff_a(ista,iv_td,it)
!                   enddo
!                 endif
                endif
              endif

!             stuck wind direction
              if(obs_std(ista,iv_dir) .ne. r_missing_data)then
                if(obs_std(ista,iv_dir)  .lt. 3.0 .and. 
     1             obs_mean(ista,iv_spd) .gt. 0.)then
                  if(dd_s(mm) .ne. badflag)then
                    obs(mm)%dd_ea_deg = 180.                  
                    dd_ea(mm)         = 180.                 
                    write(6,138)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),obs_mean(ista,iv_dir) 
138		    format(' stuck wind direction for   '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,f8.1,'deg')
                    do it = 1,mt              
                      write(6,*)' stuck time ',it 
     1                         ,obs_a(ista,iv_dir,it)
     1                         ,obs_a(ista,iv_spd,it)
                    enddo
                  endif
                endif
              endif

!             stuck wind speed
              if(obs_mean(ista,iv_spd) .ne. r_missing_data)then
                if(obs_mean(ista,iv_spd) .le. 0.01)then   
                  if(ff_s(mm) .ne. badflag)then
                    obs(mm)%ff_ea_kt = 50.                  
                    ff_ea(mm)        = 50.                 
                    write(6,139)mm,ista,obs(mm)%i,obs(mm)%j
     1                       ,stn_a(ista),obs_mean(ista,iv_spd)               
139		    format(' stuck wind speed for       '
     1                    ,i6,i6,4x,2i5,1x,a5,1x,f8.1,'kt')
                    do it = 1,mt              
                      write(6,*)' stuck time ',it 
     1                         ,obs_a(ista,iv_dir,it)
     1                         ,obs_a(ista,iv_spd,it)
                    enddo
                  endif
                endif
              endif

150         enddo ! mm

            i4_elapsed = ishow_timer()

        endif
c
c.....  qc the backgrounds.  if td > t, set td = t...temporary fix until lgb
c.....  can do this....
c
	do j=1,nj
	do i=1,ni
	   if(td_bk(i,j) .gt. t_bk(i,j)) td_bk(i,j) = t_bk(i,j)
	enddo !i
	enddo !j
c
c.....  set up arrays for the verify routine.
c
	do i=1,ni
	   x1a(i) = float(i)
	enddo !i
	do j=1,nj
	   x2a(j) = float(j)
	enddo !j
c
c
cx....	now call mdatlaps to put the data on the grid.
c
	call mdat_laps(i4time,atime,ni,nj,mxstn,laps_cycle_time,lat,
     &     lon,topo,x1a,x2a,y2a,redp_lvl,
     &     lon_s, elev_s, t_s, td_s, dd_s, ff_s, pstn_s, pmsl_s, 
     &     alt_s, pred_s, vis_s, stations, rii, rjj, ii, jj, n_obs_b, 
     &     n_sao_g, obs,
     &     u_bk, v_bk, t_bk, td_bk, rp_bk, mslp_bk, stnp_bk, vis_bk,
     &     wt_u, wt_v, wt_rp, wt_mslp, ilaps_bk, 
     &     u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1,
     &     back_t,back_td,back_uv,back_sp,back_rp,back_mp,back_vis,
     &     jstatus) 
c
	if(jstatus(1) .ne. 1) then
	   print *,' from mdat_laps:  error return.  stop.'
	   stop 
	endif

        i4_elapsed = ishow_timer()
c
c
c.....	call lapsvanl to do the actual variational analysis, and calculate
c.....	derived variables, etc.  the output file goes to the lapsprd 
c.....	directory (machine dependent) and has the extension '.lsx'.
c
c       units at this stage are: p  in mb
c                                t  in f
c                                td in f
c                                u  in kt
c                                v  in kt
c                                dd in kt
	call laps_vanl(i4time,filename,ni,nj,nk,mxstn,
     &     itheta,redp_lvl,sfc_nl_parms,laps_cycle_time,
     &     dt,del,gam,ak,lat,lon,topo,ldf,grid_spacing, 
     &     lat_s, lon_s, elev_s, t_s, td_s, dd_s, ff_s, pstn_s, pmsl_s,       
     &     pred_s,
     &     vis_s, stations, n_obs_b, n_sao_b, n_sao_g, obs,
     &     u_bk,v_bk,t_bk,td_bk,rp_bk,mslp_bk,stnp_bk,vis_bk,tgd_bk_f,   
     &     wt_u, wt_v, wt_t, wt_td, wt_rp, wt_mslp, wt_vis, ilaps_bk, 
     &     back_t,back_td,back_uv,back_sp,back_rp,back_mp,back_vis,
     &     bad_t,bad_td,bad_u,bad_v,bad_p,bad_mp,bad_th,bad_the,
     &     bad_vis,bad_tb8,
     &     u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1,
     &     x1a,x2a,y2a,ii,jj,jstatus)
c
	if(jstatus(3) .ne. 1) then
	  print *,' from laps_vanl: error return.' 
	endif
c
c.....	that's about it...let's go home.
c
	print *,' end of laps surface analysis'
	return
c
	end


