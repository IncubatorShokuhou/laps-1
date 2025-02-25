cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
c
c
	subroutine mdat_laps(i4time,atime,ni,nj,mxstn,laps_cycle_time,
     &     lat,lon,topo,x1a,x2a,y2a,redp_lvl,
     &     lon_s, elev_s, t_s, td_s, dd_s, ff_s, pstn_s, pmsl_s, alt_s, 
     &     pred_s,
     &     vis_s, stn, rii, rjj, ii, jj, n_obs_b, n_sao_g, obs,
     &     u_bk, v_bk, t_bk, td_bk, rp_bk, mslp_bk, stnp_bk, vis_bk, 
     &     wt_u, wt_v, wt_rp, wt_mslp, ilaps_bk, 
     &     u1, v1, rp1, t1, td1, sp1, tb81, mslp1, vis1, elev1,
     &     back_t,back_td,back_uv,back_sp,back_rp,back_mp,back_vis,
     &     jstatus)
c
c*******************************************************************************
c
c	rewritten version of the mcginley mesodat program
c	-- rewritten again for the laps surface analysis...1-12-88
c
c	changes:
c 	p.a. stamus	06-27-88  changed laps grid to new dimensions.
c			07-13-88  restructured for new data formats
c				                  and maps first-guess.
c			07-26-88  finished new stuff.
c			08-05-88  rotate sao winds to the projection.
c			08-23-88  changes for laps library routines.
c			09-22-88  make 1st guess optional.
c			10-11-88  make filenames time dependent.
c-------------------------------------------------------------------------------
c			12-15-88  rewritten.
c			01-09-89  changed header, added meanpres calc.
c			03-29-89  corrected for staggered grid. removed
c					dependence on file for correct times.
c			04-12-89  changed to do only 1 time per run.
c			05-12-89  fixed i4time error.
c			06-01-89  new grid -- add nest6grid.
c			11-07-89  add nummeso/sao to the output header.
c			--------------------------------------------------------
c			03-12-90  subroutine version.
c			04-06-90  pass in ihr,del,gam,ak for header.
c			04-11-90  pass in laps lat/lon and topography.
c			04-16-90  bag cloud stuff except ceiling.
c			04-17-90  add msl pressure.
c			06-18-90  new cloud check routine.
c			06-19-90  add topo.
c			10-03-90  changes for new vas data setup.
c			10-30-90  put barnes anl on boundaries.
c			02-15-91  add solar radiation code.
c			05-01-91  add hartsough qc code.
c			11-01-91  changes for new grids.
c			01-15-92  add visibility analysis.
c			01-22-93  changes for new lso/lvd.
c			07-29-93  changes for new barnes_wide routine.
c                       02-24-94  remove ceiling ht stuff.
c                       04-14-94  changes for cray port.
c                       07-20-94  add include file.
c                       09-31-94  change to lga from lma for bkgs.
c                       02-03-95  move background calls to driver routine.
c                       02-24-95  add code to check bkgs for bad winds.
c                       08-08-95  changes for new verify code.
c                       03-22-96  changes for 30 min cycle.
c                       10-09-96  grid stn elevs for use in temp anl.
c                       11-18-96  ck num obs.
c                       12-13-96  more porting changes...common for
c                                   sfc data, lgs grids. bag stations.in
c                       08-27-97  changes for dynamic laps.
c                       09-24-98  if missing background, do a smooth barnes
c                                   so something is there.
c                       09-30-98  housekeeping.
c                       01-28-99  rm smooth barnes for missing bkgs.
c                       07-08-99  removed east/west/north/south from call.
c                                   changed character array 'stn'.  rm *4
c                                   from all declarations.  use barnes for
c                                   missing backgrounds. turn off for now 
c                                   msl p calc at stns that don't report it.
c	                07-26-99  set back_mp off until can check.
c                       08-17-99  change vis for obs of 10+ miles, turn off bkg.
c                                   if lgb bkgs good, calculate a bkg red_p
c                                   instead of using previous lsx.
c                       09-17-99  add check of red_p bkg calc.  if badflags get
c                                   in, don't use it.
c
c	notes:
c
c******************************************************************************
c
	include 'laps_sfc.inc'
c
c..... stuff for the sfc data and other station info (lso +)

        include 'sfcob.inc'
        type (sfcob) obs(mxstn)
c
	real lon_s(mxstn), elev_s(mxstn)
	real t_s(mxstn), td_s(mxstn), dd_s(mxstn), ff_s(mxstn)
	real pstn_s(mxstn), pmsl_s(mxstn), alt_s(mxstn)
	real vis_s(mxstn)
	real rii(mxstn), rjj(mxstn)
        real zero_ea(mxstn) ! dummy
c
	integer ii(mxstn), jj(mxstn)
c
	character stn(mxstn)*20 
c
c.....	arrays for derived variables from obs data
c
	real uu(mxstn), vv(mxstn), pred_s(mxstn)
c
c.....	stuff for satellite data.
c
	integer lvl_v(1)
	character var_v(1)*3, units_v(1)*10
	character comment_v(1)*125, ext_v*31
c
c.....  stuff for intermediate grids (old lgs file)
c
	real u1(ni,nj), v1(ni,nj)
	real t1(ni,nj), td1(ni,nj), tb81(ni,nj)
	real rp1(ni,nj), sp1(ni,nj), mslp1(ni,nj)
	real vis1(ni,nj), elev1(ni,nj)
	real dum1(ni,nj), dum2(ni,nj)
c
c..... other arrays for intermediate grids 
c
        real wwu(ni,nj), wwv(ni,nj)
	real wp(ni,nj), wsp(ni,nj), wmslp(ni,nj)
	real wt(ni,nj), wtd(ni,nj), welev(ni,nj), wvis(ni,nj)
c
        real fnorm(0:ni-1,0:nj-1)
	real x1a(ni), x2a(nj), y2a(ni,nj)    !interp routine
	real d1(ni,nj)   ! work array
c
c..... laps lat/lon grids.
c
	real lat(ni,nj),lon(ni,nj), topo(ni,nj)
c
	real lapse_t, lapse_td
	character atime*24
c
c.....	grids for the background fields...use if not enough sao data.
c
        real u_bk(ni,nj), v_bk(ni,nj), t_bk(ni,nj), td_bk(ni,nj)
        real wt_u(ni,nj), wt_v(ni,nj)
        real rp_bk(ni,nj), mslp_bk(ni,nj), stnp_bk(ni,nj)
        real wt_rp(ni,nj), wt_mslp(ni,nj) 
        real vis_bk(ni,nj) 
        integer back_t, back_td, back_rp, back_uv, back_vis, back_sp
        integer back_mp
c
c.....  stuff for checking the background fields.
c
	real interp_spd(mxstn), bk_speed(ni,nj)
	parameter(threshold = 2.)  ! factor for diff check
	parameter(spdt      = 20.) ! spd min for diff check
	character stn_mx*5, stn_mn*5, amax_stn_id*5
c       
	integer jstatus(20)
        logical l_bilinear_interp 
c
c.....	start.  set up constants.
c
        i4_elapsed = ishow_timer()

	call tagit('mdatlaps', 19990917)

!       this will interpolate the grid to the stations using a significantly
!       faster bilinear interpolation instead of the splines.
        l_bilinear_interp = .true.

	jstatus(1) = -1		! put something in the status
	jstatus(2) = -1
	ibt = 0
	imax = ni
	jmax = nj
	icnt = 0
	delt = 0.035
        fill_val = 1.e37
        smsng = 1.e37
        zero_ea = 0.
c
c.....  zero out the sparse obs arrays.
c
	call zero(u1,    imax,jmax)
	call zero(v1,    imax,jmax)
	call zero(t1,    imax,jmax)
	call zero(td1,   imax,jmax)
	call zero(rp1,   imax,jmax)
	call zero(sp1,   imax,jmax)
	call zero(mslp1, imax,jmax)
	call zero(vis1,  imax,jmax)
	call zero(elev1, imax,jmax)
	do i=1,mxstn
	   uu(i) = 0.
	   vv(i) = 0.
	enddo !i
c
c.....  stuff for checking the background windspeed.
c
	if(ilaps_bk.ne.1 .or. back_uv.ne.1) then
	   call constant(bk_speed,badflag, imax,jmax)
	else
	   call windspeed(u_bk,v_bk,bk_speed, imax,jmax)
	endif
c
c.....	rotate sao winds to the projection grid, then change dd,fff to u,v
c
	do i=1,n_obs_b
	   if(dd_s(i).eq.badflag .or. ff_s(i).eq.badflag) then
	      uu(i) = badflag
	      vv(i) = badflag
           elseif(obs(i)%dd_ea_deg .gt. 30. .or. ! withhold based on accuracy
     1            obs(i)%ff_ea_kt  .gt. 10.)then
	      uu(i) = badflag
	      vv(i) = badflag
	   else ! dd_s is true n and dd_rot is grid n
              rlat_dum = -999.
	      dd_rot = dd_s(i) - projrot_latlon( rlat_dum , lon_s(i)
     1                                                    ,istatus )

	      dd_rot = mod( (dd_rot + 360.), 360.)
	      call decompwind_gm(dd_rot,ff_s(i),uu(i),vv(i),istatus)     
	      if(uu(i).lt.-150. .or. uu(i).gt.150.) uu(i) = badflag
	      if(vv(i).lt.-150. .or. vv(i).gt.150.) vv(i) = badflag
	   endif
	enddo !i
c
c.....  before continuing, use the sao data to check the backgrounds.
c.....  find the background at each station location, then compare
c.....  to the current observation.  if the difference is larger than the
c.....  threshold, zero out the background weights for that variable.
c
c
c.....  first find the max in the background wind speed field.
c
	print *,' '
	print *,' checking background windspeeds...'
	print *,' '
	if(ilaps_bk.ne.1 .or. back_uv.ne.1) then
	   print *,' no background wind fields availible...skipping...'
	   go to 415
	endif
c
	do j=1,jmax
	do i=1,imax
	   if(bk_speed(i,j) .gt. bksp_mx) then
	      bksp_mx = bk_speed(i,j)
	      ibksp = i
	      jbksp = j
	   endif
	enddo !i
	enddo !j

        i4_elapsed = ishow_timer()
c
c.....  find the 2nd derivative table for use by the splines later.
c
        if(.not. l_bilinear_interp)then
	    call splie2(x1a,x2a,bk_speed,imax,jmax,y2a)
        endif
c
c.....  now call the spline routine for each station in the grid.
c
	ithresh = 0
	ibkthresh = 0
	diff_mx = -1.e30
	diff_mn = 1.e30
	amax_stn = -1.e30
	do i=1,n_obs_b
	   if(ii(i).lt.1 .or. ii(i).gt.imax) go to 330
	   if(jj(i).lt.1 .or. jj(i).gt.jmax) go to 330
	   aii = float(ii(i))
	   ajj = float(jj(i))

           if(l_bilinear_interp)then
               call bilinear_laps(aii,ajj,imax,jmax,bk_speed
     1                           ,interp_spd(i))
           else
	       call splin2(x1a,x2a,bk_speed,y2a,imax,jmax,aii,ajj,
     &                 interp_spd(i))
           endif

	   if(ff_s(i) .le. badflag) then
	      diff = badflag
	   else
	      diff = interp_spd(i) - ff_s(i)
	      if(ff_s(i) .lt. 1.) then
		 percent = -1.
	      else
	         percent = ( abs(diff) / ff_s(i) ) * 100.
	      endif
	   endif
	   write(6,400) 
     &      i,stn(i)(1:5),ii(i),jj(i),interp_spd(i),ff_s(i),diff,percent
	   if(diff .eq. badflag) go to 330
	   diff = abs( diff )         ! only really care about magnitude
	   if(diff .gt. diff_mx) then
	      diff_mx = diff
	      stn_mx = stn(i)(1:5)
	   endif
	   if(diff .lt. diff_mn) then
	      diff_mn = diff
	      stn_mn = stn(i)(1:5)
	   endif
	   if(diff.gt.(threshold * ff_s(i)).and.ff_s(i).gt.spdt) then
	      ithresh = ithresh+1
	   endif
	   if(ff_s(i) .gt. amax_stn) then
	      amax_stn = ff_s(i)
	      amax_stn_id = stn(i)(1:5)
	   endif
 330	enddo !i
 400	format(1x,i5,':',1x,a5,' at i,j ',2i5,':',3f12.2,f12.0)
	write(6,405) diff_mx, stn_mx
 405	format(1x,' max difference of ',f12.2,'  at ',a5)
	write(6,406) diff_mn, stn_mn
 406	format(1x,' min difference of ',f12.2,'  at ',a5)
	write(6,410) ithresh, threshold, spdt
 410	format(1x,' there were ',i4,
     &            ' locations exceeding threshold of ',f6.3,
     &            ' at speeds greater than ',f6.1,' kts.')

        i4_elapsed = ishow_timer()

c
c.....  if too many stations exceed threshold, or if the max in the 
c.....  background is too much larger than the max in the obs, backgrounds 
c.....  probably bad.  zero out the wt arrays so they won't be used.
c
	print *,' '
	write(6,420) bksp_mx, ibksp, jbksp
 420	format(1x,' background field max: ',f12.2,' at ',i4,',',i4)
	write(6,421) amax_stn, amax_stn_id
 421	format(1x,' max speed at station: ',f12.2,' at ',a5)
c
	if(bksp_mx .ge. 60.) then
	   if(bksp_mx .gt. amax_stn*2.66) ibkthresh = 1
	endif
c
	if(ithresh.gt.2 .or. ibkthresh.gt.0) then
	   write(6,412)
 412	   format(1x,
     &      '  possible bad wind/pressure backgrounds...skipping.')
	   call zero(wt_u, imax,jmax)
	   call zero(wt_v, imax,jmax)
	   call zero(wt_rp, imax,jmax)
	   call zero(wt_mslp, imax,jmax)
	endif
	print *,' '
c
c.....	now reduce station pressures to standard levels...1500 m (for co) 
c.....  and msl.  use background 700 mb and 850 mb data from lga (or equiv).
c
cc	call mean_lapse(n_obs_b,elev_s,t_s,td_s,a_t,lapse_t,a_td,
cc     &                    lapse_td,hbar,badflag)
c
c.....  set standard lapse rates in deg f.
c
        lapse_t = -.01167
        lapse_td = -.007
c
c.....  if we have good background t, td, and station p fields from lgb (or 
c.....  fsf), and there is no reduced pressure background, use them
c.....  to calculate a reduced pressure background.  otherwise, skip this
c.....  section to use what we've already found above.
c
        write(6,*)' back_rp = ',back_rp

	if(back_t.eq.1 .and. back_td.eq.1 .and. back_sp.eq.1
     1                 .and. back_rp.ne.1) then
	   print *
     1          ,' have good backgrounds...calculating rp_bk from them.'
	   do j=1,jmax
	   do i=1,imax
	      if(stnp_bk(i,j).le.badflag .or. t_bk(i,j).le.badflag 
     &                                  .or. td_bk(i,j).le.badflag) then
		 rp_bk(i,j) = badflag
	      else
		 call reduce_p(t_bk(i,j),td_bk(i,j),stnp_bk(i,j),
     &                         topo(i,j),lapse_t,lapse_td,rp_bk(i,j),
     &                         redp_lvl,badflag) ! 1500 m for co
	      endif
	   enddo !i
	   enddo !j
c
c.....  check for badflags in the field.  don't use if we find one.
c
	   back_rp = 1
	   do j=1,jmax
	   do i=1,imax
	      if(rp_bk(i,j) .le. badflag) back_rp = 0
	   enddo !i
	   enddo !j
c
	   call check_field_2d(rp_bk,imax,jmax,fill_val,istatus)
	   if(istatus .ne. 1) back_rp = 0
c
	   print *,' done.  back_rp = ', back_rp
	endif
c
c
c.....  now, back to the analysis.
c.....	convert altimeters to station pressure. original sta pres is superceded
c
        write(6,*)' original station pressures utilized'
 415	do j=1,n_obs_b
	  if(alt_s(j) .gt. badflag) then
	    pstn_s(j) = alt_2_sfc_press(alt_s(j), elev_s(j)) !conv alt to sp
          elseif(pstn_s(j) .gt. badflag)then
            write(6,983) j, stn(j)(1:5), pstn_s(j)
          endif
	enddo !j
c
c       we are now leaving all the original pmsl_s obs in place, regardless
c       of whether pstn_s obs are present or not.
c
	sum_diffp = 0.
	num_diffp = 0
	print *,' '
	print *,' calculating reduced pressures'
	print *,'----------------------------------'
	do k=1,n_obs_b
	  if(pstn_s(k).le.badflag .or. t_s(k).le.badflag 
     &                           .or. td_s(k).le.badflag    
!    &                           .or. obs(i)%t_ea_f  .gt. 8.0
!    &                           .or. obs(i)%td_ea_f .gt. 8.0
     &                                                       )then
	    pred_s(k) = badflag
!           pmsl_s(k) = badflag
	  else
	    call reduce_p(t_s(k),td_s(k),pstn_s(k),elev_s(k),lapse_t,
     &                       lapse_td,pred_s(k),redp_lvl,badflag)  ! 1500 m for co
	    call reduce_p(t_s(k),td_s(k),pstn_s(k),elev_s(k),lapse_t,
     &                       lapse_td,p_msl,0.,badflag)        ! msl
	    if(pmsl_s(k).gt.900. .and. pmsl_s(k).lt.1100.) then
              if(p_msl .ne. badflag) then
		 diff_ps = p_msl - pmsl_s(k)
		 write(6,983) k, stn(k)(1:5), p_msl, pmsl_s(k), diff_ps
		 sum_diffp = sum_diffp + diff_ps
		 num_diffp = num_diffp + 1
	      endif
	    else
cc	       pmsl_s(k) = p_msl
	    endif
	  endif

          if(pred_s(k) .ne. badflag .and. back_rp .eq. 1)then ! compare to bkg
            i = ii(k)
            j = jj(k)
            if(i.ge.1 .and. i.le.ni .and. j.ge.1 .and. j.le.nj)then
              diff = pred_s(k) - rp_bk(ii(k),jj(k))
              if(abs(diff) .gt. 10.)then
                  write(6,981)k,stn(k)(1:5),pred_s(k),rp_bk(ii(k),jj(k))
     1                       ,diff,elev_s(k)
 981		  format(' reduced pressure flagged for ',i5,2x,a,4f9.1)
                  pred_s(k) = badflag
              endif
            endif
          endif

        enddo !k
	print *,' '
	if(num_diffp .le. 0) then
	   print *,' bad num_diffp'
	else
	   bias = sum_diffp / float(num_diffp)
	   print *,'num: ', num_diffp,'   msl pressure bias = ', bias
	endif
 983    format(1x,i5,2x,a8,':',3f12.2)
	print *,' '
c
c.....	change vis observations that are more than 10 miles to 11 miles,
c.....  so the high end of the analysis will be "greater than 10 miles".
c.....  then convert visibility to log( vis ) for the analysis.
c
	do i=1,n_obs_b
	   if(vis_s(i) .ge. 11.) vis_s(i) = 11.
	enddo !i
	call viss2log(vis_s,mxstn,n_obs_b,badflag)
c
c.....	read in the band 8 brightness temps (deg k)
c
	ext_v = 'lvd'
	lvl_v(1) = 0
	var_v(1) = 's8w'	! satellite...band 8, warm pixel (k)
c
	call get_laps_2dvar(i4time,970,i4time_nearest,lat,lon,
     &                      dum1,dum2,
     &                      ext_v,var_v,
     &        units_v,comment_v,imax,jmax,tb81,lvl_v,istatus)
c
	if(istatus .ne. 1) then
	   write(6,962) atime
 962	   format(1x,' +++ satellite data not available for the ',a24,
     &           ' analysis. +++')
	   call zero(tb81, imax,jmax)
	   go to 800
	endif
	ibt = 1
c
c.....  read in any other data here
c
 800	continue
c
c.....	put the data on the grids - accumulate sums.
c
c.....	winds:
c
	call put_winds(uu,vv,mxstn,n_obs_b,u1,v1,wwu,wwv,icnt,
     &                 imax,jmax,rii,rjj,ii,jj,badflag)
	icnt_t = icnt

        ea_thr = 999. ! lower this to enable qc for each variable
c
c.....	temperatures:
c
	print *,' put_thermo for t:'
	call put_thermo(t_s,obs(:)%t_ea_f,mxstn,n_obs_b,t1,wt,icnt,
     &                  imax,jmax,ii,jj,ea_thr,badflag)
c
c.....	dew points: 
c
	print *,' put_thermo for td:'
	call put_thermo(td_s,obs(:)%td_ea_f,mxstn,n_obs_b,td1,wtd,icnt,
     &                  imax,jmax,ii,jj,ea_thr,badflag)
c
c.....	put the reduced pressure on the grid
c
	print *,' put_thermo for p:'
	call put_thermo(pred_s,zero_ea,mxstn,n_obs_b,rp1,wp,icnt,
     &                  imax,jmax,ii,jj,ea_thr,badflag)
c
c.....	put the station pressure on the grid
c
	print *,' put_thermo for ps:'
	call put_thermo(pstn_s,zero_ea,mxstn,n_obs_b,sp1,wsp,icnt,
     &                  imax,jmax,ii,jj,ea_thr,badflag)
c
c.....	put the msl pressure on the grid
c
	print *,' put_thermo for msl:'
	call put_thermo(pmsl_s,zero_ea,mxstn,n_obs_b,mslp1,wmslp,icnt,
     &                  imax,jmax,ii,jj,ea_thr,badflag)
c
c.....	visibility:
c
	print *,' put_thermo for vis:'
	call put_thermo(vis_s,zero_ea,mxstn,n_obs_b,vis1,wvis,icnt,
     &                  imax,jmax,ii,jj,ea_thr,badflag)
c
c.....	station elevation:
c
	print *,' put_thermo for elev:'
	call put_thermo(elev_s,zero_ea,mxstn,n_obs_b,elev1,welev,icnt,
     &                  imax,jmax,ii,jj,ea_thr,badflag)
c
c.....	now find the values at the gridpts.
c       divide by the weights to average multiple obs on a gridpoint if needed
c
        write(6,1010) icnt_t
1010	format(1x,'data set 1 initialized with ',i6,' observations')
        call procar(u1,imax,jmax,wwu,imax,jmax,-1)
        call procar(v1,imax,jmax,wwv,imax,jmax,-1)
        call procar(t1,imax,jmax,wt,imax,jmax,-1)
        call procar(td1,imax,jmax,wtd,imax,jmax,-1)
        call procar(rp1,imax,jmax,wp,imax,jmax,-1)
        call procar(sp1,imax,jmax,wsp,imax,jmax,-1)
        call procar(mslp1,imax,jmax,wmslp,imax,jmax,-1)
        call procar(vis1,imax,jmax,wvis,imax,jmax,-1)
        call procar(elev1,imax,jmax,welev,imax,jmax,-1)
c
c.....  now that the data is ready, check the backgrounds.  if they
c.....  are missing, fill the background field for the variable
c.....  with a smooth barnes analysis of the obs.  this will allow
c.....  us to cold start the analysis, or run the analysis in a 
c.....  stand-alone mode.
c
        i4_elapsed = ishow_timer()
c
	n_obs_var = 0
	npass = 1
	rom2 = 0.005
	if(back_t .ne. 1) then
	   print *,' '
	   print *,
     & ' **warning. no t background. using smooth barnes anl of obs'
	   call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(t_bk,imax,jmax,t1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(t_bk,imax,jmax,fill_val,istatus)
	endif
c
	if(back_td .ne. 1) then
	   print *,' '
	   print *,
     & ' **warning. no td background. using smooth barnes anl of obs'
	   call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(td_bk,imax,jmax,td1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(td_bk,imax,jmax,fill_val,istatus)
	endif
c
	if(back_uv .ne. 1) then  !u1/v1 already rotated to grid north
	   print *,' '
	   print *,
     & ' **warning. no wind background. using smooth barnes anl of obs'
           call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(u_bk,imax,jmax,u1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(u_bk,imax,jmax,fill_val,istatus)
	   call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(v_bk,imax,jmax,v1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(v_bk,imax,jmax,fill_val,istatus)
	endif
c
	if(back_sp .ne. 1) then
	   print *,' '
	   print *,
     & ' **warning. no sfc p background. using smooth barnes anl of obs'
	   call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(stnp_bk,imax,jmax,sp1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(stnp_bk,imax,jmax,fill_val,istatus)
	endif
c
	if(back_rp .ne. 1) then
	   print *,' '
	   print *, ' **warning. no reduced p background.',
     &              ' using smooth barnes anl of obs'
	   call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(rp_bk,imax,jmax,rp1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(rp_bk,imax,jmax,fill_val,istatus)
	endif
c
	if(back_mp .ne. 1) then
	   print *,' '
	   print *,
     & ' **warning. no msl p background. using smooth barnes anl of obs'
	   call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(mslp_bk,imax,jmax,mslp1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(mslp_bk,imax,jmax,fill_val,istatus)
	endif
c
	back_vis = 0
	if(back_vis .ne. 1) then
	   print *,' '
	   print *,
     & ' **warning. no vis background. using smooth barnes anl of obs.'
	   call dynamic_wts(imax,jmax,n_obs_var,rom2,d,fnorm)
	   call barnes2(vis_bk,imax,jmax,vis1,smsng,mxstn,npass,fnorm)
	   call check_field_2d(vis_bk,imax,jmax,fill_val,istatus)
	endif
c
c.....	fill in the boundary of each field.
c
c
        i4_elapsed = ishow_timer()

        i_boundary = -1 ! 0 ! 2

	print *,' boundary for    u:'
	call set_boundary(u1,imax,jmax,ii,jj,uu,u_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for    v:'
	call set_boundary(v1,imax,jmax,ii,jj,vv,v_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for    t:'
	call set_boundary(t1,imax,jmax,ii,jj,t_s,t_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for   td:'
	call set_boundary(td1,imax,jmax,ii,jj,td_s,td_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for redp:'
	call set_boundary(rp1,imax,jmax,ii,jj,pred_s,rp_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for stnp:'
	call set_boundary(sp1,imax,jmax,ii,jj,pstn_s,stnp_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for mslp:'
	call set_boundary(mslp1,imax,jmax,ii,jj,pmsl_s,mslp_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for  vis:'
	call set_boundary(vis1,imax,jmax,ii,jj,vis_s,vis_bk,n_obs_b,
     &                    badflag,mxstn,i_boundary)
	print *,' boundary for elev:'
	call set_boundary(elev1,imax,jmax,ii,jj,elev_s,topo,n_obs_b,
     &                    badflag,mxstn,i_boundary)
c
        i4_elapsed = ishow_timer()
c
c.....	check the brightness temperatures for clouds.
c
	if(ilaps_bk.eq.0 .or. back_t.eq.0) then
	   print *,' ++ no previous temperature est for cloud routine ++'
	   go to 720
	endif
	call zero(d1,imax,jmax)
	call conv_f2k(t_bk,d1,imax,jmax)
	call clouds(imax,jmax,topo,d1,badflag,tb81,i4time,
     &              laps_cycle_time,lat,lon,1.e37)
c
c.....  convert tb8 from k to f...watch for 0.0's where clds removed.
c
720	continue
	do j=1,jmax
	do i=1,imax
	   if(tb81(i,j) .gt. 400.) tb81(i,j) = 0.
	   if(tb81(i,j) .ne. 0.) then
	      tb81(i,j) = (1.8 * (tb81(i,j) - 273.15)) + 32.
	   endif
	enddo !i
	enddo !j
c
c..... that's it here....
c
	jstatus(1) = 1		! everything's ok...
	print *,' normal completion of mdatlaps'
c
	return
	end
c
c
	subroutine put_winds(u_in,v_in,max_stn,num_sta,u,v,wwu,wwv,icnt,
     &                       ni,nj,rii,rjj,ii,jj,badflag)
c
c*******************************************************************************
c
c	routine to put the u and v wind components on the laps grid...properly
c       located on the staggered u- and v- grids.
c
c	changes:
c		p.a. stamus	12-01-88  original (cut from old mdatlaps)
c				03-29-89  fix for staggered grid.
c                               02-24-94  change method to use rii,rjj locatns.
c                               08-27-97  pass in bad flag value.
c
c	inputs/outputs:
c	   variable     var type   i/o   description
c	  ----------   ---------- ----- -------------
c	   u_in            ra       i    array of u wind components at stations
c	   v_in            ra       i      "      v  "       "       "     "
c	   max_stn         i        i    max number of stations (for dimension)
c	   num_sta         i        i    number of stations in input file
c	   u, v            ra       o    u and v component grids.
c	   wwu             ra       o    weight grid for u.
c	   wwv             ra       o    weight grid for v.
c	   icnt            i        o    number of stations put on grid.
c          rii, rjj        ra       i    i,j locations of stations (real)
c          ii, jj          ia       i    i,j locations of stations (integer)
c          badflag         r        i    bad flag value.
c
c	user notes:
c
c*******************************************************************************
c
	real u_in(max_stn), v_in(max_stn), u(ni,nj), v(ni,nj)
	real wwu(ni,nj), wwv(ni,nj)
        real rii(max_stn), rjj(max_stn)
        integer ii(max_stn), jj(max_stn)
c
	zeros = 1.e-30
	call zero(wwu, ni,nj)
	call zero(wwv, ni,nj)
c
	do 10 ista=1,num_sta
c
c.....	find ixx, iyy to put data at proper location at the grid square
c
	  ixxu = ii(ista)
	  iyyu = rjj(ista) + 0.5   ! grid offset for u-grid from major grid
	  ixxv = rii(ista) + 0.5   ! grid offset for v-grid from major grid
	  iyyv = jj(ista)
	  icnt = icnt + 1
c
c.....	put wind components on the u and v grids
c
	  if(u_in(ista).eq.badflag .or. v_in(ista).eq.badflag) go to 10
	  if(u_in(ista) .eq. 0.) u_in(ista) = zeros
	  if(v_in(ista) .eq. 0.) v_in(ista) = zeros
	  if(ixxu.lt.1 .or. ixxu.gt.ni) go to 15
	  if(iyyu.lt.1 .or. iyyu.gt.nj) go to 15
	  u(ixxu,iyyu) = u_in(ista) + u(ixxu,iyyu)
	  wwu(ixxu,iyyu) = wwu(ixxu,iyyu) + 1.
15	  if(ixxv.lt.1 .or. ixxv.gt.ni) go to 10
	  if(iyyv.lt.1 .or. iyyv.gt.nj) go to 10
	  v(ixxv,iyyv) = v_in(ista) + v(ixxv,iyyv)
	  wwv(ixxv,iyyv) = wwv(ixxv,iyyv) + 1.
10	continue
c
	return
	end
c
c
	subroutine put_thermo(var_in,ea_in,max_stn,num_sta,x,w,icnt,
     &                        ni,nj,ii,jj,ea_thr,badflag)
c
c*******************************************************************************
c
c	routine to put non-wind variables on the 'major' laps grid.
c
c	changes:
c		p.a. stamus	12-01-88  original (cut from old mdatlaps)
c				03-29-89  fix for staggered grid.
c				04-19-89  added ii,jj for qc routine.
c				10-30-90  ii,jj now from 'find_ij'.
c                               02-24-94  new ii,jj arrays.
c                               08-27-97  pass in bad flag value.
c
c	inputs/outputs:
c	   variable     var type   i/o   description
c	  ----------   ---------- ----- -------------
c	   var_in          ra       i    array of the station ob. 
c	   max_stn         i        i    max number of stations (for dimension)
c	   num_sta         i        i    number of stations in input file
c	   x               ra       o    grid for the variable. 
c	   w               ra       o    weight grid.
c	   icnt            i        o    number of stations put on grid.
c          ii, jj          ia       i    i,j locations of the stations (integer)
c          badflag         r        i    bad flag value.
c
c	user notes:
c
c*******************************************************************************
c
	real var_in(max_stn), x(ni,nj), w(ni,nj)
        real ea_in(max_stn)
        integer ii(max_stn), jj(max_stn)
c
	zeros = 1.e-30
        call zero(w,ni,nj)

        icnt_var = 0
c
	do 10 ista=1,num_sta
c
	  ixx = ii(ista)
	  iyy = jj(ista)
	  icnt = icnt + 1
c
c.....	put variable on the laps grid
c
          if(ixx.lt.1 .or. ixx.gt.ni) go to 10
          if(iyy.lt.1 .or. iyy.gt.nj) go to 10
	  if(var_in(ista) .eq. badflag) go to 10
	  if(ea_in(ista) .gt. ea_thr) go to 10
	  if(var_in(ista) .eq. 0.) var_in(ista) = zeros
	  x(ixx,iyy) = var_in(ista) + x(ixx,iyy)
	  w(ixx,iyy) = w(ixx,iyy) + 1.
          icnt_var = icnt_var + 1
10	continue
c
        write(6,*)' number of obs for this field (in put_thermo) = '        
     1           ,icnt_var

	return
	end
c
c
        subroutine procar(a,imax,jmax,b,imax1,jmax1,iproc)
        real a(imax,jmax),b(imax1,jmax1)
        do 2 j=1,jmax
	   jj=j
	   if(jmax.gt.jmax1) jj=jmax1
	   do 2 i=1,imax
	      ii=i
	      if(imax.gt.imax1) ii=imax1
c
	      if(b(ii,jj) .ne. 0.) then !3,4,3
		 if(iproc .eq. -1) then
 		    a(i,j) = a(i,j) / b(ii,jj)
		 elseif(iproc .eq. 1) then
 3		    a(i,j) = a(i,j) * b(ii,jj)
		 endif
	      endif
c
	      go to 2
!    4   a(i,j)=0.
 4	      continue
 2	continue 
c   
	return  
	end
c
c
c
	subroutine set_boundary(x,imax,jmax,ii,jj,x_ob,x_bk,n_obs_b,
     &                          badflag,mxstn,i_boundary)
c
c======================================================================
c
c       routine to fill the boundaries of the grid for a variable.
c       first, we try a barnes analysis of observations near the
c       boundary, both inside and out.  if there aren't enough obs,
c       we fall back to using the background field.  this routine
c       fills 2 points all the way around the grid (the spline needs
c       this).
c
c       orginal:  p. stamus  noaa/fsl  7 jan 1999 (from fill_bounds
c                                                   and back_bounds).
c       changes:  
c
c======================================================================
c
	real x(imax,jmax)                ! array to fill boundarys
	real x_bk(imax,jmax)             ! array with background field
	real x_ob(mxstn)                 ! array of observations
	real fnorm(0:imax-1,0:jmax-1)    ! barnes weights
	real dum(imax,jmax)              ! work array
c
	integer ii(mxstn), jj(mxstn)       ! obs location (in gridpts)
c
	npass = 1
	call zero(dum,imax,jmax)
c
c.....	first, try the wide-area barnes, filling the dum array.
c
        if(i_boundary .eq. -1)then
           write(6,*)' skipping dyn wts, barnes_wide & boundary fill'
           return
        endif

	rom2 = 0.005
	call dynamic_wts(imax,jmax,0,rom2,d,fnorm)
	call barnes_wide(dum,imax,jmax,ii,jj,x_ob,n_obs_b,badflag,
     &                   mxstn,npass,fnorm,istatus)

        if(i_boundary .eq. 0)then
           write(6,*)' skipping boundary fill'
           return
        endif
c
c.....	copy the boundaries from the dummy array to the main array.
c
	if(istatus .eq. 1) then
	   print *,' filling boundary using barnes.'
c
	   do i=1,imax
	      do j=1,i_boundary
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = dum(i,j)
	      enddo !j
	      do j=jmax+1-i_boundary,jmax
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = dum(i,j)
	      enddo !j
	   enddo !i
	   do j=1,jmax
	      do i=1,i_boundary
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = dum(i,j)
	      enddo !i
	      do i=imax+1-i_boundary,imax
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = dum(i,j)
	      enddo !i
	   enddo !j
c
	else
	   print *,' problem calculating boundary. using background.'
c
	   do i=1,imax
	      do j=1,i_boundary
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = x_bk(i,j)
	      enddo !j
	      do j=jmax+1-i_boundary,jmax
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = x_bk(i,j)
	      enddo !j
	   enddo !i
	   do j=1,jmax
	      do i=1,i_boundary
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = x_bk(i,j)
	      enddo !i
	      do i=imax+1-i_boundary,imax
		 if(x(i,j).eq.0. .or. x(i,j).eq.badflag) x(i,j) = x_bk(i,j)
	      enddo !i
	   enddo !j
c
	endif
c
c....   that's it.
c
	return
	end
