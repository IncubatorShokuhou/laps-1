       subroutine lso_reader_meso(m,nvar,ncycles,laps_cycle_time,
     &        badflag,o,w,olat,olon,otime,istarttime,nx,ny,bkgd)
c
c*********************************************************************
c     subroutine reads surface observations for the advanced meso
c     analysis system for recursive and wavelet approaches.
 
c     original: john mcginley, noaa/fsl  spring 2004 
c     changes:  yuanfu xie,    noaa/fsl  spring 2004
c
c     notes: units conversion of u and v; (yuanfu xie)
c            read in background fields.   (yuanfu xie)
c
c*********************************************************************
c
c
c
c
c
c
c ....lso reader arrays
c
        real lat(m), lon(m), elev(m)
        real t(m), t_ea(m), max24t(m), min24t(m)
        real td(m), td_ea(m), rh(m), rh_ea(m)
        real dd(m), ddg(m), dd_ea(m)
        real ff(m), ffg(m), ff_ea(m)
        real alt(m), pmsl(m), pstn(m)
        real alt_ea(m), delp(m), p_ea(m)
        real vis(m), vis_ea(m)
        real solar(m), solar_ea(m)
        real sfct(m), sfct_ea(m)
        real sfcm(m), sfcm_ea(m)
        real pcp1(m), pcp3(m), pcp6(m), pcp24(m)
        real snow(m), snow_ea(m), pcp_ea(m)
        real store_cldht(m,5),pctmin
c
        integer i4time, wmoid(m), jstatus
        integer time(m), delpch(m), kkk_s(m)
        integer index(m), indexa(m), indexb(m)
        integer minobs,ihr,maxsta,istarttime
c
        character store_cldamt(m,5)*4
        character stations(m)*20, provider(m)*11
        character stations_out(m)*20
        character reptype(m)*6, autostntype(m)*6
c
c.....  lso write arrays
c
        real  o(nvar,m*ncycles),olat(m*ncycles),olon(m*ncycles),
     &   otime(m*ncycles),w(m*ncycles)
        real  utrue,vtrue,ugrid,vgrid

	real  badsfc  ! yuanfu xie modified

c arrays for data input
        integer  obstime(m),kloud(m),idp3(m)
c
c switches 
        integer on,off
c
        character providera(m)*11, providerb(m)*11
        character reptypea(m)*6, reptypeb(m)*6
c
        logical exists, flagstart
        data exists/.false./
        data flagstart/.false./
c        
c character arrys for file names, station names, time, wx symbols
c
	integer     i4prev(ncycles), laps_cycle_time, cycle
	character*9 filename, fname1, fname2
        character   atime*24, atime_cur*24, stn(m)*5, wx(m)*25
        character   stna(m)*5, stnb(m)*5
        character   dir_mon*256, dir_out*256, ext_out*31,dir_qcr*256
        character   dir_s*256,ext_s*31,units*10,comment*125,var_s*3
        character   nanvar*10
c
	! background:
	character*31 bkg_ext
 	integer      nx,ny,bkg_time
        real	     bkgd(nx,ny,ncycles,6)

	! lapse rates: see laps mdatlap.f under sfc
	real         lapse_t,lapse_td

	lapse_t = -.01167
        lapse_td = -.007
c
c
c..... start here.  first get the time from the user or scheduler.
c
        narg = iargc()
c
           call get_systime(i4time,filename,istatus)
           call i4time_fname_lp(filename,i4time,status)

	print*,'i4time: ',i4time,mod(i4time,86400)
c
c.....  set the data cycle in seconds and figure out the 
c.....  previous time variables and filenames.
c
        cycle = laps_cycle_time
        cycler=float(cycle)/3600.
        do n=1,ncycles
	 i4prev(n) = i4time -(n-1)* cycle
c	 i4prev(n) = i4time -(n+1)* cycle
        enddo
        istarttime=i4prev(ncycles)
        do i=1,m
           index(i)  = 0
           indexa(i) = 0
           indexb(i) = 0
        enddo !i
c
c
c
c
c read ncycles worth of data  populate the o,w,otime,olat,olon arrays
c
  
        ! get the value for bad surface data:
        call get_sfc_badflag(badsfc,istatus)
c
	umx = -1000.0
	umm = 1000.0
	dmx = -1000.0
	dmm = 1000.0
        nobs=0 
        do n=1,ncycles
        call make_fnam_lp(i4prev(n), fname1, istatus)

c       background: time order is reverse of laps reading order - yuanfu.
        call get_background_sfc(i4prev(n),'temp',bkg_ext,bkg_time,
     &       bkgd(1,1,ncycles-n+1,1),laps_cycle_time,nx,ny,jstatus)
        if (jstatus .eq. 0) then
	   print*,'lso_reader_meso: error in reading temp background'
	   stop
	endif
	print*,'temp bkgd: ',bkgd(1,1,ncycles-n+1,1),bkg_time
        call get_background_sfc(i4prev(n),'visb',bkg_ext,bkg_time,
     &       bkgd(1,1,ncycles-n+1,4),laps_cycle_time,nx,ny,jstatus)
        if (jstatus .eq. 0) then
	   print*,'lso_reader_meso: no background for visb'
	   bkgd(1:nx,1:ny,ncycles-n+1,4) = 0.0
	endif
	print*,'visb bkgd: ',bkgd(1,1,ncycles-n+1,4),bkg_time
        call get_bkgwind_sfc(i4prev(n),bkg_ext,bkg_time,
     &       bkgd(1,1,ncycles-n+1,2),bkgd(1,1,ncycles-n+1,3),
     &       laps_cycle_time,nx,ny,jstatus)
        if (jstatus .eq. 0) then
	   print*,'lso_reader_meso: error in reading wind background'
	   stop
	endif
	print*,'wind bkgd: ',bkgd(1,1,ncycles-n+1,2:3),bkg_time
        call get_background_sfc(i4prev(n),'dewp',bkg_ext,bkg_time,
     &       bkgd(1,1,ncycles-n+1,5),laps_cycle_time,nx,ny,jstatus)
        if (jstatus .eq. 0) then
	   print*,'lso_reader_meso: error in reading dewpoint background'
	   stop
	endif
	print*,'dewp bkgd: ',bkgd(1,1,ncycles-n+1,5),bkg_time
        call get_background_sfc(i4prev(n),'redp',bkg_ext,bkg_time,
     &       bkgd(1,1,ncycles-n+1,6),laps_cycle_time,nx,ny,jstatus)
        if (jstatus .eq. 0) then
	   print*,'lso_reader_meso: error in reading redp pressure bkgrd'
	   stop
	endif
	print*,'redp bkgd: ',bkgd(1,1,ncycles-n+1,6),bkg_time
c
        call read_surface_data(i4prev(n),atime_cur,n_obs_g,n_obs_b,time,
     &     wmoid,stations,provider,wx,reptype,autostntype,lat,lon,elev,
     &     t,td,rh,dd,ff,ddg,ffg,alt,pstn,pmsl,delpch,delp,vis,solar,
     &     sfct,sfcm,pcp1,pcp3,pcp6,pcp24,snow,kkk_s,max24t,min24t,t_ea,
     &     td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,sfct_ea,
     &     sfcm_ea,pcp_ea,snow_ea,store_cldamt,store_cldht,
     &     m,jstatus)
        iflag=0
        if(jstatus .ne. 1) then
           print *,' no current lso data for ', filename
     &      
        else
           print *,' found lso data (current) at ', atime_cur
           print *,'filename: ',fname1,n_obs_g,n_obs_b
           maxsta=n_obs_b
           ! call convuv(dd,ff,u,v,maxsta,m,badflag)

	otmn = 86000.0
	otmx = 0.0
           do k=1,maxsta
            nobs=nobs+1
            otime(nobs)=time(k)
            olat(nobs)=lat(k)
            olon(nobs)=lon(k)
            o(1,nobs)=t(k)
	if (otmn .gt. time(k)) otmn = time(k)
	if (otmx .lt. time(k)) otmx = time(k)
            ! if (dd(k).ne.badflag) then
            if ((dd(k).ne.badflag) .and. (dd(k).ne.badsfc)) then 
	    ! yuanfu xie modified
             call disp_to_uv(dd(k),ff(k),utrue,vtrue)
             call uvtrue_to_uvgrid(utrue,vtrue,ugrid,vgrid,lon(k))
             o(2,nobs)=ugrid*0.5277777778     ! m/s
             o(3,nobs)=vgrid*0.5277777778     ! m/s
	if (umm .gt. o(2,nobs)) umm = o(2,nobs)
	if (umx .lt. o(2,nobs)) umx = o(2,nobs)
            else
             o(2,nobs)=badflag
             o(3,nobs)=badflag
            endif

	    ! use either altimeter or station pressure: yuanfu
	    ! this station pressure is needed for reduced pressure:
            if ((alt(k).ne.badflag) .and. (alt(k).ne.badsfc)) then
               o(4,nobs)=alt_2_sfc_press(alt(k),elev(k))
 	    else 
	       o(4,nobs)=pstn(k)
	    endif

            o(5,nobs)=td(k)
	if ((o(5,nobs) .ne. badsfc) .and. (o(5,nobs) .ne. badflag) 
     1	   .and. (dmm .gt. o(5,nobs))) dmm = o(5,nobs)
	if ((o(5,nobs) .ne. badsfc) .and. (o(5,nobs) .ne. badflag) 
     1     .and. (dmx .lt. o(5,nobs))) dmx = o(5,nobs)

	    ! 6. reduced pressure:
	    call reduce_p(t(k),td(k),o(4,nobs),elev(k),lapse_t,
     1	                  lapse_td,o(6,nobs),0.0,badflag)

	! print*,'reduced: ',o(4,nobs),elev(k),o(6,nobs)

	    ! use 4th variable for visibility:
	    o(4,nobs) = vis(k)

           enddo !k
	print*,'time range: ',otmn,otmx
        endif
           print*, nobs,maxsta,'obs read for cycle ',n
        enddo ! on n ...ob scycles
c
c...  that's it.
c     
	print*,'um : ',umm,umx
	print*,'dp : ',dmm,dmx
	w(1:10) = 0.0

        return
        end
      
