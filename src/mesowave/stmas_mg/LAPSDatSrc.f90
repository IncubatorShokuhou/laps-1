!dis    forecast systems laboratory
!dis    noaa/oar/erl/fsl
!dis    325 broadway
!dis    boulder, co     80303
!dis
!dis    forecast research division
!dis    local analysis and prediction branch
!dis    laps
!dis
!dis    this software and its documentation are in the public domain and
!dis    are furnished "as is."  the united states government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  they assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  all modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  if significant modifications or enhancements
!dis    are made to this software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

module lapsdatsrc

!==========================================================
!  this module defines laps data retrieval functionality.
!
!  history:
!	creation: yuanfu xie	8-2005
!==========================================================

  use definition
  use mem_namelist

contains

subroutine lapsinfo

!==========================================================
!  this routine configures the stmas for its analyses.
!  tasks:
!	1. 
!
!  note: three letter variables are local; six global;
!
!  history: 
! 	creation: yuanfu xie	6-2005
!==========================================================

  implicit none

  ! local variables:
  character*9 :: fnm
  integer :: err	! error indicator
  integer :: i

  ! variables for accessing wind parameters:
  character*150 :: static,filenm
  integer       :: length

  !*********************
  ! laps configuration:
  !*********************

  ! get number of gridpoints:
  call get_grid_dim_xy(numgrd(1),numgrd(2),err)
  if (err .ne. 1) print*, 'stmas>lapsinfo: error getting numgrd'

  ! get laps cycle time:
  call get_laps_cycle_time(lapsdt,err)
  if (err .ne. 1) print*, 'stmas>lapsinfo: error getting cycle time'

  ! get current system time:
  call get_systime(i4time,fnm,err)
  write(6,1) fnm(1:9), 'lso'	! assume stmas using lso instead of lso_qc
1 format('stmas>lapsinfo: getting surface data at: ',a9,' from ',a6)

  ! get a flag for missing data:
  call get_r_missing_data(mising,err)
  if (err .ne. 1) print*, 'stmas>lapsinfo: error getting mising flag'

  ! get a flag for bad surface data:
  call get_sfc_badflag(badsfc,err)
  if (err .ne. 1) print*, 'stmas>lapsinfo: error getting badsfc flag'

  call get_maxstns(mxstts,err)
  if (err .ne. 1) print*, 'stmas>lapsinfo: error getting maxstations'

  ! use a new laps wind parameter scheme:
  call get_directory('static',static,length)
  filenm = static(1:length)//'/wind.nl'
  call read_namelist_laps ('wind',filenm)

  ! maximum number of stations includes both surface and snd surface:
  mxstts = mxstts+max_pr

  ! check:
  if (verbal .eq. 1) then
    write(*,2) numgrd(1:3),lapsdt,mxstts,mxstts-max_pr,max_pr,mising,badsfc
  endif
2 format('stmas>lapsinfo: num  gridpoints: ',3i6,/, &
	 'stmas>lapsinfo: laps cycle time: ',i6,/, &
	 'stmas>lapsinfo: maxnumber sites: ',i6,/, &
	 'stmas>lapsinfo: maxnumber surfs: ',i6,/, &
	 'stmas>lapsinfo: maxnumber sonde: ',i6,/, &
	 'stmas>lapsinfo: missing/bad ids: ',e16.8,f16.5)

end subroutine lapsinfo

subroutine lapsconf

!==========================================================
!  this routine reads in necessary configuration grid data.
!
!  history: 
!	creation: 8-2005 by yuanfu xie.
!==========================================================

  implicit none

  ! local variables:
  character*9 :: fnm,ext,var
  character :: unt*60,com*60
  integer :: i,j,err	! error indicator

  ! physical grid points and spacing info:
  ext = 'nest7grid'
  var = 'lat'
  call rd_laps_static(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,latgrd,phydxy(1),err)
  if (err .ne. 1) print*,'stmas>lapsinfo: error getting lat'
  var = 'lon'
  call rd_laps_static(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,longrd,phydxy(2),err)
  if (err .ne. 1) print*,'stmas>lapsinfo: error getting lon'

  ! topography:
  call read_static_grid(numgrd(1),numgrd(2),'avg',topogr,err)
  if (err .ne. 1) then
    write(6,*) 'lapsconf: error get laps avg'
    stop
  endif

  var = 'ldf'
  call rd_laps_static(dirstc,ext,numgrd(1),numgrd(2),1, &
		      var,unt,com,lndfac,grdspc(2),err)
  if (err .ne. 1) print*,'stmas>lapsinfo: error getting lon'
  ! removing meaningless land factors:
  do j=1,numgrd(2)
    do i=1,numgrd(1)
      if (lndfac(i,j) .gt. 1.0) lndfac(i,j) = 1.0
      if (lndfac(i,j) .lt. 0.0) lndfac(i,j) = 0.0
    enddo
  enddo

  ! get reduced pressure level:
  call get_laps_redp(rdplvl,err)

  ! time step:
  grdspc(3) = lapsdt

  !*********************
  ! stmas configuration:
  !*********************

  ! analysis grid numbers:
  numgrd = numgrd+2*numfic

  ! analysis domain:
  domain(1,1:2) = 1.0-numfic(1:2)	! x/y: use grid numbers
  domain(2,1:2) = float(numgrd(1:2)-numfic(1:2))
  ! domain(1,3) = mod(i4time-(numtmf-1)*lapsdt,86400) 
  ! domain(2,3) = domain(1,3)+(numtmf-1)*lapsdt
  ! new time window using i4time instead of hhmm:
  i4wndw(1) = i4time-(numtmf-1)*lapsdt
  i4wndw(2) = i4time
  domain(1,3) = 0.0
  domain(2,3) = i4wndw(2)-i4wndw(1)

  grdspc(1:2) = 1.0			! based the domain setting

end subroutine lapsconf

subroutine lapsbkgd

!==========================================================
!  this routine reads into lga background fields from laps.
!
!  history: 
! 	creation: yuanfu xie	8-2005
!==========================================================

  implicit none

  character*31 :: ext	! extension of bkg file used (laps)
  integer :: i,j,k,err	! err: 1 normal; 0 file not found
  integer :: tim	! time of bkg file used (laps)
  integer :: iwv,id,it	! index for v component of wind

  ! set laps circle time frames:
  do i=0,numtmf-1
    i4prev(numtmf-i) = i4time-i*lapsdt
  enddo

  ! check laps time frames:
  if (verbal .eq. 1) then
    do i=1,numtmf
      write(*,11) i,i4prev(i),mod(i4prev(i),86400)
    enddo
  endif
11 format('stmas>lapsbkgd: laps time stamp',i3,':',i11,i7)

  ! read background fields:
  do j=1,numvar

    ! wind:
    if (varnam(j) .eq. 'wndu') then
      ! search index for v component of wind:
      iwv = 0
      do i=1,numvar
        if (varnam(i) .eq. 'wndv') iwv = i
      enddo
      if (iwv .eq. 0) then
	write(*,12)
	stop
      endif
12 format('stmas>lapsbkgd: v component of wind is missing!')

      ! get wind fields:
      if (needbk(j) .eq. 1) then
        do i=1,numtmf
          call get_bkgwind_sfc(i4prev(i),ext,tim, &
            bkgrnd(1,1,i,j),bkgrnd(1,1,i,iwv),lapsdt, &
            numgrd(1),numgrd(2),err)
          if (err .eq. 0) write(*,13) varnam(j),i4prev(i),i,j
        enddo
      else
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
        bkgrnd(1:numgrd(1),1:numgrd(2),i,iwv) = 0.0
      endif
    !pcp1:
    else if (varnam(j) .eq. 'pcp1') then        ! precip 1hr bkg        added by min-ken hsieh
      ! get pcp1hr fields:
      if (needbk(j) .eq. 1) then
        do i=1,numtmf
          call get_modelfg_2d(i4prev(i),'pcp',numgrd(1),numgrd(2),bkgrnd(1,1,i,j),err)
          if(err .ne. 1)then
            write(6,*)' no model first guess preicp, using zero field'
            bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
          endif
        enddo
      else
        bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,j) = 0.0
      endif
    else if (varnam(j) .eq. 'pcp3') then        ! precip 3hr bkg        added by min-ken hsieh
      ! assign zero fields:
      do i=1,numtmf
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
      enddo
    else if (varnam(j) .eq. 'pcp6') then        ! precip 6hr bkg        added by min-ken hsieh
      ! assign zero fields:
      do i=1,numtmf
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
      enddo
    else if (varnam(j) .eq. 'pc24') then        ! precip 24hr bkg       added by min-ken hsieh
      ! assign zero fields:
      do i=1,numtmf
        bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = 0.0
      enddo
    else if (varnam(j) .eq. 'tgd ') then	! ground/skin temperature added by yuanfu
      if (needbk(j) .eq. 1) then
        do i=1,numtmf
          call get_modelfg_2d(i4prev(i),varnam(j),numgrd(1),numgrd(2),bkgrnd(1,1,i,j),err)
          ! temporarily add a test of tgd. if it is small (<1.0), assume
          ! laps does not provide good tgd data. when paula or steve fix lga/lgb
          ! we can change back to test err only!!!
          ! if (err .ne. 1) then
          if (err .ne. 1 .or. bkgrnd(1,1,i,j) .lt. 1.0) then
            ! use t and td and landfactor as used in laps sfc:
            it = 0	! search temp read in already
            id = 0	! search dewp read in already
            do k=1,j-1
              if (varnam(k) .eq. 'temp') it = k
              if (varnam(k) .eq. 'dewp') id = k
            enddo
            if (it .eq. 0 .and. id .eq. 0) then
              print*,'lapsbkgd: place skin/ground temperature analysis after temp/dewp in stmas_mg.vr'
              stop
            endif
            print*,'lapsbkgd: no background for skin/ground temp; use sfc temp and dewp'
            bkgrnd(1:numgrd(1),1:numgrd(2),i,j) = &
              bkgrnd(1:numgrd(1),1:numgrd(2),i,it)*lndfac(1:numgrd(1),1:numgrd(2))+ &
              bkgrnd(1:numgrd(1),1:numgrd(2),i,id)*(1.0-lndfac(1:numgrd(1),1:numgrd(2)))
          endif
        enddo
      else
        bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,j) = 0.0
      endif
    !other fields:
    else if ((varnam(j) .ne. 'wndv') .and. &    ! v in with u
             (varnam(j) .ne. 'ceil')) then      ! no ceiling bkg
      if (needbk(j) .eq. 1) then
        do i=1,numtmf
          call get_background_sfc(i4prev(i),varnam(j),ext,tim, &
            bkgrnd(1,1,i,j),lapsdt,numgrd(1),numgrd(2),err)
          if (err .eq. 0) then
            write(*,13) varnam(j),i4prev(i),i,j
            stop
          endif
        enddo
      else
        bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,j) = 0.0
      endif
    endif
  enddo
13 format('stmas>lapsbkgd: background is not found for: ',a4,i16,2i3)

  ! save background if requested for debugging:
  if (savdat .eq. 1) then
    open(unit=10,file='stmas_bkg.dat',form='formatted')
    write(10,*) numgrd(1:2),numtmf,domain
    write(10,*) bkgrnd(1:numgrd(1),1:numgrd(2),1:numtmf,saveid)
    close(10)
  endif

end subroutine lapsbkgd

subroutine lapsobsv(m)

!==========================================================
!  this routine reads in lso observation data from laps.
!
!  history:
! 	creation: yuanfu xie	8-2005
!	modified: 6-2006 by yuanfu xie: add sfc pressure.
!==========================================================

  implicit none

  integer, intent(in) :: m	! maximum number of sites

  ! local variables:
  character*24 :: tim		! obs file time
  character :: stn(m)*20	! station names
  character :: prd(m)*11	! provider names
  character :: pwx(m)*25	! present weather
  character :: rtp(m)*6		! report type
  character :: stp(m)*6		! station type (manual/auto)
  character :: amt(m,5)*4	! cloud amount

  integer :: nog,nob		! number obs over grid/box
  integer :: nss		! sfc sonde data start number = nob+1
  integer :: nsg,nsb		! number sfc sonde obs over grid/box
  integer :: wid(m)		! wmo id
  integer :: otm(m)		! observation time
  integer :: cld(m)		! number of cloud layers

  real :: lat(m),lon(m), &	! lat/lon
	    elv(m)		! elevation
  real :: tmp(m),tmpea(m), &	! temperature/expected accuracy
	    dew(m),dewea(m), &	! dewpoint/ea
	    rhd(m),rhdea(m), &	! relative humidity/ea
	    wdi(m),wdiea(m), &	! wind direction/ea
	    spd(m),spdea(m), &	! wind speed/ea
	    gdi(m),	     &	! gust wind direction
	    gsp(m),          &	! gust wind speed
	    alt(m),altea(m), &	! altimeter/ea
	    spr(m),prsea(m), &	! station pressure/ea
	    msp(m), &		! mean sea level pressure/ea
	    pcc(m),pccea(m), &	! 3-hour pressure change character
	    pch(m),pchea(m), &	! 3-hour pressure change
	    vis(m),visea(m), &	! visibility/ea
	    sol(m),solea(m), &	! solar/ea
	    slt(m),sltea(m), &	! soil/water temperature/ea
	    slm(m),slmea(m), &	! soil moist/ea
	    pc1(m),pcpea(m), &	! 1-hour precipitation/ea
	    pc3(m),pc6(m), &
	    p24(m), &		! 3,6,24-hour precipitation
	    snw(m),snwea(m), &	! snow depth/ea
	    mxt(m),mnt(m)	! 24-hour maximum/minimum temperature
  real :: cht(m,5)		! cloud layer heights
  
  integer :: i,j,k,err,iwv,kmax ! kmax need for reading sonde sfc data
  ! old obs time treatment:
  ! integer :: hrs,mns,nit	! time: hours, minutes and mid-night

  ! no matter what sts value, get_sfc_obtime converts the obstime:
  integer :: ot,sts		! obstime conversion status 

  integer :: jmin,jmax          ! for checking the min/max obs
  real :: vmin,vmax

  real :: xyt(3)		! x, y and t
  real :: prs,alt_2_sfc_press

  real :: t_c,d_c,f_to_c,c_to_f,dwpt	! for converting rh to td

  ! terrain interpolation variables:
  integer :: ix,ix1,jy,jy1
  real :: alpha,beta,topo

  ! read observation data by laps time frames:
  numobs = 0
  do i=1,numtmf

    ! frame by frame: read_surface_dataqc or read_surface_data
    call read_surface_data(i4prev(i),tim,nog,nob, &
	otm,wid,stn,prd,pwx,rtp,stp,lat,lon,elv, &
	tmp,dew,rhd,wdi,spd,gdi,gsp,alt,spr,msp,pcc,pch, &
	vis,sol,slt,slm,pc1,pc3,pc6,p24,snw,cld,mxt,mnt, &
	tmpea,dewea,rhdea,wdiea,spdea,altea,prsea,visea, &
	solea,sltea,slmea,pcpea,snwea,amt,cht,mxstts,err)

    if (err .ne. 1) then
      print*,'lapsobsv: error in reading lso data, check!',i
    endif

    ! read in surface data from sondes: 
    ! (laps surface data is lso and snd both)
    nss = nob+1
    nsg = 0	! read_sfc_snd does not initialize nsg
    nsb = 0	! read_sfc_snd does not initialize nsb
    call get_laps_dimensions(kmax,err)
    call read_sfc_snd(i4prev(i),tim,nsg,nsb, &
        otm(nss),wid(nss),stn(nss),prd(nss),pwx(nss), &
        rtp(nss),stp(nss),lat(nss),lon(nss),elv(nss), &
        tmp(nss),dew(nss),rhd(nss),wdi(nss),spd(nss), &
        gdi(nss),gsp(nss),alt(nss),spr(nss),msp(nss),pcc(nss),pch(nss), &
        vis(nss),sol(nss),slt(nss),slm(nss),pc1(nss), &
        pc3(nss),pc6(nss),p24(nss),snw(nss),cld(nss),mxt(nss),mnt(nss), &
        tmpea(nss),dewea(nss),rhdea(nss),wdiea(nss),spdea(nss),altea(nss), &
        prsea(nss),visea(nss),solea(nss),sltea(nss),slmea(nss),pcpea(nss), &
        snwea(nss),amt(nss,1),cht(nss,1),mxstts,latgrd,longrd, &
        numgrd(1),numgrd(2),kmax,max_pr,max_pr_levels,topogr,sts)
    ! combine surface data together:
    nob = nob+nsb	! add snd data to the total

    if (nob .eq. 0) then
      print*,'lapsobsv: no sfc obs found!'
      ! stop ! cover it up for allowing analysis to proceed without obs
    else
        ! convert laps surface obs time to i4time:
      do j=1,nob
        ot = otm(j)
        ! otm from lso is in hhmm:
        if (ot .ge. 0 .and. ot .le. 2400) then
          call get_sfc_obtime(ot,i4prev(i),otm(j),sts)
        else
          ! obs time error:
          otm(j) = -100
          if (verbal .eq. 1) &
            print*,'lapsobs: invalid obs time: set to bad value: ',ot,j,i
        endif
      enddo

    endif

    if (err .ne. 1 .and. sts .ne. 1) then
      ! lso data cannot be read in:
      write(*,21) i
    else
      ! assign lso data to the corresponding arrays:

      ! retrieve location and time sequence:
      do j=1,nob	! through all obs sites
	! x and y:
	call latlon_to_rlapsgrid(lat(j),lon(j), &
		latgrd,longrd,numgrd(1),numgrd(2), &
		xyt(1),xyt(2),err)

        ! good station locations:
        if (err .eq. 1) then

          ! laps latlon_to_rlapsgrid returns xyt(1:2) ranging from:
          ! (0.5,numgrd(1:2)+0.5). stmas treats them on the edge of
          ! the domain:
          if (xyt(1) .lt. 1) xyt(1) = 1.0
          if (xyt(1) .gt. numgrd(1)) xyt(1) = numgrd(1) 
          if (xyt(2) .lt. 1) xyt(2) = 1.0
          if (xyt(2) .gt. numgrd(2)) xyt(2) = numgrd(2) 

	  ! t: from laps time form: hhmm to seconds
	  ! hrs = otm(j)/100
	  ! mns = otm(j)-hrs*100
	  ! xyt(3) = hrs*3600+mns*60
	  ! if (otm(j) .lt. 0) xyt(3) = 2*86400	! void: bad data

          ! use i4time to handle obs times:
          xyt(3) = otm(j)-i4wndw(1)

	  ! adjust the time when crossing the midnight:
	  if ((xyt(3)+86400 .ge. domain(1,3)) .and. &
	      (xyt(3)+86400 .le. domain(2,3)) ) &
	    xyt(3) = xyt(3)+86400

	  ! pass the location/time to obs arrays:
          do k=1,numvar
	    rawobs(2:4,j+numobs(k),k) = xyt(1:3)
	  enddo

        else  ! bad station location
          do k=1,numvar
            rawobs(1:4,j+numobs(k),k) = badsfc
          enddo
        endif

      enddo

      ! place the observations into right variables:
      do j=1,numvar
	stanam(1+numobs(j):nob+numobs(j),j) = stn(1:nob)	!added by min-ken,hsieh:stanam for stmasver
	
	select case (varnam(j))
	case ("temp")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = tmp(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = tmpea(1:nob)
	case ("dewp")
          ! use td obs or rh obs if td missing:
          do k=1,nob
            if (dew(k) .ne. mising .and. &
                dew(k) .ne. badsfc) then
              rawobs(1,numobs(j)+k,j) = dew(k)
	      weight(numobs(j)+k,j) = dewea(k)
            elseif (rhd(k) .ne. mising .and. &
                    rhd(k) .ne. badsfc .and. &
                    tmp(k) .ne. mising .and. &
                    tmp(k) .ne. badsfc) then
                ! obs value:
                t_c = f_to_c(tmp(k))
                d_c = dwpt(t_c,rhd(k))
                rawobs(1,numobs(j)+k,j) = c_to_f(d_c)
                dew(k) = rawobs(1,numobs(j)+k,j) ! replace dew with one from rh
                ! obs expect accuracy:
                t_c = f_to_c(tmpea(k))
                d_c = dwpt(t_c,rhdea(k))
                weight(numobs(j)+k,j) = c_to_f(d_c)
            else
              rawobs(1,numobs(j)+k,j) = badsfc
	      weight(numobs(j)+k,j) = badsfc
            endif
          enddo
	case ("visb")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = vis(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = visea(1:nob)
        case ("ceil")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = cht(1:nob,1)
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0
        case ("mslp")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = msp(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = prsea(1:nob)
        case ("tgd ")
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = slt(1:nob)	! soil temp
	  weight(1+numobs(j):nob+numobs(j),j) = sltea(1:nob)
	case ("redp")
	  do k=1,nob
	    ! collect either station pressure or altimeter:
	    if ((spr(k) .ne. mising) .and. &
		(spr(k) .ne. badsfc)) then
	      prs = spr(k)
	    elseif ((alt(k) .ne. mising) .and. &
		    (alt(k) .ne. badsfc) .and. &
                    (rawobs(2,numobs(j)+k,j) .ne. mising) .and. &
                    (rawobs(3,numobs(j)+k,j) .ne. badsfc)) then

              ! grid indices and interpolation coefficiences
              ix = int(rawobs(2,numobs(j)+k,j))
              ix1 = min(ix+1,numgrd(1))
              jy = int(rawobs(3,numobs(j)+k,j))
              jy1 = min(jy+1,numgrd(2))

              alpha = rawobs(2,numobs(j)+k,j)-ix
              beta  = rawobs(3,numobs(j)+k,j)-jy

              topo = (1.0-alpha)*(1.0-beta)*topogr(ix ,jy )+ &
                          alpha *(1.0-beta)*topogr(ix1,jy )+ &
                     (1.0-alpha)*     beta *topogr(ix ,jy1)+ &
                          alpha *     beta *topogr(ix1,jy1)

	      prs = alt_2_sfc_press(alt(k),topo)
	    else
	      prs = badsfc
	    endif

	    ! convert to reduced pressure at the reduced pressure level, rdplvl:
            if (prs .ne. badsfc) then
              if (abs(elv(k)-rdplvl) .le. 10) then 
                ! elevation is close (10m hardcoded) to the reduced level, 
                ! use pressure as reduced pressure:
                rawobs(1,numobs(j)+k,j) = prs
              else
                ! use reduce_p to reduce prs to the reduced pressure level:
	        if (tmp(k) .ne. badsfc .and. dew(k) .ne. badsfc) then
	          call reduce_p(tmp(k),dew(k),prs,topo, &
		                lapses(1),lapses(2), &
		                rawobs(1,numobs(j)+k,j),rdplvl,badsfc)
                else
	          rawobs(1,numobs(j)+k,j) = badsfc
                endif
              endif
	    else
	      rawobs(1,numobs(j)+k,j) = badsfc
            endif
	  enddo
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0 ! altea(1:nob)
	case ("sfcp")
	  do k=1,nob
	    ! collect either station pressure or altimeter:
	    if ((spr(k) .ne. mising) .and. &
		(spr(k) .ne. badsfc)) then
	      rawobs(1,numobs(j)+k,j) = spr(k)
	    elseif ((alt(k) .ne. mising) .and. &
		    (alt(k) .ne. badsfc) .and. &
                    (rawobs(2,numobs(j)+k,j) .ne. mising) .and. &
                    (rawobs(3,numobs(j)+k,j) .ne. badsfc) ) then

              ! grid indices and interpolation coefficiences
              ix = int(rawobs(2,numobs(j)+k,j))
              ix1 = min(ix+1,numgrd(1))
              jy = int(rawobs(3,numobs(j)+k,j))
              jy1 = min(jy+1,numgrd(2))

              alpha = rawobs(2,numobs(j)+k,j)-ix
              beta  = rawobs(3,numobs(j)+k,j)-jy

              topo = (1.0-alpha)*(1.0-beta)*topogr(ix ,jy )+ &
                          alpha *(1.0-beta)*topogr(ix1,jy )+ &
                     (1.0-alpha)*     beta *topogr(ix ,jy1)+ &
                          alpha *     beta *topogr(ix1,jy1)

	      rawobs(1,numobs(j)+k,j) = alt_2_sfc_press(alt(k),topo)

	    else
	      rawobs(1,numobs(j)+k,j) = badsfc
	    endif
	  enddo
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0 ! altea(1:nob)
	case ("wndu")
	  ! find the index for v component:
	  iwv = 0
	  do k=1,numvar
	    if (varnam(k) .eq. "wndv") iwv = k
	  enddo
	  if (iwv .eq. 0) then
	    write(*,24)
	  endif
24 format('stmas>lapsobs: warning: no v component wind analysis!')

	  ! convert wind from direction/speed to u/v:
	  do k=1,nob
	    if ((wdi(k) .eq. mising) .or. &
	        (wdi(k) .eq. badsfc) .or. &
	        (spd(k) .eq. mising) .or. &
		(spd(k) .eq. badsfc)) then
	      rawobs(1,numobs(j)+k,j) = badsfc
	      rawobs(1,numobs(j)+k,iwv) = badsfc
	    else
	      ! conversion:
	      call disp_to_uv(wdi(k),spd(k),xyt(1),xyt(2))
	      call uvtrue_to_uvgrid(xyt(1),xyt(2), &
		rawobs(1,numobs(j)+k,j), &
		rawobs(1,numobs(j)+k,iwv),lon(k))
	    endif
	  enddo
	  weight(1+numobs(j):nob+numobs(j),j) = 1.0
	  weight(1+numobs(j):nob+numobs(j),iwv) = 1.0
	case ("wndv")
	  ! v should be in already. see case ("wndu").
	case ("pcp1")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = pc1(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	case ("pcp3")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = pc3(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	case ("pcp6")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = pc6(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	case ("pc24")						!added by min-ken hsieh
	  rawobs(1,1+numobs(j):nob+numobs(j),j) = p24(1:nob)
	  weight(1+numobs(j):nob+numobs(j),j) = pcpea(1:nob)
	case default
	  write(*,22) varnam(j)
	end select
      enddo
    endif

    ! update frm:
    numobs(1:numvar) = numobs(1:numvar)+nob
  enddo
21 format('stmas>lapsobs: cannot read in lso data and snd data: ',i8)
22 format('stmas>lapsobs: no such var in lso data: ',a4)

  ! check number of obs:
  do i=1,numvar
    write(*,23) varnam(i),numobs(i)
    ! if (verbal .eq. 1) write(*,23) varnam(i),numobs(i)
  enddo
23 format('stmas>lapsobsv: numobs of (raw) ',a4,': ',i8)

  ! remove invalid data:
  call rmvinvld

  ! remove redundant observations:
  ! call rmvdupls

  ! check the data ranges:
  if (verbal .eq. 1) then
    do i=1,numvar
      vmin = 1000.0
      vmax = -1000.0
      do j=1,numobs(i)
        if (rawobs(1,j,i) .lt. vmin) then
          vmin = rawobs(1,j,i)
          jmin = j
        endif
        if (rawobs(1,j,i) .gt. vmax) then
          vmax = rawobs(1,j,i)
          jmax = j
        endif
      enddo
       
      write(*,26) varnam(i),vmin,jmin,stanam(jmin,i),vmax,jmax,stanam(jmax,i)
!		  minval(rawobs(1,1:numobs(i),i)), &
!		  maxval(rawobs(1,1:numobs(i),i))
    enddo
  endif
26 format('stmas>lapsobsv: ',a4,' min/max values: ', f11.2,i4,1x,a6,f11.2,i4,1x,a6)

  ! write out requested obs for testing:
  if (savdat .eq. 1) then
    open(unit=10,file='stmas_ob1.dat',form='formatted')
    write(10,*) numobs(saveid),numtmf,domain,grdspc(3)
    write(10,*) rawobs(1:4,1:numobs(saveid),saveid)
    close(10)
  endif

end subroutine lapsobsv

subroutine rmvinvld

!==========================================================
!  this routine removes the invalid data (missing/badsfc)
!  data from observations.
!
!  history:
!	creation: 8-2005 by yuanfu xie
!==========================================================

  implicit none

  integer :: i,j,num

  ! for every variables:
  do i=1,numvar

    num = numobs(i)
    numobs(i) = 0
    ! for every obs:
    do j=1,num
      if ((rawobs(1,j,i) .ne. mising) .and. &
	  (rawobs(1,j,i) .ne. badsfc)) then
	! valid data:
	numobs(i) = numobs(i)+1
	rawobs(1:4,numobs(i),i) = rawobs(1:4,j,i)
	weight(numobs(i),i) = weight(j,i)
	stanam(numobs(i),i) = stanam(j,i)		!added by min-ken,hsieh:stanam for stmasver
      endif
    enddo
  enddo
 
  ! check numbers of obs left:
  do i=1,numvar
    if (verbal .eq. 1) write(*,31) varnam(i),numobs(i)
  enddo
31 format('stmas>laps_qcs: numobs of (vld) ',a4,': ',i8)

end subroutine rmvinvld

subroutine rmvdupls

!==========================================================
!  this routine removes the duplicated observations.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: i,j,k,l,nob
  real :: dis

  obsspc = 1.0e10
  do i=1,numvar

    ! use weight as flag:
    do j=1,numobs(i)
      do k=j+1,numobs(i)
	dis = 0.0
	do l=2,4
	  dis = dis+(rawobs(l,j,i)-rawobs(l,k,i))* &
		    (rawobs(l,j,i)-rawobs(l,k,i))
	enddo

	! mark those redundants:
	if (dis .lt. epsiln) then
	  weight(k,i) = 0.0
	  if (verbal .eq. 1) then
	    write(*,11) rawobs(1:4,j,i),rawobs(1:4,k,i), &
		varnam(i),j,k
	  endif
	else
	  ! minimal observation spacing:
          if ((abs(rawobs(2,j,i)-rawobs(2,k,i)) .gt. 0.0) .and. &
              (abs(rawobs(2,j,i)-rawobs(2,k,i)) .lt. obsspc(1,i))) &
            obsspc(1,i) = abs(rawobs(2,j,i)-rawobs(2,k,i))
          if ((abs(rawobs(3,j,i)-rawobs(3,k,i)) .gt. 0.0) .and. &
              (abs(rawobs(3,j,i)-rawobs(3,k,i)) .lt. obsspc(2,i))) &
            obsspc(2,i) = abs(rawobs(3,j,i)-rawobs(3,k,i))
          if ((abs(rawobs(4,j,i)-rawobs(4,k,i)) .gt. 0.0) .and. &
              (abs(rawobs(4,j,i)-rawobs(4,k,i)) .lt. obsspc(3,i))) &
            obsspc(3,i) = abs(rawobs(4,j,i)-rawobs(4,k,i))
	endif
      enddo
11 format('stmas>rmvdupls: redundant data: ',/,4f14.4,/, &
	4f14.4,/,a6,2i8)

    enddo
    write(*,12) varnam(i),obsspc(1:3,i)
12 format('stmas>rmvdupls: minimal obs (',a4,') spacing: ',3e12.4)

    ! remove redundants:
    nob = 0
    do j=1,numobs(i)
      if (weight(j,i) .gt. 0.0) then
        nob = nob+1
        rawobs(1:4,nob,i) = rawobs(1:4,j,i)
	weight(nob,i) = weight(j,i)
	stanam(nob,i) = stanam(j,i)		!added by yuanfu: stanam for stmasver
      endif
    enddo
    numobs(i) = nob

  enddo

end subroutine rmvdupls

subroutine lapsunit

!==========================================================
!  this routine converts lso observation units into a unit
!  consistent with the background.
!
!  history:
!	creation: 8-2005 by yuanfu xie.
!	modified: 6-2006 by yuanfu xie: add sfc pressure.
!==========================================================

  implicit none

  integer :: i

  ! check all variables:
  do i=1,numvar

    ! find necessary conversion:
    select case (varnam(i))
    case ("temp")
      ! convert to kelvin from fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    case ("dewp")
      ! convert to kelvin from fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    case ("tgd ")
      ! convert to kelvin from fahrenheit:
      rawobs(1,1:numobs(i),i) = &
	(rawobs(1,1:numobs(i),i)-32.0)*5.0/9.0+temp_0
    case ("wndu")
      ! convert to m/s from knots:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*knt2ms
    case ("wndv")
      ! convert to m/s from knots:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*knt2ms
    case ("visb")
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mile2m
    case ("redp")
      ! convert to pascal from mb:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mb2pas
    case ("sfcp")
      ! convert to pascal from mb:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mb2pas
    case ("mslp")
      ! convert to pascal from mb:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*mb2pas
    case ("pcp1")
      ! convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    case ("pcp3")
      ! convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    case ("pcp6")
      ! convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    case ("pc24")
      ! convert to meter from inch:
      rawobs(1,1:numobs(i),i) = &
	rawobs(1,1:numobs(i),i)*inch2m
    end select
  enddo

end subroutine lapsunit

subroutine laps_qcs

!==========================================================
!  this routine runs quality control over data by threshold
!  values and standard deviation.
!
!  history:
! 	creation: yuanfu xie	6-2005
!==========================================================

  implicit none

  ! interpolation indices and coefficients:
  call lapsintp

  ! optional qcs:
  if (qc_val .eq. 1) call thrshold

  ! save qced obs:
  call cpyqcobs

end subroutine laps_qcs

subroutine cpyqcobs

!==========================================================
!  this routine copies qced observation data from rawobs to
!  qc_obs after all qc is done. this routine can be avoid
!  if the qc_obs array keeps the same structure as rawobs.
!
!  history:
!	creation: 8-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: i,j

  ! copy:
  do i=1,numvar
    do j=1,numobs(i)
      qc_obs(1:4,j,i) = rawobs(1:4,j,i)

      ! save innovations:
      if (needbk(i) .eq. 1) &
        qc_obs(1,j,i) = qc_obs(1,j,i)-bkgobs(j,i)

    enddo
  enddo

end subroutine cpyqcobs

subroutine thrshold

!==========================================================
!  this routine does the threshold value qc checks.
!
!  history:
!	creation: 8-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: i,j,num

  ! check:
  do i=1,numvar
    if (needbk(i) .eq. 1) then
      num = numobs(i)
      numobs(i) = 0
      do j=1,num

        ! qc check: avoid bkg = mising with roundoff error:
        if (abs(rawobs(1,j,i)-bkgobs(j,i)) .le. thresh(i)) then
	  numobs(i) = numobs(i)+1
	  rawobs(1:4,numobs(i),i) = rawobs(1:4,j,i)
	  weight(numobs(i),i) = weight(j,i)
	  stanam(numobs(i),i) = stanam(j,i)		!added by min-ken,hsieh:stanam for stmasver
	  indice(1:6,numobs(i),i) = indice(1:6,j,i)
	  coeffs(1:6,numobs(i),i) = coeffs(1:6,j,i)
	  bkgobs(numobs(i),i) = bkgobs(j,i)
        else
          print*,'thresh out: ',rawobs(1,j,i),bkgobs(j,i),thresh(i),j,i,stanam(j,i)
        endif
      enddo
    endif
  enddo
 
  ! check numbers of obs left:
  do i=1,numvar
    if (verbal .eq. 1) write(*,31) varnam(i),numobs(i)
  enddo
31 format('stmas>laps_qcs: numobs of (vlu) ',a4,': ',i8)

end subroutine thrshold

subroutine lapsintp

!==========================================================
!  this routine interpolates gridpoints to observation site
!  and saves the indices and coefficients.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!       modification:
!                 25-08-2008 by min-ken hsieh
!                 pass stanam into grid2obs to map each obs its stn name for stmasver
!==========================================================

  implicit none

  integer :: i,j,ix,iy,it

  do i=1,numvar
    call grid2obs(indice(1,1,i),coeffs(1,1,i), &
	rawobs(1,1,i),numobs(i),weight(1,i),stanam(1,i),numgrd, &
	grdspc,domain)

    ! compute background values at observation sites:
    do j=1,numobs(i)

      ! interpolate background to the obs site:
      bkgobs(j,i) = 0.0
      do it=3,6,3
        do iy=2,5,3
          do ix=1,4,3
            bkgobs(j,i) = bkgobs(j,i) + &
	      bkgrnd(indice(ix,j,i), &
		     indice(iy,j,i), &
		     indice(it,j,i),i)* &
	      coeffs(ix,j,i)*coeffs(iy,j,i)*coeffs(it,j,i)
          enddo
        enddo
      enddo
    enddo

  enddo

end subroutine lapsintp

subroutine stmasver
!==========================================================
!  this routine prepare all parameters and pass them to
!  verify subroutine in src/lib/laps_routine.f.
!  it will output verify file in log/qc directory, just
!  like what laps_sfc.x dose.
!  history:
!       creation: 22-8-2008 by min-ken hsieh.
!==========================================================

  implicit none

  !local vaiables
  integer :: i, j, len, istatus, err
  integer :: iunit     		 !log file handle
  integer :: ii(mxstts),jj(mxstts)

  ! since we use bilinear interpolation, these arrays are dummy..
  real :: x1a(numgrd(1)), x2a(numgrd(2)), y2a(numgrd(1),numgrd(2)) 
  real :: ea

  character*60  :: title
  character*256 :: ver_file
  character*9   :: a9time

  !for time loop
  integer :: kt,nn
  integer :: obstime
  real :: obsofthistime(mxstts)
  character*20 :: staofthistime(mxstts)
  character*1 :: tmtag

  ! open log file
  iunit = 11
  call make_fnam_lp(i4time,a9time,istatus)
  if(istatus .eq. 0) goto 999
  call get_directory('log', ver_file, len)
  ver_file = ver_file(1:len)//'qc/stmas.ver.'//a9time(6:9)
  call s_len(ver_file, len)
  !print*, "min-ken",ver_file
  open(iunit,file=ver_file(1:len),status='unknown',err=999)

  ! variable loop
  do i=1,numvar
    !because we only care about real obs
    !we do not apply verify to those obs made by bkgrnd

    ! qc_obs array actually store obs - obsbkg field (done by cpyqcobs)
    ! we need to add obsbkg back to qc_obs
    if (needbk(i) .eq. 1) &
      qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)+bkgobs(1:numobs(i),i)



    select case (varnam(i))
      case ("temp")

	!time loop
	do kt = 1,numtmf
	  write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
	      obsofthistime(nn) = qc_obs(1,j,i)
	      staofthistime(nn) = stanam(j,i)
	      ii(nn) = qc_obs(2,j,i)
	      jj(nn) = qc_obs(3,j,i)
            endif
          enddo
	
          title = 'temperature background verification of tmf = '//tmtag//' (deg c)'
          ea = 1.50*5.0/9.0
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
 
          title = 'temperature verification of tmf = '//tmtag//' (deg c)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo

      case ("wndu")

        !time loop
        do kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

          title = 'u wind component background verification of tmf = '//tmtag//' (m/s)'
          ea = 2.00*knt2ms
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'u wind component verification of tmf = '//tmtag//' (m/s)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo

      case ("wndv")

        !time loop
        do kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

          title = 'v wind component background verification of tmf = '//tmtag//' (m/s)'
          ea = 2.00*knt2ms
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'v wind component verification of tmf = '//tmtag//' (m/s)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo

      case ("dewp")

        !time loop
        do kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

          title = 'dew point background verification of tmf = '//tmtag//' (deg c)'
          ea = 2.00*5.0/9.0
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
 
          title = 'dew point verification of tmf = '//tmtag//' (deg c)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo

      case ("tgd ")

        !time loop
        do kt = 1,numtmf
          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

          title = 'tgd background verification of tmf = '//tmtag//' (deg c)'
          ea = 2.00*5.0/9.0
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
 
          title = 'tdg verification of tmf = '//tmtag//' (deg c)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo

      case ("redp")
	!convert pascal to mb
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)/mb2pas

        !time loop
        do kt = 1,numtmf
 	  !convert pascal to mb
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

          title = 'reduced pressure background verification of tmf = '//tmtag//' (mb)'
          ea = 0.68
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'reduced pressure verification of tmf = '//tmtag//' (mb)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo


      case ("sfcp")
	!convert pascal to mb
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)/mb2pas

        !time loop
        do kt = 1,numtmf
	  !convert pascal to mb
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

          title = 'sfc pressure background verification of tmf = '//tmtag//' (mb)'
          ea = 0.68
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'sfc pressure verification of tmf = '//tmtag//' (mb)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo


      case ("mslp")
	!convert pascal to mb
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)/mb2pas

        !time loop
        do kt = 1,numtmf
	  !convert pascal to mb
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)/mb2pas

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

          title = 'msl pressure background verification of tmf = '//tmtag//' (mb)'
          ea = 0.68
          call verify(bkgrnd(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)

          title = 'msl pressure verification of tmf = '//tmtag//' (mb)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo


      case ("pcp1")
	!convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        do kt = 1,numtmf
	  !convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

	  !on precipitation we only verify obs and analysis
          title = '1hr precipitation verification of tmf = '//tmtag//' (meter)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo


      case ("pcp3")
	!convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        do kt = 1,numtmf
	  !convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

	  !on precipitation we only verify obs and analysis
          title = '3hr precipitation verification of tmf = '//tmtag//' (meter)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo


      case ("pcp6")
	!convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        do kt = 1,numtmf
	  !convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

	  !on precipitation we only verify obs and analysis
          title = '6hr precipitation verification of tmf = '//tmtag//' (meter)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo


      case ("pc24")
	!convert meter to mm 
	qc_obs(1,1:numobs(i),i) = qc_obs(1,1:numobs(i),i)*1000.0

        !time loop
        do kt = 1,numtmf
	  !convert meter to mm
	  analys(1:numgrd(1),1:numgrd(2),kt,i) = analys(1:numgrd(1),1:numgrd(2),kt,i)*1000.0
	  bkgrnd(1:numgrd(1),1:numgrd(2),kt,i) = bkgrnd(1:numgrd(1),1:numgrd(2),kt,i)*1000.0

          write(tmtag,'(i1)') kt
          nn= 0
          obstime = domain(1,3)+(kt-1)*lapsdt
          do j=1,numobs(i)
            if(int(qc_obs(4,j,i)).eq.obstime) then
              nn = nn + 1
              obsofthistime(nn) = qc_obs(1,j,i)
              staofthistime(nn) = stanam(j,i)
              ii(nn) = qc_obs(2,j,i)
              jj(nn) = qc_obs(3,j,i)
            endif
          enddo

	  !on precipitation we only verify obs and analysis
          title = '24hr precipitation verification of tmf = '//tmtag//' (meter)'
          call verify(analys(1:numgrd(1),1:numgrd(2),kt,i),obsofthistime(1:nn),		&
      	    	    staofthistime(1:nn),nn,title,iunit,					&
                    numgrd(1),numgrd(2),mxstts,x1a,x2a,y2a,ii,jj,ea,badsfc)
	enddo


      case default
        write(*,22) varnam(i)

    end select
  enddo	! end of variable loop

  close(iunit)

  !remember to free memory of bkgobs
  !
  deallocate(bkgobs,stat=err)
  if (err .ne. 0) then
    print*,'stmas>stmasver: cannot deallocate bkgobs memory!'
    stop
  endif
  deallocate(stanam,stat=err)
  if (err .ne. 0) then
    print*,'stmas>stmasver: cannot deallocate stanam memory!'
    stop
  endif


  !everything is done.
  print *,' normal completion of stmasver'
  return 

999  print *,'error opening ',ver_file(1:len)
  return

22 format('stmas>stmasver: do not know to verify such var: ',a4)

end subroutine stmasver

subroutine addbkgrd
!==========================================================
!  this routine mask out obs covered area defined by radius 
!  and then add background grid data into obs vector as if 
!  we have obs in those areas. 
!
!  history:
!       creation: 26-08-2008 by min-ken hsieh
!       modification:
!                 11-2008 by min-ken hsieh
!                            using jb term in stmasana instead of adding bkg to obs here.
!                            we only mark uncovered areas here.
!==========================================================

  implicit none

  !local variable
  integer :: i,j,kx,ky,kt,nn
  integer :: firstcoveredgrid(2),lastcoveredgrid(2)
  integer :: obstime
  logical :: land(numgrd(1),numgrd(2))		!land mask
  logical :: sameasstn(numgrd(1),numgrd(2))	!grid land/sea is the same as stn
  logical :: stnoverland

  character*1 :: tmtag

  !initialize land array
  do j = 1,numgrd(2)
    do i = 1,numgrd(1)
      land(i,j) = (lndfac(i,j).gt.0.0)
    enddo
  enddo

  !variable loop
  do i=1,numvar
    ! find out areas have been covered by obs for each time frame
    do kt=1,numtmf
      uncovr(:,:,kt,i) = .true.
      nn= 0
      obstime = domain(1,3)+(kt-1)*lapsdt
      do j=1,numobs(i)
        if(int(qc_obs(4,j,i)).eq.obstime) then
          nn = nn + 1
          firstcoveredgrid(1:2) = max0(1,floor(qc_obs(2:3,j,i))-radius(i))
          lastcoveredgrid(1:2) = min0(numgrd(1:2),floor(qc_obs(2:3,j,i))+radius(i)+1)

          !land/water process
          if(lndsea(i).eq.1) then
            stnoverland = (lndfac(int(qc_obs(2,j,i)),int(qc_obs(3,j,i))).gt.0.0)
            do ky=firstcoveredgrid(2),lastcoveredgrid(2)
	      do kx=firstcoveredgrid(1),lastcoveredgrid(1)
                sameasstn(kx,ky)= (stnoverland .eqv. land(kx,ky))
              enddo
            enddo
          else
            sameasstn = .true.
          endif
    
          ! mask out
          do ky=firstcoveredgrid(2),lastcoveredgrid(2)
            do kx=firstcoveredgrid(1),lastcoveredgrid(1)
	      uncovr(kx,ky,kt,i) = .not.sameasstn(kx,ky)
            enddo
          enddo
 
        endif
      enddo

    enddo

  enddo ! end of variable loop
  return

end subroutine addbkgrd

subroutine jbgridpt

!==========================================================
!  this routine identifies those gridpoints over which the
!  background are needed for j_b term in the cost function.
!  note. 
!  (1) this intends to improve the efficient of min-ken's
!  addbkgrd routine which takes 5 minutes over conus domain.
!  (2) consider spatial only, i.e., no temporal variation
!  as the namelist, stmas_mg.vr, now provides spatial radius
!  of data coverage.
!
!  history:
!       creation: nov, 2008 by yuanfu xie.
!==========================================================

  ! local variables:
  integer :: i,j,ic,jc,kc,iv,jo,is,ie,js,je,km,kp
  logical :: o,g
  real    :: r2,ol

  uncovr = .true.

  !variable loop:
  do iv=1,numvar
    ! time frames:
    r2 = radius(iv)*radius(iv)
    do jo=1,numobs(iv)

      ! center gridpoint:
      ic = int(qc_obs(2,jo,iv))
      jc = int(qc_obs(3,jo,iv))
      kc = (int(qc_obs(4,jo,iv))-domain(1,3))/lapsdt+1
      km = max0(kc-1,1)
      kp = min0(kc+1,numgrd(3))
      if ((verbal .eq. 1) .and. ((kc .lt. 1) .or. (kc .gt. numgrd(3)))) then
        print*,'obs out of range: ',iv,jo,qc_obs(4,jo,iv),domain(1,3)
        cycle
      endif
      if (kc .eq. 1) kc = 2
      if (kc .eq. numgrd(3)) kc = kc-1

      ! landfactor:
      ol = lndfac(ic,jc)+lndfac(ic+1,jc)+lndfac(ic,jc+1)+lndfac(ic+1,jc+1)
      o = .false.
      if ((ol .gt. 0.0) .or. (lndsea(iv) .eq. 0)) o = .true.

      ! covered circle:
      is = max0(ic-radius(iv),1)
      ie = min0(ic+radius(iv),numgrd(1))
      js = max0(jc-radius(iv),1)
      je = min0(jc+radius(iv),numgrd(2))

      do j=js,je
        do i=is,ie
          g = .false.
          if ((lndfac(i,j) .gt. 0.0) .or. (lndsea(iv) .eq. 0)) g = .true.
          ! within the influence radius and landfactors the same:
          if (((i-ic)*(i-ic)+(j-jc)*(j-jc) .lt. r2) .and. (o .eqv. g)) &
            uncovr(i,j,km:kp,iv) = .false.
        enddo
      enddo

    enddo
  enddo
end subroutine jbgridpt



subroutine jbingaus

!==========================================================
!  this routine identifies those gridpoints over which the
!  background are needed for j_b term in the cost function.
!  note. 
!  this is a modified version of jbgridpt.
!
!  history:
!       creation: nov, 2008 by yuanfu xie.
!==========================================================

  ! local variables:
  integer :: i,j,ic,jc,kc,iv,jo,is,ie,js,je,iex
  logical :: o,g
  real    :: r2,ol

  diagnl = 1.0

  iex = 10		! uncover gridpoint with gaussian decay

  !variable loop:
  do iv=1,numvar
    ! time frames:
    r2 = radius(iv)*radius(iv)
    do jo=1,numobs(iv)

      ! center gridpoint:
      ic = int(qc_obs(2,jo,iv))
      jc = int(qc_obs(3,jo,iv))
      kc = (int(qc_obs(4,jo,iv))-domain(1,3))/lapsdt+1
      if ((kc .lt. 1) .or. (kc .gt. numgrd(3))) then
        print*,'obs out of range: ',iv,jo,qc_obs(4,jo,iv),domain(1,3)
        cycle
      endif
      if (kc .eq. 1) kc = 2
      if (kc .eq. numgrd(3)) kc = kc-1

      ! landfactor:
      ol = lndfac(ic,jc)*lndfac(ic+1,jc)*lndfac(ic,jc+1)*lndfac(ic+1,jc+1)
      o = .false.
      if ((ol .gt. 0.0) .or. (lndsea(iv) .eq. 0)) o = .true.

      ! covered circle:
      is = max0(ic-radius(iv)-iex,1)
      ie = min0(ic+radius(iv)+iex,numgrd(1))
      js = max0(jc-radius(iv)-iex,1)
      je = min0(jc+radius(iv)+iex,numgrd(2))

      do j=js,je
        do i=is,ie
          g = .false.
          if ((lndfac(i,j) .gt. 0.0) .or. (lndsea(iv) .eq. 0)) g = .true.
          ! within the influence radius and landfactors the same:
          if (((i-ic)*(i-ic)+(j-jc)*(j-jc) .lt. r2) .and. (o .eqv. g)) then
            diagnl(i,j,iv) = 0.0
          else
            if (o .eqv. g) diagnl(i,j,iv) = &
              1.0-exp(r2-(i-ic)*(i-ic)-(j-jc)*(j-jc))
          endif
        enddo
      enddo

    enddo
  enddo
end subroutine jbingaus

end module lapsdatsrc
