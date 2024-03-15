!dis forecast systems laboratory
!dis noaa/oar/erl/fsl
!dis 325 broadway
!dis boulder, co 80303
!dis
!dis forecast research division
!dis local analysis and prediction branch
!dis laps
!dis
!dis this software and its documentation are in the public domain and
!dis are furnished "as is." the united states government, its
!dis instrumentalities, officers, employees, and agents make no
!dis warranty, express or implied, as to the usefulness of the software
!dis and documentation for any purpose. they assume no responsibility
!dis (1) for the use of the software and documentation; or (2) to provide
!dis technical support to users.
!dis
!dis permission to use, copy, modify, and distribute this software is
!dis hereby granted, provided that the entire disclaimer notice appears
!dis in all copies. all modifications to this software must be clearly
!dis documented, and are solely the responsibility of the agent making
!dis the modifications. if significant modifications or enhancements
!dis are made to this software, the fsl software policy manager
!dis (softwaremgr@fsl.noaa.gov) should be notified.
!dis

module prepostprc

!==========================================================
! this module sets up stmas pre/post processes.
!
! history:
! creation: yuanfu xie 8-2005
!==========================================================

  use definition

contains

subroutine prpstnls

!==========================================================
! this routine reads in namelists for stmas multigrid.
!
! history:
! creation: yuanfu xie 8-2005
!==========================================================

  implicit none

  character*200 :: dir
  integer :: nam,ios

  ! get namelist for stmas:
  call get_directory('static',dirstc,dirlen)
  dir = dirstc(1:dirlen)//'stmas_mg.nl'
  open(unit=11,file=dir(1:dirlen+12),form='formatted')
  read(11,nml=stmas,iostat=ios)
  close(11)
  ! check the numtmf consistent to numgrd
  if (numtmf .ne. numgrd(3)) then
    print*,'stmas>prpstnls: error in stmas_mg.nl, numtmf does not match numgrd'
    stop
  endif

  ! check the limit for number of variables to be analyzed:
  if (numvar .gt. maxvar) then
    write(*,*) 'stmas>prpstnls: error: too many variables!'
    write(*,*) 'stmas>prpstnls: increase maxvar and rerun!'
    stop
  endif

  ! get names of analyzing variables:
  dir = dirstc(1:dirlen)//'stmas_mg.vr'
  open(unit=11,file=dir(1:dirlen+12),form='formatted')
  do nam=1,numvar
    read(11,*) varnam(nam),thresh(nam),needbk(nam),bounds(nam),radius(nam),pnlt_v(nam),lndsea(nam),slevel(nam)
  enddo
  close(11)

end subroutine prpstnls

subroutine prpstlsx

!==========================================================
! this routine writes out the analyses into lsx gridded
! data in netcdf format.
!
! history:
! creation: yuanfu xie 6-2005
! modified: yuanfu xie 6-2006: write out more time frame.
! modified: yuanfu xie 6-2006: add sfc pressure and change
!			       the equations for theta and
!			       thetae from redp to sfcp.
!==========================================================

  use definition

  implicit none

  ! local variables:
  character :: ext*3,dir*200
  character*125 :: cmt(lsxvar)
  character*3 :: vnm(lsxvar)		! variable names
  character*3 :: vun(lsxvar)		! units
  character*3 :: crd(lsxvar)		! coordinates

  integer :: i,j,ncm
  integer :: ngd(2)	! actual grid points
  integer :: lvl(lsxvar) ! number of levels of each field
  integer :: itm	! time frame to write out
  integer :: i4t	! i4 time corresponding to itm
  integer :: idx(maxvar)! indices for derived variables;
  integer :: iwv	! index of v; 
  integer :: nvr,mvr	! number of variables;
  integer :: len
  integer :: sts	! return status

  real, external :: ept
  real :: tmp,dew,td_c,pr_m
  real :: gdx(numgrd(1),numgrd(2))	! grid spacing (x)
  real :: gdy(numgrd(1),numgrd(2))	! grid spacing (y)
  real :: gdp(numgrd(1),numgrd(2))	! pressure in mb
  real :: gdt(numgrd(1),numgrd(2))	! theta
  real :: dum(numgrd(1),numgrd(2))	! unused value
  real :: dmm(numgrd(1),numgrd(2))	! unused value
  real :: mrc(numgrd(1),numgrd(2))
  real :: dat(numgrd(1)-2*numfic(1),numgrd(2)-2*numfic(2),lsxvar)

  ! added by min-ken hsieh
  ! declare l1s array to store 1 hr and 24 hr precip.
  ! then output to l1s file via write_laps_data
  real :: l1s(numgrd(1)-2*numfic(1),numgrd(2)-2*numfic(2),2)
  character*125 :: l1s_cmt(2)
  character*3 :: l1s_vnm(2)		! variable names
  character*3 :: l1s_vun(2)		! units
  character*3 :: l1s_crd(2)		! coordinates
  character*3 :: lwm_vnm(2)		! variable names

  integer :: l1s_lvl(2)			! number of levels of each field

  ! mixing ratio function from td and p:
  real :: wmr

  ! time frame to write out:
  !do itm = numgrd(3)-numfic(3),max0(1,numgrd(3)-numfic(3)-2),-1	! time frame
  ! yuanfu: change output to current time frame only for forecast:
  do itm = numgrd(3)-numfic(3),max0(1,numgrd(3)-numfic(3)),-1	! time frame

  i4t = i4time-(numgrd(3)-numfic(3)-itm)*lapsdt	! i4time corresponding to itm
  mvr = lsxvar
  lvl = 0
  crd = 'agl'

  ! number of gridpoints without fictitous points:
  ngd = numgrd(1:2)-2*numfic(1:2)
  gdx = phydxy(1)
  gdy = phydxy(2)

  ! parsing the variable names:
  nvr = 0
  do i=1,numvar

    select case (varnam(i))
    case ("temp")	! temperature
      nvr = nvr+1
      if (nvr .gt. lsxvar) then
	write(*,2)
        stop
      endif
      vnm(nvr) = 't  '
      vun(nvr) = 'k  '
      cmt(nvr) = 'surface temperature'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("visb")	! visibility
      nvr = nvr+1
      if (nvr .gt. lsxvar) then
	write(*,2)
        stop
      endif
      vnm(nvr) = 'vis'
      vun(nvr) = 'm  '
      cmt(nvr) = 'visibilty'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("ceil")	! cloud ceiling
      nvr = nvr+1
      if (nvr .gt. lsxvar) then
	write(*,2)
        stop
      endif
      vnm(nvr) = 'cc'
      vun(nvr) = 'm  '
      cmt(nvr) = 'cloud ceiling'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("dewp")	! dew point temperature
      nvr = nvr+1
      if (nvr .gt. lsxvar) then
	write(*,2)
        stop
      endif
      vnm(nvr) = 'td '
      vun(nvr) = 'k  '
      cmt(nvr) = 'dewpoint'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("tgd ")	! ground temperature
      nvr = nvr+1
      if (nvr .gt. lsxvar) then
	write(*,2)
        stop
      endif
      vnm(nvr) = 'tgd'
      vun(nvr) = 'k  '
      cmt(nvr) = 'ground temperature'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("mslp")	! mean sea level pressure
      nvr = nvr+1
      if (nvr .gt. lsxvar) then
	write(*,2)
        stop
      endif
      vnm(nvr) = 'msl'
      vun(nvr) = 'pa '
      cmt(nvr) = 'msl pressure'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("redp")	! reduced pressure
      if (nvr+2 .gt. lsxvar) then
	write(*,2)
        stop
      endif
      nvr = nvr+1
      vnm(nvr) = 'p  '
      vun(nvr) = 'pa '
      write(cmt(nvr),11) int(rdplvl),' m reduced pressure'
 11   format(i5,a19)
      ! cmt(nvr) = '0  m reduced pressure'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
      if (press_pert .eq. 1) then
        nvr = nvr+1
        vnm(nvr) = 'pp '
        vun(nvr) = 'pa '
        cmt(nvr) = '0  m reduced pressure change'
        call preschng(dat(1,1,nvr-1),ngd,dat(1,1,nvr))
      endif
    case ("sfcp")	! surface pressure
      if (nvr+1 .gt. lsxvar) then
	write(*,2)
        stop
      endif
      nvr = nvr+1
      vnm(nvr) = 'ps '
      vun(nvr) = 'pa '
      cmt(nvr) = 'surface pressure'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("wndu") 	! u wind binded with wndv
      iwv = 0
      do j=1,numvar
        if (varnam(j) .eq. 'wndv') iwv = j
      enddo
      if (iwv .eq. 0) then
	write(*,3)
        stop
      endif
      if (nvr+4 .gt. lsxvar) then
	write(*,2)
        stop
      endif
      nvr = nvr+1
      vnm(nvr) = 'u  '
      vun(nvr) = 'm/s'
      cmt(nvr) = 'u component'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,i)
      nvr = nvr+1
      vnm(nvr) = 'v  '
      vun(nvr) = 'm/s'
      cmt(nvr) = 'v component'
      dat(1:ngd(1),1:ngd(2),nvr) = &
	analys(numfic(1)+1:numgrd(1)-numfic(1), &
	       numfic(2)+1:numgrd(2)-numfic(2),itm,iwv)

      ! write surface wind to lwm for plot purpose:
      ! lwm_vnm(1) = 'su '
      ! lwm_vnm(2) = 'sv '
      ! call put_laps_multi_2d(i4t,'lwm',lwm_vnm,vun(nvr-1),cmt(nvr-1), &
      !                        dat(1,1,nvr-1),ngd(1),ngd(2),2,sts)

      nvr = nvr+1
      vnm(nvr) = 'vor'
      vun(nvr) = '/s '
      cmt(nvr) = 'vorticity'
      ! interior:
      dat(2:ngd(1)-1,2:ngd(2)-1,nvr) = 0.5/phydxy(2)*( &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+1:numgrd(2)-numfic(2)-2,itm,i)- &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+3:numgrd(2)-numfic(2)  ,itm,i))+ &
		 		       0.5/phydxy(1)*( &
	analys(numfic(1)+3:numgrd(1)-numfic(1)  , &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,iwv)- &
	analys(numfic(1)+1:numgrd(1)-numfic(1)-2, &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,iwv))
      ! extrapolate boundary values:
      call extraplt(dat(1,1,nvr),ngd)
      nvr = nvr+1
      vnm(nvr) = 'div'
      vun(nvr) = '/s '
      cmt(nvr) = 'divergence'
      ! interior:
      dat(2:ngd(1)-1,2:ngd(2)-1,nvr) = -0.5/phydxy(2)*( &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+1:numgrd(2)-numfic(2)-2,itm,iwv)- &
	analys(numfic(1)+2:numgrd(1)-numfic(1)-1, &
	       numfic(2)+3:numgrd(2)-numfic(2)  ,itm,iwv))+ &
		 		       0.5/phydxy(1)*( &
	analys(numfic(1)+3:numgrd(1)-numfic(1)  , &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,i)- &
	analys(numfic(1)+1:numgrd(1)-numfic(1)-2, &
	       numfic(2)+2:numgrd(2)-numfic(2)-1,itm,i))
      ! extrapolate boundary values:
      call extraplt(dat(1,1,nvr),ngd)
    case ("wndv")	! v wind dealt with u
    case ("pcp1")	! 1 hour precip accumulation		added by min-ken hsieh
      ! make a copy in l1s array
      l1s(1:ngd(1),1:ngd(2),1) = &
        analys(numfic(1)+1:numgrd(1)-numfic(1), &
               numfic(2)+1:numgrd(2)-numfic(2),itm,i)

    case ("pc24")	! 24 hour precip accumulation		added by min-ken hsieh
      ! make a copy in l1s array
      l1s(1:ngd(1),1:ngd(2),2) = &
        analys(numfic(1)+1:numgrd(1)-numfic(1), &
               numfic(2)+1:numgrd(2)-numfic(2),itm,i)
    case ("pcp3")
    case ("pcp6")
    case default
      write(*,1) varnam(i)
    end select
  enddo

1 format('prpstprc>prpstlsx: no such variable to write: ',a4)
2 format('prpstprc>prpstlsx: too much variables to write; aborts')
3 format('prpstprc>prpstlsx: v wind is missing')

  !======================
  ! other derived fields:
  !======================
    
  ! theta and theta_e:	potential/equivalent potential temperature:
  ! search necessary basic variables:
  ncm = 0
  idx = 0
  do j=1,numvar
    if (varnam(j) .eq. "temp") then
      ncm = ncm+1
      idx(1) = j
    endif
    if (varnam(j) .eq. "dewp") then
      ncm = ncm+1
      idx(2) = j
    endif
    if (varnam(j) .eq. "sfcp") then
      ncm = ncm+1
      idx(3) = j
    endif
  enddo
  ! found necessary variables for potential temperature:
  if ((idx(1) .ne. 0) .and. (idx(3) .ne. 0)) then
    nvr = nvr+1
    if (nvr .gt. lsxvar) then
      write(*,2)
      stop
    endif
    vnm(nvr) = 'th '
    vun(nvr) = 'k  '
    cmt(nvr) = 'potential temperature'
    ! theta:
    gdt(1:numgrd(1),1:numgrd(2)) = &
      analys(1:numgrd(1),1:numgrd(2),itm,idx(1))* &
	(1.0e5/analys(1:numgrd(1),1:numgrd(2),itm,idx(3))) &
	**(gascnt/spheat)
    ! remove fictitous point and convert to r4:
    dat(1:ngd(1),1:ngd(2),nvr) = &
      gdt(numfic(1)+1:numgrd(1)-numfic(1), &
	  numfic(2)+1:numgrd(2)-numfic(2))
  endif
  ! found necessary variables for mixing ratio:
  if ((idx(2) .ne. 0) .and. (idx(3) .ne. 0)) then
    nvr = nvr+1
    if (nvr .gt. lsxvar) then
      write(*,2)
      stop
    endif
    vnm(nvr) = 'mr '
    vun(nvr) = 'g/kg'
    cmt(nvr) = 'mixing ratio'
    ! mixing ratio from td and p:
    do j=1,ngd(2)
      do i=1,ngd(1)
        td_c = analys(i,j,itm,idx(2))-temp_0	! dewpoint in celsius
        pr_m = analys(i,j,itm,idx(3))/100.0	! pressure in mb
        dat(i,j,nvr) = wmr(pr_m,td_c)
      enddo
    enddo
  endif
  ! found necessary variables for relative humidity:
  if ((idx(1) .ne. 0) .and. (idx(2) .ne. 0)) then
    nvr = nvr+1
    if (nvr .gt. lsxvar) then
      write(*,2)
      stop
    endif
    vnm(nvr) = 'rh '
    vun(nvr) = 'm  '
    cmt(nvr) = 'relative humidity'
    ! rh from t and td: hum returns rh in [0,1]
    call hum(analys(1,1,itm,idx(1)),analys(1,1,itm,idx(2)), &
              dat(1,1,nvr),numgrd(1),numgrd(2),dum,dmm)
    dat(1:ngd(1),1:ngd(2),nvr) = dat(1:ngd(1),1:ngd(2),nvr)*100		! rh in [0,100]
  endif
  ! found necessary variables for equivalent potential temperature:
  if (ncm .eq. 3) then
    nvr = nvr+1
    if (nvr .gt. lsxvar) then
      write(*,2)
      stop
    endif
    vnm(nvr) = 'the'
    vun(nvr) = 'k  '
    cmt(nvr) = 'equivalent potential temperature'
    ! use laps function ept to compute theta_e:
    gdp = analys(1:numgrd(1),1:numgrd(2),itm,idx(3))/mb2pas! pressure(mb)
    do j=1,ngd(2)
      do i=1,ngd(1)
	tmp = analys(numfic(1)+i,numfic(2)+j,itm,idx(1))-temp_0
	dew = analys(numfic(1)+i,numfic(2)+j,itm,idx(2))-temp_0
        dat(i,j,nvr) = ept(tmp,dew,gdp(i+numfic(1),j+numfic(2)))+temp_0
      enddo
    enddo
  endif

  ! moisture convergence:
  ! search necessary basic variables:
  ncm = 0
  idx = 0
  do j=1,numvar
    if (varnam(j) .eq. "temp") then
      ncm = ncm+1
      idx(1) = j
    endif
    if (varnam(j) .eq. "dewp") then
      ncm = ncm+1
      idx(2) = j
    endif
    if (varnam(j) .eq. "redp") then
      ncm = ncm+1
      idx(3) = j
    endif
    if (varnam(j) .eq. "wndu") then
      ncm = ncm+1
      idx(4) = j
    endif
    if (varnam(j) .eq. "wndv") then
      ncm = ncm+1
      idx(5) = j
    endif
  enddo
  ! found necessary variables:
  if (ncm .eq. 5) then
    nvr = nvr+1
    if (nvr .gt. lsxvar) then
      write(*,2)
      stop
    endif
    vnm(nvr) = 'mrc'
    vun(nvr) = '/s '
    cmt(nvr) = 'moisture convergence'
    call meso_anl(analys(1,1,itm,idx(4)),analys(1,1,itm,idx(5)), &
      gdp,analys(1,1,itm,idx(1)),analys(1,1,itm,idx(2)), &
      gdt,gdx,gdy,dum,mrc,dum,dum,dum,numgrd(1),numgrd(2))
    dat(1:ngd(1),1:ngd(2),nvr) = &
      mrc(numfic(1)+1:numgrd(1)-numfic(1), &
	  numfic(2)+1:numgrd(2)-numfic(2))
  endif

  !====================
  !  write to lsx file:
  !====================

  ! get the directory for lsx file:
  ext = 'lsx'
  call get_directory(ext,dir,len)

  ! write data to a lsx file:
  call write_laps_data(i4t,dir,ext,ngd(1),ngd(2),nvr,nvr, &
     		       vnm,lvl,crd,vun,cmt,dat,sts)

  ! write data to a lsx file under balance:
  dir(len-3:len+8) = 'balance/lsx/'
  call write_laps_data(i4t,dir(1:len+8),ext,ngd(1),ngd(2),nvr,nvr, &
     		       vnm,lvl,crd,vun,cmt,dat,sts)

  !====================
  !  write to l1s file:
  !====================

  ! get the directory for l1s file:
  ext = 'l1s'
  call get_directory(ext,dir,len)
  l1s_cmt(1) = 'laps 60 minute snow accumulation'
  l1s_cmt(2) = 'storm total precip. accum.'
  l1s_vnm(1) = 'r01'
  l1s_vnm(2) = 'rto'
  l1s_vun = 'm'
  l1s_crd = 'msl'
  l1s_lvl = 0

  ! write data to a lsx file:
  ! temporarily turn off l1s output as soil analysis needs snow as well
  ! however, mile's work is only on rain.
  ! call write_laps_data(i4t,dir,ext,ngd(1),ngd(2),2,2, &
  !   		       l1s_vnm,l1s_lvl,l1s_crd,l1s_vun,l1s_cmt,l1s,sts)
  ! end of writing a series of time frames:
  enddo

end subroutine prpstlsx

subroutine preschng(pres,ngrd,chng)

!==========================================================
!  this routine computes pressure change from pressure.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(2)
  real, intent(in) :: pres(ngrd(1),ngrd(2))
  real, intent(out) :: chng(ngrd(1),ngrd(2))

  ! local variables:
  integer :: i,j,nbi,nbj,nbs,ips
  real :: gam,kap,val,wgt,wsm
  real :: prs(ngrd(1),ngrd(2))

  ! barnes parameters:
  gam = 0.2
  kap = 2.0e3/gam	! kapa_0
  nbs = 60		! neighbors to analyzed
  
  ! initial:
  chng = 0.0
  prs = 0.0

  ! barnes analysis:
  do ips=1,1

    kap = kap*gam	! adjust kappa

    ! every gridpoint:
    do j=1,ngrd(2)
      do i=1,ngrd(1)
	
	! for all considered neighbors: every other one to save time
	val = 0.0
	wsm = 0.0
	do nbj=-nbs,nbs,4
	  do nbi=-nbs,nbs,4

	    wgt = exp(-(float(nbi)**2+float(nbj)**2)/kap)
	    if ((i+nbi .ge. 1) .and. (i+nbi .le. ngrd(1)) .and. &
	        (j+nbj .ge. 1) .and. (j+nbj .le. ngrd(2))) then
	      val = val+(pres(i+nbi,j+nbj)-prs(i+nbi,j+nbj))*wgt
	      wsm = wsm + wgt
	    endif
	  
	  enddo
	enddo

        ! update grid value:
        chng(i,j) = chng(i,j)+val/wsm

      enddo
    enddo

    ! save iterated result:
    prs = chng

  enddo

  ! compute pressure changes:
  chng = pres-prs

end subroutine preschng

subroutine extraplt(grid,ngrd)

!==========================================================
!  this routine extrapolates interior points to its boundary
!  values assuming ngrd > 3.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(2)
  real, intent(inout) :: grid(ngrd(1),ngrd(2))

  ! x:
  grid(1,2:ngrd(2)-1) = &
    2.0*grid(2,2:ngrd(2)-1)-grid(3,2:ngrd(2)-1)
  grid(ngrd(1),2:ngrd(2)-1) = &
    2.0*grid(ngrd(1)-1,2:ngrd(2)-1)-grid(ngrd(1)-2,2:ngrd(2)-1)

  ! y:
  grid(1:ngrd(1),1) = &
    2.0*grid(1:ngrd(1),2)-grid(1:ngrd(1),3)
  grid(1:ngrd(1),ngrd(2)) = &
    2.0*grid(1:ngrd(1),ngrd(2)-1)-grid(1:ngrd(1),ngrd(2)-2)

end subroutine extraplt

end module prepostprc
