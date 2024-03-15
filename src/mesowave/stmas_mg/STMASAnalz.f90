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

module stmasanalz

!==========================================================
!  this module contains stmas variational analysis based
!  a multigrid fitting technique.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!       modification:
!                 10-2008 by min-ken hsieh
!                 11-2008 by min-ken hsieh
!                 stmasana
!                 functn3d
!                 gradnt3d
!==========================================================

  use definition

contains

subroutine stmasana(anal,ngrd,dxyt,domn,bkgd,nfrm, &
		    obsv,nobs,wght,stna,ospc,indx, &
		    coef, bund,ipar,rpar,vnam,pnlt,&
                    slvl,ucvr,diag)
		    

!==========================================================
!  this routine reads into necessary namelists for stmas
!  analysis through a multigrid technique.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!       modification:
!                 08-2008 by min-ken hsieh
!                 add parameter stna to map each obs its stn name for stmasver
!       modification:
!                 10-2008 by min-ken hsieh
!                 bound option only apply to last level
!                 pass in penalty for each var
!                 pass in slevel for each var
!       modification:
!                 11-2008 by min-ken hsieh
!                 pass in uncovr for each var
!                 call intplbkg to interpolate bkg to multigrid
!                 and pass interpolated bkg and uncover (hbg/huc) to minimize
!       modification:
!                 12-2008 by yuanfu xie
!                 pass in diag(nol) array for j_b term.
!       modification:
!                 01-2009 by yuanfu xie
!                 calculate number of multigrid levels using slvl
!                 assuming level 1 with gridpoint 3 3 3.
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(3)	! grid numbers
  integer, intent(in) :: nfrm		! bkgd time frames
  integer, intent(inout) :: nobs	! number of obs
  integer, intent(in) :: indx(6,nobs)	! indices (intepolate)
  integer, intent(in) :: bund		! bound constraints
  integer, intent(in) :: ipar(1)	! integer parameter
  integer, intent(in) :: slvl           ! level to start analysis. added by min-ken
  real, intent(in) :: dxyt(3)		! grid spacing
  real, intent(in) :: domn(2,3)		! domain
  real, intent(in) :: bkgd(ngrd(1),ngrd(2),nfrm)
  real, intent(inout) :: obsv(4,nobs)
  real, intent(inout) :: wght(nobs)	! obs weightings
  real, intent(in) :: pnlt		! penalty for each variable
  character*20, intent(inout):: stna(nobs)
					! obs station name by min-ken hsieh
  character*4, intent(in):: vnam        ! variable name(used in lvl addition and output each level result)
  real, intent(in) :: ospc(3)		! obs spacing
  real, intent(in) :: coef(6,nobs)	! coeffients
  real, intent(in) :: rpar(1)		! real parameter
  real, intent(in) :: diag(ngrd(1),ngrd(2))! diagnol array for j_b

  real, intent(inout) :: anal(ngrd(1),ngrd(2),ngrd(3))

  logical, intent(in) :: ucvr(ngrd(1),ngrd(2),ngrd(3)) !uncovered array by min-ken

  ! local variables:
  integer :: lvl(3)			! number of multigrid levels
  integer :: lvc(3)			! account of the levels
  integer :: mgd(3)			! multigrid points
  integer :: mld(3)			! leading dimensions
  integer :: inc(3)			! grid change increment
  integer :: mcl			! number of multigrid cycles
  integer :: mlv			! maximum number levels
  integer :: idx(6,nobs)		! indices at a multigrid
  integer :: i,j,k,l,ier,ii,ij,ik
  real :: dis,rsz			! distance and resizes
  real :: dgd(3)			! grid spacing at a multigrid
  real :: coe(6,nobs)			! coefficients at a multigrid
  real, allocatable, dimension(:,:,:) &
       :: sln				! multigrid solutions
  integer, allocatable, dimension(:,:,:) &
       :: huc				! interpolated uncover by min-ken


  ! number of multigrid v-cycles:
  mcl = 0

  ! calculate total multigrid levels needed:
  do i=1,3
    lvl(i) = 0
    ! initial grid distance 3 gridpoints:
    dis = (domn(2,i)-domn(1,i))/2.0
1   continue
    if (dis .gt. ospc(i)) then
      lvl(i) = lvl(i)+1
      dis = dis*0.5
      goto 1
    endif
    ! lvl counts addition multigrid levels except the current 3x3x3 level
    lvl(i) = min0(lvl(i),int(alog(float(ngrd(i)-1))/alog(2.0))-1)
  enddo
  if (verbal .eq. 1) write(*,2) lvl(1:3)
2 format('stmasana: number of levels: ',3i3)

  ! leading dimensions according to the total levels:
  mld = 2**(lvl+1)+1

  ! allocate memory for mutligrid solutions:
  allocate(sln(mld(1),mld(2),mld(3)), stat=ier)
  
  ! allocate interpolated bkg and uncover
  allocate(huc(mld(1),mld(2),mld(3)), stat=ier)
  
  ! start multigrid analysis:

  ! modified by min-ken hsieh
  ! mgd should be determinated by slvl,
  ! but it cannot be bigger than mld
  do i=1,3
    mgd(i) = min0(2**(slvl-1)+1,mld(i))
  enddo
  ! maximum multigrid levels:
  mlv = maxval(lvl(1:3))-slvl+1

  ! multigrid cycles:
  do l=1,mcl*2+1

    ! down/up cycle:
    if (mod(l,2) .eq. 1) then
      rsz = 2.0
    else
      rsz = 0.5
    endif

    ! through all possible levels:
    ! modified by min-ken hsieh
    ! because we start from slvl
    ! lvc are no longer starting from 0
    !lvc = 0
    lvc = slvl - 1

    do k=1,mlv

      ! redefine number of gridpoints:
      do j=1,3
        if (lvc(j) .le. lvl(j)) then
	  mgd(j) = (mgd(j)-1)*rsz+1
	  inc(j) = 2			! grid change
	  lvc(j) = lvc(j)+1		! count levels
	else
	  inc(j) = 1			! grid unchange
	endif
      enddo

      ! interpolate a multigrid to observations:
      dgd = (domn(2,1:3)-domn(1,1:3))/float(mgd-1)
      call grid2obs(idx,coe,obsv,nobs,wght,stna,mgd,dgd,domn)

      ! added by min-ken hsieh
      ! interpolate a background to multigrid
      call intplbkg(ucvr,ngrd,huc,mgd,mld,dgd,dxyt,domn)

      ! initial guesses:
      if ((l .eq. 1) .and. (k .eq. 1)) &
	!sln = 0.5*(maxval(obsv(1,1:nobs))-minval(obsv(1,1:nobs)))
	sln = 0.0

      ! down and up cycle:
      if (mod(l,2) .eq. 0) then
	call projectn(sln,mld,mgd,inc)	! up cycle
      else
	call interpln(sln,mld,mgd,inc)	! down cycle
      endif

      ! added by min-ken hsieh
      ! smooth:
      !call smoother(sln,mld,mgd)

      ! analyzing:
      !bound the pcp analysis on finest level
      if (k.eq.mlv) then
        call minimize(sln,mld,mgd,obsv,nobs,idx,coe,wght,bund,pnlt,huc)
      else
        call minimize(sln,mld,mgd,obsv,nobs,idx,coe,wght,0,pnlt,huc)
      endif

    enddo
  enddo

  ! map the multigrid solution to the grid requested:
  call mul2grid(sln,mld,mgd,anal,ngrd,domn)

  ! deallocate multigrid:
  deallocate(sln,stat=ier)
  deallocate(huc,stat=ier)

end subroutine stmasana

subroutine smoother(sltn,mlds,ngrd)

!==========================================================
!  this routine smooth a grid with its neighbor values 
!
!  history:
!	creation: 9-2008 by min-ken hsieh
!==========================================================

  implicit none

  integer, intent(in) :: mlds(3),ngrd(3)
  real, intent(inout) :: sltn(mlds(1),mlds(2),mlds(3))


  ! smooth:
  ! x-direction
  sltn(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) = &
    0.5*sltn(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3))+ &
    0.5*sltn(3:ngrd(1),1:ngrd(2),1:ngrd(3)) 
  ! y-direction
  sltn(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) = &
    0.5*sltn(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3))+ &
    0.5*sltn(1:ngrd(1),3:ngrd(2),1:ngrd(3)) 
  ! t-direction
  sltn(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) = &
    0.5*sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2)+ &
    0.5*sltn(1:ngrd(1),1:ngrd(2),3:ngrd(3)) 

end subroutine smoother

subroutine projectn(sltn,mlds,ngrd,incr)

!==========================================================
!  this routine projects a grid value function to a coarse
!  resolution whose numbers of gridpoints are ngrd.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer, intent(in) :: mlds(3),incr(3),ngrd(3)
  real, intent(inout) :: sltn(mlds(1),mlds(2),mlds(3))

  ! local variables:
  integer :: mgd(3)

  mgd = incr*(ngrd-1)+1

  ! projection based grid increments:
  sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)) = &
    sltn(1:mgd(1):incr(1),1:mgd(2):incr(2),1:mgd(3):incr(3))

end subroutine projectn

subroutine interpln(sltn,mlds,ngrd,incr)

!==========================================================
!  this routine interpolates a grid value function to a 
!  fine resolution with ngrd gridpoints.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer, intent(in) :: mlds(3),ngrd(3),incr(3)
  real, intent(inout) :: sltn(mlds(1),mlds(2),mlds(3))

  ! local variables:
  integer :: mgd(3)

  mgd = (ngrd-1)/incr+1

  ! projection based grid increments:
  sltn(1:ngrd(1):incr(1),1:ngrd(2):incr(2),1:ngrd(3):incr(3)) = &
    sltn(1:mgd(1),1:mgd(2),1:mgd(3))

  ! interpolate:
  
  ! x direction:
  if (incr(1) .eq. 2) &
    sltn(2:ngrd(1)-1:2,1:ngrd(2):incr(2),1:ngrd(3):incr(3)) = &
      0.5*( &
      sltn(1:ngrd(1)-2:2,1:ngrd(2):incr(2),1:ngrd(3):incr(3))+ &
      sltn(3:ngrd(1)  :2,1:ngrd(2):incr(2),1:ngrd(3):incr(3)) )

  ! y direction:
  if (incr(2) .eq. 2) &
    sltn(1:ngrd(1),2:ngrd(2)-1:2,1:ngrd(3):incr(3)) = 0.5*( &
      sltn(1:ngrd(1),1:ngrd(2)-2:2,1:ngrd(3):incr(3))+ &
      sltn(1:ngrd(1),3:ngrd(2)  :2,1:ngrd(3):incr(3)) )

  ! t direction:
  if (incr(3) .eq. 2) &
    sltn(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1:2) = 0.5*( &
      sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2:2)+ &
      sltn(1:ngrd(1),1:ngrd(2),3:ngrd(3)  :2) )

end subroutine interpln

subroutine mul2grid(sltn,mled,mgrd,anal,ngrd,domn)

!==========================================================
!  this routine projects a grid value function to a coarse
!  resolution.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer, intent(in) :: mled(3),mgrd(3),ngrd(3)
  real, intent(in) :: sltn(mled(1),mled(2),mled(3)),domn(2,3)
  real, intent(out) :: anal(ngrd(1),ngrd(2),ngrd(3))

  ! local variables:
  integer :: i,j,k,idx(6),ix,iy,it,ier
  real :: xyt(3),dsm(3),dsn(3),coe(6)

  ! grid spacings:
  dsm = (domn(2,1:3)-domn(1,1:3))/float(mgrd-1)
  dsn = (domn(2,1:3)-domn(1,1:3))/float(ngrd-1)

  ! interpolate:
  do k=1,ngrd(3)
    xyt(3) = domn(1,3)+(k-1)*dsn(3)
    do j=1,ngrd(2)
      xyt(2) = domn(1,2)+(j-1)*dsn(2)
      do i=1,ngrd(1)
	xyt(1) = domn(1,1)+(i-1)*dsn(1)

	call intplt3d(xyt,mgrd,dsm,domn,idx,coe,ier)

	! evaluate:
	anal(i,j,k) = 0.0
	do it=3,6,3
	  do iy=2,5,3
	    do ix=1,4,3
	      anal(i,j,k) = anal(i,j,k)+ &
		sltn(idx(ix),idx(iy),idx(it))* &
		coe(ix)*coe(iy)*coe(it)
	    enddo
	  enddo
	enddo

      enddo
    enddo
  enddo

end subroutine mul2grid

subroutine minimize(sltn,mlds,ngrd,obsv,nobs,indx,coef, &
		    wght,bund,pnlt,huc)

!==========================================================
!  this routine minimizes the stmas cost function.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!       modified: 10-2008 by min-ken hsieh.
!                 pass in penalty for each var
!       modified: 11-2008 by min-ken hsieh.
!                 pass in hbg,huc to calcualte jb
!       modified: 12-2013 by yuanfu xie:
!                 new lbfgsb 3.0 fixes machine eps calculation
!                 yuanfu adjust lbfgsb 3.0 for solving a super
!                 large minimization.
!                 a) switch to lbfgsb.3.0 with yuanfu's
!                    modification for super large minimization
!                 b) adjust wkspc dimension for new lbfgsb and
!                    add a new parameter pg for new lbfgsb
!                 c) change automatic arrays of wkspc, wrka, 
!                    bdlow, bdupp, nbund to allocatables
!==========================================================

  implicit none

  integer, intent(in) :: mlds(3),ngrd(3),nobs,indx(6,nobs)
  integer, intent(in) :: bund
  real, intent(in) :: obsv(4,nobs),coef(6,nobs),wght(nobs)
  real, intent(inout) :: sltn(mlds(1),mlds(2),mlds(3))
  real, intent(in) :: pnlt	!penalty for each variable

  integer, intent(in) :: huc(mlds(1),mlds(2),mlds(3))	! interpolated uncover array by min-ken

  !** lbfgs_b variables:
  integer, parameter :: msave=7		! max iter save

  character*60 :: ctask,csave		! evaluation flag

  real, allocatable :: wkspc(:),bdlow(:),bdupp(:)
  real :: factr,pg,dsave(29)            ! yuanfu added pg for using setulb

  integer :: iprnt,isbmn,isave(44)
  integer, allocatable :: nbund(:) ! bound flags
  integer, allocatable :: iwrka(:)

  logical :: lsave(4)

  !** end of lbfgs_b declarations.

  ! local variables:
  integer :: itr,sts,nvr,i,j,k,istatus
  real :: fcn,fnp,fnm,ggg,eps
  real :: gdt(ngrd(1),ngrd(2),ngrd(3))	! gradients
  real :: grd(ngrd(1),ngrd(2),ngrd(3))  ! grid function
  real :: egd(ngrd(1),ngrd(2),ngrd(3))  ! grid function
  
  ! dec 2013: yuanfu changed the wkspc dimension from (2*msave+4)
  ! to (2*msave+5) for both old and new versions of lbfgsb:
  allocate(wkspc(ngrd(1)*ngrd(2)*ngrd(3)*(2*msave+5)+ &
                12*msave*msave+12*msave), &
           iwrka(3*ngrd(1)*ngrd(2)*ngrd(3)), stat=istatus)
  allocate(bdlow(ngrd(1)*ngrd(2)*ngrd(3)), &
           bdupp(ngrd(1)*ngrd(2)*ngrd(3)), &
           nbund(ngrd(1)*ngrd(2)*ngrd(3)), stat=istatus)
           

  ! start lbfgs_b:
  ctask = 'start'
  nvr = ngrd(1)*ngrd(2)*ngrd(3)		! number of controls
  factr = 1.0d2
  pg = 1.0e-4        ! yuanfu added for using setulb
  iprnt = 1
  isbmn = 1

  ! simple bound constraints for controls (only low bound used here):
  nbund = bund
  if (bund .eq. 1) bdlow = 0.0

  ! initial:
  itr = 0
  grd = sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3))

  ! iterations:
1 continue

  call setulb(nvr,msave,grd,bdlow,bdupp,nbund,fcn,gdt,factr, &
	      pg,wkspc,iwrka,ctask,iprnt,csave,lsave,isave,dsave)

  ! exit iteration if succeed:
  if (ctask(1:11) .eq. 'convergence') goto 2

  ! function and gradient values are needed:
  if (ctask(1:2) .eq. 'fg') then

    ! function value:
    call functn3d(fcn,grd,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

    ! gradient values:
    call gradnt3d(gdt,grd,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

    ! check gradients:
    ! i=int(0.5*ngrd(1))
    ! j=int(0.3*ngrd(2))
    ! k=3
    ! eps = 1.0e-1
    ! egd = grd
    ! egd(i,j,k) = egd(i,j,k)+eps
    ! call functn3d(fnp,egd,ngrd,obsv,nobs,indx,coef,wght)
    ! egd(i,j,k) = egd(i,j,k)-2.0*eps
    ! call functn3d(fnm,egd,ngrd,obsv,nobs,indx,coef,wght)
    ! ggg = (fnp-fnm)*0.5/eps
    ! write(*,9) ggg,gdt(i,j,k),i,j,k

  endif
9 format('stmas>stmasmin: gradient check: ',2e12.4,3i6)

  ! exit if irregularity is encountered:
  if ((ctask(1:2) .ne. 'fg') .and. (ctask(1:5) .ne. 'new_x')) then
    write(*,*) 'stmas>stmasmin: irregularity termination of lbfgsb'
    goto 2
  endif

  ! a new iteration point:
  if (ctask(1:5) .eq. 'new_x') itr = itr+1
  if (itr .lt. maxitr) goto 1

  ! exit of lbfgsb iteration:
2 continue

  ! save the solution:
  sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)) = grd

end subroutine minimize

subroutine functn3d(fctn,grid,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

!==========================================================
!  this routine evaluates the stmas cost function.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!       modified: 10-2008 by min-ken hsieh.
!                 pass in penalty for each var
!       modified: 11-2008 by min-ken hsieh.
!                 pass in bkg params to calcualte jb
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(3),mlds(3),nobs,indx(6,nobs)
  real, intent(in) :: grid(ngrd(1),ngrd(2),ngrd(3))
  real, intent(in) :: coef(6,nobs),obsv(4,nobs),wght(nobs)
  real, intent(in) :: pnlt
  real, intent(out) :: fctn

  integer, intent(in) :: huc(mlds(1),mlds(2),mlds(3))	! uncover by min-ken

  ! local variables:
  integer :: i,j,k,io
  real    :: gvl,pnl

  ! initial:
  fctn = 0.0
  penalt = pnlt

  ! jb term:
  do k=1,ngrd(3)
    do j=1,ngrd(2)
      do i=1,ngrd(1)
        fctn = fctn+huc(i,j,k)*(grid(i,j,k)**2)
      enddo
    enddo
  enddo

  ! jo term:
  do io=1,nobs

    ! grid value at obs position:
    gvl = 0.0
    do k=3,6,3
      do j=2,5,3
	do i=1,4,3
	  gvl = gvl+grid(indx(i,io), &
			 indx(j,io), &
			 indx(k,io))* &
		coef(i,io)*coef(j,io)*coef(k,io)
	enddo
      enddo
    enddo

    ! residual:
    fctn = fctn+(obsv(1,io)-gvl)**2

  enddo

  fctn = fctn*0.5

  ! penalty term: penalize the laplacian:

  ! x:
  pnl = 0.0
  do k=1,ngrd(3)
    do j=1,ngrd(2)
      do i=2,ngrd(1)-1
        pnl = pnl + &
          (grid(i-1,j,k)-2.0*grid(i,j,k)+grid(i+1,j,k))**2
      enddo
    enddo
  enddo

  ! y:
  do k=1,ngrd(3)
    do j=2,ngrd(2)-1
      do i=1,ngrd(1)
        pnl = pnl + &
          (grid(i,j-1,k)-2.0*grid(i,j,k)+grid(i,j+1,k))**2
      enddo
    enddo
  enddo

  ! t:
  do k=2,ngrd(3)-1
    do j=1,ngrd(2)
      do i=1,ngrd(1)
        pnl = pnl + &
          (grid(i,j,k-1)-2.0*grid(i,j,k)+grid(i,j,k+1))**2
      enddo
    enddo
  enddo

  fctn = fctn+penalt*pnl

end subroutine functn3d

subroutine gradnt3d(grdt,grid,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

!==========================================================
!  this routine evaluates gradients of stmas cost function.
!
!  history:
!	creation: jun. 2005 by yuanfu xie.
!       modified: 10-2008 by min-ken hsieh.
!                 pass in penalty for each var
!       modified: 11-2008 by min-ken hsieh.
!                 pass in hbg, huc to calcualte jb
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(3),mlds(3),nobs,indx(6,nobs)
  real, intent(in) :: grid(ngrd(1),ngrd(2),ngrd(3))
  real, intent(in) :: coef(6,nobs),obsv(4,nobs),wght(nobs)
  real, intent(in) :: pnlt
  real, intent(out) :: grdt(ngrd(1),ngrd(2),ngrd(3))

  integer, intent(in) :: huc(mlds(1),mlds(2),mlds(3))	! uncover array by min-ken

  ! local variables:
  integer :: i,j,k,io
  real :: gvl,pnl
  real :: lpl(ngrd(1),ngrd(2),ngrd(3))
  real :: gll(ngrd(1),ngrd(2),ngrd(3))

  ! initial:
  grdt = 0.0
  penalt = pnlt

  ! jb term:
  
  do k=1,ngrd(3)
    do j=1,ngrd(2)
      do i=1,ngrd(1)
        grdt(i,j,k) = huc(i,j,k)*(grid(i,j,k))
      enddo
    enddo
  enddo


  ! jo term:
  do io=1,nobs

    ! grid value at obs position:
    gvl = 0.0
    do k=3,6,3
      do j=2,5,3
	do i=1,4,3
	  gvl = gvl+grid(indx(i,io), &
			 indx(j,io), &
			 indx(k,io))* &
		coef(i,io)*coef(j,io)*coef(k,io)
	enddo
      enddo
    enddo

    ! gradient:
    do k=3,6,3
      do j=2,5,3
        do i=1,4,3
          grdt(indx(i,io),indx(j,io),indx(k,io)) = &
            grdt(indx(i,io),indx(j,io),indx(k,io))- &
            (obsv(1,io)-gvl)*coef(i,io)*coef(j,io)*coef(k,io)
        enddo
      enddo
    enddo

  enddo

  ! penalty:

  gll = 0.0
  ! x:
  lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) = &
         grid(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3)) &
    -2.0*grid(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) &
        +grid(3:ngrd(1)  ,1:ngrd(2),1:ngrd(3))
  gll(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3)) = &
  gll(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3)) + &
        lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3))
  gll(3:ngrd(1)  ,1:ngrd(2),1:ngrd(3)) = &
  gll(3:ngrd(1)  ,1:ngrd(2),1:ngrd(3)) + &
        lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3))
  gll(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) = &
  gll(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) - &
    2.0*lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3))

  ! y:
  lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) = &
         grid(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3)) &
    -2.0*grid(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) &
        +grid(1:ngrd(1),3:ngrd(2)  ,1:ngrd(3))
  gll(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3)) = &
  gll(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3)) + &
        lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3))
  gll(1:ngrd(1),3:ngrd(2)  ,1:ngrd(3)) = &
  gll(1:ngrd(1),3:ngrd(2)  ,1:ngrd(3)) + &
        lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3))
  gll(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) = &
  gll(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) - &
    2.0*lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3))

  ! t:
  lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) = &
         grid(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2) &
    -2.0*grid(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) &
        +grid(1:ngrd(1),1:ngrd(2),3:ngrd(3)  )
  gll(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2) = &
  gll(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2) + &
        lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1)
  gll(1:ngrd(1),1:ngrd(2),3:ngrd(3)  ) = &
  gll(1:ngrd(1),1:ngrd(2),3:ngrd(3)  ) + &
        lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1)
  gll(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) = &
  gll(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) - &
    2.0*lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1)

  ! add gll to grdt:
  grdt(1:ngrd(1),1:ngrd(2),1:ngrd(3)) = &
    grdt(1:ngrd(1),1:ngrd(2),1:ngrd(3)) + &
    2.0*penalt*gll(1:ngrd(1),1:ngrd(2),1:ngrd(3))

end subroutine gradnt3d

subroutine stmasinc

!==========================================================
!  this routine adds the analysis increments to background
!  fields using the laps land/water factors.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  ! local variables:
  integer i,j,k,l

  do l=1,numvar
    if (needbk(l) .eq. 1) then
      do k=1,numtmf

        ! add increment to the background:
	analys(1:numgrd(1),1:numgrd(2),k,l) = &
	  analys(1:numgrd(1),1:numgrd(2),k,l)+ &
	  bkgrnd(1:numgrd(1),1:numgrd(2),k,l)

	! treat land factor:
!	analys(1:numgrd(1),1:numgrd(2),k,l) = &
!	  lndfac(1:numgrd(1),1:numgrd(2))* &
!	  analys(1:numgrd(1),1:numgrd(2),k,l)+ &
!	  (1.0-lndfac(1:numgrd(1),1:numgrd(2)))* &
!	  bkgrnd(1:numgrd(1),1:numgrd(2),k,l)
      enddo
    endif
  enddo

end subroutine stmasinc

subroutine intplbkg(ucvr,ngrd,huc,mgd,mld,dgd,dxyt,domn)

!==========================================================
!  this routine interpolate background field and uncover to
!  multigrid.
!
!  history:
!	creation: 11-2008 by min-ken hsieh
!
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(3),mgd(3),mld(3)
  real, intent(in) :: dgd(3),dxyt(3),domn(2,3)
  logical, intent(in) :: ucvr(ngrd(1),ngrd(2),ngrd(3))

  integer, intent(out) :: huc(mld(1),mld(2),mld(3))

  ! local variables:
  integer :: ier,i,j,k,ix,iy,it
  real :: xyt(3)
  integer :: idx(6,mgd(1),mgd(2),mgd(3))
  real :: coe(6,mgd(1),mgd(2),mgd(3))
  logical :: lhuc(mld(1),mld(2),mld(3))

  ! interpolation for each background grid:
  do k=1,mgd(3)
    xyt(3) = domn(1,3)+(k-1)*dgd(3)
    do j=1,mgd(2)
      xyt(2) = domn(1,2)+(j-1)*dgd(2)
      do i=1,mgd(1)
        xyt(1) = domn(1,1)+(i-1)*dgd(1)
        call intplt3d(xyt(1:3),ngrd,dxyt,domn, &
	    	      idx(1,i,j,k),coe(1,i,j,k),ier)

        ! check:
        if (ier .ne. 0) then
          print*, 'background interpolation error!!',i,j,k,xyt(1),xyt(2),xyt(3)
        else
          ! evaluate:
          huc(i,j,k) = 0
	  lhuc(i,j,k) = .true.
	  do it=3,6,3
	    do iy=2,5,3
	      do ix=1,4,3
		! see if this grid is uncovered
		lhuc(i,j,k) = lhuc(i,j,k) .and. &
		  ucvr(idx(ix,i,j,k),idx(iy,i,j,k),idx(it,i,j,k))
		   
	      enddo
	    enddo
	  enddo

          if(lhuc(i,j,k)) then
            huc(i,j,k) = 1
          endif

        endif
        
      enddo
    enddo
  enddo

end subroutine intplbkg

end module stmasanalz
