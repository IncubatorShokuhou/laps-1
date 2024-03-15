subroutine minimize(id,ds)

!*********************************************************
!  this routine minimizes the cost function derived from
!  a surface analysis using recursive filter.
!
!  history: jan. 2004 by yuanfu xie.
!           sep. 2004 by yuanfu xie penalizing div.
!*********************************************************

  implicit none

  integer, intent(in) :: id
  real,    intent(in) :: ds(3)

  ! lbfgs_b variables:
  include 'lbfgsb.f90'

  ! local variables:
  integer          :: itr,nvar,idp,istatus
  ! adjoint variables:
  double precision :: f,adf,dv(mvar),g(mvar)
  real             :: ada(mx,my,mt,2),rv(mvar)
  
  ! unified u/v analysis:
  idp = id
  if (id .eq. 201) idp = id+1

  ! start lbfgs_b
  ctask = 'start'
  factr = 1.0d+2
  iprnt = 1
  isbmn = 1

  ! allocate memory:
  allocate(bdlow(mvar),bdupp(mvar),nbund(mvar),iwrka(3*mvar), &
	wk(mvar*(2*msave+4)+12*msave*msave+12*msave), &
	stat=istatus)
  if (istatus .ne. 0) then
     print*,'minimize: no space for lbfgsb workspace'
     stop
  endif

  nbund = 0

  ! initial:
  itr = 0
  nvar = n(1)*n(2)*n(3)*(idp-id+1)
  rv(1:nvar) = reshape(a(1:n(1),1:n(2),1:n(3),id:idp),(/nvar/))
  dv(1:nvar) = rv(1:nvar)

  ! looping:
1 continue
  call lbfgsb(nvar,msave,dv,bdlow,bdupp,nbund,f,g,factr,wk, &
              iwrka,ctask,iprnt,isbmn,csave,lsave,isave,dsave)

  ! transfer a vector to a grid:
  rv(1:nvar) = dv(1:nvar)
  a(1:n(1),1:n(2),1:n(3),id:idp) = reshape(rv(1:nvar), &
                               (/n(1),n(2),n(3),idp-id+1/))
  ! exit if succeed:
  if (ctask(1:11) .eq. 'convergence') goto 2

  ! function and gradient are needed:
  if (ctask(1:2) .eq. 'fg') then

     ! function value:
     if ((id .ne. 201) .and. (id .ne. 301)) then
        call functn(f,a(1,1,1,id),l,n,id,np,al)
     else
        ! call functndiv(f,a(1,1,1,id),l,n,ds,id,np,al)
     endif

     ! gradient value:
     adf = 1.0d0
     ada = 0.0
     if (id .ne. 201) then
        call adfunctn( a(1,1,1,id), l, n, id, np, al, adf, ada )
     else
	! call adfunctndiv( a(1,1,1,id), l, n, ds, id, np, al, adf, ada )
     endif

     rv(1:nvar) = reshape(ada(1:n(1),1:n(2),1:n(3),1:idp-id+1), &
	(/ nvar /))
     g(1:nvar) = rv(1:nvar)

  endif

  ! exit if irregularity is encountered:
  if ((ctask(1:2) .ne. 'fg') .and. (ctask(1:5) .ne. 'new_x')) then
     print*,'error in lbfgs_b'
     goto 2
  endif

  ! count number of iteration and compute relative error:
  if (ctask(1:5) .eq. 'new_x') then
     itr = itr+1
     ! print*,''
  endif

  ! if number of iterations of lbfgs_b exceeds the limit:
  if (itr .lt. maxitr) goto 1

  ! when an error in lbfgs_b occurs: exit
2 continue

  ! convert to analysis:
  call rf3d(a(1,1,1,id),l,n,al(1,id),np(1,id))
  if (id .eq. 201) &
     call rf3d(a(1,1,1,idp),l,n,al(1,idp),np(1,idp))

  ! deallocate memory:
  deallocate(bdlow,bdupp,nbund,iwrka,wk,stat=istatus)
  
end subroutine minimize
