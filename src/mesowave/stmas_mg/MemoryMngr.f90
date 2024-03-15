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

module memorymngr

!==========================================================
!  this module manages memory allocation and deallocation.
!
!  history:
!	creation: yuanfu xie	8-2005
!==========================================================

  use definition

contains

subroutine lapsmemo

!==========================================================
!  this routine allocates necessary memory for stmas usage.
!
!  history: may. 2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: err

  allocate(i4prev(numtmf), &
	   latgrd(numgrd(1),numgrd(2)), &
	   longrd(numgrd(1),numgrd(2)), &
	   topogr(numgrd(1),numgrd(2)), &
	   rawobs(4,numtmf*mxstts,numvar), &
	   bkgobs(numtmf*mxstts,numvar), &
	   stanam(numtmf*mxstts,numvar), &		!added by min-ken.hsieh: stanam for stmasver
	   stat=err)
  if (err .ne. 0) then
    print*,'stmas>getmemos: cannot allocate enough memory!'
    stop
  endif

end subroutine lapsmemo

subroutine lapsrels

!==========================================================
!  this routine releases necessary memory for stmas usage.
!
!  history: may. 2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: err

!
! modified by min-ken hsieh
! because we still need bkgobs for verify
! we deallocate it in stmasver subroutine
!
  deallocate(i4prev, &
	     latgrd, &
	     longrd, &
	     topogr, &
	     rawobs, &
	     stat=err)

  if (err .ne. 0) then
    print*,'stmas>relsmemo: cannot delete allocated memory!'
    stop
  endif

end subroutine lapsrels

subroutine intrmemo

!==========================================================
!  this routine allocates necessary memory for interactive
!  variables.
!
!  history: aug. 2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: err

  allocate(bkgrnd(numgrd(1),numgrd(2),numtmf,numvar), &
           uncovr(numgrd(1),numgrd(2),numtmf,numvar), &	! add by min-ken, uncovered used in 
	   diagnl(numgrd(1),numgrd(2),numvar), &	! add b's diagonal array for j_b term
	   lndfac(numgrd(1),numgrd(2)), &		! addbkgrd and stmasana
	   qc_obs(4,numtmf*mxstts,numvar), &
	   weight(numtmf*mxstts,numvar), &
	   indice(6,numtmf*mxstts,numvar), &
	   coeffs(6,numtmf*mxstts,numvar), &
	   stat=err)
  if (err .ne. 0) then
    print*,'stmas>intrmemo: cannot allocate memory!'
    stop
  endif

end subroutine intrmemo

subroutine intrrels

!==========================================================
!  this routine releases necessary memory for interactive
!  variables.
!
!  history: aug. 2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: err

  deallocate(bkgrnd, &
	     uncovr, &
             diagnl, &
	     lndfac, &
	     qc_obs, &
	     weight, &
	     indice, &
	     coeffs, &
	     stat=err)
  if (err .ne. 0) then
    print*,'stmas>intrmemo: cannot release memory!'
    stop
  endif

end subroutine intrrels

subroutine stmasmem

!==========================================================
!  this routine allocates necessary memory for stmas 
!  analysis variables.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: err

  allocate(analys(numgrd(1),numgrd(2),numgrd(3),numvar), &
	   stat=err)
  if (err .ne. 0) then
    print*,'stmas>stmasmem: cannot allocate memory!'
    stop
  endif

end subroutine stmasmem

subroutine stmarels

!==========================================================
!  this routine allocates necessary memory for stmas 
!  analysis variables.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer :: err

  deallocate(analys, &
	     stat=err)
  if (err .ne. 0) then
    print*,'stmas>stmasmem: cannot release memory!'
    stop
  endif

end subroutine stmarels

end module memorymngr
