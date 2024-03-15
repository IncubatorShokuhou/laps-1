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
!dis     technical support to users.
!dis
!dis    permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  all modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  if significant modifications or enhancements
!dis    are made to this software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

module definition

!==========================================================
!  this module defines all parameters and variables needed
!  in stmas multigrid analysis.
!
!  history:
!	creation: yuanfu xie 8-2005
!==========================================================

!==========================================================
!  define parameters.
!==========================================================

  integer, parameter :: maxvar=15
  integer, parameter :: lsxvar=21

  !****************
  ! laps constants:
  !****************

  ! lapse rates: see laps mdatlaps.f under sfc:
  ! 1: temperature; 2: dewpoint
  real, parameter :: lapses(2) = (/-0.01167, -0.007/)

  ! unit conversions:
  real, parameter :: temp_0 = 273.16		! absolute 0
  real, parameter :: knt2ms = 0.51444444	! knot 2 m/s
  real, parameter :: mile2m = 1609.0		! mile to meter
  real, parameter :: mb2pas = 100.0		! mb 2 pascal
  real, parameter :: inch2m = 0.0254		! inch to meter
  real, parameter :: gascnt = 287.0		! gas constant 
						! for dry air
  real, parameter :: spheat = 1004.0		! specific heat 
						! at constant pres

  !****************************
  ! machine related parameters:
  !****************************

  real, parameter :: epsiln = 1.0e-18

!==========================================================
!  this header file defines all necessary variables.
!
!  note: all global variables are six letter long.
!
!  history: 
!	creation: yuanfu xie	6-2005
!       modified: yuanfu xie    12-2008 for adding diagnl.
!==========================================================

  !----------------
  ! laps variables:
  !----------------
  character*4 :: varnam(maxvar)
  character*256 :: dirstc

  character*8 :: date0,date1
  character*10 :: time0,time1
  character*5 :: zone0,zone1
  integer :: timing0(8),timing1(8) 

  integer :: dirlen

  integer :: lapsdt		! laps cycle time (secs)
  integer :: mxstts		! maximum number obs sites
  integer :: i4time		! laps analysis time 
  integer :: i4wndw(2)		! analysis time window
  integer :: qc_val		! qc flag for threshold check
  integer :: qc_std		! qc flag for standard dev-check
  integer, allocatable, dimension(:) &
	  :: i4prev		! lapstimes

  real   :: mising		! laps bad data flag
  real   :: badsfc		! laps bad sfc flag
  real   :: thresh(maxvar)	! threshold for obs agaist bkg
  real   :: pnlt_v(maxvar)	! penalty of each variable !added by min-ken hsieh, used in stmasana
  real   :: rdplvl		! reduced pressure level definted in laps
  real, allocatable, dimension(:,:) &
	 :: latgrd,longrd,topogr! grid lat/lon/topo
  real, allocatable, dimension(:,:,:) &
	 :: rawobs		! raw observations
  real, allocatable, dimension(:,:) &
	 :: bkgobs		! background at obs sites

  !-----------------------------------------------------
  ! interactive variables between laps and minimization:
  !-----------------------------------------------------
  integer :: numfic(3)		! number fictitious points
  integer :: numgrd(3)		! analysis grid numbers
  integer :: numtmf		! number of laps bkgd frames
  integer :: numvar		! number of analysis vars
  integer :: needbk(maxvar)	! add background to incements
  integer :: bounds(maxvar)	! bound constraints
  integer :: radius(maxvar)     ! obs cover grids ! added by min-ken hsieh, used in addbkgrd
  integer :: lndsea(maxvar)     ! land/sea process! added by min-ken hsieh, used in addbkgrd
  integer :: slevel(maxvar)     ! starting level of analysis! added by min-ken hsieh, used in stmasana
  integer :: verbal		! print message option
  integer :: press_pert		! 1: compute pressure perturbation
  integer :: savdat		! save background and obs
  integer :: saveid		! index for saving variable
  integer :: numobs(maxvar)	! maximum number of obs
  integer, allocatable, dimension(:,:,:) &
	:: indice		! interpolation indices

  real :: grdspc(3)		! gridspacing in x, y, t
  real :: domain(2,3)		! analysis domain
  real :: obsspc(3,maxvar)	! observation spacing
  real, allocatable, dimension(:,:,:,:) &
       :: bkgrnd		! background fields
  real, allocatable, dimension(:,:) &
       :: lndfac		! land factor
  real, allocatable, dimension(:,:,:) &
       :: qc_obs		! qced observations
  real, allocatable, dimension(:,:) &
       :: weight		! observations weights
  real, allocatable, dimension(:,:,:) &
       :: coeffs		! interpolation coefficients
  real, allocatable, dimension(:,:,:) &
       :: diagnl		! diagonal array of b for j_b term
  logical, allocatable, dimension(:,:,:,:) &
       :: uncovr		! uncovered grids !added by min-ken hsieh, used in addbkgrd and stmasana

  !------------------------
  ! minimization variables:
  !------------------------

  integer :: stmasi(1)		! stmas integer parameters
  integer :: maxitr		! maximum iterations

  real :: phydxy(2)		! physical spacing for derivatives
  real :: stmasr(1)		! stmas real parameters
  real :: penalt		! penalty parameter
  real, allocatable, dimension(:,:,:,:) &
	:: analys		! analyzed fields
  !------------------------
  !verification variables:
  !                             !modified by min-ken hsieh
  !------------------------
  character*20,allocatable, dimension(:,:) &
	 :: stanam 		!store stn name for each variable

  ! namelists:
  namelist /stmas/numfic,numtmf,numgrd,numvar,savdat,saveid,verbal,press_pert
  namelist /stmas/qc_val,qc_std
  namelist /stmas/maxitr,stmasi,penalt,stmasr

contains			! common utility routines

subroutine grid2obs(indx,coef,obsv,nobs,wght,stna,ngrd,dxyt,domn)

!==========================================================
!  this routine finds the indices and coefficients for an
!  interpolation scheme from a grid function to observation
!  sites.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!       modification:
!                 25-08-2008 by min-ken hsieh
!                 add parameter stna to map each obs its stn name for stmasver
!
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(3)
  real, intent(in) :: dxyt(3),domn(2,3)

  integer, intent(out) :: indx(6,nobs)
  integer, intent(inout) :: nobs
  real, intent(out) :: coef(6,nobs)
  real, intent(inout) :: obsv(4,nobs),wght(nobs)
  character*20, intent(inout) :: stna(nobs)		!by min-ken hsieh

  ! local variables:
  integer :: i,ier
  integer :: nib		! number of obs in box

  ! count obs in box:
  nib = 0
  ! interpolation for each obs:
  do i=1,nobs
    call intplt3d(obsv(2:4,i),ngrd,dxyt,domn, &
		  indx(1,i),coef(1,i),ier)

    ! check:
    if (ier .eq. 0) then
      nib = nib+1

      ! save the obs and its weight:
      obsv(1:4,nib) = obsv(1:4,i)
      wght(nib) = wght(i)
      stna(nib) = stna(i)

      indx(1:6,nib) = indx(1:6,i)
      coef(1:6,nib) = coef(1:6,i)

    else
	if (verbal .eq. 1) &
         print*,'grid2obs: obs out of the analysi domain ',obsv(2:4,i),i,domn

    endif
  enddo

  ! count inbox obs:
  nobs = nib

end subroutine grid2obs

subroutine intplt3d(pstn,ngrd,gspc,domn,indx,coef,ierr)

!==========================================================
!  this routine returns interpolation coefficients and
!  indices of a given location from a grid in 3d space.
!
!  history:
!	creation: 9-2005 by yuanfu xie.
!==========================================================

  implicit none

  integer, intent(in) :: ngrd(3)	! numbers gridpoint
  real,    intent(in) :: pstn(3)	! interpolate point
  real,    intent(in) :: gspc(3)	! grid spacing
  real,    intent(in) :: domn(2,3)	! domain

  integer, intent(out) :: indx(6)	! indices
  integer, intent(out) :: ierr		! error flag
  real,    intent(out) :: coef(6)	! coefficients

  ! local variables:
  integer :: i

  ierr = 0
  ! check if it is in box?
  do i=1,3
     if (pstn(i) .lt. domn(1,i)) ierr = -i
     if (pstn(i) .gt. domn(2,i)) ierr = i
  enddo

  ! indices:
  indx(1:3) = (pstn-domn(1,1:3))/gspc

  ! coefficients:
  coef(1:3) = (pstn-indx(1:3)*gspc-domn(1,1:3))/gspc

  indx(1:3) = indx(1:3)+1

  indx(4:6) = min(indx(1:3)+1,ngrd)
  coef(4:6) = coef(1:3)
  coef(1:3) = 1.0-coef(4:6)

end subroutine intplt3d

end module definition
