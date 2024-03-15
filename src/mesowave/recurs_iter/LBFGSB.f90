
!====================================================
!  this header file declares necessary variables for
!  the optimization package: lbfgsb.
!
!  author: yuanfu xie   noaa/fsl
!  date:   august 2000
!
!====================================================
      
!.....variables required by lbfgs_b:

integer, parameter :: msave=7
integer, parameter :: mvar = mx*my*mt
character*60       :: ctask,csave
double precision,allocatable,dimension(:) :: &
	wk,bdlow,bdupp
double precision   :: factr,dsave(29)
integer            :: iprnt,isbmn,isave(44)
integer,allocatable,dimension(:) :: nbund,iwrka
logical            :: lsave(4)

!.....end of lbfgs_b declarations.
