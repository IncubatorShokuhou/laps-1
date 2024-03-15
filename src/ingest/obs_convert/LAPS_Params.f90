module laps_params

!==============================================================================
!doc  this routine defines all parameters required for laps.
!doc
!doc  history:
!doc	creation:	yuanfu xie	may 2007
!doc    modified:       yuanfu xie	mar 2008 adding radar variables
!==============================================================================

  implicit none

  ! constants:
  integer,parameter :: output_channel=20,bfrtbl_channel=21,points_channel=22
  ! for lapsplot of observations:
  integer,parameter :: tmgout_channel=23,pigout_channel=24,prgout_channel=25,&
                       sagout_channel=26
  real,parameter :: inches_conv2mm=2.54		! from inch to milli-meter

  ! bufr encoding:
  character*80, parameter :: header_prebufr= &
    'sid typ t29 rpt dhr yob xob elv itp sqn procn'
  integer,      parameter :: header_numitem=11
  character*80, parameter :: obsdat_prebufr= &
    'zob pob tob uob vob qob pmo prss pwo'
  integer,      parameter :: obsdat_numitem=9
  character*80, parameter :: obserr_prebufr= &
    'zoe poe toe woe qoe pwe'
  integer,      parameter :: obserr_numitem=6
  character*80, parameter :: obsqms_prebufr= &
    'zqm pqm tqm wqm qqm pmq pwq'
  integer,      parameter :: obsqms_numitem=7

  character*80, parameter :: sfcdat_prebufr= &
    'prss ??? tob ??? wdir1 nwspd1'
  integer,      parameter :: sfcdat_numitem=9
  character*80, parameter :: sfcerr_prebufr= &
   'sfcerr prebufr unknown'
  integer,      parameter :: sfcerr_numitem=6
  character*80, parameter :: sfcqms_prebufr= &
   'sfcqms prebufr unknown'
  integer,      parameter :: sfcqms_numitem=7

  ! missing data:
  real*8,       parameter :: missng_prebufr=10.0e10

  ! reference value for ice/water relative humidity:
  real,         parameter :: absolu_tmpzero=273.15
  real,         parameter :: temptr_referen=-5.0 !-132.0

  ! general vars:
  character*8 :: format_request		! data format requested to convert to
  character*9 :: system_asctime 	! ascii system time
  character   :: analys_hourcha*2, &	! analysis hour
                 analys_minutes*2, &	! analysis minute
                 yyyear_julians*5	! yy year and julian days
  character   :: string_asctime*16, &	! time in a string
                 number_asctime*14      ! time in yyyymmddhhmm form
  integer     :: number_gridpts(4)	! number grid points in x, y, z and t
  integer     :: length_anatime		! anal time window
  integer     :: yyyymm_ddhhmin(5)	! time in integers
  integer     :: system_in4time		! i4 system time
  integer     :: radars_timetol         ! radar time tolerance
  integer     :: ivalue_missing		! missing value in integer
  real        :: rvalue_missing		! missing value in real
  real        :: sfcobs_invalid 	! bad surface observation
  real        :: grid4d_spacing(4)	! grid distance of 4 dimemsion
  real, allocatable, dimension(:) :: &
		pressr_grid1dm		! pressure levels
  real, allocatable, dimension(:,:) :: &
	        domain_latitde,        & ! grid latitudes
		domain_longitd,        & ! grid longitude
		domain_topogrp           ! grid topography
  real, allocatable,dimension(:,:,:) :: &
                height_grid3dm,&	! laps 3d height background field
                temptr_grid3dm,&	! laps 3d height background field
                uuwind_grid3dm,&	! laps 3d height background field
                vvwind_grid3dm	! laps 3d height background field

  ! use a new laps wind parameter scheme and so the followings are covered
  ! except those used:

  ! for wind data:
  ! logical     :: useobs_raobdat	! use raob data
  ! logical     :: useobs_cldrfwd	! use cloud drift wind
  ! logical     :: useobs_radarwd	! use radial wind
  ! integer     :: thresh_radarob(3)	! threshold numbers doppler obs/level
					! by factor 2,4,9
  integer     :: maxnum_proflrs	! max num profilers
  integer     :: maxlvl_proflrs	! max level of prof

  integer     :: maxnum_sondes	! max num sondes
  integer     :: maxlvl_sondes 	! max level of sondes

  ! integer     :: weight_options	! laps weight scheme
  ! real        :: weight_bkgwind	! weight for background wind
  ! real        :: weight_radarwd	! weight for radar wind obs
  ! real        :: thresh_rmswind	! threshold for rms fit of analysis to
					! the observations

  ! name list for wind:
  ! namelist /wind_nl/ useobs_raobdat,useobs_cldrfwd, &
  !		     useobs_radarwd,thresh_radarob, &
  !		     weight_bkgwind,weight_radarwd, &
  !		     thresh_rmswind,maxnum_proflrs, &
  !		     maxlvl_proflrs,weight_options

end module laps_params
