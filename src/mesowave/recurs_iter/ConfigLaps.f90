module configlaps

!==========================================================
!  this module sets up laps configuration and access to its
!  lso observation data.
!
!  history: may. 2004 by yuanfu xie 
!           adapted from dr. mcginley's codes.
!==========================================================

  use definition

  implicit none

  ! laps variables:
  integer, parameter :: nvlaps = 6 	!number of variables
  integer, parameter :: ncycles = 6     !number of cycles to carry
  integer   :: nx,ny,nfic 			! grid dimensions
  integer   :: laps_cycle_time,len,i1,i2,i3,i4,istatus,ilaps,ivar
  integer   :: maxstations 		! maximum number of ob stations
  integer   :: istarttime
  character :: ext_s*30,dir_s*200      
  character :: units*60,comment*60,name*100
  character :: maproj*6,nest7grid*9,var_s*3 ! map projection character string
  real      :: stanlat,stanlat2,stanlon,badflag,grid_spacingx,grid_spacingy
  ! gridded lat, lon , observation arrary o (variables, stations*time)
  real, allocatable, dimension (:,:) :: lat,lon,ldf,olaps
  ! observation lat, lon, ob time and weight
  real, allocatable, dimension (:) :: olat,olon,otime,wght

  ! background fields:
  real,   allocatable, dimension (:,:,:,:) :: bkgd

contains

  include 'gridbarnes.f90'
  include 'lapsinfo.f90'
  include 'lso_data_qc.f90'
  include 'writeanalysis.f90'

end module
