subroutine laps_config

!==============================================================================
!doc  this routine configures laps in order to ingest data by
!doc  laps libraries.
!doc
!doc  history:
!doc	creation:	yuanfu xie	may 2007
!==============================================================================

  use laps_params
  use mem_namelist

  implicit none

  ! local variables:
  character*4 :: varnam
  integer :: status

  ! variables for accessing wind parameters:
  character*150 :: static,filenm
  integer       :: length

  integer :: nyear,nmonth,nday,nhour,nmin,nsec

  format_request = 'bufr'	! testing!!!

  ! spatial dim: get_grid_dim_xy also sets laps_config:
  call get_grid_dim_xy(number_gridpts(1),number_gridpts(2),status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading xy dimensions', status,number_gridpts(1:2)
    stop
  endif
  call get_laps_dimensions(number_gridpts(3),status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading vertical dimension'
    stop
  endif

  ! system time:
  call get_systime(system_in4time,system_asctime,status)
  !call get_systime_all(system_in4time,system_asctime,analys_hourcha, &
  !                     analys_minutes,string_asctime,yyyear_julians, status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in getting system time'
    stop
  endif
  ! convert i4time to yyyymm_ddhhmin:
  call cv_i4tim_int_lp(system_in4time,nyear,nmonth,nday,nhour, &
                        nmin,nsec,status)
  ! call read_systim(string_asctime,number_asctime,yyyymm_ddhhmin,status)
  yyyymm_ddhhmin(1) = 1900+nyear
  yyyymm_ddhhmin(2) = nmonth
  yyyymm_ddhhmin(3) = nday
  yyyymm_ddhhmin(4) = nhour
  yyyymm_ddhhmin(5) = nmin

  ! laps cycle time:
  call get_laps_cycle_time(length_anatime,status)
  length_anatime = 0.5*length_anatime
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading analysis time window'
    stop
  endif

  ! missing data:
  call get_r_missing_data(rvalue_missing,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading real missing value'
    stop
  endif
  call get_i2_missing_data(ivalue_missing,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading integer missing value'
    stop
  endif
  call get_sfc_badflag(sfcobs_invalid,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading bad value for surface'
    stop
  endif

  ! domain configuration:
  call laps_malloc
  call read_static_grid(number_gridpts(1),number_gridpts(2),'lat',&
                        domain_latitde,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading latitude'
    stop
  endif
  call read_static_grid(number_gridpts(1),number_gridpts(2),'lon',&
                        domain_longitd,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading longitude'
    stop
  endif
  call read_static_grid(number_gridpts(1),number_gridpts(2),'avg',&
                        domain_topogrp,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading longitude'
    stop
  endif
  ! note: grid spacing currently shows the same distance in x and y.
  ! thus grid4d_spacing's first element used in get_grid_spacing_actual:
  call get_grid_spacing_actual( &
    domain_latitde((number_gridpts(1)-1)/2+1, &
                   (number_gridpts(2)-1)/2+1), &
    domain_longitd((number_gridpts(1)-1)/2+1, &
                   (number_gridpts(2)-1)/2+1), grid4d_spacing,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading grid x-y spacing'
    stop
  endif
  grid4d_spacing(2) = grid4d_spacing(1)
  ! height field:
  varnam = 'ht'
  call get_modelfg_3d(system_in4time,varnam,number_gridpts(1), &
                      number_gridpts(2),number_gridpts(3), &
                      height_grid3dm,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error in reading height field'
    stop
  endif
  ! pressure levels:
  call get_pres_1d(system_in4time,number_gridpts(3),pressr_grid1dm,status)
  ! temperature:
  varnam = 't3'
  call get_modelfg_3d(system_in4time,varnam,number_gridpts(1), &
                      number_gridpts(2),number_gridpts(3), &
                temptr_grid3dm,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error retrieving temperature'
    stop
  endif
  ! u wind:
  varnam = 'u3'
  call get_modelfg_3d(system_in4time,varnam,number_gridpts(1), &
                      number_gridpts(2),number_gridpts(3), &
                uuwind_grid3dm,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error retrieving u wind'
    stop
  endif
  ! v wind:
  varnam = 'v3'
  call get_modelfg_3d(system_in4time,varnam,number_gridpts(1), &
                      number_gridpts(2),number_gridpts(3), &
                vvwind_grid3dm,status)
  if (status .ne. 1) then
    write(6,*) 'laps_config: error retrieving v wind'
    stop
  endif
  ! write background temporarily:
  ! write(11,*) number_gridpts(1:3)
  ! write(11,*) height_grid3dm,temptr_grid3dm,uuwind_grid3dm,vvwind_grid3dm
  ! write(11,*) domain_latitde,domain_longitd,domain_topogrp,pressr_grid1dm

  ! write(12,*) height_grid3dm,pressr_grid1dm
  
  ! get wind parameters:
  ! call get_wind_parms(useobs_raobdat,useobs_cldrfwd,useobs_radarwd, &
  !		      thresh_radarob(1),thresh_radarob(2),thresh_radarob(3), &
  !		      weight_bkgwind,weight_radarwd,thresh_rmswind, &
  !                    maxnum_proflrs,maxlvl_proflrs,weight_options,status)
  ! if (status .ne. 1) then
  !   write(6,*) 'laps_config: error in reading wind parameters'
  ! endif
 
  ! use a new laps wind parameter scheme:
  call get_directory('static',static,length)
  filenm = static(1:length)//'/wind.nl'
  call read_namelist_laps ('wind',filenm)
  maxnum_proflrs = max_pr		! laps wind parameter max_pr;
  maxlvl_proflrs = max_pr_levels		! laps wind parameter max_pr_lvls;

  maxnum_sondes = max_snd_grid
  maxlvl_sondes = max_snd_levels

end subroutine laps_config


subroutine laps_malloc

!==============================================================================
!doc  this routine allocates memory for laps variables.
!doc
!doc  history:
!doc	creation:	yuanfu xie	may 2007
!==============================================================================

  use laps_params

  implicit none

  ! local variable:
  integer :: status

  ! two dimension arrays:
  allocate(domain_latitde(number_gridpts(1),number_gridpts(2)), &
           domain_longitd(number_gridpts(1),number_gridpts(2)), &
	   domain_topogrp(number_gridpts(1),number_gridpts(2)), &
	   pressr_grid1dm(number_gridpts(3)), &
           height_grid3dm(number_gridpts(1),number_gridpts(2),  &
                          number_gridpts(3)), &
           temptr_grid3dm(number_gridpts(1),number_gridpts(2),  &
                          number_gridpts(3)), &
           uuwind_grid3dm(number_gridpts(1),number_gridpts(2),  &
                          number_gridpts(3)), &
           vvwind_grid3dm(number_gridpts(1),number_gridpts(2),  &
                          number_gridpts(3)), stat=status)
  if (status .ne. 0) then
    write(6,*) 'laps_malloc: error in allocating memory for lat/lon/topo'
    stop
  endif

end subroutine laps_malloc


subroutine laps_dalloc

!==============================================================================
!doc  this routine allocates memory for laps variables.
!doc
!doc  history:
!doc	creation:	yuanfu xie	may 2007
!==============================================================================

  use laps_params

  implicit none

  ! local variable:
  integer :: status

  ! two dimension arrays:
  deallocate(domain_latitde, domain_longitd, domain_topogrp, pressr_grid1dm, &
		height_grid3dm, temptr_grid3dm, uuwind_grid3dm, &
		vvwind_grid3dm,stat=status)
  if (status .ne. 0) then
    write(6,*) 'laps_dalloc: error in deallocating memory for lat/lon/topo'
    stop
  endif

end subroutine laps_dalloc

