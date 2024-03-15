subroutine read_rass_data(i4time_sys,cycle_time,max_rass,max_levels,rvalue_missing, &
                            i4time_obs,num_rass,lvl_rass,c5_names,c8_type, &
                            lat_rass,lon_rass,elv_rass,err_rass, &
                            hgt_rass,tmp_rass,istatus)

!==============================================================================
!doc  this routine reads rass raw data from laps lrs files.
!doc  it is the same as the first part of read_tsnd.f
!doc
!doc  history:
!doc	creation:	yuanfu xie	jan. 2009
!==============================================================================

  implicit none

  integer, intent(in) :: i4time_sys,cycle_time,max_rass,max_levels
  real,    intent(in) :: rvalue_missing

  character, intent(out) :: c5_names(max_rass)*5,c8_type(max_rass)*8
  integer,   intent(out) :: num_rass,lvl_rass(max_rass),i4time_obs(max_rass)
  real,      intent(out) :: lat_rass(max_rass),lon_rass(max_rass), &
                            elv_rass(max_rass),err_rass(max_rass), &
                            hgt_rass(max_rass,max_levels), &		! m
                            tmp_rass(max_rass,max_levels)		! c

  ! local variables:
  character :: ext*3, c_filespec*255,a9time*9
  integer   :: i4time_file,lag_time,i4time_rass_offset, &
               sttnid,i_qc,i_rass,i_level,n_level,istatus
  logical   :: l_string_contains
  real      :: rcycles
  real,parameter :: surface_rass_buffer = 30.0
  
  ext = 'lrs'
  call get_filespec(ext,2,c_filespec,istatus)
  call get_file_time(c_filespec,i4time_sys,i4time_file)

  lag_time = 0 ! middle of rass hourly sampling period
  i4time_rass_offset = i4time_sys - (i4time_file + lag_time)
  rcycles = float(i4time_rass_offset) / float(cycle_time)

  write(6,*)' i4time_rass_offset/rcycles = ',i4time_rass_offset,rcycles

  ! initialization:
  num_rass = 0
  lvl_rass = 0

  if (abs(i4time_rass_offset) .gt. cycle_time) then
    write(6,*)' rass file/offset is > laps cycle time'
    write(6,*)' skipping the use of rass'
    return
  endif

  call open_lapsprd_file_read(12,i4time_file,ext,istatus)
  if (istatus .ne. 1) then
    print*,'read_rass_data: cannot open data file: ',ext
    return
  endif

  do i_rass = 1,max_rass

    num_rass = num_rass+1	! count

340 continue

    read(12,401,end=500) sttnid,n_level, &
	lat_rass(num_rass),lon_rass(num_rass),elv_rass(num_rass), &
        c5_names(num_rass),a9time,c8_type(num_rass)
401 format(i12,i12,f11.0,f15.0,f15.0,5x,a5,3x,a9,1x,a8)

    ! convert a9time to i4time:
    call cv_asc_i4time(a9time,i4time_obs(num_rass),istatus)

    if (n_level .gt. max_levels) then
      write(6,*)'read_rass_data error: too many levels in the lrs file: ', &
                i_rass,n_level,max_levels
      istatus = 0
      return
    endif

    if (l_string_contains(c8_type(num_rass),'sat',istatus)) then
       err_rass(num_rass) = 5.0
    else
       err_rass(num_rass) = 1.0
    endif

    do i_level = 1,n_level

      read(12,*,err=340) hgt_rass(num_rass,i_level),tmp_rass(num_rass,i_level),i_qc

      ! suppose lrs files contain temperature in kelvin:
      if (i_qc .eq. 1 .and. &
          hgt_rass(num_rass,i_level) .lt. rvalue_missing .and. &
          tmp_rass(num_rass,i_level) .gt. 200.0 .and. &
          tmp_rass(num_rass,i_level) .lt. 450.0 .and. &
          hgt_rass(num_rass,i_level) .gt. elv_rass(num_rass)+surface_rass_buffer .and. &
          i_level .le. max_levels ) then

        lvl_rass(num_rass) = lvl_rass(num_rass) + 1		! count levels

        hgt_rass(num_rass,lvl_rass(num_rass)) = hgt_rass(num_rass,i_level)
        tmp_rass(num_rass,lvl_rass(num_rass)) = tmp_rass(num_rass,i_level)-273.15	! in centigrade

      endif

    enddo ! level

    if (lvl_rass(num_rass) .eq. 0) num_rass=num_rass-1

  enddo ! rass

500 continue
  num_rass=num_rass-1	! when reaching the end of file, num_rass has been added

  close(12)	! close lrs file

end subroutine read_rass_data
