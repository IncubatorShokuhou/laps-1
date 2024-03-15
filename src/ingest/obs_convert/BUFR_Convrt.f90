subroutine bufr_proflr(numproflr,numlevels,stationid,i4obstime, &
                        latitudes,longitude,elevation,obsvntype, &
                        maxproflr,heightobs,uuwindobs,vvwindobs, &
                        wndrmserr,prssfcobs,tmpsfcobs,rehsfcobs, &
                        uwdsfcobs,vwdsfcobs)

!==============================================================================
!doc  this routine converts laps profiler data into prepbufr format.
!doc
!doc  history:
!doc	creation:	yuanfu xie/shiow-ming deng	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character, intent(in) :: stationid(*)*5		! station id
  character, intent(in) :: obsvntype(*)*8		! obs type
  integer,   intent(in) :: numproflr,numlevels(*)	! number of profilers/levels
  integer,   intent(in) :: maxproflr			! maximum number profilrs
                                                	! instead of maxnum_proflrs
                                                	! avoid confusion on memory.
  integer,   intent(in) :: i4obstime(*)			! i4 obs times

  ! observations:
  real,      intent(in) :: latitudes(*),longitude(*),elevation(*)
  real,      intent(in) :: heightobs(maxproflr,*), &	! upair (m)
                           uuwindobs(maxproflr,*), &	! upair (m/s)
                           vvwindobs(maxproflr,*), &	! upair (m/s)
                           wndrmserr(maxproflr,*)	! wind rms error (m/s)
  real,      intent(in) :: prssfcobs(*),tmpsfcobs(*), &	! surface
                           rehsfcobs(*), &		! surface
                           uwdsfcobs(*),vwdsfcobs(*)	! surface

  ! local variables:
  character :: sttnid*8,subset*8
  integer   :: i,j,k,indate,zeroch,status
  real      :: ri,rj,rk,di,sp,height_to_pressure	! obs gridpoint locatioin
  real      :: make_ssh					! laps routine from rh to sh
  real*8    :: header(header_numitem),obsdat(obsdat_numitem,225)
  real*8    :: obserr(obserr_numitem,225),obsqms(obsqms_numitem,225)
  equivalence(sttnid,header(1))

  print*,'number of pofilers: ',numproflr,numlevels(1:numproflr)

  ! obs date: year/month/day
  zeroch = ichar('0')
  indate = yyyymm_ddhhmin(1)*1000000+yyyymm_ddhhmin(2)*10000+ &
           yyyymm_ddhhmin(3)*100+yyyymm_ddhhmin(4)

  ! write data:
  do i=1,numproflr

    ! station id:
    sttnid = stationid(i)

    ! data type:
    if (obsvntype(i) .eq. 'profiler') then
      subset = 'proflr'
      header(2) = 223		! profiler code
      header(3) = 71		! input report type
    else if (obsvntype(i) .eq. 'vad') then
      subset = 'vadwnd'
      header(2) = 224		! profiler code
      header(3) = 72		! input report type
    else
      print*,'bufr_proflr: error: unknown profiler data type!'
      stop
    endif

    ! time:
    if (abs(i4obstime(i)-system_in4time) .gt. length_anatime) cycle	! out of time

    ! laps cycle time:
    header(4) = yyyymm_ddhhmin(4)
    ! dhr: obs time different from the cycle time:
    header(5) = (i4obstime(i)-system_in4time)/3600.0+yyyymm_ddhhmin(5)/60.0

    ! lat/lon/elevation:
    header(6) = latitudes(i)
    header(7) = longitude(i)
    header(8) = elevation(i)

    ! ignore obs outside the analysis domain:
    call latlon_to_rlapsgrid(latitudes(i),longitude(i), &
                              domain_latitde,domain_longitd, &
                              number_gridpts(1),number_gridpts(2), &
                              ri,rj,status)
    if (ri .lt. 1 .or. ri .gt. number_gridpts(1) .or. &
        rj .lt. 1 .or. rj .gt. number_gridpts(2) .or. &
        status .ne. 1) cycle

    header(9) = 99			! instrument type: common code table c-2
    header(10) = i			! report sequence number
    header(11) = 0			! mpi process number

    ! upair observations:
    ! missing data conversion:
    obsdat = missng_prebufr
    obserr = missng_prebufr
    obsqms = missng_prebufr

    do j=1,numlevels(i)
      if (heightobs(i,j) .ne. rvalue_missing .and. &
          (uuwindobs(i,j) .ne. rvalue_missing) .and. &
          (vvwindobs(i,j) .ne. rvalue_missing)) then

        obsdat(1,j) = heightobs(i,j)
        obsdat(4,j) = uuwindobs(i,j)
        obsdat(5,j) = vvwindobs(i,j)
        ! assuming 1 m/s error as default:
        obserr(4,j) = 1.0
        if (wndrmserr(i,j) .ne. rvalue_missing) &
          obserr(4,j) = wndrmserr(i,j)

        ! save data into prg:
        rk = height_to_pressure(heightobs(i,j),height_grid3dm,pressr_grid1dm, &
                       number_gridpts(1),number_gridpts(2),number_gridpts(3), &
                       nint(ri),nint(rj))

        ! find the grid level for this pressure value:
        do k=1,number_gridpts(3)
          if (rk .ge. pressr_grid1dm(k)) exit
        enddo
        if (k .gt. number_gridpts(3)) then
          rk = number_gridpts(3)+1			! out of height
        elseif (k .le. 1) then
          rk = 1
          if (rk .gt. pressr_grid1dm(1)) rk = 0	! out of height
        else
          rk = k-(rk-pressr_grid1dm(k))/(pressr_grid1dm(k-1)-pressr_grid1dm(k))
        endif

        ! convert to di and sp:
        call uv_to_disp(uuwindobs(i,j),vvwindobs(i,j),di,sp)

        write(prgout_channel,11) ri,rj,rk,di,sp,obsvntype(i)
11	format(3f8.1,2f10.3,3x,a8)

      endif
    enddo
    obserr(1,1:numlevels(i)) = 10.0	! assume 10 meter error for height
    obsqms(1,1:numlevels(i)) = 0	! quality mark - bufr code table: 
    obsqms(4,1:numlevels(i)) = 0	! 0 always assimilated

    ! write to bufr file:
    call openmb(output_channel,subset,indate)
    call ufbint(output_channel,header,header_numitem,1,status,header_prebufr)
    call ufbint(output_channel,obsdat,obsdat_numitem, &
                 numlevels(i),status,obsdat_prebufr)
    call ufbint(output_channel,obserr,obserr_numitem, &
                 numlevels(i),status,obserr_prebufr)
    call ufbint(output_channel,obsqms,obsqms_numitem, &
                 numlevels(i),status,obsqms_prebufr)
    call writsb(output_channel)

    ! write surface data: future development debug needed
    ! sfcdat(1) = prssfcobs(i)
    ! sfcdat(2) = missng_prebufr
    ! sfcdat(3) = tmpsfcobs(i)
    ! sfcdat(4) = missng_prebufr
    ! use -132 as temp_reference:
    ! temperature >= temp_reference: rh is water rh;
    ! temperature <  temp_reference: rh is ice rh;
    ! assume the surface obs of rh is water rh here:
    ! sfcdat(5) = make_ssh(prssfcobs(i),tmpsfcobs(i),rehsfcobs(i)/100.0,&
    !              temptr_referen)*0.001 ! kg/kg
    ! call openmb(output_channel,'adpsfc',indate)
    ! header is needed to reflect the observation code!!!
    ! call ufbint(output_channel,header,header_numitem,1,1,header_prebufr)
    ! call ufbint(output_channel,sfcdat,surfac_numitem,1,1,surfac_prebufr)
    ! call ufbint(output_channel,sfcqms,sfcqms_numitem,1,1,surfac_prebufr)
    ! call writsb(output_channel)

  enddo

end subroutine bufr_proflr

subroutine bufr_rass(numproflr,numlevels,stationid,i4obstime, &
                        latitudes,longitude,elevation,obsvntype, &
                        maxproflr,heightobs,temptrobs,temptrerr)

!==============================================================================
!doc  this routine converts laps profiler data into prepbufr format.
!doc
!doc  history:
!doc	creation:	yuanfu xie/shiow-ming deng	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character, intent(in) :: stationid(*)*5		! station id
  character, intent(in) :: obsvntype(*)*8		! obs type
  integer,   intent(in) :: numproflr,numlevels(*)	! number of profilers/levels
  integer,   intent(in) :: maxproflr			! maximum number profilrs
                                                	! instead of maxnum_proflrs
                                                	! avoid confusion on memory.
  integer,   intent(in) :: i4obstime(*)			! i4 obs times

  ! observations:
  real,      intent(in) :: latitudes(*),longitude(*),elevation(*)
  real,      intent(in) :: heightobs(maxproflr,*), &	! upair (m)
                           temptrobs(maxproflr,*), &	! upair (c)
                           temptrerr(maxproflr)		! temperature error (c)

  ! local variables:
  character :: sttnid*8,subset*8
  integer   :: i,j,k,indate,zeroch,status
  real      :: ri,rj,rk,di,sp,height_to_pressure	! obs gridpoint locatioin
  real      :: make_ssh					! laps routine from rh to sh
  real*8    :: header(header_numitem),obsdat(obsdat_numitem,225)
  real*8    :: obserr(obserr_numitem,225),obsqms(obsqms_numitem,225)
  equivalence(sttnid,header(1))

  print*,'number of pofilers: ',numproflr,numlevels(1:numproflr)

  ! obs date: year/month/day
  zeroch = ichar('0')
  indate = yyyymm_ddhhmin(1)*1000000+yyyymm_ddhhmin(2)*10000+ &
           yyyymm_ddhhmin(3)*100+yyyymm_ddhhmin(4)

  ! write data:
  do i=1,numproflr

    ! station id:
    sttnid = stationid(i)

    ! data type:
    subset = 'adpupa'
    header(2) = 120		! profiler code
    header(3) = 77		! input report type

    ! time:
    if (abs(i4obstime(i)-system_in4time) .gt. length_anatime) cycle	! out of time

    ! laps cycle time:
    header(4) = yyyymm_ddhhmin(4)
    ! dhr: obs time different from the cycle time:
    header(5) = (i4obstime(i)-system_in4time)/3600.0+yyyymm_ddhhmin(5)/60.0

    ! lat/lon/elevation:
    header(6) = latitudes(i)
    header(7) = longitude(i)
    header(8) = elevation(i)

    ! ignore obs outside the analysis domain:
    call latlon_to_rlapsgrid(latitudes(i),longitude(i), &
                              domain_latitde,domain_longitd, &
                              number_gridpts(1),number_gridpts(2), &
                              ri,rj,status)
    if (ri .lt. 1 .or. ri .gt. number_gridpts(1) .or. &
        rj .lt. 1 .or. rj .gt. number_gridpts(2) .or. &
        status .ne. 1) cycle

    header(9) = 99			! instrument type: common code table c-2
    header(10) = i			! report sequence number
    header(11) = 0			! mpi process number

    ! upair observations:
    ! missing data conversion:
    obsdat = missng_prebufr
    obserr = missng_prebufr
    obsqms = missng_prebufr

    do j=1,numlevels(i)
      if (heightobs(i,j) .ne. rvalue_missing .and. &
          (temptrobs(i,j) .ne. rvalue_missing) ) then

        obsdat(1,j) = heightobs(i,j)
        obsdat(3,j) = temptrobs(i,j)
        ! assuming 1 degree error as default:
        obserr(3,j) = 1.0
        if (temptrerr(i) .ne. rvalue_missing) &
          obserr(3,j) = temptrerr(i)

        ! save data into tmg:
        rk = height_to_pressure(heightobs(i,j),height_grid3dm,pressr_grid1dm, &
                       number_gridpts(1),number_gridpts(2),number_gridpts(3), &
                       nint(ri),nint(rj))
        ! find the grid level for this pressure value:
        do k=1,number_gridpts(3)
          if (rk .ge. pressr_grid1dm(k)) exit
        enddo
        if (k .gt. number_gridpts(3)) then
          rk = number_gridpts(3)+1			! out of height
        elseif (k .le. 1) then
          rk = 1
          if (rk .gt. pressr_grid1dm(1)) rk = 0	! out of height
        else
          rk = k-(rk-pressr_grid1dm(k))/(pressr_grid1dm(k-1)-pressr_grid1dm(k))
        endif

        write(tmgout_channel,11) ri,rj,rk,temptrobs(i,j)+237.15,obsvntype(i)
11	format(3f10.4,f10.3,3x,a8)

      endif
    enddo
    obserr(1,1:numlevels(i)) = 10.0	! assume 10 meter error for height
    obsqms(1,1:numlevels(i)) = 0	! quality mark - bufr code table: 
    obsqms(4,1:numlevels(i)) = 0	! 0 always assimilated

    ! write to bufr file:
    call openmb(output_channel,subset,indate)
    call ufbint(output_channel,header,header_numitem,1,status,header_prebufr)
    call ufbint(output_channel,obsdat,obsdat_numitem, &
                 numlevels(i),status,obsdat_prebufr)
    call ufbint(output_channel,obserr,obserr_numitem, &
                 numlevels(i),status,obserr_prebufr)
    call ufbint(output_channel,obsqms,obsqms_numitem, &
                 numlevels(i),status,obsqms_prebufr)
    call writsb(output_channel)

    ! write surface data: future development debug needed
    ! sfcdat(1) = prssfcobs(i)
    ! sfcdat(2) = missng_prebufr
    ! sfcdat(3) = tmpsfcobs(i)
    ! sfcdat(4) = missng_prebufr
    ! use -132 as temp_reference:
    ! temperature >= temp_reference: rh is water rh;
    ! temperature <  temp_reference: rh is ice rh;
    ! assume the surface obs of rh is water rh here:
    ! sfcdat(5) = make_ssh(prssfcobs(i),tmpsfcobs(i),rehsfcobs(i)/100.0,&
    !              temptr_referen)*0.001 ! kg/kg
    ! call openmb(output_channel,'adpsfc',indate)
    ! header is needed to reflect the observation code!!!
    ! call ufbint(output_channel,header,header_numitem,1,1,header_prebufr)
    ! call ufbint(output_channel,sfcdat,surfac_numitem,1,1,surfac_prebufr)
    ! call ufbint(output_channel,sfcqms,sfcqms_numitem,1,1,surfac_prebufr)
    ! call writsb(output_channel)

  enddo

end subroutine bufr_rass

subroutine bufr_sondes(numsondes,numlevels,stationid,i4obstime, &
                         latitudes,longitude,elevation,obsvntype,heightobs, &
                         pressrobs,temptrobs,dewpntobs,uuwindobs,vvwindobs)

!==============================================================================
!doc  this routine converts laps sonde data into prepbufr format.
!doc
!doc  history:
!doc	creation:	yuanfu xie/shiow-ming deng	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character, intent(in) :: stationid(*)*5		! station id
  character, intent(in) :: obsvntype(*)*8		! obs type
  integer,   intent(in) :: numsondes                    ! number of profilers
  integer               :: numlevels(*)			! number of levels
  integer,   intent(in) :: i4obstime(*)			! i4 obs times

  ! observations:
  real,      intent(in) :: latitudes(*),longitude(*),elevation(*)
  real            :: heightobs(maxnum_sondes,*), &	! upair (m)
                           pressrobs(maxnum_sondes,*), &	! upair (mb)
                           temptrobs(maxnum_sondes,*), &	! upair (c)
                           dewpntobs(maxnum_sondes,*), &	! upair (c)
                           uuwindobs(maxnum_sondes,*), &	! upair (m/s)
                           vvwindobs(maxnum_sondes,*)		! upair (m/s)

  ! local variables:
  character :: sttnid*8,subset*8
  integer   :: i,j,k,indate,zeroch,status
  real      :: ssh2,ri,rj,rk,di,sp,p,height_to_pressure	! laps function for specific humidity
  real*8    :: header(header_numitem),obsdat(obsdat_numitem,225)
  real*8    :: obserr(obserr_numitem,225),obsqms(obsqms_numitem,225)
  equivalence(sttnid,header(1))

  ! obs date: year/month/day
  zeroch = ichar('0')
  indate = yyyymm_ddhhmin(1)*1000000+yyyymm_ddhhmin(2)*10000+ &
           yyyymm_ddhhmin(3)*100+yyyymm_ddhhmin(4)

  ! write data:
  do i=1,numsondes

    ! no valid level data: in aug. 2010, yuanfu added 150 upper limit
    if (numlevels(i) .le. 0) cycle

    ! set maximum vertical levels limit of 150: burf has problem for too many levels (yuanfu)
    if (numlevels(i) .gt. 150) numlevels(i) = 150

    ! station id:
    sttnid = stationid(i)

    ! data type:
    subset = 'adpupa'
    select case (obsvntype(i))
    case ('radiomtr')
      header(2) = 120		! prepbufr report type: table 2
      header(3) = 11		! input report type: table 6
      obserr(1,1:numlevels(i)) = 0.0	! zero error for height
      obserr(3,1:numlevels(i)) = 2.0	! assume 2 deg error
      obserr(4,1:numlevels(i)) = 0.5	! assume 0.5 m/s error
      obserr(6,1:numlevels(i)) = 2.0	! assume 2 mm error
    case ('raob','tower')
      header(2) = 120           ! prepbufr report type: table 2
      header(3) = 11            ! input report type: table 6
      obserr(1,1:numlevels(i)) = 0.0    ! zero error for height
      obserr(3,1:numlevels(i)) = 2.0    ! assume 2 deg error
      obserr(4,1:numlevels(i)) = 0.5    ! assume 0.5 m/s error
      obserr(6,1:numlevels(i)) = 2.0    ! assume 2 mm error

      ! test setting raob height obs off:
      heightobs(i,2:numlevels(i)) = rvalue_missing	! assume no actual height obs above ground
    case ('poessnd')
      header(2) = 257		! prepbufr report type: table 2
      header(3) = 63		! input report type: table 6: satellite-derived wind
      obserr(1,1:numlevels(i)) = 0.0	! zero error for height
      obserr(3,1:numlevels(i)) = 2.0	! assume 2 deg error
      obserr(4,1:numlevels(i)) = 1.0	! 2 m/s error assumed for satwind
      obserr(6,1:numlevels(i)) = 2.0	! assume 2 mm error
    case ('goes11','goes12','goes13')
      header(2) = 151		! prepbufr report type: table 2
      header(3) = 63		! input report type: table 6: satellite-derived wind
      obserr(1,1:numlevels(i)) = 0.0	! zero error for height
      obserr(3,1:numlevels(i)) = 2.0	! assume 2 deg error
      obserr(4,1:numlevels(i)) = 1.0	! 2 m/s error assumed for satwind
      obserr(6,1:numlevels(i)) = 2.0	! assume 2 mm error
    case ('dropsnd')
      header(2) = 132           ! prepbufr report type: table 2
      header(3) = 31            ! input report type: table 6: dropsonde
      obserr(1,1:numlevels(i)) = 0.0    ! zero error for height
      obserr(3,1:numlevels(i)) = 2.0    ! assume 2 deg error
      obserr(4,1:numlevels(i)) = 1.0    ! 2 m/s error assumed for satwind
      obserr(6,1:numlevels(i)) = 2.0    ! assume 2 mm error
    case default
      print*,'bufr_sondes: unknown observation data type! ',obsvntype(i),i
      stop
    end select

    ! time:
    if (abs(i4obstime(i)-system_in4time) .gt. length_anatime) cycle	! out of time

    header(4) = yyyymm_ddhhmin(4)
    header(5) = (i4obstime(i)-system_in4time)/3600.0+yyyymm_ddhhmin(5)/60.0

    ! lat/lon/elevation:
    header(6) = latitudes(i)
    header(7) = longitude(i)
    header(8) = elevation(i)

    ! ignore obs outside the analysis domain:
    call latlon_to_rlapsgrid(latitudes(i),longitude(i), &
                              domain_latitde,domain_longitd, &
                              number_gridpts(1),number_gridpts(2), &
                              ri,rj,status)
    if (ri .lt. 1 .or. ri .gt. number_gridpts(1) .or. &
        rj .lt. 1 .or. rj .gt. number_gridpts(2) .or. &
        status .ne. 1) cycle

    header(9) = 90			! instrument type: common code table c-2
    header(10) = i			! report sequence number
    header(11) = 0			! mpi process number

    ! upair observations:
    ! missing data conversion:
    obsdat = missng_prebufr
    obsqms = missng_prebufr
    do j=1,numlevels(i)

      ! no height or pressure info:
      p = pressrobs(i,j)*100	! pascal when finding grid level
      if (heightobs(i,j) .ge. rvalue_missing .and. &
           pressrobs(i,j) .ge. rvalue_missing) cycle	! invalid data

      if (heightobs(i,j) .ne. rvalue_missing) then
	obsdat(1,j) = heightobs(i,j)
        ! when pressure is missing, use height to convert pressure
        if (p .ge. rvalue_missing) &
          p = height_to_pressure(heightobs(i,j),height_grid3dm, &
                         pressr_grid1dm,number_gridpts(1),number_gridpts(2), &
                         number_gridpts(3),nint(ri),nint(rj))
        obserr(1,j) = 10.0	! 10 meter error for height obs
      endif

      ! find the grid level for this pressure value:
      do k=1,number_gridpts(3)
        if (p .ge. pressr_grid1dm(k)) exit
      enddo
      if (k .gt. number_gridpts(3)) then
        rk = number_gridpts(3)+1		! out of height
      elseif (k .le. 1) then
        rk = 1
        if (p .gt. pressr_grid1dm(1)) rk = 0	! out of height
      else
        rk = k-(p-pressr_grid1dm(k))/(pressr_grid1dm(k-1)-pressr_grid1dm(k))
      endif

      ! other obs:
      obsdat(2,j) = p/100.0	! bufr pressure in mb
      obserr(2,j) = 10.0

      if (temptrobs(i,j) .ne. rvalue_missing) then
	obsdat(3,j) = temptrobs(i,j)
	obserr(3,j) = 1.0	! 1 degree error for temperature

        ! save temperature data into tmg file:
        write(tmgout_channel,11) ri,rj,rk,temptrobs(i,j)+273.15,obsvntype(i)
11      format(3f10.4,f10.3,3x,a8)
      endif
      if (uuwindobs(i,j) .ne. rvalue_missing .and. &
           vvwindobs(i,j) .ne. rvalue_missing) then
        obsdat(4,j) = uuwindobs(i,j)
        obsdat(5,j) = vvwindobs(i,j)
        obserr(4,j) = 0.1

        ! save wind obs into prg file:
        call uv_to_disp(uuwindobs(i,j),vvwindobs(i,j),di,sp)

        write(prgout_channel,12) ri,rj,rk,di,sp,obsvntype(i)
12	format(3f8.1,2f10.3,3x,a8)
      endif

      ! specific humidity:
      if ((dewpntobs(i,j) .ne. rvalue_missing) .and. &
          (temptrobs(i,j) .ne. rvalue_missing) .and. &
           dewpntobs(i,j) .le. temptrobs(i,j)) then
        obsdat(6,j) = ssh2(p/100.0,temptrobs(i,j), &
                          dewpntobs(i,j),temptr_referen)*1000.0 !mg/kg
print*,'sh: ',p/100.0,obsdat(6,j)/1000.0,2.0*10.0**(-24.0/log10(p/100.0)+9.0), &
        0.5*ssh2(p/100.0,temptrobs(i,j),temptrobs(i,j),temptr_referen),j,numlevels(i)
        obserr(5,j) = 1.0	! assume 1mg/kg error by yuanfu
      endif
    enddo
    obsqms(1:6,1:numlevels(i)) = 1	! quality mark - bufr code table: 
					! good assumed.

    ! write to bufr file:
    call openmb(output_channel,subset,indate)
    call ufbint(output_channel,header,header_numitem,1,status,header_prebufr)
    call ufbint(output_channel,obsdat,obsdat_numitem, &
                 numlevels(i),status,obsdat_prebufr)
    call ufbint(output_channel,obserr,obserr_numitem, &
                 numlevels(i),status,obserr_prebufr)
    call ufbint(output_channel,obsqms,obsqms_numitem, &
                 numlevels(i),status,obsqms_prebufr)
    call writsb(output_channel)

  enddo

end subroutine bufr_sondes

subroutine bufr_sfcobs(numberobs,hhmintime,latitudes,longitude,stationid, &
                         obsvntype,providers,elevation,mslprsobs,mslprserr, &
                         stnprsobs,stnprserr,temptrobs,temptrerr,wind2dobs, &
                         wind2derr,relhumobs,relhumerr,sfcprsobs,precp1obs, &
                         precp1err)


!==============================================================================
!doc  this routine converts laps sonde data into prepbufr format.
!doc
!doc  history:
!doc	creation:	yuanfu xie/shiow-ming deng	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character, intent(in) :: stationid(*)*20		! station id
  character, intent(in) :: providers(*)*11		! provider's name
  character, intent(in) :: obsvntype(*)*6		! obs type
  integer,   intent(in) :: numberobs                    ! number of profilers
  integer,   intent(in) :: hhmintime(*)			! i4 obs times

  ! observations:
  real,      intent(in) :: latitudes(*),longitude(*),elevation(*)
  ! pressure in mb, temp in c, wind in m/s
  real,      intent(in) :: mslprsobs(*), &		! mean sea level pressure
                           mslprserr(*), &		! error
                           stnprsobs(*), &		! station pressure
                           stnprserr(*), &		! station pressure error
                           temptrobs(*), &		! temperature
                           temptrerr(*), &		! temperature error
                           wind2dobs(2,*), &  		! 2d wind obs
                           wind2derr(2,*), &  		! 2d windobs errkr
                           relhumobs(*), &  		! relative humidity
                           relhumerr(*), &  		! humidity error
                           sfcprsobs(*), &     		! surface pressure
                           precp1obs(*), &   		! precipitation
                           precp1err(*)			! precipitation error

  ! local variables:
  character :: sttnid*8,subset*8
  integer   :: i,k,indate,zeroch,status
  integer   :: i4time
  real      :: make_ssh		! laps function for specific humidity from rh
  real      :: ri,rj,rk,di,sp,height_to_pressure
  real*8    :: header(header_numitem),obsdat(obsdat_numitem)
  real*8    :: obserr(obserr_numitem),obsqms(obsqms_numitem)
  equivalence(sttnid,header(1))

  print*,'number of surface obs: ',numberobs

  ! obs date: year/month/day
  zeroch = ichar('0')
  indate = yyyymm_ddhhmin(1)*1000000+yyyymm_ddhhmin(2)*10000+ &
           yyyymm_ddhhmin(3)*100+yyyymm_ddhhmin(4)

  ! write data:
  do i=1,numberobs

    ! station id:
    sttnid = stationid(i)

    ! data type:
    subset = 'adpsfc'
    select case (obsvntype(i))
    case ('martim','synop')	! marine and synop
      header(2) = 281		! bufr report type:  table 2
      header(3) = 511		! input report type: table 6
    case ('metar','speci','ldad') ! speci: special metar data
      header(2) = 181		! bufr report type:  table 2
      header(3) = 512		! input report type: table 6
    case ('dropsn')
      header(2) = 132           ! prepbufr report type: table 2
      header(3) = 31            ! input report type: table 6: dropsonde
    case default
      print*,'bufr_sfcobs: unkown observation data type! ',obsvntype(i),' skip ',i
      cycle
      ! close(output_channel)
      ! stop
    end select

    ! time:
    call get_sfc_obtime(hhmintime(i),system_in4time,i4time,status)
    if (abs(i4time-system_in4time) .gt. length_anatime) cycle	! out of time

    header(4) = yyyymm_ddhhmin(4)
    header(5) = (i4time-system_in4time)/3600.0+yyyymm_ddhhmin(5)/60.0

    ! lat/lon/elevation:
    header(6) = latitudes(i)
    header(7) = longitude(i)
    header(8) = elevation(i)

    ! ignore obs outside the analysis domain:
    call latlon_to_rlapsgrid(latitudes(i),longitude(i), &
                              domain_latitde,domain_longitd, &
                              number_gridpts(1),number_gridpts(2), &
                              ri,rj,status)
    if (ri .lt. 1 .or. ri .gt. number_gridpts(1) .or. &
        rj .lt. 1 .or. rj .gt. number_gridpts(2) .or. &
        status .ne. 1) cycle

    header(9) = 90			! instrument type: common code table c-2
    header(10) = i			! report sequence number
    header(11) = 0			! mpi process number

    ! upair observations:
    ! missing data conversion:
    obsdat = missng_prebufr
    obserr = missng_prebufr
    obsqms = missng_prebufr

    ! surface obs: zob is the elevation height:
    obsdat(1) = elevation(i)
    obserr(1) = 10.0		! 10 meter error assumed

    ! get pressure from height first:
    rk = height_to_pressure(elevation(i),height_grid3dm, &
                         pressr_grid1dm,number_gridpts(1),number_gridpts(2), &
                         number_gridpts(3),nint(ri),nint(rj))

    ! surface pressure:
    if ((stnprsobs(i) .ne. rvalue_missing) .and. &
        (stnprsobs(i) .ne. sfcobs_invalid)) then
      obsdat(2) = stnprsobs(i)
      rk = stnprsobs(i)*100	! use pascal to find grid level
    endif

    ! find the grid level for this pressure value:
    do k=1,number_gridpts(3)
      if (rk .ge. pressr_grid1dm(k)) exit
    enddo
    if (k .gt. number_gridpts(3)) then
      rk = number_gridpts(3)+1			! out of height
    elseif (k .le. 1) then
      rk = 1
      if (rk .gt. pressr_grid1dm(1)) rk = 0	! out of height
    else
      rk = k-(rk-pressr_grid1dm(k))/(pressr_grid1dm(k-1)-pressr_grid1dm(k))
    endif

    ! other pressure obs:
    if ((mslprsobs(i) .ne. rvalue_missing) .and. &
        (mslprsobs(i) .ne. sfcobs_invalid)) obsdat(7) = mslprsobs(i)
    if ((sfcprsobs(i) .ne. rvalue_missing) .and. &
        (sfcprsobs(i) .ne. sfcobs_invalid)) obsdat(8) = sfcprsobs(i)

    ! temperature obs:
    if ((temptrobs(i) .ne. rvalue_missing) .and. &
        (temptrobs(i) .ne. sfcobs_invalid)) then
	obsdat(3) = temptrobs(i)

        ! temp error:
        obserr(3) = 1.0		! assume 1 degree default error
        if ((temptrerr(i) .ne. rvalue_missing) .and. &
            (temptrerr(i) .ne. sfcobs_invalid)) obserr(3) = temptrerr(i)

        ! save temperature into tmg file:
        write(tmgout_channel,11) ri,rj,rk,obsdat(3)+273.15,obsvntype(i)
11	format(3f10.4,f10.3,3x,a8)
    endif

    ! wind obs:
    if ((wind2dobs(1,i) .ne. rvalue_missing) .and. &
         (wind2dobs(1,i) .ne. sfcobs_invalid) .and. &
         (wind2dobs(2,i) .ne. rvalue_missing) .and. &
         (wind2dobs(2,i) .ne. sfcobs_invalid)) then
      obsdat(4) = wind2dobs(1,i)
      obsdat(5) = wind2dobs(2,i)

      ! save wind obs into prg file:
       call uv_to_disp(wind2dobs(1,i),wind2dobs(2,i),di,sp)

       write(prgout_channel,12) ri,rj,rk,di,sp,obsvntype(i)
12     format(3f8.1,2f10.3,3x,a8)
    endif

    obserr(4) = 0.1	! assume 0.1 m/s wind error as default
    if ((wind2derr(1,i) .ne. rvalue_missing) .and. &
        (wind2derr(1,i) .ne. sfcobs_invalid) .and. &
        (wind2derr(2,i) .ne. rvalue_missing) .and. &
        (wind2derr(2,i) .ne. sfcobs_invalid)) &
	obserr(4) = sqrt(wind2derr(1,i)**2+wind2derr(2,i)**2)
    ! specific humidity:
    if ((sfcprsobs(i) .ne. rvalue_missing) .and. &
        (sfcprsobs(i) .ne. sfcobs_invalid) .and. &
        (temptrobs(i) .ne. rvalue_missing) .and. &
        (temptrobs(i) .ne. sfcobs_invalid) .and. &
        (relhumobs(i) .ne. rvalue_missing) .and. &
        (relhumobs(i) .ne. sfcobs_invalid)) then
      obsdat(6) = make_ssh(sfcprsobs(i),temptrobs(i),relhumobs(i)/100.0,&
                           temptr_referen)*1000.0 ! mg/kg
      obserr(5) = 1.0	! assume 1mg/kg error by yuanfu
      !if ((stnprserr(i) .ne. rvalue_missing) .and. &
      !  (stnprserr(i) .ne. sfcobs_invalid) .and. &
      !  (temptrerr(i) .ne. rvalue_missing) .and. &
      !  (temptrerr(i) .ne. sfcobs_invalid) .and. &
      !  (relhumerr(i) .ne. rvalue_missing) .and. &
      !  (relhumerr(i) .ne. sfcobs_invalid)) &
      !obserr(5) = make_ssh(stnprserr(i),obserr(3),relhumerr(i)/100.0,&
      !                     temptr_referen)*1000.0 ! mg/kg
    endif
    if ((precp1obs(i) .ne. rvalue_missing) .and. &
        (precp1obs(i) .ne. sfcobs_invalid)) obsdat(9) = precp1obs(i)*inches_conv2mm
    if ((precp1err(i) .ne. rvalue_missing) .and. &
        (precp1err(i) .ne. sfcobs_invalid)) obserr(6) = precp1err(i)
    obsqms(1:5) = 0	! quality mark - bufr code table: 
			! 0 always assimilated.

    ! write to bufr file:
    call openmb(output_channel,subset,indate)
    call ufbint(output_channel,header,header_numitem,1,status,header_prebufr)
    call ufbint(output_channel,obsdat,obsdat_numitem,1,status,obsdat_prebufr)
    call ufbint(output_channel,obserr,obserr_numitem,1,status,obserr_prebufr)
    call ufbint(output_channel,obsqms,obsqms_numitem,1,status,obsqms_prebufr)
    call writsb(output_channel)

  enddo

end subroutine bufr_sfcobs

subroutine bufr_cdwaca(numberobs,obsvarray,obi4array)

!==============================================================================
!doc  this routine converts laps cloud drift wind and acars data into prepbufr
!doc  format.
!doc
!doc  history:
!doc	creation:	yuanfu xie/shiow-ming deng	jun 2007
!==============================================================================

  use laps_params

  implicit none

  integer, intent(in) :: numberobs,obi4array(3,*)
  real*8,  intent(in) :: obsvarray(7,*)

  ! local variables:
  character :: sttnid*8,subset*8,pigname*8
  integer   :: i,k,indate,zeroch,status
  real      :: ri,rj,rk,u,v,di,sp,height_to_pressure
  real*8    :: header(header_numitem),obsdat(obsdat_numitem)
  real*8    :: obserr(obserr_numitem),obsqms(obsqms_numitem)
  equivalence(sttnid,header(1))

  print*,'number of cloud drift wind and acar obs: ',numberobs

  ! obs date: year/month/day
  zeroch = ichar('0')
  indate = yyyymm_ddhhmin(1)*1000000+yyyymm_ddhhmin(2)*10000+ &
           yyyymm_ddhhmin(3)*100+yyyymm_ddhhmin(4)

  ! write data:
  do i=1,numberobs

    header(2:3) = obi4array(1:2,i)	! code and report type
    select case (obi4array(2,i))
    case (241)
      ! station id:
      sttnid = 'cdw'
      pigname = 'cdw'
      ! data type:
      subset = 'satwnd'
    case (130,230,666)	! temporarily use 666 for wisdom data
      ! station id:
      sttnid = 'acar'
      pigname = 'pin'
      if (obi4array(2,i) .eq. 666) pigname = 'wis'
      ! data type:
      subset = 'aircar'
    case default
      print*,'bufr_cdwaca: unknown observation data type! ',obi4array(2,i)
      stop
    end select

    ! time:
    if (abs(obi4array(3,i)-system_in4time) .gt. length_anatime) cycle	! out of time

    header(4) = yyyymm_ddhhmin(4)
    header(5) = (obi4array(3,i)-system_in4time)/3600.0+yyyymm_ddhhmin(5)/60.0

    ! lat/lon/elevation:
    header(6) = obsvarray(1,i)
    header(7) = obsvarray(2,i)
    header(8) = obsvarray(3,i)

    ! ignore obs outside the analysis domain:
    di = obsvarray(1,i)		! from real*8 to real
    sp = obsvarray(2,i)
    call latlon_to_rlapsgrid(di,sp, &
                              domain_latitde,domain_longitd, &
                              number_gridpts(1),number_gridpts(2), &
                              ri,rj,status)
    if (ri .lt. 1 .or. ri .gt. number_gridpts(1) .or. &
        rj .lt. 1 .or. rj .gt. number_gridpts(2) .or. &
        status .ne. 1) cycle

    header(9) = 90			! instrument type: common code table c-2
					! cannot find code table for acars instrument
    header(10) = i			! report sequence number
    header(11) = 0			! mpi process number

    ! upair observations:
    ! missing data conversion:
    obsdat = missng_prebufr
    obserr = missng_prebufr
    obsqms = missng_prebufr

    ! height obs:
    if (obsvarray(3,i) .ne. rvalue_missing) then
      obsdat(1) = obsvarray(3,i)	! height
      ! get pressure from height first:
      di = obsvarray(3,i)		! real*8 to real
      rk = height_to_pressure(di,height_grid3dm, &
                         pressr_grid1dm,number_gridpts(1),number_gridpts(2), &
                         number_gridpts(3),nint(ri),nint(rj))
    endif

    ! pressure obs:
    if (obsvarray(4,i) .ne. rvalue_missing) then
      obsdat(2) = obsvarray(4,i)	! pressure in mb
      rk = obsvarray(4,i)*100.0		! pascal
    endif
    ! find the grid level for this pressure value:
    do k=1,number_gridpts(3)
      if (rk .ge. pressr_grid1dm(k)) exit
    enddo
    if (k .gt. number_gridpts(3)) then
      rk = number_gridpts(3)+1			! out of height
    elseif (k .le. 1) then
      rk = 1
      if (rk .gt. pressr_grid1dm(1)) rk = 0	! out of height
    else
      rk = k-(rk-pressr_grid1dm(k))/(pressr_grid1dm(k-1)-pressr_grid1dm(k))
    endif

    ! temperature obs:
    if (obsvarray(7,i) .ne. rvalue_missing) then
      obsdat(3) = obsvarray(7,i)	! laps_ingest(c) uses read_acars_ob (f)

      ! save temperature obs into tmg file:
      write(tmgout_channel,11) ri,rj,rk,obsvarray(7,i)+273.15,sttnid
11    format(3f10.4,f10.3,3x,a8)
    endif

    ! wind obs:
    if (obsvarray(5,i) .ne. rvalue_missing .and. &
         obsvarray(6,i) .ne. rvalue_missing) then
      obsdat(4) = obsvarray(5,i)	! u
      obsdat(5) = obsvarray(6,i)	! v

      ! save wind obs into prg file:
      u = obsvarray(5,i)	! convert to real from real*8
      v = obsvarray(6,i)
      call uv_to_disp(u,v,di,sp)

      write(pigout_channel,12) ri,rj,rk,di,sp,pigname(1:3)
12    format(1x,3f8.1,2f8.1,x,a3)
    endif

    obserr(1) = 0.0	! zero error for height
    obserr(3) = 0.0	! assume 0 deg error
    obserr(4) = 0.0	! assume 0 m/s error
    obsqms(1:6) = 1	! quality mark - bufr code table: 
					! good assumed.

    ! write to bufr file:
    call openmb(output_channel,subset,indate)
    call ufbint(output_channel,header,header_numitem,1,status,header_prebufr)
    call ufbint(output_channel,obsdat,obsdat_numitem,1,status,obsdat_prebufr)
    call ufbint(output_channel,obserr,obserr_numitem,1,status,obserr_prebufr)
    call ufbint(output_channel,obsqms,obsqms_numitem,1,status,obsqms_prebufr)
    call writsb(output_channel)
  enddo

end subroutine bufr_cdwaca
