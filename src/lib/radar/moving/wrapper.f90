!====================================================================
!  this wrapper is designed to scan all possible airborne radar data
!  files and to read all superobs data.
!
!  history:
!
!       creation: yuanfu xie august 2009.
!====================================================================

subroutine wrapper(i4time,halftime,path,maxdat,numdat,obsi4t, &
                   lat,lon,hgh,azi,elv,dis,wnd)

  implicit none

  character*256, intent(in) :: path	! path to the airborne data
  integer, intent(in) :: i4time,halftime,maxdat
  integer, intent(out) ::numdat
  integer, intent(out) :: obsi4t(maxdat)	! obs i4time
  real,    intent(out) :: lat(maxdat),lon(maxdat),hgh(maxdat), &
                          azi(maxdat),elv(maxdat) ! azimuth/elevation angles
  real,    intent(out) :: dis(maxdat),wnd(maxdat) ! radius and radial wind

  ! local variables:
  integer, parameter :: max_files = 2000
  integer :: ifile,nfiles,l1,l2,istatus
  integer :: begini4,endngi4
  character*256 :: c_filenames(max_files),filter
  character*13 :: begintime,endngtime,obstime
  character*9 :: begin9,endng9,wfo_fname13_to_fname9

  ! initialize numdat:
  numdat = 0

  ! list of all available files:
  ! call get_file_time(path,nfiles,c_filenames,max_files,istatus)
  filter = '*'
  call getfilenames_c(path,c_filenames,nfiles,filter,istatus)

  ! for file in the time window:
  ! hardcopy for handling time after year 2000 only: !!!!
  begintime(1:2) = '20'
  endngtime(1:2) = '20'
  begintime(9:9) = '_'
  endngtime(9:9) = '_'
  do ifile=1,nfiles
    ! yymmdd:
    begintime(3:8) = c_filenames(ifile)(1:6)
    endngtime(3:8) = c_filenames(ifile)(1:6)
    ! hhmm:
    begintime(10:13) = c_filenames(ifile)(10:13)
    endngtime(10:13) = c_filenames(ifile)(15:18)

    ! convert to fname9:
    begin9 = wfo_fname13_to_fname9(begintime)
    endng9 = wfo_fname13_to_fname9(endngtime)

    ! convert to i4 time:
    call cv_asc_i4time(begin9,begini4,istatus)
    call cv_asc_i4time(endng9,endngi4,istatus)

    ! find a overlap:
    if (i4time+halftime .ge. begini4 .and. i4time-halftime .le. endngi4) then

      ! open data file:
      call s_len(path,l1)
      call s_len(c_filenames(ifile),l2)
      open(13,file=path(1:l1)//c_filenames(ifile)(1:l2),status='old')

      ! read a record:
      numdat = numdat+1
      if (numdat .gt. maxdat) then
        print*,'error: too many radial wind data: ',numdat,' > ',maxdat
        stop
      endif

 1    read(13,*,end=2) obstime(1:12),lat(numdat),lon(numdat),hgh(numdat), &
                       azi(numdat),elv(numdat),dis(numdat),wnd(numdat)

      obstime(10:13) = obstime(9:12)
      obstime(9:9) = '_'
      obstime(1:9) = wfo_fname13_to_fname9(obstime(1:13))

      ! obs i4time:
      call cv_asc_i4time(obstime(1:9),obsi4t(numdat),istatus)
      ! skip invalid obs time:
      if (obsi4t(numdat) .lt. i4time-halftime .or. &
          obsi4t(numdat) .gt. i4time+halftime) then
        numdat = numdat-1
      endif
      ! next record:
      goto 1
 2    close(13)

    endif
  enddo

end subroutine wrapper
