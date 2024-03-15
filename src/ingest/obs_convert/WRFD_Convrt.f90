subroutine wrfd_proflr(numprfls,numlevel,sttnames,obstimes,obsvnlat,obsvnlon, &
                       obsvnelv,reportyp,maxprflr,heightob,uwindobs,vwindobs,pssfcobs, &
                       tmpsfcob,rhsfcobs,usfcobsv,vsfcobsv)

!==============================================================================
!doc  this routine converts profiler data into wrf data format.
!doc
!doc  history:
!doc	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character,intent(in) :: reportyp(*)*8,sttnames(*)*5
  integer,  intent(in) :: numprfls,numlevel(*)
  integer,  intent(in) :: obstimes(maxprflr),maxprflr
  real,     intent(in) :: obsvnlat(maxprflr), &
                          obsvnlon(maxprflr), &
                          obsvnelv(maxprflr)
  real                 :: pssfcobs(*),tmpsfcob(*), &
                          rhsfcobs(*),usfcobsv(*),vsfcobsv(*)
  real                 :: heightob(maxprflr,*), &
                          uwindobs(maxprflr,*), &
                          vwindobs(maxprflr,*)

  ! local variables:
  character :: obtime*14,varch1*40,varch2*40
  integer   :: i,j,status
  real,parameter :: missng=-888888.000

  ! write out each surface station:
  do i=1,numprfls
    obtime = number_asctime
    call make_fnam_lp(obstimes(i),varch1,status)
    write(obtime(9:14),2) varch1(6:9),'00'
2   format(a4,a2)
    if (obtime(9:9) .eq. ' ') obtime(9:9) = '0'

    ! write obs time:
    write(output_channel,101) obtime

    ! write lat/lon:
    write(output_channel,102) obsvnlat(i),obsvnlon(i)

    ! write station name and data name:
    varch1 = '                                         '
    varch1(1:20) = sttnames(i)
    write(output_channel,1021) varch1,'profiler data from laps data ingest'

    ! write elevation etc:
    varch1 = '                                         '
    varch2 = '                                         '
    varch1 = reportyp(i)
    write(output_channel,103) varch1,varch2,obsvnelv(i),.true.,.false.,numlevel(i)

    ! write data at all levels:
    do j=1,numlevel(i)
      if ((pssfcobs(i) .eq. rvalue_missing) .or. &
          (pssfcobs(i) .eq. sfcobs_invalid) .or. j .gt. 1) pssfcobs(i) = missng
      if (heightob(i,j) .eq. rvalue_missing) heightob(i,j) = missng
      if ((tmpsfcob(i) .eq. rvalue_missing) .or. &
          (tmpsfcob(i) .eq. sfcobs_invalid) .or. j .gt. 1) tmpsfcob(i) = missng
      if (uwindobs(i,j) .eq. rvalue_missing) uwindobs(i,j) = missng
      if (vwindobs(i,j) .eq. rvalue_missing) vwindobs(i,j) = missng
      if ((rhsfcobs(i) .eq. rvalue_missing) .or. &
          (rhsfcobs(i) .eq. sfcobs_invalid) .or. j .gt. 1) rhsfcobs(i) = missng

      write(output_channel,104) pssfcobs(i),	0.0, &
                                heightob(i,j),  0.0, &
                                tmpsfcob(i),    0.0, &
                                uwindobs(i,j),  0.0, &
                                vwindobs(i,j),  0.0, &
                                rhsfcobs(i),    0.0
    enddo

  enddo

  ! wrf data formats:
  include 'wrfd_format.inc'

end subroutine wrfd_proflr

subroutine wrfd_rass(numprfls,numlevel,sttnames,obstimes,obsvnlat,obsvnlon, &
                       obsvnelv,reportyp,maxprflr,heightob,temptrob)

!==============================================================================
!doc  this routine converts profiler data into wrf data format.
!doc
!doc  history:
!doc	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character,intent(in) :: reportyp(*)*8,sttnames(*)*5
  integer,  intent(in) :: numprfls,numlevel(*)
  integer,  intent(in) :: obstimes(maxprflr),maxprflr
  real,     intent(in) :: obsvnlat(maxprflr), &
                          obsvnlon(maxprflr), &
                          obsvnelv(maxprflr)
  real                 :: heightob(maxprflr,*), &
                          temptrob(maxprflr,*)

  ! local variables:
  character :: obtime*14,varch1*40,varch2*40
  integer   :: i,j,status
  real,parameter :: missng=-888888.000

  ! write out each surface station:
  do i=1,numprfls
    obtime = number_asctime
    call make_fnam_lp(obstimes(i),varch1,status)
    write(obtime(9:14),2) varch1(6:9),'00'
2   format(a4,a2)
    if (obtime(9:9) .eq. ' ') obtime(9:9) = '0'

    ! write obs time:
    write(output_channel,101) obtime

    ! write lat/lon:
    write(output_channel,102) obsvnlat(i),obsvnlon(i)

    ! write station name and data name:
    varch1 = '                                         '
    varch1(1:20) = sttnames(i)
    write(output_channel,1021) varch1,'profiler data from laps data ingest'

    ! write elevation etc:
    varch1 = '                                         '
    varch2 = '                                         '
    varch1 = reportyp(i)
    write(output_channel,103) varch1,varch2,obsvnelv(i),.true.,.false.,numlevel(i)

    ! write data at all levels:
    do j=1,numlevel(i)
      if (heightob(i,j) .eq. rvalue_missing) heightob(i,j) = missng
      if (temptrob(i,j) .eq. rvalue_missing) temptrob(i,j) = missng

      write(output_channel,104) heightob(i,j),  0.0, &
                                temptrob(i,j),  0.0
    enddo

  enddo

  ! wrf data formats:
  include 'wrfd_format.inc'

end subroutine wrfd_rass

subroutine wrfd_sondes(numprfls,numlevel,sttnames,obstimes,obsvnlat,obsvnlon, &
                       obsvnelv,reportyp,heightob,pressobs,temptobs,dewtdobs, &
                       uwindobs,vwindobs)

!==============================================================================
!doc  this routine converts sonde data into wrf data format.
!doc
!doc  history:
!doc	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params
  
  implicit none

  character,intent(in) :: reportyp(*)*8,sttnames(*)*5
  integer,  intent(in) :: numprfls,numlevel(*)
  integer,  intent(in) :: obstimes(maxnum_proflrs,*)
  real,     intent(in) :: obsvnlat(maxnum_proflrs,*), &
                          obsvnlon(maxnum_proflrs,*), &
                          heightob(maxnum_proflrs,*), &
                          obsvnelv(maxnum_proflrs)
  real                 :: pressobs(maxnum_proflrs,*), &
                          temptobs(maxnum_proflrs,*), &
                          dewtdobs(maxnum_proflrs,*), &
                          uwindobs(maxnum_proflrs,*), &
                          vwindobs(maxnum_proflrs,*)

  ! local variables:
  character :: obtime*14,varch1*40,varch2*40
  integer   :: i,j,status
  real      :: humidity
  real,parameter :: missng=-888888.000

  ! write out each surface station:
  do i=1,numprfls
    do j=1,numlevel(i)

      obtime = number_asctime
      call make_fnam_lp(obstimes(i,j),varch1,status)
      write(obtime(9:14),2) varch1(6:9),'00'
2     format(a4,a2)
      if (obtime(9:9) .eq. ' ') obtime(9:9) = '0'

      ! write obs time:
      write(output_channel,101) obtime

      ! write lat/lon:
      write(output_channel,102) obsvnlat(i,j),obsvnlon(i,j)

      ! write station name and data name:
      varch1 = '                                         '
      varch2 = '                                         '
      varch1(1:20) = sttnames(i)
      write(output_channel,1021) varch1,'sonde data from laps data ingest'

      ! write elevation etc:
      varch1 = '                                         '
      varch2 = '                                         '
      varch1 = reportyp(i)
      write(output_channel,103) varch1,varch2,heightob(i,j),.true.,.false.,1

      ! convert laps missing to wrf missing:
      if (pressobs(i,j) .eq. rvalue_missing) then
        pressobs(i,j) = missng
      else
        pressobs(i,j) = pressobs(i,j)*100.0	! pascals
      endif
      if (temptobs(i,j) .eq. rvalue_missing) then
        temptobs(i,j) = missng
      else
        if (dewtdobs(i,j) .eq. rvalue_missing) then
          dewtdobs(i,j) = missng
        else
          ! save rh in dewtdobs:
          dewtdobs(i,j) = humidity(temptobs(i,j),dewtdobs(i,j))
        endif
        temptobs(i,j) = (temptobs(i,j) - 32.0)*5.0/9.0+absolu_tmpzero
      endif
      if (uwindobs(i,j) .eq. rvalue_missing) then
        uwindobs(i,j) = missng
      endif
      if (vwindobs(i,j) .eq. rvalue_missing) then
        vwindobs(i,j) = missng
      endif

      ! write obs:
      write(output_channel,104) pressobs(i,j),	0.0, &
                                heightob(i,j),  0.0, &
                                temptobs(i,j),  0.0, &
                                uwindobs(i,j),  0.0, &
                                vwindobs(i,j),  0.0, &
                                dewtdobs(i,j),  0.0

    enddo

  enddo

  ! wrf data formats:
  include 'wrfd_format.inc'

end subroutine wrfd_sondes

subroutine wrfd_sfcobs(numbrobs,obstimes,obsvnlat,obsvnlon, &
                       sttnames,reportyp,prvdname,obsvnelv,mslpress,mslpserr, &
                       refpress,refpserr,temptobs,tmperror,windinuv,uvwnderr, &
                       obsvnrhs,errorrhs,sfcpress,preciptn,errprecp)

!==============================================================================
!doc  this routine converts surface obs into wrf data format.
!doc
!doc  history:
!doc	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character, intent(in) :: sttnames(*)*20, &
                           reportyp(*)*6,prvdname(*)*11
  integer,   intent(in) :: numbrobs
  integer,   intent(in) :: obstimes(*)
  real                  :: obsvnlat(*),obsvnlon(*),obsvnelv(*), &
                           mslpress(*),mslpserr(*), &
                           refpress(*),refpserr(*), &
                           temptobs(*),tmperror(*), &
                           windinuv(2,*),uvwnderr(2,*), &
                           obsvnrhs(*),errorrhs(*), &
                           sfcpress(*),preciptn(*),errprecp(*)

  ! local variables:
  character :: obtime*14,varch1*40,varch2*40
  integer   :: i,status
  real      :: prserr(2)
  real,parameter :: missng=-888888.000

  varch1 = '                                         '
  varch2 = '                                         '

  ! write out each surface station:
  do i=1,numbrobs
    obtime = number_asctime
    write(obtime(9:14),2) obstimes(i),'00'
2   format(i4,a2)
    if (obtime(9:9) .eq. ' ') obtime(9:9) = '0'

    ! write obs time:
    write(output_channel,101) obtime

    ! write lat/lon:
    write(output_channel,102) obsvnlat(i),obsvnlon(i)

    ! write station name and data name:
    varch1(1:20) = sttnames(i)
    write(output_channel,1021) varch1,'all-sfc from laps data ingest'

    ! write elevation etc:
    varch1 = '                                         '
    varch2 = '                                         '
    varch1 = reportyp(i)
    write(output_channel,103) varch1,varch2,obsvnelv(i),.false.,.false.,1

    ! convert laps missing to wrf missing:
    if ((mslpress(i) .eq. rvalue_missing) .or. &
        (mslpress(i) .eq. sfcobs_invalid)) then
      mslpress(i) = missng
      prserr(1) = missng
    else
      mslpress(i) = mslpress(i)*100.0		! to pascals
      prserr(1) = mslpserr(i)*100.0		! to pascals
    endif
    if ((refpress(i) .eq. rvalue_missing) .or. &
        (refpress(i) .eq. sfcobs_invalid)) then
      refpress(i) = missng
      prserr(2) = missng
    else
      refpress(i) = refpress(i)*100.0		! to pascals
      prserr(2) = refpserr(i)*100.0		! to pascals
    endif
    if ((temptobs(i) .eq. rvalue_missing) .or. &
        (temptobs(i) .eq. sfcobs_invalid)) then
      temptobs(i) = missng
      tmperror(i) = missng
    else
      temptobs(i) = (temptobs(i) - 32.0)*5.0/9.0+absolu_tmpzero
      tmperror(i) = tmperror(i)*5.0/9.0
    endif
    if ((windinuv(1,i) .eq. rvalue_missing) .or. &
        (windinuv(1,i) .eq. sfcobs_invalid)) then
      windinuv(1,i) = missng
      uvwnderr(1,i) = missng
    endif
    if ((windinuv(2,i) .eq. rvalue_missing) .or. &
        (windinuv(2,i) .eq. sfcobs_invalid)) then
      windinuv(2,i) = missng
      uvwnderr(2,i) = missng
    endif
    if ((obsvnrhs(i) .eq. rvalue_missing) .or. &
        (obsvnrhs(i) .eq. sfcobs_invalid)) then
      obsvnrhs(i) = missng
      errorrhs(i) = missng
    endif
    if ((sfcpress(i) .eq. rvalue_missing) .or. &
        (sfcpress(i) .eq. sfcobs_invalid)) then
      sfcpress(i) = missng
      mslpserr(i) = missng
    endif
    if ((preciptn(i) .eq. rvalue_missing) .or. &
        (preciptn(i) .eq. sfcobs_invalid)) then
      preciptn(i) = missng
      errprecp(i) = missng
    endif

    ! write obs:
    write(output_channel,105) mslpress(i),	prserr(1), &
                              refpress(i),	prserr(2), &
                              missng,		missng, &
                              temptobs(i),	tmperror(i), &
                              windinuv(1,i),	uvwnderr(1,i), &
                              windinuv(2,i),	uvwnderr(2,i), &
                              obsvnrhs(i),	errorrhs(i), &
                              sfcpress(i),	mslpserr(i), &
                              preciptn(i),	errprecp(i)
  enddo

  ! wrf data formats:
  include 'wrfd_format.inc'

end subroutine wrfd_sfcobs

subroutine wrfd_cdwaca(numberob,obsarray,obsi4dat)
  print*,'please complete this routine!'
  stop
end subroutine wrfd_cdwaca
