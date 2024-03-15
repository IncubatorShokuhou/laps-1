subroutine laps_ingest

!==============================================================================
!doc  this routine ingests observation data by laps libraries.
!doc
!doc  history:
!doc	creation:	yuanfu xie	may 2007
!==============================================================================

  use laps_params

  implicit none

  integer :: istatus

  ! open channels for lapsplot output: tmg, pig, prg, sag etc:
  call open_lapsprd_file(pigout_channel,system_in4time,'pig',istatus)
  call open_lapsprd_file(prgout_channel,system_in4time,'prg',istatus)
  call open_lapsprd_file(sagout_channel,system_in4time,'sag',istatus)
  call open_lapsprd_file(tmgout_channel,system_in4time,'tmg',istatus)

  ! radar:
  ! call read_radar		! prefer to use laps gridded radar data for now yuanfu

  ! profiler:
  call conv_proflr

  ! rass:
  call conv_rass

  ! sonde:
  call conv_sondes

  ! surface obs:
  call conv_sfcobs

  ! cdw and acars:
  call conv_cdwaca

  close(pigout_channel)
  close(prgout_channel)
  close(tmgout_channel)
  close(sagout_channel)

end subroutine laps_ingest

subroutine conv_proflr

!==============================================================================
!doc  this routine reads and converts profiler observation data into data format
!doc  requested.
!doc
!doc  history:
!doc	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character :: fspecs*225, extnsn*31,c5name(maxnum_proflrs)*5, &
               obtype(maxnum_proflrs)*8
  integer   :: iofile,status,nprflr,nlevel(maxnum_proflrs)
  integer   :: i4time,obtime(maxnum_proflrs)
  real	    :: prflat(maxnum_proflrs),prflon(maxnum_proflrs),&
               prfelv(maxnum_proflrs)
  real      :: hghtob(maxnum_proflrs,maxlvl_proflrs)	! height obs
  real      :: uwndob(maxnum_proflrs,maxlvl_proflrs)	! u wind obs
  real      :: vwndob(maxnum_proflrs,maxlvl_proflrs)	! v wind obs
  real      :: rmsobs(maxnum_proflrs,maxlvl_proflrs)	! rms
  real      :: tsfcob(maxnum_proflrs)			! sfc t obs
  real      :: psfcob(maxnum_proflrs)			! sfc p obs
  real      :: rhsfco(maxnum_proflrs)			! sfc rh obs
  real      :: usfcob(maxnum_proflrs)			! sfc u obs
  real      :: vsfcob(maxnum_proflrs)			! sfc v obs

  iofile = 12
  nlevel = 0

  ! open nearest profiler file to the laps analysis time:
  extnsn = 'pro'
  call get_filespec(extnsn,2,fspecs,status)

  call get_file_time(fspecs,system_in4time,i4time)

  ! check if data file is close to system time:
  if (abs(system_in4time-i4time) .gt. length_anatime) then
    print*,'conv_proflr: no recent profiler data files'
    return
  endif

  ! read profiler data: note read_pro_data returns m/s wind:
  call read_pro_data(iofile,i4time,extnsn,maxnum_proflrs,maxlvl_proflrs, & ! i
                     nprflr,nlevel,prflat,prflon,prfelv,c5name,obtime,   & ! o
                     obtype,hghtob,uwndob,vwndob,rmsobs,tsfcob,psfcob,   & ! o
                     rhsfco,usfcob,vsfcob,status)                          ! o
  call laps_divider
  write(6,*) 'conv_proflr: number of profilers data read: ',nprflr
  call laps_divider

  ! convert to the requested format:
  if (nprflr .gt. 0) then
    write(6,*) 'conv_proflr: levels at each profiler: ',nlevel(1:nprflr)
    if      (format_request .eq. 'bufr') then
      call bufr_proflr(nprflr,nlevel,c5name,obtime,prflat,prflon,prfelv, &
                       obtype,maxnum_proflrs,hghtob,uwndob,vwndob,rmsobs, &
                       psfcob,tsfcob,rhsfco,usfcob,vsfcob)
    else if (format_request .ne. 'wrf') then
      call wrfd_proflr(nprflr,nlevel,c5name,obtime,prflat,prflon,prfelv, &
                       obtype,maxnum_proflrs,hghtob,uwndob,vwndob,psfcob, &
                       tsfcob,rhsfco,usfcob,vsfcob)
    endif
  endif

end subroutine conv_proflr

subroutine conv_rass

!==============================================================================
!doc  this routine reads and converts rass observation data into data format
!doc  requested.
!doc
!doc  history:
!doc	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character :: c5name(maxnum_sondes)*5,obtype(maxnum_sondes)*8
  integer   :: status,nrass,nlevel(maxnum_sondes)
  integer   :: i4time,obtime(maxnum_sondes)
  real	    :: prflat(maxnum_sondes),prflon(maxnum_sondes),&
               prfelv(maxnum_sondes)
  real      :: hghtob(maxnum_sondes,maxlvl_sondes)	! height obs
  real      :: tempob(maxnum_sondes,maxlvl_sondes)	! temperature obs
  real      :: errobs(maxnum_sondes)			! temperature error

  ! read profiler data:
  call read_rass_data(system_in4time,length_anatime,maxnum_sondes, &
                     maxlvl_sondes,rvalue_missing,obtime, & ! i
                     nrass,nlevel,c5name,obtype,prflat,prflon,prfelv,errobs, &
                     hghtob,tempob,status)                          ! o
  call laps_divider
  write(6,*) 'conv_rass: number of rass data read: ',nrass
  call laps_divider

  ! convert to the requested format:
  if (nrass .gt. 0) then
    write(6,*) 'conv_rass: levels at each profiler: ',nlevel(1:nrass)
    if      (format_request .eq. 'bufr') then
      call bufr_rass(nrass,nlevel,c5name,obtime,prflat,prflon,prfelv, &
                       obtype,maxnum_sondes,hghtob,tempob,errobs)
    else if (format_request .ne. 'wrf') then
      call wrfd_rass(nrass,nlevel,c5name,obtime,prflat,prflon,prfelv, &
                       obtype,maxnum_sondes,hghtob,tempob)
    endif
  endif

end subroutine conv_rass


subroutine conv_sondes

!==============================================================================
!  this routine reads and converts sonding observation data into data format
!  requested.
!
!  history:
!	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character :: fspecs*225, extnsn*3,c5name(maxnum_sondes)*5, &
               obtype(maxnum_sondes)*8
  integer   :: iofile,status,nprflr,nlevel(maxnum_sondes), &
               window,modeob
  integer   :: i4time,prftim(maxnum_sondes,maxlvl_sondes),i
  real	    :: prflat(maxnum_sondes,maxlvl_sondes), &
               prflon(maxnum_sondes,maxlvl_sondes), &
               prfelv(maxnum_sondes)
  real      :: hghtob(maxnum_sondes,maxlvl_sondes)	! height obs
  real      :: prsobs(maxnum_sondes,maxlvl_sondes)	! pressure obs
  real      :: uwndob(maxnum_sondes,maxlvl_sondes)	! u wind obs
  real      :: vwndob(maxnum_sondes,maxlvl_sondes)	! v wind obs
  real      :: tempob(maxnum_sondes,maxlvl_sondes)	! t obs
  real      :: dewobs(maxnum_sondes,maxlvl_sondes)	! dew point obs

  iofile = 12
  nlevel = 0
  window = 0			! i4time window mimic read_profiles.f
  modeob = 3			! key levels off of wind data

  ! open nearest profiler file to the laps analysis time:
  extnsn = 'snd'
  call get_filespec(extnsn,2,fspecs,status)

  call get_file_time(fspecs,system_in4time,i4time)
  if (abs(system_in4time-i4time) .gt. length_anatime) then
    write(6,*) 'conv_sondes: warning: nearest sonde file is outside window'
  else

    ! read sonde data:
    call read_snd_data2(iofile,i4time,extnsn,maxnum_sondes,maxlvl_sondes, & ! i
                       domain_latitde,domain_longitd,domain_topogrp,number_gridpts(1),     & ! i
                       number_gridpts(2),number_gridpts(3),height_grid3dm,  & ! i
                       .true.,modeob,                                       & ! i
                       nprflr,prfelv,nlevel,c5name,obtype,hghtob,           & ! o
                       prsobs,uwndob,vwndob,tempob,dewobs,prflat,prflon,    & ! o
                       prftim,status)                                         ! o
    call laps_divider
    write(6,*) 'conv_sondes: number of sonde data read: ',nprflr
    call laps_divider

    ! convert to the requested format:
    if (status .eq. 1 .and. nprflr .gt. 0) then
      write(6,*) 'conv_sondes: levels at each sonde: '
      do i=1,nprflr
        if (nlevel(i) .gt. 0) write(6,*) i,': ',nlevel(i)
      enddo
      if      (format_request .eq. 'bufr') then
        call bufr_sondes(nprflr,nlevel,c5name,prftim,prflat,prflon,prfelv, &
                         obtype,hghtob,prsobs,tempob,dewobs,uwndob,vwndob)
      else if (format_request .ne. 'wrf') then
        call wrfd_sondes(nprflr,nlevel,c5name,prftim,prflat,prflon,prfelv, &
                         obtype,hghtob,prsobs,tempob,dewobs,uwndob,vwndob)
      endif
    endif

  endif

end subroutine conv_sondes

subroutine conv_sfcobs

!==============================================================================
!  this routine reads and converts surface observations into requested data
!  format.
!
!  history:
!	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  character*24 :: filetm		! obs file time
  integer :: maxstn,status,i		! maximum stations
  integer :: nobgrd,nobbox		! number of obs over grid and box
  integer :: nss,nsngrd,nsnbox		! number of snd obs over grid and box

  character*20,allocatable,dimension(:) :: stname	! station names
  character*11,allocatable,dimension(:) :: pvname	! provider names
  character*25,allocatable,dimension(:) :: prstwx	! present weather
  character*6, allocatable,dimension(:) :: rptype	! report type
  character*6, allocatable,dimension(:) :: stntyp	! station type
  integer,     allocatable,dimension(:) :: obtime,wmoids! obs time/wmo id
  integer,     allocatable,dimension(:) :: cldlyr,prschc! cloud layer/prs chg
  real,        allocatable,dimension(:) :: &
    obslat,obslon,obselv,obstmp,errtmp,obsdew,errdew,obsrhs,errrhs, &
    obsdir,errdir,obsspd,errspd,gusdir,gusspd,obsalt,erralt,stnprs, &
    mslprs,prsch3,errprs,obsvis,errvis,obssol,errsol,sfctmp,errsft, &
    sfcmoi,errsfm,precp1,precp3,precp6,prec24,errpcp,snowcv,errsnw, &
    maxtmp,mintmp !,igrid,jgrid

  character*4, allocatable,dimension(:,:) :: cldamt     ! cloud amount

  real,        allocatable,dimension(:,:) :: cldhgt	! cloud heights
  real,        allocatable,dimension(:,:) :: obswnd, &  ! obs wind
                                             errwnd	! obs wind error

  ! get maximum number of surface stations:
  call get_maxstns(maxstn,status)
  if (status .ne. 1) then
    write(6,*) 'conv_sondes: error in reading maximum sfc stations'
    stop
  endif

  ! add number of sonde surface:
  maxstn = maxstn+maxnum_sondes

  ! allocatable memory for surface variables:
  allocate(obtime(maxstn),wmoids(maxstn),stname(maxstn), &
           pvname(maxstn),prstwx(maxstn),rptype(maxstn), &
           stntyp(maxstn),obslat(maxstn),obslon(maxstn), &
           obselv(maxstn),obstmp(maxstn),obsdew(maxstn), &
           obsrhs(maxstn),obsdir(maxstn),obsspd(maxstn), &
           gusdir(maxstn),gusspd(maxstn),obsalt(maxstn), &
           stnprs(maxstn),mslprs(maxstn),prschc(maxstn), &
           prsch3(maxstn),obsvis(maxstn),obssol(maxstn), &
           sfctmp(maxstn),sfcmoi(maxstn),precp1(maxstn), &
           precp3(maxstn),precp6(maxstn),prec24(maxstn), &
           snowcv(maxstn),cldlyr(maxstn),maxtmp(maxstn), &
           mintmp(maxstn),errtmp(maxstn),errdew(maxstn), &
           errrhs(maxstn),errdir(maxstn),errspd(maxstn), &
           erralt(maxstn),errprs(maxstn),errvis(maxstn), &
           errsol(maxstn),errsft(maxstn),errsfm(maxstn), &
           errpcp(maxstn),errsnw(maxstn), & ! igrid(maxstn),jgrid(maxstn), &
           cldamt(maxstn,5),cldhgt(maxstn,5), &
           obswnd(2,maxstn),errwnd(2,maxstn), stat=status)
  if (status .ne. 0) then
    write(6,*) 'conv_sfcobs: error in allocating memory for surface data'
    stop
  endif

  ! read sfc obs: read_surface_data returns knots for wind:
  call read_surface_data(system_in4time,filetm,nobgrd,nobbox,obtime,wmoids,&
                         stname,pvname,prstwx,rptype,stntyp,obslat,obslon, &
                         obselv,obstmp,obsdew,obsrhs,obsdir,obsspd,gusdir, &
                         gusspd,obsalt,stnprs,mslprs,prschc,prsch3,obsvis, &
                         obssol,sfctmp,sfcmoi,precp1,precp3,precp6,prec24, &
                         snowcv,cldlyr,maxtmp,mintmp,                      &
                         errtmp,errdew,errrhs,errdir,errspd,erralt,errprs, &
                         errvis,errsol,errsft,errsfm,errpcp,errsnw,cldamt, &
                         cldhgt,maxstn,status)

  if (status .ne. 1) then
    write(6,*) 'conv_sfcobs: error in reading surface obs'
    stop
  endif

  ! read sonde sfc obs: read_surface_data returns knots for wind:
  nss = nobbox+1
  nsngrd = 0		! note: read_sfc_snd does not initialize nsngrd
  nsnbox = 0		! note: read_sfc_snd does not initialize nsnbox
  ! the following call command allows snd surface data ingested but these data
  ! are ingested from sonde reader already. 
  ! call read_sfc_snd(system_in4time,filetm,nsngrd,nsnbox,obtime(nss),wmoids(nss),&
  !                      stname(nss),pvname(nss),prstwx(nss),rptype(nss), &
  !                      stntyp(nss),obslat(nss),obslon(nss),obselv(nss), & 
  !                      obstmp(nss),obsdew(nss),obsrhs(nss),obsdir(nss), &
  !                      obsspd(nss),gusdir(nss),gusspd(nss),obsalt(nss), & 
  !                      stnprs(nss),mslprs(nss),prschc(nss),prsch3(nss), &
  !                      obsvis(nss),obssol(nss),sfctmp(nss),sfcmoi(nss), &
  !                      precp1(nss),precp3(nss),precp6(nss),prec24(nss), &
  !                      snowcv(nss),cldlyr(nss),maxtmp(nss),mintmp(nss), &
  !                      errtmp(nss),errdew(nss),errrhs(nss),errdir(nss), &
  !                      errspd(nss),erralt(nss),errprs(nss),errvis(nss), &
  !                      errsol(nss),errsft(nss),errsfm(nss),errpcp(nss), &
  !                      errsnw(nss),cldamt(nss,1),cldhgt(nss,1),maxstn, &
  !                      domain_latitde,domain_longitd,number_gridpts(1), &
  !                      number_gridpts(2),number_gridpts(3),maxnum_sondes, &
  !                      maxlvl_sondes,domain_topogrp,status)
  !
  ! if (nsnbox .lt. 1) then
  !   write(6,*) 'conv_sfcobs: no sonde sfc data available'
  ! endif
  ! account total sfc data:
  ! the following command allows snd surface data ingested but these data
  ! are ingested from sonde reader already. 
  ! nobbox = nobbox+nsnbox

  ! convert wind direction and speed into u and v:
  do i=1,nobbox
    if ((obsdir(i) .ne. rvalue_missing) .and. &
        (obsdir(i) .ne. sfcobs_invalid) .and. &
        (obsspd(i) .ne. rvalue_missing) .and. &
        (obsspd(i) .ne. sfcobs_invalid)) then
      call disp_to_uv(obsdir(i),obsspd(i),obswnd(1,i),obswnd(2,i))
      call disp_to_uv(errdir(i),errspd(i),errwnd(1,i),errwnd(2,i))

      ! convert wind into m/s:
      obswnd(1:2,i) = obswnd(1:2,i)*0.514791
    else
      obswnd(1:2,i) = rvalue_missing
      errwnd(1:2,i) = rvalue_missing
    endif

    ! temperature and dew:
    if ((obstmp(i) .ne. rvalue_missing) .and. &
        (obstmp(i) .ne. sfcobs_invalid) ) obstmp(i) = (obstmp(i)-32.0)*5.0/9.0
    if ((obsdew(i) .ne. rvalue_missing) .and. &
        (obsdew(i) .ne. sfcobs_invalid) ) obsdew(i) = (obsdew(i)-32.0)*5.0/9.0
  enddo

  call laps_divider
  write(6,*) 'conv_sfcobs: number of surface obs data read: ',nobbox,nobgrd
  call laps_divider

  ! convert to the requested data format:
  if      (format_request .eq. 'bufr') then
    call bufr_sfcobs(nobbox,obtime, &
                     obslat,obslon,stname,rptype,pvname,obselv, &
                     mslprs,errprs,stnprs,errprs,obstmp,errtmp, &
                     obswnd,errwnd,obsrhs,errrhs,stnprs,precp1,errpcp)
  else if (format_request .eq. 'wrf' ) then
    call wrfd_sfcobs(nobbox,obtime, &
                     obslat,obslon,stname,rptype,pvname,obselv, &
                     mslprs,errprs,stnprs,errprs,obstmp,errtmp, &
                     obswnd,errwnd,obsrhs,errrhs,stnprs,precp1,errpcp)
  endif

  ! deallocatable memory for surface variables:
  deallocate(obtime,wmoids,stname, &
             pvname,prstwx,rptype, &
             stntyp,obslat,obslon, &
             obselv,obstmp,obsdew, &
             obsrhs,obsdir,obsspd, &
             gusdir,gusspd,obsalt, &
             stnprs,mslprs,prschc, &
             prsch3,obsvis,obssol, &
             sfctmp,sfcmoi,precp1, &
             precp3,precp6,prec24, &
             snowcv,cldlyr,maxtmp, &
             mintmp,errtmp,errdew, &
             errrhs,errdir,errspd, &
             erralt,errprs,errvis, &
             errsol,errsft,errsfm, &
             errpcp,errsnw, &
             cldamt,cldhgt, &
             obswnd,errwnd, stat=status)
  if (status .ne. 0) then
    write(6,*) 'conv_sfcobs: error in deallocating memory for surface obs'
    stop
  endif

end subroutine conv_sfcobs

subroutine conv_cdwaca

!==============================================================================
!  this routine reads and converts cloud drift wind and acars (pirep) data into
!  requested data format
!
!  history:
!	creation:	yuanfu xie	jun 2007
!==============================================================================

  use laps_params

  implicit none

  ! local variables:
  integer, parameter :: maxobs=300000	! locally defined,change it if needed
  character :: extend*3,obstyp*4,asctim*9
  integer   :: numobs,status
  integer   :: obsint(3,maxobs)
  logical   :: geoalt
  real      :: obslat,obslon,obselv,obsprs,obsdir,obsspd,obstmp,obsnon,ztopsa
  real      :: ucmpnt,vcmpnt,ri,rj,height_to_pressure
  real*8    :: obarry(7,maxobs)

  ! cloud drift wind:
  numobs = 0
  ! acar (pirep):
  extend = 'pin'

  ! 1. temp:
  obstyp = 'temp'
  call open_lapsprd_file_read(points_channel,system_in4time,extend,status)
  if (status .ne. 1) then
    print*,'conv_cdwaca: no acar (pirep) temp data'
  else
    status = 0
    do
      call read_acars_ob(points_channel,obstyp,obslat,obslon,obselv,obstmp,obsnon, &
                               asctim,0,geoalt,status)
      if (status .ne. 0) exit
      numobs = numobs+1
      if (numobs .gt. maxobs) then
        print*,'conv_cdw_aca: data array is too small, maxobs needs enlarge ',numobs,maxobs
        stop
      endif
      obarry(1:7,numobs) = rvalue_missing
      obarry(1,numobs) = obslat
      obarry(2,numobs) = obslon
      if (geoalt) then
        ! geometric height:
        obarry(3,numobs) = obselv
      else
        ! pressure height in standard atmosphere:
        obarry(4,numobs) = ztopsa(obselv)	! ztopsa returns pressure in mb
      endif
      obarry(7,numobs) = obstmp-273.15	! read_acars_ob returns kelvin
      obsint(1,numobs) = 41		! use bufr report type code: 
                                        ! 41 acars temperature
      if (geoalt) then 
        obsint(2,numobs) = 666		! temporarily use 666 for widsom data
      else
        obsint(2,numobs) = 130		! pilor report: temperature
      endif
      call cv_asc_i4time(asctim,obsint(3,numobs))
      !print*,'acar temp: ',obslat,obslon,obselv,obstmp,obsint(3,numobs)
    enddo
  endif
  close(points_channel)
  print*,'conv_cdwaca: total of cdw + acar temp obs: ',numobs

  ! 1. wind:
  obstyp = 'wind'
  call open_lapsprd_file_read(points_channel,system_in4time,extend,status)
  if (status .ne. 1) then
    print*,'conv_cdwaca: no acar (pirep) wind data'
  else
    status = 0
    do
      call read_acars_ob(points_channel,obstyp,obslat,obslon,obselv,obsdir,obsspd, &
                               asctim,0,geoalt,status)
      if (status .ne. 0) exit
      numobs = numobs+1
      if (numobs .gt. maxobs) then
        print*,'conv_cdwaca: data array is too small, maxobs needs enlarge ',numobs,maxobs
        stop
      endif
      obarry(1:7,numobs) = rvalue_missing
      obarry(1,numobs) = obslat
      obarry(2,numobs) = obslon
      if (geoalt) then 
        ! geometric height:
        obarry(3,numobs) = obselv
      else
        ! pressure height in standard atmosphere:
        obarry(4,numobs) = ztopsa(obselv)	! ztopsa returns pressure in mb
      endif
      call disp_to_uv(obsdir,obsspd,obstmp,obsnon)
      obarry(5,numobs) = obstmp !obsdir	! when calling uvtrue_to_uvgrid
      obarry(6,numobs) = obsnon !obsspd
      obsint(1,numobs) = 41		! use bufr report type code: 
                                        ! acars wind
      if (geoalt) then
        obsint(2,numobs) = 666		! temporarily use 666 for wisdom data
      else
        obsint(2,numobs) = 230		! pilor report: wind
      endif
      call cv_asc_i4time(asctim,obsint(3,numobs))
      ! print*,'acar wind: ',obslat,obslon,obselv,obsdir,obsspd,obsint(3,numobs)
    enddo
  endif
  close(points_channel)

  extend = 'cdw'
  call open_lapsprd_file_read(points_channel,system_in4time,extend,status)
  if (status .ne. 1) then
    print*,'conv_cdwaca: no cloud drift wind data'
  else
    status = 0
    do
      call read_laps_cdw_wind(points_channel,obslat,obslon,obsprs,obsdir,obsspd, &
                               asctim,status)
      if (status .ne. 0) exit
      numobs = numobs+1
      if (numobs .gt. maxobs) then
        print*,'conv_cdw_aca: data array is too small, maxobs needs enlarge ',numobs,maxobs
        stop
      endif
      ! print*,'conv_cdwaca: found cdw data, complete this code'

      obarry(1:7,numobs) = rvalue_missing
      obarry(1,numobs) = obslat
      obarry(2,numobs) = obslon
      obarry(4,numobs) = obsprs/100.0	! bufr pob in mb

      ! convert dir/spd into u/v for bufr
      call disp_to_uv(obsdir,obsspd,ucmpnt,vcmpnt)
      obarry(5,numobs) = ucmpnt
      obarry(6,numobs) = vcmpnt

      obsint(1,numobs) = 63		! use bufr report type code: 
                                        ! 63 satellite derived wind
      obsint(2,numobs) = 241		! rerort type: satwind
      call cv_asc_i4time(asctim,obsint(3,numobs))
      ! print*,'cdw: ',obslat,obslon,obsprs,obsdir,obsspd,obsint(3,numobs)
    enddo
  endif
  close(points_channel)
  print*,'conv_cdwaca: total of cdw obs: ',numobs

  print*,'conv_cdwaca: total of cdw + acar temp + acar wind obs: ',numobs

  if (numobs .eq. 0) return

  if (format_request .eq. 'bufr') then
    call bufr_cdwaca(numobs,obarry,obsint)
  else if (format_request .eq. 'wrf' ) then
    call wrfd_cdwaca(numobs,obarry,obsint)
  endif

end subroutine conv_cdwaca


subroutine read_radar

!==============================================================================
!doc  this routine reads multiple radar radial wind velocity using laps'
!doc  get_multiradar_vel.
!doc
!doc  history:
!doc	creation:	yuanfu xie	mar 2008
!==============================================================================

  use laps_params
  use mem_namelist		! laps wind parameter module

  implicit none

  ! local variables:
  character*31 :: radext(max_radars)	! possible radar name extensions
  character*4  :: radnam(max_radars)	! radar station names
  integer      :: nradar		! number of radar available
  integer      :: lcycle		! laps cycle time
  integer      :: sttrad,sttnqy 	! radar and its nyquist status
  integer      :: ngrdrd(max_radars) 	! number of gridpoints with measurable vel
  integer      :: ngrdrd_old(max_radars) 	! number of gridpoints with measurable vel
  integer      :: radtim(max_radars)	! radar observation time
  integer      :: radids(max_radars) 	! radar ids
  integer      :: i,j,k,l
  logical      :: cluttr		! .true. -- remove 3d radar clutter
  real         :: radvel_old(number_gridpts(1),number_gridpts(2),number_gridpts(3),max_radars)
  real         :: radvel(number_gridpts(1),number_gridpts(2),number_gridpts(3),max_radars)
                  ! radar 4d velocity grid
  real         :: radnqy(number_gridpts(1),number_gridpts(2),number_gridpts(3),max_radars)
                  ! radar 4d nyquist velocity
  real         :: uvzero(number_gridpts(1),number_gridpts(2),number_gridpts(3),2)
                  ! zero uv grids used for calling laps qc_radar_obs
  real         :: uvbkgd(number_gridpts(1),number_gridpts(2),number_gridpts(3),2)
                  ! uv background grids used for calling laps qc_radar_obs
  real         :: uv4dml(number_gridpts(1),number_gridpts(2),number_gridpts(3),-1:1,2)
                  ! uv background grids used for calling laps qc_radar_obs
  real         :: volnqy(max_radars)		! volume nyquist velocity
  real         :: radlat(max_radars),radlon(max_radars),radhgt(max_radars)
  real         :: uvgrid(2)

  include 'main_sub.inc'

  cluttr = .true.		! true. -- remove 3d radar clutter
  radars_timetol = 900
  call get_multiradar_vel(system_in4time,radars_timetol,radtim,max_radars, &
                           nradar,radext,rvalue_missing,cluttr, &
                           number_gridpts(1),number_gridpts(2),number_gridpts(3), &
                           radvel,radnqy,radids,volnqy,ngrdrd,radlat,radlon,radhgt, &
                           radnam,sttrad,sttnqy)

  ! set uv zero grids:
  uvzero = 0.0

  ! get laps background:
  call get_laps_cycle_time(lcycle,sttrad)
  call get_fg_wind_new(system_in4time,lcycle,number_gridpts(1),number_gridpts(2), &
                        number_gridpts(3),-1,1,uv4dml(1,1,1,-1,1),uv4dml(1,1,1,-1,2), &
                        uvbkgd(1,1,1,1),uvbkgd(1,1,1,2),sttrad)

  ! convert to grid north from true north:
  if (( .not. l_grid_north_bkg) .and. l_grid_north_anal) then
    write(6,*) ' rotating first guess (background) to grid north'

    do k=1,number_gridpts(3)
      do j=1,number_gridpts(2)
        do i=1,number_gridpts(1)
          call uvtrue_to_uvgrid(uvbkgd(1,1,1,1),uvbkgd(1,1,1,2), &
                                 uvgrid(1),uvgrid(2),domain_longitd(i,j))
          uvbkgd(i,j,k,1:2) = uvgrid(1:2)
        enddo
      enddo
    enddo
  endif ! end of conversion to grid north


  radvel_old = radvel
  ngrdrd_old = ngrdrd

  ! nyquist unfolding:
!!! unfinished: need more work apr. 2008
  do l=1,nradar
    ! qc and unfolding radar nyquist:
!    call qc_radar_obs(number_gridpts(1),number_gridpts(2),number_gridpts(3), &
!                       rvalue_missing,radvel(1,1,1,l),radnqy(1,1,1,l),ngrdrd(l), &
!                       domain_latitde,domain_longitd,radlat(l),radlon(l),radhgt(l), &
!                       uvzero(1,1,1,1),uvzero(1,1,1,2),uvbkgd(1,1,1,1),uvbkgd(1,1,1,2), &
!                       volnqy(l),l_correct_unfolding,l_grid_north,sttrad)
    print*,'status of qc_radar_obs: ',sttrad

  if (ngrdrd_old(l) .ne. ngrdrd(l)) print*,'number grid radar change: ',ngrdrd_old(l),ngrdrd(l)
    do k=1,number_gridpts(3)
      do j=1,number_gridpts(2)
        do i=1,number_gridpts(1)
          if (radvel(i,j,k,l) .ne. rvalue_missing) print*,'radial: ', &
	      radvel(i,j,k,l),radnqy(i,j,k,l),i,j,k,l,ngrdrd(l),volnqy(l)
          if (radvel_old(i,j,k,l) .ne. radvel(i,j,k,l)) print*,'unfolded: ', &
              radvel_old(i,j,k,l),radvel(i,j,k,l)
        enddo
      enddo
    enddo

  enddo

  stop

end subroutine read_radar
