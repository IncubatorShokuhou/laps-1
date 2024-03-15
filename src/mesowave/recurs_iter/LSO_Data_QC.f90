subroutine lso_data

!==========================================================
!  this routine reads the lso from laps.
!
!  history: may. 2004 by yuanfu xie.
!==========================================================

  implicit none

  real, parameter :: temp0 = 273.16
  real :: badsfc

  integer :: maxtime,mintime,hour,minute,seconds,midnight
  integer :: imx,imm,mem_error,i,j,nsts,iobs
  integer, allocatable, dimension(:) :: idsts
  integer, allocatable, dimension(:,:) :: naccn
  real    :: vmx,vmm
  integer, allocatable, dimension(:,:,:) :: numob
  real,    allocatable, dimension(:,:,:) :: stobs

  ! get the value for bad surface data:
  call get_sfc_badflag(badsfc,istatus)

  call get_laps_info(nx,ny,stanlat,stanlat2,stanlon, &
                     laps_cycle_time,badflag,maxstations,maproj)

  print*,'start to allocate space'

  allocate(lat(nx,ny), lon(nx,ny), ldf(nx,ny),stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: cannot allocate space for lat/lon/ldf'
     stop
  endif
  allocate (olaps(nvlaps,maxstations*ncycles), &
            olat(maxstations*ncycles),olon(maxstations*ncycles), &
            otime(maxstations*ncycles),wght(maxstations*ncycles), &
            stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: cannot allocate space for olaps,olat,olon,otime,wght'
     stop
  endif

  ! space for background fields:
  allocate(bkgd(nx,ny,ncycles,nvlaps),stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data: cannot allocate space for background fields'
     stop
  endif

  print*,'space allocated'
	
  ext_s = 'nest7grid'
  var_s = 'lat'
  call rd_laps_static(dir_s,ext_s,nx,ny,1,var_s,units, &
                      comment,lat ,grid_spacingy,istatus)
  ! debug: 
  print*,'gridspace y: ',grid_spacingy

  var_s = 'lon'
  call rd_laps_static(dir_s,ext_s,nx,ny,1,var_s,units, &
                      comment,lon ,grid_spacingx,istatus)
  ! debug: 
  print*,'gridspace x: ',grid_spacingx

  var_s = 'ldf'
  call rd_laps_static(dir_s,ext_s,nx,ny,1,var_s,units, &
                      comment,ldf ,grid_spacingx,istatus)
  do j=1,ny
     do i=1,nx
	if (ldf(i,j) .gt. 1.0) ldf(i,j) = 1.0
	if (ldf(i,j) .lt. 0.0) ldf(i,j) = 0.0
     enddo
  enddo

  ! clear unused portion of the array:
  olaps = badflag
  otime = -1.0   ! to search the maximum hour
  bkgd = 0.0
  call lso_reader_meso(maxstations,nvlaps,ncycles, &
       laps_cycle_time,badflag,olaps,wght,olat,olon, &
       otime,istarttime,nx,ny,bkgd)

  ! debug: print*,'max: ',maxstations,laps_cycle_time
  ! debug: print*,'obs: ',olaps(1:6,1),lat(1,1),lon(1,1)
  ! debug: print*,'obs: ',olaps(1:6,2),lat(nx,ny),lon(nx,ny)

  ! convert john's data to iterative recursive filter data:
  dm(1,1:2) = 1.0-nfic
  dm(2,1) = float(nx)+nfic
  dm(2,2) = float(ny)+nfic

  ! midnight:
  midnight = 0
  if (maxval(otime) .ge. 2300.00) midnight = 1

  ! space for qc analysis: 
  allocate(idsts(ncycles*maxstations*nvlaps), stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error in allocating space for idsts'
     stop
  endif

  nobs = 0
  maxtime = -100
  mintime = 100000000
  do ilaps=1,maxstations*ncycles
     do ivar=1,nvlaps
	if ((olaps(ivar,ilaps) .ne. badflag) .and. &
            (olaps(ivar,ilaps) .ne. badsfc)) then

 	   ! valid obs:
	   nobs = nobs+1
	   if (nobs .gt. mobs) then
	      print*,'lso_data: too many obs!'
	      stop
	   endif

	   ! separating hours and minutes:
	   hour = otime(ilaps)/100
	   minute = otime(ilaps)-hour*100

	   ! handle data cross midnight:
 	   if ((midnight .eq. 1) .and. (hour .le. 2)) hour = 24+hour

	   ! debug: print*,'oll: ',ivar,olaps(ivar,ilaps), &
           ! olat(ilaps),olon(ilaps),otime(ilaps),hour,minute

	   seconds = hour*3600.0+minute*60.0

	   if (maxtime .lt. seconds) maxtime = seconds
	   if (mintime .gt. seconds) mintime = seconds

           vid(nobs) = ivar
	   o(1,nobs) = olaps(ivar,ilaps)

	   ! convert f to k:
	   if ((ivar .eq. 1) .or. (ivar .eq. 5)) &
	      o(1,nobs) = (o(1,nobs)-32.0)*5.0/9.0+temp0

           ! convert altimeter to pascal:
           if ((ivar .eq. 4) .or. (ivar .eq. 6)) &
              o(1,nobs) = o(1,nobs)*100.0

           call latlon_to_rlapsgrid(olat(ilaps),olon(ilaps),lat,lon, &
                                 nx,ny,o(2,nobs),o(3,nobs),istatus)
	   ! debug: print*,'grid loc: ',o(2,nobs),o(3,nobs),nobs
           o(4,nobs) = seconds

	   ! weight:
	   w(nobs) = 1.0
	endif
     enddo
  enddo
  print*,'max time: ',maxtime/3600,mod(maxtime,3600)/60
  print*,'min time: ',mintime/3600,mod(mintime,3600)/60
  print*,'istarttime: ',istarttime,mod(istarttime,86400)/3600.0
  if ((maxtime-mintime)/3600 .gt. 3) then
     print*,'lso_data_qc: it is hard coded to run'
     print*,' the analysis less than 3 hour!'
     stop
  endif

  print*,'total number of obs: ',nobs

  ! qc:
  ! 1. find out number of different stations and identical ones:
  nsts = 0
  do iobs=1,nobs
     do i=1,iobs-1
        if ((abs(o(2,i)-o(2,iobs)) .lt. 1.0e-3) .and. &
            (abs(o(3,i)-o(3,iobs)) .lt. 1.0e-3)) then
           idsts(iobs) = idsts(i)
           goto 11
        endif
     enddo
     nsts = nsts+1
     idsts(iobs) = nsts
11   continue
  enddo
  print*,'number of obs stations: ',nsts
  
  allocate(stobs(ncycles*15,nvlaps,nsts), stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error to allocate space for stobs'
     stop
  endif
  allocate(numob(ncycles*15,nvlaps,nsts), stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error to allocate space for numob'
     stop
  endif
  allocate(naccn(nvlaps,nsts), stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error to allocate space for naccn'
     stop
  endif
  naccn = 0
  do iobs=1,nobs
     naccn(vid(iobs),idsts(iobs)) = naccn(vid(iobs),idsts(iobs))+1
     numob(naccn(vid(iobs),idsts(iobs)),vid(iobs),idsts(iobs)) = iobs
     stobs(naccn(vid(iobs),idsts(iobs)),vid(iobs),idsts(iobs)) = o(1,iobs)
  enddo

  print*,'naccn(1,1) = ',naccn(1,1),stobs(1:naccn(1,1),1,1)
  do iobs=1,naccn(1,1)
     print*,'first: ',o(2:4,numob(iobs,1,1)),olat(numob(iobs,1,1)), &
                                             olon(numob(iobs,1,1))
  enddo

  ! time interval:
  dm(1,3) = mod(istarttime,86400)  	! second of the time
  dm(2,3) = dm(1,3)+(ncycles-1)*laps_cycle_time

  ! spacings: x, y and t
  d(1:3) = (dm(2,1:3)-dm(1,1:3))/float(n(1:3)-1)

  print*,'analysis start time: ',int(dm(1,3))/3600,mod(int(dm(1,3)),3600)/60
  print*,'analysis endng time: ',int(dm(2,3))/3600,mod(int(dm(2,3)),3600)/60

  ! check dimension for variables:
  if (nvlaps .gt. mv) then
     print*,'readobsn: too much variables'
     stop
  endif
  ! check if number of obs is fit into the array:
  if (nobs .gt. mobs) then
     print*,'readobsn: too many obs'
     stop
  endif
  if (nobs .le. 0) then
     print*,'readobsn: no obs to analyze'
     stop
  endif

  print*,'spacing: ',d

  ! write surface.dat for testing:
  ! open(unit=10,file='surface.dat',form='formatted')

  ! domain info:
  print*,'domain: ', dm(1,1),dm(2,1),dm(1,2),dm(2,2),dm(1,3),dm(2,3)
  vmx = -1000.0
  vmm = 1000.0
  do ilaps=1,nobs

     if (vid(ilaps) .eq. 2) then
	if (vmx .lt. o(1,ilaps)) then
	   vmx = o(1,ilaps)
	   imx = ilaps
        endif
        if (vmm .gt. o(1,ilaps)) then
	   vmm = o(1,ilaps)
	   imm = ilaps
        endif
     endif
     ! write(10,1) vid(ilaps),o(1:4,ilaps),w(ilaps)

  enddo
1 format(i2,5e14.6)
  ! close(10)
  print*,'max/min u obs: ',vmx,vmm,imx,imm,w(imm)

  print*,'time interval: ',minval(o(4,1:nobs)),maxval(o(4,1:nobs))

  ! deallocate the memory:
  deallocate(lat,lon,stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error in deallocating lat/lon'
     stop
  endif
  deallocate(olaps,olat,olon,otime,wght,stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error in deallocating olaps,olat,olon,otime,wght'
     stop
  endif
  
  deallocate(stobs, stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error to deallocate space for stobs'
     stop
  endif
  deallocate(numob, stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error to deallocate space for numob'
     stop
  endif
  deallocate(naccn, stat=mem_error)
  if (mem_error .ne. 0) then
     print*,'lso_data_qc: error to deallocate space for naccn'
     stop
  endif
  deallocate(idsts,stat=mem_error)

end subroutine lso_data
