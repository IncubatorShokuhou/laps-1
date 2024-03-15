subroutine lapsinfo

!==========================================================
!  this routine gets the necesary laps information.
!
!  history: may. 2004 by yuanfu xie.
!==========================================================

  call get_laps_info(nx,ny,stanlat,stanlat2,stanlon, &
                     laps_cycle_time,badflag ,maxstations,maproj)

  ! debug: print*,'laps: ',nx,ny,stanlat,stanlat2,stanlon, &
  !                   laps_cycle_time,badflag ,maxstations,maproj(1:6)

  l(1:4) = (/mx,my,mt,mv/)

  ! set number of fictitious grid points surrounding the grid:
  nfic = 20

  n(1) = nx+2*nfic
  n(2) = ny+2*nfic
  n(3) = ncycles
  n(4) = nvlaps

  print*,'n = ',n
  print*,'nxy without fictitous: ',nx,ny

  ! qc: threshold value check:
  qc_cons(1) = 10.0		! t
  qc_cons(5) = 10.0		! dt
  qc_cons(2) = 5.0		! u
  qc_cons(3) = 5.0		! v
  qc_cons(4) = 400.0		! msl p
  qc_cons(6) = 400.0		! red p
  
  ! check:
  if ((n(1) .gt. mx) .or. (n(2) .gt. my) .or. &
      (n(3) .gt. mt) .or. (n(4) .gt. mv)) then
     print*,'lapsinfo: analysis array is too small!'
     print*,'x: ',n(1),mx
     print*,'y: ',n(2),my
     print*,'t: ',n(3),mt
     print*,'v: ',n(4),mv
     stop
  endif

  ! read in recursive filter parameters:
  call get_directory('static', dir_s, len)
  name(1:len) = dir_s(1:len)
  name(len+1:len+12) = 'namelist.txt'

  print*,'path to namelist: ',name(1:len+12)

  open(unit=13,file=name(1:len+12),form='formatted')
  do i1=1,4
     read(13,*)
  enddo
  do i1=1,n(4)
     read(13,*) al(1:3,i1),np(1:3,i1)
  enddo
  ! number of minimization iterations:
  read(13,*) maxitr

  ! number of recursive filter iterations:
  read(13,*) nrf(1:n(4))
  close(13)

end subroutine
