!dis    forecast systems laboratory
!dis    noaa/oar/erl/fsl
!dis    325 broadway
!dis    boulder, co     80303
!dis
!dis    forecast research division
!dis    local analysis and prediction branch
!dis    laps
!dis
!dis    this software and its documentation are in the public domain and
!dis    are furnished "as is."  the united states government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  they assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  all modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  if significant modifications or enhancements
!dis    are made to this software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

program stmas_mg

!==========================================================
!  this program is to analyze surface data using time and
!  space observation through a multigrid technique.
!
!  history: 
! 	creation: yuanfu xie	6-2005
!==========================================================

  use definition
  use lapsdatsrc
  use memorymngr
  use prepostprc
  use stmasanalz

  implicit none

  integer :: i,j,k,mi,mj,mk,nnn
  real :: a

  ! namelist:
  call prpstnls

  ! laps parameters:
  call lapsinfo

  ! allocate memory for laps and interactive vars:
  call lapsmemo		! for laps usage
  call intrmemo		! for interactive variables

  ! laps grid configuration:
  call lapsconf

  ! background fields:
  call lapsbkgd

  ! surface lso obs:
  call lapsobsv(mxstts)

  ! convert to units consistent with background:
  call lapsunit

  ! stmas qc:
  call date_and_time(date0,time0,zone0,timing0)
  call laps_qcs
  call date_and_time(date1,time1,zone1,timing1)
  print*,'time used for laps_qcs: ',(timing1(6)-timing0(6))+(timing1(7)-timing0(7))/60.0
  print*,'times used for laps_qcs: ',timing1(6),timing0(6),timing1(7),timing0(7)

  ! add background to obs where obs is spare !by min-ken hseih
  call date_and_time(date0,time0,zone0,timing0)
  ! call addbkgrd
  call jbgridpt
  ! uncovr = .false.
  call date_and_time(date1,time1,zone1,timing1)
  print*,'time used for addbkgrd: ',(timing1(6)-timing0(6))+(timing1(7)-timing0(7))/60.0
  print*,'times used for addbkgrd: ',timing1(6),timing0(6),timing1(7),timing0(7)

  ! release memory for laps:
  call lapsrels

  ! allocate memory for stmas:
  call stmasmem

  ! test analytic function:
  ! call stmastst

  ! stmas analyses:
  do i=1,numvar
    write(*,*) 'stmas_main: start analyzing ',varnam(i)
    ! check if there is any obs for analysis:
    if (numobs(i) .gt. 0) then
      ! modified by min-ken hsieh
      ! pass stanam, varnam into stmasana
      !
      call date_and_time(date0,time0,zone0,timing0)
      call stmasana(analys(1,1,1,i),numgrd,grdspc, &
	domain,bkgrnd(1,1,1,i),numtmf, &
	qc_obs(1,1,i),numobs(i),weight(1,i), stanam(1,i),&
	obsspc(1,i),indice(1,1,i),coeffs(1,1,i), &
	bounds(i),stmasi,stmasr,varnam(i),pnlt_v(i), &
        slevel(i),uncovr(1,1,1,i),diagnl(1,1,i))
      call date_and_time(date1,time1,zone1,timing1)
      print*,'time used for stmasana ',i,(timing1(6)-timing0(6))+(timing1(7)-timing0(7))/60.0
      print*,'times used for stmasana: ',timing1(6),timing0(6),timing1(7),timing0(7)
    else
      ! no analysis:
      analys(1:numgrd(1),1:numgrd(2),1:numgrd(3),i) = 0.0
    endif
  enddo

  ! add increment to the background:

  call date_and_time(date0,time0,zone0,timing0)
  call stmasinc
  call date_and_time(date1,time1,zone1,timing1)
  print*,'time used for stmasinc: ',(timing1(6)-timing0(6))+(timing1(7)-timing0(7))/60.0
  print*,'times used for stmasinc: ',timing1(6),timing0(6),timing1(7),timing0(7)

  ! write analyses to lsx
  if (savdat .eq. 1) then
    write(*,*) numgrd,numvar
    write(11,*) analys
  endif
  call date_and_time(date0,time0,zone0,timing0)
  call prpstlsx
  call date_and_time(date1,time1,zone1,timing1)
  print*,'time used for prpstlsx: ',(timing1(6)-timing0(6))+(timing1(7)-timing0(7))/60.0
  print*,'times used for prpstlsx: ',timing1(6),timing0(6),timing1(7),timing0(7)

  ! modified by min-ken hsieh
  ! verify stmas analysis
  call date_and_time(date0,time0,zone0,timing0)
  call stmasver
  call date_and_time(date1,time1,zone1,timing1)
  print*,'time used for stmasver: ',(timing1(6)-timing0(6))+(timing1(7)-timing0(7))/60.0
  print*,'times used for stmasver: ',timing1(6),timing0(6),timing1(7),timing0(7)

  ! release dynamic memory:
  call intrrels
  call stmarels

  ! end of analysis:
  write(*,*) 'stmas analysis succeeds'

end program stmas_mg
