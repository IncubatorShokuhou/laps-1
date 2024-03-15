module grib
  
! module containing routines to allow output of grib data
! requires the ncep w3fi library routines and io routines developed
! for ruc. 
!
! reference:ncep office note 388, grib (edition 1)
  use map_utils
  implicit none


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_igds(proj,igds)
    
    implicit none
    type(proj_info),intent(in)           :: proj
    integer, intent(out)                 :: igds(18)
    real, parameter                      :: pi = 3.1415927
    real, parameter                      :: deg_per_rad = 180./pi
    real, parameter                      :: rad_per_deg = pi / 180.
    real                                 :: scale , latne, lonne
    igds(1) = 0                          ! number of vertical (not used)
    igds(2) = 255                        ! pv,pl or 255 (not used)
    
    ! set data representation type (gds octet 6, table 6 from reference)
    select case (proj%code)
      case (proj_ps)
         igds(3) = 5
         igds(4) = proj%nx        ! e-w dimension
         igds(5) = proj%ny        ! n-s dimension
         igds(6) = nint(proj%lat1*1000.)   ! sw lat in millidegrees
         igds(7) = nint(proj%lon1*1000.)   ! sw lon in millidegrees
         igds(8) = 8                       ! 
         igds(9) =nint(proj%stdlon*1000.)
 
         ! adjust grid lengths to be true at 60 deg latitude
         scale = (1.+sin(60.*rad_per_deg)) / &
                 (1.+proj%hemi*sin(proj%truelat1*rad_per_deg))
         igds(10) = nint(proj%dx) * scale
         igds(11) = nint(proj%dx) * scale
         if (proj%hemi .eq. -1 ) then   
           igds(12) = 1
         else
           igds(12) = 0
         endif
         igds(13) = 64
         igds(14) = 0
         igds(15) = 0
         igds(16) = 0
         igds(17) = 0
         igds(18) = 0
      case (proj_lc)
         igds(3) = 3
         igds(4) = proj%nx        ! e-w dimension
         igds(5) = proj%ny        ! n-s dimension
         igds(6) = nint(proj%lat1*1000.)   ! sw lat in millidegrees
         igds(7) = nint(proj%lon1*1000.)   ! sw lon in millidegrees
         igds(8) = 8                       ! 
         igds(9) = nint(proj%stdlon*1000.)
         igds(10) = nint(proj%dx) 
         igds(11) = nint(proj%dx) 
         if (proj%hemi .eq. -1 ) then
           igds(12) = 1
         else
           igds(12) = 0
         endif
         igds(13) = 64
         igds(14) = 0
         igds(15) = nint(proj%truelat1*1000.)
         igds(16) = nint(proj%truelat2*1000.)
         igds(17) = -90000  ! latitude of southern pole?
         igds(18) = nint(proj%stdlon*1000.)  ! longitude of southern pole?
      case (proj_merc)
         igds(3) = 1
         igds(4) = proj%nx        ! e-w dimension
         igds(5) = proj%ny        ! n-s dimension
         igds(6) = nint(proj%lat1*1000.)   ! sw lat in millidegrees
         igds(7) = nint(proj%lon1*1000.)   ! sw lon in millidegrees
         igds(8) = 0                       !  
         ! compute lat/lon at nx/ny
         call ij_to_latlon(proj,float(proj%nx),float(proj%ny),latne,lonne)
         igds(9) = nint(latne*1000.)
         igds(10) = nint(lonne*1000.)
         igds(11) = nint(proj%dx)
         igds(12) = nint(proj%dx)
         igds(13) = 64
         igds(14:18) = 0
      case default 
         print *, 'projection code not supported: ', proj%code
         stop
    end select
    return
  end subroutine make_igds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_id(table_version,center_id,subcenter_id,process_id, &
                     param,leveltype,level1,level2,yyyyr,mmr,ddr,hhr,minr, &
                     timeunit, timerange, timeperiod1, timeperiod2,&
                     scalep10,id)

    ! routine to make the id integer array passed into w3fi72.  the id
    ! array is used to generate the grib message pds section.

    implicit none
    integer, intent(in)         :: table_version 
    integer, intent(in)         :: center_id
    integer, intent(in)         :: subcenter_id
    integer, intent(in)         :: process_id
    integer, intent(in)         :: param       ! grib code for variable
    integer, intent(in)         :: leveltype   ! grib code for level type
    integer, intent(in)         :: level1      ! value of level1
    integer, intent(in)         :: level2      ! value of level2
    integer, intent(in)         :: yyyyr       ! ref (initial) time year
    integer, intent(in)         :: mmr         ! ref (initial) time month
    integer, intent(in)         :: ddr         ! ref (initial) day of month
    integer, intent(in)         :: hhr         ! ref (initial) hour (utc)
    integer, intent(in)         :: minr        ! ref (initial) minute
    integer, intent(in)         :: timeperiod1 ! values used for forecast 
    integer, intent(in)         :: timeperiod2 !   time wrt reference time
    integer, intent(in)         :: timeunit    ! grib time unit (table 4)
    integer, intent(in)         :: timerange   ! grib time range (table 5)
    integer, intent(in)         :: scalep10    ! scaling power of 10 
    integer, intent(out)        :: id(27)      
    ! some stuff from the grib setup that are for now supplied by the
    ! include file at compile time, but should be moved to a namelist 
    ! eventually
    id(1) =  28
    id(2) = table_version
    id(3) = center_id
    id(4) = process_id
    id(5) = 255    ! we use the gds section to define grid
    id(6) = 1      ! gds included
    id(7) = 0      ! no bms or bitmask

    ! stuff from arguments
    id(8) = param       ! see table 2 of reference document
    id(9) = leveltype   ! see table 3 of reference document

    ! level stuff, dependent upon value of leveltype
    if ( ((leveltype.ge.1).and.(leveltype.le.100)) .or. &
          (leveltype.eq.102).or.(leveltype.eq.103).or.&
          (leveltype.eq.105).or.(leveltype.eq.107).or.&
          (leveltype.eq.109).or.(leveltype.eq.109).or.&
          (leveltype.eq.111).or.(leveltype.eq.113).or.&
          (leveltype.eq.115).or.(leveltype.eq.117).or.&
          (leveltype.eq.119).or.(leveltype.eq.125).or.&
          (leveltype.eq.160).or.(leveltype.eq.200).or.&
          (leveltype.eq.201) ) then
      id(10) = 0
      id(11) = level1
    else
      id(10) = level1
      id(11) = level2
    endif

    ! set reference time, which is the valid time for analyses and the
    ! initial time for forecasts
    id(12) = mod(yyyyr,100)  ! year of century
    id(13) = mmr             ! month
    id(14) = ddr             ! day
    id(15) = hhr             ! hour (utc)
    id(16) = minr            ! minute 
    id(17) = timeunit        ! unit indicator from table 4
    id(18) = timeperiod1
    id(19) = timeperiod2
    id(20) = timerange 
 
    ! flags to describe averaging (not yet used for our application)
    id(21) = 0               ! number included in average
    id(22) = 0               ! number missing from average

    ! miscellaneous stuff
    id(23) =(yyyyr/100) + 1  ! integer math to get century
    id(24) = subcenter_id   
    id(25) = scalep10        ! scaling power of 10 for precision preservation
    id(26) = 0           
    id(27) = 0               ! not used  
  end subroutine make_id 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_grib(itype,fld,id,igds,funit,startb,itot,istatus)

    ! subroutine to generate the grib message and write it to a file which
    ! must already be opened using open_grib routine.

    implicit none
    integer,intent(in)            :: itype ! (0 for real data, 1 for int data)
    real,intent(in)               :: fld(:)  ! input real data array
      ! note: if integer output desired, then pass the data in as real
      ! anyway, and this routine will convert it back to integer.

    integer,intent(in)            :: id(27)   ! input id array for pds
    integer,intent(in)            :: igds(18) ! input igds array for gds
    integer,intent(in)            :: funit   ! unit # for output
    integer,intent(in)            :: startb  ! starting record number
    integer,intent(out)           :: itot    ! number of bytes written
    integer,intent(out)           :: istatus ! status (0=ok)

    integer,parameter             :: maxbuf=512000
    integer,allocatable           :: ifld(:)
    integer,allocatable           :: ibmap(:)
    character(len=1)              :: kbuf(maxbuf)
    integer                       :: ibitl
    integer                       :: ipflag
    integer                       :: igflag
    integer                       :: igrid
    integer                       :: icomp
    integer                       :: ibflag
    integer                       :: iblen
    integer                       :: ibdsfl(9)
    character(len=1)              :: pds(28)
    integer                       :: npts,jerr
    integer                       :: nxny
    integer, external             :: c_write_g
    integer                       :: iwrite
    integer                       :: i,j,ibuf 
    istatus = 0

    nxny = igds(4)*igds(5)    ! nx*ny
    ! allocate the integer data array.  if itype is 1, then 
    ! also populate it.  otherwise, it is just used as a dummy
    ! argument
    allocate(ifld(nxny))
    if (itype .eq. 1) then 
      do j=0,igds(5)-1
        do i=1,igds(4)
          ifld( j*igds(4)+i ) = nint(fld( j*igds(4)+i ))
        enddo
      enddo
    endif
    ! allocate ibmap (dummy for now)
    allocate(ibmap(nxny))
    
    kbuf(:) = char(0)
    ! for now, there are a lot of hard coded options that we can
    ! make more flexible later on (e.g., use of bit maps, etc.)

    ibitl = 0   ! let compute pick optimum packgin length
    ipflag = 0  ! create pds from user supplied id array
    igflag = 1  ! create gds from supplied igds array
    igrid = 255 ! using gds from supplied igds array
    icomp = 1   ! grid-oriented winds (always the case in laps)
    ibflag = 0
    iblen = nxny
    ibdsfl(1) = 0  ! grid point data
    ibdsfl(2) = 0  ! simple packing
    ibdsfl(3) = itype
    ibdsfl(4) = 0  ! no additional flags at octet 14
    ibdsfl(5) = 0  ! always set to 0 (reserved)
    ibdsfl(6) = 0  ! single datum at each gridpoint
    ibdsfl(7) = 0  ! no secondary bit maps present
    ibdsfl(8) = 0  ! second order values have constant widths
    ibdsfl(9) =0
 
    ! make the grib message, which will be put into kbuf with information
    ! on its exact length provided in itot
    call w3fi72(itype,fld,ifld,ibitl,ipflag,id,pds,igflag,igrid,igds, &
                icomp,ibflag,ibmap,iblen,ibdsfl,npts,kbuf,itot,jerr)
    deallocate(ibmap)
    deallocate(ifld)
    ! check error status
    if (jerr .ne. 0) then
      print *, 'error creating grib message...jerr = ', jerr
      print *, 'npts/itot = ',npts,itot
      istatus = 1
      return
    endif
    if (itot .gt. maxbuf) then
      print *, 'message size larger than buffer allocation!'
      istatus = 1
      return
    endif
    ! ready to write the message
    !print *,'writing ',itot,'bytes starting at ', startb

    ! fortran method
    !do i = 1,itot
    !  write(funit,rec=i+startb-1) kbuf(i)
    !enddo
    ! c method
    iwrite = c_write_g(0,itot,kbuf,funit)
    if (iwrite .ne. 0) then
      print *, 'error writing grib message -- ', iwrite
      istatus = 1
      return
    endif
    !print *,'write completed.'
    return
  end subroutine write_grib
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine open_grib_f(fname,funit)
   
    ! opens a grib file for writing and returns the funit
    implicit none
    character(len=*),intent(in)      :: fname
    integer,intent(out)                :: funit
    logical                            :: opened 
    unitloop: do funit=7,1023
      inquire(unit=funit,opened=opened) 
      if (.not.opened) exit unitloop
    enddo unitloop
    open(file=fname,unit=funit,access='direct',form='unformatted',&
         recl=1) 
    return
  end subroutine open_grib_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine open_grib_c(fname,funit)
  
    ! opens a grib file for writing and returns the funit
    implicit none
    character(len=*),intent(in)      :: fname
    integer,intent(out)                :: funit
    integer, external                  :: c_open_g
    integer                            :: length
    logical                            :: opened
    call s_len(fname,length)
    funit = -1
    funit = c_open_g(fname(1:length)//char(0),'w'//char(0))
    return
  end subroutine open_grib_c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine close_grib_f(funit)

    implicit none
    integer, intent(in)       :: funit

    close(funit) 
 
    return
  end subroutine close_grib_f    
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine close_grib_c(funit)

    implicit none
    integer, intent(in)       :: funit
    integer, external         :: c_close_g
    integer                   :: iretc

    iretc = -1
    iretc = c_close_g(funit)
    print *, 'grib file closed with iretc = ', iretc
    return
  end subroutine close_grib_c
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module grib
