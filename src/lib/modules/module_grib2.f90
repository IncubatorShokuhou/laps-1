module grib2
  
! module containing routines to allow output of grib2 data
! requires the nws/mdl pk_grib2 library.

  use map_utils
!  use constants
  implicit none
  
  private 
  integer, parameter                :: kfildo = 6 
  ! variables to contain grib sections
  integer, parameter                :: nidat = 0
  integer, parameter                :: nrdat = 0
  integer                           :: idat ! normally would have nidat elements
  real                              :: rdat ! normally would have nrdat elements
  integer, parameter                :: ns0 = 16
  integer, parameter                :: ns1 = 21
  integer, parameter                :: ns3 = 96
  integer, parameter                :: ns4 = 60
  integer, parameter                :: ns5 = 49
  integer, parameter                :: ns6 = 6
  integer, parameter                :: ns7 = 8
  integer, parameter                :: ndjer = 30
  integer, parameter                :: l3264b = 32
  integer                           :: idum
  real                              :: rdum
  integer, parameter                :: master_table = 1
  integer, parameter                :: local_table = 1
  integer, parameter                :: grib_edition = 2                                                                                                                                                      
  integer                           :: nd5
  integer                           :: is0 (ns0)
  integer                           :: is1 (ns1)
  integer                           :: is3 (ns3)
  integer                           :: is4 (ns4)
  integer                           :: is5 (ns5)
  integer                           :: is6 (ns6)
  integer                           :: is7 (ns7)
  integer, allocatable              :: ib(:,:)
  integer, allocatable              :: ipack(:)
  integer                           :: jer(ndjer,2)
  integer                           :: kjer
  integer,parameter                 :: minpk=14
  logical                           :: big_endian
  integer                           :: ig2status
  logical                           :: made_sec1
  logical                           :: made_sec3
 
  public open_binary
  public init_grib2_file
  public close_grib2_file
  public close_binary
  public write_grib2_template0
  public write_grib2_template8
contains

  subroutine init_grib2_file(fname,grid,center_id,subcenter_id, &
                             reftime_sig,year,month,day,hour, &
                             minute,second,prod_status,data_type,&
                             funit,istatus)

    ! opens a new grib output file, initializes some of the grib
    ! headers
    implicit none
    character(len=*),intent(in)  :: fname
    type(proj_info)              :: grid
    integer,intent(in)           :: center_id,subcenter_id
    integer,intent(in)           :: reftime_sig
    integer,intent(in)           :: year,month,day,hour,minute,second
    integer,intent(in)           :: prod_status,data_type
    integer,intent(out)          :: funit
    integer,intent(out)          :: istatus

    istatus = 0
    ig2status = 0
    call open_binary(fname,funit)
    call make_g2_sec1(center_id, subcenter_id, reftime_sig,&
                     year,month,day,hour,minute,second,prod_status,data_type)
    call make_g2_sec3(grid)

    if (ig2status.ne.0) then
      print *, "error initializing grib2 file"
      istatus = 1
    endif
  
  end subroutine init_grib2_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine close_grib2_file(funit)
  
    implicit none
    integer, intent(in) :: funit
    
    call close_binary(funit)
    if (allocated(ib)) deallocate(ib)
    if (allocated(ipack)) deallocate(ipack)
    return
  end subroutine close_grib2_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_g2_sec0(discipline)

    implicit none
    integer, intent(in)  :: discipline  ! discipline table number
    
    is0(:) = 0
    is0(7) = discipline
    is0(8) = grib_edition
    return
  end subroutine make_g2_sec0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_g2_sec1(center_id, subcenter_id, reftime_sig,&
                          year,month,day,hour,minute,second,prod_status,data_type)

    implicit none
    integer, intent(in)    :: center_id            ! common code table c-1
    integer, intent(in)    :: subcenter_id         ! 
    integer, intent(in)    :: reftime_sig          ! signficance of reference time, code table 1.2, where 0=analysis, 1=fcst
    integer, intent(in)    :: year                 ! reference time year (4-digits)
    integer, intent(in)    :: month                ! reference time month
    integer, intent(in)    :: day                  ! reference time day 
    integer, intent(in)    :: hour                 ! reference time hour
    integer, intent(in)    :: minute               ! reference time minute
    integer, intent(in)    :: second               ! reference time second
    integer, intent(in)    :: prod_status          ! production status
    integer, intent(in)    :: data_type            ! code table 1.4, 0=analysis,1=forecat,2=analysis and forecast
                                                   ! 7=radar observations 

    is1(:) = 255
    is1(5) = 1  ! section id
    is1(6) = center_id
    is1(8) = subcenter_id
    is1(10) = master_table
    is1(11) = local_table
    is1(12) = reftime_sig
    is1(13) = year
    is1(15) = month
    is1(16) = day
    is1(17) = hour
    is1(18) = minute
    is1(19) = second
    is1(20) = prod_status
    is1(21) = data_type 
    made_sec1 = .true.
    return
  end subroutine make_g2_sec1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_g2_sec3(grid)

    implicit none
    type(proj_info),intent(in)           :: grid
    real                                 :: lat1,lat2,lon1,lon2,lov
    real, parameter                      :: llscale = 1.e6
    real, parameter                      :: dxscale = 1000.
    real, parameter                      :: eradscale = 1
    ig2status = 0 
    ! initialize is3 

    is3(:) = 0 

    ! set section number
    is3(5) = 3
    is3(6) = 0  ! method of specifying grid (0 = uses templates)
    is3(7) = grid%nx * grid%ny
    is3(11) = 0 ! regular grid
    is3(12) = 0 ! regular grid

    ! starting lat/lon
    lat1 = grid%lat1
    lon1 = grid%lon1
    call ij_to_latlon(grid,float(grid%nx),float(grid%ny),lat2,lon2)
    lov = grid%stdlon 
    if (lon1 .lt. 0) lon1 = lon1 + 360.
    if (lon2 .lt. 0) lon2 = lon2 + 360.
    if (lov .lt. 0) lov = lov + 360.

    if (grid%code .eq. proj_merc) then
      is3(13) = 10
      is3(15) = 1  ! spherical earth, radius provided by user
      is3(16) = 0  ! scale factor for spherical earth radius
      is3(17) = nint(earth_radius_m*eradscale)
      is3(31) = grid%nx
      is3(35) = grid%ny
      is3(39) = nint(lat1 * llscale)
      is3(43) = nint(lon1 * llscale)
      is3(47) = 56 ! resolution/component flag
      is3(48) = nint(grid%truelat1 * llscale)
      is3(52) = nint(lat2 * llscale)
      is3(56) = nint(lon2 * llscale)
      is3(60) = 64
      is3(61) = 0
      is3(65) = nint(grid%dx * dxscale)
      is3(69) = is3(65)
    elseif (grid%code .eq. proj_lc) then
      is3(13) = 30
      is3(15) = 1  ! spherical earth, radius provided by user
      is3(16) = 0  ! scale factor for spherical earth radius
      is3(17) = nint(earth_radius_m*eradscale)
      is3(31) = grid%nx
      is3(35) = grid%ny
      is3(39) = nint(lat1 * llscale)
      is3(43) = nint(lon1 * llscale)
      is3(47) = 8 ! resolution/component flag
      is3(48) = nint(grid%truelat1 * llscale)
      is3(52) = nint(lov * llscale)
      is3(56) = nint(grid%dx * dxscale)
      is3(60) = is3(56)
      if (grid%truelat1 .gt. 0) then
        is3(64) = 0
      else
        is3(64) = 128
      endif
      is3(65) = 64
      is3(66) = nint(grid%truelat1 * llscale)
      is3(70) = nint(grid%truelat2 * llscale)
      ! i still don't understand the definition of "lat/lon"
      ! of southern pole of project...hope this is correct
      is3(74) = nint(-90 * llscale)
      is3(78) = nint(lov * llscale)
                                                                                                                                                        
    elseif (grid%code .eq. proj_ps) then
      is3(13) = 20
      is3(15) = 1  ! spherical earth, radius provided by user
      is3(16) = 0  ! scale factor for spherical earth radius
      is3(17) = nint(earth_radius_m*eradscale)
      is3(31) = grid%nx
      is3(35) = grid%ny
      is3(39) = nint(lat1 * llscale)
      is3(43) = nint(lon1 * llscale)
      is3(47) = 8 ! resolution/component flag
      is3(48) = nint(grid%truelat1 * llscale)
      is3(52) = nint(lov * llscale)
      is3(56) = nint(grid%dx * dxscale)
      is3(60) = is3(56)
      if (grid%truelat1 .gt. 0) then
        is3(64) = 0
      else
        is3(64) = 128
      endif
      is3(65) = 64
    else
      print *, "output grib2:  projection not supported!"
      ig2status = 1
    endif
   
    ! allocate space for the packed data
    nd5 = (250 + (grid%nx * grid%ny) + (grid%nx * grid%ny) / 8 + &
         nidat + nrdat )
    if (allocated(ipack)) deallocate(ipack)
    if (allocated(ib)) deallocate(ib)
    allocate(ipack(nd5))
    allocate(ib(grid%nx,grid%ny))
    made_sec3 = .true. 
    return
  end subroutine make_g2_sec3 
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine open_binary(fname,funit)
  
    ! opens a grib file for writing and returns the funit
    implicit none
    character(len=*),intent(in)      :: fname
    integer,intent(out)                :: funit
    integer, external                  :: c_open_g
    integer                            :: length
    logical, external                  :: pk_endian
    length = len_trim(fname)
    funit = -1
    funit = c_open_g(fname(1:length)//char(0),'w'//char(0))
    big_endian = pk_endian()
    print *, "opened grib2 output file: ", trim(fname)
    print *, "  unit number: ",funit
    print *, "  big_endian: ",big_endian
    return
  end subroutine open_binary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine close_binary(funit)

    implicit none
    integer, intent(in)       :: funit
    integer, external         :: c_close_g
    integer                   :: iretc

    iretc = -1
    iretc = c_close_g(funit)
    print *, 'grib file closed with iretc = ', iretc
    return
  end subroutine close_binary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_grib2_template0(funit,discipline,param_category,param_number, &
                                   process_type,bg_process_id, &
                                   af_process_id,cutoff_hr, &
                                   cutoff_min,time_unit_indicator, &
                                   ftime,lev1_type,lev1_scale, lev1_value,&
                                   lev2_type,lev2_scale,lev2_value, &
                                   pack_method,data_scale,miss_mgmt,&
                                   nx,ny,newrec,inomiss,xmissp,xmisss,data_value)
 
    ! this subroutine outputs grib2 data to an already opened grib2 file, which
    ! was opened by the open_grib2 routine and specified via funit.  
    ! it assumes you have already correctly populated is1 and is3 via
    ! prior calls to make_g2_sec1 and make_g2_sec3.  
    ! this particular routine writes data using product definition template 4.0
    ! and data respresentation template 5.2 (gridpoint complex)
    implicit none

    integer,intent(in)   :: funit     ! file specifier of open grib2 file
    integer,intent(in)   :: discipline ! see code table 0.0
    integer,intent(in)   :: param_category  ! table 4.1
    integer,intent(in)   :: param_number    ! table 4.2
    integer,intent(in)   :: process_type    ! table 4.3
    integer,intent(in)   :: bg_process_id   ! center defined
    integer,intent(in)   :: af_process_id   ! center defined
    integer,intent(in)   :: cutoff_hr
    integer,intent(in)   :: cutoff_min
    integer,intent(in)   :: time_unit_indicator ! table 4.4
    integer,intent(in)   :: ftime  ! forecast time
    integer,intent(in)   :: lev1_type,lev1_scale,lev1_value
    integer,intent(in)   :: lev2_type,lev2_scale,lev2_value
    integer,intent(in)   :: pack_method,data_scale,miss_mgmt
    integer,intent(in)   :: nx,ny,newrec,inomiss
    real,intent(in)      :: xmissp,xmisss
    real,intent(in)      :: data_value(nx,ny)
    integer              :: dummy1
    integer  :: ibitmap 
    integer  :: imissp, imisss

    imissp = nint(xmissp)
    imisss = nint(xmisss)
    ig2status = 0   

    if (.not. made_sec3) then
      ig2status = 1
      print *,"called write_grib2_template0 before section 3 was made!"
      return
    endif
    if (.not. made_sec1) then
      ig2status = 1
      print *,"called write_grib2_template0 before section 1 was made!"
      return
    endif
  
    ! sections 1 and 3 should already be set up, we ignore section 2,
    ! so finish section 0, 4, 5, and 6 
    ! set up secion 0 (is0)
    call make_g2_sec0(discipline)

    ! set up section 5 (is5)
    call make_g2_sec5(pack_method,data_scale,miss_mgmt)
    if (ig2status.ne.0) then
      print *, "problem creating section 5"
      return
    endif

    ! set up section 6 (bitmap)...currently not using
    is6(:) = 0
    is6(5) = 6
    is6(6) = 255
    ibitmap = 0
  
    ! set up section 4 ..pdss
    is4(:) = 0
    is4(5) = 4
    is4(6) = 0
    is4(8) = 0
    is4(10) = param_category
    is4(11) = param_number
    is4(12) = process_type  
    is4(13) = bg_process_id
    is4(14) = af_process_id
    is4(15) = cutoff_hr
    is4(17) = cutoff_min
    is4(18) = time_unit_indicator
    is4(19) = ftime
    is4(23) = lev1_type
    is4(24) = lev1_scale
    is4(25) = lev1_value
    is4(29) = lev2_type
    is4(30) = lev2_scale
    is4(31) = lev2_value

    ! set up section 7
    is7(:) = 0
    is7(5) = 7
    call pk_grib2(kfildo,data_value,dummy1,nx,ny,idat,nidat,rdat,nrdat,&
                  is0,ns0,is1,ns1,is3,ns3,is4,ns4,is5,ns5,is6,ns6, &
                  is7,ns7,ib,ibitmap,ipack,nd5,imissp,xmissp,imisss,xmisss,&
                  newrec,minpk,inomiss,l3264b,jer,ndjer,kjer)

    call error_check
    if (.not. big_endian) then
      call c_swap4(nd5*4,ipack)
    endif
    call c_write_g(0,is0(9),ipack,funit)

  end subroutine write_grib2_template0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_g2_sec5(pack_method,data_scale,miss_mgmt)

    ! sets up section 5
    implicit none

    integer, intent(in) :: pack_method,data_scale,miss_mgmt

    is5(:) = 0
    is5(5) = 5

    if ((pack_method .ne. 0) .and. (pack_method .ne. 2))then

      print *, "currently, module_grib2 only supports pack_method 0 or 2"
      ig2status = 1
      return
    endif
    is5(10) = pack_method
    is5(18) = data_scale
    if (pack_method .eq. 2) then
      is5(22) = 1 ! general splitting
      is5(23) = miss_mgmt
    endif
    return
  end subroutine make_g2_sec5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_grib2_template8(funit,discipline,param_category,param_number, &
                                   process_type,bg_process_id, &
                                   af_process_id,cutoff_hr, &
                                   cutoff_min,time_unit_indicator, &
                                   ftime,lev1_type,lev1_scale, lev1_value,&
                                   lev2_type,lev2_scale,lev2_value, &
                                   eyear,emon,eday,ehour,emin,esec,&
                                   ntimes,ntimes_miss,stattype,periodtype,&
                                   etime_unit,etime_value,&
                                   pack_method,data_scale,miss_mgmt,&
                                   nx,ny,newrec,inomiss,xmissp,xmisss,data_value)
                                                                                                                                                                                                                                  
    ! this subroutine outputs grib2 data to an already opened grib2 file, which
    ! was opened by the open_grib2 routine and specified via funit.
    ! it assumes you have already correctly populated is1 and is3 via
    ! prior calls to make_g2_sec1 and make_g2_sec3.
    ! this particular routine writes data using product definition template 4.8

    implicit none
                                                                                                                                                                                                                                  
    integer,intent(in)   :: funit     ! file specifier of open grib2 file
    integer,intent(in)   :: discipline ! see code table 0.0
    integer,intent(in)   :: param_category  ! table 4.1
    integer,intent(in)   :: param_number    ! table 4.2
    integer,intent(in)   :: process_type    ! table 4.3
    integer,intent(in)   :: bg_process_id   ! center defined
    integer,intent(in)   :: af_process_id   ! center defined
    integer,intent(in)   :: cutoff_hr
    integer,intent(in)   :: cutoff_min
    integer,intent(in)   :: time_unit_indicator ! table 4.4
    integer,intent(in)   :: ftime  ! forecast time
    integer,intent(in)   :: lev1_type,lev1_scale,lev1_value
    integer,intent(in)   :: lev2_type,lev2_scale,lev2_value
    integer,intent(in)   :: eyear,emon,eday,ehour,emin,esec
    integer,intent(in)   :: ntimes,ntimes_miss,stattype,periodtype
    integer,intent(in)   :: etime_unit,etime_value
    integer,intent(in)   :: pack_method,data_scale,miss_mgmt
    integer,intent(in)   :: nx,ny,newrec,inomiss
    real, intent(in)     :: xmissp, xmisss
    real,intent(in)      :: data_value(nx,ny)
    integer              :: dummy1
    integer  :: ibitmap
    integer  :: imissp, imisss  

    imissp = nint(xmissp)
    imisss = nint(xmisss)                                                                                                                                                                                       
    ig2status = 0
                                                                                                                                                                                                                                  
    if (.not. made_sec3) then
      ig2status = 1
      print *,"called write_grib2_template8 before section 3 was made!"
      return
    endif
    if (.not. made_sec1) then
      ig2status = 1
      print *,"called write_grib2_template8 before section 1 was made!"
      return
    endif
                                                                                                                                                                                                                                  
    ! sections 1 and 3 should already be set up, we ignore section 2,
    ! so finish section 0, 4, 5, and 6
    ! set up secion 0 (is0)
    call make_g2_sec0(discipline)
                                                                                                                                                                                                                                  
    ! set up section 5 (is5)
    call make_g2_sec5(pack_method,data_scale,miss_mgmt)
    if (ig2status.ne.0) then
      print *, "problem creating section 5"
      return
    endif
                                                                                                                                                                                                                                  
    ! set up section 6 (bitmap)...currently not using
    is6(:) = 0
    is6(5) = 6
    is6(6) = 255
    ibitmap = 0

    ! set up section 4 ..pdss
    is4(:) = 0
    is4(5) = 4
    is4(6) = 0
    is4(8) = 8
    is4(10) = param_category
    is4(11) = param_number
    is4(12) = process_type
    is4(13) = bg_process_id
    is4(14) = af_process_id
    is4(15) = cutoff_hr
    is4(17) = cutoff_min
    is4(18) = time_unit_indicator
    is4(19) = ftime
    is4(23) = lev1_type
    is4(24) = lev1_scale
    is4(25) = lev1_value
    is4(29) = lev2_type
    is4(30) = lev2_scale
    is4(31) = lev2_value
    ! time of end of overall time interval
    is4(35) = eyear
    is4(37) = emon
    is4(38) = eday
    is4(39) = ehour
    is4(40) = emin
    is4(41) = esec
    ! no. time range specs desc time intvls used to cal statistically-processed field
    is4(42) = ntimes
    is4(43) = ntimes_miss
    is4(47) = stattype
    is4(48) = periodtype
    ! unit of time for time range over which statistical processing is done (table 4.4)
    is4(49) = etime_unit
    is4(50) = etime_value
    is4(54) = 255
    is4(55:58) = 0

    ! set up section 7
    is7(:) = 0
    is7(5) = 7
    call pk_grib2(kfildo,data_value,dummy1,nx,ny,idat,nidat,rdat,nrdat,&
                  is0,ns0,is1,ns1,is3,ns3,is4,ns4,is5,ns5,is6,ns6, &
                  is7,ns7,ib,ibitmap,ipack,nd5,imissp,xmissp,imisss,xmisss,&
                  newrec,minpk,inomiss,l3264b,jer,ndjer,kjer)
                                                                                                                                                                                                                                  
    if (.not. big_endian) then
      call c_swap4(nd5*4,ipack)
    endif
    call c_write_g(0,is0(9),ipack,funit)
                                                                                                                                                                                                                                  
  end subroutine write_grib2_template8
 
!#####################################################
   subroutine error_check
 
     integer  :: j

     
     do j = 1, ndjer

       if (jer(j,2) .ne. 0) then 
         print '("warning: grib2 write error code: ",i4,3x,i2)', jer(j,1),jer(j,2)
       endif
     enddo
     return
        


    
   end subroutine error_check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module grib2
