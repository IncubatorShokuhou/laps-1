!dis
!dis    open source license/disclaimer, forecast systems laboratory
!dis    noaa/oar/fsl, 325 broadway boulder, co 80305
!dis
!dis    this software is distributed under the open source definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    in particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - all modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - if significant modifications or enhancements are made to this
!dis    software, the fsl software policy manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    this software and its documentation are in the public domain
!dis    and are furnished "as is."  the authors, the united states
!dis    government, its instrumentalities, officers, employees, and
!dis    agents make no warranty, express or implied, as to the usefulness
!dis    of the software and documentation for any purpose.  they assume
!dis    no responsibility (1) for the use of the software and
!dis    documentation; or (2) to provide technical support to users.
!dis
!dis
!  module to handle i/o from a wrf v2 netcdf file
!

module wrf_netcdf

! f90 module to deal with reading wrf output files (wrfv1 netcdf)

  use grid_utils
  include 'netcdf.inc'
  private
    
  public open_wrfnc
  public close_wrfnc
  public get_wrfsi_map
  public get_wrf2_map
  public get_wrf2_timeinfo
  public get_wrf_misc
  public get_wrf2_misc
  public get_wrf_scalar
  public get_wrf_1d
  public get_wrfnc_2d
  public get_wrfnc_3d
  public make_wrf_file_name
  public make_wrf2_file_name
  public wrfio_wait

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine open_wrfnc(fname, lun, status)

    ! opens a wrf netcdf file, returning the cdf id of the file

    implicit none
    include 'netcdf.inc'  
    ! arguments

    character(len=255),intent(in)          :: fname
    integer, intent(out)                   :: lun
    integer, intent(out)                   :: status

    status = 0
    lun = ncopn(fname, ncnowrit,   status)

    if (status .ne. 0) then
       print *, 'error opening netcdf file: ',trim(fname)
    endif

    return
  end subroutine open_wrfnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine close_wrfnc(cdfid)

    implicit none
    integer, intent(in)  :: cdfid
    integer :: status 
    include 'netcdf.inc'

    status = nf_close(cdfid)
    if (status .ne. nf_noerr) then
      print *, 'problem closing the netcdf file.'
    endif
    return
 end subroutine close_wrfnc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_map(lun, projection, lat1, lon1, stdlon, &
                         truelat1, truelat2, dx_m, dy_m, &
                         nx, ny, status)
  
    ! reads from a wrfsi static file to get mapping information for a wrf domain
    ! note:  the parameters returned are the non-staggered grid parameters,
    ! as defined in the wrfsi.nl file.  the dimensions found in the wrfout
    ! netcdf files require have one less element in each direction.

    ! assumes that the static file has been opened with the open_wrfnc routine
    ! and the cdf id is the "lun". 

    implicit none
    include 'netcdf.inc' 
    ! args
    integer, intent(in)                    :: lun  
    character (len=32), intent(out)        :: projection
    real, intent(out)                      :: lat1
    real, intent(out)                      :: lon1 
    real, intent(out)                      :: stdlon
    real, intent(out)                      :: truelat1
    real, intent(out)                      :: truelat2
    real, intent(out)                      :: dx_m
    real, intent(out)                      :: dy_m
    integer, intent(out)                   :: nx
    integer, intent(out)                   :: ny
    integer, intent(out)                   :: status
  
    !locals
    integer :: vid, rcode
    character(len=132) :: dum
    character(len=132) :: grid_type
    integer, dimension(2) :: startc, countc
    integer, dimension(4) :: start, count
    status = 0

    ! get x-y dimensions

    vid = ncvid(lun, 'x', rcode)
    call ncdinq(lun,vid,dum,nx,rcode)
    if (rcode .ne. 0) then
       print *, 'error getting nx.'
       status = 1
    endif

    vid = ncvid(lun, 'y', rcode)
    call ncdinq(lun,vid,dum,ny,rcode)
    if (rcode .ne. 0) then
       print *, 'error getting nx.'
       status = 1
    endif

    !  get projection
    startc = (/1,1/)
    countc = (/132,1/)
    vid = ncvid(lun, 'grid_type', rcode)
    call ncvgtc(lun,vid,startc,countc,grid_type, 132, rcode)

    select case(grid_type)
      case ('mercator') 
        projection = 'mercator                        '
      case ('secant lambert conformal')
        projection = 'lambert conformal               '
      case ('tangential lambert conformal')
        projection = 'lambert conformal               '
      case ('polar stereographic')
        projection = 'polar stereographic             '
      case default
        print *, 'unrecognized grid type: ', grid_type
        status = 1
    end select

    ! get sw corner point
    vid = ncvid(lun, 'lat', rcode)
    call ncvgt1(lun, vid, 1, lat1, rcode)
    vid = ncvid(lun, 'lon', rcode)
    call ncvgt1(lun, vid, 1, lon1, rcode) 

    ! get truelat1/trulat2
    vid = ncvid(lun, 'latin1', rcode)
    call ncvgt1(lun, vid, 1, truelat1, rcode)
    vid = ncvid(lun, 'latin2', rcode)
    call ncvgt1(lun, vid, 1, truelat2, rcode)

    ! get standard longitude
    vid = ncvid(lun, 'lov', rcode)
    call ncvgt1(lun, vid, 1, stdlon, rcode)

    ! get dx/dy
    vid = ncvid(lun, 'dx', rcode)
    call ncvgt1(lun, vid, 1, dx_m, rcode)
    vid = ncvid(lun, 'dy', rcode)
    call ncvgt1(lun, vid, 1, dy_m, rcode)

    ! clean up
    if (lon1 .gt. 180.) lon1 = lon1 - 360.
    if (stdlon .gt. 180.) stdlon = stdlon - 360.
 
    return

  end subroutine get_wrfsi_map 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrf2_timeinfo(lun,reftime,dt,itimestep,tau_hr,tau_min,tau_sec, &
                               status)

    implicit none
    include 'netcdf.inc'

    integer, intent(in)                  :: lun
    character(len=19), intent(out)       :: reftime
    real, intent(out)                    :: dt
    real                                 :: xtime
    integer, intent(out)                 :: itimestep
    integer, intent(out)                 :: tau_hr, tau_min, tau_sec, status

    integer   :: rcode, vid

    reftime = '????-??-??_??:??:??' 
    tau_hr = 0
    tau_min = 0
    tau_sec = 0
    status = 0

    rcode = nf_get_att_text(lun,0,"simulation_start_date",reftime)
    rcode = nf_get_att_real(lun,0,"dt",dt)
    rcode = nf_inq_varid(lun,"itimestep",vid)
    rcode = nf_get_var_int(lun,vid,itimestep)

    !!!!! modified by wei-ting (20130307) :            !!!!!
    !!!!!     use xtime instead of itimestep for wrfv3 !!!!!
    if (rcode .lt. 0) then
      write(6,*)' warning: itimestep not found, looking for xtime in get_wrf2_timeinfo'
      rcode = nf_inq_varid(lun,"xtime",vid)
      rcode = nf_get_var_real(lun,vid,xtime)
      itimestep = int(xtime*60/dt)
    endif
    !!!!! end of modifying !!!!!
    
    tau_hr = int(float(itimestep)*dt)/3600
    tau_sec = mod(int(float(itimestep)*dt),3600)
    tau_min = tau_sec / 60
    tau_sec = mod(tau_sec,60) 

    return
end subroutine get_wrf2_timeinfo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrf2_map(lun, stag, projcode, lat1, lon1, stdlon, &
                          truelat1, truelat2, dx_m, dy_m, &
                          nx, ny, nz, status)

    implicit none
    include 'netcdf.inc'
    ! args
    integer, intent(in)                    :: lun
    character (len=1)                      :: stag
    integer, intent(out)                   :: projcode  
    real, intent(out)                      :: lat1
    real, intent(out)                      :: lon1
    real, intent(out)                      :: stdlon
    real, intent(out)                      :: truelat1
    real, intent(out)                      :: truelat2
    real, intent(out)                      :: dx_m
    real, intent(out)                      :: dy_m
    integer, intent(out)                   :: nx
    integer, intent(out)                   :: ny
    integer, intent(out)                   :: nz
    integer, intent(out)                   :: status

  
    ! local variables
    integer   :: rcode,dimid,vid
    integer   :: projcode_wrf
    integer   :: wrf_version
    real, allocatable, dimension(:,:) :: dum_2d

    ! initialization
    lat1 = -999.
    lon1 = -999.
    projcode = 0     
    projcode_wrf = 0
    stdlon = -999.
    truelat1 = -999.
    truelat2 = -999.
    dx_m  = -999.
    dy_m  = -999.
    nx = 0
    ny = 0
    nz =0
    status = 0

    rcode = nf_get_att_int(lun, 0, "map_proj", projcode_wrf)
    select case (projcode_wrf)
      case(1)
        projcode = 3
      case(2)
        projcode = 5
      case(3)
        projcode = 1
      case default
        projcode = 99
    end select 
    rcode = nf_get_att_real(lun, 0, "stand_lon",stdlon)
    rcode = nf_get_att_real(lun, 0, "truelat1",truelat1)
    rcode = nf_get_att_real(lun, 0, "truelat2",truelat2)
    rcode = nf_get_att_real(lun, 0, "dx", dx_m)
    rcode = nf_get_att_real(lun, 0, "dy", dy_m)
    
    ! dimenions
    rcode = nf_inq_dimid(lun, "west_east", dimid)
    rcode = nf_inq_dimlen(lun, dimid, nx)
    rcode = nf_inq_dimid(lun, "south_north", dimid)
    rcode = nf_inq_dimlen(lun, dimid, ny)

!   allocate 
    allocate(dum_2d (nx,ny))

    select case (stag)
      case('t')
        rcode = nf_inq_varid(lun,"lat_ll_t",vid)
        if(rcode .eq. -49)then ! assume version 3 wrf
          write(6,*)' warning: lat_ll_t not found, looking for xlat'
          rcode = nf_inq_varid(lun,"xlat",vid)
          rcode = nf_get_var_real(lun,vid,dum_2d)
          lat1 = dum_2d(1,1)
          wrf_version = 3
        else
          rcode = nf_get_var_real(lun,vid,lat1)
          wrf_version = 2
        endif

        rcode = nf_inq_varid(lun,"lon_ll_t",vid)
        if(rcode .eq. -49)then ! assume version 3 wrf
          write(6,*)' warning: lon_ll_t not found, looking for xlong'
          rcode = nf_inq_varid(lun,"xlong",vid)
          rcode = nf_get_var_real(lun,vid,dum_2d)
          lon1 = dum_2d(1,1)
        else
          rcode = nf_get_var_real(lun,vid,lon1)
        endif

      case('u')
        rcode = nf_inq_varid(lun,"lat_ll_u",vid)
        if(rcode .eq. -49)then ! assume version 3 wrf
          write(6,*)' warning: lat_ll_u not found, looking for xlat_u'
          rcode = nf_inq_varid(lun,"xlat_u",vid)
          rcode = nf_get_var_real(lun,vid,dum_2d)
          lat1 = dum_2d(1,1)
          wrf_version = 3
        else
          rcode = nf_get_var_real(lun,vid,lat1)
          wrf_version = 2
        endif

        rcode = nf_inq_varid(lun,"lon_ll_u",vid)
        if(rcode .eq. -49)then ! assume version 3 wrf
          write(6,*)' warning: lon_ll_u not found, looking for xlong_u'
          rcode = nf_inq_varid(lun,"xlong_u",vid)
          rcode = nf_get_var_real(lun,vid,dum_2d)
          lon1 = dum_2d(1,1)
        else
          rcode = nf_get_var_real(lun,vid,lon1)
        endif
        nx = nx + 1
 
      case('v')
        rcode = nf_inq_varid(lun,"lat_ll_v",vid)
        if(rcode .eq. -49)then ! assume version 3 wrf
          write(6,*)' warning: lat_ll_v not found, looking for xlat_v'
          rcode = nf_inq_varid(lun,"xlat_v",vid)
          rcode = nf_get_var_real(lun,vid,dum_2d)
          lat1 = dum_2d(1,1)
          wrf_version = 3
        else
          rcode = nf_get_var_real(lun,vid,lat1)
          wrf_version = 2
        endif

        rcode = nf_inq_varid(lun,"lon_ll_v",vid)
        if(rcode .eq. -49)then ! assume version 3 wrf
          write(6,*)' warning: lon_ll_v not found, looking for xlong_v'
          rcode = nf_inq_varid(lun,"xlong_v",vid)
          rcode = nf_get_var_real(lun,vid,dum_2d)
          lon1 = dum_2d(1,1)
        else
          rcode = nf_get_var_real(lun,vid,lon1)
        endif
        ny = ny + 1

    end select

    deallocate(dum_2d)

!   if(wrf_version .eq. 2)then
!       write(6,*)' version 2, looking for bottom_top for nz '
        rcode = nf_inq_dimid(lun, "bottom_top", dimid)
        rcode = nf_inq_dimlen(lun, dimid, nz)
!   else
!       write(6,*)' version 3, looking for bottom_top_stag for nz '
!       rcode = nf_inq_dimid(lun, "bottom_top_stag", dimid)
!       rcode = nf_inq_dimlen(lun, dimid, nz)
!   endif

    write(6,*)' wrf nx,ny,nz = ',nx,ny,nz

    return
  end subroutine get_wrf2_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrf_misc(cdfid, nzh, nzf, ptop, clwflag, iceflag, &
                          graupelflag)

    ! subroutine to get a few key parameters from the wrf model for the
    ! model post processor.

    implicit none
    include 'netcdf.inc'
    integer, intent(in)                :: cdfid   ! netcdf file handle
    integer, intent(out)               :: nzh     ! number of half-levels
    integer, intent(out)               :: nzf     ! number of full-levels
    real,    intent(out)               :: ptop    ! top pressure in pa
    logical, intent(out)               :: clwflag ! cloud liquid fields avail
    logical, intent(out)               :: iceflag ! ice species avail
    logical, intent(out)               :: graupelflag  ! grauple included

    ! local variables

    integer :: vid, rcode,mp_level
    character(len=132) :: dum

    ! get ptop, which is not really needed for anything in the em version
    vid = ncvid(cdfid, 'p_top', rcode)
    rcode  = nf_get_var_real(cdfid,vid,ptop)

    ! get the vertical dimensions
    rcode = nf_get_att_int(cdfid, nf_global, 'bottom-top_grid_dimension',nzh)
    nzf = nzh + 1

    ! determine which microphysics package is used and set flags
    ! accordingly

    rcode = nf_get_att_int(cdfid, nf_global, 'mp_physics',mp_level)

    select case (mp_level)

      case(0)  ! no microphysics
        clwflag = .false.
        iceflag = .false.
        graupelflag = .false.

      case(1) ! kessler warm rain
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.

      case(2) ! lin et al.
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.

      case(3)   ! ncep 3-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.
      case(4)   ! ncep 5-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.

      case(5)   ! eta ferrier 2-class
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.

      case default
        print *, 'warning:  cannot determine microphysics option!'
        print *, '          assuming all species present.'
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.

    end select

    return

  end subroutine get_wrf_misc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrf2_misc(cdfid, nzh, nzf, ptop, clwflag, iceflag, &
                          graupelflag)

    ! subroutine to get a few key parameters from the wrf model for the 
    ! model post processor.

    implicit none
    include 'netcdf.inc' 
    integer, intent(in)                :: cdfid   ! netcdf file handle
    integer, intent(out)               :: nzh     ! number of half-levels
    integer, intent(out)               :: nzf     ! number of full-levels
    real,    intent(out)               :: ptop    ! top pressure in pa
    logical, intent(out)               :: clwflag ! cloud liquid fields avail
    logical, intent(out)               :: iceflag ! ice species avail
    logical, intent(out)               :: graupelflag  ! grauple included

    ! local variables
     
    integer :: vid, rcode,mp_level
    character(len=132) :: dum

    ! get ptop, which is not really needed for anything in the em version
    vid = ncvid(cdfid, 'p_top', rcode)
    rcode  = nf_get_var_real(cdfid,vid,ptop)

    ! get the vertical dimensions
    rcode = nf_get_att_int(cdfid, nf_global, 'bottom-top_grid_dimension',nzf)
    nzh = nzf - 1

    ! determine which microphysics package is used and set flags 
    ! accordingly
 
    rcode = nf_get_att_int(cdfid, nf_global, 'mp_physics',mp_level) 
   
    select case (mp_level)

      case(0)  ! no microphysics
        clwflag = .false.
        iceflag = .false. 
        graupelflag = .false.

      case(1) ! kessler warm rain
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.

      case(2) ! lin et al.
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.
  
      case(3)   ! wsm 3-class simple
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.
  
      case(4)   ! wsm 5-class simple
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.

      case(5)   ! eta ferrier 2-class
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.
      
      case(6)  ! wsm 6-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.
  
      case(98) ! ncep 3-class
        clwflag = .true.
        iceflag = .false.
        graupelflag = .false.
 
      case(99) ! ncep 5-class
        clwflag = .true.
        iceflag = .true.
        graupelflag = .false.

      case default
        print *, 'warning:  cannot determine microphysics option!'
        print *, '          assuming all species present.'
        clwflag = .true.
        iceflag = .true.
        graupelflag = .true.

    end select

    return

  end subroutine get_wrf2_misc 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrf_1d(cdfid, varname, data1d, status)

     implicit none
     include 'netcdf.inc'
     integer, intent(in)                :: cdfid
     character(len=*),intent(in)        :: varname
     real, intent(inout)                :: data1d(:)
     integer, intent(out)               :: status
     
     integer                            :: vid
    
     status = nf_inq_varid(cdfid,trim(varname),vid)
     status = nf_get_var_real(cdfid,vid,data1d)
     if (status.ne.nf_noerr) then
       print *, 'netcdf error getting ', trim(varname),status
       status = 1
     else  
       status = 0
     endif
     return
  end subroutine get_wrf_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfnc_3d(cdfid, varname, output_stagger, nx,ny,nz,time, &
                          data3d, status)

    ! subroutine to return a 3d grid from wrf output.

    implicit none
    include 'netcdf.inc'
    
    integer, intent(in)                    :: cdfid
    character(len=*),intent(in)            :: varname
    character(len=1),intent(in)            :: output_stagger
    integer,intent(in)                     :: nx,ny,nz
    integer,intent(in)                     :: time
    real,intent(out)                       :: data3d(nx,ny,nz)
    integer, intent(out)                   :: status

    ! local vars
    integer                                :: vid,rcode,attid,ndims,natts
    integer                                :: dims(10),dimids(10)
    integer                                :: istart(10),iend(10)
    character(len=10)                      :: varstagger
    character(len=80)                      :: vname
    real, allocatable                      :: dum3d(:,:,:)
    integer                                :: i,nx_in,ny_in,nz_in,ivtype
    character(len=2)                       :: conv
    status = 0 
    rcode = nf_inq_varid(cdfid,trim(varname),vid)
    if (rcode .ne. nf_noerr) then
       print *, 'problem getting varid for ', trim(varname), ': ',rcode
       status = 1
       return
    endif
   
    ! get the dimensions of this variable
    rcode = nf_inq_var(cdfid,vid,vname,ivtype,ndims,dimids,natts)
    if (ndims .ne. 4) then
      print *, 'dimension problem in get_wrfnc_3d.'
      print *, 'data in file should be 4d, but has ', ndims, ' dims.'
      status = 1
      return
    endif

    do i = 1,ndims
      rcode = nf_inq_dimlen(cdfid,dimids(i),dims(i) )
    enddo
  
    ! check the dimensions 
    if (time .gt. dims(4) ) then
      print *, 'requested time index exceeds input times: ', time, dims(4)
      status = 1
      return
    endif
    nx_in = dims(1)
    ny_in = dims(2)
    nz_in = dims(3)
   
    istart(1) = 1
    iend(1)   = nx_in
    istart(2) = 1
    iend(2)   = ny_in
    istart(3) = 1
    iend(3)   = nz_in
    istart(4) = time
    iend(4)   = 1
       

    rcode = nf_get_att_text(cdfid,vid,'stagger',varstagger)
    if (rcode .ne. nf_noerr) then
      print *, 'problem getting stagger: ',rcode
      status = 1
      return
    endif
    if (varstagger(1:1).eq. "x") then
      varstagger = "u"
    elseif (varstagger(1:1).eq. "y") then
      varstagger = "v"
    else
      varstagger = "t"
    endif

    ! check input vs. output dimensions

    if (varstagger(1:1) .eq. output_stagger) then
       ! in/out dimensions must match
       if ( (nx_in .ne. nx) .or. (ny_in .ne. ny) .or. &
            (nz_in .ne. nz) ) then
          status = 1
       else 
         conv = "  "
       endif
    else
       select case (output_stagger)
       
         case ('a')
           if (varstagger(1:1) .eq. "t") then
             if ( (nx_in+1 .ne. nx) .or. (ny_in+1 .ne. ny) .or. &
                  (nz_in .ne. nz) ) then
                status = 1
             else 
               conv = 'ta'
             endif
           elseif(varstagger(1:1) .eq. "u") then
             if ( (nx_in .ne. nx) .or. (ny_in+1 .ne. ny) .or. &
                  (nz_in .ne. nz) ) then
                status = 1
             else  
               conv = "ua"
             endif
           elseif(varstagger(1:1) .eq. "v") then
             if ( (nx_in+1 .ne. nx) .or. (ny_in .ne. ny) .or. &
                  (nz_in .ne. nz) ) then
                status = 1
             else  
                conv = "va"
             endif
           else 
             print *, 'stagger conversion not supported.'
             status = 1
           endif
         case('t')
           if(varstagger(1:1) .eq. "u") then
             if ( (nx_in-1 .ne. nx) .or. (ny_in .ne. ny) .or. &
                  (nz_in .ne. nz) ) then
                status = 1
             else  
                conv = "ut"
             endif
           elseif(varstagger(1:1) .eq. "v") then
             if ( (nx_in .ne. nx) .or. (ny_in-1 .ne. ny) .or. &
                  (nz_in .ne. nz) ) then
                status = 1
             else 
                conv = "vt"
             endif
           else
             print *, 'stagger conversion not supported.'
           endif
         case default
           print *, 'requested stagger conversion not supported.'
           status = 1
       end select
    endif
    if (status .ne. 0) then
      print *, 'stagger_in/stagger_req: ',varstagger(1:1),output_stagger 
      print *, 'mismatch in dimensions: '
      print *, 'nx_in / nx_req : ', nx_in, nx
      print *, 'ny_in / ny_req : ', ny_in, ny
      print *, 'nz_in / nz_req : ', nz_in, nz
      return
    endif

    ! if no destaggering is required, we can just populate the output
    ! array directly
  
    if (conv .eq. "  " ) then
      call ncvgt( cdfid, vid, istart, iend, data3d, rcode)
      if (rcode .ne. nf_noerr) then
        print *, 'problem getting ', trim(varname)
        status = 1 
        return
      endif
    else
      allocate(dum3d(nx_in,ny_in,nz_in)) 
      call ncvgt( cdfid, vid, istart, iend, dum3d, rcode)
      if (rcode .ne. nf_noerr) then
        print *, 'problem getting ', trim(varname)
        status = 1
        return
      endif
      select case (conv)
        case ('ta') 
          call arakawa_c_t2n(dum3d,nx_in,ny_in,nz_in,data3d)
        case ('ua')
          call arakawa_c_u2n(dum3d,nx_in,ny_in,nz_in,data3d) 
        case ('va')
          call arakawa_c_v2n(dum3d,nx_in,ny_in,nz_in,data3d)
        case ('ut')
          call arakawa_c_u2t(dum3d,nx_in,ny_in,nz_in,data3d) 
        case ('vt')
          call arakawa_c_v2t(dum3d,nx_in,ny_in,nz_in,data3d)
        case default
          print *, 'no stagger case found.  should not be here.'
          status = 1
      end select
      deallocate(dum3d)
    endif 

    return 
  end subroutine get_wrfnc_3d  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfnc_2d(cdfid, varname, output_stagger, nx,ny,time, &
                          data2d, status)

    ! subroutine to return a 2d grid from wrf output.

    implicit none
    include 'netcdf.inc'

    integer, intent(in)                    :: cdfid
    character(len=*),intent(in)            :: varname
    character(len=1),intent(in)            :: output_stagger
    integer,intent(in)                     :: nx,ny
    integer,intent(in)                     :: time
    real,intent(out)                       :: data2d(nx,ny)
    integer, intent(out)                   :: status

    ! local vars
    integer                                :: vid,rcode,attid,ndims,natts
    integer                                :: dims(10),dimids(10)
    integer                                :: istart(10),iend(10)
    character(len=10)                      :: varstagger
    character(len=80)                      :: vname
    real, allocatable                      :: dum2d(:,:)
    integer                                :: i,nx_in,ny_in,ivtype
    character(len=2)                       :: conv
    status = 0
    rcode = nf_inq_varid(cdfid,trim(varname),vid)
    if (rcode .ne. nf_noerr) then
       print *, 'problem getting varid for ', trim(varname), ': ',rcode
       status = 1
       return
    endif

    ! get the dimensions of this variable
    rcode = nf_inq_var(cdfid,vid,vname,ivtype,ndims,dimids,natts)
    if (ndims .ne. 3) then
      print *, 'dimension problem in get_wrfnc_2d.'
      print *, 'data in file should be 3d, but has ', ndims, ' dims.'
      status = 1
      return
    endif

    do i = 1,ndims
      rcode = nf_inq_dimlen(cdfid,dimids(i),dims(i) )
    enddo

    ! check the dimensions
    if (time .gt. dims(3) ) then
      print *, 'requested time index exceeds input times: ', time, dims(3)
      status = 1
      return
    endif
    nx_in = dims(1)
    ny_in = dims(2)

    istart(1) = 1
    iend(1)   = nx_in
    istart(2) = 1
    iend(2)   = ny_in
    istart(3) = time
    iend(3)   = 1


    rcode = nf_get_att_text(cdfid,vid,'stagger',varstagger)
    if (rcode .ne. nf_noerr) then
      print *, 'problem getting stagger: ',rcode
      status = 1
      return
    endif
    if (varstagger(1:1).eq. "x") then
      varstagger = "u"
    elseif (varstagger(1:1).eq. "y") then
      varstagger = "v"
    else
      varstagger = "t"
    endif

    ! check input vs. output dimensions

    if (varstagger(1:1) .eq. output_stagger) then
       ! in/out dimensions must match
       if ( (nx_in .ne. nx) .or. (ny_in .ne. ny)) then
          status = 1
       else
         conv = "  "
       endif
    else
       select case (output_stagger)

         case ('a')
           if (varstagger(1:1) .eq. "t") then
             if ( (nx_in+1 .ne. nx) .or. (ny_in+1 .ne. ny))then
                status = 1
             else
               conv = 'ta'
             endif
           elseif(varstagger(1:1) .eq. "u") then
             if ( (nx_in .ne. nx) .or. (ny_in+1 .ne. ny)) then  
                status = 1
             else
               conv = "ua"
             endif
           elseif(varstagger(1:1) .eq. "v") then
             if ( (nx_in+1 .ne. nx) .or. (ny_in .ne. ny)) then  
                status = 1
             else
                conv = "va"
             endif
           else
             print *, 'stagger conversion not supported.'
             status = 1
           endif
         case('t')
           if(varstagger(1:1) .eq. "u") then
             if ( (nx_in-1 .ne. nx) .or. (ny_in .ne. ny) ) then 
                status = 1
             else
                conv = "ut"
             endif
           elseif(varstagger(1:1) .eq. "v") then
             if ( (nx_in .ne. nx) .or. (ny_in-1 .ne. ny)) then  
                status = 1
             else
                conv = "vt"
             endif
           else
             print *, 'stagger conversion not supported.'
           endif
         case default
           print *, 'requested stagger conversion not supported.'
           status = 1
       end select
    endif
    if (status .ne. 0) then
      print *, 'stagger_in/stagger_req: ',varstagger(1:1),output_stagger 
      print *, 'mismatch in dimensions: '
      print *, 'nx_in / nx_req : ', nx_in, nx
      print *, 'ny_in / ny_req : ', ny_in, ny
      return
    endif

    ! if no destaggering is required, we can just populate the output
    ! array directly

    if (conv .eq. "  " ) then
      call ncvgt( cdfid, vid, istart, iend, data2d, rcode)
      if (rcode .ne. nf_noerr) then
        print *, 'problem getting ', trim(varname)
        status = 1
        return
      endif
    else
      allocate(dum2d(nx_in,ny_in))
      call ncvgt( cdfid, vid, istart, iend, dum2d, rcode)
 print *, 'var:',varname,' pre-stagger min/max:',minval(dum2d),maxval(dum2d)
      if (rcode .ne. nf_noerr) then
        print *, 'problem getting ', trim(varname)
        status = 1
        return
      endif
      select case (conv)
        case ('ta')
          call arakawa_c_t2n(dum2d,nx_in,ny_in,1,data2d)
        case ('ua')
          call arakawa_c_u2n(dum2d,nx_in,ny_in,1,data2d)
        case ('va')
          call arakawa_c_v2n(dum2d,nx_in,ny_in,1,data2d)
        case ('ut')
          call arakawa_c_u2t(dum2d,nx_in,ny_in,1,data2d)
        case ('vt') 
          call arakawa_c_v2t(dum2d,nx_in,ny_in,1,data2d)
        case default
          print *, 'no stagger case found.  should not be here.'
          status = 1
      end select
      deallocate(dum2d) 
    if (varname .eq. 'gsw') then
      where(data2d .lt. 0) data2d = 0.
    endif
    endif

    return
  end subroutine get_wrfnc_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_wrf_file_name(prddir,domnum,time_min,fname)

! creates the wrf output file name based on the working directory,
! domain number, and number of minutes into simulation

    implicit none

    character(len=*), intent(in)               :: prddir
    integer, intent(in)                        :: domnum
    integer, intent(in)                        :: time_min
    character(len=255), intent(out)            :: fname

    character(len=2)                           :: domstr
    character(len=6)                           :: timestr

    write(domstr, '(i2.2)') domnum
    write(timestr, '(i6.6)') time_min

    fname = trim(prddir) // '/wrfout_d' // domstr // '_' // timestr
    return
  end subroutine make_wrf_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_wrf2_file_name(prddir,domnum,timestr,fname)

    ! creates the wrf output file name based on the working directory, 
    ! domain number, and number of minutes into simulation

    implicit none

    character(len=*), intent(in)               :: prddir
    integer, intent(in)                        :: domnum
    character(len=24), intent(in)              :: timestr
    character(len=255), intent(out)            :: fname

    character(len=2)                           :: domstr

    write(domstr, '(i2.2)') domnum

    fname = trim(prddir) // '/wrfout_d' // domstr // '_' // timestr(1:19)
    return
  end subroutine make_wrf2_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine wrfio_wait(filename,max_wait_sec)

    implicit none
    character(len=255)  :: filename
    character(len=8)    :: date_ready
    character(len=10)   :: time_ready
    logical             :: file_exists
    logical             :: file_ready
    integer             :: num_checks
    integer             :: max_wait_sec
    integer, parameter  :: pause_sec = 30
    integer             :: secs_waited
    integer             :: cdf, rcode, dimid, ntimes,status,vid,itimestep
    real                :: dt,xtime
    file_ready = .false.
    file_exists = .false.
    num_checks = 0
    print *, "checking status of ",trim(filename)
    do while (.not.file_ready)
      inquire(file=trim(filename), exist=file_exists)
      if (.not. file_exists) then
        print *, 'file not ready: ', trim(filename)
        print '(a,i3,a)', 'sleeping for ', pause_sec, ' seconds'
        call sleep(pause_sec)
        num_checks = num_checks + 1
        secs_waited = num_checks * pause_sec
        print '(a,i5,a)', 'total sleep time now ', secs_waited, ' seconds'
        if (secs_waited .ge. max_wait_sec) then
          print *, 'io_wait:  timeout waiting for file: ', trim(filename)
          print '(a,i5,a)', '    maximum wait time set to ', max_wait_sec, 's'
          stop 'io_wait'
        endif
      else 
        num_checks = 0
        ! make sure it has been populated with at least 1 time period    
        do while (.not. file_ready)
       
          call open_wrfnc(filename,cdf,status)  
          rcode = nf_inq_dimid(cdf, "time", dimid)
          rcode = nf_inq_dimlen(cdf, dimid, ntimes)
           
          print *,"rcode/ntimes =", rcode,ntimes
          if (ntimes.gt.0) then
            rcode = nf_inq_varid(cdf,"itimestep",vid)
            rcode = nf_get_var_int(cdf,vid,itimestep)
            
            !!!!! modified by wei-ting (20130307) :            !!!!!
            !!!!!     use xtime instead of itimestep for wrfv3 !!!!!
            if (rcode .lt. 0) then
              write(6,*)' warning: itimestep not found, looking for xtime in wrfio_wait'
              rcode = nf_get_att_real(cdf,0,"dt",dt)
              rcode = nf_inq_varid(cdf,"xtime",vid)
              rcode = nf_get_var_real(cdf,vid,xtime)
              itimestep = int(xtime*60/dt)
            endif
            !!!!! end of modifying !!!!!

            print *,"rcode/itimestep ", rcode,itimestep
            if ((rcode .eq. nf_noerr) .and. (itimestep .ge. 0))then
              call date_and_time(date_ready,time_ready)
              print *, trim(filename), ' ready at ', date_ready, '/',time_ready
              file_ready = .true.
              !call sleep(pause_sec)
            else 
              call sleep(pause_sec)
            endif
          else 
            call sleep(pause_sec)
            num_checks = num_checks + 1
            secs_waited = num_checks * pause_sec
            if (secs_waited .ge. max_wait_sec) then
              print *, 'io_wait:  timeout waiting for file: ', trim(filename)
              print '(a,i5,a)', '    maximum wait time set to ', max_wait_sec, 's'
              stop 'io_wait'
            endif
          endif
          call close_wrfnc(cdf)
        enddo 
      endif
    enddo 
    return
  end subroutine wrfio_wait

end module wrf_netcdf
