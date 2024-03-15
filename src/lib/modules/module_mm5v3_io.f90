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

module mm5v3_io

! f90 module to deal with reading data from mm5 v3 data files.

  private
  ! mm5v3 header fields

  type mm5_big_header
    integer               :: bhi (50,20)
    real                  :: bhr (20,20)
    character (len=80)    :: bhic(50,20)
    character (len=80)    :: bhrc(20,20)
  end type mm5_big_header

  type mm5_sub_header
    integer               :: ndim
    integer               :: start_index  (4)
    integer               :: end_index    (4)
    real                  :: xtime
    character(len=4)      :: staggering
    character(len=4)      :: ordering
    character(len=24)     :: current_date
    character(len=9)      :: name
    character(len=25)     :: units
    character(len=46)     :: description
  end type mm5_sub_header

  integer, parameter      :: bh_flag = 0
  integer, parameter      :: sh_flag = 1
  integer, parameter      :: eot_flag = 2
  character(len=4),parameter       :: crs_test = 'c   '
  character(len=4),parameter       :: dot_text = 'd   '
  character(len=4),parameter       :: press3d_text = 'yxp '
  character(len=4),parameter       :: sigma3d_text = 'yxs '
  character(len=4),parameter       :: wsigma_text = 'yxw '
  character(len=4),parameter       :: field2d_text = 'yx  '
  character(len=4),parameter       :: landuse_text = 'ca  '
  character(len=4),parameter       :: boundns_text = 'xsb '
  character(len=4),parameter       :: boundwe_text = 'ysb '
  character(len=4),parameter       :: boundns_w_text = 'xwb '
  character(len=4),parameter       :: boundwe_w_text = 'ywb '
  character(len=4),parameter       :: presslev_text  = 'p   '
  character(len=4),parameter       :: sigmalev_text  = 's   '
  integer, parameter               :: lamcon_flag = 1
  integer, parameter               :: polstr_flag = 2
  integer, parameter               :: merctr_flag = 3

  public open_mm5v3
  public get_mm5_map
  public get_mm5_time_info
  public get_mm5_misc
  public get_mm5_scalar
  public get_mm5_1d
  public get_mm5_2d
  public get_mm5_3d
  public make_data_file_name
  public make_terrain_file_name
  public io_wait

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine open_mm5v3 (fname, lun, status)
    
    ! opens an mm5v3 file

    implicit none
   
    ! arguments

    character(len=255),intent(in)          :: fname
    integer, intent(out)                   :: lun
    integer, intent(out)                   :: status

    ! local variables

    integer                                :: i
    logical                                :: used

    status = 0
    
    lun = -1
    ! find an available unit number
    find_lun: do i = 7,1023
      inquire (unit=i, opened=used)
      if(used) then
        cycle find_lun
      else
        lun = i
        exit find_lun
      endif
    enddo find_lun
    if (lun .lt. 0) then
      print '(a)', 'no available unit numbers!'
      status = 1
    else
      open ( file=trim(fname), unit=lun, access='sequential', &
             form='unformatted', status='unknown', iostat=status, &
             position='rewind')
    endif
  end subroutine open_mm5v3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_big_header(lun, bh)
    implicit none
    integer, intent(in)               :: lun
    type(mm5_big_header), intent(out) :: bh 
    logical                           :: opened
    integer                           :: hdr_flag
    inquire (unit=lun, opened=opened)
    if (.not.opened) then
      print *, 'get_mm5_big_header: file not opened yet!'
    else
      rewind (lun)
      read (lun) hdr_flag
      if (hdr_flag .eq. bh_flag) then
        read(lun) bh%bhi, bh%bhr, bh%bhic, bh%bhrc
      else
        print *, 'get_mm5_big_header: bh not found!'
      endif
    endif
    return
  end subroutine get_mm5_big_header
         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_map(lun, projection, proj_cent_lat, proj_cent_lon, &
                         truelat1, truelat2, cone_factor, pole_point, &
                         grid_spacing, nx, ny, &
                         status)

    implicit none
    
    ! arguments
    integer, intent(in)                    :: lun  ! logical unit number
    character (len=32), intent(out)        :: projection
    real, intent(out)                      :: proj_cent_lat
    real, intent(out)                      :: proj_cent_lon
    real, intent(out)                      :: truelat1
    real, intent(out)                      :: truelat2
    real, intent(out)                      :: cone_factor
    real, intent(out)                      :: pole_point
    real, intent(out)                      :: grid_spacing
    integer, intent(out)                   :: nx
    integer, intent(out)                   :: ny
    integer, intent(out)                   :: status
    ! local variables

    logical                                :: file_opened
    type (mm5_big_header)                  :: bh
    type (mm5_sub_header)                  :: sh
    integer                                :: header_flag

    status = 0
    ! inquire to see if file is already open. 
    inquire (unit=lun, opened=file_opened)
    if (.not.file_opened) then
      print '(a,i4)', 'file unit not opened: ', lun
      print '(a)', 'call open_mm5v3 first.'
      status=1
    endif
    if (status .eq. 0) then
      call get_mm5_big_header(lun,bh)
      ! set projection
      select case (bh%bhi(7,1))
        case (lamcon_flag)
          projection = 'lambert conformal               '
        case (polstr_flag) 
          projection = 'polar stereographic             '
        case (merctr_flag) 
          projection = 'mercator                        '
      end select
      ! get projection center lat/lon
      proj_cent_lat = bh%bhr(2,1)
      proj_cent_lon = bh%bhr(3,1)
      ! get cone factor and true latitudes
      truelat1 = bh%bhr(5,1)
      truelat2 = bh%bhr(6,1)
      cone_factor = bh%bhr(4,1)
      pole_point = bh%bhr(7,1)
      ! get grid spacing and dimensions
      grid_spacing = bh%bhr(9,1)
      nx = bh%bhi(17,1)
      ny = bh%bhi(16,1)

      ! print out diagnostics
      print '(a)', 'mm5 map projection parameters'
      print '(a)', '-----------------------------'
      print '(2a)','proj: ', trim(projection)
      print '(a,f10.3)', 'proj center lat: ', proj_cent_lat
      print '(a,f10.3)', 'proj center lon: ', proj_cent_lon
      print '(a,f10.3)', 'true latitude 1: ', truelat1
      print '(a,f10.3)', 'true latitude 2: ', truelat2
      print '(a,f6.3)', 'cone factor: ', cone_factor
      print '(a,f10.3)', 'pole point: ', pole_point
      print '(a,f10.1)', 'grid spacing: ', grid_spacing
      print '(a,i5)', 'x dimension: ', nx
      print '(a,i5)', 'y dimension: ', ny
    endif
  end subroutine get_mm5_map   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine find_mm5_var(lun, varname, valtime, sh, status)
 
    ! searches an mm5v3 file for a variable and returns its subheader

    implicit none
    
    ! arguments
    integer, intent(in)                 :: lun
    character(len=9), intent(in)        :: varname
    type(mm5_sub_header), intent(out)   :: sh
    character(len=24), intent(in)       :: valtime
    integer, intent(out)                 :: status

    ! local variables
    logical                             :: file_opened
    integer                             :: nx,ny,nz
    integer                             :: header_flag
    logical                             :: var_found
  
    status = 0
    var_found = .false.
    inquire(unit=lun, opened=file_opened)
    if (.not.file_opened) then 
      print *, 'find_mm5_var: file not open in search for: ', trim(varname)
      status = 1
    else
      rewind (lun)
      search_loop: do while (.not. var_found)
        read (lun, end=999, iostat=status) header_flag
        select case (header_flag)
          case(bh_flag) 
            read(lun)
          case(eot_flag) 
            cycle search_loop
          case(sh_flag)
            read (lun) sh%ndim, sh%start_index, sh%end_index, sh%xtime, &
                       sh%staggering, sh%ordering, sh%current_date, &
                       sh%name, sh%units, sh%description
            if (sh%name .eq. varname) then
              if (valtime(1:10) .lt. '1960-01-01') then
                var_found = .true.
              else
                if (valtime(1:16).eq.sh%current_date(1:16))then
                  var_found = .true.
                endif
              endif
            else
              read (lun)      ! skip the next data item
            endif
          case default
            print '(a,i4)', 'bad header flag:', header_flag
        end select
      enddo search_loop
 999  if (.not. var_found) then
        print *, 'find_mm5_var: unable to find ', trim(varname), &
                 ' at time ', valtime(1:16)
        status = 1
      endif
    endif
  end subroutine find_mm5_var
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_scalar(lun, varname, valtime, vardata, status)
  
    ! gets a scalar data value from an mm5 data file
 
    implicit none

    ! arguments
    integer, intent(in)                :: lun
    character(len=9), intent(in)       :: varname
    character(len=24),intent(in)       :: valtime
    real, intent(out)                :: vardata
    integer, intent(out)               :: status

    ! local variables
    type (mm5_sub_header)              :: sh
    logical                            :: file_opened

    status = 0
    inquire(unit=lun, opened=file_opened)
    if (.not. file_opened) then
      print '(a)', 'you need to call open_mm5v3 first!'
      status = 1
    else
      call find_mm5_var(lun, varname, valtime, sh, status)
      if (status .ne. 0) then 
        print '(3a,i4)', 'variable ', trim(varname), ' not found in unit ',lun
      else
        if (sh%ndim .eq. 0) then
          read (lun) vardata
        else
          print '(3a)', 'variable ', trim(varname),' is not a scalar!'
          print '(a,i4)', 'number of dimensions found: ', sh%ndim
          status = 1
        endif
      endif
    endif
  end subroutine get_mm5_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_1d(lun,varname,valtime,vardata,status)
    
    ! gets a 1d array variable from the mm5v3 file
   implicit none

    ! arguments
    integer, intent(in)                :: lun
    character(len=9), intent(in)       :: varname
    character(len=24),intent(in)       :: valtime
    real, intent(out)                  :: vardata(:)
    integer, intent(out)               :: status

    ! local variables
    type (mm5_sub_header)              :: sh
    logical                            :: file_opened

    status = 0
    inquire(unit=lun, opened=file_opened)
    if (.not. file_opened) then
      print '(a)', 'you need to call open_mm5v3 first!'
      status = 1
    else
      call find_mm5_var(lun, varname, valtime, sh, status)
      if (status .ne. 0) then 
        print '(3a,i4)', 'variable ', trim(varname), ' not found in unit ',lun
      else
        if (sh%ndim .eq. 1) then
          read (lun) vardata
        else
          print '(3a)', 'variable ', trim(varname),' is not 1-d!'
          print '(a,i4)', 'number of dimensions found: ', sh%ndim
          status = 1
        endif
      endif
    endif
  end subroutine get_mm5_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_2d(lun,varname,valtime,vardata,output_stagger,status)
    
    ! gets a 2d array variable from the mm5v3 file, reforms the array
    ! to be (nx,ny) (mm5 uses ny,nx) and optionally transforms the 
    ! data to the desired output grid (cross or dot points) if needed

    implicit none

    ! arguments
    integer, intent(in)                :: lun
    character(len=9), intent(in)       :: varname
    character(len=24),intent(in)       :: valtime
    real, intent(out)                  :: vardata( : , : )
    character(len=4), intent(in)       :: output_stagger
    integer, intent(out)               :: status

    ! local variables
    type (mm5_sub_header)              :: sh
    logical                            :: file_opened
    integer                            :: xdim, ydim, i, j
    real, allocatable                  :: temp_array ( : , : )

    status = 0
    inquire(unit=lun, opened=file_opened)
    if (.not. file_opened) then
      print '(a)', 'you need to call open_mm5v3 first!'
      status = 1
    else
      call find_mm5_var(lun, varname, valtime, sh, status)
      if (status .ne. 0) then 
        print '(3a,i4)', 'variable ', trim(varname), ' not found in unit ',lun
      else
        if (sh%ndim .eq. 2) then
          xdim = sh%end_index(2) - sh%start_index(2) + 1
          ydim = sh%end_index(1) - sh%start_index(1) + 1
          allocate(temp_array (ydim,xdim))
          read (lun) temp_array
          ! re-order the array to by (nx,ny)
          do j = 1,ydim
            do i = 1,xdim
              vardata(i,j) = temp_array(j,i)
            enddo
          enddo
          deallocate (temp_array)
          ! determine if we need to do crs2dot or dot2crs
          if ((sh%staggering .eq. 'c   ').and.&
              (output_stagger .eq. 'd   ') ) then
            call crs2dot(vardata,xdim,ydim)
          else if ((sh%staggering .eq. 'd   ').and.&
              (output_stagger .eq. 'c   ') ) then
            call dot2crs(vardata,xdim,ydim)
          endif
        else
          print '(3a)', 'variable ', trim(varname),' is not 2-d!'
          print '(a,i4)', 'number of dimensions found: ', sh%ndim
          status = 1
        endif
      endif
    endif
  end subroutine get_mm5_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_3d(lun,varname,valtime,vardata,output_stagger,status)
    
    ! gets a 3d array variable from the mm5v3 file, reforms the array
    ! to be (nx,ny,nz) (mm5 uses ny,nx) and optionally transforms the 
    ! data to the desired output grid (cross or dot points) if needed

    implicit none

    ! arguments
    integer, intent(in)                :: lun
    character(len=9), intent(in)       :: varname
    character(len=24),intent(in)       :: valtime
    real, intent(out)                  :: vardata( : , : , : )
    character(len=4), intent(in)       :: output_stagger
    integer, intent(out)               :: status

    ! local variables
    type (mm5_sub_header)              :: sh
    logical                            :: file_opened
    integer                            :: xdim, ydim, zdim, i, j, k
    real, allocatable                  :: temp_array ( : , : , : )

    status = 0
    inquire(unit=lun, opened=file_opened)
    if (.not. file_opened) then
      print '(a)', 'you need to call open_mm5v3 first!'
      status = 1
    else
      call find_mm5_var(lun, varname, valtime, sh, status)
      if (status .ne. 0) then 
        print '(3a,i4)', 'variable ', trim(varname), ' not found in unit ',lun
      else
        if (sh%ndim .eq. 3) then
          xdim = sh%end_index(2) - sh%start_index(2) + 1
          ydim = sh%end_index(1) - sh%start_index(1) + 1
          zdim = sh%end_index(3) - sh%start_index(3) + 1
          allocate(temp_array (ydim,xdim, zdim))
          read (lun) temp_array
          ! re-order the array to be (nx,ny,nz) and make it go from bottom of
          ! atmosphere to top of atmosphere
          do k = 1,zdim
            do j = 1,ydim
              do i = 1,xdim
                vardata(i,j,k) = temp_array(j,i,zdim-k+1)
              enddo
            enddo
          enddo
          deallocate(temp_array)

          ! determine if we need to do crs2dot or dot2crs
          if ((sh%staggering .eq. 'c   ').and.&
              (output_stagger .eq. 'd   ') ) then
            do k=1,zdim
              call crs2dot(vardata(:,:,k),xdim,ydim)
            enddo
          else if ((sh%staggering .eq. 'd   ').and.&
                   (output_stagger .eq. 'c   ') ) then
            do k=1,zdim
              call dot2crs(vardata(:,:,k),xdim,ydim)
            enddo
          endif
        else
          print '(3a)', 'variable ', trim(varname),' is not 3-d!'
          print '(a,i4)', 'number of dimensions found: ', sh%ndim
          status = 1
        endif
      endif
    endif
  end subroutine get_mm5_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_time_info(lun, sim_start_year, sim_start_month, &
                               sim_start_day, sim_start_hour, sim_start_min, &
                               sim_start_sec, sim_start_frac, sim_stop_min, &
                               dom_start_min, dom_stop_min, &
                               fdda_on, fdda_start_min, fdda_stop_min, &
                               tapfrq, buffrq, status)

  ! gets all of the time variables related to domain start/stop, simulation
  ! start/stop, output frequency, etc.
  
    implicit none
   
    ! arguments
  
    integer, intent(in)                :: lun
    integer, intent(out)               :: sim_start_year
    integer, intent(out)               :: sim_start_month
    integer, intent(out)               :: sim_start_day
    integer, intent(out)               :: sim_start_hour
    integer, intent(out)               :: sim_start_min
    integer, intent(out)               :: sim_start_sec
    integer, intent(out)               :: sim_start_frac
    real, intent(out)                  :: sim_stop_min
    real, intent(out)                  :: dom_start_min
    real, intent(out)                  :: dom_stop_min
    logical,intent(out)                :: fdda_on
    real, intent(out)                  :: fdda_start_min
    real, intent(out)                  :: fdda_stop_min
    real, intent(out)                  :: tapfrq
    real, intent(out)                  :: buffrq
    integer, intent(out)               :: status

    ! locals
    logical                            :: file_opened
    type(mm5_big_header)               :: bh

    inquire(unit=lun, opened=file_opened)
    if (.not.file_opened) then
      print *, 'get_mm5_time_info: data file not open!'
      status = 1
    else
      ! find the big header by rewinding to the beginning
      call get_mm5_big_header(lun,bh)
      sim_start_year = bh%bhi(5,11)
      sim_start_month = bh%bhi(6,11)
      sim_start_day = bh%bhi(7,11)
      sim_start_hour = bh%bhi(8,11)
      sim_start_min = bh%bhi(9,11)
      sim_start_sec = bh%bhi(10,11)
      sim_start_frac = bh%bhi(11,11)
      sim_stop_min = bh%bhr(1,12)
      tapfrq = bh%bhr(4,12)
      buffrq = bh%bhr(5,12)
      dom_start_min = bh%bhr(1,14)
      dom_stop_min = min(bh%bhr(2,14),sim_stop_min)
      print *, 'dom_start_min/dom_stop_min=',dom_start_min,dom_stop_min
      if ((bh%bhi(1,16).eq.1).or.(bh%bhi(6,16).eq.1).or. &
         (bh%bhi(14,16).eq.1)) then
         fdda_on = .true.
         fdda_start_min = bh%bhr(1,16)
         fdda_stop_min = bh%bhr(2,16)
      else
         fdda_on = .false.
         fdda_start_min = 0.
         fdda_stop_min = 0.
      endif
    endif
  end subroutine get_mm5_time_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine crs2dot(field,dim1,dim2)
  
    ! transforms grid on cross-points to grid on dot points
    
    implicit none
    integer :: dim1 , dim2
    real , dimension(dim1,dim2) :: field,dummy
    integer :: i , j 
   
    dummy(2:dim1-1,2:dim2-1)  = ( field(1:dim1-2,1:dim2-2) +  & 
                                  field(1:dim1-2,2:dim2-1) +  &  
                                  field(2:dim1-1,1:dim2-2) +  &
                                  field(2:dim1-1,2:dim2-1))*0.25
 
    dummy(2:dim1-1,1:dim2:dim2-1)=( field(1:dim1-2,1:dim2-1:dim2-2) + & 
                                    field(2:dim1-1,1:dim2-1:dim2-2)) * 0.5

    dummy(1:dim1:dim1-1,2:dim2-1)=(field(1:dim1-1:dim1-2,1:dim2-2) + & 
                                   field(1:dim1-1:dim1-2,2:dim2-1) )*0.5

    dummy(1:dim1:dim1-1,1:dim2:dim2-1) = field(1:dim1-1:dim1-2,1:dim2-1:dim2-2)

    field = dummy

  end subroutine crs2dot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dot2crs(vardot, dim1, dim2)

    ! transforms variable on dot points to cross points
    implicit none

    integer, intent(in)               :: dim1
    integer, intent(in)               :: dim2
    real, intent(inout)               :: vardot(dim1,dim2)
    real                              :: dummy(dim1,dim2)  
    integer                           :: i
    integer                           :: j
  
!-----------------------------------------------------------------------
!--   interpolate dot point values to cross points using 
!--   four-point interpolation.
!-----------------------------------------------------------------------

    do j = 1, dim2-1
      do i = 1, dim1-1 

        dummy(i,j) = 0.25 * ( vardot(i,j) + vardot(i+1,j) + &
                              vardot(i,j+1) + vardot(i+1,j+1) )

      enddo
    enddo
    dummy(dim1,:) = vardot(dim1,:)
    dummy(:,dim2) = vardot(:,dim2)
    vardot = dummy

  end subroutine dot2crs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_mm5_misc(lun, ksig, ptop, pmslbase, tmslbase, &
                          dtdlnpbase, tisobase, &
                          clwflag, iceflag,graupelflag)

    implicit none
    integer, intent(in)       :: lun
    integer, intent(out)      :: ksig
    real,    intent(out)      :: ptop
    real,    intent(out)      :: pmslbase
    real,    intent(out)      :: tmslbase
    real,    intent(out)      :: dtdlnpbase
    real,    intent(out)      :: tisobase
    logical, intent(out)      :: clwflag
    logical, intent(out)      :: iceflag
    logical, intent(out)      :: graupelflag
    logical                   :: opened
    type (mm5_big_header)     :: bh
    inquire(unit=lun, opened=opened)

    clwflag = .false.
    iceflag = .false.
    graupelflag = .false.
    if (.not.opened) then 
      print *, 'get_mm5_misc: file not opened.  call mm5v3_open first.'
    else
      call get_mm5_big_header(lun,bh)
      ksig = bh%bhi(12,11)
      ptop = bh%bhr(2,2)
      pmslbase = bh%bhr(2,5)
      tmslbase = bh%bhr(3,5)
      dtdlnpbase = bh%bhr(4,5)
      tisobase = bh%bhr(5,5)
      if (bh%bhi(20,11).eq.1) clwflag = .true.
      if (bh%bhi(18,11).eq.1) iceflag = .true.
      if (bh%bhi(19,11).eq.1) graupelflag = .true.
    endif
    return
  end subroutine get_mm5_misc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine make_data_file_name(mm5_data_root, domain_num_str, &
                                 split_output, time_period,file_name, &
                                 file_num3)
    implicit none
    character(len=255),intent(in)    :: mm5_data_root
    character(len=1),intent(in)      :: domain_num_str
    logical, intent(in)              :: split_output
    integer, intent(in)              :: time_period
    character(len=255),intent(out)   :: file_name
    logical, intent(in)              :: file_num3
    character(len=2)                 :: time_period_str
    character(len=3)                 :: time_period_str3

    file_name = trim(mm5_data_root)//'/mm5prd/raw/mmout_domain'//domain_num_str
    if (split_output) then
      if (.not.file_num3) then      
        write(time_period_str, '(i2.2)') time_period
        file_name = trim(file_name) // '_' // time_period_str
      else
        write(time_period_str3, '(i3.3)') time_period
        file_name = trim(file_name) // '_' // time_period_str3
      endif
    endif
    return

  end subroutine make_data_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine make_terrain_file_name(mm5_data_root, domain_num_str,file_name)
    implicit none
    character(len=255),intent(in)    :: mm5_data_root
    character(len=1),intent(in)      :: domain_num_str
    character(len=255),intent(out)   :: file_name

    file_name = trim(mm5_data_root)//'/static/terrain_domain'//domain_num_str
    return
  end subroutine make_terrain_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine io_wait(filename,max_wait_sec)

    implicit none
    character(len=255)  :: filename
    character(len=255)  :: flagfile
    character(len=8)    :: date_ready
    character(len=10)    :: time_ready
    logical             :: file_ready
    integer             :: num_checks
    integer             :: max_wait_sec
    integer, parameter  :: pause_sec = 30
    integer             :: secs_waited
    file_ready = .false.
    num_checks = 0
    flagfile = trim(filename) // '.done'
    do while (.not.file_ready)
      inquire(file=trim(flagfile), exist=file_ready)
      ! in case this file was just created, wait to 
      ! give the file a chance to be completely written.  also, this
      ! keeps us from banging on the disk unnecessarily.
      if (.not. file_ready) then
        print *, 'file not ready: ', trim(filename)
        print '(a,i3,a)', 'sleeping for ', pause_sec, ' seconds'
        call sleep(pause_sec)
        num_checks = num_checks + 1
        secs_waited = num_checks * pause_sec
        print '(a,i5,a)', 'total sleep time now ', secs_waited, ' seconds'
        if (secs_waited .ge. max_wait_sec) then
          print *, 'io_wait:  timeout waiting for file: ', trim(filename)
          print '(a,i5,a)', '    maximum wait time set to ', max_wait_sec, 's'
          stop 'io_timeout'
        endif
      else 
        call date_and_time(date_ready,time_ready) 
        print *, trim(filename), ' ready at ', date_ready, '/',time_ready
      endif
    enddo 
    return
  end subroutine io_wait
  end module mm5v3_io
