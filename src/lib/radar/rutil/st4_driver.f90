  !==============================================================================================!
  ! subroutine st4_driver(infile,time,debug,err) by brad beechler 2010/03                        !
  !==============================================================================================!
  !                                                                                              !
  !      purpose  - reads in a netcdf file called <filename> that contains stage iv data         !
  !                                                                                              !
  !      inputs   - infile  : a text string with the input file's path and name                  !
  !                 time    : the laps i4time of the analysis (needed for header info)           ! 
  !                 debug   : the level of detail to output                                      !
  !                          0 - no output (except for errors)                                   !
  !                          1 - minimal output                                                  !
  !                          2 - detailed output                                                 !
  !                                                                                              !
  !      outputs  - err: a code that gives the exit status of the routine                        !
  !                     0 - successful (passed all internal tests)                               !
  !                     1 - couln't open the stage iv netcdf file (error from netcdf follows)    !
  !                     2 - couldn't define laps domain                                          !
  !                     3 - problem writing the output file                                      !
  !                                                                                              !
  !      requires - libraries  :                                                                 !
  !                     netcdf                                                                   !
  !                 subroutines:                                                                 !
  !                     cdf(status, err)       error handling                                    !
  !                     alloc_err(errorcode)   error handling                                    !
  !----------------------------------------------------------------------------------------------!
  subroutine st4_driver(infile,lapsi4t,debug,err)
    use netcdf
    implicit none
    real, parameter       :: st4_missing = 1e20                ! missing value for stage iv
    integer, intent(in)   :: debug                             ! input
    integer*4, intent(in) :: lapsi4t                           ! input
    integer, intent(out)  :: err                               ! output
    character*(*) :: infile                                    ! i/o strings
    integer       :: ncid, dimxid, dimyid                      ! netcdf ids
    integer       :: latvarid, lonvarid, datavarid             ! netcdf ids
    integer       :: xcount, ycount                            ! looping
    integer       :: num_good                                  ! counting
    integer       :: inxdim, inydim, outxdim, outydim, cdfout  ! holders
    real          :: laps_missing                              ! missing value for laps
    real          :: data_min, data_max, data_mean             ! holders
    real          :: la1, lo1, lov, in_dx, in_dy               ! holders
    real          :: lat_corners(4), lon_corners(4)            ! holders
    real, allocatable, dimension(:,:) :: indata, outdata       ! data arrays
    real, allocatable, dimension(:,:) :: inlat, outlat         ! data arrays
    real, allocatable, dimension(:,:) :: inlon, outlon         ! data arrays
    ! these are specific to laps routines
    character*150 :: directory
    character*125 :: comment
    character*31  :: extension
    character*10  :: units
    character*3   :: variable
    character*4   :: lvl_coord
    integer       :: lvl
    character*132 :: cmodel
    real, allocatable, dimension(:,:) :: grx, gry
    integer nxc,nyc,nzc
    real sw(2),ne(2),rota,lat0,lon0
    real xmin,ymin,dx,dy
    common /psgrid/nxc,nyc,nzc,lat0,lon0,rota,sw,ne
    common /pscorner/xmin,ymin,dx,dy

    lvl = 0
    err = 0
    if (debug.ne.0) write(6,'(a)') 'st4_driver v1.0 by brad beechler'
   
    ! this part is the netcdf side of the story where we open the stage iv
    ! file and read in all of the domain info and data
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    call cdf(nf90_open(infile, nf90_nowrite, ncid), err) ! open the file read-only.
    if (err.ne.0) return
    ! get the x and y dimensions (for some reason ncl has reversed these, no big deal)
    call cdf(nf90_inq_dimid(ncid, 'g5_x_0', dimyid), err) ! get the dimid of the data
    call cdf(nf90_inq_dimid(ncid, 'g5_y_1', dimxid), err)
    call cdf(nf90_inquire_dimension(ncid, dimxid, len = inxdim), err)
    call cdf(nf90_inquire_dimension(ncid, dimyid, len = inydim), err)
    if (err.ne.0) return ! a strange error occurred, maybe corrupt file

    allocate(indata(inxdim,inydim), stat=err) ; if (err.ne.0) call alloc_error(err)
    allocate(inlat(inxdim,inydim), stat=err) ; if (err.ne.0) call alloc_error(err)
    allocate(inlon(inxdim,inydim), stat=err) ; if (err.ne.0) call alloc_error(err)
    
    ! get the domain information
    call cdf(nf90_inq_varid(ncid, 'g5_lat_0', latvarid), err)
    call cdf(nf90_inq_varid(ncid, 'g5_lon_1', lonvarid), err)
    call cdf(nf90_get_att(ncid, latvarid, 'la1', la1), err)
    call cdf(nf90_get_att(ncid, latvarid, 'lo1', lo1), err)
    call cdf(nf90_get_att(ncid, latvarid, 'lov', lov), err)
    call cdf(nf90_get_att(ncid, latvarid, 'dx',  in_dx), err)
    call cdf(nf90_get_att(ncid, latvarid, 'dy',  in_dy), err)
    call cdf(nf90_get_att(ncid, latvarid, 'corners',  lat_corners), err)
    call cdf(nf90_get_att(ncid, lonvarid, 'corners',  lon_corners), err)
    if (err.ne.0) return ! a strange error occurred, maybe corrupt file
    
    ! get the data information this is different for different accumulation times so
    ! we will look for all three kinds, only one can be there
    cdfout = nf90_inq_varid(ncid, 'a_pcp_gds5_sfc_acc24h', datavarid)
    if (cdfout.ne.0) cdfout = nf90_inq_varid(ncid, 'a_pcp_gds5_sfc_acc6h',  datavarid)
    if (cdfout.ne.0) cdfout = nf90_inq_varid(ncid, 'a_pcp_gds5_sfc_acc1h',  datavarid)
    if (cdfout.ne.0) then ! something is wrong and we need to abort
      err = 1 ! set output err value to the netcdf problem code
      if (debug.ne.0) write(6,'(a)') 'aborting from netcdf: '//nf90_strerror(cdfout)
      return
    endif
    call cdf(nf90_get_var(ncid, datavarid, indata, start = (/ 1, 1 /) ), err)

    ! get latitudes and longitudes
    call cdf(nf90_inq_varid(ncid, 'g5_lat_0', datavarid), err)    
    call cdf(nf90_get_var(ncid, datavarid, inlat, start = (/ 1, 1 /) ), err)
    call cdf(nf90_inq_varid(ncid, 'g5_lon_1', datavarid), err)    
    call cdf(nf90_get_var(ncid, datavarid, inlon, start = (/ 1, 1 /) ), err)
    
    ! convert to meters and set missing values to laps values
    call get_r_missing_data(laps_missing,err)
    if (err.ne.1) then
       err = 2 ! set output err value to the laps problem code
       write(6,'(a)') 'error getting the laps missing data value'
       return
    endif
    data_min = st4_missing
    data_max = -st4_missing
    num_good = 0
    data_mean = 0
    do xcount = 1,inxdim ! loop over x and y
    do ycount = 1,inydim
      if (indata(xcount,ycount).ge.st4_missing) then
        indata(xcount,ycount) = laps_missing
      else if (debug.gt.1) then ! only need to do this for debug output
        indata(xcount,ycount) = indata(xcount,ycount) / 1000.0
        num_good = num_good + 1
        if (indata(xcount,ycount).gt.data_max) data_max = indata(xcount,ycount)
        if (indata(xcount,ycount).lt.data_min) data_min = indata(xcount,ycount)
        data_mean = ((data_mean * (num_good-1)) + indata(xcount,ycount)) / num_good
      endif
    enddo
    enddo ! loop over x and y
    
    call cdf(nf90_close(ncid), err) ! close the stage iv file
    
    if (debug.gt.1) then ! output grid information
      write(6,'(a)') ''
      write(6,'(a)') 'input information  '
      write(6,'(a)') '^^^^^^^^^^^^^^^^^^ '
      write(6,'(a,i4,a,i4)') 'dimensions : ',inxdim,' x ',inydim
      !write(6,'(a,f7.2,a,f7.2,a,f7.2)') 'la1/lo1 lov: ',la1,' / ', lo1, '   ', lov
      !write(6,'(a,f9.2,a,f9.2)') 'dx/dy      : ',in_dx,' / ', in_dy
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'lat corners: ',lat_corners(1),', ', &
                       lat_corners(2),', ',lat_corners(3),', ',lat_corners(4)
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'lon corners: ',lon_corners(1),', ', &
                       lon_corners(2),', ',lon_corners(3),', ',lon_corners(4)
      write(6,'(a,i7,a,i7,a,f6.2,a)') 'number of good values/total (%) : ', num_good, &
              ' / ',inxdim*inydim,'  (', (num_good*100.0)/(inxdim*inydim),'% )'
      write(6,'(a,f6.3,a,f6.3,a,f7.3)') 'data in mm (w/o missing) min/mean/max : ', &
                       1000.0*data_min,' / ',1000.0*data_mean,' / ',1000.0*data_max
    endif

    ! this part is the laps side of things that gets the output information and writes
    ! the data to the laps grid defined in $laps_data_root
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    call get_grid_dim_xy(outxdim,outydim,err)
    if (err.ne.1) then
      err = 2 ! set output err value to the laps problem code
      if (debug.ne.0) write(6,'(a)') 'from laps: error reading nest7grid.parms'
      return
    endif

    allocate(outdata(outxdim,outydim),stat=err) ; if (err.ne.0) call alloc_error(err)
    allocate(outlat(outxdim,outydim),stat=err) ; if (err.ne.0) call alloc_error(err)
    allocate(outlon(outxdim,outydim),stat=err) ; if (err.ne.0) call alloc_error(err)
    allocate(grx(outxdim,outydim),stat=err) ; if (err.ne.0) call alloc_error(err)
    allocate(gry(outxdim,outydim),stat=err) ; if (err.ne.0) call alloc_error(err)
    
    call read_static_grid(outxdim,outydim,'lat',outlat,err)
    if (err.ne.1) then
      err = 2
      if (debug.ne.0) write(6,'(a)') 'from laps: error reading latitude'
      return
    endif
    call read_static_grid(outxdim,outydim,'lon',outlon,err)
    if (err.ne.1) then
      err = 2
      if (debug.ne.0) write(6,'(a)') 'from laps: error reading longitude'
      return
    endif
    
    ! set global values for laps
    nxc   = inxdim
    nyc   = inydim
    nzc   = 1
    lat0  = 90.0
    lon0  = lov
    rota  = 0.0
    sw(1) = lat_corners(1) ! southwest
    sw(2) = lon_corners(1)
    ne(1) = lat_corners(3) ! northeast
    ne(2) = lon_corners(3)
    xmin  = 1
    ymin  = 1
    dx    = in_dx
    dy    = in_dy

    ! remap the data to the laps domain
    call init_hinterp(inxdim,inydim,outxdim,outydim,'ps',outlat,outlon,grx,gry,0,cmodel,.false.)
    call hinterp_field(inxdim,inydim,outxdim,outydim,1,grx,gry,indata,outdata,.false.)

    ! quality control negative values produced by the cubic interpolation 
    do xcount = 1,outxdim ! loop over x and y 
    do ycount = 1,outydim
      if (outdata(xcount,ycount).lt.0.0) outdata(xcount,ycount) = 0.0
    enddo
    enddo ! loop over x and y

    data_min = laps_missing
    data_max = -laps_missing
    num_good = 0
    data_mean = 0
    do xcount = 1,outxdim ! loop over x and y
    do ycount = 1,outydim
      if (outdata(xcount,ycount).lt.laps_missing) then
        num_good = num_good + 1
        if (outdata(xcount,ycount).gt.data_max) data_max = outdata(xcount,ycount)
        if (outdata(xcount,ycount).lt.data_min) data_min = outdata(xcount,ycount)
        data_mean = ((data_mean * (num_good-1)) + outdata(xcount,ycount)) / num_good
      endif
    enddo
    enddo ! loop over x and y
    
    if (debug.gt.1) then ! output grid information
      write(6,'(a)') ''
      write(6,'(a)') 'output information '
      write(6,'(a)') '^^^^^^^^^^^^^^^^^^ '
      write(6,'(a,i4,a,i4)') 'dimensions : ',outxdim,' x ',outydim
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'lat corners: ',outlat(1,1),', ', &
           outlat(outxdim,1),', ',outlat(1,outydim),', ',outlat(outxdim,outydim)
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'lat corners: ',outlon(1,1),', ', &
           outlon(outxdim,1),', ',outlon(1,outydim),', ',outlon(outxdim,outydim)
      write(6,'(a,i7,a,i7,a,f6.2,a)') 'number of good values/total (%) : ', num_good, &
              ' / ',outxdim*outydim,'  (', (num_good*100.0)/(outxdim*outydim),'% )'
      write(6,'(a,f6.3,a,f6.3,a,f7.3)') 'data in mm (w/o missing) min/mean/max : ', &
                       1000.0*data_min,' / ',1000.0*data_mean,' / ',1000.0*data_max
    endif

    if (data_mean.eq.0.0) write(6,'(a)') ''
    if (data_mean.eq.0.0) write(6,'(a)') '***warning***'
    if (data_mean.eq.0.0) write(6,'(a)') '   the output data mean is zero! check that your'
    if (data_mean.eq.0.0) write(6,'(a)') '   check that your laps_data_root is set up correctly'
    if (data_mean.eq.0.0) write(6,'(a)') '***warning***'
    if (data_mean.eq.0.0) write(6,'(a)') ''

    ! output the files using the laps standards
    directory = './'
    extension = 'st4'
    variable  = 'ppt'
    units     = 'mm'
    comment   = 'stage iv precipitation data (mm) remapped to laps domain'
    call write_laps_data(lapsi4t,directory,extension,outxdim,outydim,1,1, &
                         variable,lvl,lvl_coord,units,comment,outdata,err)
    err = 0
    return
  end subroutine st4_driver

  ! the following are used by the st4_driver

  !==============================================================================================!
  ! subroutine cdf(status, err)                                                                  !
  !==============================================================================================!
  !---   purpose - error handling framework for netcdf                                           !
  !----------------------------------------------------------------------------------------------!
  subroutine cdf(status, err)
    use netcdf
    integer, intent (in)  :: status
    integer, intent (out) :: err
      err = 0
    if (status.ne.nf90_noerr) then 
      err = 1
      write(6,'(a)') 'from netcdf: '//nf90_strerror(status)
      return
    endif
  end subroutine cdf  

  !==============================================================================================!
  !  subroutine | alloc_error(errcode)                                                           !
  !==============================================================================================!
  !---   purpose - displays information about a memory allocation error                          !
  !---   input - error code                                                                      !
  !----------------------------------------------------------------------------------------------!
  subroutine alloc_error(errcode)
    integer, intent(in) :: errcode
    write(6,'(a)') 'an error occurred while allocating memory.'
    if (errcode.eq.1) write(6,'(a)') 'error in system routine attempting to do allocation'
    if (errcode.eq.2) write(6,'(a)') 'an invalid data object has been specified for allocation'
    if (errcode.eq.3) write(6,'(a)') 'error in system routine attempting to do allocation'
    if (errcode.eq.4) write(6,'(a)') 'an invalid data object has been specified for allocation'
    stop
  end subroutine alloc_error
