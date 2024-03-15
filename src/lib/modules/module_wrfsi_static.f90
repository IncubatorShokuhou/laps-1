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

module wrfsi_static

  use time_utils
  ! f90 module to interact with the wrfsi static file, which is found
  ! in moad_dataroot/static.  
  !
  implicit none

  ! declare an allocatable array to use for fields which may have a variable
  ! 3rd dimension (such as categorical land use fraction, soil category fraction, etc.) 
  ! that can be allocated and returned to the calling routine
 
   real, allocatable                 :: static_landusef(:,:,:)
   real, allocatable                 :: static_soilbot(:,:,:)
   real, allocatable                 :: static_soiltop(:,:,:)
  


  ! delcare some variable name ids.  these are variables that
  ! map the netcdf names based on the cdl.  the naming convention is
  ! as follows:
  !
  ! name_var_n:  non-staggered variable
  ! name_var_m:  variable on the mass points of the staggered wrf grid
  ! name_var_u:  variable on the same points as the u-wind 
  ! name_var_v:     "     "   "   "     "     "  "  v-wind 


  ! grid spacing stuff
  character(len=2), parameter  :: name_dx = 'dx'
  character(len=2), parameter  :: name_dy = 'dy'
  
  ! projection info stuff
  character(len=6), parameter  :: name_truelat1 = 'latin1'
  character(len=6), parameter  :: name_truelat2 = 'latin2'
  character(len=3), parameter  :: name_stdlon = 'lov'
  character(len=3), parameter  :: name_lat1 = 'la1'
  character(len=3), parameter  :: name_lon1 = 'lo1'
  character(len=9), parameter  :: name_type = 'grid_type'
  character(len=2), parameter  :: name_nx = 'nx'
  character(len=2), parameter  :: name_ny = 'ny'
  ! latitude and longitude arrays

  character(len=3), parameter  :: name_lat_n = 'lat'
  character(len=3), parameter  :: name_lon_n = 'lon'
  character(len=3), parameter  :: name_lat_t = 'lac'
  character(len=3), parameter  :: name_lon_t = 'loc'
  character(len=3), parameter  :: name_lat_u = 'lab'
  character(len=3), parameter  :: name_lon_u = 'lob'
  character(len=3), parameter  :: name_lat_v = 'laa'
  character(len=3), parameter  :: name_lon_v = 'loa'
!mp-bls
  character(len=3), parameter  :: name_lat_h= 'lah'
  character(len=3), parameter  :: name_lon_h = 'loh'
  character(len=3), parameter  :: name_lat_w = 'lav'
  character(len=3), parameter  :: name_lon_w = 'lov'
!mp-bls

  ! stuff related to topography
  character(len=3), parameter  :: name_ter_n = 'avg'
  character(len=3), parameter  :: name_ter_t = 'avc'
  character(len=3), parameter  :: name_tersd_n = 'std'
  character(len=3), parameter  :: name_terenv_n = 'env'
  character(len=3), parameter  :: name_tergradln_n = 'sln'
  character(len=3), parameter  :: name_tergradlt_n = 'slt'

  ! coriolis
  character(len=3), parameter  :: name_hcor_t = 'cph'
  character(len=3), parameter  :: name_vcor_t = 'cpv'

  ! map factors
  character(len=3), parameter  :: name_mapfac_n = 'mfl'
  character(len=3), parameter  :: name_mapfac_t = 'mfc'
  character(len=3), parameter  :: name_mapfac_u = 'mfb'
  character(len=3), parameter  :: name_mapfac_v = 'mfa'

  ! sin of alpha angles (angle between longitude and stdlon)
  character(len=3), parameter  :: name_sinalpha_t = 'spr'
  character(len=3), parameter  :: name_cosalpha_t = 'cpr'

  ! land use categories
  character(len=3), parameter  :: name_landuse_t = 'use'

  ! land mask field
  character(len=3), parameter  :: name_lwmask_t = 'lnd'
  ! albedo climatology
  character(len=3), parameter  :: name_albedo_n = 'alb'

  ! mean annual deep soil temperature on mass grid.
  character(len=3), parameter  :: name_amt_t = 'tmp'

  include "netcdf.inc"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine open_wrfsi_static(dataroot,nestid,cdfid)
  
    implicit none
    character(len=*), intent(in)   :: dataroot
    integer, intent(in)            :: nestid
    integer, intent(out)           :: cdfid
    character(len=255)            :: staticfile
    logical                       :: static_exists
    integer                       :: status
    character(len=2)              :: nestid_str

    write(nestid_str, '(i2.2)') nestid
    staticfile = trim(dataroot) // '/static/static.wrfsi.d' // nestid_str
    inquire(file=staticfile, exist=static_exists)
    if (static_exists) then
      status = nf_open(trim(staticfile),nf_nowrite,cdfid)
      if (status .ne. nf_noerr) then
        print '(a,i5)', 'netcdf error opening wrf static file: ',status
        stop 'open_wrfsi_static'
      end if 
    else
!mp-bls
!       search for rotlat version??
!      print '(a)', 'static file not found ', staticfile
!      print '(a)', 'look for nmm version'
      staticfile = trim(dataroot) // '/static/static.wrfsi.rotlat'
      inquire(file=staticfile, exist=static_exists)
      if (static_exists) then
        status = nf_open(trim(staticfile),nf_nowrite,cdfid)
        if (status .ne. nf_noerr) then
          print '(a,i5)', 'netcdf error opening wrf static file: ',status
          stop 'open_wrfsi_static'
        end if
      else

        print '(a)', 'rotlat static file not found, either: ', staticfile
        stop 'open_wrfsi_static'
      endif

    endif
    return
  end subroutine open_wrfsi_static      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_dims(dataroot, nestid, nx, ny)
  
    ! subroutine to return the horizontal dimensions of wrf static file
    ! contained in the input dataroot

    implicit none
    character(len=*), intent(in)  :: dataroot
    integer         , intent(in)  :: nestid
    integer         , intent(out) :: nx
    integer         , intent(out) :: ny

    integer                       :: cdfid,vid, status

    call open_wrfsi_static(dataroot,nestid,cdfid)
    status = nf_inq_dimid(cdfid, 'x', vid)
    status = nf_inq_dimlen(cdfid, vid, nx)
    status = nf_inq_dimid(cdfid, 'y', vid)
    status = nf_inq_dimlen(cdfid, vid, ny) 
    status = nf_close(cdfid)  
    return
  end subroutine get_wrfsi_static_dims     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_proj(dataroot,nestid,proj_type, lat1, lon1, dx, dy, &
                                   stdlon, truelat1, truelat2)

    ! returns basic projection information from the wrf static file found
    ! in dataroot

    implicit none
    character(len=*) , intent(in)      :: dataroot
    integer, intent(in)                :: nestid
    character(len=32), intent(out)     :: proj_type
    real             , intent(out)     :: lat1
    real             , intent(out)     :: lon1
    real             , intent(out)     :: dx
    real             , intent(out)     :: dy
    real             , intent(out)     :: stdlon
    real             , intent(out)     :: truelat1
    real             , intent(out)     :: truelat2
   
    integer                            :: cdfid, vid,status
    integer                            :: nx,ny
    real,allocatable                   :: lats(:,:),lons(:,:)
    character(len=132)                 :: grid_type
    
    call get_wrfsi_static_dims(dataroot,nestid,nx,ny)
    allocate(lats(nx,ny))
    allocate(lons(nx,ny))
    call open_wrfsi_static(dataroot,nestid,cdfid)   
    status = nf_inq_varid( cdfid , name_type, vid )   
    status = nf_get_var_text( cdfid , vid , grid_type )  
    if (grid_type(1:19) .eq. 'polar stereographic') then 
      proj_type = 'polar stereographic             ' 
    else if ( grid_type(1:24) .eq. 'secant lambert conformal' ) then
      proj_type = 'lambert conformal               '   
    else if ( grid_type(1:28) .eq. 'tangential lambert conformal') then
      proj_type = 'lambert conformal               '  
    else if (grid_type(1:8)  .eq. 'mercator'                 ) then
      proj_type = 'mercator                        '    
!mp-bls
    else if (grid_type(1:15)  .eq. 'rotated lat-lon'                 ) then
      write(6,*) 'setting proj_type to rotated latlon'
      proj_type = 'rotated latlon                  '    
!mp-bls
    else
      print '(a,a)', 'unrecognized projection:', proj_type
      stop 'get_wrfsi_static_proj'
    end if                                                       

    ! get sw corner lat/lon of projection
    call get_wrfsi_static_2d(dataroot,nestid,name_lat_n,lats)
    call get_wrfsi_static_2d(dataroot,nestid,name_lon_n,lons)
    lat1 = lats(1,1)
    lon1 = lons(1,1)
    if (lon1 .lt. -180.) lon1 = lon1 + 360.
    if (lon1 .gt. +180.) lon1 = lon1 - 360.
    print '(a,f10.2,a,f10.2)', 'wrf lat1 = ',lat1, &
        ' wrf lon1 = ', lon1
    ! get dx and dy, convert to meters from kilometers

    status =  nf_inq_varid( cdfid, name_dx, vid )
    status = nf_get_var_real(cdfid,vid,dx)
    status = nf_inq_varid( cdfid, name_dy, vid )
    status = nf_get_var_real(cdfid,vid,dy)
    print '(a,f10.2,a,f10.2)', 'wrf delta-x = ',dx, &
        ' wrf delta-y = ', dy  
    ! get standard longitude
    status = nf_inq_varid ( cdfid , name_stdlon, vid )
    status = nf_get_var_real(cdfid , vid , stdlon)    
    if (stdlon .lt. -180.) stdlon = stdlon + 360.
    if (stdlon .gt. 180.) stdlon = stdlon - 360.
    print '(a,f10.3)', 'wrf standard lon = ', stdlon
    ! get true latitudes
    status = nf_inq_varid(cdfid , name_truelat1 , vid)
    status = nf_get_var_real( cdfid , vid , truelat1)
    status = nf_inq_varid(cdfid , name_truelat2, vid)
    status = nf_get_var_real( cdfid , vid , truelat2 )
    print '(a,2f10.3)', 'wrf standard lats = ', truelat1, truelat2
    status = nf_close(cdfid)
    deallocate(lats)
    deallocate(lons)
    return
  end subroutine get_wrfsi_static_proj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_latlon(dataroot, nestid,stagger, lat, lon)

    ! subroutine to get lat/lon arrays for desired grid stagger

    implicit none
    character(len=*), intent(in)        :: dataroot
    integer , intent(in)                :: nestid
    character(len=1), intent(in)        :: stagger
    real                                :: lat(:,:)
    real                                :: lon(:,:)

    integer                             :: status, cdfid, vid
    character(len=3)                    :: varname_lat, varname_lon

    call open_wrfsi_static(dataroot, nestid,cdfid)    

    if ( (stagger .eq. 'n').or.(stagger .eq. 'n').or.(stagger.eq. ' '))then
      varname_lat = name_lat_n
      varname_lon = name_lon_n
    else if ( (stagger .eq. 't').or.(stagger.eq.'t'))then
      varname_lat = name_lat_t
      varname_lon = name_lon_t  
    else if ( (stagger .eq. 'u').or.(stagger.eq.'u'))then
      varname_lat = name_lat_u
      varname_lon = name_lon_u
    else if ( (stagger .eq. 'v').or.(stagger.eq.'v'))then 
      varname_lat = name_lat_v
      varname_lon = name_lon_v
!mp-bls
    else if ( (stagger .eq. 'h').or.(stagger.eq.'h'))then 
      varname_lat = name_lat_h
      varname_lon = name_lon_h
    else if ( (stagger .eq. 'w').or.(stagger.eq.'w'))then 
      varname_lat = name_lat_w
      varname_lon = name_lon_w
!mp-bls
    else
      print '(2a)', 'unrecongized stagger code: ', stagger
      stop 'get_wrfsi_static_latlon'
    endif
  
    status = nf_inq_varid(cdfid, varname_lat, vid)
    status = nf_get_var_real(cdfid,vid,lat)    
    status = nf_inq_varid(cdfid, varname_lon, vid)
    status = nf_get_var_real(cdfid,vid,lon) 
  
    status = nf_close(cdfid) 
    return

  end subroutine get_wrfsi_static_latlon 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_mapfac(dataroot, nestid,stagger, mapfac) 

    ! subroutine to get alpha arrays for desired grid stagger

    implicit none
    character(len=*), intent(in)        :: dataroot
    integer, intent(in)                 :: nestid
    character(len=1), intent(in)        :: stagger
    real, intent(out)                   :: mapfac(:,:)

    integer                             :: cdfid, vid, status
    character(len=3)                    :: varname_map

    call open_wrfsi_static(dataroot,nestid,cdfid)
    if ( (stagger .eq. 'n').or.(stagger .eq. 'n').or.(stagger.eq. ' '))then
      varname_map = name_mapfac_n
    else if ( (stagger .eq. 'u').or.(stagger.eq.'u'))then
      varname_map = name_mapfac_u
    else if ( (stagger .eq. 'v').or.(stagger.eq.'v'))then
      varname_map = name_mapfac_v
    else if ( (stagger .eq. 't').or.(stagger.eq.'t')) then
      varname_map = name_mapfac_t
    else
      print '(2a)', 'unrecongized stagger code: ', stagger
      stop 'get_wrfsi_static_mapfac'
    endif

    status = nf_inq_varid(cdfid, varname_map, vid)
    status = nf_get_var_real(cdfid,vid,mapfac)

    status = nf_close(cdfid)
    return

  end subroutine get_wrfsi_static_mapfac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_landmask(dataroot, nestid,landmask)
 
    ! subroutine to get the 2d land mask field from the wrfsi static
    ! file
   
    implicit none
    character(len=*), intent(in)        :: dataroot
    integer, intent(in)                 :: nestid
    real            , intent(out)       :: landmask(:,:)

    integer                             :: cdfid, vid, status
    integer                             :: lu_water
    integer                             :: lu_ice
    character(len=4)                    :: lu_source
    integer                             :: n_cat, dim3,nx,ny
    real, allocatable                   :: landuse(:,:,:)
    integer                             :: i,j

    call open_wrfsi_static(dataroot,nestid,cdfid)  
    status = nf_inq_varid(cdfid, name_lwmask_t, vid)
    print '(a,2i4)','landmask inq_varid status and vid = ' ,status,vid
    if (status .ne. nf_enotvar) then
      status = nf_get_var_real(cdfid,vid,landmask)
      landmask = float(nint(landmask))
      if ( (minval(landmask) .lt. 0.) .or. &
           (maxval(landmask) .gt. 1.) ) then
        print '(a,2f12.8)', 'landmask values lt 0 or gt 1 found!:', &
          minval(landmask), maxval(landmask)
          stop 'get_wrfsi_static_landmask'
      endif
      landmask = float(nint(landmask)) ! prevents precision problems
    else
      print '(a)','problem getting wrf landmask from static.'
      print '(a)', '(old version of static file?)'
      print '(a)', 'attempting to use land use...'
      call get_wrfsi_static_landuse_info(dataroot,nestid,lu_source,n_cat,dim3, &
                                         lu_water,lu_ice)
      call get_wrfsi_static_dims(dataroot,nestid, nx,ny)
      allocate (landuse(nx,ny,dim3))
      call get_wrfsi_static_landuse(dataroot, nestid,landuse)
      if (dim3 .eq. 1) then
        print '(a)', 'using dominant category version...'
        where(nint(landuse(:,:,1)) .ne. lu_water) landmask = 1.
        where(nint(landuse(:,:,1)) .eq. lu_water) landmask = 0.
      else if (dim3 .eq. n_cat) then
        print '(a)','using fractional category version with a 50% water threshold...'
        landmask = 0.
        do j = 1,ny
          do i=1,nx
            if (landuse(i,j,lu_water).lt. 0.5) landmask(i,j) = 1.
          enddo
        enddo
      endif
      deallocate(landuse)
    endif
    status = nf_close(cdfid)
    return
  end subroutine get_wrfsi_static_landmask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_landuse_info(dataroot, nestid, source, n_categories, dim3, &
                                           water_ind,ice_ind )
 
    ! subroutine to query what type of land use we have in the static file
    !  source is the source of data (only 'usgs' is supported for now)
    !  n_categories tells how many different categories
    !  dim3 tells the third dimension of the landuse data array:
    !     landuse(nx,ny,dim3), where dim3 = 1 if it is a dominant
    !      category value or dim3 = n_categories if it is 
    !      fraction by category
    !  water_ind is the index value for the water type
    !  ice_ind is the index value for the ice type
   
    implicit none
    character(len=*), intent(in)        :: dataroot
    integer, intent(in)                 :: nestid
    character(len=4), intent(out)       :: source
    integer         , intent(out)       :: n_categories
    integer         , intent(out)       :: dim3
    integer         , intent(out)       :: water_ind
    integer         , intent(out)       :: ice_ind

    integer                             :: cdfid, vid, status
    integer                             :: dimid(4)
   

    call open_wrfsi_static(dataroot,nestid,cdfid)
    ! get land use source (hard coded for usgs 24 cat)
    source = 'usgs'

    ! get number of categories
    n_categories = 24
    
    ! get value that represents water (hard coded for usgs 24 cat)
    water_ind = 16
    
    ! get value that represents ice (hard coded for usgs 24 cat)
    ice_ind = 24

    ! query for dimension info
    status = nf_inq_varid(cdfid,name_landuse_t,vid)
    status = nf_inq_vardimid(cdfid,vid,dimid)
    status = nf_inq_dimlen(cdfid,dimid(3),dim3)
    
    print '(a)', 'static file landuse info:'
    print '(2a)', 'source: ', source
    print '(a,i4)', 'number of categories: ',n_categories
    print '(a,i4)', 'number of levels in array: ', dim3
    print '(a,i4)', 'water value: ', water_ind
    print '(a,i4)', 'ice value: ', ice_ind
    status = nf_close(cdfid)
    return
  end subroutine get_wrfsi_static_landuse_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_landuse(dataroot,nestid,landuse)
    
    ! gets the landuse data
    implicit none
    character(len=*), intent(in)  :: dataroot
    integer, intent(in)           :: nestid
    real, intent(inout)             :: landuse(:,:,:)
 
    integer                             :: cdfid, vid, status
   
    call open_wrfsi_static(dataroot,nestid,cdfid)
    status = nf_inq_varid(cdfid,name_landuse_t,vid)
    status = nf_get_var_real(cdfid,vid,landuse)
    if (status .ne. nf_noerr) then
      print '(a)', 'problem getting land use data.'
    endif  
    status = nf_close(cdfid)
    return
  end subroutine get_wrfsi_static_landuse  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_landusef(dataroot, nestid,source,ncat, water, ice)
    
    ! reads the individual 2d categorical landuse fraction arrays from
    ! the static file and populates the 3d variable delcared at the top of this
    ! module (static_landusef)

    implicit none
    character(len=*), intent(in)         :: dataroot
    integer, intent(in)                  :: nestid
    character(len=4), intent(out)        :: source
    integer         , intent(out)        :: ncat
    integer         , intent(out)        :: water
    integer         , intent(out)        :: ice
    integer                              :: dom_nx,dom_ny
    integer                              :: cdfid,vid,status
    integer                              :: cat
    character(len=3)                     :: vname
 
 
    ! for now, we only support 16-category wmo/fao data set, so hard code the
    ! info
    source = 'usgs'
    ncat = 24
    water = 16
    ice = 24
    print '(a)', 'attempting to read categorical landuse fractions...'
    print '(2a)','source = ', source
    print '(a,i2)','number of categories = ', ncat
    print '(a,i2)','index used for water = ', water
    print '(a,i2)','index used for ice = ', ice
    call get_wrfsi_static_dims(dataroot, nestid,dom_nx, dom_ny)
    if (allocated(static_landusef)) deallocate(static_landusef)
    allocate(static_landusef(dom_nx,dom_ny,ncat))
    call open_wrfsi_static(dataroot,nestid,cdfid)      
    do cat = 1, ncat
      write(vname,'("u",i2.2)') cat 
      status = nf_inq_varid(cdfid, vname, vid) 
      status = nf_get_var_real(cdfid,vid,static_landusef(:,:,cat))
      if (status .ne. nf_noerr) then
        print '(a)', 'problem getting landuse category #',cat
      endif
    enddo
    status = nf_close(cdfid)
    return
  end subroutine get_wrfsi_static_landusef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  subroutine get_wrfsi_static_soil(dataroot, nestid,source,ncat, st_water)
    
    ! reads the 2-layer soil categorical fractions and populates the arrays
    ! static_soiltop and static_soilbot

    implicit none
    character(len=*), intent(in)         :: dataroot 
    integer, intent(in)                  :: nestid
    character(len=4), intent(out)        :: source
    integer         , intent(out)        :: ncat
    integer         , intent(out)        :: st_water

    integer                              :: dom_nx,dom_ny
    integer                              :: cdfid,vid,status
    integer                              :: cat
    character(len=3)                     :: vname
 
 
    ! for now, we only support 16-category wmo/fao data set, so hard code the
    ! info
    source = 'fao '
    ncat = 16
    st_water = 14
    print '(a)', 'attempting to read categorical soil type fractions...'
    print '(2a)','source = ', source
    print '(a,i2)','number of categories = ', ncat
    call get_wrfsi_static_dims(dataroot, nestid,dom_nx, dom_ny)
    if (allocated(static_soiltop)) deallocate(static_soiltop)
    if (allocated(static_soilbot)) deallocate(static_soilbot)
    allocate(static_soiltop(dom_nx,dom_ny,ncat))
    allocate(static_soilbot(dom_nx,dom_ny,ncat))
    call open_wrfsi_static(dataroot,nestid,cdfid)      
    do cat = 1, ncat
      write(vname,'("b",i2.2)') cat 
      status = nf_inq_varid(cdfid, vname, vid) 
      status = nf_get_var_real(cdfid,vid,static_soilbot(:,:,cat))
      if (status .ne. nf_noerr) then
        print '(a)', 'problem getting bottom soil category #',cat
      endif
      write(vname,'("t",i2.2)') cat 
      status = nf_inq_varid(cdfid, vname, vid) 
      status = nf_get_var_real(cdfid,vid,static_soiltop(:,:,cat))
      if (status .ne. nf_noerr) then
        print '(a)', 'problem getting top soil category #',cat
      endif
    enddo
    status = nf_close(cdfid)
    return
  end subroutine get_wrfsi_static_soil  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_monthly(dataroot,nestid, dtypechar, time, data)
   
    ! returns a time-interpolated (valid for time) 2d array for either
    ! greenness ("g") or albedo ("a"), based on user-supplied type character.
    ! the monthly values are 
    ! valid on the 15th day of each month.  this routine only interpolates to the
    ! nearest day and does not account for leap years, but this should not be any
    ! big deal.  

    implicit none
    character(len=*) , intent(in)               :: dataroot
    integer, intent(in)                         :: nestid
    character(len=19), intent(in)               :: time    ! yyyy-mm-dd_hh:mm:ss
    character(len=1) , intent(in)               :: dtypechar
    real             , intent(out)              :: data(:,:)
    
    integer                                     :: midmonth_day(12)
    integer                                     :: valid_day
    integer                                     :: yyyyjjj
    real                                        :: sss
    integer                                     :: m, d1, d2, m1, m2
    integer                                     :: nx,ny
    character(len=3)                            :: varname
    real, allocatable                           :: data1(:,:), data2(:,:)
    real                                        :: w1, w2

    ! midmonth_day is the julian day of the year corresponding to the 15th day
    ! of each month for a standard (non-leap) year
    data midmonth_day / 15, 43, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349 /

    ! check data type character to make sure it is either greenness or albedo.
    if ((dtypechar .ne. "a").and.(dtypechar.ne."g")) then
      print *, 'unknown data type character passed into get_wrfsi_static_monthly:', &
               dtypechar
      print *,'current supported values are a (albedo) and g (greenness fraction).'
      stop 'get_wrfsi_monthly'
    else if (dtypechar .eq. 'a') then
      print *, 'getting time-interpolated albedo.'
    else if (dtypechar .eq. 'g') then
      print *, 'getting time-interpolated greenness.'
    endif

    ! convert date string into integer yyyyjjj and sss
    call mm5_to_wrf_date(time, yyyyjjj, sss)
    valid_day = mod(yyyyjjj,1000)
    print *, 'time-interpolating to day # ', valid_day
    ! find bounding months
    if ((valid_day .lt. midmonth_day(1)) .or. (valid_day .gt. midmonth_day(12))) then
      ! december and january are bounding months
      d1 = midmonth_day(12)
      d2 = midmonth_day(1)
      m1 = 12
      m2 = 1
    else
      find_bounds: do m = 1, 11
        d1 = midmonth_day(m)
        d2 = midmonth_day(m+1)
        if (valid_day .eq. d1) then
           d2 = d1
           m1 = m
           m2 = m1
           exit find_bounds
        else if (valid_day .eq. d2) then
           d1 = d2
           m1 = m + 1
           m2 = m1
           exit find_bounds
        else if ((valid_day .gt. d1).and.(valid_day .lt. d2)) then
           m1 = m
           m2 = m + 1
           exit find_bounds
        endif
      enddo find_bounds
    endif

    ! if d1 = d2, then we don't need any interpolation, just get that month's 
    ! data values
    if ( d1 .eq. d2) then
      write(varname, '(a1,i2.2)') dtypechar,m1
      call get_wrfsi_static_2d(dataroot,nestid,varname,data)
    else
      ! we need to get the two months of bounding data and time interpolate
      call get_wrfsi_static_dims(dataroot,nestid, nx, ny)
      allocate(data1 (nx,ny))
      allocate(data2 (nx,ny))
      write(varname, '(a1,i2.2)') dtypechar,m1
      call get_wrfsi_static_2d(dataroot,nestid,varname,data1)
      write(varname, '(a1,i2.2)') dtypechar,m2
      call get_wrfsi_static_2d(dataroot,nestid,varname,data2)

      ! compute weights
      if (d2 .gt. d1) then
        w1 = ( float(d2) - float(valid_day) ) / float(d2-d1)
      else ! we must be between dec 15 and jan 15
        if (valid_day .lt. midmonth_day(1)) then ! we are in january
           w1 = ( float(d2) - float(valid_day) ) / 31.
        else ! we are in december
           w1 = ( 366. - float(valid_day) + float(midmonth_day(1)) ) / 31.
        endif
      endif
      w2 = 1. - w1
      data = w1*data1 + w2*data2
      deallocate(data1)
      deallocate(data2)
    endif
    return
  end subroutine get_wrfsi_static_monthly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_monthly_all(dataroot, nestid, dtypechar, data)
    implicit none
    character(len=*) , intent(in)               :: dataroot
    integer, intent(in)                         :: nestid
    character(len=1) , intent(in)               :: dtypechar
    real             , intent(out)              :: data(:,:,:)
    integer                                     :: m,nx,ny
    real, allocatable                           :: month_slab(:,:)
    character(len=3)                            :: varname

    ! allocate the dummy array to hold each months data
    call get_wrfsi_static_dims(dataroot,nestid, nx, ny)
    allocate (month_slab(nx,ny))
    month_slab(:,:) = 0.
    month_loop: do m = 1, 12  
      write (varname, '(a1,i2.2)') dtypechar, m
      call get_wrfsi_static_2d(dataroot,nestid,varname,month_slab)
      data(:,:,m) = month_slab(:,:)
      month_slab(:,:) = 0.
    enddo month_loop
    deallocate(month_slab)
    return
  end subroutine get_wrfsi_static_monthly_all
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_wrfsi_static_2d(dataroot, nestid,varname, data)

    ! gets any 2d variable from the static file
    character(len=*), intent(in)  :: dataroot
    integer, intent(in)           :: nestid
    character(len=*), intent(in)  :: varname
    real, intent(out)             :: data(:,:)
 
    integer                             :: cdfid, vid, status
   
    call open_wrfsi_static(dataroot,nestid,cdfid)
    status = nf_inq_varid(cdfid,varname,vid)
    status = nf_get_var_real(cdfid,vid,data)
    if (status .ne. nf_noerr) then
      print '(a)', 'problem getting 2d data.'
    endif 
    status = nf_close(cdfid) 
    return
  end subroutine get_wrfsi_static_2d    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module wrfsi_static
