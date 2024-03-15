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



module mem_static

  type static_type
  integer                     :: nx             ! x grid points
  integer                     :: ny             ! y grid points
  integer                     :: nz             ! z grid points
  real                        :: swlat          ! sw corner lat
  real                        :: swlon          ! se corner lon
  real                        :: nelat          ! ne corner lat
  real                        :: nelon          ! ne corner lon
  real                        :: dx             ! x-direction grid spacing
  real                        :: dy             ! y-direction grid spacing
  real                        :: cenlon         ! center longitude
  real                        :: cenlat         ! center latitude
  real                        :: stdlon         ! orientation longitude (polelon)
  real                        :: stdlat1        ! standard lat 1
  real                        :: stdlat2        ! standard lat 2 (polelat)
  real, pointer               :: topo(:,:)      ! laps topographic height
  real, pointer               :: glat(:,:)      ! laps latitude array
  real, pointer               :: glon(:,:)      ! laps longitude array
  real, pointer               :: ldf(:,:)       ! land fraction
  real, pointer               :: pbl(:,:)       ! pbl height
  real, pointer               :: akk(:,:)       ! drag coefficient
  character (len=132)         :: grid_type      ! map projection type
  end type

  type(static_type) :: stgrid

  real, pointer, dimension(:,:) :: lat, lon, topo, ldf

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc_static_arrays (nx, ny)

implicit none
integer :: nx, ny

    allocate ( stgrid%glat (nx,ny) )
    allocate ( stgrid%glon (nx,ny) )
    allocate ( stgrid%topo (nx,ny) )
    allocate ( stgrid%ldf  (nx,ny) )
    allocate ( stgrid%pbl (nx,ny) )
    allocate ( stgrid%akk  (nx,ny) )

return
end subroutine alloc_static_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine point_static_arrays ()

implicit none

lat =>  stgrid%glat
lon =>  stgrid%glon
topo  =>  stgrid%topo
ldf  =>  stgrid%ldf

return
end subroutine point_static_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloc_static_arrays ()

implicit none

    deallocate (stgrid%glat)
    deallocate (stgrid%glon)
    deallocate (stgrid%topo)
    deallocate (stgrid%ldf )
    deallocate (stgrid%pbl)
    deallocate (stgrid%akk )

return
end subroutine dealloc_static_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fill_static      (imax, jmax                           &
                            ,n_flds, dx, dy, lov, latin1         &
                            ,latin2, origin, var, comment         &
                            ,data, model, grid_spacing,map_proj   &
                            ,la1,lo1,la2,lo2                      &
                            ,center_lat, center_lon,lli,llj,uri   &
                            ,urj,parent_id,ratio_2_parent         &
                            ,status)

implicit none

integer     ::  imax, jmax, n_flds,status  &
               ,lli,llj,uri,urj,parent_id,ratio_2_parent
real        ::  dx, dy, lov, latin1, latin2
character(len=*) :: origin,var(n_flds),comment(n_flds),model,map_proj
real        :: data(imax,jmax,n_flds)
real        :: grid_spacing,la1,lo1,la2,lo2,center_lat,center_lon

stgrid%nx               = imax
stgrid%ny               = jmax
stgrid%swlat            = la1
stgrid%swlon            = lo1
stgrid%nelat            = la2
stgrid%nelon            = lo2
stgrid%dx               = grid_spacing
stgrid%dy               = grid_spacing
stgrid%cenlon           = center_lon
stgrid%cenlat           = center_lat
stgrid%stdlon           = lov
stgrid%stdlat1          = latin1
stgrid%stdlat2          = latin2
stgrid%topo(:,:)        = data(:,:,3)
stgrid%glat(:,:)        = data(:,:,1)
stgrid%glon(:,:)        = data(:,:,2)
stgrid%ldf(:,:)         = data(:,:,4)
stgrid%grid_type        = map_proj

!print*,'-----------:',imax,jmax,la1,lo1,la2,lo2,grid_spacing,grid_spacing  &
!,center_lon,center_lat,lov,latin1,latin2,data(1,1,3),data(1,1,1),data(1,1,2) &
!,map_proj


!         var(1)    = 'lat'
!         var(2)    = 'lon'
!         var(3)    = 'avg'
!         var(4)    = 'ldf'
!         var(5)    = 'use'
!         var(6)    = 'alb'  !now used for max snow alb 2-20-03 js.
!         var(7)    = 'std'
!         var(8)    = 'sln'
!         var(9)    = 'slt'
!         var(10)   = 'stl'
!         var(11)   = 'sbl'
!         var(12)   = 'lnd'  ! land-water mask based on usgs landuse
!         i=12
!         do j=1,12
!            write(cat,'(i2.2)')j
!            if(cat(1:1).eq.' ')cat(1:1)='0'
!            if(cat(2:2).eq.' ')cat(2:2)='0'
!            var(i+j)= 'g'//cat   ! vegetation greenness fraction
!         enddo
!
!         var(25)='tmp'
!         i=25
!         do j=1,12
!            write(cat,'(i2.2)')j
!            if(cat(1:1).eq.' ')cat(1:1)='0'
!            if(cat(2:2).eq.' ')cat(2:2)='0'
!            var(i+j)= 'a'//cat   ! monthly albedo
!         enddo
!
!         var(ngrids)   = 'zin'
!
!         comment(1) = 'lat: from model by j. smart/ s. albers 2-03\0'
!         comment(2) = 'lon: from model by j. smart/ s. albers 2-03\0'
!         comment(3) = 'average terrain elevation (m) \0'
!         comment(4) = 'land fraction: derived from usgs land use \0'
!         comment(5) = 'land use dominant category (usgs 24 category) \0'
!         comment(6) = 'maximum snow albedo; defined over land only \0'
!         comment(7) = 'standard deviation of elevation data (m)\0'
!         comment(8) = 'mean longitudinal terrain slope (m/m)\0'
!         comment(9) = 'mean latitudinal terrain slope (m/m)\0'
!         comment(10)= 'top layer (0-30cm) dominant category soiltype\0'
!         comment(11)= 'bot layer (30-90cm) dominant category soiltype\0'
!         comment(12)= 'land-water mask (0=water; 1 otherwise) \0'

return
end subroutine fill_static

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module mem_static
