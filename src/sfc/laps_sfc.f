cdis    forecast systems laboratory
cdis    noaa/oar/erl/fsl
cdis    325 broadway
cdis    boulder, co     80303
cdis 
cdis    forecast research division
cdis    local analysis and prediction branch
cdis    laps 
cdis 
cdis    this software and its documentation are in the public domain and 
cdis    are furnished "as is."  the united states government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  they assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  all modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  if significant modifications or enhancements 
cdis    are made to this software, the fsl software policy manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
c
c
	program laps_sfc

!       we are using just a subset of the global parameters in this driver
        use mem_namelist, only: read_namelist_laps
        use mem_namelist, only: nx_l,ny_l,nk_laps,maxstns,
     &                          laps_cycle_time,grid_spacing_m

        use mem_grid, only: lat,lon,topo,ldf

        use mem_sfcanl

        character*150 static_dir,filename
        character*9 a9time
c
!       read global parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/nest7grid.parms'
        call read_namelist_laps('lapsparms',filename)

!       read surface parameters into module memory structure
        filename = static_dir(1:len_dir)//'/surface_analysis.nl'
        call read_namelist_laps('sfc_anal',filename)

!       allocate static arrays (lat, lon, topo, ldf)
        allocate( lat(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate lat'
            stop
        endif

        allocate( lon(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate lon'
            stop
        endif

        allocate( topo(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate topo'
            stop
        endif

        allocate( ldf(nx_l,ny_l), stat=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' error: could not allocate ldf'
            stop
        endif

!       read static arrays
        call read_static_grid(nx_l,ny_l,'lat',lat,istatus)
        if(istatus .ne. 1)then
            stop
        endif

        call read_static_grid(nx_l,ny_l,'lon',lon,istatus)
        if(istatus .ne. 1)then
            stop
        endif

        call read_static_grid(nx_l,ny_l,'avg',topo,istatus)
        if(istatus .ne. 1)then
            stop
        endif

        call read_static_grid(nx_l,ny_l,'ldf',ldf,istatus)
        if(istatus .ne. 1)then
            stop
        endif

!       get system analysis time
        call get_systime(i4time,a9time,istatus)

!       note these are being allocated later on just before they are filled
!       call alloc_sfc_arrays(nx_l,ny_l)
!       call point_sfcanl_arrays(nx_l,ny_l)

	call laps_sfc_sub(i4time)

        call deallocate_sfcanl_arrays()
c
	end
c
