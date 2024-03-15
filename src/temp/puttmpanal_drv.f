cdis   
cdis    open source license/disclaimer, forecast systems laboratory
cdis    noaa/oar/fsl, 325 broadway boulder, co 80305
cdis    
cdis    this software is distributed under the open source definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    in particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - all modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - if significant modifications or enhancements are made to this
cdis    software, the fsl software policy manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    this software and its documentation are in the public domain
cdis    and are furnished "as is."  the authors, the united states
cdis    government, its instrumentalities, officers, employees, and
cdis    agents make no warranty, express or implied, as to the usefulness
cdis    of the software and documentation for any purpose.  they assume
cdis    no responsibility (1) for the use of the software and
cdis    documentation; or (2) to provide technical support to users.
cdis   
cdis
cdis
cdis   
cdis

        use mem_namelist, only: read_namelist_laps
        use mem_namelist, only: nx_l,ny_l,nk_laps
 
        use mem_grid, only: lat,lon,topo

        use mem_temp

        character*9 a9time        
        character*150 static_dir,filename
c
!       read global parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/nest7grid.parms'
        call read_namelist_laps('lapsparms',filename)

!       read temp parameters into module memory structure
        filename = static_dir(1:len_dir)//'/temp.nl'
        call read_namelist_laps('temp_anal',filename)

!       allocate static arrays (lat, lon, topo)
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

!       read static arrays (lat, lon, topo)
        call read_static_grid(nx_l,ny_l,'lat',lat,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps lat'
            stop
        endif

        call read_static_grid(nx_l,ny_l,'lon',lon,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps lon'
            stop
        endif

        call read_static_grid(nx_l,ny_l,'avg',topo,istatus)
        if(istatus .ne. 1)then
            write(6,*)' error getting laps topo'
            stop
        endif

!       get system time
        call get_systime(i4time,a9time,istatus)
        if(istatus .ne. 1)go to 999
        write(6,*)' systime = ',a9time

        call alloc_temp_arrays(nx_l,ny_l,nk_laps)
        call point_temp_arrays()

        call laps_temp      (i4time)

        call deallocate_temp_arrays()

999     continue

        end

