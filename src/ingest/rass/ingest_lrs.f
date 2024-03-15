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

!       1997 jul      ken dritz        added call to get_grid_dim_xy.
!       1997 jul      ken dritz        pass nx_l, ny_l to ingest_lrs.

        character*9 a9_time
        character*8 c8_project, c8_blpfmt, c8_blpfmt_in

        call get_systime(i4time,a9_time,istatus)
        if(istatus .ne. 1)go to 999

        call get_grid_dim_xy(nx_l,ny_l,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'error getting horizontal domain dimensions'
           go to 999
        endif

        write(6,*)
        write(6,*)' running wpdn (nimbus) rass ingest'
        call ingest_lrs(i4time,nx_l,ny_l,j_status)
        write(6,*)' return from wpdn (nimbus) rass ingest'

        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)goto999

        call get_c8_blpfmt(c8_blpfmt_in,istatus)
        if(istatus .ne. 1)goto999

        if(c8_project .eq. 'rsa')then
            write(6,*)
            write(6,*)' running rsa/ldad local rass ingest '
            call ingest_rsalrs(i4time,nx_l,ny_l,j_status)
            write(6,*)' return from rsa/ldad local rass ingest'
        elseif(c8_project   .eq. 'wfo' .or. 
     1         c8_blpfmt_in .eq. 'madis'    )then
            write(6,*)
            write(6,*)' running madis (wfo) multi-agency profile ingest'       
            call ingest_madis_map(i4time,nx_l,ny_l,'lrs',j_status)
            write(6,*)' return from madis (wfo) map ingest'
        else
            write(6,*)
            write(6,*)' running blp (nimbus) rass ingest'
            call ingest_blplrs(i4time,nx_l,ny_l,j_status)
            write(6,*)' return from blp (nimbus) rass ingest'
        endif

999     continue
        end

